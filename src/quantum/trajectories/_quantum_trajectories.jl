module QuantumTrajectories

using Distributions: Uniform
using Random: AbstractRNG, default_rng
using LinearAlgebra
using OrdinaryDiffEqLinear
using OrdinaryDiffEqTsit5
using TestItems
using ForwardDiff

using ..QuantumSystems:
    AbstractQuantumSystem, QuantumSystem, CompositeQuantumSystem, OpenQuantumSystem
using ..Pulses:
    AbstractPulse,
    AbstractSplinePulse,
    ZeroOrderPulse,
    LinearSplinePulse,
    CubicSplinePulse,
    GaussianPulse,
    ErfPulse,
    CompositePulse
using ..Pulses: duration, drive_name, n_drives, sample
using ..Pulses: get_knot_times, get_knot_count, get_knot_values, get_knot_derivatives
import ..Pulses: duration, drive_name, interpolation_basis
import ..Rollouts
import ..Rollouts: rollout!
using ..Rollouts:
    UnitaryODEProblem,
    UnitaryOperatorODEProblem,
    KetODEProblem,
    KetOperatorODEProblem,
    DensityODEProblem
using ..Rollouts: unitary_fidelity, _sim_index
using ..Rollouts: PiecewiseExponential, piecewise_exponential_solve, _time_independent
import SciMLBase
using ..EmbeddedOperators: AbstractPiccoloOperator, EmbeddedOperator, unembed
using ..Isomorphisms: operator_to_iso_vec, ket_to_iso, iso_to_ket, iso_vec_to_operator
using ..Isomorphisms: density_to_compact_iso, compact_iso_to_density

import NamedTrajectories: NamedTrajectory, get_times

import ..QuantumSystems: default_algorithm

"""
    default_algorithm(sys::QuantumSystem)
    default_algorithm(sys::CompositeQuantumSystem)
    default_algorithm(sys::OpenQuantumSystem)

Return the system's default ODE algorithm for trajectory rollouts.
Uses `Tsit5()` for non-Hermitian systems (where Magnus adaptive error
control fails), `MagnusAdapt4()` for Hermitian systems.

For `CompositeQuantumSystem`, Hermiticity is derived from the subsystems —
Hermitian iff every subsystem is Hermitian. The coupling drift and drives are
assumed Hermitian (the typical closed-system case).

For `OpenQuantumSystem`, the Lindblad generator is non-Hermitian by
construction, so `Tsit5()` is the only sensible default.

On a system with no explicit time dependence, a rollout with `algorithm = nothing` does
not reach this for spline or ZOH pulses (#379): a `ZeroOrderPulse` is propagated exactly
by [`PiecewiseExponential`](@ref Piccolo.Quantum.Rollouts.PiecewiseExponential), and a
`LinearSplinePulse` / `CubicSplinePulse` integrates with `Tsit5()`, whose error tracks
`abstol` / `reltol`. `MagnusAdapt4`'s tolerance is inert there: OrdinaryDiffEqLinear's
embedded error estimate for it reduces, for a state-independent generator, to a
commutator with the generator's second difference over the step, which vanishes
identically wherever `H(t)` is affine in `t` (every linear-spline and ZOH interval), so
the step size is set by the `tstops` alone (measured: a `1.6e-4` error floor on a
9-level, 200 ns linear-spline fixture at every tolerance from `1e-6` to `1e-12`). See
`docs/experiments/379-rollout-accuracy.md`.
"""
function default_algorithm(sys::QuantumSystem)
    return sys.hermitian ? MagnusAdapt4() : Tsit5()
end

function default_algorithm(sys::CompositeQuantumSystem)
    return all(s -> s.hermitian, sys.subsystems) ? MagnusAdapt4() : Tsit5()
end

function default_algorithm(::OpenQuantumSystem)
    return Tsit5()
end

# ---------------------------------------------------------------------------- #
# Rollout dispatch: exact piecewise exponential for ZOH, ODE otherwise (#379)
# ---------------------------------------------------------------------------- #

# A zero-order-hold pulse on a time-independent system has a piecewise-constant
# generator, so its exact propagator is a product of matrix exponentials.
_piecewise_exponential_eligible(system, pulse) =
    pulse isa ZeroOrderPulse && _time_independent(system)

function _use_piecewise_exponential(algorithm, system, pulse)
    eligible = _piecewise_exponential_eligible(system, pulse)
    algorithm isa PiecewiseExponential &&
        !eligible &&
        throw(
            ArgumentError(
                "PiecewiseExponential() needs a ZeroOrderPulse on a system with no " *
                "explicit time dependence; got a $(nameof(typeof(pulse))). Pass an ODE " *
                "algorithm (e.g. Tsit5()) instead.",
            ),
        )
    return algorithm isa PiecewiseExponential || (isnothing(algorithm) && eligible)
end

# The ODE algorithm for `algorithm = nothing`. A spline pulse on a time-independent
# system has a piecewise-polynomial generator, on which `MagnusAdapt4`'s step-size
# control is inert (see `default_algorithm`); `Tsit5` tracks the tolerance there.
_default_rollout_algorithm(system, pulse) =
    pulse isa Union{LinearSplinePulse,CubicSplinePulse} && _time_independent(system) ?
    Tsit5() : default_algorithm(system)

# Solve one rollout ODE problem, or propagate it exactly when the pulse is ZOH and no
# algorithm is requested. `algorithm = nothing` otherwise means `_default_rollout_algorithm`.
function _solve_rollout(
    prob,
    system,
    pulse,
    save_times,
    algorithm;
    abstol,
    reltol,
    kwargs...,
)
    _use_piecewise_exponential(algorithm, system, pulse) &&
        return piecewise_exponential_solve(
            prob,
            system,
            pulse,
            get_knot_times(pulse),
            save_times,
        )
    alg = isnothing(algorithm) ? _default_rollout_algorithm(system, pulse) : algorithm
    return solve(prob, alg; saveat = save_times, abstol, reltol, kwargs...)
end

# The ensemble counterpart of `_solve_rollout`, one member per initial state.
function _solve_ensemble_rollout(
    base_prob,
    initials,
    system,
    pulse,
    save_times,
    algorithm;
    abstol,
    reltol,
    kwargs...,
)
    _use_piecewise_exponential(algorithm, system, pulse) &&
        return piecewise_exponential_solve(
            base_prob,
            system,
            pulse,
            get_knot_times(pulse),
            save_times;
            initials,
        )
    prob_func(prob, i_or_ctx, _repeat = nothing) =
        remake(prob, u0 = initials[_sim_index(i_or_ctx)])
    alg = isnothing(algorithm) ? _default_rollout_algorithm(system, pulse) : algorithm
    return solve(
        EnsembleProblem(base_prob; prob_func = prob_func),
        alg;
        trajectories = length(initials),
        saveat = save_times,
        abstol,
        reltol,
        kwargs...,
    )
end

# The algorithm a stored solution was produced with. In-place `rollout!` with
# `algorithm = nothing` reuses it: the trajectory's `solution` field is typed by its
# first solve, so a different algorithm could not be stored back into it.
_stored_algorithm(sol::SciMLBase.AbstractODESolution) = sol.alg
_stored_algorithm(sol::SciMLBase.EnsembleSolution) =
    isempty(sol.u) ? nothing : _stored_algorithm(first(sol.u))
_stored_algorithm(_) = nothing

@testitem "default_algorithm — CompositeQuantumSystem dispatches via subsystems" begin
    using OrdinaryDiffEqLinear: MagnusAdapt4
    using OrdinaryDiffEqTsit5: Tsit5
    using Piccolo.Quantum.QuantumSystems: default_algorithm

    sys_h = QuantumSystem(PAULIS[:Z], [PAULIS[:X]], [(-1.0, 1.0)])
    sys_nh = QuantumSystem(PAULIS[:Z], [PAULIS[:X]], [(-1.0, 1.0)]; hermitian = false)

    csys_h = CompositeQuantumSystem(
        zeros(ComplexF64, 4, 4),
        Matrix{ComplexF64}[],
        [sys_h, sys_h],
        Float64[],
    )
    csys_nh = CompositeQuantumSystem(
        zeros(ComplexF64, 4, 4),
        Matrix{ComplexF64}[],
        [sys_h, sys_nh],
        Float64[],
    )

    @test default_algorithm(csys_h) isa MagnusAdapt4
    @test default_algorithm(csys_nh) isa Tsit5
end

@testitem "default_algorithm — OpenQuantumSystem returns Tsit5" begin
    using OrdinaryDiffEqTsit5: Tsit5
    using Piccolo.Quantum.QuantumSystems: default_algorithm

    L = ComplexF64[0.1 0.0; 0.0 0.0]  # decay operator
    sys = OpenQuantumSystem(PAULIS[:Z], [PAULIS[:X]], [1.0]; dissipation_operators = [L])
    @test default_algorithm(sys) isa Tsit5
end

@testitem "rollout(::DensityTrajectory) does not require explicit algorithm" begin
    using LinearAlgebra

    L = ComplexF64[0.1 0.0; 0.0 0.0]
    sys = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρg = ComplexF64[0.0 0.0; 0.0 1.0]
    qtraj = DensityTrajectory(sys, ρ0, ρg, 1.0)

    # Smoke-test: rollout! / rollout should pick a default algorithm via
    # default_algorithm(::OpenQuantumSystem) without the caller providing one.
    qtraj_new = rollout(qtraj)
    @test qtraj_new isa DensityTrajectory
end

# ---------------------------------------------------------------------------- #
# Module layout convention
#
# `AbstractQuantumTrajectory` behaviors (accessors, name helpers, display)
# live in `trajectory_interface.jl`. User-facing rollouts, fidelity, and
# system swaps live in `rollouts_extensions.jl` so they are easy to locate.
# Concrete trajectory types each have their own file.
# ---------------------------------------------------------------------------- #

# Abstract type
include("abstract_trajectory.jl")

# Concrete trajectory types
include("unitary_trajectory.jl")
include("ket_trajectory.jl")
include("ensemble_trajectory.jl")
include("density_trajectory.jl")
include("multi_density_trajectory.jl")
include("sampling_trajectory.jl")

# AbstractQuantumTrajectory interface (getters, name accessors, display)
include("trajectory_interface.jl")

# Rollouts / fidelity — user-facing solution-testing entry points
include("rollouts_extensions.jl")

# Extract pulse from optimized trajectories
include("extract_pulse.jl")

# NamedTrajectory conversion
include("named_trajectory_conversion.jl")

end  # module
