module QuantumTrajectories

using Distributions: Uniform
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
import ..Pulses: duration, drive_name
import ..Rollouts
import ..Rollouts: rollout!
using ..Rollouts:
    UnitaryODEProblem,
    UnitaryOperatorODEProblem,
    KetODEProblem,
    KetOperatorODEProblem,
    DensityODEProblem
using ..Rollouts: unitary_fidelity, _sim_index
using ..EmbeddedOperators: AbstractPiccoloOperator, EmbeddedOperator, unembed
using ..Isomorphisms: operator_to_iso_vec, ket_to_iso, iso_to_ket, iso_vec_to_operator
using ..Isomorphisms: density_to_compact_iso, compact_iso_to_density

import NamedTrajectories: NamedTrajectory, get_times

import ..QuantumSystems: default_algorithm

"""
    default_algorithm(sys::QuantumSystem)
    default_algorithm(sys::CompositeQuantumSystem)
    default_algorithm(sys::OpenQuantumSystem)

Return the default ODE algorithm for trajectory rollouts.
Uses `Tsit5()` for non-Hermitian systems (where Magnus adaptive error
control fails), `MagnusAdapt4()` for Hermitian systems.

For `CompositeQuantumSystem`, Hermiticity is derived from the subsystems —
Hermitian iff every subsystem is Hermitian. The coupling drift and drives are
assumed Hermitian (the typical closed-system case).

For `OpenQuantumSystem`, the Lindblad generator is non-Hermitian by
construction, so `Tsit5()` is the only sensible default — and it
matches the explicit `algorithm = Tsit5()` default in the
`DensityTrajectory` / `MultiDensityTrajectory` constructors.
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

# Abstract type and common interface
include("abstract_trajectory.jl")

# Concrete trajectory types
include("unitary_trajectory.jl")
include("ket_trajectory.jl")
include("ensemble_trajectory.jl")
include("density_trajectory.jl")
include("multi_density_trajectory.jl")
include("sampling_trajectory.jl")

# Interface methods (getters, accessors, fidelity)
include("interface.jl")

# Extract pulse from optimized trajectories
include("extract_pulse.jl")

# NamedTrajectory conversion
include("named_trajectory_conversion.jl")

end  # module
