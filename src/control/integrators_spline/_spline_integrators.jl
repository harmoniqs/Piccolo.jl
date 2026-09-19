module SplineIntegrators

# ============================================================================
# SplineIntegrator — shared struct + DENSE cells (Piccolo home).
#
# Open-core slice 3b (harmoniqs/Piccolissimo.jl#430). This module owns the
# `SplineIntegrator` struct, the spline type trait, the dense algorithm types
# and data (Tsit5 / MagnusGL4 / MagnusAdapt4 / Rodas5P), the dense cells for
# every trajectory kind (residual call operators, analytic sensitivity ODE
# builders, Jacobian/Hessian assembly), and the shared interval-coefficient
# kernel (`spline_interval_coeffs.jl`).
#
# The MATRIX-FREE cells/kernels, the MF layout machinery, the GPU variant and
# the Chebyshev/Magnus MF cell data construction stay in Piccolissimo as
# method extensions on the shared types declared here (the 3a extension-seam
# pattern; see the `matrix_free_step_scratch` / `_build_magnus_ket_data` /
# `_chebyshev_forward!` seams below and harmoniqs/Piccolo.jl#314's
# `matrix_free_jacobian_op` precedent).
# ============================================================================ #

export AbstractSplineIntegrator
export SplineIntegrator
export du_name, ddu_name
export build_sensitivity_ode, build_sensitivity_problems, extract_sensitivity_solution!
export PropagatorResult, get_propagator, get_sensitivities, get_sensitivities_flat
export IntegrationAlgorithm, Tsit5Alg, MagnusGL4Alg, MagnusAdapt4Alg
export Tsit5Data, MagnusGL4Data, MagnusAdapt4Data
export ChebyshevAlg
export refresh_sensitivities!
export SplineType, LinearSpline, CubicSpline
export SplineIntervalCoeffs, interval_coeff!, interval_coeff_dir!, interval_vjp_scatter!
export interval_hvp_scatter!
export LindbladDuhamelTape, DensityLindbladData, compact_iso_hs_weights
# Knot-point sensitivity / VJP-JVP surface (moved with the dense cells; the
# top-level KnotPointSensitivity module in Piccolissimo re-exports these).
export gauss_legendre_01,
    compute_sensitivity_kick_first_order,
    compute_sensitivity_kick_exact,
    KnotPointPropagationData,
    setup_knot_point_propagation,
    ket_vjp,
    ket_jvp
export unitary_rollout_trajectory

using OrdinaryDiffEqTsit5
# #342: the SciML integrator interface, for the preallocated per-interval HVP
# solvers. `init` / `solve!` / `reinit!` are exported by a dozen loaded packages
# (and `solve!` by DirectTrajOpt too), so bring them in under explicit names
# rather than leaving the bare bindings ambiguous.
using OrdinaryDiffEqTsit5: init as ode_init, solve! as ode_solve!, reinit! as ode_reinit!
using SciMLBase: ODEProblem, solve, remake
using OrdinaryDiffEqLinear: MagnusGL4, MagnusAdapt4
using LinearAlgebra
using SparseArrays
using NamedTrajectories
using NamedTrajectories: KnotPoint
using DirectTrajOpt: AbstractIntegrator
using DirectTrajOpt.Integrators: AbstractBilinearIntegrator, test_integrator
import DirectTrajOpt.Integrators:
    get_jacobian_structure, get_hessian_of_lagrangian_structure
using DirectTrajOpt.CommonInterface
import DirectTrajOpt.CommonInterface:
    evaluate!,
    jacobian_structure,
    jacobian!,
    hessian_structure,
    hessian_of_lagrangian,
    hessian_of_lagrangian!
import DirectTrajOpt.CommonInterface: eval_jacobian, eval_hessian_of_lagrangian
using TestItemRunner
using TrajectoryIndexingUtils
using DataInterpolations
using DataInterpolations: ExtrapolationType
using ..Quantum
using ..Quantum:
    QuantumSystem,
    KetTrajectory,
    UnitaryTrajectory,
    MultiKetTrajectory,
    AbstractQuantumTrajectory
using ..Quantum: DensityTrajectory, MultiDensityTrajectory, OpenQuantumSystem
using ..Quantum: SamplingTrajectory, sampling_member_states
import ..Quantum: get_system, drive_name, state_name, state_names, get_pulse, duration
using ..Quantum: Isomorphisms, compact_lindbladian_parts
using ..Quantum:
    AbstractDrive,
    LinearDrive,
    NonlinearDrive,
    drive_coeff,
    drive_coeff_jac,
    drive_coeff_hess,
    drive_matrix,
    has_nonlinear_drives,
    active_controls
# The shared operator seam (slice 3b): the abstract dynamics-operator layer and
# the MatrixOperator bridge moved from Piccolissimo's Operators module so the
# dense sensitivity ODE builders run here without a proprietary dependency.
using ..Quantum:
    AbstractDynamicsOperator,
    MatrixOperator,
    apply!,
    state_dim,
    materialize,
    _to_operator,
    dynamics_operator
import SciMLOperators
const SciMLMatrixOperator = SciMLOperators.MatrixOperator

const ⊗ = kron

abstract type AbstractSplineIntegrator <: AbstractIntegrator end

# ── Extension seams for Piccolissimo's matrix-free cells (slice 3b) ──────── #
#
# Declared here (empty) so the moved dense code can dispatch to the
# matrix-free machinery that remains in Piccolissimo, exactly like 3a's
# `matrix_free_jacobian_op` hook. Piccolissimo's module loads AFTER this
# submodule and attaches its methods; the dependency arrow stays strictly
# Piccolissimo → Piccolo.

"""
    _build_magnus_ket_data(alg, sys, u_dim, control_dim, ketdim, order, traj, u, global_dim, global_names)

Hook building the Ket Magnus cell's algorithm data (a `ChebyshevData` — the
matrix-free midpoint-Magnus cores are shared with the Chebyshev cell). Declared
empty here; the concrete method stays proprietary in Piccolissimo
(`integrators/shared/alg_data.jl`, ADR-0003).
"""
function _build_magnus_ket_data end

"""
    _chebyshev_forward!(alg_data::ChebyshevData, spline_type, k, pₖ, ψₖ)

Hook for the matrix-free exp-action forward through the standalone cores
(`magnus_forward!` / `chebyshev_expv!`) — the propagator Φ is never
materialized. Declared empty here; the concrete method stays proprietary in
Piccolissimo. Serves both ChebyshevAlg and the Magnus spline cell on
KetTrajectory.
"""
function _chebyshev_forward! end

"""
    _stiff_rodas5p_solve(prob, tol::Float64; saveat = nothing, save_everystep = true)

Hook for the STIFF forward/sensitivity solves of the `Rodas5PAlg` cells.
Declared empty here (slice 3b, director de-scope 2026-08-30: the Rosenbrock
hard dep is resolver-unsatisfiable against ExponentialUtilities 1.35.1 —
ExponentialUtilities 1.35.1 needs LinearSolve 5.x, every Rosenbrock version
caps LinearSolve below it). `Rodas5PAlg`/`Rodas5PData` — the dataless tag and
the parametric workspace container — stay in Piccolo exactly like
`ChebyshevData`; the CONCRETE `Rodas5P` solver and its construction live in
Piccolissimo, which attaches this method. The method is
`solve(prob, Rodas5P(autodiff = AutoFiniteDiff()); abstol = tol, reltol = tol,
saveat = saveat, save_everystep = save_everystep)` — byte-identical to the
pre-split call sites.
"""
function _stiff_rodas5p_solve end

"""
    matrix_free_layout(𝒮, traj) -> SplineMatrixFreeLayout

Hook returning the cached layout-keyed index tables for the matrix-free
products. Declared empty here because the moved dense Multiket `eval_jacobian`
dispatches to it (the `#205` matrix-free hook); the layout machinery —
`SplineMatrixFreeLayout`, the value-keyed cache and the memo — stays proprietary
in Piccolissimo, whose module attaches the concrete method (the 3a seam
pattern).
"""
function matrix_free_layout end

"""
    _multiket_directional_hvp_probs(𝒮) -> Union{Nothing,Vector{ODEProblem}}

Hook probed by the moved dense MultiKet `eval_hessian_of_lagrangian`: whether
this cell carries matrix-free directional HVP problems. Declared empty here;
the concrete method stays proprietary in Piccolissimo (the #205/#335 MF seam).
"""
function _multiket_directional_hvp_probs end

include("spline_types.jl")
include("propagator_result.jl")
include("algorithms.jl")
include("alg_data.jl")
include("spline_interval_coeffs.jl")
include("complex_real_interface.jl")
include("knot_point_sensitivity.jl")
include("vjp_jvp.jl")
include("spline_integrator_type.jl")
include("spline_integrator_unitary.jl")
include("spline_integrator_ket.jl")
include("spline_integrator_multiket.jl")
include("spline_integrator_sampling.jl")
include("spline_integrator_density.jl")
include("spline_integrator_multidensity.jl")

function unitary_rollout_trajectory(
    u_fn::Function,
    G::Function,
    T::Float64;
    samples::Int = 100,
    kwargs...,
)
    ketdim = size(G(u_fn(0.0), 0.0), 1) ÷ 2

    Id = I(ketdim)

    Ũ⃗_init = operator_to_iso_vec(1.0I(ketdim))

    f! = (dx, x, p, t) -> mul!(dx, Id ⊗ G(u_fn(t), t), x)

    prob = ODEProblem(f!, Ũ⃗_init, (0.0, T))

    times = collect(range(0.0, T, samples))

    Ũ⃗_traj = stack(solve(prob, Tsit5(); abstol = 1e-12, reltol = 1e-12, saveat = times).u)

    return _unitary_trajectory(Ũ⃗_traj, stack([u_fn(t) for t ∈ times]), times; kwargs...)
end

# Helper function for unitary_rollout_trajectory
function _unitary_trajectory(
    Ũ⃗_traj::AbstractMatrix,
    controls::AbstractMatrix,
    times::AbstractVector;
    U_goal = nothing,
    control_bounds = nothing,
)
    u_dim = size(controls, 1)
    ketdim = Int(sqrt(size(Ũ⃗_traj, 1) ÷ 2))
    Δt = diff(times)
    Δt = [Δt; Δt[end]]

    data = (Ũ⃗ = Ũ⃗_traj, u = controls, Δt = Δt, t = times)
    initial = (Ũ⃗ = Isomorphisms.operator_to_iso_vec(1.0I(ketdim)), u = zeros(u_dim))
    final = (u = zeros(u_dim),)
    goal = isnothing(U_goal) ? (;) : (Ũ⃗ = Isomorphisms.operator_to_iso_vec(U_goal),)
    bounds = (
        Ũ⃗ = (-ones(size(Ũ⃗_traj, 1)), ones(size(Ũ⃗_traj, 1))),
        Δt = (1e-3minimum(data.Δt), 2maximum(data.Δt)),
    )
    if !isnothing(control_bounds)
        bounds = merge(bounds, (u = control_bounds,))
    end

    return NamedTrajectory(
        data;
        controls = (:u,),
        timestep = :Δt,
        bounds = bounds,
        initial = initial,
        final = final,
        goal = goal,
    )
end

end
