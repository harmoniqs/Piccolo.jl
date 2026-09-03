# ============================================================================ #
# AdjointRobustness — the adjoint-robustness objective family (Piccolo home).  #
#                                                                              #
# Open-core slice 3d (harmoniqs/Piccolissimo.jl#432). This module owns the      #
# published adjoint (variational) robustness objectives                         #
# (`AdjointRobustnessObjective`, `KetAdjointRobustnessObjective`), the          #
# propagation methods (`ExponentialPropagation`, `SplineODEPropagation`,        #
# `ZKickPropagation`), the Gauss-Newton HVP/VJP/JVP primitives defined ON the   #
# objectives, the `RobustControlProblem` robust wrapper, and their colocated    #
# testitems — moved from Piccolissimo per docs/moved-file-review.md row 28 of   #
# the open-core split's manifest.                                               #
#                                                                              #
# What STAYED in Piccolissimo (proprietary): the solver-backend HVP machinery   #
# itself — the Altissimo backend paths, `knot_hvp` carriers, loss recognition   #
# and the device objective walkers. They consume this surface exactly as        #
# external consumers would (the re-export seam plus `build_objective_callbacks` #
# calling `Objectives.get_full_hessian` / `Objectives.gradient!` on the moved   #
# objectives), so the dependency arrow ends Piccolissimo -> Piccolo and never   #
# the reverse (the slice's Key Decision).                                      #
# ============================================================================ #

module AdjointRobustness

export AdjointRobustnessObjective, KetAdjointRobustnessObjective, RobustControlProblem
export ZKickPropagation
export test_robustness_objective
export objective_hvp!

using LinearAlgebra
using SparseArrays
using ForwardDiff
using FiniteDiff
using NamedTrajectories
using TrajectoryIndexingUtils
using DirectTrajOpt
using DirectTrajOpt: AbstractObjective, objective_value, gradient!, get_full_hessian
using DirectTrajOpt: DirectTrajOptProblem, AbstractConstraint
using DirectTrajOpt.Objectives: hessian_structure, test_objective
using ..Quantum.Isomorphisms: Isomorphisms, iso_vec_to_operator, operator_to_iso_vec
using ..Quantum.QuantumTrajectories: UnitaryTrajectory
using ..Quantum.Isomorphisms: iso_to_ket, ket_to_iso
using ..Quantum.QuantumTrajectories: state_name as piccolo_state_name
using ..Quantum.Pulses: drive_name as piccolo_drive_name
using ..Quantum.QuantumTrajectories: state_names as piccolo_state_names
using ..Quantum.QuantumSystems: VariationalQuantumSystem
using ..QuantumControlProblems: QuantumControlProblem, get_trajectory
using ..Quantum.QuantumTrajectories: AbstractQuantumTrajectory
using ..QuantumConstraints: FinalUnitaryFidelityConstraint
using ..Quantum.QuantumTrajectories: KetTrajectory, MultiKetTrajectory
using ..QuantumConstraints: FinalKetFidelityConstraint, FinalCoherentKetFidelityConstraint
using NamedTrajectories.StructKnotPoint: KnotPoint
using ..QuantumIntegrators.ExponentialIntegrators: HermitianExponentialIntegrator
using ..QuantumIntegrators.ExponentialIntegrators: x_name, exp_eigen
import ..QuantumIntegrators.SplineIntegrators
using ..QuantumIntegrators.SplineIntegrators: SplineIntegrator
using OrdinaryDiffEqTsit5: ODEProblem, solve, Tsit5
using ExponentialAction: expv
using TestItemRunner: @testitem
import Test

# ============================================================================ #
#                         Propagation Methods                                  #
# ============================================================================ #

"""
    PropagationMethod

Abstract type for how the adjoint robustness objective propagates V between knot points.
"""
abstract type PropagationMethod end

"""
    ExponentialPropagation <: PropagationMethod

Use block matrix exponential (for piecewise-constant / zero-order-hold controls).
Compatible with: HermitianExponentialIntegrator.
"""
struct ExponentialPropagation <: PropagationMethod
    G::Function     # u -> 2d×2d iso generator
end

"""
    SplineODEPropagation <: PropagationMethod

Use augmented ODE integration (for spline-interpolated controls).
Compatible with: SplineIntegrator (all algorithm variants).
"""
struct SplineODEPropagation <: PropagationMethod
    G::Function     # u -> 2d×2d iso generator
    order::Int      # 1=linear, 3=cubic
    tol::Float64    # ODE solver tolerance
    derivative_name::Union{Nothing,Symbol}  # :du for cubic, nothing for linear
end

"""
    ZKickPropagation <: PropagationMethod

Propagation for Z-kick + free-evolution dynamics: P_k = Rz(ϕ_k) · exp(Δt_k · G_H).

The error coupling occurs during free evolution only (the Z-kick is instantaneous).
The block exponential computes exp(Δt · [G_H, 0; G_E, G_H]) for the free-evolution
part, then the Z-kick rotation is applied after.
"""
struct ZKickPropagation <: PropagationMethod
    G_H::Matrix{Float64}   # iso(-iH) — free-evolution generator
end

# ============================================================================ #
#                       AdjointRobustnessObjective                             #
# ============================================================================ #

"""
    AdjointRobustnessObjective <: AbstractObjective

Objective for minimizing first-order error susceptibility using the adjoint
(variational) method from Kamen et al., arXiv:2602.10349.

Computes:
    E_V = (Q / t_f²) Σ_m ||U_N† V_N^(m)||_F² / d

where V_N^(m) is the adjoint state propagated forward via block matrix exponential
(for `ExponentialPropagation`) or augmented ODE (for `SplineODEPropagation`).

# Fields
- `G_errors::Vector{Matrix{Float64}}`: Error operators in iso form, G_E = iso(-iE)
- `propagation::PropagationMethod`: Propagation method (exponential or spline ODE)
- `state_name::Symbol`: Name of unitary state variable (e.g., :Ũ⃗)
- `control_name::Symbol`: Name of control variable (e.g., :a or :u)
- `ketdim::Int`: Hilbert space dimension d
- `iso_dim::Int`: 2d (iso matrix row dimension)
- `state_dim::Int`: 2d² (iso-vectorized state dimension)
- `Q::Float64`: Overall weight
- `exact_hessian::Bool`: If `true`, `get_full_hessian` returns the exact Hessian
  (via finite differences over `gradient!`) rather than the Gauss–Newton
  approximation `J^T J`. Defaults to `false`. Set this when the system has
  `NonlinearDrive`s — see `Piccolissimo.jl/gauss_newton_explanation.md`: the
  GN approximation drops a `Σ_d (∂²c_d/∂u²) G_d` term that is identically
  zero for `LinearDrive` but nonzero for nonlinear drive coefficients.
"""
struct AdjointRobustnessObjective <: AbstractObjective
    G_errors::Vector{Matrix{Float64}}
    propagation::PropagationMethod
    state_name::Symbol
    control_name::Symbol
    ketdim::Int
    iso_dim::Int
    state_dim::Int
    Q::Float64
    exact_hessian::Bool
end

# Backward-compat outer ctor: old positional callers (no exact_hessian) get
# `false` by default.
AdjointRobustnessObjective(
    G_errors::Vector{<:AbstractMatrix},
    propagation::PropagationMethod,
    state_name::Symbol,
    control_name::Symbol,
    ketdim::Int,
    iso_dim::Int,
    state_dim::Int,
    Q::Float64,
) = AdjointRobustnessObjective(
    G_errors,
    propagation,
    state_name,
    control_name,
    ketdim,
    iso_dim,
    state_dim,
    Q,
    false,
)

function Base.show(io::IO, obj::AdjointRobustnessObjective)
    n_err = length(obj.G_errors)
    print(io, "AdjointRobustnessObjective(d=$(obj.ketdim), $(n_err) error ops, Q=$(obj.Q))")
end

# ============================================================================ #
#                            Constructors                                       #
# ============================================================================ #

"""
    AdjointRobustnessObjective(
        integrator::HermitianExponentialIntegrator{UnitaryTrajectory},
        error_operators::Vector{<:AbstractMatrix},
        traj::NamedTrajectory;
        Q::Float64 = 1.0,
    )

Construct from a HermitianExponentialIntegrator and error operators.

# Arguments
- `integrator`: The dynamics integrator (provides G(u) closure and dimensions)
- `error_operators`: Hermitian matrices [E1, E2, ...] representing error channels
- `traj`: Trajectory (used for dimension inference)
- `Q`: Weight on the robustness objective
"""
function AdjointRobustnessObjective(
    integrator::HermitianExponentialIntegrator{UnitaryTrajectory},
    error_operators::Vector{<:AbstractMatrix},
    traj::NamedTrajectory;
    Q::Float64 = 1.0,
    exact_hessian::Bool = false,
)
    ketdim = integrator.ketdim
    iso_dim = 2 * ketdim
    state_dim = iso_dim * ketdim  # 2d²

    G_errors = [Matrix{Float64}(Isomorphisms.G(ComplexF64.(E))) for E in error_operators]

    return AdjointRobustnessObjective(
        G_errors,
        ExponentialPropagation(integrator.G),
        x_name(integrator),
        integrator.u_name,
        ketdim,
        iso_dim,
        state_dim,
        Q,
        exact_hessian,
    )
end

"""
    AdjointRobustnessObjective(
        sys::QuantumSystem,
        error_operators::Vector{<:AbstractMatrix},
        state_name::Symbol,
        control_name::Symbol,
        traj::NamedTrajectory;
        Q::Float64 = 1.0,
    )

Construct from a QuantumSystem directly.
"""
function AdjointRobustnessObjective(
    sys,  # QuantumSystem
    error_operators::Vector{<:AbstractMatrix},
    state_nm::Symbol,
    control_nm::Symbol,
    traj::NamedTrajectory;
    Q::Float64 = 1.0,
    exact_hessian::Bool = false,
)
    ketdim = size(error_operators[1], 1)
    iso_dim = 2 * ketdim
    state_dim = iso_dim * ketdim

    G_errors = [Matrix{Float64}(Isomorphisms.G(ComplexF64.(E))) for E in error_operators]
    G = u -> sys.G(u, 0.0)

    return AdjointRobustnessObjective(
        G_errors,
        ExponentialPropagation(G),
        state_nm,
        control_nm,
        ketdim,
        iso_dim,
        state_dim,
        Q,
        exact_hessian,
    )
end

"""
    AdjointRobustnessObjective(
        integrator::SplineIntegrator{UnitaryTrajectory},
        sys,
        error_operators::Vector{<:AbstractMatrix},
        traj::NamedTrajectory;
        Q::Float64 = 1.0,
    )

Construct from a SplineIntegrator and QuantumSystem. Uses ODE-based propagation
that respects spline-interpolated controls (linear or cubic Hermite).
"""
function AdjointRobustnessObjective(
    integrator::SplineIntegrator{UnitaryTrajectory},
    sys,  # AbstractQuantumSystem
    error_operators::Vector{<:AbstractMatrix},
    traj::NamedTrajectory;
    Q::Float64 = 1.0,
    exact_hessian::Bool = false,
)
    ketdim = integrator.ketdim
    iso_dim = 2 * ketdim
    state_dim = iso_dim * ketdim

    G_errors = [Matrix{Float64}(Isomorphisms.G(ComplexF64.(E))) for E in error_operators]
    G = u -> sys.G(u, 0.0)

    order = SplineIntegrators.spline_order(integrator)
    u_name = integrator.u_name
    derivative_name = order == 3 ? Symbol("d" * String(u_name)) : nothing

    prop = SplineODEPropagation(G, order, integrator.tol, derivative_name)

    return AdjointRobustnessObjective(
        G_errors,
        prop,
        integrator.x_names[1],
        u_name,
        ketdim,
        iso_dim,
        state_dim,
        Q,
        exact_hessian,
    )
end

"""
    AdjointRobustnessObjective(
        varsys::VariationalQuantumSystem,
        state_name::Symbol,
        control_name::Symbol,
        traj::NamedTrajectory;
        Q::Float64 = 1.0,
    )

Construct from a VariationalQuantumSystem. The variational generators `G_vars`
are used as error operators.
"""
function AdjointRobustnessObjective(
    varsys::VariationalQuantumSystem,
    state_nm::Symbol,
    control_nm::Symbol,
    traj::NamedTrajectory;
    Q::Float64 = 1.0,
    exact_hessian::Bool = false,
)
    ketdim = varsys.levels
    iso_dim = 2 * ketdim
    state_dim = iso_dim * ketdim

    # G_vars[m](a) = iso(-i * H_var_m), which is exactly G_E
    a_zero = zeros(varsys.n_drives)
    G_errors = [Matrix{Float64}(gv(a_zero)) for gv in varsys.G_vars]

    G = u -> varsys.G(u)

    return AdjointRobustnessObjective(
        G_errors,
        ExponentialPropagation(G),
        state_nm,
        control_nm,
        ketdim,
        iso_dim,
        state_dim,
        Q,
        exact_hessian,
    )
end

"""
    AdjointRobustnessObjective(
        H::AbstractMatrix,
        error_operators::Vector{<:AbstractMatrix},
        state_name::Symbol,
        control_name::Symbol,
        traj::NamedTrajectory;
        Q::Float64 = 1.0,
    )

Construct from a Hamiltonian matrix for Z-kick dynamics. Uses `ZKickPropagation`
which computes P_k = Rz(ϕ_k) · exp(Δt_k · G_H) at each step.

The control variable contains Z-kick angles ϕ_k. The error operators act during
free evolution only (the Z-kick is instantaneous).
"""
function AdjointRobustnessObjective(
    H::AbstractMatrix{<:Complex},
    error_operators::Vector{<:AbstractMatrix},
    state_nm::Symbol,
    control_nm::Symbol,
    traj::NamedTrajectory;
    Q::Float64 = 1.0,
    exact_hessian::Bool = false,
)
    ketdim = size(error_operators[1], 1)
    iso_dim = 2 * ketdim
    state_dim = iso_dim * ketdim

    G_errors = [Matrix{Float64}(Isomorphisms.G(ComplexF64.(E))) for E in error_operators]
    G_H = Matrix{Float64}(Isomorphisms.G(ComplexF64.(H)))

    return AdjointRobustnessObjective(
        G_errors,
        ZKickPropagation(G_H),
        state_nm,
        control_nm,
        ketdim,
        iso_dim,
        state_dim,
        Q,
        exact_hessian,
    )
end

# ============================================================================ #
#                     KetAdjointRobustnessObjective                            #
# ============================================================================ #

"""
    KetAdjointRobustnessObjective <: AbstractObjective

Objective for minimizing first-order error susceptibility of ket trajectories
using the adjoint (variational) method from Kamen et al., arXiv:2602.10349.

Computes:
    E_V = Q Σ_m Σ_i |⟨ψ_goal_i | v_N^(m,i)⟩|² / n_kets

where V_N^(m) is the adjoint state (iso_dim × n_kets) propagated forward via
block matrix exponential, and v_N^(m,i) is the i-th column converted to a ket.

# Fields
- `G_errors::Vector{Matrix{Float64}}`: Error operators in iso form, G_E = iso(-iE)
- `propagation::PropagationMethod`: Propagation method (exponential or spline ODE)
- `state_names::Vector{Symbol}`: Names of ket state variables (e.g., [:ψ̃1, :ψ̃2, ...])
- `goals::Vector{Vector{ComplexF64}}`: Target ket states
- `control_name::Symbol`: Name of control variable (e.g., :u)
- `ketdim::Int`: Hilbert space dimension d
- `iso_dim::Int`: 2d (iso ket dimension)
- `n_kets::Int`: Number of kets in ensemble
- `Q::Float64`: Overall weight
- `exact_hessian::Bool`: see `AdjointRobustnessObjective`
"""
struct KetAdjointRobustnessObjective <: AbstractObjective
    G_errors::Vector{Matrix{Float64}}
    propagation::PropagationMethod
    state_names::Vector{Symbol}
    goals::Vector{Vector{ComplexF64}}
    control_name::Symbol
    ketdim::Int
    iso_dim::Int
    n_kets::Int
    Q::Float64
    exact_hessian::Bool
end

KetAdjointRobustnessObjective(
    G_errors::Vector{<:AbstractMatrix},
    propagation::PropagationMethod,
    state_names::Vector{Symbol},
    goals::Vector{<:AbstractVector},
    control_name::Symbol,
    ketdim::Int,
    iso_dim::Int,
    n_kets::Int,
    Q::Float64,
) = KetAdjointRobustnessObjective(
    G_errors,
    propagation,
    state_names,
    [Vector{ComplexF64}(g) for g in goals],
    control_name,
    ketdim,
    iso_dim,
    n_kets,
    Q,
    false,
)

function Base.show(io::IO, obj::KetAdjointRobustnessObjective)
    n_err = length(obj.G_errors)
    print(
        io,
        "KetAdjointRobustnessObjective(d=$(obj.ketdim), $(obj.n_kets) kets, $(n_err) error ops, Q=$(obj.Q))",
    )
end

# ============================================================================ #
#                  KetAdjointRobustnessObjective Constructors                   #
# ============================================================================ #

"""
    KetAdjointRobustnessObjective(
        integrator::HermitianExponentialIntegrator{MultiKetTrajectory},
        error_operators::Vector{<:AbstractMatrix},
        goals::Vector{<:AbstractVector},
        traj::NamedTrajectory;
        Q::Float64 = 1.0,
    )

Construct from a HermitianExponentialIntegrator for MultiKetTrajectory.
"""
function KetAdjointRobustnessObjective(
    integrator::HermitianExponentialIntegrator{MultiKetTrajectory},
    error_operators::Vector{<:AbstractMatrix},
    goals::Vector{<:AbstractVector},
    traj::NamedTrajectory;
    Q::Float64 = 1.0,
    exact_hessian::Bool = false,
)
    ketdim = integrator.ketdim
    iso_dim = 2 * ketdim
    n_kets = length(integrator.x_names)

    G_errors = [Matrix{Float64}(Isomorphisms.G(ComplexF64.(E))) for E in error_operators]

    return KetAdjointRobustnessObjective(
        G_errors,
        ExponentialPropagation(integrator.G),
        collect(integrator.x_names),
        [Vector{ComplexF64}(g) for g in goals],
        integrator.u_name,
        ketdim,
        iso_dim,
        n_kets,
        Q,
        exact_hessian,
    )
end

"""
    KetAdjointRobustnessObjective(
        integrator::HermitianExponentialIntegrator{KetTrajectory},
        error_operators::Vector{<:AbstractMatrix},
        goal::AbstractVector,
        traj::NamedTrajectory;
        Q::Float64 = 1.0,
    )

Construct from a HermitianExponentialIntegrator for KetTrajectory (single ket).
"""
function KetAdjointRobustnessObjective(
    integrator::HermitianExponentialIntegrator{KetTrajectory},
    error_operators::Vector{<:AbstractMatrix},
    goal::AbstractVector,
    traj::NamedTrajectory;
    Q::Float64 = 1.0,
    exact_hessian::Bool = false,
)
    ketdim = integrator.ketdim
    iso_dim = 2 * ketdim

    G_errors = [Matrix{Float64}(Isomorphisms.G(ComplexF64.(E))) for E in error_operators]

    return KetAdjointRobustnessObjective(
        G_errors,
        ExponentialPropagation(integrator.G),
        collect(integrator.x_names),
        [Vector{ComplexF64}(goal)],
        integrator.u_name,
        ketdim,
        iso_dim,
        1,
        Q,
        exact_hessian,
    )
end

# ============================================================================ #
#                        Forward Propagation                                    #
# ============================================================================ #

"""
    ForwardPassResult

Stores all intermediates from the forward propagation of the adjoint state,
needed for both objective evaluation and backward gradient pass.
"""
struct ForwardPassResult
    # Per-error-operator adjoint states at each knot: V[m][k] is iso_dim × ketdim
    Vs::Vector{Vector{Matrix{Float64}}}
    # Propagators at each knot point: P[k] is iso_dim × iso_dim
    Ps::Vector{Matrix{Float64}}
    # Cross-coupling blocks: Φ21[m][k] is iso_dim × iso_dim
    Φ21s::Vector{Vector{Matrix{Float64}}}
    # Total time
    t_f::Float64
end

"""
    forward_pass(obj::AdjointRobustnessObjective, traj::NamedTrajectory) -> ForwardPassResult

Propagate the adjoint state V forward through all knot points for each error operator.
Dispatches on `obj.propagation` to use block matrix exponential or ODE integration.
"""
function forward_pass(obj::AdjointRobustnessObjective, traj::NamedTrajectory)
    return _forward_pass(obj.propagation, obj, traj)
end

"""
    forward_pass(obj::KetAdjointRobustnessObjective, traj::NamedTrajectory) -> ForwardPassResult

Propagate the adjoint state V forward for ket trajectories. Gathers kets from
multiple state variables into a concatenated matrix Ψ_k (iso_dim × n_kets).
"""
function forward_pass(obj::KetAdjointRobustnessObjective, traj::NamedTrajectory)
    return _forward_pass_ket(obj.propagation, obj, traj)
end

# --- Exponential propagation (zero-order hold) ---

function _forward_pass(prop::ExponentialPropagation, obj, traj)
    N = traj.N
    d2 = obj.iso_dim
    d = obj.ketdim
    n_errors = length(obj.G_errors)

    t_f = sum(traj[k].timestep for k = 1:(N-1))

    Ps = Vector{Matrix{Float64}}(undef, N - 1)
    Φ21s = [Vector{Matrix{Float64}}(undef, N - 1) for _ = 1:n_errors]
    Vs = [Vector{Matrix{Float64}}(undef, N) for _ = 1:n_errors]

    for m = 1:n_errors
        Vs[m][1] = zeros(d2, d)
    end

    for k = 1:(N-1)
        zk = traj[k]
        u_k = zk[obj.control_name]
        Δt_k = zk.timestep
        Ũ_k = reshape(zk[obj.state_name], d2, d)

        G_k = prop.G(u_k)

        for m = 1:n_errors
            G_E = obj.G_errors[m]

            aug_G = zeros(2d2, 2d2)
            aug_G[1:d2, 1:d2] .= G_k
            aug_G[(d2+1):2d2, 1:d2] .= G_E
            aug_G[(d2+1):2d2, (d2+1):2d2] .= G_k

            aug_exp = exp(Δt_k * aug_G)

            if m == 1
                Ps[k] = aug_exp[1:d2, 1:d2]
            end
            Φ21s[m][k] = aug_exp[(d2+1):2d2, 1:d2]

            Vs[m][k+1] = Ps[k] * Vs[m][k] + Φ21s[m][k] * Ũ_k
        end
    end

    return ForwardPassResult(Vs, Ps, Φ21s, t_f)
end

# --- Z-kick propagation ---

"""
Build Rz(ϕ) in iso form (2d × 2d real matrix).
"""
function _iso_Rz(ϕ)
    # Generic in eltype(ϕ) so ForwardDiff.Dual works through _add_control_gradient!
    c = exp(-im * ϕ / 2)
    Rz = Diagonal([c, conj(c)])
    return Matrix(Isomorphisms.iso(Rz))
end

function _forward_pass(prop::ZKickPropagation, obj, traj)
    N = traj.N
    d2 = obj.iso_dim
    d = obj.ketdim
    n_errors = length(obj.G_errors)
    G_H = prop.G_H

    t_f = sum(traj[k].timestep for k = 1:(N-1))

    Ps = Vector{Matrix{Float64}}(undef, N - 1)
    Φ21s = [Vector{Matrix{Float64}}(undef, N - 1) for _ = 1:n_errors]
    Vs = [Vector{Matrix{Float64}}(undef, N) for _ = 1:n_errors]

    for m = 1:n_errors
        Vs[m][1] = zeros(d2, d)
    end

    for k = 1:(N-1)
        zk = traj[k]
        ϕ_k = zk[obj.control_name][1]
        Δt_k = zk.timestep
        Ũ_k = reshape(zk[obj.state_name], d2, d)

        R_k = _iso_Rz(ϕ_k)

        for m = 1:n_errors
            G_E = obj.G_errors[m]

            # Block exponential for free evolution + error coupling
            aug_G = zeros(2d2, 2d2)
            aug_G[1:d2, 1:d2] .= G_H
            aug_G[(d2+1):2d2, 1:d2] .= G_E
            aug_G[(d2+1):2d2, (d2+1):2d2] .= G_H

            aug_exp = exp(Δt_k * aug_G)
            F_k = aug_exp[1:d2, 1:d2]
            Φ21_free = aug_exp[(d2+1):2d2, 1:d2]

            # Apply Z-kick after free evolution
            if m == 1
                Ps[k] = R_k * F_k
            end
            Φ21s[m][k] = R_k * Φ21_free

            Vs[m][k+1] = Ps[k] * Vs[m][k] + Φ21s[m][k] * Ũ_k
        end
    end

    return ForwardPassResult(Vs, Ps, Φ21s, t_f)
end

# --- Spline ODE propagation ---

"""
Interpolate controls at normalized time τ ∈ [0,1] for spline propagation.
"""
function _interpolate_controls(τ, u_k, u_k1, Δt_k, order, du_k = nothing, du_k1 = nothing)
    if order == 1
        return (1 - τ) .* u_k .+ τ .* u_k1
    else  # order == 3
        h00 = 2τ^3 - 3τ^2 + 1
        h10 = (τ^3 - 2τ^2 + τ) * Δt_k
        h01 = -2τ^3 + 3τ^2
        h11 = (τ^3 - τ^2) * Δt_k
        return h00 .* u_k .+ h10 .* du_k .+ h01 .* u_k1 .+ h11 .* du_k1
    end
end

function _forward_pass(prop::SplineODEPropagation, obj, traj)
    N = traj.N
    d2 = obj.iso_dim
    d = obj.ketdim
    n_errors = length(obj.G_errors)
    Φ_dim = d2^2

    t_f = sum(traj[k].timestep for k = 1:(N-1))

    Ps = Vector{Matrix{Float64}}(undef, N - 1)
    Φ21s = [Vector{Matrix{Float64}}(undef, N - 1) for _ = 1:n_errors]
    Vs = [Vector{Matrix{Float64}}(undef, N) for _ = 1:n_errors]

    for m = 1:n_errors
        Vs[m][1] = zeros(d2, d)
    end

    # Augmented state: [vec(Φ); vec(W₁); vec(W₂); ...] where Wₘ is per-error cross-coupling
    aug_dim = Φ_dim * (1 + n_errors)

    for k = 1:(N-1)
        u_k = traj[k][obj.control_name]
        u_k1 = traj[k+1][obj.control_name]
        Δt_k = traj[k].timestep
        Ũ_k = reshape(traj[k][obj.state_name], d2, d)

        # Get derivatives for cubic spline
        du_k = nothing
        du_k1 = nothing
        if prop.order == 3 && !isnothing(prop.derivative_name)
            du_k = traj[k][prop.derivative_name]
            du_k1 = traj[k+1][prop.derivative_name]
        end

        # Build augmented ODE
        G_errors_local = obj.G_errors
        G_func = prop.G
        order = prop.order

        function aug_ode!(dx, x, _p, τ)
            u_τ = _interpolate_controls(τ, u_k, u_k1, Δt_k, order, du_k, du_k1)
            G_τ = G_func(u_τ)

            Φ_mat = reshape(@view(x[1:Φ_dim]), d2, d2)
            dΦ_mat = reshape(@view(dx[1:Φ_dim]), d2, d2)

            # dΦ/dτ = Δt * G(u(τ)) * Φ
            mul!(dΦ_mat, G_τ, Φ_mat, Δt_k, false)

            for mm = 1:n_errors
                offset = mm * Φ_dim
                W_mat = reshape(@view(x[(offset+1):(offset+Φ_dim)]), d2, d2)
                dW_mat = reshape(@view(dx[(offset+1):(offset+Φ_dim)]), d2, d2)

                # dW/dτ = Δt * G(u(τ)) * W + Δt * G_E * Φ
                mul!(dW_mat, G_τ, W_mat, Δt_k, false)
                mul!(dW_mat, G_errors_local[mm], Φ_mat, Δt_k, true)
            end
        end

        # Initial conditions: Φ(0) = I, W_m(0) = 0
        x0 = zeros(aug_dim)
        for i = 1:d2
            x0[(i-1)*d2+i] = 1.0  # Identity for Φ
        end

        prob = ODEProblem(aug_ode!, x0, (0.0, 1.0))
        sol = solve(
            prob,
            Tsit5();
            abstol = prop.tol,
            reltol = prop.tol,
            save_everystep = false,
        )
        final_state = sol.u[end]

        # Extract propagator P_k
        Ps[k] = reshape(final_state[1:Φ_dim], d2, d2)

        # Extract cross-couplings Ψ_k for each error operator
        for m = 1:n_errors
            offset = m * Φ_dim
            Φ21s[m][k] = reshape(final_state[(offset+1):(offset+Φ_dim)], d2, d2)

            # Chain: V_{k+1} = P_k * V_k + Ψ_k * Ũ_k
            Vs[m][k+1] = Ps[k] * Vs[m][k] + Φ21s[m][k] * Ũ_k
        end
    end

    return ForwardPassResult(Vs, Ps, Φ21s, t_f)
end

# --- Ket-specific forward pass (exponential propagation) ---

function _forward_pass_ket(
    prop::ExponentialPropagation,
    obj::KetAdjointRobustnessObjective,
    traj,
)
    N = traj.N
    d2 = obj.iso_dim
    n = obj.n_kets
    n_errors = length(obj.G_errors)

    t_f = sum(traj[k].timestep for k = 1:(N-1))

    Ps = Vector{Matrix{Float64}}(undef, N - 1)
    Φ21s = [Vector{Matrix{Float64}}(undef, N - 1) for _ = 1:n_errors]
    Vs = [Vector{Matrix{Float64}}(undef, N) for _ = 1:n_errors]

    for m = 1:n_errors
        Vs[m][1] = zeros(d2, n)
    end

    for k = 1:(N-1)
        zk = traj[k]
        u_k = zk[obj.control_name]
        Δt_k = zk.timestep

        # Gather kets into matrix Ψ_k (iso_dim × n_kets)
        Ψ_k = zeros(d2, n)
        for (i, sname) in enumerate(obj.state_names)
            Ψ_k[:, i] = zk[sname]
        end

        G_k = prop.G(u_k)

        for m = 1:n_errors
            G_E = obj.G_errors[m]

            aug_G = zeros(2d2, 2d2)
            aug_G[1:d2, 1:d2] .= G_k
            aug_G[(d2+1):2d2, 1:d2] .= G_E
            aug_G[(d2+1):2d2, (d2+1):2d2] .= G_k

            aug_exp = exp(Δt_k * aug_G)

            if m == 1
                Ps[k] = aug_exp[1:d2, 1:d2]
            end
            Φ21s[m][k] = aug_exp[(d2+1):2d2, 1:d2]

            Vs[m][k+1] = Ps[k] * Vs[m][k] + Φ21s[m][k] * Ψ_k
        end
    end

    return ForwardPassResult(Vs, Ps, Φ21s, t_f)
end

# ============================================================================ #
#                  Residual computation (for GN Hessian)                        #
# ============================================================================ #
#
# The objective f = Q/d * Σ_m ‖A_m‖²_F where A_m = U_N† V_N^(m), so we can write
# f = ½ ‖r(Z)‖² with residual
#
#     r_m = √(2Q/d) · [vec(Re A_m); vec(Im A_m)]   ∈ ℝ^{2d²}
#
# stacking all error operators gives r ∈ ℝ^{2 d² n_err}. Then the analytic
# Gauss–Newton Hessian is H_GN = J(Z)^T J(Z), Z = traj decision variables.
#
# The functions below compute r(Z) in an element-type-generic way so that
# `ForwardDiff.jacobian(r, Z)` yields J analytically. This replaces the previous
# O(|Z|²) finite-difference Hessian with O(|Z|) ForwardDiff cost (one column of J
# per Dual partial), assembled via J^T J.
#
# These are a parallel, simpler forward pass than `_forward_pass` — they cache
# only the final V_N^(m) (no Ps or Φ21s), and they accept Dual element types.
#
# All `_propagate_final_V` and `_compute_residual` methods take an explicit
# `data::AbstractVector{T}` instead of reading from `traj.datavec`. This avoids
# `NamedTrajectory(traj; datavec=Dual_vec)` reconstruction (which fails because
# bounds/initial/final/goal NamedTuples retain Float64 element type and break
# `inspect_dims_pairs` dispatch). We construct `KnotPoint`s manually from
# `data`, reusing `traj`'s component metadata.

# --- Iso ↔ complex helpers (zero-allocation views where possible) ---

"""
    _iso_to_minus_i_H(G_iso) -> Matrix{Complex{T}}

Convert an iso-form generator `G_iso = iso(-iH)` (real `2d × 2d`) to its complex
form `-iH` (`d × d`). For `H` Hermitian (the standard case), `iso(-iH) = [Im(H)
Re(H); -Re(H) Im(H)]`, so `-iH = G_iso[1:d, 1:d] + i · G_iso[(d+1):2d, 1:d]`.
"""
@inline function _iso_to_minus_i_H(G_iso::AbstractMatrix{T}) where {T<:Real}
    d = size(G_iso, 1) ÷ 2
    return view(G_iso, 1:d, 1:d) .+ im .* view(G_iso, (d+1):(2d), 1:d)
end

"""
    _iso_view_to_complex_matrix(Ũ_iso_view) -> Matrix{Complex{T}}

Convert a `2d × n` real iso view (top half = real part, bottom half = imag part
of each column) to a `d × n` complex matrix.
"""
@inline function _iso_view_to_complex_matrix(Ũ_iso_view::AbstractMatrix{T}) where {T<:Real}
    d = size(Ũ_iso_view, 1) ÷ 2
    return view(Ũ_iso_view, 1:d, :) .+ im .* view(Ũ_iso_view, (d+1):(2d), :)
end

"""
    _make_knot_view(traj::NamedTrajectory, data::AbstractVector{T}, k::Int) -> KnotPoint{T}

Build a `KnotPoint` view at knot `k` over an arbitrary-eltype `data` vector,
using `traj`'s components / names / control_names. Bypasses the typed
`NamedTrajectory` constructor so that `T` may be `ForwardDiff.Dual`.
"""
@inline function _make_knot_view(
    traj::NamedTrajectory,
    data::AbstractVector{T},
    k::Int,
) where {T<:Real}
    base = (k - 1) * traj.dim
    data_view = view(data, (base+1):(base+traj.dim))
    ts_idx = first(traj.components[traj.timestep])
    return KnotPoint(
        k,
        data_view,
        data_view[ts_idx],
        traj.components,
        traj.names,
        traj.control_names,
    )
end

# --- Exponential propagation ---

function _propagate_final_V(
    prop::ExponentialPropagation,
    obj::AdjointRobustnessObjective,
    traj::NamedTrajectory,
    data::AbstractVector{T},
) where {T<:Real}
    N = traj.N
    d2 = obj.iso_dim
    d = obj.ketdim
    n_errors = length(obj.G_errors)

    # Vs in COMPLEX d × d (was 2d × d real iso). Saves a factor of ~2 in
    # `expv` matvec cost and avoids a final iso → complex conversion in
    # `_compute_residual`.
    Vs = [zeros(Complex{T}, d, d) for _ = 1:n_errors]

    # Pre-deiso the static error operators once per call.
    miEs = [_iso_to_minus_i_H(G_E) for G_E in obj.G_errors]

    for k = 1:(N-1)
        zk = _make_knot_view(traj, data, k)
        u_k = zk[obj.control_name]
        Δt_k = zk.timestep
        Ũ_k_iso = reshape(zk[obj.state_name], d2, d)
        U_k = _iso_view_to_complex_matrix(Ũ_k_iso)  # d × d complex

        G_k_iso = prop.G(u_k)
        miH_k = _iso_to_minus_i_H(G_k_iso)  # d × d complex

        # Augmented complex matrix `aug_C = [-iH_k 0; -iE_m -iH_k]` of size
        # `2d × 2d`. expv on `[U_col; V_col]` produces `[-iH_k * U_col; ...]`
        # at the top and `(action of P*V + Φ21*U) at the bottom`, which is
        # exactly the new V column — same trick as the iso version, just
        # smaller and complex.
        for m = 1:n_errors
            miE_m = miEs[m]
            aug_C = zeros(Complex{T}, 2d, 2d)
            aug_C[1:d, 1:d] .= miH_k
            aug_C[(d+1):(2d), 1:d] .= miE_m
            aug_C[(d+1):(2d), (d+1):(2d)] .= miH_k

            new_Vm = zeros(Complex{T}, d, d)
            v_aug = Vector{Complex{T}}(undef, 2d)
            for j = 1:d
                v_aug[1:d] .= view(U_k, :, j)
                v_aug[(d+1):(2d)] .= view(Vs[m], :, j)
                result = expv(Δt_k, aug_C, v_aug)
                new_Vm[:, j] .= view(result, (d+1):(2d))
            end
            Vs[m] = new_Vm
        end
    end
    return Vs
end

# --- Z-kick propagation ---

function _propagate_final_V(
    prop::ZKickPropagation,
    obj::AdjointRobustnessObjective,
    traj::NamedTrajectory,
    data::AbstractVector{T},
) where {T<:Real}
    N = traj.N
    d2 = obj.iso_dim
    d = obj.ketdim
    n_errors = length(obj.G_errors)

    # Complex internals (see ExponentialPropagation note).
    Vs = [zeros(Complex{T}, d, d) for _ = 1:n_errors]

    # Pre-deiso static operators
    miH = _iso_to_minus_i_H(prop.G_H)
    miEs = [_iso_to_minus_i_H(G_E) for G_E in obj.G_errors]

    for k = 1:(N-1)
        zk = _make_knot_view(traj, data, k)
        ϕ_k = zk[obj.control_name][1]
        Δt_k = zk.timestep
        Ũ_k_iso = reshape(zk[obj.state_name], d2, d)
        U_k = _iso_view_to_complex_matrix(Ũ_k_iso)

        # Z-kick rotation as 2×2 complex diag (assumes d=2 single-qubit semantics
        # in the existing ZKick model). Generalize via _iso_Rz path: it returns
        # a 2d × 2d real iso, deiso to d × d complex.
        R_k_iso = _iso_Rz(ϕ_k)
        R_k = _iso_view_to_complex_matrix(R_k_iso[:, 1:d])

        # Free-evolution augmented complex matrix [-iH 0; -iE -iH], 2d × 2d
        for m = 1:n_errors
            miE_m = miEs[m]
            aug_C = zeros(Complex{T}, 2d, 2d)
            aug_C[1:d, 1:d] .= miH
            aug_C[(d+1):(2d), 1:d] .= miE_m
            aug_C[(d+1):(2d), (d+1):(2d)] .= miH

            new_Vm = zeros(Complex{T}, d, d)
            v_aug = Vector{Complex{T}}(undef, 2d)
            for j = 1:d
                v_aug[1:d] .= view(U_k, :, j)
                v_aug[(d+1):(2d)] .= view(Vs[m], :, j)
                result = expv(Δt_k, aug_C, v_aug)
                # Bottom block = Φ21_free*U + F*V. Apply Z-kick R_k after.
                new_Vm[:, j] .= R_k * view(result, (d+1):(2d))
            end
            Vs[m] = new_Vm
        end
    end
    return Vs
end

# --- Spline ODE propagation ---

function _propagate_final_V(
    prop::SplineODEPropagation,
    obj::AdjointRobustnessObjective,
    traj::NamedTrajectory,
    data::AbstractVector{T},
) where {T<:Real}
    # SplineODE keeps iso-form internals: the augmented ODE state is
    # (Φ, W₁, …, W_{n_err}) — each Φ is `d2 × d2` real iso, the natural shape
    # for the ODE RHS `dΦ = G(u(τ)) Φ`. Converting to complex internally would
    # require also rewriting the RHS `mul!`s for complex arithmetic. Defer that
    # until we tackle the SplineODE adjoint costate refactor (the spec's PR2).
    #
    # Output Vs is converted to complex at the end so it matches the other
    # _propagate_final_V methods' return type.
    N = traj.N
    d2 = obj.iso_dim
    d = obj.ketdim
    n_errors = length(obj.G_errors)
    Φ_dim = d2^2

    Vs_iso = [zeros(T, d2, d) for _ = 1:n_errors]

    aug_dim = Φ_dim * (1 + n_errors)
    G_func = prop.G
    order = prop.order
    G_errors_local = obj.G_errors

    for k = 1:(N-1)
        zk = _make_knot_view(traj, data, k)
        zk1 = _make_knot_view(traj, data, k + 1)
        u_k = zk[obj.control_name]
        u_k1 = zk1[obj.control_name]
        Δt_k = zk.timestep
        Ũ_k = reshape(zk[obj.state_name], d2, d)

        du_k = nothing
        du_k1 = nothing
        if order == 3 && !isnothing(prop.derivative_name)
            du_k = zk[prop.derivative_name]
            du_k1 = zk1[prop.derivative_name]
        end

        function aug_ode!(dx, x, _p, τ)
            u_τ = _interpolate_controls(τ, u_k, u_k1, Δt_k, order, du_k, du_k1)
            G_τ = G_func(u_τ)

            Φ_mat = reshape(@view(x[1:Φ_dim]), d2, d2)
            dΦ_mat = reshape(@view(dx[1:Φ_dim]), d2, d2)
            mul!(dΦ_mat, G_τ, Φ_mat, Δt_k, false)

            for mm = 1:n_errors
                offset = mm * Φ_dim
                W_mat = reshape(@view(x[(offset+1):(offset+Φ_dim)]), d2, d2)
                dW_mat = reshape(@view(dx[(offset+1):(offset+Φ_dim)]), d2, d2)
                mul!(dW_mat, G_τ, W_mat, Δt_k, false)
                mul!(dW_mat, G_errors_local[mm], Φ_mat, Δt_k, true)
            end
        end

        x0 = zeros(T, aug_dim)
        for i = 1:d2
            x0[(i-1)*d2+i] = one(T)
        end

        prob = ODEProblem(aug_ode!, x0, (zero(T), one(T)))
        sol = solve(
            prob,
            Tsit5();
            abstol = prop.tol,
            reltol = prop.tol,
            save_everystep = false,
        )
        final_state = sol.u[end]

        P_k = reshape(final_state[1:Φ_dim], d2, d2)
        for m = 1:n_errors
            offset = m * Φ_dim
            Φ21_k = reshape(final_state[(offset+1):(offset+Φ_dim)], d2, d2)
            Vs_iso[m] = P_k * Vs_iso[m] + Φ21_k * Ũ_k
        end
    end

    # Convert iso real `Vs_iso` (2d × d) → complex `Vs` (d × d) for the
    # uniform return shape across all _propagate_final_V methods.
    return [_iso_view_to_complex_matrix(V_iso) for V_iso in Vs_iso]
end

# --- Ket version, exponential propagation ---

function _propagate_final_V(
    prop::ExponentialPropagation,
    obj::KetAdjointRobustnessObjective,
    traj::NamedTrajectory,
    data::AbstractVector{T},
) where {T<:Real}
    N = traj.N
    d = obj.ketdim
    n = obj.n_kets
    n_errors = length(obj.G_errors)

    # Complex internals (see ExponentialPropagation note).
    Vs = [zeros(Complex{T}, d, n) for _ = 1:n_errors]
    miEs = [_iso_to_minus_i_H(G_E) for G_E in obj.G_errors]

    for k = 1:(N-1)
        zk = _make_knot_view(traj, data, k)
        u_k = zk[obj.control_name]
        Δt_k = zk.timestep

        # Gather kets into d × n complex matrix Ψ_k_complex
        Ψ_k = Matrix{Complex{T}}(undef, d, n)
        for (i, sname) in enumerate(obj.state_names)
            ψ̃_iso = zk[sname]  # length 2d real
            Ψ_k[:, i] .= view(ψ̃_iso, 1:d) .+ im .* view(ψ̃_iso, (d+1):(2d))
        end

        G_k_iso = prop.G(u_k)
        miH_k = _iso_to_minus_i_H(G_k_iso)

        for m = 1:n_errors
            miE_m = miEs[m]
            aug_C = zeros(Complex{T}, 2d, 2d)
            aug_C[1:d, 1:d] .= miH_k
            aug_C[(d+1):(2d), 1:d] .= miE_m
            aug_C[(d+1):(2d), (d+1):(2d)] .= miH_k

            new_Vm = zeros(Complex{T}, d, n)
            v_aug = Vector{Complex{T}}(undef, 2d)
            for i = 1:n
                v_aug[1:d] .= view(Ψ_k, :, i)
                v_aug[(d+1):(2d)] .= view(Vs[m], :, i)
                result = expv(Δt_k, aug_C, v_aug)
                new_Vm[:, i] .= view(result, (d+1):(2d))
            end
            Vs[m] = new_Vm
        end
    end
    return Vs
end

# --- Pack residual: real vector r so that ‖r‖²/2 == objective_value(obj, traj) ---

function _compute_residual(
    obj::AdjointRobustnessObjective,
    traj::NamedTrajectory,
    data::AbstractVector{T},
) where {T<:Real}
    d = obj.ketdim
    n_errors = length(obj.G_errors)

    zk_N = _make_knot_view(traj, data, traj.N)
    Ũ⃗_N = collect(zk_N[obj.state_name])  # collect to materialize from view
    U_N = iso_vec_to_operator(Ũ⃗_N)  # d × d complex

    # Vs[m] is now d × d complex (was 2d × d real iso) — skip the
    # `iso_vec_to_operator(vec(Vs[m]))` conversion that previously round-tripped
    # through iso form.
    Vs = _propagate_final_V(obj.propagation, obj, traj, data)

    scale = sqrt(2.0 * obj.Q / d)
    r = Vector{T}(undef, 2 * d^2 * n_errors)
    idx = 0
    for m = 1:n_errors
        A = U_N' * Vs[m]
        for j = 1:d, i = 1:d
            idx += 1
            r[idx] = scale * real(A[i, j])
            idx += 1
            r[idx] = scale * imag(A[i, j])
        end
    end
    return r
end

function _compute_residual(
    obj::KetAdjointRobustnessObjective,
    traj::NamedTrajectory,
    data::AbstractVector{T},
) where {T<:Real}
    n = obj.n_kets
    n_errors = length(obj.G_errors)

    # Vs[m] is now d × n complex (was 2d × n real iso) — direct overlap
    # `goals[i]' * Vs[m][:, i]` instead of iso_to_ket conversion per column.
    Vs = _propagate_final_V(obj.propagation, obj, traj, data)

    scale = sqrt(2.0 * obj.Q / n)
    r = Vector{T}(undef, 2 * n * n_errors)
    idx = 0
    for m = 1:n_errors
        V_m = Vs[m]
        for i = 1:n
            o = obj.goals[i]' * view(V_m, :, i)
            idx += 1
            r[idx] = scale * real(o)
            idx += 1
            r[idx] = scale * imag(o)
        end
    end
    return r
end

# Convenience wrappers: read from traj.datavec
_compute_residual(obj::AdjointRobustnessObjective, traj::NamedTrajectory) =
    _compute_residual(obj, traj, traj.datavec)
_compute_residual(obj::KetAdjointRobustnessObjective, traj::NamedTrajectory) =
    _compute_residual(obj, traj, traj.datavec)

# ============================================================================ #
#                          Objective Value                                      #
# ============================================================================ #

"""
    objective_value(obj::AdjointRobustnessObjective, traj::NamedTrajectory)

Compute E_V = Q Σ_m ||U_N† V_N^(m)||_F² / d
"""
function DirectTrajOpt.objective_value(
    obj::AdjointRobustnessObjective,
    traj::NamedTrajectory,
)
    fwd = forward_pass(obj, traj)
    N = traj.N

    # Get U_N (complex operator)
    Ũ⃗_N = traj[N][obj.state_name]
    U_N = iso_vec_to_operator(Ũ⃗_N)

    E_V = 0.0
    for m in eachindex(obj.G_errors)
        V_N = iso_vec_to_operator(vec(fwd.Vs[m][N]))
        A = U_N' * V_N
        E_V += real(tr(A' * A)) / obj.ketdim
    end
    return obj.Q * E_V
end

"""
    objective_value(obj::KetAdjointRobustnessObjective, traj::NamedTrajectory)

Compute E_V = Q Σ_m Σ_i |⟨ψ_goal_i | v_N^(m,i)⟩|² / n_kets
"""
function DirectTrajOpt.objective_value(
    obj::KetAdjointRobustnessObjective,
    traj::NamedTrajectory,
)
    fwd = forward_pass(obj, traj)
    N = traj.N

    E_V = 0.0
    for m in eachindex(obj.G_errors)
        V_N = fwd.Vs[m][N]  # iso_dim × n_kets
        for i = 1:obj.n_kets
            v_i = iso_to_ket(V_N[:, i])
            E_V += abs2(obj.goals[i]' * v_i) / obj.n_kets
        end
    end
    return obj.Q * E_V
end

# ============================================================================ #
#                              Gradient                                         #
# ============================================================================ #

"""
    gradient!(∇::AbstractVector, obj::AdjointRobustnessObjective, traj::NamedTrajectory)

Compute the analytic gradient. Since `f = ½ ‖r‖²` for residual `r` (see
`_compute_residual`), the gradient is `∇f = J^T r` — exactly the residual
VJP applied to `r` itself. This is a thin wrapper over `_residual_vjp!`.

The backward-pass machinery, terminal-seed math, and per-knot control
contributions all live in `_residual_vjp!`. Cost: one forward pass + one
backward pass = O(N).
"""
function DirectTrajOpt.gradient!(
    ∇::AbstractVector,
    obj::AdjointRobustnessObjective,
    traj::NamedTrajectory,
)
    n_data = length(traj.datavec)
    r = _compute_residual(obj, traj)
    if length(∇) == n_data
        _residual_vjp!(∇, obj, traj, r)
    else
        # ∇ may include trailing global_data slot; objective doesn't depend on
        # globals, so those entries stay zero.
        fill!(∇, 0.0)
        ∇_view = view(∇, 1:n_data)
        _residual_vjp!(∇_view, obj, traj, r)
    end
    return nothing
end

"""
    gradient!(∇::AbstractVector, obj::KetAdjointRobustnessObjective, traj::NamedTrajectory)

Ket variant. Same structure: gradient is residual VJP with cotangent `r`.
"""
function DirectTrajOpt.gradient!(
    ∇::AbstractVector,
    obj::KetAdjointRobustnessObjective,
    traj::NamedTrajectory,
)
    n_data = length(traj.datavec)
    r = _compute_residual(obj, traj)
    if length(∇) == n_data
        _residual_vjp!(∇, obj, traj, r)
    else
        fill!(∇, 0.0)
        ∇_view = view(∇, 1:n_data)
        _residual_vjp!(∇_view, obj, traj, r)
    end

    return nothing
end

# ============================================================================ #
#                  JVP / VJP / HVP — Altissimo-aligned primitives               #
# ============================================================================ #
#
# The residual is r ∈ ℝ^{n_res} where n_res = 2 d² n_err (operator) or
# 2 n_kets n_err (ket). The Jacobian J = ∂r/∂Z (Z = traj.datavec) gives
#
#     ∇f(Z) = J^T r            (gradient of f = ½ ‖r‖²)
#     H_GN(Z) = J^T J          (Gauss–Newton Hessian)
#     H_GN(Z) v = J^T (J v)    (matrix-free GN HVP — Altissimo primitive)
#
# We never materialize J for HVP. The pieces:
#   _residual_jvp!(δr, obj, traj, v)  : δr = J v, via ForwardDiff Dual{1}
#   _residual_vjp!(δZ, obj, traj, w)  : δZ = J^T w, via the analytic adjoint
#                                       backward pass with terminal seeds
#                                       built from the cotangent w
#   gradient!(∇, obj, traj)           : VJP with w = r (cf. _residual_vjp!)
#
# `get_full_hessian` materializes H_GN by VJP-on-each-basis-vector of the
# residual (cost: O(n_res) backward passes, vs. O(|Z|) for ForwardDiff over
# the full Jacobian). That's the right choice for n_res ≪ |Z|, which is
# typical (n_res = 2 d² n_err, |Z| = N · dim).

# --- Cotangent → complex W_m unpack (operator) ---

function _cotangent_to_W(obj::AdjointRobustnessObjective, w::AbstractVector{T}) where {T}
    d = obj.ketdim
    n_errors = length(obj.G_errors)
    Ws = Vector{Matrix{Complex{T}}}(undef, n_errors)
    for m = 1:n_errors
        block_offset = (m - 1) * 2 * d^2
        W = Matrix{Complex{T}}(undef, d, d)
        idx = block_offset
        for j = 1:d, i = 1:d
            re = w[idx+1]
            im = w[idx+2]
            W[i, j] = complex(re, im)
            idx += 2
        end
        Ws[m] = W
    end
    return Ws
end

# --- VJP (operator): δZ = J^T w ---

"""
    _residual_vjp!(δZ, obj::AdjointRobustnessObjective, traj, w)

Pull the real residual cotangent `w` back to the trajectory cotangent
`δZ = J^T w`, where `J = ∂r/∂Z` and `r = _compute_residual(obj, traj)`.
Same backward-pass machinery as `gradient!` but with the terminal seeds
built from `w` instead of from the residual itself.

Recovers `gradient!` exactly when `w == r`.
"""
function _residual_vjp!(
    δZ::AbstractVector,
    obj::AdjointRobustnessObjective,
    traj::NamedTrajectory,
    w::AbstractVector,
)
    # Convenience: compute forward pass internally.
    fwd = forward_pass(obj, traj)
    return _residual_vjp!(δZ, obj, traj, fwd, w)
end

function _residual_vjp!(
    δZ::AbstractVector,
    obj::AdjointRobustnessObjective,
    traj::NamedTrajectory,
    fwd::ForwardPassResult,
    w::AbstractVector,
)
    fill!(δZ, 0.0)
    N = traj.N
    d2 = obj.iso_dim
    d = obj.ketdim
    n_errors = length(obj.G_errors)
    scale_r = sqrt(2.0 * obj.Q / d)

    state_comps = collect(traj.components[obj.state_name])
    control_comps = collect(traj.components[obj.control_name])
    Δt_comps = collect(traj.components[traj.timestep])
    z_dim = traj.dim

    Ũ⃗_N = traj[N][obj.state_name]
    U_N = iso_vec_to_operator(Ũ⃗_N)

    Ws = _cotangent_to_W(obj, w)

    for m = 1:n_errors
        V_N_iso = fwd.Vs[m][N]
        V_N = iso_vec_to_operator(vec(V_N_iso))
        W_m = Ws[m]

        # Terminal V cotangent: scale_r · iso_vec(U_N W_m)
        Λ_complex = U_N * W_m
        Λ = reshape(scale_r * operator_to_iso_vec(Λ_complex), d2, d)

        # Terminal U_N gradient contribution: scale_r · iso_vec(V_N W_m^†)
        grad_U_complex = V_N * W_m'
        δZ[slice(N, state_comps, z_dim)] .+= scale_r * operator_to_iso_vec(grad_U_complex)

        # Backward pass — identical to `gradient!` once Λ and ∂U seeded
        for k = (N-1):-1:1
            Ũ⃗_k = traj[k][obj.state_name]
            Ũ_k = reshape(Ũ⃗_k, d2, d)

            ∂Ũ_k = fwd.Φ21s[m][k]' * Λ
            δZ[slice(k, state_comps, z_dim)] .+= vec(∂Ũ_k)

            V_k = fwd.Vs[m][k]
            G_E = obj.G_errors[m]

            _add_control_gradient!(
                obj.propagation,
                δZ,
                Λ,
                V_k,
                Ũ_k,
                obj,
                traj,
                k,
                G_E,
                d2,
                state_comps,
                control_comps,
                Δt_comps,
                z_dim,
            )

            Λ = fwd.Ps[k]' * Λ
        end
    end
    return δZ
end

# --- VJP (ket) ---

function _residual_vjp!(
    δZ::AbstractVector,
    obj::KetAdjointRobustnessObjective,
    traj::NamedTrajectory,
    w::AbstractVector,
)
    fwd = forward_pass(obj, traj)
    return _residual_vjp!(δZ, obj, traj, fwd, w)
end

function _residual_vjp!(
    δZ::AbstractVector,
    obj::KetAdjointRobustnessObjective,
    traj::NamedTrajectory,
    fwd::ForwardPassResult,
    w::AbstractVector,
)
    fill!(δZ, 0.0)
    N = traj.N
    d2 = obj.iso_dim
    n = obj.n_kets
    n_errors = length(obj.G_errors)
    scale_r = sqrt(2.0 * obj.Q / n)

    control_comps = collect(traj.components[obj.control_name])
    Δt_comps = collect(traj.components[traj.timestep])
    z_dim = traj.dim
    all_state_comps = [collect(traj.components[sname]) for sname in obj.state_names]

    for m = 1:n_errors
        # Terminal Λ_N: column i is scale_r · ket_to_iso(z_{m,i} · ψ_goal_i)
        Λ = zeros(d2, n)
        block_offset = (m - 1) * 2 * n
        for i = 1:n
            re = w[block_offset+2i-1]
            im = w[block_offset+2i]
            z = complex(re, im)
            Λ[:, i] = scale_r * ket_to_iso(z * obj.goals[i])
        end

        # Backward pass (no terminal U_N — ket states aren't decision-var
        # endpoints in this objective)
        for k = (N-1):-1:1
            Φ21_Λ = fwd.Φ21s[m][k]' * Λ
            for i = 1:n
                δZ[slice(k, all_state_comps[i], z_dim)] .+= Φ21_Λ[:, i]
            end

            Ψ_k = zeros(d2, n)
            for (i, sname) in enumerate(obj.state_names)
                Ψ_k[:, i] = traj[k][sname]
            end

            V_k = fwd.Vs[m][k]
            G_E = obj.G_errors[m]

            _add_control_gradient!(
                obj.propagation,
                δZ,
                Λ,
                V_k,
                Ψ_k,
                obj,
                traj,
                k,
                G_E,
                d2,
                nothing,
                control_comps,
                Δt_comps,
                z_dim,
            )

            Λ = fwd.Ps[k]' * Λ
        end
    end
    return δZ
end

# --- Batched VJP (operator): process K cotangents in one backward pass ---
#
# The per-knot ForwardDiff.jacobian is cotangent-independent, so it's computed
# ONCE per (knot, error op) and contracted against all K cotangents'
# Λ-matrices via a single matmul. Cost on the Hessian-assembly path:
#
#   unbatched: n_res · (fwd + bwd + n_err · (N-1) · jac)
#   batched:         fwd + bwd + n_err · (N-1) · jac     (× small K-scaling on
#                                                         matmul side)
#
# Savings ≈ n_res = 2d²·n_err. For d=4, n_err=1: ~32×.

function _residual_vjp_batch!(
    δZ_batch::AbstractMatrix{Float64},    # (n_data, K)
    obj::AdjointRobustnessObjective,
    traj::NamedTrajectory,
    fwd::ForwardPassResult,
    W_batch::AbstractMatrix{Float64},     # (n_res, K)
)
    fill!(δZ_batch, 0.0)
    K = size(W_batch, 2)
    N = traj.N
    d2 = obj.iso_dim
    d = obj.ketdim
    n_errors = length(obj.G_errors)
    scale_r = sqrt(2.0 * obj.Q / d)

    state_comps = collect(traj.components[obj.state_name])
    control_comps = collect(traj.components[obj.control_name])
    Δt_comps = collect(traj.components[traj.timestep])
    z_dim = traj.dim

    Ũ⃗_N = traj[N][obj.state_name]
    U_N = iso_vec_to_operator(Ũ⃗_N)

    for m = 1:n_errors
        V_N_iso = fwd.Vs[m][N]
        V_N = iso_vec_to_operator(vec(V_N_iso))

        # Build batched terminal seeds: for each cotangent k_cot, extract
        # complex W_{m,k_cot} from W_batch's m-block and form Λ and ∂U_N.
        # Λ_batch layout: (d2, d*K) — each (d2, d) slice per cotangent.
        Λ_batch = zeros(d2, d * K)
        ∂UN_batch = zeros(d2 * d, K)
        block_offset = (m - 1) * 2 * d^2
        for k_cot = 1:K
            W_mk = Matrix{ComplexF64}(undef, d, d)
            idx = 0
            for j = 1:d, i = 1:d
                W_mk[i, j] = complex(
                    W_batch[block_offset+idx+1, k_cot],
                    W_batch[block_offset+idx+2, k_cot],
                )
                idx += 2
            end
            Λ_k = scale_r * operator_to_iso_vec(U_N * W_mk)
            Λ_batch[:, ((k_cot-1)*d+1):(k_cot*d)] .= reshape(Λ_k, d2, d)
            ∂UN_batch[:, k_cot] .= scale_r * operator_to_iso_vec(V_N * W_mk')
        end

        # Accumulate terminal U_N contribution
        δZ_batch[slice(N, state_comps, z_dim), :] .+= ∂UN_batch

        # Backward pass — batched
        for k = (N-1):-1:1
            Ũ⃗_k = traj[k][obj.state_name]
            Ũ_k = reshape(Ũ⃗_k, d2, d)

            # ∂Ũ_k contribution: Φ21^T * Λ_batch — one BLAS-3 matmul for all K
            ∂Ũ_k_batch = fwd.Φ21s[m][k]' * Λ_batch   # (d2, d*K)
            # Reshape each column block to vec form and accumulate per-cotangent
            ∂Ũ_k_vecd = reshape(∂Ũ_k_batch, d2 * d, K)
            δZ_batch[slice(k, state_comps, z_dim), :] .+= ∂Ũ_k_vecd

            # Control gradient — batched
            V_k = fwd.Vs[m][k]
            G_E = obj.G_errors[m]
            _add_control_gradient_batched!(
                obj.propagation,
                δZ_batch,
                Λ_batch,
                V_k,
                Ũ_k,
                obj,
                traj,
                k,
                G_E,
                d2,
                d,
                K,
                state_comps,
                control_comps,
                Δt_comps,
                z_dim,
            )

            # Backward propagate Λ_batch — one BLAS-3 matmul
            Λ_batch = fwd.Ps[k]' * Λ_batch
        end
    end
    return δZ_batch
end

# --- JVP: δr = J v via ForwardDiff Dual{1} ---

"""
    _residual_jvp!(δr, obj, traj, v)

Push the real trajectory tangent `v` forward through the residual to get
`δr = J v`. Implemented via `ForwardDiff.derivative` along the line
`datavec + t * v` — one Dual{1} pass through the forward propagation.
"""
function _residual_jvp!(
    δr::AbstractVector,
    obj::Union{AdjointRobustnessObjective,KetAdjointRobustnessObjective},
    traj::NamedTrajectory,
    v::AbstractVector,
)
    res = ForwardDiff.derivative(
        t -> _compute_residual(obj, traj, traj.datavec .+ t .* v),
        zero(eltype(traj.datavec)),
    )
    copyto!(δr, res)
    return δr
end

# --- Objective HVP: H_GN v = J^T (J v) — Altissimo's `compute_objective_hvp!` ---

"""
    objective_hvp!(Hv, obj, traj, v)

Matrix-free Gauss–Newton Hessian-vector product:

    H_GN v = J^T (J v) = vjp(jvp(v))

Cost: one JVP (forward Dual{1}) + one VJP (analytic adjoint backward) —
roughly two gradient-evaluation costs, **independent of |Z|**. Suitable
for plugging directly into Altissimo's `compute_objective_hvp!` slot.

Note: drops the second-order residual term `Σ_k r_k ∇²r_k` (vanishes at
the optimum, but is non-zero away from it — see `gauss_newton_explanation.md`
for the NonlinearDrive failure mode).
"""
function objective_hvp!(
    Hv::AbstractVector,
    obj::Union{AdjointRobustnessObjective,KetAdjointRobustnessObjective},
    traj::NamedTrajectory,
    v::AbstractVector,
)
    n_res = length(_compute_residual(obj, traj))
    δr = Vector{Float64}(undef, n_res)
    _residual_jvp!(δr, obj, traj, v)
    _residual_vjp!(Hv, obj, traj, δr)
    return Hv
end

# ============================================================================ #
#                  Batched per-knot control gradient                            #
# ============================================================================ #
#
# Batched counterparts of `_add_control_gradient!` — compute the per-knot
# ForwardDiff.jacobian ONCE and contract against all K cotangents' Λ-matrices
# via one BLAS-3 matmul. Drop-in replacement for the scalar loop inside
# `_residual_vjp_batch!`.

function _add_control_gradient_batched!(
    prop::ExponentialPropagation,
    δZ_batch::AbstractMatrix{Float64},    # (n_data, K)
    Λ_batch::AbstractMatrix{Float64},     # (d2, d*K)
    V_k::Matrix{Float64},
    Ũ_k::AbstractMatrix,
    obj::AbstractObjective,
    traj::NamedTrajectory,
    k::Int,
    G_E::Matrix{Float64},
    d2::Int,
    d::Int,
    K::Int,
    state_comps,
    control_comps,
    Δt_comps,
    z_dim,
)
    u_k = traj[k][obj.control_name]
    Δt_k = traj[k].timestep
    u_indices = slice(k, control_comps, z_dim)
    Δt_index = slice(k, Δt_comps, z_dim)
    G_func = prop.G
    n_u = length(u_k)

    # Same ForwardDiff.jacobian closure as the scalar path — shared across all
    # K cotangents because it depends only on (u_k, Δt_k), not on Λ.
    function block_exp_action(params)
        u = @view params[1:n_u]
        Δt = params[n_u+1]
        G_k = G_func(u)
        TT = eltype(params)
        aug_G = zeros(TT, 2d2, 2d2)
        aug_G[1:d2, 1:d2] .= G_k
        aug_G[(d2+1):2d2, 1:d2] .= G_E
        aug_G[(d2+1):2d2, (d2+1):2d2] .= G_k

        out = Vector{TT}(undef, d2 * d)
        v_aug = Vector{TT}(undef, 2d2)
        for j = 1:d
            v_aug[1:d2] .= view(Ũ_k, :, j)
            v_aug[(d2+1):2d2] .= view(V_k, :, j)
            result = expv(Δt, aug_G, v_aug)
            out[((j-1)*d2+1):(j*d2)] .= view(result, (d2+1):(2d2))
        end
        return out
    end

    params = [u_k; Δt_k]
    J = ForwardDiff.jacobian(block_exp_action, params)  # (d2·d, n_u+1)

    # Contract J^T against Λ_batch for all K cotangents in one matmul.
    # Λ_batch (d2, d*K) reshaped to (d2·d, K); grad_batch = (n_u+1, K).
    Λ_flat = reshape(Λ_batch, d2 * d, K)
    grad_batch = J' * Λ_flat  # (n_u+1, K)

    δZ_batch[u_indices, :] .+= view(grad_batch, 1:n_u, :)
    δZ_batch[Δt_index, :] .+= view(grad_batch, (n_u+1):(n_u+1), :)

    return nothing
end

function _add_control_gradient_batched!(
    prop::ZKickPropagation,
    δZ_batch::AbstractMatrix{Float64},
    Λ_batch::AbstractMatrix{Float64},
    V_k::Matrix{Float64},
    Ũ_k::AbstractMatrix,
    obj::AbstractObjective,
    traj::NamedTrajectory,
    k::Int,
    G_E::Matrix{Float64},
    d2::Int,
    d::Int,
    K::Int,
    state_comps,
    control_comps,
    Δt_comps,
    z_dim,
)
    ϕ_k = traj[k][obj.control_name][1]
    Δt_k = traj[k].timestep
    u_indices = slice(k, control_comps, z_dim)
    Δt_index = slice(k, Δt_comps, z_dim)
    G_H = prop.G_H

    function zkick_block_action(params)
        ϕ = params[1]
        Δt = params[2]
        R = _iso_Rz(ϕ)
        TT = eltype(params)

        aug_G = zeros(TT, 2d2, 2d2)
        aug_G[1:d2, 1:d2] .= G_H
        aug_G[(d2+1):2d2, 1:d2] .= G_E
        aug_G[(d2+1):2d2, (d2+1):2d2] .= G_H

        out = Vector{TT}(undef, d2 * d)
        v_aug = Vector{TT}(undef, 2d2)
        for j = 1:d
            v_aug[1:d2] .= view(Ũ_k, :, j)
            v_aug[(d2+1):2d2] .= view(V_k, :, j)
            result = expv(Δt, aug_G, v_aug)
            out[((j-1)*d2+1):(j*d2)] .= R * view(result, (d2+1):(2d2))
        end
        return out
    end

    J = ForwardDiff.jacobian(zkick_block_action, [ϕ_k, Δt_k])
    Λ_flat = reshape(Λ_batch, d2 * d, K)
    grad_batch = J' * Λ_flat  # (2, K)

    δZ_batch[u_indices, :] .+= view(grad_batch, 1:1, :)
    δZ_batch[Δt_index, :] .+= view(grad_batch, 2:2, :)

    return nothing
end

function _add_control_gradient_batched!(
    prop::SplineODEPropagation,
    δZ_batch::AbstractMatrix{Float64},
    Λ_batch::AbstractMatrix{Float64},
    V_k::Matrix{Float64},
    Ũ_k::AbstractMatrix,
    obj::AdjointRobustnessObjective,
    traj::NamedTrajectory,
    k::Int,
    G_E::Matrix{Float64},
    d2::Int,
    d::Int,
    K::Int,
    state_comps,
    control_comps,
    Δt_comps,
    z_dim,
)
    N = traj.N
    u_k = traj[k][obj.control_name]
    u_k1 = traj[k+1][obj.control_name]
    Δt_k = traj[k].timestep
    n_u = length(u_k)
    Φ_dim = d2^2
    G_func = prop.G
    order = prop.order
    tol = prop.tol

    du_k = nothing
    du_k1 = nothing
    if order == 3 && !isnothing(prop.derivative_name)
        du_k = traj[k][prop.derivative_name]
        du_k1 = traj[k+1][prop.derivative_name]
    end

    function ode_action(params)
        if order == 1
            uu_k = params[1:n_u]
            uu_k1 = params[(n_u+1):(2n_u)]
            ΔΔt = params[2n_u+1]
            ddu_k = nothing
            ddu_k1 = nothing
        else  # order == 3
            uu_k = params[1:n_u]
            uu_k1 = params[(n_u+1):(2n_u)]
            ddu_k = params[(2n_u+1):(3n_u)]
            ddu_k1 = params[(3n_u+1):(4n_u)]
            ΔΔt = params[4n_u+1]
        end

        TT = eltype(params)
        aug_dim = 2 * Φ_dim
        function aug_ode!(dx, x, _p, τ)
            u_τ = _interpolate_controls(τ, uu_k, uu_k1, ΔΔt, order, ddu_k, ddu_k1)
            G_τ = G_func(u_τ)

            Φ_mat = reshape(@view(x[1:Φ_dim]), d2, d2)
            dΦ_mat = reshape(@view(dx[1:Φ_dim]), d2, d2)
            W_mat = reshape(@view(x[(Φ_dim+1):(2Φ_dim)]), d2, d2)
            dW_mat = reshape(@view(dx[(Φ_dim+1):(2Φ_dim)]), d2, d2)

            mul!(dΦ_mat, G_τ, Φ_mat, ΔΔt, false)
            mul!(dW_mat, G_τ, W_mat, ΔΔt, false)
            mul!(dW_mat, G_E, Φ_mat, ΔΔt, true)
        end

        x0 = zeros(TT, aug_dim)
        for i = 1:d2
            x0[(i-1)*d2+i] = one(TT)
        end

        prob = ODEProblem(aug_ode!, x0, (zero(TT), one(TT)))
        sol = solve(prob, Tsit5(); abstol = tol, reltol = tol, save_everystep = false)
        fs = sol.u[end]

        P = reshape(fs[1:Φ_dim], d2, d2)
        Ψ = reshape(fs[(Φ_dim+1):(2Φ_dim)], d2, d2)
        return vec(P * V_k + Ψ * Ũ_k)
    end

    if order == 1
        params = [u_k; u_k1; Δt_k]
    else
        params = [u_k; u_k1; du_k; du_k1; Δt_k]
    end

    J = ForwardDiff.jacobian(ode_action, params)
    Λ_flat = reshape(Λ_batch, d2 * d, K)
    grad_batch = J' * Λ_flat  # (n_params, K)

    u_indices_k = slice(k, control_comps, z_dim)
    u_indices_k1 = slice(k + 1, control_comps, z_dim)
    Δt_index = slice(k, Δt_comps, z_dim)

    δZ_batch[u_indices_k, :] .+= view(grad_batch, 1:n_u, :)
    δZ_batch[u_indices_k1, :] .+= view(grad_batch, (n_u+1):(2n_u), :)

    if order == 3 && !isnothing(prop.derivative_name)
        du_comps = collect(traj.components[prop.derivative_name])
        du_indices_k = slice(k, du_comps, z_dim)
        du_indices_k1 = slice(k + 1, du_comps, z_dim)
        δZ_batch[du_indices_k, :] .+= view(grad_batch, (2n_u+1):(3n_u), :)
        δZ_batch[du_indices_k1, :] .+= view(grad_batch, (3n_u+1):(4n_u), :)
        δZ_batch[Δt_index, :] .+= view(grad_batch, (4n_u+1):(4n_u+1), :)
    else
        δZ_batch[Δt_index, :] .+= view(grad_batch, (2n_u+1):(2n_u+1), :)
    end

    return nothing
end

# ============================================================================ #
#                          Per-knot control gradient                            #
# ============================================================================ #

"""
    _add_control_gradient!(prop, ∇, Λ, V_k, Ũ_k, obj, traj, k, G_E, d2, ...)

Add gradient contributions w.r.t. controls and timestep at knot k.
Dispatches on propagation type.
"""
function _add_control_gradient!(
    prop::ExponentialPropagation,
    ∇::AbstractVector,
    Λ::Matrix{Float64},
    V_k::Matrix{Float64},
    Ũ_k::AbstractMatrix,
    obj::AbstractObjective,
    traj::NamedTrajectory,
    k::Int,
    G_E::Matrix{Float64},
    d2::Int,
    state_comps,
    control_comps,
    Δt_comps,
    z_dim,
)
    u_k = traj[k][obj.control_name]
    Δt_k = traj[k].timestep
    u_indices = slice(k, control_comps, z_dim)
    Δt_index = slice(k, Δt_comps, z_dim)
    G_func = prop.G
    n_u = length(u_k)
    d = size(V_k, 2)

    # Action of the augmented exp on `[Ũ_col; V_col]`: bottom block is
    # `Φ21*Ũ_col + P*V_col`, which is exactly the new V column. Use `expv`
    # (Krylov, Dual-friendly) instead of forming exp(Δt*aug_G) — Julia 1.12's
    # `LinearAlgebra.exp!` only handles BlasFloats so it breaks under
    # ForwardDiff. ForwardDiff over expv gives an analytic per-knot Jacobian.
    function block_exp_action(params)
        u = @view params[1:n_u]
        Δt = params[n_u+1]
        G_k = G_func(u)
        TT = eltype(params)
        aug_G = zeros(TT, 2d2, 2d2)
        aug_G[1:d2, 1:d2] .= G_k
        aug_G[(d2+1):2d2, 1:d2] .= G_E
        aug_G[(d2+1):2d2, (d2+1):2d2] .= G_k

        out = Vector{TT}(undef, d2 * d)
        v_aug = Vector{TT}(undef, 2d2)
        for j = 1:d
            v_aug[1:d2] .= view(Ũ_k, :, j)
            v_aug[(d2+1):2d2] .= view(V_k, :, j)
            result = expv(Δt, aug_G, v_aug)
            out[((j-1)*d2+1):(j*d2)] .= view(result, (d2+1):2d2)
        end
        return out
    end

    params = [u_k; Δt_k]
    J = ForwardDiff.jacobian(block_exp_action, params)

    Λ_vec = vec(Λ)
    grad_params = J' * Λ_vec

    ∇[u_indices] .+= grad_params[1:n_u]
    ∇[Δt_index] .+= grad_params[n_u+1]

    return nothing
end

function _add_control_gradient!(
    prop::ZKickPropagation,
    ∇::AbstractVector,
    Λ::Matrix{Float64},
    V_k::Matrix{Float64},
    Ũ_k::AbstractMatrix,
    obj::AbstractObjective,
    traj::NamedTrajectory,
    k::Int,
    G_E::Matrix{Float64},
    d2::Int,
    state_comps,
    control_comps,
    Δt_comps,
    z_dim,
)
    ϕ_k = traj[k][obj.control_name][1]
    Δt_k = traj[k].timestep
    u_indices = slice(k, control_comps, z_dim)
    Δt_index = slice(k, Δt_comps, z_dim)
    G_H = prop.G_H

    d = size(V_k, 2)

    # Same expv-columnwise pattern as ExponentialPropagation, plus the Z-kick
    # rotation R_k applied to the bottom-block result of the free-evolution
    # action. ForwardDiff over expv gives an analytic per-knot Jacobian.
    function zkick_block_action(params)
        ϕ = params[1]
        Δt = params[2]
        R = _iso_Rz(ϕ)
        TT = eltype(params)

        aug_G = zeros(TT, 2d2, 2d2)
        aug_G[1:d2, 1:d2] .= G_H
        aug_G[(d2+1):2d2, 1:d2] .= G_E
        aug_G[(d2+1):2d2, (d2+1):2d2] .= G_H

        out = Vector{TT}(undef, d2 * d)
        v_aug = Vector{TT}(undef, 2d2)
        for j = 1:d
            v_aug[1:d2] .= view(Ũ_k, :, j)
            v_aug[(d2+1):2d2] .= view(V_k, :, j)
            result = expv(Δt, aug_G, v_aug)
            out[((j-1)*d2+1):(j*d2)] .= R * view(result, (d2+1):2d2)
        end
        return out
    end

    J = ForwardDiff.jacobian(zkick_block_action, [ϕ_k, Δt_k])
    Λ_vec = vec(Λ)
    grad_params = J' * Λ_vec

    ∇[u_indices] .+= grad_params[1:1]
    ∇[Δt_index] .+= grad_params[2:2]

    return nothing
end

function _add_control_gradient!(
    prop::SplineODEPropagation,
    ∇::AbstractVector,
    Λ::Matrix{Float64},
    V_k::Matrix{Float64},
    Ũ_k::AbstractMatrix,
    obj::AdjointRobustnessObjective,
    traj::NamedTrajectory,
    k::Int,
    G_E::Matrix{Float64},
    d2::Int,
    state_comps,
    control_comps,
    Δt_comps,
    z_dim,
)
    N = traj.N
    u_k = traj[k][obj.control_name]
    u_k1 = traj[k+1][obj.control_name]
    Δt_k = traj[k].timestep
    n_u = length(u_k)
    Φ_dim = d2^2
    G_func = prop.G
    order = prop.order
    tol = prop.tol

    # Get derivative components for cubic spline
    du_k = nothing
    du_k1 = nothing
    if order == 3 && !isnothing(prop.derivative_name)
        du_k = traj[k][prop.derivative_name]
        du_k1 = traj[k+1][prop.derivative_name]
    end

    # Build parameter vector: [u_k, u_{k+1}, (du_k, du_{k+1} for cubic), Δt]
    # ForwardDiff propagates Duals through the augmented ODE solve (Tsit5
    # supports Dual element types), giving an analytic per-knot Jacobian.
    function ode_action(params)
        if order == 1
            uu_k = params[1:n_u]
            uu_k1 = params[(n_u+1):2n_u]
            ΔΔt = params[2n_u+1]
            ddu_k = nothing
            ddu_k1 = nothing
        else  # order == 3
            uu_k = params[1:n_u]
            uu_k1 = params[(n_u+1):2n_u]
            ddu_k = params[(2n_u+1):3n_u]
            ddu_k1 = params[(3n_u+1):4n_u]
            ΔΔt = params[4n_u+1]
        end

        TT = eltype(params)
        aug_dim = 2 * Φ_dim
        function aug_ode!(dx, x, _p, τ)
            u_τ = _interpolate_controls(τ, uu_k, uu_k1, ΔΔt, order, ddu_k, ddu_k1)
            G_τ = G_func(u_τ)

            Φ_mat = reshape(@view(x[1:Φ_dim]), d2, d2)
            dΦ_mat = reshape(@view(dx[1:Φ_dim]), d2, d2)
            W_mat = reshape(@view(x[(Φ_dim+1):2Φ_dim]), d2, d2)
            dW_mat = reshape(@view(dx[(Φ_dim+1):2Φ_dim]), d2, d2)

            mul!(dΦ_mat, G_τ, Φ_mat, ΔΔt, false)
            mul!(dW_mat, G_τ, W_mat, ΔΔt, false)
            mul!(dW_mat, G_E, Φ_mat, ΔΔt, true)
        end

        x0 = zeros(TT, aug_dim)
        for i = 1:d2
            x0[(i-1)*d2+i] = one(TT)
        end

        prob = ODEProblem(aug_ode!, x0, (zero(TT), one(TT)))
        sol = solve(prob, Tsit5(); abstol = tol, reltol = tol, save_everystep = false)
        fs = sol.u[end]

        P = reshape(fs[1:Φ_dim], d2, d2)
        Ψ = reshape(fs[(Φ_dim+1):2Φ_dim], d2, d2)
        return vec(P * V_k + Ψ * Ũ_k)
    end

    if order == 1
        params = [u_k; u_k1; Δt_k]
    else
        params = [u_k; u_k1; du_k; du_k1; Δt_k]
    end

    J = ForwardDiff.jacobian(ode_action, params)
    Λ_vec = vec(Λ)
    grad_params = J' * Λ_vec

    # Distribute gradients to correct indices
    u_indices_k = slice(k, control_comps, z_dim)
    u_indices_k1 = slice(k + 1, control_comps, z_dim)
    Δt_index = slice(k, Δt_comps, z_dim)

    ∇[u_indices_k] .+= grad_params[1:n_u]
    ∇[u_indices_k1] .+= grad_params[(n_u+1):2n_u]

    if order == 3 && !isnothing(prop.derivative_name)
        du_comps = collect(traj.components[prop.derivative_name])
        du_indices_k = slice(k, du_comps, z_dim)
        du_indices_k1 = slice(k + 1, du_comps, z_dim)
        ∇[du_indices_k] .+= grad_params[(2n_u+1):3n_u]
        ∇[du_indices_k1] .+= grad_params[(3n_u+1):4n_u]
        ∇[Δt_index] .+= grad_params[4n_u+1]
    else
        ∇[Δt_index] .+= grad_params[2n_u+1]
    end

    return nothing
end

# ============================================================================ #
#                          Hessian Structure                                    #
# ============================================================================ #

"""
    hessian_structure(obj::AdjointRobustnessObjective, traj::NamedTrajectory)

Return sparsity structure of the Hessian. The adjoint robustness objective couples
all state and control variables across all knot points (dense structure).
"""
function DirectTrajOpt.Objectives.hessian_structure(
    obj::AdjointRobustnessObjective,
    traj::NamedTrajectory,
)
    Z_dim = traj.dim * traj.N + traj.global_dim
    structure = spzeros(Z_dim, Z_dim)

    state_comps = collect(traj.components[obj.state_name])
    control_comps = collect(traj.components[obj.control_name])
    Δt_comps = collect(traj.components[traj.timestep])
    z_dim = traj.dim
    N = traj.N

    # Collect all relevant indices
    all_indices = Int[]
    for k = 1:N
        append!(all_indices, slice(k, state_comps, z_dim))
    end
    for k = 1:(N-1)
        append!(all_indices, slice(k, control_comps, z_dim))
        append!(all_indices, slice(k, Δt_comps, z_dim))
    end
    # Include derivative components for spline propagation
    if obj.propagation isa SplineODEPropagation &&
       !isnothing(obj.propagation.derivative_name)
        du_comps = collect(traj.components[obj.propagation.derivative_name])
        for k = 1:N
            append!(all_indices, slice(k, du_comps, z_dim))
        end
    end

    # Mark all cross-terms (dense block). Use COO construction — assigning
    # element-by-element to a SparseMatrixCSC is O(nnz) per element, which
    # collapsed ~1000²=1M assignments into ~8 minutes on the CZ demo.
    n = length(all_indices)
    I_idx = repeat(all_indices; outer = n)
    J_idx = repeat(all_indices; inner = n)
    V_idx = ones(Float64, n * n)
    structure = sparse(I_idx, J_idx, V_idx, Z_dim, Z_dim)

    return structure
end

# ============================================================================ #
#                          Gauss-Newton Hessian                                #
# ============================================================================ #

"""
    get_full_hessian(obj::AdjointRobustnessObjective, traj::NamedTrajectory)

Analytic Gauss–Newton approximation of the Hessian: H ≈ J^T J

where J = ∂r/∂Z is the Jacobian of the residual

    r_m = √(2Q/d) · [vec(Re(U_N† V_N^(m))); vec(Im(U_N† V_N^(m)))]

stacked over all error operators m. Since f = ½ ‖r‖², we have ∇²f = J^T J +
Σ_k r_k ∇² r_k; GN drops the second term, which vanishes at the optimum.

J is computed analytically via ForwardDiff through the residual function. Cost:
O(|Z|) Dual-passes through the forward propagation, vs. O(|Z|²) for the previous
finite-difference Hessian.

Note: GN is an approximation away from the optimum. For NonlinearDrive systems
the dropped second-order term has a constant component (∂²H/∂u² ≠ 0); see
[gauss_newton_explanation.md](https://github.com/harmoniqs/Piccolissimo.jl/blob/main/gauss_newton_explanation.md).
"""
function DirectTrajOpt.Objectives.get_full_hessian(
    obj::AdjointRobustnessObjective,
    traj::NamedTrajectory,
)
    if obj.exact_hessian
        return _assemble_exact_hessian_via_fd_gradient(obj, traj)
    else
        # Build J row-by-row via VJP on each residual basis vector — cost is
        # O(n_res · backward-pass) ≪ O(|Z| · forward-pass) for typical n_res ≪ |Z|.
        # Then H_GN = J^T J. Matches the legacy `length(traj.datavec)`² shape
        # (objective is independent of global_data, so global block of H_GN is 0).
        return _assemble_gn_hessian_via_vjp(obj, traj)
    end
end

function _assemble_gn_hessian_via_vjp(
    obj::AdjointRobustnessObjective,
    traj::NamedTrajectory,
)
    n_data = length(traj.datavec)
    n_res = length(_compute_residual(obj, traj))
    Z_dim = traj.dim * traj.N + traj.global_dim

    fwd = forward_pass(obj, traj)

    # Batched VJP: n_res basis cotangents in one backward pass. `J_T[:, j]` is
    # row j of J (= gradient of residual component j). Then H_GN = J^T J =
    # J_T * J_T^T.
    J_T = zeros(n_data, n_res)
    W_batch = Matrix{Float64}(I, n_res, n_res)
    _residual_vjp_batch!(J_T, obj, traj, fwd, W_batch)
    H_data = J_T * J_T'  # (n_data × n_data)
    # Embed into (Z_dim × Z_dim). Objective does not depend on global_data,
    # so the global block is zero. CompositeObjective.get_full_hessian expects
    # this Z_dim shape to add with peers.
    if Z_dim == n_data
        return sparse(Matrix(H_data))
    else
        H_full = spzeros(Z_dim, Z_dim)
        H_full[1:n_data, 1:n_data] = H_data
        return H_full
    end
end

# Ket path still uses the unbatched per-cotangent VJP because the ket
# `_residual_vjp!` takes a different internal shape — batching the ket variant
# is tracked as a follow-up. Cost: cached forward pass is reused across VJPs,
# so each VJP pays only backward-pass cost.
function _assemble_gn_hessian_via_vjp(
    obj::KetAdjointRobustnessObjective,
    traj::NamedTrajectory,
)
    n_data = length(traj.datavec)
    n_res = length(_compute_residual(obj, traj))
    Z_dim = traj.dim * traj.N + traj.global_dim
    fwd = forward_pass(obj, traj)
    J = zeros(n_res, n_data)
    e = zeros(n_res)
    row_buf = Vector{Float64}(undef, n_data)
    for j = 1:n_res
        e[j] = 1.0
        _residual_vjp!(row_buf, obj, traj, fwd, e)
        view(J, j, :) .= row_buf
        e[j] = 0.0
    end
    H_data = J' * J
    if Z_dim == n_data
        return sparse(Matrix(H_data))
    else
        H_full = spzeros(Z_dim, Z_dim)
        H_full[1:n_data, 1:n_data] = H_data
        return H_full
    end
end

# Exact Hessian fallback: H = ∂(∇f)/∂Z via FD on `gradient!`. Use this when
# `exact_hessian=true` (typically because the system has NonlinearDrives, where
# the GN approximation drops a non-vanishing curvature term — see
# `gauss_newton_explanation.md`).
#
# Cost: O(|Z| · gradient-cost) — much more expensive than GN's O(n_res · gradient).
# Symmetrized to remove FD asymmetry noise.
function _assemble_exact_hessian_via_fd_gradient(
    obj::Union{AdjointRobustnessObjective,KetAdjointRobustnessObjective},
    traj::NamedTrajectory,
)
    n_data = length(traj.datavec)
    Z_dim = traj.dim * traj.N + traj.global_dim
    Z = collect(traj.datavec)
    g_buf = Vector{Float64}(undef, n_data)

    grad_at(z) = begin
        traj_z = NamedTrajectory(traj; datavec = z)
        gradient!(g_buf, obj, traj_z)
        return copy(view(g_buf, 1:n_data))
    end

    H = FiniteDiff.finite_difference_jacobian(grad_at, Z)
    H_sym = (H + H') / 2
    if Z_dim == n_data
        return sparse(Matrix(H_sym))
    else
        H_full = spzeros(Z_dim, Z_dim)
        H_full[1:n_data, 1:n_data] = H_sym
        return H_full
    end
end

# --- Ket hessian structure ---

function DirectTrajOpt.Objectives.hessian_structure(
    obj::KetAdjointRobustnessObjective,
    traj::NamedTrajectory,
)
    Z_dim = traj.dim * traj.N + traj.global_dim
    structure = spzeros(Z_dim, Z_dim)

    control_comps = collect(traj.components[obj.control_name])
    Δt_comps = collect(traj.components[traj.timestep])
    z_dim = traj.dim
    N = traj.N

    all_indices = Int[]
    for sname in obj.state_names
        scomps = collect(traj.components[sname])
        for k = 1:N
            append!(all_indices, slice(k, scomps, z_dim))
        end
    end
    for k = 1:(N-1)
        append!(all_indices, slice(k, control_comps, z_dim))
        append!(all_indices, slice(k, Δt_comps, z_dim))
    end

    # Dense block via COO construction (see operator variant for why).
    n = length(all_indices)
    I_idx = repeat(all_indices; outer = n)
    J_idx = repeat(all_indices; inner = n)
    V_idx = ones(Float64, n * n)
    structure = sparse(I_idx, J_idx, V_idx, Z_dim, Z_dim)

    return structure
end

"""
    get_full_hessian(obj::KetAdjointRobustnessObjective, traj::NamedTrajectory)

Analytic Gauss–Newton Hessian for the ket variant; same structure as the operator
case with residual

    r_{m,i} = √(2Q/n) · [Re ⟨ψ_goal_i|v_i^(m)⟩; Im ⟨ψ_goal_i|v_i^(m)⟩]
"""
function DirectTrajOpt.Objectives.get_full_hessian(
    obj::KetAdjointRobustnessObjective,
    traj::NamedTrajectory,
)
    if obj.exact_hessian
        return _assemble_exact_hessian_via_fd_gradient(obj, traj)
    else
        return _assemble_gn_hessian_via_vjp(obj, traj)
    end
end

# ============================================================================ #
#                       RobustControlProblem Wrapper                            #
# ============================================================================ #

"""
    _extract_regularizers(objective, state_sym::Symbol)

Extract regularizer terms from a (possibly composite) objective, dropping
any terms that depend on the state variable (i.e., infidelity objectives).
"""
function _extract_regularizers(objective, state_sym::Symbol)
    objs = hasproperty(objective, :objectives) ? objective.objectives : [objective]

    regularizers = filter(objs) do term
        term_syms = if hasproperty(term, :syms)
            term.syms
        elseif hasproperty(term, :var_names)
            term.var_names
        elseif hasproperty(term, :name) && term.name isa Symbol
            (term.name,)
        else
            Symbol[]
        end
        state_sym ∉ term_syms
    end

    return isempty(regularizers) ? nothing : reduce(+, regularizers)
end

"""
    _rob_fidelity_constraint(qtraj, final_fidelity, traj)

Create a fidelity constraint dispatched on trajectory type.
"""
function _rob_fidelity_constraint(qtraj::UnitaryTrajectory, final_fidelity, traj)
    U_goal = qtraj.goal
    state_sym = piccolo_state_name(qtraj)
    return FinalUnitaryFidelityConstraint(U_goal, state_sym, final_fidelity, traj)
end

function _rob_fidelity_constraint(qtraj::KetTrajectory, final_fidelity, traj)
    ψ_goal = qtraj.goal
    state_sym = piccolo_state_name(qtraj)
    return FinalKetFidelityConstraint(ψ_goal, state_sym, final_fidelity, traj)
end

function _rob_fidelity_constraint(qtraj::MultiKetTrajectory, final_fidelity, traj)
    goals = qtraj.goals
    snames = piccolo_state_names(qtraj)
    return FinalCoherentKetFidelityConstraint(goals, snames, final_fidelity, traj)
end

"""
    _find_dynamics_integrator(integrators)

Find the dynamics integrator from a list of integrators.
Returns the first HermitianExponentialIntegrator or SplineIntegrator found.
"""
function _find_dynamics_integrator(integrators)
    for ℰ in integrators
        if ℰ isa HermitianExponentialIntegrator
            return ℰ
        elseif ℰ isa SplineIntegrator
            return ℰ
        end
    end
    error("No supported dynamics integrator found for robustness objective")
end

"""
    _build_robustness_objective(integrator, qtraj, error_operators, traj, sys, Q)

Build the appropriate robustness objective based on integrator type.
"""
function _build_robustness_objective(
    integrator::HermitianExponentialIntegrator{UnitaryTrajectory},
    qtraj,
    error_operators,
    traj,
    sys,
    Q,
)
    return AdjointRobustnessObjective(integrator, error_operators, traj; Q = Q)
end

function _build_robustness_objective(
    integrator::HermitianExponentialIntegrator{MultiKetTrajectory},
    qtraj::MultiKetTrajectory,
    error_operators,
    traj,
    sys,
    Q,
)
    return KetAdjointRobustnessObjective(
        integrator,
        error_operators,
        qtraj.goals,
        traj;
        Q = Q,
    )
end

function _build_robustness_objective(
    integrator::HermitianExponentialIntegrator{KetTrajectory},
    qtraj::KetTrajectory,
    error_operators,
    traj,
    sys,
    Q,
)
    return KetAdjointRobustnessObjective(
        integrator,
        error_operators,
        qtraj.goal,
        traj;
        Q = Q,
    )
end

function _build_robustness_objective(
    integrator::SplineIntegrator,
    qtraj,
    error_operators,
    traj,
    sys,
    Q,
)
    @assert !isnothing(sys) "QuantumSystem required for SplineIntegrator robustness"
    return AdjointRobustnessObjective(integrator, sys, error_operators, traj; Q = Q)
end

"""
    _extract_regularizers_multi(objective, state_syms::Vector{Symbol})

Extract regularizer terms from a (possibly composite) objective, dropping
any terms that depend on ANY of the state variables.
"""
function _extract_regularizers_multi(objective, state_syms::AbstractVector{Symbol})
    objs = hasproperty(objective, :objectives) ? objective.objectives : [objective]

    regularizers = filter(objs) do term
        term_syms = if hasproperty(term, :syms)
            term.syms
        elseif hasproperty(term, :var_names)
            term.var_names
        elseif hasproperty(term, :name) && term.name isa Symbol
            (term.name,)
        else
            Symbol[]
        end
        all(s ∉ term_syms for s in state_syms)
    end

    return isempty(regularizers) ? nothing : reduce(+, regularizers)
end

"""
    RobustControlProblem(
        qcp::QuantumControlProblem;
        error_operators,
        sys = nothing,
        final_fidelity = 0.9999,
        Q_robustness = 1.0,
        keep_infidelity_objective = false,
    )

Wrap an existing `QuantumControlProblem` to minimize error susceptibility
subject to a fidelity floor constraint. Follows the `MinimumTimeProblem`
composition pattern.

# Arguments
- `qcp`: Base problem (from `SmoothPulseProblem`, `SplinePulseProblem`, etc.)
- `error_operators`: Vector of Hermitian error matrices [E1, E2, ...]
- `sys`: QuantumSystem (required for SplineIntegrator, optional for exponential)
- `final_fidelity`: Minimum fidelity constraint (default 0.9999)
- `Q_robustness`: Weight on the robustness objective
- `keep_infidelity_objective`: If true, keep infidelity in objective alongside robustness
"""
function RobustControlProblem(
    qcp::QuantumControlProblem{<:Any,QT};
    error_operators::Vector{<:AbstractMatrix} = Matrix{ComplexF64}[],
    sys = nothing,
    final_fidelity::Float64 = 0.9999,
    Q_robustness::Float64 = 1.0,
    keep_infidelity_objective::Bool = false,
    robustness_objective::Union{Nothing,AbstractObjective} = nothing,
) where {QT<:AbstractQuantumTrajectory}
    # 1. Deep-copy trajectory and constraints
    traj = deepcopy(qcp.prob.trajectory)
    constraints = deepcopy(qcp.prob.constraints)

    # 2. Build robustness objective — use provided one or auto-construct from integrator
    rob_obj = if !isnothing(robustness_objective)
        robustness_objective
    else
        integrator = _find_dynamics_integrator(qcp.prob.integrators)
        _build_robustness_objective(
            integrator,
            qcp.qtraj,
            error_operators,
            traj,
            sys,
            Q_robustness,
        )
    end

    # 3. Compose objectives
    if keep_infidelity_objective
        J = qcp.prob.objective + rob_obj
    else
        # For MultiKet/Ket, need to filter out objectives referencing ANY state name
        state_syms = if qcp.qtraj isa MultiKetTrajectory
            piccolo_state_names(qcp.qtraj)
        else
            [piccolo_state_name(qcp.qtraj)]
        end
        reg = _extract_regularizers_multi(qcp.prob.objective, state_syms)
        J = isnothing(reg) ? rob_obj : reg + rob_obj
    end

    # 4. Add fidelity floor constraint
    fidelity_constraint = _rob_fidelity_constraint(qcp.qtraj, final_fidelity, traj)
    if fidelity_constraint isa AbstractVector
        append!(constraints, fidelity_constraint)
    else
        push!(constraints, fidelity_constraint)
    end

    # 5. Build new problem, reusing integrators
    # Use positional (inner) constructor to avoid re-adding trajectory constraints
    # (deepcopied constraints already include them from the original problem)
    new_prob = DirectTrajOptProblem(traj, J, qcp.prob.integrators, constraints)
    return QuantumControlProblem(deepcopy(qcp.qtraj), new_prob)
end

# ============================================================================ #
#                            Test Helper                                        #
# ============================================================================ #

"""
    test_robustness_objective(obj, traj; atol=1e-9, rtol=1e-9, show_diff=false)

Validate an `AdjointRobustnessObjective` or `KetAdjointRobustnessObjective`:

1. `objective_value` returns a real scalar.
2. ‖r‖²/2 ≈ objective value (residual packing is correct).
3. `gradient!` matches `ForwardDiff.gradient` of the objective `½‖r‖²`.
4. `get_full_hessian` (analytic GN, `J_AD^T J_AD`) matches the independent assembly
   `J_ref^T J_ref` where `J_ref = ForwardDiff.jacobian(_compute_residual)`. This
   validates the GN implementation **without** asserting that GN equals the exact
   Hessian (which is only true at the optimum — see `gauss_newton_explanation.md`).

The reference is ForwardDiff (machine-precision), so these comparisons hold to
~1e-9 rather than FiniteDiff's ~1e-4.
5. The GN Hessian is symmetric and PSD (modulo numerical noise).

Replaces `DirectTrajOpt.Objectives.test_objective` for these objectives, which
assumed the analytic Hessian equals the FD-of-objective Hessian (it does not
under the GN approximation, except at the optimum).
"""
function test_robustness_objective(
    obj::Union{AdjointRobustnessObjective,KetAdjointRobustnessObjective},
    traj::NamedTrajectory;
    atol = 1e-9,
    rtol = 1e-9,
    show_diff = false,
)
    n_data = length(traj.datavec)

    # 1. Scalar objective
    f = DirectTrajOpt.objective_value(obj, traj)
    Test.@test f isa Real

    # 2. ‖r‖²/2 == f (residual packing matches the objective)
    r = _compute_residual(obj, traj)
    Test.@test isapprox(0.5 * sum(abs2, r), f; rtol = 1e-12, atol = 1e-14)

    # 3. Gradient vs ForwardDiff of ½‖r‖² (over the data part of Z; objective is
    # independent of global_data, so we only test that block here). We
    # differentiate ½‖r‖² rather than `objective_value` directly: the two are
    # equal (asserted in step 2) and the residual path is the Dual-generic one.
    Z⃗_vec = collect(vec(traj))
    Z_dim = length(Z⃗_vec)
    ∇ = zeros(Z_dim)
    DirectTrajOpt.gradient!(∇, obj, traj)
    ∇_ad = ForwardDiff.gradient(traj.datavec) do td
        return 0.5 * sum(abs2, _compute_residual(obj, traj, td))
    end
    if show_diff
        for i = 1:n_data
            if abs(∇[i] - ∇_ad[i]) > atol ||
               abs((∇[i] - ∇_ad[i]) / (∇_ad[i] + 1e-10)) > rtol
                println("  ∇[$i]: analytic=$(∇[i]) ad=$(∇_ad[i])")
            end
        end
    end
    Test.@test all(isapprox.(view(∇, 1:n_data), ∇_ad; atol = atol, rtol = rtol))

    # 4. GN Hessian matches independently-assembled J_ref^T J_ref, where the
    # residual Jacobian J_ref is computed by ForwardDiff (machine precision).
    H_gn = Matrix(DirectTrajOpt.Objectives.get_full_hessian(obj, traj))
    J_ref = ForwardDiff.jacobian(traj.datavec) do td
        return _compute_residual(obj, traj, td)
    end
    H_check = J_ref' * J_ref
    if show_diff
        println("  Hessian max abs diff: $(maximum(abs, H_gn - H_check))")
    end
    Test.@test isapprox(H_gn, H_check; atol = atol, rtol = rtol)

    # 5. Symmetry + PSD
    Test.@test isapprox(H_gn, H_gn'; atol = 1e-10)
    λmin = minimum(real, eigvals(Hermitian((H_gn + H_gn') / 2)))
    Test.@test λmin > -1e-8

    return nothing
end

# ============================================================================ #
#                              Tests                                            #
# ============================================================================ #

@testitem "AdjointRobustnessObjective exact_hessian flag" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra
    using ForwardDiff
    import Piccolo.Control.AdjointRobustness: _compute_residual

    # Linear-drive system — GN and exact Hessian both well-defined; exact may
    # still differ from GN away from optimum (because of the residual-times-
    # second-order-sensitivity term, which doesn't vanish for general U_N).
    T = 5.0
    N = 6
    sys = OpenQuantumSystem(GATES.Z / 2, [GATES.X / 2, GATES.Y / 2], [1.0, 1.0])
    qtraj = UnitaryTrajectory(sys, GATES.H, T)
    traj = NamedTrajectory(qtraj, N)
    integrator = HermitianExponentialIntegrator(qtraj, N)
    E_ops = [Matrix{ComplexF64}(GATES.Z)]

    rob_gn = AdjointRobustnessObjective(integrator, E_ops, traj; Q = 1.0)
    rob_exact =
        AdjointRobustnessObjective(integrator, E_ops, traj; Q = 1.0, exact_hessian = true)

    @test rob_gn.exact_hessian == false
    @test rob_exact.exact_hessian == true

    H_gn = Matrix(get_full_hessian(rob_gn, traj))
    H_exact = Matrix(get_full_hessian(rob_exact, traj))

    # Both symmetric
    @test isapprox(H_gn, H_gn'; atol = 1e-10)
    @test isapprox(H_exact, H_exact'; atol = 1e-8)

    # GN is PSD by construction; exact can have tiny negative eigenvalues from
    # FD noise away from optimum, but should be approximately PSD here
    @test minimum(real, eigvals(Hermitian((H_gn + H_gn') / 2))) > -1e-8

    # GN and exact agree at the **objective-value-zero** trajectory (V_N is
    # identically 0 when error operators are zero, so the dropped term r·∇²r is
    # 0 trivially). Here, test the non-trivial case: the exact Hessian must match
    # the TRUE exact Hessian of ½‖r‖², assembled by nested-Dual ForwardDiff on the
    # residual (= Jᵀ J + Σ_k r_k ∇²r_k). The exact_hessian impl is FiniteDiff over
    # `gradient!`, so the residual is the impl's FD error (~3e-8), not the
    # reference — far tighter than the old FD-vs-FD 1e-3.
    n_data = length(traj.datavec)
    H_ad_obj = ForwardDiff.hessian(traj.datavec) do z
        return 0.5 * sum(abs2, _compute_residual(rob_gn, traj, z))
    end
    @test isapprox(H_exact, H_ad_obj; atol = 1e-6, rtol = 1e-6)

    println("✓ exact_hessian flag test passed")
end

@testitem "AdjointRobustnessObjective JVP / VJP / HVP primitives" begin
    using Piccolo
    using Piccolo.Control.AdjointRobustness:
        _residual_jvp!, _residual_vjp!, _compute_residual
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra
    using ForwardDiff

    # Small operator problem
    T = 5.0
    N = 6
    sys = OpenQuantumSystem(GATES.Z / 2, [GATES.X / 2, GATES.Y / 2], [1.0, 1.0])
    qtraj = UnitaryTrajectory(sys, GATES.H, T)
    traj = NamedTrajectory(qtraj, N)
    integrator = HermitianExponentialIntegrator(qtraj, N)
    E_ops = [Matrix{ComplexF64}(GATES.Z), Matrix{ComplexF64}(GATES.X)]
    rob_obj = AdjointRobustnessObjective(integrator, E_ops, traj; Q = 1.0)

    n_data = length(traj.datavec)
    r = _compute_residual(rob_obj, traj)
    n_res = length(r)

    # Independent ForwardDiff Jacobian of the residual (machine precision). The
    # analytic VJP (and hence gradient/HVP) is validated against J_refᵀ, so this
    # is a genuine cross-check even though the JVP itself is built from Duals.
    J_ref = ForwardDiff.jacobian(z -> _compute_residual(rob_obj, traj, z), traj.datavec)

    # --- JVP: δr = J v ---
    v = randn(n_data)
    δr = zeros(n_res)
    _residual_jvp!(δr, rob_obj, traj, v)
    @test isapprox(δr, J_ref * v; atol = 1e-10, rtol = 1e-10)

    # JVP linearity
    α = 1.7
    δr_α = zeros(n_res)
    _residual_jvp!(δr_α, rob_obj, traj, α .* v)
    @test isapprox(δr_α, α .* δr; atol = 1e-12)

    # --- VJP: δZ = J^T w ---
    w = randn(n_res)
    δZ = zeros(n_data)
    _residual_vjp!(δZ, rob_obj, traj, w)
    @test isapprox(δZ, J_ref' * w; atol = 1e-9, rtol = 1e-9)

    # VJP linearity
    β = -2.3
    δZ_β = zeros(n_data)
    _residual_vjp!(δZ_β, rob_obj, traj, β .* w)
    @test isapprox(δZ_β, β .* δZ; atol = 1e-10)

    # --- VJP with w = r recovers gradient! ---
    ∇_grad = zeros(traj.dim * traj.N + traj.global_dim)
    gradient!(∇_grad, rob_obj, traj)
    ∇_vjp = zeros(n_data)
    _residual_vjp!(∇_vjp, rob_obj, traj, r)
    @test isapprox(view(∇_grad, 1:n_data), ∇_vjp; atol = 1e-8)

    # --- HVP: H_GN v = J^T (J v) = vjp(jvp(v)) ; also vs J_refᵀ(J_ref v) ---
    Hv = zeros(n_data)
    objective_hvp!(Hv, rob_obj, traj, v)
    H_gn = Matrix(get_full_hessian(rob_obj, traj))
    @test isapprox(Hv, H_gn * v; atol = 1e-9, rtol = 1e-9)
    @test isapprox(Hv, J_ref' * (J_ref * v); atol = 1e-9, rtol = 1e-9)

    # HVP symmetry: vᵀ H w' == w'ᵀ H v for v, w' both in input space
    w_input = randn(n_data)
    Hw = zeros(n_data)
    objective_hvp!(Hw, rob_obj, traj, w_input)
    @test isapprox(dot(v, Hw), dot(w_input, Hv); atol = 1e-9, rtol = 1e-9)

    println("✓ JVP / VJP / HVP primitives test passed")
end

@testitem "AdjointRobustnessObjective basic test" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    # Single qubit system: H = ω Z + u_x X
    T = 5.0
    N = 10
    sys = OpenQuantumSystem(GATES.Z / 2, [GATES.X / 2], [1.0])

    # Create unitary trajectory targeting Hadamard gate
    qtraj = UnitaryTrajectory(sys, GATES.H, T)
    traj = NamedTrajectory(qtraj, N)

    # Create integrator
    integrator = HermitianExponentialIntegrator(qtraj, N)

    # Error operator: Z dephasing
    E_ops = [Matrix{ComplexF64}(GATES.Z)]

    # Create objective
    rob_obj = AdjointRobustnessObjective(integrator, E_ops, traj; Q = 1.0)

    # Test objective value is real and non-negative
    J = objective_value(rob_obj, traj)
    @test J isa Float64
    @test J >= 0.0

    # Test using DirectTrajOpt's test_objective (gradient + Hessian vs finite diff)
    test_robustness_objective(rob_obj, traj; atol = 1e-9, rtol = 1e-9)

    println("✓ AdjointRobustnessObjective basic test passed")
end

@testitem "AdjointRobustnessObjective multiple error operators" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    # Single qubit with X and Y drives
    T = 5.0
    N = 8
    sys = OpenQuantumSystem(GATES.Z / 2, [GATES.X / 2, GATES.Y / 2], [1.0, 1.0])

    qtraj = UnitaryTrajectory(sys, GATES.H, T)
    traj = NamedTrajectory(qtraj, N)
    integrator = HermitianExponentialIntegrator(qtraj, N)

    # Multiple error operators: X and Z
    E_ops = [Matrix{ComplexF64}(GATES.X), Matrix{ComplexF64}(GATES.Z)]

    rob_obj = AdjointRobustnessObjective(integrator, E_ops, traj; Q = 1.0)

    J = objective_value(rob_obj, traj)
    @test J >= 0.0

    # Verify multi-error matches sum of single-error objectives
    rob_X = AdjointRobustnessObjective(integrator, [E_ops[1]], traj; Q = 1.0)
    rob_Z = AdjointRobustnessObjective(integrator, [E_ops[2]], traj; Q = 1.0)
    @test objective_value(rob_obj, traj) ≈
          objective_value(rob_X, traj) + objective_value(rob_Z, traj)

    test_robustness_objective(rob_obj, traj; atol = 1e-9, rtol = 1e-9)

    println("✓ AdjointRobustnessObjective multiple error operators test passed")
end

@testitem "AdjointRobustnessObjective V=0 for zero error" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    T = 5.0
    N = 6
    sys = OpenQuantumSystem(GATES.Z / 2, [GATES.X / 2], [1.0])

    qtraj = UnitaryTrajectory(sys, GATES.H, T)
    traj = NamedTrajectory(qtraj, N)
    integrator = HermitianExponentialIntegrator(qtraj, N)

    # Zero error operator → E_V should be 0
    E_zero = [zeros(ComplexF64, 2, 2)]
    rob_obj = AdjointRobustnessObjective(integrator, E_zero, traj; Q = 1.0)
    @test objective_value(rob_obj, traj) ≈ 0.0 atol = 1e-12

    println("✓ AdjointRobustnessObjective V=0 for zero error test passed")
end

# ============================================================================ #
# Q and error amplitude scaling
# ============================================================================ #

@testitem "AdjointRobustnessObjective Q and error scaling" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    T = 5.0
    N = 8
    sys = OpenQuantumSystem(GATES.Z / 2, [GATES.X / 2], [1.0])
    qtraj = UnitaryTrajectory(sys, GATES.H, T)
    traj = NamedTrajectory(qtraj, N)
    integrator = HermitianExponentialIntegrator(qtraj, N)
    E_ops = [Matrix{ComplexF64}(GATES.Z)]

    # Q scaling: J(Q=α) = α * J(Q=1)
    rob_Q1 = AdjointRobustnessObjective(integrator, E_ops, traj; Q = 1.0)
    rob_Q3 = AdjointRobustnessObjective(integrator, E_ops, traj; Q = 3.0)
    rob_Q05 = AdjointRobustnessObjective(integrator, E_ops, traj; Q = 0.5)

    J1 = objective_value(rob_Q1, traj)
    @test objective_value(rob_Q3, traj) ≈ 3.0 * J1 rtol = 1e-12
    @test objective_value(rob_Q05, traj) ≈ 0.5 * J1 rtol = 1e-12

    # Error amplitude scaling: J(αE) = α² * J(E) (quadratic in error)
    α = 2.5
    E_scaled = [α * Matrix{ComplexF64}(GATES.Z)]
    rob_scaled = AdjointRobustnessObjective(integrator, E_scaled, traj; Q = 1.0)
    @test objective_value(rob_scaled, traj) ≈ α^2 * J1 rtol = 1e-10

    # Scaling with small α
    β = 0.1
    E_small = [β * Matrix{ComplexF64}(GATES.Z)]
    rob_small = AdjointRobustnessObjective(integrator, E_small, traj; Q = 1.0)
    @test objective_value(rob_small, traj) ≈ β^2 * J1 rtol = 1e-10

    println("✓ AdjointRobustnessObjective Q and error scaling test passed")
end

# ============================================================================ #
# Gradient properties: non-zero, descent direction, tight FD validation
# ============================================================================ #

@testitem "AdjointRobustnessObjective gradient properties" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    T = 5.0
    N = 8
    sys = OpenQuantumSystem(GATES.Z / 2, [GATES.X / 2, GATES.Y / 2], [1.0, 1.0])
    qtraj = UnitaryTrajectory(sys, GATES.H, T)
    traj = NamedTrajectory(qtraj, N)
    integrator = HermitianExponentialIntegrator(qtraj, N)
    E_ops = [Matrix{ComplexF64}(GATES.Z)]
    rob_obj = AdjointRobustnessObjective(integrator, E_ops, traj; Q = 1.0)

    # Gradient is non-zero for non-trivial case
    Z_dim = traj.dim * traj.N + traj.global_dim
    ∇ = zeros(Z_dim)
    gradient!(∇, rob_obj, traj)
    @test norm(∇) > 0.0

    # Small step in -∇ direction decreases objective (validates sign)
    J0 = objective_value(rob_obj, traj)
    ε = 1e-7
    Z = collect(vec(traj))
    Z_new = Z - ε * ∇
    traj_new = NamedTrajectory(traj; datavec = Z_new)
    J_new = objective_value(rob_obj, traj_new)
    @test J_new < J0

    # Tight tolerance gradient/Hessian validation (2-drive system)
    test_robustness_objective(rob_obj, traj; atol = 1e-9, rtol = 1e-9)

    println("✓ AdjointRobustnessObjective gradient properties test passed")
end

# ============================================================================ #
# Forward pass validation: recurrence, orthogonality, consistency
# ============================================================================ #

@testitem "AdjointRobustnessObjective forward pass validation" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    T = 5.0
    N = 8
    sys = OpenQuantumSystem(GATES.Z / 2, [GATES.X / 2], [1.0])
    qtraj = UnitaryTrajectory(sys, GATES.H, T)
    traj = NamedTrajectory(qtraj, N)
    integrator = HermitianExponentialIntegrator(qtraj, N)
    E_ops = [Matrix{ComplexF64}(GATES.X), Matrix{ComplexF64}(GATES.Z)]
    rob_obj = AdjointRobustnessObjective(integrator, E_ops, traj; Q = 1.0)

    fwd = Piccolo.Control.AdjointRobustness.forward_pass(rob_obj, traj)

    # V_1 = 0 for each error operator
    for m in eachindex(E_ops)
        @test fwd.Vs[m][1] ≈ zeros(rob_obj.iso_dim, rob_obj.ketdim) atol = 1e-15
    end

    # P_k matrices preserve orthogonality (iso form of unitary → P^T P ≈ I)
    for k = 1:(N-1)
        @test fwd.Ps[k]' * fwd.Ps[k] ≈ I(rob_obj.iso_dim) atol = 1e-10
        @test fwd.Ps[k] * fwd.Ps[k]' ≈ I(rob_obj.iso_dim) atol = 1e-10
    end

    # Total time matches sum of Δt_k
    dt_sum = sum(traj[k].timestep for k = 1:(N-1))
    @test fwd.t_f ≈ dt_sum

    # Recurrence holds: V_{k+1} = P_k V_k + Φ21_k Ũ_k
    for m in eachindex(E_ops)
        for k = 1:(N-1)
            Ũ⃗_k = traj[k][rob_obj.state_name]
            Ũ_k = reshape(Ũ⃗_k, rob_obj.iso_dim, rob_obj.ketdim)
            expected = fwd.Ps[k] * fwd.Vs[m][k] + fwd.Φ21s[m][k] * Ũ_k
            @test fwd.Vs[m][k+1] ≈ expected atol = 1e-10
        end
    end

    # V_N from forward pass gives consistent objective value
    U_N = iso_vec_to_operator(traj[N][rob_obj.state_name])
    E_V_manual = let ev = 0.0
        for m in eachindex(E_ops)
            V_N = iso_vec_to_operator(vec(fwd.Vs[m][N]))
            A = U_N' * V_N
            ev += real(tr(A' * A)) / rob_obj.ketdim
        end
        ev
    end
    @test rob_obj.Q * E_V_manual ≈ objective_value(rob_obj, traj) rtol = 1e-12

    println("✓ AdjointRobustnessObjective forward pass validation test passed")
end

# ============================================================================ #
# Analytical identity evolution: exact closed-form E_V = Q * tr(E†E) / d
# ============================================================================ #

@testitem "AdjointRobustnessObjective analytical identity evolution" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    # System with zero drift: H(u) = u * σ_x
    # At u=0: G(0) = 0, P_k = I, Φ21_k = Δt_k * G_E
    # V_N = T * G_E, so ||U†V||² = T² * ||E||²
    # Exact result: E_V = Q * T² * tr(E†E) / d
    sys = OpenQuantumSystem([GATES.X], [1.0])
    d = 2

    N = 5
    T = 2.0
    dt = T / (N - 1)

    # Manual trajectory: all unitaries = identity, zero controls
    Ũ⃗_I = operator_to_iso_vec(Matrix{ComplexF64}(I, d, d))
    Ũ⃗_data = repeat(Ũ⃗_I, 1, N)
    u_data = zeros(1, N)
    Δt_data = fill(dt, 1, N)

    traj = NamedTrajectory(
        (Ũ⃗ = Ũ⃗_data, u = u_data, Δt = Δt_data);
        timestep = :Δt,
        controls = (:u, :Δt),
    )

    # E = σ_z: T²*tr(σ_z†σ_z)/d = 4*2/2 = 4.0
    E_z = [Matrix{ComplexF64}(GATES.Z)]
    rob_z = AdjointRobustnessObjective(sys, E_z, :Ũ⃗, :u, traj; Q = 1.0)
    @test objective_value(rob_z, traj) ≈ T^2 * 1.0 atol = 1e-8

    # E = σ_x: same (also unitary), E_V = T² * 1.0
    E_x = [Matrix{ComplexF64}(GATES.X)]
    rob_x = AdjointRobustnessObjective(sys, E_x, :Ũ⃗, :u, traj; Q = 1.0)
    @test objective_value(rob_x, traj) ≈ T^2 * 1.0 atol = 1e-8

    # E = α * σ_z: E_V = T² * α² * tr(σ_z†σ_z)/d = T² * α²
    α = 0.7
    E_scaled = [α * Matrix{ComplexF64}(GATES.Z)]
    rob_scaled = AdjointRobustnessObjective(sys, E_scaled, :Ũ⃗, :u, traj; Q = 1.0)
    @test objective_value(rob_scaled, traj) ≈ T^2 * α^2 atol = 1e-8

    # Q scaling on analytical result
    rob_Q5 = AdjointRobustnessObjective(sys, E_z, :Ũ⃗, :u, traj; Q = 5.0)
    @test objective_value(rob_Q5, traj) ≈ 5.0 * T^2 atol = 1e-8

    # Multiple errors: X + Z → T² * (1.0 + 1.0) = T² * 2.0
    E_both = [Matrix{ComplexF64}(GATES.X), Matrix{ComplexF64}(GATES.Z)]
    rob_both = AdjointRobustnessObjective(sys, E_both, :Ũ⃗, :u, traj; Q = 1.0)
    @test objective_value(rob_both, traj) ≈ T^2 * 2.0 atol = 1e-8

    # Different N (more knot points): same analytical answer
    N2 = 20
    dt2 = T / (N2 - 1)
    traj2 = NamedTrajectory(
        (Ũ⃗ = repeat(Ũ⃗_I, 1, N2), u = zeros(1, N2), Δt = fill(dt2, 1, N2));
        timestep = :Δt,
        controls = (:u, :Δt),
    )
    rob_z2 = AdjointRobustnessObjective(sys, E_z, :Ũ⃗, :u, traj2; Q = 1.0)
    @test objective_value(rob_z2, traj2) ≈ T^2 * 1.0 atol = 1e-8

    # Gradient/Hessian validation on analytical trajectory
    test_robustness_objective(rob_z, traj; atol = 1e-9, rtol = 1e-9)

    println("✓ AdjointRobustnessObjective analytical identity evolution test passed")
end

# ============================================================================ #
# 3-level system
# ============================================================================ #

@testitem "AdjointRobustnessObjective 3-level system" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    # 3-level system with two drives
    H_drift = diagm(ComplexF64[1.0, 0.0, -1.0])
    H_drive1 = zeros(ComplexF64, 3, 3)
    H_drive1[1, 2] = H_drive1[2, 1] = 1.0
    H_drive2 = zeros(ComplexF64, 3, 3)
    H_drive2[2, 3] = H_drive2[3, 2] = 1.0

    sys = OpenQuantumSystem(H_drift, [H_drive1, H_drive2], [1.0, 1.0])

    T = 3.0
    N = 8
    goal = Matrix{ComplexF64}(I, 3, 3)  # Identity gate
    qtraj = UnitaryTrajectory(sys, goal, T)
    traj = NamedTrajectory(qtraj, N)
    integrator = HermitianExponentialIntegrator(qtraj, N)

    # Error: dephasing
    E_dephase = diagm(ComplexF64[1.0, 0.0, -1.0])
    E_ops = [E_dephase]

    rob_obj = AdjointRobustnessObjective(integrator, E_ops, traj; Q = 1.0)

    # Basic type and value checks
    J = objective_value(rob_obj, traj)
    @test J isa Float64
    @test J >= 0.0

    # Dimensions are correct for d=3
    @test rob_obj.ketdim == 3
    @test rob_obj.iso_dim == 6
    @test rob_obj.state_dim == 18

    # Gradient is non-zero
    Z_dim = traj.dim * traj.N + traj.global_dim
    ∇ = zeros(Z_dim)
    gradient!(∇, rob_obj, traj)
    @test norm(∇) > 0.0

    # Forward pass V_1 = 0
    fwd = Piccolo.Control.AdjointRobustness.forward_pass(rob_obj, traj)
    @test fwd.Vs[1][1] ≈ zeros(6, 3) atol = 1e-15

    # P_k orthogonality for d=3
    for k = 1:(N-1)
        @test fwd.Ps[k]' * fwd.Ps[k] ≈ I(6) atol = 1e-10
    end

    # Gradient/Hessian validation
    test_robustness_objective(rob_obj, traj; atol = 1e-9, rtol = 1e-9)

    println("✓ AdjointRobustnessObjective 3-level system test passed")
end

# ============================================================================ #
# Composite objective: combine with QuadraticRegularizer
# ============================================================================ #

@testitem "AdjointRobustnessObjective composite objective" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    T = 5.0
    N = 8
    sys = OpenQuantumSystem(GATES.Z / 2, [GATES.X / 2], [1.0])
    qtraj = UnitaryTrajectory(sys, GATES.H, T)
    traj = NamedTrajectory(qtraj, N)
    integrator = HermitianExponentialIntegrator(qtraj, N)

    E_ops = [Matrix{ComplexF64}(GATES.Z)]
    rob_obj = AdjointRobustnessObjective(integrator, E_ops, traj; Q = 1.0)

    # Combine with QuadraticRegularizer via + operator
    reg_obj = QuadraticRegularizer(:u, traj, 0.1)
    composite = rob_obj + reg_obj

    # Composite value = sum of individual values
    J_rob = objective_value(rob_obj, traj)
    J_reg = objective_value(reg_obj, traj)
    @test objective_value(composite, traj) ≈ J_rob + J_reg

    # Scaled composite
    scaled = 2.0 * rob_obj
    @test objective_value(scaled, traj) ≈ 2.0 * J_rob

    # Verify composite gradient = sum of individual gradients
    Z_dim = traj.dim * traj.N + traj.global_dim
    ∇_composite = zeros(Z_dim)
    ∇_rob = zeros(Z_dim)
    ∇_reg = zeros(Z_dim)
    gradient!(∇_composite, composite, traj)
    gradient!(∇_rob, rob_obj, traj)
    gradient!(∇_reg, reg_obj, traj)
    @test ∇_composite ≈ ∇_rob + ∇_reg

    println("✓ AdjointRobustnessObjective composite objective test passed")
end

# ============================================================================ #
# Integration test: solve a robust Hadamard gate (paper Sec. V)
# ============================================================================ #

# ============================================================================ #
# SplineIntegrator tests
# ============================================================================ #

@testitem "AdjointRobustnessObjective SplineIntegrator linear" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    T = 1.0
    N = 6
    sys = OpenQuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])

    times = collect(range(0.0, T, length = N))
    pulse = LinearSplinePulse(fill(0.5, 2, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, GATES.H)
    traj = NamedTrajectory(qtraj, N)

    integrator = SplineIntegrator(qtraj, N)

    E_ops = [Matrix{ComplexF64}(GATES.Z)]
    rob_obj = AdjointRobustnessObjective(integrator, sys, E_ops, traj; Q = 1.0)

    J = objective_value(rob_obj, traj)
    @test J isa Float64
    @test J >= 0.0

    # SplineIntegrator uses an adaptive Tsit5 ODE solve: ForwardDiff propagates
    # Duals through the solver, which selects slightly different steps than the
    # primal, capping analytic-vs-AD agreement near the solver tolerance (~1e-7)
    # rather than machine precision. 1e-6 is still 100× tighter than the prior
    # FiniteDiff bound; the analytic math is exact (HermitianExponential ≈ 1e-15).
    test_robustness_objective(rob_obj, traj; atol = 1e-6, rtol = 1e-6)

    println("✓ AdjointRobustnessObjective SplineIntegrator linear test passed")
end

@testitem "AdjointRobustnessObjective SplineIntegrator cubic" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    T = 1.0
    N = 6
    sys = OpenQuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])

    times = collect(range(0.0, T, length = N))
    pulse = CubicSplinePulse(fill(0.5, 2, N), fill(0.0, 2, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, GATES.H)
    traj = NamedTrajectory(qtraj, N)

    integrator = SplineIntegrator(qtraj, N; spline_order = 3)

    E_ops = [Matrix{ComplexF64}(GATES.Z)]
    rob_obj = AdjointRobustnessObjective(integrator, sys, E_ops, traj; Q = 1.0)

    J = objective_value(rob_obj, traj)
    @test J isa Float64
    @test J >= 0.0

    # Adaptive Tsit5 ODE solve — see the SplineIntegrator-linear test for why
    # ForwardDiff agreement caps near solver tolerance (~1e-7), not 1e-9.
    test_robustness_objective(rob_obj, traj; atol = 1e-6, rtol = 1e-6)

    println("✓ AdjointRobustnessObjective SplineIntegrator cubic test passed")
end

@testitem "AdjointRobustnessObjective SplineIntegrator V=0 for zero error" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    T = 1.0
    N = 5
    sys = OpenQuantumSystem(GATES.Z, [GATES.X], [1.0])

    times = collect(range(0.0, T, length = N))
    pulse = LinearSplinePulse(fill(0.5, 1, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, GATES.H)
    traj = NamedTrajectory(qtraj, N)

    integrator = SplineIntegrator(qtraj, N)

    E_zero = [zeros(ComplexF64, 2, 2)]
    rob_obj = AdjointRobustnessObjective(integrator, sys, E_zero, traj; Q = 1.0)
    @test objective_value(rob_obj, traj) ≈ 0.0 atol = 1e-10

    println("✓ AdjointRobustnessObjective SplineIntegrator V=0 for zero error test passed")
end

@testitem "AdjointRobustnessObjective SplineIntegrator multiple errors" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    T = 1.0
    N = 6
    sys = OpenQuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])

    times = collect(range(0.0, T, length = N))
    pulse = LinearSplinePulse(fill(0.5, 2, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, GATES.H)
    traj = NamedTrajectory(qtraj, N)

    integrator = SplineIntegrator(qtraj, N)

    E_X = [Matrix{ComplexF64}(GATES.X)]
    E_Z = [Matrix{ComplexF64}(GATES.Z)]
    E_both = [Matrix{ComplexF64}(GATES.X), Matrix{ComplexF64}(GATES.Z)]

    rob_X = AdjointRobustnessObjective(integrator, sys, E_X, traj; Q = 1.0)
    rob_Z = AdjointRobustnessObjective(integrator, sys, E_Z, traj; Q = 1.0)
    rob_both = AdjointRobustnessObjective(integrator, sys, E_both, traj; Q = 1.0)

    @test objective_value(rob_both, traj) ≈
          objective_value(rob_X, traj) + objective_value(rob_Z, traj)

    # Adaptive Tsit5 ODE solve — see the SplineIntegrator-linear test for why
    # ForwardDiff agreement caps near solver tolerance (~1e-7), not 1e-9.
    test_robustness_objective(rob_both, traj; atol = 1e-6, rtol = 1e-6)

    println("✓ AdjointRobustnessObjective SplineIntegrator multiple errors test passed")
end

# ============================================================================ #
# VariationalQuantumSystem test
# ============================================================================ #

@testitem "AdjointRobustnessObjective VariationalQuantumSystem" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    T = 5.0
    N = 8

    # Regular system for trajectory creation
    sys = OpenQuantumSystem(GATES.Z / 2, [GATES.X / 2, GATES.Y / 2], [1.0, 1.0])

    # Variational system: same drift/drives, Z dephasing as variational perturbation
    varsys = VariationalQuantumSystem(
        GATES.Z / 2,
        [GATES.X / 2, GATES.Y / 2],
        [GATES.Z],
        [1.0, 1.0],
    )

    qtraj = UnitaryTrajectory(sys, GATES.H, T)
    traj = NamedTrajectory(qtraj, N)
    integrator = HermitianExponentialIntegrator(qtraj, N)

    # Explicit error operator constructor
    E_ops = [Matrix{ComplexF64}(GATES.Z)]
    rob_explicit = AdjointRobustnessObjective(integrator, E_ops, traj; Q = 1.0)

    # VariationalQuantumSystem constructor (uses G_vars as error operators)
    rob_varsys = AdjointRobustnessObjective(
        varsys,
        rob_explicit.state_name,
        rob_explicit.control_name,
        traj;
        Q = 1.0,
    )

    # Both should produce the same objective value
    @test objective_value(rob_varsys, traj) ≈ objective_value(rob_explicit, traj)

    # Gradient check
    test_robustness_objective(rob_varsys, traj; atol = 1e-9, rtol = 1e-9)

    println("✓ AdjointRobustnessObjective VariationalQuantumSystem test passed")
end

# ============================================================================ #
# RobustControlProblem wrapper tests
# ============================================================================ #

@testitem "RobustControlProblem wrapper" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    T = 10.0
    N = 15
    sys = OpenQuantumSystem(GATES.Z / 2, [GATES.X / 2, GATES.Y / 2], [1.0, 1.0])

    qtraj = UnitaryTrajectory(sys, GATES.H, T)
    integrator = HermitianExponentialIntegrator(qtraj, N)
    qcp = SmoothPulseProblem(
        qtraj,
        N;
        integrator = integrator,
        Q = 100.0,
        R = 1e-2,
        piccolo_options = PiccoloOptions(verbose = false),
    )

    E_ops = [Matrix{ComplexF64}(GATES.Z)]

    rcp = RobustControlProblem(
        qcp;
        error_operators = E_ops,
        final_fidelity = 0.99,
        Q_robustness = 10.0,
    )

    @test rcp isa QuantumControlProblem

    # Should have at least one constraint (the fidelity floor)
    @test length(rcp.prob.constraints) >= 1

    # Solve a few iterations with L-BFGS (no Hessian needed)
    solve!(rcp; max_iter = 10, print_level = 0, eval_hessian = false)

    final_traj = get_trajectory(rcp)
    @test final_traj isa NamedTrajectory

    println("✓ RobustControlProblem wrapper test passed")
end

@testitem "RobustControlProblem with SplinePulseProblem" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    T = 1.0
    N = 10
    sys = OpenQuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])

    times = collect(range(0.0, T, length = N))
    pulse = LinearSplinePulse(fill(0.5, 2, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, GATES.H)
    integrator = SplineIntegrator(qtraj, N; spline_order = 1)

    qcp = SplinePulseProblem(qtraj, N; Q = 100.0, R = 1e-2, integrator = integrator)

    E_ops = [Matrix{ComplexF64}(GATES.Z)]

    rcp = RobustControlProblem(
        qcp;
        error_operators = E_ops,
        sys = sys,
        final_fidelity = 0.99,
        Q_robustness = 10.0,
    )

    @test rcp isa QuantumControlProblem

    # Solve a few iterations
    solve!(rcp; max_iter = 10, print_level = 0, eval_hessian = false)

    final_traj = get_trajectory(rcp)
    @test final_traj isa NamedTrajectory

    println("✓ RobustControlProblem with SplinePulseProblem test passed")
end

# ============================================================================ #
#                  KetAdjointRobustnessObjective Tests                         #
# ============================================================================ #

@testitem "KetAdjointRobustnessObjective MultiKetTrajectory basic" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    T = 5.0
    N = 8
    sys = OpenQuantumSystem(GATES.Z / 2, [GATES.X / 2, GATES.Y / 2], [1.0, 1.0])

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    initials = [ψ0, ψ1]
    goals = [ψ1, ψ0]

    qtraj = MultiKetTrajectory(sys, initials, goals, T)
    traj = NamedTrajectory(qtraj, N)
    integrator = HermitianExponentialIntegrator(qtraj, N)

    E_ops = [Matrix{ComplexF64}(GATES.Z)]

    rob_obj = KetAdjointRobustnessObjective(integrator, E_ops, goals, traj; Q = 1.0)

    # Test objective value is real and non-negative
    J = objective_value(rob_obj, traj)
    @test J isa Float64
    @test J >= 0.0

    # Test gradient + Hessian vs finite diff
    test_robustness_objective(rob_obj, traj; atol = 1e-9, rtol = 1e-9)

    println("✓ KetAdjointRobustnessObjective MultiKetTrajectory basic test passed")
end

@testitem "KetAdjointRobustnessObjective MultiKetTrajectory multiple errors" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    T = 5.0
    N = 8
    sys = OpenQuantumSystem(GATES.Z / 2, [GATES.X / 2, GATES.Y / 2], [1.0, 1.0])

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    initials = [ψ0, ψ1]
    goals = [ψ1, ψ0]

    qtraj = MultiKetTrajectory(sys, initials, goals, T)
    traj = NamedTrajectory(qtraj, N)
    integrator = HermitianExponentialIntegrator(qtraj, N)

    E_X = [Matrix{ComplexF64}(GATES.X)]
    E_Z = [Matrix{ComplexF64}(GATES.Z)]
    E_both = [Matrix{ComplexF64}(GATES.X), Matrix{ComplexF64}(GATES.Z)]

    rob_X = KetAdjointRobustnessObjective(integrator, E_X, goals, traj; Q = 1.0)
    rob_Z = KetAdjointRobustnessObjective(integrator, E_Z, goals, traj; Q = 1.0)
    rob_both = KetAdjointRobustnessObjective(integrator, E_both, goals, traj; Q = 1.0)

    @test objective_value(rob_both, traj) ≈
          objective_value(rob_X, traj) + objective_value(rob_Z, traj)

    test_robustness_objective(rob_both, traj; atol = 1e-9, rtol = 1e-9)

    println(
        "✓ KetAdjointRobustnessObjective MultiKetTrajectory multiple errors test passed",
    )
end

@testitem "KetAdjointRobustnessObjective V=0 for zero error" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    T = 5.0
    N = 6
    sys = OpenQuantumSystem(GATES.Z / 2, [GATES.X / 2], [1.0])

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    initials = [ψ0, ψ1]
    goals = [ψ1, ψ0]

    qtraj = MultiKetTrajectory(sys, initials, goals, T)
    traj = NamedTrajectory(qtraj, N)
    integrator = HermitianExponentialIntegrator(qtraj, N)

    E_zero = [zeros(ComplexF64, 2, 2)]
    rob_obj = KetAdjointRobustnessObjective(integrator, E_zero, goals, traj; Q = 1.0)
    @test objective_value(rob_obj, traj) ≈ 0.0 atol = 1e-12

    println("✓ KetAdjointRobustnessObjective V=0 for zero error test passed")
end

@testitem "KetAdjointRobustnessObjective KetTrajectory" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    T = 5.0
    N = 8
    sys = OpenQuantumSystem(GATES.Z / 2, [GATES.X / 2, GATES.Y / 2], [1.0, 1.0])

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]

    qtraj = KetTrajectory(sys, ψ0, ψ1, T)
    traj = NamedTrajectory(qtraj, N)
    integrator = HermitianExponentialIntegrator(qtraj, N)

    E_ops = [Matrix{ComplexF64}(GATES.Z)]
    rob_obj = KetAdjointRobustnessObjective(integrator, E_ops, ψ1, traj; Q = 1.0)

    J = objective_value(rob_obj, traj)
    @test J isa Float64
    @test J >= 0.0

    test_robustness_objective(rob_obj, traj; atol = 1e-9, rtol = 1e-9)

    println("✓ KetAdjointRobustnessObjective KetTrajectory test passed")
end

@testitem "KetAdjointRobustnessObjective Q and error scaling" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    T = 5.0
    N = 8
    sys = OpenQuantumSystem(GATES.Z / 2, [GATES.X / 2], [1.0])

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    initials = [ψ0, ψ1]
    goals = [ψ1, ψ0]

    qtraj = MultiKetTrajectory(sys, initials, goals, T)
    traj = NamedTrajectory(qtraj, N)
    integrator = HermitianExponentialIntegrator(qtraj, N)
    E_ops = [Matrix{ComplexF64}(GATES.Z)]

    # Q scaling
    rob_Q1 = KetAdjointRobustnessObjective(integrator, E_ops, goals, traj; Q = 1.0)
    rob_Q3 = KetAdjointRobustnessObjective(integrator, E_ops, goals, traj; Q = 3.0)

    J1 = objective_value(rob_Q1, traj)
    @test objective_value(rob_Q3, traj) ≈ 3.0 * J1 rtol = 1e-12

    # Error amplitude scaling: J(αE) = α² * J(E)
    α = 2.5
    E_scaled = [α * Matrix{ComplexF64}(GATES.Z)]
    rob_scaled = KetAdjointRobustnessObjective(integrator, E_scaled, goals, traj; Q = 1.0)
    @test objective_value(rob_scaled, traj) ≈ α^2 * J1 rtol = 1e-10

    println("✓ KetAdjointRobustnessObjective Q and error scaling test passed")
end

@testitem "KetAdjointRobustnessObjective 3-level system" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    # 3-level system with two drives
    H_drift = diagm(ComplexF64[1.0, 0.0, -1.0])
    H_drive1 = zeros(ComplexF64, 3, 3)
    H_drive1[1, 2] = H_drive1[2, 1] = 1.0
    H_drive2 = zeros(ComplexF64, 3, 3)
    H_drive2[2, 3] = H_drive2[3, 2] = 1.0

    sys = OpenQuantumSystem(H_drift, [H_drive1, H_drive2], [1.0, 1.0])

    T = 3.0
    N = 8

    e = i -> (v = zeros(ComplexF64, 3); v[i] = 1.0; v)
    ψ0 = e(1)
    ψ1 = e(2)
    ψ2 = e(3)
    initials = [ψ0, ψ1]
    goals = [ψ1, ψ0]

    qtraj = MultiKetTrajectory(sys, initials, goals, T)
    traj = NamedTrajectory(qtraj, N)
    integrator = HermitianExponentialIntegrator(qtraj, N)

    E_dephase = diagm(ComplexF64[1.0, 0.0, -1.0])
    E_ops = [E_dephase]

    rob_obj = KetAdjointRobustnessObjective(integrator, E_ops, goals, traj; Q = 1.0)

    J = objective_value(rob_obj, traj)
    @test J isa Float64
    @test J >= 0.0
    @test rob_obj.ketdim == 3
    @test rob_obj.iso_dim == 6
    @test rob_obj.n_kets == 2

    test_robustness_objective(rob_obj, traj; atol = 1e-9, rtol = 1e-9)

    println("✓ KetAdjointRobustnessObjective 3-level system test passed")
end

@testitem "RobustControlProblem with MultiKetTrajectory" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    T = 10.0
    N = 15
    sys = OpenQuantumSystem(GATES.Z / 2, [GATES.X / 2, GATES.Y / 2], [1.0, 1.0])

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    initials = [ψ0, ψ1]
    goals = [ψ1, ψ0]

    qtraj = MultiKetTrajectory(sys, initials, goals, T)
    integrator = HermitianExponentialIntegrator(qtraj, N)

    qcp = SmoothPulseProblem(
        qtraj,
        N;
        integrator = integrator,
        Q = 100.0,
        R = 1e-2,
        piccolo_options = PiccoloOptions(verbose = false),
    )

    E_ops = [Matrix{ComplexF64}(GATES.Z)]

    rcp = RobustControlProblem(
        qcp;
        error_operators = E_ops,
        final_fidelity = 0.90,
        Q_robustness = 10.0,
    )

    @test rcp isa QuantumControlProblem

    # Should have fidelity constraint(s)
    @test length(rcp.prob.constraints) >= 1

    # Solve a few iterations with L-BFGS
    solve!(rcp; max_iter = 10, print_level = 0, eval_hessian = false)

    final_traj = get_trajectory(rcp)
    @test final_traj isa NamedTrajectory

    println("✓ RobustControlProblem with MultiKetTrajectory test passed")
end

end # module AdjointRobustness
