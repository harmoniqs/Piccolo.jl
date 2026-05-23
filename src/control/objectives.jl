module QuantumObjectives

export KetInfidelityObjective
export CoherentKetInfidelityObjective
export CoherentKetFreePhaseInfidelityObjective
export KetFreePhaseInfidelityObjective
export UnitaryInfidelityObjective
export DensityMatrixInfidelityObjective
export DensityMatrixPureStateInfidelityObjective
export UnitarySensitivityObjective
export UnitaryFreePhaseInfidelityObjective
export LeakageObjective
export BendingEnergyObjective
export build_bending_energy_blocks

using LinearAlgebra
using SparseArrays
using NamedTrajectories
using TrajectoryIndexingUtils
using ...Quantum
using DirectTrajOpt
using TestItems

# --------------------------------------------------------- 
#                       Kets
# ---------------------------------------------------------

function ket_fidelity_loss(ψ̃::AbstractVector, ψ_goal::AbstractVector{<:Complex})
    ψ = iso_to_ket(ψ̃)
    return abs2(ψ_goal' * ψ)
end

"""
    KetInfidelityObjective(ψ̃_name, traj; Q=100.0)

Create a terminal objective for ket state infidelity, using the goal from `traj.goal[ψ̃_name]`.
"""
function KetInfidelityObjective(ψ̃_name::Symbol, traj::NamedTrajectory; Q = 100.0)
    ψ_goal = iso_to_ket(traj.goal[ψ̃_name])
    ℓ = ψ̃ -> abs(1 - ket_fidelity_loss(ψ̃, ψ_goal))
    return TerminalObjective(ℓ, ψ̃_name, traj; Q = Q)
end

"""
    KetInfidelityObjective(ψ_goal, ψ̃_name, traj; Q=100.0)

Create a terminal objective for ket state infidelity with an explicit goal state.

This variant is useful for SamplingProblem and EnsembleTrajectory where the goal
is shared across multiple state variables that don't have individual goals in `traj.goal`.

# Arguments
- `ψ_goal::AbstractVector{<:Complex}`: The target ket state (complex vector)
- `ψ̃_name::Symbol`: Name of the isomorphic state variable in the trajectory
- `traj::NamedTrajectory`: The trajectory

# Keyword Arguments
- `Q::Float64=100.0`: Weight on the infidelity objective
"""
function KetInfidelityObjective(
    ψ_goal::AbstractVector{<:Complex},
    ψ̃_name::Symbol,
    traj::NamedTrajectory;
    Q = 100.0,
)
    ℓ = ψ̃ -> abs(1 - ket_fidelity_loss(ψ̃, ComplexF64.(ψ_goal)))
    return TerminalObjective(ℓ, ψ̃_name, traj; Q = Q)
end

# ---------------------------------------------------------
#                  Coherent Ket Fidelity
# ---------------------------------------------------------

"""
    coherent_ket_fidelity(ψ̃s, ψ_goals)

Compute coherent fidelity across multiple ket states:

    F_coherent = |1/n ∑ᵢ ⟨ψᵢ_goal|ψᵢ⟩|²

This requires all overlaps to have consistent phases (global phase alignment),
which is necessary for implementing gates via state transfer.

# Arguments
- `ψ̃s::Vector{<:AbstractVector}`: List of isomorphic state vectors
- `ψ_goals::Vector{<:AbstractVector{<:Complex}}`: List of goal states
"""
function coherent_ket_fidelity(ψ̃s, ψ_goals::AbstractVector{<:AbstractVector{<:Complex}})
    n = length(ψ̃s)
    @assert n == length(ψ_goals) "Number of states must match number of goals"

    # Sum of overlaps (complex)
    overlap_sum = sum(ψ_goals[i]' * iso_to_ket(ψ̃s[i]) for i = 1:n)

    # Coherent fidelity: |⟨sum⟩/n|²
    return abs2(overlap_sum / n)
end

"""
    CoherentKetInfidelityObjective(ψ_goals, ψ̃_names, traj; Q=100.0)

Create a terminal objective for coherent ket state infidelity across multiple states.

Coherent fidelity is defined as:
    F_coherent = |1/n ∑ᵢ ⟨ψᵢ_goal|ψᵢ⟩|²

Unlike incoherent fidelity (average of individual |⟨ψᵢ_goal|ψᵢ⟩|²), coherent fidelity 
requires all state overlaps to have aligned phases. This is essential when implementing
a gate via multiple state transfers - the gate should have a single global phase,
not independent phases per state.

# Arguments
- `ψ_goals::Vector{<:AbstractVector{<:Complex}}`: Target ket states
- `ψ̃_names::Vector{Symbol}`: Names of isomorphic state variables in trajectory
- `traj::NamedTrajectory`: The trajectory

# Keyword Arguments
- `Q::Float64=100.0`: Weight on the infidelity objective

# Example
```julia
# For implementing X gate via |0⟩→|1⟩ and |1⟩→|0⟩
goals = [ComplexF64[0, 1], ComplexF64[1, 0]]
names = [:ψ̃1, :ψ̃2]
obj = CoherentKetInfidelityObjective(goals, names, traj; Q=100.0)
```
"""
function CoherentKetInfidelityObjective(
    ψ_goals::Vector{<:AbstractVector{<:Complex}},
    ψ̃_names::Vector{Symbol},
    traj::NamedTrajectory;
    Q::Float64 = 100.0,
)
    n_states = length(ψ_goals)
    @assert length(ψ̃_names) == n_states "Number of names must match number of goals"

    # Convert goals to ComplexF64
    goals = [ComplexF64.(g) for g in ψ_goals]

    # Get component indices for each state at terminal time
    state_comps = [traj.components[name] for name in ψ̃_names]
    state_dims = [length(comp) for comp in state_comps]

    # Loss function operating on concatenated terminal states
    function ℓ(z_terminal)
        # Extract each state from the concatenated vector
        ψ̃s = Vector{Vector{eltype(z_terminal)}}(undef, n_states)
        offset = 0
        for i = 1:n_states
            ψ̃s[i] = z_terminal[(offset+1):(offset+state_dims[i])]
            offset += state_dims[i]
        end

        # Coherent infidelity: 1 - F_coherent
        return abs(1 - coherent_ket_fidelity(ψ̃s, goals))
    end

    # Pass vector of component names for multi-component terminal objective
    return TerminalObjective(ℓ, ψ̃_names, traj; Q = Q)
end

# ---------------------------------------------------------
#              Free-Phase Ket Fidelity
# ---------------------------------------------------------

"""
    CoherentKetFreePhaseInfidelityObjective(goals_fn, ψ̃_names, θ_names, traj; Q=100.0)

Coherent ket infidelity with optimizable single-qubit Z-phase rotations.

`goals_fn(θ)` returns phase-rotated goal kets. Phase variables `θ` are stored
in the trajectory's `global_data` and optimized alongside the pulse.

The objective minimizes:
    1 - |1/n Σᵢ ⟨goal_i(θ)|ψ_i⟩|²

where `goal_i(θ) = Φ(θ) · goal_i` with `Φ(θ) = Z₁(θ₁) ⊗ Z₂(θ₂) ⊗ ⋯`.

# Arguments
- `goals_fn::Function`: Maps phase vector `θ` to phased goal kets
- `ψ̃_names::Vector{Symbol}`: Names of isomorphic state variables in trajectory
- `θ_names::AbstractVector{Symbol}`: Names of phase global variables
- `traj::NamedTrajectory`: The trajectory

# Keyword Arguments
- `Q::Float64=100.0`: Weight on the infidelity objective
"""
function CoherentKetFreePhaseInfidelityObjective(
    goals_fn::Function,
    ψ̃_names::Vector{Symbol},
    θ_names::AbstractVector{Symbol},
    traj::NamedTrajectory;
    Q::Float64 = 100.0,
)
    n_states = length(ψ̃_names)
    state_dims = [length(traj.components[name]) for name in ψ̃_names]
    total_state_dim = sum(state_dims)

    function ℓ(z)
        x = z[1:total_state_dim]
        θ = z[(total_state_dim+1):end]

        # Extract individual ket states from concatenated vector
        ψ̃s = Vector{Vector{eltype(x)}}(undef, n_states)
        offset = 0
        for i = 1:n_states
            ψ̃s[i] = x[(offset+1):(offset+state_dims[i])]
            offset += state_dims[i]
        end

        phased_goals = goals_fn(θ)
        return abs(1 - coherent_ket_fidelity(ψ̃s, phased_goals))
    end

    return GlobalKnotPointObjective(
        ℓ,
        ψ̃_names,
        collect(θ_names),
        traj;
        Qs = [Q],
        times = [traj.N],
    )
end

"""
    KetFreePhaseInfidelityObjective(goal_fn, ψ̃_name, θ_names, traj; Q=100.0)

Single-ket infidelity with optimizable Z-phase rotations.

`goal_fn(θ)` returns the phase-rotated goal ket. Minimizes `1 - |⟨goal(θ)|ψ⟩|²`.

# Arguments
- `goal_fn::Function`: Maps phase vector `θ` to phased goal ket
- `ψ̃_name::Symbol`: Name of isomorphic state variable in trajectory
- `θ_names::AbstractVector{Symbol}`: Names of phase global variables
- `traj::NamedTrajectory`: The trajectory

# Keyword Arguments
- `Q::Float64=100.0`: Weight on the infidelity objective
"""
function KetFreePhaseInfidelityObjective(
    goal_fn::Function,
    ψ̃_name::Symbol,
    θ_names::AbstractVector{Symbol},
    traj::NamedTrajectory;
    Q::Float64 = 100.0,
)
    d_state = length(traj.components[ψ̃_name])

    function ℓ(z)
        ψ̃ = z[1:d_state]
        θ = z[(d_state+1):end]
        phased_goal = goal_fn(θ)
        return abs(1 - ket_fidelity_loss(ψ̃, phased_goal))
    end

    return TerminalObjective(ℓ, ψ̃_name, θ_names, traj; Q = Q)
end

# ---------------------------------------------------------
#                       Unitaries
# ---------------------------------------------------------

function unitary_fidelity_loss(
    Ũ⃗::AbstractVector{<:Real},
    U_goal::AbstractMatrix{<:Complex{<:Real}},
)
    U = iso_vec_to_operator(Ũ⃗)
    n = size(U, 1)
    return abs2(tr(U_goal' * U)) / n^2
end

function unitary_fidelity_loss(Ũ⃗::AbstractVector{<:Real}, op::EmbeddedOperator)
    U_goal = unembed(op)
    U = iso_vec_to_operator(Ũ⃗)[op.subspace, op.subspace]
    n = length(op.subspace)
    M = U_goal'U
    return 1 / (n * (n + 1)) * (abs(tr(M'M)) + abs2(tr(M)))
end

function UnitaryInfidelityObjective(
    U_goal::AbstractPiccoloOperator,
    Ũ⃗_name::Symbol,
    traj::NamedTrajectory;
    Q = 100.0,
)
    ℓ = Ũ⃗ -> abs(1 - unitary_fidelity_loss(Ũ⃗, U_goal))
    return TerminalObjective(ℓ, Ũ⃗_name, traj; Q = Q)
end

function UnitaryFreePhaseInfidelityObjective(
    U_goal::Function,
    Ũ⃗_name::Symbol,
    θ_names::AbstractVector{Symbol},
    traj::NamedTrajectory;
    Q = 100.0,
)
    d = sum(traj.global_dims[n] for n in θ_names)
    function ℓ(z)
        Ũ⃗, θ = z[1:(end-d)], z[(end-d+1):end]
        return abs(1 - QuantumObjectives.unitary_fidelity_loss(Ũ⃗, U_goal(θ)))
    end
    return TerminalObjective(ℓ, Ũ⃗_name, θ_names, traj; Q = Q)
end

function UnitaryFreePhaseInfidelityObjective(
    U_goal::Function,
    Ũ⃗_name::Symbol,
    θ_name::Symbol,
    traj::NamedTrajectory;
    kwargs...,
)
    return UnitaryFreePhaseInfidelityObjective(U_goal, Ũ⃗_name, [θ_name], traj; kwargs...)
end

# ---------------------------------------------------------
#                       Density Matrices
# ---------------------------------------------------------

function density_matrix_infidelity_loss(
    ρ̃::AbstractVector,
    ρ_goal::AbstractMatrix{<:Complex{Float64}},
)
    ρ = compact_iso_to_density(ρ̃)
    ℱ = real(tr(ρ * ρ_goal))
    return abs(1 - ℱ)
end

"""
    DensityMatrixInfidelityObjective(ρ̃_name, ρ_goal, traj; Q=100.0)

Terminal objective for density matrix fidelity using the compact isomorphism.

Minimizes `|1 - tr(ρ * ρ_goal)|` where `ρ` is reconstructed from the compact
iso vector via `compact_iso_to_density`.
"""
function DensityMatrixInfidelityObjective(
    ρ̃_name::Symbol,
    ρ_goal::AbstractMatrix{<:Complex{Float64}},
    traj::NamedTrajectory;
    Q = 100.0,
)
    ℓ = ρ̃ -> density_matrix_infidelity_loss(ρ̃, ρ_goal)
    return TerminalObjective(ℓ, ρ̃_name, traj; Q = Q)
end

function density_matrix_pure_state_infidelity_loss(
    ρ̃::AbstractVector,
    ψ::AbstractVector{<:Complex{Float64}},
)
    ρ = compact_iso_to_density(ρ̃)
    ℱ = real(ψ' * ρ * ψ)
    return abs(1 - ℱ)
end

function DensityMatrixPureStateInfidelityObjective(
    ρ̃_name::Symbol,
    ψ_goal::AbstractVector{<:Complex{Float64}},
    traj::NamedTrajectory;
    Q = 100.0,
)
    ℓ = ρ̃ -> density_matrix_pure_state_infidelity_loss(ρ̃, ψ_goal)
    return TerminalObjective(ℓ, ρ̃_name, traj; Q = Q)
end

# ---------------------------------------------------------
#                       Sensitivity
# ---------------------------------------------------------

function unitary_fidelity_loss(Ũ⃗::AbstractVector{<:Real})
    U = iso_vec_to_operator(Ũ⃗)
    n = size(U, 1)
    return abs2(tr(U' * U)) / n^2
end

function UnitarySensitivityObjective(
    name::Symbol,
    traj::NamedTrajectory,
    times::AbstractVector{Int};
    Qs::AbstractVector{<:Float64} = fill(1.0, length(times)),
    scale::Float64 = 1.0,
)
    ℓ = Ũ⃗ -> scale^4 * unitary_fidelity_loss(Ũ⃗)

    return KnotPointObjective(ℓ, name, traj; Qs = Qs, times = times)
end

# ---------------------------------------------------------
#                       Leakage
# ---------------------------------------------------------

"""
    LeakageObjective(indices, name, traj::NamedTrajectory)

Construct a `KnotPointObjective` that penalizes leakage of `name` at the knot points specified by `times` at any `indices` that are outside the computational subspace.

"""
function LeakageObjective(
    indices::AbstractVector{Int},
    name::Symbol,
    traj::NamedTrajectory;
    times = 1:traj.N,
    Qs::AbstractVector{<:Float64} = fill(1.0, length(times)),
)
    leakage_objective(x) = sum(abs2, x[indices]) / length(indices)

    return KnotPointObjective(leakage_objective, name, traj; Qs = Qs, times = times)
end

# ---------------------------------------------------------
#                  Bending-energy regularizer
# ---------------------------------------------------------

@doc raw"""
    build_bending_energy_blocks(h::AbstractVector{<:Real})

Assemble the three sparse symmetric blocks `(K_uu, K_um, K_mm)` of the bending
energy quadratic form for a single-channel cubic Hermite spline through
`N = length(h) + 1` knots, with segment lengths `h[k] = t_{k+1} - t_k`.

The bending energy of a cubic Hermite spline through knot values
``u_1, \ldots, u_N`` with Hermite tangents ``m_1, \ldots, m_N`` admits a closed
form per segment ``[t_k, t_{k+1}]`` of width ``h_k``:

```math
\int_{t_k}^{t_{k+1}} \ddot p(t)^2 \, dt
  = \frac{12(u_{k+1}-u_k)^2}{h_k^3}
  - \frac{12(u_{k+1}-u_k)(m_k+m_{k+1})}{h_k^2}
  + \frac{4(m_k^2 + m_k m_{k+1} + m_{k+1}^2)}{h_k}
```

Summing over segments produces a sparse symmetric quadratic form

```math
E = u^T K_{uu} u + 2\, u^T K_{um} m + m^T K_{mm} m
```

with all three blocks `N × N` and tridiagonal in structure (each segment
couples adjacent knots).

# Reference

Reinsch 1967 derives the analogous Gram matrix for cubic splines parameterised
by knot values alone; this routine extends to the full Hermite ``(u, m)``
representation that `CubicSplinePulse` uses.
"""
function build_bending_energy_blocks(h::AbstractVector{<:Real})
    N = length(h) + 1
    @assert N ≥ 2 "Need at least one segment (length(h) ≥ 1, got $(length(h)))"
    @assert all(>(0), h) "All segment widths h_k must be positive"

    K_uu = spzeros(N, N)
    K_um = spzeros(N, N)
    K_mm = spzeros(N, N)

    @inbounds for k in 1:(N-1)
        hk = h[k]
        h2 = hk^2
        h3 = hk^3

        # u-u block: coefficients of u_k^2, u_{k+1}^2, u_k u_{k+1}
        K_uu[k, k]     += 12 / h3
        K_uu[k+1, k+1] += 12 / h3
        K_uu[k, k+1]   += -12 / h3  # symmetric off-diagonal stores half the cross coeff
        K_uu[k+1, k]   += -12 / h3

        # u-m block: cross term has factor of 2 in the global form
        # The per-segment cross term: -12/h^2 * (u_{k+1}-u_k)(m_k+m_{k+1})
        #   = -12/h^2 (u_{k+1} m_k + u_{k+1} m_{k+1} - u_k m_k - u_k m_{k+1})
        # We absorb the factor of 2 from "E = u^T K_{uu} u + 2 u^T K_{um} m + m^T K_{mm} m"
        # by halving the per-segment coefficients.
        K_um[k, k]     += +6 / h2   # u_k m_k coefficient -(-12/h^2)/2 = +6/h^2
        K_um[k, k+1]   += +6 / h2   # u_k m_{k+1}
        K_um[k+1, k]   += -6 / h2   # u_{k+1} m_k
        K_um[k+1, k+1] += -6 / h2   # u_{k+1} m_{k+1}

        # m-m block: 4/h * (m_k^2 + m_k m_{k+1} + m_{k+1}^2)
        K_mm[k, k]     += 4 / hk
        K_mm[k+1, k+1] += 4 / hk
        K_mm[k, k+1]   += 2 / hk   # cross coeff is 4/h, but symmetric storage halves it
        K_mm[k+1, k]   += 2 / hk
    end

    return K_uu, K_um, K_mm
end

@doc raw"""
    BendingEnergyObjective <: AbstractObjective

Penalises ``\int_0^T \ddot u_c(t)^2 \, dt`` for each control channel ``c`` of a
cubic Hermite spline, summed over channels:

```math
J_{\mathrm{bend}} = \lambda \sum_{c=1}^{n_{\text{drives}}}
    \int_0^T \ddot u_c(t)^2 \, dt
  = \lambda \sum_{c=1}^{n_{\text{drives}}}
    \left( u_c^T K_{uu} u_c + 2\, u_c^T K_{um} m_c + m_c^T K_{mm} m_c \right)
```

where ``u_c`` is the vector of knot values for channel ``c`` and ``m_c`` the
matching Hermite tangents. Unlike `QuadraticRegularizer(:u, ...)` and
`QuadraticRegularizer(:du, ...)`, which penalise pulse *amplitude* and
*tangent magnitude*, this objective penalises the spline's *second derivative*
— the integral of the bending-strain energy. For a flat constant pulse the
contribution is exactly zero regardless of amplitude; for a sinusoid at
frequency ``\omega`` it scales as ``\omega^4``.

# Knot times

The segment widths ``h_k`` used to build the blocks are baked in at
construction time from `get_times(traj)`. Treating ``h_k`` as constant
(rather than re-deriving the blocks from ``\Delta t`` at every evaluation)
keeps the objective a true quadratic — analytical gradient and Hessian — and
matches how Reinsch-style smoothing splines are usually implemented. Free-time
problems still optimise ``\Delta t``, but the bending coefficients reflect
the initial schedule. The penalty remains a soft regulariser, so this is a
fine tradeoff for the typical ~1–3% timestep adjustment seen in min-time
problems.

# Fields
- `u_name::Symbol`: control variable name (e.g. `:u`)
- `du_name::Symbol`: derivative variable name (e.g. `:du`)
- `λ::Float64`: scalar bending-energy weight
- `K_uu::SparseMatrixCSC{Float64,Int}`: ``N×N`` block coupling knot values to themselves
- `K_um::SparseMatrixCSC{Float64,Int}`: ``N×N`` cross block coupling knot values to tangents
- `K_mm::SparseMatrixCSC{Float64,Int}`: ``N×N`` block coupling tangents to themselves
- `n_drives::Int`: number of control channels (the same K blocks apply to each)

# Constructor

```julia
BendingEnergyObjective(u_name, du_name, λ, traj::NamedTrajectory)
```

reads segment widths from `get_times(traj)` and builds the blocks once.
"""
struct BendingEnergyObjective <: AbstractObjective
    u_name::Symbol
    du_name::Symbol
    λ::Float64
    K_uu::SparseMatrixCSC{Float64,Int}
    K_um::SparseMatrixCSC{Float64,Int}
    K_mm::SparseMatrixCSC{Float64,Int}
    n_drives::Int
end

function BendingEnergyObjective(
    u_name::Symbol,
    du_name::Symbol,
    λ::Real,
    traj::NamedTrajectory,
)
    @assert haskey(traj.components, u_name) "BendingEnergyObjective: trajectory has no component :$u_name"
    @assert haskey(traj.components, du_name) "BendingEnergyObjective: trajectory has no component :$du_name"
    @assert traj.dims[u_name] == traj.dims[du_name] "BendingEnergyObjective: $u_name and $du_name must have the same number of channels"

    times = get_times(traj)
    h = diff(times)
    K_uu, K_um, K_mm = build_bending_energy_blocks(h)

    return BendingEnergyObjective(
        u_name,
        du_name,
        Float64(λ),
        K_uu,
        K_um,
        K_mm,
        traj.dims[u_name],
    )
end

function Base.show(io::IO, obj::BendingEnergyObjective)
    print(
        io,
        "BendingEnergyObjective on (:$(obj.u_name), :$(obj.du_name)) ",
        "(λ = $(obj.λ), n_drives = $(obj.n_drives))",
    )
end

# Helper: gather per-channel knot-value and tangent vectors from a trajectory.
# For a control component of size (n_drives, N), traj[name] returns the
# n_drives × N matrix indexed by (channel, knot).
function _bending_extract(traj::NamedTrajectory, name::Symbol)
    return traj[name]  # n_drives × N matrix
end

function DirectTrajOpt.Objectives.objective_value(
    obj::BendingEnergyObjective,
    traj::NamedTrajectory,
)
    U = _bending_extract(traj, obj.u_name)
    M = _bending_extract(traj, obj.du_name)
    J = 0.0
    @inbounds for c in 1:obj.n_drives
        uc = view(U, c, :)
        mc = view(M, c, :)
        J += dot(uc, obj.K_uu * uc) + 2 * dot(uc, obj.K_um * mc) + dot(mc, obj.K_mm * mc)
    end
    return obj.λ * J
end

function DirectTrajOpt.Objectives.gradient!(
    ∇::AbstractVector,
    obj::BendingEnergyObjective,
    traj::NamedTrajectory,
)
    U = _bending_extract(traj, obj.u_name)
    M = _bending_extract(traj, obj.du_name)

    u_comps = traj.components[obj.u_name]
    m_comps = traj.components[obj.du_name]
    N = traj.N

    @inbounds for c in 1:obj.n_drives
        uc = U[c, :]
        mc = M[c, :]

        # ∂J/∂u_c = 2 λ (K_uu u_c + K_um m_c)
        ∂u = 2 * obj.λ * (obj.K_uu * uc + obj.K_um * mc)
        # ∂J/∂m_c = 2 λ (K_um' u_c + K_mm m_c)
        ∂m = 2 * obj.λ * (transpose(obj.K_um) * uc + obj.K_mm * mc)

        # Scatter per-channel gradients into the global ∇ vector.
        # At knot k, the index of variable :u channel c is slice(k, u_comps, traj.dim)[c].
        for k in 1:N
            u_idx = slice(k, u_comps, traj.dim)[c]
            m_idx = slice(k, m_comps, traj.dim)[c]
            ∇[u_idx] += ∂u[k]
            ∇[m_idx] += ∂m[k]
        end
    end

    return nothing
end

function DirectTrajOpt.Objectives.hessian_structure(
    obj::BendingEnergyObjective,
    traj::NamedTrajectory,
)
    Z_dim = traj.dim * traj.N + traj.global_dim
    structure = spzeros(Z_dim, Z_dim)

    N = traj.N

    # The K blocks are tridiagonal in structure, so each adjacent (k, k+1) pair
    # contributes a couple non-zeros. We populate the structure at every nonzero
    # entry of K_uu, K_um, K_mm for every channel.
    @inbounds for c in 1:obj.n_drives
        # Walk K_uu, K_um, K_mm sparse entries
        for (Kmat, name_i, name_j) in (
            (obj.K_uu, obj.u_name, obj.u_name),
            (obj.K_um, obj.u_name, obj.du_name),
            (obj.K_mm, obj.du_name, obj.du_name),
        )
            rows = rowvals(Kmat)
            comps_i = traj.components[name_i]
            comps_j = traj.components[name_j]
            for col in 1:N
                for ptr in nzrange(Kmat, col)
                    row = rows[ptr]
                    i_idx = slice(row, comps_i, traj.dim)[c]
                    j_idx = slice(col, comps_j, traj.dim)[c]
                    structure[i_idx, j_idx] = 1.0
                    structure[j_idx, i_idx] = 1.0
                end
            end
        end
    end

    return structure
end

function DirectTrajOpt.Objectives.get_full_hessian(
    obj::BendingEnergyObjective,
    traj::NamedTrajectory,
)
    Z_dim = traj.dim * traj.N + traj.global_dim
    H = spzeros(Z_dim, Z_dim)

    u_comps = traj.components[obj.u_name]
    m_comps = traj.components[obj.du_name]
    N = traj.N
    λ = obj.λ

    @inbounds for c in 1:obj.n_drives
        # ∂²J/∂u_c ∂u_c = 2λ K_uu
        for col in 1:N
            for ptr in nzrange(obj.K_uu, col)
                row = rowvals(obj.K_uu)[ptr]
                val = nonzeros(obj.K_uu)[ptr]
                i_idx = slice(row, u_comps, traj.dim)[c]
                j_idx = slice(col, u_comps, traj.dim)[c]
                H[i_idx, j_idx] += 2 * λ * val
            end
        end

        # ∂²J/∂u_c ∂m_c = 2λ K_um  (and symmetric transpose)
        for col in 1:N
            for ptr in nzrange(obj.K_um, col)
                row = rowvals(obj.K_um)[ptr]
                val = nonzeros(obj.K_um)[ptr]
                u_idx = slice(row, u_comps, traj.dim)[c]
                m_idx = slice(col, m_comps, traj.dim)[c]
                H[u_idx, m_idx] += 2 * λ * val
                H[m_idx, u_idx] += 2 * λ * val
            end
        end

        # ∂²J/∂m_c ∂m_c = 2λ K_mm
        for col in 1:N
            for ptr in nzrange(obj.K_mm, col)
                row = rowvals(obj.K_mm)[ptr]
                val = nonzeros(obj.K_mm)[ptr]
                i_idx = slice(row, m_comps, traj.dim)[c]
                j_idx = slice(col, m_comps, traj.dim)[c]
                H[i_idx, j_idx] += 2 * λ * val
            end
        end
    end

    return H
end

# ---------------------------------------------------------
#                       Tests
# ---------------------------------------------------------

using TestItems

@testitem "CoherentKetInfidelityObjective" begin
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    # Create a simple trajectory with two ket states
    N = 10
    ket_dim = 4  # iso dim for 2-level system

    # Two state variables
    ψ̃1 = normalize(randn(ket_dim, N))
    ψ̃2 = normalize(randn(ket_dim, N))
    u = randn(1, N)
    Δt = fill(0.1, N)

    traj = NamedTrajectory(
        (ψ̃1 = ψ̃1, ψ̃2 = ψ̃2, u = u, Δt = Δt);
        timestep = :Δt,
        controls = :u,
    )

    # Goal states for X gate: |0⟩→|1⟩ and |1⟩→|0⟩
    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    goals = [ψ1, ψ0]  # |0⟩→|1⟩, |1⟩→|0⟩

    # Create coherent objective
    obj = CoherentKetInfidelityObjective(goals, [:ψ̃1, :ψ̃2], traj; Q = 100.0)

    @test obj isa DirectTrajOpt.Objectives.KnotPointObjective

    # Test that objective can be evaluated
    J = objective_value(obj, traj)
    @test J isa Float64
    @test 0.0 <= J <= 100.0  # Infidelity scaled by Q

    # Test gradient computation
    ∇ = zeros(traj.dim * traj.N + traj.global_dim)
    gradient!(∇, obj, traj)
    @test !all(∇ .== 0)  # Should have non-zero gradient

    # Test coherent vs incoherent behavior:
    # Create perfect states with SAME phase
    ψ̃1_perfect = zeros(ket_dim, N)
    ψ̃2_perfect = zeros(ket_dim, N)
    for k = 1:N
        ψ̃1_perfect[:, k] = ket_to_iso(ψ1)  # |0⟩ should go to |1⟩
        ψ̃2_perfect[:, k] = ket_to_iso(ψ0)  # |1⟩ should go to |0⟩
    end

    traj_perfect = NamedTrajectory(
        (ψ̃1 = ψ̃1_perfect, ψ̃2 = ψ̃2_perfect, u = u, Δt = Δt);
        timestep = :Δt,
        controls = :u,
    )

    J_perfect = objective_value(obj, traj_perfect)
    @test J_perfect < 1e-10  # Should be ~0 for perfect coherent transfer

    # Create perfect states with OPPOSITE phases (phase mismatch)
    ψ̃1_phase = zeros(ket_dim, N)
    ψ̃2_phase = zeros(ket_dim, N)
    for k = 1:N
        ψ̃1_phase[:, k] = ket_to_iso(ψ1)       # +|1⟩
        ψ̃2_phase[:, k] = ket_to_iso(-ψ0)      # -|0⟩ (opposite phase!)
    end

    traj_phase = NamedTrajectory(
        (ψ̃1 = ψ̃1_phase, ψ̃2 = ψ̃2_phase, u = u, Δt = Δt);
        timestep = :Δt,
        controls = :u,
    )

    obj_phase = CoherentKetInfidelityObjective(goals, [:ψ̃1, :ψ̃2], traj_phase; Q = 100.0)
    J_phase = objective_value(obj_phase, traj_phase)

    # Coherent fidelity should be low due to phase mismatch!
    # overlap_sum = ⟨ψ1|ψ1⟩ + ⟨ψ0|(-ψ0)⟩ = 1 + (-1) = 0
    # F_coherent = |0/2|² = 0
    @test J_phase > 50.0  # Should be high infidelity (close to Q * 1.0)
end

@testitem "coherent_ket_fidelity accepts generic Complex types" begin
    using LinearAlgebra
    using ForwardDiff

    # Test that coherent_ket_fidelity works with Complex{Float32} (relaxed from Complex{Float64})
    ψ0_f32 = Complex{Float32}[1.0f0, 0.0f0]
    ψ1_f32 = Complex{Float32}[0.0f0, 1.0f0]
    goals_f32 = [ψ1_f32, ψ0_f32]

    ψ̃s = [ket_to_iso(ComplexF64.(ψ1_f32)), ket_to_iso(ComplexF64.(ψ0_f32))]
    F = Piccolo.QuantumObjectives.coherent_ket_fidelity(ψ̃s, goals_f32)
    @test F ≈ 1.0

    # Test that it works with Complex{ForwardDiff.Dual} (needed for autodiff)
    # This validates the type signature relaxation
    goals_f64 = [ComplexF64[0.0, 1.0], ComplexF64[1.0, 0.0]]
    F2 = Piccolo.QuantumObjectives.coherent_ket_fidelity(ψ̃s, goals_f64)
    @test F2 ≈ 1.0
end

@testitem "CoherentKetFreePhaseInfidelityObjective" begin
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    # Setup: 2-level system, 2 kets for X gate
    N = 10
    ket_dim = 4  # iso dim for 2-level system

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    goals = [ψ1, ψ0]

    # Build a goals_fn that applies Z-phase to each goal
    # For a single qubit: phase_diag = [1, exp(iθ)]
    function goals_fn(θ)
        phase_diag = [one(eltype(θ)), exp(im * θ[1])]
        return [phase_diag .* g for g in goals]
    end

    # Create trajectory with global phase variable
    ψ̃1 = zeros(ket_dim, N)
    ψ̃2 = zeros(ket_dim, N)
    for k = 1:N
        ψ̃1[:, k] = ket_to_iso(ψ1)
        ψ̃2[:, k] = ket_to_iso(ψ0)
    end

    traj = NamedTrajectory(
        (ψ̃1 = ψ̃1, ψ̃2 = ψ̃2, u = randn(1, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
        global_data = [0.0],
        global_components = (φ_1 = 1:1,),
    )

    θ_names = [:φ_1]
    obj = CoherentKetFreePhaseInfidelityObjective(
        goals_fn,
        [:ψ̃1, :ψ̃2],
        θ_names,
        traj;
        Q = 100.0,
    )

    # With θ=0, perfect states should give ~0 infidelity
    J = objective_value(obj, traj)
    @test J < 1e-10

    # Test gradient is computable
    ∇ = zeros(traj.dim * traj.N + traj.global_dim)
    gradient!(∇, obj, traj)

    # Test with arbitrary (normalized) states: objective should be in [0, Q].
    # Per-knot iso kets normalized to unit norm so the Pedersen-style fidelity
    # is bounded in [0, 1] and the assertion is well-defined regardless of seed.
    ψ̃1_rand = randn(ket_dim, N)
    ψ̃2_rand = randn(ket_dim, N)
    for k = 1:N
        ψ̃1_rand[:, k] ./= norm(ψ̃1_rand[:, k])
        ψ̃2_rand[:, k] ./= norm(ψ̃2_rand[:, k])
    end
    u_init = 0.05 * cos.(reshape(2π .* (0:(N-1)) ./ max(N - 1, 1), 1, N))

    traj_rand = NamedTrajectory(
        (ψ̃1 = ψ̃1_rand, ψ̃2 = ψ̃2_rand, u = u_init, Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
        global_data = [0.5],
        global_components = (φ_1 = 1:1,),
    )

    obj_rand = CoherentKetFreePhaseInfidelityObjective(
        goals_fn,
        [:ψ̃1, :ψ̃2],
        θ_names,
        traj_rand;
        Q = 100.0,
    )
    J_rand = objective_value(obj_rand, traj_rand)
    @test isfinite(J_rand)
    @test 0.0 <= J_rand <= 100.0
end

@testitem "KetFreePhaseInfidelityObjective" begin
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    N = 10
    ket_dim = 4  # iso dim for 2-level system
    ψ_goal = ComplexF64[0.0, 1.0]

    function goal_fn(θ)
        phase_diag = [one(eltype(θ)), exp(im * θ[1])]
        return phase_diag .* ψ_goal
    end

    # Perfect state at final time
    ψ̃ = zeros(ket_dim, N)
    for k = 1:N
        ψ̃[:, k] = ket_to_iso(ψ_goal)
    end

    traj = NamedTrajectory(
        (ψ̃ = ψ̃, u = randn(1, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
        global_data = [0.0],
        global_components = (φ_1 = 1:1,),
    )

    obj = KetFreePhaseInfidelityObjective(goal_fn, :ψ̃, [:φ_1], traj; Q = 100.0)

    # Perfect state with zero phase -> zero infidelity
    J = objective_value(obj, traj)
    @test J < 1e-10

    # Test gradient
    ∇ = zeros(traj.dim * traj.N + traj.global_dim)
    gradient!(∇, obj, traj)
end

end
