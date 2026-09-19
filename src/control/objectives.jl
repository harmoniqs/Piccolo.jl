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
export HermiteBendingEnergyRegularizer

using LinearAlgebra
using SparseArrays
using NamedTrajectories
using ...Quantum
using DirectTrajOpt
using TrajectoryIndexingUtils
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
    coherent_ket_fidelity(ψ̃s, ψ_goals; weights=nothing)

Compute coherent fidelity across multiple ket states:

    F_coherent = |1/n ∑ᵢ ⟨ψᵢ_goal|ψᵢ⟩|²

This requires all overlaps to have consistent phases (global phase alignment),
which is necessary for implementing gates via state transfer.

With per-state `weights`, the mean of overlaps becomes a *weighted* mean,
normalized by the weight sum:

    F_coherent = |∑ᵢ wᵢ ⟨ψᵢ_goal|ψᵢ⟩ / ∑ᵢ wᵢ|²

Only the ratios between weights matter — they need not be normalized. Uniform
weights are the identity case, recovering the unweighted formula exactly.

# Arguments
- `ψ̃s::Vector{<:AbstractVector}`: List of isomorphic state vectors
- `ψ_goals::Vector{<:AbstractVector{<:Complex}}`: List of goal states

# Keyword Arguments
- `weights::Union{Nothing, AbstractVector{<:Real}}=nothing`: Per-state weights.
  `nothing` (the default) means uniform weighting.
"""
function coherent_ket_fidelity(
    ψ̃s,
    ψ_goals::AbstractVector{<:AbstractVector{<:Complex}};
    weights::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    n = length(ψ̃s)
    @assert n == length(ψ_goals) "Number of states must match number of goals"
    isnothing(weights) ||
        @assert n == length(weights) "Number of states must match number of weights"

    # Uniform weights are the identity case — take the unweighted path so the
    # value is bit-for-bit what it was before weights existed
    if isnothing(weights) || allequal(weights)
        # Sum of overlaps (complex)
        overlap_sum = sum(ψ_goals[i]' * iso_to_ket(ψ̃s[i]) for i = 1:n)

        # Coherent fidelity: |⟨sum⟩/n|²
        return abs2(overlap_sum / n)
    else
        # Weighted sum of overlaps (complex)
        overlap_sum = sum(weights[i] * (ψ_goals[i]' * iso_to_ket(ψ̃s[i])) for i = 1:n)

        # Weighted coherent fidelity: |⟨weighted sum⟩/∑w|²
        return abs2(overlap_sum / sum(weights))
    end
end

"""
    coherent_fidelity_weights(weights, n_states) -> Union{Nothing, Vector{Float64}}

Canonicalize per-state weights for coherent fidelity.

Returns `nothing` — the unweighted code path — when `weights` is `nothing` or
uniform. Uniform weights are mathematically the identity case, and routing them
back to the unweighted formula makes that identity hold *bit-for-bit* rather
than up to floating-point rounding, so existing unweighted results do not move.

Otherwise returns the weights normalized to sum to one, so that a constructed
objective captures self-describing weights.
"""
coherent_fidelity_weights(::Nothing, ::Int) = nothing

function coherent_fidelity_weights(weights::AbstractVector{<:Real}, n_states::Int)
    @assert length(weights) == n_states "Number of weights must match number of states"
    @assert all(≥(0), weights) "Weights must be non-negative"
    @assert sum(weights) > 0 "Weights must not all be zero"
    allequal(weights) && return nothing
    return collect(Float64, weights) ./ sum(weights)
end

"""
    CoherentKetInfidelityObjective(ψ_goals, ψ̃_names, traj; Q=100.0)

Create a terminal objective for coherent ket state infidelity across multiple states.

Coherent fidelity is defined as:
    F_coherent = |∑ᵢ wᵢ ⟨ψᵢ_goal|ψᵢ⟩ / ∑ᵢ wᵢ|²

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
- `weights::Union{Nothing, AbstractVector{<:Real}}=nothing`: Per-state weights on the
  coherent mean of overlaps. Normalized at construction; only their ratios matter.
  `nothing` or uniform weights give the unweighted fidelity |1/n ∑ᵢ ⟨ψᵢ_goal|ψᵢ⟩|².

# Example
```julia
# For implementing X gate via |0⟩→|1⟩ and |1⟩→|0⟩
goals = [ComplexF64[0, 1], ComplexF64[1, 0]]
names = [:ψ̃1, :ψ̃2]
obj = CoherentKetInfidelityObjective(goals, names, traj; Q=100.0)

# Emphasize the first state transfer
obj = CoherentKetInfidelityObjective(goals, names, traj; Q=100.0, weights=[0.9, 0.1])
```
"""
function CoherentKetInfidelityObjective(
    ψ_goals::Vector{<:AbstractVector{<:Complex}},
    ψ̃_names::Vector{Symbol},
    traj::NamedTrajectory;
    Q::Float64 = 100.0,
    weights::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    n_states = length(ψ_goals)
    @assert length(ψ̃_names) == n_states "Number of names must match number of goals"

    # Convert goals to ComplexF64
    goals = [ComplexF64.(g) for g in ψ_goals]

    # Normalize once at construction, so the loss captures self-describing weights
    ws = coherent_fidelity_weights(weights, n_states)

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
        return abs(1 - coherent_ket_fidelity(ψ̃s, goals; weights = ws))
    end

    # Pass vector of component names for multi-component terminal objective.
    # (Matrix-free per-knot HVP carrier construction lives in Piccolissimo.)
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
    1 - |Σᵢ wᵢ ⟨goal_i(θ)|ψ_i⟩ / Σᵢ wᵢ|²

where `goal_i(θ) = Φ(θ) · goal_i` with `Φ(θ) = Z₁(θ₁) ⊗ Z₂(θ₂) ⊗ ⋯`. Weights
apply to the *phased* overlaps, so weighting composes with the phase rotation.

# Arguments
- `goals_fn::Function`: Maps phase vector `θ` to phased goal kets
- `ψ̃_names::Vector{Symbol}`: Names of isomorphic state variables in trajectory
- `θ_names::AbstractVector{Symbol}`: Names of phase global variables
- `traj::NamedTrajectory`: The trajectory

# Keyword Arguments
- `Q::Float64=100.0`: Weight on the infidelity objective
- `weights::Union{Nothing, AbstractVector{<:Real}}=nothing`: Per-state weights on the
  coherent mean of overlaps. Normalized at construction; only their ratios matter.
  `nothing` or uniform weights give the unweighted fidelity |1/n Σᵢ ⟨goal_i(θ)|ψ_i⟩|².
"""
function CoherentKetFreePhaseInfidelityObjective(
    goals_fn::Function,
    ψ̃_names::Vector{Symbol},
    θ_names::AbstractVector{Symbol},
    traj::NamedTrajectory;
    Q::Float64 = 100.0,
    weights::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    n_states = length(ψ̃_names)
    state_dims = [length(traj.components[name]) for name in ψ̃_names]
    total_state_dim = sum(state_dims)

    # Normalize once at construction, so the loss captures self-describing weights
    ws = coherent_fidelity_weights(weights, n_states)

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
        return abs(1 - coherent_ket_fidelity(ψ̃s, phased_goals; weights = ws))
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
    # Matrix-free per-knot HVP carrier construction lives in Piccolissimo.
    ℓ = Ũ⃗ -> abs(1 - unitary_fidelity_loss(Ũ⃗, U_goal))
    return TerminalObjective(ℓ, Ũ⃗_name, traj; Q = Q)
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
    times = 1:(traj.N),
    Qs::AbstractVector{<:Float64} = fill(1.0, length(times)),
)
    leakage_objective(x) = sum(abs2, x[indices]) / length(indices)

    return KnotPointObjective(leakage_objective, name, traj; Qs = Qs, times = times)
end

# ---------------------------------------------------------
#           Bending Energy (cubic-Hermite smoothness)
# ---------------------------------------------------------

# Ported from Piccolissimo (hermite_bending_energy_regularizer.jl) per #309:
# smoothness-as-chassis. Bending energy ∫u″²dt is the quantity that predicts
# how a model-optimal pulse survives contact with reality (internal P-series
# calibration studies, Aug 2026), so the closed-form regularizer lives in the
# public package and SplinePulseProblem defaults it ON for cubic splines.

"""
    HermiteBendingEnergyRegularizer <: AbstractObjective

Penalizes the integrated squared curvature (bending energy) of a cubic Hermite
spline pulse:

```
    J = (1/2) Σ_d R_d  ∫₀^T u_d''(t)² dt
```

For a Hermite cubic on segment k of length Δt_k with endpoint accelerations
`a_start(k), a_end(k)`, the per-segment integral has the closed form

```
    ∫₀^Δt_k f''(τ)² dτ = (Δt_k / 3) (a_start² + a_end² + a_start·a_end)
```

(exact: f''(τ) is linear on each segment, so Simpson at the two endpoints +
midpoint is exact).

Only meaningful for cubic-Hermite-parameterized trajectories (components
`u`, `du`, `Δt`): for a piecewise-constant or linear-spline pulse the second
derivative is distributional (Dirac deltas) and bending energy is not a
property of the pulse.

# Fields
- `R::Vector{Float64}`: Per-drive regularization weight (length `n_drives`).
- `control_name::Symbol`: Name of control variable (default `:u`).
- `derivative_name::Symbol`: Name of derivative variable (default `:du`).
- `timestep_name::Symbol`: Name of timestep variable (default `:Δt`).

# Constructor
```julia
HermiteBendingEnergyRegularizer(
    traj::NamedTrajectory;
    R::Union{Real, AbstractVector{<:Real}} = 1.0,
    control_name::Symbol = :u,
    derivative_name::Symbol = :du,
    timestep_name::Symbol = :Δt,
)
```

If `R` is a scalar, it is broadcast across all drives. If it is a vector, its
length must equal `n_drives` — useful when channels have very different units.
"""
struct HermiteBendingEnergyRegularizer <: AbstractObjective
    R::Vector{Float64}
    control_name::Symbol
    derivative_name::Symbol
    timestep_name::Symbol
end

function HermiteBendingEnergyRegularizer(
    traj::NamedTrajectory;
    R::Union{Real,AbstractVector{<:Real}} = 1.0,
    control_name::Symbol = :u,
    derivative_name::Symbol = :du,
    timestep_name::Symbol = :Δt,
)
    @assert haskey(traj.components, control_name) "Trajectory must have component :$control_name"
    @assert haskey(traj.components, derivative_name) "Trajectory must have component :$derivative_name"
    @assert haskey(traj.components, timestep_name) "Trajectory must have component :$timestep_name"

    n_drives = traj.dims[control_name]
    R_vec = R isa Real ? fill(Float64(R), n_drives) : Vector{Float64}(R)
    @assert length(R_vec) == n_drives "R must be scalar or length-$n_drives vector (got length $(length(R_vec)))"

    return HermiteBendingEnergyRegularizer(
        R_vec,
        control_name,
        derivative_name,
        timestep_name,
    )
end

# Endpoint accelerations of a Hermite cubic on segment k.
@inline function _bending_a_start(u_k, u_kp1, du_k, du_kp1, Δt_k)
    inv_Δt = inv(Δt_k)
    inv_Δt2 = inv_Δt * inv_Δt
    return 6.0 * inv_Δt2 * (u_kp1 - u_k) - 2.0 * inv_Δt * (2.0 * du_k + du_kp1)
end

@inline function _bending_a_end(u_k, u_kp1, du_k, du_kp1, Δt_k)
    inv_Δt = inv(Δt_k)
    inv_Δt2 = inv_Δt * inv_Δt
    return -6.0 * inv_Δt2 * (u_kp1 - u_k) + 2.0 * inv_Δt * (du_k + 2.0 * du_kp1)
end

# Per-segment per-drive objective kernel for ForwardDiff Hessian assembly.
# Input z = (u_k, u_{k+1}, du_k, du_{k+1}, Δt_k); returns (R_d/2) ∫ f''² dτ.
@inline function _bending_seg_kernel(z::AbstractVector{T}, R_d::Real) where {T}
    u_k, u_kp1, du_k, du_kp1, Δt_k = z[1], z[2], z[3], z[4], z[5]
    a_s = _bending_a_start(u_k, u_kp1, du_k, du_kp1, Δt_k)
    a_e = _bending_a_end(u_k, u_kp1, du_k, du_kp1, Δt_k)
    return (T(R_d) * Δt_k / 6.0) * (a_s^2 + a_e^2 + a_s * a_e)
end

function DirectTrajOpt.objective_value(
    reg::HermiteBendingEnergyRegularizer,
    traj::NamedTrajectory,
)
    u = traj[reg.control_name]
    du = traj[reg.derivative_name]
    Δt = traj[reg.timestep_name]

    n_drives, N = size(u)
    J = 0.0

    @inbounds for k = 1:(N-1)
        Δt_k = Δt[k]
        Δt_over_6 = Δt_k / 6.0   # (R_i / 2) · (Δt / 3) = R_i · Δt / 6
        for i = 1:n_drives
            a_s = _bending_a_start(u[i, k], u[i, k+1], du[i, k], du[i, k+1], Δt_k)
            a_e = _bending_a_end(u[i, k], u[i, k+1], du[i, k], du[i, k+1], Δt_k)
            J += reg.R[i] * Δt_over_6 * (a_s^2 + a_e^2 + a_s * a_e)
        end
    end

    return J
end

function DirectTrajOpt.gradient!(
    ∇J::AbstractVector,
    reg::HermiteBendingEnergyRegularizer,
    traj::NamedTrajectory,
)
    u = traj[reg.control_name]
    du = traj[reg.derivative_name]
    Δt = traj[reg.timestep_name]

    u_comps = traj.components[reg.control_name]
    du_comps = traj.components[reg.derivative_name]
    Δt_comp1 = traj.components[reg.timestep_name][1]  # scalar offset

    n_drives, N = size(u)
    D = traj.dim

    @inbounds for k = 1:(N-1)
        Δt_k = Δt[k]
        inv_Δt = 1.0 / Δt_k
        inv_Δt2 = inv_Δt * inv_Δt
        inv_Δt3 = inv_Δt2 * inv_Δt

        # Hoist segment base offsets out of the drive loop.
        base_k = D * (k - 1)
        base_kp1 = D * k
        idx_Δt = base_k + Δt_comp1

        for i = 1:n_drives
            u_k = u[i, k]
            u_kp1 = u[i, k+1]
            du_k = du[i, k]
            du_kp1 = du[i, k+1]

            a_s = _bending_a_start(u_k, u_kp1, du_k, du_kp1, Δt_k)
            a_e = _bending_a_end(u_k, u_kp1, du_k, du_kp1, Δt_k)

            # J_seg = (R_i/2)(Δt/3)(a_s² + a_e² + a_s·a_e) = (R_i Δt / 6) S.
            # ∂S/∂a_s = 2 a_s + a_e ;  ∂S/∂a_e = a_s + 2 a_e.
            coef = reg.R[i] * Δt_k / 6.0
            ws = coef * (2.0 * a_s + a_e)
            we = coef * (a_s + 2.0 * a_e)

            # Partial derivatives of a_s, a_e:
            #   ∂a_s/∂u_k = -6/Δt² ;  ∂a_s/∂u_{k+1} = +6/Δt²
            #   ∂a_s/∂du_k = -4/Δt ;  ∂a_s/∂du_{k+1} = -2/Δt
            #   ∂a_e/∂u_k = +6/Δt² ;  ∂a_e/∂u_{k+1} = -6/Δt²
            #   ∂a_e/∂du_k = +2/Δt ;  ∂a_e/∂du_{k+1} = +4/Δt
            u_i = u_comps[i]
            du_i = du_comps[i]
            idx_u_k = base_k + u_i
            idx_u_kp1 = base_kp1 + u_i
            idx_du_k = base_k + du_i
            idx_du_kp1 = base_kp1 + du_i

            six_inv_Δt2 = 6.0 * inv_Δt2
            two_inv_Δt = 2.0 * inv_Δt
            four_inv_Δt = 4.0 * inv_Δt

            ∇J[idx_u_k] += -ws * six_inv_Δt2 + we * six_inv_Δt2
            ∇J[idx_u_kp1] += ws * six_inv_Δt2 - we * six_inv_Δt2
            ∇J[idx_du_k] += -ws * four_inv_Δt + we * two_inv_Δt
            ∇J[idx_du_kp1] += -ws * two_inv_Δt + we * four_inv_Δt

            # ∂a_s/∂Δt = -12 (u_{k+1}-u_k)/Δt³ + 2 (2 du_k + du_{k+1})/Δt²
            # ∂a_e/∂Δt = +12 (u_{k+1}-u_k)/Δt³ - 2 (du_k + 2 du_{k+1})/Δt²
            twelve_du_inv_Δt3 = 12.0 * (u_kp1 - u_k) * inv_Δt3
            da_s_dΔt = -twelve_du_inv_Δt3 + 2.0 * (2.0 * du_k + du_kp1) * inv_Δt2
            da_e_dΔt = twelve_du_inv_Δt3 - 2.0 * (du_k + 2.0 * du_kp1) * inv_Δt2

            S = a_s^2 + a_e^2 + a_s * a_e
            ∇J[idx_Δt] += reg.R[i] * S / 6.0 + ws * da_s_dΔt + we * da_e_dΔt
        end
    end

    return nothing
end

# ---------------------------------------------------------------------------- #
#  Analytical per-segment Hessian block                                         #
# ---------------------------------------------------------------------------- #
# For J_seg = coef·S with coef = R·Δt/6 and S = a_s² + a_e² + a_s·a_e, the exact
# 5×5 Hessian over z = (u_k, u_{k+1}, du_k, du_{k+1}, Δt) (indices 1..5) is the
# full product-rule expansion. a_s, a_e are LINEAR in u/du (so ∂²a vanishes there)
# but NONLINEAR in Δt; coef is linear in Δt. Writing ∂a_s, ∂a_e and ∂²a_s, ∂²a_e
# (the latter only along Δt), and with cS = 2a_s+a_e, cE = a_s+2a_e:
#
#   ∂S/∂x   = cS·∂a_s/∂x + cE·∂a_e/∂x
#   ∂²S/∂xy = cS·∂²a_s/∂xy + cE·∂²a_e/∂xy
#             + 2 a_s' a_s'' + 2 a_e' a_e'' + a_s' a_e'' + a_e' a_s''   (sym in x,y)
#   H_xy(x,y≠Δt) = coef·∂²S/∂xy
#   H_x,Δt       = (R/6)·∂S/∂x + coef·∂²S/∂x∂Δt
#   H_Δt,Δt      = 2·(R/6)·∂S/∂Δt + coef·∂²S/∂Δt²
#
# This block matches the per-segment ForwardDiff kernel `_bending_seg_kernel` to
# machine precision. Pure arithmetic on gathered scalars — array-type agnostic.

const _BEND_BLK = 5

"""
    _bending_seg_block!(H5, u_k, u_kp1, du_k, du_kp1, Δt, R_d)

Fill the 5×5 analytical Hessian block `H5` for one segment / drive over
`z = (u_k, u_{k+1}, du_k, du_{k+1}, Δt)`. Pure scalar arithmetic — array-agnostic.
"""
@inline function _bending_seg_block!(H5, u_k, u_kp1, du_k, du_kp1, Δt, R_d)
    inv_Δt = 1.0 / Δt
    inv_Δt2 = inv_Δt * inv_Δt
    inv_Δt3 = inv_Δt2 * inv_Δt
    inv_Δt4 = inv_Δt3 * inv_Δt

    Δu = u_kp1 - u_k
    a_s = 6.0 * inv_Δt2 * Δu - 2.0 * inv_Δt * (2.0 * du_k + du_kp1)
    a_e = -6.0 * inv_Δt2 * Δu + 2.0 * inv_Δt * (du_k + 2.0 * du_kp1)

    coef = R_d * Δt / 6.0
    coefp = R_d / 6.0          # ∂coef/∂Δt

    cS = 2.0 * a_s + a_e       # ∂S/∂a_s
    cE = a_s + 2.0 * a_e       # ∂S/∂a_e

    # First derivatives of a_s, a_e w.r.t. the 5 vars (index order as above).
    das = (
        -6.0 * inv_Δt2,                                            # ∂a_s/∂u_k
        6.0 * inv_Δt2,                                             # ∂a_s/∂u_{k+1}
        -4.0 * inv_Δt,                                            # ∂a_s/∂du_k
        -2.0 * inv_Δt,                                            # ∂a_s/∂du_{k+1}
        -12.0 * inv_Δt3 * Δu + 2.0 * inv_Δt2 * (2.0 * du_k + du_kp1),  # ∂a_s/∂Δt
    )
    dae = (
        6.0 * inv_Δt2,                                             # ∂a_e/∂u_k
        -6.0 * inv_Δt2,                                            # ∂a_e/∂u_{k+1}
        2.0 * inv_Δt,                                             # ∂a_e/∂du_k
        4.0 * inv_Δt,                                             # ∂a_e/∂du_{k+1}
        12.0 * inv_Δt3 * Δu - 2.0 * inv_Δt2 * (du_k + 2.0 * du_kp1),   # ∂a_e/∂Δt
    )

    # Second derivatives of a_s, a_e — nonzero only when one index is Δt (=5),
    # since a is linear in u/du. Stored as the Δt-row vectors (∂²a/∂Δt∂x).
    d2as_dt = (
        12.0 * inv_Δt3,                                           # ∂²a_s/∂Δt∂u_k
        -12.0 * inv_Δt3,                                          # ∂²a_s/∂Δt∂u_{k+1}
        4.0 * inv_Δt2,                                            # ∂²a_s/∂Δt∂du_k
        2.0 * inv_Δt2,                                            # ∂²a_s/∂Δt∂du_{k+1}
        36.0 * inv_Δt4 * Δu - 4.0 * inv_Δt3 * (2.0 * du_k + du_kp1),   # ∂²a_s/∂Δt²
    )
    d2ae_dt = (
        -12.0 * inv_Δt3,                                          # ∂²a_e/∂Δt∂u_k
        12.0 * inv_Δt3,                                           # ∂²a_e/∂Δt∂u_{k+1}
        -2.0 * inv_Δt2,                                           # ∂²a_e/∂Δt∂du_k
        -4.0 * inv_Δt2,                                           # ∂²a_e/∂Δt∂du_{k+1}
        -36.0 * inv_Δt4 * Δu + 4.0 * inv_Δt3 * (du_k + 2.0 * du_kp1),  # ∂²a_e/∂Δt²
    )

    # ∂S/∂x for all x (used by the Δt cross/diagonal blocks).
    @inline dS(x) = cS * das[x] + cE * dae[x]

    @inbounds for a = 1:_BEND_BLK
        for b = a:_BEND_BLK
            # Quadratic-form part of ∂²S/∂a∂b (always present).
            quad =
                2.0 * das[a] * das[b] +
                2.0 * dae[a] * dae[b] +
                das[a] * dae[b] +
                dae[a] * das[b]
            # Second-derivative part: nonzero only if a==5 or b==5 (Δt involved).
            d2 = 0.0
            if a == _BEND_BLK
                d2 = cS * d2as_dt[b] + cE * d2ae_dt[b]
            elseif b == _BEND_BLK
                d2 = cS * d2as_dt[a] + cE * d2ae_dt[a]
            end
            d2S = quad + d2

            if a < _BEND_BLK && b < _BEND_BLK
                Hab = coef * d2S
            elseif a < _BEND_BLK && b == _BEND_BLK
                Hab = coefp * dS(a) + coef * d2S
            elseif a == _BEND_BLK && b < _BEND_BLK
                Hab = coefp * dS(b) + coef * d2S
            else  # a == b == 5 (Δt,Δt)
                Hab = 2.0 * coefp * dS(_BEND_BLK) + coef * d2S
            end

            H5[a, b] = Hab
            H5[b, a] = Hab
        end
    end

    return nothing
end

function DirectTrajOpt.Objectives.hessian_structure(
    reg::HermiteBendingEnergyRegularizer,
    traj::NamedTrajectory,
)
    Z_dim = traj.dim * traj.N + traj.global_dim
    structure = spzeros(Z_dim, Z_dim)

    u_comps = traj.components[reg.control_name]
    du_comps = traj.components[reg.derivative_name]
    Δt_comps = traj.components[reg.timestep_name]

    n_drives = length(u_comps)

    for k = 1:(traj.N-1)
        indices = Int[]
        for d = 1:n_drives
            push!(indices, slice(k, u_comps, traj.dim)[d])
            push!(indices, slice(k+1, u_comps, traj.dim)[d])
            push!(indices, slice(k, du_comps, traj.dim)[d])
            push!(indices, slice(k+1, du_comps, traj.dim)[d])
        end
        push!(indices, slice(k, Δt_comps, traj.dim)[1])

        for i in indices
            for j in indices
                structure[i, j] = 1.0
            end
        end
    end

    return structure
end

"""
    get_full_hessian(reg::HermiteBendingEnergyRegularizer, traj)

EXACT analytical per-segment Hessian. Each (segment k, drive d) contributes
a dense 5×5 block over `(u_k, u_{k+1}, du_k, du_{k+1}, Δt_k)` built by the closed-form
`_bending_seg_block!` (replacing the per-segment ForwardDiff). Blocks are
scattered into the sparse global Hessian; cost is `O(N · n_drives)` tiny blocks and
the result is exact to machine precision (matches the ForwardDiff kernel to ~1e-12).
"""
function DirectTrajOpt.Objectives.get_full_hessian(
    reg::HermiteBendingEnergyRegularizer,
    traj::NamedTrajectory,
)
    Z_dim = traj.dim * traj.N + traj.global_dim
    H = spzeros(Z_dim, Z_dim)

    u = traj[reg.control_name]
    du = traj[reg.derivative_name]
    Δt = traj[reg.timestep_name]

    u_comps = traj.components[reg.control_name]
    du_comps = traj.components[reg.derivative_name]
    Δt_comp1 = traj.components[reg.timestep_name][1]

    n_drives, N = size(u)
    D = traj.dim

    H_seg = Matrix{Float64}(undef, _BEND_BLK, _BEND_BLK)
    idx_buf = Vector{Int}(undef, _BEND_BLK)

    @inbounds for k = 1:(N-1)
        Δt_k = Δt[k]
        base_k = D * (k - 1)
        base_kp1 = D * k
        idx_Δt = base_k + Δt_comp1

        for i = 1:n_drives
            _bending_seg_block!(
                H_seg,
                u[i, k],
                u[i, k+1],
                du[i, k],
                du[i, k+1],
                Δt_k,
                reg.R[i],
            )

            u_i = u_comps[i]
            du_i = du_comps[i]
            idx_buf[1] = base_k + u_i
            idx_buf[2] = base_kp1 + u_i
            idx_buf[3] = base_k + du_i
            idx_buf[4] = base_kp1 + du_i
            idx_buf[5] = idx_Δt

            for a = 1:_BEND_BLK
                ia = idx_buf[a]
                for b = 1:_BEND_BLK
                    H[ia, idx_buf[b]] += H_seg[a, b]
                end
            end
        end
    end

    return H
end

# ---------------------------------------------------------
#                       Tests
# ---------------------------------------------------------

using TestItems

@testitem "HermiteBendingEnergyRegularizer gradient/Hessian" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using Random

    Random.seed!(309)
    T = 5
    u = rand(2, T)
    du = rand(2, T)
    Δt = fill(0.2, T)

    traj = NamedTrajectory((u = u, du = du, Δt = Δt); timestep = :Δt, controls = :u)

    reg_scalar = HermiteBendingEnergyRegularizer(traj; R = 1.0)
    ∇ = zeros(traj.dim * traj.N + traj.global_dim)
    DirectTrajOpt.gradient!(∇, reg_scalar, traj)
    H = DirectTrajOpt.Objectives.get_full_hessian(reg_scalar, traj)
    @test isfinite(DirectTrajOpt.objective_value(reg_scalar, traj))
    @test !all(∇ .== 0)
    @test isfinite(sum(abs2, H))

    reg_vec = HermiteBendingEnergyRegularizer(traj; R = [0.3, 2.7])
    @test reg_vec.R == [0.3, 2.7]
end

@testitem "HermiteBendingEnergyRegularizer matches segment formula" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt

    # Direct check of the closed-form integral on a known segment.
    # Hermite endpoints y1=0, y2=1, m1=0, m2=0 over Δt=1:
    # the cubic is f(τ) = 3τ² - 2τ³ on [0,1], so f''(τ) = 6 - 12τ.
    # ∫₀¹ (6 - 12τ)² dτ = 36 - 72 + 48 = 12. J = (R/2)·12 = 6.0.
    N = 2
    u = reshape([0.0, 1.0], 1, N)
    du = reshape([0.0, 0.0], 1, N)
    Δt = fill(1.0, N)

    traj = NamedTrajectory((; u = u, du = du, Δt = Δt); timestep = :Δt, controls = (:u,))
    reg = HermiteBendingEnergyRegularizer(traj; R = 1.0)

    @test isapprox(DirectTrajOpt.objective_value(reg, traj), 6.0; atol = 1e-12)
end

@testitem "HermiteBendingEnergyRegularizer zero for linear pulse" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt

    # u(t) linear ⇒ f''(t) = 0 ⇒ bending energy = 0
    N = 6
    Δt_val = 0.4
    times = collect(range(0.0, (N - 1) * Δt_val, length = N))

    u = reshape(times, 1, N)   # u(t) = t
    du = ones(1, N)            # slope 1
    Δt = fill(Δt_val, N)

    traj = NamedTrajectory((; u = u, du = du, Δt = Δt); timestep = :Δt, controls = (:u,))

    reg = HermiteBendingEnergyRegularizer(traj; R = 1.0)
    @test DirectTrajOpt.objective_value(reg, traj) < 1e-20
end

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

    traj =
        NamedTrajectory((ψ̃1 = ψ̃1, ψ̃2 = ψ̃2, u = u, Δt = Δt); timestep = :Δt, controls = :u)

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

    # Per-state weights must reach the objective value.
    # Asymmetric states: ⟨ψ1|ψ̃1⟩ = 1, ⟨ψ0|ψ̃2⟩ = ½
    ψ̃1_asym = zeros(ket_dim, N)
    ψ̃2_asym = zeros(ket_dim, N)
    for k = 1:N
        ψ̃1_asym[:, k] = ket_to_iso(ψ1)
        ψ̃2_asym[:, k] = ket_to_iso(0.5 * ψ0)
    end

    traj_asym = NamedTrajectory(
        (ψ̃1 = ψ̃1_asym, ψ̃2 = ψ̃2_asym, u = u, Δt = Δt);
        timestep = :Δt,
        controls = :u,
    )

    obj_w1 = CoherentKetInfidelityObjective(
        goals,
        [:ψ̃1, :ψ̃2],
        traj_asym;
        Q = 100.0,
        weights = [0.9, 0.1],
    )
    obj_w2 = CoherentKetInfidelityObjective(
        goals,
        [:ψ̃1, :ψ̃2],
        traj_asym;
        Q = 100.0,
        weights = [0.1, 0.9],
    )

    # F = |0.9·1 + 0.1·½|² = 0.9025  and  |0.1·1 + 0.9·½|² = 0.3025
    @test objective_value(obj_w1, traj_asym) ≈ 100.0 * (1 - 0.9025)
    @test objective_value(obj_w2, traj_asym) ≈ 100.0 * (1 - 0.3025)

    # Two different weight vectors give two different objective values
    @test objective_value(obj_w1, traj_asym) != objective_value(obj_w2, traj_asym)

    # Uniform weights leave the unweighted value exactly where it is
    obj_unweighted = CoherentKetInfidelityObjective(goals, [:ψ̃1, :ψ̃2], traj_asym; Q = 100.0)
    obj_uniform = CoherentKetInfidelityObjective(
        goals,
        [:ψ̃1, :ψ̃2],
        traj_asym;
        Q = 100.0,
        weights = [0.5, 0.5],
    )
    @test objective_value(obj_uniform, traj_asym) ===
          objective_value(obj_unweighted, traj_asym)

    # Weighted objectives stay differentiable
    ∇_w = zeros(traj_asym.dim * traj_asym.N + traj_asym.global_dim)
    gradient!(∇_w, obj_w1, traj_asym)
    @test !all(∇_w .== 0)
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

@testitem "weighted coherent_ket_fidelity" begin
    using LinearAlgebra
    const coherent_ket_fidelity = Piccolo.QuantumObjectives.coherent_ket_fidelity

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    goals = [ψ1, ψ0]

    # Overlaps: ⟨ψ1|ψ1⟩ = 1, ⟨ψ0|½ψ0⟩ = ½
    ψ̃s = [ket_to_iso(ψ1), ket_to_iso(0.5 * ψ0)]

    # Weighted mean of overlaps, normalized by the weight sum
    @test coherent_ket_fidelity(ψ̃s, goals; weights = [0.9, 0.1]) ≈ abs2(0.9 * 1 + 0.1 * 0.5)
    @test coherent_ket_fidelity(ψ̃s, goals; weights = [0.1, 0.9]) ≈ abs2(0.1 * 1 + 0.9 * 0.5)

    # Distinct weights give distinct fidelities
    @test coherent_ket_fidelity(ψ̃s, goals; weights = [0.9, 0.1]) !=
          coherent_ket_fidelity(ψ̃s, goals; weights = [0.1, 0.9])

    # Weights need not be normalized — only their ratio matters
    @test coherent_ket_fidelity(ψ̃s, goals; weights = [9.0, 1.0]) ≈
          coherent_ket_fidelity(ψ̃s, goals; weights = [0.9, 0.1])

    # Uniform weights are the identity case: bit-for-bit equal to the unweighted value,
    # including when 1/n is not exactly representable
    ψ̃s3 = [ket_to_iso(ψ1), ket_to_iso(0.5 * ψ0), ket_to_iso(0.25 * ψ1)]
    goals3 = [ψ1, ψ0, ψ1]
    F_plain = coherent_ket_fidelity(ψ̃s3, goals3)
    @test coherent_ket_fidelity(ψ̃s3, goals3; weights = nothing) === F_plain
    @test coherent_ket_fidelity(ψ̃s3, goals3; weights = fill(1 / 3, 3)) === F_plain
    @test coherent_ket_fidelity(ψ̃s3, goals3; weights = fill(1.0, 3)) === F_plain
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

    # Per-state weights must reach the free-phase objective too (issue #263).
    # Asymmetric states at θ=0: ⟨ψ1|ψ̃1⟩ = 1, ⟨ψ0|ψ̃2⟩ = ½
    ψ̃1_asym = zeros(ket_dim, N)
    ψ̃2_asym = zeros(ket_dim, N)
    for k = 1:N
        ψ̃1_asym[:, k] = ket_to_iso(ψ1)
        ψ̃2_asym[:, k] = ket_to_iso(0.5 * ψ0)
    end
    traj_asym = NamedTrajectory(
        (ψ̃1 = ψ̃1_asym, ψ̃2 = ψ̃2_asym, u = randn(1, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
        global_data = [0.0],
        global_components = (φ_1 = 1:1,),
    )

    free_phase_value(ws) = objective_value(
        CoherentKetFreePhaseInfidelityObjective(
            goals_fn,
            [:ψ̃1, :ψ̃2],
            θ_names,
            traj_asym;
            Q = 100.0,
            weights = ws,
        ),
        traj_asym,
    )

    # Goals are phase-rotated by θ=0 here, so the weighted mean is analytic:
    # F = |0.9·1 + 0.1·½|² = 0.9025  and  |0.1·1 + 0.9·½|² = 0.3025
    @test free_phase_value([0.9, 0.1]) ≈ 100.0 * (1 - 0.9025)
    @test free_phase_value([0.1, 0.9]) ≈ 100.0 * (1 - 0.3025)
    @test free_phase_value([0.9, 0.1]) != free_phase_value([0.1, 0.9])

    # Omitted and uniform weights leave today's value exactly where it is
    @test free_phase_value(nothing) === free_phase_value([0.5, 0.5])
    @test free_phase_value(nothing) === free_phase_value([1.0, 1.0])

    # Weighting composes with the phase rotation rather than replacing it:
    # a nonzero phase still moves a weighted objective
    traj_phase = NamedTrajectory(
        (ψ̃1 = ψ̃1_asym, ψ̃2 = ψ̃2_asym, u = randn(1, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
        global_data = [0.7],
        global_components = (φ_1 = 1:1,),
    )
    obj_phase_w = CoherentKetFreePhaseInfidelityObjective(
        goals_fn,
        [:ψ̃1, :ψ̃2],
        θ_names,
        traj_phase;
        Q = 100.0,
        weights = [0.9, 0.1],
    )
    @test objective_value(obj_phase_w, traj_phase) != free_phase_value([0.9, 0.1])

    # Weighted free-phase objectives stay differentiable
    ∇_w = zeros(traj_asym.dim * traj_asym.N + traj_asym.global_dim)
    gradient!(
        ∇_w,
        CoherentKetFreePhaseInfidelityObjective(
            goals_fn,
            [:ψ̃1, :ψ̃2],
            θ_names,
            traj_asym;
            Q = 100.0,
            weights = [0.9, 0.1],
        ),
        traj_asym,
    )
    @test !all(∇_w .== 0)
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
