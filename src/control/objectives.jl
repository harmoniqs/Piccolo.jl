module QuantumObjectives

export KetInfidelityObjective
export CoherentKetInfidelityObjective
export UnitaryInfidelityObjective
export DensityMatrixInfidelityObjective
export DensityMatrixPureStateInfidelityObjective
export UnitarySensitivityObjective
export UnitaryFreePhaseInfidelityObjective
export LeakageObjective

using LinearAlgebra
using NamedTrajectories
using ...Quantum
using DirectTrajOpt
using TestItems

# --------------------------------------------------------- 
#                       Kets
# ---------------------------------------------------------

function ket_fidelity_loss(ψ̃::AbstractVector, ψ_goal::AbstractVector{<:Complex{Float64}})
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
function coherent_ket_fidelity(ψ̃s, ψ_goals::Vector{<:AbstractVector{<:Complex{Float64}}})
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

end
