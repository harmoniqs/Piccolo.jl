module QuantumConstraints

using ..QuantumObjectives
using ..QuantumObjectives: ket_fidelity_loss, unitary_fidelity_loss, coherent_ket_fidelity

using DirectTrajOpt
using LinearAlgebra
using NamedTrajectories
using ...Quantum

export FinalKetFidelityConstraint
export FinalUnitaryFidelityConstraint
export FinalCoherentKetFidelityConstraint
export LeakageConstraint

# ---------------------------------------------------------
#                        Kets
# ---------------------------------------------------------

function FinalKetFidelityConstraint(
    ψ_goal::AbstractVector{<:Complex{Float64}},
    ψ̃_name::Symbol,
    final_fidelity::Float64,
    traj::NamedTrajectory,
)
    terminal_constraint = ψ̃ -> [final_fidelity - ket_fidelity_loss(ψ̃, ψ_goal)]

    return NonlinearKnotPointConstraint(
        terminal_constraint,
        ψ̃_name,
        traj,
        equality = false,
        times = [traj.N],
    )
end

# ---------------------------------------------------------
#                  Coherent Ket Fidelity
# ---------------------------------------------------------

"""
    FinalCoherentKetFidelityConstraint(ψ_goals, ψ̃_names, final_fidelity, traj)

Create a final fidelity constraint using coherent ket fidelity across multiple states.

Coherent fidelity: F = |1/n ∑ᵢ ⟨ψᵢ_goal|ψᵢ⟩|²

This constraint enforces that all state overlaps have aligned phases, which is 
essential when implementing a gate via multiple state transfers (e.g., MultiKetTrajectory).

# Arguments
- `ψ_goals::Vector{<:AbstractVector{<:Complex}}`: Target ket states
- `ψ̃_names::Vector{Symbol}`: Names of isomorphic state variables in trajectory
- `final_fidelity::Float64`: Minimum fidelity threshold (constraint: F ≥ final_fidelity)
- `traj::NamedTrajectory`: The trajectory

# Example
```julia
# For implementing X gate via |0⟩→|1⟩ and |1⟩→|0⟩
goals = [ComplexF64[0, 1], ComplexF64[1, 0]]
names = [:ψ̃1, :ψ̃2]
constraint = FinalCoherentKetFidelityConstraint(goals, names, 0.99, traj)
```
"""
function FinalCoherentKetFidelityConstraint(
    ψ_goals::Vector{<:AbstractVector{<:Complex}},
    ψ̃_names::Vector{Symbol},
    final_fidelity::Float64,
    traj::NamedTrajectory,
)
    n_states = length(ψ_goals)
    @assert length(ψ̃_names) == n_states "Number of names must match number of goals"

    # Convert goals to ComplexF64
    goals = [ComplexF64.(g) for g in ψ_goals]

    # Get component info for extracting states from concatenated vector
    state_dims = [traj.dims[name] for name in ψ̃_names]

    function terminal_constraint(z_terminal)
        # Extract each state from the concatenated vector
        ψ̃s = Vector{Vector{eltype(z_terminal)}}(undef, n_states)
        offset = 0
        for i = 1:n_states
            ψ̃s[i] = z_terminal[(offset+1):(offset+state_dims[i])]
            offset += state_dims[i]
        end

        # Constraint: final_fidelity - F_coherent ≤ 0
        return [final_fidelity - coherent_ket_fidelity(ψ̃s, goals)]
    end

    return NonlinearKnotPointConstraint(
        terminal_constraint,
        ψ̃_names,
        traj,
        equality = false,
        times = [traj.N],
    )
end

# ---------------------------------------------------------
#                        Unitaries
# ---------------------------------------------------------

function FinalUnitaryFidelityConstraint(
    U_goal::AbstractPiccoloOperator,
    Ũ⃗_name::Symbol,
    final_fidelity::Float64,
    traj::NamedTrajectory,
)
    terminal_constraint = Ũ⃗ -> [final_fidelity - unitary_fidelity_loss(Ũ⃗, U_goal)]

    return NonlinearKnotPointConstraint(
        terminal_constraint,
        Ũ⃗_name,
        traj,
        equality = false,
        times = [traj.N],
    )
end

function FinalUnitaryFidelityConstraint(
    U_goal::Function,
    Ũ⃗_name::Symbol,
    θ_names::AbstractVector{Symbol},
    final_fidelity::Float64,
    traj::NamedTrajectory,
)
    θ_dim = sum(traj.global_dims[n] for n in θ_names)
    function terminal_constraint(z)
        Ũ⃗, θ = z[1:(end-θ_dim)], z[(end-θ_dim+1):end]
        return [final_fidelity - unitary_fidelity_loss(Ũ⃗, U_goal(θ))]
    end

    return NonlinearGlobalKnotPointConstraint(
        terminal_constraint,
        Ũ⃗_name,
        θ_names,
        traj,
        equality = false,
        times = [traj.N],
    )
end

# ---------------------------------------------------------
# Leakage Constraint
# ---------------------------------------------------------

"""
    LeakageConstraint(value, indices, name, traj::NamedTrajectory)

Construct a `KnotPointConstraint` that bounds leakage of `name` at the knot points specified by `times` at any `indices` that are outside the computational subspace.

"""
function LeakageConstraint(
    value::Float64,
    indices::AbstractVector{Int},
    name::Symbol,
    traj::NamedTrajectory;
    times = 1:traj.N,
)
    leakage_constraint(x) = abs2.(x[indices]) .- value

    return NonlinearKnotPointConstraint(
        leakage_constraint,
        name,
        traj,
        equality = false,
        times = times,
    )
end

end
