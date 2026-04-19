module QuantumConstraints

using ..QuantumObjectives
using ..QuantumObjectives: ket_fidelity_loss, unitary_fidelity_loss, coherent_ket_fidelity

using DirectTrajOpt
using LinearAlgebra
using NamedTrajectories
using TestItems
using ...Quantum

export FinalKetFidelityConstraint
export FinalUnitaryFidelityConstraint
export FinalCoherentKetFidelityConstraint
export FinalDensityFidelityConstraint
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

"""
    FinalCoherentKetFidelityConstraint(goals_fn, ψ̃_names, θ_names, final_fidelity, traj)

Free-phase version: `goals_fn(θ)` returns phase-adjusted goal kets.
Uses `NonlinearGlobalKnotPointConstraint` to include global phase variables.
"""
function FinalCoherentKetFidelityConstraint(
    goals_fn::Function,
    ψ̃_names::Vector{Symbol},
    θ_names::AbstractVector{Symbol},
    final_fidelity::Float64,
    traj::NamedTrajectory,
)
    n_states = length(ψ̃_names)
    state_dims = [traj.dims[name] for name in ψ̃_names]
    total_state_dim = sum(state_dims)

    function terminal_constraint(z)
        x = z[1:total_state_dim]
        θ = z[(total_state_dim+1):end]

        # Extract each state from the concatenated vector
        ψ̃s = Vector{Vector{eltype(x)}}(undef, n_states)
        offset = 0
        for i = 1:n_states
            ψ̃s[i] = x[(offset+1):(offset+state_dims[i])]
            offset += state_dims[i]
        end

        phased_goals = goals_fn(θ)
        return [final_fidelity - coherent_ket_fidelity(ψ̃s, phased_goals)]
    end

    return NonlinearGlobalKnotPointConstraint(
        terminal_constraint,
        ψ̃_names,
        θ_names,
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
#                     Density Matrices
# ---------------------------------------------------------

@doc raw"""
    FinalDensityFidelityConstraint(ρ_goal, ρ̃_name, final_fidelity, traj)

Enforce a minimum state-transfer fidelity on the final density matrix of a
`DensityTrajectory` knot-point variable.

The trajectory stores the density matrix in the compact real isomorphism
representation (see `density_to_compact_iso` / `compact_iso_to_density`): a
real vector of length ``n^2`` for a Hermitian ``n \times n`` density matrix.

The constraint enforces

```math
F(\rho_{\text{final}}) = \mathrm{Re}\,\mathrm{tr}(\rho_{\text{goal}}\, \rho_{\text{final}}) \geq F_{\text{threshold}},
```

which is the exact state-transfer fidelity when ``\rho_{\text{goal}}`` is pure.
``F`` is real-linear (but not complex-linear) in the compact-iso vector, so the
Jacobian is a constant row vector; `NonlinearKnotPointConstraint` recovers it
via automatic differentiation at constraint construction.

# Arguments
- `ρ_goal::AbstractMatrix{<:Number}`: Target density matrix (``n \times n``,
  Hermitian).
- `ρ̃_name::Symbol`: Name of the compact-iso density component in the
  trajectory (e.g. `:ρ⃗̃`).
- `final_fidelity::Float64`: Minimum fidelity threshold
  (constraint: ``F \geq F_{\text{threshold}}``).
- `traj::NamedTrajectory`: The trajectory carrying the density component.
"""
function FinalDensityFidelityConstraint(
    ρ_goal::AbstractMatrix{<:Number},
    ρ̃_name::Symbol,
    final_fidelity::Float64,
    traj::NamedTrajectory,
)
    ρ_goal_c = Matrix{ComplexF64}(ρ_goal)

    function terminal_constraint(ρ̃)
        ρ = compact_iso_to_density(ρ̃)
        F = real(tr(ρ_goal_c * ρ))
        return [final_fidelity - F]
    end

    return NonlinearKnotPointConstraint(
        terminal_constraint,
        ρ̃_name,
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

# ---------------------------------------------------------
#                         Tests
# ---------------------------------------------------------

@testitem "FinalDensityFidelityConstraint via MinimumTimeProblem" begin
    using DirectTrajOpt
    using LinearAlgebra
    using NamedTrajectories

    # 2-level open system, no dissipation, σx drive: |0⟩⟨0| → |1⟩⟨1|
    sys = OpenQuantumSystem(zeros(ComplexF64, 2, 2), [PAULIS.X], [1.0])

    ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρg = ComplexF64[0.0 0.0; 0.0 1.0]

    T = 10.0
    N = 50

    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    qtraj = DensityTrajectory(sys, pulse, ρ0, ρg)

    qcp_smooth = SmoothPulseProblem(
        qtraj,
        N;
        Q = 100.0,
        R = 1e-2,
        Δt_bounds = (0.01, 0.5),
    )
    solve!(qcp_smooth; max_iter = 100, verbose = false, print_level = 1)

    # Convert to minimum-time — this is the path the dispatch stub used to block
    qcp_mintime =
        MinimumTimeProblem(qcp_smooth; final_fidelity = 0.95, D = 50.0)

    @test qcp_mintime isa QuantumControlProblem{<:DensityTrajectory}

    # Solve minimum-time problem
    solve!(qcp_mintime; max_iter = 100, verbose = false, print_level = 1)

    # Roll out with the optimized pulse and verify the fidelity constraint is respected
    traj = get_trajectory(qcp_mintime)
    ρ̃_final = traj[end][:ρ⃗̃]
    ρ_final = compact_iso_to_density(ρ̃_final)

    # Final density matrix should be Hermitian and trace-preserving
    @test ρ_final ≈ ρ_final' atol = 1e-6
    @test real(tr(ρ_final)) ≈ 1.0 atol = 1e-2

    # The constraint enforces F = Re tr(ρ_goal · ρ_final) ≥ 0.95 at the final knot
    fid = real(tr(ρg * ρ_final))
    @test fid ≥ 0.94  # small tolerance for solver feasibility slack
end

@testitem "FinalDensityFidelityConstraint direct construction" begin
    using DirectTrajOpt
    using LinearAlgebra

    # Smoke test the constructor independently of the minimum-time path.
    sys = OpenQuantumSystem(zeros(ComplexF64, 2, 2), [PAULIS.X], [1.0])

    ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρg = ComplexF64[0.5 0.5; 0.5 0.5]  # |+⟩⟨+|

    T = 2.0
    N = 10
    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    qtraj = DensityTrajectory(sys, pulse, ρ0, ρg)

    qcp = SmoothPulseProblem(qtraj, N; Q = 10.0, R = 1e-2)
    traj = get_trajectory(qcp)

    constraint = FinalDensityFidelityConstraint(ρg, :ρ⃗̃, 0.9, traj)

    @test constraint isa DirectTrajOpt.AbstractNonlinearConstraint
    # Inequality constraint, evaluated only at final knot
    @test constraint.equality == false
    @test constraint.times == [traj.N]
    @test constraint.g_dim == 1
end

end
