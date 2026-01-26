export MinimumTimeProblem

@doc raw"""
    MinimumTimeProblem(qcp::QuantumControlProblem; kwargs...)

Convert an existing quantum control problem to minimum-time optimization.

**IMPORTANT**: This function requires an existing `QuantumControlProblem` (e.g., from `SmoothPulseProblem`).
It cannot be created directly from a quantum trajectory. The workflow is:
1. Create base problem with `SmoothPulseProblem` (or similar)
2. Solve base problem to get feasible solution
3. Convert to minimum-time with `MinimumTimeProblem`

This ensures the problem starts from a good initialization and maintains solution quality
through the final fidelity constraint.

# Type Dispatch
Automatically handles different quantum trajectory types through the type parameter:
- `QuantumControlProblem{UnitaryTrajectory}` → Uses `FinalUnitaryFidelityConstraint`
- `QuantumControlProblem{KetTrajectory}` → Uses `FinalKetFidelityConstraint`
- `QuantumControlProblem{DensityTrajectory}` → Not yet implemented

The optimization problem is:

```math
\begin{aligned}
\underset{\vec{\tilde{q}}, u, \Delta t}{\text{minimize}} & \quad
J_{\text{original}}(\vec{\tilde{q}}, u) + D \sum_t \Delta t_t \\
\text{ subject to } & \quad \text{original dynamics \& constraints} \\
& F_{\text{final}} \geq F_{\text{threshold}} \\
& \quad \Delta t_{\text{min}} \leq \Delta t_t \leq \Delta t_{\text{max}} \\
\end{aligned}
```

where q represents the quantum state (unitary, ket, or density matrix).

# Arguments
- `qcp::QuantumControlProblem`: Existing quantum control problem to convert

# Keyword Arguments
- `final_fidelity::Float64=0.99`: Minimum fidelity constraint at final time
- `D::Float64=100.0`: Weight on minimum-time objective ∑Δt
- `piccolo_options::PiccoloOptions=PiccoloOptions()`: Piccolo solver options

# Returns
- `QuantumControlProblem`: New problem with minimum-time objective and fidelity constraint

# Examples
```julia
# Standard workflow
sys = QuantumSystem(H_drift, H_drives, drive_bounds)
pulse = ZeroOrderPulse(0.1 * randn(n_drives, N), collect(range(0.0, T, length=N)))
qtraj = UnitaryTrajectory(sys, pulse, U_goal)

# Step 1: Create and solve base smooth pulse problem (with Δt_bounds for free time)
qcp_smooth = SmoothPulseProblem(qtraj, N; Q=100.0, R=1e-2, Δt_bounds=(0.01, 0.5))
solve!(qcp_smooth; max_iter=100)

# Step 2: Convert to minimum-time
qcp_mintime = MinimumTimeProblem(qcp_smooth; final_fidelity=0.99, D=100.0)
solve!(qcp_mintime; max_iter=100)

# Compare durations
duration_before = sum(get_timesteps(get_trajectory(qcp_smooth)))
duration_after = sum(get_timesteps(get_trajectory(qcp_mintime)))
@assert duration_after <= duration_before

# Nested transformations also work
qcp_final = MinimumTimeProblem(
    RobustnessProblem(qcp_smooth);  # Future feature
    final_fidelity=0.95
)
```

# Convenience Constructors

You can also update the goal when creating minimum-time problem:

```julia
# Different goal for minimum-time optimization
qcp_mintime = MinimumTimeProblem(qcp_smooth; goal=U_goal_new, final_fidelity=0.98)
```
"""
function MinimumTimeProblem(
    qcp::QuantumControlProblem{QT};
    goal::Union{Nothing,AbstractPiccoloOperator,AbstractVector}=nothing,
    final_fidelity::Float64=0.99,
    D::Float64=100.0,
    piccolo_options::PiccoloOptions=PiccoloOptions()
) where {QT<:AbstractQuantumTrajectory}

    if piccolo_options.verbose
        println("    constructing MinimumTimeProblem from QuantumControlProblem{$QT}...")
        println("\tfinal fidelity constraint: $(final_fidelity)")
        println("\tminimum-time weight D: $(D)")
    end

    # Copy trajectory and constraints from original problem
    traj = deepcopy(qcp.prob.trajectory)
    constraints = deepcopy(qcp.prob.constraints)

    # Add minimum-time objective to existing objective
    J = qcp.prob.objective + MinimumTimeObjective(traj, D=D)

    # Use updated goal if provided, otherwise use original
    qtraj_for_constraint = if isnothing(goal)
        qcp.qtraj
    else
        # Create new quantum trajectory with updated goal
        _update_goal(qcp.qtraj, goal)
    end

    # Add final fidelity constraint - dispatches on QT type parameter!
    fidelity_constraint = _final_fidelity_constraint(
        qtraj_for_constraint,
        final_fidelity,
        traj
    )

    # Handle single constraint or multiple constraints
    if fidelity_constraint isa AbstractVector
        append!(constraints, fidelity_constraint)
    else
        push!(constraints, fidelity_constraint)
    end

    # Create new optimization problem with same integrators
    new_prob = DirectTrajOptProblem(
        traj,
        J,
        qcp.prob.integrators,
        constraints
    )

    # Return new QuantumControlProblem with potentially updated qtraj
    return QuantumControlProblem(qtraj_for_constraint, new_prob)
end

# ============================================================================= #
# Type-specific helper functions
# ============================================================================= #

# Helper to update goal in quantum trajectory (convenience constructor support)
function _update_goal(qtraj::UnitaryTrajectory, new_goal::AbstractPiccoloOperator)
    # Create new trajectory with updated goal, using same pulse
    return UnitaryTrajectory(
        get_system(qtraj),
        qtraj.pulse,
        new_goal
    )
end

function _update_goal(qtraj::KetTrajectory, new_goal::AbstractVector{<:Number})
    # Keep initial state and pulse, update goal
    return KetTrajectory(
        get_system(qtraj),
        qtraj.pulse,
        qtraj.initial,
        new_goal
    )
end

function _update_goal(qtraj::DensityTrajectory, new_goal::AbstractMatrix)
    return DensityTrajectory(
        get_system(qtraj),
        qtraj.pulse,
        qtraj.initial,
        new_goal
    )
end

# Fidelity constraint functions - dispatch on quantum trajectory type

function _final_fidelity_constraint(
    qtraj::UnitaryTrajectory,
    final_fidelity::Float64,
    traj::NamedTrajectory
)
    U_goal = qtraj.goal
    state_sym = state_name(qtraj)
    return FinalUnitaryFidelityConstraint(U_goal, state_sym, final_fidelity, traj)
end

function _final_fidelity_constraint(
    qtraj::KetTrajectory,
    final_fidelity::Float64,
    traj::NamedTrajectory
)
    ψ_goal = qtraj.goal
    state_sym = state_name(qtraj)
    return FinalKetFidelityConstraint(ψ_goal, state_sym, final_fidelity, traj)
end

function _final_fidelity_constraint(
    qtraj::DensityTrajectory,
    final_fidelity::Float64,
    traj::NamedTrajectory
)
    # TODO: Implement density matrix fidelity constraint when available
    throw(ArgumentError("Final fidelity constraint for DensityTrajectory not yet implemented"))
end

# ============================================================================= #
# MultiKetTrajectory Support
# ============================================================================= #

"""
    _final_fidelity_constraint(qtraj::MultiKetTrajectory, final_fidelity, traj)

Create a coherent fidelity constraint for an MultiKetTrajectory.

Uses coherent fidelity: F = |1/n ∑ᵢ ⟨ψᵢ_goal|ψᵢ⟩|²

This enforces that all state transfers have aligned global phases, which is 
essential when implementing a gate via state transfer (e.g., X gate via 
|0⟩→|1⟩ and |1⟩→|0⟩).
"""
function _final_fidelity_constraint(
    qtraj::MultiKetTrajectory,
    final_fidelity::Float64,
    traj::NamedTrajectory
)
    snames = state_names(qtraj)
    goals = qtraj.goals
    
    # Use coherent fidelity constraint for proper phase alignment
    return FinalCoherentKetFidelityConstraint(goals, snames, final_fidelity, traj)
end

function _ensemble_fidelity_constraint(
    base_qtraj::DensityTrajectory,
    goal::AbstractMatrix,
    state_sym::Symbol,
    final_fidelity::Float64,
    traj::NamedTrajectory
)
    throw(ArgumentError("Final fidelity constraint for DensityTrajectory ensemble not yet implemented"))
end

# Update goal for MultiKetTrajectory is not typically needed since goals are embedded
# in the trajectory. But we provide a fallback that errors clearly.
function _update_goal(qtraj::MultiKetTrajectory, new_goal)
    throw(ArgumentError(
        "Updating goals for MultiKetTrajectory is not directly supported. " *
        "Create a new MultiKetTrajectory with the desired goals instead."
    ))
end

# ============================================================================= #
# Tests
# ============================================================================= #

@testitem "MinimumTimeProblem from SmoothPulseProblem (Unitary)" begin
    using DirectTrajOpt
    using NamedTrajectories

    T = 1.0
    N = 50
    
    # Create and solve smooth pulse problem
    sys = QuantumSystem(0.1 * GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    pulse = ZeroOrderPulse(0.1 * randn(2, N), collect(range(0.0, T, length=N)))
    qtraj = UnitaryTrajectory(sys, pulse, GATES[:H])
    qcp_smooth = SmoothPulseProblem(qtraj, N; Q=100.0, R=1e-2, Δt_bounds=(0.01, 0.5))

    solve!(qcp_smooth; max_iter=50, verbose=false, print_level=1)
    duration_before = sum(get_timesteps(get_trajectory(qcp_smooth)))

    # Convert to minimum-time problem
    qcp_mintime = MinimumTimeProblem(qcp_smooth; final_fidelity=0.95, D=100.0)

    @test qcp_mintime isa QuantumControlProblem
    @test qcp_mintime isa QuantumControlProblem{<:UnitaryTrajectory}
    @test haskey(get_trajectory(qcp_mintime).components, :du)
    @test haskey(get_trajectory(qcp_mintime).components, :ddu)

    # Test accessors
    @test get_system(qcp_mintime) === sys
    @test qtraj.goal === GATES[:H]

    # Solve minimum-time problem
    solve!(qcp_mintime; max_iter=50, verbose=false, print_level=1)
    duration_after = sum(get_timesteps(get_trajectory(qcp_mintime)))

    # Duration should decrease (or stay same if already optimal)
    @test duration_after <= duration_before
end

@testitem "MinimumTimeProblem from SmoothPulseProblem (Ket)" begin
    using DirectTrajOpt
    using NamedTrajectories

    T = 1.0
    N = 50
    
    sys = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    pulse = ZeroOrderPulse(0.1 * randn(2, N), collect(range(0.0, T, length=N)))
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

    # Create smooth pulse problem
    qcp_smooth = SmoothPulseProblem(qtraj, N; Q=50.0, R=1e-3, Δt_bounds=(0.01, 0.5))
    solve!(qcp_smooth; max_iter=10, verbose=false, print_level=1)

    # Convert to minimum-time
    qcp_mintime = MinimumTimeProblem(qcp_smooth; final_fidelity=0.90, D=50.0)

    @test qcp_mintime isa QuantumControlProblem{<:KetTrajectory}
    @test haskey(get_trajectory(qcp_mintime).components, :du)

    # Test problem solve
    solve!(qcp_mintime; max_iter=10, print_level=1, verbose=false)
end

@testitem "MinimumTimeProblem with updated goal" begin
    using DirectTrajOpt

    T = 1.0
    N = 50
    
    sys = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    pulse = ZeroOrderPulse(0.1 * randn(2, N), collect(range(0.0, T, length=N)))
    qtraj = UnitaryTrajectory(sys, pulse, GATES[:H])

    qcp_smooth = SmoothPulseProblem(qtraj, N; Q=100.0, R=1e-2, Δt_bounds=(0.01, 0.5))
    solve!(qcp_smooth; max_iter=5, verbose=false, print_level=1)

    # Create minimum-time with different goal
    qcp_mintime = MinimumTimeProblem(
        qcp_smooth;
        goal=GATES[:X],  # Different goal!
        final_fidelity=0.95,
        D=100.0
    )

    @test qcp_mintime isa QuantumControlProblem
    @test qcp_mintime.qtraj.goal === GATES[:X]  # Goal should be updated
end

@testitem "MinimumTimeProblem type dispatch" begin
    using DirectTrajOpt

    T = 1.0
    N = 50
    
    # Test that type parameter is correct for different trajectory types
    sys = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

    # Unitary
    pulse_u = ZeroOrderPulse(0.1 * randn(2, N), collect(range(0.0, T, length=N)))
    qtraj_u = UnitaryTrajectory(sys, pulse_u, GATES[:H])
    qcp_u = SmoothPulseProblem(qtraj_u, N)
    qcp_mintime_u = MinimumTimeProblem(qcp_u)
    @test qcp_mintime_u isa QuantumControlProblem{<:UnitaryTrajectory}

    # Ket
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    pulse_k = ZeroOrderPulse(0.1 * randn(2, N), collect(range(0.0, T, length=N)))
    qtraj_k = KetTrajectory(sys, pulse_k, ψ_init, ψ_goal)
    qcp_k = SmoothPulseProblem(qtraj_k, N)
    qcp_mintime_k = MinimumTimeProblem(qcp_k)
    @test qcp_mintime_k isa QuantumControlProblem{<:KetTrajectory}
end

@testitem "MinimumTimeProblem with SamplingTrajectory" begin
    using DirectTrajOpt
    using NamedTrajectories

    T = 1.0
    N = 50
    
    # Robust minimum-time gate optimization
    sys_nominal = QuantumSystem(0.1 * GATES[:Z], [GATES[:X]], [1.0])
    sys_perturbed = QuantumSystem(0.11 * GATES[:Z], [GATES[:X]], [1.0])

    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length=N)))
    qtraj = UnitaryTrajectory(sys_nominal, pulse, GATES[:X])
    qcp = SmoothPulseProblem(qtraj, N; Q=100.0, R=1e-2, Δt_bounds=(0.01, 0.5))

    # Create sampling problem
    sampling_prob = SamplingProblem(qcp, [sys_nominal, sys_perturbed]; Q=100.0)
    solve!(sampling_prob; max_iter=50, verbose=false, print_level=1)

    duration_before = sum(get_timesteps(get_trajectory(sampling_prob)))

    # Convert to minimum-time
    mintime_prob = MinimumTimeProblem(sampling_prob; final_fidelity=0.90, D=50.0)

    @test mintime_prob isa QuantumControlProblem{<:SamplingTrajectory}
    @test mintime_prob.qtraj isa SamplingTrajectory

    # Should have fidelity constraints for each sample
    # (one per system in the sampling ensemble)

    # Solve minimum-time
    solve!(mintime_prob; max_iter=20, verbose=false, print_level=1)

    duration_after = sum(get_timesteps(get_trajectory(mintime_prob)))
    @test duration_after <= duration_before * 1.1  # Allow small tolerance
end

@testitem "MinimumTimeProblem with SamplingTrajectory (Ket)" begin
    using DirectTrajOpt

    T = 1.0
    N = 50
    
    # Robust minimum-time state transfer
    sys_nominal = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    sys_perturbed = QuantumSystem(1.1 * GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    pulse = ZeroOrderPulse(0.1 * randn(2, N), collect(range(0.0, T, length=N)))
    qtraj = KetTrajectory(sys_nominal, pulse, ψ_init, ψ_goal)

    qcp = SmoothPulseProblem(qtraj, N; Q=50.0, R=1e-3, Δt_bounds=(0.01, 0.5))

    # Create sampling problem
    sampling_prob = SamplingProblem(qcp, [sys_nominal, sys_perturbed]; Q=50.0)
    solve!(sampling_prob; max_iter=15, verbose=false, print_level=1)

    # Convert to minimum-time
    mintime_prob = MinimumTimeProblem(sampling_prob; final_fidelity=0.85, D=30.0)

    @test mintime_prob isa QuantumControlProblem{<:SamplingTrajectory}

    # Solve
    solve!(mintime_prob; max_iter=15, verbose=false, print_level=1)
end

@testitem "MinimumTimeProblem with MultiKetTrajectory" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    T = 1.0
    N = 50
    
    # Multi-state minimum-time optimization (X gate via state transfer)
    sys = QuantumSystem(0.1 * GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    
    # Create ensemble of ket states for X gate
    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    
    pulse = ZeroOrderPulse(0.1 * randn(2, N), collect(range(0.0, T, length=N)))
    ensemble_qtraj = MultiKetTrajectory(sys, pulse, [ψ0, ψ1], [ψ1, ψ0])
    
    # Create and solve smooth pulse problem
    qcp_smooth = SmoothPulseProblem(ensemble_qtraj, N; Q=100.0, R=1e-2, Δt_bounds=(0.01, 0.5))
    solve!(qcp_smooth; max_iter=30, verbose=false, print_level=1)
    
    duration_before = sum(get_timesteps(get_trajectory(qcp_smooth)))
    
    # Convert to minimum-time problem
    qcp_mintime = MinimumTimeProblem(qcp_smooth; final_fidelity=0.90, D=50.0)
    
    @test qcp_mintime isa QuantumControlProblem{<:MultiKetTrajectory}
    @test qcp_mintime.qtraj isa MultiKetTrajectory
    
    # Should have fidelity constraints for each ensemble member
    # (one per state transfer in the ensemble)
    
    # Solve minimum-time problem
    solve!(qcp_mintime; max_iter=30, verbose=false, print_level=1)
    
    duration_after = sum(get_timesteps(get_trajectory(qcp_mintime)))
    
    # Duration should not increase significantly
    @test duration_after <= duration_before * 1.1

    # Verify fidelity constraints are met for both states
    traj = get_trajectory(qcp_mintime)
    snames = state_names(qcp_mintime.qtraj)
    goals = qcp_mintime.qtraj.goals
    
    for (name, goal) in zip(snames, goals)
        ψ̃_final = traj[end][name]
        ψ_final = iso_to_ket(ψ̃_final)
        fid = fidelity(ψ_final, goal)
        @test fid >= 0.89  # Just under constraint to account for numerical tolerance
    end
end

@testitem "MinimumTimeProblem with time-dependent UnitaryTrajectory" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    # Time-dependent Hamiltonian
    ω = 2π * 5.0
    H(u, t) = GATES[:Z] + u[1] * cos(ω * t) * GATES[:X] + u[2] * sin(ω * t) * GATES[:Y]
    
    T = 5.0
    N = 50
    sys = QuantumSystem(H, [1.0, 1.0])
    
    pulse = ZeroOrderPulse(0.1 * randn(2, N), collect(range(0.0, T, length=N)))
    qtraj = UnitaryTrajectory(sys, pulse, GATES[:H])
    
    # Create and solve smooth pulse problem
    qcp_smooth = SmoothPulseProblem(qtraj, N; Q=100.0, R=1e-2, Δt_bounds=(0.01, 0.5))
    
    # TimeConsistencyConstraint is auto-applied
    @test length(qcp_smooth.prob.integrators) == 3  # dynamics + 2 derivatives
    
    solve!(qcp_smooth; max_iter=30, verbose=false, print_level=1)
    
    duration_before = sum(get_timesteps(get_trajectory(qcp_smooth)))
    
    # Convert to minimum-time
    qcp_mintime = MinimumTimeProblem(qcp_smooth; final_fidelity=0.85, D=50.0)
    
    @test qcp_mintime isa QuantumControlProblem{<:UnitaryTrajectory}
    
    # Solve minimum-time problem
    solve!(qcp_mintime; max_iter=30, verbose=false, print_level=1)
    
    duration_after = sum(get_timesteps(get_trajectory(qcp_mintime)))
    @test duration_after <= duration_before * 1.1
end

@testitem "MinimumTimeProblem with time-dependent KetTrajectory" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    # Time-dependent Hamiltonian
    ω = 2π * 5.0
    H(u, t) = GATES[:Z] + u[1] * cos(ω * t) * GATES[:X]
    
    T = 5.0
    N = 50
    sys = QuantumSystem(H, [1.0])
    
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    
    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length=N)))
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)
    
    # Create and solve smooth pulse problem
    qcp_smooth = SmoothPulseProblem(qtraj, N; Q=50.0, R=1e-3, Δt_bounds=(0.01, 0.5))
    
    # TimeConsistencyConstraint is auto-applied
    @test length(qcp_smooth.prob.integrators) == 3  # dynamics + 2 derivatives
    
    solve!(qcp_smooth; max_iter=100, verbose=false, print_level=1)
    
    duration_before = sum(get_timesteps(get_trajectory(qcp_smooth)))
    
    # Convert to minimum-time
    qcp_mintime = MinimumTimeProblem(qcp_smooth; final_fidelity=0.85, D=50.0)
    
    @test qcp_mintime isa QuantumControlProblem{<:KetTrajectory}
    
    # Solve minimum-time problem
    solve!(qcp_mintime; max_iter=30, verbose=false, print_level=1)
    
    duration_after = sum(get_timesteps(get_trajectory(qcp_mintime)))
    @test duration_after <= duration_before * 1.1
end

@testitem "MinimumTimeProblem with time-dependent MultiKetTrajectory" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    # Time-dependent Hamiltonian
    ω = 2π * 5.0
    H(u, t) = GATES[:Z] + u[1] * cos(ω * t) * GATES[:X] + u[2] * sin(ω * t) * GATES[:Y]
    
    T = 5.0
    N = 50
    sys = QuantumSystem(H, [1.0, 1.0])
    
    # Create ensemble: |0⟩ → |1⟩ and |1⟩ → |0⟩
    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    
    pulse = ZeroOrderPulse(0.1 * randn(2, N), collect(range(0.0, T, length=N)))
    qtraj = MultiKetTrajectory(sys, pulse, [ψ0, ψ1], [ψ1, ψ0])
    
    # Create and solve smooth pulse problem
    qcp_smooth = SmoothPulseProblem(qtraj, N; Q=50.0, R=1e-3, Δt_bounds=(0.01, 0.5))
    
    # TimeConsistencyConstraint is auto-applied
    # 2 dynamics + 2 derivatives = 4 integrators
    @test length(qcp_smooth.prob.integrators) == 4
    
    solve!(qcp_smooth; max_iter=30, verbose=false, print_level=1)
    
    duration_before = sum(get_timesteps(get_trajectory(qcp_smooth)))
    
    # Convert to minimum-time
    qcp_mintime = MinimumTimeProblem(qcp_smooth; final_fidelity=0.80, D=50.0)
    
    @test qcp_mintime isa QuantumControlProblem{<:MultiKetTrajectory}
    
    # Solve minimum-time problem
    solve!(qcp_mintime; max_iter=30, verbose=false, print_level=1)
    
    duration_after = sum(get_timesteps(get_trajectory(qcp_mintime)))
    @test duration_after <= duration_before * 1.1
end

@testitem "MinimumTimeProblem with time-dependent SamplingTrajectory (Unitary)" tags=[:experimental] begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    # Time-dependent Hamiltonian with oscillating drive
    ω = 2π * 5.0
    H1(u, t) = GATES[:Z] + u[1] * cos(ω * t) * GATES[:X]
    H2(u, t) = 1.1 * GATES[:Z] + u[1] * cos(ω * t) * GATES[:X]  # Perturbed
    
    T = 1.0
    N = 50
    sys_nominal = QuantumSystem(H1, [1.0])
    sys_perturbed = QuantumSystem(H2, [1.0])
    
    U_goal = GATES[:X]
    pulse = ZeroOrderPulse(randn(1, N), collect(range(0.0, T, length=N)))
    qtraj = UnitaryTrajectory(sys_nominal, pulse, U_goal)
    
    qcp = SmoothPulseProblem(qtraj, N; Q=100.0, R=1e-2, Δt_bounds=(0.01, 0.5))
    
    # Create sampling problem
    sampling_prob = SamplingProblem(qcp, [sys_nominal, sys_perturbed]; Q=100.0)
    
    @test sampling_prob isa QuantumControlProblem
    @test sampling_prob.qtraj isa SamplingTrajectory{<:AbstractPulse, <:UnitaryTrajectory}
    
    # Solve sampling problem first
    solve!(sampling_prob; max_iter=100, verbose=false, print_level=1)
    
    duration_before = sum(get_timesteps(get_trajectory(sampling_prob)))
    
    # Convert to minimum-time
    sampling_mintime = MinimumTimeProblem(sampling_prob; final_fidelity=0.60, D=50.0)
    
    @test sampling_mintime isa QuantumControlProblem{<:SamplingTrajectory}
    
    # Solve minimum-time problem
    solve!(sampling_mintime; max_iter=30, verbose=false, print_level=1)
    
    duration_after = sum(get_timesteps(get_trajectory(sampling_mintime)))
    @test duration_after <= duration_before * 1.1
end

@testitem "MinimumTimeProblem with time-dependent SamplingTrajectory (Ket)" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    # Time-dependent Hamiltonian
    ω = 2π * 5.0
    H1(u, t) = GATES[:Z] + u[1] * cos(ω * t) * GATES[:X]
    H2(u, t) = 1.1 * GATES[:Z] + u[1] * cos(ω * t) * GATES[:X]  # Perturbed
    
    T = 1.0
    N = 50
    sys_nominal = QuantumSystem(H1, [1.0])
    sys_perturbed = QuantumSystem(H2, [1.0])
    
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    
    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length=N)))
    qtraj = KetTrajectory(sys_nominal, pulse, ψ_init, ψ_goal)
    
    qcp = SmoothPulseProblem(qtraj, N; Q=50.0, R=1e-3, Δt_bounds=(0.01, 0.5))
    
    # Create sampling problem
    sampling_prob = SamplingProblem(qcp, [sys_nominal, sys_perturbed]; Q=50.0)
    
    @test sampling_prob isa QuantumControlProblem
    @test sampling_prob.qtraj isa SamplingTrajectory{<:AbstractPulse, <:KetTrajectory}
    
    # Solve sampling problem first
    solve!(sampling_prob; max_iter=100, verbose=false, print_level=1)
    
    duration_before = sum(get_timesteps(get_trajectory(sampling_prob)))
    
    # Convert to minimum-time
    sampling_mintime = MinimumTimeProblem(sampling_prob; final_fidelity=0.60, D=50.0)
    
    @test sampling_mintime isa QuantumControlProblem{<:SamplingTrajectory}
    
    # Solve minimum-time problem
    solve!(sampling_mintime; max_iter=30, verbose=false, print_level=1)
    
    duration_after = sum(get_timesteps(get_trajectory(sampling_mintime)))
    @test duration_after <= duration_before * 1.1
end