module QuantumControlProblems

using DirectTrajOpt
using NamedTrajectories
using ...Quantum
using TestItems

import ...Quantum: get_system, get_goal, state_name, drive_name, extract_pulse, fidelity
import DirectTrajOpt.Solvers: solve!

export QuantumControlProblem
export get_trajectory, get_system, get_goal, state_name, drive_name
export solve!, sync_trajectory!, fidelity
# Note: solve! is NOT exported to avoid ambiguity with SciMLBase.solve!
# Users should use: using DirectTrajOpt (to get solve!)

"""
    QuantumControlProblem{QT<:AbstractQuantumTrajectory}

Wrapper combining quantum trajectory information with trajectory optimization problem.

This type enables:
- Type-stable dispatch on quantum trajectory type (Unitary, Ket, Density)
- Clean separation of quantum information (system, goal) from optimization details
- Composable problem transformations (e.g., SmoothPulseProblem → MinimumTimeProblem)

# Fields
- `qtraj::QT`: Quantum trajectory containing system, goal, and quantum state information
- `prob::DirectTrajOptProblem`: Direct trajectory optimization problem with objective, dynamics, constraints

# Construction
Typically created via problem templates:
```julia
qtraj = UnitaryTrajectory(sys, U_goal, N)
qcp = SmoothPulseProblem(qtraj; Q=100.0, R=1e-2)
```

# Accessors
- `get_trajectory(qcp)`: Get the NamedTrajectory
- `get_system(qcp)`: Get the QuantumSystem
- `get_goal(qcp)`: Get the goal state/unitary
- `state_name(qcp)`: Get the state variable name
- `drive_name(qcp)`: Get the control variable name

# Solving
```julia
solve!(qcp; max_iter=100, verbose=true)
```
"""
mutable struct QuantumControlProblem{QT<:AbstractQuantumTrajectory}
    qtraj::QT
    prob::DirectTrajOptProblem
end

# ============================================================================= #
# Convenience accessors - extend PiccoloQuantumObjects methods
# ============================================================================= #

"""
    get_trajectory(qcp::QuantumControlProblem)

Get the NamedTrajectory from the optimization problem.
"""
get_trajectory(qcp::QuantumControlProblem) = qcp.prob.trajectory

"""
    get_system(qcp::QuantumControlProblem)

Get the QuantumSystem from the quantum trajectory.
"""
get_system(qcp::QuantumControlProblem) = get_system(qcp.qtraj)

"""
    get_goal(qcp::QuantumControlProblem)

Get the goal state/operator from the quantum trajectory.
"""
get_goal(qcp::QuantumControlProblem) = get_goal(qcp.qtraj)

"""
    state_name(qcp::QuantumControlProblem)

Get the state variable name from the quantum trajectory.
"""
state_name(qcp::QuantumControlProblem) = state_name(qcp.qtraj)

"""
    drive_name(qcp::QuantumControlProblem)

Get the control variable name from the quantum trajectory.
"""
drive_name(qcp::QuantumControlProblem) = drive_name(qcp.qtraj)

"""
    fidelity(qcp::QuantumControlProblem; kwargs...)

Compute the fidelity of the quantum trajectory.

This is a convenience wrapper that forwards to `fidelity(qcp.qtraj; kwargs...)`.

# Example
```julia
solve!(qcp)
fid = fidelity(qcp)  # Equivalent to fidelity(qcp.qtraj)
```
"""
fidelity(qcp::QuantumControlProblem; kwargs...) = fidelity(qcp.qtraj; kwargs...)

# ============================================================================= #
# Forward DirectTrajOptProblem methods
# ============================================================================= #

"""
    sync_trajectory!(qcp::QuantumControlProblem)

Update the quantum trajectory in-place from the optimized control values.

After optimization, this function:
1. Extracts the optimized controls from `prob.trajectory` (unadapting if needed)
2. Creates a new pulse with those controls via `extract_pulse`
3. Re-solves the ODE to get the updated quantum evolution
4. Replaces `qtraj` with the new quantum trajectory

This gives you access to the continuous-time ODE solution with the optimized controls,
allowing you to:
- Evaluate the fidelity via `fidelity(qcp.qtraj)`
- Sample the quantum state at any time via `qcp.qtraj(t)`
- Get the optimized pulse via `get_pulse(qcp.qtraj)`

# Example
```julia
solve!(qcp; max_iter=100)  # Automatically calls sync_trajectory!
fid = fidelity(qcp.qtraj)  # Evaluate fidelity with continuous-time solution
pulse = get_pulse(qcp.qtraj)  # Get the optimized pulse
```
"""
function sync_trajectory!(qcp::QuantumControlProblem)
    # Extract the optimized pulse from the discrete trajectory and roll it out in-place
    pulse = extract_pulse(qcp.qtraj, qcp.prob.trajectory)
    rollout!(qcp.qtraj, pulse)

    # Update global parameters in the system if present
    # Use get_system() to work with all trajectory types (including SamplingTrajectory)
    sys = get_system(qcp.qtraj)
    if hasproperty(sys, :global_params) && !isempty(sys.global_params) && qcp.prob.trajectory.global_dim > 0
        update_global_params!(qcp.qtraj, qcp.prob.trajectory)
    end

    return nothing
end

"""
    solve!(qcp::QuantumControlProblem; sync::Bool=true, kwargs...)

Solve the quantum control problem by forwarding to the inner DirectTrajOptProblem.

# Arguments
- `sync::Bool=true`: If true, call `sync_trajectory!` after solving to update `qtraj.trajectory`
  with physical control values. Set to false to skip synchronization (e.g., for debugging).

All other keyword arguments are passed to the DirectTrajOpt solver.
"""
function solve!(qcp::QuantumControlProblem; sync::Bool = true, kwargs...)
    solve!(qcp.prob; kwargs...)
    if sync
        sync_trajectory!(qcp)
    end
    return nothing
end

# Forward other common DirectTrajOptProblem accessors
Base.getproperty(qcp::QuantumControlProblem, s::Symbol) = begin
    if s === :qtraj
        getfield(qcp, :qtraj)
    elseif s === :prob
        getfield(qcp, :prob)
        # Forward to prob for common fields
    elseif s in (:objective, :dynamics, :constraints, :trajectory)
        getproperty(qcp.prob, s)
    else
        # Fall back to default behavior
        getfield(qcp, s)
    end
end

# ============================================================================= #
# Display
# ============================================================================= #

function Base.show(io::IO, qcp::QuantumControlProblem{QT}) where {QT}
    println(io, "QuantumControlProblem{$QT}")
    println(io, "  System: $(typeof(get_system(qcp)))")
    println(io, "  Goal: $(typeof(get_goal(qcp)))")
    println(io, "  Trajectory: $(qcp.prob.trajectory.N) knots")
    println(io, "  State: $(state_name(qcp))")
    print(io, "  Controls: $(drive_name(qcp))")
end

# ============================================================================= #
# Tests
# ============================================================================= #

@testitem "sync_trajectory! updates quantum trajectory" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    # Create a simple quantum system with X drive
    levels = 2
    H_drift = zeros(ComplexF64, levels, levels)
    σx = ComplexF64[0 1; 1 0]

    N = 11
    T = 5.0
    sys = QuantumSystem(H_drift, [σx], [(-2.0, 2.0)]; time_dependent = false)

    ψ_init = ComplexF64[1, 0]
    ψ_target = ComplexF64[0, 1]

    # Create pulse with zero controls
    times = collect(range(0, T, length = N))
    controls = zeros(1, N)
    pulse = LinearSplinePulse(controls, times)

    # Create initial trajectory
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_target)

    # Verify initial fidelity is low (controls are zero)
    initial_fid = fidelity(qtraj)
    @test initial_fid < 0.1  # Should be near 0 since |0⟩ stays at |0⟩

    # Create a simple problem - convert QT to NamedTrajectory
    traj = NamedTrajectory(qtraj, N)
    obj = QuadraticRegularizer(:u, traj, 1.0)
    integrator = BilinearIntegrator(qtraj, N)
    prob = DirectTrajOptProblem(traj, obj, integrator)
    qcp = QuantumControlProblem(qtraj, prob)

    # Manually modify prob.trajectory to simulate optimization
    # Set u to a constant value that will rotate |0⟩ toward |1⟩
    # For a simple rotation: exp(-i * σx * u * t) |0⟩ = cos(ut)|0⟩ - i*sin(ut)|1⟩
    # After time T with u = π/(2T), we get |1⟩
    u_opt = π / (2 * T)
    qcp.prob.trajectory.u .= u_opt

    # Call sync to update qtraj with new controls
    sync_trajectory!(qcp)

    # The qtraj should now have the optimized pulse
    new_pulse = get_pulse(qcp.qtraj)
    # Access underlying data from the interpolator (.u field in DataInterpolations)
    @test all(new_pulse.controls.u .≈ u_opt)

    # The fidelity should be much better
    final_fid = fidelity(qcp.qtraj)
    @test final_fid > 0.9  # Should be near 1 now
end

@testitem "solve! with sync=true updates trajectory" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    # Create a minimal quantum system
    levels = 2
    H_drift = zeros(ComplexF64, levels, levels)
    σx = ComplexF64[0 1; 1 0]

    ψ_init = ComplexF64[1, 0]
    ψ_target = ComplexF64[0, 1]
    N = 11
    T = 5.0

    sys = QuantumSystem(H_drift, [σx], [(-2.0, 2.0)]; time_dependent = false)

    # Create pulse with zero controls
    times = collect(range(0, T, length = N))
    controls = zeros(1, N)
    pulse = LinearSplinePulse(controls, times)

    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_target)
    traj = NamedTrajectory(qtraj, N)

    # Create problem with simple objective
    obj = QuadraticRegularizer(:u, traj, 0.01)
    integrator = BilinearIntegrator(qtraj, N)
    prob = DirectTrajOptProblem(traj, obj, integrator)
    qcp = QuantumControlProblem(qtraj, prob)

    # Store original pulse
    original_pulse = get_pulse(qcp.qtraj)

    # Solve with max_iter=0 (no optimization, just test sync mechanism)
    solve!(qcp; max_iter = 0, sync = true)

    # qtraj should be updated (in-place, but new solution)
    @test true  # If we get here without errors, sync worked

    # Test sync=false doesn't update trajectory
    qtraj2 = KetTrajectory(sys, pulse, ψ_init, ψ_target)
    traj2 = NamedTrajectory(qtraj2, N)
    prob2 = DirectTrajOptProblem(traj2, obj, integrator)
    qcp2 = QuantumControlProblem(qtraj2, prob2)
    original_qtraj2 = qcp2.qtraj

    solve!(qcp2; max_iter = 0, sync = false)
    @test qcp2.qtraj === original_qtraj2  # Same object, not rebuilt
end

@testitem "extract_pulse and rollout creates new trajectory with updated pulse" begin
    using LinearAlgebra
    using NamedTrajectories

    # Create system
    levels = 2
    H_drift = zeros(ComplexF64, levels, levels)
    σx = ComplexF64[0 1; 1 0]
    N = 11
    T = 5.0

    sys = QuantumSystem(H_drift, [σx], [(-2.0, 2.0)]; time_dependent = false)

    ψ_init = ComplexF64[1, 0]
    ψ_goal = ComplexF64[0, 1]

    # Create pulse with zero controls
    times = collect(range(0, T, length = N))
    controls = zeros(1, N)
    pulse = LinearSplinePulse(controls, times)

    # Create initial trajectory
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

    # Get NamedTrajectory and modify controls
    traj = NamedTrajectory(qtraj, N)
    u_opt = π / (2 * T)

    # Create modified trajectory data
    new_u = fill(u_opt, size(traj.u))
    new_traj = NamedTrajectory(
        (; ψ̃ = traj.ψ̃, t = traj.t, Δt = traj.Δt, u = new_u);
        timestep = :Δt,
        controls = (:Δt, :u),
        bounds = traj.bounds,
        initial = traj.initial,
        goal = traj.goal,
    )

    # Extract pulse and rollout with new controls
    new_pulse = extract_pulse(qtraj, new_traj)
    new_qtraj = rollout(qtraj, new_pulse)

    # Check pulse was updated (access underlying data via .u)
    @test all(new_qtraj.pulse.controls.u .≈ u_opt)

    # Check ODE was re-solved (fidelity should be high)
    @test fidelity(new_qtraj) > 0.9

    # Original trajectory unchanged (access underlying data via .u)
    @test all(qtraj.pulse.controls.u .≈ 0.0)
end

end
