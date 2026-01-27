# # Quickstart Guide
#
# This guide shows you how to set up and solve a quantum optimal control problem
# in Piccolo.jl. We'll synthesize a single-qubit X gate.
#
# ## The Problem
#
# We want to find control pulses that implement an X gate on a single qubit
# with system Hamiltonian:
# ```math
# H(t) = \frac{\omega}{2} \sigma_z + u_1(t) \sigma_x + u_2(t) \sigma_y
# ```

using Piccolo

# ## Step 1: Define the Quantum System
#
# First, we define our quantum system by specifying the drift Hamiltonian (always-on),
# the drive Hamiltonians (controllable), and the bounds on control amplitudes.

## Drift Hamiltonian: qubit frequency term
H_drift = 0.5 * PAULIS[:Z]

## Drive Hamiltonians: X and Y controls
H_drives = [PAULIS[:X], PAULIS[:Y]]

## Maximum control amplitudes
drive_bounds = [1.0, 1.0]

## Create the quantum system
sys = QuantumSystem(H_drift, H_drives, drive_bounds)

# ## Step 2: Create an Initial Pulse
#
# We need an initial guess for the control pulse. `ZeroOrderPulse` represents
# piecewise constant controls.

## Time parameters
T = 10.0   # Total gate duration
N = 100    # Number of timesteps

## Create time vector
times = collect(range(0, T, length = N))

## Random initial controls (scaled by drive bounds)
initial_controls = 0.1 * randn(2, N)

## Create the pulse
pulse = ZeroOrderPulse(initial_controls, times)

# ## Step 3: Define the Goal via a Trajectory
#
# A `UnitaryTrajectory` combines the system, pulse, and target gate.

## Target: X gate
U_goal = GATES[:X]

## Create the trajectory
qtraj = UnitaryTrajectory(sys, pulse, U_goal)

# ## Step 4: Set Up the Optimization Problem
#
# `SmoothPulseProblem` creates the optimization problem with:
# - Fidelity objective (weight Q)
# - Regularization for smooth controls (weight R)
# - Derivative bounds for control smoothness

qcp = SmoothPulseProblem(
    qtraj,
    N;
    Q = 100.0,       # Fidelity weight
    R = 1e-2,        # Regularization weight
    ddu_bound = 1.0,  # Control acceleration bound
)

# ## Step 5: Solve!

solve!(qcp; max_iter = 20, verbose = false, print_level = 1)

# ## Step 6: Analyze Results
#
# After solving, we can check the fidelity and examine the optimized controls.

## Check final fidelity
println("Final fidelity: ", fidelity(qcp))

## Access the trajectory
traj = get_trajectory(qcp)
println("Control dimensions: ", size(traj[:u]))

# We can convert the isomorphic unitary back to matrix form:

U_final = iso_vec_to_operator(traj[:Ũ⃗][:, end])
println("\nFinal unitary:")
display(round.(U_final, digits = 3))

# ## Visualization
#
# Piccolo provides specialized plotting functions for quantum trajectories:

using CairoMakie

## Plot the unitary evolution (state populations over time)
fig = plot_unitary_populations(traj)

# ## Minimum Time Optimization
#
# Now let's find the shortest gate duration that achieves 99% fidelity.
#
# First, we need to create a new problem with variable timesteps enabled:

## Create problem with free-time optimization
qcp_free = SmoothPulseProblem(
    qtraj,
    N;
    Q = 100.0,
    R = 1e-2,
    ddu_bound = 1.0,
    Δt_bounds = (0.01, 0.5),  # Enable variable timesteps
)
solve!(qcp_free; max_iter = 20, verbose = false, print_level = 1)

## Convert to minimum time problem
qcp_mintime = MinimumTimeProblem(qcp_free; final_fidelity = 0.99)
solve!(qcp_mintime; max_iter = 20, verbose = false, print_level = 1)

# Compare durations:

initial_duration = sum(get_timesteps(get_trajectory(qcp_free)))
minimum_duration = sum(get_timesteps(get_trajectory(qcp_mintime)))

println("\nInitial duration:  ", round(initial_duration, digits = 3))
println("Minimum duration:  ", round(minimum_duration, digits = 3))
println("Time saved:        ", round(initial_duration - minimum_duration, digits = 3))
println("Final fidelity:    ", round(fidelity(qcp_mintime), digits = 4))

# Plot the time-optimal solution:

fig_mintime = plot_unitary_populations(get_trajectory(qcp_mintime))

# ## Next Steps
#
# - Learn about different [Problem Templates](@ref) for various optimization scenarios
# - Explore [Tutorials](@ref) for more complex examples
# - See [Concepts](@ref) for detailed documentation of types and functions
