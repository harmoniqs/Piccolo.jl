# # [Your First Gate](@id first-gate-tutorial)
#
# This tutorial walks through synthesizing your first quantum gate with Piccolo.jl.
# We'll implement an X gate (NOT gate) on a single qubit.
#
# ## What We're Doing
#
# We want to find control pulses that implement:
# ```math
# X = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}
# ```
#
# Our qubit has Hamiltonian:
# ```math
# H(t) = \frac{\omega}{2}\sigma_z + u_x(t)\sigma_x + u_y(t)\sigma_y
# ```
#
# The optimizer will find ``u_x(t)`` and ``u_y(t)`` that produce the X gate.

# ## Setup
#
# First, load the required packages:

using Piccolo
using CairoMakie
using Random
Random.seed!(42)  # For reproducibility

# ## Step 1: Define the Quantum System
#
# A `QuantumSystem` needs:
# - **Drift Hamiltonian**: Always-on terms (qubit frequency)
# - **Drive Hamiltonians**: Controllable interactions
# - **Drive bounds**: Maximum control amplitudes

## The drift Hamiltonian: ω/2 σ_z (qubit frequency)
## We set ω = 1.0 for simplicity
H_drift = 0.5 * PAULIS[:Z]

## The drive Hamiltonians: σ_x and σ_y controls
H_drives = [PAULIS[:X], PAULIS[:Y]]

## Maximum amplitude for each drive (in same units as H_drift)
drive_bounds = [1.0, 1.0]

## Create the system
sys = QuantumSystem(H_drift, H_drives, drive_bounds)

# Let's check what we created:

println("System has $(sys.levels) energy levels")
println("System has $(sys.n_drives) control drives")

# ## Step 2: Create an Initial Pulse
#
# We need an initial guess for the control pulse. `ZeroOrderPulse` represents
# piecewise constant controls - the standard choice for most problems.

## Gate duration and discretization
T = 10.0   # Total time (in units where ω = 1)
N = 100    # Number of timesteps

## Time vector
times = collect(range(0, T, length=N))

## Random initial controls (small amplitude)
## Shape: (n_drives, N) = (2, 100)
initial_controls = 0.1 * randn(2, N)

## Create the pulse
pulse = ZeroOrderPulse(initial_controls, times)

# Check the pulse:

println("Pulse duration: ", duration(pulse))
println("Number of drives: ", n_drives(pulse))
println("Control at t=5: ", pulse(5.0))

# ## Step 3: Define the Goal
#
# A `UnitaryTrajectory` combines the system, pulse, and target gate.

## Our target: the X gate
U_goal = GATES[:X]

println("Target gate:")
display(U_goal)

## Create the trajectory
qtraj = UnitaryTrajectory(sys, pulse, U_goal)

# ## Step 4: Set Up the Optimization Problem
#
# `SmoothPulseProblem` creates an optimization problem with:
# - Fidelity objective (weight `Q`)
# - Control regularization (weight `R`)
# - Smoothness via derivative bounds

qcp = SmoothPulseProblem(
    qtraj, N;
    Q = 100.0,       # Fidelity weight (higher = prioritize fidelity)
    R = 1e-2,        # Regularization weight (higher = smoother controls)
    ddu_bound = 1.0  # Limit on control acceleration
)

# ## Step 5: Solve!
#
# The `solve!` function runs the optimizer:

solve!(qcp; max_iter=20, verbose=false, print_level=1)

# ## Step 6: Analyze the Results
#
# First, check the fidelity:

println("\n=== Results ===")
println("Final fidelity: ", round(fidelity(qcp), digits=6))

# Get the optimized trajectory:

traj = get_trajectory(qcp)

## Check the final unitary
U_final = iso_vec_to_operator(traj[:Ũ⃗][:, end])
println("\nFinal unitary (should be close to X gate):")
display(round.(U_final, digits=3))

# ## Step 7: Visualize
#
# Plot the optimized control pulses:

fig = Figure(size=(800, 400))

## Time axis
plot_times = cumsum([0; get_timesteps(traj)])[1:end-1]

## Control pulses
ax1 = Axis(fig[1, 1], xlabel="Time", ylabel="Control Amplitude", title="Optimized Controls")
lines!(ax1, plot_times, traj[:u][1, :], label="u_x (σ_x drive)", linewidth=2)
lines!(ax1, plot_times, traj[:u][2, :], label="u_y (σ_y drive)", linewidth=2)
axislegend(ax1, position=:rt)

fig

# ## Understanding the Solution
#
# The optimizer found control pulses that:
# 1. **Start and end smoothly** (due to derivative regularization)
# 2. **Stay within bounds** (due to drive_bounds)
# 3. **Achieve high fidelity** (due to the Q-weighted objective)
#
# The X gate rotates the qubit state around the X-axis by π radians.
# You can see the controls create the right rotation!

# ## What's Next?
#
# Now that you've synthesized your first gate, try:
#
# 1. **Different gates**: Change `U_goal` to `GATES[:H]` (Hadamard) or `GATES[:T]`
# 2. **Faster gates**: Reduce `T` and see how fidelity changes
# 3. **Smoother pulses**: Increase `R` or decrease `ddu_bound`
# 4. **Time-optimal**: Add `Δt_bounds` and use `MinimumTimeProblem`
#
# Continue to the [State Transfer](@ref state-transfer-tutorial) tutorial to learn about
# preparing specific quantum states.
