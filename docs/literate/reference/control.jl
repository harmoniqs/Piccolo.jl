# # [Control Module API](@id control-api)
#
# The Control module provides the optimal control framework: problems, objectives, constraints, integrators, and problem templates.
# This guide demonstrates the API with runnable examples.

# ## Setup
#
# First, load the required packages:

using Piccolo
using CairoMakie
using SparseArrays # hide
using Random
Random.seed!(42)

# # Quantum Control Problems
#
# ## QuantumControlProblem
#
# The main wrapper type that combines a quantum trajectory with an optimization problem.
# Let's create a complete example:

## Define a simple qubit system
H_drift = 0.5 * PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
drive_bounds = [1.0, 1.0]
sys = QuantumSystem(H_drift, H_drives, drive_bounds)

## Create initial pulse
T = 10.0
N = 50
times = collect(range(0, T, length = N))
initial_controls = 0.1 * randn(2, N)
pulse = ZeroOrderPulse(initial_controls, times)

## Define goal and trajectory
U_goal = GATES[:X]
qtraj = UnitaryTrajectory(sys, pulse, U_goal)

# # Problem Templates
#
# ## SmoothPulseProblem
#
# For piecewise constant controls (`ZeroOrderPulse`). Adds discrete derivative variables for smoothness:

qcp = SmoothPulseProblem(
    qtraj,
    N;
    Q = 100.0,        # Fidelity weight
    R = 1e-2,         # Regularization weight
    ddu_bound = 1.0,   # Control acceleration bound
)

println("Problem type: ", typeof(qcp))
println("Initial fidelity: ", round(fidelity(qcp), digits = 4))

# ## Accessing Problem Components

traj = get_trajectory(qcp)
system = get_system(qcp)
goal = get_goal(qcp)

# The trajectory has state and control variables. System info:
system.levels

# The goal operator:
goal |> sparse

# ## Solving the Problem
#
# Use the `solve!` function to optimize:

cached_solve!(qcp, "control_ref_basic"; max_iter = 20, verbose = false, print_level = 1)

# After optimization:

fidelity(qcp)


# Visualize the optimized controls:

traj = get_trajectory(qcp)
plot_times = cumsum([0; get_timesteps(traj)])[1:(end-1)]

fig = Figure(size = (800, 400))
ax = Axis(
    fig[1, 1],
    xlabel = "Time",
    ylabel = "Control Amplitude",
    title = "Optimized Controls",
)
lines!(ax, plot_times, traj[:u][1, :], label = "u_x", linewidth = 2)
lines!(ax, plot_times, traj[:u][2, :], label = "u_y", linewidth = 2)
axislegend(ax, position = :rt)
fig

# # Problem Composition
#
# ## MinimumTimeProblem
#
# Wrap an existing problem to minimize gate duration while maintaining fidelity:

## Start with a basic problem for a shorter gate
pulse_min = ZeroOrderPulse(0.1 * randn(2, 30), collect(range(0, 5.0, length = 30)))
qtraj2 = UnitaryTrajectory(sys, pulse_min, GATES[:X])
qcp_base = SmoothPulseProblem(qtraj2, 30; Q = 100.0, R = 1e-2, ddu_bound = 1.0)

## Wrap in minimum time problem
qcp_mintime = MinimumTimeProblem(
    qcp_base;
    final_fidelity = 0.99,  # Target fidelity
    D = 100.0,                # Time penalty weight
)

# ## SamplingProblem
#
# For robust optimization over multiple system variants.
# This is useful for handling parameter uncertainties:

## Create system variants with different anharmonicities
sys1 = TransmonSystem(levels = 3, δ = -0.18, drive_bounds = [0.2, 0.2])
sys2 = TransmonSystem(levels = 3, δ = -0.20, drive_bounds = [0.2, 0.2])
sys3 = TransmonSystem(levels = 3, δ = -0.22, drive_bounds = [0.2, 0.2])
system_variants = [sys1, sys2, sys3]

## For 3-level systems, use embedded X gate on qubit subspace
X_embedded = EmbeddedOperator(GATES[:X], sys2)

## Create base problem
pulse_short = ZeroOrderPulse(0.1 * randn(2, 30), collect(range(0, 5.0, length = 30)))
qtraj_robust = UnitaryTrajectory(sys2, pulse_short, X_embedded.operator)
qcp_single = SmoothPulseProblem(qtraj_robust, 30; Q = 100.0)

## Wrap in sampling problem for robust optimization
qcp_robust = SamplingProblem(qcp_single, system_variants; Q = 100.0)

# The robust optimization problem uses multiple system samples:

println("Number of system samples: ", length(system_variants)) # nothing

# # PiccoloOptions
#
# Configuration for solver behavior and additional features:

opts = PiccoloOptions(
    leakage_constraint = true,
    leakage_constraint_value = 1e-2,
    leakage_cost = 1e-2,
    verbose = true,
)

# The options configured:

println("Leakage constraint: ", opts.leakage_constraint)
println("Leakage value: ", opts.leakage_constraint_value) # nothing

# # Objectives
#
# ## UnitaryInfidelityObjective
#
# The primary objective for gate synthesis measures how far the evolved unitary is from the target.
# This is automatically included in problem templates, but you can access it:

# Example: Creating an infidelity objective manually
# (This is normally done automatically by problem templates)
manual_traj = get_trajectory(qcp)
obj = UnitaryInfidelityObjective(goal, :Ũ⃗, manual_traj; Q = 100.0)

# ## LeakageObjective
#
# Penalizes population in leakage levels for multilevel systems.
# This is added automatically when using `PiccoloOptions` with leakage suppression.

# # Constraints
#
# ## FinalUnitaryFidelityConstraint
#
# Ensures the final gate achieves a minimum fidelity:

# Example: Creating a fidelity constraint manually
# (This is normally done automatically by MinimumTimeProblem)
constraint_traj = get_trajectory(qcp)
fidelity_constraint = FinalUnitaryFidelityConstraint(goal, :Ũ⃗, 0.99, constraint_traj)

# ## FinalKetFidelityConstraint
#
# For state transfer problems, similar pattern:
ψ0 = [1.0 + 0im, 0.0 + 0im]
ψ1 = [0.0 + 0im, 1.0 + 0im]
qtraj_ket = KetTrajectory(sys, pulse, ψ0, ψ1)
qcp_ket = SmoothPulseProblem(
    qtraj_ket,
    N;
    Q = 100.0,        # Fidelity weight
    R = 1e-2,         # Regularization weight
    ddu_bound = 1.0,   # Control acceleration bound
)

ket_traj = get_trajectory(qcp_ket)
ψ_goal_constraint = ψ1
ket_constraint = FinalKetFidelityConstraint(ψ_goal_constraint, :ψ̃, 0.99, ket_traj)

# # Solving Options
#
# ## Common solve! Arguments

## Example with common options:
## solve!(qcp;
##     max_iter = 100,        # Maximum iterations
##     tol = 1e-6,            # Convergence tolerance
##     verbose = true,        # Print solver output
##     print_level = 1        # Ipopt verbosity (0-12)
## )

# # Summary
#
# The Control module provides:
#
# - **Problem Type**: `QuantumControlProblem` wraps trajectory and optimization
# - **Templates**: `SmoothPulseProblem`, `SplinePulseProblem`, `MinimumTimeProblem`, `SamplingProblem`
# - **Options**: `PiccoloOptions` for leakage suppression and other features
# - **Objectives**: Infidelity, leakage, sensitivity objectives
# - **Constraints**: Fidelity constraints, leakage constraints
# - **Solving**: `solve!` function with extensive configuration options
#
# For complete function signatures and advanced usage, see the [API Reference](@ref) overview.
