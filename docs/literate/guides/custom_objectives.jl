# # [Custom Objectives](@id custom-objectives)
#
# While Piccolo.jl provides standard objectives for fidelity and regularization,
# you may need custom objectives for specialized optimization goals. This guide
# shows how to create and use them.

# ## Setup
#
# We'll work with a simple single-qubit X gate problem:

using Piccolo
using DirectTrajOpt
using LinearAlgebra
using Random
Random.seed!(42)

## Define system
H_drift = 0.5 * PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
drive_bounds = [1.0, 1.0]
sys = QuantumSystem(H_drift, H_drives, drive_bounds)

## Create initial pulse and trajectory
T = 10.0
N = 50
times = collect(range(0, T, length = N))
initial_controls = 0.1 * randn(2, N)
pulse = ZeroOrderPulse(initial_controls, times)

U_goal = GATES[:X]
qtraj = UnitaryTrajectory(sys, pulse, U_goal)

# ## Objective Types
#
# Piccolo uses two objective types from DirectTrajOpt.jl:
#
# | Type | When Evaluated | Constructor |
# |------|----------------|-------------|
# | `TerminalObjective` | At final timestep only | `TerminalObjective(loss, name, traj; Q=weight)` |
# | `KnotPointObjective` | At specified timesteps | `KnotPointObjective(loss, name, traj; Qs=weights)` |
#
# Both take a **loss function** as their first argument — a regular Julia function
# that maps a component vector to a scalar cost.

# ## Creating a Custom Terminal Objective
#
# A terminal objective evaluates a loss function on a trajectory component at the
# final timestep. Let's create one that penalizes trace distance from the target:

## Build the problem to get the internal NamedTrajectory
qcp = SmoothPulseProblem(qtraj, N; Q = 100.0, R = 1e-2, ddu_bound = 1.0)
traj = get_trajectory(qcp)

## Define a custom loss: trace distance penalty on the isomorphic unitary vector
trace_distance_loss(Ũ⃗) = let
    U = iso_vec_to_operator(Ũ⃗)
    n = size(U, 1)
    diff = U - U_goal
    real(tr(diff' * diff)) / n
end

## Wrap it in a TerminalObjective (evaluated on :Ũ⃗ at the final timestep)
custom_terminal_obj = TerminalObjective(trace_distance_loss, :Ũ⃗, traj; Q = 50.0)

# ## Creating a Custom Knotpoint Objective
#
# A knotpoint objective evaluates at every (or selected) timesteps.
# Let's penalize control energy:

## Loss function: squared norm of control values
control_energy_loss(u) = dot(u, u)

## Applied at all timesteps on the :u component
custom_knotpoint_obj = KnotPointObjective(
    control_energy_loss, :u, traj;
    Qs = fill(0.1, N),
    times = 1:N
)

# ## Adding Custom Objectives to a Problem
#
# ### Building the objective manually
#
# Combine objective terms with `+`:

## Start with the standard infidelity objective
J = UnitaryInfidelityObjective(U_goal, :Ũ⃗, traj; Q = 100.0)

## Add standard regularization
J += QuadraticRegularizer(:u, traj, 1e-2)
J += QuadraticRegularizer(:du, traj, 1e-2)
J += QuadraticRegularizer(:ddu, traj, 1e-2)

## Add our custom objective
J += custom_knotpoint_obj

typeof(J)

# ### Adding to an existing problem

qcp = SmoothPulseProblem(qtraj, N; Q = 100.0, R = 1e-2, ddu_bound = 1.0)
qcp.prob.objective += custom_knotpoint_obj

# ## Solving and Comparing
#
# Solve with the extra control energy penalty:

cached_solve!(qcp, "custom_objectives_with_penalty"; max_iter = 50, verbose = false, print_level = 1)
fidelity(qcp)

# Compare against the standard problem:

qcp_standard = SmoothPulseProblem(qtraj, N; Q = 100.0, R = 1e-2, ddu_bound = 1.0)
cached_solve!(qcp_standard, "custom_objectives_standard"; max_iter = 50, verbose = false, print_level = 1)
fidelity(qcp_standard)

# ## Example: Leakage-Style Penalty
#
# Here's how the built-in `LeakageObjective` works under the hood —
# it's just a `KnotPointObjective` with a custom loss function:

leakage_indices = [3, 4]  # Indices of leakage states in the isomorphic vector

leakage_loss(x) = sum(abs2, x[leakage_indices]) / length(leakage_indices)

leakage_obj = KnotPointObjective(
    leakage_loss, :Ũ⃗, traj;
    Qs = fill(1.0, N),
    times = 1:N
)

# ## Tips for Custom Objectives

# ### 1. Scale Appropriately
#
# Match the scale of built-in objectives:
#
# - Built-in fidelity uses `Q ~ 100`
# - Custom objectives should use similar magnitude
#
# ```julia
# custom_obj = TerminalObjective(my_loss, :Ũ⃗, traj; Q = 50.0)
# ```

# ### 2. Ensure Smoothness
#
# Avoid discontinuities that can cause optimization issues:
#
# ```julia
# ## Bad: discontinuous
# penalty = x > 0 ? x^2 : 0
#
# ## Good: smooth approximation
# penalty = max(0, x)^2
#
# ## Also good: softplus
# penalty = log(1 + exp(k * x)) / k
# ```

# ### 3. Test Independently
#
# Verify your objective computes expected values:
#
# ```julia
# ## Evaluate the full objective on a trajectory
# value = objective_value(custom_obj, traj)
# println("Objective value: ", value)
#
# ## Check gradient
# ∇ = zeros(traj.dim * traj.N + traj.global_dim)
# gradient!(∇, custom_obj, traj)
# ```

# ### 4. Start Simple
#
# Add custom objectives incrementally:
#
# ```julia
# ## Step 1: Solve with standard objectives
# qcp = SmoothPulseProblem(qtraj, N)
# solve!(qcp)
#
# ## Step 2: Check if solution needs improvement
# println("Fidelity: ", fidelity(qcp))
#
# ## Step 3: Add custom objective with small weight
# custom_obj = KnotPointObjective(my_loss, :u, traj; Qs = fill(0.1, N))
# qcp.prob.objective += custom_obj
# solve!(qcp)
#
# ## Step 4: Increase weight if needed
# ```

# ### 5. Gradients are Automatic
#
# DirectTrajOpt uses ForwardDiff for automatic differentiation, so gradients
# are computed automatically for any loss function you provide. No need to
# implement gradient methods manually.

# ## See Also
#
# - [Objectives](@ref objectives-concept) - Built-in objectives
# - [Constraints](@ref constraints-concept) - Hard constraints (vs soft objective penalties)
# - [Problem Templates](@ref problem-templates-overview) - Using objectives in problems
