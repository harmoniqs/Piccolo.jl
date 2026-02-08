# # [Leakage Suppression](@id leakage-suppression)
#
# When working with multilevel quantum systems (like transmons), population can
# "leak" from the computational subspace to higher energy levels. This guide
# shows how to suppress leakage in Piccolo.jl.

# ## The Problem
#
# Consider a 3-level transmon where we want to implement gates only on the
# |0⟩ and |1⟩ states. During optimization, population might temporarily
# occupy |2⟩, which can:
# - Reduce gate fidelity
# - Cause errors in subsequent operations
# - Lead to non-unitary dynamics if |2⟩ decays

# ## EmbeddedOperator
#
# The key tool for handling subspace gates is `EmbeddedOperator`. It creates
# a full unitary that acts as the specified gate on the computational subspace
# and as identity on leakage levels.
#
# For a 3-level system with 2-level computational subspace:
#
# ```math
# U_{\text{embedded}} = \begin{pmatrix} U_{\text{gate}} & 0 \\ 0 & I \end{pmatrix}
# ```

# ### Basic Usage

using Piccolo
using Random
Random.seed!(42)

## Create a 3-level transmon
sys = TransmonSystem(levels = 3, δ = 0.2, drive_bounds = [0.2, 0.2])

## Define X gate in computational subspace
U_goal = EmbeddedOperator(:X, sys)

println("Full operator size: ", size(U_goal.operator))
println("Computational subspace: ", U_goal.subspace)

# ### Construction Options

## From symbol (gate name)
U_X = EmbeddedOperator(:X, sys)
U_H = EmbeddedOperator(:H, sys)

## From matrix
custom_gate = ComplexF64[1 0; 0 exp(im * π / 4)]  # T gate
U_T = EmbeddedOperator(custom_gate, sys)

# ### Accessing Subspace Information

## Computational subspace indices
println("Subspace indices: ", U_goal.subspace)

## Leakage indices in isomorphic vector space
leak_indices = get_iso_vec_leakage_indices(U_goal)
println("Leakage indices (iso vec): ", leak_indices)

# ## Leakage via PiccoloOptions
#
# The easiest way to add leakage handling is through `PiccoloOptions`.

# ### Add Leakage Objective
#
# Penalize population in leakage states with a soft cost:

opts_cost = PiccoloOptions(
    leakage_cost = 10.0,  # Weight on leakage penalty
    verbose = false,
)

# ### Add Leakage Constraint
#
# Enforce a hard bound on leakage:

opts_constraint = PiccoloOptions(
    leakage_constraint = true,
    leakage_constraint_value = 1e-3,  # Max 0.1% leakage
    verbose = false,
)

# ### Both Together

opts = PiccoloOptions(
    leakage_cost = 10.0,
    leakage_constraint = true,
    leakage_constraint_value = 1e-3,
    verbose = false,
)

# ## Complete Example: X Gate on 3-Level Transmon

## 1. Create multilevel system
sys = TransmonSystem(levels = 3, δ = 0.2, drive_bounds = [0.2, 0.2])

## 2. Define embedded gate
U_goal = EmbeddedOperator(:X, sys)

## 3. Create trajectory
T, N = 20.0, 100
times = collect(range(0, T, length = N))
pulse = ZeroOrderPulse(0.05 * randn(2, N), times)
qtraj = UnitaryTrajectory(sys, pulse, U_goal)

## 4. Configure leakage suppression
opts = PiccoloOptions(
    leakage_cost = 10.0,
    leakage_constraint = true,
    leakage_constraint_value = 1e-3,
    verbose = false,
)

## 5. Solve
qcp = SmoothPulseProblem(qtraj, N; Q = 100.0, piccolo_options = opts)
solve!(qcp; max_iter = 150, verbose = false, print_level = 1)

println("Fidelity: ", round(fidelity(qcp), digits = 6))

# ## Manual Leakage Objectives and Constraints
#
# For more control, add leakage handling manually using `LeakageObjective`
# and `LeakageConstraint`. Both take the leakage indices in the isomorphic
# vector space:

traj = get_trajectory(qcp)
leak_indices = get_iso_vec_leakage_indices(U_goal)

## Leakage objective (soft penalty)
leak_obj = LeakageObjective(leak_indices, :Ũ⃗, traj; Qs = fill(10.0, N))

## Leakage constraint (hard bound)
leak_constraint = LeakageConstraint(1e-3, leak_indices, :Ũ⃗, traj)

# ## Strategies for Difficult Problems

# ### 1. Start Without Leakage Constraints
#
# Get a working solution first, then add constraints:

## Step 1: Optimize without leakage constraints
qcp_initial = SmoothPulseProblem(qtraj, N; Q = 100.0)
solve!(qcp_initial; max_iter = 100, verbose = false, print_level = 1)

## Step 2: Add leakage suppression
opts = PiccoloOptions(
    leakage_cost = 10.0,
    leakage_constraint = true,
    leakage_constraint_value = 1e-3,
    verbose = false,
)
qcp_leakage = SmoothPulseProblem(qtraj, N; Q = 100.0, piccolo_options = opts)
solve!(qcp_leakage; max_iter = 150, verbose = false, print_level = 1)

println("Without leakage suppression: ", round(fidelity(qcp_initial), digits = 6))
println("With leakage suppression:    ", round(fidelity(qcp_leakage), digits = 6))

# ### 2. Increase Gate Time
#
# Faster gates often have more leakage. Try a longer gate time if leakage is high.

# ### 3. Use More Levels
#
# Include more levels to capture dynamics accurately:

sys_4level = TransmonSystem(levels = 4, δ = 0.2, drive_bounds = [0.2, 0.2])

# ### 4. Reduce Drive Amplitude
#
# Lower drives can reduce leakage:

sys_low_drive = TransmonSystem(levels = 3, δ = 0.2, drive_bounds = [0.1, 0.1])

# ## Analyzing Leakage
#
# After solving, visualize the state evolution including leakage levels:

using CairoMakie

traj = get_trajectory(qcp)
fig = plot_unitary_populations(traj)

# ## See Also
#
# - [Operators](@ref operators-concept) - EmbeddedOperator details
# - [System Templates](@ref system-templates) - Creating multilevel systems
# - [Objectives](@ref objectives-concept) - LeakageObjective documentation
# - [Constraints](@ref constraints-concept) - LeakageConstraint documentation
