# [Leakage Suppression](@id leakage-suppression)

When working with multilevel quantum systems (like transmons), population can "leak" from the computational subspace to higher energy levels. This guide shows how to suppress leakage in Piccolo.jl.

## The Problem

Consider a 3-level transmon where we want to implement gates only on the |0⟩ and |1⟩ states. During optimization, population might temporarily occupy |2⟩, which can:
- Reduce gate fidelity
- Cause errors in subsequent operations
- Lead to non-unitary dynamics if |2⟩ decays

## EmbeddedOperator

The key tool for handling subspace gates is `EmbeddedOperator`.

### Basic Usage

```julia
using Piccolo

# Create a 3-level transmon
sys = TransmonSystem(levels=3, δ=0.2, drive_bounds=[0.2, 0.2])

# Define X gate in computational subspace
U_goal = EmbeddedOperator(:X, sys)
```

### How It Works

`EmbeddedOperator` creates a full unitary that:
1. Acts as the specified gate on the computational subspace
2. Acts as identity on leakage levels

For a 3-level system with 2-level computational subspace:

```math
U_{\text{embedded}} = \begin{pmatrix} U_{\text{gate}} & 0 \\ 0 & I \end{pmatrix}
```

### Construction Options

```julia
# From symbol
U_X = EmbeddedOperator(:X, sys)
U_H = EmbeddedOperator(:H, sys)
U_CZ = EmbeddedOperator(:CZ, sys)  # For 2-qubit systems

# From matrix
custom_gate = [1 0; 0 exp(im*π/4)]  # T gate
U_T = EmbeddedOperator(custom_gate, sys)

# Specifying subspace explicitly
U = EmbeddedOperator(:X, sys; subspace_indices=[1, 2])
```

### Accessing Subspace Information

```julia
U_goal = EmbeddedOperator(:X, sys)

# Computational subspace indices
comp_indices = U_goal.subspace_indices  # [1, 2]

# Leakage indices
leak_indices = U_goal.leakage_indices   # [3]

# Full operator
full_U = U_goal.operator  # 3×3 matrix
```

## Leakage via PiccoloOptions

The easiest way to add leakage handling is through `PiccoloOptions`.

### Add Leakage Objective

Penalize population in leakage states:

```julia
opts = PiccoloOptions(
    leakage_cost = 10.0  # Weight on leakage penalty
)

qcp = SmoothPulseProblem(qtraj, N; piccolo_options=opts)
```

### Add Leakage Constraint

Enforce a hard bound on leakage:

```julia
opts = PiccoloOptions(
    leakage_constraint = true,
    leakage_constraint_value = 1e-3  # Max 0.1% leakage
)

qcp = SmoothPulseProblem(qtraj, N; piccolo_options=opts)
```

### Both Together

```julia
opts = PiccoloOptions(
    leakage_cost = 10.0,           # Objective penalty
    leakage_constraint = true,      # Hard constraint
    leakage_constraint_value = 1e-3
)
```

## Complete Example

### X Gate on 3-Level Transmon

```julia
using Piccolo

# 1. Create multilevel system
sys = TransmonSystem(levels=3, δ=0.2, drive_bounds=[0.2, 0.2])

# 2. Define embedded gate
U_goal = EmbeddedOperator(:X, sys)

# 3. Create trajectory
T, N = 20.0, 100
times = collect(range(0, T, length=N))
pulse = ZeroOrderPulse(0.05 * randn(2, N), times)
qtraj = UnitaryTrajectory(sys, pulse, U_goal)

# 4. Configure leakage suppression
opts = PiccoloOptions(
    leakage_cost = 10.0,
    leakage_constraint = true,
    leakage_constraint_value = 1e-3
)

# 5. Solve
qcp = SmoothPulseProblem(qtraj, N; Q=100.0, piccolo_options=opts)
solve!(qcp; max_iter=150)

println("Fidelity: ", fidelity(qcp))
```

### Hadamard Gate with Leakage Check

```julia
using Piccolo

sys = TransmonSystem(levels=4, δ=0.2, drive_bounds=[0.2, 0.2])
U_goal = EmbeddedOperator(:H, sys)

T, N = 25.0, 100
times = collect(range(0, T, length=N))
pulse = ZeroOrderPulse(0.05 * randn(2, N), times)
qtraj = UnitaryTrajectory(sys, pulse, U_goal)

# Strong leakage suppression
opts = PiccoloOptions(
    leakage_cost = 50.0,
    leakage_constraint = true,
    leakage_constraint_value = 1e-4
)

qcp = SmoothPulseProblem(qtraj, N; Q=200.0, piccolo_options=opts)
solve!(qcp; max_iter=200)

# Check results
println("Gate fidelity: ", fidelity(qcp))
```

## Manual Leakage Objectives and Constraints

For more control, add leakage handling manually:

### LeakageObjective

```julia
# Create manually
leak_obj = LeakageObjective(:Ũ⃗, U_goal; Q=10.0)

# Add to problem via constraints parameter
qcp = SmoothPulseProblem(qtraj, N; constraints=[leak_obj])
```

### LeakageConstraint

```julia
leak_constraint = LeakageConstraint(:Ũ⃗, U_goal, 1e-3)

qcp = SmoothPulseProblem(qtraj, N; constraints=[leak_constraint])
```

## Strategies for Difficult Problems

### 1. Start Without Leakage Constraints

Get a working solution first, then add constraints:

```julia
# Step 1: Optimize without leakage constraints
qcp_initial = SmoothPulseProblem(qtraj, N; Q=100.0)
solve!(qcp_initial; max_iter=100)

# Step 2: Add leakage suppression
opts = PiccoloOptions(leakage_cost=10.0, leakage_constraint=true, leakage_constraint_value=1e-3)
qcp_leakage = SmoothPulseProblem(qtraj, N; Q=100.0, piccolo_options=opts)
solve!(qcp_leakage; max_iter=150)
```

### 2. Increase Gate Time

Faster gates often have more leakage:

```julia
# If leakage is high, try longer gate time
T_long = 30.0
times_long = collect(range(0, T_long, length=N))
pulse_long = ZeroOrderPulse(0.05 * randn(2, N), times_long)
```

### 3. Use More Levels

Include more levels to capture dynamics accurately:

```julia
# 4 levels instead of 3
sys = TransmonSystem(levels=4, δ=0.2, drive_bounds=[0.2, 0.2])
```

### 4. Reduce Drive Amplitude

Lower drives can reduce leakage:

```julia
sys = TransmonSystem(levels=3, δ=0.2, drive_bounds=[0.1, 0.1])  # Reduced
```

## Analyzing Leakage

### Visualize Population Dynamics

After solving, plot the state populations:

```julia
using CairoMakie

traj = get_trajectory(qcp)

# Plot unitary elements including leakage levels
fig = plot(traj, [:Ũ⃗])
```

### Compute Final Leakage

```julia
traj = get_trajectory(qcp)
final_state = traj[:Ũ⃗][:, end]
U_final = iso_vec_to_operator(final_state)

# Compute population in leakage states
# (depends on your specific setup)
```

## See Also

- [Operators](@ref operators-concept) - EmbeddedOperator details
- [System Templates](@ref system-templates) - Creating multilevel systems
- [Objectives](@ref objectives-concept) - LeakageObjective documentation
- [Constraints](@ref constraints-concept) - LeakageConstraint documentation
