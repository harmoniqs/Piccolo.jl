# [Constraints](@id constraints-concept)

Constraints define hard requirements that optimized solutions must satisfy. Unlike objectives (which are minimized), constraints are enforced exactly.

## Overview

Piccolo.jl supports several constraint types:

| Constraint Type | Description |
|-----------------|-------------|
| Bound constraints | Limits on variable values |
| Fidelity constraints | Minimum fidelity requirements |
| Leakage constraints | Maximum allowed leakage |
| Equality constraints | Custom equality requirements |

## Bound Constraints

### Control Bounds

Control bounds are specified in the `QuantumSystem` and automatically enforced:

```julia
# Bounds specified at system creation
drive_bounds = [1.0, 0.5]  # Drive 1: ±1.0, Drive 2: ±0.5
sys = QuantumSystem(H_drift, H_drives, drive_bounds)
```

### Derivative Bounds

Derivative bounds limit how fast controls can change:

```julia
qcp = SmoothPulseProblem(
    qtraj, N;
    du_bound=0.5,    # Max control jump per timestep
    ddu_bound=0.1    # Max control acceleration
)
```

### Timestep Bounds

For free-time optimization:

```julia
qcp = SmoothPulseProblem(
    qtraj, N;
    Δt_bounds=(0.01, 0.5)  # Min and max timestep
)
```

## Fidelity Constraints

Used with `MinimumTimeProblem` to enforce minimum gate quality.

### FinalUnitaryFidelityConstraint

```julia
constraint = FinalUnitaryFidelityConstraint(
    :Ũ⃗,        # State variable name
    U_goal,    # Target unitary
    0.99       # Minimum fidelity
)
```

### FinalKetFidelityConstraint

```julia
constraint = FinalKetFidelityConstraint(
    :ψ̃,        # State variable name
    ψ_goal,    # Target state
    0.99       # Minimum fidelity
)
```

### FinalCoherentKetFidelityConstraint

For `MultiKetTrajectory`:

```julia
constraint = FinalCoherentKetFidelityConstraint(
    [:ψ̃1, :ψ̃2],      # State variable names
    [ψ_goal1, ψ_goal2], # Target states
    0.99                 # Minimum fidelity
)
```

### Automatic Setup in MinimumTimeProblem

You typically don't create fidelity constraints manually. `MinimumTimeProblem` adds them automatically:

```julia
# Automatically adds FinalUnitaryFidelityConstraint for UnitaryTrajectory
qcp_mintime = MinimumTimeProblem(qcp_base; final_fidelity=0.99)
```

## Leakage Constraints

### LeakageConstraint

Bounds population outside the computational subspace:

```julia
constraint = LeakageConstraint(
    :Ũ⃗,        # State variable name
    op,        # EmbeddedOperator defining subspace
    1e-3       # Maximum allowed leakage
)
```

### Via PiccoloOptions

The easier approach:

```julia
opts = PiccoloOptions(
    leakage_constraint=true,
    leakage_constraint_value=1e-3
)

qcp = SmoothPulseProblem(qtraj, N; piccolo_options=opts)
```

## Adding Custom Constraints

### To Problem Templates

Pass constraints via the `constraints` keyword:

```julia
my_constraint = MyCustomConstraint(...)

qcp = SmoothPulseProblem(
    qtraj, N;
    constraints=[my_constraint]
)
```

### Multiple Constraints

```julia
constraints = [
    FinalUnitaryFidelityConstraint(:Ũ⃗, U_goal, 0.99),
    LeakageConstraint(:Ũ⃗, op, 1e-3)
]

qcp = SmoothPulseProblem(qtraj, N; constraints=constraints)
```

## PiccoloOptions

`PiccoloOptions` provides a convenient way to configure common constraint settings:

```julia
opts = PiccoloOptions(
    # Leakage handling
    leakage_constraint=true,
    leakage_constraint_value=1e-3,
    leakage_cost=10.0,  # Also adds objective

    # Timestep handling
    timesteps_all_equal=false,

    # Complex control constraints
    complex_control_norm_constraint_name=:u_norm,
    complex_control_norm_constraint_radius=1.0,

    # Verbosity
    verbose=true
)

qcp = SmoothPulseProblem(qtraj, N; piccolo_options=opts)
```

### Key Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `leakage_constraint` | `Bool` | `false` | Enable leakage constraint |
| `leakage_constraint_value` | `Float64` | `1e-3` | Maximum leakage |
| `leakage_cost` | `Float64` | `0.0` | Leakage objective weight |
| `timesteps_all_equal` | `Bool` | `false` | Force uniform timesteps |
| `verbose` | `Bool` | `false` | Print solver progress |

## Constraints vs Objectives

Understanding when to use each:

| Use Case | Constraint | Objective |
|----------|------------|-----------|
| Must achieve F ≥ 0.99 | ✓ Fidelity constraint | |
| Prefer higher fidelity | | ✓ Infidelity objective |
| Control must be ≤ 1.0 | ✓ Bound constraint | |
| Prefer smaller controls | | ✓ Regularization |
| Leakage must be < 1e-3 | ✓ Leakage constraint | |
| Prefer less leakage | | ✓ Leakage objective |

### Trade-offs

**Constraints:**
- Guarantee satisfaction (if feasible)
- Can make problem harder to solve
- May be infeasible

**Objectives:**
- More flexible
- Easier optimization
- No guarantees

### Recommendation

Start with objectives, add constraints for hard requirements:

```julia
# Step 1: Optimize with objectives only
qcp = SmoothPulseProblem(qtraj, N; Q=100.0)
solve!(qcp)

# Step 2: If fidelity constraint needed, use MinimumTimeProblem
if need_guaranteed_fidelity
    qcp_constrained = MinimumTimeProblem(qcp; final_fidelity=0.99)
    solve!(qcp_constrained)
end
```

## Constraint Feasibility

### Checking Feasibility

If optimization fails, constraints may be infeasible:

```julia
# Check if constraints are satisfied
traj = get_trajectory(qcp)
# Evaluate constraint violation manually if needed
```

### Common Causes of Infeasibility

1. **Fidelity too high**: Target fidelity may not be achievable
2. **Time too short**: Insufficient time for the gate
3. **Bounds too tight**: Controls can't reach required values
4. **Conflicting constraints**: e.g., low leakage with high fidelity in short time

### Solutions

1. **Relax constraints**: Lower fidelity target, increase time
2. **Better initialization**: Start from a working solution
3. **Adjust bounds**: Allow larger controls or longer time
4. **Use objectives first**: Find a good solution, then add constraints

## Best Practices

### 1. Start Without Constraints

```julia
# First, find a good solution
qcp = SmoothPulseProblem(qtraj, N; Q=100.0)
solve!(qcp; max_iter=200)
println("Achieved fidelity: ", fidelity(qcp))
```

### 2. Add Constraints Gradually

```julia
# Then add constraints if needed
if fidelity(qcp) > 0.99
    qcp_constrained = MinimumTimeProblem(qcp; final_fidelity=0.99)
    solve!(qcp_constrained)
end
```

### 3. Use Margin for Robustness

```julia
# Target slightly higher than needed
qcp = MinimumTimeProblem(qcp_base; final_fidelity=0.995)  # Want 0.99
```

### 4. Monitor Constraint Satisfaction

```julia
solve!(qcp; max_iter=100)
println("Final fidelity: ", fidelity(qcp))
# Check other constraints as needed
```

## See Also

- [Objectives](@ref) - Soft optimization targets
- [Problem Templates](@ref) - Using constraints in practice
- [MinimumTimeProblem](@ref minimum-time) - Fidelity-constrained time optimization
