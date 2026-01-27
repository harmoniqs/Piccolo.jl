# [Custom Objectives](@id custom-objectives)

While Piccolo.jl provides standard objectives for fidelity and regularization, you may need custom objectives for specialized optimization goals.

## Objective Types

Piccolo.jl uses objectives from DirectTrajOpt.jl:

| Type | When Evaluated |
|------|----------------|
| `TerminalObjective` | At final timestep only |
| `KnotPointObjective` | At every timestep |

## Creating Custom Objectives

### Terminal Objective Example

A terminal objective adds cost based on the final state only.

```julia
using DirectTrajOpt

# Custom terminal objective: penalize distance from target state
struct MyTerminalObjective <: TerminalObjective
    target::Vector{Float64}
    weight::Float64
end

function DirectTrajOpt.evaluate(obj::MyTerminalObjective, traj, k)
    state = traj[:ψ̃][:, k]
    error = state - obj.target
    return obj.weight * dot(error, error)
end
```

### Knotpoint Objective Example

A knotpoint objective adds cost at every timestep.

```julia
# Custom objective: penalize control energy over time
struct ControlEnergyObjective <: KnotPointObjective
    name::Symbol
    weight::Float64
end

function DirectTrajOpt.evaluate(obj::ControlEnergyObjective, traj, k)
    u = traj[obj.name][:, k]
    return obj.weight * dot(u, u)
end
```

## Adding Custom Objectives

### Via constraints Parameter

```julia
my_obj = MyTerminalObjective(target_state, 10.0)

qcp = SmoothPulseProblem(
    qtraj, N;
    constraints = [my_obj]  # Can include objectives too
)
```

### Direct Problem Modification

```julia
# Create problem first
qcp = SmoothPulseProblem(qtraj, N)

# Add custom objective
my_obj = MyTerminalObjective(target_state, 10.0)
add_objective!(qcp.prob, my_obj)

# Then solve
solve!(qcp)
```

## Common Custom Objectives

### Spectral Penalty

Penalize high-frequency content in controls:

```julia
struct SpectralPenalty <: KnotPointObjective
    control_name::Symbol
    weight::Float64
end

function DirectTrajOpt.evaluate(obj::SpectralPenalty, traj, k)
    if k == 1
        return 0.0
    end
    u_k = traj[obj.control_name][:, k]
    u_km1 = traj[obj.control_name][:, k-1]
    diff = u_k - u_km1
    return obj.weight * dot(diff, diff)
end
```

### State Constraint Penalty

Soft constraint on state staying in a region:

```julia
struct StateRegionPenalty <: KnotPointObjective
    state_name::Symbol
    threshold::Float64
    weight::Float64
end

function DirectTrajOpt.evaluate(obj::StateRegionPenalty, traj, k)
    state = traj[obj.state_name][:, k]
    violation = max(0, norm(state) - obj.threshold)
    return obj.weight * violation^2
end
```

### Time-Varying Weight

Apply different weights at different times:

```julia
struct TimeVaryingRegularizer <: KnotPointObjective
    control_name::Symbol
    weights::Vector{Float64}  # One weight per timestep
end

function DirectTrajOpt.evaluate(obj::TimeVaryingRegularizer, traj, k)
    u = traj[obj.control_name][:, k]
    return obj.weights[k] * dot(u, u)
end
```

## Gradient Computation

DirectTrajOpt uses automatic differentiation (ForwardDiff), so gradients are computed automatically for most objectives. However, you can provide analytical gradients for efficiency:

```julia
function DirectTrajOpt.gradient!(grad, obj::MyObjective, traj, k)
    # Fill grad with ∂J/∂x for each variable x
    # This is optional - ForwardDiff is used if not provided
end
```

## Combining Multiple Objectives

```julia
# Create multiple custom objectives
obj1 = MyTerminalObjective(target, 10.0)
obj2 = SpectralPenalty(:u, 1.0)
obj3 = StateRegionPenalty(:ψ̃, 2.0, 5.0)

# Add all to problem
qcp = SmoothPulseProblem(qtraj, N; constraints=[obj1, obj2, obj3])
```

## Tips for Custom Objectives

### 1. Scale Appropriately

Match the scale of built-in objectives:

```julia
# Built-in fidelity uses Q ~ 100
# Custom objectives should use similar scale
my_obj = MyObjective(weight=50.0)  # Same order of magnitude
```

### 2. Ensure Smoothness

Avoid discontinuities that can cause optimization issues:

```julia
# Bad: discontinuous
penalty = x > 0 ? x^2 : 0

# Good: smooth approximation
penalty = max(0, x)^2
# Or use softplus: log(1 + exp(k*x)) / k
```

### 3. Test Independently

Verify your objective computes expected values:

```julia
# Create test trajectory
test_traj = ...

# Evaluate objective
value = evaluate(my_obj, test_traj, N)
println("Objective value: ", value)

# Check gradient (if provided)
grad = zeros(length_of_variables)
gradient!(grad, my_obj, test_traj, N)
```

### 4. Start Simple

Add custom objectives incrementally:

```julia
# Step 1: Solve with standard objectives
qcp = SmoothPulseProblem(qtraj, N)
solve!(qcp)

# Step 2: Check if solution needs improvement
# ... analyze results ...

# Step 3: Add custom objective with small weight
qcp = SmoothPulseProblem(qtraj, N; constraints=[MyObjective(1.0)])
solve!(qcp)

# Step 4: Increase weight if needed
```

## See Also

- [Objectives](@ref) - Built-in objectives
- [Constraints](@ref) - Hard constraints (vs soft objective penalties)
- [Problem Templates](@ref) - Using objectives in problems
