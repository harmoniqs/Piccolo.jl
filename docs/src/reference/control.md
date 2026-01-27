# [Control Module](@id control-api)

The Control module provides the optimal control framework: problems, objectives, constraints, integrators, and problem templates.

## Quantum Control Problems

### QuantumControlProblem

The main wrapper type that combines a quantum trajectory with an optimization problem.

```julia
qcp = SmoothPulseProblem(qtraj, N; kwargs...)  # Creates a QuantumControlProblem
```

**Key functions:**
- `get_trajectory(qcp)` - Get the NamedTrajectory
- `get_system(qcp)` - Get the QuantumSystem
- `get_goal(qcp)` - Get the goal state/operator
- `sync_trajectory!(qcp)` - Update quantum trajectory from optimized values
- `fidelity(qcp)` - Compute fidelity of the solution

## Problem Templates

### SmoothPulseProblem

For piecewise constant controls (`ZeroOrderPulse`). Adds discrete derivative variables for smoothness.

```julia
qcp = SmoothPulseProblem(qtraj, N; Q=100.0, R=1e-2, ddu_bound=1.0)
```

### SplinePulseProblem

For spline-based controls (`LinearSplinePulse`, `CubicSplinePulse`).

```julia
qcp = SplinePulseProblem(qtraj, N; Q=100.0, R=1e-2)
```

### MinimumTimeProblem

Wraps an existing problem to minimize gate duration while maintaining fidelity.

```julia
qcp_mintime = MinimumTimeProblem(qcp; final_fidelity=0.99, D=100.0)
```

### SamplingProblem

For robust optimization over multiple system variants.

```julia
qcp_robust = SamplingProblem(qcp, [sys1, sys2, sys3]; Q=100.0)
```

## Options

### PiccoloOptions

Configuration for solver behavior and additional features like leakage suppression.

```julia
opts = PiccoloOptions(
    leakage_constraint = true,
    leakage_constraint_value = 1e-2,
    leakage_cost = 1e-2,
    verbose = true
)
```

## Objectives

The following objective types are available:

- `UnitaryInfidelityObjective` - Unitary gate fidelity
- `KetInfidelityObjective` - State transfer fidelity
- `CoherentKetInfidelityObjective` - Phase-coherent multi-state fidelity
- `UnitaryFreePhaseInfidelityObjective` - Gate fidelity ignoring global phase
- `LeakageObjective` - Penalize population in leakage levels
- `UnitarySensitivityObjective` - Sensitivity analysis

Regularization is provided via `QuadraticRegularizer` from DirectTrajOpt.jl.

## Constraints

The following constraint types are available:

- `FinalUnitaryFidelityConstraint` - Ensure minimum gate fidelity
- `FinalKetFidelityConstraint` - Ensure minimum state fidelity
- `FinalCoherentKetFidelityConstraint` - Phase-coherent fidelity constraint
- `LeakageConstraint` - Limit leakage to higher levels

## Integrators

Dynamics integrators (`BilinearIntegrator`, `DerivativeIntegrator`) are provided by DirectTrajOpt.jl and used internally by the problem templates.

## Solving

The `solve!` function is provided by DirectTrajOpt.jl and reexported:

```julia
solve!(qcp::QuantumControlProblem; max_iter=100, kwargs...)
```

Solves the optimization problem and updates the trajectory. Automatically calls `sync_trajectory!` to update the quantum trajectory with the optimized solution.

### Common Keyword Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `max_iter` | `Int` | Maximum iterations |
| `tol` | `Float64` | Convergence tolerance |
| `verbose` | `Bool` | Print solver output |
| `print_level` | `Int` | Ipopt print level (0-12) |
