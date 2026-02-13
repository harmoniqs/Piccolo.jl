# [Problem Templates](@id problem-templates-overview)

Problem templates are the primary user-facing API in Piccolo.jl. They provide high-level constructors that set up quantum optimal control problems with sensible defaults while allowing fine-grained customization.

## Overview

Piccolo.jl provides four main problem templates:

| Template | Purpose | Pulse Type |
|----------|---------|------------|
| [`SmoothPulseProblem`](@ref smooth-pulse) | Piecewise constant controls with smoothness regularization | `ZeroOrderPulse` |
| [`SplinePulseProblem`](@ref spline-pulse) | Spline-based controls for smooth pulses | `LinearSplinePulse`, `CubicSplinePulse` |
| [`MinimumTimeProblem`](@ref minimum-time) | Time-optimal control (wraps another problem) | Any |
| [`SamplingProblem`](@ref sampling) | Robust optimization over parameter variations | Any |

## Choosing a Template

Use this decision flowchart to select the right template for your task:

```
        ┌──────────────────────┐
        │  What pulse type?    │
        └──────────┬───────────┘
           ┌───────┴───────┐
           ▼               ▼
    ZeroOrderPulse    SplinePulse
           │               │
           ▼               ▼
  ┌─────────────────┐ ┌─────────────────┐
  │ SmoothPulse     │ │ SplinePulse     │
  │ Problem         │ │ Problem         │
  └────────┬────────┘ └────────┬────────┘
           └─────────┬─────────┘
                     ▼
           ┌──────────────────┐
           │ Need robustness? │
           └───┬──────────┬───┘
            Yes│          │No
               ▼          │
       ┌──────────────┐   │
       │ Sampling     │   │
       │ Problem      │   │
       └──────┬───────┘   │
              └─────┬─────┘
                    ▼
           ┌────────────────┐
           │ Minimize time? │
           └───┬────────┬───┘
            Yes│        │No
               ▼        │
       ┌──────────────┐ │
       │ MinimumTime  │ │
       │ Problem      │ │
       └──────┬───────┘ │
              └────┬────┘
                   ▼
              ┌────────┐
              │ Done!  │
              └────────┘
```

### Quick Decision Guide

**Choose your base problem:**
- **`SmoothPulseProblem`**: Start here for most problems. Uses piecewise constant pulses with derivative regularization for smoothness.
- **`SplinePulseProblem`**: Use when you need inherently smooth pulses or want to warm-start from a previous solution.

**Add robustness (optional):**
- **`SamplingProblem`**: Wrap your base problem to optimize over multiple system variants (e.g., different detuning values).

**Minimize time (optional):**
- **`MinimumTimeProblem`**: Wrap any problem to find the shortest gate duration that achieves a target fidelity.

## Composability

Problem templates can be chained together:

```julia
# Step 1: Create base problem with free time enabled
qcp_base = SmoothPulseProblem(qtraj, N; Δt_bounds=(0.01, 0.5))
solve!(qcp_base; max_iter=100)

# Step 2: Add robustness
qcp_robust = SamplingProblem(qcp_base, [sys_nominal, sys_perturbed])
solve!(qcp_robust; max_iter=100)

# Step 3: Minimize time while maintaining fidelity
qcp_mintime = MinimumTimeProblem(qcp_robust; final_fidelity=0.99)
solve!(qcp_mintime; max_iter=100)
```

## Common Workflow

All problem templates follow the same workflow:

```julia
# 1. Define quantum system
sys = QuantumSystem(H_drift, H_drives, drive_bounds)

# 2. Create initial pulse
pulse = ZeroOrderPulse(initial_controls, times)

# 3. Create trajectory with goal
qtraj = UnitaryTrajectory(sys, pulse, U_goal)

# 4. Set up optimization problem
qcp = SmoothPulseProblem(qtraj, N; Q=100.0, R=1e-2)

# 5. Solve
solve!(qcp; max_iter=100)

# 6. Analyze results
println("Fidelity: ", fidelity(qcp))
optimized_pulse = get_pulse(qcp.qtraj)
```

## Key Concepts

### Pulse Type Matching

Each base problem template requires a specific pulse type:

| Problem Template | Required Pulse Type |
|-----------------|---------------------|
| `SmoothPulseProblem` | `ZeroOrderPulse` |
| `SplinePulseProblem` | `LinearSplinePulse` or `CubicSplinePulse` |

Using the wrong pulse type will result in a helpful error message.

### Trajectory Types

Problem templates work with different trajectory types depending on your goal:

| Trajectory Type | Use Case |
|-----------------|----------|
| `UnitaryTrajectory` | Gate synthesis (most common) |
| `KetTrajectory` | State transfer |
| `DensityTrajectory` | Open system evolution |
| `MultiKetTrajectory` | Multiple state transfers with coherent phases |

### Regularization Parameters

Most problems use these regularization parameters:

- **`Q`**: Weight on infidelity objective (higher = prioritize fidelity)
- **`R`**: Base regularization weight for control smoothness
- **`R_u`**, **`R_du`**, **`R_ddu`**: Per-derivative regularization weights

### Free-Time Optimization

To enable time optimization with `MinimumTimeProblem`, set `Δt_bounds` in your base problem:

```julia
qcp = SmoothPulseProblem(qtraj, N; Δt_bounds=(0.01, 0.5))
```

## Detailed Documentation

- [SmoothPulseProblem](@ref smooth-pulse) - Full parameter reference and examples
- [SplinePulseProblem](@ref spline-pulse) - Spline-based optimization
- [MinimumTimeProblem](@ref minimum-time) - Time-optimal control
- [SamplingProblem](@ref sampling) - Robust optimization
- [Composing Templates](@ref composition) - Advanced composition patterns
