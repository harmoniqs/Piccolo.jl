# [SplinePulseProblem](@id spline-pulse)

`SplinePulseProblem` sets up trajectory optimization with spline-based pulses where derivative variables represent spline slopes or tangents. This is ideal for inherently smooth control pulses and warm-starting from previous solutions.

## When to Use

Use `SplinePulseProblem` when:
- You need inherently smooth control pulses
- You're warm-starting from a previously optimized solution
- You want cubic spline smoothness without derivative regularization
- You're working with hardware that expects smooth pulse shapes

## Pulse Requirements

`SplinePulseProblem` works with spline pulse types:

| Pulse Type | Derivative Meaning |
|------------|-------------------|
| `LinearSplinePulse` | `:du` represents constrained slope between knots |
| `CubicSplinePulse` | `:du` represents independent Hermite tangents |

```julia
# Linear spline
pulse = LinearSplinePulse(controls, times)
qtraj = UnitaryTrajectory(sys, pulse, U_goal)
qcp = SplinePulseProblem(qtraj)  # ✓ Works

# Cubic spline
pulse = CubicSplinePulse(controls, tangents, times)
qtraj = UnitaryTrajectory(sys, pulse, U_goal)
qcp = SplinePulseProblem(qtraj)  # ✓ Works
```

## Constructor Variants

### Use Native Knot Times (Recommended for Warm-Starting)

```julia
SplinePulseProblem(qtraj::AbstractQuantumTrajectory{<:AbstractSplinePulse}; kwargs...)
```

Uses the pulse's native knot times without resampling. Best for warm-starting from a previous solution.

### Resample to N Timesteps

```julia
SplinePulseProblem(qtraj::AbstractQuantumTrajectory{<:AbstractSplinePulse}, N::Int; kwargs...)
```

Resamples the pulse to `N` uniformly spaced timesteps.

### Use Specific Times

```julia
SplinePulseProblem(qtraj::AbstractQuantumTrajectory{<:AbstractSplinePulse}, times::AbstractVector; kwargs...)
```

Resamples the pulse to the specified time points.

## Parameter Reference

### Objective Weights

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `Q` | `Float64` | `100.0` | Weight on infidelity objective |
| `R` | `Float64` | `1e-2` | Base regularization weight |
| `R_u` | `Float64` or `Vector{Float64}` | `R` | Regularization on control values |
| `R_du` | `Float64` or `Vector{Float64}` | `R` | Regularization on derivatives/tangents |

### Bounds

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `du_bound` | `Float64` | `Inf` | Maximum derivative/slope bound |
| `Δt_bounds` | `Tuple{Float64, Float64}` | `nothing` | Time-step bounds for free-time optimization |

### Advanced Options

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `integrator` | `AbstractIntegrator` | `nothing` | Custom integrator (uses `BilinearIntegrator` if `nothing`) |
| `global_names` | `Vector{Symbol}` | `nothing` | Global parameters to optimize |
| `global_bounds` | `Dict{Symbol, ...}` | `nothing` | Bounds on global variables |
| `constraints` | `Vector{AbstractConstraint}` | `[]` | Additional constraints |
| `piccolo_options` | `PiccoloOptions` | `PiccoloOptions()` | Solver options |

## Examples

### Basic Spline Optimization

```julia
using Piccolo

# Define system
H_drift = PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
sys = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

# Create cubic spline pulse
T, N = 10.0, 50
times = collect(range(0, T, length=N))
controls = 0.1 * randn(2, N)
tangents = zeros(2, N)  # Initial tangents
pulse = CubicSplinePulse(controls, tangents, times)

qtraj = UnitaryTrajectory(sys, pulse, GATES[:X])

# Solve using native knot times
qcp = SplinePulseProblem(qtraj; Q=100.0, du_bound=10.0)
solve!(qcp; max_iter=100)
```

### Warm-Starting from Previous Solution

```julia
# Load previously optimized pulse
using JLD2
@load "optimized_pulse.jld2" saved_pulse

# Create new trajectory with saved pulse
qtraj = UnitaryTrajectory(sys, saved_pulse, U_goal)

# Use native knot times (no resampling)
qcp = SplinePulseProblem(qtraj)
solve!(qcp; max_iter=50)  # Converges quickly from good initial guess
```

### Resampling to Different Resolution

```julia
# Original pulse has 50 knots
pulse = CubicSplinePulse(controls_50, tangents_50, times_50)
qtraj = UnitaryTrajectory(sys, pulse, U_goal)

# Resample to 100 knots for finer control
qcp = SplinePulseProblem(qtraj, 100; Q=100.0)
solve!(qcp; max_iter=100)
```

### Linear vs Cubic Splines

**Linear Splines**: The derivative `:du` represents the slope between knots. A `DerivativeIntegrator` constraint enforces `du[k] = (u[k+1] - u[k]) / Δt`.

```julia
pulse = LinearSplinePulse(controls, times)
qtraj = UnitaryTrajectory(sys, pulse, U_goal)
qcp = SplinePulseProblem(qtraj)
# du is constrained to be consistent with u
```

**Cubic Splines (Hermite)**: The derivative `:du` represents independent Hermite tangents at each knot. These are free optimization variables with no inter-knot constraint.

```julia
pulse = CubicSplinePulse(controls, tangents, times)
qtraj = UnitaryTrajectory(sys, pulse, U_goal)
qcp = SplinePulseProblem(qtraj)
# du (tangents) are independent variables
```

## Trajectory Structure

Unlike `SmoothPulseProblem` which has three derivative levels (`:u`, `:du`, `:ddu`), `SplinePulseProblem` only has one:

| Variable | Description |
|----------|-------------|
| `:u` | Control values at knot points |
| `:du` | Derivatives/tangents (meaning depends on spline type) |

The reduced number of variables can lead to faster optimization while maintaining smooth pulses through the spline interpolation.

## See Also

- [SmoothPulseProblem](@ref smooth-pulse) - For piecewise constant controls
- [MinimumTimeProblem](@ref minimum-time) - Time-optimal control
- [Pulses](@ref) - Detailed pulse type documentation
