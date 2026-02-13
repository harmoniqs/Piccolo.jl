# [Pulses](@id pulses-concept)

Pulses in Piccolo.jl parameterize how control signals vary over time. The choice of pulse type affects both the optimization problem structure and the resulting control smoothness.

## Overview

| Pulse Type | Description | Use With |
|------------|-------------|----------|
| `ZeroOrderPulse` | Piecewise constant | `SmoothPulseProblem` |
| `LinearSplinePulse` | Linear interpolation | `SplinePulseProblem` |
| `CubicSplinePulse` | Cubic Hermite splines | `SplinePulseProblem` |
| `GaussianPulse` | Parametric Gaussian | Analytical evaluation |
| `CompositePulse` | Combination of pulses | Various |

## ZeroOrderPulse

Piecewise constant (zero-order hold) controls. The most common choice for initial optimization.

### Construction

```julia
using Piccolo

n_drives = 2
N = 100
T = 10.0

# Time points (N points define N-1 intervals)
times = collect(range(0, T, length=N))

# Control values: n_drives × N matrix
controls = 0.1 * randn(n_drives, N)

pulse = ZeroOrderPulse(controls, times)
```

### Properties

```julia
# Evaluate at time t
u = pulse(t)  # Returns control vector at time t

# Duration
T = duration(pulse)

# Number of drives
n = n_drives(pulse)

# Sample at new times
new_pulse = sample(pulse, new_times)
```

### Visualization

```
Control Value
    │    ┌──┐
    │    │  │  ┌──┐
    │────┘  │  │  └───
    │       └──┘
    └─────────────────── Time
```

### Use Case

- **Primary use**: `SmoothPulseProblem`
- **Characteristics**: Simple structure, smoothness via derivative regularization
- **Best for**: Initial optimization, most quantum control problems

## LinearSplinePulse

Linear interpolation between control knots.

### Construction

```julia
times = collect(range(0, T, length=N))
controls = 0.1 * randn(n_drives, N)

pulse = LinearSplinePulse(controls, times)
```

### Properties

- Continuous control values
- Discontinuous first derivative (at knots)
- Derivative = slope between knots

### Visualization

```
Control Value
    │      /\
    │     /  \    /
    │    /    \  /
    │───/      \/
    └─────────────────── Time
```

### Use Case

- **Primary use**: `SplinePulseProblem`
- **Characteristics**: Continuous but not smooth
- **Best for**: Simple continuous controls

## CubicSplinePulse

Cubic Hermite spline interpolation with independent tangents at each knot.

### Construction

```julia
times = collect(range(0, T, length=N))
controls = 0.1 * randn(n_drives, N)
tangents = zeros(n_drives, N)  # Initial tangents (slopes)

pulse = CubicSplinePulse(controls, tangents, times)
```

### Properties

- Continuous control values AND first derivatives
- Tangents are independent optimization variables
- Smooth C¹ continuous curves

### Visualization

```
Control Value
    │      ╭─╮
    │     ╱   ╲    ╱
    │    ╱     ╲  ╱
    │───╱       ╲╱
    └─────────────────── Time
```

### Use Case

- **Primary use**: `SplinePulseProblem`
- **Characteristics**: Inherently smooth, good for warm-starting
- **Best for**: Hardware-friendly pulses, refining solutions

## GaussianPulse

Parametric Gaussian envelope with analytical form.

### Construction

```julia
# Parameters per drive
amplitudes = [0.5, 0.3]
centers = [5.0, 5.0]
widths = [1.0, 1.5]

pulse = GaussianPulse(amplitudes, centers, widths)
```

### Mathematical Form

```math
u_i(t) = A_i \exp\left(-\frac{(t - \mu_i)^2}{2\sigma_i^2}\right)
```

### Use Case

- **Analytical evaluation**: No discretization required
- **Initialization**: Generate initial guesses
- **Best for**: Simple pulse shapes, understanding pulse structure

## CompositePulse

Combine multiple pulses (sum or sequence).

### Construction

```julia
pulse1 = GaussianPulse([0.5], [2.0], [0.5])
pulse2 = GaussianPulse([0.3], [8.0], [0.5])

# Sum the pulses
composite = CompositePulse([pulse1, pulse2])
```

### Use Case

- **Complex shapes**: Build from simpler components
- **Multi-stage pulses**: Different behaviors in different time windows

## Pulse Operations

### Evaluation

```julia
# Single time point
u = pulse(t)  # Vector of length n_drives

# Multiple time points
us = [pulse(t) for t in times]
```

### Sampling

Convert between pulse types or resample:

```julia
# Resample to new times
new_times = collect(range(0, T, length=200))
resampled = sample(pulse, new_times)
```

### Duration

```julia
T = duration(pulse)
```

### Drive Names

```julia
names = drive_name(pulse)  # Returns symbol or vector of symbols
```

## Choosing a Pulse Type

### Decision Guide

```
Start
  │
  ▼
Is this your first optimization attempt?
  │
  ├── Yes → ZeroOrderPulse + SmoothPulseProblem
  │
  └── No → Do you have a previous solution?
              │
              ├── Yes → CubicSplinePulse + SplinePulseProblem (warm-start)
              │
              └── No → Do you need smooth pulses?
                          │
                          ├── Yes → CubicSplinePulse
                          │
                          └── No → ZeroOrderPulse
```

### Practical Recommendations

| Scenario | Recommended Pulse |
|----------|-------------------|
| Starting fresh | `ZeroOrderPulse` |
| Refining a solution | `CubicSplinePulse` |
| Hardware requires smooth pulses | `CubicSplinePulse` |
| Simple continuous pulses | `LinearSplinePulse` |
| Analytical pulse design | `GaussianPulse` |

## Converting Between Pulse Types

### ZeroOrderPulse → CubicSplinePulse

```julia
# After optimization with SmoothPulseProblem
optimized_zop = get_pulse(qcp.qtraj)

# Sample to get control values
controls = hcat([optimized_zop(t) for t in times]...)

# Estimate tangents (finite differences)
tangents = similar(controls)
for k in 1:N-1
    tangents[:, k] = (controls[:, k+1] - controls[:, k]) / (times[k+1] - times[k])
end
tangents[:, N] = tangents[:, N-1]

# Create cubic spline
cubic_pulse = CubicSplinePulse(controls, tangents, times)
```

### Arbitrary → ZeroOrderPulse

```julia
# Sample any pulse type
times = collect(range(0, duration(pulse), length=N))
controls = hcat([pulse(t) for t in times]...)

zop = ZeroOrderPulse(controls, times)
```

## Best Practices

### 1. Initialize Appropriately

```julia
# Scale by drive bounds
max_amp = 0.1 * maximum(drive_bounds)
controls = max_amp * randn(n_drives, N)
```

### 2. Use Enough Time Points

```julia
# Rule of thumb: ~10 points per characteristic time scale
T = 10.0  # Total time
τ = 1.0   # Shortest feature you want to capture
N = ceil(Int, 10 * T / τ)
```

### 3. Start with ZeroOrderPulse

Even if you need smooth pulses, optimize with `ZeroOrderPulse` first, then convert:

```julia
# Step 1: Optimize with ZeroOrderPulse
pulse_init = ZeroOrderPulse(randn(n_drives, N), times)
qtraj = UnitaryTrajectory(sys, pulse_init, U_goal)
qcp = SmoothPulseProblem(qtraj, N)
solve!(qcp)

# Step 2: Convert and refine with CubicSplinePulse
cubic_pulse = CubicSplinePulse(get_pulse(qcp.qtraj))
qtraj_refined = UnitaryTrajectory(sys, cubic_pulse, U_goal)
qcp_refined = SplinePulseProblem(qtraj_refined)
solve!(qcp_refined; max_iter=50)
```

## See Also

- [SmoothPulseProblem](@ref smooth-pulse) - Using `ZeroOrderPulse`
- [SplinePulseProblem](@ref spline-pulse) - Using spline pulses
- [Trajectories](@ref trajectories-concept) - Combining pulses with systems
