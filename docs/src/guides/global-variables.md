# [Global Variables](@id global-variables)

Global variables allow you to optimize system parameters alongside control pulses. This is useful for calibrating system properties or finding optimal operating points.

## When to Use Global Variables

- **System calibration**: Finding optimal qubit frequencies or coupling strengths
- **Gate design**: Optimizing hardware parameters for better gates
- **Sensitivity analysis**: Understanding parameter effects

## Basic Concept

Standard optimization only varies control pulses `u(t)`. With global variables, you can also optimize system parameters that appear in the Hamiltonian.

```math
H(u, t; \theta) = H_{\text{drift}}(\theta) + \sum_i u_i(t) H_{\text{drive},i}(\theta)
```

Where `θ` are global parameters to optimize.

## Setup Requirements

Global variable optimization requires:
1. A **time-dependent system** with parameters
2. A **custom integrator** that supports globals (e.g., from Piccolissimo)
3. **Global bounds** specification

## Example: Optimizing Detuning

### Step 1: Define Parametric System

```julia
using Piccolo

# Hamiltonian with global parameter δ (detuning)
function H(u, t)
    δ = u[end]
    return δ * PAULIS[:Z] + u[1] * PAULIS[:X] + u[2] * PAULIS[:Y]
end

# Create system with named global parameter
sys = QuantumSystem(H, [1.0, 1.0]; time_dependent=true, global_params=(δ=0.1,))
```

### Step 2: Use Custom Integrator

Global variable optimization requires an integrator that handles globals:

```julia
# Example with a compatible integrator
# (specific integrator depends on your setup)
using Piccolissimo

# Create trajectory
T, N = 10.0, 100
times = collect(range(0, T, length=N))
pulse = ZeroOrderPulse(0.1 * randn(2, N), times)
qtraj = UnitaryTrajectory(sys, pulse, GATES[:X])

# Create integrator with global support
integrator = HermitianExponentialIntegrator(qtraj, N; global_names=[:δ])
```

### Step 3: Set Up Problem with Global Bounds

```julia
qcp = SmoothPulseProblem(
    qtraj, N;
    integrator = integrator,
    global_names = [:δ],
    global_bounds = Dict(:δ => (0.05, 0.2))  # δ ∈ [0.05, 0.2]
)

solve!(qcp; max_iter=100)

# Access optimized global value
traj = get_trajectory(qcp)
optimized_δ = traj[:δ][1]  # Global variables are constant
println("Optimized detuning: ", optimized_δ)
```

## Global Bounds Format

Global bounds can be specified in several ways:

```julia
# Symmetric bounds: δ ∈ [0.1 - 0.05, 0.1 + 0.05] = [0.05, 0.15]
global_bounds = Dict(:δ => 0.05)

# Asymmetric bounds: δ ∈ [0.05, 0.2]
global_bounds = Dict(:δ => (0.05, 0.2))

# Multiple globals
global_bounds = Dict(
    :δ => (0.05, 0.2),
    :J => 0.01  # J ∈ [J₀ - 0.01, J₀ + 0.01]
)
```

## Multiple Global Variables

```julia
# Hamiltonian with multiple parameters
function H(u, t; ω=1.0, δ=0.1, J=0.05)
    return ω * PAULIS[:Z] + δ * PAULIS[:Z] ⊗ PAULIS[:Z] +
           J * (PAULIS[:X] ⊗ PAULIS[:X]) +
           u[1] * PAULIS[:X] ⊗ I(2) + u[2] * I(2) ⊗ PAULIS[:X]
end

sys = QuantumSystem(H, bounds; global_params=(ω=1.0, δ=0.1, J=0.05))

# Optimize ω and J (keep δ fixed)
qcp = SmoothPulseProblem(
    qtraj, N;
    integrator = integrator,
    global_names = [:ω, :J],
    global_bounds = Dict(:ω => (0.9, 1.1), :J => 0.02)
)
```

## Best Practices

### 1. Start with Good Initial Values

```julia
# Use physically reasonable initial parameters
sys = QuantumSystem(H, bounds; global_params=(δ=0.15,))  # Close to expected
```

### 2. Use Reasonable Bounds

```julia
# Don't make bounds too wide
global_bounds = Dict(:δ => (0.1, 0.2))  # ~50% variation, not 10x
```

### 3. Consider Sensitivity

Global parameters affect all timesteps, so they have strong influence on the solution. Start with tighter bounds and relax if needed.

### 4. Verify Physical Meaning

After optimization, check that global values are physically reasonable:

```julia
optimized_δ = traj[:δ][1]
@assert 0.05 < optimized_δ < 0.25 "Detuning outside expected range"
```

## Comparison with SamplingProblem

| Feature | Global Variables | SamplingProblem |
|---------|-----------------|-----------------|
| Parameter role | Optimized | Fixed (sampled) |
| Goal | Find optimal parameter | Robust to parameter variation |
| Use case | Calibration, design | Uncertainty handling |

For robustness to uncertainty, use `SamplingProblem`. For finding optimal parameters, use global variables.

## Limitations

- Requires custom integrator support
- More complex optimization landscape
- Convergence can be slower

## See Also

- [SamplingProblem](@ref sampling) - For parameter robustness
- [Quantum Systems](@ref) - System construction
- [Problem Templates](@ref) - Main optimization API
