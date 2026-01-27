# [Core Concepts](@id getting-started-concepts)

This page introduces the fundamental concepts in Piccolo.jl. For detailed documentation, see the [Concepts](@ref concepts-overview) section.

## The Quantum Optimal Control Problem

Quantum optimal control finds control pulses that steer a quantum system to achieve a desired outcome. In Piccolo.jl, this means:

1. **Define a system**: What physical hardware are you controlling?
2. **Specify a goal**: What gate or state do you want to achieve?
3. **Optimize controls**: Find pulse shapes that achieve the goal

## Key Objects

### QuantumSystem

A `QuantumSystem` represents the physical hardware:

```julia
using Piccolo

# Hamiltonian: H(t) = H_drift + u₁(t)·H_drive₁ + u₂(t)·H_drive₂
H_drift = PAULIS[:Z]  # Always-on term
H_drives = [PAULIS[:X], PAULIS[:Y]]  # Controllable terms
drive_bounds = [1.0, 1.0]  # Maximum control amplitudes

sys = QuantumSystem(H_drift, H_drives, drive_bounds)
```

### Pulse

A `Pulse` specifies how controls vary in time:

```julia
T = 10.0  # Total duration
N = 100   # Number of time points
times = collect(range(0, T, length=N))
controls = 0.1 * randn(2, N)  # Initial control values

pulse = ZeroOrderPulse(controls, times)
```

### Trajectory

A `Trajectory` combines system, pulse, and goal:

```julia
U_goal = GATES[:X]  # Target: X gate
qtraj = UnitaryTrajectory(sys, pulse, U_goal)
```

### QuantumControlProblem

A `QuantumControlProblem` sets up the optimization:

```julia
qcp = SmoothPulseProblem(qtraj, N; Q=100.0, R=1e-2)
```

## The Workflow

```
┌─────────────────┐
│ Define System   │  sys = QuantumSystem(...)
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ Create Pulse    │  pulse = ZeroOrderPulse(...)
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ Set Goal        │  qtraj = UnitaryTrajectory(sys, pulse, U_goal)
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ Create Problem  │  qcp = SmoothPulseProblem(qtraj, N)
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ Solve!          │  solve!(qcp; max_iter=100)
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ Analyze Results │  fidelity(qcp), plot(...)
└─────────────────┘
```

## Key Parameters

### Fidelity Weight (Q)

`Q` controls how much the optimizer prioritizes achieving high fidelity:

- **Q = 100** (default): Good balance
- **Q = 1000**: Strongly prioritize fidelity
- **Q = 10**: Allow more flexibility

### Regularization (R)

`R` controls pulse smoothness:

- **R = 1e-2** (default): Moderate smoothness
- **R = 0.1**: Very smooth pulses
- **R = 1e-4**: Allow aggressive controls

### Number of Timesteps (N)

`N` is the discretization resolution:

- **N = 50**: Coarse, fast optimization
- **N = 100**: Typical
- **N = 200+**: Fine control, slower

## Common Patterns

### Gate Synthesis

Most common use case - synthesize a quantum gate:

```julia
# X gate
qcp = SmoothPulseProblem(
    UnitaryTrajectory(sys, pulse, GATES[:X]),
    N
)
solve!(qcp)
```

### State Preparation

Prepare a specific quantum state:

```julia
ψ_init = [1.0, 0.0]  # |0⟩
ψ_goal = [0.0, 1.0]  # |1⟩
qcp = SmoothPulseProblem(
    KetTrajectory(sys, pulse, ψ_init, ψ_goal),
    N
)
solve!(qcp)
```

### Time-Optimal Control

Find the shortest duration:

```julia
# Enable variable timesteps
qcp_base = SmoothPulseProblem(qtraj, N; Δt_bounds=(0.01, 0.5))
solve!(qcp_base)

# Minimize time
qcp_fast = MinimumTimeProblem(qcp_base; final_fidelity=0.99)
solve!(qcp_fast)
```

### Robust Control

Optimize for parameter uncertainty:

```julia
# Solve for nominal system
qcp_base = SmoothPulseProblem(qtraj, N)
solve!(qcp_base)

# Add robustness
perturbed_systems = [sys_low, sys_nominal, sys_high]
qcp_robust = SamplingProblem(qcp_base, perturbed_systems)
solve!(qcp_robust)
```

## Next Steps

- **[Quickstart](@ref quickstart)**: Run your first optimization
- **[Problem Templates](@ref problem-templates-overview)**: Learn the main API
- **[Tutorials](@ref tutorials-overview)**: Step-by-step examples
