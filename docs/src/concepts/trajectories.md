# Trajectories

Trajectories in Piccolo.jl represent a complete quantum control scenario: a system, a control pulse, and a goal. They are the central object passed to problem templates.

## Overview

A trajectory encapsulates:
- **Quantum system**: The hardware being controlled
- **Pulse**: How controls vary over time
- **Goal**: The desired outcome (gate, state, etc.)
- **State evolution**: The quantum state over time

## Trajectory Types

| Type | Use Case | State Variable |
|------|----------|----------------|
| `UnitaryTrajectory` | Gate synthesis | Unitary matrix `U(t)` |
| `KetTrajectory` | State preparation | State vector `|ψ(t)⟩` |
| `DensityTrajectory` | Open systems | Density matrix `ρ(t)` |
| `MultiKetTrajectory` | Multiple state transfers | Multiple `|ψᵢ(t)⟩` |
| `SamplingTrajectory` | Robust optimization | States for multiple systems |

## UnitaryTrajectory

For synthesizing quantum gates.

### Construction

```julia
using Piccolo

# Define system
H_drift = PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
sys = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

# Create pulse
T, N = 10.0, 100
times = collect(range(0, T, length=N))
pulse = ZeroOrderPulse(0.1 * randn(2, N), times)

# Create trajectory with goal
U_goal = GATES[:X]
qtraj = UnitaryTrajectory(sys, pulse, U_goal)
```

### Properties

```julia
# Access components
sys = qtraj.system
pulse = get_pulse(qtraj)
goal = qtraj.goal

# State at a specific time
U_t = qtraj(t)

# Final state
U_final = qtraj(duration(pulse))

# Fidelity
F = fidelity(qtraj)
```

### Rollout

After optimization, the trajectory contains the full ODE solution:

```julia
# Solve optimization
qcp = SmoothPulseProblem(qtraj, N)
solve!(qcp)

# sync_trajectory! is called automatically by solve!
# Now qtraj contains the ODE solution

# Evaluate at any time
U_halfway = qcp.qtraj(T/2)
```

## KetTrajectory

For state preparation (state-to-state transfer).

### Construction

```julia
# Initial and goal states
ψ_init = ComplexF64[1, 0]  # |0⟩
ψ_goal = ComplexF64[0, 1]  # |1⟩

qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)
```

### When to Use

Use `KetTrajectory` when:
- Preparing a specific quantum state
- Global phase doesn't matter
- Single state transfer

## MultiKetTrajectory

For gates defined by multiple state mappings with coherent phases.

### Construction

```julia
# Define state pairs
ψ0 = ComplexF64[1, 0]  # |0⟩
ψ1 = ComplexF64[0, 1]  # |1⟩

# X gate: |0⟩ → |1⟩ and |1⟩ → |0⟩
initial_states = [ψ0, ψ1]
goal_states = [ψ1, ψ0]

qtraj = MultiKetTrajectory(sys, pulse, initial_states, goal_states)
```

### Coherent Fidelity

`MultiKetTrajectory` uses `CoherentKetInfidelityObjective`, which ensures:
- Each state reaches its target
- Relative phases between states are preserved

This is important for gates where phase relationships matter (unlike `UnitaryTrajectory` which tracks the full unitary).

## DensityTrajectory

For open quantum systems with dissipation.

### Construction

```julia
# Open system with collapse operators
open_sys = OpenQuantumSystem(H_drift, H_drives, bounds, c_ops)

# Initial density matrix
ρ_init = [1.0 0.0; 0.0 0.0]  # Pure |0⟩
ρ_goal = [0.0 0.0; 0.0 1.0]  # Pure |1⟩

qtraj = DensityTrajectory(open_sys, pulse, ρ_init, ρ_goal)
```

### Note

`DensityTrajectory` support for fidelity objectives is still in development. For most open-system problems, consider using `KetTrajectory` with an effective Hamiltonian or using the `Rollouts` module directly.

## SamplingTrajectory

For robust optimization over parameter variations. This is created internally by `SamplingProblem`.

### How It Works

`SamplingTrajectory` wraps multiple system variants:

```julia
# Created internally by SamplingProblem
# You typically don't construct this directly

qcp_base = SmoothPulseProblem(qtraj, N)
solve!(qcp_base)

systems = [sys_nominal, sys_high, sys_low]
qcp_robust = SamplingProblem(qcp_base, systems)
# This creates a SamplingTrajectory internally
```

The resulting trajectory has:
- Single shared control pulse
- Multiple state trajectories (one per system)

## Common Operations

### Fidelity

Compute the fidelity of a trajectory:

```julia
F = fidelity(qtraj)  # Requires trajectory to have been rolled out

# For QuantumControlProblem (automatically synced after solve!)
F = fidelity(qcp)
```

### Rollout

Simulate the trajectory with a new pulse:

```julia
# Create new trajectory with updated pulse
new_pulse = ZeroOrderPulse(new_controls, times)
rolled_out = rollout(qtraj, new_pulse)

# In-place update
rollout!(qtraj, new_pulse)
```

### Extracting the Pulse

```julia
pulse = get_pulse(qtraj)

# Get control values at specific times
u = pulse(t)

# Get all knot values
controls = pulse.controls
times = pulse.times
```

### State Names

Trajectories have named state variables:

```julia
# Default names
state_name(UnitaryTrajectory(...))  # :Ũ⃗ (isomorphic unitary)
state_name(KetTrajectory(...))      # :ψ̃ (isomorphic ket)
state_name(DensityTrajectory(...))  # :ρ̃ (isomorphic density)
```

## Internal Representation

### Isomorphic States

Piccolo.jl uses real isomorphisms internally for optimization. Complex states are converted to real vectors:

```julia
# Complex ket |ψ⟩ → real vector ψ̃
ψ = ComplexF64[1, im] / √2
ψ̃ = ket_to_iso(ψ)  # [Re(ψ); Im(ψ)]

# Unitary U → real vector Ũ⃗
U = GATES[:H]
Ũ = operator_to_iso_vec(U)
```

See [Isomorphisms](@ref) for details.

### Named Trajectory Integration

After solving, the `NamedTrajectory` stores discrete time points:

```julia
# Access the underlying NamedTrajectory
traj = get_trajectory(qcp)

# State at timestep k
state_k = traj[:Ũ⃗, k]

# Controls at timestep k
u_k = traj[:u, k]

# All timesteps
Δts = get_timesteps(traj)
```

## Best Practices

### 1. Match Pulse Type to Problem

```julia
# For SmoothPulseProblem
pulse = ZeroOrderPulse(controls, times)
qtraj = UnitaryTrajectory(sys, pulse, U_goal)
qcp = SmoothPulseProblem(qtraj, N)  # ✓

# For SplinePulseProblem
pulse = CubicSplinePulse(controls, tangents, times)
qtraj = UnitaryTrajectory(sys, pulse, U_goal)
qcp = SplinePulseProblem(qtraj)  # ✓
```

### 2. Initialize with Reasonable Controls

Random initialization often works, but scaled appropriately:

```julia
# Scale by drive bounds
max_amp = maximum(sys.drive_bounds)
initial_controls = 0.1 * max_amp * randn(n_drives, N)
```

### 3. Use sync_trajectory! After Solving

`solve!` calls `sync_trajectory!` automatically, but if you modify the underlying optimization trajectory manually:

```julia
# Manual modification requires explicit sync
qcp.prob.trajectory[:u] .= new_values
sync_trajectory!(qcp)
```

## See Also

- [Quantum Systems](@ref) - System definitions
- [Pulses](@ref) - Control parameterizations
- [Problem Templates](@ref) - Using trajectories in optimization
