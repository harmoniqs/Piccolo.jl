# [Quantum Module](@id quantum-api)

The Quantum module provides all quantum mechanical building blocks: systems, trajectories, pulses, operators, and isomorphisms.

## Quantum Systems

### Types

- `QuantumSystem` - Standard closed quantum system with drift and drive Hamiltonians
- `OpenQuantumSystem` - Open system with Lindblad collapse operators
- `CompositeQuantumSystem` - Multi-component systems with couplings
- `VariationalQuantumSystem` - Systems with optimizable parameters

### Constructor

```julia
sys = QuantumSystem(H_drift, H_drives, drive_bounds)
sys = QuantumSystem(H_drift, H_drives, drive_bounds; time_dependent=false)
```

### Functions

- `get_drift(sys)` - Get the drift Hamiltonian
- `get_drives(sys)` - Get the drive Hamiltonians
- `get_c_ops(sys)` - Get collapse operators (for open systems)

## System Templates

Pre-built system constructors for common hardware:

- `TransmonSystem` - Superconducting transmon qubit
- `MultiTransmonSystem` - Multiple coupled transmons
- `IonChainSystem` - Trapped ion chain
- `RadialMSGateSystem` - Molmer-Sorensen gate system
- `RydbergChainSystem` - Rydberg atom array
- `CatSystem` - Cat qubit in a cavity

## Quantum Trajectories

### Types

- `UnitaryTrajectory` - For unitary gate synthesis
- `KetTrajectory` - For state-to-state transfer
- `DensityTrajectory` - For open system evolution
- `MultiKetTrajectory` - For multiple state transfers (gate via states)
- `SamplingTrajectory` - For robust optimization over system variants

### Constructor Examples

```julia
# Unitary gate
qtraj = UnitaryTrajectory(sys, pulse, U_goal)

# State transfer
qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

# Multiple state mappings
qtraj = MultiKetTrajectory(sys, pulse, [ψ0, ψ1], [ψ1, ψ0])
```

### Functions

- `fidelity(qtraj)` - Compute fidelity with goal
- `rollout(qtraj, pulse)` - Create new trajectory with different pulse
- `rollout!(qtraj, pulse)` - Update trajectory in-place
- `get_pulse(qtraj)` - Get the pulse from trajectory
- `get_system(qtraj)` - Get the quantum system
- `state_name(qtraj)` - Get the state variable name (e.g., `:Ũ⃗`, `:ψ̃`)
- `drive_name(qtraj)` - Get the control variable name (e.g., `:u`)

## Pulses

### Types

- `ZeroOrderPulse` - Piecewise constant (use with `SmoothPulseProblem`)
- `LinearSplinePulse` - Linear interpolation between knots
- `CubicSplinePulse` - Cubic spline with derivative continuity
- `GaussianPulse` - Analytic Gaussian envelope
- `CompositePulse` - Combine multiple pulses

### Constructor Examples

```julia
# Piecewise constant
pulse = ZeroOrderPulse(controls, times)

# Linear spline
pulse = LinearSplinePulse(knot_values, knot_times)

# Gaussian
pulse = GaussianPulse(amplitudes, sigma, duration)
```

### Functions

- `duration(pulse)` - Get total pulse duration
- `n_drives(pulse)` - Get number of control drives
- `sample(pulse, times)` - Sample pulse at given times
- `pulse(t)` - Evaluate pulse at time t (callable)

## Operators

### EmbeddedOperator

For gates on a subspace of a larger Hilbert space:

```julia
op = EmbeddedOperator(:X, sys)  # X gate on qubit subspace
op = EmbeddedOperator(U_gate, sys, subspace_indices)
```

### Functions

- `lift_operator(op, dims)` - Lift operator to larger space
- `direct_sum(A, B)` - Direct sum of operators
- `embed(op)` - Get full-space operator
- `unembed(U, op)` - Extract subspace operator

### Constants

- `PAULIS` - Dictionary of Pauli matrices (`:X`, `:Y`, `:Z`, `:I`)
- `GATES` - Dictionary of common gates (`:X`, `:Y`, `:Z`, `:H`, `:T`, `:S`, `:CNOT`, etc.)

### Utility Functions

- `annihilate(n)` - Bosonic annihilation operator for n levels
- `create(n)` - Bosonic creation operator for n levels

## Isomorphisms

Functions for converting between complex quantum objects and real vectors (for optimization):

### Ket States

```julia
ψ̃ = ket_to_iso(ψ)      # Complex ket → real isomorphic vector
ψ = iso_to_ket(ψ̃)      # Real vector → complex ket
```

### Operators/Unitaries

```julia
Ũ⃗ = operator_to_iso_vec(U)    # Complex operator → real vector
U = iso_vec_to_operator(Ũ⃗)    # Real vector → complex operator
```

### Density Matrices

```julia
ρ̃ = density_to_iso_vec(ρ)     # Density matrix → real vector
ρ = iso_vec_to_density(ρ̃)     # Real vector → density matrix
```

## Dynamics

ODE problem types for quantum evolution:

- `KetODEProblem` - State vector evolution
- `UnitaryODEProblem` - Unitary propagator evolution
- `DensityODEProblem` - Density matrix evolution (Lindblad)

These are used internally by trajectories for rollout computations.
