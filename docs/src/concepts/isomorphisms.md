# Isomorphisms

Piccolo.jl uses real isomorphisms to convert complex quantum states and operators into real vectors suitable for numerical optimization.

## Why Isomorphisms?

Optimization algorithms work with real numbers. Quantum states and unitaries are complex, so we need to:
1. Convert complex objects to real vectors for optimization
2. Convert back to complex form for physics calculations

## State Isomorphisms

### Ket States

A complex ket state `|ψ⟩` is converted to a real vector:

```julia
using Piccolo

# Complex ket
ψ = ComplexF64[1, im] / √2

# Convert to isomorphic form
ψ̃ = ket_to_iso(ψ)  # [Re(ψ); Im(ψ)]

# Convert back
ψ_recovered = iso_to_ket(ψ̃)
```

### Mathematical Form

```math
|\psi\rangle = \begin{pmatrix} a + ib \\ c + id \end{pmatrix}
\quad\rightarrow\quad
\tilde{\psi} = \begin{pmatrix} a \\ c \\ b \\ d \end{pmatrix}
```

The isomorphism stacks real parts followed by imaginary parts.

### Density Matrices

For density matrices (used in open systems):

```julia
# Complex density matrix
ρ = [1.0 0.0; 0.0 0.0]  # |0⟩⟨0|

# Convert to isomorphic form
ρ̃ = density_to_iso(ρ)

# Convert back
ρ_recovered = iso_to_density(ρ̃)
```

## Operator Isomorphisms

### Unitary Operators

Unitary matrices are vectorized and converted to real form:

```julia
# Complex unitary
U = GATES[:H]  # Hadamard

# Convert to isomorphic vector
Ũ = operator_to_iso_vec(U)

# Convert back
U_recovered = iso_vec_to_operator(Ũ)
```

### Mathematical Form

The unitary `U` is first vectorized (column-major) then split into real and imaginary parts:

```math
U = \begin{pmatrix} U_{11} & U_{12} \\ U_{21} & U_{22} \end{pmatrix}
\rightarrow
\text{vec}(U) = \begin{pmatrix} U_{11} \\ U_{21} \\ U_{12} \\ U_{22} \end{pmatrix}
\rightarrow
\tilde{U} = \begin{pmatrix} \text{Re}(\text{vec}(U)) \\ \text{Im}(\text{vec}(U)) \end{pmatrix}
```

## Hamiltonian Isomorphisms

The Schrödinger equation becomes a real linear system under isomorphism.

### Generator Form

For the Schrödinger equation:
```math
i\frac{d|\psi\rangle}{dt} = H|\psi\rangle
```

The isomorphic form uses a real generator `G`:
```math
\frac{d\tilde{\psi}}{dt} = G(H)\tilde{\psi}
```

Where `G(H)` is the isomorphic generator:

```julia
H = PAULIS[:Z]  # Hamiltonian

# Get isomorphic generator
G = hamiltonian_to_iso_generator(H)
```

### Mathematical Form

```math
G(H) = \begin{pmatrix} -\text{Im}(H) & -\text{Re}(H) \\ \text{Re}(H) & -\text{Im}(H) \end{pmatrix}
```

This transforms the complex linear ODE into a real linear ODE.

## Using Isomorphisms in Practice

### Accessing Trajectory Data

Trajectories store isomorphic states:

```julia
# After solving
traj = get_trajectory(qcp)

# Isomorphic unitary at timestep k
Ũ_k = traj[:Ũ⃗][:, k]

# Convert to complex unitary
U_k = iso_vec_to_operator(Ũ_k)
```

### Computing Fidelity

```julia
# Get final isomorphic state
Ũ_final = traj[:Ũ⃗][:, end]

# Convert to unitary
U_final = iso_vec_to_operator(Ũ_final)

# Compute fidelity with target
U_goal = GATES[:X]
d = size(U_goal, 1)
F = abs(tr(U_goal' * U_final))^2 / d^2
```

### Custom State Analysis

```julia
# Analyze ket trajectory
for k in 1:size(traj[:ψ̃], 2)
    ψ̃_k = traj[:ψ̃][:, k]
    ψ_k = iso_to_ket(ψ̃_k)

    # Compute populations
    populations = abs2.(ψ_k)
    println("Step $k: ", populations)
end
```

## Function Reference

### Ket Conversions

| Function | Description |
|----------|-------------|
| `ket_to_iso(ψ)` | Complex ket → real vector |
| `iso_to_ket(ψ̃)` | Real vector → complex ket |

### Operator Conversions

| Function | Description |
|----------|-------------|
| `operator_to_iso_vec(U)` | Complex operator → real vector |
| `iso_vec_to_operator(Ũ)` | Real vector → complex operator |

### Density Matrix Conversions

| Function | Description |
|----------|-------------|
| `density_to_iso(ρ)` | Complex density matrix → real vector |
| `iso_to_density(ρ̃)` | Real vector → complex density matrix |

### Hamiltonian Conversions

| Function | Description |
|----------|-------------|
| `hamiltonian_to_iso_generator(H)` | Hamiltonian → real generator |

## Variable Naming Convention

Piccolo.jl uses a tilde notation to distinguish isomorphic variables:

| Physical | Isomorphic | Meaning |
|----------|------------|---------|
| `ψ` | `ψ̃` | Ket state |
| `U` | `Ũ` | Unitary operator |
| `ρ` | `ρ̃` | Density matrix |

In trajectories:
- `:ψ̃` - Isomorphic ket state
- `:Ũ⃗` - Isomorphic vectorized unitary
- `:ρ̃` - Isomorphic vectorized density matrix

## Dimension Reference

For a system with `d` levels:

| Object | Complex Dimension | Isomorphic Dimension |
|--------|-------------------|---------------------|
| Ket `\|ψ⟩` | `d` complex | `2d` real |
| Unitary `U` | `d×d` complex | `2d²` real |
| Density `ρ` | `d×d` complex | `2d²` real |

## See Also

- [Trajectories](@ref) - How isomorphisms are used in trajectories
- [Operators](@ref) - Working with quantum operators
