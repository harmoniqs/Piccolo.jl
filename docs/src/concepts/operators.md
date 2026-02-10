# [Operators](@id operators-concept)

Piccolo.jl provides specialized operator types for working with composite and multilevel quantum systems.

## Overview

| Type | Purpose |
|------|---------|
| `EmbeddedOperator` | Gates acting on subspace of larger system |
| `LiftedOperator` | Operators lifted to composite Hilbert spaces |
| `DirectSum` | Block diagonal operator combinations |

## EmbeddedOperator

`EmbeddedOperator` represents a gate that acts on a computational subspace within a larger Hilbert space. This is essential for multilevel systems like transmons.

### Construction

```julia
using Piccolo

# From a gate symbol
sys = TransmonSystem(levels=3, δ=0.2, drive_bounds=[0.2, 0.2])
U_X = EmbeddedOperator(:X, sys)

# From a matrix
T_gate = [1 0; 0 exp(im*π/4)]
U_T = EmbeddedOperator(T_gate, sys)

# Explicit subspace specification
U = EmbeddedOperator(:H, sys; subspace_indices=[1, 2])
```

### Supported Gate Symbols

Standard single-qubit gates:
- `:I`, `:X`, `:Y`, `:Z` - Pauli gates
- `:H` - Hadamard
- `:T`, `:S` - Phase gates
- `:sqrtX`, `:sqrtY` - Square root gates

Two-qubit gates (for composite systems):
- `:CX`, `:CNOT` - Controlled-NOT
- `:CZ` - Controlled-Z
- `:SWAP` - SWAP
- `:sqrtiSWAP` - √iSWAP

### Properties

```julia
U_goal = EmbeddedOperator(:X, sys)

# Access the full operator matrix
full_matrix = U_goal.operator

# Computational subspace indices
comp_indices = U_goal.subspace_indices  # e.g., [1, 2]

# Leakage indices
leak_indices = U_goal.leakage_indices   # e.g., [3]

# System reference
system = U_goal.system
```

### How Embedding Works

For a 3-level system with computational subspace {|0⟩, |1⟩}:

```math
X_{\text{embedded}} = \begin{pmatrix}
0 & 1 & 0 \\
1 & 0 & 0 \\
0 & 0 & 1
\end{pmatrix}
```

The gate acts as X on the first two levels and as identity on the third.

### Use with Trajectories

```julia
sys = TransmonSystem(levels=3, δ=0.2, drive_bounds=[0.2, 0.2])
U_goal = EmbeddedOperator(:X, sys)

# Create trajectory - EmbeddedOperator is accepted as goal
qtraj = UnitaryTrajectory(sys, pulse, U_goal)

# Solve with leakage handling
opts = PiccoloOptions(leakage_constraint=true, leakage_constraint_value=1e-3)
qcp = SmoothPulseProblem(qtraj, N; piccolo_options=opts)
```

## LiftedOperator

`LiftedOperator` lifts an operator from a subsystem to a composite Hilbert space via tensor products.

### Construction

```julia
# Operator on qubit 1 of a 2-qubit system
X_on_q1 = LiftedOperator(PAULIS[:X], 1, [2, 2])  # X ⊗ I

# Operator on qubit 2
Z_on_q2 = LiftedOperator(PAULIS[:Z], 2, [2, 2])  # I ⊗ Z

# With different subsystem dimensions
# Qubit-qutrit system: qubit (2 levels) on subsystem 1, qutrit (3 levels) on subsystem 2
dims = [2, 3]
X_on_qubit = LiftedOperator(PAULIS[:X], 1, dims)
```

### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `op` | `Matrix` | Operator to lift |
| `subsystem` | `Int` | Which subsystem (1-indexed) |
| `dims` | `Vector{Int}` | Dimensions of all subsystems |

### Example: Two-Qubit Hamiltonian

```julia
dims = [2, 2]  # Two qubits

# Individual qubit terms
Z1 = LiftedOperator(PAULIS[:Z], 1, dims)
Z2 = LiftedOperator(PAULIS[:Z], 2, dims)
X1 = LiftedOperator(PAULIS[:X], 1, dims)
X2 = LiftedOperator(PAULIS[:X], 2, dims)

# Build Hamiltonian
ω1, ω2, J = 1.0, 1.1, 0.05
H_drift = ω1 * Z1 + ω2 * Z2 + J * Z1 * Z2
H_drives = [X1, X2]
```

## DirectSum

`DirectSum` creates block diagonal operators, useful for parallel operations on independent subspaces.

### Construction

```julia
# Direct sum of two operators
A = [1 0; 0 -1]
B = [0 1; 1 0]
C = DirectSum([A, B])
# Result: [1 0 0 0; 0 -1 0 0; 0 0 0 1; 0 0 1 0]
```

### Use Case

Direct sums are useful for:
- Operations on multiple independent qubits
- Block diagonal Hamiltonians
- Composite systems with no coupling

## Common Patterns

### Multi-Qubit Gates

```julia
# Two-transmon system
sys = MultiTransmonSystem(n_qubits=2, levels=3, δs=[0.2, 0.22], J=0.01)

# CZ gate on computational subspace
U_CZ = EmbeddedOperator(:CZ, sys)

qtraj = UnitaryTrajectory(sys, pulse, U_CZ)
```

### Mixed Subsystems

```julia
# Qubit coupled to resonator (cavity)
# Qubit: 2 levels, Cavity: 10 levels
dims = [2, 10]

# Qubit operators
σ_x = LiftedOperator(PAULIS[:X], 1, dims)
σ_z = LiftedOperator(PAULIS[:Z], 1, dims)

# Cavity operators
a = LiftedOperator(annihilate(10), 2, dims)
a_dag = LiftedOperator(create(10), 2, dims)

# Jaynes-Cummings coupling
g = 0.1
H_coupling = g * (σ_x * (a + a_dag))
```

### Custom Subspace Selection

```julia
# 4-level system, use levels 0 and 2 as computational subspace
sys = QuantumSystem(H_drift, H_drives, bounds)  # 4 levels

# Custom subspace (skip level 1)
U_goal = EmbeddedOperator(:X, sys; subspace_indices=[1, 3])
```

## Utility Functions

### Creating Standard Operators

```julia
# Annihilation and creation operators
levels = 5
a = annihilate(levels)
a_dag = create(levels)
n = number(levels)  # a† a

# Pauli matrices
I, X, Y, Z = PAULIS[:I], PAULIS[:X], PAULIS[:Y], PAULIS[:Z]
```

### Tensor Products

```julia
# Manual tensor product
H_ZZ = kron(PAULIS[:Z], PAULIS[:Z])

# Or use LiftedOperator for clarity
Z1 = LiftedOperator(PAULIS[:Z], 1, [2, 2])
Z2 = LiftedOperator(PAULIS[:Z], 2, [2, 2])
H_ZZ = Z1 * Z2  # Equivalent
```

## See Also

- [Leakage Suppression](@ref leakage-suppression) - Using EmbeddedOperator for leakage control
- [Quantum Systems](@ref quantum-systems) - Building Hamiltonians with operators
- [Isomorphisms](@ref isomorphisms-concept) - Converting operators to optimization-friendly forms
