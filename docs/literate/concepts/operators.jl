# # [Operators](@id operators-concept)
#
# Piccolo.jl provides specialized operator types for working with composite and
# multilevel quantum systems.
#
# ## Overview
#
# | Type / Function | Purpose |
# |------|---------|
# | `EmbeddedOperator` | Gates acting on subspace of larger system |
# | `lift_operator` | Operators lifted to composite Hilbert spaces |
# | `direct_sum` | Block diagonal operator combinations |
#
# ## EmbeddedOperator
#
# `EmbeddedOperator` represents a gate that acts on a computational subspace
# within a larger Hilbert space. This is essential for multilevel systems like
# transmons.
#
# ### Construction

using Piccolo

## From a gate symbol
sys = TransmonSystem(levels = 3, δ = 0.2, drive_bounds = [0.2, 0.2])
U_X = EmbeddedOperator(:X, sys)

## From a matrix
T_gate = [1 0; 0 exp(im * π / 4)]
U_T = EmbeddedOperator(T_gate, sys)

# ### Supported Gate Symbols
#
# Standard single-qubit gates:
# - `:I`, `:X`, `:Y`, `:Z` - Pauli gates
# - `:H` - Hadamard
# - `:T`, `:S` - Phase gates
# - `:sqrtX`, `:sqrtY` - Square root gates
#
# Two-qubit gates (for composite systems):
# - `:CX`, `:CNOT` - Controlled-NOT
# - `:CZ` - Controlled-Z
# - `:SWAP` - SWAP
# - `:sqrtiSWAP` - √iSWAP
#
# ### Properties

U_goal = EmbeddedOperator(:X, sys)

## Access the full operator matrix
size(U_goal.operator)

## Computational subspace
U_goal.subspace

# ### How Embedding Works
#
# For a 3-level system with computational subspace {|0⟩, |1⟩}:
#
# ```math
# X_{\text{embedded}} = \begin{pmatrix}
# 0 & 1 & 0 \\
# 1 & 0 & 0 \\
# 0 & 0 & 1
# \end{pmatrix}
# ```
#
# The gate acts as X on the first two levels and as identity on the third.
U_goal.operator

# ### Use with Trajectories

T_dur, N = 10.0, 100
times = collect(range(0, T_dur, length = N))
pulse = ZeroOrderPulse(0.01 * randn(2, N), times)

## EmbeddedOperator is accepted as goal
qtraj = UnitaryTrajectory(sys, pulse, U_goal)

## Solve with leakage handling
opts = PiccoloOptions(leakage_constraint = true, leakage_constraint_value = 1e-3)
qcp = SmoothPulseProblem(qtraj, N; piccolo_options = opts)
solve!(qcp; max_iter = 100)
println("Embedded gate fidelity: ", fidelity(qcp))

# ## lift_operator
#
# `lift_operator` lifts an operator from a subsystem to a composite Hilbert
# space via tensor products.
#
# ### Construction

## Operator on qubit 1 of a 2-qubit system
X_on_q1 = lift_operator(PAULIS[:X], 1, [2, 2])  # X ⊗ I
X_on_q1

## Operator on qubit 2
Z_on_q2 = lift_operator(PAULIS[:Z], 2, [2, 2])  # I ⊗ Z
Z_on_q2

# ### Arguments
#
# | Argument | Type | Description |
# |----------|------|-------------|
# | `op` | `Matrix` | Operator to lift |
# | `subsystem` | `Int` | Which subsystem (1-indexed) |
# | `dims` | `Vector{Int}` | Dimensions of all subsystems |
#
# ### Example: Two-Qubit Hamiltonian

dims = [2, 2]  # Two qubits

## Individual qubit terms
Z1 = lift_operator(PAULIS[:Z], 1, dims)
Z2 = lift_operator(PAULIS[:Z], 2, dims)
X1 = lift_operator(PAULIS[:X], 1, dims)
X2 = lift_operator(PAULIS[:X], 2, dims) # nothing

## Build Hamiltonian: H = ω1*Z1 + ω2*Z2 + J*Z1*Z2
ω1, ω2, J = 1.0, 1.1, 0.05
H_drift_2q = ω1 * Z1 + ω2 * Z2 + J * Z1 * Z2
H_drift_2q

# ## direct_sum
#
# `direct_sum` creates block diagonal operators, useful for parallel operations
# on independent subspaces.
#
# ### Construction

A = [1 0; 0 -1]
B = [0 1; 1 0]
C = direct_sum(A, B)
C

# ### Use Case
#
# Direct sums are useful for:
# - Operations on multiple independent qubits
# - Block diagonal Hamiltonians
# - Composite systems with no coupling
#
# ## Common Patterns
#
# ### Mixed Subsystems
#
# ```julia
# # Qubit coupled to resonator (cavity)
# # Qubit: 2 levels, Cavity: 10 levels
# dims = [2, 10]
#
# # Qubit operators
# σ_x = lift_operator(PAULIS[:X], 1, dims)
# σ_z = lift_operator(PAULIS[:Z], 1, dims)
#
# # Cavity operators
# a = lift_operator(annihilate(10), 2, dims)
# a_dag = lift_operator(create(10), 2, dims)
#
# # Jaynes-Cummings coupling
# g = 0.1
# H_coupling = g * (σ_x * (a + a_dag))
# ```
#
# ## Utility Functions
#
# ### Creating Standard Operators

## Annihilation and creation operators
levels = 5
a = annihilate(levels)
a_dag = create(levels)
n_op = a_dag * a  # Number operator
n_op

## Pauli matrices
I2, X, Y, Z = PAULIS[:I], PAULIS[:X], PAULIS[:Y], PAULIS[:Z]

# ### Tensor Products

## Manual tensor product
H_ZZ = kron(PAULIS[:Z], PAULIS[:Z])

## Or use lift_operator for clarity
Z1 = lift_operator(PAULIS[:Z], 1, [2, 2])
Z2 = lift_operator(PAULIS[:Z], 2, [2, 2])
H_ZZ_lifted = Z1 * Z2

H_ZZ ≈ H_ZZ_lifted

# ## See Also
#
# - [Leakage Suppression](@ref leakage-suppression) - Using EmbeddedOperator for leakage control
# - [Quantum Systems](@ref quantum-systems) - Building Hamiltonians with operators
# - [Isomorphisms](@ref isomorphisms-concept) - Converting operators to optimization-friendly forms
