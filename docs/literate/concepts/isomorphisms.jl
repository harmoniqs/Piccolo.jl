# # [Isomorphisms](@id isomorphisms-concept)
#
# Piccolo.jl uses real isomorphisms to convert complex quantum states and
# operators into real vectors suitable for numerical optimization.
#
# ## Why Isomorphisms?
#
# Optimization algorithms work with real numbers. Quantum states and unitaries
# are complex, so we need to:
# 1. Convert complex objects to real vectors for optimization
# 2. Convert back to complex form for physics calculations
#
# ## State Isomorphisms
#
# ### Ket States
#
# A complex ket state `|ψ⟩` is converted to a real vector:

using Piccolo
using LinearAlgebra

## Complex ket
ψ = ComplexF64[1, im] / √2

## Convert to isomorphic form
ψ̃ = ket_to_iso(ψ)

## Convert back
ψ_recovered = iso_to_ket(ψ̃)
ψ ≈ ψ_recovered

# ### Mathematical Form
#
# ```math
# |\psi\rangle = \begin{pmatrix} a + ib \\ c + id \end{pmatrix}
# \quad\rightarrow\quad
# \tilde{\psi} = \begin{pmatrix} a \\ c \\ b \\ d \end{pmatrix}
# ```
#
# The isomorphism stacks real parts followed by imaginary parts.
#
# ## Operator Isomorphisms
#
# ### Unitary Operators
#
# Unitary matrices are vectorized and converted to real form:

## Complex unitary
U = GATES[:H]  # Hadamard
U

## Convert to isomorphic vector
Ũ = operator_to_iso_vec(U)
length(Ũ)

## Convert back
U_recovered = iso_vec_to_operator(Ũ)
U ≈ U_recovered

# ### Mathematical Form
#
# The unitary `U` is first vectorized (column-major) then split into real and
# imaginary parts:
#
# ```math
# U = \begin{pmatrix} U_{11} & U_{12} \\ U_{21} & U_{22} \end{pmatrix}
# \rightarrow
# \text{vec}(U) = \begin{pmatrix} U_{11} \\ U_{21} \\ U_{12} \\ U_{22} \end{pmatrix}
# \rightarrow
# \tilde{U} = \begin{pmatrix} \text{Re}(\text{vec}(U)) \\ \text{Im}(\text{vec}(U)) \end{pmatrix}
# ```
#
# ## Using Isomorphisms in Practice
#
# ### Accessing Trajectory Data
#
# Trajectories store isomorphic states. Let's solve a problem and inspect:

H_drift = PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
sys = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

T, N = 10.0, 100
times = collect(range(0, T, length = N))
pulse = ZeroOrderPulse(0.1 * randn(2, N), times)
U_goal = GATES[:X]

qtraj = UnitaryTrajectory(sys, pulse, U_goal)
qcp = SmoothPulseProblem(qtraj, N; Q = 100.0)
solve!(qcp; max_iter = 50)
fidelity(qcp)

# ### Inspecting Isomorphic State

traj = get_trajectory(qcp)

## Isomorphic unitary at final timestep
Ũ_final = traj[:Ũ⃗][:, end]
length(Ũ_final)

## Convert to complex unitary
U_final = iso_vec_to_operator(Ũ_final)
size(U_final)

# ### Computing Fidelity Manually

d = size(U_goal, 1)
F = abs(tr(U_goal' * U_final))^2 / d^2
F

# ## Function Reference
#
# ### Ket Conversions
#
# | Function | Description |
# |----------|-------------|
# | `ket_to_iso(ψ)` | Complex ket → real vector |
# | `iso_to_ket(ψ̃)` | Real vector → complex ket |
#
# ### Operator Conversions
#
# | Function | Description |
# |----------|-------------|
# | `operator_to_iso_vec(U)` | Complex operator → real vector |
# | `iso_vec_to_operator(Ũ)` | Real vector → complex operator |
#
# ## Variable Naming Convention
#
# Piccolo.jl uses a tilde notation to distinguish isomorphic variables:
#
# | Physical | Isomorphic | Meaning |
# |----------|------------|---------|
# | `ψ` | `ψ̃` | Ket state |
# | `U` | `Ũ` | Unitary operator |
# | `ρ` | `ρ̃` | Density matrix |
#
# In trajectories:
# - `:ψ̃` - Isomorphic ket state
# - `:Ũ⃗` - Isomorphic vectorized unitary
# - `:ρ̃` - Isomorphic vectorized density matrix
#
# ## Dimension Reference
#
# For a system with `d` levels:
#
# | Object | Complex Dimension | Isomorphic Dimension |
# |--------|-------------------|---------------------|
# | Ket `\|ψ⟩` | `d` complex | `2d` real |
# | Unitary `U` | `d×d` complex | `2d²` real |
# | Density `ρ` | `d×d` complex | `2d²` real |

## Verify dimensions for a 2-level system
d = 2
length(ket_to_iso(zeros(ComplexF64, d)))           ## 2d real elements

#-

length(operator_to_iso_vec(zeros(ComplexF64, d, d))) ## 2d² real elements

# ## See Also
#
# - [Trajectories](@ref trajectories-concept) - How isomorphisms are used in trajectories
# - [Operators](@ref operators-concept) - Working with quantum operators
