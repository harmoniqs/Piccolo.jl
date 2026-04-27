# ```@copybutton
# literate/concepts/isomorphisms.jl
# ```
#
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
# ## Density Matrix Isomorphisms
#
# ### The Compact Representation
#
# Density matrices `ρ` are Hermitian (`ρ = ρ†`), so a `d×d` density matrix has
# only `d²` independent real parameters — not `2d²`. Piccolo.jl exploits this
# with a **compact isomorphism** that halves the state dimension compared to
# the naive `[Re(vec(ρ)); Im(vec(ρ))]` approach.
#
# The compact vector `x ∈ ℝᵈ²` packs the upper triangle of `ρ`:
# 1. **Real parts of the upper triangle** (column-major, `j ≤ k`): `d(d+1)/2` entries
# 2. **Imaginary parts of the strict upper triangle** (`j < k`): `d(d-1)/2` entries
#
# Total: `d(d+1)/2 + d(d-1)/2 = d²`

using SparseArrays  #hide

## Example: compact isomorphism for a 3-level density matrix
ρ = ComplexF64[
    0.5 0.1+0.2im 0.0+0.1im
    0.1-0.2im 0.3 0.05+0.0im
    0.0-0.1im 0.05 0.2
]

x = density_to_compact_iso(ρ)
length(x)  ## d² = 9 (not 2d² = 18)

#-

## Round-trip: compact iso back to density matrix
ρ_recovered = compact_iso_to_density(x)
ρ ≈ ρ_recovered

# ### Lift and Projection Matrices
#
# The compact representation is related to the full `2d²`-dimensional iso_vec
# by a pair of sparse matrices:
#
# - **Lift matrix** `L` (`2d² × d²`): maps compact → full, reconstructing the
#   lower triangle from Hermiticity
# - **Projection matrix** `P` (`d² × 2d²`): maps full → compact, extracting
#   the upper triangle
#
# These satisfy `P * L = I(d²)` (left inverse).
#
# ```math
# \underbrace{\tilde{\rho}_{\text{full}}}_{2d^2} = L \underbrace{x}_{d^2},
# \qquad
# x = P \, \tilde{\rho}_{\text{full}}
# ```

d_ex = 3
L = density_lift_matrix(d_ex)
P = density_projection_matrix(d_ex)
(size(L), size(P))

#-

## P * L is the identity
P * L ≈ I(d_ex^2)

# ### Why This Matters
#
# For open quantum system optimization, the Lindbladian generators `𝒢` are
# reduced from `(2d²)² → (d²)²` entries via `𝒢_compact = P * 𝒢 * L`. For a
# cat qubit with `d = 14` (cat ⊗ buffer), this cuts the state vector from 392
# to 196 entries and the generator from ~150k to ~38k entries — roughly a 4×
# speedup in integration.
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
cached_solve!(qcp, "isomorphisms_example"; max_iter = 50)
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
# ### Density Matrix Conversions
#
# | Function | Description |
# |----------|-------------|
# | `density_to_compact_iso(ρ)` | Hermitian matrix → compact real vector (`d²`) |
# | `compact_iso_to_density(x)` | Compact real vector → Hermitian matrix |
# | `density_lift_matrix(d)` | Sparse `L` matrix: compact → full (`2d² × d²`) |
# | `density_projection_matrix(d)` | Sparse `P` matrix: full → compact (`d² × 2d²`) |
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
# | Density `ρ` (compact) | `d×d` Hermitian | `d²` real |

## Verify dimensions for a 2-level system
d = 2
length(ket_to_iso(zeros(ComplexF64, d)))           ## 2d real elements

#-

length(operator_to_iso_vec(zeros(ComplexF64, d, d))) ## 2d² real elements

#-

ρ_test = ComplexF64[1 0; 0 0]
length(density_to_compact_iso(ρ_test))              ## d² real elements

# ## See Also
#
# - [Trajectories](@ref trajectories-concept) - How isomorphisms are used in trajectories
# - [Operators](@ref operators-concept) - Working with quantum operators
