# # [Quantum Systems](@id quantum-systems)
#
# A quantum system in Piccolo.jl specifies the Hamiltonian and the hardware
# constraints. It is the first ingredient of every optimization problem.
#
# ## The Hamiltonian Model
#
# The most common model is the **bilinear control Hamiltonian**:
#
# ```math
# H(\boldsymbol{u}, t) = H_{\text{drift}} + \sum_{i=1}^{m} u_i(t)\, H_{\text{drive},i}
# ```
#
# where ``H_{\text{drift}} \in \mathbb{C}^{d \times d}`` is the always-on
# (free-evolution) Hamiltonian, ``\{H_{\text{drive},i}\}_{i=1}^{m}`` are the
# controllable interactions, and ``\boldsymbol{u}(t) \in \mathbb{R}^m`` are the
# time-dependent control amplitudes.
#
# Piccolo.jl also supports a **generalized drive model** where each drive term
# has an arbitrary scalar coefficient function of the control vector:
#
# ```math
# H(\boldsymbol{u}, t) = H_{\text{drift}} + \sum_{d} c_d(\boldsymbol{u})\, H_d
# ```
#
# This includes the standard linear case (``c_d(\boldsymbol{u}) = u_i``) as well
# as nonlinear drives (``c_d(\boldsymbol{u}) = u_1^2 + u_2^2``, etc.) that arise
# in displaced frames, cross-Kerr couplings, and other physics.
#
# For optimization the Hamiltonian is converted to a **generator** of the
# Schrödinger equation:
#
# ```math
# G(\boldsymbol{u}) = -i\!\left( H_{\text{drift}} + \sum_{d} c_d(\boldsymbol{u})\, H_d \right)
# ```
#
# so that ``\dot{\psi} = G\,\psi`` (ket evolution) or ``\dot{U} = G\,U``
# (unitary propagator). The generator is stored internally in isomorphic (real) form.
# See [Isomorphisms](@ref isomorphisms-concept) for how complex matrices become
# real generators.
#
# ## QuantumSystem
#
# `QuantumSystem` is the primary type for closed quantum systems.
#
# ### Construction

using Piccolo

## Single qubit: H(u) = σ_z + u₁ σ_x + u₂ σ_y
H_drift = PAULIS[:Z]                    # ωq/2 σz
H_drives = [PAULIS[:X], PAULIS[:Y]]     # Controllable terms
drive_bounds = [1.0, 1.0]               # Maximum amplitude for each drive

sys = QuantumSystem(H_drift, H_drives, drive_bounds)

# ### Constructor Variants
#
# ```julia
# # Full specification (linear drives)
# sys = QuantumSystem(H_drift, H_drives, drive_bounds)
#
# # No drift (pure control)
# sys = QuantumSystem(H_drives, drive_bounds)
#
# # No drives (free evolution)
# sys = QuantumSystem(H_drift)
#
# # Typed drives (supports nonlinear coefficients)
# sys = QuantumSystem(H_drift, drives, drive_bounds)
# ```
#
# ## Drive Bounds
#
# Drive bounds set the box constraints ``u_{\min} \leq u_i \leq u_{\max}``
# for each channel:
## Vector: per-drive bounds
sys_vector = QuantumSystem(H_drift, H_drives, [0.5, 1.0])

# ## Accessing System Properties

## Number of energy levels
sys.levels

## Number of control drives
sys.n_drives

## Get Hamiltonian components
H_d = get_drift(sys)
H_dr = get_drives(sys)
H_d

# ## Nonlinear Drives
#
# When the Hamiltonian contains terms whose coefficient is a nonlinear function of
# the controls — for example ``|\alpha|^2`` terms in a displaced frame or products
# of control amplitudes — you can use **typed drive terms** via the
# `QuantumSystem(H_drift, drives, drive_bounds)` constructor.
#
# ### Drive Types
#
# Each drive pairs a Hermitian operator ``H_d`` with a scalar coefficient
# ``c_d(\boldsymbol{u})``:
#
# | Type | Coefficient | Use case |
# |------|-------------|----------|
# | `LinearDrive(H, i)` | ``u_i`` | Standard bilinear control |
# | `NonlinearDrive(H, f)` | ``f(\boldsymbol{u})`` | Displaced frames, cross-Kerr, products |
# | `NonlinearDrive(H, f, ∂f)` | ``f(\boldsymbol{u})`` | Same, with hand-written Jacobian |
#
# The **2-argument** `NonlinearDrive(H, f)` constructor is recommended — it
# auto-generates the Jacobian via ForwardDiff:

using SparseArrays

## Auto-Jacobian (recommended):
d_auto = NonlinearDrive(PAULIS[:Z], u -> u[1]^2 + u[2]^2)

## Check: ∂(u₁² + u₂²)/∂u₁ at u = [3, 4] should be 6
drive_coeff_jac(d_auto, [3.0, 4.0, 0.0], 1)

# The **3-argument** form `NonlinearDrive(H, f, ∂f)` is available when you want
# full control over the Jacobian (e.g., for performance-critical inner loops):
#
# ```julia
# NonlinearDrive(
#     PAULIS[:Z],
#     u -> u[1]^2 + u[2]^2,                             # coefficient
#     (u, j) -> j == 1 ? 2u[1] : j == 2 ? 2u[2] : 0.0  # hand-written Jacobian
# )
# ```
#
# Hand-written Jacobians are **validated against ForwardDiff** during
# `QuantumSystem` construction — sign errors or off-by-one bugs are caught
# immediately.
#
# ### Structural Sparsity
#
# When a nonlinear drive depends on only a few control indices (e.g., ``u_3^2 + u_4^2``
# in a 10-control system), you can declare `active_controls` to skip unnecessary
# Jacobian evaluations in the sensitivity equations:

d_sparse = NonlinearDrive(
    PAULIS[:Z],
    u -> u[1]^2 + u[2]^2;
    active_controls = [1, 2],   ## only u₁ and u₂ affect this drive
)
active_controls(d_sparse)

# For `LinearDrive`, `active_controls` is automatically `[index]`.
#
# ### Example: Mixed Linear and Nonlinear Drives
#
# Consider a qubit in a displaced frame with controls ``\boldsymbol{u} = (u_1, u_2)``
# and a Hamiltonian:
#
# ```math
# H(\boldsymbol{u}) = u_1\,\sigma_x + u_2\,\sigma_y + (u_1^2 + u_2^2)\,\sigma_z
# ```

drives = AbstractDrive[
    LinearDrive(sparse(ComplexF64.(PAULIS[:X])), 1),     ## u₁ σx
    LinearDrive(sparse(ComplexF64.(PAULIS[:Y])), 2),     ## u₂ σy
    NonlinearDrive(PAULIS[:Z], u -> u[1]^2 + u[2]^2),   ## (u₁² + u₂²) σz
]

sys_nl = QuantumSystem(zeros(ComplexF64, 2, 2), drives, [1.0, 1.0])

# The resulting system works with all the same trajectories and problem templates:

sys_nl.n_drives   ## control dimension (2, not 3 — drives share controls)

#-

length(sys_nl.drives)  ## number of drive terms (3)

#-

has_nonlinear_drives(sys_nl.drives)

# Note that `n_drives` is the control dimension (length of ``\boldsymbol{u}``),
# which may differ from the number of drive terms when nonlinear drives combine
# multiple control channels.
#
# ### Accessing Drive Terms
#
# Use `get_drive_terms(sys)` to access the full `AbstractDrive` objects with
# coefficient functions, Jacobians, and active control indices:

get_drive_terms(sys_nl)

# Use `get_drives(sys)` for just the Hamiltonian matrices (one per drive term):

get_drives(sys_nl)

# ### Backward Compatibility
#
# The standard matrix-based constructor `QuantumSystem(H_drift, H_drives, bounds)`
# still works exactly as before — it automatically creates `LinearDrive` objects
# internally:

sys_linear = QuantumSystem(PAULIS[:Z], [PAULIS[:X], PAULIS[:Y]], [1.0, 1.0])
sys_linear.drives  ## auto-populated LinearDrives

# ## OpenQuantumSystem
#
# For systems with dissipation the dynamics follow the Lindblad master equation:
#
# ```math
# \dot{\rho} = \underbrace{-i[H,\,\rho]}_{\text{unitary}} + \underbrace{\sum_k \left( L_k \rho L_k^\dagger - \tfrac{1}{2}\{L_k^\dagger L_k,\, \rho\} \right)}_{\text{dissipation}}
# ```
#
# Vectorizing ``\rho`` turns this into a linear ODE ``\dot{\vec\rho} = \mathcal{G}\,\vec\rho``
# with the Lindbladian superoperator
#
# ```math
# \mathcal{G}(\boldsymbol{u}) = \mathcal{G}_{\text{drift}} + \sum_{d} c_d(\boldsymbol{u})\, \mathcal{G}_{d}
# ```
#
# which has the same generalized drive structure as the closed-system generator.
#
# ```julia
# L_ops = [
#     sqrt(γ1) * annihilate(levels),  # Energy relaxation (T1)
#     sqrt(γ2) * PAULIS[:Z]           # Pure dephasing (T2)
# ]
#
# # Matrix-based (linear drives only)
# open_sys = OpenQuantumSystem(
#     H_drift, H_drives, drive_bounds;
#     dissipation_operators=L_ops
# )
#
# # Typed drives (supports nonlinear coefficients)
# open_sys = OpenQuantumSystem(
#     H_drift, drives, drive_bounds;
#     dissipation_operators=L_ops
# )
#
# # From an existing QuantumSystem (preserves nonlinear drives)
# open_sys = OpenQuantumSystem(qsys; dissipation_operators=L_ops)
# ```
#
# ### Compact Lindbladian Generators
#
# Since density matrices are Hermitian (``\rho = \rho^\dagger``), only
# ``d^2`` real parameters are independent — not ``2d^2``. The function
# `compact_lindbladian_generators(sys)` returns generators
# ``\mathcal{G}_c = P\,\mathcal{G}\,L`` of size ``d^2 \times d^2`` (instead
# of ``2d^2 \times 2d^2``), where ``L`` and ``P`` are the lift and projection
# matrices from the compact isomorphism.
# See [Isomorphisms](@ref isomorphisms-concept) for the full construction.
#
# ## CompositeQuantumSystem
#
# For multi-qubit or multi-subsystem setups the total Hilbert space is a
# tensor product ``\mathcal{H} = \mathcal{H}_1 \otimes \mathcal{H}_2``:
#
# ```julia
# sys1 = QuantumSystem(H1_drift, H1_drives, bounds1)
# sys2 = QuantumSystem(H2_drift, H2_drives, bounds2)
# H_coupling = J * kron(PAULIS[:Z], PAULIS[:Z])
#
# composite_sys = CompositeQuantumSystem([sys1, sys2], H_coupling)
# ```
#
# ## Common Gates and Operators
#
# ### Pauli Matrices

PAULIS[:I]  ## Identity

#-

PAULIS[:X]  ## Pauli-X (NOT)

#-

PAULIS[:Y]  ## Pauli-Y

#-

PAULIS[:Z]  ## Pauli-Z

# ### Common Gates
(:I, :X, :Y, :Z, :H, :T, :S, :CX, :CZ)

# ### Creation/Annihilation Operators
#
# For a ``d``-level system the annihilation operator is
# ``a = \sum_{n=1}^{d-1} \sqrt{n}\,|n{-}1\rangle\langle n|``:

levels = 5
a = annihilate(levels)
a_dag = create(levels)
n_op = a_dag * a  # Number operator

# Number operator (5 levels):
n_op

# ## System Templates
#
# Pre-built templates for common physical platforms.
#
# ### Transmon Qubits
#
# The transmon Hamiltonian (in the rotating frame, truncated to ``d`` levels):
#
# ```math
# H_{\text{transmon}} = \sum_{n=0}^{d-1} \left(n\omega_q + \tfrac{n(n-1)}{2}\,\delta\right)|n\rangle\langle n|
# ```
#
# where ``\delta`` is the anharmonicity.

sys_transmon = TransmonSystem(levels = 3, δ = 0.2, drive_bounds = [0.2, 0.2])
sys_transmon.levels, sys_transmon.n_drives

# ### Other Templates
#
# ```julia
# # Trapped ions
sys = IonChainSystem(N_ions = 2, ωq = 1.0, η = 0.1)
#
# # Rydberg atoms
sys = RydbergChainSystem(N = 3, drive_bounds = [1.0])
# ```
#
# See [System Templates](@ref system-templates) in the How-To Guides for detailed usage.
#
# ## Best Practices
#
# ### 1. Check Hermiticity
#
# All Hamiltonians must be Hermitian (``H = H^\dagger``). Piccolo.jl
# validates this at construction.
#
# ### 2. Normalize Units
#
# Use consistent units throughout. A common choice:
# - Energy/frequency: GHz (or ``2\pi \times`` GHz)
# - Time: nanoseconds
# - Control amplitudes: GHz
#
# ### 3. Include All Relevant Levels
#
# For transmon qubits, include at least one level above the computational
# subspace to capture leakage:

## Good: 3 levels for qubit (captures leakage to |2⟩)
sys_3level = TransmonSystem(levels = 3, δ = 0.2)
sys_3level.levels

## Better for high-fidelity: 4+ levels
sys_4level = TransmonSystem(levels = 4, δ = 0.2)
sys_4level.levels

# ## See Also
#
# - [Trajectories](@ref trajectories-concept) - Combining systems with pulses and goals
# - [Problem Templates](@ref problem-templates-overview) - Setting up optimization problems
# - [System Templates](@ref system-templates) - Pre-built physical system models
