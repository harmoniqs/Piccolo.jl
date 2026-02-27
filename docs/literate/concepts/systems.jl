# # [Quantum Systems](@id quantum-systems)
#
# A quantum system in Piccolo.jl specifies the Hamiltonian and the hardware
# constraints. It is the first ingredient of every optimization problem.
#
# ## The Hamiltonian Model
#
# Piccolo.jl uses the standard bilinear control Hamiltonian:
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
# For optimization the Hamiltonian is converted to a **generator** of the
# Schrödinger equation:
#
# ```math
# G(\boldsymbol{u}) = -i\!\left( H_{\text{drift}} + \sum_{i=1}^{m} u_i\, H_{\text{drive},i} \right)
# ```
#
# so that ``\dot{\psi} = G\,\psi`` (ket evolution) or ``\dot{U} = G\,U``
# (unitary propagator).  The generator inherits the bilinear control structure:
#
# ```math
# G(\boldsymbol{u}) = G_{\text{drift}} + \sum_{i=1}^{m} u_i\, G_{\text{drive},i}
# ```
#
# and is stored internally in isomorphic (real) form.
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
# # Full specification
# sys = QuantumSystem(H_drift, H_drives, drive_bounds)
#
# # No drift (pure control)
# sys = QuantumSystem(H_drives, drive_bounds)
#
# # No drives (free evolution)
# sys = QuantumSystem(H_drift)
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
# \mathcal{G}(\boldsymbol{u}) = \mathcal{G}_{\text{drift}} + \sum_{i=1}^{m} u_i\, \mathcal{G}_{\text{drive},i}
# ```
#
# which has the same bilinear control structure as the closed-system generator.
#
# ```julia
# L_ops = [
#     sqrt(γ1) * annihilate(levels),  # Energy relaxation (T1)
#     sqrt(γ2) * PAULIS[:Z]           # Pure dephasing (T2)
# ]
#
# open_sys = OpenQuantumSystem(
#     H_drift, H_drives, drive_bounds;
#     dissipation_operators=L_ops
# )
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
