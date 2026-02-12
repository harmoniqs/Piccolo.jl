# # [Quantum Systems](@id quantum-systems)
#
# Quantum systems in Piccolo.jl represent the physical hardware you're
# controlling. They encapsulate the Hamiltonian structure and control bounds.
#
# ## The Hamiltonian Model
#
# Piccolo.jl uses the standard quantum control Hamiltonian:
#
# ```math
# H(u, t) = H_{\text{drift}} + \sum_{i=1}^{n} u_i(t) H_{\text{drive},i}
# ```
#
# Where:
# - `H_drift`: The always-on system Hamiltonian (e.g., qubit frequencies, couplings)
# - `H_drive,i`: The i-th controllable interaction (e.g., microwave drives)
# - `u_i(t)`: The control amplitude for drive `i` at time `t`
#
# ## QuantumSystem
#
# `QuantumSystem` is the primary type for closed quantum systems.
#
# ### Matrix-Based Construction
#
# The most common way to create a system:

using Piccolo

## Single qubit with X and Y drives
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
# Drive bounds specify the maximum control amplitude for each drive channel:
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
# For systems with dissipation, use `OpenQuantumSystem`:
#
# ```julia
# # Collapse operators (Lindblad form)
# c_ops = [
#     sqrt(γ1) * annihilate(levels),  # Energy relaxation (T1)
#     sqrt(γ2) * PAULIS[:Z]           # Pure dephasing (T2)
# ]
#
# open_sys = OpenQuantumSystem(H_drift, H_drives, drive_bounds, c_ops)
# ```
#
# The dynamics follow the Lindblad master equation:
#
# ```math
# \dot{\rho} = -i[H, \rho] + \sum_k \left( L_k \rho L_k^\dagger - \frac{1}{2}\{L_k^\dagger L_k, \rho\} \right)
# ```
#
# ## CompositeQuantumSystem
#
# For multi-qubit or multi-subsystem setups:
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
# Piccolo.jl provides standard quantum operators:
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

levels = 5
a = annihilate(levels)
a_dag = create(levels)
n_op = a_dag * a  # Number operator

# Number operator (5 levels):
n_op

# ## System Templates
#
# Piccolo.jl provides pre-built templates for common physical systems:
#
# ### Transmon Qubits

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
# Hamiltonians should be Hermitian. Piccolo.jl validates this:
#
# ```julia
# # This will warn if H is not Hermitian
# sys = QuantumSystem(H_drift, H_drives, bounds)
# ```
#
# ### 2. Normalize Units
#
# Use consistent units throughout. A common choice:
# - Energy/frequency: GHz (or 2π × GHz)
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
