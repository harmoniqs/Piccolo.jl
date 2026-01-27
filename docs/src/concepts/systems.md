# Quantum Systems

Quantum systems in Piccolo.jl represent the physical hardware you're controlling. They encapsulate the Hamiltonian structure and control bounds.

## The Hamiltonian Model

Piccolo.jl uses the standard quantum control Hamiltonian:

```math
H(u, t) = H_{\text{drift}} + \sum_{i=1}^{n} u_i(t) H_{\text{drive},i}
```

Where:
- `H_drift`: The always-on system Hamiltonian (e.g., qubit frequencies, couplings)
- `H_drive,i`: The i-th controllable interaction (e.g., microwave drives)
- `u_i(t)`: The control amplitude for drive `i` at time `t`

## QuantumSystem

`QuantumSystem` is the primary type for closed quantum systems.

### Matrix-Based Construction

The most common way to create a system:

```julia
using Piccolo

# Single qubit with X and Y drives
H_drift = PAULIS[:Z]  # ωq/2 σz
H_drives = [PAULIS[:X], PAULIS[:Y]]
drive_bounds = [1.0, 1.0]  # Maximum amplitude for each drive

sys = QuantumSystem(H_drift, H_drives, drive_bounds)
```

### Function-Based Construction

For time-dependent or parametric Hamiltonians:

```julia
# Time-dependent Hamiltonian
function H(u, t)
    ω = 1.0 + 0.1 * sin(2π * t)  # Time-varying frequency
    return ω * PAULIS[:Z] + u[1] * PAULIS[:X] + u[2] * PAULIS[:Y]
end

sys = QuantumSystem(H, [1.0, 1.0]; time_dependent=true)
```

### Constructor Variants

```julia
# Full specification
sys = QuantumSystem(H_drift, H_drives, drive_bounds)

# No drift (pure control)
sys = QuantumSystem(H_drives, drive_bounds)

# No drives (free evolution)
sys = QuantumSystem(H_drift)

# Function-based
sys = QuantumSystem(H::Function, drive_bounds)
```

## Drive Bounds

Drive bounds specify the maximum control amplitude for each drive channel:

```julia
# Scalar: same bound for all drives
sys = QuantumSystem(H_drift, H_drives, 1.0)

# Vector: per-drive bounds
sys = QuantumSystem(H_drift, H_drives, [0.5, 1.0])  # Drive 1: ±0.5, Drive 2: ±1.0

# Tuple: asymmetric bounds
sys = QuantumSystem(H_drift, H_drives, [(0.0, 1.0), (-0.5, 0.5)])
```

## Accessing System Properties

```julia
# Number of energy levels
sys.levels

# Number of control drives
sys.n_drives

# Get Hamiltonian components
H_d = get_drift(sys)
H_drives = get_drives(sys)

# Evaluate Hamiltonian at specific controls
H = sys.H(u, t)

# Get superoperator form (for density matrix evolution)
G = sys.G  # G(ρ) = -i[H, ρ]
```

## OpenQuantumSystem

For systems with dissipation, use `OpenQuantumSystem`:

```julia
# Collapse operators (Lindblad form)
c_ops = [
    sqrt(γ1) * annihilate(levels),  # Energy relaxation (T1)
    sqrt(γ2) * PAULIS[:Z]           # Pure dephasing (T2)
]

open_sys = OpenQuantumSystem(H_drift, H_drives, drive_bounds, c_ops)
```

The dynamics follow the Lindblad master equation:

```math
\dot{\rho} = -i[H, \rho] + \sum_k \left( L_k \rho L_k^\dagger - \frac{1}{2}\{L_k^\dagger L_k, \rho\} \right)
```

## CompositeQuantumSystem

For multi-qubit or multi-subsystem setups:

```julia
# Two-qubit system
sys1 = QuantumSystem(H1_drift, H1_drives, bounds1)
sys2 = QuantumSystem(H2_drift, H2_drives, bounds2)
H_coupling = J * kron(PAULIS[:Z], PAULIS[:Z])

composite_sys = CompositeQuantumSystem([sys1, sys2], H_coupling)
```

## VariationalQuantumSystem

For systems with optimizable parameters:

```julia
# System with variational parameters
sys = QuantumSystem(H_drift, H_drives, drive_bounds)
var_sys = VariationalQuantumSystem(sys, [:δ, :J])  # Optimize δ and J
```

## Common Gates and Operators

Piccolo.jl provides standard quantum operators:

### Pauli Matrices

```julia
using Piccolo

PAULIS[:I]  # Identity
PAULIS[:X]  # Pauli-X (NOT)
PAULIS[:Y]  # Pauli-Y
PAULIS[:Z]  # Pauli-Z
```

### Common Gates

```julia
GATES[:I]        # Identity
GATES[:X]        # Pauli-X
GATES[:Y]        # Pauli-Y
GATES[:Z]        # Pauli-Z
GATES[:H]        # Hadamard
GATES[:T]        # T gate
GATES[:S]        # S gate
GATES[:CX]       # CNOT
GATES[:CZ]       # Controlled-Z
GATES[:SWAP]     # SWAP
GATES[:sqrtiSWAP] # √iSWAP
```

### Creation/Annihilation Operators

```julia
# For multilevel systems
a = annihilate(levels)  # Annihilation operator
a_dag = create(levels)  # Creation operator
n = number(levels)      # Number operator
```

## System Templates

Piccolo.jl provides pre-built templates for common physical systems:

### Transmon Qubits

```julia
sys = TransmonSystem(
    levels=3,           # Number of energy levels
    δ=0.2,              # Anharmonicity (GHz)
    drive_bounds=[0.2, 0.2]  # X and Y drive bounds
)
```

### Trapped Ions

```julia
sys = IonChainSystem(
    n_ions=2,
    ω_motional=1.0,     # Motional mode frequency
    η=0.1,              # Lamb-Dicke parameter
    drive_bounds=0.5
)
```

### Rydberg Atoms

```julia
sys = RydbergChainSystem(
    n_atoms=3,
    Ω_max=1.0,          # Maximum Rabi frequency
    Δ_range=(-5.0, 5.0) # Detuning range
)
```

See [System Templates](@ref) in the How-To Guides for detailed usage.

## Best Practices

### 1. Check Hermiticity

Hamiltonians should be Hermitian. Piccolo.jl validates this:

```julia
# This will warn if H is not Hermitian
sys = QuantumSystem(H_drift, H_drives, bounds)
```

### 2. Normalize Units

Use consistent units throughout. A common choice:
- Energy/frequency: GHz (or 2π × GHz)
- Time: nanoseconds
- Control amplitudes: GHz

### 3. Include All Relevant Levels

For transmon qubits, include at least one level above the computational subspace to capture leakage:

```julia
# Good: 3 levels for qubit (captures leakage to |2⟩)
sys = TransmonSystem(levels=3, δ=0.2, drive_bounds=0.2)

# Better for high-fidelity: 4+ levels
sys = TransmonSystem(levels=4, δ=0.2, drive_bounds=0.2)
```

## See Also

- [Trajectories](@ref) - Combining systems with pulses and goals
- [Problem Templates](@ref) - Setting up optimization problems
- [System Templates](@ref) - Pre-built physical system models
