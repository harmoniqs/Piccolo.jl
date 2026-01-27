# [System Templates](@id system-templates)

Piccolo.jl provides pre-built system templates for common physical platforms. These templates handle the Hamiltonian construction with physically meaningful parameters.

## Available Templates

| Template | Physical System |
|----------|-----------------|
| `TransmonSystem` | Superconducting transmon qubits |
| `MultiTransmonSystem` | Multiple coupled transmons |
| `IonChainSystem` | Linear chain of trapped ions |
| `RadialMSGateSystem` | Mølmer-Sørensen gates for ions |
| `RydbergChainSystem` | Rydberg atom arrays |
| `CatSystem` | Bosonic cat qubits in cavities |

## TransmonSystem

For superconducting transmon qubits with anharmonicity.

### Basic Usage

```julia
using Piccolo

# 3-level transmon with X and Y drives
sys = TransmonSystem(
    levels = 3,              # Number of energy levels
    δ = 0.2,                 # Anharmonicity (GHz)
    drive_bounds = [0.2, 0.2]  # Max amplitude for X, Y drives
)
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `levels` | `Int` | required | Number of transmon levels (≥2) |
| `δ` | `Float64` | required | Anharmonicity (typically negative, ~-0.2 GHz) |
| `drive_bounds` | `Vector{Float64}` | required | Bounds for X and Y drives |
| `ω` | `Float64` | `0.0` | Qubit frequency (often set to 0 in rotating frame) |

### Hamiltonian Structure

```math
H = \omega a^\dagger a + \frac{\delta}{2} a^\dagger a (a^\dagger a - 1) + u_x(t) (a + a^\dagger) + u_y(t) i(a^\dagger - a)
```

### Example: X Gate on 3-Level Transmon

```julia
using Piccolo

# Create system
sys = TransmonSystem(levels=3, δ=0.2, drive_bounds=[0.2, 0.2])

# Goal: X gate embedded in computational subspace
U_goal = EmbeddedOperator(:X, sys)

# Setup and solve
T, N = 20.0, 100
times = collect(range(0, T, length=N))
pulse = ZeroOrderPulse(0.05 * randn(2, N), times)
qtraj = UnitaryTrajectory(sys, pulse, U_goal)

qcp = SmoothPulseProblem(qtraj, N; Q=100.0)
solve!(qcp; max_iter=100)
```

## MultiTransmonSystem

For multiple coupled transmon qubits.

### Usage

```julia
sys = MultiTransmonSystem(
    n_qubits = 2,
    levels = 3,
    δs = [0.2, 0.22],        # Anharmonicities (can differ)
    J = 0.01,                 # Coupling strength
    drive_bounds = 0.2
)
```

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `n_qubits` | `Int` | Number of qubits |
| `levels` | `Int` | Levels per qubit |
| `δs` | `Vector{Float64}` | Anharmonicities |
| `J` | `Float64` | ZZ coupling strength |
| `drive_bounds` | `Float64` or `Vector` | Control bounds |

## IonChainSystem

For trapped ion quantum computing.

### Usage

```julia
sys = IonChainSystem(
    n_ions = 2,
    ω_motional = 1.0,    # Motional mode frequency (MHz)
    η = 0.1,              # Lamb-Dicke parameter
    drive_bounds = 0.5
)
```

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `n_ions` | `Int` | Number of ions |
| `ω_motional` | `Float64` | Motional mode frequency |
| `η` | `Float64` | Lamb-Dicke parameter |
| `drive_bounds` | `Float64` | Laser drive bounds |

## RadialMSGateSystem

Specialized for Mølmer-Sørensen entangling gates.

### Usage

```julia
sys = RadialMSGateSystem(
    n_ions = 2,
    ω_modes = [1.0, 1.1],  # Radial mode frequencies
    ηs = [0.1, 0.08],      # Lamb-Dicke parameters
    drive_bounds = 0.5
)
```

## RydbergChainSystem

For Rydberg atom arrays.

### Usage

```julia
sys = RydbergChainSystem(
    n_atoms = 3,
    Ω_max = 1.0,          # Maximum Rabi frequency
    Δ_range = (-5.0, 5.0), # Detuning range
    positions = [0.0, 1.0, 2.0]  # Atom positions (for interaction)
)
```

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `n_atoms` | `Int` | Number of atoms |
| `Ω_max` | `Float64` | Maximum Rabi frequency |
| `Δ_range` | `Tuple` | Detuning bounds |
| `positions` | `Vector` | Atom positions (affects Rydberg interactions) |

### Hamiltonian Structure

```math
H = \sum_i \frac{\Omega(t)}{2} \sigma_x^i - \Delta(t) n_i + \sum_{i<j} \frac{C_6}{|r_i - r_j|^6} n_i n_j
```

## CatSystem

For bosonic cat qubits in superconducting cavities.

### Usage

```julia
sys = CatSystem(
    n_photons = 10,       # Photon number cutoff
    α = 2.0,              # Cat state amplitude
    κ = 0.01,             # Two-photon dissipation rate
    drive_bounds = 0.5
)
```

## Creating Custom Templates

You can create your own system templates:

```julia
function MyCustomSystem(; param1, param2, drive_bounds)
    levels = compute_levels(param1)

    # Build Hamiltonian components
    H_drift = build_drift(param1, param2, levels)
    H_drives = build_drives(param1, levels)

    return QuantumSystem(H_drift, H_drives, drive_bounds)
end
```

## Best Practices

### 1. Include Enough Levels

For transmons, always include at least one level above the computational space:

```julia
# For single-qubit gates: 3 levels minimum
sys = TransmonSystem(levels=3, δ=0.2, drive_bounds=0.2)

# For two-qubit gates: 3 levels per qubit
sys = MultiTransmonSystem(n_qubits=2, levels=3, ...)
```

### 2. Use Realistic Parameters

Match your physical system:

```julia
# Typical transmon parameters
sys = TransmonSystem(
    levels = 4,
    δ = -0.2,      # Negative anharmonicity, ~200 MHz
    drive_bounds = [0.05, 0.05]  # ~50 MHz max drive
)
```

### 3. Combine with EmbeddedOperator

For multilevel systems, define gates in the computational subspace:

```julia
sys = TransmonSystem(levels=3, δ=0.2, drive_bounds=0.2)
U_goal = EmbeddedOperator(:X, sys)  # X gate on |0⟩, |1⟩ subspace
```

## See Also

- [Quantum Systems](@ref quantum-systems) - General system documentation
- [Leakage Suppression](@ref leakage-suppression) - Handling higher levels
- [Operators](@ref operators-concept) - EmbeddedOperator details
