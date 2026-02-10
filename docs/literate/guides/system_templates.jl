# # [System Templates](@id system-templates)
#
# Piccolo.jl provides pre-built system templates for common physical platforms.
# These templates handle the Hamiltonian construction with physically meaningful parameters.

# ## Available Templates
#
# | Template | Physical System |
# |----------|-----------------|
# | `TransmonSystem` | Superconducting transmon qubits |
# | `MultiTransmonSystem` | Multiple coupled transmons |
# | `IonChainSystem` | Linear chain of trapped ions |
# | `RadialMSGateSystem` | Molmer-Sorensen gates for ions |
# | `RydbergChainSystem` | Rydberg atom arrays |
# | `CatSystem` | Bosonic cat qubits in cavities |

# ## TransmonSystem
#
# For superconducting transmon qubits with anharmonicity.

# ### Basic Usage

using Piccolo
using LinearAlgebra
using Random
Random.seed!(42)

## 3-level transmon with X and Y drives
sys = TransmonSystem(
    levels = 3,
    δ = 0.2,                   # Anharmonicity (GHz)
    drive_bounds = [0.2, 0.2]  # Max amplitude for X, Y drives
)

sys.levels, sys.n_drives

# ### Parameters
#
# | Parameter | Type | Default | Description |
# |-----------|------|---------|-------------|
# | `levels` | `Int` | `3` | Number of transmon levels (≥2) |
# | `δ` | `Float64` | `0.2` | Anharmonicity (typically ~0.2 GHz) |
# | `drive_bounds` | `Vector{Float64}` | `[1.0, 1.0]` | Bounds for X and Y drives |
# | `ω` | `Float64` | `4.0` | Qubit frequency (often 0 in rotating frame) |
# | `lab_frame` | `Bool` | `false` | Whether to use lab frame |
#
# ### Hamiltonian Structure
#
# ```math
# H = \omega a^\dagger a + \frac{\delta}{2} a^\dagger a (a^\dagger a - 1) + u_x(t) (a + a^\dagger) + u_y(t) i(a^\dagger - a)
# ```

# ### Example: X Gate on 3-Level Transmon

sys = TransmonSystem(levels = 3, δ = 0.2, drive_bounds = [0.2, 0.2])

## Goal: X gate embedded in computational subspace
U_goal = EmbeddedOperator(:X, sys)

## Setup and solve
T, N = 20.0, 100
times = collect(range(0, T, length = N))
pulse = ZeroOrderPulse(0.05 * randn(2, N), times)
qtraj = UnitaryTrajectory(sys, pulse, U_goal)

qcp = SmoothPulseProblem(qtraj, N; Q = 100.0)
cached_solve!(qcp, "system_templates_transmon"; max_iter = 100, verbose = false, print_level = 1)

fidelity(qcp)

# ## MultiTransmonSystem
#
# For multiple coupled transmon qubits.

# ### Usage

sys_multi = MultiTransmonSystem(
    [4.0, 4.1],               # Qubit frequencies (GHz)
    [0.2, 0.22],              # Anharmonicities (can differ)
    [0.0 0.01; 0.01 0.0];    # Coupling matrix
    drive_bounds = 0.2,
    levels_per_transmon = 3,
)

# ### Parameters
#
# | Parameter | Type | Description |
# |-----------|------|-------------|
# | `ωs` | `Vector{Float64}` | Qubit frequencies |
# | `δs` | `Vector{Float64}` | Anharmonicities |
# | `gs` | `Matrix{Float64}` | Coupling strength matrix |
# | `levels_per_transmon` | `Int` | Levels per qubit |
# | `drive_bounds` | `Float64` or `Vector` | Control bounds |

# ## IonChainSystem
#
# For trapped ion quantum computing.

# ### Usage

sys_ion = IonChainSystem(
    N_ions = 2,
    mode_levels = 5,
    ωq = 1.0,            # Qubit frequency
    ωm = 0.1,            # Motional mode frequency
    η = 0.1,             # Lamb-Dicke parameter
    drive_bounds = fill(0.5, 4),
)

# ### Parameters
#
# | Parameter | Type | Description |
# |-----------|------|-------------|
# | `N_ions` | `Int` | Number of ions |
# | `mode_levels` | `Int` | Motional mode Fock space truncation |
# | `ωq` | `Float64` | Qubit frequency |
# | `ωm` | `Float64` | Motional mode frequency |
# | `η` | `Float64` | Lamb-Dicke parameter |
# | `drive_bounds` | `Vector` | Laser drive bounds |

# ## RadialMSGateSystem
#
# Specialized for Molmer-Sorensen entangling gates.

# ### Usage

sys_ms = RadialMSGateSystem(
    N_ions = 2,
    mode_levels = 5,
    ωm_radial = [5.0, 5.0, 5.1, 5.1],  # Radial mode frequencies
    δ = 0.2,
    η = 0.1,
    drive_bounds = fill(1.0, 2),
)

# ## RydbergChainSystem
#
# For Rydberg atom arrays.

# ### Usage

sys_rydberg = RydbergChainSystem(
    N = 3,                # Number of atoms
    C = 862690 * 2π,      # Rydberg interaction coefficient
    distance = 8.7,       # Atom spacing (μm)
    cutoff_order = 1,     # Nearest-neighbor interactions
    drive_bounds = [1.0, 1.0, 1.0],
)

# ### Parameters
#
# | Parameter | Type | Description |
# |-----------|------|-------------|
# | `N` | `Int` | Number of atoms |
# | `C` | `Float64` | Rydberg interaction coefficient |
# | `distance` | `Float64` | Atom spacing (μm) |
# | `cutoff_order` | `Int` | Interaction range (1 = nearest neighbor) |
#
# ### Hamiltonian Structure
#
# ```math
# H = \sum_i \frac{\Omega(t)}{2} \sigma_x^i - \Delta(t) n_i + \sum_{i<j} \frac{C_6}{|r_i - r_j|^6} n_i n_j
# ```

# ## CatSystem
#
# For bosonic cat qubits in superconducting cavities:
#
# ```julia
# sys_cat = CatSystem(
#     cat_levels = 13,       # Photon number cutoff
#     buffer_levels = 3,     # Buffer mode levels
#     g2 = 0.36,            # Two-photon drive strength
#     drive_bounds = [1.0, 1.0],
# )
# ```
#
# !!! note
#     `CatSystem` returns an `OpenQuantumSystem` for Lindbladian dynamics.

# ## Creating Custom Templates
#
# You can create your own system templates by wrapping `QuantumSystem`:

function MyCustomSystem(; ω, δ, drive_bounds)
    levels = 3
    ## Build Hamiltonian using Piccolo's annihilation operator
    a = annihilate(levels)
    n = a' * a

    H_drift = ω * n + (δ / 2) * n * (n - I)
    H_drives = [a + a', 1.0im * (a' - a)]

    return QuantumSystem(H_drift, H_drives, drive_bounds)
end

my_sys = MyCustomSystem(ω = 4.0, δ = 0.2, drive_bounds = [0.2, 0.2])
my_sys.levels

# ## Best Practices

# ### 1. Include Enough Levels
#
# For transmons, always include at least one level above the computational space:

## For single-qubit gates: 3 levels minimum
sys = TransmonSystem(levels = 3, δ = 0.2, drive_bounds = [0.2, 0.2])

# ### 2. Use Realistic Parameters
#
# Match your physical system:

## Typical transmon parameters
sys = TransmonSystem(
    levels = 4,
    δ = 0.2,                       # ~200 MHz anharmonicity
    drive_bounds = [0.05, 0.05],   # ~50 MHz max drive
)

# ### 3. Combine with EmbeddedOperator
#
# For multilevel systems, define gates in the computational subspace:

sys = TransmonSystem(levels = 3, δ = 0.2, drive_bounds = [0.2, 0.2])
U_goal = EmbeddedOperator(:X, sys)  # X gate on |0⟩, |1⟩ subspace

# ## See Also
#
# - [Quantum Systems](@ref quantum-systems) - General system documentation
# - [Leakage Suppression](@ref leakage-suppression) - Handling higher levels
# - [Operators](@ref operators-concept) - EmbeddedOperator details
