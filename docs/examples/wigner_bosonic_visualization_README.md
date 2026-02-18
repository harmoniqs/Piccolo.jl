# Wigner Function Visualization for Bosonic Qubits

This example demonstrates publication-quality Wigner function visualizations for bosonic quantum states in Piccolo.jl.

## Overview

This addresses issue #57: Create visual designs of bosonic qubits represented in Wigner function space. The Wigner function provides a quasi-probability distribution in phase space that beautifully visualizes quantum states of harmonic oscillators.

## Features

- **Coherent States**: Gaussian Wigner functions (classical-like states)
- **Fock States**: Exhibits quantum negativity  
- **Cat States**: Shows interference fringes from quantum superposition
- **3D Visualizations**: Dramatic surface plots highlighting quantum features
- **Animations**: Time evolution of rotating cat states
- **CatSystem Integration**: Uses Piccolo's built-in templates

## Quick Start

```julia
using Piccolo
using QuantumToolbox
using CairoMakie  # For static plots
using NamedTrajectories

# Create even cat state
N = 30  # Fock space cutoff
α = 2.0  # Cat state parameter

ψ_plus = coherent(N, α)
ψ_minus = coherent(N, -α)
cat_state = (ψ_plus.data + ψ_minus.data) / norm(ψ_plus.data + ψ_minus.data)

# Plot Wigner function
traj = NamedTrajectory(
    (ψ̃ = hcat(ket_to_iso(cat_state)), Δt = [1.0]),
    timestep = :Δt
)

fig = plot_wigner(traj, 1; xvec=-4:0.08:4, yvec=-4:0.08:4)
fig.figure[1, 1, Top()] = Label(
    fig.figure,
    "Even Cat State (|α⟩ + |-α⟩)/√2",
    fontsize = 20,
    font = :bold
)
display(fig)

# Save publication-quality figure
save("cat_state_wigner.png", fig, px_per_unit = 2)
```

## Full Documentation

See the complete guide at: `docs/literate/guides/wigner_bosonic_qubits.jl`

This file contains:
- Coherent state visualizations
- Fock state quantum negativity
- Cat state interference patterns
- Multi-panel comparison plots
- 3D surface renderings
- Animation examples
- CatSystem template integration

## Building the Examples

To generate all figures from the literate guide:

```bash
cd docs
julia --project -e 'include("make.jl")'
```

The rendered markdown and figures will appear in `docs/build/generated/guides/wigner_bosonic_qubits/`.

## Publication Tips

- Use `save(filename, fig, px_per_unit = 2)` for high-DPI displays
- PDF format: `save("figure.pdf", fig)` - best for LaTeX
- SVG format: `save("figure.svg", fig)` - ideal for web/presentations
- Adjust colormaps (`:RdBu`, `:viridis`) to match publication style

## Key Insights

**Wigner Function Properties:**
- Classical states: W(α) ≥ 0 everywhere (non-negative)
- Quantum states: Can have negative regions (quantum interference)
- Coherent states: Gaussian peaks centered at α  
- Cat states: Two peaks with interference fringes
- Fock states: Negative regions indicating quantum behavior

**Why It Matters:**
- Gold standard for visualizing bosonic quantum states
- Reveals quantum-classical boundary
- Essential for cat qubit and GKP qubit research
- Shows quantum error correction capabilities

## References

1. Wigner, E. P. (1932). "On the Quantum Correction For Thermodynamic Equilibrium"
2. Mirrahimi, M. et al. (2014). "Dynamically protected cat-qubits"
3. Ofek, N. et al. (2016). "Extending the lifetime of a quantum bit with error correction"

## Issue Resolution

This implementation resolves issue #57 by providing:
- ✅ Wigner function visualization (2D heatmap and 3D surface)
- ✅ Demonstration with cat states using CatSystem
- ✅ Publication-quality figures suitable for website
- ✅ Examples using existing `coherent_ket` and `plot_wigner` helpers
- ✅ Integration with Piccolo's Makie-based visualization system

## License

MIT License - See LICENSE file in repository root.
