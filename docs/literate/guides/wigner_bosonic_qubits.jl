# # [Wigner Function Visualization of Bosonic Qubits](@id wigner_bosonic_qubits)
#
# This guide demonstrates publication-quality Wigner function visualizations for bosonic 
# quantum states. Bosonic qubits encode quantum information in the continuous-variable 
# states of a harmonic oscillator (e.g., a superconducting cavity). The Wigner function 
# provides a quasi-probability distribution in phase space that beautifully reveals the 
# quantum nature of these states.
#
#md # !!! tip "What You'll Learn"
#md #     - Create Wigner function plots for coherent states, Fock states, and cat states
#md #     - Understand how Wigner negativity reveals quantum features
#md #     - Generate publication-quality figures for bosonic quantum control
#md #     - Animate Wigner function evolution over time

# ## Setup

using Piccolo
using QuantumToolbox
using CairoMakie  # For static plots
using NamedTrajectories

# ## Mathematical Background
#
# The Wigner function W(α) maps a quantum state |ψ⟩ to phase space coordinates 
# (Re(α), Im(α)). For a density matrix ρ:
#
# ```math
# W(\alpha) = \frac{2}{\pi} \text{Tr}[\rho \hat{D}(\alpha) \hat{\Pi} \hat{D}^\dagger(\alpha)]
# ```
#
# where D̂(α) is the displacement operator and Π̂ is the parity operator.
#
# **Key Properties:**
# - Classical states: W(α) ≥ 0 everywhere
# - Quantum states: Can have negative regions (quantum interference)
# - Coherent states |α⟩: Gaussian peaks centered at α
# - Cat states (|α⟩ + |-α⟩)/√2: Show interference fringes

# ## Coherent States
#
# Coherent states are the most classical quantum states. Their Wigner functions are 
# Gaussian peaks.

## Parameters
N = 30  # Fock space cutoff
α = 2.0 + 1.5im  # Coherent state parameter

## Create coherent state
ψ_coherent = coherent(N, α)

## Create trajectory (single timestep for static plot)
traj_coherent = NamedTrajectory(
    (ψ̃ = hcat(ket_to_iso(ψ_coherent.data)), Δt = [1.0]),
    timestep = :Δt
)

## Plot Wigner function
fig_coherent = plot_wigner(
    traj_coherent, 
    1; 
    state_name = :ψ̃,
    xvec = -4:0.1:4,
    yvec = -4:0.1:4
)

## Customize for publication quality
fig_coherent.figure[1, 1, Top()] = Label(
    fig_coherent.figure,
    "Coherent State |$(round(real(α), digits=1)) + $(round(imag(α), digits=1))i⟩",
    fontsize = 20,
    font = :bold
)

fig_coherent

#md # !!! info "Coherent States"
#md #     The Wigner function shows a single Gaussian peak centered at α = 2.0 + 1.5i 
#md #     in phase space. These states minimize the Heisenberg uncertainty principle 
#md #     and are produced by ideal lasers.

# ## Fock States
#
# Fock states |n⟩ are energy eigenstates with exactly n photons. Higher Fock states 
# exhibit pronounced quantum features in their Wigner functions.

## Create Fock state |3⟩
n_photons = 3
ψ_fock = fock(N, n_photons)

traj_fock = NamedTrajectory(
    (ψ̃ = hcat(ket_to_iso(ψ_fock.data)), Δt = [1.0]),
    timestep = :Δt
)

## Plot with higher resolution to capture negative regions
fig_fock = plot_wigner(
    traj_fock,
    1;
    state_name = :ψ̃,
    xvec = -4:0.05:4,
    yvec = -4:0.05:4
)

fig_fock.figure[1, 1, Top()] = Label(
    fig_fock.figure,
    "Fock State |$n_photons⟩ (Quantum Negativity)",
    fontsize = 20,
    font = :bold
)

fig_fock

#md # !!! warning "Quantum Negativity"
#md #     Notice the negative regions (blue) in the Wigner function! These cannot 
#md #     occur in classical physics and are a signature of quantum behavior. Fock 
#md #     states are highly non-classical and useful for quantum sensing.

# ## Cat States
#
# Cat states (Schrödinger cat states) are superpositions of coherent states: 
# |cat⟩ = (|α⟩ + |-α⟩)/√2. They exhibit beautiful interference patterns.

## Create even cat state
α_cat = 2.0  # Real alpha for symmetry
ψ_plus = coherent(N, α_cat)
ψ_minus = coherent(N, -α_cat)

## Superposition (unnormalized initially)
cat_state = ψ_plus.data + ψ_minus.data
cat_state = cat_state / norm(cat_state)  # Normalize

traj_cat = NamedTrajectory(
    (ψ̃ = hcat(ket_to_iso(cat_state)), Δt = [1.0]),
    timestep = :Δt
)

fig_cat = plot_wigner(
    traj_cat,
    1;
    state_name = :ψ̃,
    xvec = -4:0.08:4,
    yvec = -4:0.08:4
)

fig_cat.figure[1, 1, Top()] = Label(
    fig_cat.figure,
    "Even Cat State (|α⟩ + |-α⟩)/√2 with α = $α_cat",
    fontsize = 20,
    font = :bold
)

fig_cat

#md # !!! tip "Interference Fringes"
#md #     The cat state Wigner function shows two Gaussian peaks at ±α connected by 
#md #     interference fringes. These fringes are quantum interference between the 
#md #     two coherent state components. Odd cat states (|α⟩ - |-α⟩)/√2 have different 
#md #     fringe patterns.

# ## Comparison Panel
#
# Create a multi-panel figure comparing different bosonic states.

fig_comparison = Figure(size = (1400, 450), fontsize = 14)

## Coherent state
ax1 = Axis(
    fig_comparison[1, 1],
    xlabel = "Re(α)",
    ylabel = "Im(α)",
    title = "Coherent State",
    aspect = DataAspect()
)
xvec = -4:0.1:4
yvec = -4:0.1:4
W_coherent = transpose(wigner(ψ_coherent, xvec, yvec))
heatmap!(ax1, xvec, yvec, W_coherent, colormap = :RdBu)
Colorbar(fig_comparison[1, 1, Right()], limits = extrema(W_coherent), colormap = :RdBu)

## Fock state
ax2 = Axis(
    fig_comparison[1, 2],
    xlabel = "Re(α)",
    ylabel = "Im(α)",
    title = "Fock State |3⟩",
    aspect = DataAspect()
)
W_fock = transpose(wigner(ψ_fock, xvec, yvec))
heatmap!(ax2, xvec, yvec, W_fock, colormap = :RdBu)
Colorbar(fig_comparison[1, 2, Right()], limits = extrema(W_fock), colormap = :RdBu)

## Cat state
ax3 = Axis(
    fig_comparison[1, 3],
    xlabel = "Re(α)",
    ylabel = "Im(α)",
    title = "Even Cat State",
    aspect = DataAspect()
)
ψ_cat_qobj = QuantumObject(cat_state)
W_cat = transpose(wigner(ψ_cat_qobj, xvec, yvec))
heatmap!(ax3, xvec, yvec, W_cat, colormap = :RdBu)
Colorbar(fig_comparison[1, 3, Right()], limits = extrema(W_cat), colormap = :RdBu)

## Overall title
Label(
    fig_comparison[0, :],
    "Wigner Functions of Bosonic Quantum States",
    fontsize = 22,
    font = :bold
)

fig_comparison

# ## 3D Surface Plot
#
# For presentations, 3D surface plots can be more visually striking.

using GLMakie  # Better for 3D rendering

## Create cat state Wigner data
xvec_3d = -3:0.15:3
yvec_3d = -3:0.15:3
W_cat_3d = transpose(wigner(ψ_cat_qobj, xvec_3d, yvec_3d))

## Create 3D surface
fig_3d = Figure(size = (800, 700))
ax_3d = Axis3(
    fig_3d[1, 1],
    xlabel = "Re(α)",
    ylabel = "Im(α)",
    zlabel = "W(α)",
    title = "Cat State Wigner Function (3D)",
    titlesize = 20,
    aspect = (1, 1, 0.5)
)

surface!(ax_3d, xvec_3d, yvec_3d, W_cat_3d, colormap = :RdBu, shading = true)
Colorbar(fig_3d[1, 2], limits = extrema(W_cat_3d), colormap = :RdBu, label = "W(α)")

fig_3d

#md # !!! tip "3D Visualization"
#md #     The 3D surface makes the quantum interference pattern more dramatic. 
#md #     The negative "valleys" in the Wigner function are clearly visible as 
#md #     dips below the zero plane.

# ## Animation: Rotating Cat State
#
# Animate a coherent state rotating in phase space to create a "rotating cat".

using GLMakie  # Required for animations

## Create time-evolving trajectory
T = 60
ω = 2π
times = range(0, 2π/ω, length = T)
α_amplitude = 2.5

## Generate rotating cat states
kets_rotating = map(times) do t
    α_t = α_amplitude * exp(im * ω * t)
    ψ_plus_t = coherent(N, α_t).data
    ψ_minus_t = coherent(N, -α_t).data
    cat_t = ψ_plus_t + ψ_minus_t
    cat_t / norm(cat_t)
end

traj_rotating = NamedTrajectory(
    (ψ̃ = hcat(ket_to_iso.(kets_rotating)...), Δt = fill(times[2] - times[1], T)),
    timestep = :Δt
)

## Animate (uncomment to run - requires GLMakie and interactive display)
# fig_anim = animate_wigner(
#     traj_rotating;
#     fps = 20,
#     xvec = -4:0.15:4,
#     yvec = -4:0.15:4,
#     mode = :inline  # Use :record to save to file
# )

#md # !!! info "Animation"
#md #     When running interactively with GLMakie, this creates a real-time animation 
#md #     showing the cat state rotating in phase space. The interference fringes 
#md #     rotate with the cat state components. Uncomment the animation code to see it!
#md #     Use `mode = :record, filename = "cat_rotation.mp4"` to save to file.

# ## Using CatSystem Template
#
# Piccolo provides a `CatSystem` template for cat qubit dynamics with realistic 
# Hamiltonian parameters.

using Piccolo

## Create cat system
sys_cat = CatSystem(
    g2 = 0.36,
    χ_aa = -7e-3,
    χ_bb = -32,
    χ_ab = 0.79,
    κa = 53e-3,
    κb = 13,
    cat_levels = 20,
    buffer_levels = 3
)

## Create target cat state using coherent_ket helper
α_target = 2.0
cat_levels = sys_cat.params[:cat_levels]
buffer_levels = sys_cat.params[:buffer_levels]

## Construct even cat state in system Hilbert space
ψ_plus_sys = kron(coherent_ket(α_target, cat_levels), [1.0; zeros(buffer_levels - 1)])
ψ_minus_sys = kron(coherent_ket(-α_target, cat_levels), [1.0; zeros(buffer_levels - 1)])
cat_sys = (ψ_plus_sys + ψ_minus_sys) / norm(ψ_plus_sys + ψ_minus_sys)

## Visualize (reduced subspace for cat cavity only)
## Project to cat subspace
cat_subspace = cat_sys[1:cat_levels]  # First cat_levels dimensions

traj_cat_sys = NamedTrajectory(
    (ψ̃ = hcat(ket_to_iso(cat_subspace)), Δt = [1.0]),
    timestep = :Δt
)

fig_cat_sys = plot_wigner(
    traj_cat_sys,
    1;
    state_name = :ψ̃,
    xvec = -4:0.1:4,
    yvec = -4:0.1:4
)

fig_cat_sys.figure[1, 1, Top()] = Label(
    fig_cat_sys.figure,
    "Cat State from CatSystem Template",
    fontsize = 20,
    font = :bold
)

fig_cat_sys

#md # !!! info "CatSystem"
#md #     The `CatSystem` template provides a realistic model of cat qubit hardware 
#md #     with buffer mode for two-photon processes. It includes Kerr nonlinearities 
#md #     and dissipation. See the system templates guide for optimal control examples.

# ## Saving Figures
#
# Export publication-ready figures in various formats.

## Save as PNG (high resolution)
# save("wigner_cat_state.png", fig_cat, px_per_unit = 2)

## Save as PDF (vector graphics for papers)
# save("wigner_cat_state.pdf", fig_cat)

## Save as SVG (web/presentations)
# save("wigner_cat_state.svg", fig_cat)

## Save comparison panel
# save("wigner_comparison.png", fig_comparison, px_per_unit = 2)

#md # !!! tip "Publication Tips"
#md #     - Use `px_per_unit = 2` or higher for high-DPI displays
#md #     - PDF format is best for LaTeX documents
#md #     - SVG is ideal for web presentations and blogs
#md #     - Adjust `colormap` to match your publication style (`:RdBu`, `:viridis`, etc.)

# ## Summary
#
# This guide demonstrated:
#
# - **Coherent states**: Gaussian Wigner functions (most classical)
# - **Fock states**: Negative regions indicating quantum behavior
# - **Cat states**: Interference fringes from quantum superposition
# - **Comparison plots**: Multi-panel publication figures
# - **3D visualization**: Dramatic rendering of quantum features
# - **Animation**: Time evolution of quantum states
# - **CatSystem integration**: Using Piccolo's templates for realistic dynamics
#
# The Wigner function is a powerful tool for visualizing and understanding bosonic 
# quantum states. Negative regions and interference patterns directly reveal the 
# quantum nature of these states, making them invaluable for quantum control and 
# quantum information processing.

#md # ## Next Steps
#md #
#md # - Explore [System Templates](@ref system_templates) for cat qubit control
#md # - See [Visualization Guide](@ref visualization) for more plotting options
#md # - Check [Getting Started](@ref getting_started) for optimal control examples

# ## References
#
# 1. Wigner, E. P. (1932). "On the Quantum Correction For Thermodynamic Equilibrium"
# 2. Mirrahimi, M. et al. (2014). "Dynamically protected cat-qubits: a new paradigm"
# 3. Ofek, N. et al. (2016). "Extending the lifetime of a quantum bit with error correction"
