#!/usr/bin/env julia
#
# Standalone script to generate Wigner function visualizations for bosonic qubits
# Issue #57: Visual Design of a Bosonic Qubit in Wigner Function Space
#
# Usage: julia wigner_visualization_demo.jl
#
# This will create several PNG files showing different bosonic quantum states:
#   - coherent_state_wigner.png
#   - fock_state_wigner.png  
#   - cat_state_wigner.png
#   - bosonic_states_comparison.png
#   - cat_state_3d.png

using Piccolo
using QuantumToolbox
using CairoMakie
using NamedTrajectories

println("=== Wigner Function Visualization for Bosonic Qubits ===\n")

# Configuration
N = 30  # Fock space cutoff
OUTPUT_DIR = "wigner_figures"
mkpath(OUTPUT_DIR)

println("1. Generating Coherent State Wigner Function...")

# Coherent state |α⟩
α = 2.0 + 1.5im
ψ_coherent = coherent(N, α)

traj_coherent = NamedTrajectory(
    (ψ̃ = hcat(ket_to_iso(ψ_coherent.data)), Δt = [1.0]),
    timestep = :Δt
)

fig_coherent = plot_wigner(
    traj_coherent, 
    1; 
    state_name = :ψ̃,
    xvec = -4:0.1:4,
    yvec = -4:0.1:4
)

fig_coherent.figure[1, 1, Top()] = Label(
    fig_coherent.figure,
    "Coherent State |$(round(real(α), digits=1)) + $(round(imag(α), digits=1))i⟩",
    fontsize = 20,
    font = :bold
)

save(joinpath(OUTPUT_DIR, "coherent_state_wigner.png"), fig_coherent, px_per_unit = 2)
println("   ✓ Saved: $(OUTPUT_DIR)/coherent_state_wigner.png")

println("\n2. Generating Fock State Wigner Function...")

# Fock state |3⟩
n_photons = 3
ψ_fock = fock(N, n_photons)

traj_fock = NamedTrajectory(
    (ψ̃ = hcat(ket_to_iso(ψ_fock.data)), Δt = [1.0]),
    timestep = :Δt
)

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

save(joinpath(OUTPUT_DIR, "fock_state_wigner.png"), fig_fock, px_per_unit = 2)
println("   ✓ Saved: $(OUTPUT_DIR)/fock_state_wigner.png")

println("\n3. Generating Cat State Wigner Function...")

# Even cat state: (|α⟩ + |-α⟩)/√2
α_cat = 2.0
ψ_plus = coherent(N, α_cat)
ψ_minus = coherent(N, -α_cat)
cat_state = (ψ_plus.data + ψ_minus.data)
cat_state = cat_state / norm(cat_state)

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

save(joinpath(OUTPUT_DIR, "cat_state_wigner.png"), fig_cat, px_per_unit = 2)
println("   ✓ Saved: $(OUTPUT_DIR)/cat_state_wigner.png")

println("\n4. Generating Comparison Panel...")

# Create comparison figure
fig_comparison = Figure(size = (1400, 450), fontsize = 14)

xvec = -4:0.1:4
yvec = -4:0.1:4

# Coherent state
ax1 = Axis(
    fig_comparison[1, 1],
    xlabel = "Re(α)",
    ylabel = "Im(α)",
    title = "Coherent State",
    aspect = DataAspect()
)
W_coherent = transpose(wigner(ψ_coherent, xvec, yvec))
heatmap!(ax1, xvec, yvec, W_coherent, colormap = :RdBu)
Colorbar(fig_comparison[1, 1, Right()], limits = extrema(W_coherent), colormap = :RdBu)

# Fock state
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

# Cat state
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

Label(
    fig_comparison[0, :],
    "Wigner Functions of Bosonic Quantum States",
    fontsize = 22,
    font = :bold
)

save(joinpath(OUTPUT_DIR, "bosonic_states_comparison.png"), fig_comparison, px_per_unit = 2)
println("   ✓ Saved: $(OUTPUT_DIR)/bosonic_states_comparison.png")

println("\n5. Generating 3D Surface Plot...")

# Note: CairoMakie doesn't support 3D well, so we use a 2D contour with elevation
# For true 3D, users should use GLMakie
fig_3d_alt = Figure(size = (900, 700), fontsize = 14)

ax_contour = Axis(
    fig_3d_alt[1, 1],
    xlabel = "Re(α)",
    ylabel = "Im(α)",
    title = "Cat State Wigner Function (Contour with Elevation)",
    aspect = DataAspect()
)

xvec_3d = -3:0.12:3
yvec_3d = -3:0.12:3
W_cat_3d = transpose(wigner(ψ_cat_qobj, xvec_3d, yvec_3d))

# Contour plot with many levels to simulate 3D
contourf!(ax_contour, xvec_3d, yvec_3d, W_cat_3d, levels = 30, colormap = :RdBu)
Colorbar(fig_3d_alt[1, 2], limits = extrema(W_cat_3d), colormap = :RdBu, label = "W(α)")

Label(
    fig_3d_alt[0, :],
    "Cat State Quantum Interference Pattern",
    fontsize = 20,
    font = :bold
)

save(joinpath(OUTPUT_DIR, "cat_state_contour.png"), fig_3d_alt, px_per_unit = 2)
println("   ✓ Saved: $(OUTPUT_DIR)/cat_state_contour.png")

println("\n6. Demonstrating CatSystem Template...")

# Create cat system using Piccolo template
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

# Create target cat state using coherent_ket helper
α_target = 2.0
cat_levels = sys_cat.params[:cat_levels]
buffer_levels = sys_cat.params[:buffer_levels]

# Construct even cat state in system Hilbert space
ψ_plus_sys = kron(coherent_ket(α_target, cat_levels), [1.0; zeros(buffer_levels - 1)])
ψ_minus_sys = kron(coherent_ket(-α_target, cat_levels), [1.0; zeros(buffer_levels - 1)])
cat_sys = (ψ_plus_sys + ψ_minus_sys) / norm(ψ_plus_sys + ψ_minus_sys)

# Project to cat subspace for visualization
cat_subspace = cat_sys[1:cat_levels]

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

save(joinpath(OUTPUT_DIR, "cat_system_wigner.png"), fig_cat_sys, px_per_unit = 2)
println("   ✓ Saved: $(OUTPUT_DIR)/cat_system_wigner.png")

println("\n" * "="^70)
println("✓ All visualizations generated successfully!")
println("  Output directory: $(OUTPUT_DIR)/")
println("\nGenerated files:")
println("  • coherent_state_wigner.png - Classical Gaussian state")
println("  • fock_state_wigner.png - Quantum negativity example")
println("  • cat_state_wigner.png - Interference fringes")
println("  • bosonic_states_comparison.png - Side-by-side comparison")
println("  • cat_state_contour.png - Elevated contour plot")
println("  • cat_system_wigner.png - Using CatSystem template")
println("\nThese figures are publication-quality and ready for:")
println("  - Research papers (save as PDF)")
println("  - Presentations (existing PNG)")
println("  - Web/blog posts (convert to SVG)")
println("  - Harmoniqs website")
println("="^70)
