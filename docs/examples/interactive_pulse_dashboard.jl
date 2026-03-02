#!/usr/bin/env julia
#
# Interactive pulse parameter dashboard with real-time sliders
# Issue #59: Animation of a Pulse Evolving with Parameters
#
# Usage: julia interactive_pulse_dashboard.jl
#
# Requires: GLMakie (for interactive displays)
# Install: using Pkg; Pkg.add("GLMakie")

using Piccolo
using GLMakie
using Random
Random.seed!(42)

println("=== Interactive Pulse Parameter Dashboard ===\n")
println("Loading... This will open an interactive window.")
println("Use sliders to adjust pulse parameters in real-time!\n")

# Create quantum system
H_drift = 0.5 * PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
sys = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

# Time parameters
T = 10.0
N = 100
times = collect(range(0, T, length=N))

# Create interactive figure
fig = Figure(size=(1400, 1000), fontsize=14)

# Title
Label(fig[0, :], "Interactive Pulse Parameter Dashboard", 
      fontsize=24, font=:bold, color=:navy)

# Sliders for pulse parameters
sg_amplitude = SliderGrid(
    fig[1, 1],
    (label="Amplitude", range=0.1:0.05:2.0, startvalue=0.8),
    (label="Frequency", range=0.5:0.1:5.0, startvalue=1.5),
    (label="Phase", range=0:0.1:2π, startvalue=0.0),
    (label="Chirp Rate", range=-1.0:0.1:1.0, startvalue=0.0),
    width=400,
    tellheight=true
)

# Extract slider values
amplitude_obs = sg_amplitude.sliders[1].value
frequency_obs = sg_amplitude.sliders[2].value
phase_obs = sg_amplitude.sliders[3].value
chirp_obs = sg_amplitude.sliders[4].value

# Pulse plot
ax_pulse = Axis(
    fig[2, 1],
    xlabel="Time (arb. units)",
    ylabel="Control Amplitude",
    title="Control Pulses"
)

# Population plot
ax_pop = Axis(
    fig[3, 1],
    xlabel="Time (arb. units)",
    ylabel="Population |Uᵢⱼ|²",
    title="Unitary Populations (Column 1)"
)

ylims!(ax_pop, 0, 1)

# Fidelity display
fidelity_box = Box(fig[1, 2], color=(:lightblue, 0.3))
fidelity_label = Label(
    fig[1, 2], 
    "Fidelity: Calculating...", 
    fontsize=20, 
    font=:bold,
    padding=(10, 10, 10, 10)
)

# Info panel
info_text = Label(
    fig[2:3, 2],
    """
    Interactive Controls:
    
    • Amplitude: Scales pulse strength
    • Frequency: Oscillation rate
    • Phase: Initial phase offset
    • Chirp Rate: Frequency sweep
    
    The plots update in real-time
    as you move the sliders.
    
    Target: X gate (σₓ)
    
    Tip: Try amplitude ≈ 0.8,
    frequency ≈ 1.5 for high
    fidelity!
    """,
    fontsize=14,
    halign=:left,
    valign=:top,
    padding=(20, 20, 20, 20)
)

# Reactive pulse computation
controls_obs = lift(amplitude_obs, frequency_obs, phase_obs, chirp_obs) do amp, freq, ph, chirp
    # Generate time-varying frequency if chirp is active
    if abs(chirp) > 0.01
        ω_t = @. freq + chirp * times / T
        phase_integrated = cumsum(ω_t) * (T / N)
        u_x = amp * cos.(phase_integrated .+ ph)
        u_y = amp * sin.(phase_integrated .+ ph)
    else
        u_x = amp * sin.(2π * freq * times / T .+ ph)
        u_y = amp * cos.(2π * freq * times / T .+ ph)
    end
    
    return [u_x'; u_y']
end

# Reactive trajectory and fidelity
traj_obs = lift(controls_obs) do controls
    pulse = ZeroOrderPulse(controls, times)
    UnitaryTrajectory(sys, pulse, GATES[:X])
end

fidelity_obs = lift(traj_obs) do traj
    unitary_fidelity(traj)
end

# Update fidelity label
on(fidelity_obs) do fid
    fidelity_label.text[] = "Fidelity: $(round(fid, digits=6))"
    
    # Color code by fidelity
    if fid > 0.99
        fidelity_label.color[] = :green
    elseif fid > 0.95
        fidelity_label.color[] = :orange
    else
        fidelity_label.color[] = :red
    end
end

# Plot controls (reactive)
u_x_obs = lift(controls_obs) do controls
    controls[1, :]
end

u_y_obs = lift(controls_obs) do controls
    controls[2, :]
end

lines!(ax_pulse, times, u_x_obs, label="uₓ (I channel)", linewidth=2, color=:blue)
lines!(ax_pulse, times, u_y_obs, label="uᵧ (Q channel)", linewidth=2, color=:red)
axislegend(ax_pulse, position=:rt)

# Plot populations (reactive)
populations_obs = lift(traj_obs) do traj
    # Get populations from trajectory
    Ũ⃗ = traj.data[:Ũ⃗]
    N_traj = size(Ũ⃗, 2)
    
    pop_00 = Float64[]
    pop_10 = Float64[]
    
    for i in 1:N_traj
        U = iso_vec_to_operator(Ũ⃗[:, i])
        push!(pop_00, abs2(U[1, 1]))
        push!(pop_10, abs2(U[2, 1]))
    end
    
    (pop_00, pop_10)
end

pop_00_obs = lift(p -> p[1], populations_obs)
pop_10_obs = lift(p -> p[2], populations_obs)

lines!(ax_pop, times, pop_00_obs, label="|U₁₁|² (stay in |0⟩)", linewidth=2, color=:green)
lines!(ax_pop, times, pop_10_obs, label="|U₂₁|² (flip to |1⟩)", linewidth=2, color=:orange)
axislegend(ax_pop, position=:rt)

# Add guide lines
hlines!(ax_pop, [0.0, 1.0], color=:gray, linestyle=:dash, linewidth=1)

# Instructions
Label(
    fig[4, :],
    "Move the sliders to see real-time updates. Watch how pulse shape affects quantum populations!",
    fontsize=16,
    color=:gray
)

# Display the figure
display(fig)

println("\n" * "="^70)
println("✓ Interactive dashboard is now running!")
println("\nInstructions:")
println("  • Drag sliders to adjust parameters")
println("  • Observe real-time pulse and population updates")
println("  • Fidelity changes color: green (>0.99), orange (>0.95), red (<0.95)")
println("  • Close window to exit")
println("\nExperiment suggestions:")
println("  1. Start with default values and observe baseline")
println("  2. Increase amplitude slowly - watch fidelity improve then degrade")
println("  3. Change frequency - see oscillation rate in pulse")
println("  4. Add chirp - frequency sweeps create curved patterns")
println("  5. Adjust phase - shifts the pulse timing")
println("="^70)

# Keep running until window is closed
println("\nPress Ctrl+C to exit...")
try
    wait(fig.scene)
catch e
    println("\nDashboard closed.")
end
