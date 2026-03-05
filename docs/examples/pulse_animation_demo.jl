#!/usr/bin/env julia
#
# Standalone script to generate pulse animation demonstrations
# Issue #59: Animation of a Pulse Evolving with Parameters
#
# Usage: julia pulse_animation_demo.jl
#
# This will create several animation files:
#   - pulse_progressive.mp4 - Basic pulse drawing itself
#   - pulse_with_populations.mp4 - Pulse + quantum state evolution
#   - amplitude_sweep.mp4 - Varying pulse amplitude
#   - chirp_pulse.mp4 - Frequency-swept pulse
#   - pulse_sequence.mp4 - Multi-gate sequence

using Piccolo
using CairoMakie
using Random
Random.seed!(42)

println("=== Pulse Parameter Animation Demo ===\n")

# Configuration
OUTPUT_DIR = "pulse_animations"
mkpath(OUTPUT_DIR)

# Create base quantum system
println("Setting up quantum system...")
H_drift = 0.5 * PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
sys = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

# Optimize an X gate
T = 10.0
N = 100
times = collect(range(0, T, length=N))
initial_controls = 0.1 * randn(2, N)
pulse = ZeroOrderPulse(initial_controls, times)
qtraj = UnitaryTrajectory(sys, pulse, GATES[:X])

qcp = SmoothPulseProblem(qtraj, N; Q=100.0, R=1e-2, ddu_bound=1.0)
println("Optimizing pulse...")
solve!(qcp; max_iter=50, verbose=false, print_level=1)

traj = get_trajectory(qcp)
plot_times = cumsum([0; get_timesteps(traj)])[1:(end-1)]
println("✓ Pulse optimized (fidelity = $(round(fidelity(qcp), digits=4)))\n")

# Animation 1: Progressive pulse reveal
println("1. Creating progressive pulse animation...")
fig = animate_name(
    traj, :u; 
    fps=24, 
    mode=:record, 
    filename=joinpath(OUTPUT_DIR, "pulse_progressive.mp4")
)
println("   ✓ Saved: $(OUTPUT_DIR)/pulse_progressive.mp4")

# Animation 2: Pulse with state evolution
println("\n2. Creating pulse + populations animation...")

fig2 = Figure(size=(1200, 800))

ax_control = Axis(
    fig2[1, 1],
    xlabel="Time",
    ylabel="Control Amplitude",
    title="Control Pulses"
)

ax_pop = Axis(
    fig2[2, 1],
    xlabel="Time",
    ylabel="Population |U_{i,j}|²",
    title="Unitary Populations (Column 1)"
)

xlims!(ax_control, plot_times[1], plot_times[end])
xlims!(ax_pop, plot_times[1], plot_times[end])
ylims!(ax_pop, 0, 1)

line_u1 = lines!(ax_control, Float64[], Float64[], label="uₓ", linewidth=2, color=:blue)
line_u2 = lines!(ax_control, Float64[], Float64[], label="uᵧ", linewidth=2, color=:red)
axislegend(ax_control, position=:rt)

line_p1 = lines!(ax_pop, Float64[], Float64[], label="|U₁₁|²", linewidth=2, color=:green)
line_p2 = lines!(ax_pop, Float64[], Float64[], label="|U₂₁|²", linewidth=2, color=:orange)
axislegend(ax_pop, position=:rt)

function update_pulse_and_state!(k)
    line_u1[1][] = plot_times[1:k]
    line_u1[2][] = traj[:u][1, 1:k]
    line_u2[1][] = plot_times[1:k]
    line_u2[2][] = traj[:u][2, 1:k]
    
    Ũ⃗ = traj[:Ũ⃗]
    pops = [abs2.(iso_vec_to_operator(Ũ⃗[:, i])[:, 1]) for i in 1:k]
    
    line_p1[1][] = plot_times[1:k]
    line_p1[2][] = [p[1] for p in pops]
    line_p2[1][] = plot_times[1:k]
    line_p2[2][] = [p[2] for p in pops]
end

animate_figure(
    fig2, 1:traj.N, update_pulse_and_state!; 
    fps=24, 
    mode=:record, 
    filename=joinpath(OUTPUT_DIR, "pulse_with_populations.mp4")
)
println("   ✓ Saved: $(OUTPUT_DIR)/pulse_with_populations.mp4")

# Animation 3: Amplitude sweep
println("\n3. Creating amplitude sweep animation...")

amplitudes = range(0.5, 1.5, length=30)
fidelities = Float64[]
pulse_data = []

base_controls = traj[:u]

for amp in amplitudes
    scaled_controls = amp * base_controls
    scaled_pulse = ZeroOrderPulse(scaled_controls, times)
    scaled_traj = UnitaryTrajectory(sys, scaled_pulse, GATES[:X])
    
    fid = unitary_fidelity(scaled_traj)
    push!(fidelities, fid)
    push!(pulse_data, scaled_controls)
end

fig3 = Figure(size=(1200, 900))

ax_pulse = Axis(
    fig3[1, 1],
    xlabel="Time",
    ylabel="Control Amplitude",
    title="Control Pulses"
)

ax_fid = Axis(
    fig3[2, 1],
    xlabel="Amplitude Scaling Factor",
    ylabel="Gate Fidelity",
    title="Fidelity vs. Amplitude"
)

lines!(ax_fid, amplitudes, fidelities, linewidth=2, color=:black)
ylims!(ax_fid, 0, 1)

scatter_current = scatter!(ax_fid, [amplitudes[1]], [fidelities[1]], 
                           markersize=20, color=:red)

line_u1_sweep = lines!(ax_pulse, plot_times, pulse_data[1][1, :], 
                       label="uₓ", linewidth=2, color=:blue)
line_u2_sweep = lines!(ax_pulse, plot_times, pulse_data[1][2, :], 
                       label="uᵧ", linewidth=2, color=:red)
axislegend(ax_pulse, position=:rt)

amp_label = Label(fig3[0, 1], "Amplitude: $(round(amplitudes[1], digits=2))", 
                  fontsize=20, font=:bold)

function update_amplitude_sweep!(k)
    line_u1_sweep[2][] = pulse_data[k][1, :]
    line_u2_sweep[2][] = pulse_data[k][2, :]
    
    scatter_current[1][] = [amplitudes[k]]
    scatter_current[2][] = [fidelities[k]]
    
    amp_label.text[] = "Amplitude: $(round(amplitudes[k], digits=2)) | Fidelity: $(round(fidelities[k], digits=4))"
end

animate_figure(
    fig3, 1:length(amplitudes), update_amplitude_sweep!; 
    fps=10, 
    mode=:record, 
    filename=joinpath(OUTPUT_DIR, "amplitude_sweep.mp4")
)
println("   ✓ Saved: $(OUTPUT_DIR)/amplitude_sweep.mp4")

# Animation 4: Chirped pulse
println("\n4. Creating frequency-swept (chirp) animation...")

function create_chirp_pulse(ω₀, ω₁, T, N)
    times = range(0, T, length=N)
    Δt = T / N
    
    ω_t = [ω₀ + (ω₁ - ω₀) * t / T for t in times]
    phase = cumsum(ω_t .* Δt)
    
    amplitude = 0.5
    u_x = amplitude * cos.(phase)
    u_y = amplitude * sin.(phase)
    
    return [u_x'; u_y'], collect(times)
end

ω₀ = 1.0
ω₁ = 3.0
chirp_controls, chirp_times = create_chirp_pulse(ω₀, ω₁, T, N)

fig4 = Figure(size=(1200, 600))

ax_chirp = Axis(
    fig4[1, 1],
    xlabel="Time",
    ylabel="Control Amplitude",
    title="Frequency-Swept Pulse (Chirp)"
)

line_chirp_u1 = lines!(ax_chirp, Float64[], Float64[], 
                       label="uₓ (I)", linewidth=2, color=:blue)
line_chirp_u2 = lines!(ax_chirp, Float64[], Float64[], 
                       label="uᵧ (Q)", linewidth=2, color=:red)
axislegend(ax_chirp, position=:rt)

freq_label = Label(fig4[0, 1], "Current Frequency: $(round(ω₀, digits=2)) rad/s", 
                   fontsize=20, font=:bold)

function update_chirp!(k)
    line_chirp_u1[1][] = chirp_times[1:k]
    line_chirp_u1[2][] = chirp_controls[1, 1:k]
    line_chirp_u2[1][] = chirp_times[1:k]
    line_chirp_u2[2][] = chirp_controls[2, 1:k]
    
    t = chirp_times[k]
    ω_current = ω₀ + (ω₁ - ω₀) * t / T
    freq_label.text[] = "Current Frequency: $(round(ω_current, digits=2)) rad/s | Time: $(round(t, digits=2)) s"
end

animate_figure(
    fig4, 1:N, update_chirp!; 
    fps=24, 
    mode=:record, 
    filename=joinpath(OUTPUT_DIR, "chirp_pulse.mp4")
)
println("   ✓ Saved: $(OUTPUT_DIR)/chirp_pulse.mp4")

# Animation 5: Pulse sequence
println("\n5. Creating multi-pulse sequence animation...")

rotation_sequence = [GATES[:X], GATES[:Y], GATES[:X]]
traj_sequence = []

for target_gate in rotation_sequence
    T_gate = 5.0
    N_gate = 50
    times_gate = collect(range(0, T_gate, length=N_gate))
    
    initial_controls_seq = 0.1 * randn(2, N_gate)
    pulse_seq = ZeroOrderPulse(initial_controls_seq, times_gate)
    qtraj_seq = UnitaryTrajectory(sys, pulse_seq, target_gate)
    
    qcp_seq = SmoothPulseProblem(qtraj_seq, N_gate; Q=100.0, R=1e-2, ddu_bound=1.0)
    solve!(qcp_seq; max_iter=50, verbose=false, print_level=1)
    
    push!(traj_sequence, get_trajectory(qcp_seq))
end

fig5 = Figure(size=(1200, 600))

ax_seq = Axis(
    fig5[1, 1],
    xlabel="Time",
    ylabel="Control Amplitude",
    title="Pulse Sequence Animation"
)

gate_names = ["X Gate", "Y Gate", "X Gate"]
gate_label = Label(fig5[0, 1], "Current: $(gate_names[1])", 
                   fontsize=22, font=:bold, color=:blue)

line_seq_u1 = lines!(ax_seq, Float64[], Float64[], 
                     label="uₓ", linewidth=2, color=:blue)
line_seq_u2 = lines!(ax_seq, Float64[], Float64[], 
                     label="uᵧ", linewidth=2, color=:red)
axislegend(ax_seq, position=:rt)

total_frames = sum([t.N for t in traj_sequence])
frame_to_gate = []
frame_to_local_idx = []

for (gate_idx, traj_gate) in enumerate(traj_sequence)
    for local_idx in 1:traj_gate.N
        push!(frame_to_gate, gate_idx)
        push!(frame_to_local_idx, local_idx)
    end
end

function update_sequence!(frame)
    gate_idx = frame_to_gate[frame]
    local_idx = frame_to_local_idx[frame]
    
    traj_current = traj_sequence[gate_idx]
    times_current = cumsum([0; get_timesteps(traj_current)])[1:(end-1)]
    
    line_seq_u1[1][] = times_current[1:local_idx]
    line_seq_u1[2][] = traj_current[:u][1, 1:local_idx]
    line_seq_u2[1][] = times_current[1:local_idx]
    line_seq_u2[2][] = traj_current[:u][2, 1:local_idx]
    
    gate_label.text[] = "Current: $(gate_names[gate_idx]) (Pulse $(gate_idx)/3)"
end

animate_figure(
    fig5, 1:total_frames, update_sequence!; 
    fps=24, 
    mode=:record, 
    filename=joinpath(OUTPUT_DIR, "pulse_sequence.mp4")
)
println("   ✓ Saved: $(OUTPUT_DIR)/pulse_sequence.mp4")

println("\n" * "="^70)
println("✓ All animations generated successfully!")
println("  Output directory: $(OUTPUT_DIR)/")
println("\nGenerated files:")
println("  • pulse_progressive.mp4 - Basic progressive reveal")
println("  • pulse_with_populations.mp4 - Pulse + quantum state evolution")
println("  • amplitude_sweep.mp4 - Parameter sweep showing amplitude effect")
println("  • chirp_pulse.mp4 - Frequency-modulated pulse")
println("  • pulse_sequence.mp4 - Multi-gate composite sequence")
println("\nThese animations are suitable for:")
println("  - Presentations and talks")
println("  - Website and blog posts")
println("  - Educational materials")
println("  - Understanding optimal control")
println("="^70)
