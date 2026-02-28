#!/usr/bin/env julia
#
# Pulse Parameter Animation Demo
# Issue #59: Animation of a Pulse Evolving with Parameters
#
# Generates several animation files demonstrating pulse evolution:
#   - pulse_progressive.mp4       — Pulse drawing itself progressively
#   - pulse_with_populations.mp4  — Pulse + quantum state evolution
#   - amplitude_sweep.mp4         — Varying pulse amplitude
#   - chirp_pulse.mp4             — Frequency-swept (chirped) pulse
#   - pulse_sequence.mp4          — Multi-gate composite sequence
#
# Usage:
#   julia pulse_animation_demo.jl
#
# Optional: pass an output directory as the first argument
#   julia pulse_animation_demo.jl my_output_dir
#
# Dependencies: Piccolo, CairoMakie
#   using Pkg; Pkg.add(["Piccolo", "CairoMakie"])

using Piccolo
using CairoMakie
using Random

Random.seed!(42)

# ─── Configuration ────────────────────────────────────────────────────────────

const OUTPUT_DIR = length(ARGS) > 0 ? ARGS[1] : "pulse_animations"
const FPS_FAST   = 24
const FPS_SLOW   = 10

mkpath(OUTPUT_DIR)

function outpath(name)
    joinpath(OUTPUT_DIR, name)
end

println("=== Pulse Parameter Animation Demo ===")
println("Output directory: $(OUTPUT_DIR)\n")

# ─── Quantum System Setup ─────────────────────────────────────────────────────

println("Setting up quantum system...")

H_drift  = 0.5 * PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
sys      = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

const T_GATE = 10.0
const N_GATE = 100
const TIMES  = collect(range(0, T_GATE; length = N_GATE))

# ─── Optimize a Reference X-Gate Pulse ───────────────────────────────────────

println("Optimizing reference X-gate pulse (this may take ~30 s)...")

initial_controls = 0.1 * randn(2, N_GATE)
pulse  = ZeroOrderPulse(initial_controls, TIMES)
qtraj  = UnitaryTrajectory(sys, pulse, GATES[:X])
qcp    = SmoothPulseProblem(qtraj, N_GATE; Q = 100.0, R = 1e-2, ddu_bound = 1.0)

solve!(qcp; max_iter = 50, verbose = false, print_level = 1)

traj       = get_trajectory(qcp)
plot_times = cumsum([0; get_timesteps(traj)])[1:(end - 1)]

@info "Reference pulse ready" fidelity = round(fidelity(qcp); digits = 4)

# ─── Helpers ─────────────────────────────────────────────────────────────────

"""Return matrix of |U_{i,j}|² for each time step (columns of Ũ⃗)."""
function extract_populations(traj, col = 1)
    Ũ⃗ = traj[:Ũ⃗]
    n  = size(Ũ⃗, 2)
    pop = Matrix{Float64}(undef, 2, n)
    for i in 1:n
        U = iso_vec_to_operator(Ũ⃗[:, i])
        pop[1, i] = abs2(U[1, col])
        pop[2, i] = abs2(U[2, col])
    end
    return pop
end

"""Style an axis for control-pulse plots."""
function style_control_axis!(ax)
    ax.xlabel = "Time (arb. units)"
    ax.ylabel = "Control amplitude"
    hlines!(ax, [0.0]; color = (:black, 0.2), linewidth = 1)
end

"""Style an axis for population plots."""
function style_population_axis!(ax)
    ax.xlabel = "Time (arb. units)"
    ax.ylabel = "Population |Uᵢⱼ|²"
    ylims!(ax, 0, 1)
    hlines!(ax, [0.0, 1.0]; color = (:gray, 0.4), linestyle = :dash, linewidth = 1)
end

# ─── Animation 1: Progressive Pulse Reveal ───────────────────────────────────

println("\n[1/5] Progressive pulse reveal...")

animate_name(
    traj, :u;
    fps      = FPS_FAST,
    mode     = :record,
    filename = outpath("pulse_progressive.mp4"),
)

println("      ✓  pulse_progressive.mp4")

# ─── Animation 2: Pulse + Quantum Populations ────────────────────────────────

println("[2/5] Pulse + quantum populations...")

pops = extract_populations(traj)

fig2 = Figure(size = (1100, 750), fontsize = 13)
Label(fig2[0, 1], "Control Pulse & Unitary Population Evolution"; fontsize = 18, font = :bold)

ax_ctrl = Axis(fig2[1, 1]; title = "Control Pulses")
ax_pop  = Axis(fig2[2, 1]; title = "Unitary Populations (column 1)")

style_control_axis!(ax_ctrl)
style_population_axis!(ax_pop)

xlims!(ax_ctrl, plot_times[1], plot_times[end])
xlims!(ax_pop,  plot_times[1], plot_times[end])

line_ux  = lines!(ax_ctrl, Float64[], Float64[]; label = "uₓ (I)",    linewidth = 2, color = :royalblue)
line_uy  = lines!(ax_ctrl, Float64[], Float64[]; label = "uᵧ (Q)",    linewidth = 2, color = :crimson)
line_p00 = lines!(ax_pop,  Float64[], Float64[]; label = "|U₁₁|² (stay |0⟩)", linewidth = 2, color = :seagreen)
line_p10 = lines!(ax_pop,  Float64[], Float64[]; label = "|U₂₁|² (flip |1⟩)", linewidth = 2, color = :darkorange)

axislegend(ax_ctrl; position = :rt)
axislegend(ax_pop;  position = :rt)

function update2!(k)
    ts = plot_times[1:k]
    line_ux[1][]  = ts;  line_ux[2][]  = traj[:u][1, 1:k]
    line_uy[1][]  = ts;  line_uy[2][]  = traj[:u][2, 1:k]
    line_p00[1][] = ts;  line_p00[2][] = pops[1, 1:k]
    line_p10[1][] = ts;  line_p10[2][] = pops[2, 1:k]
end

animate_figure(fig2, 1:traj.N, update2!;
    fps      = FPS_FAST,
    mode     = :record,
    filename = outpath("pulse_with_populations.mp4"),
)

println("      ✓  pulse_with_populations.mp4")

# ─── Animation 3: Amplitude Sweep ─────────────────────────────────────────────

println("[3/5] Amplitude sweep...")

n_amp       = 40
amplitudes  = range(0.3, 1.8; length = n_amp)
base_ctrl   = traj[:u]
amp_fids    = Vector{Float64}(undef, n_amp)
amp_ctrls   = Vector{Matrix{Float64}}(undef, n_amp)

for (i, amp) in enumerate(amplitudes)
    sc = amp .* base_ctrl
    t  = UnitaryTrajectory(sys, ZeroOrderPulse(sc, TIMES), GATES[:X])
    amp_fids[i]  = unitary_fidelity(t)
    amp_ctrls[i] = sc
end

fig3 = Figure(size = (1100, 850), fontsize = 13)
Label(fig3[0, 1], "Amplitude Sweep"; fontsize = 18, font = :bold)

ax_p3  = Axis(fig3[1, 1]; title = "Control Pulses")
ax_f3  = Axis(fig3[2, 1]; title = "Gate Fidelity vs. Amplitude Scaling",
              xlabel = "Amplitude scaling factor", ylabel = "Fidelity")

style_control_axis!(ax_p3)
xlims!(ax_p3, plot_times[1], plot_times[end])
ylims!(ax_f3, 0, 1.05)

# Static fidelity curve
lines!(ax_f3, collect(amplitudes), amp_fids; linewidth = 2, color = :black)
hlines!(ax_f3, [0.99]; color = (:green, 0.5), linestyle = :dash, linewidth = 1, label = "0.99 threshold")
axislegend(ax_f3; position = :lb)

line_ux3   = lines!(ax_p3, plot_times, amp_ctrls[1][1, :]; label = "uₓ", linewidth = 2, color = :royalblue)
line_uy3   = lines!(ax_p3, plot_times, amp_ctrls[1][2, :]; label = "uᵧ", linewidth = 2, color = :crimson)
scat_amp   = scatter!(ax_f3, [amplitudes[1]], [amp_fids[1]]; markersize = 16, color = :red)
amp_lbl    = Label(fig3[1, 1, Top()], ""; fontsize = 14, halign = :right, padding = (0, 10, 5, 0))

axislegend(ax_p3; position = :rt)

function update3!(k)
    line_ux3[2][]  = amp_ctrls[k][1, :]
    line_uy3[2][]  = amp_ctrls[k][2, :]
    scat_amp[1][]  = [amplitudes[k]]
    scat_amp[2][]  = [amp_fids[k]]
    amp_lbl.text[] = "scale = $(round(amplitudes[k]; digits=2))   fidelity = $(round(amp_fids[k]; digits=4))"
end

animate_figure(fig3, 1:n_amp, update3!;
    fps      = FPS_SLOW,
    mode     = :record,
    filename = outpath("amplitude_sweep.mp4"),
)

println("      ✓  amplitude_sweep.mp4")

# ─── Animation 4: Chirped (Frequency-Swept) Pulse ────────────────────────────

println("[4/5] Chirped pulse...")

function make_chirp_pulse(ω_start, ω_end, T, N; amplitude = 0.5)
    ts  = collect(range(0, T; length = N))
    Δt  = T / N
    ω_t = @. ω_start + (ω_end - ω_start) * ts / T
    φ   = cumsum(ω_t .* Δt)
    u_x = amplitude .* cos.(φ)
    u_y = amplitude .* sin.(φ)
    return [u_x'; u_y'], ts
end

chirp_ctrl, chirp_ts = make_chirp_pulse(0.5, 4.0, T_GATE, N_GATE)
ω_inst(t) = 0.5 + (4.0 - 0.5) * t / T_GATE   # instantaneous frequency

fig4 = Figure(size = (1100, 600), fontsize = 13)
Label(fig4[0, 1], "Frequency-Swept (Chirped) Pulse"; fontsize = 18, font = :bold)

ax4 = Axis(fig4[1, 1])
style_control_axis!(ax4)
xlims!(ax4, chirp_ts[1], chirp_ts[end])

line_cx = lines!(ax4, Float64[], Float64[]; label = "uₓ (I)", linewidth = 2, color = :royalblue)
line_cy = lines!(ax4, Float64[], Float64[]; label = "uᵧ (Q)", linewidth = 2, color = :crimson)
freq_lbl = Label(fig4[1, 1, Top()], ""; fontsize = 14, halign = :right, padding = (0, 10, 5, 0))

axislegend(ax4; position = :rt)

function update4!(k)
    ts = chirp_ts[1:k]
    line_cx[1][] = ts;  line_cx[2][] = chirp_ctrl[1, 1:k]
    line_cy[1][] = ts;  line_cy[2][] = chirp_ctrl[2, 1:k]
    ω = ω_inst(chirp_ts[k])
    freq_lbl.text[] = "t = $(round(chirp_ts[k]; digits=1)) s   ω(t) = $(round(ω; digits=2)) rad/s"
end

animate_figure(fig4, 1:N_GATE, update4!;
    fps      = FPS_FAST,
    mode     = :record,
    filename = outpath("chirp_pulse.mp4"),
)

println("      ✓  chirp_pulse.mp4")

# ─── Animation 5: Multi-Gate Pulse Sequence ───────────────────────────────────

println("[5/5] Multi-gate pulse sequence (X → Y → X)...")

const SEQ_TARGETS    = [GATES[:X], GATES[:Y], GATES[:X]]
const SEQ_GATE_NAMES = ["X gate", "Y gate", "X gate (echo)"]
const T_SEQ  = 5.0
const N_SEQ  = 50
const TS_SEQ = collect(range(0, T_SEQ; length = N_SEQ))

seq_trajs = map(SEQ_TARGETS) do target
    ic   = 0.1 * randn(2, N_SEQ)
    pl   = ZeroOrderPulse(ic, TS_SEQ)
    qt   = UnitaryTrajectory(sys, pl, target)
    prob = SmoothPulseProblem(qt, N_SEQ; Q = 100.0, R = 1e-2, ddu_bound = 1.0)
    solve!(prob; max_iter = 50, verbose = false, print_level = 1)
    get_trajectory(prob)
end

# Build frame index
seq_gate_idx  = Int[]
seq_local_idx = Int[]
for (g, t) in enumerate(seq_trajs), k in 1:t.N
    push!(seq_gate_idx,  g)
    push!(seq_local_idx, k)
end

gate_colors = [:royalblue, :seagreen, :darkorange]

# Pre-compute per-gate time vectors and x-axis limits (avoid repeated allocations
# inside the hot animation loop).
seq_ts    = [cumsum([0; get_timesteps(tr)])[1:(end - 1)] for tr in seq_trajs]
seq_xlims = [maximum(cumsum([0; get_timesteps(tr)])) for tr in seq_trajs]

fig5 = Figure(size = (1100, 600), fontsize = 13)
Label(fig5[0, 1], "Multi-Gate Pulse Sequence"; fontsize = 18, font = :bold)

ax5 = Axis(fig5[1, 1])
style_control_axis!(ax5)

# Use an Observable for color so Makie tracks the attribute correctly.
line_color_obs = Observable{Symbol}(:royalblue)
line_sx  = lines!(ax5, Float64[], Float64[]; label = "uₓ", linewidth = 2, color = line_color_obs)
line_sy  = lines!(ax5, Float64[], Float64[]; label = "uᵧ", linewidth = 2, color = :crimson)
gate_lbl = Label(fig5[1, 1, Top()], ""; fontsize = 15, font = :bold,
                 halign = :left, padding = (10, 0, 5, 0))

axislegend(ax5; position = :rt)

function update5!(frame)
    g  = seq_gate_idx[frame]
    k  = seq_local_idx[frame]
    tr = seq_trajs[g]
    ts = seq_ts[g]

    line_sx[1][] = ts[1:k];  line_sx[2][] = tr[:u][1, 1:k]
    line_sy[1][] = ts[1:k];  line_sy[2][] = tr[:u][2, 1:k]

    # Update color via its Observable — the correct Makie pattern.
    line_color_obs[] = gate_colors[g]

    gate_lbl.text[]  = "$(SEQ_GATE_NAMES[g])  (gate $(g)/$(length(SEQ_TARGETS)))"
    gate_lbl.color[] = gate_colors[g]

    xlims!(ax5, 0, seq_xlims[g])
end

animate_figure(fig5, 1:length(seq_gate_idx), update5!;
    fps      = FPS_FAST,
    mode     = :record,
    filename = outpath("pulse_sequence.mp4"),
)

println("      ✓  pulse_sequence.mp4")

# ─── Summary ──────────────────────────────────────────────────────────────────

println()
println("=" ^ 68)
println("✓ All animations saved to: $(OUTPUT_DIR)/")
println()
println("  pulse_progressive.mp4       — Pulse drawing itself step-by-step")
println("  pulse_with_populations.mp4  — Pulse + unitary population evolution")
println("  amplitude_sweep.mp4         — Fidelity landscape across amplitudes")
println("  chirp_pulse.mp4             — Frequency-modulated (chirped) pulse")
println("  pulse_sequence.mp4          — X → Y → X composite gate sequence")
println()
println("Suggested uses:")
println("  • Harmoniqs website / blog posts")
println("  • Conference talks and tutorial notebooks")
println("  • README GIFs (convert with: ffmpeg -i file.mp4 file.gif)")
println("=" ^ 68)
