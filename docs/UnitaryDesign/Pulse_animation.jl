#!/usr/bin/env julia
#
# docs/examples/pulse_animation_demo.jl
# 
#
# Usage:  julia --project docs/examples/pulse_animation_demo.jl
#
# Output (written to ./pulse_animations/)
#   pulse_static.png            static plot_unitary_populations
#   pulse_progressive.mp4       pulse drawing itself (animate_name)
#   pulse_with_populations.mp4  pulse + unitary populations
#   amplitude_sweep.mp4         fidelity vs. amplitude (rollout API)
#   chirp_pulse.mp4             frequency-swept pulse reveal
#
#
# ═══════════════════════════════════════════════════════════════════

using Piccolo
using CairoMakie
using Random

Random.seed!(42)

const OUTPUT_DIR = "pulse_animations"
mkpath(OUTPUT_DIR)
out(name) = joinpath(OUTPUT_DIR, name)

println("=== Pulse Animation Demo ===\n")

# ── 1. Quantum system + optimised trajectory ─

println("[1/5] Building system and optimising X gate...")

H_drift  = 0.5 * PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
sys      = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

T      = 10.0
N      = 100
# length = N 
times  = collect(range(0, T; length = N))

pulse_init = ZeroOrderPulse(0.1 * randn(2, N), times)
qtraj      = UnitaryTrajectory(sys, pulse_init, GATES[:X])
qcp        = SmoothPulseProblem(qtraj, N; Q = 100.0, R = 1e-2, ddu_bound = 1.0)
solve!(qcp; max_iter = 50, verbose = false, print_level = 1)

traj = get_trajectory(qcp)

# 
# plot_times = cumsum([0; get_timesteps(traj)])[1:(end-1)]
plot_times = cumsum([0.0; get_timesteps(traj)])[1:end-1]
controls   = traj[:u]     # (2 × N)

@info "X gate fidelity" fidelity(qcp)

# ── 2. Static snapshot — plot_unitary_populations ───
#
# plot_unitary_populations(traj) returns a Figure.
# Used for standalone static output .

println("\n[2/5] Static snapshot via plot_unitary_populations...")

fig_static = plot_unitary_populations(traj)
save(out("pulse_static.png"), fig_static)
println("      Saved: $(out("pulse_static.png"))")

# ── 3. Animation A — progressive reveal via animate_name ─
#
# 
#   fig = animate_name(traj, :u; fps=24, mode=:record, filename=...)

println("\n[3/5] Animation A — progressive pulse reveal (animate_name)...")

animate_name(traj, :u; fps = 24, mode = :record,
    filename = out("pulse_progressive.mp4"))
println("      Saved: $(out("pulse_progressive.mp4"))")

# ── 4. Animation B — pulse + populations co-evolution ──
#
# Population extraction uses the exact pattern 
#
#   Ũ⃗ = traj[:Ũ⃗]                                          
#   pops = [abs2.(iso_vec_to_operator(Ũ⃗[:, i])[:, 1])    
#            for i in 1:k]
#

println("\n[4/5] Animation B — pulse + populations (animate_figure)...")

# Pre-compute all population arrays 
Utilde = traj[:Ũ⃗]   
# :
pops_full = [abs2.(iso_vec_to_operator(Utilde[:, i])[:, 1]) for i in 1:traj.N]
pop0_full = [p[1] for p in pops_full]
pop1_full = [p[2] for p in pops_full]

fig_b     = Figure(size = (1200, 800))
ax_ctrl_b = Axis(fig_b[1, 1];
    xlabel = "Time (arb. units)", ylabel = "Control amplitude",
    title  = "Control Pulses")
xlims!(ax_ctrl_b, plot_times[1], plot_times[end])

ax_pop_b  = Axis(fig_b[2, 1];
    xlabel = "Time (arb. units)", ylabel = "|Uij|^2",
    title  = "Unitary Populations — column 1")
xlims!(ax_pop_b, plot_times[1], plot_times[end])
ylims!(ax_pop_b, 0, 1)
hlines!(ax_pop_b, [0.0, 1.0]; color = :gray, linestyle = :dash, linewidth = 1)

ln_ux_b = lines!(ax_ctrl_b, Float64[], Float64[]; label = "ux (I)", linewidth = 2, color = :steelblue)
ln_uy_b = lines!(ax_ctrl_b, Float64[], Float64[]; label = "uy (Q)", linewidth = 2, color = :firebrick)
axislegend(ax_ctrl_b; position = :rt)

ln_p0_b = lines!(ax_pop_b, Float64[], Float64[]; label = "|U_11|^2", linewidth = 2, color = :mediumseagreen)
ln_p1_b = lines!(ax_pop_b, Float64[], Float64[]; label = "|U_21|^2", linewidth = 2, color = :darkorange)
axislegend(ax_pop_b; position = :rt)

# Callback pattern mirrors reference demo update_pulse_and_state
function update_b!(k)
    ts = plot_times[1:k]
    ln_ux_b[1][] = ts
    ln_ux_b[2][] = controls[1, 1:k]
    ln_uy_b[1][] = ts
    ln_uy_b[2][] = controls[2, 1:k]
    ln_p0_b[1][] = ts
    ln_p0_b[2][] = pop0_full[1:k]
    ln_p1_b[1][] = ts
    ln_p1_b[2][] = pop1_full[1:k]
end

# 1:traj.N 
animate_figure(fig_b, 1:traj.N, update_b!;
    fps = 24, mode = :record, filename = out("pulse_with_populations.mp4"))
println("      Saved: $(out("pulse_with_populations.mp4"))")

# ─amplitude sweep via rollout ──
#
# (rollouts API docs):
#   qtraj_new = rollout(qtraj, new_pulse)
#   fid = fidelity(qtraj_new)

println("\n[5/5] Animation C — amplitude sweep (rollout API)...")

alphas      = collect(range(0.4, 1.6; length = 40))
sweep_fids  = Float64[]
sweep_ctrls = Vector{Matrix{Float64}}()

for alpha in alphas
    scaled_pulse = ZeroOrderPulse(alpha .* controls, times)
    new_traj     = rollout(qtraj, scaled_pulse)   #  rollouts API docs
    push!(sweep_fids,  fidelity(new_traj))         #  rollouts API docs
    push!(sweep_ctrls, alpha .* controls)
end

fig_c     = Figure(size = (1200, 900))
Label(fig_c[0, 1], "Amplitude Sweep — rollout(qtraj, pulse)";
    fontsize = 18, font = :bold)

ax_ctrl_c = Axis(fig_c[1, 1];
    xlabel = "Time (arb. units)", ylabel = "Amplitude", title = "Scaled Pulses")
xlims!(ax_ctrl_c, plot_times[1], plot_times[end])

ax_fid_c  = Axis(fig_c[2, 1];
    xlabel = "Scale factor alpha", ylabel = "Gate fidelity",
    title = "Fidelity vs. Amplitude")
ylims!(ax_fid_c, 0, 1)
lines!(ax_fid_c, alphas, sweep_fids; linewidth = 2, color = :black)

ln_ux_c  = lines!(ax_ctrl_c, plot_times, sweep_ctrls[1][1, :];
    label = "ux", linewidth = 2, color = :steelblue)
ln_uy_c  = lines!(ax_ctrl_c, plot_times, sweep_ctrls[1][2, :];
    label = "uy", linewidth = 2, color = :firebrick)
axislegend(ax_ctrl_c; position = :rt)

dot_c  = scatter!(ax_fid_c, [alphas[1]], [sweep_fids[1]]; markersize = 18, color = :red)
lbl_c  = Label(fig_c[3, 1],
    "alpha=$(round(alphas[1]; digits=2))  fidelity=$(round(sweep_fids[1]; digits=4))";
    fontsize = 15)

function update_c!(k)
    ln_ux_c[2][]  = sweep_ctrls[k][1, :]
    ln_uy_c[2][]  = sweep_ctrls[k][2, :]
    dot_c[1][]    = [alphas[k]]
    dot_c[2][]    = [sweep_fids[k]]
    lbl_c.text[]  = "alpha=$(round(alphas[k]; digits=2))" *
                    "  fidelity=$(round(sweep_fids[k]; digits=4))"
end

animate_figure(fig_c, 1:length(alphas), update_c!;
    fps = 12, mode = :record, filename = out("amplitude_sweep.mp4"))
println("      Saved: $(out("amplitude_sweep.mp4"))")

# ─chirped pulse reveal 

println("\nchirped pulse reveal...")

omega0    = 1.0
omega1    = 3.0
omega_t   = @. omega0 + (omega1 - omega0) * times / T
phi       = cumsum(omega_t) .* (T / N)
chirp_amp = 0.5
chirp_ux  = chirp_amp .* cos.(phi)
chirp_uy  = chirp_amp .* sin.(phi)

fig_d    = Figure(size = (1200, 600))
ax_chirp = Axis(fig_d[1, 1];
    xlabel = "Time (arb. units)", ylabel = "Amplitude",
    title  = "Frequency-Swept (Chirped) Pulse")
xlims!(ax_chirp, times[1], times[end])
ylims!(ax_chirp, -chirp_amp * 1.2, chirp_amp * 1.2)

freq_lbl = Label(fig_d[0, 1],
    "omega(t) = $(round(omega0; digits=2)) rad/arb.unit";
    fontsize = 18, font = :bold)

ln_cux = lines!(ax_chirp, Float64[], Float64[]; label = "ux (I)", linewidth = 2, color = :steelblue)
ln_cuy = lines!(ax_chirp, Float64[], Float64[]; label = "uy (Q)", linewidth = 2, color = :firebrick)
axislegend(ax_chirp; position = :rt)

function update_d!(k)
    ln_cux[1][] = times[1:k]
    ln_cux[2][] = chirp_ux[1:k]
    ln_cuy[1][] = times[1:k]
    ln_cuy[2][] = chirp_uy[1:k]
    omega_now   = omega0 + (omega1 - omega0) * times[k] / T
    freq_lbl.text[] = "omega(t)=$(round(omega_now; digits=2))" *
                      "   t=$(round(times[k]; digits=1))"
end

animate_figure(fig_d, 1:traj.N, update_d!;    # 1:traj.N 
    fps = 24, mode = :record, filename = out("chirp_pulse.mp4"))
println("      Saved: $(out("chirp_pulse.mp4"))")

println("""

============================================================
  All outputs: $(OUTPUT_DIR)/

  pulse_static.png           plot_unitary_populations (standalone)
  pulse_progressive.mp4      animate_name(traj, :u)
  pulse_with_populations.mp4 animate_figure + iso_vec_to_operator
  amplitude_sweep.mp4        rollout(qtraj, pulse) sweep
  chirp_pulse.mp4            frequency-swept pulse reveal

  Piccolo APIs 
    plot_unitary_populations(traj)         
    animate_name(traj, :u; ...)            
    animate_figure(fig, 1:traj.N, fn; ...) 
    rollout(qtraj, pulse)                  rollouts API docs
    fidelity(qtraj)                        rollouts API docs
    traj[:Ũ⃗]                              
    iso_vec_to_operator(Ũ⃗[:,i])[:, 1]    
    get_trajectory(qcp)                    
    get_timesteps(traj)                    
============================================================
""")
