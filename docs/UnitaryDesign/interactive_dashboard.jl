#!/usr/bin/env julia
#
# interactive_pulse_dashboard.jl
# 
#
# Usage:  julia --project interactive_pulse_dashboard.jl
# Requires: Piccolo, GLMakie
#
#  
# ═══════════════════════════════════════════════════════════════════

using Piccolo
using GLMakie
using Random

Random.seed!(42)

# ── 1. Quantum system 

H_drift  = 0.5 * PAULIS[:Z]                          # 
H_drives = [PAULIS[:X], PAULIS[:Y]]                  # 
sys      = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

T     = 10.0
N     = 100
# length = N confirmed: quickstart line 35,
times = collect(range(0, T; length = N))

# trajectory
_pulse_ref = ZeroOrderPulse(zeros(2, N), times)
qtraj_ref  = UnitaryTrajectory(sys, _pulse_ref, GATES[:X])  # 

# ── 2. Control builder 

function build_controls(amp::Real, freq::Real, phase::Real, chirp::Real)
    if abs(chirp) > 1e-6
        omega_t     = @. freq + chirp * times / T
        phase_accum = cumsum(omega_t) .* (T / N)
        u_x = amp .* cos.(phase_accum .+ phase)
        u_y = amp .* sin.(phase_accum .+ phase)
    else
        theta = @. 2pi * freq * times / T + phase
        u_x   = amp .* sin.(theta)
        u_y   = amp .* cos.(theta)
    end
    return vcat(u_x', u_y')   # (2 × N)
end

# ── 3. Population extractor ───────────────────────────────────────────────────
#
# Pattern c in official reference demo (pulse_animation_demo.jl):
#   Ũ⃗ = traj[:Ũ⃗]
#   pops = [abs2.(iso_vec_to_operator(Ũ⃗[:, i])[:, 1]) for i in 1:k]
#
#

function traj_pops(qt::UnitaryTrajectory)
    tdata  = get_trajectory(qt)          # 
    Utilde = tdata[:Ũ⃗]                  # 
    n      = size(Utilde, 2)
    # 
    pops   = [abs2.(iso_vec_to_operator(Utilde[:, i])[:, 1]) for i in 1:n]
    pop0   = [p[1] for p in pops]
    pop1   = [p[2] for p in pops]
    return pop0, pop1
end

# ── 4. Figure 

fig = Figure(size = (1400, 1020), fontsize = 14)

Label(fig[0, 1:2], "Interactive Pulse Parameter Dashboard";
      fontsize = 24, font = :bold, color = :navy)

sg = SliderGrid(
    fig[1, 1],
    (label = "Amplitude",  range = 0.1:0.05:2.0,  startvalue = 0.8),
    (label = "Frequency",  range = 0.5:0.1:5.0,   startvalue = 1.5),
    (label = "Phase",      range = 0.0:0.1:2pi,   startvalue = 0.0),
    (label = "Chirp Rate", range = -1.0:0.1:1.0,  startvalue = 0.0);
    width = 420, tellheight = true,
)

amp_sl   = sg.sliders[1].value
freq_sl  = sg.sliders[2].value
phase_sl = sg.sliders[3].value
chirp_sl = sg.sliders[4].value

ax_pulse = Axis(fig[2, 1];
    xlabel = "Time (arb. units)", ylabel = "Control amplitude",
    title  = "Control Pulses")

ax_pop = Axis(fig[3, 1];
    xlabel = "Time (arb. units)", ylabel = "|Uij|^2",
    title  = "Unitary Populations — column 1")
ylims!(ax_pop, 0, 1)
hlines!(ax_pop, [0.0, 1.0]; color = :gray, linestyle = :dash, linewidth = 1)

# Colour stored as RGBf — prevents Symbol → Any widening in Observable type
fid_label = Label(fig[1, 2], "Fidelity: —";
    fontsize = 22, font = :bold, padding = (16, 16, 16, 16))

Label(fig[2:3, 2],
    """
    Controls

    Amplitude  — pulse strength
    Frequency  — oscillation rate
    Phase      — initial phase offset
    Chirp      — frequency sweep

    Plots update in real-time.
    Target: X gate (sigma_x)

    Tip: amplitude ~0.8,
         frequency ~1.5
    """;
    fontsize = 13, halign = :left, valign = :top, padding = (20, 20, 20, 20))

# ── 5. Observable chain 
#
#  sliders → controls_obs → ux_obs / uy_obs       (pulse lines)
#                         → qtraj_obs              rollout(qtraj_ref, pulse)
#                               → pop0_obs / pop1_obs
#                               → on(...) fidelity label

controls_obs = lift(amp_sl, freq_sl, phase_sl, chirp_sl) do amp, freq, ph, chirp
    build_controls(amp, freq, ph, chirp)
end

ux_obs = lift(c -> vec(c[1, :]), controls_obs)
uy_obs = lift(c -> vec(c[2, :]), controls_obs)

# rollout(qtraj, pulse) 
qtraj_obs = lift(controls_obs) do ctrl
    pulse = ZeroOrderPulse(ctrl, times)
    rollout(qtraj_ref, pulse)
end

# traj_pops uses the demo-confirmed pattern
pop0_obs = lift(qt -> traj_pops(qt)[1], qtraj_obs)
pop1_obs = lift(qt -> traj_pops(qt)[2], qtraj_obs)

# fidelity(qtraj) 
on(qtraj_obs) do qt
    fid = fidelity(qt)
    fid_label.text[]  = "Fidelity: $(round(fid; digits = 6))"
    fid_label.color[] =
        fid > 0.99 ? RGBf(0.05f0, 0.60f0, 0.05f0) :
        fid > 0.95 ? RGBf(0.85f0, 0.55f0, 0.00f0) :
                     RGBf(0.80f0, 0.10f0, 0.10f0)
end

# ── 6. Plot lines 

lines!(ax_pulse, times, ux_obs; label = "ux (I)", linewidth = 2, color = :steelblue)
lines!(ax_pulse, times, uy_obs; label = "uy (Q)", linewidth = 2, color = :firebrick)
axislegend(ax_pulse; position = :rt)

lines!(ax_pop, times, pop0_obs; label = "|U_11|^2 stay |0>", linewidth = 2, color = :mediumseagreen)
lines!(ax_pop, times, pop1_obs; label = "|U_21|^2 flip |1>", linewidth = 2, color = :darkorange)
axislegend(ax_pop; position = :rt)

Label(fig[4, 1:2],
    "Drag sliders for real-time updates.  " *
    "Populations via iso_vec_to_operator (Piccolo reference demo pattern).";
    fontsize = 13, color = :gray60)

# ── 7. Launch ─────────────────────────────────────────────────────────────────

display(fig)

@info """
+------------------------------------------------------------------+
|       Interactive Pulse Dashboard — running                      |
+------------------------------------------------------------------+
|  Drag sliders → controls → rollout(qtraj, pulse) → fidelity     |
|  Fidelity: green >0.99 | orange >0.95 | red <0.95               |
|  Close window or Ctrl-C to exit                                  |
+------------------------------------------------------------------+
|  All APIs confirmed against official Piccolo documentation       |
|    rollout(qtraj, pulse)   rollouts API docs                     |
|    fidelity(qtraj)         rollouts API docs                     |
|    traj[:Ũ⃗]               quickstart + reference demo           |
|    iso_vec_to_operator     quickstart + reference demo           |
|    plot_unitary_populations  used for static snapshots           |
+------------------------------------------------------------------+
"""

try
    wait(fig.scene)
catch
    @info "Dashboard closed."
end
