#!/usr/bin/env julia
#
# Interactive Pulse Parameter Dashboard
# Issue #59: Animation of a Pulse Evolving with Parameters
#
# Opens a live GLMakie window with four sliders (amplitude, frequency,
# phase, chirp rate). All plots update in real-time as you drag them.
#
# Usage:
#   julia interactive_pulse_dashboard.jl
#
# Dependencies: Piccolo, GLMakie
#   using Pkg; Pkg.add(["Piccolo", "GLMakie"])

using Piccolo
using GLMakie
using Random

Random.seed!(42)

# ─── Quantum System ───────────────────────────────────────────────────────────

H_drift  = 0.5 * PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
sys      = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

const T_GATE = 10.0
const N_GATE = 100
const TIMES  = collect(range(0, T_GATE; length = N_GATE))

# ─── Layout ───────────────────────────────────────────────────────────────────

fig = Figure(size = (1450, 1000), fontsize = 14)

Label(fig[0, 1:2], "Interactive Pulse Parameter Dashboard";
      fontsize = 24, font = :bold, color = :navy)

# Slider panel (left column, row 1)
sg = SliderGrid(
    fig[1, 1],
    (label = "Amplitude",   range = 0.1:0.05:2.0, startvalue = 0.8),
    (label = "Frequency",   range = 0.5:0.1:5.0,  startvalue = 1.5),
    (label = "Phase (rad)", range = 0.0:0.1:2π,   startvalue = 0.0),
    (label = "Chirp rate",  range = -1.0:0.1:1.0, startvalue = 0.0),
    width = 420,
    tellheight = true,
)

amp_obs   = sg.sliders[1].value
freq_obs  = sg.sliders[2].value
phase_obs = sg.sliders[3].value
chirp_obs = sg.sliders[4].value

# Pulse axes (left column, rows 2-3)
ax_pulse = Axis(fig[2, 1];
    title  = "Control Pulses",
    xlabel = "Time (arb. units)",
    ylabel = "Amplitude")

ax_pop = Axis(fig[3, 1];
    title  = "Unitary Populations (column 1)",
    xlabel = "Time (arb. units)",
    ylabel = "|Uij|^2")

ylims!(ax_pop, 0, 1)
hlines!(ax_pop, [0.0, 1.0]; color = (:gray, 0.4), linestyle = :dash, linewidth = 1)

# Right panel: fidelity readout + help text
fidelity_lbl = Label(fig[1, 2], "Fidelity\nCalculating...";
    fontsize = 20, font = :bold,
    padding  = (20, 20, 20, 20),
    halign   = :center, valign = :center)

Label(fig[2:3, 2],
    "Slider guide\n\n" *
    "Amplitude  - pulse strength\n" *
    "Frequency  - oscillation rate\n" *
    "Phase      - initial phase offset\n" *
    "Chirp rate - linear frequency sweep\n\n" *
    "Fidelity colour key\n\n" *
    "Green  > 0.99   excellent\n" *
    "Orange > 0.95   good\n" *
    "Red   <= 0.95   needs tuning\n\n" *
    "Target gate: X\n\n" *
    "Tip: amp=0.8, freq=1.5\n" *
    "gives high fidelity.";
    fontsize = 13,
    halign   = :left,
    valign   = :top,
    padding  = (20, 20, 20, 20),
)

Label(fig[4, 1:2],
    "Drag any slider to see real-time updates of pulse shape, populations, and gate fidelity.";
    fontsize = 15, color = :gray40)

# ─── Reactive Observables ─────────────────────────────────────────────────────

# 1. Controls matrix (2 x N)
controls_obs = lift(amp_obs, freq_obs, phase_obs, chirp_obs) do amp, freq, ph, chirp
    if abs(chirp) > 0.01
        omega_t = @. freq + chirp * TIMES / T_GATE
        # Integrate ω(t) to get phase; include ph as initial phase offset.
        phi = cumsum(omega_t) .* (T_GATE / N_GATE) .+ ph
        u_x = @. amp * cos(phi)
        u_y = @. amp * sin(phi)
    else
        phi = @. 2pi * freq * TIMES / T_GATE + ph
        u_x = @. amp * sin(phi)
        u_y = @. amp * cos(phi)
    end
    [u_x'; u_y']
end

# 2. Trajectory
traj_obs = lift(controls_obs) do ctrl
    UnitaryTrajectory(sys, ZeroOrderPulse(ctrl, TIMES), GATES[:X])
end

# 3. Fidelity scalar
fidelity_obs = lift(unitary_fidelity, traj_obs)

# 4. Populations — field name is :Ũ⃗ (iso-vectorised unitary), consistent with
#    Piccolo's NamedTrajectory convention and the animation demo helper.
pop_obs = lift(traj_obs) do traj
    Ũ⃗ = traj.data[:Ũ⃗]
    n  = size(Ũ⃗, 2)
    p00 = Vector{Float64}(undef, n)
    p10 = Vector{Float64}(undef, n)
    for i in 1:n
        U      = iso_vec_to_operator(Ũ⃗[:, i])
        p00[i] = abs2(U[1, 1])
        p10[i] = abs2(U[2, 1])
    end
    (p00, p10)
end

# ─── Fidelity label update ────────────────────────────────────────────────────

on(fidelity_obs) do fid
    fidelity_lbl.text[]  = "Fidelity\n$(round(fid; digits = 6))"
    fidelity_lbl.color[] = fid > 0.99 ? :green :
                           fid > 0.95 ? :darkorange : :red
end

# ─── Pulse Plot ───────────────────────────────────────────────────────────────

ux_obs = lift(c -> c[1, :], controls_obs)
uy_obs = lift(c -> c[2, :], controls_obs)

lines!(ax_pulse, TIMES, ux_obs; label = "ux (I)", linewidth = 2, color = :royalblue)
lines!(ax_pulse, TIMES, uy_obs; label = "uy (Q)", linewidth = 2, color = :crimson)
hlines!(ax_pulse, [0.0]; color = (:black, 0.15), linewidth = 1)
axislegend(ax_pulse; position = :rt)

# ─── Population Plot ──────────────────────────────────────────────────────────

p00_obs = lift(p -> p[1], pop_obs)
p10_obs = lift(p -> p[2], pop_obs)

lines!(ax_pop, TIMES, p00_obs; label = "|U11|^2 (stay |0>)", linewidth = 2, color = :seagreen)
lines!(ax_pop, TIMES, p10_obs; label = "|U21|^2 (flip |1>)", linewidth = 2, color = :darkorange)
axislegend(ax_pop; position = :rt)

# ─── Launch ───────────────────────────────────────────────────────────────────

display(fig)

println()
println("=" ^ 68)
println("Interactive dashboard is running!")
println()
println("  Drag any slider to update pulse shape and fidelity in real-time.")
println("  Close the window (or press Ctrl+C here) to exit.")
println("=" ^ 68)
println()

try
    wait(fig.scene)
catch _
    println("Dashboard closed.")
end
