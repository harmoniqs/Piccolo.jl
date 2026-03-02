#!/usr/bin/env julia
#
# Interactive Pulse Parameter Dashboard
# Issue #59: Animation of a Pulse Evolving with Parameters
#
# Opens a live GLMakie window with four sliders.
# All plots update in real time as you drag them.
#
# Usage:
#   julia interactive_pulse_dashboard.jl
#
# Dependencies:
#   using Pkg; Pkg.add(["Piccolo", "GLMakie"])
#
# API used (verified against QuantumCollocation.jl docs, Jan 2026):
#   PAULIS[:X], GATES[:X]
#   UnitaryTrajectory(sys, U_goal, T)
#   SmoothPulseProblem(qtraj, N; Q, R)
#   solve!(qcp; options=IpoptOptions(max_iter=N))
#   fidelity(qcp)
#   get_trajectory(qcp) → NamedTrajectory
#   traj[:u], traj[:Ũ⃗], traj[:Δt]
#   iso_vec_to_operator(v)

using Piccolo
using GLMakie
using Random
using LinearAlgebra

Random.seed!(42)

# ─── Quantum system ───────────────────────────────────────────────────────────

H_drift  = 0.1 * PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
sys      = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

const T_GATE = 10.0
const N_GATE = 51

# Fixed time axis for display (N_GATE equally-spaced points over T_GATE).
# The actual traj[:Δt] may vary slightly after optimisation, but for the
# parametric slider view we always use this uniform grid.
const TIMES = collect(range(0, T_GATE; length=N_GATE))

# ─── Pure physics helpers (no Piccolo calls needed) ───────────────────────────

"""
Build 2×N parametric controls from the four slider values.
Row 1 = uₓ (I channel), Row 2 = uᵧ (Q channel).
"""
function make_controls(amp, freq, ph, chirp)
    if abs(chirp) > 0.01
        ω_t = @. freq + chirp * TIMES / T_GATE
        φ   = cumsum(ω_t) .* (T_GATE / N_GATE) .+ ph
        u_x = @. amp * cos(φ)
        u_y = @. amp * sin(φ)
    else
        φ   = @. 2π * freq * TIMES / T_GATE + ph
        u_x = @. amp * sin(φ)
        u_y = @. amp * cos(φ)
    end
    return [u_x'; u_y']
end

"""
Propagate 2×N controls through the quantum system.
Returns (p00, p10, fidelity) where pij = |U_{i,1}|² at each timestep.

We use a simple matrix-exponential Magnus step (first-order) which is fast
enough for real-time slider updates and avoids a full Ipopt round-trip.
The fidelity is computed using the standard unitary fidelity formula.
"""
function propagate(ctrl::Matrix{Float64}, U_goal)
    N    = size(ctrl, 2)
    dt   = T_GATE / N
    U    = Matrix{ComplexF64}(I, 2, 2)

    p00 = Vector{Float64}(undef, N)
    p10 = Vector{Float64}(undef, N)

    for k in 1:N
        H_k = H_drift + ctrl[1, k] * H_drives[1] + ctrl[2, k] * H_drives[2]
        U   = exp(-im * H_k * dt) * U
        p00[k] = abs2(U[1, 1])
        p10[k] = abs2(U[2, 1])
    end

    d   = size(U, 1)
    fid = abs2(tr(U_goal' * U)) / d^2
    return p00, p10, fid
end

fid_color(f) = f > 0.99 ? :green : f > 0.95 ? :darkorange : :red

# ─── Figure layout ────────────────────────────────────────────────────────────

fig = Figure(size=(1450, 1000), fontsize=14)

Label(fig[0, 1:2], "Interactive Pulse Parameter Dashboard";
      fontsize=24, font=:bold, color=:navy)

sg = SliderGrid(
    fig[1, 1],
    (label="Amplitude",    range=0.1:0.05:2.0, startvalue=0.5),
    (label="Frequency",    range=0.5:0.1:5.0,  startvalue=1.5),
    (label="Phase  (rad)", range=0.0:0.1:2π,   startvalue=0.0),
    (label="Chirp rate",   range=-1.0:0.1:1.0, startvalue=0.0),
    width=420, tellheight=true,
)

amp_obs   = sg.sliders[1].value
freq_obs  = sg.sliders[2].value
phase_obs = sg.sliders[3].value
chirp_obs = sg.sliders[4].value

ax_pulse = Axis(fig[2, 1];
    title="Control Pulses", xlabel="Time (arb. units)", ylabel="Amplitude")
hlines!(ax_pulse, [0.0]; color=(:black, 0.15), linewidth=1)

ax_pop = Axis(fig[3, 1];
    title="Unitary Populations (column 1)",
    xlabel="Time (arb. units)", ylabel="|Uᵢⱼ|²")
ylims!(ax_pop, -0.05, 1.05)
hlines!(ax_pop, [0.0, 1.0]; color=(:gray, 0.4), linestyle=:dash)

fid_lbl = Label(fig[1, 2], "Fidelity\n—";
    fontsize=20, font=:bold, padding=(20, 20, 20, 20),
    halign=:center, valign=:center)

Label(fig[2:3, 2],
    "Slider guide\n\n" *
    "Amplitude   — pulse strength\n" *
    "Frequency   — oscillation rate\n" *
    "Phase       — initial phase offset\n" *
    "Chirp rate  — linear frequency sweep\n\n" *
    "Fidelity colour key\n\n" *
    "  Green  > 0.99   excellent\n" *
    "  Orange > 0.95   good\n" *
    "  Red   ≤ 0.95   needs tuning\n\n" *
    "Target: X gate\n\n" *
    "Note: fidelity computed via\n" *
    "first-order Magnus propagation\n" *
    "(fast; use [Optimize!] in\n" *
    "pulse_gui.jl for exact GRAPE).";
    fontsize=13, halign=:left, valign=:top, padding=(20, 20, 20, 20))

Label(fig[4, 1:2],
    "Drag any slider to update pulse shape, populations, and fidelity in real time.";
    fontsize=15, color=:gray40)

# ─── Reactive graph ───────────────────────────────────────────────────────────

# 1. Controls (pure math — instant)
ctrl_obs = lift(amp_obs, freq_obs, phase_obs, chirp_obs) do amp, freq, ph, chirp
    make_controls(amp, freq, ph, chirp)
end

# 2. Physics result (Magnus propagation — ~1 ms, imperceptible)
result_obs = lift(ctrl_obs) do ctrl
    p00, p10, fid = propagate(ctrl, GATES[:X])
    (p00=p00, p10=p10, fid=fid)
end

fid_obs = lift(r -> r.fid, result_obs)
p00_obs = lift(r -> r.p00, result_obs)
p10_obs = lift(r -> r.p10, result_obs)

# ─── Fidelity label ──────────────────────────────────────────────────────────

on(fid_obs) do fid
    fid_lbl.text[]  = "Fidelity\n$(round(fid; digits=6))"
    fid_lbl.color[] = fid_color(fid)
end

# ─── Pulse plot ───────────────────────────────────────────────────────────────

ux_obs = lift(c -> c[1, :], ctrl_obs)
uy_obs = lift(c -> c[2, :], ctrl_obs)

lines!(ax_pulse, TIMES, ux_obs; label="uₓ (I)", linewidth=2, color=:royalblue)
lines!(ax_pulse, TIMES, uy_obs; label="uᵧ (Q)", linewidth=2, color=:crimson)
axislegend(ax_pulse; position=:rt)

# ─── Population plot ─────────────────────────────────────────────────────────

lines!(ax_pop, TIMES, p00_obs; label="|U₁₁|²  (stay |0⟩)", linewidth=2, color=:seagreen)
lines!(ax_pop, TIMES, p10_obs; label="|U₂₁|²  (flip |1⟩)", linewidth=2, color=:darkorange)
axislegend(ax_pop; position=:rt)

# ─── Launch ───────────────────────────────────────────────────────────────────

display(fig)

println()
println("=" ^ 60)
println("Interactive dashboard running!")
println("  Drag sliders → real-time pulse + population updates.")
println("  Close the window or Ctrl-C to exit.")
println("=" ^ 60)

try
    wait(fig.scene)
catch _
    println("Dashboard closed.")
end
