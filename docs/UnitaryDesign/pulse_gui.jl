#!/usr/bin/env julia
#
# Pulse Designer GUI
# Issue #59: Animation of a Pulse Evolving with Parameters
#
# Full GUI built on GLMakie:
#   • Gate selector drop-down (X, Y, H)
#   • 4 slider parameters — real-time pulse + physics update
#   • Fidelity readout (colour-coded)
#   • 4 tabbed views:
#       Tab 1 — Pulse & Populations (+ optimised overlay after [Optimize!])
#       Tab 2 — IQ Plane (phasor trajectory)
#       Tab 3 — Bloch Sphere ⟨X⟩, ⟨Y⟩, ⟨Z⟩ vs time
#       Tab 4 — Optimisation convergence plot
#   • [Optimize!] — runs GRAPE warm-started from current sliders
#   • [Save PNG]   — snapshots the window
#
# Usage:
#   julia pulse_gui.jl
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

# ═══════════════════════════════════════════════════════════════════════════════
# 1.  Quantum system
# ═══════════════════════════════════════════════════════════════════════════════

H_drift  = 0.1 * PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
sys      = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

const T_GATE = 10.0
const N_GATE = 51
const TIMES  = collect(range(0, T_GATE; length=N_GATE))

# Only gates confirmed to exist in GATES (X, Y, H are standard)
const GATE_LABELS = ["X  (NOT)", "Y", "H  (Hadamard)"]
const GATE_MAP    = Dict(
    "X  (NOT)"      => GATES[:X],
    "Y"             => GATES[:Y],
    "H  (Hadamard)" => GATES[:H],
)

# ═══════════════════════════════════════════════════════════════════════════════
# 2.  Physics helpers — Magnus propagation (fast, no solver overhead)
# ═══════════════════════════════════════════════════════════════════════════════

"""Build 2×N controls from the four slider parameters."""
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
Propagate 2×N controls through the quantum system using first-order Magnus.
Returns (p00, p10, bloch_x, bloch_y, bloch_z, fidelity).
All vectors are length N.
"""
function propagate(ctrl::Matrix{Float64}, U_goal)
    N   = size(ctrl, 2)
    dt  = T_GATE / N
    U   = Matrix{ComplexF64}(I, 2, 2)

    p00 = Vector{Float64}(undef, N)
    p10 = Vector{Float64}(undef, N)
    bx  = Vector{Float64}(undef, N)
    by  = Vector{Float64}(undef, N)
    bz  = Vector{Float64}(undef, N)

    for k in 1:N
        H_k = H_drift + ctrl[1, k] * H_drives[1] + ctrl[2, k] * H_drives[2]
        U   = exp(-im * H_k * dt) * U

        # Column 1 of U is the state |0⟩ → U|0⟩
        α, β   = U[1, 1], U[2, 1]
        p00[k] = abs2(α)
        p10[k] = abs2(β)
        bx[k]  = 2real(conj(α) * β)
        by[k]  = 2imag(conj(α) * β)
        bz[k]  = abs2(α) - abs2(β)
    end

    d   = size(U, 1)
    fid = abs2(tr(U_goal' * U)) / d^2
    return p00, p10, bx, by, bz, fid
end

fid_color(f) = f > 0.99 ? :green : f > 0.95 ? :darkorange : :red

# ═══════════════════════════════════════════════════════════════════════════════
# 3.  GUI state
# ═══════════════════════════════════════════════════════════════════════════════

mutable struct GUIState
    busy         :: Bool
    opt_run      :: Int
    opt_fids     :: Vector{Float64}
    opt_runs     :: Vector{Float64}
end
state = GUIState(false, 0, Float64[], Float64[])

# ═══════════════════════════════════════════════════════════════════════════════
# 4.  Figure skeleton
# ═══════════════════════════════════════════════════════════════════════════════

fig = Figure(size=(1650, 950), fontsize=13)

Label(fig[0, 1:2], "Piccolo Pulse Designer";
      fontsize=26, font=:bold, color=:navy)

# ═══════════════════════════════════════════════════════════════════════════════
# 5.  Left control panel
# ═══════════════════════════════════════════════════════════════════════════════

lp = fig[1, 1] = GridLayout()

Label(lp[1, 1], "Target Gate"; font=:bold, halign=:left)
gate_menu = Menu(lp[2, 1]; options=GATE_LABELS, default=GATE_LABELS[1], width=210)

Label(lp[3, 1], "Pulse Parameters"; font=:bold, halign=:left)
sg = SliderGrid(
    lp[4, 1],
    (label="Amplitude",    range=0.05:0.05:2.0, startvalue=0.5),
    (label="Frequency",    range=0.1:0.1:5.0,   startvalue=1.5),
    (label="Phase  (rad)", range=0.0:0.1:2π,    startvalue=0.0),
    (label="Chirp rate",   range=-1.0:0.1:1.0,  startvalue=0.0),
    width=310, tellheight=true,
)

amp_obs   = sg.sliders[1].value
freq_obs  = sg.sliders[2].value
phase_obs = sg.sliders[3].value
chirp_obs = sg.sliders[4].value

Label(lp[5, 1], "Parametric Fidelity"; font=:bold, halign=:left)
fid_lbl = Label(lp[6, 1], "—";
    fontsize=22, font=:bold, halign=:center, padding=(8, 8, 4, 4))

Label(lp[7, 1], "Actions"; font=:bold, halign=:left)
btn_opt = Button(lp[8, 1]; label="▶  Optimize Pulse", width=230)
btn_png = Button(lp[9, 1]; label="📷  Save PNG",        width=230)

status_lbl = Label(lp[10, 1], "Ready.";
    fontsize=11, color=:gray40, halign=:left, padding=(2, 2, 6, 6))

rowgap!(lp, 5)

# ═══════════════════════════════════════════════════════════════════════════════
# 6.  Tab bar (row -1, right column)
# ═══════════════════════════════════════════════════════════════════════════════

TAB_NAMES  = ["Pulse & Populations", "IQ Plane", "Bloch Sphere", "Optimisation"]
tab_grid   = fig[-1, 2] = GridLayout()
tab_btns   = Vector{Button}(undef, length(TAB_NAMES))
active_tab = Observable(1)

for (i, name) in enumerate(TAB_NAMES)
    b = Button(tab_grid[1, i]; label=name,
               buttoncolor=i == 1 ? :steelblue : :lightgray,
               labelcolor =i == 1 ? :white     : :black)
    tab_btns[i] = b
    on(b.clicks) do _
        active_tab[] = i
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# 7.  Tab content (all four share fig[1, 2])
# ═══════════════════════════════════════════════════════════════════════════════

# ── Tab 1 ──────────────────────────────────────────────────────────────────────
tab1 = fig[1, 2] = GridLayout()

ax_pulse = Axis(tab1[1, 1];
    title="Control Pulses  (parametric)",
    xlabel="Time (arb. units)", ylabel="Amplitude")
hlines!(ax_pulse, [0.0]; color=(:black, 0.12))

ax_pop = Axis(tab1[2, 1];
    title="Unitary Populations  |column 1|²",
    xlabel="Time (arb. units)", ylabel="Population")
ylims!(ax_pop, -0.05, 1.05)
hlines!(ax_pop, [0.0, 1.0]; color=(:gray, 0.35), linestyle=:dash)

l_ux     = lines!(ax_pulse, TIMES, zeros(N_GATE);
    label="uₓ (I)", linewidth=2, color=:royalblue)
l_uy     = lines!(ax_pulse, TIMES, zeros(N_GATE);
    label="uᵧ (Q)", linewidth=2, color=:crimson)
axislegend(ax_pulse; position=:rt)

l_p00     = lines!(ax_pop, TIMES, fill(0.5, N_GATE);
    label="|U₁₁|²  parametric", linewidth=2, color=:seagreen)
l_p10     = lines!(ax_pop, TIMES, fill(0.5, N_GATE);
    label="|U₂₁|²  parametric", linewidth=2, color=:darkorange)
l_opt_p00 = lines!(ax_pop, TIMES, fill(NaN, N_GATE);
    label="|U₁₁|²  optimised", linewidth=2, color=:seagreen,   linestyle=:dash)
l_opt_p10 = lines!(ax_pop, TIMES, fill(NaN, N_GATE);
    label="|U₂₁|²  optimised", linewidth=2, color=:darkorange, linestyle=:dash)
axislegend(ax_pop; position=:rt)

# ── Tab 2 ──────────────────────────────────────────────────────────────────────
tab2 = fig[1, 2] = GridLayout()

ax_iq = Axis(tab2[1, 1];
    title="IQ Plane  (phasor trajectory)",
    xlabel="In-phase uₓ", ylabel="Quadrature uᵧ",
    aspect=DataAspect())

θ_c = range(0, 2π; length=200)
lines!(ax_iq, cos.(θ_c), sin.(θ_c); color=(:gray, 0.25), linewidth=1, linestyle=:dash)

iq_line  = lines!(ax_iq, zeros(N_GATE), zeros(N_GATE);
    linewidth=2, color=:royalblue, label="path")
iq_start = scatter!(ax_iq, [0.0], [0.0]; markersize=14, color=:green,  marker=:circle, label="t=0")
iq_end   = scatter!(ax_iq, [0.0], [0.0]; markersize=14, color=:red,    marker=:rect,   label="t=T")
axislegend(ax_iq; position=:rt)

# ── Tab 3 ──────────────────────────────────────────────────────────────────────
tab3 = fig[1, 2] = GridLayout()

ax_bloch = Axis(tab3[1, 1];
    title="Bloch Sphere Projections  (U(t)|0⟩)",
    xlabel="Time (arb. units)", ylabel="Expectation value")
ylims!(ax_bloch, -1.15, 1.15)
hlines!(ax_bloch, [-1.0, 0.0, 1.0]; color=(:gray, 0.3), linestyle=:dash)

bl_x = lines!(ax_bloch, TIMES, zeros(N_GATE); label="⟨X⟩", linewidth=2, color=:royalblue)
bl_y = lines!(ax_bloch, TIMES, zeros(N_GATE); label="⟨Y⟩", linewidth=2, color=:crimson)
bl_z = lines!(ax_bloch, TIMES, ones(N_GATE);  label="⟨Z⟩", linewidth=2, color=:seagreen)
axislegend(ax_bloch; position=:rt)

# ── Tab 4 ──────────────────────────────────────────────────────────────────────
tab4 = fig[1, 2] = GridLayout()

ax_conv = Axis(tab4[1, 1];
    title="Optimisation Convergence  (fidelity per run)",
    xlabel="Run index", ylabel="Gate Fidelity")
ylims!(ax_conv, 0, 1.05)
hlines!(ax_conv, [0.99]; color=(:green, 0.5), linestyle=:dash, linewidth=1)

conv_line    = lines!(ax_conv, Float64[], Float64[]; linewidth=2, color=:steelblue)
conv_scatter = scatter!(ax_conv, Float64[], Float64[]; markersize=10, color=:steelblue)
opt_info_lbl = Label(tab4[2, 1],
    "Press [Optimize!] to run GRAPE and see convergence.";
    fontsize=13, color=:gray40, halign=:left)

# ── Tab visibility — hide inactive tabs via Fixed(0) sizes ───────────────────

ALL_TABS  = [tab1, tab2, tab3, tab4]
TAB_NROWS = [2, 1, 1, 2]

function set_tab_visible!(tab, nrows, visible)
    sz = visible ? Auto() : Fixed(0)
    colsize!(tab, 1, sz)
    for r in 1:nrows
        rowsize!(tab, r, sz)
    end
end

on(active_tab) do t
    for (i, (tab, nr)) in enumerate(zip(ALL_TABS, TAB_NROWS))
        set_tab_visible!(tab, nr, i == t)
    end
    for (j, b) in enumerate(tab_btns)
        b.buttoncolor[] = j == t ? :steelblue : :lightgray
        b.labelcolor[]  = j == t ? :white     : :black
    end
end

notify(active_tab)   # show tab 1, hide others at startup

# ═══════════════════════════════════════════════════════════════════════════════
# 8.  Reactive graph
# ═══════════════════════════════════════════════════════════════════════════════

current_target = Observable{Matrix{ComplexF64}}(GATE_MAP[GATE_LABELS[1]])

on(gate_menu.selection) do sel
    isnothing(sel) && return
    current_target[] = GATE_MAP[sel]
end

ctrl_obs = lift(amp_obs, freq_obs, phase_obs, chirp_obs) do amp, freq, ph, chirp
    make_controls(amp, freq, ph, chirp)
end

result_obs = lift(ctrl_obs, current_target) do ctrl, target
    p00, p10, bx, by, bz, fid = propagate(ctrl, target)
    (p00=p00, p10=p10, bx=bx, by=by, bz=bz, fid=fid)
end

fid_obs = lift(r -> r.fid,  result_obs)
p00_obs = lift(r -> r.p00,  result_obs)
p10_obs = lift(r -> r.p10,  result_obs)
bx_obs  = lift(r -> r.bx,   result_obs)
by_obs  = lift(r -> r.by,   result_obs)
bz_obs  = lift(r -> r.bz,   result_obs)

# ── Tab 1 wires ───────────────────────────────────────────────────────────────

on(ctrl_obs) do ctrl
    l_ux[2][] = ctrl[1, :]
    l_uy[2][] = ctrl[2, :]
end

on(p00_obs) do p; l_p00[2][] = p end
on(p10_obs) do p; l_p10[2][] = p end

on(fid_obs) do fid
    fid_lbl.text[]  = "$(round(fid * 100; digits=3)) %"
    fid_lbl.color[] = fid_color(fid)
end

# ── Tab 2 wires ───────────────────────────────────────────────────────────────

on(ctrl_obs) do ctrl
    ux, uy = ctrl[1, :], ctrl[2, :]
    iq_line[1][]  = ux;  iq_line[2][]  = uy
    iq_start[1][] = [ux[1]];   iq_start[2][] = [uy[1]]
    iq_end[1][]   = [ux[end]]; iq_end[2][]   = [uy[end]]
end

# ── Tab 3 wires ───────────────────────────────────────────────────────────────

on(bx_obs) do v; bl_x[2][] = v end
on(by_obs) do v; bl_y[2][] = v end
on(bz_obs) do v; bl_z[2][] = v end

# ═══════════════════════════════════════════════════════════════════════════════
# 9.  Optimize button
# ═══════════════════════════════════════════════════════════════════════════════

on(btn_opt.clicks) do _
    state.busy && return

    state.busy         = true
    status_lbl.text[]  = "⏳  Optimising…  (~30 s)"
    status_lbl.color[] = :darkorange

    init_ctrl = ctrl_obs[]
    target    = current_target[]

    try
        # Build problem warm-started from current slider controls.
        # UnitaryTrajectory initialises controls to zero; we overwrite via
        # the NamedTrajectory returned by get_trajectory before solving.
        qt  = UnitaryTrajectory(sys, target, T_GATE)
        qcp = SmoothPulseProblem(qt, N_GATE; Q=100.0, R=1e-2)

        traj_init = get_trajectory(qcp)
        traj_init[:u] .= init_ctrl    # warm-start controls

        solve!(qcp; options=IpoptOptions(max_iter=150))

        opt_traj = get_trajectory(qcp)
        opt_fid  = fidelity(qcp)

        # ── Convergence plot (Tab 4) ──────────────────────────────────────────
        state.opt_run += 1
        push!(state.opt_fids, opt_fid)
        push!(state.opt_runs, Float64(state.opt_run))

        conv_line[1][]    = state.opt_runs
        conv_line[2][]    = state.opt_fids
        conv_scatter[1][] = state.opt_runs
        conv_scatter[2][] = state.opt_fids
        opt_info_lbl.text[] =
            "Run $(state.opt_run)  →  fidelity = $(round(opt_fid * 100; digits=4)) %"

        # ── Optimised population overlay on Tab 1 ────────────────────────────
        # Extract populations from the optimised NamedTrajectory.
        U⃗ = opt_traj[:Ũ⃗]
        T_opt = size(U⃗, 2)
        opt_p00 = [abs2(iso_vec_to_operator(U⃗[:, k])[1, 1]) for k in 1:T_opt]
        opt_p10 = [abs2(iso_vec_to_operator(U⃗[:, k])[2, 1]) for k in 1:T_opt]

        if T_opt == N_GATE
            l_opt_p00[2][] = opt_p00
            l_opt_p10[2][] = opt_p10
        else
            ts_opt = collect(range(0, T_GATE; length=T_opt))
            l_opt_p00[1][] = ts_opt;  l_opt_p00[2][] = opt_p00
            l_opt_p10[1][] = ts_opt;  l_opt_p10[2][] = opt_p10
        end

        status_lbl.text[]  =
            "✓  fidelity = $(round(opt_fid * 100; digits=4)) %  (run $(state.opt_run))"
        status_lbl.color[] = fid_color(opt_fid)

    catch err
        status_lbl.text[]  = "✗  $(sprint(showerror, err))"
        status_lbl.color[] = :red
        @error "Optimisation failed" exception=(err, catch_backtrace())
    finally
        state.busy = false
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# 10.  Save PNG button
# ═══════════════════════════════════════════════════════════════════════════════

on(btn_png.clicks) do _
    try
        outdir  = mkpath("pulse_gui_exports")
        stamp   = string(round(Int, time()))
        outfile = joinpath(outdir, "pulse_$(stamp).png")
        save(outfile, fig)
        status_lbl.text[]  = "✓  PNG → $(outfile)"
        status_lbl.color[] = :seagreen
    catch err
        status_lbl.text[]  = "✗  PNG failed: $(sprint(showerror, err))"
        status_lbl.color[] = :red
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# 11.  Layout sizing and initial data push
# ═══════════════════════════════════════════════════════════════════════════════

colsize!(fig.layout, 1, Fixed(350))
colsize!(fig.layout, 2, Relative(0.73))
rowsize!(fig.layout, -1, Fixed(36))
rowsize!(fig.layout, 0,  Fixed(44))
rowsize!(fig.layout, 1,  Auto())

notify(amp_obs)   # populate all plots before window appears

# ═══════════════════════════════════════════════════════════════════════════════
# 12.  Launch
# ═══════════════════════════════════════════════════════════════════════════════

display(fig)

println()
println("=" ^ 65)
println("Piccolo Pulse Designer — running")
println()
println("  Gate drop-down : X, Y, H")
println("  Sliders        : reshape the parametric pulse in real time")
println("  Fidelity       : colour-coded (green ≥ 0.99, orange ≥ 0.95)")
println("  [Optimize!]    : run GRAPE warm-started from current sliders")
println("                   (optimised populations overlaid on Tab 1)")
println("  [Save PNG]     : snapshot the current window")
println()
println("  Tabs: Pulse & Populations | IQ Plane | Bloch Sphere | Optimisation")
println("  Exports → ./pulse_gui_exports/")
println("  Close the window or Ctrl-C to quit.")
println("=" ^ 65)

try
    wait(fig.scene)
catch _
    println("GUI closed.")
end
