#!/usr/bin/env julia
#
# Pulse Designer GUI
# Issue #59: Animation of a Pulse Evolving with Parameters
#
# A full GUI application built on GLMakie that adds on top of the two existing
# scripts (pulse_animation_demo.jl, interactive_pulse_dashboard.jl) with:
#
#   • Gate selector drop-down  (X, Y, H, S, T)
#   • Real-time parametric pulse editor  (4 sliders)
#   • Fidelity readout with colour coding
#   • Four tabbed views:
#       Tab 1 — Pulse & Populations (same as dashboard, plus optimised overlay)
#       Tab 2 — IQ Plane  (uₓ vs uᵧ phasor trajectory)
#       Tab 3 — Bloch Sphere projections  (⟨X⟩, ⟨Y⟩, ⟨Z⟩ vs time)
#       Tab 4 — Optimisation convergence plot
#   • [Optimize!] button — runs GRAPE from current pulse as warm start
#   • [Export MP4] button — saves progressive-reveal animation via CairoMakie
#   • [Save PNG]   button — snapshots the live window
#
# Usage:
#   julia pulse_gui.jl
#
# Dependencies:
#   using Pkg; Pkg.add(["Piccolo", "GLMakie", "CairoMakie"])
#
# Layout sketch:
#
#   row -1 ║  [Tab 1]  [Tab 2]  [Tab 3]  [Tab 4]          (tab bar)
#   ───────╫──────────────────────────────────────────────────────────
#   row  0 ║  Piccolo Pulse Designer                        (title)
#   ───────╫────────────────────────╥─────────────────────────────────
#          ║  Gate: [X ▼]           ║                                │
#          ║  ── Sliders ──         ║   (active tab content)         │
#   row  1 ║  Fidelity: 99.xx %     ║                                │
#          ║  [▶ Optimize]          ║                                │
#          ║  [🎬 Export MP4]       ║                                │
#          ║  [📷 Save PNG]         ║                                │
#          ║  status...             ║                                │
#   ───────╫────────────────────────╨─────────────────────────────────

using Piccolo
using GLMakie
using CairoMakie          # needed for MP4 export — imported at top level
using Random
using LinearAlgebra

# Restore GLMakie as the active backend after CairoMakie import
GLMakie.activate!()

Random.seed!(42)

# ═══════════════════════════════════════════════════════════════════════════════
# 1.  Quantum system
# ═══════════════════════════════════════════════════════════════════════════════

const H_DRIFT  = 0.5 * PAULIS[:Z]
const H_DRIVES = [PAULIS[:X], PAULIS[:Y]]
const SYS      = QuantumSystem(H_DRIFT, H_DRIVES, [1.0, 1.0])

const T_GATE = 10.0
const N_GATE = 100
const TIMES  = collect(range(0, T_GATE; length = N_GATE))

const GATE_LABELS = ["X  (NOT)", "Y", "H  (Hadamard)", "S  (Phase)", "T"]
const GATE_MAP    = Dict(
    "X  (NOT)"      => GATES[:X],
    "Y"             => GATES[:Y],
    "H  (Hadamard)" => GATES[:H],
    "S  (Phase)"    => GATES[:S],
    "T"             => GATES[:T],
)

# ═══════════════════════════════════════════════════════════════════════════════
# 2.  Pure helper functions  (no Makie, no globals)
# ═══════════════════════════════════════════════════════════════════════════════

"""Build a 2×N control matrix from the four pulse parameters."""
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

"""Propagate controls and return the unitary trajectory."""
function propagate(ctrl, target)
    UnitaryTrajectory(SYS, ZeroOrderPulse(ctrl, TIMES), target)
end

"""
Extract (p00, p10) = (|U₁₁|², |U₂₁|²) from a trajectory.
Uses traj.data[:Ũ⃗], the iso-vectorised unitary stored by Piccolo.
"""
function populations(traj; col = 1)
    Ũ⃗ = traj.data[:Ũ⃗]
    n  = size(Ũ⃗, 2)
    p00 = Vector{Float64}(undef, n)
    p10 = Vector{Float64}(undef, n)
    for i in 1:n
        U      = iso_vec_to_operator(Ũ⃗[:, i])
        p00[i] = abs2(U[1, col])
        p10[i] = abs2(U[2, col])
    end
    return p00, p10
end

"""
Compute Bloch-sphere components ⟨X⟩, ⟨Y⟩, ⟨Z⟩ of the state U(t)|col⟩.
For a 2×2 unitary column α=U[1,col], β=U[2,col]:
  ⟨X⟩ = 2Re(ᾱβ),  ⟨Y⟩ = 2Im(ᾱβ),  ⟨Z⟩ = |α|²-|β|²
"""
function bloch(traj; col = 1)
    Ũ⃗ = traj.data[:Ũ⃗]
    n  = size(Ũ⃗, 2)
    bx = Vector{Float64}(undef, n)
    by = Vector{Float64}(undef, n)
    bz = Vector{Float64}(undef, n)
    for i in 1:n
        U      = iso_vec_to_operator(Ũ⃗[:, i])
        α, β   = U[1, col], U[2, col]
        bx[i]  = 2real(conj(α) * β)
        by[i]  = 2imag(conj(α) * β)
        bz[i]  = abs2(α) - abs2(β)
    end
    return bx, by, bz
end

fid_color(f) = f > 0.99 ? :green : f > 0.95 ? :darkorange : :red

# ═══════════════════════════════════════════════════════════════════════════════
# 3.  GUI state (mutable, lives outside Observable graph)
# ═══════════════════════════════════════════════════════════════════════════════

mutable struct GUIState
    optimizing    :: Bool
    opt_run       :: Int               # counts completed optimisation runs
    opt_fidelities:: Vector{Float64}
    opt_runs      :: Vector{Int}       # x-axis: run index
end

state = GUIState(false, 0, Float64[], Int[])

# ═══════════════════════════════════════════════════════════════════════════════
# 4.  Figure skeleton
# ═══════════════════════════════════════════════════════════════════════════════

fig = Figure(size = (1650, 950), fontsize = 13)

# Row -1: tab bar (spans right column only)
# Row  0: title  (spans both columns)
# Row  1: left controls | right tab content
Label(fig[0, 1:2], "Piccolo Pulse Designer";
      fontsize = 26, font = :bold, color = :navy)

# ═══════════════════════════════════════════════════════════════════════════════
# 5.  Left control panel
# ═══════════════════════════════════════════════════════════════════════════════

lp = fig[1, 1] = GridLayout()

Label(lp[1, 1], "Target Gate"; font = :bold, halign = :left)
gate_menu = Menu(lp[2, 1]; options = GATE_LABELS, default = GATE_LABELS[1], width = 210)

Label(lp[3, 1], "Pulse Parameters"; font = :bold, halign = :left)
sg = SliderGrid(
    lp[4, 1],
    (label = "Amplitude",    range = 0.05:0.05:2.0, startvalue = 0.8),
    (label = "Frequency",    range = 0.1:0.1:5.0,   startvalue = 1.5),
    (label = "Phase  (rad)", range = 0.0:0.1:2π,    startvalue = 0.0),
    (label = "Chirp rate",   range = -1.0:0.1:1.0,  startvalue = 0.0),
    (label = "Timesteps  N", range = 20:10:200,      startvalue = N_GATE),
    width = 310, tellheight = true,
)

amp_obs   = sg.sliders[1].value
freq_obs  = sg.sliders[2].value
phase_obs = sg.sliders[3].value
chirp_obs = sg.sliders[4].value
# Slider 5 (N) is read imperatively inside the optimize callback, not reactively.

Label(lp[5, 1], "Parametric Fidelity"; font = :bold, halign = :left)
fid_lbl = Label(lp[6, 1], "—";
    fontsize = 22, font = :bold, halign = :center, padding = (8, 8, 4, 4))

Label(lp[7, 1], "Actions"; font = :bold, halign = :left)
btn_opt    = Button(lp[8,  1]; label = "▶  Optimize Pulse", width = 230)
btn_mp4    = Button(lp[9,  1]; label = "🎬  Export MP4",     width = 230)
btn_png    = Button(lp[10, 1]; label = "📷  Save PNG",        width = 230)

status_lbl = Label(lp[11, 1], "Ready.";
    fontsize = 11, color = :gray40, halign = :left, padding = (2, 2, 6, 6))

rowgap!(lp, 5)

# ═══════════════════════════════════════════════════════════════════════════════
# 6.  Tab bar  (row -1, right column)
# ═══════════════════════════════════════════════════════════════════════════════

TAB_NAMES  = ["Pulse & Populations", "IQ Plane", "Bloch Sphere", "Optimisation"]
n_tabs     = length(TAB_NAMES)
tab_grid   = fig[-1, 2] = GridLayout()
tab_btns   = Vector{Button}(undef, n_tabs)
active_tab = Observable(1)

for (i, name) in enumerate(TAB_NAMES)
    b = Button(tab_grid[1, i]; label = name,
               buttoncolor = i == 1 ? :steelblue : :lightgray,
               labelcolor  = i == 1 ? :white     : :black)
    tab_btns[i] = b
    on(b.clicks) do _
        active_tab[] = i
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# 7.  Right panel — four tab content GridLayouts in the same cell [1, 2]
#
#     GLMakie has no native visibility toggle on GridLayouts, so we simulate
#     tabs by assigning all four to fig[1, 2] and then hiding the inactive
#     ones by setting their row/colsizes to Fixed(0).  Each tab's axes are
#     placed at [1,1] inside their own sub-GridLayout so the Fixed(0) trick
#     works cleanly on exactly one row and one column.
# ═══════════════════════════════════════════════════════════════════════════════

# ─── Tab 1: Pulse & Populations ───────────────────────────────────────────────
tab1 = fig[1, 2] = GridLayout()

ax_pulse = Axis(tab1[1, 1];
    title  = "Control Pulses  (parametric)",
    xlabel = "Time (arb. units)", ylabel = "Amplitude")
hlines!(ax_pulse, [0.0]; color = (:black, 0.12), linewidth = 1)

ax_pop = Axis(tab1[2, 1];
    title  = "Unitary Populations  |column 1|²",
    xlabel = "Time (arb. units)", ylabel = "Population")
ylims!(ax_pop, -0.05, 1.05)
hlines!(ax_pop, [0.0, 1.0]; color = (:gray, 0.35), linestyle = :dash)

line_ux  = lines!(ax_pulse, TIMES, zeros(N_GATE); label = "uₓ (I)", linewidth = 2, color = :royalblue)
line_uy  = lines!(ax_pulse, TIMES, zeros(N_GATE); label = "uᵧ (Q)", linewidth = 2, color = :crimson)
axislegend(ax_pulse; position = :rt)

line_p00     = lines!(ax_pop, TIMES, zeros(N_GATE);
    label = "|U₁₁|²  parametric", linewidth = 2, color = :seagreen)
line_p10     = lines!(ax_pop, TIMES, zeros(N_GATE);
    label = "|U₂₁|²  parametric", linewidth = 2, color = :darkorange)
# Optimised overlay — starts full of NaN so nothing renders until opt runs.
line_opt_p00 = lines!(ax_pop, TIMES, fill(NaN, N_GATE);
    label = "|U₁₁|²  optimised", linewidth = 2, color = :seagreen,   linestyle = :dash)
line_opt_p10 = lines!(ax_pop, TIMES, fill(NaN, N_GATE);
    label = "|U₂₁|²  optimised", linewidth = 2, color = :darkorange, linestyle = :dash)
axislegend(ax_pop; position = :rt)

# ─── Tab 2: IQ Plane ──────────────────────────────────────────────────────────
tab2 = fig[1, 2] = GridLayout()

ax_iq = Axis(tab2[1, 1];
    title  = "IQ Plane  (complex control phasor trajectory)",
    xlabel = "In-phase  uₓ", ylabel = "Quadrature  uᵧ",
    aspect = DataAspect())

θ_c = range(0, 2π; length = 200)
lines!(ax_iq, cos.(θ_c), sin.(θ_c);
    color = (:gray, 0.25), linewidth = 1, linestyle = :dash)   # unit-circle guide

iq_traj  = lines!(ax_iq, zeros(N_GATE), zeros(N_GATE);
    linewidth = 2, color = :royalblue, label = "phasor path")
iq_start = scatter!(ax_iq, [0.0], [0.0]; markersize = 14, color = :green, marker = :circle, label = "t=0")
iq_end   = scatter!(ax_iq, [0.0], [0.0]; markersize = 14, color = :red,   marker = :rect,   label = "t=T")
axislegend(ax_iq; position = :rt)

# ─── Tab 3: Bloch Sphere ──────────────────────────────────────────────────────
tab3 = fig[1, 2] = GridLayout()

ax_bloch = Axis(tab3[1, 1];
    title  = "Bloch Sphere Projection  (U(t)|0⟩)",
    xlabel = "Time (arb. units)", ylabel = "Expectation value")
ylims!(ax_bloch, -1.15, 1.15)
hlines!(ax_bloch, [-1.0, 0.0, 1.0]; color = (:gray, 0.3), linestyle = :dash)

bl_x = lines!(ax_bloch, TIMES, zeros(N_GATE); label = "⟨X⟩", linewidth = 2, color = :royalblue)
bl_y = lines!(ax_bloch, TIMES, zeros(N_GATE); label = "⟨Y⟩", linewidth = 2, color = :crimson)
bl_z = lines!(ax_bloch, TIMES, zeros(N_GATE); label = "⟨Z⟩", linewidth = 2, color = :seagreen)
axislegend(ax_bloch; position = :rt)

# ─── Tab 4: Optimisation ──────────────────────────────────────────────────────
tab4 = fig[1, 2] = GridLayout()

ax_conv = Axis(tab4[1, 1];
    title  = "Optimisation Convergence  (fidelity per run)",
    xlabel = "Run index", ylabel = "Gate Fidelity")
ylims!(ax_conv, 0.0, 1.05)
hlines!(ax_conv, [0.99]; color = (:green, 0.5), linestyle = :dash, linewidth = 1)

# Both x and y are Float64 to avoid type confusion when updating.
conv_line    = lines!(ax_conv, Float64[], Float64[]; linewidth = 2, color = :steelblue)
conv_scatter = scatter!(ax_conv, Float64[], Float64[]; markersize = 10, color = :steelblue)
opt_info_lbl = Label(tab4[2, 1],
    "Press [Optimize!] on the left to run GRAPE and see fidelity here.";
    fontsize = 13, color = :gray40, halign = :left)

# ── Tab switching — hide inactive tabs via Fixed(0) row/col sizes ─────────────

ALL_TABS = [tab1, tab2, tab3, tab4]

# Each tab has exactly 1 or 2 rows and 1 column.
TAB_NROWS = [2, 1, 1, 2]   # tab1 has ax_pulse + ax_pop; others have 1 axis + optional label

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
    # Update tab button colours
    for (j, b) in enumerate(tab_btns)
        b.buttoncolor[] = j == t ? :steelblue : :lightgray
        b.labelcolor[]  = j == t ? :white     : :black
    end
end

# Show tab 1, hide 2-4 on startup
notify(active_tab)

# ═══════════════════════════════════════════════════════════════════════════════
# 8.  Reactive Observable graph
# ═══════════════════════════════════════════════════════════════════════════════

current_target = Observable{Matrix{ComplexF64}}(GATE_MAP[GATE_LABELS[1]])

on(gate_menu.selection) do sel
    isnothing(sel) && return
    current_target[] = GATE_MAP[sel]
end

ctrl_obs = lift(amp_obs, freq_obs, phase_obs, chirp_obs) do amp, freq, ph, chirp
    make_controls(amp, freq, ph, chirp)
end

traj_obs = lift(ctrl_obs, current_target) do ctrl, target
    propagate(ctrl, target)
end

fid_obs  = lift(unitary_fidelity, traj_obs)
pop_obs  = lift(traj -> populations(traj), traj_obs)
blch_obs = lift(traj -> bloch(traj),       traj_obs)

# ── Tab 1 updates ─────────────────────────────────────────────────────────────

on(ctrl_obs) do ctrl
    line_ux[2][] = ctrl[1, :]
    line_uy[2][] = ctrl[2, :]
end

on(pop_obs) do (p00, p10)
    line_p00[2][] = p00
    line_p10[2][] = p10
end

on(fid_obs) do fid
    fid_lbl.text[]  = "$(round(fid * 100; digits = 3)) %"
    fid_lbl.color[] = fid_color(fid)
end

# ── Tab 2 updates ─────────────────────────────────────────────────────────────

on(ctrl_obs) do ctrl
    ux, uy = ctrl[1, :], ctrl[2, :]
    iq_traj[1][]  = ux
    iq_traj[2][]  = uy
    iq_start[1][] = [ux[1]];   iq_start[2][] = [uy[1]]
    iq_end[1][]   = [ux[end]]; iq_end[2][]   = [uy[end]]
end

# ── Tab 3 updates ─────────────────────────────────────────────────────────────

on(blch_obs) do (bx, by, bz)
    bl_x[2][] = bx
    bl_y[2][] = by
    bl_z[2][] = bz
end

# ═══════════════════════════════════════════════════════════════════════════════
# 9.  Optimize button callback
# ═══════════════════════════════════════════════════════════════════════════════

on(btn_opt.clicks) do _
    state.optimizing && return          # prevent re-entrant clicks

    state.optimizing   = true
    status_lbl.text[]  = "⏳  Optimising…  (up to ~30 s)"
    status_lbl.color[] = :darkorange

    init_ctrl = ctrl_obs[]             # warm start from current slider values
    target    = current_target[]
    N_opt     = sg.sliders[5].value[]  # from the N slider

    try
        pl    = ZeroOrderPulse(init_ctrl, TIMES)
        qt    = UnitaryTrajectory(SYS, pl, target)
        qcp   = SmoothPulseProblem(qt, N_opt; Q = 100.0, R = 1e-2, ddu_bound = 1.0)

        solve!(qcp; max_iter = 150, verbose = false, print_level = 1)

        opt_traj = get_trajectory(qcp)
        opt_fid  = fidelity(qcp)

        # ── Convergence plot (Tab 4) ───────────────────────────────────────────
        # Piccolo's solver doesn't expose a per-iteration callback by default.
        # We record the terminal fidelity per button press.  If your build of
        # Piccolo stores a history in qcp (e.g. qcp.data[:fidelity_history]),
        # unpack that vector here instead.
        state.opt_run += 1
        push!(state.opt_fidelities, opt_fid)
        push!(state.opt_runs,       Float64(state.opt_run))

        conv_line[1][]    = state.opt_runs
        conv_line[2][]    = state.opt_fidelities
        conv_scatter[1][] = state.opt_runs
        conv_scatter[2][] = state.opt_fidelities
        opt_info_lbl.text[] = "Run $(state.opt_run)  →  fidelity = $(round(opt_fid * 100; digits = 4)) %"

        # ── Optimised population overlay (Tab 1) ──────────────────────────────
        opt_p00, opt_p10 = populations(opt_traj)
        n_opt = length(opt_p00)

        if n_opt == N_GATE
            # Same grid as parametric — update y only
            line_opt_p00[2][] = opt_p00
            line_opt_p10[2][] = opt_p10
        else
            # Different N — update both x and y
            ts_opt = collect(range(0, T_GATE; length = n_opt))
            line_opt_p00[1][] = ts_opt;  line_opt_p00[2][] = opt_p00
            line_opt_p10[1][] = ts_opt;  line_opt_p10[2][] = opt_p10
        end

        status_lbl.text[]  = "✓  fidelity = $(round(opt_fid * 100; digits = 4)) %   (run $(state.opt_run))"
        status_lbl.color[] = fid_color(opt_fid)

    catch err
        status_lbl.text[]  = "✗  $(sprint(showerror, err))"
        status_lbl.color[] = :red
        @error "Optimisation failed" exception = (err, catch_backtrace())
    finally
        state.optimizing = false
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# 10.  Export MP4 button callback
# ═══════════════════════════════════════════════════════════════════════════════

on(btn_mp4.clicks) do _
    status_lbl.text[]  = "🎬  Exporting MP4…"
    status_lbl.color[] = :steelblue

    try
        outdir  = mkpath("pulse_gui_exports")
        stamp   = string(round(Int, time()))
        outfile = joinpath(outdir, "pulse_$(stamp).mp4")

        # Snapshot current state before spawning the animation loop
        ctrl_snap = ctrl_obs[]
        pop_snap  = pop_obs[]           # (p00, p10) tuple

        # Build a CairoMakie figure for reliable headless recording
        fig_exp = CairoMakie.Figure(size = (900, 600), fontsize = 13)
        gate_name = something(gate_menu.selection[], GATE_LABELS[1])
        Label(fig_exp[0, 1], "Pulse Export — $(gate_name)";
              fontsize = 16, font = :bold)

        ax_e1 = CairoMakie.Axis(fig_exp[1, 1];
            title = "Control Pulses", xlabel = "Time", ylabel = "Amplitude")
        ax_e2 = CairoMakie.Axis(fig_exp[2, 1];
            title = "Populations",   xlabel = "Time", ylabel = "|U|²")
        CairoMakie.ylims!(ax_e2, 0, 1)
        CairoMakie.xlims!(ax_e1, TIMES[1], TIMES[end])
        CairoMakie.xlims!(ax_e2, TIMES[1], TIMES[end])

        le_ux  = CairoMakie.lines!(ax_e1, Float64[], Float64[];
            label = "uₓ", color = :royalblue, linewidth = 2)
        le_uy  = CairoMakie.lines!(ax_e1, Float64[], Float64[];
            label = "uᵧ", color = :crimson,   linewidth = 2)
        le_p00 = CairoMakie.lines!(ax_e2, Float64[], Float64[];
            label = "|U₁₁|²", color = :seagreen,   linewidth = 2)
        le_p10 = CairoMakie.lines!(ax_e2, Float64[], Float64[];
            label = "|U₂₁|²", color = :darkorange, linewidth = 2)
        CairoMakie.axislegend(ax_e1; position = :rt)
        CairoMakie.axislegend(ax_e2; position = :rt)

        function update_export!(k)
            ts_k = TIMES[1:k]
            le_ux[1][]  = ts_k;  le_ux[2][]  = ctrl_snap[1, 1:k]
            le_uy[1][]  = ts_k;  le_uy[2][]  = ctrl_snap[2, 1:k]
            le_p00[1][] = ts_k;  le_p00[2][] = pop_snap[1][1:k]
            le_p10[1][] = ts_k;  le_p10[2][] = pop_snap[2][1:k]
        end

        animate_figure(fig_exp, 1:N_GATE, update_export!;
            fps = 24, mode = :record, filename = outfile)

        status_lbl.text[]  = "✓  Saved → $(outfile)"
        status_lbl.color[] = :seagreen

    catch err
        status_lbl.text[]  = "✗  Export failed: $(sprint(showerror, err))"
        status_lbl.color[] = :red
        @error "MP4 export failed" exception = (err, catch_backtrace())
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# 11.  Save PNG button callback
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
# 12.  Column / row sizing and initial data push
# ═══════════════════════════════════════════════════════════════════════════════

colsize!(fig.layout, 1, Fixed(350))   # left control panel
colsize!(fig.layout, 2, Relative(0.73))

rowsize!(fig.layout, -1, Fixed(36))   # tab bar
rowsize!(fig.layout, 0,  Fixed(44))   # title
rowsize!(fig.layout, 1,  Auto())      # main content

# Populate all plots before window appears
notify(amp_obs)

# ═══════════════════════════════════════════════════════════════════════════════
# 13.  Launch
# ═══════════════════════════════════════════════════════════════════════════════

display(fig)

println()
println("=" ^ 72)
println("Piccolo Pulse Designer  —  running")
println()
println("  Left panel")
println("    Gate drop-down  : select target unitary (X, Y, H, S, T)")
println("    Sliders         : reshape the parametric pulse in real time")
println("    Fidelity        : colour-coded gate fidelity (green ≥ 0.99)")
println("    [Optimize!]     : run GRAPE from the current pulse as warm start")
println("    [Export MP4]    : save a progressive-reveal animation (CairoMakie)")
println("    [Save PNG]      : snapshot the current window")
println()
println("  Right tabs")
println("    Pulse & Populations  — control signals and |U|² evolution")
println("    IQ Plane             — complex phasor trajectory in the (uₓ, uᵧ) plane")
println("    Bloch Sphere         — ⟨X⟩, ⟨Y⟩, ⟨Z⟩ of U(t)|0⟩ vs time")
println("    Optimisation         — fidelity per optimization run")
println()
println("  Exports → ./pulse_gui_exports/")
println("  Close the window or Ctrl-C to quit.")
println("=" ^ 72)
println()

try
    wait(fig.scene)
catch _
    println("GUI closed.")
end
