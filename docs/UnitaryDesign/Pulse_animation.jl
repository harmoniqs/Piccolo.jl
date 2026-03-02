#!/usr/bin/env julia
#
# Pulse Parameter Animation Demo
# Issue #59: Animation of a Pulse Evolving with Parameters
#
# Generates five .mp4 files:
#   pulse_progressive.mp4       — optimised pulse drawing itself step-by-step
#   pulse_with_populations.mp4  — pulse + unitary population evolution
#   amplitude_sweep.mp4         — fidelity vs amplitude scaling
#   chirp_pulse.mp4             — frequency-swept (chirped) IQ pulse
#   pulse_sequence.mp4          — X → Y → X composite gate sequence
#
# Usage:
#   julia pulse_animation_demo.jl                  # output → ./pulse_animations/
#   julia pulse_animation_demo.jl my_output_dir
#
# Dependencies:
#   using Pkg; Pkg.add(["Piccolo", "CairoMakie"])
#
# API used (verified against QuantumCollocation.jl docs, Jan 2026):
#   PAULIS[:X]  / PAULIS.X  — both work
#   GATES[:X]   / GATES.X   — both work
#   UnitaryTrajectory(sys, U_goal, T)
#   SmoothPulseProblem(qtraj, N; Q, R)
#   solve!(qcp; options=IpoptOptions(max_iter=N))
#   fidelity(qcp)
#   get_trajectory(qcp)          → NamedTrajectory
#   traj[:u]                     → 2×N control matrix
#   traj[:Ũ⃗]                    → iso-vec unitary, columns = timesteps
#   traj[:Δt]                    → length-N timestep vector
#   iso_vec_to_operator(v)       → 2×2 complex matrix
#   operator_to_iso_vec(U)       → iso-vec column

using Piccolo
using CairoMakie
using Random

Random.seed!(42)

# ─── Output dir ───────────────────────────────────────────────────────────────

const OUTPUT_DIR = length(ARGS) > 0 ? ARGS[1] : "pulse_animations"
mkpath(OUTPUT_DIR)
outpath(name) = joinpath(OUTPUT_DIR, name)

println("=== Pulse Parameter Animation Demo ===")
println("Output → $(OUTPUT_DIR)/\n")

# ─── Quantum system ───────────────────────────────────────────────────────────

H_drift  = 0.1 * PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
sys      = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

const T_GATE = 10.0   # total gate duration
const N_GATE = 51     # collocation timesteps

# ─── Helpers ──────────────────────────────────────────────────────────────────

"""
Build a cumulative time axis from a solved trajectory's Δt vector.
Returns a length-N Float64 vector starting at 0.
"""
function time_axis(traj)
    Δt = traj[:Δt]
    cumsum([0.0; Δt[1:end-1]])   # length = traj.T
end

"""
Extract (p00, p10) = (|U_{1,col}|², |U_{2,col}|²) from a NamedTrajectory.
Uses traj[:Ũ⃗] — the iso-vectorised unitary stored by Piccolo.
"""
function extract_pops(traj; col=1)
    Umat = traj[:Ũ⃗]              # dims: iso_vec_length × T
    T    = size(Umat, 2)
    p00  = Vector{Float64}(undef, T)
    p10  = Vector{Float64}(undef, T)
    for k in 1:T
        U       = iso_vec_to_operator(Umat[:, k])
        p00[k]  = abs2(U[1, col])
        p10[k]  = abs2(U[2, col])
    end
    return p00, p10
end

"""Apply standard styling to a control-pulse Axis."""
function style_ctrl!(ax)
    ax.xlabel = "Time (arb. units)"
    ax.ylabel = "Control amplitude"
    hlines!(ax, [0.0]; color=(:black, 0.15), linewidth=1)
end

"""Apply standard styling to a population Axis."""
function style_pop!(ax)
    ax.xlabel = "Time (arb. units)"
    ax.ylabel = "Population |Uᵢⱼ|²"
    ylims!(ax, -0.05, 1.05)
    hlines!(ax, [0.0, 1.0]; color=(:gray, 0.4), linestyle=:dash)
end

# ─── Solve reference X-gate ───────────────────────────────────────────────────

println("Optimising reference X-gate…")

qtraj = UnitaryTrajectory(sys, GATES[:X], T_GATE)
qcp   = SmoothPulseProblem(qtraj, N_GATE; Q=100.0, R=1e-2)
solve!(qcp; options=IpoptOptions(max_iter=100))

traj    = get_trajectory(qcp)
ts      = time_axis(traj)
opt_fid = fidelity(qcp)

println("  fidelity = $(round(opt_fid; digits=4))\n")

# ─── Animation 1: Progressive pulse reveal ────────────────────────────────────

println("[1/5] Progressive pulse reveal…")

let
    u   = traj[:u]    # 2×T
    T   = size(u, 2)
    fps = 24

    fig = Figure(size=(900, 400), fontsize=13)
    Label(fig[0, 1], "Progressive Pulse Reveal"; fontsize=16, font=:bold)
    ax  = Axis(fig[1, 1])
    style_ctrl!(ax)
    xlims!(ax, ts[1], ts[end])

    l_ux = lines!(ax, Float64[], Float64[]; label="uₓ (I)", linewidth=2, color=:royalblue)
    l_uy = lines!(ax, Float64[], Float64[]; label="uᵧ (Q)", linewidth=2, color=:crimson)
    axislegend(ax; position=:rt)

    record(fig, outpath("pulse_progressive.mp4"); framerate=fps) do io
        for k in 1:T
            l_ux[1][] = ts[1:k];  l_ux[2][] = u[1, 1:k]
            l_uy[1][] = ts[1:k];  l_uy[2][] = u[2, 1:k]
            recordframe!(io)
        end
    end
end

println("      ✓  pulse_progressive.mp4")

# ─── Animation 2: Pulse + populations ────────────────────────────────────────

println("[2/5] Pulse + populations…")

let
    u         = traj[:u]
    p00, p10  = extract_pops(traj)
    T         = size(u, 2)
    fps       = 24

    fig = Figure(size=(1100, 750), fontsize=13)
    Label(fig[0, 1], "Control Pulse & Population Evolution"; fontsize=18, font=:bold)

    ax_c = Axis(fig[1, 1]; title="Control Pulses")
    ax_p = Axis(fig[2, 1]; title="Unitary Populations (column 1)")
    style_ctrl!(ax_c);  style_pop!(ax_p)
    xlims!(ax_c, ts[1], ts[end]);  xlims!(ax_p, ts[1], ts[end])

    l_ux  = lines!(ax_c, Float64[], Float64[]; label="uₓ (I)", linewidth=2, color=:royalblue)
    l_uy  = lines!(ax_c, Float64[], Float64[]; label="uᵧ (Q)", linewidth=2, color=:crimson)
    l_p00 = lines!(ax_p, Float64[], Float64[]; label="|U₁₁|²", linewidth=2, color=:seagreen)
    l_p10 = lines!(ax_p, Float64[], Float64[]; label="|U₂₁|²", linewidth=2, color=:darkorange)
    axislegend(ax_c; position=:rt);  axislegend(ax_p; position=:rt)

    record(fig, outpath("pulse_with_populations.mp4"); framerate=fps) do io
        for k in 1:T
            t_k = ts[1:k]
            l_ux[1][]  = t_k;  l_ux[2][]  = u[1, 1:k]
            l_uy[1][]  = t_k;  l_uy[2][]  = u[2, 1:k]
            l_p00[1][] = t_k;  l_p00[2][] = p00[1:k]
            l_p10[1][] = t_k;  l_p10[2][] = p10[1:k]
            recordframe!(io)
        end
    end
end

println("      ✓  pulse_with_populations.mp4")

# ─── Animation 3: Amplitude sweep ────────────────────────────────────────────

println("[3/5] Amplitude sweep…")

let
    base_u  = traj[:u]           # 2×T from the optimised trajectory
    n_amp   = 30
    amps    = range(0.3, 1.8; length=n_amp)

    amp_fids  = Vector{Float64}(undef, n_amp)
    amp_ctrls = Vector{Matrix{Float64}}(undef, n_amp)

    # For each scaling: build a fresh UnitaryTrajectory from scaled controls
    # by constructing a NamedTrajectory manually (the only public way to hold
    # arbitrary controls in the current API without an optimizer round-trip).
    for (i, a) in enumerate(amps)
        sc   = a .* base_u
        Δt_v = traj[:Δt]
        n    = size(sc, 2)

        # Propagate: compute unitaries using Piccolo's iso formalism
        # The standard way is operator_to_iso_vec + Pade propagator, but the
        # simplest public-API approach is to create a NamedTrajectory and use
        # unitary_fidelity(traj, U_goal) if available, or propagate manually.
        # We propagate with a simple first-order Magnus step since we only need
        # the fidelity approximation for plotting (exact fidelity would require
        # a full solver round-trip per frame which is too slow for 30 frames).
        U = Matrix{ComplexF64}(I, 2, 2)
        for k in 1:n
            H_k = H_drift + sc[1, k] * H_drives[1] + sc[2, k] * H_drives[2]
            U   = exp(-im * H_k * Δt_v[k]) * U
        end
        goal = GATES[:X]
        d    = size(U, 1)
        amp_fids[i]  = abs2(tr(goal' * U)) / d^2   # unitary fidelity formula
        amp_ctrls[i] = sc
    end

    fig = Figure(size=(1100, 850), fontsize=13)
    Label(fig[0, 1], "Amplitude Sweep"; fontsize=18, font=:bold)
    ax_p = Axis(fig[1, 1]; title="Control Pulses")
    ax_f = Axis(fig[2, 1];
        title="Gate Fidelity vs Amplitude",
        xlabel="Amplitude scaling factor", ylabel="Fidelity")
    style_ctrl!(ax_p)
    xlims!(ax_p, ts[1], ts[end]);  ylims!(ax_f, 0, 1.05)

    lines!(ax_f, collect(amps), amp_fids; linewidth=2, color=:black)
    hlines!(ax_f, [0.99]; color=(:green, 0.5), linestyle=:dash, linewidth=1)

    l_ux3 = lines!(ax_p, ts, amp_ctrls[1][1, :]; label="uₓ", linewidth=2, color=:royalblue)
    l_uy3 = lines!(ax_p, ts, amp_ctrls[1][2, :]; label="uᵧ", linewidth=2, color=:crimson)
    s_dot = scatter!(ax_f, [amps[1]], [amp_fids[1]]; markersize=16, color=:red)
    lbl3  = Label(fig[1, 1, Top()], ""; fontsize=14, halign=:right, padding=(0, 10, 5, 0))
    axislegend(ax_p; position=:rt)

    record(fig, outpath("amplitude_sweep.mp4"); framerate=10) do io
        for k in 1:n_amp
            l_ux3[2][]  = amp_ctrls[k][1, :]
            l_uy3[2][]  = amp_ctrls[k][2, :]
            s_dot[1][]  = [amps[k]]
            s_dot[2][]  = [amp_fids[k]]
            lbl3.text[] = "scale=$(round(amps[k];digits=2))  fid=$(round(amp_fids[k];digits=4))"
            recordframe!(io)
        end
    end
end

println("      ✓  amplitude_sweep.mp4")

# ─── Animation 4: Chirped pulse ───────────────────────────────────────────────

println("[4/5] Chirped pulse…")

let
    N_c  = 100
    T_c  = T_GATE
    ts_c = collect(range(0, T_c; length=N_c))
    Δt_c = T_c / N_c

    ω_start, ω_end = 0.5, 4.0
    ω_t = @. ω_start + (ω_end - ω_start) * ts_c / T_c
    φ   = cumsum(ω_t .* Δt_c)
    u_x = 0.5 .* cos.(φ)
    u_y = 0.5 .* sin.(φ)

    ω_inst(t) = ω_start + (ω_end - ω_start) * t / T_c

    fig = Figure(size=(1100, 500), fontsize=13)
    Label(fig[0, 1], "Frequency-Swept (Chirped) Pulse"; fontsize=18, font=:bold)
    ax  = Axis(fig[1, 1])
    style_ctrl!(ax)
    xlims!(ax, ts_c[1], ts_c[end])

    l_cx = lines!(ax, Float64[], Float64[]; label="uₓ (I)", linewidth=2, color=:royalblue)
    l_cy = lines!(ax, Float64[], Float64[]; label="uᵧ (Q)", linewidth=2, color=:crimson)
    lbl4 = Label(fig[1, 1, Top()], ""; fontsize=14, halign=:right, padding=(0, 10, 5, 0))
    axislegend(ax; position=:rt)

    record(fig, outpath("chirp_pulse.mp4"); framerate=24) do io
        for k in 1:N_c
            t_k = ts_c[1:k]
            l_cx[1][] = t_k;  l_cx[2][] = u_x[1:k]
            l_cy[1][] = t_k;  l_cy[2][] = u_y[1:k]
            ω = ω_inst(ts_c[k])
            lbl4.text[] = "t=$(round(ts_c[k];digits=1))  ω(t)=$(round(ω;digits=2)) rad/s"
            recordframe!(io)
        end
    end
end

println("      ✓  chirp_pulse.mp4")

# ─── Animation 5: Multi-gate sequence (X → Y → X echo) ──────────────────────

println("[5/5] Gate sequence X → Y → X…")

let
    targets    = [GATES[:X], GATES[:Y], GATES[:X]]
    gate_names = ["X gate", "Y gate", "X echo"]
    gate_cols  = [:royalblue, :seagreen, :darkorange]
    T_seq      = 5.0
    N_seq      = 31

    seq_trajs = map(targets) do target
        qt   = UnitaryTrajectory(sys, target, T_seq)
        prob = SmoothPulseProblem(qt, N_seq; Q=100.0, R=1e-2)
        solve!(prob; options=IpoptOptions(max_iter=100))
        get_trajectory(prob)
    end

    seq_ts    = [time_axis(t) for t in seq_trajs]
    seq_xlims = [ts_g[end] for ts_g in seq_ts]

    # Build flat frame list
    gate_idx  = Int[]
    local_idx = Int[]
    for (g, tr) in enumerate(seq_trajs), k in 1:tr.T
        push!(gate_idx, g)
        push!(local_idx, k)
    end

    fig = Figure(size=(1100, 500), fontsize=13)
    Label(fig[0, 1], "Multi-Gate Pulse Sequence"; fontsize=18, font=:bold)
    ax  = Axis(fig[1, 1])
    style_ctrl!(ax)

    # Observable for reactive line colour
    col_obs = Observable{Symbol}(:royalblue)
    l_sx = lines!(ax, Float64[], Float64[]; label="uₓ", linewidth=2, color=col_obs)
    l_sy = lines!(ax, Float64[], Float64[]; label="uᵧ", linewidth=2, color=:crimson)
    lbl5 = Label(fig[1, 1, Top()], ""; fontsize=15, font=:bold,
                 halign=:left, padding=(10, 0, 5, 0))
    axislegend(ax; position=:rt)

    record(fig, outpath("pulse_sequence.mp4"); framerate=24) do io
        for frame in eachindex(gate_idx)
            g  = gate_idx[frame]
            k  = local_idx[frame]
            tr = seq_trajs[g]
            ts_g = seq_ts[g]

            l_sx[1][] = ts_g[1:k];  l_sx[2][] = tr[:u][1, 1:k]
            l_sy[1][] = ts_g[1:k];  l_sy[2][] = tr[:u][2, 1:k]
            col_obs[]      = gate_cols[g]
            lbl5.text[]    = "$(gate_names[g])  ($(g)/$(length(targets)))"
            lbl5.color[]   = gate_cols[g]
            xlims!(ax, 0, seq_xlims[g])
            recordframe!(io)
        end
    end
end

println("      ✓  pulse_sequence.mp4")

# ─── Summary ──────────────────────────────────────────────────────────────────

println()
println("=" ^ 60)
println("✓ All animations saved to: $(OUTPUT_DIR)/")
println("  Convert to GIF: ffmpeg -i file.mp4 file.gif")
println("=" ^ 60)
