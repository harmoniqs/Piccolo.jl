# # Pulse Animation Guide
#
# This guide walks through animating quantum control pulses in Piccolo.jl —
# from a static snapshot all the way to an interactive real-time dashboard.
#
# ## Why animate pulses?
#
# Pulse visualisation is central to quantum optimal control.  Watching a
# pulse draw itself alongside the evolving quantum populations gives physical
# intuition that a static plot cannot provide.
#
# ## Prerequisites
#
# ```julia
# using Pkg
# Pkg.add(["Piccolo", "CairoMakie"])   # headless / file output
# # or
# Pkg.add(["Piccolo", "GLMakie"])      # interactive windows
# ```
#

# ## Step 1 — Setup

using Piccolo
using CairoMakie   # swap for GLMakie for interactive windows
using Random

Random.seed!(42)

# ## Step 2 — Define the system and optimise

H_drift  = 0.5 * PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
sys      = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

T      = 10.0
N      = 100
# length = N: 
times  = collect(range(0, T; length = N))

pulse  = ZeroOrderPulse(0.1 * randn(2, N), times)
qtraj  = UnitaryTrajectory(sys, pulse, GATES[:X])
qcp    = SmoothPulseProblem(qtraj, N; Q = 100.0, R = 1e-2, ddu_bound = 1.0)
solve!(qcp; max_iter = 50, verbose = false, print_level = 1)

traj = get_trajectory(qcp)   # 

## 
plot_times = cumsum([0.0; get_timesteps(traj)])[1:end-1]
controls   = traj[:u]

@info "Fidelity" fidelity(qcp)

# ## Step 3 — Static snapshot: `plot_unitary_populations`
#
# `plot_unitary_populations(traj)` returns a standalone `Figure` containing
# all |U_ij|² population curves.  Confirmed in quickstart:
#
# ```julia
# fig = plot_unitary_populations(traj)          # 
# fig_mintime = plot_unitary_populations(...)   # 
# ```
#
# **Important:** because it returns a `Figure`, it cannot be embedded as a
# sub-axis inside another figure.  For co-animated layouts the population
# data must be extracted via `traj[:Ũ⃗]` — see Step 5.

fig_static = plot_unitary_populations(traj)
## save("pulse_static.png", fig_static)

# ## Step 4 — Progressive reveal: `animate_name`
#
# `animate_name` draws successive prefixes of a named trajectory variable.
# 
#
# ```julia
# fig = animate_name(traj, :u; fps=24, mode=:record, filename=...)
# ```

animate_name(traj, :u; fps = 24, mode = :record,
    filename = "pulse_progressive.mp4")

# ## Step 5 — Co-animate pulse and populations: `animate_figure`
#
# For a layout that places control pulses and unitary populations side by
# side we need per-axis population data.
#
# ```julia
# Ũ⃗ = traj[:Ũ⃗]                                          # 
# pops = [abs2.(iso_vec_to_operator(Ũ⃗[:, i])[:, 1])    # 
#          for i in 1:k]
# ```
#
# 
#
# The frame range `1:traj.N` is also confirmed in the reference demo
# (line 308: `animate_figure(fig2, 1:traj.N, update_pulse_and_state!; ...)`).

Utilde = traj[:Ũ⃗]   # c

## :
pops_all = [abs2.(iso_vec_to_operator(Utilde[:, i])[:, 1]) for i in 1:traj.N]
pop0_all = [p[1] for p in pops_all]
pop1_all = [p[2] for p in pops_all]

fig_anim  = Figure(size = (1200, 800))

ax_ctrl   = Axis(fig_anim[1, 1]; xlabel = "Time", ylabel = "Amplitude",
    title = "Control Pulses")
xlims!(ax_ctrl, plot_times[1], plot_times[end])

ax_pop    = Axis(fig_anim[2, 1]; xlabel = "Time", ylabel = "|Uij|^2",
    title = "Unitary Populations — column 1")
xlims!(ax_pop, plot_times[1], plot_times[end])
ylims!(ax_pop, 0, 1)
hlines!(ax_pop, [0.0, 1.0]; color = :gray, linestyle = :dash, linewidth = 1)

ln_ux = lines!(ax_ctrl, Float64[], Float64[]; label = "ux (I)", linewidth = 2, color = :steelblue)
ln_uy = lines!(ax_ctrl, Float64[], Float64[]; label = "uy (Q)", linewidth = 2, color = :firebrick)
axislegend(ax_ctrl; position = :rt)

ln_p0 = lines!(ax_pop, Float64[], Float64[]; label = "|U_11|^2", linewidth = 2, color = :mediumseagreen)
ln_p1 = lines!(ax_pop, Float64[], Float64[]; label = "|U_21|^2", linewidth = 2, color = :darkorange)
axislegend(ax_pop; position = :rt)

## Callback structure mirrors reference demo update_pulse_and_state
function update_pulse_and_pop!(k)
    ts      = plot_times[1:k]
    ln_ux[1][] = ts;  ln_ux[2][] = controls[1, 1:k]
    ln_uy[1][] = ts;  ln_uy[2][] = controls[2, 1:k]
    ln_p0[1][] = ts;  ln_p0[2][] = pop0_all[1:k]
    ln_p1[1][] = ts;  ln_p1[2][] = pop1_all[1:k]
end

## 1:traj.N confirmed: 
animate_figure(fig_anim, 1:traj.N, update_pulse_and_pop!;
    fps = 24, mode = :record, filename = "pulse_with_populations.mp4")

# ## Step 6 — Parameter sweep: `rollout(qtraj, pulse)`
#
#
#
# ```julia
# qtraj_new = rollout(qtraj, new_pulse)   # rollouts API docs
# fid = fidelity(qtraj_new)              # rollouts API docs
# ```
#
# This re-uses the system and goal from `qtraj` so the caller never has to
# reconstruct them, and delegates all ODE propagation to Piccolo internals.

alphas = collect(range(0.4, 1.6; length = 40))
fids   = Float64[]
ctrls_list = Vector{Matrix{Float64}}()

for alpha in alphas
    scaled_pulse = ZeroOrderPulse(alpha .* controls, times)
    new_traj     = rollout(qtraj, scaled_pulse)   # rollouts API docs
    push!(fids,       fidelity(new_traj))          # rollouts API docs
    push!(ctrls_list, alpha .* controls)
end

fig_sweep = Figure(size = (1200, 900))
Label(fig_sweep[0, 1], "Amplitude Sweep — rollout(qtraj, pulse)";
    fontsize = 18, font = :bold)

ax_ctrl_s = Axis(fig_sweep[1, 1]; xlabel = "Time", ylabel = "Amplitude",
    title = "Scaled Pulses")
xlims!(ax_ctrl_s, plot_times[1], plot_times[end])

ax_fid_s  = Axis(fig_sweep[2, 1]; xlabel = "alpha", ylabel = "Fidelity",
    title = "Gate Fidelity vs. Amplitude")
ylims!(ax_fid_s, 0, 1)
lines!(ax_fid_s, alphas, fids; linewidth = 2, color = :black)

ln_ux_s = lines!(ax_ctrl_s, plot_times, ctrls_list[1][1, :];
    label = "ux", linewidth = 2, color = :steelblue)
ln_uy_s = lines!(ax_ctrl_s, plot_times, ctrls_list[1][2, :];
    label = "uy", linewidth = 2, color = :firebrick)
axislegend(ax_ctrl_s; position = :rt)

dot_s  = scatter!(ax_fid_s, [alphas[1]], [fids[1]]; markersize = 18, color = :red)
lbl_s  = Label(fig_sweep[3, 1],
    "alpha=$(round(alphas[1]; digits=2))  fidelity=$(round(fids[1]; digits=4))";
    fontsize = 15)

function update_sweep!(k)
    ln_ux_s[2][]  = ctrls_list[k][1, :]
    ln_uy_s[2][]  = ctrls_list[k][2, :]
    dot_s[1][]    = [alphas[k]]
    dot_s[2][]    = [fids[k]]
    lbl_s.text[]  = "alpha=$(round(alphas[k]; digits=2))" *
                    "  fidelity=$(round(fids[k]; digits=4))"
end

animate_figure(fig_sweep, 1:length(alphas), update_sweep!;
    fps = 12, mode = :record, filename = "amplitude_sweep.mp4")

# ## Step 7 — Real-time interactive dashboard (GLMakie)
#
# For interactive exploration swap CairoMakie for GLMakie and use
# Makie Observables.  The full implementation lives in
# `interactive_pulse_dashboard.jl`.  Key design decisions:
#
# ### Observable chain
#
# ```
# amp_sl, freq_sl, phase_sl, chirp_sl   GLMakie slider Observables
#         │
#         ▼
#   controls_obs   build_controls(amp, freq, phase, chirp)
#         │
#         ├──▶ ux_obs / uy_obs          pulse plot lines
#         │
#         └──▶ qtraj_obs               rollout(qtraj_ref, pulse)
#                   │
#                   ├──▶ pop0_obs      iso_vec_to_operator → |U_11|²
#                   ├──▶ pop1_obs      iso_vec_to_operator → |U_21|²
#                   └──▶ on(...)       fidelity(qtraj) → update label
# ```
#
