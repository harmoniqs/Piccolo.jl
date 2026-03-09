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
# 

# ## Step 1 — Setup

using Piccolo
using CairoMakie   # 
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
# all |U_ij|² population curves.  
#
# ```julia
# fig = plot_unitary_populations(traj)          
# fig_mintime = plot_unitary_populations(...)   
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
# side we need per-axis population data.  The correct way to extract it is
# 
#
# ```julia
# Ũ⃗ = traj[:Ũ⃗]                                          
# pops = [abs2.(iso_vec_to_operator(Ũ⃗[:, i])[:, 1])    
#          for i in 1:k]
# ```
#
# 

Utilde = traj[:Ũ⃗]   # 

## 
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

# Explicit Observables — the only pattern that is stable across all Makie
# versions and backends.  Makie constructs Point2f internally from the two
# Observables; .
ux_x = Observable(Float64[]);  ux_y = Observable(Float64[])
uy_x = Observable(Float64[]);  uy_y = Observable(Float64[])
lines!(ax_ctrl, ux_x, ux_y; label = "ux (I)", linewidth = 2, color = :steelblue)
lines!(ax_ctrl, uy_x, uy_y; label = "uy (Q)", linewidth = 2, color = :firebrick)
axislegend(ax_ctrl; position = :rt)

p0_x = Observable(Float64[]);  p0_y = Observable(Float64[])
p1_x = Observable(Float64[]);  p1_y = Observable(Float64[])
lines!(ax_pop, p0_x, p0_y; label = "|U_11|^2", linewidth = 2, color = :mediumseagreen)
lines!(ax_pop, p1_x, p1_y; label = "|U_21|^2", linewidth = 2, color = :darkorange)
axislegend(ax_pop; position = :rt)

## Callback — update Observables directly, never index into the line object
function update_pulse_and_pop!(k)
    ts      = plot_times[1:k]
    ux_x[]  = ts;  ux_y[]  = controls[1, 1:k]
    uy_x[]  = ts;  uy_y[]  = controls[2, 1:k]
    p0_x[]  = ts;  p0_y[]  = pop0_all[1:k]
    p1_x[]  = ts;  p1_y[]  = pop1_all[1:k]
end

## 1:traj.N 
animate_figure(fig_anim, 1:traj.N, update_pulse_and_pop!;
    fps = 24, mode = :record, filename = "pulse_with_populations.mp4")

# ## Step 6 — Parameter sweep: `rollout(qtraj, pulse)`
#
# Piccolo's 
#
# ```julia
# qtraj_new = rollout(qtraj, new_pulse)   #
# fid = fidelity(qtraj_new)              # 
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

# Explicit Observables for sweep lines and scatter dot
ux_s_x = Observable(plot_times);        ux_s_y = Observable(ctrls_list[1][1, :])
uy_s_x = Observable(plot_times);        uy_s_y = Observable(ctrls_list[1][2, :])
lines!(ax_ctrl_s, ux_s_x, ux_s_y; label = "ux", linewidth = 2, color = :steelblue)
lines!(ax_ctrl_s, uy_s_x, uy_s_y; label = "uy", linewidth = 2, color = :firebrick)
axislegend(ax_ctrl_s; position = :rt)

dot_x = Observable([alphas[1]]);  dot_y = Observable([fids[1]])
scatter!(ax_fid_s, dot_x, dot_y; markersize = 18, color = :red)

lbl_s = Label(fig_sweep[3, 1],
    "alpha=$(round(alphas[1]; digits=2))  fidelity=$(round(fids[1]; digits=4))";
    fontsize = 15)

function update_sweep!(k)
    ux_s_x[]  = plot_times;   ux_s_y[]  = ctrls_list[k][1, :]
    uy_s_x[]  = plot_times;   uy_s_y[]  = ctrls_list[k][2, :]
    dot_x[]   = [alphas[k]];  dot_y[]   = [fids[k]]
    lbl_s.text[] = "alpha=$(round(alphas[k]; digits=2))" *
                   "  fidelity=$(round(fids[k]; digits=4))"
end

animate_figure(fig_sweep, 1:length(alphas), update_sweep!;
    fps = 12, mode = :record, filename = "amplitude_sweep.mp4")

# 
# ## Further reading
#
# - Piccolo quickstart guide
# - [`rollout` and `fidelity` API](../../api/rollouts.md)
# - [`plot_unitary_populations`](../../api/visualizations.md)
# - [`interactive_pulse_dashboard.jl`](../interactive_pulse_dashboard.jl)

# - [GLMakie Observables](https://docs.makie.org/stable/documentation/observables/)
