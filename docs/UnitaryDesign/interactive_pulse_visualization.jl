using Piccolo
using GLMakie
GLMakie.activate!()
using Random

Random.seed!(42)

# ==============================
# Step 1 — Define system
# ==============================

H_drift  = 0.5 * PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
sys      = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

T = 10.0
N = 100
times = collect(range(0, T; length = N))

pulse  = ZeroOrderPulse(0.1 .* randn(2, N), times)
qtraj  = UnitaryTrajectory(sys, pulse, GATES[:X])
qcp    = SmoothPulseProblem(qtraj, N; Q = 100.0, R = 1e-2, ddu_bound = 1.0)

solve!(qcp; max_iter = 50, verbose = false)

traj = get_trajectory(qcp)

@info "Optimised fidelity:" fidelity(qcp)

plot_times = cumsum([0.0; get_timesteps(traj)])[1:end-1]
controls   = traj[:u]

# ==============================
# Step 2 — Static interactive plot
# ==============================

fig_static = plot_unitary_populations(traj)
display(fig_static)

# ==============================
# Step 3 — Progressive pulse animation
# ==============================

display(
    animate_name(traj, :u;
        fps = 24,
        mode = :inline)   # FIXED
)

# ==============================
# Step 4 — Co-animate pulse + populations
# ==============================

Utilde = traj[:Ũ⃗]

pops_all = [abs2.(iso_vec_to_operator(Utilde[:, i])[:, 1]) for i in 1:traj.N]
pop0_all = [p[1] for p in pops_all]
pop1_all = [p[2] for p in pops_all]

fig_anim = Figure(size = (1000, 800))

ax_ctrl = Axis(fig_anim[1, 1],
    xlabel = "Time",
    ylabel = "Amplitude",
    title  = "Control Pulses")

xlims!(ax_ctrl, plot_times[1], plot_times[end])

ax_pop = Axis(fig_anim[2, 1],
    xlabel = "Time",
    ylabel = "|Uij|²",
    title  = "Unitary Populations — column 1")

xlims!(ax_pop, plot_times[1], plot_times[end])
ylims!(ax_pop, 0, 1)

# Observables
ux_x = Observable(Float64[])
ux_y = Observable(Float64[])
uy_x = Observable(Float64[])
uy_y = Observable(Float64[])

lines!(ax_ctrl, ux_x, ux_y; linewidth = 2, label = "ux")
lines!(ax_ctrl, uy_x, uy_y; linewidth = 2, label = "uy")
axislegend(ax_ctrl)

p0_x = Observable(Float64[])
p0_y = Observable(Float64[])
p1_x = Observable(Float64[])
p1_y = Observable(Float64[])

lines!(ax_pop, p0_x, p0_y; linewidth = 2, label = "|U₁₁|²")
lines!(ax_pop, p1_x, p1_y; linewidth = 2, label = "|U₂₁|²")
axislegend(ax_pop)

display(fig_anim)

function update_pulse_and_pop!(k)
    ts = plot_times[1:k]
    ux_x[] = ts
    ux_y[] = controls[1, 1:k]
    uy_x[] = ts
    uy_y[] = controls[2, 1:k]
    p0_x[] = ts
    p0_y[] = pop0_all[1:k]
    p1_x[] = ts
    p1_y[] = pop1_all[1:k]
end

animate_figure(fig_anim, 1:traj.N, update_pulse_and_pop!;
    fps = 24,
    mode = :inline)   # FIXED

# ==============================
# Step 5 — Amplitude sweep
# ==============================

alphas = collect(range(0.4, 1.6; length = 40))
fids = Float64[]
ctrls_list = Matrix{Float64}[]

for alpha in alphas
    scaled_pulse = ZeroOrderPulse(alpha .* controls, times)
    new_traj = rollout(qtraj, scaled_pulse)
    push!(fids, fidelity(new_traj))
    push!(ctrls_list, alpha .* controls)
end

fig_sweep = Figure(size = (1000, 800))

ax_ctrl_s = Axis(fig_sweep[1, 1],
    xlabel = "Time",
    ylabel = "Amplitude",
    title  = "Scaled Pulses")

xlims!(ax_ctrl_s, plot_times[1], plot_times[end])

ax_fid_s = Axis(fig_sweep[2, 1],
    xlabel = "α",
    ylabel = "Fidelity",
    title  = "Gate Fidelity vs Amplitude")

ylims!(ax_fid_s, 0, 1)
lines!(ax_fid_s, alphas, fids; linewidth = 2)

ux_s_x = Observable(plot_times)
ux_s_y = Observable(ctrls_list[1][1, :])
uy_s_x = Observable(plot_times)
uy_s_y = Observable(ctrls_list[1][2, :])

lines!(ax_ctrl_s, ux_s_x, ux_s_y; linewidth = 2, label = "ux")
lines!(ax_ctrl_s, uy_s_x, uy_s_y; linewidth = 2, label = "uy")
axislegend(ax_ctrl_s)

dot_x = Observable([alphas[1]])
dot_y = Observable([fids[1]])
scatter!(ax_fid_s, dot_x, dot_y; markersize = 20)

display(fig_sweep)

function update_sweep!(k)
    ux_s_y[] = ctrls_list[k][1, :]
    uy_s_y[] = ctrls_list[k][2, :]
    dot_x[]  = [alphas[k]]
    dot_y[]  = [fids[k]]
end

animate_figure(fig_sweep, 1:length(alphas), update_sweep!;
    fps = 12,
    mode = :inline)   # FIXED

# ==============================
# Step 6 — Real-time slider dashboard
# ==============================

fig_dash = Figure(size = (1100, 750))

ax_dash = Axis(fig_dash[1, 1],
    xlabel = "Time",
    ylabel = "Amplitude",
    title  = "Interactive Pulse")

xlims!(ax_dash, plot_times[1], plot_times[end])

amp_sl   = Slider(fig_dash[2, 1], range = 0.1:0.05:2.0, startvalue = 1.0)
freq_sl  = Slider(fig_dash[3, 1], range = 0.1:0.1:3.0, startvalue = 1.0)
phase_sl = Slider(fig_dash[4, 1], range = 0:0.1:2π, startvalue = 0.0)

ux_obs = Observable(controls[1, :])
uy_obs = Observable(controls[2, :])

lines!(ax_dash, plot_times, ux_obs; linewidth = 2, label = "ux")
lines!(ax_dash, plot_times, uy_obs; linewidth = 2, label = "uy")
axislegend(ax_dash)

fid_label = Label(fig_dash[5, 1],
    "Fidelity: $(round(fidelity(qcp), digits=4))",
    fontsize = 18)

display(fig_dash)

function rebuild_controls(amp, freq, phase)
    u1 = amp .* sin.(freq .* plot_times .+ phase)
    u2 = amp .* cos.(freq .* plot_times .+ phase)
    return vcat(u1', u2')
end

onany(amp_sl.value, freq_sl.value, phase_sl.value) do amp, freq, phase
    new_ctrls = rebuild_controls(amp, freq, phase)

    ux_obs[] = new_ctrls[1, :]
    uy_obs[] = new_ctrls[2, :]

    new_pulse = ZeroOrderPulse(new_ctrls, times)
    new_traj  = rollout(qtraj, new_pulse)

    fid = fidelity(new_traj)
    fid_label.text[] = "Fidelity: $(round(fid, digits=4))"
end

println("Interactive dashboard running.")
wait(Condition())