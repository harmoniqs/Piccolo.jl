#!/usr/bin/env julia

using Piccolo
using PiccoloMakieExt
using GLMakie
using Observables
using LinearAlgebra

GLMakie.activate!()

# ------------------------------------------------------------
# Quantum System Definition (NO custom propagation)
# ------------------------------------------------------------

H_drift = 0.5 * PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]

sys = QuantumSystem(
    H_drift,
    H_drives,
    [1.0, 1.0]
)

target = GATES[:X]

# ------------------------------------------------------------
# Time Grid
# ------------------------------------------------------------

T = 10.0
N = 200
times = range(0, T, length=N)

# ------------------------------------------------------------
# GUI Figure
# ------------------------------------------------------------

fig = Figure(size = (1200, 800))
Label(fig[0, :], "Interactive Pulse Animation Dashboard",
      fontsize = 24)

# ------------------------------------------------------------
# Sliders
# ------------------------------------------------------------

sg = SliderGrid(
    fig[1, 1],
    (label = "Amplitude", range = 0.0:0.05:2.0, startvalue = 0.8),
    (label = "Frequency", range = 0.5:0.1:5.0, startvalue = 1.5),
    (label = "Phase", range = 0:0.1:2π, startvalue = 0.0),
    (label = "Chirp", range = -1.0:0.1:1.0, startvalue = 0.0),
    width = 350
)

amp_obs   = sg.sliders[1].value
freq_obs  = sg.sliders[2].value
phase_obs = sg.sliders[3].value
chirp_obs = sg.sliders[4].value

# ------------------------------------------------------------
# Pulse Construction (reactive, lightweight only)
# ------------------------------------------------------------

controls_obs = lift(amp_obs, freq_obs, phase_obs, chirp_obs) do A, f, ϕ, χ

    ω_t = @. f + χ * times / T
    phase = cumsum(ω_t) * (T / N)

    u₁ = A .* cos.(phase .+ ϕ)
    u₂ = A .* sin.(phase .+ ϕ)

    hcat(u₁, u₂)'   # match Piccolo control convention
end

# ------------------------------------------------------------
# Rollout via Piccolo internal API (NO reimplementation)
# ------------------------------------------------------------

traj_obs = lift(controls_obs) do controls

    pulse = ZeroOrderPulse(controls, times)

    rollout(sys;
        pulse = pulse,
        target = target,
        save_trajectory = true
    )
end

# ------------------------------------------------------------
# Fidelity (reactive but cheap)
# ------------------------------------------------------------

fid_obs = lift(traj_obs) do traj
    unitary_fidelity(traj)
end

fid_label = Label(fig[1, 2], "Fidelity: 0.0", fontsize = 20)

on(fid_obs) do fid
    fid_label.text[] = "Fidelity: $(round(fid, digits=6))"
end

# ------------------------------------------------------------
# Pulse Plot (clean reactive update)
# ------------------------------------------------------------

ax_pulse = Axis(fig[2, 1],
    xlabel = "Time",
    ylabel = "Amplitude",
    title  = "Control Pulses"
)

u1_obs = lift(c -> c[1, :], controls_obs)
u2_obs = lift(c -> c[2, :], controls_obs)

lines!(ax_pulse, times, u1_obs, linewidth = 2)
lines!(ax_pulse, times, u2_obs, linewidth = 2)

axislegend(ax_pulse, ["uₓ", "uᵧ"])

# ------------------------------------------------------------
# Populations via Piccolo built-in visualization
# (NO manual unitary extraction)
# ------------------------------------------------------------

ax_pop = Axis(fig[3, 1],
    xlabel = "Time",
    ylabel = "Population",
    title  = "Unitary Populations"
)

plot_unitary_populations!(ax_pop, traj_obs)

# ------------------------------------------------------------
# Animate Pulse Drawing Itself Over Time
# Reuse animate_figure instead of custom animation
# ------------------------------------------------------------

animate_button = Button(fig[2, 2], label = "Play Animation")

on(animate_button.clicks) do _

    animate_figure(
        fig,
        traj_obs[],
        :controls;       # animate control drawing
        framerate = 30,
        mode = :inline   # can switch to :record for MP4
    )
end

# ------------------------------------------------------------
# Display
# ------------------------------------------------------------

display(fig)

println("Dashboard running...")
println("Adjust sliders to update pulse and populations.")
println("Click 'Play Animation' to animate pulse drawing.")