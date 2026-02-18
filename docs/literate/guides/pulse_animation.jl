# # [Pulse Animation with Parameter Evolution](@id pulse_animation)
#
# This guide demonstrates how to create engaging animations of control pulses that evolve
# as parameters change. We'll show pulses "drawing" themselves over time, quantum state
# evolution alongside the pulse, and interactive parameter control with sliders.
#
#md # !!! tip "What You'll Learn"
#md #     - Animate pulses progressively revealing over time
#md #     - Show quantum state evolution synchronized with pulse
#md #     - Create parameter sweeps (amplitude, frequency, duration)
#md #     - Build interactive dashboards with GLMakie sliders

# ## Setup

using Piccolo
using CairoMakie  # For recorded animations
using Random
Random.seed!(42)

# ## Basic Pulse Animation
#
# Let's start by animating a pulse "drawing" itself from left to right, showing how
# the control evolves over time.

## Create a simple qubit system
H_drift = 0.5 * PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
sys = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

## Design an X gate pulse
T = 10.0
N = 100
times = collect(range(0, T, length = N))
initial_controls = 0.1 * randn(2, N)
pulse = ZeroOrderPulse(initial_controls, times)
qtraj = UnitaryTrajectory(sys, pulse, GATES[:X])

qcp = SmoothPulseProblem(qtraj, N; Q = 100.0, R = 1e-2, ddu_bound = 1.0)
solve!(qcp; max_iter = 50, verbose = false, print_level = 1)

traj = get_trajectory(qcp)
plot_times = cumsum([0; get_timesteps(traj)])[1:(end-1)]

#md # !!! info "Progressive Reveal Animation"
#md #     Using `animate_name`, we can show the pulse appearing progressively.
#md #     This creates a "drawing" effect that's visually engaging.

## Uncomment to run animation (requires GLMakie for :inline mode)
# using GLMakie
# fig_pulse = animate_name(traj, :u; fps=30, mode=:inline)

## For recorded animation:
fig_pulse = animate_name(traj, :u; fps=24, mode=:record, filename="pulse_evolution.mp4")

#md # The animation shows each control amplitude appearing from left to right.

# ## Pulse with State Evolution
#
# More informative: show the pulse AND the quantum state populations evolving together.

fig = Figure(size = (1200, 800))

## Control subplot (top)
ax_control = Axis(
    fig[1, 1],
    xlabel = "Time",
    ylabel = "Control Amplitude",
    title = "Control Pulses"
)

## Population subplot (bottom)
ax_pop = Axis(
    fig[2, 1],
    xlabel = "Time",
    ylabel = "Population |U_{i,j}|²",
    title = "Unitary Populations"
)

## Set x-axis limits
xlims!(ax_control, plot_times[1], plot_times[end])
xlims!(ax_pop, plot_times[1], plot_times[end])
ylims!(ax_pop, 0, 1)

## Initialize empty plots
line_u1 = lines!(ax_control, Float64[], Float64[], label="u_x", linewidth=2, color=:blue)
line_u2 = lines!(ax_control, Float64[], Float64[], label="u_y", linewidth=2, color=:red)
axislegend(ax_control, position=:rt)

line_p1 = lines!(ax_pop, Float64[], Float64[], label="|U₁₁|²", linewidth=2, color=:green)
line_p2 = lines!(ax_pop, Float64[], Float64[], label="|U₂₁|²", linewidth=2, color=:orange)
axislegend(ax_pop, position=:rt)

## Animation function
function update_pulse_and_state!(k)
    ## Update controls
    line_u1[1][] = plot_times[1:k]
    line_u1[2][] = traj[:u][1, 1:k]
    line_u2[1][] = plot_times[1:k]
    line_u2[2][] = traj[:u][2, 1:k]
    
    ## Update populations (first column of unitary)
    Ũ⃗ = traj[:Ũ⃗]
    U_dim = sys.isodim
    pops = [abs2.(iso_vec_to_operator(Ũ⃗[:, i])[:, 1]) for i in 1:k]
    
    line_p1[1][] = plot_times[1:k]
    line_p1[2][] = [p[1] for p in pops]
    line_p2[1][] = plot_times[1:k]
    line_p2[2][] = [p[2] for p in pops]
end

## Animate
animate_figure(fig, 1:traj.N, update_pulse_and_state!; fps=24, mode=:record, filename="pulse_with_states.mp4")

#md # !!! tip "Synchronized Animation"
#md #     The top panel shows control pulses, bottom shows quantum populations.
#md #     You can see how the pulse shapes drive the quantum state evolution!

# ## Parameter Sweep: Amplitude Scaling
#
# Animate how changing pulse amplitude affects the gate fidelity.

## Create pulses with different amplitudes
amplitudes = range(0.5, 1.5, length=30)
fidelities = Float64[]
pulse_data = []

base_controls = traj[:u]

for amp in amplitudes
    ## Scale pulse
    scaled_controls = amp * base_controls
    scaled_pulse = ZeroOrderPulse(scaled_controls, times)
    scaled_traj = UnitaryTrajectory(sys, scaled_pulse, GATES[:X])
    
    ## Compute fidelity
    fid = unitary_fidelity(scaled_traj)
    push!(fidelities, fid)
    push!(pulse_data, scaled_controls)
end

#md # Now create animation showing pulse changing with amplitude parameter:

fig_sweep = Figure(size = (1200, 900))

## Pulse plot
ax_pulse = Axis(
    fig_sweep[1, 1],
    xlabel = "Time",
    ylabel = "Control Amplitude",
    title = "Control Pulses"
)

## Fidelity vs amplitude plot
ax_fid = Axis(
    fig_sweep[2, 1],
    xlabel = "Amplitude Scaling Factor",
    ylabel = "Gate Fidelity",
    title = "Fidelity vs. Amplitude"
)

## Plot fidelity curve
lines!(ax_fid, amplitudes, fidelities, linewidth=2, color=:black)
ylims!(ax_fid, 0, 1)

## Current point marker
scatter_current = scatter!(ax_fid, [amplitudes[1]], [fidelities[1]], 
                           markersize=20, color=:red)

## Pulse lines
line_sweep_u1 = lines!(ax_pulse, plot_times, pulse_data[1][1, :], 
                       label="u_x", linewidth=2, color=:blue)
line_sweep_u2 = lines!(ax_pulse, plot_times, pulse_data[1][2, :], 
                       label="u_y", linewidth=2, color=:red)
axislegend(ax_pulse, position=:rt)

## Amplitude label
amp_label = Label(fig_sweep[0, 1], "Amplitude: $(round(amplitudes[1], digits=2))", 
                  fontsize=20, font=:bold)

## Update function
function update_amplitude_sweep!(k)
    ## Update pulse
    line_sweep_u1[2][] = pulse_data[k][1, :]
    line_sweep_u2[2][] = pulse_data[k][2, :]
    
    ## Update marker
    scatter_current[1][] = [amplitudes[k]]
    scatter_current[2][] = [fidelities[k]]
    
    ## Update label
    amp_label.text[] = "Amplitude: $(round(amplitudes[k], digits=2)) | Fidelity: $(round(fidelities[k], digits=4))"
end

## Animate
animate_figure(fig_sweep, 1:length(amplitudes), update_amplitude_sweep!; 
               fps=10, mode=:record, filename="amplitude_sweep.mp4")

#md # !!! info "Parameter Sweep"
#md #     This animation shows how pulse amplitude affects gate fidelity. The red dot
#md #     moves along the fidelity curve while the pulse scales up and down. There's
#md #     typically an optimal amplitude for maximum fidelity!

# ## Frequency Modulation Animation
#
# Show a pulse with time-varying frequency (chirp).

## Create chirped pulse
function create_chirp_pulse(ω₀, ω₁, T, N)
    times = range(0, T, length=N)
    Δt = T / N
    
    ## Instantaneous frequency sweeps linearly
    ω_t = [ω₀ + (ω₁ - ω₀) * t / T for t in times]
    
    ## Integrate to get phase
    phase = cumsum(ω_t .* Δt)
    
    ## Create amplitude-modulated pulse
    amplitude = 0.5
    u_x = amplitude * cos.(phase)
    u_y = amplitude * sin.(phase)
    
    return [u_x'; u_y'], times
end

ω₀ = 1.0  # Start frequency
ω₁ = 3.0  # End frequency
chirp_controls, chirp_times = create_chirp_pulse(ω₀, ω₁, T, N)

## Visualize chirp with spectrogram-style animation
fig_chirp = Figure(size = (1200, 600))

ax_chirp = Axis(
    fig_chirp[1, 1],
    xlabel = "Time",
    ylabel = "Control Amplitude",
    title = "Frequency-Swept Pulse (Chirp)"
)

line_chirp_u1 = lines!(ax_chirp, Float64[], Float64[], 
                       label="u_x (I)", linewidth=2, color=:blue)
line_chirp_u2 = lines!(ax_chirp, Float64[], Float64[], 
                       label="u_y (Q)", linewidth=2, color=:red)
axislegend(ax_chirp, position=:rt)

## Current time marker
vlines!(ax_chirp, [0.0], color=:black, linestyle=:dash, linewidth=2)

freq_label = Label(fig_chirp[0, 1], "Current Frequency: $(round(ω₀, digits=2)) rad/s", 
                   fontsize=20, font=:bold)

function update_chirp!(k)
    line_chirp_u1[1][] = chirp_times[1:k]
    line_chirp_u1[2][] = chirp_controls[1, 1:k]
    line_chirp_u2[1][] = chirp_times[1:k]
    line_chirp_u2[2][] = chirp_controls[2, 1:k]
    
    ## Update frequency label
    t = chirp_times[k]
    ω_current = ω₀ + (ω₁ - ω₀) * t / T
    freq_label.text[] = "Current Frequency: $(round(ω_current, digits=2)) rad/s | Time: $(round(t, digits=2)) s"
end

animate_figure(fig_chirp, 1:N, update_chirp!; 
               fps=24, mode=:record, filename="chirp_pulse.mp4")

#md # !!! tip "Frequency Sweeps"
#md #     Chirped pulses are useful for robust control and adiabatic passage. The
#md #     frequency increases linearly, visible as increasing oscillation frequency.

# ## Multi-Pulse Sequence Animation
#
# Animate a sequence of pulses with different parameters.

## Create sequence: X, Y, Z rotations
rotation_sequence = [GATES[:X], GATES[:Y], GATES[:X]]
pulse_sequence = []
traj_sequence = []

for (i, target_gate) in enumerate(rotation_sequence)
    T_gate = 5.0
    N_gate = 50
    times_gate = collect(range(0, T_gate, length=N_gate))
    
    ## Random initial controls
    initial_controls_seq = 0.1 * randn(2, N_gate)
    pulse_seq = ZeroOrderPulse(initial_controls_seq, times_gate)
    qtraj_seq = UnitaryTrajectory(sys, pulse_seq, target_gate)
    
    qcp_seq = SmoothPulseProblem(qtraj_seq, N_gate; Q = 100.0, R = 1e-2, ddu_bound = 1.0)
    solve!(qcp_seq; max_iter = 50, verbose = false, print_level = 1)
    
    push!(traj_sequence, get_trajectory(qcp_seq))
end

## Animate sequence
fig_seq = Figure(size = (1200, 600))

ax_seq = Axis(
    fig_seq[1, 1],
    xlabel = "Time",
    ylabel = "Control Amplitude",
    title = "Pulse Sequence Animation"
)

gate_names = ["X Gate", "Y Gate", "X Gate"]
gate_label = Label(fig_seq[0, 1], "Current: $(gate_names[1])", 
                   fontsize=22, font=:bold, color=:blue)

line_seq_u1 = lines!(ax_seq, Float64[], Float64[], 
                     label="u_x", linewidth=2, color=:blue)
line_seq_u2 = lines!(ax_seq, Float64[], Float64[], 
                     label="u_y", linewidth=2, color=:red)
axislegend(ax_seq, position=:rt)

## Build combined animation
total_frames = sum([t.N for t in traj_sequence])
frame_to_gate = []
frame_to_local_idx = []

for (gate_idx, traj_gate) in enumerate(traj_sequence)
    for local_idx in 1:traj_gate.N
        push!(frame_to_gate, gate_idx)
        push!(frame_to_local_idx, local_idx)
    end
end

function update_sequence!(frame)
    gate_idx = frame_to_gate[frame]
    local_idx = frame_to_local_idx[frame]
    
    traj_current = traj_sequence[gate_idx]
    times_current = cumsum([0; get_timesteps(traj_current)])[1:(end-1)]
    
    line_seq_u1[1][] = times_current[1:local_idx]
    line_seq_u1[2][] = traj_current[:u][1, 1:local_idx]
    line_seq_u2[1][] = times_current[1:local_idx]
    line_seq_u2[2][] = traj_current[:u][2, 1:local_idx]
    
    gate_label.text[] = "Current: $(gate_names[gate_idx]) (Pulse $(gate_idx)/3)"
end

animate_figure(fig_seq, 1:total_frames, update_sequence!; 
               fps=24, mode=:record, filename="pulse_sequence.mp4")

#md # !!! info "Pulse Sequences"
#md #     Shows consecutive gates being drawn. Useful for visualizing composite pulses
#md #     or dynamical decoupling sequences.

# ## Interactive Dashboard with Sliders
#
# This requires GLMakie for interactivity. Here's a template that users can run
# with GLMakie to get real-time parameter control.

#md # ```julia
#md # using GLMakie
#md # using Piccolo
#md # 
#md # # Create base system and trajectory
#md # H_drift = 0.5 * PAULIS[:Z]
#md # H_drives = [PAULIS[:X], PAULIS[:Y]]
#md # sys = QuantumSystem(H_drift, H_drives, [1.0, 1.0])
#md # 
#md # # Initial pulse
#md # T = 10.0
#md # N = 100
#md # times = collect(range(0, T, length=N))
#md # base_controls = 0.5 * [sin.(2π * times / T)'; cos.(2π * times / T)']
#md # 
#md # # Create interactive figure
#md # fig = Figure(size=(1200, 800))
#md # 
#md # # Sliders
#md # amplitude_slider = Slider(fig[1, 1], range=0.1:0.1:2.0, startvalue=1.0)
#md # frequency_slider = Slider(fig[2, 1], range=0.5:0.1:5.0, startvalue=1.0)
#md # phase_slider = Slider(fig[3, 1], range=0:0.1:2π, startvalue=0.0)
#md # 
#md # Label(fig[1, 1, Left()], "Amplitude:", width=100)
#md # Label(fig[2, 1, Left()], "Frequency:", width=100)
#md # Label(fig[3, 1, Left()], "Phase:", width=100)
#md # 
#md # # Plot area
#md # ax = Axis(fig[4, 1], xlabel="Time", ylabel="Control Amplitude")
#md # 
#md # # Observables for reactive updates
#md # amplitude = amplitude_slider.value
#md # frequency = frequency_slider.value
#md # phase = phase_slider.value
#md # 
#md # # Computed controls (reactive to slider changes)
#md # controls_reactive = lift(amplitude, frequency, phase) do amp, freq, ph
#md #     u_x = amp * sin.(2π * freq * times / T .+ ph)
#md #     u_y = amp * cos.(2π * freq * times / T .+ ph)
#md #     return [u_x'; u_y']
#md # end
#md # 
#md # # Plot that updates automatically
#md # u_plot = lift(controls_reactive) do controls
#md #     [times, controls[1, :]]
#md # end
#md # 
#md # lines!(ax, u_plot, label="u_x", linewidth=2, color=:blue)
#md # 
#md # # Fidelity display
#md # fidelity_label = lift(controls_reactive) do controls
#md #     pulse = ZeroOrderPulse(controls, times)
#md #     traj = UnitaryTrajectory(sys, pulse, GATES[:X])
#md #     fid = unitary_fidelity(traj)
#md #     "Fidelity: $(round(fid, digits=4))"
#md # end
#md # 
#md # Label(fig[0, 1], fidelity_label, fontsize=20, font=:bold)
#md # 
#md # display(fig)
#md # ```

#md # !!! tip "Interactive Control"
#md #     The code above creates sliders that update the pulse in real-time. Moving
#md #     the sliders instantly redraws the pulse and updates the fidelity. This is
#md #     powerful for intuition building and parameter tuning!

# ## Summary
#
# This guide demonstrated:
#
# - **Progressive pulse animation**: Pulses "drawing" themselves over time
# - **Synchronized state evolution**: Showing controls and populations together
# - **Parameter sweeps**: Amplitude, frequency, and other parameter variations
# - **Multi-pulse sequences**: Animating composite pulse sequences
# - **Interactive dashboards**: Real-time control with GLMakie sliders (template provided)
#
# These animations are useful for:
# - Understanding optimal control results
# - Presentations and publications
# - Website and blog content
# - Building intuition about quantum control
# - Debugging control problems

#md # ## Saving Animations
#md # 
#md # Use `mode=:record` with a filename:
#md # ```julia
#md # animate_figure(fig, frames, update!; mode=:record, filename="my_animation.mp4")
#md # ```
#md # 
#md # Supported formats: `.mp4`, `.gif`, `.webm`, `.mkv`

#md # ## Next Steps
#md #
#md # - Explore [Visualization Guide](@ref visualization) for more plotting options
#md # - See [System Templates](@ref system_templates) for different quantum systems
#md # - Check [Problem Templates](@ref problem_templates) for optimization setup

# ## References
#
# 1. GLMakie.jl documentation: https://docs.makie.org/stable/
# 2. Beautiful Makie: https://beautiful.makie.org/
# 3. Piccolo visualization infrastructure: PiccoloMakieExt.jl
