# # [Visualizations API](@id viz-api)
#
# The Visualizations module provides plotting utilities for analyzing quantum control solutions.
# This guide demonstrates visualization capabilities with runnable examples.

# ## Setup
#
# Load the required packages:

using Piccolo
using CairoMakie
using Random
Random.seed!(42)

# Create a simple example problem to visualize:

## System
H_drift = 0.5 * PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
sys = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

## Pulse and trajectory
T = 10.0
N = 100
times = collect(range(0, T, length = N))
pulse = ZeroOrderPulse(0.1 * randn(2, N), times)
qtraj = UnitaryTrajectory(sys, pulse, GATES[:X])

## Solve
qcp = SmoothPulseProblem(qtraj, N; Q = 100.0, R = 1e-2, ddu_bound = 1.0)
cached_solve!(qcp, "viz_ref_unitary"; max_iter = 50, verbose = false, print_level = 1)

# Optimization complete:

fidelity(qcp)

# # Basic Trajectory Plotting
#
# ## Plotting Trajectory Components
#
# The primary plotting interface is provided by NamedTrajectories.jl.
# You can plot specific components of the trajectory:

traj = get_trajectory(qcp)

## Plot controls only
fig1 = plot(traj, [:u])
fig1

# ## Plotting Multiple Components
#
# Plot controls and their derivatives:

fig2 = plot(traj, [:u, :du, :ddu])
fig2

# # Quantum-Specific Plots
#
# ## plot_unitary_populations
#
# Plot the population evolution for a unitary trajectory.
# This shows how the quantum state populations evolve during the gate:

fig3 = plot_unitary_populations(traj)
fig3

# The plot shows how populations in different basis states change over time.
# For an X gate starting from |0⟩, we expect to see population transfer to |1⟩.

# ## plot_state_populations
#
# For state (ket) trajectories, use `plot_state_populations`:

## Create a state transfer problem
ψ_init = [1.0 + 0im, 0.0 + 0im]  # |0⟩
ψ_goal = [0.0 + 0im, 1.0 + 0im]  # |1⟩

pulse_ket = ZeroOrderPulse(0.1 * randn(2, N), times)
qtraj_ket = KetTrajectory(sys, pulse_ket, ψ_init, ψ_goal)
qcp_ket = SmoothPulseProblem(qtraj_ket, N; Q = 100.0, R = 1e-2, ddu_bound = 1.0)
cached_solve!(qcp_ket, "viz_ref_ket"; max_iter = 50, verbose = false, print_level = 1)

traj_ket = get_trajectory(qcp_ket)
fig4 = plot_state_populations(traj_ket)
fig4

# # Custom Plotting
#
# ## Manual Control Plotting
#
# For custom visualizations, extract trajectory data and use Makie directly:

plot_times = cumsum([0; get_timesteps(traj)])[1:(end-1)]
controls = traj[:u]

fig5 = Figure(size = (800, 400))
ax = Axis(
    fig5[1, 1],
    xlabel = "Time",
    ylabel = "Control Amplitude",
    title = "Custom Control Plot",
)

for i = 1:size(controls, 1)
    lines!(ax, plot_times, controls[i, :], label = "Drive $i", linewidth = 2)
end
axislegend(ax, position = :rt)

fig5

# ## Subplot Layouts
#
# Create multi-panel figures:

fig6 = Figure(size = (1000, 600))

## Controls
ax1 = Axis(fig6[1, 1], xlabel = "Time", ylabel = "Amplitude", title = "Controls")
lines!(ax1, plot_times, traj[:u][1, :], label = "u_x", linewidth = 2)
lines!(ax1, plot_times, traj[:u][2, :], label = "u_y", linewidth = 2)
axislegend(ax1, position = :rt)

## Control derivatives
ax2 =
    Axis(fig6[1, 2], xlabel = "Time", ylabel = "Derivative", title = "Control Derivatives")
lines!(ax2, plot_times, traj[:du][1, :], label = "du_x", linewidth = 2)
lines!(ax2, plot_times, traj[:du][2, :], label = "du_y", linewidth = 2)
axislegend(ax2, position = :rt)

fig6

# ## Phase Space Plots
#
# Visualize control in phase space:

fig7 = Figure(size = (600, 600))
ax = Axis(fig7[1, 1], xlabel = "u_x", ylabel = "u_y", title = "Control Phase Space")
lines!(ax, traj[:u][1, :], traj[:u][2, :], linewidth = 2)
scatter!(
    ax,
    [traj[:u][1, 1]],
    [traj[:u][2, 1]],
    color = :green,
    markersize = 15,
    label = "Start",
)
scatter!(
    ax,
    [traj[:u][1, end]],
    [traj[:u][2, end]],
    color = :red,
    markersize = 15,
    label = "End",
)
axislegend(ax, position = :rt)
fig7

# # Saving Figures
#
# ## Saving to Files
#
# Save figures in various formats:

## PNG (raster)
save("control_plot.png", fig1)

## PDF (vector graphics)
save("control_plot.pdf", fig1)

## SVG (vector graphics)
save("control_plot.svg", fig1)

# # Animation
#
# ## Animating Trajectories
#
# Create animations to visualize time evolution.
# Note: Animation functions like `animate_bloch` and `animate_wigner`
# are available for specific visualization types.

## Example: Save an animation of control evolution
## animate_figure(fig, traj; filename="animation.gif")

# # Figure Customization
#
# ## Custom Figure Size

fig_custom = Figure(size = (1200, 800))

# ## Custom Colors and Styles

fig8 = Figure(size = (800, 400))
ax = Axis(fig8[1, 1], xlabel = "Time", ylabel = "Amplitude", title = "Styled Plot")

lines!(
    ax,
    plot_times,
    traj[:u][1, :],
    color = :blue,
    linewidth = 3,
    linestyle = :solid,
    label = "u_x",
)
lines!(
    ax,
    plot_times,
    traj[:u][2, :],
    color = :red,
    linewidth = 3,
    linestyle = :dash,
    label = "u_y",
)

axislegend(ax, position = :rt, backgroundcolor = (:white, 0.8))
fig8

# # Makie Backends
#
# Choose a backend based on your needs:

## Other backends available:
## using GLMakie     # Interactive, 3D visualizations
## using WGLMakie    # Web-based, for Jupyter notebooks
## Currently using CairoMakie for static figures

# # Tips and Common Patterns
#
# ## Extracting Time and Control Data

timesteps = get_timesteps(traj)
cumulative_times = cumsum([0; timesteps])
control_times = cumulative_times[1:(end-1)]
control_values = traj[:u]

# Trajectory data:

println("Number of timesteps: ", length(timesteps))
println("Total duration: ", sum(timesteps))
println("Control array size: ", size(control_values)) # nothing

# ## Plotting Multiple Trajectories

## Create a second trajectory for comparison
pulse2 = ZeroOrderPulse(0.15 * randn(2, N), times)
qtraj2 = UnitaryTrajectory(sys, pulse2, GATES[:X])
qcp2 = SmoothPulseProblem(qtraj2, N; Q = 100.0, R = 2e-2, ddu_bound = 1.0)
cached_solve!(qcp2, "viz_ref_comparison"; max_iter = 50, verbose = false, print_level = 1)
traj2 = get_trajectory(qcp2)

fig9 = Figure(size = (800, 400))
ax = Axis(fig9[1, 1], xlabel = "Time", ylabel = "u_x", title = "Comparing Two Solutions")
lines!(ax, plot_times, traj[:u][1, :], label = "R=1e-2", linewidth = 2)
lines!(ax, plot_times, traj2[:u][1, :], label = "R=2e-2", linewidth = 2)
axislegend(ax, position = :rt)
fig9

# # Summary
#
# The Visualizations module provides:
#
# - **Basic plotting**: `plot(trajectory, components)` for any trajectory variables
# - **Quantum plots**: `plot_unitary_populations`, `plot_state_populations`
# - **Animation**: `animate_figure`, `animate_bloch`, `animate_wigner`
# - **Custom plotting**: Full access to Makie's powerful plotting capabilities
# - **Multiple backends**: CairoMakie (static), GLMakie (interactive), WGLMakie (web)
# - **Export options**: PNG, PDF, SVG formats
#
# For complete function signatures and advanced usage, see the [API Reference](@ref) overview.
