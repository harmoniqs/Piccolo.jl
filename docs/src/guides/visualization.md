# [Visualization](@id visualization)

Piccolo.jl provides visualization tools for analyzing optimization results. This guide covers plotting controls, states, and populations.

## Setup

Visualization requires a Makie backend:

```julia
using Piccolo
using CairoMakie  # or GLMakie for interactive plots
```

## Basic Trajectory Plotting

### Plot Controls

```julia
# After solving
qcp = SmoothPulseProblem(qtraj, N)
solve!(qcp)

traj = get_trajectory(qcp)

# Plot controls
fig = plot(traj, [:u])
```

### Plot Controls and State

```julia
# Plot multiple components
fig = plot(traj, [:u, :Ũ⃗])
```

### Plot Specific Variables

```julia
# Plot first-order derivatives
fig = plot(traj, [:du])

# Plot all derivatives
fig = plot(traj, [:u, :du, :ddu])
```

## Customizing Plots

### Figure Size

```julia
fig = Figure(size=(800, 600))
ax = Axis(fig[1, 1])
# Add custom plotting...
```

### Manual Control Plots

```julia
using CairoMakie

traj = get_trajectory(qcp)
times = cumsum([0; get_timesteps(traj)])

fig = Figure(size=(800, 400))
ax = Axis(fig[1, 1], xlabel="Time", ylabel="Control Amplitude")

# Plot each drive
for i in 1:size(traj[:u], 1)
    lines!(ax, times[1:end-1], traj[:u][i, :], label="Drive $i")
end

axislegend(ax)
fig
```

## State Population Plots

### Unitary Populations

For `UnitaryTrajectory`, visualize how unitary elements evolve:

```julia
# After solving
traj = get_trajectory(qcp)

# Plot unitary isomorphic components
fig = plot(traj, [:Ũ⃗])
```

### Ket State Populations

For `KetTrajectory`:

```julia
# Plot ket state components
fig = plot(traj, [:ψ̃])
```

### Custom Population Plot

```julia
using CairoMakie

traj = get_trajectory(qcp)
times = cumsum([0; get_timesteps(traj)])

# Convert isomorphic state to physical state at each timestep
populations = zeros(sys.levels, size(traj[:ψ̃], 2))
for k in 1:size(traj[:ψ̃], 2)
    ψ = iso_to_ket(traj[:ψ̃][:, k])
    populations[:, k] = abs2.(ψ)
end

fig = Figure(size=(800, 400))
ax = Axis(fig[1, 1], xlabel="Time", ylabel="Population")

for i in 1:sys.levels
    lines!(ax, times[1:end-1], populations[i, :], label="|$i⟩")
end

axislegend(ax)
fig
```

## Fidelity Evolution

Track fidelity during the pulse:

```julia
using CairoMakie

# Compute fidelity at each timestep
fidelities = Float64[]
for k in 1:size(traj[:Ũ⃗], 2)
    U_k = iso_vec_to_operator(traj[:Ũ⃗][:, k])
    F_k = abs(tr(U_goal' * U_k))^2 / sys.levels^2
    push!(fidelities, F_k)
end

fig = Figure(size=(800, 300))
ax = Axis(fig[1, 1], xlabel="Timestep", ylabel="Fidelity")
lines!(ax, 1:length(fidelities), fidelities)
hlines!(ax, [0.99], color=:red, linestyle=:dash, label="Target")
axislegend(ax)
fig
```

## Comparing Solutions

### Before and After Optimization

```julia
using CairoMakie

fig = Figure(size=(800, 600))

# Initial controls
ax1 = Axis(fig[1, 1], title="Initial", ylabel="Control")
lines!(ax1, initial_controls[1, :], label="Drive 1")
lines!(ax1, initial_controls[2, :], label="Drive 2")

# Optimized controls
traj = get_trajectory(qcp)
ax2 = Axis(fig[2, 1], title="Optimized", xlabel="Timestep", ylabel="Control")
lines!(ax2, traj[:u][1, :], label="Drive 1")
lines!(ax2, traj[:u][2, :], label="Drive 2")

fig
```

### Multiple Solutions

```julia
using CairoMakie

fig = Figure(size=(800, 400))
ax = Axis(fig[1, 1], xlabel="Time", ylabel="Control (Drive 1)")

# Plot solutions with different Q values
for (Q, color) in [(10, :blue), (100, :green), (1000, :red)]
    qcp = SmoothPulseProblem(qtraj, N; Q=Q)
    solve!(qcp)
    traj = get_trajectory(qcp)
    lines!(ax, traj[:u][1, :], color=color, label="Q=$Q")
end

axislegend(ax)
fig
```

## Animation

Create animations of state evolution:

```julia
using CairoMakie

traj = get_trajectory(qcp)
n_frames = size(traj[:Ũ⃗], 2)

# Create observable for animation
frame = Observable(1)

fig = Figure(size=(600, 400))
ax = Axis(fig[1, 1], title=@lift("Timestep $($frame)"))

# Plot that updates with frame
U_current = @lift(iso_vec_to_operator(traj[:Ũ⃗][:, $frame]))
heatmap!(ax, @lift(abs.($U_current)))

# Record animation
record(fig, "evolution.mp4", 1:n_frames; framerate=30) do f
    frame[] = f
end
```

## Saving Figures

```julia
# Save as PNG
save("controls.png", fig)

# Save as PDF (vector graphics)
save("controls.pdf", fig)

# Save as SVG
save("controls.svg", fig)
```

## Plotting Tips

### 1. Use Appropriate Resolution

```julia
fig = Figure(size=(1200, 800), fontsize=14)  # High-res for papers
```

### 2. Add Labels and Legends

```julia
ax = Axis(fig[1, 1],
    xlabel="Time (ns)",
    ylabel="Amplitude (GHz)",
    title="Optimized Control Pulse"
)
```

### 3. Use Consistent Colors

```julia
colors = Makie.wong_colors()
lines!(ax, data1, color=colors[1])
lines!(ax, data2, color=colors[2])
```

### 4. Show Bounds

```julia
# Show drive bounds as shaded region
band!(ax, times, -bound * ones(length(times)), bound * ones(length(times)),
      color=(:gray, 0.2), label="Bounds")
```

## See Also

- [Problem Templates](@ref problem-templates-overview) - Generating solutions to plot
- [Trajectories](@ref trajectories-concept) - Understanding trajectory structure
