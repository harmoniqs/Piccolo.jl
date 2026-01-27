# [Visualizations](@id viz-api)

The Visualizations module provides plotting utilities for analyzing quantum control solutions.

## Trajectory Plotting

The primary plotting interface is provided by NamedTrajectories.jl:

```julia
using CairoMakie
using Piccolo

# Plot specific components
fig = plot(trajectory, [:u, :Ũ⃗])

# Plot controls only
fig = plot(trajectory, [:u])

# Plot with derivatives
fig = plot(trajectory, [:u, :du, :ddu])
```

## Quantum-Specific Plots

### plot_unitary_populations

Plot the population evolution of a unitary trajectory.

```julia
fig = plot_unitary_populations(trajectory)
fig = plot_unitary_populations(trajectory; kwargs...)
```

### plot_state_populations

Plot the population evolution of a state (ket) trajectory.

```julia
fig = plot_state_populations(trajectory)
fig = plot_state_populations(trajectory; kwargs...)
```

## Animation

### animate_figure

Create an animation from a sequence of figures.

```julia
animate_figure(fig, trajectory; filename="animation.gif")
```

### animate_bloch

Animate a qubit state trajectory on the Bloch sphere.

```julia
animate_bloch(trajectory; filename="bloch.gif")
```

### animate_wigner

Animate Wigner function evolution for a bosonic trajectory.

```julia
animate_wigner(trajectory; filename="wigner.gif")
```

## Plot Customization

### Figure Size

```julia
fig = Figure(size=(800, 600))
```

### Adding to Existing Figure

```julia
fig = Figure()
ax = Axis(fig[1, 1])
# Add custom content
```

### Saving Figures

```julia
save("figure.png", fig)
save("figure.pdf", fig)  # Vector graphics
save("figure.svg", fig)
```

## Tips

### Makie Backends

Choose a backend based on your needs:

```julia
using CairoMakie  # Static figures, PDFs
using GLMakie     # Interactive, 3D
using WGLMakie    # Web/notebooks
```

### Common Patterns

```julia
# Extract trajectory data
traj = get_trajectory(qcp)
times = cumsum([0; get_timesteps(traj)])
controls = traj[:u]

# Manual plotting
fig = Figure()
ax = Axis(fig[1, 1], xlabel="Time", ylabel="Control")
for i in 1:size(controls, 1)
    lines!(ax, times[1:end-1], controls[i, :], label="Drive $i")
end
axislegend(ax)
fig
```
