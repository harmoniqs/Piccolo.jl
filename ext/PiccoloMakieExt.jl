module PiccoloMakieExt

using Piccolo
using Makie
using NamedTrajectories
using TestItems

# Animation implementations - extend Piccolo stubs

"""
    animate_figure(
        fig::Figure,
        frames::AbstractVector{Int},
        update_frame!::Function;
        mode::Symbol = :inline,
        fps::Int = 24,
        filename::String = "animation.mp4"
    ) -> Figure

Animate a Makie figure by updating frames with a custom update function.

This is a low-level animation function that provides flexible control over how each frame
is rendered. For trajectory-specific animations, see [`animate_name`](@ref), [`animate_bloch`](@ref),
or [`animate_wigner`](@ref).

# Arguments
- `fig::Figure`: The Makie figure to animate. Should already contain initial plot elements.
- `frames::AbstractVector{Int}`: Vector of frame indices to iterate through (e.g., `1:100` or `[1, 5, 10]`).
- `update_frame!::Function`: Callback function that takes a frame index and updates the figure.
  Should modify `fig` in place and return nothing.

# Keyword Arguments
- `mode::Symbol`: Animation mode. Options:
  - `:inline` - Display animation in an interactive window (requires `GLMakie`). Loops continuously.
  - `:record` - Save animation to file (works with `CairoMakie` or `GLMakie`).
  Default is `:inline`.
- `fps::Int`: Frames per second for animation playback. Default is `24`.
- `filename::String`: Output filename when `mode = :record`. Default is `"animation.mp4"`.
  Supported formats: `.mp4`, `.gif`, `.webm`, `.mkv`.

# Returns
- The input `Figure` object (potentially modified).

# Backend Compatibility
- **GLMakie**: Supports both `:inline` and `:record` modes.
- **CairoMakie**: Only supports `:record` mode. Will warn if `:inline` is attempted.
- **WGLMakie**: Supports `:inline` in Jupyter/web environments.

# Examples

## Basic line animation
```julia
using GLMakie  # Required for :inline mode

# Create initial figure
fig = Figure()
ax = Axis(fig[1, 1], xlabel="x", ylabel="sin(x + t)")
xs = range(0, 2π, length=100)

# Initial plot
lines!(ax, xs, sin.(xs))

# Animate with phase shift
function update_frame!(t)
    empty!(ax)
    lines!(ax, xs, sin.(xs .+ t/10))
end

animate_figure(fig, 1:100, update_frame!; fps=30)
```

## Save animation to file
```julia
using CairoMakie  # Works with CairoMakie for recording

fig = Figure()
ax = Axis(fig[1, 1])

# Animate growing circle
function update_frame!(i)
    empty!(ax)
    θ = range(0, 2π, length=100)
    r = i / 50  # Grow from 0 to 2
    lines!(ax, r * cos.(θ), r * sin.(θ))
end

animate_figure(fig, 1:100, update_frame!; mode=:record, filename="circle.mp4")
```

## Animate trajectory data
```julia
using GLMakie
using NamedTrajectories

# Assume `traj` is a NamedTrajectory with :u controls
traj = # ... your trajectory
times = get_times(traj)

fig = Figure()
ax = Axis(fig[1, 1], xlabel="Time", ylabel="Control")

# Initialize empty plot
lines!(ax, Float64[], Float64[])

function update_frame!(k)
    empty!(ax)
    lines!(ax, times[1:k], traj[:u][1, 1:k])
end

animate_figure(fig, 1:traj.N, update_frame!)
```

See also: [`animate_name`](@ref), [`animate_bloch`](@ref), [`animate_wigner`](@ref)
"""
function Piccolo.animate_figure(
    fig::Figure,
    frames::AbstractVector{Int},
    update_frame!::Function;
    mode::Symbol = :inline,
    fps::Int = 24,
    filename::String = "animation.mp4",
)
    if mode == :inline
        display(fig) # open the scene
        if !isopen(fig.scene)
            @warn "Unable to open :inline animation for the current backend. This is a known limitation of CairoMakie. Consider setting mode = :record."
        end
        @async begin # don't block
            while isopen(fig.scene)
                for i in frames
                    update_frame!(i)
                    sleep(1 / fps)
                end
                if !isopen(fig.scene)
                    break # exit after close
                else
                    sleep(10 / fps)
                end
            end
        end
    elseif mode == :record
        record(fig, filename, frames; framerate = fps) do i
            update_frame!(i)
        end
    else
        throw(ArgumentError("Unsupported mode: $mode. Use :inline or :record."))
    end

    return fig
end


"""
    animate_name(
        traj::NamedTrajectory,
        name::Symbol;
        fps::Int = 24,
        mode::Symbol = :inline,
        filename::String = "name_animation.mp4",
        kwargs...
    ) -> Figure

Animate the time evolution of a trajectory component, progressively revealing data over time.

This function creates an animation showing how a variable in a `NamedTrajectory` evolves,
with data appearing progressively from left to right. Each frame displays data from the
start up to the current time index.

# Arguments
- `traj::NamedTrajectory`: The trajectory containing the data to animate.
- `name::Symbol`: The component name to animate (e.g., `:u` for controls, `:ψ̃` for states).

# Keyword Arguments
- `fps::Int`: Frames per second for animation. Default is `24`.
- `mode::Symbol`: Animation mode:
  - `:inline` - Display in interactive window (requires `GLMakie`).
  - `:record` - Save to file.
  Default is `:inline`.
- `filename::String`: Output filename when `mode = :record`. Default is `"name_animation.mp4"`.
- `kwargs...`: Additional arguments passed to the internal plotting function (currently unused).

# Returns
- A Makie `Figure` object containing the animation.

# Notes
- **Backend limitation**: CairoMakie cannot re-render in `:inline` mode. Use `:record` or switch to GLMakie.
- The animation loops continuously in `:inline` mode until the window is closed.
- For multi-dimensional data (e.g., multiple controls), all dimensions are plotted together.

# Examples

## Animate control evolution
```julia
using GLMakie
using NamedTrajectories

# Create trajectory with control pulses
N = 100
times = range(0, 10, length=N)
u = 0.5 * sin.(range(0, 4π, length=N))  # Sinusoidal control

traj = NamedTrajectory(
    (u = reshape(u, 1, N), Δt = diff([times..., times[end]])),
    controls = :u,
    timestep = :Δt
)

# Animate controls appearing over time
animate_name(traj, :u; fps=30)
```

## Save multi-control animation
```julia
using CairoMakie

# Multiple control signals
N = 150
times = range(0, 15, length=N)
u1 = 0.3 * sin.(range(0, 6π, length=N))
u2 = 0.4 * cos.(range(0, 4π, length=N))

traj = NamedTrajectory(
    (u = [u1'; u2'], Δt = fill(times[2] - times[1], N)),
    controls = :u,
    timestep = :Δt
)

animate_name(traj, :u; mode=:record, filename="controls.mp4", fps=24)
```

## Animate quantum state evolution
```julia
using GLMakie

# Trajectory with quantum states (e.g., from a solved problem)
traj = get_trajectory(qcp)  # Assume qcp is a solved QuantumControlProblem

# Animate state populations appearing
animate_name(traj, :ψ̃; fps=20)
```

See also: [`animate_figure`](@ref), [`animate_bloch`](@ref), [`NamedTrajectories.plot`](https://docs.harmoniqs.co/NamedTrajectories/dev/generated/plotting/)
"""
function Piccolo.animate_name(
    traj::NamedTrajectory,
    name::Symbol;
    fps::Int = 24,
    mode::Symbol = :inline,
    filename = "name_animation.mp4",
    kwargs...,
)
    fig = Figure()
    ax = Axis(fig[1, 1])

    # Set limits
    times = get_times(traj)
    xlims = [times[1], times[end]]
    ylims = collect(extrema(traj[name]))

    scatter!(ax, xlims, ylims, markersize = 0.0)
    QuantumObjectPlots.plot_name!(ax, traj, name, indices = 1:1, kwargs...)

    # TODO: Unclear how to set observables via indices, so redraw
    # WARNING: Known bug with CairoMakie, which cannot re-render
    function update_frame!(i)
        empty!(ax.scene)
        scatter!(ax, xlims, ylims, markersize = 0.0)
        QuantumObjectPlots.plot_name!(ax, traj, name, indices = 1:i)
    end

    return Piccolo.animate_figure(
        fig,
        1:traj.N,
        update_frame!,
        mode = mode,
        fps = fps,
        filename = filename,
    )
end


@testitem "Test animate_name for NamedTrajectory" begin
    using QuantumToolbox
    using NamedTrajectories
    using CairoMakie

    x = ComplexF64[1.0; 0.0]
    y = ComplexF64[0.0, 1.0]

    comps = (ψ̃ = hcat(ket_to_iso(x), ket_to_iso(y)), Δt = [1.0; 1.0])
    traj = NamedTrajectory(comps)

    fig = animate_name(traj, :ψ̃, mode = :inline, fps = 10)
    @test fig isa Figure
end


end
