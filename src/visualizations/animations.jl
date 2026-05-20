module AnimatedPlots

export animate_figure
export animate_name

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

Low-level animation primitive providing flexible control over how each frame is rendered.
For trajectory-specific animations, see `animate_name`, `animate_bloch`, or `animate_wigner`.

Defined as a stub here; the implementation is provided by the `PiccoloMakieExt` extension
and is loaded when a Makie backend (`CairoMakie`, `GLMakie`, `WGLMakie`) is available.

# Arguments
- `fig::Figure`: The Makie figure to animate. Should already contain initial plot elements.
- `frames::AbstractVector{Int}`: Frame indices to iterate through (e.g., `1:100`).
- `update_frame!::Function`: Callback `(i) -> ...` that mutates `fig` for frame `i`.

# Keyword Arguments
- `mode::Symbol`: `:inline` (interactive window, requires `GLMakie`) or `:record`
  (save to file, works with `CairoMakie` or `GLMakie`). Default `:inline`.
- `fps::Int`: Frames per second. Default `24`.
- `filename::String`: Output path when `mode = :record`. Supported: `.mp4`, `.gif`,
  `.webm`, `.mkv`. Default `"animation.mp4"`.

# Returns
- The input `Figure` (potentially mutated).

# Backend Compatibility
- **GLMakie**: both `:inline` and `:record`.
- **CairoMakie**: `:record` only.
- **WGLMakie**: `:inline` in Jupyter / web.

# Examples

```julia
using CairoMakie  # for :record mode

fig = Figure()
ax = Axis(fig[1, 1])

function update_frame!(i)
    empty!(ax)
    θ = range(0, 2π, length=100)
    r = i / 50
    lines!(ax, r * cos.(θ), r * sin.(θ))
end

animate_figure(fig, 1:100, update_frame!; mode=:record, filename="circle.mp4")
```

See also: `animate_name`, `animate_bloch`, `animate_wigner`.
"""
function animate_figure end

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

Creates an animation showing how a variable in a `NamedTrajectory` evolves, with data
appearing left-to-right. Each frame shows data from `1:i` for the current index `i`.

Defined as a stub here; the implementation is provided by the `PiccoloMakieExt` extension.

# Arguments
- `traj::NamedTrajectory`: The trajectory containing the data to animate.
- `name::Symbol`: Component name (e.g., `:u`, `:ψ̃`).

# Keyword Arguments
- `fps::Int`: Frames per second. Default `24`.
- `mode::Symbol`: `:inline` (requires `GLMakie`) or `:record`. Default `:inline`.
- `filename::String`: Output path when `mode = :record`. Default `"name_animation.mp4"`.
- `kwargs...`: Forwarded to the internal plotting call.

# Returns
- A Makie `Figure` containing the animation.

# Notes
- **CairoMakie cannot re-render in `:inline` mode** — use `:record` or switch to GLMakie.
- The animation loops continuously in `:inline` mode until the window is closed.

# Examples

```julia
using CairoMakie

# Save controls animation to file
animate_name(traj, :u; mode=:record, filename="controls.mp4", fps=24)
```

```julia
using GLMakie

# Interactive controls animation
animate_name(traj, :u; fps=30)
```

See also: `animate_figure`, `animate_bloch`, and
[`NamedTrajectories.plot`](https://docs.harmoniqs.co/NamedTrajectories/dev/generated/plotting/).
"""
function animate_name end

end
