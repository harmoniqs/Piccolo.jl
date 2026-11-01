module QuantumToolboxPlots

export animate_bloch
export animate_wigner
export plot_bloch!
export plot_wigner!

using TestItems

"""
    animate_bloch(
        traj::NamedTrajectory;
        fps::Int = 24,
        mode::Symbol = :inline,
        filename::String = "bloch_animation.mp4",
        state_name::Symbol = :ψ̃,
        state_type::Symbol = :ket,
        subspace::AbstractVector{Int} = 1:2,
        kwargs...
    ) -> Figure

Animate the evolution of a quantum state on the Bloch sphere.

Creates an animation showing a qubit state's trajectory on the Bloch sphere, with a
moving vector arrow tracking the current state.

Defined as a stub here; the implementation is provided by the `PiccoloQuantumToolboxExt`
extension and loaded when both `QuantumToolbox` and a Makie backend are present.

# Arguments
- `traj::NamedTrajectory`: The trajectory containing quantum states to animate.

# Keyword Arguments
- `fps::Int`: Frames per second. Default `24`.
- `mode::Symbol`: `:inline` (requires `GLMakie`) or `:record`. Default `:inline`.
- `filename::String`: Output path when `mode = :record`. Default `"bloch_animation.mp4"`.
- `state_name::Symbol`: Component name for the quantum state. Default `:ψ̃`.
- `state_type::Symbol`: `:ket` or `:density`. Default `:ket`.
- `subspace::AbstractVector{Int}`: Subspace indices for multi-level systems. Default `1:2`.
- `kwargs...`: Forwarded to `plot_bloch`.

# Returns
- A Makie `Figure` containing the animation.

# Examples

```julia
using CairoMakie

# Save Bloch animation to file
animate_bloch(traj_ket; mode=:record, filename="bloch.mp4", fps=24)
```

See also: `plot_bloch`, `animate_figure`, `animate_wigner`.
"""
function animate_bloch end

"""
    plot_bloch!(fig::Figure, traj::NamedTrajectory, idx::Int; kwargs...) -> Figure

In-place update of the Bloch vector observable on `fig` for frame `idx`. Used internally
by `animate_bloch`; requires a figure previously produced by `plot_bloch(...; index=...)`.
"""
function plot_bloch! end

"""
    animate_wigner(
        traj::NamedTrajectory;
        mode::Symbol = :inline,
        fps::Int = 24,
        filename::String = "wigner_animation.mp4",
        state_name::Symbol = :ψ̃,
        state_type::Symbol = :ket,
        kwargs...
    ) -> Figure

Animate the time evolution of the Wigner function for a quantum state trajectory.

Creates an animation showing how the Wigner quasi-probability distribution evolves in
phase space over time.

Defined as a stub here; the implementation is provided by the `PiccoloQuantumToolboxExt`
extension.

# Arguments
- `traj::NamedTrajectory`: The trajectory containing quantum states to animate.

# Keyword Arguments
- `mode::Symbol`: `:inline` (requires `GLMakie`) or `:record`. Default `:inline`.
- `fps::Int`: Frames per second. Default `24`.
- `filename::String`: Output path when `mode = :record`. Default `"wigner_animation.mp4"`.
- `state_name::Symbol`: Component name for the quantum state. Default `:ψ̃`.
- `state_type::Symbol`: `:ket` or `:density`. Default `:ket`.
- `kwargs...`: Forwarded to `plot_wigner` (e.g., `xvec`, `yvec`, `colorbar`).

# Returns
- A Makie `Figure` containing the animation.

# Examples

```julia
using CairoMakie

animate_wigner(traj; mode=:record, filename="wigner.mp4", fps=24,
               xvec=-4:0.1:4, yvec=-4:0.1:4)
```

See also: `plot_wigner`, `animate_bloch`, `animate_figure`.
"""
function animate_wigner end

"""
    plot_wigner!(fig::Figure, traj::NamedTrajectory, idx::Int) -> Figure

In-place update of the Wigner heatmap on `fig` for frame `idx`. Used internally by
`animate_wigner`; requires a figure previously produced by `plot_wigner(traj, 1; ...)`.
"""
function plot_wigner! end

end
