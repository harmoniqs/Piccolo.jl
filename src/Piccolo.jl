module Piccolo

using Reexport

# Foundation packages (stay separate - they're generic)
@reexport using TrajectoryIndexingUtils
@reexport using NamedTrajectories
@reexport using DirectTrajOpt

# Quantum objects: systems, gates, pulses, trajectories
include("quantum/_quantum.jl")
@reexport using .Quantum

# Optimal control: objectives, constraints, problem templates
include("control/_control.jl")
@reexport using .Control

# Stub functions for plotting extensions
# These are implemented in ext/PiccoloMakieExt.jl when Makie is loaded
export animate_figure, animate_name

"""
    animate_figure(fig, frames, update_frame!; mode=:inline, fps=24, filename="animation.mp4")

Animate a Makie figure. Requires CairoMakie or GLMakie to be loaded.
"""
function animate_figure end

"""
    animate_name(traj, name; fps=24, mode=:inline, filename="name_animation.mp4", kwargs...)

Animate the evolution of a variable in a NamedTrajectory. Requires Makie to be loaded.
"""
function animate_name end

end
