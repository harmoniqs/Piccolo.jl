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

# Plotting
include("pplots/_pplots.jl")
@reexport using .PPlots

end
