module PPlots

using Reexport

include("animated_plots.jl")
@reexport using .AnimatedPlots

include("quantum_object_plots/_quantum_object_plots.jl")
@reexport using .QuantumObjectPlots

include("quantum_toolbox_plots.jl")
@reexport using .QuantumToolboxPlots

end
