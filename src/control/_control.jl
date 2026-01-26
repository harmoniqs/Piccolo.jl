module Control

using Reexport

# Core dependencies
using LinearAlgebra
using SparseArrays
using Random
using Distributions

# Piccolo foundation
using DirectTrajOpt
using NamedTrajectories
using TrajectoryIndexingUtils

# Quantum module (sibling)
using ..Quantum

# Dependencies
using ExponentialAction
using TestItems

include("piccolo_options.jl")
@reexport using .Options

include("quantum_control_problem.jl")
@reexport using .QuantumControlProblems

include("quantum_objectives.jl")
@reexport using .QuantumObjectives

include("quantum_constraints.jl")
@reexport using .QuantumConstraints

include("quantum_integrators.jl")
@reexport using .QuantumIntegrators

include("problem_templates/_problem_templates.jl")
@reexport using .ProblemTemplates

end
