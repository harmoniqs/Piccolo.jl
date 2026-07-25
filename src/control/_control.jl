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

include("options.jl")
@reexport using .Options

include("problems.jl")
@reexport using .QuantumControlProblems

include("objectives.jl")
@reexport using .QuantumObjectives

# After objectives.jl: `verify` needs the `*_fidelity_loss` functions defined there to score
# the optimizer's own collocation state.
include("verification.jl")
@reexport using .Verification

include("constraints.jl")
@reexport using .QuantumConstraints

include("integrators.jl")
@reexport using .QuantumIntegrators

include("display/_display.jl")
@reexport using .ProblemDisplay

include("templates/_problem_templates.jl")
@reexport using .ProblemTemplates

end
