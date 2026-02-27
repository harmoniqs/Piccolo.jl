module Quantum

using Reexport

# Core dependencies
using LinearAlgebra
using SparseArrays
using ForwardDiff
using NamedTrajectories

# ODE / integration dependencies
using DataInterpolations
using ExponentialAction
using OrdinaryDiffEqLinear
using OrdinaryDiffEqTsit5
using SciMLBase
using SymbolicIndexingInterface

# Other dependencies
using ProgressMeter
using SpecialFunctions
using TestItems

# Primitives
include("primitives/gates.jl")
@reexport using .Gates

include("primitives/isomorphisms.jl")
@reexport using .Isomorphisms

include("primitives/pulses.jl")
@reexport using .Pulses

# Object utils (depends on gates)
include("object_utils.jl")
@reexport using .QuantumObjectUtils

# Operators - lifted (no dependencies)
include("operators/lifted_operators.jl")
@reexport using .LiftedOperators

# Systems
include("systems/_quantum_systems.jl")
@reexport using .QuantumSystems

# Operators - embedded and direct_sums (depend on systems)
include("operators/embedded_operators.jl")
@reexport using .EmbeddedOperators

include("operators/direct_sums.jl")
@reexport using .DirectSums

# System utils (depends on embedded_operators, systems, etc.)
include("system_utils.jl")
@reexport using .QuantumSystemUtils

# Dynamics
include("dynamics.jl")
@reexport using .Rollouts

# Trajectories
include("trajectories/_quantum_trajectories.jl")
@reexport using .QuantumTrajectories

# Templates
include("templates/_quantum_system_templates.jl")
@reexport using .QuantumSystemTemplates

end
