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

# `duration` must resolve to the Pulses function in this namespace: with
# NamedTrajectories >= 0.9.3 co-resolved (`using NamedTrajectories` above),
# TimeWarp also exports `duration`, which left the name ambiguous here —
# UndefVarError at `Quantum.duration` and at downstream `import ..Quantum:
# duration` sites (#323). Bind explicitly from the defining module.
using .Pulses: duration

# Object utils (depends on gates)
include("object_utils.jl")
@reexport using .QuantumObjectUtils

# Operators - lifted (no dependencies)
include("operators/lifted_operators.jl")
@reexport using .LiftedOperators

# Systems
include("systems/_quantum_systems.jl")
# Shared operator seam (slice 3b, #430): abstract dynamics-operator layer +
# MatrixOperator bridge moved from Piccolissimo — the dense spline cells and
# the interval-coefficient kernel build on them; structured operators extend
# them from Piccolissimo.
include("operators/abstract_dynamics_operator.jl")
include("operators/matrix_operator.jl")

@reexport using .QuantumSystems

# Encodings (depend on gates; used by embedded operators)
include("encodings/dual_rail.jl")
@reexport using .DualRailEncodings

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
