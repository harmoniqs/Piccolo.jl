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

# ─────────────────────────────────────────────────────────────────────────── #
# Adjoint-robustness objective family (open-core slice 3d, #432).              #
#                                                                             #
# Moved from Piccolissimo (moved-file manifest row 28): the published adjoint  #
# (variational) robustness objectives, the propagation methods, the GN         #
# HVP/VJP/JVP primitives defined on the objectives, the RobustControlProblem   #
# wrapper, and their colocated testitems. Included AFTER the integrator        #
# modules (the family's spline-ODE propagation and exponential-integrator      #
# constructors bind `QuantumIntegrators` names); the re-export block below is  #
# the top-level seam, kept complete by the drift guard in                       #
# test/test_robustness_reexport_seam.jl.                                        #
#                                                                             #
# What STAYED in Piccolissimo: the solver-backend HVP machinery (Altissimo      #
# backend paths, knot_hvp carriers) — it consumes this surface through the      #
# re-export seam exactly as an external consumer would, so the dependency       #
# arrow ends Piccolissimo → Piccolo and never the reverse.                      #
# ─────────────────────────────────────────────────────────────────────────── #
include("objectives_robustness/_adjoint_robustness.jl")
using .AdjointRobustness

# The submodule name itself rides the seam (the SplineConstraints convention):
# Piccolissimo's seam and the drift guard reach the module through this export.
# DirectTrajOpt-owned generics the family extends (objective_value, gradient!,
# hessian_structure, get_full_hessian) are deliberately NOT re-exported here:
# interface functions come from DirectTrajOpt, which `using Piccolo` already
# re-exports.
export AdjointRobustness
export AdjointRobustnessObjective,
    KetAdjointRobustnessObjective,
    RobustControlProblem,
    ZKickPropagation,
    test_robustness_objective,
    objective_hvp!

include("display/_display.jl")
@reexport using .ProblemDisplay

include("templates/_problem_templates.jl")
@reexport using .ProblemTemplates

end
