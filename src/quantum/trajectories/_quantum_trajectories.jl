module QuantumTrajectories

using LinearAlgebra
using OrdinaryDiffEqLinear
using OrdinaryDiffEqTsit5
using TestItems
using ForwardDiff

using ..QuantumSystems: AbstractQuantumSystem, QuantumSystem, OpenQuantumSystem
using ..Pulses:
    AbstractPulse,
    AbstractSplinePulse,
    ZeroOrderPulse,
    LinearSplinePulse,
    CubicSplinePulse,
    GaussianPulse,
    ErfPulse,
    CompositePulse
using ..Pulses: duration, drive_name, n_drives
using ..Pulses: get_knot_times, get_knot_count, get_knot_values, get_knot_derivatives
import ..Pulses: duration, drive_name
import ..Rollouts
import ..Rollouts: rollout!
using ..Rollouts:
    UnitaryODEProblem,
    UnitaryOperatorODEProblem,
    KetODEProblem,
    KetOperatorODEProblem,
    DensityODEProblem
using ..Rollouts: unitary_fidelity
using ..EmbeddedOperators: AbstractPiccoloOperator, EmbeddedOperator
using ..Isomorphisms: operator_to_iso_vec, ket_to_iso, iso_to_ket, iso_vec_to_operator

import NamedTrajectories: NamedTrajectory, get_times

# Abstract type and common interface
include("abstract_trajectory.jl")

# Concrete trajectory types
include("unitary_trajectory.jl")
include("ket_trajectory.jl")
include("ensemble_trajectory.jl")
include("density_trajectory.jl")
include("sampling_trajectory.jl")

# Interface methods (getters, accessors, fidelity)
include("interface.jl")

# Extract pulse from optimized trajectories
include("extract_pulse.jl")

# NamedTrajectory conversion
include("named_trajectory_conversion.jl")

end  # module
