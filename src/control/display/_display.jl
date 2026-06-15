module ProblemDisplay

using LinearAlgebra
using Printf
using UnicodePlots

using NamedTrajectories
using DirectTrajOpt

using ..QuantumControlProblems: QuantumControlProblem
using ..QuantumControlProblems
using ..QuantumObjectives
using ...Quantum
using ...Quantum:
    AbstractQuantumTrajectory,
    UnitaryTrajectory,
    KetTrajectory,
    DensityTrajectory,
    MultiKetTrajectory,
    SamplingTrajectory,
    state_name,
    state_names,
    drive_name,
    get_goal,
    fidelity,
    iso_to_ket,
    EmbeddedOperator,
    unembed

export inspect, ProblemInspection
export show_problem
export pulse_lineplot

include("inspect.jl")
include("show.jl")
include("plot.jl")

end  # module
