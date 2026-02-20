module QuantumObjectPlots

export plot_unitary_populations
export plot_state_populations
export plot_weyl_trajectory
export plot_name!

using LaTeXStrings
using Makie
import Makie: plot
using LinearAlgebra
using NamedTrajectories
using Piccolo
using TestItems

# Utility functions
function get_layout(index::Int, n::Int)
    √n = isqrt(n) + 1
    return ((index - 1) ÷ √n + 1, (index - 1) % √n + 1)
end

# Stub for plot_name! - to be implemented or delegated to NamedTrajectories
function plot_name!(ax, traj::NamedTrajectory, name::Symbol; indices = 1:traj.N, kwargs...)
    times = get_times(traj)[indices]
    data = traj[name][:, indices]
    for (i, row) in enumerate(eachrow(data))
        lines!(ax, times, collect(row); kwargs...)
    end
end

# Include submodules
include("unitary_populations.jl")
include("state_populations.jl")
include("weyl_trajectory.jl")

end
