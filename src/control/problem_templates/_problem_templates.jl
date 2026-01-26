module ProblemTemplates

using ..QuantumObjectives
using ..QuantumConstraints
using ..QuantumIntegrators
using ..QuantumControlProblems: QuantumControlProblem, get_trajectory, get_system
using ..Options

using TrajectoryIndexingUtils
using NamedTrajectories
using DirectTrajOpt
using ...Quantum
using ...Quantum: SamplingTrajectory, MultiKetTrajectory, 
    state_names, get_weights, AbstractSplinePulse, AbstractPulse, ZeroOrderPulse,
    LinearSplinePulse, CubicSplinePulse

using ExponentialAction
using LinearAlgebra
using SparseArrays
using TestItems

const ⊗ = kron

include("smooth_pulse_problem.jl")
include("spline_pulse_problem.jl")
include("minimum_time_problem.jl")
include("sampling_problem.jl")

function apply_piccolo_options!(
    piccolo_options::PiccoloOptions,
    constraints::AbstractVector{<:AbstractConstraint},
    traj::NamedTrajectory;
    state_names::Union{Nothing,Symbol,AbstractVector{Symbol}}=nothing,
    state_leakage_indices::Union{Nothing,AbstractVector{Int},AbstractVector{<:AbstractVector{Int}}}=nothing,
)
    J = NullObjective(traj)

    if piccolo_options.leakage_constraint
        val = piccolo_options.leakage_constraint_value
        if piccolo_options.verbose
            println("\tapplying leakage suppression: $(state_names) < $(val)")
        end

        if isnothing(state_leakage_indices)
            throw(ValueError("Leakage indices are required for leakage suppression."))
        end

        if state_names isa Symbol
            state_names = [state_names]
            state_leakage_indices = [state_leakage_indices]
        end

        for (name, indices) ∈ zip(state_names, state_leakage_indices)
            J += LeakageObjective(indices, name, traj, Qs=fill(piccolo_options.leakage_cost, traj.N))
            push!(constraints, LeakageConstraint(val, indices, name, traj))
        end
    end

    if piccolo_options.timesteps_all_equal
        if piccolo_options.verbose
            println("\tapplying timesteps_all_equal constraint: $(traj.timestep)")
        end
        push!(
            constraints,
            TimeStepsAllEqualConstraint()
        )
    end

    if !isnothing(piccolo_options.complex_control_norm_constraint_name)
        if piccolo_options.verbose
            println("\tapplying complex control norm constraint: $(piccolo_options.complex_control_norm_constraint_name)")
        end
        norm_con = NonlinearKnotPointConstraint(
            u -> [norm(u)^2 - piccolo_options.complex_control_norm_constraint_radius^2],
            piccolo_options.complex_control_norm_constraint_name,
            traj;
            equality=false,
        )
        push!(constraints, norm_con)
    end

    return J
end

"""
    add_global_bounds_constraints!(constraints, global_bounds, traj; verbose=false)

Add GlobalBoundsConstraint entries for each global variable specified in `global_bounds`.

Converts bounds from user-friendly formats to the format expected by GlobalBoundsConstraint:
- `Float64`: Symmetric scalar bounds (applied symmetrically to all dimensions)
- `Tuple{Float64, Float64}`: Asymmetric scalar bounds (expanded to vectors)
- `Vector` or `Tuple{Vector, Vector}`: Already in correct format (passed through)

Modifies `constraints` in place.
"""
function add_global_bounds_constraints!(
    constraints::AbstractVector{<:AbstractConstraint},
    global_bounds,
    traj::NamedTrajectory;
    verbose::Bool=false
)
    if isnothing(global_bounds)
        return
    end
    
    for (name, bounds) in global_bounds
        if !haskey(traj.global_components, name)
            error("Global variable :$name not found in trajectory. Available: $(keys(traj.global_components))")
        end
        global_dim = length(traj.global_components[name])
        # Convert bounds to format expected by GlobalBoundsConstraint
        if bounds isa Float64
            # Symmetric scalar bounds
            bounds_value = bounds
        elseif bounds isa Tuple{Float64, Float64}
            # Asymmetric scalar bounds -> convert to vector tuple
            bounds_value = (fill(bounds[1], global_dim), fill(bounds[2], global_dim))
        else
            # Already in correct format (Vector or Tuple of Vectors)
            bounds_value = bounds
        end
        push!(constraints, GlobalBoundsConstraint(name, bounds_value))
        if verbose
            println("    added GlobalBoundsConstraint for :$name with bounds $bounds_value")
        end
    end
end

@testitem "add_global_bounds_constraints! helper function" begin
    using NamedTrajectories
    using DirectTrajOpt
    using .ProblemTemplates: add_global_bounds_constraints!
    
    # Create a trajectory with global components for testing
    # global_data is a flat vector, global_components maps names to index ranges
    N = 5
    traj = NamedTrajectory(
        (x = rand(2, N), u = rand(1, N), Δt = fill(0.1, N));
        timestep=:Δt,
        controls=:u,
        global_data=[0.1, 0.5, 0.3],  # flat vector
        global_components=(δ = 1:1, ω = 2:3)  # δ is scalar, ω is 2D
    )
    
    # Test 1: nothing global_bounds is a no-op
    constraints1 = AbstractConstraint[]
    add_global_bounds_constraints!(
        constraints1, nothing, traj
    )
    @test isempty(constraints1)
    
    # Test 2: Float64 symmetric scalar bounds
    constraints2 = AbstractConstraint[]
    add_global_bounds_constraints!(
        constraints2, Dict(:δ => 0.5), traj
    )
    @test length(constraints2) == 1
    @test constraints2[1] isa BoundsConstraint
    @test constraints2[1].is_global
    
    # Test 3: Tuple{Float64, Float64} asymmetric scalar bounds (expanded to vectors)
    constraints3 = AbstractConstraint[]
    add_global_bounds_constraints!(
        constraints3, Dict(:ω => (-0.2, 0.8)), traj
    )
    @test length(constraints3) == 1
    @test constraints3[1] isa BoundsConstraint
    @test constraints3[1].is_global
    
    # Test 4: Multiple globals with mixed bound types
    constraints4 = AbstractConstraint[]
    global_bounds = Dict{Symbol, Union{Float64, Tuple{Float64, Float64}}}(
        :δ => 0.5,           # symmetric
        :ω => (-0.2, 0.8)    # asymmetric
    )
    add_global_bounds_constraints!(
        constraints4, global_bounds, traj
    )
    @test length(constraints4) == 2
    @test all(c -> c isa BoundsConstraint && c.is_global, constraints4)
    
    # Test 5: Error when global variable doesn't exist
    constraints5 = AbstractConstraint[]
    @test_throws "Global variable :nonexistent not found" begin
        add_global_bounds_constraints!(
            constraints5, Dict(:nonexistent => 0.5), traj
        )
    end
    
    # Test 6: Verbose output (just ensure it doesn't error)
    constraints6 = AbstractConstraint[]
    add_global_bounds_constraints!(
        constraints6, Dict(:δ => 0.5), traj; verbose=true
    )
    @test length(constraints6) == 1
end

end
