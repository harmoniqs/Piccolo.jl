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
using ...Quantum:
    SamplingTrajectory,
    MultiKetTrajectory,
    state_names,
    get_weights,
    AbstractSplinePulse,
    AbstractPulse,
    ZeroOrderPulse,
    LinearSplinePulse,
    CubicSplinePulse

using ExponentialAction
using LinearAlgebra
using SparseArrays
using TestItems

const ⊗ = kron

include("smooth_pulse_problem.jl")
include("bang_bang_pulse_problem.jl")
include("spline_pulse_problem.jl")
include("minimum_time_problem.jl")
include("sampling_problem.jl")

function apply_piccolo_options!(
    piccolo_options::PiccoloOptions,
    constraints::AbstractVector{<:AbstractConstraint},
    traj::NamedTrajectory;
    state_names::Union{Nothing,Symbol,AbstractVector{Symbol}} = nothing,
    state_leakage_indices::Union{
        Nothing,
        AbstractVector{Int},
        AbstractVector{<:AbstractVector{Int}},
    } = nothing,
)
    J = NullObjective(traj)

    if piccolo_options.leakage_constraint
        val = piccolo_options.leakage_constraint_value
        if piccolo_options.verbose
            println("\tapplying leakage suppression: $(state_names) < $(val)")
        end

        if isnothing(state_leakage_indices)
            throw(ArgumentError("Leakage indices are required for leakage suppression."))
        end

        if state_names isa Symbol
            state_names = [state_names]
            state_leakage_indices = [state_leakage_indices]
        end

        for (name, indices) ∈ zip(state_names, state_leakage_indices)
            J += LeakageObjective(
                indices,
                name,
                traj,
                Qs = fill(piccolo_options.leakage_cost, traj.N),
            )
            push!(constraints, LeakageConstraint(val, indices, name, traj))
        end
    end

    if piccolo_options.timesteps_all_equal
        if piccolo_options.verbose
            println("\tapplying timesteps_all_equal constraint: $(traj.timestep)")
        end
        push!(constraints, TimeStepsAllEqualConstraint())
    end

    let cname = piccolo_options.complex_control_norm_constraint_name
        # Bind into a local so the !isnothing guard narrows the Union for JET.
        if cname !== nothing
            if piccolo_options.verbose
                println("\tapplying complex control norm constraint: $cname")
            end
            norm_con = NonlinearKnotPointConstraint(
                u -> [norm(u)^2 - piccolo_options.complex_control_norm_constraint_radius^2],
                cname,
                traj;
                equality = false,
            )
            push!(constraints, norm_con)
        end
    end

    return J
end

"""
    setup_free_phase_globals!(n_qubits, global_data, global_bounds; initial_phases=nothing, verbose=false)

Add per-qubit Z-phase variables (`φ_1`, `φ_2`, …) to `global_data` and
`global_bounds` dicts. Returns the vector of phase variable names and the
(possibly initialized) dicts.

Mutates `global_data` and `global_bounds` in place if they are not `nothing`;
otherwise creates new dicts.

# Keyword Arguments
- `initial_phases::Union{Nothing,Vector{Float64}}`: Initial values for the phase variables.
  If `nothing`, all phases are initialized to 0. If provided, must have length `n_qubits`.
- `verbose::Bool`: Print diagnostic information.
"""
function setup_free_phase_globals!(
    n_qubits::Int,
    global_data::Union{Nothing,Dict{Symbol,Vector{Float64}}},
    global_bounds::Union{Nothing,Dict{Symbol,<:Union{Float64,Tuple{Float64,Float64}}}};
    initial_phases::Union{Nothing,Vector{Float64}} = nothing,
    verbose::Bool = false,
)
    θ_names = [Symbol(:φ_, i) for i = 1:n_qubits]

    if !isnothing(initial_phases)
        @assert length(initial_phases) == n_qubits "initial_phases must have length $n_qubits"
    end

    if isnothing(global_data)
        global_data = Dict{Symbol,Vector{Float64}}()
    end
    for (i, name) in enumerate(θ_names)
        global_data[name] = [isnothing(initial_phases) ? 0.0 : initial_phases[i]]
    end

    if isnothing(global_bounds)
        global_bounds = Dict{Symbol,Union{Float64,Tuple{Float64,Float64}}}()
    end
    for name in θ_names
        if !haskey(global_bounds, name)
            global_bounds[name] = (-2π, 2π)
        end
    end

    if verbose
        println("    free_phase=true: added phase variables $θ_names")
    end

    return θ_names, global_data, global_bounds
end

"""
    _make_free_phase_goal(op::EmbeddedOperator)

Build a function `θ -> EmbeddedOperator` that applies single-qubit Z-phase rotations
to the goal gate. For an N-qubit gate, `θ` has N elements (one phase per qubit).

The phase-adjusted gate is `(Z(θ₁) ⊗ Z(θ₂) ⊗ ⋯) ⋅ U_goal`, where `Z(θ) = diag(1, e^{iθ})`.
"""
function _make_free_phase_goal(op::EmbeddedOperator)
    U_base = unembed(op)
    subspace = op.subspace
    levels = op.subsystem_levels
    n_qubits = length(levels)
    n_sub = size(U_base, 1)

    # Type-generic for ForwardDiff compatibility (θ may contain Dual numbers)
    function U_goal_fn(θ)
        phase_diag = map(1:n_sub) do i
            bits = i - 1
            phase = sum(
                θ[j] for j = 1:n_qubits if (bits >> (n_qubits - j)) & 1 == 1;
                init = zero(eltype(θ)),
            )
            return exp(im * phase)
        end
        phased = Diagonal(phase_diag) * U_base
        return EmbeddedOperator(Matrix(phased), subspace, levels)
    end
    return U_goal_fn
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
    global_bounds::Union{Nothing,AbstractDict{Symbol,<:Any}},
    traj::NamedTrajectory;
    verbose::Bool = false,
)
    if isnothing(global_bounds)
        return
    end

    for pair in global_bounds
        name::Symbol = pair.first
        bounds = pair.second
        if !haskey(traj.global_components, name)
            error(
                "Global variable :$name not found in trajectory. Available: $(keys(traj.global_components))",
            )
        end
        global_dim = length(traj.global_components[name])
        # Convert bounds to format expected by GlobalBoundsConstraint
        bounds_value = if bounds isa Float64
            # Symmetric scalar bounds
            bounds
        elseif bounds isa Tuple{Float64,Float64}
            # Asymmetric scalar bounds -> convert to vector tuple
            (fill(bounds[1], global_dim), fill(bounds[2], global_dim))
        else
            # Already in correct format (Vector or Tuple of Vectors)
            bounds
        end
        push!(constraints, GlobalBoundsConstraint(name, bounds_value))
        if verbose
            println("    added GlobalBoundsConstraint for :$name with bounds $bounds_value")
        end
    end
end

"""
    apply_calibration_targets!(constraints, calibration_targets, traj; verbose=false)

For each global variable name in `calibration_targets`, pin it at its current
value in `traj.global_data` via `fix_global_variable!`. Removes any existing
`GlobalBoundsConstraint`/`GlobalEqualityConstraint` on the name to avoid MOI
conflicts (so calling sites can pass `global_bounds` and `calibration_targets`
together — the latter wins, which is the intended semantics).

This implements the "calibration target" semantics: a global declared here is
a knob managed by an external calibration step, not a free NLP variable. It is
held at nominal during the QCP solve so the optimizer cannot drift it as a slack
variable. Downstream tools (e.g. QILC) selectively unpin individual targets at
their own discretion.

Modifies `constraints` in place. No-op when `calibration_targets` is empty.
"""
function apply_calibration_targets!(
    constraints::AbstractVector{<:AbstractConstraint},
    calibration_targets::Vector{Symbol},
    traj::NamedTrajectory;
    verbose::Bool = false,
)
    isempty(calibration_targets) && return
    for name in calibration_targets
        if !haskey(traj.global_components, name)
            error(
                "calibration_targets entry :$name not found in trajectory globals. " *
                "Available: $(keys(traj.global_components))",
            )
        end
        value = collect(Float64, traj.global_data[traj.global_components[name]])
        fix_global_variable!(constraints, name, value)
        if verbose
            println("    pinned :$name as calibration_target at value $value")
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
        timestep = :Δt,
        controls = :u,
        global_data = [0.1, 0.5, 0.3],  # flat vector
        global_components = (δ = 1:1, ω = 2:3),  # δ is scalar, ω is 2D
    )

    # Test 1: nothing global_bounds is a no-op
    constraints1 = AbstractConstraint[]
    add_global_bounds_constraints!(constraints1, nothing, traj)
    @test isempty(constraints1)

    # Test 2: Float64 symmetric scalar bounds
    constraints2 = AbstractConstraint[]
    add_global_bounds_constraints!(constraints2, Dict(:δ => 0.5), traj)
    @test length(constraints2) == 1
    @test constraints2[1] isa BoundsConstraint
    @test constraints2[1].is_global

    # Test 3: Tuple{Float64, Float64} asymmetric scalar bounds (expanded to vectors)
    constraints3 = AbstractConstraint[]
    add_global_bounds_constraints!(constraints3, Dict(:ω => (-0.2, 0.8)), traj)
    @test length(constraints3) == 1
    @test constraints3[1] isa BoundsConstraint
    @test constraints3[1].is_global

    # Test 4: Multiple globals with mixed bound types
    constraints4 = AbstractConstraint[]
    global_bounds = Dict{Symbol,Union{Float64,Tuple{Float64,Float64}}}(
        :δ => 0.5,           # symmetric
        :ω => (-0.2, 0.8),    # asymmetric
    )
    add_global_bounds_constraints!(constraints4, global_bounds, traj)
    @test length(constraints4) == 2
    @test all(c -> c isa BoundsConstraint && c.is_global, constraints4)

    # Test 5: Error when global variable doesn't exist
    constraints5 = AbstractConstraint[]
    @test_throws "Global variable :nonexistent not found" begin
        add_global_bounds_constraints!(constraints5, Dict(:nonexistent => 0.5), traj)
    end

    # Test 6: Verbose output (just ensure it doesn't error)
    constraints6 = AbstractConstraint[]
    add_global_bounds_constraints!(constraints6, Dict(:δ => 0.5), traj; verbose = true)
    @test length(constraints6) == 1
end

@testitem "apply_calibration_targets! pins globals at nominal" begin
    using NamedTrajectories
    using DirectTrajOpt
    using Piccolo.Control.ProblemTemplates:
        apply_calibration_targets!, add_global_bounds_constraints!

    N = 5
    traj = NamedTrajectory(
        (x = rand(2, N), u = rand(1, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
        global_data = [0.7, 1.2, 0.3],
        global_components = (a = 1:1, b = 2:2, c = 3:3),
    )

    # Test 1: empty → no-op
    constraints = AbstractConstraint[]
    apply_calibration_targets!(constraints, Symbol[], traj)
    @test isempty(constraints)

    # Test 2: pin :b at its nominal value 1.2
    apply_calibration_targets!(constraints, [:b], traj)
    eq_cons = filter(
        c -> c isa EqualityConstraint && c.var_names == :b && c.is_global,
        constraints,
    )
    @test length(eq_cons) == 1
    @test eq_cons[1].values ≈ [1.2]

    # Test 3: replaces an existing GlobalBoundsConstraint on the same name
    constraints2 = AbstractConstraint[]
    add_global_bounds_constraints!(constraints2, Dict(:b => (-2.0, 2.0)), traj)
    @test any(c -> c isa BoundsConstraint && c.var_names == :b && c.is_global, constraints2)
    apply_calibration_targets!(constraints2, [:b], traj)
    @test !any(
        c -> c isa BoundsConstraint && c.var_names == :b && c.is_global,
        constraints2,
    )
    @test count(
        c -> c isa EqualityConstraint && c.var_names == :b && c.is_global,
        constraints2,
    ) == 1

    # Test 4: error when name is not a global
    constraints3 = AbstractConstraint[]
    @test_throws ErrorException apply_calibration_targets!(
        constraints3,
        [:nonexistent],
        traj,
    )

    # Test 5: pin multiple targets at once
    constraints4 = AbstractConstraint[]
    apply_calibration_targets!(constraints4, [:a, :c], traj)
    eq_a = filter(
        c -> c isa EqualityConstraint && c.var_names == :a && c.is_global,
        constraints4,
    )
    eq_c = filter(
        c -> c isa EqualityConstraint && c.var_names == :c && c.is_global,
        constraints4,
    )
    @test length(eq_a) == 1 && eq_a[1].values ≈ [0.7]
    @test length(eq_c) == 1 && eq_c[1].values ≈ [0.3]
end

@testitem "apply_piccolo_options! throws ArgumentError for missing leakage indices" begin
    using NamedTrajectories
    using DirectTrajOpt
    using .ProblemTemplates: apply_piccolo_options!

    N = 5
    traj = NamedTrajectory(
        (x = rand(2, N), u = rand(1, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
    )

    piccolo_opts = PiccoloOptions(leakage_constraint = true)
    constraints = AbstractConstraint[]

    # Should throw ArgumentError (not ValueError) when leakage indices are missing
    @test_throws ArgumentError apply_piccolo_options!(
        piccolo_opts,
        constraints,
        traj;
        state_names = :x,
        state_leakage_indices = nothing,
    )
end

end
