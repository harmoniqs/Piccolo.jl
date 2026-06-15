module ProblemTemplates

using ..QuantumObjectives
using ..QuantumConstraints
using ..QuantumIntegrators
using ..QuantumControlProblems: QuantumControlProblem, get_trajectory, get_system
using ..Options
using ..ProblemDisplay: show_problem

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

# ============================================================================ #
# Verbose-output helpers
# ============================================================================ #

# Display-level gates. `:silent` shows nothing during construction.
# `:compact` shows the outer "constructing X" header line.
# `:standard` (default) additionally triggers `display(qcp)` at the end of the
# outer constructor (handled by template return path).
# `:detailed` additionally shows inner template details + plot.
_show_header(opts::PiccoloOptions) = display_level(opts.display) >= DISPLAY_COMPACT
_show_details(opts::PiccoloOptions) = display_level(opts.display) >= DISPLAY_DETAILED
_should_inspect(opts::PiccoloOptions) = display_level(opts.display) >= DISPLAY_STANDARD

"""
    _maybe_display(qcp, opts) -> qcp

Auto-render `qcp` after outer-constructor return when `opts.display` is
`:standard` or `:detailed`. `:standard` shows the rich tree view;
`:detailed` adds sparsity + a terminal pulse plot. Returns `qcp` so it can
be used as `return _maybe_display(qcp, opts)`.
"""
function _maybe_display(qcp, opts::PiccoloOptions)
    lvl = display_level(opts.display)
    if lvl >= DISPLAY_STANDARD
        detail = lvl >= DISPLAY_DETAILED ? :full : :standard
        show_problem(stdout, qcp; detail = detail)
        println()
    end
    return qcp
end

# Strip type parameters for human-readable logs:
#   KetTrajectory{CubicSplinePulse{...}}  → "KetTrajectory"
#   CubicSplinePulse{DataInterpolations…} → "CubicSplinePulse"
_typename(T::Type) = string(nameof(T))
_typename(x) = _typename(typeof(x))

# Compact bound formatting for verbose construction logs. Handles the four
# concrete shapes used by GlobalBoundsConstraint / update_bound!:
#   Float64               → "±0.0195"
#   (lo::Float64, hi)     → "±2π" when symmetric, else "[lo, hi]"
#   Vector{Float64}       → "±[…]"
#   (lo::Vec, hi::Vec)    → "±[…]" when symmetric, else "[lo, hi]"
function _fmt_bounds(b)
    if b isa Real
        return string("±", round(b; sigdigits = 4))
    elseif b isa Tuple{<:Real,<:Real}
        lo, hi = b
        return lo == -hi ? string("±", round(hi; sigdigits = 4)) :
               string("[", round(lo; sigdigits = 4), ", ", round(hi; sigdigits = 4), "]")
    elseif b isa AbstractVector{<:Real}
        return string("±", round.(b; sigdigits = 4))
    elseif b isa Tuple{<:AbstractVector,<:AbstractVector}
        lo, hi = b
        return all(lo .== .-hi) ? string("±", round.(hi; sigdigits = 4)) :
               string("[", round.(lo; sigdigits = 4), ", ", round.(hi; sigdigits = 4), "]")
    else
        return string(b)
    end
end

include("smooth_pulse_problem.jl")
include("bang_bang_pulse_problem.jl")
include("spline_pulse_problem.jl")
include("minimum_time_problem.jl")
include("sampling_problem.jl")

"""
    _unbind_state!(traj::NamedTrajectory, name::Symbol)

Widen the bounds for variable `name` to `(-Inf, Inf)`, effectively removing the
box constraint while keeping the NamedTuple type stable (cannot delete keys from
a parametric NamedTuple). No-op if `name` is not in `traj.bounds`.
"""
function _unbind_state!(traj::NamedTrajectory, name::Symbol)
    name ∈ keys(traj.bounds) || return nothing
    d = traj.dims[name]
    update_bound!(traj, name, (-fill(Inf, d), fill(Inf, d)))
    return nothing
end

"""
    _safe_bound_times(name::Symbol, traj::NamedTrajectory) -> Vector{Int}

Compute time indices where it is safe to add bounds on variable `name`,
excluding timesteps already pinned by initial or final equality constraints.
Follows the same pattern as `get_trajectory_constraints` in DirectTrajOpt
(`problems.jl:174-187`).
"""
function _safe_bound_times(name::Symbol, traj::NamedTrajectory)
    has_initial = name ∈ keys(traj.initial)
    has_final = name ∈ keys(traj.final)
    if has_initial && has_final
        return collect(2:(traj.N-1))
    elseif has_initial
        return collect(2:traj.N)
    elseif has_final
        return collect(1:(traj.N-1))
    else
        return collect(1:traj.N)
    end
end

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
    iso_layout::Symbol = :block,
)
    J = NullObjective(traj)

    if piccolo_options.leakage_constraint
        val = piccolo_options.leakage_constraint_value
        if _show_details(piccolo_options)
            println("    applying leakage suppression: $(state_names) < $(val)")
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
        if _show_details(piccolo_options)
            println("    applying timesteps_all_equal constraint: $(traj.timestep)")
        end
        push!(constraints, TimeStepsAllEqualConstraint())
    end

    let cname = piccolo_options.complex_control_norm_constraint_name
        # Bind into a local so the !isnothing guard narrows the Union for JET.
        if cname !== nothing
            if _show_details(piccolo_options)
                println("    applying complex control norm constraint: $cname")
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

    if !piccolo_options.bound_state
        # Widen default [-1, 1] state bounds to (-Inf, Inf), effectively removing
        # the box constraint. The NamedTuple type is preserved (cannot delete keys).
        _names = if state_names isa Symbol
            [state_names]
        elseif isnothing(state_names)
            Symbol[]
        else
            state_names
        end
        for name in _names
            if name ∈ keys(traj.bounds)
                if _show_details(piccolo_options)
                    println("    unbinding state :$name (bound_state=false)")
                end
                _unbind_state!(traj, name)
            end
        end
    end

    if piccolo_options.bound_state_l2
        if isnothing(state_names)
            throw(ArgumentError("state_names required for bound_state_l2 constraint."))
        end
        if _show_details(piccolo_options)
            println("    applying bound_state_l2 constraint: $(state_names), |z|² ≤ 1")
        end
        _names = state_names isa Symbol ? [state_names] : state_names
        for name in _names
            ts = _safe_bound_times(name, traj)
            isempty(ts) && continue
            push!(constraints, BoundStateL2Constraint(name, traj, iso_layout; times = ts))
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
        _θ_list = join(θ_names, ", ")
        println("    free_phase: added $_θ_list")
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
            println("    global bound :$name = $(_fmt_bounds(bounds_value))")
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
            println(
                "    pinned :$name as calibration_target = $(round.(value; sigdigits=4))",
            )
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

@testitem "bound_state=true preserves state bounds" begin
    using NamedTrajectories
    using DirectTrajOpt

    apply_piccolo_options! = Piccolo.Control.ProblemTemplates.apply_piccolo_options!

    N = 5
    traj = NamedTrajectory(
        (ψ̃ = rand(4, N), u = rand(1, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
        bounds = (ψ̃ = (-ones(4), ones(4)), u = 1.0),
    )

    piccolo_opts = PiccoloOptions(bound_state = true)
    constraints = AbstractConstraint[]

    apply_piccolo_options!(piccolo_opts, constraints, traj; state_names = :ψ̃)

    # bound_state=true (default) keeps existing traj.bounds on ψ̃
    @test :ψ̃ ∈ keys(traj.bounds)
    # No extra BoundsConstraint added to the constraints vector
    bc = filter(c -> c isa BoundsConstraint, constraints)
    @test isempty(bc)
end

@testitem "bound_state=false widens state bounds to ±Inf" begin
    using NamedTrajectories
    using DirectTrajOpt

    apply_piccolo_options! = Piccolo.Control.ProblemTemplates.apply_piccolo_options!

    N = 5
    traj = NamedTrajectory(
        (ψ̃ = rand(4, N), u = rand(1, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
        bounds = (ψ̃ = (-ones(4), ones(4)), u = 1.0),
    )

    piccolo_opts = PiccoloOptions(bound_state = false)
    constraints = AbstractConstraint[]

    apply_piccolo_options!(piccolo_opts, constraints, traj; state_names = :ψ̃)

    # State bounds widened to ±Inf (effectively no constraint)
    lb, ub = traj.bounds[:ψ̃]
    @test all(lb .== -Inf)
    @test all(ub .== Inf)
    # Control bounds preserved
    @test all(abs.(traj.bounds[:u][1]) .< Inf)
end

@testitem "bound_state=false with multiple state names" begin
    using NamedTrajectories
    using DirectTrajOpt

    apply_piccolo_options! = Piccolo.Control.ProblemTemplates.apply_piccolo_options!

    N = 5
    traj = NamedTrajectory(
        (ψ̃1 = rand(4, N), ψ̃2 = rand(4, N), u = rand(1, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
        bounds = (ψ̃1 = 1.0, ψ̃2 = 1.0, u = 1.0),
    )

    piccolo_opts = PiccoloOptions(bound_state = false)
    constraints = AbstractConstraint[]

    apply_piccolo_options!(piccolo_opts, constraints, traj; state_names = [:ψ̃1, :ψ̃2])

    lb1, ub1 = traj.bounds[:ψ̃1]
    @test all(lb1 .== -Inf) && all(ub1 .== Inf)
    lb2, ub2 = traj.bounds[:ψ̃2]
    @test all(lb2 .== -Inf) && all(ub2 .== Inf)
    # Control bounds preserved
    @test all(abs.(traj.bounds[:u][1]) .< Inf)
end

@testitem "bound_state_l2 adds NonlinearKnotPointConstraint (block layout)" begin
    using NamedTrajectories
    using DirectTrajOpt

    apply_piccolo_options! = Piccolo.Control.ProblemTemplates.apply_piccolo_options!

    N = 5
    # 4-dim iso-vec = 2 complex components
    traj = NamedTrajectory(
        (ψ̃ = rand(4, N), u = rand(1, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
    )

    piccolo_opts = PiccoloOptions(bound_state_l2 = true)
    constraints = AbstractConstraint[]

    apply_piccolo_options!(
        piccolo_opts,
        constraints,
        traj;
        state_names = :ψ̃,
        iso_layout = :block,
    )

    nlc = filter(c -> c isa AbstractNonlinearConstraint, constraints)
    @test length(nlc) == 1
    @test nlc[1].equality == false
    # 2 complex components → g_dim = 2 per knot point
    @test nlc[1].g_dim == 2
end

@testitem "bound_state_l2 throws when state_names is nothing" begin
    using NamedTrajectories
    using DirectTrajOpt

    apply_piccolo_options! = Piccolo.Control.ProblemTemplates.apply_piccolo_options!

    N = 5
    traj = NamedTrajectory(
        (x = rand(2, N), u = rand(1, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
    )

    piccolo_opts = PiccoloOptions(bound_state_l2 = true)
    constraints = AbstractConstraint[]

    @test_throws ArgumentError apply_piccolo_options!(piccolo_opts, constraints, traj;)
end

@testitem "_safe_bound_times respects initial/final constraints" begin
    using NamedTrajectories

    _safe_bound_times = Piccolo.Control.ProblemTemplates._safe_bound_times

    N = 5
    # No initial/final → all times
    traj1 = NamedTrajectory(
        (x = rand(2, N), u = rand(1, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
    )
    @test _safe_bound_times(:x, traj1) == collect(1:N)

    # Initial only → skip time 1
    traj2 = NamedTrajectory(
        (x = rand(2, N), u = rand(1, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
        initial = (x = rand(2),),
    )
    @test _safe_bound_times(:x, traj2) == collect(2:N)

    # Both initial and final → skip endpoints
    traj3 = NamedTrajectory(
        (x = rand(2, N), u = rand(1, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
        initial = (x = rand(2),),
        final = (x = rand(2),),
    )
    @test _safe_bound_times(:x, traj3) == collect(2:(N-1))
end

end
