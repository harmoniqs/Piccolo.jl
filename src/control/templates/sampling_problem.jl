export SamplingProblem, SamplingParams

# ----------------------------------------------------------------------------- #
# The wrapper type
# ----------------------------------------------------------------------------- #
#
# Wrappers are hand-written rather than declared through a macro. `@problem_template`
# exists because a *template's* alias bound encodes the pulse x trajectory
# compatibility matrix that `Specs.emit_schema` must mirror -- one declaration, one
# matrix, no drift. A wrapper has no such matrix: it is polymorphic in the problem
# it wraps, so a `@problem_wrapper` macro would generate four lines and buy no
# anti-drift property. Wrappers therefore stay hand-written, with their registry
# entry reflected from `SamplingParams` through the same `reflect_params` the
# templates use (see `Specs.register_all!`).

@doc raw"""
    SamplingParams <: AbstractTemplateParams

The typed keyword surface of [`SamplingProblem`](@ref). `weights` is empty when the
variants are equally weighted.
"""
Base.@kwdef struct SamplingParams <: AbstractTemplateParams
    Q::Float64 = 100.0
    weights::Vector{Float64} = Float64[]
    calibration_targets::Vector{Symbol} = Symbol[]
end

@doc raw"""
    SamplingProblem{P<:AbstractQuantumControlProblem, QT<:SamplingTrajectory} <: AbstractProblemWrapper

Robust-by-ensemble wrapper: one shared control acting on several systems, each
evolving under its own dynamics, with a weighted sum of fidelity objectives.

The wrap history lives in the **type**:

```julia
SamplingProblem{QuantumControlProblem{SmoothPulseTemplate, UnitaryTrajectory{ZeroOrderPulse,…}}}
```

is a sampled *smooth* unitary problem, and says so — the type-level mirror of the
ordered `[[wrappers]]` list in a spec. Accessors forward through
[`inner`](@ref)/[`quantum_trajectory`](@ref); see [`AbstractProblemWrapper`](@ref)
for the field contract.
"""
mutable struct SamplingProblem{P<:AbstractQuantumControlProblem,QT<:SamplingTrajectory} <:
               AbstractProblemWrapper
    inner::P
    qtraj::QT
    prob::DirectTrajOptProblem
    params::SamplingParams
    spec::Union{Nothing,AbstractProblemSpec}
end

wrapper_kind(::SamplingProblem) = :sampling

rewrap(w::SamplingProblem, qtraj::SamplingTrajectory, prob::DirectTrajOptProblem) =
    SamplingProblem(inner(w), qtraj, prob, template_params(w), retained_spec(w))

# Note: SamplingTrajectory is now exported from PiccoloQuantumObjects

# ============================================================================= #
# Sampling Objective Utilities
# ============================================================================= #

"""
    extract_regularization(objective, state_sym::Symbol, new_traj::NamedTrajectory) -> AbstractObjective

Extract regularization terms (non-state-dependent objectives) from a composite objective,
filtering to only include terms for variables that exist in the new trajectory.

Used by `SamplingProblem` to extract shared regularizers (e.g., control penalty) from
the base problem while excluding regularizers for variables that don't exist in the
sampling trajectory (e.g., `:du`, `:ddu` which are added by `SmoothPulseProblem`).
"""
function extract_regularization(objective, state_sym::Symbol, new_traj::NamedTrajectory)
    objs = hasproperty(objective, :objectives) ? objective.objectives : [objective]

    regularizers = filter(objs) do term
        # Get variable names this term depends on
        term_syms = if hasproperty(term, :syms)
            term.syms
        elseif hasproperty(term, :var_names)
            term.var_names
        elseif hasproperty(term, :name) && term.name isa Symbol
            # QuadraticRegularizer has a single :name field
            (term.name,)
        else
            Symbol[]
        end
        # Only include if:
        # 1. It doesn't depend on the state symbol
        # 2. All its variables exist in the new trajectory
        state_sym ∉ term_syms && all(s -> s ∈ new_traj.names, term_syms)
    end

    return isempty(regularizers) ? NullObjective() : reduce(+, regularizers)
end

# ============================================================================= #
# Sampling State Objective (dispatch-based)
# ============================================================================= #

"""
    sampling_state_objective(qtraj, traj, state_sym, Q)

Create the state-dependent objective for a sampling member.
Dispatches on quantum trajectory type.
"""
function sampling_state_objective(
    qtraj::UnitaryTrajectory,
    traj::NamedTrajectory,
    state_sym::Symbol,
    Q::Float64,
)
    goal = get_goal(qtraj)
    # Use free-phase objective when globals exist and goal is EmbeddedOperator
    if traj.global_dim > 0 && goal isa EmbeddedOperator
        # `Symbol[traj.global_names...]` rather than `collect(traj.global_names)` so
        # the eltype is `Symbol` even when the source tuple is empty (JET infers
        # `Vector{Union{}}` from `collect` on a possibly-empty tuple, which then
        # fails to match `θ_names::AbstractVector{Symbol}` downstream).
        θ_names = Symbol[traj.global_names...]
        U_goal_fn = _make_free_phase_goal(goal)
        return UnitaryFreePhaseInfidelityObjective(
            U_goal_fn,
            state_sym,
            θ_names,
            traj;
            Q = Q,
        )
    else
        return UnitaryInfidelityObjective(goal, state_sym, traj; Q = Q)
    end
end

function sampling_state_objective(
    qtraj::KetTrajectory,
    traj::NamedTrajectory,
    state_sym::Symbol,
    Q::Float64,
)
    ψ_goal = get_goal(qtraj)
    return KetInfidelityObjective(ψ_goal, state_sym, traj; Q = Q)
end

function sampling_state_objective(
    qtraj::MultiKetTrajectory,
    traj::NamedTrajectory,
    state_syms::Vector{Symbol},
    Q::Float64,
)
    # Per-member coherent ensemble objective — reuses the MultiKetTrajectory
    # machinery: one coherent infidelity over the member's state transfers,
    # weighted by the base trajectory's per-state weights.
    return CoherentKetInfidelityObjective(
        get_goal(qtraj),
        state_syms,
        traj;
        Q = Q,
        weights = qtraj.weights,
    )
end

# Density bases: Piccolo ships no public density fidelity objective, so the
# density sampling-objective cells are intentionally loud — construction fails
# here with an error naming the extension point rather than silently installing
# a null objective that would optimize nothing. A downstream package (e.g.
# Piccolissimo) registers the density sampling objective by extending
# `sampling_state_objective` for the density trajectory types.
const _DENSITY_SAMPLING_OBJECTIVE_ERROR = """
SamplingProblem cannot build a state objective for a density-matrix base trajectory: \
Piccolo does not provide a public density fidelity objective. The sampling-trajectory \
conversion and the compact-Lindbladian bilinear integrators are available, but solving \
requires a density sampling objective from a downstream package.

Extension point: extend `sampling_state_objective` for the density trajectory types, e.g.

    Piccolo.ProblemTemplates.sampling_state_objective(
        ::DensityTrajectory, traj::NamedTrajectory, state_sym::Symbol, Q::Float64,
    )

(Piccolissimo registers its density sampling objective through this hook.) No null \
objective was installed; this loud error is intentional."""

# Generic fallback (no density-specific methods): keeps the density cells loud
# while letting a downstream package EXTEND `sampling_state_objective` with typed
# density methods at module top level — no method overwrite, no `__init__` eval.
function sampling_state_objective(
    qtraj::Union{DensityTrajectory,MultiDensityTrajectory},
    traj::NamedTrajectory,
    state_sym,
    Q::Float64,
)
    error(_DENSITY_SAMPLING_OBJECTIVE_ERROR)
end

# ============================================================================= #
# Integrator keyword resolution
# ============================================================================= #

"""
    _validate_sampling_integrator_vector(integrators, n_slots; source) -> Vector{AbstractIntegrator}

Validate a user-supplied integrator vector for a sampling problem: every element
must be an `AbstractIntegrator`, and there must be one integrator per ensemble
dynamics slot (one per member for single-state bases, one per (member, sub-state)
for multi-state bases). Errors are actionable and name the offending shape.
"""
function _validate_sampling_integrator_vector(
    integrators::AbstractVector,
    n_slots::Int;
    source::String,
)
    bad = findfirst(i -> !(i isa AbstractIntegrator), integrators)
    if !isnothing(bad)
        error(
            "SamplingProblem $source: every integrator must be an AbstractIntegrator, " *
            "but element $bad is a $(typeof(integrators[bad])). " *
            "Pass integrators built against the sampling trajectory, e.g. " *
            "`BilinearIntegrator(SamplingTrajectory(qtraj, systems), N)`.",
        )
    end
    if length(integrators) != n_slots
        error(
            "SamplingProblem $source: the ensemble has $n_slots dynamics slot(s) " *
            "(one per member, times sub-states for multi-state bases), but received " *
            "$(length(integrators)) integrator(s). Pass one integrator per slot, or a " *
            "factory `(sampling_qtraj, N) -> integrators`.",
        )
    end
    return AbstractIntegrator[integrators...]
end

"""
    _resolve_sampling_integrators(integrator, sampling_qtraj, N, n_slots) -> Vector{AbstractIntegrator}

Resolve the `integrator` keyword of `SamplingProblem` into the dynamics integrator
vector. Three call shapes, aligned with the other problem templates:

- `nothing` — default `BilinearIntegrator(sampling_qtraj, N)`.
- an `AbstractIntegrator` instance — valid only for single-slot ensembles.
- a vector of integrators — one per ensemble dynamics slot.
- a factory `Function` — called as `integrator(sampling_qtraj, N)`, returning an
  integrator or vector of integrators (validated like the shapes above).
"""
function _resolve_sampling_integrators(integrator, sampling_qtraj, N::Int, n_slots::Int)
    if isnothing(integrator)
        default_int = BilinearIntegrator(sampling_qtraj, N)
        return AbstractIntegrator[(default_int isa AbstractVector ? default_int :
                                   [default_int])...,]
    elseif integrator isa AbstractIntegrator
        if n_slots != 1
            error(
                "SamplingProblem integrator keyword: a single integrator instance was " *
                "provided, but the ensemble has $n_slots dynamics slots (one per member). " *
                "Pass a per-member vector of integrators or a factory " *
                "`(sampling_qtraj, N) -> integrators`.",
            )
        end
        return AbstractIntegrator[integrator]
    elseif integrator isa AbstractVector
        return _validate_sampling_integrator_vector(
            integrator,
            n_slots;
            source = "integrator keyword",
        )
    else
        # Factory: called as integrator(sampling_qtraj, N)
        result = integrator(sampling_qtraj, N)
        if result isa AbstractIntegrator
            if n_slots != 1
                error(
                    "SamplingProblem integrator factory returned a single integrator, " *
                    "but the ensemble has $n_slots dynamics slots (one per member). " *
                    "Return one integrator per slot, e.g. " *
                    "`(sq, n) -> BilinearIntegrator(sq, n)`.",
                )
            end
            return AbstractIntegrator[result]
        elseif result isa AbstractVector
            return _validate_sampling_integrator_vector(
                result,
                n_slots;
                source = "integrator factory",
            )
        else
            error(
                "SamplingProblem integrator factory must return an AbstractIntegrator " *
                "or a vector of integrators; got $(typeof(result)).",
            )
        end
    end
end

# ============================================================================= #
# SamplingProblem Constructor
# ============================================================================= #

@doc raw"""
    SamplingProblem(qcp::AbstractQuantumControlProblem, systems::Vector{<:AbstractQuantumSystem}; kwargs...)

Construct a `SamplingProblem` from an existing `QuantumControlProblem` and a list of systems.

This creates a robust optimization problem where the controls are shared across all systems,
but each system evolves according to its own dynamics. The objective is the weighted sum of
fidelity objectives for each system.

# Arguments
- `qcp::AbstractQuantumControlProblem`: The base problem (defines nominal trajectory, objective, etc.)
- `systems::Vector{<:AbstractQuantumSystem}`: List of systems to optimize over

# Keyword Arguments
- `weights::Vector{Float64}=fill(1.0, length(systems))`: Weights for each system
- `Q::Float64=100.0`: Weight on infidelity objective (explicit, not extracted from base problem)
- `integrator::Union{Nothing, Function}=nothing`: Optional integrator factory function. When
  provided, it is called as `integrator(sampling_qtraj, N)` and must return an integrator or
  vector of integrators. When `nothing` (default), `BilinearIntegrator` is used.
- `calibration_targets::Vector{Symbol}=Symbol[]`: Names of globals declared as **calibration targets** — knobs an external calibration step manages, not free NLP variables. SamplingProblem builds a fresh constraint list (rather than inheriting from the base `qcp`), so calibration_target pins set on the base `qcp` are *not* automatically carried over — pass them here explicitly. Default empty: globals stay free.
- `piccolo_options::PiccoloOptions=PiccoloOptions()`: Options for the solver

# Returns
- `SamplingProblem{typeof(qcp)}`: a wrapper problem carrying the sampling trajectory and
  the wrapped problem in its type
"""
function SamplingProblem(
    qcp::AbstractQuantumControlProblem,
    systems::Vector{<:AbstractQuantumSystem};
    weights::Vector{Float64} = fill(1.0, length(systems)),
    Q::Float64 = 100.0,
    integrator::Union{Nothing,Function} = nothing,
    calibration_targets::Vector{Symbol} = Symbol[],
    piccolo_options::PiccoloOptions = PiccoloOptions(),
)
    if _show_header(piccolo_options)
        println("constructing SamplingProblem")
        println("    systems: $(length(systems))")
    end

    base_qtraj = quantum_trajectory(qcp)
    state_sym = state_name(base_qtraj)
    base_traj = get_trajectory(qcp)

    # 1. Create SamplingTrajectory wrapper (new API: no stored trajectory)
    sampling_qtraj = SamplingTrajectory(base_qtraj, systems; weights)

    # 2. Build trajectory from sampling trajectory (this creates duplicated states)
    #    Propagate Δt bounds and global variables from base problem
    N = base_traj.N
    Δt_bounds = if haskey(base_traj.bounds, :Δt)
        (base_traj.bounds[:Δt][1][1], base_traj.bounds[:Δt][2][1])
    else
        nothing
    end

    new_traj = NamedTrajectory(sampling_qtraj, N; Δt_bounds = Δt_bounds)

    # Propagate global variables (e.g., free-phase φ_1, φ_2) from base trajectory.
    # Directly mutate the struct fields — NamedTrajectory doesn't have a public API
    # for adding globals after construction.
    if base_traj.global_dim > 0
        new_traj = NamedTrajectory(
            new_traj.datavec,
            new_traj.components,
            new_traj.N;
            timestep = new_traj.timestep,
            controls = new_traj.control_names,
            bounds = new_traj.bounds,
            initial = new_traj.initial,
            final = isnothing(new_traj.final_) ? NamedTuple() : new_traj.final_,
            goal = new_traj.goal,
            global_data = copy(base_traj.global_data),
            global_components = base_traj.global_components,
        )
    end
    snames = state_names(sampling_qtraj)

    # 2b. Preserve the base problem's control-derivative structure. A smooth base
    #     problem carries :du, :ddu components, DerivativeIntegrators, and
    #     regularizers on them; the sampling rebuild keeps that chain (shared
    #     across members, like the control itself) so the robust problem is the
    #     same problem as the base. The chain is detected by the `d`-prefix
    #     naming convention of `add_control_derivatives` (:du, :ddu, …).
    #     Components the sampling conversion already carries from the pulse data
    #     (e.g. a CubicSplinePulse's :du Hermite tangents) are NOT re-added.
    control_sym = drive_name(base_qtraj)
    deriv_names = Symbol[]
    while Symbol("d"^(length(deriv_names) + 1) * string(control_sym)) ∈ base_traj.names
        push!(deriv_names, Symbol("d"^(length(deriv_names) + 1) * string(control_sym)))
    end

    n_existing = 0
    while n_existing < length(deriv_names) && deriv_names[n_existing+1] ∈ new_traj.names
        n_existing += 1
    end

    if 0 < n_existing < length(deriv_names)
        error(
            "SamplingProblem cannot preserve the base problem's derivative chain " *
            "$(deriv_names): the sampling trajectory already carries " *
            "$(deriv_names[1:n_existing]) from the pulse data, and partially " *
            "extending a derivative chain is not supported. Rebuild the base " *
            "problem with a matching chain, or open an issue.",
        )
    elseif n_existing == 0 && !isempty(deriv_names)
        derivative_bounds = if all(haskey(base_traj.bounds, dn) for dn in deriv_names)
            Tuple(base_traj.bounds[dn] for dn in deriv_names)
        else
            nothing
        end
        new_traj = add_control_derivatives(
            new_traj,
            length(deriv_names);
            control_name = control_sym,
            derivative_bounds = derivative_bounds,
        )
    end

    # 3. Build objective: weighted per-member state objectives + shared regularization.
    #    member_states is one Symbol per member for single-state bases, one
    #    Vector{Symbol} per member for multi-state bases (multi-ket, multi-density).
    member_states = sampling_member_states(sampling_qtraj)
    J_state = sum(
        sampling_state_objective(base_qtraj, new_traj, names, w * Q) for
        (names, w) in zip(member_states, weights)
    )
    J_reg = extract_regularization(direct_problem(qcp).objective, state_sym, new_traj)
    J_total = J_state + J_reg

    # 4. Build integrators: dynamics for each system (instance / per-slot vector /
    #    factory — validated at construction), then the base problem's derivative
    #    chain (u → du → ddu …), preserved from step 2b.
    n_slots = sum(sampling_member_states(sampling_qtraj)) do ms
        ms isa Symbol ? 1 : length(ms)
    end
    all_integrators = _resolve_sampling_integrators(integrator, sampling_qtraj, N, n_slots)

    # Replicate the base problem's derivative integrators exactly (smooth:
    # u→du→ddu; bang-bang / linear spline: u→du; cubic spline: none — its :du
    # Hermite tangents are free DOFs, not a constrained chain).
    for base_int in direct_problem(qcp).integrators
        base_int isa DerivativeIntegrator || continue
        if base_int.x_name ∈ new_traj.names && base_int.ẋ_name ∈ new_traj.names
            push!(
                all_integrators,
                DerivativeIntegrator(base_int.x_name, base_int.ẋ_name, new_traj),
            )
        else
            error(
                "SamplingProblem cannot preserve the base problem's " *
                "DerivativeIntegrator ($(base_int.x_name) → $(base_int.ẋ_name)): " *
                "those components do not exist in the sampling trajectory. " *
                "Only control-derivative chains (:du, :ddu, …) are preserved.",
            )
        end
    end

    # 5. Construct problem (TimeConsistencyConstraint auto-applied)
    constraints = AbstractConstraint[]

    # Pin calibration targets at nominal — knobs the user has declared an
    # external calibration step will manage, not free NLP variables. SamplingProblem
    # builds a fresh constraint list (cf. base qcp's), so this is not inherited
    # automatically; the user passes the list explicitly.
    apply_calibration_targets!(
        constraints,
        calibration_targets,
        new_traj;
        verbose = _show_details(piccolo_options),
    )

    prob =
        DirectTrajOptProblem(new_traj, J_total, all_integrators; constraints = constraints)

    return _maybe_display(
        SamplingProblem(
            qcp,
            sampling_qtraj,
            prob,
            SamplingParams(;
                Q = Q,
                weights = weights,
                calibration_targets = calibration_targets,
            ),
            retained_spec(qcp),
        ),
        piccolo_options,
    )
end

# ============================================================================= #
# Composability with MinimumTimeProblem
# ============================================================================= #

function _update_goal(qtraj::SamplingTrajectory, new_goal)
    new_base = _update_goal(qtraj.base_trajectory, new_goal)
    # Reconstruct the SamplingTrajectory around the updated base — a thin
    # SamplingTrajectory only owns (base_trajectory, systems, weights).
    return SamplingTrajectory(new_base, qtraj.systems; weights = qtraj.weights)
end

function _final_fidelity_constraint(
    qtraj::SamplingTrajectory,
    final_fidelity::Float64,
    traj::NamedTrajectory;
    subsystem_levels::Union{Nothing,Vector{Int}} = nothing,
)
    constraints = [
        _sampling_fidelity_constraint(qtraj.base_trajectory, name, final_fidelity, traj) for name in state_names(qtraj)
    ]
    return constraints
end

# Dispatch on trajectory type for fidelity constraint
function _sampling_fidelity_constraint(
    qtraj::UnitaryTrajectory,
    state_sym::Symbol,
    final_fidelity::Float64,
    traj::NamedTrajectory,
)
    return FinalUnitaryFidelityConstraint(qtraj.goal, state_sym, final_fidelity, traj)
end

function _sampling_fidelity_constraint(
    qtraj::KetTrajectory,
    state_sym::Symbol,
    final_fidelity::Float64,
    traj::NamedTrajectory,
)
    return FinalKetFidelityConstraint(qtraj.goal, state_sym, final_fidelity, traj)
end

# Tests
@testitem "SamplingProblem Construction" begin
    using DirectTrajOpt

    T = 10.0
    N = 50

    # Define system
    sys = QuantumSystem(GATES[:Z], [GATES[:X]], [1.0])

    # Create base problem
    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    qtraj = UnitaryTrajectory(sys, pulse, GATES[:H])
    qcp = SmoothPulseProblem(qtraj, N; Q = 100.0)

    # Create sampling problem with 2 systems
    systems = [sys, sys] # Identical systems for testing
    sampling_prob = SamplingProblem(qcp, systems)

    @test sampling_prob isa SamplingProblem
    @test sampling_prob isa AbstractQuantumControlProblem
    @test inner(sampling_prob) isa QuantumControlProblem
    @test sampling_prob.qtraj isa SamplingTrajectory
    @test length(sampling_prob.qtraj.systems) == 2

    # Check trajectory components (now use numbered suffix :Ũ⃗1, :Ũ⃗2, etc.)
    traj = get_trajectory(sampling_prob)
    @test haskey(traj.components, :Ũ⃗1)
    @test haskey(traj.components, :Ũ⃗2)
    @test haskey(traj.components, :u)

    # Check integrators
    # 2 dynamics integrators (one per system) + 2 derivative integrators
    # carried over from the smooth base problem's :du/:ddu chain
    @test length(sampling_prob.prob.integrators) == 4
end

@testitem "SamplingProblem preserves smooth base structure" begin
    using DirectTrajOpt
    using LinearAlgebra

    # Regression pin for the silent structure drop: a SamplingProblem built from
    # a smooth base problem must retain the base's derivative components (:du,
    # :ddu), its DerivativeIntegrators, and its regularizer objectives — the
    # robust problem is otherwise not the same problem as the base.

    T = 10.0
    N = 50

    sys = QuantumSystem(GATES[:Z], [GATES[:X]], [1.0])

    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    qtraj = UnitaryTrajectory(sys, pulse, GATES[:H])
    qcp = SmoothPulseProblem(qtraj, N; Q = 100.0)

    sampling_prob = SamplingProblem(qcp, [sys, sys])
    traj = get_trajectory(sampling_prob)

    # Derivative components carried into the sampling trajectory
    @test haskey(traj.components, :du)
    @test haskey(traj.components, :ddu)

    # Derivative integrators carried: 2 dynamics (one per member) + 2 derivative
    n_derivative = count(i -> i isa DerivativeIntegrator, sampling_prob.prob.integrators)
    @test n_derivative == 2
    @test length(sampling_prob.prob.integrators) == 4

    # Regularizer objectives carried: the smooth base's quadratic regularizers
    # on :u, :du, :ddu survive the rebuild (recursively flatten composites).
    function _regularizer_names(obj)
        terms = hasproperty(obj, :objectives) ? obj.objectives : (obj,)
        names = Symbol[]
        for t in terms
            if hasproperty(t, :objectives)
                append!(names, _regularizer_names(t))
            elseif t isa QuadraticRegularizer
                push!(names, t.name)
            end
        end
        return names
    end
    reg_names = _regularizer_names(sampling_prob.prob.objective)
    @test :u ∈ reg_names
    @test :du ∈ reg_names
    @test :ddu ∈ reg_names

    # The carried structure is live: objective evaluates finite, and the
    # derivative constraints are consistent at the initial trajectory.
    @test isfinite(sampling_prob.prob.objective(sampling_prob.trajectory))
    for integrator in sampling_prob.prob.integrators
        if integrator isa DerivativeIntegrator
            δ = zeros(integrator.dim)
            DirectTrajOpt.evaluate!(δ, integrator, traj)
            @test norm(δ, Inf) < 1e-8
        end
    end

    # --- Spline bases: the pulse carries :du, so the rebuild must not
    #     duplicate it; which DerivativeIntegrators survive follows the base. ---

    times = collect(range(0.0, T, length = N))

    # Cubic spline: :du is Hermite tangents (free DOFs) — no DerivativeIntegrator
    cubic_pulse = CubicSplinePulse(0.1 * randn(1, N), 0.1 * randn(1, N), times)
    cubic_qtraj = UnitaryTrajectory(sys, cubic_pulse, GATES[:H])
    cubic_qcp = SplinePulseProblem(cubic_qtraj, N; Q = 100.0, integrator_type = :pwc)

    cubic_sampling = SamplingProblem(cubic_qcp, [sys, sys])
    cubic_traj = get_trajectory(cubic_sampling)
    @test haskey(cubic_traj.components, :du)
    @test count(==(:du), cubic_traj.names) == 1
    @test count(i -> i isa DerivativeIntegrator, cubic_sampling.prob.integrators) == 0
    @test :du ∈ _regularizer_names(cubic_sampling.prob.objective)
    @test isfinite(cubic_sampling.prob.objective(cubic_sampling.trajectory))

    # Linear spline: :du is a constrained derivative — one DerivativeIntegrator
    linear_pulse = LinearSplinePulse(0.1 * randn(1, N), times)
    linear_qtraj = UnitaryTrajectory(sys, linear_pulse, GATES[:H])
    linear_qcp = SplinePulseProblem(linear_qtraj, N; Q = 100.0)

    linear_sampling = SamplingProblem(linear_qcp, [sys, sys])
    linear_traj = get_trajectory(linear_sampling)
    @test haskey(linear_traj.components, :du)
    @test count(==(:du), linear_traj.names) == 1
    @test count(i -> i isa DerivativeIntegrator, linear_sampling.prob.integrators) == 1
    @test :du ∈ _regularizer_names(linear_sampling.prob.objective)
    @test isfinite(linear_sampling.prob.objective(linear_sampling.trajectory))
end

@testitem "SamplingProblem Solving" tags = [:sampling_problem] begin
    using DirectTrajOpt

    T = 10.0
    N = 50

    # Simple robust optimization
    # System with uncertainty in drift
    sys_nominal = QuantumSystem(GATES[:Z], [GATES[:X]], [1.0])
    sys_perturbed = QuantumSystem(1.1 * GATES[:Z], [GATES[:X]], [1.0])

    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    qtraj = UnitaryTrajectory(sys_nominal, pulse, GATES[:X])
    qcp = SmoothPulseProblem(qtraj, N; Q = 100.0)

    sampling_prob = SamplingProblem(qcp, [sys_nominal, sys_perturbed])

    # Solve
    solve!(sampling_prob; max_iter = 10, verbose = false)

    # Check that we have a solution
    @test sampling_prob.prob.objective(sampling_prob.trajectory) < 1e10 # Just check it didn't blow up
end

@testitem "SamplingProblem with KetTrajectory" begin
    using DirectTrajOpt

    T = 1.0
    N = 50

    # Robust state transfer over parameter uncertainty
    sys_nominal = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    sys_perturbed = QuantumSystem(1.1 * GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    pulse = ZeroOrderPulse(0.1 * randn(2, N), collect(range(0.0, T, length = N)))
    qtraj = KetTrajectory(sys_nominal, pulse, ψ_init, ψ_goal)

    qcp = SmoothPulseProblem(qtraj, N; Q = 50.0, R = 1e-3)

    # Create sampling problem
    sampling_prob = SamplingProblem(qcp, [sys_nominal, sys_perturbed]; Q = 50.0)

    @test sampling_prob isa SamplingProblem
    @test sampling_prob isa AbstractQuantumControlProblem
    @test inner(sampling_prob) isa QuantumControlProblem
    @test sampling_prob.qtraj isa SamplingTrajectory

    # Check trajectory has sample states (now use numbered suffix)
    traj = get_trajectory(sampling_prob)
    @test haskey(traj.components, :ψ̃1)
    @test haskey(traj.components, :ψ̃2)

    # Solve
    solve!(sampling_prob; max_iter = 10, verbose = false, print_level = 1)
end

@testitem "SamplingProblem with custom weights" begin
    using DirectTrajOpt

    T = 10.0
    N = 50

    sys1 = QuantumSystem(GATES[:Z], [GATES[:X]], [1.0])
    sys2 = QuantumSystem(1.1 * GATES[:Z], [GATES[:X]], [1.0])
    sys3 = QuantumSystem(0.9 * GATES[:Z], [GATES[:X]], [1.0])

    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    qtraj = UnitaryTrajectory(sys1, pulse, GATES[:X])
    qcp = SmoothPulseProblem(qtraj, N; Q = 100.0)

    # Non-uniform weights - emphasize nominal system
    weights = [0.6, 0.2, 0.2]
    sampling_prob = SamplingProblem(qcp, [sys1, sys2, sys3]; weights = weights, Q = 100.0)

    @test sampling_prob.qtraj.weights == weights
    @test length(sampling_prob.qtraj.systems) == 3

    # Should have 3 sample states (numbered suffix)
    traj = get_trajectory(sampling_prob)
    @test haskey(traj.components, :Ũ⃗1)
    @test haskey(traj.components, :Ũ⃗2)
    @test haskey(traj.components, :Ũ⃗3)

    solve!(sampling_prob; max_iter = 5, verbose = false, print_level = 1)
end

@testitem "SamplingProblem + MinimumTimeProblem composition" begin
    using DirectTrajOpt

    T = 1.0
    N = 50

    # Robust minimum-time optimization
    sys_nominal = QuantumSystem(0.1 * GATES[:Z], [GATES[:X]], [1.0])
    sys_perturbed = QuantumSystem(0.11 * GATES[:Z], [GATES[:X]], [1.0])

    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    qtraj = UnitaryTrajectory(sys_nominal, pulse, GATES[:X])
    qcp = SmoothPulseProblem(qtraj, N; Q = 100.0, R = 1e-2, Δt_bounds = (0.01, 0.5))

    # First create sampling problem
    sampling_prob = SamplingProblem(qcp, [sys_nominal, sys_perturbed]; Q = 100.0)
    solve!(sampling_prob; max_iter = 20, verbose = false, print_level = 1)

    # Then convert to minimum-time
    mintime_prob = MinimumTimeProblem(sampling_prob; final_fidelity = 0.90, D = 50.0)

    @test mintime_prob isa SamplingProblem
    @test mintime_prob isa AbstractQuantumControlProblem
    @test inner(mintime_prob) isa QuantumControlProblem
    @test mintime_prob.qtraj isa SamplingTrajectory

    # Solve minimum-time
    solve!(mintime_prob; max_iter = 20, verbose = false, print_level = 1)
end

@testitem "SamplingProblem with EmbeddedOperator" begin
    using DirectTrajOpt

    # Minimal setup (reproducing the bug from main.jl)
    T = 1.0
    N = 10
    times = collect(range(0, T, length = N))
    initial_controls = zeros(2, N)
    pulse = ZeroOrderPulse(initial_controls, times)

    # Create systems
    sys1 = TransmonSystem(levels = 2, δ = -0.18)
    sys2 = TransmonSystem(levels = 2, δ = -0.20)

    # Create embedded operator (this is what caused the bug)
    X_embedded = EmbeddedOperator(GATES[:X], sys2)

    # Create trajectory with embedded operator as goal
    qtraj = UnitaryTrajectory(sys2, pulse, X_embedded)
    qcp = SmoothPulseProblem(qtraj, N; Q = 100.0)

    # This should not fail - it's the bug we're fixing
    sampling_prob = SamplingProblem(qcp, [sys1, sys2])

    @test sampling_prob isa SamplingProblem
    @test sampling_prob isa AbstractQuantumControlProblem
    @test inner(sampling_prob) isa QuantumControlProblem
    @test sampling_prob.qtraj isa SamplingTrajectory
    @test length(sampling_prob.qtraj.systems) == 2
end

@testitem "SamplingProblem with DensityTrajectory" tags = [:density, :skip] begin
    # TODO: DensityTrajectory support for SamplingProblem is not yet complete
    # Needs: BilinearIntegrator dispatch, SamplingTrajectory NamedTrajectory conversion
    @test_skip "DensityTrajectory support not yet implemented"
end

@testitem "SamplingProblem with custom integrator factory" begin
    using DirectTrajOpt
    using LinearAlgebra

    T = 10.0
    N = 50

    sys_nominal = QuantumSystem(GATES[:Z], [GATES[:X]], [1.0])
    sys_perturbed = QuantumSystem(1.1 * GATES[:Z], [GATES[:X]], [1.0])

    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    qtraj = UnitaryTrajectory(sys_nominal, pulse, GATES[:X])
    qcp = SmoothPulseProblem(qtraj, N; Q = 100.0)

    # Custom integrator factory — reimplements default BilinearIntegrator logic
    custom_factory(sqtraj, n) = BilinearIntegrator(sqtraj, n)

    sampling_prob =
        SamplingProblem(qcp, [sys_nominal, sys_perturbed]; integrator = custom_factory)

    @test sampling_prob isa SamplingProblem
    @test sampling_prob isa AbstractQuantumControlProblem
    @test inner(sampling_prob) isa QuantumControlProblem
    # 2 factory-provided dynamics integrators + the 2 DerivativeIntegrators
    # (u→du, du→ddu) preserved from the smooth base problem.
    @test length(sampling_prob.prob.integrators) == 4

    solve!(sampling_prob; max_iter = 5, verbose = false, print_level = 1)
end

@testitem "SamplingProblem is a parametric wrapper: wrap history lives in the type" begin
    using DirectTrajOpt
    using LinearAlgebra

    T = 1.0
    N = 8
    sys_a = QuantumSystem(0.1 * GATES[:Z], [GATES[:X]], [1.0])
    sys_b = QuantumSystem(0.11 * GATES[:Z], [GATES[:X]], [1.0])
    opts = PiccoloOptions(display = :silent)

    times = collect(range(0.0, T, length = N))
    qtraj = UnitaryTrajectory(sys_a, ZeroOrderPulse(0.1 * randn(1, N), times), GATES[:X])
    base = SmoothPulseProblem(
        qtraj,
        N;
        Q = 100.0,
        R = 1e-2,
        Δt_bounds = (0.01, 0.5),
        piccolo_options = opts,
    )
    sp = SamplingProblem(base, [sys_a, sys_b]; Q = 100.0, piccolo_options = opts)

    # the type IS the ordered wrapper list
    @test sp isa SamplingProblem
    @test sp isa AbstractProblemWrapper
    @test sp isa AbstractQuantumControlProblem
    @test !(sp isa QuantumControlProblem)          # a wrapper is NOT the tagged type
    @test inner(sp) === base
    @test inner(sp) isa SmoothPulseProblem
    @test base_problem(sp) === base

    # wrapper identity + params
    @test wrapper_kind(sp) === :sampling
    @test template_params(sp) isa SamplingParams
    @test template_params(sp).Q == 100.0
    @test length(template_params(sp).weights) == 2
    @test retained_spec(sp) === nothing

    # the template tag of the wrapped problem shows through
    @test template_tag(sp) === template_tag(base)

    # accessors forward through the wrapper's own trajectory, not the inner one
    @test quantum_trajectory(sp) isa SamplingTrajectory
    @test get_system(sp) isa AbstractQuantumSystem
    @test get_trajectory(sp) === direct_problem(sp).trajectory
    @test drive_name(sp) === drive_name(quantum_trajectory(sp))
    @test sp.objective === direct_problem(sp).objective   # property forwarding

    # min-time through a wrapper returns the SAME wrapper type (no flattening)
    mt = MinimumTimeProblem(sp; final_fidelity = 0.5, D = 10.0, piccolo_options = opts)
    @test mt isa SamplingProblem
    @test inner(mt) === base
    @test wrapper_kind(mt) === :sampling
    @test template_params(mt) === template_params(sp)
    @test length(direct_problem(mt).constraints) > length(direct_problem(sp).constraints)   # gained the fidelity constraint

    # min-time on a tagged problem preserves the tag and the retained params —
    # it is a recipe over the composition axes, not a wrapper
    mt_base =
        MinimumTimeProblem(base; final_fidelity = 0.5, D = 10.0, piccolo_options = opts)
    @test mt_base isa SmoothPulseProblem{<:UnitaryTrajectory}
    @test !(mt_base isa AbstractProblemWrapper)
    @test template_tag(mt_base) === SmoothPulseTemplate()
    @test template_params(mt_base) === template_params(base)
end
