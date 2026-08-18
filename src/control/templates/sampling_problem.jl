export SamplingProblem

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
    SamplingProblem(qcp::QuantumControlProblem, systems::Vector{<:AbstractQuantumSystem}; kwargs...)

Construct a `SamplingProblem` from an existing `QuantumControlProblem` and a list of systems.

This creates a robust optimization problem where the controls are shared across all systems,
but each system evolves according to its own dynamics. The objective is the weighted sum of
fidelity objectives for each system.

# Arguments
- `qcp::QuantumControlProblem`: The base problem (defines nominal trajectory, objective, etc.)
- `systems::Vector{<:AbstractQuantumSystem}`: List of systems to optimize over

# Keyword Arguments
- `weights::Vector{Float64}=fill(1.0, length(systems))`: Weights for each system
- `Q::Float64=100.0`: Weight on infidelity objective (explicit, not extracted from base problem)
- `integrator::Union{Nothing, AbstractIntegrator, Vector{<:AbstractIntegrator}, Function}=nothing`:
  Optional integrator specification, aligned with the other problem templates. Accepts an
  integrator instance (single-slot ensembles only), a vector with one integrator per
  ensemble dynamics slot, or a factory called as `integrator(sampling_qtraj, N)` returning
  an integrator or vector of integrators. When `nothing` (default), `BilinearIntegrator`
  is used. Malformed shapes fail at construction with an actionable error.
- `calibration_targets::Vector{Symbol}=Symbol[]`: Names of globals declared as **calibration targets** — knobs an external calibration step manages, not free NLP variables. SamplingProblem builds a fresh constraint list (rather than inheriting from the base `qcp`), so calibration_target pins set on the base `qcp` are *not* automatically carried over — pass them here explicitly. Default empty: globals stay free.
- `piccolo_options::PiccoloOptions=PiccoloOptions()`: Options for the solver

# Structure preservation

The rebuild preserves the base problem's control-derivative structure: derivative
components (`:du`, `:ddu`, …) are added to the sampling trajectory (shared across
members, like the control itself), the base's `DerivativeIntegrator` chain is
carried over, and regularizer objectives from the base whose variables exist in
the sampling trajectory are retained. A robust version of a smooth problem is
therefore the same problem, ensemble-replicated.

**Known limitation:** bang-bang L1/slack structure does not survive the rebuild.
The `:s_du` slack component, its `LinearRegularizer`, and the `L1SlackConstraint`
are tied to the base problem's constraint list, which `SamplingProblem` rebuilds
fresh (cf. `calibration_targets`). A bang-bang base yields a one-derivative
smooth sampling problem — the `:du` chain and control regularizers survive, the
L1 sparsity machinery does not.

# Returns
- `QuantumControlProblem{SamplingTrajectory}`: A new problem with the sampling trajectory
"""
function SamplingProblem(
    qcp::QuantumControlProblem,
    systems::Vector{<:AbstractQuantumSystem};
    weights::Vector{Float64} = fill(1.0, length(systems)),
    Q::Float64 = 100.0,
    integrator::Union{Nothing,AbstractIntegrator,AbstractVector,Function} = nothing,
    calibration_targets::Vector{Symbol} = Symbol[],
    piccolo_options::PiccoloOptions = PiccoloOptions(),
)
    if _show_header(piccolo_options)
        println("constructing SamplingProblem")
        println("    systems: $(length(systems))")
    end

    base_qtraj = qcp.qtraj
    state_sym = state_name(base_qtraj)
    control_sym = drive_name(base_qtraj)
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
    # 2b. Preserve the base problem's control-derivative structure. A smooth base
    #     problem carries :du, :ddu components, DerivativeIntegrators, and
    #     regularizers on them; the sampling rebuild keeps that chain (shared
    #     across members, like the control itself) so the robust problem is the
    #     same problem as the base. The chain is detected by the `d`-prefix
    #     naming convention of `add_control_derivatives` (:du, :ddu, …).
    #     Components the sampling conversion already carries from the pulse data
    #     (e.g. a CubicSplinePulse's :du Hermite tangents) are NOT re-added.
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
    J_reg = extract_regularization(qcp.prob.objective, state_sym, new_traj)
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
    for base_int in qcp.prob.integrators
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

    return _maybe_display(QuantumControlProblem(sampling_qtraj, prob), piccolo_options)
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
    # Per-member constraints: one name per member for single-state bases, one
    # name-vector per member for multi-state bases. A method may return a single
    # constraint or a vector of them (multi-density) — flatten either way.
    constraints = AbstractConstraint[]
    for member_states in sampling_member_states(qtraj)
        c = _sampling_fidelity_constraint(
            qtraj.base_trajectory,
            member_states,
            final_fidelity,
            traj,
        )
        c isa AbstractVector ? append!(constraints, c) : push!(constraints, c)
    end
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

function _sampling_fidelity_constraint(
    qtraj::MultiKetTrajectory,
    state_syms::Vector{Symbol},
    final_fidelity::Float64,
    traj::NamedTrajectory,
)
    # Per-member coherent fidelity constraint over the member's ket sub-states
    return FinalCoherentKetFidelityConstraint(
        qtraj.goals,
        state_syms,
        final_fidelity,
        traj;
        weights = qtraj.weights,
    )
end

# Density bases: the min-time constraint cell IS supported publicly (unlike the
# objective cell) — FinalDensityFidelityConstraint is public machinery. It only
# becomes reachable in a full solve once a downstream density sampling objective
# registers (see the sampling_state_objective extension point).
function _sampling_fidelity_constraint(
    qtraj::DensityTrajectory,
    state_sym::Symbol,
    final_fidelity::Float64,
    traj::NamedTrajectory,
)
    return FinalDensityFidelityConstraint(qtraj.goal, state_sym, final_fidelity, traj)
end

function _sampling_fidelity_constraint(
    qtraj::MultiDensityTrajectory,
    state_syms::Vector{Symbol},
    final_fidelity::Float64,
    traj::NamedTrajectory,
)
    # One density fidelity constraint per sub-state of this member
    return [
        FinalDensityFidelityConstraint(goal, name, final_fidelity, traj) for
        (goal, name) in zip(qtraj.goals, state_syms)
    ]
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

    @test sampling_prob isa QuantumControlProblem
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

    @test sampling_prob isa QuantumControlProblem
    @test sampling_prob.qtraj isa SamplingTrajectory

    # Check trajectory has sample states (now use numbered suffix)
    traj = get_trajectory(sampling_prob)
    @test haskey(traj.components, :ψ̃1)
    @test haskey(traj.components, :ψ̃2)

    # Solve
    solve!(sampling_prob; max_iter = 10, verbose = false, print_level = 1)
end

@testitem "SamplingProblem with MultiKetTrajectory" begin
    using DirectTrajOpt

    T = 1.0
    N = 21

    # Systems with a global parameter (to pin globals propagation) and drift variation
    sys_nom = QuantumSystem(
        GATES[:Z],
        [GATES[:X], GATES[:Y]],
        [1.0, 1.0];
        global_params = (δ = 0.5,),
    )
    sys_var = QuantumSystem(
        1.1 * GATES[:Z],
        [GATES[:X], GATES[:Y]],
        [1.0, 1.0];
        global_params = (δ = 0.5,),
    )

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    pulse = ZeroOrderPulse(0.1 * randn(2, N), collect(range(0.0, T, length = N)))
    qtraj = MultiKetTrajectory(sys_nom, pulse, [ψ0, ψ1], [ψ1, ψ0])

    qcp = SmoothPulseProblem(qtraj, N; Q = 50.0, R = 1e-3)

    sampling_prob = SamplingProblem(qcp, [sys_nom, sys_var]; Q = 50.0)

    @test sampling_prob isa QuantumControlProblem
    @test sampling_prob.qtraj isa SamplingTrajectory

    traj = get_trajectory(sampling_prob)

    # 2 systems × 2 kets = 4 per-member state components, with bounds propagated
    for sn in [:ψ̃1, :ψ̃2, :ψ̃3, :ψ̃4]
        @test haskey(traj.components, sn)
        @test haskey(traj.bounds, sn)
    end

    # One shared control and one shared Δt
    @test haskey(traj.components, :u)
    @test :u ∈ traj.control_names
    @test :Δt ∈ traj.control_names

    # Globals propagated from the base problem's trajectory
    @test traj.global_dim == get_trajectory(qcp).global_dim
    @test haskey(traj.global_components, :δ)
    @test traj.global_data == get_trajectory(qcp).global_data

    # 4 dynamics integrators (one per (member, ket)) + 2 derivative integrators
    # carried over from the smooth base problem's :du/:ddu chain
    @test length(sampling_prob.prob.integrators) == 6

    # Per-member coherent objectives: evaluating must be finite
    @test isfinite(sampling_prob.prob.objective(sampling_prob.trajectory))
end

@testitem "SamplingProblem Solving with MultiKetTrajectory" tags = [:sampling_problem] begin
    using DirectTrajOpt

    T = 1.0
    N = 21

    # Robust multi-state transfer over drift uncertainty: X gate via |0⟩→|1⟩, |1⟩→|0⟩
    sys_nominal = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    sys_perturbed = QuantumSystem(1.1 * GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    pulse = ZeroOrderPulse(0.1 * randn(2, N), collect(range(0.0, T, length = N)))
    qtraj = MultiKetTrajectory(sys_nominal, pulse, [ψ0, ψ1], [ψ1, ψ0])

    qcp = SmoothPulseProblem(qtraj, N; Q = 50.0, R = 1e-3)

    sampling_prob = SamplingProblem(qcp, [sys_nominal, sys_perturbed]; Q = 50.0)

    # Short end-to-end optimization on CPU
    solve!(sampling_prob; max_iter = 10, verbose = false, print_level = 1)

    # Didn't blow up
    @test sampling_prob.prob.objective(sampling_prob.trajectory) < 1e10
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

    @test mintime_prob isa QuantumControlProblem
    @test mintime_prob.qtraj isa SamplingTrajectory

    # Solve minimum-time
    solve!(mintime_prob; max_iter = 20, verbose = false, print_level = 1)
end

@testitem "SamplingProblem + MinimumTimeProblem composition (MultiKet)" begin
    using DirectTrajOpt

    T = 1.0
    N = 21

    sys_nominal = QuantumSystem(0.1 * GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    sys_perturbed = QuantumSystem(0.11 * GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    pulse = ZeroOrderPulse(0.1 * randn(2, N), collect(range(0.0, T, length = N)))
    qtraj = MultiKetTrajectory(sys_nominal, pulse, [ψ0, ψ1], [ψ1, ψ0])

    qcp = SmoothPulseProblem(qtraj, N; Q = 100.0, R = 1e-2, Δt_bounds = (0.01, 0.5))

    sampling_prob = SamplingProblem(qcp, [sys_nominal, sys_perturbed]; Q = 100.0)
    solve!(sampling_prob; max_iter = 10, verbose = false, print_level = 1)

    # Per-member final-fidelity constraints: one coherent constraint per member
    # (2 members), each over the member's 2 ket components
    cons = Piccolo.ProblemTemplates._final_fidelity_constraint(
        sampling_prob.qtraj,
        0.80,
        get_trajectory(sampling_prob),
    )
    @test length(cons) == 2
    @test all(c -> c isa NonlinearKnotPointConstraint, cons)

    mintime_prob = MinimumTimeProblem(sampling_prob; final_fidelity = 0.80, D = 50.0)

    @test mintime_prob isa QuantumControlProblem
    @test mintime_prob.qtraj isa SamplingTrajectory

    solve!(mintime_prob; max_iter = 10, verbose = false, print_level = 1)
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

    @test sampling_prob isa QuantumControlProblem
    @test sampling_prob.qtraj isa SamplingTrajectory
    @test length(sampling_prob.qtraj.systems) == 2
end

@testitem "SamplingProblem with DensityTrajectory fails loudly (extension point)" tags =
    [:density] begin
    using DirectTrajOpt

    T = 1.0
    N = 11

    # Open systems with dissipation and a drift variation
    L = ComplexF64[0.0 0.1; 0.0 0.0]
    sys_nom = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    sys_var =
        OpenQuantumSystem(0.95 * PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])

    ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρg = ComplexF64[0.0 0.0; 0.0 1.0]
    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    qtraj = DensityTrajectory(sys_nom, pulse, ρ0, ρg)
    qcp = SmoothPulseProblem(qtraj, N; Q = 100.0)

    # Piccolo ships no public density fidelity objective: construction must fail
    # loudly, naming the extension point — never silently install a null objective.
    err = try
        SamplingProblem(qcp, [sys_nom, sys_var])
        nothing
    catch e
        e
    end

    @test err isa Exception
    @test occursin("sampling_state_objective", sprint(showerror, err))
    @test occursin("Piccolissimo", sprint(showerror, err))
end

@testitem "SamplingTrajectory (Density) min-time fidelity constraint" tags = [:density] begin
    using DirectTrajOpt

    T = 1.0
    N = 11

    L = ComplexF64[0.0 0.1; 0.0 0.0]
    sys_nom = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    sys_var =
        OpenQuantumSystem(0.95 * PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])

    ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρg = ComplexF64[0.0 0.0; 0.0 1.0]
    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    base_qtraj = DensityTrajectory(sys_nom, pulse, ρ0, ρg)

    # Built directly: SamplingProblem construction errors loudly on the density
    # objective cell (pinned above), so the min-time machinery is exercised at
    # the dispatch level. Behavior pinned: per-member FinalDensityFidelityConstraint
    # (public machinery), usable as-is once a downstream objective registers.
    sampling_qtraj = SamplingTrajectory(base_qtraj, [sys_nom, sys_var])
    traj = NamedTrajectory(sampling_qtraj, N)

    cons = Piccolo.ProblemTemplates._final_fidelity_constraint(sampling_qtraj, 0.9, traj)
    @test length(cons) == 2
    @test all(c -> c isa NonlinearKnotPointConstraint, cons)
end

@testitem "SamplingProblem with MultiDensityTrajectory fails loudly (extension point)" tags =
    [:density] begin
    using DirectTrajOpt

    T = 1.0
    N = 11

    L = ComplexF64[0.0 0.1; 0.0 0.0]
    sys_nom = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    sys_var =
        OpenQuantumSystem(0.95 * PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])

    ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρ1 = ComplexF64[0.0 0.0; 0.0 1.0]
    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    qtraj = MultiDensityTrajectory(sys_nom, pulse, [ρ0, ρ1], [ρ1, ρ0])

    # No SmoothPulseProblem method exists for MultiDensityTrajectory; a minimal
    # base problem suffices — SamplingProblem only reads the trajectory and the
    # objective before the (loud) density objective cell is reached.
    base_prob = DirectTrajOptProblem(
        NamedTrajectory(qtraj, N),
        NullObjective(),
        AbstractIntegrator[],
    )
    qcp = QuantumControlProblem(qtraj, base_prob)

    # Same loud cell as the single-density base, reached through the per-member
    # (Vector{Symbol}) objective dispatch
    err = try
        SamplingProblem(qcp, [sys_nom, sys_var])
        nothing
    catch e
        e
    end

    @test err isa Exception
    @test occursin("sampling_state_objective", sprint(showerror, err))
    @test occursin("Piccolissimo", sprint(showerror, err))
end

@testitem "SamplingTrajectory (MultiDensity) min-time fidelity constraints" tags =
    [:density] begin
    using DirectTrajOpt

    T = 1.0
    N = 11

    L = ComplexF64[0.0 0.1; 0.0 0.0]
    sys_nom = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    sys_var =
        OpenQuantumSystem(0.95 * PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])

    ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρ1 = ComplexF64[0.0 0.0; 0.0 1.0]
    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    base_qtraj = MultiDensityTrajectory(sys_nom, pulse, [ρ0, ρ1], [ρ1, ρ0])

    sampling_qtraj = SamplingTrajectory(base_qtraj, [sys_nom, sys_var])
    traj = NamedTrajectory(sampling_qtraj, N)

    # One FinalDensityFidelityConstraint per (member, density) — 2 × 2 = 4
    cons = Piccolo.ProblemTemplates._final_fidelity_constraint(sampling_qtraj, 0.9, traj)
    @test length(cons) == 4
    @test all(c -> c isa NonlinearKnotPointConstraint, cons)
end

@testitem "SamplingProblem integrator keyword call shapes" begin
    using DirectTrajOpt
    using LinearAlgebra

    T = 10.0
    N = 50

    sys_nominal = QuantumSystem(GATES[:Z], [GATES[:X]], [1.0])
    sys_perturbed = QuantumSystem(1.1 * GATES[:Z], [GATES[:X]], [1.0])
    systems = [sys_nominal, sys_perturbed]

    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    qtraj = UnitaryTrajectory(sys_nominal, pulse, GATES[:X])
    qcp = SmoothPulseProblem(qtraj, N; Q = 100.0)

    # The three call shapes for the bilinear case — all equivalent to the
    # default. Integrators are name-based, so user-built integrators against
    # NamedTrajectory(sampling_qtraj, N) are valid inside the rebuilt problem.
    sampling_qtraj = SamplingTrajectory(qtraj, systems)
    default_vector = BilinearIntegrator(sampling_qtraj, N)
    @test default_vector isa AbstractVector && length(default_vector) == 2

    prob_default = SamplingProblem(qcp, systems)
    prob_factory =
        SamplingProblem(qcp, systems; integrator = (sq, n) -> BilinearIntegrator(sq, n))
    prob_vector = SamplingProblem(qcp, systems; integrator = default_vector)

    for p in (prob_default, prob_factory, prob_vector)
        # 2 bilinear dynamics + 2 derivative integrators from the smooth base
        @test count(i -> i isa BilinearIntegrator, p.prob.integrators) == 2
        @test length(p.prob.integrators) == 4
    end

    # Equivalent behavior: identical objective values at the initial trajectory
    J0 = prob_default.prob.objective(prob_default.trajectory)
    @test prob_factory.prob.objective(prob_factory.trajectory) ≈ J0
    @test prob_vector.prob.objective(prob_vector.trajectory) ≈ J0

    # Instance shape: a single integrator instance covers a single-member ensemble
    one_system = [sys_nominal]
    instance = BilinearIntegrator(SamplingTrajectory(qtraj, one_system), N)[1]
    prob_instance = SamplingProblem(qcp, one_system; integrator = instance)
    prob_default_one = SamplingProblem(qcp, one_system)
    @test count(i -> i isa BilinearIntegrator, prob_instance.prob.integrators) == 1
    @test prob_instance.prob.objective(prob_instance.trajectory) ≈
          prob_default_one.prob.objective(prob_default_one.trajectory)

    # Malformed call shapes fail at construction with actionable errors

    # Per-member vector with the wrong count (1 integrator for 2 members)
    err_count = try
        SamplingProblem(qcp, systems; integrator = [default_vector[1]])
        nothing
    catch e
        e
    end
    @test err_count isa Exception
    @test occursin("2", sprint(showerror, err_count))
    @test occursin("integrator", lowercase(sprint(showerror, err_count)))

    # Control: a well-formed per-member vector typed as AbstractIntegrator
    # constructs fine
    prob_typed = SamplingProblem(
        qcp,
        systems;
        integrator = AbstractIntegrator[default_vector[1], default_vector[1]],
    )
    @test prob_typed isa QuantumControlProblem

    # Vector with non-integrator contents
    err_junk = try
        SamplingProblem(qcp, systems; integrator = Any[default_vector[1], :not_an_integrator])
        nothing
    catch e
        e
    end
    @test err_junk isa Exception
    @test occursin("integrator", lowercase(sprint(showerror, err_junk)))

    # Factory returning something that is not integrator(s)
    err_factory = try
        SamplingProblem(qcp, systems; integrator = (sq, n) -> 42)
        nothing
    catch e
        e
    end
    @test err_factory isa Exception
    @test occursin("factory", lowercase(sprint(showerror, err_factory)))

    # Single instance for a multi-member ensemble leaves slots uncovered
    err_slots = try
        SamplingProblem(qcp, systems; integrator = default_vector[1])
        nothing
    catch e
        e
    end
    @test err_slots isa Exception
    @test occursin("integrator", lowercase(sprint(showerror, err_slots)))
end

@testitem "SamplingProblem from bang-bang base: L1 slack limitation" begin
    using DirectTrajOpt
    using LinearAlgebra

    # Pins the documented limitation (issue #267): bang-bang L1/slack structure
    # does NOT survive the sampling rebuild. The slack component :s_du, its
    # LinearRegularizer, and the L1SlackConstraint are bang-bang machinery tied
    # to the base problem's constraint list, which SamplingProblem deliberately
    # rebuilds fresh. The shared derivative chain (:du + DerivativeIntegrator)
    # and control regularizers DO survive — the result is a one-derivative
    # smooth sampling problem, not an L1-sparsified one.

    T = 10.0
    N = 50

    sys = QuantumSystem(GATES[:Z], [GATES[:X]], [1.0])
    sys_perturbed = QuantumSystem(1.1 * GATES[:Z], [GATES[:X]], [1.0])

    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    qtraj = UnitaryTrajectory(sys, pulse, GATES[:H])
    qcp = BangBangPulseProblem(qtraj, N; Q = 100.0)

    # Sanity: the bang-bang base really carries the slack structure
    @test haskey(get_trajectory(qcp).components, :s_du)

    sampling_prob = SamplingProblem(qcp, [sys, sys_perturbed])
    traj = get_trajectory(sampling_prob)

    # What survives: the 1-derivative chain and its integrator
    @test haskey(traj.components, :du)
    @test count(i -> i isa DerivativeIntegrator, sampling_prob.prob.integrators) == 1

    # What does not: no slack component, no L1 slack regularizer
    @test !haskey(traj.components, :s_du)
    function _has_linear_regularizer(obj)
        terms = hasproperty(obj, :objectives) ? obj.objectives : (obj,)
        return any(
            t ->
                hasproperty(t, :objectives) ? _has_linear_regularizer(t) :
                t isa LinearRegularizer,
            terms,
        )
    end
    @test !_has_linear_regularizer(sampling_prob.prob.objective)

    # The rebuilt problem is still well-posed
    @test isfinite(sampling_prob.prob.objective(sampling_prob.trajectory))
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

    @test sampling_prob isa QuantumControlProblem
    # 2 factory-supplied dynamics integrators + 2 derivative integrators carried
    # over from the smooth base problem's :du/:ddu chain
    @test length(sampling_prob.prob.integrators) == 4

    solve!(sampling_prob; max_iter = 5, verbose = false, print_level = 1)
end
