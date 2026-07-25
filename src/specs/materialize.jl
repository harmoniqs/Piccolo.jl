export materialize

using Random: MersenneTwister
using DirectTrajOpt:
    DirectTrajOptProblem,
    CompositeObjective,
    KnotPointObjective,
    NullObjective,
    QuadraticRegularizer,
    MinimumTimeObjective,
    AbstractObjective,
    AbstractConstraint

# ===========================================================================
# materialize — the function barrier between untyped spec data and the typed
# Piccolo world. For `control` specs it looks up registry factories, builds
# system → goal → pulse → trajectory, calls the existing template function, then
# applies the composition axes (Task 8) and wrappers (Task 9). All validation is
# a pre-pass returning structured `SpecError`s — never a materialization crash.
# ===========================================================================

"""
    materialize(spec::ProblemSpec) -> QuantumControlProblem

Materialize a `control` [`ProblemSpec`](@ref) into a live `QuantumControlProblem`
by looking up the registered system/template factories and calling the existing
Piccolo template function with mapped keyword arguments.

Phase-1 scope (wire-format-first): `system.kind ∈ {:template, :raw}` (`:composite`
deferred); solver backend execution is `ipopt`-only (backend dispatch is the
runner's concern, Task 12); the `:robust` wrapper is schema-only and returns a
structured "deferred" [`SpecError`](@ref). All trait/compatibility violations are
collected and thrown as a [`SpecValidationError`](@ref) *before* any Piccolo
object is constructed.
"""
function materialize(spec::ProblemSpec; piccolo_options = nothing)
    errs = SpecError[]
    _validate!(spec, errs)
    isempty(errs) || throw(SpecValidationError(errs))

    sys = _build_system(spec.system)
    goal = _build_goal(spec.goal, sys)
    qtraj = _build_pulse_trajectory(spec, sys, goal)
    base_qcp = _call_template(spec, qtraj; piccolo_options = piccolo_options)
    qcp = _apply_composition(spec, base_qcp)
    qcp = _apply_wrappers(spec, qcp; piccolo_options = piccolo_options)
    return qcp
end

# ---------------------------------------------------------------------------
# System
# ---------------------------------------------------------------------------

# Narrow a Vector{Any} of homogeneous numbers to a concrete Vector so template
# kwargs like `drive_bounds::Vector{Float64}` typecheck. TOML/JSON parsing boxes
# arrays as Vector{Any}; leave anything non-numeric untouched.
function _concretize(x)
    if x isa AbstractVector && !isempty(x)
        if all(v -> v isa Integer && !(v isa Bool), x)
            return Int[v for v in x]
        elseif all(v -> v isa Real && !(v isa Bool), x)
            return Float64[v for v in x]
        end
    end
    return x
end

_concretize_params(params::AbstractDict) =
    Dict{Symbol,Any}(k => _concretize(v) for (k, v) in params)

function _build_system(s::SystemSpec)
    if s.kind === :template
        entry = lookup_system(s.template)
        return entry.factory(; _concretize_params(s.params)...)
    elseif s.kind === :raw
        # Minimal Phase-1 raw path: build a QuantumSystem from explicit matrices.
        H_drift = _to_matrix(s.H_drift)
        H_drives = [_to_matrix(h) for h in s.H_drives]
        drive_bounds =
            _concretize(get(s.params, :drive_bounds, fill(1.0, length(H_drives))))
        return Quantum.QuantumSystem(H_drift, H_drives, drive_bounds)
    else
        # :composite deferred in Phase 1 (validation should have flagged it).
        throw(
            SpecValidationError([
                SpecError(
                    "system.kind",
                    "composite systems are deferred in Phase 1 (template + raw only)",
                    string(s.kind),
                ),
            ]),
        )
    end
end

# nested [re, im] pairs → ComplexF64 matrix (used by :raw H_drift/H_drives and
# inline unitary goal matrices).
function _to_matrix(x)
    x isa AbstractMatrix && return ComplexF64.(x)
    rows = collect(x)
    n = length(rows)
    M = Matrix{ComplexF64}(undef, n, length(rows[1]))
    for (i, row) in enumerate(rows), (j, entry) in enumerate(row)
        M[i, j] =
            entry isa AbstractVector ? ComplexF64(entry[1], entry[2]) : ComplexF64(entry)
    end
    return M
end

# ---------------------------------------------------------------------------
# Goal
# ---------------------------------------------------------------------------

function _build_goal(goal::GoalSpec, sys)
    if goal.kind === :unitary
        op = goal.gate !== nothing ? GATES[goal.gate] : _to_matrix(goal.matrix)
        if sys isa Quantum.CompositeQuantumSystem
            return EmbeddedOperator(op, sys)
        elseif goal.subspace !== nothing && length(goal.subspace) == 1
            return EmbeddedOperator(
                op,
                sys;
                subspace = goal.subspace[1],
                levels = sys.levels,
            )
        else
            return EmbeddedOperator(op, sys)
        end
    elseif goal.kind === :ket
        # Phase-1 ket goals: `target`/`initial` are Julia-parseable complex vectors.
        return _parse_ket(goal.target)
    else
        throw(
            SpecValidationError([
                SpecError(
                    "goal.kind",
                    "unknown goal kind",
                    string(goal.kind),
                    ["unitary", "ket"],
                ),
            ]),
        )
    end
end

_parse_ket(::Nothing) =
    throw(SpecValidationError([SpecError("goal.target", "ket goal requires a `target`")]))
_parse_ket(s::AbstractString) = ComplexF64.(eval(Meta.parse(s)))

# ---------------------------------------------------------------------------
# Pulse + trajectory (review correction 2: the pulse TYPE is fixed here, not by
# the template — build the concrete pulse ourselves, then wrap it).
# ---------------------------------------------------------------------------

function _controls(p::PulseSpec, n_drives::Int)
    N_seed = 11
    if p.init === :random
        rng = MersenneTwister(p.seed)
        return 0.1 .* randn(rng, n_drives, N_seed)
    else
        return zeros(n_drives, N_seed)
    end
end

function _build_pulse_trajectory(spec::ProblemSpec, sys, goal)
    p = spec.pulse
    N = spec.problem.N
    n_drives = sys.n_drives
    times = collect(range(0.0, p.T, N))
    controls =
        p.init === :random ? 0.1 .* randn(MersenneTwister(p.seed), n_drives, N) :
        zeros(n_drives, N)
    pulse = if p.kind === :zero_order
        ZeroOrderPulse(controls, times)
    elseif p.kind === :linear_spline
        LinearSplinePulse(controls, times)
    elseif p.kind === :cubic_spline
        CubicSplinePulse(controls, zeros(n_drives, N), times)
    else
        throw(
            SpecValidationError([
                SpecError(
                    "pulse.kind",
                    "unknown pulse kind",
                    string(p.kind),
                    ["zero_order", "linear_spline", "cubic_spline"],
                ),
            ]),
        )
    end

    gk = spec.goal.kind
    if gk === :unitary
        return UnitaryTrajectory(sys, pulse, goal)
    elseif gk === :ket
        ψ0 =
            spec.goal.initial === nothing ?
            ComplexF64[i == 1 ? 1.0 : 0.0 for i = 1:sys.levels] :
            _parse_ket(spec.goal.initial)
        return KetTrajectory(sys, pulse, ψ0, goal)
    else
        throw(
            SpecValidationError([
                SpecError("goal.kind", "unknown goal kind", string(gk), ["unitary", "ket"]),
            ]),
        )
    end
end

# ---------------------------------------------------------------------------
# Integrator block (Task 9): bilinear/omitted → nothing (template default path);
# exponential/spline → Piccolissimo factory (not loaded in this worktree, so
# validation blocks those before we get here).
# ---------------------------------------------------------------------------

function _build_integrator(spec::ProblemSpec, qtraj)
    spec.integrator === nothing && return nothing
    ik = spec.integrator.kind
    ik === :bilinear && return nothing
    entry = lookup_integrator(ik)
    return entry.factory(qtraj, spec.problem.N; alg = spec.integrator.alg)
end

# ---------------------------------------------------------------------------
# Template call (Task 7 happy path). Coercions per review correction 11.
# ---------------------------------------------------------------------------

# global_bounds::Dict{Symbol,Any} → the template's typed bounds dict.
function _coerce_global_bounds(gb::AbstractDict)
    out = Dict{Symbol,Union{Float64,Tuple{Float64,Float64}}}()
    for (k, v) in gb
        if v isa AbstractVector && length(v) == 2
            out[k] = (Float64(v[1]), Float64(v[2]))
        elseif v isa Tuple
            out[k] = (Float64(v[1]), Float64(v[2]))
        else
            out[k] = Float64(v)
        end
    end
    return out
end

function _call_template(spec::ProblemSpec, qtraj; piccolo_options = nothing)
    p = spec.problem
    fac = lookup_template(p.template).factory
    kwargs = Dict{Symbol,Any}(:Q => p.Q, :R => p.R, :free_phase => p.free_phase)
    p.R_u === nothing || (kwargs[:R_u] = p.R_u)
    p.R_du === nothing || (kwargs[:R_du] = p.R_du)
    p.R_ddu === nothing || (kwargs[:R_ddu] = p.R_ddu)
    isfinite(p.du_bound) && (kwargs[:du_bound] = p.du_bound)
    p.ddu_bound === nothing || (kwargs[:ddu_bound] = Float64(p.ddu_bound))
    p.initial_phases === nothing || (kwargs[:initial_phases] = p.initial_phases)
    isempty(p.calibration_targets) || (kwargs[:calibration_targets] = p.calibration_targets)
    isempty(p.global_names) || (kwargs[:global_names] = p.global_names)
    isempty(p.global_bounds) ||
        (kwargs[:global_bounds] = _coerce_global_bounds(p.global_bounds))
    p.free_dt isa Free && (kwargs[:Δt_bounds] = (p.free_dt.lo, p.free_dt.hi))
    intg = _build_integrator(spec, qtraj)
    intg === nothing || (kwargs[:integrator] = intg)
    piccolo_options === nothing || (kwargs[:piccolo_options] = piccolo_options)
    return fac(qtraj, p.N; kwargs...)
end

# ---------------------------------------------------------------------------
# Composition algebra (Task 8): goal_treatment, extra [[objectives]], and the
# min_time recipe. `free_dt → Δt_bounds` is already threaded into the template
# call (`_call_template`); here we shape the objective + fidelity constraint.
# ---------------------------------------------------------------------------

# A single-term objective is a bare objective (no `.objectives`); normalize to a
# vector so callers never assume `.objectives` exists (review nit 10).
_obj_terms(J) = J isa CompositeObjective ? collect(J.objectives) : AbstractObjective[J]

# The infidelity objective is the terminal-knot state objective the template adds
# first (`TerminalObjective` → a `KnotPointObjective` at `times == [N]` on the
# state var). Leakage (times `1:N`) and regularizers do not match, so this drops
# only the infidelity term under `goal_treatment=:constraint`.
_is_infidelity(o, state_sym::Symbol, N::Int) =
    o isa KnotPointObjective && state_sym in o.var_names && o.times == [N]

function _build_objective_term(o::ObjectiveTermSpec, qtraj, traj)
    k, w = o.kind, o.weight
    dsym = drive_name(qtraj)
    if k === :time
        return MinimumTimeObjective(traj; D = w)
    elseif k === :reg_u
        return QuadraticRegularizer(dsym, traj, w)
    elseif k === :reg_du
        return QuadraticRegularizer(Symbol(:d, dsym), traj, w)
    elseif k === :reg_ddu
        return QuadraticRegularizer(Symbol(:d, :d, dsym), traj, w)
    elseif k === :leakage
        idxs = get_iso_vec_leakage_indices(qtraj.goal)
        return LeakageObjective(idxs, state_name(qtraj), traj; Qs = fill(w, traj.N))
    elseif k === :sensitivity
        return UnitarySensitivityObjective(state_name(qtraj), traj, [traj.N]; Qs = [w])
    else
        throw(
            SpecValidationError([
                SpecError(
                    "problem.objectives",
                    "unsupported objective term kind",
                    string(k),
                ),
            ]),
        )
    end
end

function _apply_composition(spec::ProblemSpec, qcp)
    p = spec.problem
    treatment = p.goal_treatment
    extra = p.objectives
    (treatment === :objective && isempty(extra)) && return qcp

    prob = qcp.prob
    traj = prob.trajectory
    qtraj = qcp.qtraj
    state_sym = state_name(qtraj)

    terms = _obj_terms(prob.objective)
    if treatment === :constraint
        terms = filter(t -> !_is_infidelity(t, state_sym, traj.N), terms)
    end
    for o in extra
        push!(terms, _build_objective_term(o, qtraj, traj))
    end
    J = isempty(terms) ? NullObjective(traj) : reduce(+, terms)

    constraints = copy(prob.constraints)
    if treatment in (:constraint, :both)
        ff = p.final_fidelity === nothing ? 0.99 : p.final_fidelity
        fc = Control.ProblemTemplates._final_fidelity_constraint(qtraj, ff, traj)
        fc isa AbstractVector ? append!(constraints, fc) : push!(constraints, fc)
    end

    new_prob = DirectTrajOptProblem(traj, J, prob.integrators; constraints = constraints)
    return QuantumControlProblem(qtraj, new_prob)
end

# ---------------------------------------------------------------------------
# Trait validations (Task 9) — every check is a registry `.compat`/param query
# that pushes a structured `SpecError`; `materialize` throws them as one
# `SpecValidationError` before building anything (never a materialization crash).
# ---------------------------------------------------------------------------

function _validate!(spec::ProblemSpec, errs::Vector{SpecError})
    p = spec.problem
    if p === nothing
        push!(errs, SpecError("problem", "missing [problem] block for a control spec"))
        return errs
    end

    te = lookup_template(p.template)
    if te === nothing
        push!(
            errs,
            SpecError(
                "problem.template",
                "unknown problem template",
                string(p.template),
                sort!(String[string(k) for k in keys(TEMPLATES)]),
            ),
        )
        return errs   # compat checks below need the entry
    end

    # system availability / composite deferral
    if spec.system.kind === :template
        if lookup_system(spec.system.template) === nothing
            push!(
                errs,
                SpecError(
                    "system.template",
                    "unknown system template",
                    string(spec.system.template),
                    sort!(String[string(k) for k in keys(SYSTEMS)]),
                ),
            )
        end
    elseif spec.system.kind === :composite
        push!(
            errs,
            SpecError(
                "system.kind",
                "composite systems are deferred in Phase 1 (template + raw only)",
                string(spec.system.kind),
            ),
        )
    end

    # integrator block
    ik = spec.integrator === nothing ? :bilinear : spec.integrator.kind
    if ik !== :bilinear && lookup_integrator(ik) === nothing
        push!(
            errs,
            SpecError(
                "integrator.kind",
                "integrator :$(ik) is not registered (Piccolissimo not loaded)",
                string(ik),
                ["bilinear"],
            ),
        )
    end

    # free_phase / globals require a non-bilinear integrator
    if (p.free_phase || !isempty(p.global_names)) && ik === :bilinear
        push!(
            errs,
            SpecError(
                "integrator",
                "free_phase/globals require an exponential or spline integrator, not bilinear",
            ),
        )
    end

    # per-template kwarg validity: R_ddu / ddu_bound only for SmoothPulseProblem
    if p.template !== :SmoothPulseProblem
        p.R_ddu === nothing || push!(
            errs,
            SpecError(
                "problem.R_ddu",
                "R_ddu is only valid for SmoothPulseProblem",
                p.R_ddu,
            ),
        )
        p.ddu_bound === nothing || push!(
            errs,
            SpecError(
                "problem.ddu_bound",
                "ddu_bound is only valid for SmoothPulseProblem",
                p.ddu_bound,
            ),
        )
    end

    # pulse/template + trajectory/template compatibility (from registry .compat)
    pulse_kinds = get(te.compat, :pulse_kinds, Symbol[])
    if spec.pulse !== nothing && !isempty(pulse_kinds) && !(spec.pulse.kind in pulse_kinds)
        push!(
            errs,
            SpecError(
                "pulse.kind",
                "pulse kind incompatible with $(p.template)",
                string(spec.pulse.kind),
                String[string(k) for k in pulse_kinds],
            ),
        )
    end
    traj_kinds = get(te.compat, :trajectory_kinds, Symbol[])
    if spec.goal !== nothing && !isempty(traj_kinds) && !(spec.goal.kind in traj_kinds)
        push!(
            errs,
            SpecError(
                "goal.kind",
                "goal kind incompatible with $(p.template)",
                string(spec.goal.kind),
                String[string(k) for k in traj_kinds],
            ),
        )
    end

    # ket + free_phase only where compat[:ket_free_phase] == true
    if spec.goal !== nothing &&
       spec.goal.kind === :ket &&
       p.free_phase &&
       !get(te.compat, :ket_free_phase, false)
        push!(
            errs,
            SpecError(
                "problem.free_phase",
                "$(p.template) does not support free_phase for ket goals",
            ),
        )
    end

    # objective terms: registered + template-specific requirements
    has_time = false
    for o in p.objectives
        has_time |= (o.kind === :time)
        if lookup_objective_term(o.kind) === nothing
            push!(
                errs,
                SpecError(
                    "problem.objectives",
                    "unknown/unavailable objective term :$(o.kind) (Piccolissimo may be required)",
                    string(o.kind),
                ),
            )
        end
        if o.kind === :reg_ddu && p.template !== :SmoothPulseProblem
            push!(
                errs,
                SpecError(
                    "problem.objectives",
                    "reg_ddu requires a ddu-carrying template (SmoothPulseProblem)",
                    string(o.kind),
                ),
            )
        end
    end

    # a `time` objective requires free Δt (fixed-dt time term is inert). The
    # reverse (free Δt without a time term) is valid Piccolo free-time behavior,
    # so only this direction is enforced (documented divergence from the plan's
    # strict "⟺").
    if has_time && !(p.free_dt isa Free)
        push!(
            errs,
            SpecError(
                "problem.free_dt",
                "a `time` objective requires free_dt = [lo, hi] (Δt must be free)",
            ),
        )
    end

    # calibration_targets ⊆ declared globals
    for ct in p.calibration_targets
        ct in p.global_names || push!(
            errs,
            SpecError(
                "problem.calibration_targets",
                "calibration target :$(ct) is not in global_names",
                string(ct),
            ),
        )
    end

    # wrappers: sampling OK; robust is schema-only (deferred structured error)
    for (i, w) in enumerate(spec.wrappers)
        if w.kind === :robust
            push!(
                errs,
                SpecError(
                    "wrappers[$i]",
                    "the robust wrapper is deferred in Phase 1 (schema-only)",
                ),
            )
        elseif w.kind !== :sampling
            push!(
                errs,
                SpecError(
                    "wrappers[$i]",
                    "unknown wrapper kind",
                    string(w.kind),
                    ["sampling"],
                ),
            )
        end
    end

    return errs
end

# ---------------------------------------------------------------------------
# Wrappers (Task 9). Phase-1 reachable wrapper is `sampling`; `robust` is blocked
# in validation. The sampling G-arm re-derives per variant by building a fresh
# system for each variant's `[system].params` overrides.
# ---------------------------------------------------------------------------

function _apply_wrappers(spec::ProblemSpec, qcp; piccolo_options = nothing)
    for w in spec.wrappers
        w.kind === :sampling &&
            (qcp = _apply_sampling(spec, qcp, w; piccolo_options = piccolo_options))
    end
    return qcp
end

function _apply_sampling(spec::ProblemSpec, qcp, w::WrapperSpec; piccolo_options = nothing)
    entry = lookup_system(spec.system.template)
    systems = [
        entry.factory(; _concretize_params(merge(spec.system.params, v))...) for
        v in w.variants
    ]
    weights = w.weights === nothing ? fill(1.0, length(systems)) : w.weights
    kwargs = Dict{Symbol,Any}(:weights => weights, :Q => spec.problem.Q)
    piccolo_options === nothing || (kwargs[:piccolo_options] = piccolo_options)
    isempty(spec.problem.calibration_targets) ||
        (kwargs[:calibration_targets] = spec.problem.calibration_targets)
    # SamplingProblem takes an integrator *factory* (not an instance); bilinear →
    # nothing (default). Non-bilinear factories come from Piccolissimo (not here).
    fac = _sampling_integrator_factory(spec)
    fac === nothing || (kwargs[:integrator] = fac)
    return SamplingProblem(qcp, systems; kwargs...)
end

function _sampling_integrator_factory(spec::ProblemSpec)
    spec.integrator === nothing && return nothing
    ik = spec.integrator.kind
    ik === :bilinear && return nothing
    entry = lookup_integrator(ik)
    alg = spec.integrator.alg
    return (sqtraj, N) -> entry.factory(sqtraj, N; alg = alg)
end

# ---------------------------------------------------------------------------
# Accessors
# ---------------------------------------------------------------------------

export get_variants

"""
    get_variants(qcp::QuantumControlProblem) -> Vector

The system variants a materialized problem optimizes over: the sampled `systems`
for a `SamplingProblem`-backed problem, else the single nominal system.
"""
get_variants(qcp::QuantumControlProblem) =
    qcp.qtraj isa Quantum.SamplingTrajectory ? qcp.qtraj.systems : [get_system(qcp)]

@testitem "materialize: cubic-spline X gate matches hand-built problem" begin
    using Piccolo, Piccolo.Specs
    MAT_CONTROL_TOML = """
    schema_version = 1
    kind = "control"
    [system]
    kind = "template"
    template = "TransmonSystem"
    params = { levels = 3, drive_bounds = [0.02, 0.02] }
    [goal]
    kind = "unitary"
    gate = "X"
    [pulse]
    kind = "cubic_spline"
    T = 10.0
    [problem]
    template = "SplinePulseProblem"
    N = 11
    """
    spec = Specs.parse_spec(MAT_CONTROL_TOML; format = :toml)
    qcp = Specs.materialize(spec)
    @test qcp isa QuantumControlProblem

    # hand-built equivalent — the pulse TYPE is fixed at trajectory construction,
    # NOT by the template (see plan review correction 2).
    sys = TransmonSystem(; levels = 3, drive_bounds = [0.02, 0.02])
    goal = EmbeddedOperator(GATES[:X], sys)
    times = collect(range(0.0, 10.0, 11))
    pulse = CubicSplinePulse(zeros(2, 11), zeros(2, 11), times)
    traj = UnitaryTrajectory(sys, pulse, goal)
    ref = SplinePulseProblem(traj, 11)

    # A single-term objective is a bare objective (no `.objectives`); normalize.
    nterms(J) = J isa DirectTrajOpt.CompositeObjective ? length(J.objectives) : 1
    @test nterms(qcp.prob.objective) == nterms(ref.prob.objective)
    @test get_trajectory(qcp).dims == get_trajectory(ref).dims
end

@testitem "composition: EnergyPolish, min_time recipe, constraint-only" begin
    using Piccolo, Piccolo.Specs

    _obj_terms(J) = J isa DirectTrajOpt.CompositeObjective ? collect(J.objectives) : [J]

    # EnergyPolish: goal_treatment=constraint + reg_du, fixed T → NO
    # MinimumTimeObjective in the built problem (spec success criterion #4).
    ENERGY_POLISH_TOML = """
    schema_version = 1
    kind = "control"
    [system]
    kind = "template"
    template = "TransmonSystem"
    params = { levels = 3, drive_bounds = [0.02, 0.02] }
    [goal]
    kind = "unitary"
    gate = "X"
    [pulse]
    kind = "cubic_spline"
    T = 10.0
    [problem]
    template = "SplinePulseProblem"
    N = 11
    goal_treatment = "constraint"
    final_fidelity = 0.999
    [[problem.objectives]]
    kind = "reg_du"
    weight = 0.01
    """
    ep = Specs.materialize(Specs.parse_spec(ENERGY_POLISH_TOML; format = :toml))
    @test !any(o -> o isa Piccolo.MinimumTimeObjective, _obj_terms(ep.prob.objective))

    # min_time recipe desugars to the same canonical spec as the hand-factored
    # both+time+free_dt form (spec success criterion #5). Phase-1 has no dedicated
    # `min_time` sugar field, so the "recipe" IS the explicit
    # goal_treatment=both + [[objectives]] time + free_dt form; the two fixtures
    # differ only by default-omission and key order, and must canonicalize equal.
    MINTIME_RECIPE_TOML = """
    schema_version = 1
    kind = "control"
    [system]
    kind = "template"
    template = "TransmonSystem"
    params = { levels = 3, drive_bounds = [0.02, 0.02] }
    [goal]
    kind = "unitary"
    gate = "X"
    [pulse]
    kind = "cubic_spline"
    T = 10.0
    [problem]
    template = "SplinePulseProblem"
    N = 11
    goal_treatment = "both"
    final_fidelity = 0.999
    free_dt = [0.05, 0.5]
    [[problem.objectives]]
    kind = "time"
    weight = 100.0
    """
    MINTIME_FACTORED_TOML = """
    schema_version = 1
    kind = "control"
    [goal]
    kind = "unitary"
    gate = "X"
    [system]
    kind = "template"
    template = "TransmonSystem"
    params = { levels = 3, drive_bounds = [0.02, 0.02] }
    [problem]
    template = "SplinePulseProblem"
    N = 11
    Q = 100.0
    R = 0.01
    goal_treatment = "both"
    final_fidelity = 0.999
    free_dt = [0.05, 0.5]
    [[problem.objectives]]
    kind = "time"
    weight = 100.0
    [pulse]
    kind = "cubic_spline"
    T = 10.0
    """
    a = Specs.parse_spec(MINTIME_RECIPE_TOML; format = :toml)
    b = Specs.parse_spec(MINTIME_FACTORED_TOML; format = :toml)
    @test Specs.canonical_json(Specs.full_dict(a)) ==
          Specs.canonical_json(Specs.full_dict(b))
    @test Specs.structure_hash(a) == Specs.structure_hash(b)

    # The min_time recipe materializes to the legacy MinimumTimeProblem NLP:
    # infidelity objective retained + MinimumTimeObjective + final-fidelity constraint.
    mt = Specs.materialize(a)
    @test any(o -> o isa Piccolo.MinimumTimeObjective, _obj_terms(mt.prob.objective))
    @test any(c -> occursin("Fidelity", string(typeof(c))), mt.prob.constraints)

    # constraint-only builds a Final*FidelityConstraint, no infidelity objective.
    CONSTRAINT_ONLY_TOML = """
    schema_version = 1
    kind = "control"
    [system]
    kind = "template"
    template = "TransmonSystem"
    params = { levels = 3, drive_bounds = [0.02, 0.02] }
    [goal]
    kind = "unitary"
    gate = "X"
    [pulse]
    kind = "cubic_spline"
    T = 10.0
    [problem]
    template = "SplinePulseProblem"
    N = 11
    goal_treatment = "constraint"
    final_fidelity = 0.99
    """
    c = Specs.materialize(Specs.parse_spec(CONSTRAINT_ONLY_TOML; format = :toml))
    @test any(cn -> occursin("Fidelity", string(typeof(cn))), c.prob.constraints)
    # infidelity objective (the terminal-knot state KnotPointObjective) is dropped.
    @test !any(
        o -> o isa DirectTrajOpt.KnotPointObjective && o.times == [11],
        _obj_terms(c.prob.objective),
    )
end

@testitem "materialize validations + sampling wrapper" begin
    using Piccolo, Piccolo.Specs

    _SYS_GOAL_PULSE = """
    schema_version = 1
    kind = "control"
    [system]
    kind = "template"
    template = "TransmonSystem"
    params = { levels = 3, drive_bounds = [0.02, 0.02] }
    [goal]
    kind = "unitary"
    gate = "X"
    [pulse]
    kind = "cubic_spline"
    T = 10.0
    """

    # free_phase without a non-bilinear integrator → structured validation error.
    # (exponential/spline integrators live in Piccolissimo, not loaded here.)
    FREE_PHASE_BILINEAR_TOML = _SYS_GOAL_PULSE * """
    [problem]
    template = "SplinePulseProblem"
    N = 11
    free_phase = true
    """
    @test_throws Specs.SpecValidationError Specs.materialize(
        Specs.parse_spec(FREE_PHASE_BILINEAR_TOML; format = :toml),
    )

    # R_ddu on SplinePulseProblem (which has no ddu) → structured error.
    SPLINE_WITH_RDDU_TOML = _SYS_GOAL_PULSE * """
    [problem]
    template = "SplinePulseProblem"
    N = 11
    R_ddu = 0.01
    """
    @test_throws Specs.SpecValidationError Specs.materialize(
        Specs.parse_spec(SPLINE_WITH_RDDU_TOML; format = :toml),
    )

    # sampling wrapper builds a SamplingProblem-backed qcp over 2 system variants.
    SAMPLING_TOML = """
    schema_version = 1
    kind = "control"
    [system]
    kind = "template"
    template = "TransmonSystem"
    params = { levels = 3, drive_bounds = [0.02, 0.02], "δ" = 0.2 }
    [goal]
    kind = "unitary"
    gate = "X"
    [pulse]
    kind = "cubic_spline"
    T = 10.0
    [problem]
    template = "SplinePulseProblem"
    N = 11
    [[wrappers]]
    kind = "sampling"
    variants = [ { "δ" = 0.19 }, { "δ" = 0.21 } ]
    """
    s = Specs.materialize(Specs.parse_spec(SAMPLING_TOML; format = :toml))
    @test s isa QuantumControlProblem
    @test length(Specs.get_variants(s)) == 2   # 2 system variants

    # robust wrapper is schema-only in Phase 1 → structured "deferred" error.
    ROBUST_TOML = _SYS_GOAL_PULSE * """
    [problem]
    template = "SplinePulseProblem"
    N = 11
    [[wrappers]]
    kind = "robust"
    """
    @test_throws Specs.SpecValidationError Specs.materialize(
        Specs.parse_spec(ROBUST_TOML; format = :toml),
    )
end
