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
    # Retain the originating spec on the materialized object. `extract_spec` returns
    # it only after verifying it against the object, which is what turns the spec
    # round-trip into a real consistency check on this function.
    Control.retain_spec!(qcp, spec)
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

"""
    _template_kwargs(spec::ProblemSpec) -> Dict{Symbol,Any}

The template keyword arguments a `control` spec maps to — the single mapping from
wire fields to template keywords. `_call_template` splats it; `extract_spec` replays
it to verify the retained params against the spec.

Runtime-object keywords (`integrator`, `piccolo_options`) are added by
`_call_template`, not here: they are never spec-expressible and never retained.
"""
function _template_kwargs(spec::ProblemSpec)
    p = spec.problem
    kwargs = Dict{Symbol,Any}(:Q => p.Q, :R => p.R, :free_phase => p.free_phase)
    p.R_u === nothing || (kwargs[:R_u] = p.R_u)
    p.R_du === nothing || (kwargs[:R_du] = p.R_du)
    p.R_ddu === nothing || (kwargs[:R_ddu] = p.R_ddu)
    isfinite(p.du_bound) && (kwargs[:du_bound] = p.du_bound)
    p.ddu_bound === nothing || (kwargs[:ddu_bound] = Float64(p.ddu_bound))
    p.initial_phases === nothing || (kwargs[:initial_phases] = p.initial_phases)
    # free_phase for ket goals: the templates assert `subsystem_levels` is set
    # (a KetTrajectory goal has no EmbeddedOperator to infer them from), so
    # thread the goal's declaration into the template call (#297).
    spec.goal.subsystem_levels === nothing ||
        (kwargs[:subsystem_levels] = spec.goal.subsystem_levels)
    isempty(p.calibration_targets) || (kwargs[:calibration_targets] = p.calibration_targets)
    isempty(p.global_names) || (kwargs[:global_names] = p.global_names)
    isempty(p.global_bounds) ||
        (kwargs[:global_bounds] = _coerce_global_bounds(p.global_bounds))
    p.free_dt isa Free && (kwargs[:Δt_bounds] = (p.free_dt.lo, p.free_dt.hi))
    return kwargs
end

function _call_template(spec::ProblemSpec, qtraj; piccolo_options = nothing)
    p = spec.problem
    fac = lookup_template(p.template).factory
    kwargs = _template_kwargs(spec)
    intg = _build_integrator(spec, qtraj)
    intg === nothing || (kwargs[:integrator] = intg)
    # The spec always names its integrator (schema-level field), so a bilinear
    # integrator is a DECLARED choice, never a silent template default — pass it
    # explicitly so spline-template pulse-type guards (#275) treat it as such.
    if intg === nothing && p.template === :SplinePulseProblem
        kwargs[:integrator_type] = :pwc
    end
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
    # `with_problem`, not a bare `QuantumControlProblem(...)`: the composition axes
    # rewrite the NLP but must not flatten the problem back to the untagged flavor —
    # the template tag and retained params are part of its identity.
    return Control.with_problem(qcp, qtraj, new_prob)
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

    # `goal` and `pulse` default to `nothing` on ProblemSpec (and the JSON Schema's
    # control branch requires only `system`), but materializing a control spec needs
    # both: `_build_goal` has no `Nothing` method, and `_build_pulse_trajectory`
    # reads `spec.pulse.T` directly. Omitting either used to surface as a
    # MethodError/FieldError — an internal error with no field path, from the one
    # code path whose stated contract is "never a materialization crash".
    if spec.goal === nothing
        push!(errs, SpecError("goal", "missing [goal] block for a control spec"))
    end
    if spec.pulse === nothing
        push!(errs, SpecError("pulse", "missing [pulse] block for a control spec"))
    end

    # system availability / composite deferral
    if spec.system.kind === :template
        # `template` is Union{Nothing,Symbol}, and `lookup_system` only accepts a
        # Symbol — so the omitted-field case has to be caught before the lookup,
        # not by it.
        if spec.system.template === nothing
            push!(
                errs,
                SpecError(
                    "system.template",
                    "missing `template` for a system with `kind = \"template\"`",
                    nothing,
                    sort!(String[string(k) for k in keys(SYSTEMS)]),
                ),
            )
        elseif lookup_system(spec.system.template) === nothing
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
    get_variants(qcp::AbstractQuantumControlProblem) -> Vector

The system variants a materialized problem optimizes over: the sampled `systems`
for a `SamplingProblem`-backed problem, else the single nominal system.
"""
get_variants(qcp::AbstractQuantumControlProblem) =
    Control.quantum_trajectory(qcp) isa Quantum.SamplingTrajectory ?
    Control.quantum_trajectory(qcp).systems : [get_system(qcp)]

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
    ref = SplinePulseProblem(traj, 11; integrator_type = :pwc)

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
    # Assert the constraint's TYPE, not its type *name*. `FinalUnitaryFidelityConstraint` is a
    # function returning a `NonlinearKnotPointConstraint` over a closure, and whether that closure's
    # type name carries the enclosing function's name is a Julia-version detail: 1.12 produces
    # `var"#FinalUnitaryFidelityConstraint##0#..."`, while 1.10/1.11 produce `var"#41#42"`. An
    # `occursin("Fidelity", string(typeof(c)))` check therefore passed locally on 1.12 and failed on
    # CI's 1.10/1.11 — it was testing the closure-naming convention, not the materialization.
    @test any(c -> c isa DirectTrajOpt.NonlinearKnotPointConstraint, mt.prob.constraints)

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
    # type, not type name — see the note above
    @test any(cn -> cn isa DirectTrajOpt.NonlinearKnotPointConstraint, c.prob.constraints)
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
    # the [[wrappers]] list is mirrored in the TYPE: a sampled spline problem
    @test s isa SamplingProblem
    @test s isa AbstractQuantumControlProblem
    @test inner(s) isa SplinePulseProblem
    @test template_tag(s) === SplinePulseTemplate()
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

@testitem "materialize: omitted required blocks are structured errors, not crashes" begin
    using Piccolo, Piccolo.Specs

    # `_validate!`'s contract is that materialize is "never a materialization
    # crash". Each omission below reaches code accepting only a non-nothing value,
    # so without a pre-pass guard it surfaced as a MethodError/FieldError — an
    # internal error with no field path, from the one entry point whose job is to
    # hand a caller something it can correct. `goal`/`system.template` were the two
    # JET flagged (`_build_goal(::Nothing, …)`, `lookup_system(::Nothing)`); `pulse`
    # is the same defect one field over, which JET did not reach.

    # Each case below isolates ONE omission (supplying the others) so a failure
    # names the guard that regressed; the last case checks they collect.
    _GOAL = """
    [goal]
    kind = "unitary"
    gate = "X"
    """
    _PULSE = """
    [pulse]
    kind = "linear_spline"
    T = 10.0
    """
    _PROBLEM = """
    [problem]
    template = "SplinePulseProblem"
    N = 11
    """
    _SYS = """
    schema_version = 1
    kind = "control"
    [system]
    kind = "template"
    template = "TransmonSystem"
    """
    _SYS_NO_TEMPLATE = """
    schema_version = 1
    kind = "control"
    [system]
    kind = "template"
    """
    paths_of(toml) = begin
        err = try
            Specs.materialize(Specs.parse_spec(toml; format = :toml))
            nothing
        catch e
            e
        end
        @test err isa Specs.SpecValidationError
        err
    end

    # [goal] omitted. It parses fine — the schema's control branch requires only
    # `system` — so this is purely a materialize-time contract.
    e_goal = paths_of(_SYS * _PULSE * _PROBLEM)
    @test "goal" in [e.path for e in e_goal.errors]

    # [pulse] omitted: `_build_pulse_trajectory` reads `spec.pulse.T` directly.
    e_pulse = paths_of(_SYS * _GOAL * _PROBLEM)
    @test "pulse" in [e.path for e in e_pulse.errors]

    # `kind = "template"` with no `template` key. This one had to be caught AHEAD of
    # `lookup_system`, not by it — `lookup_system` accepts only a Symbol, so the
    # existing `lookup_system(...) === nothing` check was itself the crash site.
    e_tmpl = paths_of(_SYS_NO_TEMPLATE * _GOAL * _PULSE * _PROBLEM)
    st = only(filter(e -> e.path == "system.template", e_tmpl.errors))
    @test occursin("missing", st.msg)
    @test st.allowed !== nothing          # still lists the registered systems
    @test "TransmonSystem" in st.allowed

    # A known system template still validates — the guard must not swallow the
    # ordinary path (and an unknown one still reports as unknown, not missing).
    e_unknown = paths_of(
        replace(_SYS, "TransmonSystem" => "NoSuchSystem") * _GOAL * _PULSE * _PROBLEM,
    )
    su = only(filter(e -> e.path == "system.template", e_unknown.errors))
    @test occursin("unknown", su.msg)

    # All three omissions at once are COLLECTED, not short-circuited — the
    # documented behaviour is that a caller fixes every field from one payload.
    e_all = paths_of(_SYS_NO_TEMPLATE * _PROBLEM)
    paths = [e.path for e in e_all.errors]
    @test "goal" in paths
    @test "pulse" in paths
    @test "system.template" in paths
end

@testitem "materialize: raw system + SmoothPulseProblem + zero-order pulse" begin
    using Piccolo, Piccolo.Specs

    # kind = "raw": explicit H_drift/H_drives as nested [re, im] pairs, default
    # and explicit drive_bounds, zero-order pulse on the SmoothPulseProblem path.
    RAW_TOML = """
    schema_version = 1
    kind = "control"
    [system]
    kind = "raw"
    H_drift = [[[0.0, 0.0], [0.0, 0.0]], [[0.0, 0.0], [0.01, 0.0]]]
    H_drives = [[[[0.0, 0.0], [1.0, 0.0]], [[1.0, 0.0], [0.0, 0.0]]]]
    params = { drive_bounds = [0.5] }
    [goal]
    kind = "unitary"
    gate = "X"
    [pulse]
    kind = "zero_order"
    T = 10.0
    [problem]
    template = "SmoothPulseProblem"
    N = 11
    """
    qcp = Specs.materialize(Specs.parse_spec(RAW_TOML; format = :toml))
    @test qcp isa QuantumControlProblem
    @test get_system(qcp).H(zeros(1), 0.0) ≈ ComplexF64[0.0 0.0; 0.0 0.01]
    @test get_system(qcp).n_drives == 1

    # omitted params → drive_bounds defaults to fill(1.0, n_drives)
    qcp2 = Specs.materialize(
        Specs.parse_spec(
            replace(RAW_TOML, "params = { drive_bounds = [0.5] }\n" => "");
            format = :toml,
        ),
    )
    @test qcp2 isa QuantumControlProblem

    # kind = "composite" is deferred — a structured validation error, never a crash
    e = try
        Specs.materialize(
            Specs.parse_spec(
                replace(RAW_TOML, "kind = \"raw\"" => "kind = \"composite\"");
                format = :toml,
            ),
        )
        nothing
    catch e
        e
    end
    @test e isa Specs.SpecValidationError
    @test "system.kind" in [x.path for x in e.errors]
end

@testitem "materialize: MultiTransmonSystem spec builds a composite problem" begin
    using Piccolo, Piccolo.Specs

    # MultiTransmonSystem takes positional (ωs, δs, gs); the registry entry wraps
    # them into a kwargs-only factory (#297). TOML carries `gs` as nested rows.
    MULTI_TOML = """
    schema_version = 1
    kind = "control"
    [system]
    kind = "template"
    template = "MultiTransmonSystem"
    params = { "ωs" = [4.0, 4.1], "δs" = [0.2, 0.21], "gs" = [[0.0, 0.01], [0.01, 0.0]], subsystem_levels = [2, 2] }
    [goal]
    kind = "unitary"
    gate = "CZ"
    [pulse]
    kind = "linear_spline"
    T = 10.0
    [problem]
    template = "SplinePulseProblem"
    N = 11
    """
    spec = Specs.parse_spec(MULTI_TOML; format = :toml)

    # the factory path alone: a CompositeQuantumSystem of 2 transmons × 2 levels
    sys = Specs._build_system(spec.system)
    @test sys isa CompositeQuantumSystem
    @test sys.subsystem_levels == [2, 2]
    @test sys.levels == 4
    @test sys.n_drives == 4            # 2 drive quadratures per transmon

    # and the full problem materializes: CZ embeds on the qubit subspace of the
    # composite system (the `_build_goal` CompositeQuantumSystem branch).
    qcp = Specs.materialize(spec)
    @test qcp isa QuantumControlProblem
    @test get_system(qcp).levels == 4
    @test qcp.qtraj.goal isa EmbeddedOperator
end

@testitem "materialize: ket goals, subspace goals, random pulse init" begin
    using Piccolo, Piccolo.Specs

    KET_TOML = """
    schema_version = 1
    kind = "control"
    [system]
    kind = "template"
    template = "TransmonSystem"
    params = { levels = 3, drive_bounds = [0.02, 0.02] }
    [goal]
    kind = "ket"
    target = "[0.0, 1.0, 0.0]"
    [pulse]
    kind = "linear_spline"
    T = 10.0
    [problem]
    template = "SplinePulseProblem"
    N = 11
    """

    # default initial state (first basis vector), parsed target ket
    qcp = Specs.materialize(Specs.parse_spec(KET_TOML; format = :toml))
    @test qcp.qtraj isa KetTrajectory
    @test qcp.qtraj.initial ≈ [1.0, 0.0, 0.0]
    @test qcp.qtraj.goal ≈ [0.0, 1.0, 0.0]

    # explicit initial: also a Julia-parseable complex vector
    toml = replace(
        KET_TOML,
        "target = \"[0.0, 1.0, 0.0]\"" => "target = \"[0.0, 1.0, 0.0]\"\ninitial = \"[0.5, 0.5, 0.0]\"",
    )
    qcp2 = Specs.materialize(Specs.parse_spec(toml; format = :toml))
    @test qcp2.qtraj.initial ≈ [0.5, 0.5, 0.0]

    # ket goal without a target: structured error naming goal.target
    e = try
        Specs.materialize(
            Specs.parse_spec(
                replace(KET_TOML, "target = \"[0.0, 1.0, 0.0]\"" => "");
                format = :toml,
            ),
        )
        nothing
    catch e
        e
    end
    @test e isa Specs.SpecValidationError
    @test "goal.target" in [x.path for x in e.errors]

    # unitary goal with a subspace: EmbeddedOperator on the computational subspace
    SUBSPACE_TOML = """
    schema_version = 1
    kind = "control"
    [system]
    kind = "template"
    template = "TransmonSystem"
    params = { levels = 3, drive_bounds = [0.02, 0.02] }
    [goal]
    kind = "unitary"
    gate = "X"
    subspace = [[1, 2]]
    [pulse]
    kind = "linear_spline"
    T = 10.0
    init = "random"
    seed = 42
    [problem]
    template = "SplinePulseProblem"
    N = 11
    """
    qcp3 = Specs.materialize(Specs.parse_spec(SUBSPACE_TOML; format = :toml))
    @test qcp3.qtraj.goal isa EmbeddedOperator
    @test qcp3.qtraj.goal.subspace == [1, 2]

    # init = "random" seeds a deterministic, nonzero control guess
    @test any(!=(0.0), get_trajectory(qcp3).u)
    qcp3_again = Specs.materialize(Specs.parse_spec(SUBSPACE_TOML; format = :toml))
    @test get_trajectory(qcp3_again).u == get_trajectory(qcp3).u  # same seed ⇒ same guess
end

@testitem "materialize: _concretize, _coerce_global_bounds, composite _build_system" begin
    using Piccolo, Piccolo.Specs

    # _concretize narrows boxed TOML/JSON vectors to concrete element types,
    # leaving non-numeric and empty vectors untouched.
    @test Specs._concretize([1, 2, 3]) == [1, 2, 3]
    @test Specs._concretize([1, 2, 3]) isa Vector{Int}
    @test Specs._concretize([0.02, 0.02]) isa Vector{Float64}
    @test Specs._concretize(["a", "b"]) == ["a", "b"]
    @test Specs._concretize([]) == []
    @test Specs._concretize(3.0) == 3.0
    @test Specs._concretize([1, 2.5]) == [1, 2.5]        # mixed → untouched
    @test Specs._concretize_params(Dict(:b => [1, 2]))[:b] isa Vector{Int}

    # _coerce_global_bounds: vector / scalar / tuple value forms
    out = Specs._coerce_global_bounds(
        Dict(:vec => [0.1, 0.9], :scl => 0.5, :tup => (0.2, 0.8)),
    )
    @test out[:vec] == (0.1, 0.9)
    @test out[:scl] == 0.5
    @test out[:tup] == (0.2, 0.8)

    # :composite reaches _build_system only via a direct call (validation blocks
    # it first on the materialize path) — the defensive throw is still SpecValidationError.
    e = try
        Specs._build_system(Specs.SystemSpec(; kind = :composite))
        nothing
    catch e
        e
    end
    @test e isa Specs.SpecValidationError

    # _build_goal's CompositeQuantumSystem branch: the goal embeds into each
    # subsystem's computational subspace (unreachable via registered spec
    # templates — MultiTransmonSystem takes positional args, so the specs
    # kwargs-only factory cannot build one).
    σx = ComplexF64[0 1; 1 0]
    σz = ComplexF64[1 0; 0 -1]
    sub1 = QuantumSystem(0.01 * σz, [σx], [1.0])
    sub2 = QuantumSystem(0.02 * σz, [σx], [1.0])
    comp = CompositeQuantumSystem(
        0.01 * kron(σx, σx),
        Matrix{ComplexF64}[],
        [sub1, sub2],
        [1.0, 1.0],
    )
    goal = Specs._build_goal(Specs.GoalSpec(; kind = :unitary, gate = :CZ), comp)
    @test goal isa EmbeddedOperator

    # _build_goal / _build_pulse_trajectory defensive throws for an unknown
    # goal kind — validation's trajectory-compat check blocks these on the
    # materialize path (every registered template rejects :banana), so they are
    # exercised directly.
    e = try
        Specs._build_goal(Specs.GoalSpec(; kind = :banana), sub1)
        nothing
    catch e
        e
    end
    @test e isa Specs.SpecValidationError

    spec = Specs.ProblemSpec(;
        system = Specs.SystemSpec(; kind = :template, template = :TransmonSystem),
        goal = Specs.GoalSpec(; kind = :banana),
        pulse = Specs.PulseSpec(; kind = :linear_spline, T = 10.0),
        problem = Specs.TemplateBlock(; template = :SplinePulseProblem, N = 11),
    )
    e = try
        Specs._build_pulse_trajectory(spec, sub1, nothing)
        nothing
    catch e
        e
    end
    @test e isa Specs.SpecValidationError
end

@testitem "materialize: integrator registration fallbacks via injected registry entry" begin
    using Piccolo, Piccolo.Specs
    using Piccolo: BilinearIntegrator

    # `exponential`/`spline` integrators live in Piccolissimo; without it the
    # registry only has the `bilinear` sentinel. Inject a stand-in entry to
    # exercise the non-bilinear branches (_build_integrator, the free_phase
    # gating pass, _sampling_integrator_factory, and the template integrator
    # kwarg pass-through), then remove it — the registries are process-global
    # and shared across test items on a worker (see registries.jl testitem).
    Specs.register_integrator!(
        :probe_integrator,
        Specs.RegistryEntry(;
            factory = (qtraj, N; alg = :probe) -> BilinearIntegrator(qtraj, N),
        ),
    )

    try
        BASE = """
        schema_version = 1
        kind = "control"
        [system]
        kind = "template"
        template = "TransmonSystem"
        params = { levels = 3, drive_bounds = [0.02, 0.02] }
        [goal]
        kind = "unitary"
        gate = "X"
        subspace = [[1, 2]]
        [pulse]
        kind = "linear_spline"
        T = 10.0
        [problem]
        template = "SplinePulseProblem"
        N = 11
        [integrator]
        kind = "probe_integrator"
        """

        # the declared integrator is built and passed to the template explicitly
        qcp = Specs.materialize(Specs.parse_spec(BASE; format = :toml))
        @test qcp isa QuantumControlProblem
        @test any(i -> i isa BilinearIntegrator, qcp.prob.integrators)

        # free_phase is gated on a non-bilinear integrator kind: with the
        # injected entry it passes validation, and global_bounds coerce into the
        # free-phase φ variable's bounds.
        for (gb_form, expected) in (
            ("global_bounds = { \"φ_1\" = [0.1, 0.9] }" => (0.1, 0.9)),
            ("global_bounds = { \"φ_1\" = 0.5 }" => (-0.5, 0.5)),
        )
            toml = replace(BASE, "[problem]" => "[problem]\nfree_phase = true\n$gb_form")
            qcp = Specs.materialize(Specs.parse_spec(toml; format = :toml))
            @test haskey(get_trajectory(qcp).global_components, :φ_1)
        end

        # the sampling wrapper resolves the same declared integrator into a
        # factory for SamplingProblem (bilinear → nothing is the other branch).
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
        kind = "linear_spline"
        T = 10.0
        [problem]
        template = "SplinePulseProblem"
        N = 11
        [integrator]
        kind = "probe_integrator"
        [[wrappers]]
        kind = "sampling"
        variants = [ { "δ" = 0.19 }, { "δ" = 0.21 } ]
        """
        s = Specs.materialize(Specs.parse_spec(SAMPLING_TOML; format = :toml))
        @test length(Specs.get_variants(s)) == 2

        # _sampling_integrator_factory: nothing for bilinear / absent kinds
        bilinear_spec = Specs.parse_spec(
            replace(SAMPLING_TOML, "kind = \"probe_integrator\"" => "kind = \"bilinear\"");
            format = :toml,
        )
        @test Specs._sampling_integrator_factory(bilinear_spec) === nothing
        no_int_spec = Specs.parse_spec(
            replace(SAMPLING_TOML, "[integrator]\nkind = \"probe_integrator\"\n" => "");
            format = :toml,
        )
        @test Specs._sampling_integrator_factory(no_int_spec) === nothing

        # ket + free_phase on SmoothPulseProblem (compat ket_free_phase = false):
        # reachable only with a non-bilinear integrator kind declared, and then
        # still a structured validation error.
        KET_SMOOTH_TOML = """
        schema_version = 1
        kind = "control"
        [system]
        kind = "template"
        template = "TransmonSystem"
        params = { levels = 3, drive_bounds = [0.02, 0.02] }
        [goal]
        kind = "ket"
        target = "[0.0, 1.0, 0.0]"
        [pulse]
        kind = "zero_order"
        T = 10.0
        [problem]
        template = "SmoothPulseProblem"
        N = 11
        free_phase = true
        [integrator]
        kind = "probe_integrator"
        """
        e = try
            Specs.materialize(Specs.parse_spec(KET_SMOOTH_TOML; format = :toml))
            nothing
        catch e
            e
        end
        @test e isa Specs.SpecValidationError
        @test "problem.free_phase" in [x.path for x in e.errors]
    finally
        delete!(Specs.INTEGRATORS, :probe_integrator)
    end
end

@testitem "materialize: ket + free_phase threads goal.subsystem_levels" begin
    using Piccolo, Piccolo.Specs
    using Piccolo: BilinearIntegrator

    # free_phase validation requires a non-bilinear integrator KIND;
    # exponential/spline live in Piccolissimo, so inject the same
    # BilinearIntegrator-backed stand-in the integrator-fallback test uses.
    Specs.register_integrator!(
        :probe_integrator,
        Specs.RegistryEntry(;
            factory = (qtraj, N; alg = :probe) -> BilinearIntegrator(qtraj, N),
        ),
    )
    try
        # Before #297 this spec passed every compat check (SplinePulseProblem,
        # ket goal, free_phase, non-bilinear integrator kind) and then crashed
        # on the template's `@assert !isnothing(subsystem_levels)` —
        # `_call_template` never threaded the goal's subsystem_levels through.
        KET_FREEPHASE_TOML = """
        schema_version = 1
        kind = "control"
        [system]
        kind = "template"
        template = "TransmonSystem"
        params = { levels = 3, drive_bounds = [0.02, 0.02] }
        [goal]
        kind = "ket"
        target = "[0.0, 1.0, 0.0]"
        subsystem_levels = [3]
        [pulse]
        kind = "linear_spline"
        T = 10.0
        [problem]
        template = "SplinePulseProblem"
        N = 11
        free_phase = true
        [integrator]
        kind = "probe_integrator"
        """
        qcp = Specs.materialize(Specs.parse_spec(KET_FREEPHASE_TOML; format = :toml))
        @test qcp isa QuantumControlProblem
        @test qcp.qtraj isa KetTrajectory
        # one free-phase global per subsystem (φ_1 for the single [3] subsystem)
        @test haskey(get_trajectory(qcp).global_components, :φ_1)
    finally
        delete!(Specs.INTEGRATORS, :probe_integrator)
    end
end

@testitem "materialize: objective term kinds build into the composed problem" begin
    using Piccolo, Piccolo.Specs

    _obj_terms(J) =
        J isa DirectTrajOpt.CompositeObjective ? collect(J.objectives) :
        AbstractObjective[J]

    BASE = """
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
    kind = "linear_spline"
    T = 10.0
    [problem]
    template = "SplinePulseProblem"
    N = 11
    """

    # reg_u: an extra QuadraticRegularizer on the drive variable
    qcp = Specs.materialize(
        Specs.parse_spec(
            BASE * "\n[[problem.objectives]]\nkind = \"reg_u\"\nweight = 0.1\n";
            format = :toml,
        ),
    )
    @test count(
        o -> o isa DirectTrajOpt.QuadraticRegularizer && o.name == :u,
        _obj_terms(qcp.prob.objective),
    ) == 2

    # reg_ddu: only valid on SmoothPulseProblem (the ddu-carrying template)
    smooth = replace(
        BASE,
        "linear_spline" => "zero_order",
        "SplinePulseProblem" => "SmoothPulseProblem",
    )
    qcp2 = Specs.materialize(
        Specs.parse_spec(
            smooth * "\n[[problem.objectives]]\nkind = \"reg_ddu\"\nweight = 0.1\n";
            format = :toml,
        ),
    )
    @test qcp2 isa QuantumControlProblem

    # leakage: a LeakageObjective is a full-horizon (times 1:N) KnotPointObjective
    # on the state variable — structurally distinct from the terminal-knot
    # infidelity objective and from the QuadraticRegularizer terms.
    qcp3 = Specs.materialize(
        Specs.parse_spec(
            BASE * "\n[[problem.objectives]]\nkind = \"leakage\"\nweight = 10.0\n";
            format = :toml,
        ),
    )
    leaky = filter(_obj_terms(qcp3.prob.objective)) do o
        o isa DirectTrajOpt.KnotPointObjective &&
            o.var_names == [state_name(qcp3.qtraj)] &&
            o.times == 1:11
    end
    @test length(leaky) == 1

    # sensitivity: a UnitarySensitivityObjective is a terminal-knot (times [N])
    # KnotPointObjective — a second one beyond the built-in infidelity term.
    qcp4 = Specs.materialize(
        Specs.parse_spec(
            BASE * "\n[[problem.objectives]]\nkind = \"sensitivity\"\nweight = 10.0\n";
            format = :toml,
        ),
    )
    terminal = count(
        o ->
            o isa DirectTrajOpt.KnotPointObjective &&
            o.var_names == [state_name(qcp4.qtraj)] &&
            o.times == [11],
        _obj_terms(qcp4.prob.objective),
    )
    @test terminal == 2

    # a REGISTERED but unsupported term kind is a structured materialize error,
    # not a silent skip or a MethodError (inject one, then remove it).
    Specs.register_objective_term!(
        :probe_term,
        Specs.RegistryEntry(; factory = Specs.ConstFactory(:probe_term)),
    )
    try
        e = try
            Specs.materialize(
                Specs.parse_spec(
                    BASE *
                    "\n[[problem.objectives]]\nkind = \"probe_term\"\nweight = 1.0\n";
                    format = :toml,
                ),
            )
            nothing
        catch e
            e
        end
        @test e isa Specs.SpecValidationError
        @test "problem.objectives" in [x.path for x in e.errors]
    finally
        delete!(Specs.OBJECTIVE_TERMS, :probe_term)
    end
end

@testitem "materialize: validation branch matrix" begin
    using Piccolo, Piccolo.Specs

    BASE = """
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
    kind = "linear_spline"
    T = 10.0
    [problem]
    template = "SplinePulseProblem"
    N = 11
    """

    # (name, toml, expected error path)
    cases = [
        ("missing [problem] block", replace(BASE, r"\[problem\][^\0]*\z" => ""), "problem"),
        (
            "unknown template",
            replace(BASE, "SplinePulseProblem" => "NoSuchTemplate"),
            "problem.template",
        ),
        (
            "unregistered integrator kind",
            BASE * "\n[integrator]\nkind = \"exponential\"\n",
            "integrator.kind",
        ),
        (
            "unknown goal kind",
            replace(BASE, "kind = \"unitary\"" => "kind = \"banana\""),
            "goal.kind",
        ),
        (
            "pulse kind incompatible",
            replace(BASE, "linear_spline" => "banana_pulse"),
            "pulse.kind",
        ),
        (
            "unknown objective term kind",
            BASE * "\n[[problem.objectives]]\nkind = \"banana\"\n",
            "problem.objectives",
        ),
        (
            "reg_ddu requires SmoothPulseProblem",
            BASE * "\n[[problem.objectives]]\nkind = \"reg_ddu\"\n",
            "problem.objectives",
        ),
        (
            "time objective requires free_dt",
            BASE * "\n[[problem.objectives]]\nkind = \"time\"\n",
            "problem.free_dt",
        ),
        (
            "calibration target not in global_names",
            BASE * "\ncalibration_targets = [\"δ\"]\n",
            "problem.calibration_targets",
        ),
        (
            "unknown wrapper kind",
            BASE * "\n[[wrappers]]\nkind = \"banana\"\n",
            "wrappers[1]",
        ),
    ]
    for (name, toml, path) in cases
        e = try
            Specs.materialize(Specs.parse_spec(toml; format = :toml))
            nothing
        catch e
            e
        end
        @test e isa Specs.SpecValidationError
        @test path in [x.path for x in e.errors]
    end
end
