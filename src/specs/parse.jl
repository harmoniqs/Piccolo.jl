export parse_spec

# ---------------------------------------------------------------------------
# Allowed key sets — mirror the struct field names in spec_structs.jl. Any key
# outside these sets is a strict-mode error (with its dotted field path).
# ---------------------------------------------------------------------------
const _CONTROL_KEYS = Set(["schema_version", "kind", "system", "goal", "pulse",
    "trajectory", "problem", "integrator", "wrappers", "solver", "warm_start"])
const _ROLLOUT_KEYS = Set(["schema_version", "kind", "input_pulse", "system",
    "rollout_kind", "alg", "n_samples", "times", "initial", "report", "referee"])
const _SYSTEM_KEYS = Set(["kind", "template", "params", "global_params",
    "components", "H_drift", "H_drives"])
const _GOAL_KEYS = Set(["kind", "gate", "matrix", "target", "initial",
    "subsystem_levels", "subspace"])
const _PULSE_KEYS = Set(["kind", "T", "init", "seed"])
const _TRAJECTORY_KEYS = Set(["kind"])
const _PROBLEM_KEYS = Set(["template", "N", "goal_treatment", "final_fidelity",
    "free_dt", "Q", "R", "R_u", "R_du", "R_ddu", "du_bound", "ddu_bound",
    "free_phase", "initial_phases", "calibration_targets", "global_names",
    "global_bounds", "objectives", "options"])
const _INTEGRATOR_KEYS = Set(["kind", "alg"])
const _SOLVER_KEYS = Set(["backend", "device", "precision", "max_iter", "tol", "strategy"])
const _WARM_START_KEYS = Set(["catalog_ref", "pulse_hash"])
const _OBJECTIVE_TERM_KEYS = Set(["kind", "weight"])
const _WRAPPER_KEYS = Set(["kind", "variants", "weights"])
const _REPORT_KEYS = Set(["fidelity", "goal", "populations", "leakage"])
const _REFEREE_KEYS = Set(["run", "solve_knots", "solve_integrator", "fidelity_reported"])

# ---------------------------------------------------------------------------
# Normalization + coercion helpers
# ---------------------------------------------------------------------------

# Recursively convert any AbstractDict/AbstractVector (TOML.jl Dicts, JSON3
# Objects/Arrays) into plain Dict{String,Any}/Vector{Any} so downstream builders
# are format-agnostic.
_to_plain(x::AbstractDict) = Dict{String,Any}(string(k) => _to_plain(v) for (k, v) in x)
_to_plain(x::AbstractVector) = Any[_to_plain(v) for v in x]
_to_plain(x) = x

_symdict(d::AbstractDict) = Dict{Symbol,Any}(Symbol(k) => v for (k, v) in d)
_symvec(v) = Symbol[Symbol(x) for x in v]

function _int(x, path, errs)
    if x isa Integer && !(x isa Bool)
        return Int(x)
    elseif x isa AbstractFloat && isinteger(x)
        return Int(x)
    else
        push!(errs, SpecError(path, "expected an integer", x))
        return 0
    end
end

function _float(x, path, errs)
    if x isa Real && !(x isa Bool)
        return Float64(x)
    else
        push!(errs, SpecError(path, "expected a number", x))
        return 0.0
    end
end

function _bool(x, path, errs)
    x isa Bool && return x
    push!(errs, SpecError(path, "expected a boolean", x))
    return false
end

_intvec(v, path, errs) = Int[_int(x, path, errs) for x in v]
_floatvec(v, path, errs) = Float64[_float(x, path, errs) for x in v]
_parse_subspace(v) = Vector{Int}[Int[Int(y) for y in row] for row in v]

function _strict_fields(raw::AbstractDict, allowed, path, errs)
    for k in keys(raw)
        ks = string(k)
        if !(ks in allowed)
            p = isempty(path) ? ks : "$path.$ks"
            push!(errs, SpecError(p, "unknown field", ks, sort!(collect(allowed))))
        end
    end
end

# ---------------------------------------------------------------------------
# Public entrypoints
# ---------------------------------------------------------------------------

"""
    parse_spec(s::AbstractString; format::Symbol=:toml)

Parse a declarative spec from a `:toml` or `:json` string into a
[`ProblemSpec`](@ref) (kind `control`) or [`RolloutSpec`](@ref) (kind
`rollout`). Strict: unknown fields anywhere are rejected. Every validation
problem is collected (not just the first) and reported as a
[`SpecValidationError`](@ref) carrying field-path [`SpecError`](@ref)s.
"""
function parse_spec(s::AbstractString; format::Symbol=:toml)
    raw = if format === :toml
        TOML.parse(s)
    elseif format === :json
        JSON3.read(s)
    else
        throw(SpecValidationError([SpecError("format", "unknown format", format, ["toml", "json"])]))
    end
    return parse_spec(raw)
end

"""
    parse_spec(raw::AbstractDict)

Build a spec from an already-parsed dictionary. Dispatches on the top-level
`kind` (default `"control"`).
"""
function parse_spec(raw::AbstractDict)
    plain = _to_plain(raw)
    errs = SpecError[]
    kind = Symbol(get(plain, "kind", "control"))
    spec = if kind === :rollout
        _parse_rollout(plain, errs)
    elseif kind === :control
        _parse_control(plain, errs)
    else
        push!(errs, SpecError("kind", "unknown spec kind", string(kind), ["control", "rollout"]))
        nothing
    end
    isempty(errs) || throw(SpecValidationError(errs))
    return spec
end

# ---------------------------------------------------------------------------
# control kind
# ---------------------------------------------------------------------------

function _parse_control(raw, errs)
    _strict_fields(raw, _CONTROL_KEYS, "", errs)
    schema_version = _int(get(raw, "schema_version", 1), "schema_version", errs)
    system = if haskey(raw, "system")
        _parse_system(raw["system"], "system", errs)
    else
        push!(errs, SpecError("system", "missing required [system] block"))
        nothing
    end
    goal = haskey(raw, "goal") ? _parse_goal(raw["goal"], "goal", errs) : nothing
    pulse = haskey(raw, "pulse") ? _parse_pulse(raw["pulse"], "pulse", errs) : nothing
    trajectory = haskey(raw, "trajectory") ? _parse_trajectory(raw["trajectory"], "trajectory", errs) : TrajectorySpec()
    problem = haskey(raw, "problem") ? _parse_problem(raw["problem"], "problem", errs) : nothing
    integrator = haskey(raw, "integrator") ? _parse_integrator(raw["integrator"], "integrator", errs) : nothing
    wrappers = haskey(raw, "wrappers") ? _parse_wrappers(raw["wrappers"], "wrappers", errs) : WrapperSpec[]
    solver = haskey(raw, "solver") ? _parse_solver(raw["solver"], "solver", errs) : SolverSpec()
    warm_start = haskey(raw, "warm_start") ? _parse_warm_start(raw["warm_start"], "warm_start", errs) : nothing
    # Never construct with a nothing where a value is required — if anything went
    # wrong, the collected errors will be thrown by the caller.
    (system === nothing || !isempty(errs)) && return nothing
    return ProblemSpec(; schema_version, kind=:control, system, goal, pulse, trajectory,
        problem, integrator, wrappers, solver, warm_start)
end

function _parse_system(raw::AbstractDict, path, errs)
    _strict_fields(raw, _SYSTEM_KEYS, path, errs)
    kind = Symbol(get(raw, "kind", "template"))
    template = haskey(raw, "template") ? Symbol(raw["template"]) : nothing
    params = _symdict(get(raw, "params", Dict{String,Any}()))
    global_params = _symdict(get(raw, "global_params", Dict{String,Any}()))
    components = get(raw, "components", nothing)
    H_drift = get(raw, "H_drift", nothing)
    H_drives = get(raw, "H_drives", nothing)
    return SystemSpec(; kind, template, params, global_params, components, H_drift, H_drives)
end
_parse_system(raw, path, errs) =
    (push!(errs, SpecError(path, "expected a table", raw)); nothing)

function _parse_goal(raw::AbstractDict, path, errs)
    _strict_fields(raw, _GOAL_KEYS, path, errs)
    if !haskey(raw, "kind")
        push!(errs, SpecError("$path.kind", "missing required field"))
        return nothing
    end
    kind = Symbol(raw["kind"])
    gate = haskey(raw, "gate") ? Symbol(raw["gate"]) : nothing
    matrix = get(raw, "matrix", nothing)
    target = haskey(raw, "target") ? String(raw["target"]) : nothing
    initial = haskey(raw, "initial") ? String(raw["initial"]) : nothing
    subsystem_levels = haskey(raw, "subsystem_levels") ?
                       _intvec(raw["subsystem_levels"], "$path.subsystem_levels", errs) : nothing
    subspace = haskey(raw, "subspace") ? _parse_subspace(raw["subspace"]) : nothing
    return GoalSpec(; kind, gate, matrix, target, initial, subsystem_levels, subspace)
end
_parse_goal(raw, path, errs) =
    (push!(errs, SpecError(path, "expected a table", raw)); nothing)

function _parse_pulse(raw::AbstractDict, path, errs)
    _strict_fields(raw, _PULSE_KEYS, path, errs)
    haskey(raw, "kind") || push!(errs, SpecError("$path.kind", "missing required field"))
    haskey(raw, "T") || push!(errs, SpecError("$path.T", "missing required field"))
    (haskey(raw, "kind") && haskey(raw, "T")) || return nothing
    kind = Symbol(raw["kind"])
    T = _float(raw["T"], "$path.T", errs)
    init = Symbol(get(raw, "init", "default"))
    seed = _int(get(raw, "seed", 0), "$path.seed", errs)
    return PulseSpec(; kind, T, init, seed)
end
_parse_pulse(raw, path, errs) =
    (push!(errs, SpecError(path, "expected a table", raw)); nothing)

function _parse_trajectory(raw::AbstractDict, path, errs)
    _strict_fields(raw, _TRAJECTORY_KEYS, path, errs)
    kind = haskey(raw, "kind") ? Symbol(raw["kind"]) : nothing
    return TrajectorySpec(; kind)
end
_parse_trajectory(raw, path, errs) =
    (push!(errs, SpecError(path, "expected a table", raw)); nothing)

function _parse_free_dt(x, path, errs)
    if x === false
        return Fixed()
    elseif x === true
        push!(errs, SpecError(path, "free_dt must be false or [lo, hi], not true", true))
        return Fixed()
    elseif x isa AbstractVector && length(x) == 2 && all(e -> e isa Real && !(e isa Bool), x)
        lo = Float64(x[1])
        hi = Float64(x[2])
        try
            return Free(lo, hi)
        catch
            push!(errs, SpecError(path, "free_dt requires 0 < lo < hi", [lo, hi]))
            return Fixed()
        end
    else
        push!(errs, SpecError(path, "free_dt must be false or [lo, hi]", x))
        return Fixed()
    end
end

function _parse_problem(raw::AbstractDict, path, errs)
    _strict_fields(raw, _PROBLEM_KEYS, path, errs)
    haskey(raw, "template") || push!(errs, SpecError("$path.template", "missing required field"))
    haskey(raw, "N") || push!(errs, SpecError("$path.N", "missing required field"))
    (haskey(raw, "template") && haskey(raw, "N")) || return nothing
    template = Symbol(raw["template"])
    N = _int(raw["N"], "$path.N", errs)
    goal_treatment = Symbol(get(raw, "goal_treatment", "objective"))
    final_fidelity = haskey(raw, "final_fidelity") ? _float(raw["final_fidelity"], "$path.final_fidelity", errs) : nothing
    free_dt = _parse_free_dt(get(raw, "free_dt", false), "$path.free_dt", errs)
    Q = _float(get(raw, "Q", 100.0), "$path.Q", errs)
    R = _float(get(raw, "R", 1e-2), "$path.R", errs)
    R_u = haskey(raw, "R_u") ? _float(raw["R_u"], "$path.R_u", errs) : nothing
    R_du = haskey(raw, "R_du") ? _float(raw["R_du"], "$path.R_du", errs) : nothing
    R_ddu = haskey(raw, "R_ddu") ? _float(raw["R_ddu"], "$path.R_ddu", errs) : nothing
    du_bound = _float(get(raw, "du_bound", Inf), "$path.du_bound", errs)
    ddu_bound = haskey(raw, "ddu_bound") ? _float(raw["ddu_bound"], "$path.ddu_bound", errs) : nothing
    free_phase = _bool(get(raw, "free_phase", false), "$path.free_phase", errs)
    initial_phases = haskey(raw, "initial_phases") ? _floatvec(raw["initial_phases"], "$path.initial_phases", errs) : nothing
    calibration_targets = _symvec(get(raw, "calibration_targets", String[]))
    global_names = _symvec(get(raw, "global_names", String[]))
    global_bounds = _symdict(get(raw, "global_bounds", Dict{String,Any}()))
    objectives = _parse_objectives(get(raw, "objectives", Any[]), "$path.objectives", errs)
    options = _symdict(get(raw, "options", Dict{String,Any}()))
    return TemplateBlock(; template, N, goal_treatment, final_fidelity, free_dt, Q, R,
        R_u, R_du, R_ddu, du_bound, ddu_bound, free_phase, initial_phases,
        calibration_targets, global_names, global_bounds, objectives, options)
end
_parse_problem(raw, path, errs) =
    (push!(errs, SpecError(path, "expected a table", raw)); nothing)

function _parse_objectives(v, path, errs)
    terms = ObjectiveTermSpec[]
    v isa AbstractVector || (push!(errs, SpecError(path, "expected a list of objective terms", v)); return terms)
    for (i, item) in enumerate(v)
        ip = "$path[$i]"
        if !(item isa AbstractDict)
            push!(errs, SpecError(ip, "objective term must be a table", item))
            continue
        end
        _strict_fields(item, _OBJECTIVE_TERM_KEYS, ip, errs)
        if !haskey(item, "kind")
            push!(errs, SpecError("$ip.kind", "missing required field"))
            continue
        end
        kind = Symbol(item["kind"])
        weight = _float(get(item, "weight", 1.0), "$ip.weight", errs)
        push!(terms, ObjectiveTermSpec(; kind, weight))
    end
    return terms
end

function _parse_integrator(raw::AbstractDict, path, errs)
    _strict_fields(raw, _INTEGRATOR_KEYS, path, errs)
    kind = Symbol(get(raw, "kind", "bilinear"))
    alg = Symbol(get(raw, "alg", "tsit5"))
    return IntegratorSpec(; kind, alg)
end
_parse_integrator(raw, path, errs) =
    (push!(errs, SpecError(path, "expected a table", raw)); nothing)

function _parse_wrappers(v, path, errs)
    wraps = WrapperSpec[]
    v isa AbstractVector || (push!(errs, SpecError(path, "expected a list of wrappers", v)); return wraps)
    for (i, item) in enumerate(v)
        ip = "$path[$i]"
        if !(item isa AbstractDict)
            push!(errs, SpecError(ip, "wrapper must be a table", item))
            continue
        end
        _strict_fields(item, _WRAPPER_KEYS, ip, errs)
        if !haskey(item, "kind")
            push!(errs, SpecError("$ip.kind", "missing required field"))
            continue
        end
        kind = Symbol(item["kind"])
        variants = Dict{Symbol,Any}[_symdict(x) for x in get(item, "variants", Any[])]
        weights = haskey(item, "weights") ? _floatvec(item["weights"], "$ip.weights", errs) : nothing
        push!(wraps, WrapperSpec(; kind, variants, weights))
    end
    return wraps
end

function _parse_solver(raw::AbstractDict, path, errs)
    _strict_fields(raw, _SOLVER_KEYS, path, errs)
    backend = Symbol(get(raw, "backend", "ipopt"))
    device = Symbol(get(raw, "device", "cpu"))
    precision = Symbol(get(raw, "precision", "f64"))
    max_iter = _int(get(raw, "max_iter", 500), "$path.max_iter", errs)
    tol = haskey(raw, "tol") ? _float(raw["tol"], "$path.tol", errs) : nothing
    strategy = Symbol(get(raw, "strategy", "direct"))
    return SolverSpec(; backend, device, precision, max_iter, tol, strategy)
end
_parse_solver(raw, path, errs) =
    (push!(errs, SpecError(path, "expected a table", raw)); nothing)

function _parse_warm_start(raw::AbstractDict, path, errs)
    _strict_fields(raw, _WARM_START_KEYS, path, errs)
    catalog_ref = haskey(raw, "catalog_ref") ? String(raw["catalog_ref"]) : nothing
    pulse_hash = haskey(raw, "pulse_hash") ? String(raw["pulse_hash"]) : nothing
    return WarmStartSpec(; catalog_ref, pulse_hash)
end
_parse_warm_start(raw, path, errs) =
    (push!(errs, SpecError(path, "expected a table", raw)); nothing)

# ---------------------------------------------------------------------------
# rollout kind
# ---------------------------------------------------------------------------

function _parse_rollout(raw, errs)
    _strict_fields(raw, _ROLLOUT_KEYS, "", errs)
    schema_version = _int(get(raw, "schema_version", 1), "schema_version", errs)
    haskey(raw, "input_pulse") || push!(errs, SpecError("input_pulse", "missing required field"))
    haskey(raw, "rollout_kind") || push!(errs, SpecError("rollout_kind", "missing required field"))
    system = if haskey(raw, "system")
        _parse_system(raw["system"], "system", errs)
    else
        push!(errs, SpecError("system", "missing required [system] block"))
        nothing
    end
    alg = Symbol(get(raw, "alg", "tsit5"))
    n_samples = haskey(raw, "n_samples") ? _int(raw["n_samples"], "n_samples", errs) : nothing
    times = haskey(raw, "times") ? _floatvec(raw["times"], "times", errs) : nothing
    initial = haskey(raw, "initial") ? String(raw["initial"]) : nothing
    report = haskey(raw, "report") ? _parse_report(raw["report"], "report", errs) : RolloutReportSpec()
    referee = haskey(raw, "referee") ? _parse_referee(raw["referee"], "referee", errs) : nothing
    (system === nothing || !haskey(raw, "input_pulse") || !haskey(raw, "rollout_kind") || !isempty(errs)) && return nothing
    input_pulse = String(raw["input_pulse"])
    rollout_kind = Symbol(raw["rollout_kind"])
    return RolloutSpec(; schema_version, input_pulse, system, rollout_kind, alg,
        n_samples, times, initial, report, referee)
end

function _parse_report(raw::AbstractDict, path, errs)
    _strict_fields(raw, _REPORT_KEYS, path, errs)
    fidelity = _bool(get(raw, "fidelity", true), "$path.fidelity", errs)
    goal = haskey(raw, "goal") ? _parse_goal(raw["goal"], "$path.goal", errs) : nothing
    populations = _bool(get(raw, "populations", false), "$path.populations", errs)
    leakage = _bool(get(raw, "leakage", false), "$path.leakage", errs)
    return RolloutReportSpec(; fidelity, goal, populations, leakage)
end
_parse_report(raw, path, errs) =
    (push!(errs, SpecError(path, "expected a table", raw)); nothing)

function _parse_referee(raw::AbstractDict, path, errs)
    _strict_fields(raw, _REFEREE_KEYS, path, errs)
    required = ("run", "solve_knots", "solve_integrator", "fidelity_reported")
    for f in required
        haskey(raw, f) || push!(errs, SpecError("$path.$f", "missing required field"))
    end
    all(f -> haskey(raw, f), required) || return nothing
    run = String(raw["run"])
    solve_knots = _int(raw["solve_knots"], "$path.solve_knots", errs)
    solve_integrator = Symbol(raw["solve_integrator"])
    fidelity_reported = _float(raw["fidelity_reported"], "$path.fidelity_reported", errs)
    return RefereeSpec(; run, solve_knots, solve_integrator, fidelity_reported)
end
_parse_referee(raw, path, errs) =
    (push!(errs, SpecError(path, "expected a table", raw)); nothing)

@testitem "parse: valid control spec, strict unknown-field rejection" begin
    using Piccolo.Specs
    toml = """
    schema_version = 1
    kind = "control"
    [system]
    kind = "template"
    template = "TransmonSystem"
    params = { levels = 3, drive_bounds = [0.02, 0.02] }
    [goal]
    kind = "unitary"
    gate = "CZ"
    subsystem_levels = [3, 3]
    [pulse]
    kind = "cubic_spline"
    T = 100.0
    [problem]
    template = "SplinePulseProblem"
    N = 100
    """
    spec = Specs.parse_spec(toml; format=:toml)
    @test spec.kind == :control
    @test spec.system.template == :TransmonSystem
    @test spec.goal.gate == :CZ
    @test spec.problem.N == 100

    bad = toml * "\nnonsense_field = true\n"
    @test_throws Specs.SpecValidationError Specs.parse_spec(bad; format=:toml)
end
