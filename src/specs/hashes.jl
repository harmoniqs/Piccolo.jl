export canonical_json, structure_hash, problem_hash, full_dict

# ===========================================================================
# Canonical JSON
# ===========================================================================
#
# Hand-rolled, deterministic JSON serializer designed to produce **byte-identical**
# output to a TypeScript `hashing.ts` (Plan 2 Task 5), so the two languages agree
# on `structure_hash`/`problem_hash`. The rules are pinned here — do NOT swap in
# Julia's native `repr`/`string`/Ryu float formatting (its `100.0`/`1.0e-5` style
# does not match JavaScript's `Number.prototype.toString()`).
#
# Pinned rules
# ------------
# * Objects: keys sorted ascending by Unicode code point. All spec keys are ASCII
#   identifiers, for which code-point order == UTF-8 byte order == JS UTF-16
#   `Array.sort()` order, so Julia `sort` and JS `.sort()` agree.
# * Strings: JSON-escaped, double-quoted (`"`, `\`, and control chars U+0000–
#   U+001F escaped; `\b \t \n \f \r` use their short forms).
# * `true`/`false`/`null` literally.
# * Numbers are **int/float-agnostic**: an integer-valued number renders as a bare
#   integer regardless of whether it arrived as `100` or `100.0` (JS cannot tell
#   `100.0` from `100`; TOML.jl and smol-toml disagree on int-vs-float parsing).
#   Non-integer numbers use the ECMAScript Number-to-String algorithm
#   (ECMA-262 §Number::toString, the same algorithm `JSON.stringify` uses):
#     - shortest round-tripping decimal significand `digits` with decimal point
#       position `n` (value = digits × 10^(n − length(digits)));
#     - `length(digits) ≤ n ≤ 21`  → `digits` then `n−len` zeros;
#     - `0 < n ≤ 21`               → `digits[1:n] "." digits[n+1:end]`;
#     - `−6 < n ≤ 0`               → `"0." (−n zeros) digits`;
#     - otherwise                  → exponential `d[.ddd] "e" ("+"|"-") |n−1|`.
#   Non-finite values are rejected (the spec layer omits `Inf`/`NaN` before here).

"""
    canonical_json(x) -> String

Serialize `x` (nested `Dict`/`Vector`/`Symbol`/`String`/`Bool`/`Real`/`Nothing`)
to canonical JSON per the pinned rules above. Deterministic and intended to be
byte-identical across Julia and TypeScript.
"""
function canonical_json(x)::String
    io = IOBuffer()
    _canon(io, x)
    return String(take!(io))
end

_canon(io::IO, ::Nothing) = print(io, "null")
_canon(io::IO, x::Bool) = print(io, x ? "true" : "false")   # before Integer
_canon(io::IO, x::Symbol) = _canon_string(io, string(x))
_canon(io::IO, x::AbstractString) = _canon_string(io, x)
_canon(io::IO, x::Integer) = print(io, string(x))
_canon(io::IO, x::AbstractFloat) = print(io, _es_number(Float64(x)))

function _canon(io::IO, x::AbstractDict)
    print(io, '{')
    ks = sort!(collect(string(k) for k in keys(x)))
    first = true
    for k in ks
        first || print(io, ',')
        first = false
        _canon_string(io, k)
        print(io, ':')
        _canon(io, x[k])
    end
    print(io, '}')
end

function _canon(io::IO, x::AbstractVector)
    print(io, '[')
    first = true
    for v in x
        first || print(io, ',')
        first = false
        _canon(io, v)
    end
    print(io, ']')
end

function _canon_string(io::IO, s::AbstractString)
    print(io, '"')
    for c in s
        if c == '"'
            print(io, "\\\"")
        elseif c == '\\'
            print(io, "\\\\")
        elseif c == '\b'
            print(io, "\\b")
        elseif c == '\t'
            print(io, "\\t")
        elseif c == '\n'
            print(io, "\\n")
        elseif c == '\f'
            print(io, "\\f")
        elseif c == '\r'
            print(io, "\\r")
        elseif c < '\x20'
            print(io, "\\u", lpad(string(UInt16(c); base = 16), 4, '0'))
        else
            print(io, c)
        end
    end
    print(io, '"')
end

# ECMAScript Number::toString for a finite Float64. Returns integer-valued numbers
# without a decimal point (so 100.0 -> "100"), matching the int/float-agnostic rule.
function _es_number(x::Float64)::String
    isfinite(x) ||
        throw(ArgumentError("canonical_json cannot serialize non-finite number $x"))
    x == 0.0 && return "0"                    # also handles -0.0
    x < 0 && return "-" * _es_number(-x)
    digits, n = _shortest_digits(x)
    k = length(digits)
    if k <= n <= 21
        return digits * "0"^(n - k)
    elseif 0 < n <= 21
        return digits[1:n] * "." * digits[(n+1):end]
    elseif -6 < n <= 0
        return "0." * "0"^(-n) * digits
    else
        mant = k == 1 ? digits : string(digits[1], ".", digits[2:end])
        e = n - 1
        return string(mant, "e", e >= 0 ? "+" : "-", abs(e))
    end
end

# Shortest round-tripping significant digits of a positive finite Float64.
# Returns (digits, n) with value == parse(BigInt, digits) × 10^(n − length(digits)).
# Uses a minimal-precision %e loop for the digit generation (not Julia's Ryu string
# formatting) so the *formatting* is fully controlled by `_es_number` above.
function _shortest_digits(x::Float64)
    local s::String
    for p = 0:17
        s = @sprintf("%.*e", p, x)
        parse(Float64, s) == x && break
    end
    mantissa, expstr = split(s, 'e')
    E = parse(Int, expstr)
    intpart, _, fracpart = partition_mantissa(mantissa)
    digits = intpart * fracpart
    # strip trailing zeros (the minimal-p loop can leave one, e.g. "1.0e+00")
    digits = rstrip_zeros(digits)
    isempty(digits) && (digits = "0")
    n = E + 1
    return digits, n
end

# split "d.ddd" (or "d") into ("d", ".", "ddd")
function partition_mantissa(m::AbstractString)
    if occursin('.', m)
        i = findfirst('.', m)
        return m[1:(i-1)], ".", m[(i+1):end]
    else
        return String(m), "", ""
    end
end

function rstrip_zeros(d::AbstractString)
    i = lastindex(d)
    while i > firstindex(d) && d[i] == '0'
        i = prevind(d, i)
    end
    return String(d[firstindex(d):i])
end

# ===========================================================================
# Wire-form projection (full_dict) — the canonical WIRE layout parse_spec re-accepts
# ===========================================================================

# Recursively map an internal value to its wire form (Symbols -> strings, nested
# Dicts/Vectors preserved with string keys). Scalars pass through.
_wireval(x::Symbol) = string(x)
_wireval(x::AbstractString) = String(x)
_wireval(x::Bool) = x
_wireval(x::Integer) = x
_wireval(x::AbstractFloat) = x
_wireval(x::AbstractVector) = Any[_wireval(v) for v in x]
_wireval(x::AbstractDict) = Dict{String,Any}(string(k) => _wireval(v) for (k, v) in x)
_wireval(x) = x

"""
    full_dict(spec::ProblemSpec) -> Dict{String,Any}

Project a spec into its canonical **wire** form: the nested `[system]`/`[goal]`/…
layout `parse_spec` re-accepts, with Symbols rendered as strings, `FreeDt` as
`false` or `[lo, hi]`, `nothing`/`Inf` fields omitted (parse restores defaults).
`problem_hash` is `sha256(canonical_json(full_dict(spec)))`; the identity
`parse_spec(canonical_json(full_dict(spec)); format=:json)` round-trips.
"""
function full_dict(spec::ProblemSpec)
    d = Dict{String,Any}()
    d["schema_version"] = spec.schema_version
    d["kind"] = "control"
    d["system"] = _system_dict(spec.system)
    spec.goal === nothing || (d["goal"] = _goal_dict(spec.goal))
    spec.pulse === nothing || (d["pulse"] = _pulse_dict(spec.pulse))
    spec.trajectory.kind === nothing ||
        (d["trajectory"] = Dict{String,Any}("kind" => string(spec.trajectory.kind)))
    spec.problem === nothing || (d["problem"] = _problem_dict(spec.problem))
    spec.integrator === nothing || (d["integrator"] = _integrator_dict(spec.integrator))
    isempty(spec.wrappers) || (d["wrappers"] = Any[_wrapper_dict(w) for w in spec.wrappers])
    d["solver"] = _solver_dict(spec.solver)
    spec.warm_start === nothing || (d["warm_start"] = _warm_start_dict(spec.warm_start))
    return d
end

function _system_dict(s::SystemSpec)
    d = Dict{String,Any}("kind" => string(s.kind))
    s.template === nothing || (d["template"] = string(s.template))
    isempty(s.params) || (d["params"] = _wireval(s.params))
    isempty(s.global_params) || (d["global_params"] = _wireval(s.global_params))
    s.components === nothing || (d["components"] = _wireval(s.components))
    s.H_drift === nothing || (d["H_drift"] = _wireval(s.H_drift))
    s.H_drives === nothing || (d["H_drives"] = _wireval(s.H_drives))
    return d
end

function _goal_dict(g::GoalSpec)
    d = Dict{String,Any}("kind" => string(g.kind))
    g.gate === nothing || (d["gate"] = string(g.gate))
    g.matrix === nothing || (d["matrix"] = _wireval(g.matrix))
    g.target === nothing || (d["target"] = g.target)
    g.initial === nothing || (d["initial"] = g.initial)
    g.subsystem_levels === nothing ||
        (d["subsystem_levels"] = Any[Int(x) for x in g.subsystem_levels])
    g.subspace === nothing ||
        (d["subspace"] = Any[Any[Int(y) for y in row] for row in g.subspace])
    return d
end

_pulse_dict(p::PulseSpec) = Dict{String,Any}(
    "kind" => string(p.kind),
    "T" => p.T,
    "init" => string(p.init),
    "seed" => p.seed,
)

function _problem_dict(p::TemplateBlock)
    d = Dict{String,Any}(
        "template" => string(p.template),
        "N" => p.N,
        "goal_treatment" => string(p.goal_treatment),
        "free_dt" => (p.free_dt isa Free ? Any[p.free_dt.lo, p.free_dt.hi] : false),
        "Q" => p.Q,
        "R" => p.R,
        "free_phase" => p.free_phase,
    )
    p.final_fidelity === nothing || (d["final_fidelity"] = p.final_fidelity)
    p.R_u === nothing || (d["R_u"] = p.R_u)
    p.R_du === nothing || (d["R_du"] = p.R_du)
    p.R_ddu === nothing || (d["R_ddu"] = p.R_ddu)
    isfinite(p.du_bound) && (d["du_bound"] = p.du_bound)   # omit Inf default (JSON has no Inf)
    p.ddu_bound === nothing || (d["ddu_bound"] = p.ddu_bound)
    p.initial_phases === nothing ||
        (d["initial_phases"] = Any[Float64(x) for x in p.initial_phases])
    isempty(p.calibration_targets) ||
        (d["calibration_targets"] = Any[string(s) for s in p.calibration_targets])
    isempty(p.global_names) || (d["global_names"] = Any[string(s) for s in p.global_names])
    isempty(p.global_bounds) || (d["global_bounds"] = _wireval(p.global_bounds))
    isempty(p.objectives) || (
        d["objectives"] = Any[
            Dict{String,Any}("kind" => string(o.kind), "weight" => o.weight) for
            o in p.objectives
        ]
    )
    isempty(p.options) || (d["options"] = _wireval(p.options))
    return d
end

_integrator_dict(i::IntegratorSpec) =
    Dict{String,Any}("kind" => string(i.kind), "alg" => string(i.alg))

function _wrapper_dict(w::WrapperSpec)
    d = Dict{String,Any}("kind" => string(w.kind))
    isempty(w.variants) || (d["variants"] = Any[_wireval(v) for v in w.variants])
    w.weights === nothing || (d["weights"] = Any[Float64(x) for x in w.weights])
    return d
end

function _solver_dict(s::SolverSpec)
    d = Dict{String,Any}(
        "backend" => string(s.backend),
        "device" => string(s.device),
        "precision" => string(s.precision),
        "max_iter" => s.max_iter,
        "strategy" => string(s.strategy),
    )
    s.tol === nothing || (d["tol"] = s.tol)
    return d
end

function _warm_start_dict(ws::WarmStartSpec)
    d = Dict{String,Any}()
    ws.catalog_ref === nothing || (d["catalog_ref"] = ws.catalog_ref)
    ws.pulse_hash === nothing || (d["pulse_hash"] = ws.pulse_hash)
    return d
end

# ===========================================================================
# Structure projection — type-determining fields only (the structure_hash carve-out)
# ===========================================================================

"""
    structure_fields(spec::ProblemSpec) -> Dict{String,Any}

The type-determining subset used for `structure_hash`: system kind/template/levels,
trajectory/pulse kinds, problem template + `goal_treatment` + `free_dt` **arm only**
(free vs fixed, not the window) + `free_phase` + sorted objective kinds +
`options.leakage_constraint`, integrator kind/alg, wrapper kinds, and solver
backend/device/precision/strategy. Deliberately excludes N/T/Q/regularizer
weights so that a resize/reweight preserves `structure_hash` while changing
`problem_hash`.
"""
function structure_fields(spec::ProblemSpec)
    d = Dict{String,Any}("kind" => "control")
    sys = Dict{String,Any}("kind" => string(spec.system.kind))
    spec.system.template === nothing || (sys["template"] = string(spec.system.template))
    haskey(spec.system.params, :levels) &&
        (sys["levels"] = _wireval(spec.system.params[:levels]))
    d["system"] = sys
    if spec.trajectory.kind !== nothing
        d["trajectory_kind"] = string(spec.trajectory.kind)
    elseif spec.goal !== nothing
        d["trajectory_kind"] = string(spec.goal.kind)
    end
    spec.pulse === nothing || (d["pulse_kind"] = string(spec.pulse.kind))
    if spec.problem !== nothing
        p = spec.problem
        prob = Dict{String,Any}(
            "template" => string(p.template),
            "goal_treatment" => string(p.goal_treatment),
            "free_dt" => (p.free_dt isa Free ? "free" : "fixed"),
            "free_phase" => p.free_phase,
            "objective_kinds" => sort!(String[string(o.kind) for o in p.objectives]),
        )
        haskey(p.options, :leakage_constraint) &&
            (prob["leakage_constraint"] = _wireval(p.options[:leakage_constraint]))
        d["problem"] = prob
    end
    spec.integrator === nothing || (
        d["integrator"] = Dict{String,Any}(
            "kind" => string(spec.integrator.kind),
            "alg" => string(spec.integrator.alg),
        )
    )
    d["wrapper_kinds"] = sort!(String[string(w.kind) for w in spec.wrappers])
    d["solver"] = Dict{String,Any}(
        "backend" => string(spec.solver.backend),
        "device" => string(spec.solver.device),
        "precision" => string(spec.solver.precision),
        "strategy" => string(spec.solver.strategy),
    )
    return d
end

"""
    structure_hash(spec::ProblemSpec) -> String

`sha256` hex digest of `canonical_json(structure_fields(spec))` — the identity key
for a problem's *shape* (stable across N/T/Q/weight changes).
"""
structure_hash(spec::ProblemSpec) =
    bytes2hex(sha256(canonical_json(structure_fields(spec))))

"""
    problem_hash(spec::ProblemSpec) -> String

`sha256` hex digest of `canonical_json(full_dict(spec))` — the identity key for the
*full* problem instance (changes when any wire field changes).
"""
problem_hash(spec::ProblemSpec) = bytes2hex(sha256(canonical_json(full_dict(spec))))

@testitem "hashes: structure vs problem hash carve-outs" begin
    using Piccolo.Specs
    CONTROL_TOML = """
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
    Q = 100.0
    [solver]
    backend = "ipopt"
    device = "cpu"
    max_iter = 500
    """
    base = Specs.parse_spec(CONTROL_TOML; format = :toml)
    # changing solver.device changes structure_hash
    dev = Specs.parse_spec(
        replace(CONTROL_TOML, "device = \"cpu\"" => "device = \"gpu\"");
        format = :toml,
    )
    @test Specs.structure_hash(base) != Specs.structure_hash(dev)
    # changing only N/T/Q changes problem_hash but NOT structure_hash
    nq = Specs.parse_spec(replace(CONTROL_TOML, "N = 100" => "N = 120"); format = :toml)
    @test Specs.structure_hash(base) == Specs.structure_hash(nq)
    @test Specs.problem_hash(base) != Specs.problem_hash(nq)
    # determinism
    @test Specs.problem_hash(base) == Specs.problem_hash(base)
end

@testitem "canonical_json: pinned int/float-agnostic numeric rule" begin
    using Piccolo.Specs
    cj = Specs.canonical_json
    # integer-valued numbers render as bare integers regardless of Int vs Float
    @test cj(100) == "100"
    @test cj(100.0) == "100"
    @test cj(0) == "0"
    @test cj(-0.0) == "0"
    @test cj(-5) == "-5"
    @test cj(2.0) == "2"
    # non-integers: ECMAScript Number::toString (matches JS `.toString()`)
    @test cj(0.02) == "0.02"
    @test cj(0.01) == "0.01"
    @test cj(0.0001) == "0.0001"
    @test cj(1e-5) == "0.00001"
    @test cj(1e-6) == "0.000001"
    @test cj(1e-7) == "1e-7"
    @test cj(0.999) == "0.999"
    @test cj(1.5) == "1.5"
    @test cj(123.45) == "123.45"
    # containers: sorted keys, no spaces, bool/null literals
    @test cj(Dict("b" => 1, "a" => 2)) == "{\"a\":2,\"b\":1}"
    @test cj(Any[1, 2.0, true, nothing]) == "[1,2,true,null]"
    @test cj(:CZ) == "\"CZ\""
end
