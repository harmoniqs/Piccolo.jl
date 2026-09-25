# ===========================================================================
# emit_schema — JSON Schema (draft-07) reflected from the registries (Task 13)
#
# The wire-format-first anti-drift source: the schema is emitted FROM the live
# registry globals + the parser's allowed-key sets, so any registration change
# (or a parser field-set change) moves the schema. Plan 2 checks the emitted
# schema into the repo and adds the CI diff gate + ajv/JSONSchema.jl fixture runs.
#
# Structure
# ---------
# * Top-level discriminated union: `oneOf` on `kind` — a self-contained `control`
#   branch (`kind = "control"`) and a `rollout` branch (`kind = "rollout"`).
# * Enums are reflected from the registries: template/system/integrator/wrapper/
#   objective-term/solver/strategy names come straight from the registry `keys`.
#   With Piccolissimo NOT loaded, this is the OSS variant (Piccolo-only names:
#   `bilinear`, `sampling`, `ipopt`, `direct`, the six objective terms, seven
#   systems, three templates). Loading Piccolissimo augments the same registries
#   and the schema grows the exponential/spline integrators, hermite_* terms, the
#   robust wrapper, and the altissimo solver — no code change here.
# * Per-block `properties` mirror the parser's allowed-key sets (`_CONTROL_KEYS`,
#   `_PROBLEM_KEYS`, …) so the schema's `additionalProperties:false` whitelist is
#   byte-for-byte the strict field set `parse_spec` accepts.
# * Conditional compatibility branches come from `RegistryEntry.compat`:
#     - `template → pulse.kind`  (e.g. SmoothPulseProblem ⇒ pulse.kind ∈ {zero_order})
#     - `free_phase ∨ global_names → integrator.kind ∈ {exponential, spline}`
#   emitted under the branch's `allOf`.
#
# Draft-07 footgun avoidance
# --------------------------
# A single top-level `additionalProperties:false` does NOT see properties that a
# sibling `if/then`/`oneOf` subschema introduces, so combining strict-unknown-
# field rejection with conditional-compat branches would reject *valid* specs.
# We therefore keep `additionalProperties:false` **per branch / per nested object**
# and make every conditional branch *only narrow enums on already-declared
# properties* (`pulse.kind`, `integrator.kind`) — never introduce a new property.
# ===========================================================================

export emit_schema, precompile_workload, type_universe

const _DRAFT07 = "http://json-schema.org/draft-07/schema#"

# Non-bilinear integrator families that support free_phase/globals. The bilinear
# sentinel means `integrator=nothing` (the template's default path), which cannot
# carry free-phase globals; the ODE-backed exponential/spline integrators can.
# This mirrors the materialize-time rule (see materialize.jl trait validations)
# and the `SpecError` message "free_phase/globals require exponential or spline
# integrator". Reflected verbatim so the schema states the rule even in the OSS
# variant where those integrators are not yet registered.
const _FREE_PHASE_INTEGRATOR_KINDS = ["exponential", "spline"]

# ---------------------------------------------------------------------------
# Small schema-fragment builders
# ---------------------------------------------------------------------------

_registry_names(reg::AbstractDict) = sort!(String[string(k) for k in keys(reg)])

_str_enum(syms) =
    Dict{String,Any}("type" => "string", "enum" => sort!(String[string(s) for s in syms]))

# a `{ "type": "string", "enum": [...] }` whose enum is reflected from a registry
_registry_enum(reg::AbstractDict) =
    Dict{String,Any}("type" => "string", "enum" => _registry_names(reg))

_typed(t::AbstractString) = Dict{String,Any}("type" => t)
_array_of(items) = Dict{String,Any}("type" => "array", "items" => items)
_any() = Dict{String,Any}()   # unconstrained (reserved / free-form) fields

function _obj(props::AbstractDict; required = String[], additional::Bool = false)
    d = Dict{String,Any}(
        "type" => "object",
        "properties" => props,
        "additionalProperties" => additional,
    )
    isempty(required) || (d["required"] = required)
    return d
end

# ---------------------------------------------------------------------------
# Per-block object schemas (property keys mirror the parser's allowed-key sets)
# ---------------------------------------------------------------------------

function _system_schema()
    props = Dict{String,Any}(
        "kind" => _str_enum([:template, :raw]),   # composite deferred (Phase 1)
        "template" => _registry_enum(SYSTEMS),
        "params" => _typed("object"),
        "global_params" => _typed("object"),
        "components" => _any(),
        "H_drift" => _any(),
        "H_drives" => _any(),
    )
    _assert_keys(props, _SYSTEM_KEYS, "system")
    return _obj(props)
end

function _goal_schema()
    props = Dict{String,Any}(
        "kind" => _str_enum([:unitary, :ket]),
        "gate" => _typed("string"),
        "matrix" => _any(),
        "target" => _typed("string"),
        "initial" => _typed("string"),
        "subsystem_levels" => _array_of(_typed("integer")),
        "subspace" => _typed("array"),
    )
    _assert_keys(props, _GOAL_KEYS, "goal")
    return _obj(props; required = ["kind"])
end

function _pulse_schema()
    props = Dict{String,Any}(
        "kind" => _str_enum([:zero_order, :linear_spline, :cubic_spline]),
        "T" => _typed("number"),
        "init" => _str_enum([:default, :zeros, :random]),
        "seed" => _typed("integer"),
    )
    _assert_keys(props, _PULSE_KEYS, "pulse")
    return _obj(props; required = ["kind", "T"])
end

function _trajectory_schema()
    props = Dict{String,Any}("kind" => _str_enum([:unitary, :ket]))
    _assert_keys(props, _TRAJECTORY_KEYS, "trajectory")
    return _obj(props)
end

function _objective_term_schema()
    props = Dict{String,Any}(
        "kind" => _registry_enum(OBJECTIVE_TERMS),
        "weight" => _typed("number"),
    )
    _assert_keys(props, _OBJECTIVE_TERM_KEYS, "objectives[]")
    return _obj(props; required = ["kind"])
end

function _problem_schema()
    # free_dt poka-yoke: `false` or a `[lo, hi]` window — never a bare `true`.
    free_dt = Dict{String,Any}(
        "oneOf" => Any[
            Dict{String,Any}("const" => false),
            Dict{String,Any}(
                "type" => "array",
                "items" => _typed("number"),
                "minItems" => 2,
                "maxItems" => 2,
            ),
        ],
    )
    props = Dict{String,Any}(
        "template" => _registry_enum(TEMPLATES),
        "N" => _typed("integer"),
        "goal_treatment" => _str_enum([:objective, :constraint, :both]),
        "final_fidelity" => _typed("number"),
        "free_dt" => free_dt,
        "Q" => _typed("number"),
        "R" => _typed("number"),
        "R_u" => _typed("number"),
        "R_du" => _typed("number"),
        "R_ddu" => _typed("number"),
        "du_bound" => _typed("number"),
        "ddu_bound" => _typed("number"),
        "free_phase" => _typed("boolean"),
        "initial_phases" => _array_of(_typed("number")),
        "calibration_targets" => _array_of(_typed("string")),
        "global_names" => _array_of(_typed("string")),
        "global_bounds" => _typed("object"),
        "objectives" => _array_of(_objective_term_schema()),
        "options" => _typed("object"),
    )
    _assert_keys(props, _PROBLEM_KEYS, "problem")
    return _obj(props; required = ["template", "N"])
end

function _integrator_schema()
    props = Dict{String,Any}(
        "kind" => _registry_enum(INTEGRATORS),
        "alg" => _str_enum([:tsit5, :magnus_gl4, :magnus_adapt4]),
    )
    _assert_keys(props, _INTEGRATOR_KEYS, "integrator")
    return _obj(props)
end

function _wrapper_schema()
    props = Dict{String,Any}(
        "kind" => _registry_enum(WRAPPERS),
        "variants" => _array_of(_typed("object")),
        "weights" => _array_of(_typed("number")),
    )
    _assert_keys(props, _WRAPPER_KEYS, "wrappers[]")
    return _obj(props; required = ["kind"])
end

function _solver_schema()
    props = Dict{String,Any}(
        "backend" => _registry_enum(SOLVERS),
        "device" => _str_enum([:cpu, :gpu]),
        "precision" => _str_enum([:f64, :f32]),
        "max_iter" => _typed("integer"),
        "tol" => _typed("number"),
        "strategy" => _registry_enum(SOLVE_STRATEGIES),
    )
    _assert_keys(props, _SOLVER_KEYS, "solver")
    return _obj(props)
end

function _warm_start_schema()
    props = Dict{String,Any}(
        "catalog_ref" => _typed("string"),
        "pulse_hash" => _typed("string"),
    )
    _assert_keys(props, _WARM_START_KEYS, "warm_start")
    return _obj(props)
end

function _report_schema()
    props = Dict{String,Any}(
        "fidelity" => _typed("boolean"),
        "goal" => _goal_schema(),
        "populations" => _typed("boolean"),
        "leakage" => _typed("boolean"),
    )
    _assert_keys(props, _REPORT_KEYS, "report")
    return _obj(props)
end

function _referee_schema()
    props = Dict{String,Any}(
        "run" => _typed("string"),
        "solve_knots" => _typed("integer"),
        "solve_integrator" => _typed("string"),
        "fidelity_reported" => _typed("number"),
    )
    _assert_keys(props, _REFEREE_KEYS, "referee")
    return _obj(
        props;
        required = ["run", "solve_knots", "solve_integrator", "fidelity_reported"],
    )
end

# Internal consistency guard: the schema's declared property set for a block must
# equal the parser's strict allowed-key set for that block, so the schema
# `additionalProperties:false` whitelist and `parse_spec` never drift apart.
function _assert_keys(props::AbstractDict, allowed::AbstractSet, block::AbstractString)
    have = Set(String.(keys(props)))
    want = Set(String.(allowed))
    have == want || error(
        "schema property set for `$block` = $(sort!(collect(have))) " *
        "≠ parser allowed keys $(sort!(collect(want)))",
    )
    return nothing
end

# ---------------------------------------------------------------------------
# Conditional compatibility branches (reflected from RegistryEntry.compat)
# ---------------------------------------------------------------------------

# For each template carrying `compat[:pulse_kinds]`, emit
#   if  problem.template == <name>
#   then pulse.kind ∈ <pulse_kinds>
# Only narrows the enum of the already-declared `pulse.kind` (no new property).
# A one-property nested constraint `{ required: [name], properties: { name: sub } }`.
_prop_constraint(name::AbstractString, sub) =
    Dict{String,Any}("required" => [name], "properties" => Dict{String,Any}(name => sub))
# The same, but without the `required` (constrain-if-present).
_prop_only(name::AbstractString, sub) =
    Dict{String,Any}("properties" => Dict{String,Any}(name => sub))

function _template_pulse_conditionals()
    conds = Any[]
    for name in _registry_names(TEMPLATES)
        entry = lookup_template(Symbol(name))
        entry === nothing && continue
        pk = get(entry.compat, :pulse_kinds, nothing)
        pk === nothing && continue
        if_part = _prop_constraint(
            "problem",
            _prop_constraint("template", Dict{String,Any}("const" => name)),
        )
        then_part = _prop_only(
            "pulse",
            _prop_only(
                "kind",
                Dict{String,Any}("enum" => sort!(String[string(x) for x in pk])),
            ),
        )
        push!(conds, Dict{String,Any}("if" => if_part, "then" => then_part))
    end
    return conds
end

# ---------------------------------------------------------------------------
# Per-template parameter schema, REFLECTED from the params structs (Phase 1b)
# ---------------------------------------------------------------------------
#
# `RegistryEntry.params` is generated by `@problem_template` from the template's
# `@kwdef` params struct (`reflect_params`: fieldnames × fieldtypes × defaults).
# Phase 1 could only hand-write the `[problem]` block's property types, so
# per-template keyword validity — "`R_ddu` exists only on `SmoothPulseParams`" —
# was invisible to the schema and had to be discovered by a bounced tool call.
#
# Now each template emits a conditional branch:
#
#   if   problem.template == <name>
#   then problem.properties.<own param>  narrowed to the field's reflected type
#        (with its declared default), AND
#        problem must NOT carry any `[problem]`-expressible parameter that belongs
#        to a *different* template's params struct.
#
# so the LLM sees the per-template keyword matrix inside the tool schema and mostly
# never generates an invalid combination.

# Julia field type → JSON Schema type. `nothing`-able and multi-typed fields are
# left unconstrained (the parser and the params struct still enforce them); only
# unambiguous scalar/array shapes are narrowed.
#
# `RegistryEntry.params` is normally *reflected* (`reflect_params`), so `type` is a
# real `Type`. But the registry API also accepts a hand-declared params schema whose
# `type` is a type NAME (`:Float64`) — Phase 1's documented shape, still used by
# third-party/ad-hoc registrations. Both are handled; anything else is left
# unconstrained rather than crashing the emitter.
_json_type_of(::Any) = nothing

const _JSON_TYPE_BY_NAME = Dict{Symbol,String}(
    :Bool => "boolean",
    :Int => "integer",
    :Int64 => "integer",
    :Integer => "integer",
    :Float64 => "number",
    :Real => "number",
    :Number => "number",
    :String => "string",
    :Symbol => "string",
)

function _json_type_of(name::Symbol)
    t = get(_JSON_TYPE_BY_NAME, name, nothing)
    return t === nothing ? nothing : Dict{String,Any}("type" => t)
end

function _json_type_of(T::Type)
    T = T === Any ? Any : T
    T <: Bool && return Dict{String,Any}("type" => "boolean")
    T <: Integer && return Dict{String,Any}("type" => "integer")
    T <: Real && return Dict{String,Any}("type" => "number")
    T <: AbstractString && return Dict{String,Any}("type" => "string")
    T === Symbol && return Dict{String,Any}("type" => "string")
    T <: AbstractVector{Symbol} && return _array_of(_typed("string"))
    T <: AbstractVector{Int} && return _array_of(_typed("integer"))
    T <: AbstractVector{<:Real} && return _array_of(_typed("number"))
    T <: AbstractDict && return _typed("object")
    if T isa Union
        # Union{Nothing,X} → X's shape; anything wider stays unconstrained.
        nn = filter(!=(Nothing), Base.uniontypes(T))
        length(nn) == 1 && return _json_type_of(only(nn))
    end
    return nothing
end

# A JSON-representable default, or `nothing` (meaning: do not advertise one).
function _json_default(v)
    v isa Bool && return v
    v isa Integer && return v
    v isa Real && return isfinite(v) ? Float64(v) : nothing
    v isa Symbol && return string(v)
    v isa AbstractVector{Symbol} && return String[string(x) for x in v]
    v isa AbstractVector{<:Real} && return collect(Float64, v)
    return nothing
end

# The params fields a template declares that the `[problem]` block can actually
# carry (i.e. that are also parser-allowed `[problem]` keys).
function _expressible_params(entry)
    return Set{String}(String(k) for k in keys(entry.params) if String(k) in _PROBLEM_KEYS)
end

# `_json_default` must tolerate anything a hand-declared entry may carry.
_json_default(::Nothing) = nothing

function _template_param_conditionals()
    conds = Any[]
    names = _registry_names(TEMPLATES)
    entries = Dict{String,Any}()
    for n in names
        e = lookup_template(Symbol(n))
        e === nothing && continue
        entries[n] = e
    end
    # union of every template's `[problem]`-expressible params, so "foreign" is
    # computable per template
    all_expressible = Set{String}()
    for (_, e) in entries
        union!(all_expressible, _expressible_params(e))
    end
    for n in sort!(collect(keys(entries)))
        e = entries[n]
        own = _expressible_params(e)
        isempty(own) && continue
        props = Dict{String,Any}()
        for k in sort!(collect(own))
            spec = e.params[Symbol(k)]
            hasproperty(spec, :type) || continue
            t = _json_type_of(spec.type)
            t === nothing && continue
            sub = copy(t)
            d = hasproperty(spec, :default) ? _json_default(spec.default) : nothing
            d === nothing || (sub["default"] = d)
            props[k] = sub
        end
        foreign = sort!(collect(setdiff(all_expressible, own)))
        then_problem = Dict{String,Any}()
        isempty(props) || (then_problem["properties"] = props)
        if !isempty(foreign)
            then_problem["allOf"] = Any[
                Dict{String,Any}("not" => Dict{String,Any}("required" => [f])) for
                f in foreign
            ]
        end
        isempty(then_problem) && continue
        if_part = _prop_constraint(
            "problem",
            _prop_constraint("template", Dict{String,Any}("const" => n)),
        )
        push!(
            conds,
            Dict{String,Any}(
                "if" => if_part,
                "then" => _prop_only("problem", then_problem),
            ),
        )
    end
    return conds
end

# free_phase ∨ non-empty global_names ⟹ integrator.kind ∈ {exponential, spline}.
function _free_phase_conditional()
    if_part = _prop_constraint(
        "problem",
        Dict{String,Any}(
            "anyOf" => Any[
                _prop_constraint("free_phase", Dict{String,Any}("const" => true)),
                _prop_constraint(
                    "global_names",
                    Dict{String,Any}("type" => "array", "minItems" => 1),
                ),
            ],
        ),
    )
    then_part = _prop_constraint(
        "integrator",
        _prop_constraint(
            "kind",
            Dict{String,Any}("enum" => copy(_FREE_PHASE_INTEGRATOR_KINDS)),
        ),
    )
    return Dict{String,Any}("if" => if_part, "then" => then_part)
end

# ---------------------------------------------------------------------------
# The two discriminated-union branches
# ---------------------------------------------------------------------------

function _control_schema()
    props = Dict{String,Any}(
        "schema_version" => Dict{String,Any}("type" => "integer", "enum" => [1]),
        "kind" => Dict{String,Any}("const" => "control"),
        "system" => _system_schema(),
        "goal" => _goal_schema(),
        "pulse" => _pulse_schema(),
        "trajectory" => _trajectory_schema(),
        "problem" => _problem_schema(),
        "integrator" => _integrator_schema(),
        "wrappers" => _array_of(_wrapper_schema()),
        "solver" => _solver_schema(),
        "warm_start" => _warm_start_schema(),
    )
    _assert_keys(props, _CONTROL_KEYS, "control")
    d = _obj(props; required = ["system"])
    allof = Any[]
    append!(allof, _template_pulse_conditionals())
    append!(allof, _template_param_conditionals())
    push!(allof, _free_phase_conditional())
    d["allOf"] = allof
    return d
end

function _rollout_schema()
    props = Dict{String,Any}(
        "schema_version" => Dict{String,Any}("type" => "integer", "enum" => [1]),
        "kind" => Dict{String,Any}("const" => "rollout"),
        "input_pulse" => _typed("string"),
        "system" => _system_schema(),
        "rollout_kind" => _str_enum([:unitary, :ket, :density]),
        "alg" => _str_enum([:tsit5, :magnus_gl4, :magnus_adapt4]),
        "n_samples" => _typed("integer"),
        "times" => _array_of(_typed("number")),
        "initial" => _typed("string"),
        "report" => _report_schema(),
        "referee" => _referee_schema(),
    )
    _assert_keys(props, _ROLLOUT_KEYS, "rollout")
    return _obj(props; required = ["input_pulse", "rollout_kind", "system"])
end

"""
    emit_schema() -> String

Emit the Piccolo `ProblemSpec` JSON Schema (draft-07) as a string, reflected from
the live registries and the parser's allowed-key sets. The schema is a
discriminated union (`oneOf` on `kind`) over a `control` branch and a `rollout`
branch; enums for templates/systems/integrators/wrappers/objective-terms/solvers/
    strategies come from the registry `keys`, and conditional compatibility branches
    (`template → pulse.kind`; `template → per-template parameter types and the
parameters it does NOT accept`; `free_phase ∨ globals → integrator.kind`) come from
each `RegistryEntry.compat` and `RegistryEntry.params` — the latter reflected by
`@problem_template` from the template's `@kwdef` params struct, so the schema's
per-template keyword matrix and the Julia constructors cannot drift apart. `additionalProperties:false` is applied per branch / per
nested object (never once at the top level) to avoid the draft-07
`additionalProperties` + `if/then` footgun. `schema_version` is an INTEGER enum
`[1]`. With Piccolissimo unloaded this is the OSS variant; loading it augments the
same registries and the schema grows automatically.
"""
function emit_schema()
    schema = Dict{String,Any}(
        "\$schema" => _DRAFT07,
        "\$id" => "https://harmoniqs.co/schemas/problemspec.schema.json",
        "title" => "Piccolo ProblemSpec (Phase 1, wire-format-first)",
        "description" =>
            "Declarative, versioned Piccolo problem specification. " *
            "Reflected from the Piccolo.Specs registries; " *
            "`oneOf` discriminates on `kind`.",
        "oneOf" => Any[_control_schema(), _rollout_schema()],
    )
    return JSON3.write(schema)
end

# ---------------------------------------------------------------------------
# Schema introspection helpers (used by tests + Plan 2 fixture tooling)
# ---------------------------------------------------------------------------

# Recursively collect every value stored under `key` at any depth of a parsed
# JSON value (works on JSON3.Object/JSON3.Array and Dict/Vector alike).
function _collect_under_key(x, key::Symbol, acc = Any[])
    if x isa AbstractDict
        for k in keys(x)
            v = x[k]
            Symbol(k) === key && push!(acc, v)
            _collect_under_key(v, key, acc)
        end
    elseif x isa AbstractVector
        for v in x
            _collect_under_key(v, key, acc)
        end
    end
    return acc
end

# Recursively collect every subschema that is a conditional (has both `if`+`then`).
function _conditionals(x, acc = Any[])
    if x isa AbstractDict
        (haskey(x, :if) && haskey(x, :then)) && push!(acc, x)
        for k in keys(x)
            _conditionals(x[k], acc)
        end
    elseif x isa AbstractVector
        for v in x
            _conditionals(v, acc)
        end
    end
    return acc
end

_as_strings(v) = String[string(x) for x in v]

"""
    schema_template_names(schema) -> Vector{String}

Given a parsed `emit_schema()` value, return the problem-template names advertised
in the schema (the `problem.template` enum).
"""
function schema_template_names(schema)
    names = String[]
    for tv in _collect_under_key(schema, :template)
        tv isa AbstractDict || continue
        e = get(tv, :enum, nothing)
        e === nothing || append!(names, _as_strings(e))
    end
    return unique!(names)
end

# does the `if` part select `key` == `value` (via `const` or `enum` membership)?
function _if_selects(cond, key::Symbol, value::AbstractString)
    for tv in _collect_under_key(cond[:if], key)
        tv isa AbstractDict || continue
        c = get(tv, :const, nothing)
        c !== nothing && string(c) == value && return true
        e = get(tv, :enum, nothing)
        e !== nothing && value in _as_strings(e) && return true
    end
    return false
end

# does the `then` part constrain `outer.inner` to an `enum` containing `value`?
function _then_enum_contains(cond, outer::Symbol, inner::Symbol, value::AbstractString)
    for ov in _collect_under_key(cond[:then], outer)
        for iv in _collect_under_key(ov, inner)
            iv isa AbstractDict || continue
            e = get(iv, :enum, nothing)
            e !== nothing && value in _as_strings(e) && return true
        end
    end
    return false
end

"""
    schema_has_conditional(schema, template, pulse_kind) -> Bool

True iff the schema carries a conditional branch selecting
`problem.template == template` and constraining `pulse.kind` to an enum containing
`pulse_kind` (the `template → pulse.kind` compatibility branch).
"""
function schema_has_conditional(
    schema,
    template::AbstractString,
    pulse_kind::AbstractString,
)
    for cond in _conditionals(schema)
        _if_selects(cond, :template, template) || continue
        _then_enum_contains(cond, :pulse, :kind, pulse_kind) && return true
    end
    return false
end

"""
    schema_free_phase_requires_nonbilinear(schema) -> Bool

True iff the schema carries the `free_phase ∨ globals → integrator.kind ∈
{exponential, spline}` compatibility branch.
"""
function schema_free_phase_requires_nonbilinear(schema)
    for cond in _conditionals(schema)
        (
            isempty(_collect_under_key(cond[:if], :free_phase)) &&
            isempty(_collect_under_key(cond[:if], :global_names))
        ) && continue
        (
            _then_enum_contains(cond, :integrator, :kind, "exponential") ||
            _then_enum_contains(cond, :integrator, :kind, "spline")
        ) && return true
    end
    return false
end

# ---------------------------------------------------------------------------
# precompile_workload — tier-1 warm-up sweep (callable; NOT wired to @compile_workload)
# ---------------------------------------------------------------------------

# Tiny tier-1 specs: one per base template/pulse pairing on the smallest usable
# TransmonSystem instance. Kept minimal so a single Ipopt step is cheap.
function _tier1_specs()
    return String[
        # cubic_spline + SplinePulseProblem
        """
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
        """,
        # zero_order + SmoothPulseProblem
        """
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
        kind = "zero_order"
        T = 10.0
        [problem]
        template = "SmoothPulseProblem"
        N = 11
        """,
    ]
end

# ---------------------------------------------------------------------------
# The registry-enumerated TYPE universe
# ---------------------------------------------------------------------------
#
# This is the payoff of the parametric-type rewrite. Because a template is a
# constrained alias of `QuantumControlProblem{Tag, QT}` and `QT` is parametric in
# the pulse, the set of concrete problem types a spec can name is *computable from
# the registry*: `TEMPLATE_DECLARATIONS` × `pulse_kinds` × `trajectory_kinds`. So
# the precompile workload is literally "enumerate the registered concrete types" —
# and since `structure_hash` covers exactly the type-determining spec fields, a warm
# worker routed by `structure_hash` is guaranteed to have the concrete types
# compiled: same `structure_hash` ⇒ same concrete types ⇒ no JIT.
#
# Piccolissimo's own declarations land in the same dict, so its templates join the
# sweep with no change here.

const _WORKLOAD_LEVELS = 2
const _WORKLOAD_N = 6
const _WORKLOAD_T = 1.0

function _workload_pulse(kind::Symbol, n_drives::Int, N::Int, times)
    kind === :zero_order && return ZeroOrderPulse(zeros(n_drives, N), times)
    kind === :linear_spline && return LinearSplinePulse(zeros(n_drives, N), times)
    kind === :cubic_spline &&
        return CubicSplinePulse(zeros(n_drives, N), zeros(n_drives, N), times)
    return nothing
end

function _workload_trajectory(kind::Symbol, sys, pulse)
    if kind === :unitary
        return UnitaryTrajectory(sys, pulse, EmbeddedOperator(GATES[:X], sys))
    elseif kind === :ket
        ψ0 = ComplexF64[i == 1 ? 1.0 : 0.0 for i = 1:sys.levels]
        ψ1 = ComplexF64[i == 2 ? 1.0 : 0.0 for i = 1:sys.levels]
        return KetTrajectory(sys, pulse, ψ0, ψ1)
    end
    return nothing
end

"""
    type_universe() -> Vector{Tuple{Symbol,Symbol,Symbol}}

Every `(template, pulse_kind, trajectory_kind)` triple the registered templates
admit — the enumeration of the concrete `QuantumControlProblem{Tag, QT}` types a
spec can name. Reflected from `TEMPLATE_DECLARATIONS`, so it grows with
Piccolissimo's templates automatically.
"""
function type_universe()
    out = Tuple{Symbol,Symbol,Symbol}[]
    for (name, d) in Control.ProblemTemplates.TEMPLATE_DECLARATIONS
        for pk in d.pulse_kinds, tk in d.trajectory_kinds
            push!(out, (name, pk, tk))
        end
    end
    sort!(out)
    return out
end

"""
    precompile_workload(; tier1_only=true)

Compile the reachable problem-construction surface at tiny size.

Two sweeps:

1. **The registry-enumerated type universe** ([`type_universe`](@ref)): construct one
   tiny problem per `(template, pulse_kind, trajectory_kind)` triple, plus the
   `sampling` wrapper over one of them, so every concrete
   `QuantumControlProblem{Tag, QT}` (and `SamplingProblem{…}`) a spec can name is
   compiled. This is the "same `structure_hash` ⇒ same concrete types ⇒ no JIT on a
   warm worker" guarantee.
2. **The tier-1 spec path**: `parse_spec` → `materialize` → one Ipopt step, so the
   declarative entry point is warm too.

Wired to `PrecompileTools.@compile_workload` in `Piccolo.jl`. Best-effort: each
combination is guarded so one failure never aborts the sweep. Returns `nothing`.
"""
function precompile_workload(; tier1_only::Bool = true)
    # Silence construction-time `@warn`s (e.g. DirectTrajOpt's default-Δt-bound
    # notice) so precompilation produces no output.
    return Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
        _precompile_workload(tier1_only)
    end
end

function _precompile_workload(tier1_only::Bool)
    register_all!()
    opts = PiccoloOptions(; display = :silent)

    # (1) the type universe — construction only; no solve needed to specialize the
    #     constructors, objectives, integrators and display path.
    sys = TransmonSystem(; levels = _WORKLOAD_LEVELS, drive_bounds = [0.1, 0.1])
    times = collect(range(0.0, _WORKLOAD_T, length = _WORKLOAD_N))
    first_qcp = nothing
    for (tname, pk, tk) in type_universe()
        try
            decl = Control.ProblemTemplates.TEMPLATE_DECLARATIONS[tname]
            pulse = _workload_pulse(pk, sys.n_drives, _WORKLOAD_N, times)
            pulse === nothing && continue
            qtraj = _workload_trajectory(tk, sys, pulse)
            qtraj === nothing && continue
            qcp = decl.alias_type(
                qtraj,
                _WORKLOAD_N;
                Δt_bounds = (1e-3, 0.5),
                piccolo_options = opts,
            )
            first_qcp === nothing && (first_qcp = qcp)
        catch err
            @debug "precompile_workload: a type-universe combination failed" tname pk tk exception =
                (err, catch_backtrace())
        end
    end
    # the wrapper types are part of the universe too
    if first_qcp !== nothing
        try
            SamplingProblem(first_qcp, [sys, sys]; piccolo_options = opts)
        catch err
            @debug "precompile_workload: the sampling wrapper failed" exception =
                (err, catch_backtrace())
        end
    end

    # (2) the declarative path
    specs = _tier1_specs()
    # (tier2 combinations reserved; tier1_only is the default and only tier today)
    tier1_only || nothing
    for src in specs
        try
            spec = parse_spec(src; format = :toml)
            qcp = materialize(spec; piccolo_options = opts)
            # Ipopt prints its EPL banner unconditionally on first solve; keep
            # precompilation output clean.
            redirect_stdout(devnull) do
                solve!(qcp; max_iter = 1, print_level = 0, verbose = false)
            end
        catch err
            @debug "precompile_workload: a tier-1 combination failed" exception =
                (err, catch_backtrace())
        end
    end
    return nothing
end

@testitem "emit_schema: discriminated union + conditional branches" begin
    using Piccolo.Specs, JSON3
    Specs.register_all!()
    sch = JSON3.read(Specs.emit_schema())
    @test haskey(sch, :oneOf)   # kind discriminator
    names = Specs.schema_template_names(sch)
    @test "SmoothPulseProblem" in names && "SplinePulseProblem" in names
    # conditional branch present: SmoothPulseProblem → pulse zero_order
    @test Specs.schema_has_conditional(sch, "SmoothPulseProblem", "zero_order")
end

@testitem "emit_schema: integer schema_version enum + per-branch additionalProperties" begin
    using Piccolo.Specs, JSON3
    Specs.register_all!()
    sch = JSON3.read(Specs.emit_schema())
    branches = sch[:oneOf]
    @test length(branches) == 2
    # each branch is self-contained with additionalProperties:false — per-branch, to
    # avoid the draft-07 top-level additionalProperties + if/then footgun.
    for b in branches
        @test b[:type] == "object"
        @test b[:additionalProperties] == false
        sv = b[:properties][:schema_version]
        @test collect(sv[:enum]) == [1]      # INTEGER enum, not the string "1"
        @test sv[:enum][1] === 1
    end
    # OSS variant: Piccolo-only integrator kinds (bilinear); Piccolissimo not loaded.
    ctrl = branches[1][:properties][:kind][:const] == "control" ? branches[1] : branches[2]
    intk = ctrl[:properties][:integrator][:properties][:kind][:enum]
    @test "bilinear" in intk
    @test !("exponential" in intk)
    # free_phase/globals conditional requires a non-bilinear integrator.
    @test Specs.schema_free_phase_requires_nonbilinear(sch)
end

@testitem "emit_schema: schema whitelist == parser allowed-key sets (anti-drift)" begin
    using Piccolo.Specs, JSON3
    Specs.register_all!()
    sch = JSON3.read(Specs.emit_schema())
    ctrl =
        sch[:oneOf][1][:properties][:kind][:const] == "control" ? sch[:oneOf][1] :
        sch[:oneOf][2]
    # control-branch property set is exactly the parser's control field whitelist
    @test Set(String.(keys(ctrl[:properties]))) == Set([
        "schema_version",
        "kind",
        "system",
        "goal",
        "pulse",
        "trajectory",
        "problem",
        "integrator",
        "wrappers",
        "solver",
        "warm_start",
    ])
    # nested problem block mirrors the parser's problem field whitelist + closed
    prob = ctrl[:properties][:problem]
    @test prob[:additionalProperties] == false
    @test "free_dt" in String.(keys(prob[:properties]))
    @test "objectives" in String.(keys(prob[:properties]))
    # free_dt is the false-or-[lo,hi] poka-yoke (never a bare true)
    @test haskey(prob[:properties][:free_dt], :oneOf)
end

@testitem "precompile_workload is callable and returns nothing" begin
    using Piccolo.Specs
    @test Specs.precompile_workload isa Function
    @test Specs.precompile_workload(; tier1_only = true) === nothing
end

@testitem "emit_schema: per-template parameter fragment is reflected from the params structs" begin
    using Piccolo, Piccolo.Specs, JSON3
    Specs.register_all!()
    sch = JSON3.read(Specs.emit_schema())
    ctrl = sch.oneOf[1].properties.kind.const == "control" ? sch.oneOf[1] : sch.oneOf[2]

    # collect the per-template parameter branches out of the control branch's allOf
    branches = Dict{String,Any}()
    for b in ctrl.allOf
        ip = get(get(b[:if], :properties, Dict()), :problem, nothing)
        ip === nothing && continue
        tp = get(get(ip, :properties, Dict()), :template, nothing)
        tp === nothing && continue
        t = get(tp, :const, nothing)
        t === nothing && continue
        th = get(get(b[:then], :properties, Dict()), :problem, nothing)
        th === nothing && continue
        branches[String(t)] = th
    end
    @test sort!(collect(keys(branches))) ==
          ["BangBangPulseProblem", "SmoothPulseProblem", "SplinePulseProblem"]

    decls = Piccolo.Control.ProblemTemplates.TEMPLATE_DECLARATIONS
    forbidden(t) =
        haskey(branches[t], :allOf) ?
        sort!(String[String(x.not.required[1]) for x in branches[t].allOf]) : String[]

    # the "does NOT accept" list is exactly the foreign `[problem]`-expressible
    # params of the *other* templates — generated from the same source as the
    # constructors' struct-field check, so the two cannot disagree
    for (name, decl) in decls
        t = string(name)
        haskey(branches, t) || continue
        own = fieldnames(decl.params_type)
        for f in forbidden(t)
            @test !(Symbol(f) in own)
            # and it really is some other template's parameter
            @test any(d -> Symbol(f) in fieldnames(d.params_type), values(decls))
        end
    end

    # R_ddu / ddu_bound exist ONLY on SmoothPulseParams, so the schema forbids them
    # for the other two — the same fact the constructors enforce
    @test "R_ddu" in forbidden("BangBangPulseProblem")
    @test "ddu_bound" in forbidden("BangBangPulseProblem")
    @test "R_ddu" in forbidden("SplinePulseProblem")
    @test "ddu_bound" in forbidden("SplinePulseProblem")
    @test isempty(forbidden("SmoothPulseProblem"))
    @test :R_ddu in fieldnames(decls[:SmoothPulseProblem].params_type)
    @test !(:R_ddu in fieldnames(decls[:BangBangPulseProblem].params_type))

    # advertised defaults come from the params structs (bang-bang's R_u default is
    # 0.0 -- bang-bang needs no amplitude regularization -- vs smooth's R)
    smooth_props = branches["SmoothPulseProblem"].properties
    @test smooth_props.Q.default == Piccolo.Control.ProblemTemplates.SmoothPulseParams().Q
    @test smooth_props.R.default == Piccolo.Control.ProblemTemplates.SmoothPulseParams().R

    # every advertised property is a real field of that template's params struct
    for (t, th) in branches
        haskey(th, :properties) || continue
        name = Symbol(t)
        own = fieldnames(decls[name].params_type)
        for k in keys(th.properties)
            @test Symbol(k) in own
        end
    end
end

@testitem "type_universe enumerates the registered concrete-problem-type universe" begin
    using Piccolo, Piccolo.Specs
    Specs.register_all!()
    u = Specs.type_universe()

    # every triple is (registered template, one of its pulse kinds, one of its
    # trajectory kinds) — reflected from TEMPLATE_DECLARATIONS, so a new template
    # (Piccolissimo's) joins the sweep with no change here
    decls = Piccolo.Control.ProblemTemplates.TEMPLATE_DECLARATIONS
    @test !isempty(u)
    for (t, pk, tk) in u
        @test haskey(decls, t)
        @test pk in decls[t].pulse_kinds
        @test tk in decls[t].trajectory_kinds
    end
    # and it is the FULL cross product, so nothing reachable is left uncompiled
    expected =
        sum(length(d.pulse_kinds) * length(d.trajectory_kinds) for d in values(decls))
    @test length(u) == expected

    # the three base templates and both spline pulse kinds are in there
    @test (:SmoothPulseProblem, :zero_order, :unitary) in u
    @test (:BangBangPulseProblem, :zero_order, :ket) in u
    @test (:SplinePulseProblem, :cubic_spline, :unitary) in u
    @test (:SplinePulseProblem, :linear_spline, :ket) in u

    # each triple names a constructible concrete type: the alias applied to the
    # trajectory type is a type, and the spec-side structure fields determine it
    @test all(decls[t].alias_type isa Type for (t, _, _) in u)
end

@testitem "emit_schema tolerates a hand-declared params schema (Symbol type names)" begin
    using Piccolo, Piccolo.Specs, JSON3
    Specs.register_all!()
    # The registry API accepts a hand-declared parameter schema whose `type` is a
    # type NAME rather than a `Type` (Phase 1's documented shape, still reachable for
    # ad-hoc/third-party registrations). The reflected per-template fragment must
    # handle it — and must not crash on an entry it cannot interpret at all.
    Specs.register_template!(
        :HandDeclaredTmpl,
        Specs.RegistryEntry(;
            factory = Specs.ConstFactory(:hand),
            params = Dict{Symbol,Any}(
                :Q => (; type = :Float64, default = 1.0),        # type NAME
                :free_phase => (; type = :Bool, default = false),
                :du_bound => (; type = :mystery),               # uninterpretable + no default
            ),
            compat = Dict{Symbol,Any}(:pulse_kinds => [:zero_order]),
        ),
    )
    try
        sch = JSON3.read(Specs.emit_schema())      # must not throw
        @test "HandDeclaredTmpl" in Specs.schema_template_names(sch)
        ctrl = sch.oneOf[1].properties.kind.const == "control" ? sch.oneOf[1] : sch.oneOf[2]
        branch = nothing
        for b in ctrl.allOf
            ip = get(get(b[:if], :properties, Dict()), :problem, nothing)
            ip === nothing && continue
            tp = get(get(ip, :properties, Dict()), :template, nothing)
            tp === nothing && continue
            get(tp, :const, nothing) == "HandDeclaredTmpl" || continue
            th = get(get(b[:then], :properties, Dict()), :problem, nothing)
            th === nothing && continue
            branch = th
        end
        @test branch !== nothing
        # the interpretable names were narrowed …
        @test branch.properties.Q.type == "number"
        @test branch.properties.Q.default == 1.0
        @test branch.properties.free_phase.type == "boolean"
        # … and the uninterpretable one was skipped rather than crashing the emitter
        @test !haskey(branch.properties, :du_bound)
    finally
        delete!(Specs.TEMPLATES, :HandDeclaredTmpl)
    end
    @test !haskey(Specs.TEMPLATES, :HandDeclaredTmpl)
end
