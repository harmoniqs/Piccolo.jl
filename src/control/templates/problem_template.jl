# ============================================================================= #
# `@problem_template` — the single declaration point for a problem template
# ============================================================================= #
#
# Compatibility knowledge used to live in method signatures
# (`SmoothPulseProblem(qtraj::AbstractQuantumTrajectory{<:ZeroOrderPulse}, …)`)
# where neither the schema emitter nor downstream dispatch could reach it, and a
# constructed problem forgot which template built it — templates were *functions*.
#
# `@problem_template` promotes that knowledge to the type system. One declaration
# emits, in order:
#
#   1. the singleton tag type (`SmoothPulseTemplate <: AbstractProblemTemplate`),
#   2. the **constrained type alias** — the `Vector{T} = Array{T,1}` pattern, so
#      `SmoothPulseProblem` becomes a *type*, not a function, and its bound *is*
#      the compatibility matrix,
#   3. constructor methods on the alias (kwargs splat into the `@kwdef` params
#      struct — per-template keyword validity is now a struct-field check — with
#      the body delegating to the existing construction logic) plus the catch-all
#      carrying the "use SplinePulseProblem for spline pulses" hint,
#   4. Holy-trait methods (`pulse_family`, `supported_trajectories`,
#      `supports_free_phase`, `template_params_type`, `state_dependent_terms`),
#   5. a `TemplateDeclaration` record from which `Specs.register_all!` builds the
#      `TEMPLATES` `RegistryEntry` (params schema + compat) — so the registry, the
#      emitted JSON Schema, and the Julia types cannot drift apart.
#
# `Live`/`Frozen` term typing rides on (4): a template declares which physics
# terms are state-dependent and the macro emits `Live`-only handling for them.

export @problem_template
export Live, Frozen, build_term

# The module that owns the trait generic functions and the tagged problem type.
# Interpolated (as a Module object) into generated code so a template declared in
# *another* package (Piccolissimo's `TweezerTransportProblem`) extends the same
# functions rather than shadowing them.
const _TRAIT_HOME = QuantumControlProblems

# The module that owns the declaration registry + generated-constructor helpers.
const _TEMPLATE_HOME = @__MODULE__

"""
    TemplateDeclaration

The `@problem_template` declaration record: everything the type system, the
validator, and the schema emitter need about one problem template, in one place.

`Specs.register_all!` turns each record into the `TEMPLATES` `RegistryEntry`
(`factory` = the constrained alias, `params` = reflected from `params_type`,
`compat` = the pulse/trajectory matrix), which is what makes the registry and the
emitted JSON Schema *derived from* the Julia declaration instead of a parallel
hand-maintained copy.

# Fields
- `name::Symbol`: registry / public name (e.g. `:SmoothPulseProblem`)
- `tag::Symbol`, `tag_type::Type`: the template tag singleton
- `alias_type::Type`: the constrained alias (a `UnionAll` over the trajectory)
- `pulse::Type`: the accepted pulse family
- `trajectories::Tuple`: accepted quantum-trajectory types
- `pulse_kinds::Vector{Symbol}`, `trajectory_kinds::Vector{Symbol}`: the wire-format
  enum values corresponding to the above
- `ket_free_phase::Bool`: whether `free_phase` works with a `KetTrajectory`
- `params_type::Type`: the `@kwdef` params struct
- `passthrough::Vector{Symbol}`: keywords that are *runtime objects*
  (`integrator`, `constraints`, `piccolo_options`, …) — accepted, never retained
  in `params`, never spec-expressible
- `state_dependent::Vector{Symbol}`: state-dependent physics terms (see `Live`)
- `requires_N::Bool`: whether the knot-count positional argument is mandatory
- `hint::String`: the catch-all error message
"""
struct TemplateDeclaration
    name::Symbol
    tag::Symbol
    tag_type::Type
    alias_type::Type
    pulse::Type
    trajectories::Tuple
    pulse_kinds::Vector{Symbol}
    trajectory_kinds::Vector{Symbol}
    ket_free_phase::Bool
    params_type::Type
    passthrough::Vector{Symbol}
    state_dependent::Vector{Symbol}
    requires_N::Bool
    hint::String
end

"""
    TEMPLATE_DECLARATIONS::Dict{Symbol,TemplateDeclaration}

Every `@problem_template` declaration, keyed by public template name. Piccolo
fills it for the three base templates at load; Piccolissimo adds its own when it
loads. `Specs.register_all!` reads it — this dict is the anti-drift seam.
"""
const TEMPLATE_DECLARATIONS = Dict{Symbol,TemplateDeclaration}()

# ----------------------------------------------------------------------------- #
# Generated-constructor helpers (called from macro expansions)
# ----------------------------------------------------------------------------- #

"""
    _template_params(PT, passthrough, tname, kwargs) -> PT

Split a template call's keyword arguments into the typed params struct `PT` and
the runtime `passthrough` names, and *reject anything else*. A keyword that is
neither a field of `PT` nor declared passthrough is an `ArgumentError` naming the
valid keywords — this is where "`R_ddu` exists only on `SmoothPulseParams`" stops
being a convention and becomes a check.
"""
function _template_params(PT::Type, passthrough::Tuple, tname::Symbol, kwargs)
    fields = fieldnames(PT)
    pkw = Dict{Symbol,Any}()
    for (k, v) in kwargs
        if k in fields
            pkw[k] = v
        elseif k in passthrough
            continue    # runtime object: accepted, not spec-expressible, not retained
        else
            throw(
                ArgumentError(
                    "$(tname): unknown keyword argument `$(k)`.\n" *
                    "  template parameters: $(join(sort!(String[string(f) for f in fields]), ", "))\n" *
                    "  runtime passthrough: $(join(sort!(String[string(f) for f in passthrough]), ", "))",
                ),
            )
        end
    end
    return PT(; pkw...)
end

"""
    _retag(qcp, tag, params) -> QuantumControlProblem{typeof(tag), QT}

Stamp a problem built by the (untagged) construction logic with its template tag
and retained params. The retained spec, if any, is carried across.
"""
_retag(
    qcp::QuantumControlProblem,
    tag::AbstractProblemTemplate,
    params::AbstractTemplateParams,
) = QuantumControlProblem{typeof(tag)}(
    qcp.qtraj,
    qcp.prob;
    params = params,
    spec = retained_spec(qcp),
)

# ----------------------------------------------------------------------------- #
# `Live` / `Frozen` — state-dependent-term poka-yoke (spec §2b)
# ----------------------------------------------------------------------------- #

@doc raw"""
    Live(term)

A **state-dependent** physics term evaluated against the *current* iterate: the
closure `term` is called during `build`, so its adjoint chain
$\partial\,\text{term}/\partial x \cdot \partial x/\partial u$ is part of the
solve.

`Live` is the only wrapper any build path accepts (see [`build_term`](@ref)).

See also [`Frozen`](@ref).
"""
struct Live{F}
    term::F
end

(l::Live)(args...; kwargs...) = l.term(args...; kwargs...)

@doc raw"""
    Frozen(values)

A state-dependent term **precomputed along a reference trajectory** and held
fixed. Legitimate for building an initial guess; *never* legitimate inside a
solve — a state-dependent term evaluated on a reference rather than the current
iterate is the recurring "frozen-drift" bug (the solve optimizes against a term
that no longer matches the state it is producing).

`Frozen` therefore stays constructible but has **no spelling** in the build path:
no method of [`build_term`](@ref) accepts it, so a `Frozen` term reaching `build`
is a `MethodError` by construction rather than a runtime check that someone can
forget to write. Freezing becomes an explicit, greppable act with no path into the
solve.

See also [`Live`](@ref).
"""
struct Frozen{T}
    values::T
end

@doc raw"""
    build_term(tag::AbstractProblemTemplate, ::Val{name}, t::Live) -> term

Unwrap a state-dependent term on its way into `build`.

**There is deliberately no `Frozen` method** — not for any tag, not for any term
name. `build_term(tag, Val(:W), Frozen(w))` is a `MethodError`, which is the
poka-yoke: the frozen-drift bug class has no spelling in the type system rather
than being caught by a validation pass.

`@problem_template` emits one tag-specialized method per name listed in the
declaration's `state_dependent` tuple; the generic method below covers templates
that carry state-dependent terms without specializing per-name behavior.
"""
build_term(::AbstractProblemTemplate, ::Val, t::Live) = t.term

"""
    _check_state_dependent_fields(PT, state_dependent, tname)

Declaration-time enforcement of **`Live`-only field types**: if a template's params
struct carries a field named in its `state_dependent` tuple, that field's declared
type must be [`Live`](@ref). A `state_dependent` term with a `::Any` (or otherwise
`Frozen`-admitting) field would give the frozen-drift bug a place to live, so
`@problem_template` refuses the declaration rather than trusting a convention.

A state-dependent term need not appear on the params struct at all (it may be
constructed inside the template); the check only fires for fields that do exist.
"""
function _check_state_dependent_fields(PT::Type, state_dependent, tname::Symbol)
    for n in state_dependent
        n in fieldnames(PT) || continue
        T = fieldtype(PT, n)
        T <: Live || error(
            "@problem_template $(tname): state-dependent term `$(n)` has field type " *
            "$(T) on $(PT), but state-dependent term fields must be typed `Live`. " *
            "A wider field type would let a `Frozen` (reference-trajectory) value " *
            "into the solve — the frozen-drift bug class. Declare `$(n)::Live`.",
        )
    end
    return nothing
end

# ----------------------------------------------------------------------------- #
# The macro
# ----------------------------------------------------------------------------- #

# `key = value` lines inside the declaration block → Dict{Symbol,Expr}.
function _parse_template_block(block)
    Meta.isexpr(block, :block) ||
        throw(ArgumentError("@problem_template expects a `begin … end` block"))
    out = Dict{Symbol,Any}()
    for line in block.args
        line isa LineNumberNode && continue
        if Meta.isexpr(line, :(=), 2) && line.args[1] isa Symbol
            out[line.args[1]] = line.args[2]
        else
            throw(
                ArgumentError(
                    "@problem_template: expected `key = value` declarations, got `$(line)`",
                ),
            )
        end
    end
    return out
end

const _TEMPLATE_KEYS = Set{Symbol}([
    :julia_name,
    :pulse,
    :trajectories,
    :pulse_kinds,
    :trajectory_kinds,
    :ket_free_phase,
    :params,
    :passthrough,
    :state_dependent,
    :builder,
    :requires_N,
    :hint,
])

@doc raw"""
    @problem_template Tag begin
        julia_name       = SmoothPulseProblem
        pulse            = ZeroOrderPulse
        trajectories     = (UnitaryTrajectory, KetTrajectory, MultiKetTrajectory)
        pulse_kinds      = (:zero_order,)
        trajectory_kinds = (:unitary, :ket)
        ket_free_phase   = false
        params           = SmoothPulseParams
        passthrough      = (:integrator, :constraints, :piccolo_options)
        builder          = _smooth_pulse_problem
        requires_N       = true
        hint             = "…"
    end

Declare a problem template once and get its type, dispatch, validation, registry
entry, and schema fragment from that single declaration.

# Generated

1. **The tag type** `Tag <: AbstractProblemTemplate` (exported).
2. **The constrained alias**
   `const julia_name{QT<:AbstractQuantumTrajectory{<:pulse}} = QuantumControlProblem{Tag, QT}`
   (exported). The bound is the compatibility matrix: with
   `pulse = ZeroOrderPulse`, `SmoothPulseProblem{UnitaryTrajectory{CubicSplinePulse}}`
   is not a type at all (`TypeError`). Julia cannot couple struct parameters, so
   the raw spelling `QuantumControlProblem{SmoothPulseTemplate, UnitaryTrajectory{CubicSplinePulse}}`
   remains a legal *type expression* — but no constructor produces it and it fails
   `isa SmoothPulseProblem`, so enforcement is at the alias/constructor surface and
   propagation is through dispatch.
3. **Constructor methods on the alias**: `julia_name(qtraj, N; kwargs...)`
   (`requires_N = false` also emits the `julia_name(qtraj; kwargs...)` arity).
   Keywords splat into `params` (unknown keyword ⇒ `ArgumentError` listing the
   valid ones); the body delegates to `builder` and stamps the result with the tag
   and retained params. Plus a **catch-all** for the wrong pulse family carrying
   `hint`.
4. **Holy traits**: [`pulse_family`](@ref), [`supported_trajectories`](@ref),
   [`template_params_type`](@ref), [`state_dependent_terms`](@ref), and — when
   `ket_free_phase = false` — `supports_free_phase(::Tag, ::Type{<:KetTrajectory}) = false`.
   The same matrix as (2), queryable by the validator and the schema emitter.
5. **The declaration record** in [`TEMPLATE_DECLARATIONS`](@ref), from which
   `Specs.register_all!` builds the `TEMPLATES` `RegistryEntry` with the params
   struct as the parameter schema.

# Keys

| key | meaning |
|---|---|
| `julia_name` | public template name; becomes the alias (required) |
| `pulse` | accepted pulse family (required) |
| `trajectories` | accepted quantum-trajectory types (required) |
| `pulse_kinds` | wire-format `pulse.kind` values (required) |
| `trajectory_kinds` | wire-format trajectory kinds (required) |
| `params` | the `@kwdef` `AbstractTemplateParams` struct (required) |
| `builder` | the construction function the constructor delegates to (required) |
| `hint` | catch-all error text (required) |
| `ket_free_phase` | default `true` |
| `passthrough` | runtime-object keywords, default `()` |
| `state_dependent` | state-dependent term names (see [`Live`](@ref)), default `()` |
| `requires_N` | default `true` |
"""
macro problem_template(tag, block)
    tag isa Symbol ||
        throw(ArgumentError("@problem_template: tag must be a symbol, got `$(tag)`"))
    d = _parse_template_block(block)
    for k in keys(d)
        k in _TEMPLATE_KEYS || throw(
            ArgumentError(
                "@problem_template: unknown declaration key `$(k)`; " *
                "valid keys: $(join(sort!(String[string(x) for x in _TEMPLATE_KEYS]), ", "))",
            ),
        )
    end
    for k in (
        :julia_name,
        :pulse,
        :trajectories,
        :pulse_kinds,
        :trajectory_kinds,
        :params,
        :builder,
        :hint,
    )
        haskey(d, k) ||
            throw(ArgumentError("@problem_template: missing required key `$(k)`"))
    end

    name = d[:julia_name]
    name isa Symbol ||
        throw(ArgumentError("@problem_template: `julia_name` must be a symbol"))
    pulse = d[:pulse]
    trajectories = d[:trajectories]
    pulse_kinds = d[:pulse_kinds]
    trajectory_kinds = d[:trajectory_kinds]
    params = d[:params]
    builder = d[:builder]
    hint = d[:hint]
    ket_free_phase = get(d, :ket_free_phase, true)
    passthrough = get(d, :passthrough, :(()))
    state_dependent = get(d, :state_dependent, :(()))
    requires_N = get(d, :requires_N, true)

    TH = _TRAIT_HOME
    MH = _TEMPLATE_HOME
    AQT = AbstractQuantumTrajectory     # interpolated as a Type, not a name to resolve
    KT = KetTrajectory
    etag = esc(tag)
    ename = esc(name)
    epulse = esc(pulse)
    eparams = esc(params)
    ebuilder = esc(builder)
    traj_bound = :($AQT{<:$epulse})

    # `state_dependent = (:W, :σ)` is a literal tuple of symbols; read it at
    # expansion time so one `build_term` method can be emitted per name.
    sd_names = Symbol[]
    if Meta.isexpr(state_dependent, :tuple)
        for a in state_dependent.args
            a isa QuoteNode && a.value isa Symbol ? push!(sd_names, a.value) :
            throw(
                ArgumentError(
                    "@problem_template: `state_dependent` must be a literal tuple of symbols",
                ),
            )
        end
    end

    ctors = Expr[]
    push!(
        ctors,
        quote
            function (::Type{P})(qtraj::$traj_bound, N; kwargs...) where {P<:$ename}
                p = $MH._template_params(
                    $eparams,
                    $(esc(passthrough)),
                    $(QuoteNode(name)),
                    kwargs,
                )
                return $MH._retag($ebuilder(qtraj, N; kwargs...), $etag(), p)
            end
            # Catch-all: wrong pulse family for this template.
            function (::Type{P})(qtraj::$AQT, N; kwargs...) where {P<:$ename}
                error(
                    $MH._wrong_pulse_message(
                        $(QuoteNode(name)),
                        $epulse,
                        qtraj,
                        $(esc(hint)),
                    ),
                )
            end
        end,
    )
    if requires_N === false
        push!(
            ctors,
            quote
                function (::Type{P})(qtraj::$traj_bound; kwargs...) where {P<:$ename}
                    p = $MH._template_params(
                        $eparams,
                        $(esc(passthrough)),
                        $(QuoteNode(name)),
                        kwargs,
                    )
                    return $MH._retag($ebuilder(qtraj; kwargs...), $etag(), p)
                end
                function (::Type{P})(qtraj::$AQT; kwargs...) where {P<:$ename}
                    error(
                        $MH._wrong_pulse_message(
                            $(QuoteNode(name)),
                            $epulse,
                            qtraj,
                            $(esc(hint)),
                        ),
                    )
                end
            end,
        )
    end

    free_phase_trait = if ket_free_phase === false
        quote
            $TH.supports_free_phase(::$etag, ::Type{<:$KT}) = false
        end
    else
        quote end
    end

    # One `build_term` method per declared state-dependent term. `Live` only:
    # no method of `build_term` accepts `Frozen`, for any tag or name.
    sd_methods = Expr[
        :($MH.build_term(::$etag, ::Val{$(QuoteNode(n))}, t::$MH.Live) = t.term) for
        n in sd_names
    ]

    return quote
        # `Live`-only field types for the declared state-dependent terms — checked
        # first, so a bad declaration defines nothing at all.
        $MH._check_state_dependent_fields(
            $eparams,
            $(esc(state_dependent)),
            $(QuoteNode(tag)),
        )

        struct $etag <: $TH.AbstractProblemTemplate end

        const $ename{QT<:$traj_bound} = $TH.QuantumControlProblem{$etag,QT}

        $TH.pulse_family(::$etag) = $epulse
        $TH.supported_trajectories(::$etag) = $(esc(trajectories))
        $TH.template_params_type(::$etag) = $eparams
        $TH.state_dependent_terms(::$etag) = $(esc(state_dependent))
        $free_phase_trait

        $(ctors...)
        $(sd_methods...)

        $MH.TEMPLATE_DECLARATIONS[$(QuoteNode(name))] = $MH.TemplateDeclaration(
            $(QuoteNode(name)),
            $(QuoteNode(tag)),
            $etag,
            $ename,
            $epulse,
            $(esc(trajectories)),
            Symbol[$(esc(pulse_kinds))...],
            Symbol[$(esc(trajectory_kinds))...],
            $(esc(ket_free_phase)),
            $eparams,
            Symbol[$(esc(passthrough))...],
            Symbol[$(esc(state_dependent))...],
            $(esc(requires_N)),
            $(esc(hint)),
        )

        export $etag, $ename
        nothing
    end
end

"""
    _wrong_pulse_message(name, pulse, qtraj, hint) -> String

The catch-all error text: which template was called, which pulse family it
accepts, which one it got, and the `hint` pointing at the right template. This is
the generated form of the hand-written "Fallback Error Method" each template used
to carry.
"""
function _wrong_pulse_message(name::Symbol, pulse::Type, qtraj, hint::AbstractString)
    return """

  $(name) is only for $(nameof(pulse)) pulses.

  You provided a trajectory with pulse type: $(_pulse_type_of(qtraj))

  $(hint)
  """
end

_pulse_type_of(::AbstractQuantumTrajectory{P}) where {P} = P

# ============================================================================= #
# Tests
# ============================================================================= #

@testitem "@problem_template generates tag, alias, constructors, traits, registry record" setup=[PiccoloTemplateHelpers] begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra
    using Piccolo.Control.ProblemTemplates: TEMPLATE_DECLARATIONS, TemplateDeclaration

    Base.@kwdef struct ProbeParams <: AbstractTemplateParams
        Q::Float64 = 100.0
        R::Float64 = 1e-2
        widget::Symbol = :off
    end

    function _probe_build(
        qtraj,
        N;
        Q::Float64 = 100.0,
        R::Float64 = 1e-2,
        widget::Symbol = :off,
        piccolo_options::PiccoloOptions = PiccoloOptions(),
    )
        traj = NamedTrajectory(qtraj, N)
        obj = QuadraticRegularizer(:u, traj, R)
        prob = DirectTrajOptProblem(traj, obj, BilinearIntegrator(qtraj, N))
        return QuantumControlProblem(qtraj, prob)
    end

    @problem_template ProbeTemplate begin
        julia_name = ProbeProblem
        pulse = ZeroOrderPulse
        trajectories = (UnitaryTrajectory, KetTrajectory)
        pulse_kinds = (:zero_order,)
        trajectory_kinds = (:unitary, :ket)
        ket_free_phase = false
        params = ProbeParams
        passthrough = (:piccolo_options,)
        builder = _probe_build
        state_dependent = (:W,)
        hint = "For spline-based pulses use SplinePulseProblem instead."
    end

    # (5) declaration record — asserted first, then removed from the process-global
    # dict so it cannot leak into another test item's registry/schema view.
    decl = TEMPLATE_DECLARATIONS[:ProbeProblem]
    @test decl isa TemplateDeclaration
    @test decl.tag === :ProbeTemplate
    @test decl.tag_type === ProbeTemplate
    @test decl.alias_type === ProbeProblem
    @test decl.pulse === ZeroOrderPulse
    @test decl.pulse_kinds == [:zero_order]
    @test decl.trajectory_kinds == [:unitary, :ket]
    @test decl.ket_free_phase == false
    @test decl.params_type === ProbeParams
    @test decl.passthrough == [:piccolo_options]
    @test decl.state_dependent == [:W]
    @test decl.requires_N == true
    delete!(TEMPLATE_DECLARATIONS, :ProbeProblem)

    # (1) the tag type
    @test ProbeTemplate <: AbstractProblemTemplate
    @test ProbeTemplate() isa AbstractProblemTemplate

    # (2) the constrained alias — a TYPE, not a function
    @test ProbeProblem isa Type
    @test !(ProbeProblem isa Function)

    # (4) Holy traits
    @test pulse_family(ProbeTemplate()) === ZeroOrderPulse
    @test supported_trajectories(ProbeTemplate()) === (UnitaryTrajectory, KetTrajectory)
    @test template_params_type(ProbeTemplate()) === ProbeParams
    @test state_dependent_terms(ProbeTemplate()) === (:W,)
    @test supports_free_phase(ProbeTemplate(), KetTrajectory) == false
    @test supports_free_phase(ProbeTemplate(), UnitaryTrajectory) == true

    # (3) constructor methods on the alias
    N = 6
    T = 2.0
    sys = OpenQuantumSystem(zeros(ComplexF64, 2, 2), [ComplexF64[0 1; 1 0]], [1.0])
    times = collect(range(0.0, T, length = N))
    zo_qtraj = UnitaryTrajectory(sys, ZeroOrderPulse(zeros(1, N), times), GATES[:X])

    qcp = ProbeProblem(
        zo_qtraj,
        N;
        Q = 42.0,
        widget = :on,
        piccolo_options = PiccoloOptions(),
    )
    @test qcp isa ProbeProblem
    @test qcp isa QuantumControlProblem
    @test qcp isa AbstractQuantumControlProblem
    @test template_tag(qcp) === ProbeTemplate()
    @test template_params(qcp) isa ProbeParams
    @test template_params(qcp).Q == 42.0
    @test template_params(qcp).widget === :on
    @test template_params(qcp).R == 1e-2          # declared default retained
    @test retained_spec(qcp) === nothing

    # per-template keyword validity is a struct-field check
    @test_throws ArgumentError ProbeProblem(zo_qtraj, N; R_ddu = 1e-3)

    # the catch-all carries the hint
    sp_qtraj = UnitaryTrajectory(sys, LinearSplinePulse(zeros(1, N), times), GATES[:X])
    err = try
        ProbeProblem(sp_qtraj, N)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("only for ZeroOrderPulse", err.msg)
    @test occursin("SplinePulseProblem", err.msg)

    # the alias bound IS the compatibility matrix: the invalid pairing is not a type
    @test_throws TypeError ProbeProblem{typeof(sp_qtraj)}
    @test ProbeProblem{typeof(zo_qtraj)} <: ProbeProblem
end

@testitem "Live/Frozen: a Frozen term has no spelling in the build path" setup=[PiccoloTemplateHelpers] begin
    # The frozen-drift bug class (a state-dependent term evaluated along a
    # reference instead of the current iterate) is made unrepresentable rather
    # than validated against: NO method of `build_term` accepts `Frozen`.
    struct ProbeLiveTag <: AbstractProblemTemplate end

    f = x -> 2x
    @test build_term(ProbeLiveTag(), Val(:W), Live(f)) === f
    @test Live(f)(3) == 6

    frozen = Frozen([1.0, 2.0])
    @test frozen.values == [1.0, 2.0]              # constructible: reference guesses are legitimate
    @test_throws MethodError build_term(ProbeLiveTag(), Val(:W), frozen)
    @test_throws MethodError build_term(ProbeLiveTag(), Val(:anything), frozen)

    # no `build_term` method anywhere admits a Frozen argument
    @test !any(
        m -> any(t -> t <: Frozen, Base.unwrap_unionall(m.sig).parameters),
        methods(build_term),
    )

    # templates default to declaring no state-dependent terms
    @test state_dependent_terms(ProbeLiveTag()) === ()
end

@testitem "@problem_template state_dependent: Live-only field types + per-tag build_term" setup=[PiccoloTemplateHelpers] begin
    using DirectTrajOpt
    using NamedTrajectories

    # A template whose state-dependent term rides on the params struct must type it
    # `Live` — anything wider would give a Frozen (reference-trajectory) value a
    # place to live, which is the frozen-drift bug class.
    Base.@kwdef struct LiveOKParams <: AbstractTemplateParams
        Q::Float64 = 100.0
        W::Live = Live(identity)
    end

    Base.@kwdef struct LiveBadParams <: AbstractTemplateParams
        Q::Float64 = 100.0
        W::Any = nothing
    end

    _sd_build(qtraj, N; kwargs...) = error("not called in this test")

    @problem_template SDGoodTemplate begin
        julia_name = SDGoodProblem
        pulse = ZeroOrderPulse
        trajectories = (UnitaryTrajectory,)
        pulse_kinds = (:zero_order,)
        trajectory_kinds = (:unitary,)
        params = LiveOKParams
        builder = _sd_build
        state_dependent = (:W,)
        hint = "n/a"
    end
    delete!(Piccolo.Control.ProblemTemplates.TEMPLATE_DECLARATIONS, :SDGoodProblem)

    # the trait records the declaration
    @test state_dependent_terms(SDGoodTemplate()) === (:W,)

    # the macro emitted a tag-specialized `build_term` for the declared name…
    f = x -> 3x
    @test build_term(SDGoodTemplate(), Val(:W), Live(f)) === f
    @test any(m -> m.sig <: Tuple{Any,SDGoodTemplate,Val{:W},Live}, methods(build_term))
    # …and still no method anywhere admits a Frozen
    @test_throws MethodError build_term(SDGoodTemplate(), Val(:W), Frozen([1.0, 2.0]))

    # a widened field type for a declared state-dependent term is refused at
    # DECLARATION time, not at solve time
    err = try
        @eval @problem_template SDBadTemplate begin
            julia_name = SDBadProblem
            pulse = ZeroOrderPulse
            trajectories = (UnitaryTrajectory,)
            pulse_kinds = (:zero_order,)
            trajectory_kinds = (:unitary,)
            params = LiveBadParams
            builder = _sd_build
            state_dependent = (:W,)
            hint = "n/a"
        end
        nothing
    catch e
        e
    end
    @test err !== nothing
    msg = sprint(showerror, err)
    @test occursin("must be typed `Live`", msg)
    @test !haskey(Piccolo.Control.ProblemTemplates.TEMPLATE_DECLARATIONS, :SDBadProblem)

    # the base templates declare no state-dependent terms (the mechanism's first
    # real customer is transport, in Piccolissimo)
    @test state_dependent_terms(SmoothPulseTemplate()) === ()
    @test state_dependent_terms(SplinePulseTemplate()) === ()
    @test state_dependent_terms(BangBangPulseTemplate()) === ()
end
