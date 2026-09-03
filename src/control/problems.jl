module QuantumControlProblems

using DirectTrajOpt
using LinearAlgebra
using NamedTrajectories
using ...Quantum
using TestItems

import ...Quantum: get_system, get_goal, state_name, drive_name, extract_pulse, fidelity
import DirectTrajOpt.Solvers: solve!

export QuantumControlProblem
export AbstractQuantumControlProblem, AbstractProblemTemplate, AbstractTemplateParams
export AbstractProblemSpec
export NoTemplate, NoParams
export AbstractProblemWrapper, inner, base_problem, wrapper_kind, rewrap, with_problem
export quantum_trajectory, direct_problem
export template_tag, template_params, retained_spec, retain_spec!
export pulse_family, supported_trajectories, supports_free_phase, template_params_type
export state_dependent_terms
export get_trajectory, get_system, get_goal, state_name, drive_name
export solve!, sync_trajectory!, fidelity
export rollout_divergence, ROLLOUT_DIVERGENCE_RTOL, stored_phases
# Note: solve! is NOT exported to avoid ambiguity with SciMLBase.solve!
# Users should use: using DirectTrajOpt (to get solve!)

# ============================================================================= #
# Abstract interface
# ============================================================================= #

@doc raw"""
    AbstractQuantumControlProblem

Supertype of every Piccolo control problem: the template-tagged
[`QuantumControlProblem`](@ref) and the parametric *wrapper* problems
(`MinimumTimeProblem`, `SamplingProblem`) that nest one problem inside another.

Downstream code that only needs "some Piccolo problem" should annotate
`::AbstractQuantumControlProblem` rather than `::QuantumControlProblem`, so
wrapper types are accepted too.

The interface every subtype provides:

| function | meaning |
|---|---|
| `get_trajectory(p)` | the `NamedTrajectory` being optimized |
| `get_system(p)` / `get_goal(p)` | the quantum system / goal |
| `state_name(p)` / `drive_name(p)` | trajectory variable names |
| `template_tag(p)` | the [`AbstractProblemTemplate`](@ref) singleton that built it |
| `template_params(p)` | the retained [`AbstractTemplateParams`](@ref) |
| `retained_spec(p)` | the originating `ProblemSpec`, or `nothing` |
| `inner(p)` | the wrapped problem (wrappers only) |
"""
abstract type AbstractQuantumControlProblem end

@doc raw"""
    AbstractProblemTemplate

Supertype of the singleton *template tags* (`SmoothPulseTemplate`,
`SplinePulseTemplate`, `BangBangPulseTemplate`, …) that parametrize
[`QuantumControlProblem`](@ref).

The tag is what makes a constructed problem remember which template built it:
`SmoothPulseProblem` is a constrained type alias of
`QuantumControlProblem{SmoothPulseTemplate, QT}`, so `qcp isa SmoothPulseProblem`
holds after construction and `solve!`/plot/analysis/wrapper methods can dispatch
on problem *flavor*.

Tags are declared by `@problem_template`, which also emits the Holy-trait methods
[`pulse_family`](@ref), [`supported_trajectories`](@ref) and
[`supports_free_phase`](@ref) for them.
"""
abstract type AbstractProblemTemplate end

@doc raw"""
    AbstractTemplateParams

Supertype of the per-template `@kwdef` parameter structs (`SmoothPulseParams`,
`SplinePulseParams`, …) generated/declared alongside a `@problem_template`.

These are plain data: one field per spec-expressible template keyword, with the
template's own default. They are the *per-template keyword truth* — `R_ddu`
exists as a field only on `SmoothPulseParams`, so passing it to
`BangBangPulseProblem` is a construction error rather than a silently ignored
keyword — and they are what `Specs.emit_schema` reflects the per-template
parameter schema from.

A constructed problem retains its params (`template_params(qcp)`).
"""
abstract type AbstractTemplateParams end

"""
    stored_phases(traj::NamedTrajectory) -> Union{Nothing,Vector{Float64}}

Free-phase globals `φ_*` held by a trajectory, concatenated in sorted name order, or
`nothing` when the trajectory carries none.

Returns an **empty vector** — distinct from `nothing` — when `φ_*` components are declared but
hold no data. That state is corrupt input rather than "no phases", and callers are expected
to refuse it rather than silently fall back to the fixed-phase answer.
"""
function stored_phases(traj)
    names = sort!([n for n in keys(traj.global_components) if startswith(string(n), "φ_")])
    isempty(names) && return nothing
    return reduce(
        vcat,
        (Vector{Float64}(traj.global_data[traj.global_components[n]]) for n in names);
        init = Float64[],
    )
end

@doc raw"""
    AbstractProblemSpec

Supertype of the declarative wire-format spec structs (`Specs.ProblemSpec`).

Declared here, in `Piccolo.Control`, rather than in `Piccolo.Specs` so that
[`QuantumControlProblem`](@ref) can carry a *typed* `spec` field even though
`Specs` is loaded after `Control`.
"""
abstract type AbstractProblemSpec end

"""
    NoTemplate() <: AbstractProblemTemplate

Tag for a hand-built problem: `QuantumControlProblem(qtraj, prob)` produced it
directly, not a registered template. Carries no compatibility claims.
"""
struct NoTemplate <: AbstractProblemTemplate end

"""
    NoParams() <: AbstractTemplateParams

Params placeholder for hand-built problems (see [`NoTemplate`](@ref)).
"""
struct NoParams <: AbstractTemplateParams end

# The three base-template tags (`SmoothPulseTemplate`, `SplinePulseTemplate`,
# `BangBangPulseTemplate`) are emitted by `@problem_template` in
# `Control.ProblemTemplates`, not declared here — the declaration block is their
# single source of truth.

# ----------------------------------------------------------------------------- #
# Holy traits — `@problem_template` adds one method per declared template.
# ----------------------------------------------------------------------------- #

"""
    pulse_family(tag::AbstractProblemTemplate) -> Type{<:AbstractPulse}

The pulse family a template accepts (e.g. `ZeroOrderPulse` for
`SmoothPulseTemplate`). This is the same fact as the template alias' trajectory
bound, queryable by the validator and the schema emitter.
"""
function pulse_family end

"""
    supported_trajectories(tag::AbstractProblemTemplate) -> Tuple{Vararg{Type}}

The quantum-trajectory types a template can build from.
"""
function supported_trajectories end

"""
    supports_free_phase(tag::AbstractProblemTemplate, ::Type{<:AbstractQuantumTrajectory}) -> Bool

Whether a template supports `free_phase = true` for a given trajectory type.
Defaults to `true`; `@problem_template` emits `false` methods for the declared
exceptions (e.g. `SmoothPulseTemplate` + `KetTrajectory`).
"""
function supports_free_phase end

"""
    template_params_type(tag::AbstractProblemTemplate) -> Type{<:AbstractTemplateParams}

The params struct a template's keywords splat into.
"""
function template_params_type end

"""
    state_dependent_terms(tag::AbstractProblemTemplate) -> Tuple{Vararg{Symbol}}

Names of the template's *state-dependent* physics terms — terms whose value
depends on the current iterate. Declared in the `@problem_template` block; the
macro emits `Live`-only handling for them (see `Live`/`Frozen`). Defaults to `()`.
"""
function state_dependent_terms end

supports_free_phase(::AbstractProblemTemplate, ::Type) = true
state_dependent_terms(::AbstractProblemTemplate) = ()

# Field-wise equality for the params structs. Julia's struct fallback is `===`,
# which is identity-based for the `Vector`/`Dict` fields these carry — so two
# separately built but value-identical params would compare unequal, and
# `extract_spec`'s params verification (and the registry's idempotency check) would
# see phantom differences.
function Base.:(==)(a::AbstractTemplateParams, b::AbstractTemplateParams)
    typeof(a) === typeof(b) || return false
    return all(isequal(getfield(a, f), getfield(b, f)) for f in fieldnames(typeof(a)))
end
function Base.hash(p::AbstractTemplateParams, h::UInt)
    h = hash(typeof(p), h)
    for f in fieldnames(typeof(p))
        h = hash(getfield(p, f), h)
    end
    return h
end

@doc raw"""
    QuantumControlProblem{T<:AbstractProblemTemplate, QT<:AbstractQuantumTrajectory}

Template-tagged quantum control problem: quantum trajectory information plus the
`DirectTrajOptProblem` that optimizes it.

The two type parameters are deliberate and ordered *template first*:

- `T` — the [`AbstractProblemTemplate`](@ref) tag identifying which template built
  the problem. This is what makes `SmoothPulseProblem` a **type** rather than a
  function: `const SmoothPulseProblem{QT<:AbstractQuantumTrajectory{<:ZeroOrderPulse}}
  = QuantumControlProblem{SmoothPulseTemplate, QT}` — the alias bound *is* the
  compatibility matrix, so `SmoothPulseProblem{UnitaryTrajectory{CubicSplinePulse}}`
  is not a type at all.
- `QT` — the quantum trajectory type (itself parametric in the pulse type), so
  dispatch on Unitary/Ket/MultiKet/Density and on the pulse family is free.

Deliberately *not* type parameters: the integrator (it is already its own
parametric type and carries the specialization where it matters) and
`PiccoloOptions` (a runtime value). Type parameters are spent only where dispatch
pays rent.

# Fields
- `qtraj::QT`: quantum trajectory: system, goal, pulse, quantum state
- `prob::DirectTrajOptProblem`: objective, dynamics, constraints
- `params::AbstractTemplateParams`: the retained, typed template keywords
- `spec::Union{Nothing,AbstractProblemSpec}`: the originating declarative spec
  (set by `Specs.materialize`), or `nothing` for hand-built problems

# Construction
Via a problem template (the public call surface, unchanged):
```julia
qtraj = UnitaryTrajectory(sys, pulse, U_goal)
qcp = SmoothPulseProblem(qtraj, N; Q=100.0, R=1e-2)   # ::SmoothPulseProblem
```
or directly, which yields the untagged `NoTemplate` flavor:
```julia
qcp = QuantumControlProblem(qtraj, prob)              # ::QuantumControlProblem{NoTemplate, …}
```

# Accessors
- `get_trajectory(qcp)`: Get the NamedTrajectory
- `get_system(qcp)`: Get the QuantumSystem
- `get_goal(qcp)`: Get the goal state/unitary
- `state_name(qcp)`: Get the state variable name
- `drive_name(qcp)`: Get the control variable name
- `template_tag(qcp)` / `template_params(qcp)` / `retained_spec(qcp)`

# Solving
```julia
solve!(qcp; max_iter=100, verbose=true)
```
"""
mutable struct QuantumControlProblem{
    T<:AbstractProblemTemplate,
    QT<:AbstractQuantumTrajectory,
} <: AbstractQuantumControlProblem
    qtraj::QT
    prob::DirectTrajOptProblem
    params::AbstractTemplateParams
    spec::Union{Nothing,AbstractProblemSpec}
end

"""
    QuantumControlProblem(qtraj, prob)

Build an untagged ([`NoTemplate`](@ref)) problem directly from a quantum
trajectory and a `DirectTrajOptProblem`. The hand-built path; template
constructors produce tagged problems instead.
"""
QuantumControlProblem(qtraj::QT, prob::DirectTrajOptProblem) where {QT} =
    QuantumControlProblem{NoTemplate,QT}(qtraj, prob, NoParams(), nothing)

"""
    QuantumControlProblem{T}(qtraj, prob; params=NoParams(), spec=nothing)

Build a problem tagged with template `T`. Used by the `@problem_template`-generated
constructors and by `Specs.materialize`.
"""
QuantumControlProblem{T}(
    qtraj::QT,
    prob::DirectTrajOptProblem;
    params::AbstractTemplateParams = NoParams(),
    spec::Union{Nothing,AbstractProblemSpec} = nothing,
) where {T<:AbstractProblemTemplate,QT} =
    QuantumControlProblem{T,QT}(qtraj, prob, params, spec)

# ============================================================================= #
# Wrapper problems
# ============================================================================= #

@doc raw"""
    AbstractProblemWrapper <: AbstractQuantumControlProblem

Supertype of the *parametric wrapper problems* — problems built **around** another
problem, e.g.

```julia
SamplingProblem{QuantumControlProblem{SmoothPulseTemplate, UnitaryTrajectory{ZeroOrderPulse}}}
```

which is the type-level mirror of the ordered `[[wrappers]]` list in a spec: a
sampled smooth CZ *says so in its type*. Before Phase 1b the wrappers were
functions returning plain problems, so the wrap history evaporated at construction.

# Field contract

Every subtype carries five fields (this is the interface — the generic accessors
below read them with `getfield`):

| field | meaning |
|---|---|
| `inner` | the wrapped [`AbstractQuantumControlProblem`](@ref) |
| `qtraj` | the wrapper's own quantum trajectory (e.g. the `SamplingTrajectory`) |
| `prob` | the wrapper's `DirectTrajOptProblem` (the NLP actually solved) |
| `params` | the wrapper's typed params |
| `spec` | retained originating spec, or `nothing` |

and implements [`wrapper_kind`](@ref) and [`rewrap`](@ref).

Everything else — `get_trajectory`, `get_system`, `get_goal`, `state_name`,
`drive_name`, `fidelity`, `solve!`, `sync_trajectory!`, property forwarding — is
inherited from the [`AbstractQuantumControlProblem`](@ref) methods.
"""
abstract type AbstractProblemWrapper <: AbstractQuantumControlProblem end

"""
    inner(w::AbstractProblemWrapper) -> AbstractQuantumControlProblem

The problem a wrapper wraps.
"""
inner(w::AbstractProblemWrapper) = getfield(w, :inner)

"""
    base_problem(p) -> QuantumControlProblem

Peel every wrapper layer and return the innermost template-tagged problem.
"""
base_problem(qcp::QuantumControlProblem) = qcp
base_problem(w::AbstractProblemWrapper) = base_problem(inner(w))

"""
    wrapper_kind(w::AbstractProblemWrapper) -> Symbol

The wrapper's registry name (`:sampling`, …) — the `[[wrappers]] kind` this layer
corresponds to.
"""
function wrapper_kind end

"""
    rewrap(w::AbstractProblemWrapper, qtraj, prob) -> typeof(w).name.wrapper

Rebuild a wrapper around a *new* quantum trajectory and NLP, preserving `inner`,
`params` and `spec`. This is how recipes that rewrite the NLP (e.g.
`MinimumTimeProblem`) keep the wrap history instead of flattening it.
"""
function rewrap end

@doc raw"""
    with_problem(p::AbstractQuantumControlProblem, qtraj, prob) -> typeof(p)

Rebuild a problem around a **new** quantum trajectory and NLP while preserving
everything that identifies it: the template tag, the retained params, the retained
spec, and (for wrappers) the wrapped problem.

This is what keeps recipes that rewrite the NLP — `MinimumTimeProblem`, the spec
layer's composition axes — from silently flattening a problem back to the untagged
flavor. Dispatches to [`rewrap`](@ref) for wrappers.
"""
with_problem(qcp::QuantumControlProblem{T}, qtraj, prob::DirectTrajOptProblem) where {T} =
    QuantumControlProblem{T}(
        qtraj,
        prob;
        params = template_params(qcp),
        spec = retained_spec(qcp),
    )
with_problem(w::AbstractProblemWrapper, qtraj, prob::DirectTrajOptProblem) =
    rewrap(w, qtraj, prob)

# ----------------------------------------------------------------------------- #
# The two structural accessors every problem provides
# ----------------------------------------------------------------------------- #

"""
    quantum_trajectory(p::AbstractQuantumControlProblem) -> AbstractQuantumTrajectory

The problem's own quantum trajectory. For a wrapper this is the *wrapper's*
trajectory (e.g. the `SamplingTrajectory`), not the inner problem's.
"""
quantum_trajectory(p::AbstractQuantumControlProblem) = getfield(p, :qtraj)

"""
    direct_problem(p::AbstractQuantumControlProblem) -> DirectTrajOptProblem

The NLP actually handed to the solver.
"""
direct_problem(p::AbstractQuantumControlProblem) = getfield(p, :prob)

# ============================================================================= #
# Tag / params / spec accessors
# ============================================================================= #

"""
    template_tag(qcp) -> AbstractProblemTemplate

The template tag singleton of a problem (`SmoothPulseTemplate()`, …). For a
wrapper, the tag of the problem it wraps.
"""
template_tag(::QuantumControlProblem{T}) where {T} = T()
template_tag(w::AbstractProblemWrapper) = template_tag(inner(w))

"""
    template_tag(::Type{<:QuantumControlProblem}) -> Type{<:AbstractProblemTemplate}

Type-level tag extraction.
"""
template_tag(::Type{<:QuantumControlProblem{T}}) where {T} = T

"""
    template_params(qcp) -> AbstractTemplateParams

The retained, typed template keywords the problem was built with.
"""
template_params(p::AbstractQuantumControlProblem) = getfield(p, :params)

"""
    retained_spec(p) -> Union{Nothing,AbstractProblemSpec}

The originating declarative spec, if the problem was materialized from one.
`nothing` for a problem built by calling a template directly.
"""
retained_spec(p::AbstractQuantumControlProblem) = getfield(p, :spec)

"""
    retain_spec!(p, spec) -> p

Attach the originating declarative spec to a problem. Called by
`Specs.materialize`; `Specs.extract_spec` reads it back *after verifying it against
the materialized object*, which is what makes the spec round-trip a real
consistency check on the materializer rather than a tautology.
"""
function retain_spec!(p::AbstractQuantumControlProblem, spec::AbstractProblemSpec)
    setfield!(p, :spec, spec)
    return p
end

# ============================================================================= #
# Convenience accessors - extend PiccoloQuantumObjects methods
# ============================================================================= #

# Defined once on `AbstractQuantumControlProblem`, so the wrapper types inherit
# them through `quantum_trajectory`/`direct_problem`.

"""
    get_trajectory(p::AbstractQuantumControlProblem)

Get the NamedTrajectory from the optimization problem.
"""
get_trajectory(p::AbstractQuantumControlProblem) = direct_problem(p).trajectory

"""
    get_system(p::AbstractQuantumControlProblem)

Get the QuantumSystem from the quantum trajectory.
"""
get_system(p::AbstractQuantumControlProblem) = get_system(quantum_trajectory(p))

"""
    get_goal(p::AbstractQuantumControlProblem)

Get the goal state/operator from the quantum trajectory.
"""
get_goal(p::AbstractQuantumControlProblem) = get_goal(quantum_trajectory(p))

"""
    state_name(p::AbstractQuantumControlProblem)

Get the state variable name from the quantum trajectory.
"""
state_name(p::AbstractQuantumControlProblem) = state_name(quantum_trajectory(p))

"""
    drive_name(p::AbstractQuantumControlProblem)

Get the control variable name from the quantum trajectory.
"""
drive_name(p::AbstractQuantumControlProblem) = drive_name(quantum_trajectory(p))

"""
    fidelity(p::AbstractQuantumControlProblem; kwargs...)

Compute the fidelity of the quantum trajectory, **applying the problem's optimized free
phases** when it has any.

A `free_phase = true` problem is optimized against the *phase-rotated* goal, so this reads the
`φ_*` globals off `p.trajectory` and forwards them. An explicit `phases = ...` always
wins; pass `phases = zeros(n)` for the fixed-phase number deliberately.

!!! warning "Behaviour changed 2026-07-25"
    This previously forwarded with no phase handling, so every free-phase problem silently
    reported its **fixed-phase** fidelity — a number its objective never optimized. Any
    free-phase figure obtained through this path before that date needs re-deriving.

# Example
```julia
solve!(qcp)
fid = fidelity(qcp)                     # φ-aware when the problem has φ globals
raw = fidelity(qcp; phases = zeros(2))  # fixed-phase, explicitly
```
"""
function fidelity(p::AbstractQuantumControlProblem; phases = nothing, kwargs...)
    if !isnothing(phases)
        return fidelity(quantum_trajectory(p); phases = phases, kwargs...)
    end

    φ = stored_phases(get_trajectory(p))
    isnothing(φ) && return fidelity(quantum_trajectory(p); kwargs...)

    if isempty(φ)
        error(
            "this problem declares free-phase globals (φ_*) but stores no phase values, so a " *
            "free-phase fidelity cannot be computed. Reporting the fixed-phase number instead " *
            "would be the very bug this check exists to prevent. Re-solve the problem, or pass " *
            "`phases = ...` explicitly if you know them.",
        )
    end

    return fidelity(quantum_trajectory(p); phases = φ, kwargs...)
end

# ============================================================================= #
# Forward DirectTrajOptProblem methods
# ============================================================================= #

"""
    sync_trajectory!(p::AbstractQuantumControlProblem)

Update the quantum trajectory in-place from the optimized control values.

After optimization, this function:
1. Extracts the optimized controls from `prob.trajectory` (unadapting if needed)
2. Creates a new pulse with those controls via `extract_pulse`
3. Re-solves the ODE to get the updated quantum evolution
4. Replaces `qtraj` with the new quantum trajectory

This gives you access to the continuous-time ODE solution with the optimized controls,
allowing you to:
- Evaluate the fidelity via `fidelity(qcp.qtraj)`
- Sample the quantum state at any time via `qcp.qtraj(t)`
- Get the optimized pulse via `get_pulse(qcp.qtraj)`

# Example
```julia
solve!(qcp; max_iter=100)  # Automatically calls sync_trajectory!
fid = fidelity(qcp.qtraj)  # Evaluate fidelity with continuous-time solution
pulse = get_pulse(qcp.qtraj)  # Get the optimized pulse
```
"""
function sync_trajectory!(
    p::AbstractQuantumControlProblem;
    check_divergence::Bool = true,
    divergence_rtol::Real = ROLLOUT_DIVERGENCE_RTOL[],
)
    qtraj = quantum_trajectory(p)
    traj = direct_problem(p).trajectory
    # Update global parameters in the system BEFORE rollout so the ODE uses optimized values
    sys = get_system(qtraj)
    if hasproperty(sys, :global_params) &&
       !isempty(sys.global_params) &&
       traj.global_dim > 0
        update_global_params!(qtraj, traj)
    end

    # Extract the optimized pulse from the discrete trajectory and roll it out in-place
    pulse = extract_pulse(qtraj, traj)
    rollout!(qtraj, pulse)

    # Both terminal states are now in hand — the collocation one the optimizer actually
    # minimized against, and the ODE one just produced. Nothing compared them before.
    check_divergence && _warn_on_rollout_divergence(p; rtol = divergence_rtol)

    return nothing
end

# ============================================================================= #
# Rollout divergence diagnostic
# ============================================================================= #

"""
    ROLLOUT_DIVERGENCE_RTOL

Relative tolerance above which [`sync_trajectory!`](@ref) warns that the optimizer's
collocation solution and an actual ODE re-rollout of the same pulse disagree at the
final time. Read as `ROLLOUT_DIVERGENCE_RTOL[]`; set as `ROLLOUT_DIVERGENCE_RTOL[] = 0.05`.

Raise it to quiet the check globally, or pass `check_divergence = false` to silence a
single call.

The default is measured, not guessed: on a 3-level X gate with everything but the
pulse–integrator pairing held fixed, a correct spline pairing gives 3.5e-7 … 6.4e-4 while a
spline under a piecewise-constant integrator gives 5.9e-2 … 1.9e-1. `5e-3` sits near the
geometric centre of that ~90× gap. Tighten it to catch subtler model error, at the cost of
firing on coarse-but-correct problems.
"""
const ROLLOUT_DIVERGENCE_RTOL = Ref(5.0e-3)

# Terminal state(s) of a rolled-out quantum trajectory as `name => iso-vector` pairs, in the
# SAME iso-vector representation the collocation `NamedTrajectory` stores them under, so the
# two are directly comparable without a fidelity in between. Comparing states rather than
# fidelities is deliberate: it needs no goal, no `EmbeddedOperator`, and no phase convention,
# and it is strictly stronger — two different states can share a fidelity. The conversions
# mirror `NamedTrajectory(::AbstractQuantumTrajectory, …)` exactly; note density uses the
# *compact* iso (n² real params), NOT `density_to_iso_vec`.
#
# Multi-state types return one pair per sub-trajectory. Their `qtraj.solution.u` is indexed by
# sub-trajectory rather than by time, so the terminal state of sub-trajectory `i` is
# `solution.u[i].u[end]` — the same access `fidelity(::MultiKetTrajectory)` uses.
#
# `nothing` means the check is NOT APPLICABLE, never that there is no divergence.
_terminal_iso_states(qtraj::UnitaryTrajectory) =
    [state_name(qtraj) => operator_to_iso_vec(qtraj.solution.u[end])]
_terminal_iso_states(qtraj::KetTrajectory) =
    [state_name(qtraj) => ket_to_iso(qtraj.solution.u[end])]
_terminal_iso_states(qtraj::DensityTrajectory) =
    [state_name(qtraj) => density_to_compact_iso(qtraj.solution.u[end])]

_terminal_iso_states(qtraj::MultiKetTrajectory) = [
    name => ket_to_iso(qtraj.solution.u[i].u[end]) for
    (i, name) in enumerate(state_names(qtraj))
]
_terminal_iso_states(qtraj::MultiDensityTrajectory) = [
    name => density_to_compact_iso(qtraj.solution.u[i].u[end]) for
    (i, name) in enumerate(state_names(qtraj))
]

# `rollout = :none` placeholder states are the tiled initial-ket guess, not solved dynamics, so
# a divergence computed against them would be meaningless rather than merely inaccurate.
_terminal_iso_states(qtraj::MultiKetTrajectory{<:AbstractPulse,RolloutStates}) =
    qtraj.solution.real_states ?
    [
        name => ket_to_iso(qtraj.solution.u[i].u[end]) for
        (i, name) in enumerate(state_names(qtraj))
    ] : nothing

_terminal_iso_states(::AbstractQuantumTrajectory) = nothing

@doc raw"""
    rollout_divergence(qcp::QuantumControlProblem) -> Union{Nothing,Float64}

Relative disagreement at the final time between the optimizer's **collocation** solution
and an **ODE re-rollout** of the same extracted pulse:

```math
\varepsilon = \frac{\lVert x^{\text{rollout}}_N - x^{\text{collocation}}_N\rVert_2}
                   {\max\!\left(\lVert x^{\text{collocation}}_N\rVert_2,\, 1\right)}
```

Multi-state trajectories (`MultiKetTrajectory`, `MultiDensityTrajectory`) are covered: the
result is the 2-norm of the stacked difference over the 2-norm of the stacked collocation
state, which reduces exactly to the single-state formula when there is one sub-state, so the
numbers are on a common scale.

Returns `nothing` — meaning *not applicable*, never *no divergence* — when the check cannot be
made: a state component absent from the optimizer's trajectory, a shape mismatch, a trajectory
type without a defined terminal state, or a `rollout = :none` `MultiKetTrajectory` still
holding placeholder states rather than solved dynamics.

Call this only after [`sync_trajectory!`](@ref) (or `solve!` with `sync = true`);
before that, `qcp.qtraj` holds a stale rollout and the number is meaningless.

When this is large, the optimizer's objective and `fidelity(qcp.qtraj)` describe **different
waveforms** and only the rollout one is physical. Usual causes: a pulse/integrator mismatch, or
a collocation grid too coarse for the dynamics.
"""
function rollout_divergence(qcp::AbstractQuantumControlProblem)
    pairs = _terminal_iso_states(qcp.qtraj)
    isnothing(pairs) && return nothing

    traj = qcp.prob.trajectory
    sq_diff = 0.0
    sq_collocation = 0.0

    for (name, x_rollout) in pairs
        haskey(traj.components, name) || return nothing
        x_collocation = traj[name][:, end]
        length(x_collocation) == length(x_rollout) || return nothing
        sq_diff += sum(abs2, x_rollout - x_collocation)
        sq_collocation += sum(abs2, x_collocation)
    end

    # Aggregated over sub-states, this is the 2-norm of the stacked difference over the 2-norm
    # of the stacked collocation state — identical to the single-state formula when there is
    # only one pair, so multi-state numbers are on the same scale as single-state ones.
    return sqrt(sq_diff) / max(sqrt(sq_collocation), 1.0)
end

function _warn_on_rollout_divergence(
    qcp::AbstractQuantumControlProblem;
    rtol::Real = ROLLOUT_DIVERGENCE_RTOL[],
)
    ε = rollout_divergence(qcp)
    (isnothing(ε) || !isfinite(ε) || ε ≤ rtol) && return nothing

    @warn """
          Optimizer and rollout disagree at the final time.

          relative divergence = $(round(ε; sigdigits = 3))  (tolerance $rtol)

          `fidelity(qcp.qtraj)` — the ODE re-rollout — is the trustworthy number. The
          optimizer's own objective was evaluated on the collocation solution, which is a
          different waveform, so any fidelity taken from it is not what this pulse does.

          Two usual causes:
            • Pulse/integrator mismatch. A spline pulse under a piecewise-constant
              integrator is the common case: the interpolation the optimizer constrained is
              not the one being integrated. Pair a spline pulse with
              `Piccolissimo.SplineIntegrator`.
            • Collocation grid too coarse for the dynamics — increase the number of knots.

          Silence: `sync_trajectory!(qcp; check_divergence = false)`, or raise
          `ROLLOUT_DIVERGENCE_RTOL[]`.
          """
    return nothing
end

"""
    solve!(p::AbstractQuantumControlProblem; sync::Bool=true, verbose::Bool=false, kwargs...)

Solve the quantum control problem by forwarding to the inner DirectTrajOptProblem.

# Arguments
- `sync::Bool=true`: If true, call `sync_trajectory!` after solving to update `qtraj.trajectory`
  with physical control values. Set to false to skip synchronization (e.g., for debugging).
- `verbose::Bool=false`: Controls the solver setup trace (evaluator construction, jacobian/hessian
  structure, NLP block assembly). Defaults to `false` so the log stays clean. Ipopt's iteration
  log is controlled separately by `print_level` (passed through to Ipopt).
- `check_divergence::Bool=true`: Warn if the optimizer's collocation solution and the ODE
  re-rollout disagree at the final time — see [`rollout_divergence`](@ref). Only has an
  effect when `sync = true`.
- `divergence_rtol::Real=ROLLOUT_DIVERGENCE_RTOL[]`: Relative tolerance for that check.

All other keyword arguments are passed to the DirectTrajOpt solver.

Returns `DirectTrajOpt.Solvers.SolveStats` (termination status, raw status string,
NLP objective, IPM iterations, solve wall time, solver symbol) —
computed after the trajectory sync, describing the trajectory you hold.
"""
function solve!(
    p::AbstractQuantumControlProblem;
    sync::Bool = true,
    verbose::Bool = false,
    check_divergence::Bool = true,
    divergence_rtol::Real = ROLLOUT_DIVERGENCE_RTOL[],
    kwargs...,
)
    # SolveStats (status, iterations, objective, wall time) from the backend —
    # returned AFTER the trajectory sync so callers get data describing the
    # trajectory they actually hold.
    stats = solve!(direct_problem(p); verbose = verbose, kwargs...)
    if sync
        sync_trajectory!(
            p;
            check_divergence = check_divergence,
            divergence_rtol = divergence_rtol,
        )
    end
    return stats
end

# Forward other common DirectTrajOptProblem accessors
Base.getproperty(p::AbstractQuantumControlProblem, s::Symbol) = begin
    if s === :qtraj || s === :prob || s === :params || s === :spec || s === :inner
        getfield(p, s)
        # Forward to prob for common fields
    elseif s in (:objective, :dynamics, :constraints, :trajectory)
        getproperty(getfield(p, :prob), s)
    else
        # Fall back to default behavior
        getfield(p, s)
    end
end

# ============================================================================= #
# Display
# ============================================================================= #
#
# `Base.show` for QuantumControlProblem is defined in
# `control/display/show.jl`. The rich display lives in submodule
# `ProblemDisplay`, which is loaded after this file so it can override the
# default.

# ============================================================================= #
# Tests
# ============================================================================= #

# NOTE: the pulse/integrator-pairing regression test does NOT live here.
#
# Demonstrating that `rollout_divergence` separates a correct pairing from a mismatch
# requires BOTH a piecewise-constant and a spline integrator, and public Piccolo ships only
# the former (`BilinearIntegrator`). Building the comparison here would mean treating
# `BilinearIntegrator` as the object of study, which it is not — it is the public default,
# not the integrator this stack actually optimizes against.
#
# The pairing test therefore lives in Piccolissimo, where `SplineIntegrator` and
# `HermitianExponentialIntegrator` both exist. Measured there (3-level X gate, objective /
# system / goal / grid / seed held fixed, only the pairing varied):
#
#   spline pulse + SplineIntegrator (correct)   divergence 3.5e-7 .. 6.4e-4
#   spline pulse + HermitianExp     (PWC)       divergence 5.9e-2 .. 1.9e-1
#
# and, more tellingly, the mismatched run's collocation infidelity reads ~1e-8 while its
# actual rolled-out infidelity is ~4e-3: the optimizer believes it converged to machine
# precision. That ~90x divergence gap is what sets `ROLLOUT_DIVERGENCE_RTOL`.
#
# The two testitems below are integrator-agnostic on purpose: they force a divergence by
# mutating the collocation state directly, so the integrator is only scaffolding.

@testitem "fidelity(qcp) applies the problem's stored free phases" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    # Regression guard for the 2026-07-25 decision that `fidelity(qcp)` becomes φ-aware.
    # Before it, this path silently reported the FIXED-phase number for a free-phase problem —
    # a value the objective never optimized. Notably, no test in the suite combined
    # `free_phase = true` with `fidelity(qcp)`, which is why the bug survived.
    σx = ComplexF64[0 1; 1 0]
    H_drift = ComplexF64[0 0 0; 0 1 0; 0 0 2]
    H_drive = ComplexF64[0 1 0; 1 0 1; 0 1 0] / √2
    sys = OpenQuantumSystem(H_drift, [H_drive], [1.0])

    N, T = 21, 10.0
    times = collect(range(0.0, T, length = N))
    pulse = LinearSplinePulse(0.1 * randn(1, N), times)
    U_goal = EmbeddedOperator(σx, [1, 2], [3])

    qtraj = UnitaryTrajectory(sys, pulse, U_goal)
    qcp = SplinePulseProblem(qtraj, N; Q = 100.0, R = 1e-2, free_phase = true)
    traj = get_trajectory(qcp)
    @test haskey(traj.global_components, :φ_1)

    # Put a deliberately non-trivial phase in, then sync so the rollout is current.
    traj.global_data[traj.global_components[:φ_1]] .= [0.7]
    sync_trajectory!(qcp; check_divergence = false)

    φ = stored_phases(traj)
    @test φ == [0.7]

    F_default = fidelity(qcp)
    F_fixed = fidelity(qcp; phases = [0.0])

    # The default now equals the explicitly-φ-aware value, not the fixed-phase one.
    @test F_default ≈ fidelity(qcp.qtraj; phases = φ)
    @test !isapprox(F_default, F_fixed; atol = 1e-10)

    # An explicit `phases` still wins over the stored value.
    @test fidelity(qcp; phases = [0.0]) ≈ fidelity(qcp.qtraj; phases = [0.0])
end

@testitem "stored_phases distinguishes no-φ from φ-declared-but-empty" begin
    using DirectTrajOpt
    using NamedTrajectories

    sys = OpenQuantumSystem(zeros(ComplexF64, 2, 2), [ComplexF64[0 1; 1 0]], [(-2.0, 2.0)])
    N, T = 11, 5.0
    times = collect(range(0, T, length = N))
    pulse = LinearSplinePulse(zeros(1, N), times)
    qtraj = KetTrajectory(sys, pulse, ComplexF64[1, 0], ComplexF64[0, 1])
    traj = NamedTrajectory(qtraj, N; Δt_bounds = (1e-3, 1.0))
    prob = DirectTrajOptProblem(
        traj,
        QuadraticRegularizer(:u, traj, 1.0),
        BilinearIntegrator(qtraj, N),
    )
    qcp = QuantumControlProblem(qtraj, prob)
    sync_trajectory!(qcp; check_divergence = false)

    # No φ globals at all => `nothing`, and the ordinary (non-free-phase) path is unaffected.
    # This is the 23-of-32 case and must stay byte-identical to the pre-φ behaviour.
    @test isnothing(stored_phases(qcp.prob.trajectory))
    @test fidelity(qcp) ≈ fidelity(qcp.qtraj)

    # An EMPTY phase vector is the corrupt-input case — the state the free-phase catalog entries
    # that never persisted φ are in. It must be refused, not silently treated as fixed-phase,
    # because returning a plausible number there is the whole bug. `verify` shares the guard and
    # takes `phases` explicitly, so it is the drivable entry point for it.
    err = try
        verify(qcp; phases = Float64[])
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("no phase values", err.msg)
end

@testitem "rollout_divergence: nothing means not-applicable, not agreement" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    sys = OpenQuantumSystem(zeros(ComplexF64, 2, 2), [ComplexF64[0 1; 1 0]], [(-2.0, 2.0)])
    N, T = 11, 5.0
    times = collect(range(0, T, length = N))
    pulse = LinearSplinePulse(zeros(1, N), times)

    qtraj = KetTrajectory(sys, pulse, ComplexF64[1, 0], ComplexF64[0, 1])
    traj = NamedTrajectory(qtraj, N)
    prob = DirectTrajOptProblem(
        traj,
        QuadraticRegularizer(:u, traj, 1.0),
        BilinearIntegrator(qtraj, N),
    )
    qcp = QuantumControlProblem(qtraj, prob)

    # A ket trajectory IS supported, so after a sync the divergence is a real number.
    sync_trajectory!(qcp; check_divergence = false)
    ε = rollout_divergence(qcp)
    @test ε isa Float64
    @test ε ≥ 0.0

    # Zero controls: collocation and ODE both sit at |0⟩, so they must agree.
    @test ε < 1.0e-8
end

@testitem "rollout_divergence covers multi-state trajectories" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    # MultiKet/MultiDensity index `solution.u` by sub-trajectory rather than by time, so they
    # were initially skipped entirely. They are the population most exposed to the pairing
    # hazard, so being skipped meant the riskiest problems got no check at all.
    σx = ComplexF64[0 1; 1 0]
    sys = OpenQuantumSystem(0.01 * ComplexF64[1 0; 0 -1], [σx], [1.0])
    N, T = 21, 10.0
    times = collect(range(0.0, T, length = N))
    pulse = LinearSplinePulse(0.1 * randn(1, N), times)
    ψ0, ψ1 = ComplexF64[1.0, 0.0], ComplexF64[0.0, 1.0]

    qtraj = MultiKetTrajectory(sys, pulse, [ψ0, ψ1], [ψ1, ψ0])
    qcp = SplinePulseProblem(qtraj, N; Q = 100.0, R = 1e-2)
    sync_trajectory!(qcp; check_divergence = false)

    ε = rollout_divergence(qcp)
    @test ε isa Float64                     # NOT `nothing` any more
    @test ε ≥ 0.0
    @test isfinite(ε)

    # Perturbing ONE sub-state's collocation terminal value must move the aggregate.
    traj = get_trajectory(qcp)
    @test haskey(traj.components, :ψ̃1) && haskey(traj.components, :ψ̃2)
    traj.ψ̃1[:, end] .+= 0.5
    ε_perturbed = rollout_divergence(qcp)
    @test ε_perturbed > ε
end

@testitem "sync_trajectory!: divergence check is opt-out and threshold-driven" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra
    using Logging

    sys = OpenQuantumSystem(zeros(ComplexF64, 2, 2), [ComplexF64[0 1; 1 0]], [(-2.0, 2.0)])
    N, T = 11, 5.0
    times = collect(range(0, T, length = N))
    pulse = LinearSplinePulse(zeros(1, N), times)

    qtraj = KetTrajectory(sys, pulse, ComplexF64[1, 0], ComplexF64[0, 1])
    traj = NamedTrajectory(qtraj, N)
    prob = DirectTrajOptProblem(
        traj,
        QuadraticRegularizer(:u, traj, 1.0),
        BilinearIntegrator(qtraj, N),
    )
    qcp = QuantumControlProblem(qtraj, prob)

    # Force a divergence: move the collocation terminal state away from the ODE result
    # without touching the controls, so the re-rollout cannot follow it.
    qcp.prob.trajectory.ψ̃[:, end] .= Float64[0, 1, 0, 0]

    @test_logs (:warn, r"disagree at the final time") sync_trajectory!(qcp)

    # ...and the same call is silent when the check is disabled.
    qcp.prob.trajectory.ψ̃[:, end] .= Float64[0, 1, 0, 0]
    @test_logs min_level = Logging.Warn sync_trajectory!(qcp; check_divergence = false)

    # ...and silent when the tolerance is loose enough to accept it.
    qcp.prob.trajectory.ψ̃[:, end] .= Float64[0, 1, 0, 0]
    @test_logs min_level = Logging.Warn sync_trajectory!(qcp; divergence_rtol = 10.0)
end

@testitem "sync_trajectory! updates quantum trajectory" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    # Create a simple quantum system with X drive
    levels = 2
    H_drift = zeros(ComplexF64, levels, levels)
    σx = ComplexF64[0 1; 1 0]

    N = 11
    T = 5.0
    sys = OpenQuantumSystem(H_drift, [σx], [(-2.0, 2.0)]; time_dependent = false)

    ψ_init = ComplexF64[1, 0]
    ψ_target = ComplexF64[0, 1]

    # Create pulse with zero controls
    times = collect(range(0, T, length = N))
    controls = zeros(1, N)
    pulse = LinearSplinePulse(controls, times)

    # Create initial trajectory
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_target)

    # Verify initial fidelity is low (controls are zero)
    initial_fid = fidelity(qtraj)
    @test initial_fid < 0.1  # Should be near 0 since |0⟩ stays at |0⟩

    # Create a simple problem - convert QT to NamedTrajectory
    traj = NamedTrajectory(qtraj, N)
    obj = QuadraticRegularizer(:u, traj, 1.0)
    integrator = BilinearIntegrator(qtraj, N)
    prob = DirectTrajOptProblem(traj, obj, integrator)
    qcp = QuantumControlProblem(qtraj, prob)

    # Manually modify prob.trajectory to simulate optimization
    # Set u to a constant value that will rotate |0⟩ toward |1⟩
    # For a simple rotation: exp(-i * σx * u * t) |0⟩ = cos(ut)|0⟩ - i*sin(ut)|1⟩
    # After time T with u = π/(2T), we get |1⟩
    u_opt = π / (2 * T)
    qcp.prob.trajectory.u .= u_opt

    # Call sync to update qtraj with new controls
    sync_trajectory!(qcp)

    # The qtraj should now have the optimized pulse
    new_pulse = get_pulse(qcp.qtraj)
    # Access underlying data from the interpolator (.u field in DataInterpolations)
    @test all(new_pulse.controls.u .≈ u_opt)

    # The fidelity should be much better
    final_fid = fidelity(qcp.qtraj)
    @test final_fid > 0.9  # Should be near 1 now
end

@testitem "solve! with sync=true updates trajectory" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    # Create a minimal quantum system
    levels = 2
    H_drift = zeros(ComplexF64, levels, levels)
    σx = ComplexF64[0 1; 1 0]

    ψ_init = ComplexF64[1, 0]
    ψ_target = ComplexF64[0, 1]
    N = 11
    T = 5.0

    sys = OpenQuantumSystem(H_drift, [σx], [(-2.0, 2.0)]; time_dependent = false)

    # Create pulse with zero controls
    times = collect(range(0, T, length = N))
    controls = zeros(1, N)
    pulse = LinearSplinePulse(controls, times)

    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_target)
    traj = NamedTrajectory(qtraj, N)

    # Create problem with simple objective
    obj = QuadraticRegularizer(:u, traj, 0.01)
    integrator = BilinearIntegrator(qtraj, N)
    prob = DirectTrajOptProblem(traj, obj, integrator)
    qcp = QuantumControlProblem(qtraj, prob)

    # Store original pulse
    original_pulse = get_pulse(qcp.qtraj)

    # Solve with max_iter=0 (no optimization, just test sync mechanism)
    solve!(qcp; max_iter = 0, sync = true)

    # qtraj should be updated (in-place, but new solution)
    @test true  # If we get here without errors, sync worked

    # Test sync=false doesn't update trajectory
    qtraj2 = KetTrajectory(sys, pulse, ψ_init, ψ_target)
    traj2 = NamedTrajectory(qtraj2, N)
    prob2 = DirectTrajOptProblem(traj2, obj, integrator)
    qcp2 = QuantumControlProblem(qtraj2, prob2)
    original_qtraj2 = qcp2.qtraj

    solve!(qcp2; max_iter = 0, sync = false)
    @test qcp2.qtraj === original_qtraj2  # Same object, not rebuilt
end

@testitem "extract_pulse and rollout creates new trajectory with updated pulse" begin
    using LinearAlgebra
    using NamedTrajectories

    # Create system
    levels = 2
    H_drift = zeros(ComplexF64, levels, levels)
    σx = ComplexF64[0 1; 1 0]
    N = 11
    T = 5.0

    sys = OpenQuantumSystem(H_drift, [σx], [(-2.0, 2.0)]; time_dependent = false)

    ψ_init = ComplexF64[1, 0]
    ψ_goal = ComplexF64[0, 1]

    # Create pulse with zero controls
    times = collect(range(0, T, length = N))
    controls = zeros(1, N)
    pulse = LinearSplinePulse(controls, times)

    # Create initial trajectory
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

    # Get NamedTrajectory and modify controls
    traj = NamedTrajectory(qtraj, N)
    u_opt = π / (2 * T)

    # Create modified trajectory data
    new_u = fill(u_opt, size(traj.u))
    new_traj = NamedTrajectory(
        (; ψ̃ = traj.ψ̃, t = traj.t, Δt = traj.Δt, u = new_u);
        timestep = :Δt,
        controls = (:Δt, :u),
        bounds = traj.bounds,
        initial = traj.initial,
        goal = traj.goal,
    )

    # Extract pulse and rollout with new controls
    new_pulse = extract_pulse(qtraj, new_traj)
    new_qtraj = rollout(qtraj, new_pulse)

    # Check pulse was updated (access underlying data via .u)
    @test all(new_qtraj.pulse.controls.u .≈ u_opt)

    # Check ODE was re-solved (fidelity should be high)
    @test fidelity(new_qtraj) > 0.9

    # Original trajectory unchanged (access underlying data via .u)
    @test all(qtraj.pulse.controls.u .≈ 0.0)
end

end
