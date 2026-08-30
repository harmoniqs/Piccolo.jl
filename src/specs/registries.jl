module SpecRegistries

using ..SpecStructs
using ...Quantum          # QuantumSystemTemplates, Gates, etc.
using ...Control          # ProblemTemplates: SmoothPulseProblem, ...
using ...Control.QuantumIntegrators.ExponentialIntegrators:
    HermitianExponentialIntegrator, NonHermitianExponentialIntegrator
using ...Control: Control

export RegistryEntry,
    ConstFactory,
    reflect_params,
    register_templates_from_declarations!,
    SYSTEMS,
    TEMPLATES,
    INTEGRATORS,
    WRAPPERS,
    OBJECTIVE_TERMS,
    SOLVERS,
    SOLVE_STRATEGIES,
    PROBLEM_KINDS,
    register_system!,
    register_template!,
    register_integrator!,
    register_wrapper!,
    register_objective_term!,
    register_solver!,
    register_strategy!,
    register_kind!,
    lookup_system,
    lookup_template,
    lookup_integrator,
    lookup_wrapper,
    lookup_objective_term,
    lookup_solver,
    lookup_strategy,
    lookup_kind,
    register_all!

"""
    ConstFactory(value) <: Function

A stable, value-comparable factory that ignores its arguments and returns
`value`. Used for the registry's sentinel/enum-echo entries (kinds, solvers,
strategies, objective-term placeholders) so that re-registering the same logical
entry across repeated `register_all!()` calls is a byte-identical no-op — unlike
a fresh anonymous closure, `ConstFactory(:x) == ConstFactory(:x)`.
"""
struct ConstFactory <: Function
    value::Any
end
(f::ConstFactory)(args...; kwargs...) = f.value
Base.:(==)(a::ConstFactory, b::ConstFactory) = a.value == b.value
Base.hash(f::ConstFactory, h::UInt) = hash(f.value, hash(ConstFactory, h))

# `bilinear` sentinel: the templates build BilinearIntegrator internally when
# `integrator=nothing`, so this factory is never actually constructed — it only
# marks "pass integrator=nothing" (we never construct BilinearIntegrator here).
_bilinear_integrator_factory(args...; kwargs...) = nothing

# MultiTransmonSystem takes POSITIONAL (ωs, δs, gs); the materializer calls
# system factories kwargs-only (`entry.factory(; params...)`), so wrap the
# constructor (#297). A named function, not a closure: `register_all!`
# idempotency compares factories with `==`, and anonymous closures are never
# `==` across calls. TOML/JSON carry the coupling matrix `gs` as nested row
# vectors — normalize to the Matrix the constructor's `@assert size(gs)` expects.
function _multi_transmon_system_factory(; ωs, δs, gs, kwargs...)
    gs_matrix = gs isa AbstractMatrix ? gs : Matrix{Float64}(permutedims(reduce(hcat, gs)))
    return MultiTransmonSystem(
        Float64.(collect(ωs)),
        Float64.(collect(δs)),
        gs_matrix;
        kwargs...,
    )
end

"""
    RegistryEntry(; factory, params, compat)

A single registry record: the `factory` (a function, or — for the problem
templates, which are now constrained *type aliases* — a `Type` that materializes
the object), a `params` schema (`name => (; type, default)`), and `compat`
metadata (e.g. `:pulse_kinds`, `:trajectory_kinds`, `:ket_free_phase`). The entry
is the single source of truth reflected by `emit_schema`.

Template entries are **generated** from `@problem_template`'s
`TEMPLATE_DECLARATIONS` records (Phase 1b), so the Julia types, the registry, and
the emitted JSON Schema cannot drift apart.
"""
Base.@kwdef struct RegistryEntry
    factory::Union{Function,Type}
    params::Dict{Symbol,Any} = Dict{Symbol,Any}()   # name → (; type, default, range, doc)
    compat::Dict{Symbol,Any} = Dict{Symbol,Any}()   # e.g. :pulse_kinds, :trajectory_kinds, :requires_nonbilinear
end

"""
    reflect_params(PT::Type) -> Dict{Symbol,Any}

Reflect a template params struct into the registry's parameter schema:
`fieldnames` × `fieldtypes` × the `@kwdef` defaults (read off a
default-constructed prototype). This replaces hand-written per-template parameter
schemas — the struct *is* the declaration.
"""
function reflect_params(PT::Type)
    out = Dict{Symbol,Any}()
    proto = PT()
    for (name, T) in zip(fieldnames(PT), fieldtypes(PT))
        out[name] = (; type = T, default = getfield(proto, name))
    end
    return out
end
const SYSTEMS = Dict{Symbol,RegistryEntry}()
const TEMPLATES = Dict{Symbol,RegistryEntry}()
const INTEGRATORS = Dict{Symbol,RegistryEntry}()
const WRAPPERS = Dict{Symbol,RegistryEntry}()
const OBJECTIVE_TERMS = Dict{Symbol,RegistryEntry}()
const SOLVERS = Dict{Symbol,RegistryEntry}()
const SOLVE_STRATEGIES = Dict{Symbol,RegistryEntry}()
const PROBLEM_KINDS = Dict{Symbol,RegistryEntry}()

# idempotent register with conflict detection
function _register!(reg::Dict, name::Symbol, entry::RegistryEntry)
    if haskey(reg, name)
        reg[name] === entry && return entry                       # same object → no-op
        _entry_equal(reg[name], entry) && return entry            # value-identical → no-op
        error("conflicting registration for $(name); a different entry already exists")
    end
    reg[name] = entry
end
_entry_equal(a::RegistryEntry, b::RegistryEntry) =
    a.factory == b.factory && a.params == b.params && a.compat == b.compat

for (fn, reg) in (
    (:register_system!, SYSTEMS),
    (:register_template!, TEMPLATES),
    (:register_integrator!, INTEGRATORS),
    (:register_wrapper!, WRAPPERS),
    (:register_objective_term!, OBJECTIVE_TERMS),
    (:register_solver!, SOLVERS),
    (:register_strategy!, SOLVE_STRATEGIES),
    (:register_kind!, PROBLEM_KINDS),
)
    @eval $fn(name::Symbol, e::RegistryEntry) = _register!($reg, name, e)
end
for (fn, reg) in (
    (:lookup_system, SYSTEMS),
    (:lookup_template, TEMPLATES),
    (:lookup_integrator, INTEGRATORS),
    (:lookup_wrapper, WRAPPERS),
    (:lookup_objective_term, OBJECTIVE_TERMS),
    (:lookup_solver, SOLVERS),
    (:lookup_strategy, SOLVE_STRATEGIES),
    (:lookup_kind, PROBLEM_KINDS),
)
    @eval $fn(name::Symbol) = get($reg, name, nothing)
end

"""
    register_templates_from_declarations!()

Turn every `@problem_template` declaration record into a `TEMPLATES`
[`RegistryEntry`](@ref): `factory` = the constrained type alias, `params` =
[`reflect_params`](@ref) of the declared params struct, `compat` = the declared
pulse/trajectory matrix. Idempotent — the record is value-identical across calls,
so re-registering is a no-op.

This is the anti-drift seam: there is exactly one place a template's
compatibility matrix is written (the `@problem_template` block), and the registry
plus the emitted JSON Schema are *derived* from it.
"""
function register_templates_from_declarations!()
    for (name, d) in Control.ProblemTemplates.TEMPLATE_DECLARATIONS
        register_template!(
            name,
            RegistryEntry(;
                factory = d.alias_type,
                params = reflect_params(d.params_type),
                compat = Dict{Symbol,Any}(
                    :pulse_kinds => copy(d.pulse_kinds),
                    :trajectory_kinds => copy(d.trajectory_kinds),
                    :ket_free_phase => d.ket_free_phase,
                ),
            ),
        )
    end
    return nothing
end

"""
    register_all!()

Register Piccolo's own systems, problem templates, the `bilinear` integrator
sentinel, the `sampling` wrapper, Piccolo-side objective terms, the `ipopt`
solver, the `direct` strategy, and the `control`/`rollout` kinds. Idempotent —
safe to call repeatedly (re-registering a value-identical entry is a no-op).
Piccolissimo augments these at load time via `register_all_piccolissimo!`.
"""
function register_all!()
    # systems (factory = the template constructor; params/compat hand-declared)
    for (name, fac) in (
        (:TransmonSystem, TransmonSystem),
        (:MultiTransmonSystem, _multi_transmon_system_factory),
        (:TransmonCavitySystem, TransmonCavitySystem),
        (:RydbergChainSystem, RydbergChainSystem),
        (:IonChainSystem, IonChainSystem),
        (:RadialMSGateSystem, RadialMSGateSystem),
        (:CatSystem, CatSystem),
    )
        register_system!(name, RegistryEntry(; factory = fac))
    end
    # templates — GENERATED from the `@problem_template` declaration records, so a
    # declaration change moves the registry (and therefore the emitted schema) with
    # it. Piccolissimo's own declarations land in the same dict when it loads.
    register_templates_from_declarations!()
    # integrators — bilinear (Piccolo's own) plus the exponential family moved
    # in from Piccolissimo (open-core slice 3a, Piccolissimo#429); :spline now
    # registers here too — the SplineIntegrator struct + dense cells moved in
    # slice 3b (Piccolissimo#430). The matrix-free cells stay proprietary; the
    # factory closure below is Piccolissimo-free (dense construction only).
    register_integrator!(:bilinear, RegistryEntry(; factory = _bilinear_integrator_factory))  # sentinel: template default path
    register_integrator!(
        :spline,
        RegistryEntry(;
            factory = (sqtraj, N; alg = nothing) -> SplineIntegrator(sqtraj, N),
        ),
    )
    register_integrator!(
        :hermitian_exponential,
        RegistryEntry(;
            factory = (sqtraj, N; alg = nothing) ->
                HermitianExponentialIntegrator(sqtraj, N),
        ),
    )
    register_integrator!(
        :nonhermitian_exponential,
        RegistryEntry(;
            factory = (sqtraj, N; alg = nothing) ->
                NonHermitianExponentialIntegrator(sqtraj, N),
        ),
    )
    # wrappers
    register_wrapper!(:sampling, RegistryEntry(; factory = SamplingProblem))
    # objective terms available in Piccolo (Piccolissimo adds hermite_* in Task 14)
    for k in (:time, :reg_u, :reg_du, :reg_ddu, :leakage, :sensitivity)
        register_objective_term!(k, RegistryEntry(; factory = ConstFactory(k)))
    end
    # solvers
    register_solver!(:ipopt, RegistryEntry(; factory = ConstFactory(:ipopt)))
    # strategies
    register_strategy!(:direct, RegistryEntry(; factory = ConstFactory(:direct)))
    # kinds
    register_kind!(:control, RegistryEntry(; factory = ConstFactory(:control)))
    register_kind!(:rollout, RegistryEntry(; factory = ConstFactory(:rollout)))
    return nothing
end

end

@testitem "registries: register, lookup, idempotency, conflict" begin
    using Piccolo.Specs
    r = Specs.RegistryEntry(;
        factory = () -> 42,
        params = Dict(:Q => (; type = :Float64, default = 100.0)),
        compat = Dict(:pulse_kinds => [:zero_order]),
    )
    Specs.register_template!(:MyTmpl, r)
    @test Specs.lookup_template(:MyTmpl).factory() == 42
    Specs.register_template!(:MyTmpl, r)          # idempotent no-op
    @test_throws Exception Specs.register_template!(
        :MyTmpl,
        Specs.RegistryEntry(; factory = () -> 0),
    )  # conflicting → error
    # Piccolo's own registrations present after register_all!
    Specs.register_all!()
    @test haskey(Specs.TEMPLATES, :SmoothPulseProblem)
    @test haskey(Specs.SYSTEMS, :TransmonSystem)
    @test Specs.lookup_kind(:control) !== nothing
    # The registries are process-global and TestItemRunner reuses workers, so the
    # ad-hoc :MyTmpl registration above must not survive this test item: a stray
    # template name leaks into `emit_schema()`'s template enum and its per-template
    # parameter branches, which would break the schema drift gate depending on which
    # test item ran first on the worker.
    delete!(Specs.TEMPLATES, :MyTmpl)
    @test !haskey(Specs.TEMPLATES, :MyTmpl)
end
