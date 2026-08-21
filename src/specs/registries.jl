module SpecRegistries

using ..SpecStructs
using ...Quantum          # QuantumSystemTemplates, Gates, etc.
using ...Control          # ProblemTemplates: SmoothPulseProblem, ...

export RegistryEntry,
    ConstFactory,
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

A single registry record: the `factory` (a function/callable that materializes
the object), a hand-declared `params` schema (`name => (; type, default, range,
doc)`), and `compat` metadata (e.g. `:pulse_kinds`, `:trajectory_kinds`,
`:ket_free_phase`). The entry is the single source of truth reflected by
`emit_schema`; when the `@problem_template` macro lands (Phase 1b) it will
*generate* this same entry.
"""
Base.@kwdef struct RegistryEntry
    factory::Function
    params::Dict{Symbol,Any} = Dict{Symbol,Any}()   # name → (; type, default, range, doc)
    compat::Dict{Symbol,Any} = Dict{Symbol,Any}()   # e.g. :pulse_kinds, :trajectory_kinds, :requires_nonbilinear
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
    # templates (factory = existing function; compat = pulse/trajectory matrix)
    register_template!(
        :SmoothPulseProblem,
        RegistryEntry(;
            factory = SmoothPulseProblem,
            compat = Dict(
                :pulse_kinds => [:zero_order],
                :trajectory_kinds => [:unitary, :ket],
                :ket_free_phase => false,
            ),
        ),
    )
    register_template!(
        :SplinePulseProblem,
        RegistryEntry(;
            factory = SplinePulseProblem,
            compat = Dict(
                :pulse_kinds => [:linear_spline, :cubic_spline],
                :trajectory_kinds => [:unitary, :ket],
                :ket_free_phase => true,
            ),
        ),
    )
    register_template!(
        :BangBangPulseProblem,
        RegistryEntry(;
            factory = BangBangPulseProblem,
            compat = Dict(
                :pulse_kinds => [:zero_order],
                :trajectory_kinds => [:unitary, :ket],
                :ket_free_phase => true,
            ),
        ),
    )
    # integrators — bilinear is Piccolo's; exponential/spline are Piccolissimo's (Task 14)
    register_integrator!(:bilinear, RegistryEntry(; factory = _bilinear_integrator_factory))  # sentinel: template default path
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
