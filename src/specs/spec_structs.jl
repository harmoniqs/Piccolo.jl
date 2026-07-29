module SpecStructs

export ProblemSpec,
    SystemSpec,
    GoalSpec,
    PulseSpec,
    TrajectorySpec,
    TemplateBlock,
    IntegratorSpec,
    WrapperSpec,
    ObjectiveTermSpec,
    SolverSpec,
    WarmStartSpec,
    RolloutSpec,
    RolloutReportSpec,
    RefereeSpec,
    FreeDt,
    Fixed,
    Free

"""
    FreeDt

Poka-yoke sum type for the timestep degree of freedom. Either [`Fixed`](@ref)
(Δt is pinned) or [`Free`](@ref)`(lo, hi)` (Δt is an optimization variable within
`[lo, hi]`). A bare `true` is deliberately unrepresentable — the parser rejects
`free_dt = true`, forcing an explicit `[lo, hi]` window.
"""
abstract type FreeDt end

"""
    Fixed() <: FreeDt

The timestep is fixed (not an optimization variable).
"""
struct Fixed <: FreeDt end

"""
    Free(lo, hi) <: FreeDt

The timestep is a free optimization variable bounded to `[lo, hi]`. Requires
`0 < lo < hi`.
"""
struct Free <: FreeDt
    lo::Float64
    hi::Float64
    function Free(lo, hi)
        (0 < lo < hi) ||
            throw(ArgumentError("Free(lo,hi) requires 0 < lo < hi; got ($lo,$hi)"))
        new(Float64(lo), Float64(hi))
    end
end

"""
    SystemSpec

Declarative description of the quantum system. `kind` is one of `:template`
(look up a registered system constructor), `:raw` (build a `QuantumSystem` from
explicit `H_drift`/`H_drives`), or `:composite` (reserved — Phase 1 supports
`:template` and `:raw` only).
"""
Base.@kwdef struct SystemSpec
    kind::Symbol = :template            # template | composite | raw
    template::Union{Nothing,Symbol} = nothing
    params::Dict{Symbol,Any} = Dict{Symbol,Any}()
    global_params::Dict{Symbol,Any} = Dict{Symbol,Any}()
    # composite/raw fields reserved (Phase-1 supports template + raw; composite reserved)
    components::Union{Nothing,Vector{Any}} = nothing
    H_drift::Any = nothing
    H_drives::Any = nothing
end

"""
    GoalSpec

Declarative description of the optimization target. `kind` is `:unitary` or
`:ket`. A `unitary` goal names a `gate` (looked up in `GATES`) or supplies an
inline `matrix` (nested `[re, im]` pairs); a `ket` goal supplies `target` (and
optionally `initial`) as ket strings.
"""
Base.@kwdef struct GoalSpec
    kind::Symbol                        # unitary | ket
    gate::Union{Nothing,Symbol} = nothing
    matrix::Any = nothing               # nested [re,im] pairs when no gate name
    target::Union{Nothing,String} = nothing   # ket string
    initial::Union{Nothing,String} = nothing
    subsystem_levels::Union{Nothing,Vector{Int}} = nothing
    subspace::Union{Nothing,Vector{Vector{Int}}} = nothing
end

"""
    PulseSpec

Declarative description of the control pulse ansatz. `kind` is one of
`:zero_order`, `:linear_spline`, or `:cubic_spline`; `T` is the scalar duration;
`init` chooses the initial guess (`:default`, `:zeros`, `:random`) seeded by
`seed`.
"""
Base.@kwdef struct PulseSpec
    kind::Symbol                        # zero_order | linear_spline | cubic_spline
    T::Float64
    init::Symbol = :default             # default | zeros | random
    seed::Int = 0
end

"""
    TrajectorySpec

Optional trajectory descriptor. In schema v1 the trajectory kind is derived from
`goal.kind`; this struct reserves an explicit override.
"""
Base.@kwdef struct TrajectorySpec
    kind::Union{Nothing,Symbol} = nothing   # optional in v1; derived from goal.kind
end

"""
    ObjectiveTermSpec

A single extra objective term with a `kind` (e.g. `:time`, `:reg_u`, `:reg_du`,
`:reg_ddu`, `:hermite_c2`, `:hermite_bending_energy`, `:sensitivity`, `:leakage`)
and a scalar `weight`.
"""
Base.@kwdef struct ObjectiveTermSpec
    kind::Symbol                        # time | reg_u | reg_du | reg_ddu | hermite_c2 | hermite_bending_energy | sensitivity | leakage
    weight::Float64 = 1.0
end

"""
    TemplateBlock

The `[problem]` block: which registered problem `template` to call, the knot
count `N`, and the composition axes (`goal_treatment`, `free_dt`, regularizer
weights, bounds, free-phase, globals, extra `objectives`, and a passthrough
`options` map for a `PiccoloOptions` subset). `global_names` defaults to empty —
globals are pinned unless explicitly listed.
"""
Base.@kwdef struct TemplateBlock
    template::Symbol                    # registry name
    N::Int
    goal_treatment::Symbol = :objective # objective | constraint | both
    final_fidelity::Union{Nothing,Float64} = nothing
    free_dt::FreeDt = Fixed()
    Q::Float64 = 100.0
    R::Float64 = 1e-2
    R_u::Union{Nothing,Float64} = nothing
    R_du::Union{Nothing,Float64} = nothing
    R_ddu::Union{Nothing,Float64} = nothing
    du_bound::Float64 = Inf
    ddu_bound::Union{Nothing,Float64} = nothing
    free_phase::Bool = false
    initial_phases::Union{Nothing,Vector{Float64}} = nothing
    calibration_targets::Vector{Symbol} = Symbol[]
    global_names::Vector{Symbol} = Symbol[]        # house default: pinned
    global_bounds::Dict{Symbol,Any} = Dict{Symbol,Any}()
    objectives::Vector{ObjectiveTermSpec} = ObjectiveTermSpec[]
    options::Dict{Symbol,Any} = Dict{Symbol,Any}()  # → PiccoloOptions subset
end

"""
    IntegratorSpec

The `[integrator]` block: `kind` is `:bilinear` (Piccolo default — a sentinel
meaning `integrator=nothing`), `:exponential`, or `:spline` (Piccolissimo);
`alg` selects the propagation algorithm for the ODE-backed integrators.
"""
Base.@kwdef struct IntegratorSpec
    kind::Symbol = :bilinear            # bilinear | exponential | spline
    alg::Symbol = :tsit5                # tsit5 | magnus_gl4 | magnus_adapt4
end

"""
    WrapperSpec

A problem wrapper (Phase 1: `:sampling`; `:robust` deferred). `variants` carries
per-variant `[system].params` overrides; `weights` optionally weights them.
"""
Base.@kwdef struct WrapperSpec
    kind::Symbol                        # sampling | robust(deferred)
    variants::Vector{Dict{Symbol,Any}} = Dict{Symbol,Any}[]
    weights::Union{Nothing,Vector{Float64}} = nothing
end

"""
    SolverSpec

The `[solver]` block. Phase 1 executes `:ipopt` only; `:altissimo` (with
`device`/`precision`) is registered for schema but its backend dispatch is
deferred.
"""
Base.@kwdef struct SolverSpec
    backend::Symbol = :ipopt            # ipopt | altissimo
    device::Symbol = :cpu               # cpu | gpu
    precision::Symbol = :f64
    max_iter::Int = 500
    tol::Union{Nothing,Float64} = nothing
    strategy::Symbol = :direct          # direct | continuation | staged
end

"""
    WarmStartSpec

Optional warm-start reference (`catalog_ref` or `pulse_hash`). Resolution is
stubbed in Phase 1.
"""
Base.@kwdef struct WarmStartSpec
    catalog_ref::Union{Nothing,String} = nothing
    pulse_hash::Union{Nothing,String} = nothing
end

"""
    RefereeSpec

The `[referee]` block recording a solve's `{run, solve_knots, solve_integrator,
fidelity_reported}` so an independent referee rollout can be verified against it.
"""
Base.@kwdef struct RefereeSpec
    run::String
    solve_knots::Int
    solve_integrator::Symbol
    fidelity_reported::Float64
end

"""
    RolloutReportSpec

Which quantities a rollout should report (`fidelity`, `populations`, `leakage`),
plus an optional override `goal`.
"""
Base.@kwdef struct RolloutReportSpec
    fidelity::Bool = true
    goal::Union{Nothing,GoalSpec} = nothing
    populations::Bool = false
    leakage::Bool = false
end

"""
    ProblemSpec

Top-level declarative spec for the `control` kind: a versioned bundle of
`system`, `goal`, `pulse`, `trajectory`, `problem`, `integrator`, `wrappers`,
`solver`, and optional `warm_start`. `schema_version` is an integer.
"""
Base.@kwdef struct ProblemSpec
    schema_version::Int = 1
    kind::Symbol = :control
    system::SystemSpec
    goal::Union{Nothing,GoalSpec} = nothing
    pulse::Union{Nothing,PulseSpec} = nothing
    trajectory::TrajectorySpec = TrajectorySpec()
    problem::Union{Nothing,TemplateBlock} = nothing
    integrator::Union{Nothing,IntegratorSpec} = nothing
    wrappers::Vector{WrapperSpec} = WrapperSpec[]
    solver::SolverSpec = SolverSpec()
    warm_start::Union{Nothing,WarmStartSpec} = nothing
end

"""
    RolloutSpec

Top-level declarative spec for the `rollout` kind: an `input_pulse` reference, a
`system`, the `rollout_kind` (`:unitary`/`:ket`/`:density`), the propagation
`alg`, sampling controls, a `report`, and an optional `referee`.
"""
Base.@kwdef struct RolloutSpec
    schema_version::Int = 1
    input_pulse::String                 # catalog:… | run_dir | jld2 path
    system::SystemSpec
    rollout_kind::Symbol                # unitary | ket | density
    alg::Symbol = :tsit5
    n_samples::Union{Nothing,Int} = nothing
    times::Union{Nothing,Vector{Float64}} = nothing
    initial::Union{Nothing,String} = nothing
    report::RolloutReportSpec = RolloutReportSpec()
    referee::Union{Nothing,RefereeSpec} = nothing
end

end

@testitem "spec structs construct and FreeDt is a sum type" begin
    using Piccolo.Specs
    fd = Specs.Free(0.5, 2.0)
    @test fd isa Specs.FreeDt
    @test Specs.Fixed() isa Specs.FreeDt
    g = Specs.GoalSpec(kind = :unitary, gate = :CZ, subsystem_levels = [3, 3])
    @test g.kind == :unitary
    s = Specs.ProblemSpec(
        schema_version = 1,
        kind = :control,
        system = Specs.SystemSpec(kind = :template, template = :TransmonSystem),
        goal = g,
    )
    @test s.kind == :control
end
