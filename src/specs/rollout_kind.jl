# ===========================================================================
# Rollout kind — materialize + run (Task 10) and referee independence + verdict
# (Task 11).
#
# A `rollout` spec re-plays a *saved pulse* through a (possibly provenance-shifted)
# system and reports fidelity — the read-only twin of the `control` kind. When it
# carries a `[referee]` block recording the solve that produced the pulse, the
# runner re-rolls at a strictly-finer resolution with a *different* integrator
# family and mints an `Agree`/`Disagree` verdict (a poka-yoke against a rollout
# that silently reuses the solve's own discretization).
# ===========================================================================

using OrdinaryDiffEqTsit5: Tsit5
using OrdinaryDiffEqLinear: MagnusGL4, MagnusAdapt4

export run_spec, RolloutResult

# ---------------------------------------------------------------------------
# Result type
# ---------------------------------------------------------------------------

"""
    RolloutResult

The outcome of `run_spec(::RolloutSpec)`: the reported `fidelity`, an optional
referee `verdict` (`nothing` unless the spec carried a `[referee]` block), and a
`report` dictionary for any extra requested quantities (populations/leakage).
"""
struct RolloutResult
    fidelity::Float64
    verdict::Union{Nothing,Any}   # narrowed to `Union{Nothing,Verdict}` semantically (Task 11)
    report::Dict{Symbol,Any}
end

# ---------------------------------------------------------------------------
# Rollout algorithm mapping. Symbols → concrete ODE algorithm instances. These
# are the same families the forward propagators use; the sensitivity equations
# are unaffected (they always use Tsit5 internally).
# ---------------------------------------------------------------------------

function _rollout_alg(a::Symbol)
    a === :tsit5 && return Tsit5()
    a === :magnus_gl4 && return MagnusGL4()
    a === :magnus_adapt4 && return MagnusAdapt4()
    throw(SpecValidationError([SpecError("alg", "unknown rollout algorithm",
        string(a), ["tsit5", "magnus_gl4", "magnus_adapt4"])]))
end

# ---------------------------------------------------------------------------
# Input pulse resolution. Phase 1 resolves only an existing `.jld2` path;
# catalog:/run_dir resolution is deferred (structured error, never a crash).
# ---------------------------------------------------------------------------

function _load_input_pulse(s::AbstractString)
    (endswith(s, ".jld2") && isfile(s)) && return load_pulse(s)
    throw(SpecValidationError([SpecError("input_pulse",
        "Phase-1 supports only an existing .jld2 pulse path (catalog:/run_dir resolution deferred)",
        s)]))
end

# ---------------------------------------------------------------------------
# Goal for the rollout comes from the [report.goal] block (a `GoalSpec`).
# ---------------------------------------------------------------------------

function _rollout_goal(spec::RolloutSpec, sys)
    spec.report.goal === nothing && throw(SpecValidationError([SpecError("report.goal",
        "rollout fidelity requires a [report.goal] block")]))
    return _build_goal(spec.report.goal, sys)
end

# ---------------------------------------------------------------------------
# Pre-pass validation (structured errors, never a materialization crash).
# ---------------------------------------------------------------------------

function _validate_rollout!(spec::RolloutSpec, errs::Vector{SpecError})
    if !(spec.rollout_kind in (:unitary, :ket, :density))
        push!(errs, SpecError("rollout_kind", "unknown rollout kind",
            string(spec.rollout_kind), ["unitary", "ket", "density"]))
    end
    if spec.system.kind === :template && lookup_system(spec.system.template) === nothing
        push!(errs, SpecError("system.template", "unknown system template",
            string(spec.system.template), sort!(String[string(k) for k in keys(SYSTEMS)])))
    elseif spec.system.kind === :composite
        push!(errs, SpecError("system.kind",
            "composite systems are deferred in Phase 1 (template + raw only)",
            string(spec.system.kind)))
    end
    return errs
end

# ---------------------------------------------------------------------------
# materialize(::RolloutSpec) — load the pulse, build the system + goal, and
# re-roll the trajectory at the requested algorithm/resolution. Returns the
# rolled-out `AbstractQuantumTrajectory` (the read-only twin of the control
# kind's `QuantumControlProblem`).
# ---------------------------------------------------------------------------

"""
    materialize(spec::RolloutSpec) -> AbstractQuantumTrajectory

Materialize a `rollout` spec into a rolled-out quantum trajectory: load the saved
pulse from `input_pulse`, build the (possibly provenance-shifted) system and the
`[report.goal]`, and re-solve the ODE with the requested `alg`/`n_samples`.

`rollout_kind = "density"` requires an `OpenQuantumSystem`; Phase-1 template
systems are closed, so it returns a structured [`SpecValidationError`](@ref).
"""
function materialize(spec::RolloutSpec)
    errs = SpecError[]
    _validate_rollout!(spec, errs)
    isempty(errs) || throw(SpecValidationError(errs))

    pulse = _load_input_pulse(spec.input_pulse)
    sys = _build_system(spec.system)
    n_save = spec.n_samples === nothing ? 101 : spec.n_samples
    alg = _rollout_alg(spec.alg)
    rk = spec.rollout_kind

    if rk === :unitary
        goal = _rollout_goal(spec, sys)
        qtraj = UnitaryTrajectory(sys, pulse, goal)
        return rollout(qtraj; algorithm=alg, n_save=n_save)
    elseif rk === :ket
        ψg = _rollout_goal(spec, sys)
        ψ0 = spec.initial === nothing ?
             throw(SpecValidationError([SpecError("initial",
                 "ket rollout requires an `initial` ket")])) :
             _parse_ket(spec.initial)
        qtraj = KetTrajectory(sys, pulse, ψ0, ψg)
        return rollout(qtraj; algorithm=alg, n_save=n_save)
    else # :density — trait check: requires an OpenQuantumSystem
        sys isa OpenQuantumSystem || throw(SpecValidationError([SpecError("rollout_kind",
            "density rollout requires an OpenQuantumSystem; Phase-1 template systems are closed",
            "density")]))
        throw(SpecValidationError([SpecError("rollout_kind",
            "density rollout is not yet wired in Phase 1", "density")]))
    end
end

# ---------------------------------------------------------------------------
# run_spec(::RolloutSpec) — solve the rollout and compute the requested report.
# The referee branch (verdict) is added in Task 11.
# ---------------------------------------------------------------------------

"""
    run_spec(spec::RolloutSpec; kwargs...) -> RolloutResult

Materialize and evaluate a `rollout` spec, returning a [`RolloutResult`](@ref).
When the spec carries a `[referee]` block, the rollout is re-verified for
independence and a `verdict` is minted (Task 11); otherwise `verdict === nothing`.
"""
function run_spec(spec::RolloutSpec; kwargs...)
    spec.referee === nothing || _verify_referee!(spec)
    qtraj = materialize(spec)
    f = spec.report.fidelity ? Float64(fidelity(qtraj)) : NaN
    verdict = spec.referee === nothing ? nothing :
              _verdict(f, spec.referee.fidelity_reported)
    return RolloutResult(f, verdict, Dict{Symbol,Any}())
end

@testitem "rollout: unitary fidelity matches a direct rollout" begin
    using Piccolo, Piccolo.Specs
    import JLD2
    using OrdinaryDiffEqTsit5: Tsit5

    # Tiny control solve (X on a single 3-level transmon) → save the pulse.
    TINY_CZ_TOML = """
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
    qcp = Specs.materialize(Specs.parse_spec(TINY_CZ_TOML; format=:toml))
    solve!(qcp; max_iter=15, print_level=0, verbose=false)
    pulse = extract_pulse(qcp.qtraj, get_trajectory(qcp))
    path = tempname() * ".jld2"
    JLD2.save(path, pulse)

    ROLLOUT_TOML = """
    schema_version = 1
    kind = "rollout"
    input_pulse = "$(path)"
    rollout_kind = "unitary"
    alg = "tsit5"
    n_samples = 51
    [system]
    kind = "template"
    template = "TransmonSystem"
    params = { levels = 3, drive_bounds = [0.02, 0.02] }
    [report]
    fidelity = true
    [report.goal]
    kind = "unitary"
    gate = "X"
    """
    rspec = Specs.parse_spec(ROLLOUT_TOML; format=:toml)
    res = Specs.run_spec(rspec)

    # Direct rollout with the SAME system/pulse/alg/resolution.
    sys = TransmonSystem(; levels=3, drive_bounds=[0.02, 0.02])
    goal = EmbeddedOperator(GATES[:X], sys)
    qtraj = UnitaryTrajectory(sys, pulse, goal)
    direct = fidelity(rollout(qtraj; algorithm=Tsit5(), n_save=51))

    @test res.verdict === nothing            # no [referee] ⇒ no verdict
    @test isapprox(res.fidelity, direct; atol=1e-8)
end
