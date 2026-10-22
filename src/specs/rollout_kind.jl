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
    throw(
        SpecValidationError([
            SpecError(
                "alg",
                "unknown rollout algorithm",
                string(a),
                ["tsit5", "magnus_gl4", "magnus_adapt4"],
            ),
        ]),
    )
end

# ---------------------------------------------------------------------------
# Input pulse resolution. Phase 1 resolves only an existing `.jld2` path;
# catalog:/run_dir resolution is deferred (structured error, never a crash).
# ---------------------------------------------------------------------------

function _load_input_pulse(s::AbstractString)
    (endswith(s, ".jld2") && isfile(s)) && return load_pulse(s)
    throw(
        SpecValidationError([
            SpecError(
                "input_pulse",
                "Phase-1 supports only an existing .jld2 pulse path (catalog:/run_dir resolution deferred)",
                s,
            ),
        ]),
    )
end

# ---------------------------------------------------------------------------
# Goal for the rollout comes from the [report.goal] block (a `GoalSpec`).
# ---------------------------------------------------------------------------

function _rollout_goal(spec::RolloutSpec, sys)
    spec.report.goal === nothing && throw(
        SpecValidationError([
            SpecError("report.goal", "rollout fidelity requires a [report.goal] block"),
        ]),
    )
    return _build_goal(spec.report.goal, sys)
end

# ---------------------------------------------------------------------------
# Pre-pass validation (structured errors, never a materialization crash).
# ---------------------------------------------------------------------------

function _validate_rollout!(spec::RolloutSpec, errs::Vector{SpecError})
    if !(spec.rollout_kind in (:unitary, :ket, :density))
        push!(
            errs,
            SpecError(
                "rollout_kind",
                "unknown rollout kind",
                string(spec.rollout_kind),
                ["unitary", "ket", "density"],
            ),
        )
    end
    # Mirrors the control path in `materialize.jl`: `template` is
    # Union{Nothing,Symbol} and `lookup_system` takes only a Symbol, so an omitted
    # `template` has to be caught ahead of the lookup rather than by it.
    if spec.system.kind === :template && spec.system.template === nothing
        push!(
            errs,
            SpecError(
                "system.template",
                "missing `template` for a system with `kind = \"template\"`",
                nothing,
                sort!(String[string(k) for k in keys(SYSTEMS)]),
            ),
        )
    elseif spec.system.kind === :template && lookup_system(spec.system.template) === nothing
        push!(
            errs,
            SpecError(
                "system.template",
                "unknown system template",
                string(spec.system.template),
                sort!(String[string(k) for k in keys(SYSTEMS)]),
            ),
        )
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
        return rollout(qtraj; algorithm = alg, n_save = n_save)
    elseif rk === :ket
        ψg = _rollout_goal(spec, sys)
        ψ0 =
            spec.initial === nothing ?
            throw(
                SpecValidationError([
                    SpecError("initial", "ket rollout requires an `initial` ket"),
                ]),
            ) : _parse_ket(spec.initial)
        qtraj = KetTrajectory(sys, pulse, ψ0, ψg)
        return rollout(qtraj; algorithm = alg, n_save = n_save)
    else # :density — trait check: requires an OpenQuantumSystem
        sys isa OpenQuantumSystem || throw(
            SpecValidationError([
                SpecError(
                    "rollout_kind",
                    "density rollout requires an OpenQuantumSystem; Phase-1 template systems are closed",
                    "density",
                ),
            ]),
        )
        throw(
            SpecValidationError([
                SpecError(
                    "rollout_kind",
                    "density rollout is not yet wired in Phase 1",
                    "density",
                ),
            ]),
        )
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
    verdict =
        spec.referee === nothing ? nothing : _verdict(f, spec.referee.fidelity_reported)
    return RolloutResult(f, verdict, Dict{Symbol,Any}())
end

# ===========================================================================
# Referee independence + verdict (Task 11)
#
# A referee rollout re-checks a solve's reported fidelity from *independent*
# discretization choices: strictly finer sampling AND a different integrator
# family. The typed `Agree`/`Disagree` verdict (both carrying the re-rolled and
# reported fidelities) is the Phase-1 subset of the merged spec's `Verdict`
# hierarchy (transport's richer set is out of scope). Forging a referee block
# whose axes are NOT independent throws rather than minting a verdict.
# ===========================================================================

export Verdict, Agree, Disagree, referee_rollout

"""
    Verdict

Abstract supertype for a referee outcome. Phase-1 concretes are [`Agree`](@ref)
and [`Disagree`](@ref); both carry `f_rerolled` (the referee's independent
re-rollout fidelity) and `f_reported` (the solve's reported fidelity).
"""
abstract type Verdict end

"""
    Agree(f_rerolled, f_reported) <: Verdict

The referee's independent re-rollout agrees with the solve's reported fidelity
(within tolerance).
"""
struct Agree <: Verdict
    f_rerolled::Float64
    f_reported::Float64
end

"""
    Disagree(f_rerolled, f_reported) <: Verdict

The referee's independent re-rollout disagrees with the solve's reported
fidelity — a flag that the solve's discretization was not converged.
"""
struct Disagree <: Verdict
    f_rerolled::Float64
    f_reported::Float64
end

# Integrator-family coarse classes. Referee independence requires the referee's
# rollout algorithm to live in a *different* class than the solve integrator.
# Collocation covers the solve-side integrators (bilinear/spline/exponential);
# tsit5 is the RK family; magnus_* the Magnus family.
function _integrator_family(s::Symbol)
    s === :tsit5 && return :rk
    (s === :magnus_gl4 || s === :magnus_adapt4) && return :magnus
    return :collocation
end

# Pick a referee rollout algorithm in a different family than the solve
# integrator: a Magnus solve/collocation → Tsit5 (RK); an RK solve → Magnus.
_referee_alg(solve_integrator::Symbol) =
    _integrator_family(solve_integrator) === :rk ? :magnus_adapt4 : :tsit5

_verdict(f_rerolled::Real, f_reported::Real; tol::Float64 = 1e-3) =
    abs(f_rerolled - f_reported) <= tol ? Agree(Float64(f_rerolled), Float64(f_reported)) :
    Disagree(Float64(f_rerolled), Float64(f_reported))

# Re-verify a referee block at run time: the rollout's own axes must be strictly
# finer than, and in a different integrator family than, the recorded solve.
# A forged/inconsistent block throws `SpecValidationError` (never mints a verdict).
function _verify_referee!(spec::RolloutSpec)
    ref = spec.referee
    n = spec.n_samples === nothing ? 101 : spec.n_samples
    errs = SpecError[]
    n > ref.solve_knots || push!(
        errs,
        SpecError(
            "referee",
            "referee resolution (n_samples=$n) must be strictly finer than the solve " *
            "($(ref.solve_knots) knots)",
            n,
        ),
    )
    _integrator_family(spec.alg) != _integrator_family(ref.solve_integrator) || push!(
        errs,
        SpecError(
            "referee",
            "referee integrator family (:$(spec.alg)) must differ from the solve " *
            "family (:$(ref.solve_integrator))",
            string(spec.alg),
        ),
    )
    isempty(errs) || throw(SpecValidationError(errs))
    return nothing
end

# Locate the single saved pulse in a completed run directory.
function _run_dir_pulse(run_dir::AbstractString)
    jld2s = filter(f -> endswith(f, ".jld2"), readdir(run_dir))
    isempty(jld2s) && throw(
        SpecValidationError([
            SpecError(
                "referee.run",
                "no .jld2 pulse found in run directory",
                String(run_dir),
            ),
        ]),
    )
    return joinpath(run_dir, first(jld2s))
end

"""
    referee_rollout(control_spec::ProblemSpec, run_dir::AbstractString) -> RolloutSpec

Construct an *independent* referee rollout for a completed control run. Reads the
run's saved pulse and `result.toml` (for the reported fidelity), then chooses
referee axes that are strictly finer (`n_samples > solve_knots`) and a different
integrator family (e.g. Tsit5 for a bilinear solve). The returned spec carries a
`[referee]` block recording the solve's `{solve_knots, solve_integrator,
fidelity_reported}`; `run_spec` re-verifies these at run time before minting a
verdict.
"""
function referee_rollout(control_spec::ProblemSpec, run_dir::AbstractString)
    control_spec.problem === nothing && throw(
        SpecValidationError([
            SpecError(
                "problem",
                "referee_rollout requires a control spec with a [problem] block",
            ),
        ]),
    )
    control_spec.goal === nothing && throw(
        SpecValidationError([
            SpecError(
                "goal",
                "referee_rollout requires a control spec with a [goal] block",
            ),
        ]),
    )

    result = TOML.parsefile(joinpath(run_dir, "result.toml"))
    fidelity_reported = Float64(result["fidelity"])
    pulse_path = _run_dir_pulse(run_dir)

    solve_knots = control_spec.problem.N
    solve_integrator =
        control_spec.integrator === nothing ? :bilinear : control_spec.integrator.kind
    n_samples = max(4 * solve_knots, solve_knots + 1)   # strictly finer
    alg = _referee_alg(solve_integrator)

    referee = RefereeSpec(;
        run = String(run_dir),
        solve_knots = solve_knots,
        solve_integrator = solve_integrator,
        fidelity_reported = fidelity_reported,
    )
    report = RolloutReportSpec(; fidelity = true, goal = control_spec.goal)

    return RolloutSpec(;
        schema_version = control_spec.schema_version,
        input_pulse = pulse_path,
        system = control_spec.system,
        rollout_kind = control_spec.goal.kind,
        alg = alg,
        n_samples = n_samples,
        report = report,
        referee = referee,
    )
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
    qcp = Specs.materialize(Specs.parse_spec(TINY_CZ_TOML; format = :toml))
    solve!(qcp; max_iter = 15, print_level = 0, verbose = false)
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
    rspec = Specs.parse_spec(ROLLOUT_TOML; format = :toml)
    res = Specs.run_spec(rspec)

    # Direct rollout with the SAME system/pulse/alg/resolution.
    sys = TransmonSystem(; levels = 3, drive_bounds = [0.02, 0.02])
    goal = EmbeddedOperator(GATES[:X], sys)
    qtraj = UnitaryTrajectory(sys, pulse, goal)
    direct = fidelity(rollout(qtraj; algorithm = Tsit5(), n_save = 51))

    @test res.verdict === nothing            # no [referee] ⇒ no verdict
    @test isapprox(res.fidelity, direct; atol = 1e-8)
end

@testitem "referee: independence enforced, verdict on referee-verified only" begin
    using Piccolo, Piccolo.Specs
    import JLD2, TOML

    CTRL_TOML = """
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
    ctrl = Specs.parse_spec(CTRL_TOML; format = :toml)

    # Build a completed run dir by hand (mimics solve_spec's artifacts): the saved
    # pulse + a result.toml carrying the reported fidelity.
    run_dir = mktempdir()
    qcp = Specs.materialize(ctrl)
    solve!(qcp; max_iter = 15, print_level = 0, verbose = false)
    fid = Float64(fidelity(qcp))
    pulse = extract_pulse(qcp.qtraj, get_trajectory(qcp))
    JLD2.save(joinpath(run_dir, "pulse-abc123.jld2"), pulse)
    open(joinpath(run_dir, "result.toml"), "w") do io
        TOML.print(
            io,
            Dict(
                "schema_version" => "1",
                "fidelity" => fid,
                "iterations" => 15,
                "wall_seconds" => 0.0,
            ),
        )
    end

    # valid referee: strictly finer resolution + a different integrator family
    # (Tsit5 for a bilinear solve).
    rr = Specs.referee_rollout(ctrl, run_dir)
    @test rr.referee !== nothing
    @test rr.n_samples > rr.referee.solve_knots
    @test rr.referee.solve_integrator === :bilinear
    res = Specs.run_spec(rr)
    @test res.verdict isa Specs.Verdict
    @test res.verdict isa Specs.Agree        # same pulse rerolls to reported fidelity

    # forged referee (resolution == solve, not strictly finer) fails run-time
    # re-verification → SpecValidationError (never mints a verdict).
    pulse_path = joinpath(run_dir, "pulse-abc123.jld2")
    FORGED_TOML = """
    schema_version = 1
    kind = "rollout"
    input_pulse = "$(pulse_path)"
    rollout_kind = "unitary"
    alg = "tsit5"
    n_samples = 11
    [system]
    kind = "template"
    template = "TransmonSystem"
    params = { levels = 3, drive_bounds = [0.02, 0.02] }
    [report]
    fidelity = true
    [report.goal]
    kind = "unitary"
    gate = "X"
    [referee]
    run = "$(run_dir)"
    solve_knots = 11
    solve_integrator = "bilinear"
    fidelity_reported = $(fid)
    """
    forged = Specs.parse_spec(FORGED_TOML; format = :toml)
    @test_throws Specs.SpecValidationError Specs.run_spec(forged)

    # rollout without a [referee] block → no verdict.
    PLAIN_TOML = """
    schema_version = 1
    kind = "rollout"
    input_pulse = "$(pulse_path)"
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
    plain = Specs.parse_spec(PLAIN_TOML; format = :toml)
    @test Specs.run_spec(plain).verdict === nothing
end

@testitem "rollout: omitted system.template is a structured error, not a crash" begin
    using Piccolo, Piccolo.Specs

    # Same guard as the control path in materialize.jl — `_validate_rollout!` called
    # `lookup_system(spec.system.template)` with template::Union{Nothing,Symbol}, so
    # an omitted key threw MethodError instead of reporting the field.
    NO_TEMPLATE = """
    schema_version = 1
    kind = "rollout"
    rollout_kind = "unitary"
    input_pulse = "pulse.jld2"
    [system]
    kind = "template"
    """
    err = try
        Specs.materialize(Specs.parse_spec(NO_TEMPLATE; format = :toml))
        nothing
    catch e
        e
    end
    @test err isa Specs.SpecValidationError
    st = only(filter(e -> e.path == "system.template", err.errors))
    @test occursin("missing", st.msg)
    @test "TransmonSystem" in st.allowed
end

@testitem "coverage: rollout alg variants, error paths, and the ket rollout lane" begin
    using Piccolo, Piccolo.Specs
    import JLD2
    using OrdinaryDiffEqTsit5: Tsit5

    # A tiny saved pulse for the happy-path lanes.
    CTRL_TOML = """
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
    qcp = Specs.materialize(Specs.parse_spec(CTRL_TOML; format = :toml))
    solve!(qcp; max_iter = 10, print_level = 0, verbose = false)
    pulse = extract_pulse(qcp.qtraj, get_trajectory(qcp))
    path = tempname() * ".jld2"
    JLD2.save(path, pulse)

    # ── alg variants: magnus_gl4 / magnus_adapt4 map; unknown alg errors ──
    for (alg, T_alg) in (:magnus_gl4 => "magnus_gl4", :magnus_adapt4 => "magnus_adapt4")
        ROLLOUT_TOML = """
        schema_version = 1
        kind = "rollout"
        input_pulse = "$(path)"
        rollout_kind = "unitary"
        alg = "$(T_alg)"
        n_samples = 21
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
        res = Specs.run_spec(Specs.parse_spec(ROLLOUT_TOML; format = :toml))
        @test isfinite(res.fidelity)
    end

    # unknown alg: the structured error names the allowed set
    BAD_ALG = """
    schema_version = 1
    kind = "rollout"
    input_pulse = "$(path)"
    rollout_kind = "unitary"
    alg = "bogus"
    [system]
    kind = "template"
    template = "TransmonSystem"
    params = { levels = 3 }
    [report]
    fidelity = true
    [report.goal]
    kind = "unitary"
    gate = "X"
    """
    err = try
        Specs.parse_spec(BAD_ALG; format = :toml) |> Specs.run_spec
        nothing
    catch e
        e
    end
    @test err isa Specs.SpecValidationError

    # ── the ket lane: goal + initial ket strings ──
    KET_TOML = """
    schema_version = 1
    kind = "rollout"
    input_pulse = "$(path)"
    rollout_kind = "ket"
    alg = "tsit5"
    initial = "ComplexF64[1.0, 0.0, 0.0]"
    [system]
    kind = "template"
    template = "TransmonSystem"
    params = { levels = 3, drive_bounds = [0.02, 0.02] }
    [report]
    fidelity = true
    [report.goal]
    kind = "ket"
    target = "ComplexF64[0.0, 1.0, 0.0]"
    """
    kres = Specs.run_spec(Specs.parse_spec(KET_TOML; format = :toml))
    @test isfinite(kres.fidelity)

    # ket without `initial`: the structured error
    KET_NOINIT = replace(KET_TOML, "initial = \"ComplexF64[1.0, 0.0, 0.0]\"" => "")
    err2 = try
        Specs.parse_spec(KET_NOINIT; format = :toml) |> Specs.run_spec
        nothing
    catch e
        e
    end
    @test err2 isa Specs.SpecValidationError

    # ── density lanes: not-yet-wired (closed system) + unknown kind ──
    DENS_TOML = """
    schema_version = 1
    kind = "rollout"
    input_pulse = "$(path)"
    rollout_kind = "density"
    alg = "tsit5"
    [system]
    kind = "template"
    template = "TransmonSystem"
    params = { levels = 3 }
    [report]
    fidelity = true
    [report.goal]
    kind = "unitary"
    gate = "X"
    """
    err3 = try
        Specs.parse_spec(DENS_TOML; format = :toml) |> Specs.run_spec
        nothing
    catch e
        e
    end
    @test err3 isa Specs.SpecValidationError

    BAD_KIND =
        replace(DENS_TOML, "rollout_kind = \"density\"" => "rollout_kind = \"hologram\"")
    err4 = try
        Specs.parse_spec(BAD_KIND; format = :toml) |> Specs.run_spec
        nothing
    catch e
        e
    end
    @test err4 isa Specs.SpecValidationError

    # non-.jld2 input_pulse: the Phase-1 structured error
    BAD_PULSE = replace(
        DENS_TOML,
        "rollout_kind = \"density\"" => "rollout_kind = \"unitary\"",
        "input_pulse = \"$(path)\"" => "input_pulse = \"catalog:xyz\"",
    )
    err5 = try
        Specs.parse_spec(BAD_PULSE; format = :toml) |> Specs.run_spec
        nothing
    catch e
        e
    end
    @test err5 isa Specs.SpecValidationError
end
