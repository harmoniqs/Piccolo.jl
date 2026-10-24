# ===========================================================================
# solve_spec / run_spec — the single vetted generic runner (Task 12)
#
# One code path takes a declarative `control` spec end-to-end: parse → validate →
# materialize → solve (the declared backend) → save the pulse (never the
# trajectory) → write the amicode `result/v1` artifact (`result.toml`).
# Per-iteration progress is streamed as `AMICODE_ITER k=v …` lines and a terminal
# `DONE fidelity=…`, the exact contract amico-run's `telemetry.classifyLine`
# parses (`AMICODE_ITER` prefix + whitespace-split k=v tokens; `^DONE(\s|$)`).
# No per-problem scripts.
# ===========================================================================

export solve_spec

import JLD2
import MadNLP
import Pkg
using DirectTrajOpt: IpoptOptions, MadNLPOptions

# Stack packages whose versions the learning-loops ledger (Plan 3) reads from
# `params.versions`. Restricted to the optimal-control stack.
const _VERSION_PKGS =
    ("Piccolo", "Piccolissimo", "DirectTrajOpt", "NamedTrajectories", "Altissimo")

# Active `Pkg` versions of the stack packages, for `result.toml`'s
# `params.versions`. The active project itself (e.g. Piccolo when running its own
# tests) is NOT listed among `Pkg.dependencies()`, so it is read separately from
# `Pkg.project()` — this is what guarantees `params.versions.Piccolo` is present.
function _stack_versions()
    versions = Dict{String,Any}()
    proj = Pkg.project()
    if proj.name !== nothing && proj.name in _VERSION_PKGS && proj.version !== nothing
        versions[proj.name] = string(proj.version)
    end
    for (_, dep) in Pkg.dependencies()
        if dep.name in _VERSION_PKGS && dep.version !== nothing
            versions[dep.name] = string(dep.version)
        end
    end
    return versions
end

# Build a PiccoloOptions from the spec's `[problem.options]` block. `display` is
# forced to `:silent` — the templates print a full problem tree at `:standard`
# which would pollute the `AMICODE_ITER`/`DONE` telemetry stream. `materialize`
# itself does not map `options`, so the runner owns this mapping.
function _piccolo_options(spec::ProblemSpec)
    o = spec.problem === nothing ? Dict{Symbol,Any}() : spec.problem.options
    kw = Dict{Symbol,Any}(:display => :silent)
    haskey(o, :leakage_constraint) &&
        (kw[:leakage_constraint] = Bool(o[:leakage_constraint]))
    haskey(o, :leakage_constraint_value) &&
        (kw[:leakage_constraint_value] = Float64(o[:leakage_constraint_value]))
    haskey(o, :leakage_cost) && (kw[:leakage_cost] = Float64(o[:leakage_cost]))
    haskey(o, :timesteps_all_equal) &&
        (kw[:timesteps_all_equal] = Bool(o[:timesteps_all_equal]))
    haskey(o, :bound_state) && (kw[:bound_state] = Bool(o[:bound_state]))
    haskey(o, :bound_state_l2) && (kw[:bound_state_l2] = Bool(o[:bound_state_l2]))
    return PiccoloOptions(; kw...)
end

# AMICODE_ITER telemetry, Ipopt arm (`solver.backend = "ipopt"`). DirectTrajOpt's
# raw `callback` (forwarded via `solve!(qcp; callback=…)`) is an Ipopt
# CallbackFunction *factory*: `callback(optimizer) -> (optimizer_state... -> Bool)`,
# where `optimizer_state` is the 11-tuple `(alg_mod, iter_count, obj_value,
# inf_pr, inf_du, …)`. We emit one `AMICODE_ITER k=v …` line per *main-loop*
# iteration (`alg_mod == 0`, i.e. skipping Ipopt's restoration phase) — the shape
# `telemetry.classifyLine` parses. The factory writes the final iteration into
# `iter_ref` so the runner can read it.
function _amicode_iter_callback(iter_ref::Base.RefValue{Int})
    factory = function (optimizer)
        return function (optimizer_state...)
            optimizer_state[1] == 0 || return true      # main IPM loop only
            it = Int(optimizer_state[2])
            obj = Float64(optimizer_state[3])
            inf_pr = Float64(optimizer_state[4])
            inf_du = Float64(optimizer_state[5])
            iter_ref[] = max(iter_ref[], it)
            @printf(
                "AMICODE_ITER iter=%d obj=%.8e inf_pr=%.6e inf_du=%.6e\n",
                it,
                obj,
                inf_pr,
                inf_du
            )
            flush(stdout)
            return true
        end
    end
    return factory
end

# AMICODE_ITER telemetry, MadNLP arm — the default backend (DirectTrajOpt's
# MadNLP-default flip since DTO 0.11, inherited by Piccolo #360). The MadNLP
# backend does not wire the raw `callback` kwarg; per-iteration emission rides
# `MadNLPOptions.intermediate_callback` as a raw `MadNLP.AbstractUserCallback`.
# The `UserCallbackRegular` filter is the MadNLP analogue of the Ipopt arm's
# `alg_mod == 0` — restoration/robust phases fire without advancing `cnt.k` and
# are skipped, keeping the emitted iters monotone (mirrors DirectTrajOpt's own
# telemetry-parity audit, DTO #155). The agnostic `(primal, iter)` callback form
# cannot carry the obj/inf_pr/inf_du columns the Inspector's stats row reads.
mutable struct _AmicodeIterCallbackMadNLP <: MadNLP.AbstractUserCallback
    iter_ref::Base.RefValue{Int}
end
function (cb::_AmicodeIterCallbackMadNLP)(solver, mode)
    mode isa MadNLP.UserCallbackRegular || return true  # main IPM loop only
    it = Int(MadNLP.get_cnt(solver).k)
    @printf(
        "AMICODE_ITER iter=%d obj=%.8e inf_pr=%.6e inf_du=%.6e\n",
        it,
        Float64(MadNLP.get_obj_val(solver)),
        Float64(MadNLP.get_inf_pr(solver)),
        Float64(MadNLP.get_inf_du(solver)),
    )
    flush(stdout)
    cb.iter_ref[] = max(cb.iter_ref[], it)
    return true
end

"""
    solve_spec(src::AbstractString; run_dir, format=:toml, max_iter=nothing, kwargs...) -> Dict

Run a `control` spec end-to-end: parse `src` (a spec string, or a path to one),
materialize it, solve with the declared solver backend — `:madnlp` (the default,
inherited from DirectTrajOpt's MadNLP-default flip in DTO 0.11) or explicitly
selected `:ipopt` — save the optimized pulse as `pulse-<problem_hash>.jld2`,
and write the amicode `result/v1` artifact `result.toml` into `run_dir`. Emits
`AMICODE_ITER …` progress lines and a terminal `DONE fidelity=…`. Returns the
`result.toml` payload as a `Dict`.

Backend scope: the spec's `solver.backend` declares the executed backend;
both `:madnlp` and `:ipopt` execute. `solver.backend = "altissimo"` is
registered for schema/emit purposes but its backend dispatch is deferred, so
it raises a structured [`SpecValidationError`](@ref) here.
"""
function solve_spec(
    src::AbstractString;
    run_dir::AbstractString,
    format::Symbol = :toml,
    max_iter::Union{Nothing,Int} = nothing,
    kwargs...,
)
    # `src` may be an inline spec string or a path to one. A multi-line string is
    # never a path; guard `isfile` so a spec body is not `stat`-ed as a filename.
    spec_src = (!occursin('\n', src) && isfile(src)) ? read(src, String) : src
    spec = parse_spec(spec_src; format = format)
    spec isa ProblemSpec || throw(
        SpecValidationError([
            SpecError(
                "kind",
                "solve_spec runs `control` specs; use `run_spec` for a rollout spec",
            ),
        ]),
    )
    return _run_control(spec, run_dir; max_iter = max_iter, kwargs...)
end

"""
    run_spec(spec::ProblemSpec; run_dir, max_iter=nothing, kwargs...) -> Dict

Generic runner dispatch for an already-parsed `control` spec: solve it and write
the run artifacts into `run_dir` (see [`solve_spec`](@ref)). The `rollout` kind
is handled by `run_spec(::RolloutSpec)` (Task 10/11).
"""
function run_spec(
    spec::ProblemSpec;
    run_dir::AbstractString,
    max_iter::Union{Nothing,Int} = nothing,
    kwargs...,
)
    return _run_control(spec, run_dir; max_iter = max_iter, kwargs...)
end

function _run_control(
    spec::ProblemSpec,
    run_dir::AbstractString;
    max_iter::Union{Nothing,Int} = nothing,
    kwargs...,
)
    # Backend gate: the spec names its backend, and the runner executes exactly
    # that one (declared == executed — the coherence the Phase-1 Ipopt pin was
    # invented for, kept as the permanent contract now that the pin is retired).
    # `:madnlp` is the inherited DTO-0.11 default (#360); `:ipopt` stays fully
    # selectable; `:altissimo` remains register-for-schema-only (its dispatch is
    # deferred with Piccolissimo's matrix-free arm).
    if spec.solver.backend === :altissimo
        throw(
            SpecValidationError([
                SpecError(
                    "solver.backend",
                    "the altissimo backend is registered for schema only (dispatch deferred); " *
                    "solve_spec executes the madnlp and ipopt backends",
                    "altissimo",
                    ["madnlp", "ipopt"],
                ),
            ]),
        )
    elseif spec.solver.backend !== :madnlp && spec.solver.backend !== :ipopt
        throw(
            SpecValidationError([
                SpecError(
                    "solver.backend",
                    "unknown solver backend",
                    string(spec.solver.backend),
                    ["madnlp", "ipopt", "altissimo"],
                ),
            ]),
        )
    end

    isdir(run_dir) || mkpath(run_dir)

    # Identity keys + the canonical wire form retained in the runner's scope (the
    # Phase-1 `extract_spec` substitute — persisted under params.spec).
    cj = canonical_json(full_dict(spec))
    shash = structure_hash(spec)
    phash = problem_hash(spec)

    # Materialize with display forced :silent so the template's problem tree does
    # not pollute the telemetry stream.
    qcp = materialize(spec; piccolo_options = _piccolo_options(spec))

    # Solve, streaming AMICODE_ITER per main-loop iteration. `options` pins the
    # executed backend to the declared one on both arms — a bare `solve!` would
    # ride DirectTrajOpt's own default instead of the spec's declaration. The
    # arms carry telemetry differently: Ipopt composes the raw `callback`
    # factory into its single CallbackFunction, while the MadNLP backend does
    # not wire `callback`, so its emission rides `intermediate_callback`.
    mi = max_iter === nothing ? spec.solver.max_iter : max_iter
    iter_ref = Ref(0)
    wall = @elapsed begin
        stats = if spec.solver.backend === :madnlp
            solve!(
                qcp;
                options = MadNLPOptions(
                    intermediate_callback = _AmicodeIterCallbackMadNLP(iter_ref),
                ),
                max_iter = mi,
                print_level = 0,
                verbose = false,
                kwargs...,
            )
        else
            solve!(
                qcp;
                options = IpoptOptions(),
                max_iter = mi,
                print_level = 0,
                verbose = false,
                callback = _amicode_iter_callback(iter_ref),
                kwargs...,
            )
        end
    end

    # Declared == executed coherence: the backend symbol the solver reported
    # must be the one the spec declared — a mismatch means the options plumbing
    # above has broken, which is exactly what this gate exists to make loud.
    stats.solver === spec.solver.backend || error(
        "backend coherence violated: spec declared $(spec.solver.backend), " *
        "solver executed $(stats.solver)",
    )

    fid = Float64(fidelity(qcp))
    iters = iter_ref[]
    # `converged` proxy: the backend's SolveStats now surfaces the termination
    # status, but the result/v1 `converged` field keeps its documented heuristic
    # — "stopped before the iteration cap" — so historical result.toml files
    # stay comparable. (The status is in `stats.status` when a caller wants it.)
    converged = iters < mi

    # Save the pulse (never the trajectory) as the single run artifact JLD2.
    pulse = extract_pulse(qcp.qtraj, get_trajectory(qcp))
    pulse_path = joinpath(run_dir, "pulse-$(phash[1:12]).jld2")
    JLD2.save(pulse_path, pulse)

    # Terminal sentinels. `DONE` is parsed by telemetry.classifyLine; the
    # AMICODE_PULSE_META line is a forward-looking pulse-preview log line for the
    # extension's run-dir reader (classifyLine treats it as a plain log line).
    n_drives = get_system(qcp).n_drives
    @printf(
        "AMICODE_PULSE_META n_drives=%d knots=%d pulse=%s\n",
        n_drives,
        spec.problem.N,
        basename(pulse_path)
    )
    @printf("DONE fidelity=%.10f iterations=%d wall_seconds=%.3f\n", fid, iters, wall)
    flush(stdout)

    # result.toml — amicode result/v1. Top-level scalars (schema_version as the
    # string "1", fidelity, iterations, wall_seconds) satisfy the base schema; the
    # hashes, converged flag, canonical spec, and the Pkg versions map go under
    # `[params]` (the only additionalProperties:true table in result/v1).
    result = Dict{String,Any}(
        "schema_version" => "1",
        "fidelity" => fid,
        "iterations" => iters,
        "wall_seconds" => wall,
        "params" => Dict{String,Any}(
            "structure_hash" => shash,
            "problem_hash" => phash,
            "converged" => converged,
            "versions" => _stack_versions(),
            "spec" => cj,
        ),
    )
    open(joinpath(run_dir, "result.toml"), "w") do io
        TOML.print(io, result)
    end
    return result
end

@testitem "solve_spec end-to-end + result.toml + round-trip" begin
    using Piccolo, Piccolo.Specs
    import TOML

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
    dir = mktempdir()
    Specs.solve_spec(TINY_CZ_TOML; run_dir = dir, format = :toml, max_iter = 15)

    @test isfile(joinpath(dir, "result.toml"))
    r = TOML.parsefile(joinpath(dir, "result.toml"))
    @test r["schema_version"] == "1"                                  # result/v1 unchanged
    @test haskey(r, "fidelity") && haskey(r, "iterations") && haskey(r, "wall_seconds")
    @test haskey(r["params"], "structure_hash") && haskey(r["params"], "problem_hash")
    @test haskey(r["params"], "converged")
    @test haskey(r["params"], "versions") && haskey(r["params"]["versions"], "Piccolo")
    @test length(filter(f -> endswith(f, ".jld2"), readdir(dir))) == 1  # pulse, not trajectory

    # spec→canonical→parse→canonical identity (Phase-1 extract_spec substitute).
    spec = Specs.parse_spec(TINY_CZ_TOML; format = :toml)
    cj = Specs.canonical_json(Specs.full_dict(spec))
    @test Specs.canonical_json(Specs.full_dict(Specs.parse_spec(cj; format = :json))) == cj
    # the persisted spec matches the canonical form.
    @test r["params"]["spec"] == cj
end

@testitem "solve_spec emits AMICODE_ITER + terminal DONE sentinels" begin
    using Piccolo, Piccolo.Specs

    TINY = """
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
    dir = mktempdir()
    # Capture stdout via a temp file: this Julia version's `redirect_stdout` does
    # not accept a bare IOBuffer (only a Pipe / IOStream / file).
    out = mktemp() do path, io
        redirect_stdout(io) do
            Specs.solve_spec(TINY; run_dir = dir, format = :toml, max_iter = 12)
        end
        flush(io)
        read(path, String)
    end
    # exact contract amico-run's telemetry.classifyLine parses:
    #   iter line: startsWith("AMICODE_ITER") with k=v tokens
    #   done line: /^DONE(\s|$)/
    @test occursin(r"(?m)^AMICODE_ITER iter=\d+", out)
    @test occursin(r"(?m)^DONE fidelity=", out)
    @test occursin(r"(?m)^AMICODE_PULSE_META ", out)
end

@testitem "solve_spec: altissimo backend is deferred (schema-only)" begin
    using Piccolo, Piccolo.Specs

    ALT_TOML = """
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
    [solver]
    backend = "altissimo"
    """
    dir = mktempdir()
    @test_throws Specs.SpecValidationError Specs.solve_spec(
        ALT_TOML;
        run_dir = dir,
        format = :toml,
    )
end

@testitem "solve_spec: explicitly declared ipopt backend executes and telemetries" begin
    using Piccolo, Piccolo.Specs

    IPOPT_TOML = """
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
    [solver]
    backend = "ipopt"
    """
    dir = mktempdir()
    out = mktemp() do path, io
        redirect_stdout(io) do
            Specs.solve_spec(IPOPT_TOML; run_dir = dir, format = :toml, max_iter = 12)
        end
        flush(io)
        read(path, String)
    end
    # The selectable Ipopt arm keeps the full telemetry contract — the raw
    # `callback` factory path, identical lines to the default MadNLP arm.
    @test occursin(r"(?m)^AMICODE_ITER iter=\d+ obj=\S+ inf_pr=\S+ inf_du=\S+", out)
    @test occursin(r"(?m)^DONE fidelity=", out)
    @test isfile(joinpath(dir, "result.toml"))
end

@testitem "solve_spec: unknown backend names the live allowed set" begin
    using Piccolo, Piccolo.Specs

    BOGUS_TOML = """
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
    [solver]
    backend = "gurobi"
    """
    dir = mktempdir()
    err = try
        Specs.solve_spec(BOGUS_TOML; run_dir = dir, format = :toml)
        nothing
    catch e
        e isa Specs.SpecValidationError ? e : rethrow()
    end
    @test err !== nothing
    # The structured error names every registered backend, madnlp first (the
    # default) — the retirement of the "ipopt only" pin (#360).
    @test err.errors[1].allowed == ["madnlp", "ipopt", "altissimo"]
end
