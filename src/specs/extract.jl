export extract_spec

# ===========================================================================
# extract_spec — the inverse of `materialize`, with teeth.
#
# Phase 1 could only assert the *wire-level* identity
# spec → canonical → parse → canonical. That never touches the materializer:
# it proves the parser and the canonicaliser agree, nothing more.
#
# Phase 1b makes the round-trip real. `materialize` retains the originating
# spec on the (mutable, now template-tagged) problem, and `extract_spec` returns
# it only after **verifying it against the materialized object**: the tag type,
# the trajectory and pulse types, the wrapper nesting, and the params. So
# `spec → materialize → extract_spec → spec` passes only if the materializer
# actually built what the spec asked for — a spec/materializer consistency check,
# not a tautology.
#
# Types and `params` alone cannot reconstruct a spec (a built `QuantumSystem` does
# not retain its registry template name, an `EmbeddedOperator` does not retain the
# `GATES` name, and the solver block never enters the problem), so for a problem
# with no retained spec `extract_spec` is explicitly best-effort: the system is
# emitted as `raw` matrices and the goal as an inline matrix.
# ===========================================================================

@doc raw"""
    extract_spec(problem::AbstractQuantumControlProblem) -> ProblemSpec

Recover the declarative spec of a live problem.

Two paths:

1. **Materialized problem** (`retained_spec(problem) !== nothing`) — the retained
   spec is returned **after verification against the object**: template tag,
   trajectory kind, pulse kind, wrapper nesting (outermost → innermost), and the
   retained `params`. A mismatch is a [`SpecValidationError`](@ref) with field
   paths, because it means the materializer and the spec disagree.

2. **Hand-built problem** (`retained_spec(problem) === nothing`) — a **best-effort,
   explicitly non-canonical** spec: `system.kind = "raw"` with the drift/drive
   matrices inlined, the goal as an inline matrix, and the `[problem]` block
   reconstructed from the template tag plus the retained params. The `[solver]`
   block is defaulted (it never entered the problem).

   A problem with no template tag at all (`NoTemplate` — assembled from a
   `DirectTrajOptProblem` by hand rather than through a template) has no spec form
   and raises a structured error: `problem.template` is required and there is
   nothing to put in it.

```julia
spec = parse_spec(toml; format=:toml)
qcp  = materialize(spec)
canonical_json(full_dict(extract_spec(qcp))) == canonical_json(full_dict(spec))   # true
```
"""
function extract_spec(problem::AbstractQuantumControlProblem)
    spec = Control.retained_spec(problem)
    spec === nothing && return _best_effort_spec(problem)
    spec isa ProblemSpec || throw(
        SpecValidationError([
            SpecError(
                "kind",
                "retained spec is not a control ProblemSpec",
                string(typeof(spec)),
            ),
        ]),
    )
    errs = SpecError[]
    _verify_spec!(problem, spec, errs)
    isempty(errs) || throw(SpecValidationError(errs))
    return spec
end

# ---------------------------------------------------------------------------
# Verification: does the object match the spec that claims to have produced it?
# ---------------------------------------------------------------------------

# pulse type ⇄ wire `pulse.kind`
const _PULSE_KIND_OF = Dict{Symbol,Symbol}(
    :ZeroOrderPulse => :zero_order,
    :LinearSplinePulse => :linear_spline,
    :CubicSplinePulse => :cubic_spline,
)

_pulse_kind(pulse) = get(_PULSE_KIND_OF, nameof(typeof(pulse)), nothing)

# trajectory type ⇄ wire trajectory kind
function _trajectory_kind(qtraj)
    n = nameof(typeof(qtraj))
    n === :UnitaryTrajectory && return :unitary
    n === :KetTrajectory && return :ket
    n === :DensityTrajectory && return :density
    n === :MultiKetTrajectory && return :multiket
    return nothing
end

"""
    _wrapper_chain(problem) -> Vector{Symbol}

The wrapper kinds of a problem, outermost first — the type-level mirror of the
spec's ordered `[[wrappers]]` list.
"""
function _wrapper_chain(problem::AbstractQuantumControlProblem)
    kinds = Symbol[]
    p = problem
    while p isa Control.AbstractProblemWrapper
        push!(kinds, Control.wrapper_kind(p))
        p = Control.inner(p)
    end
    return kinds
end

function _verify_spec!(problem, spec::ProblemSpec, errs::Vector{SpecError})
    p = spec.problem
    if p === nothing
        push!(errs, SpecError("problem", "retained spec has no [problem] block"))
        return errs
    end

    # --- wrapper nesting: the ordered [[wrappers]] list, outermost first --------
    have = _wrapper_chain(problem)
    want = Symbol[w.kind for w in spec.wrappers]
    if have != want
        push!(
            errs,
            SpecError(
                "wrappers",
                "materialized wrapper nesting does not match the spec",
                string(have),
                String[string(w) for w in want],
            ),
        )
    end

    base = Control.base_problem(problem)

    # --- template tag ----------------------------------------------------------
    decl = get(Control.ProblemTemplates.TEMPLATE_DECLARATIONS, p.template, nothing)
    tag_ok = true
    if decl === nothing
        tag_ok = false
        push!(
            errs,
            SpecError(
                "problem.template",
                "unknown problem template",
                string(p.template),
                sort!(
                    String[
                        string(k) for
                        k in keys(Control.ProblemTemplates.TEMPLATE_DECLARATIONS)
                    ],
                ),
            ),
        )
    elseif !(base isa decl.alias_type)
        tag_ok = false
        push!(
            errs,
            SpecError(
                "problem.template",
                "materialized problem is not an instance of this template's type",
                string(nameof(typeof(Control.template_tag(base)))),
            ),
        )
    end

    # --- trajectory / pulse types ---------------------------------------------
    qtraj = Control.quantum_trajectory(base)
    tk = _trajectory_kind(qtraj)
    want_tk =
        spec.trajectory.kind === nothing ?
        (spec.goal === nothing ? nothing : spec.goal.kind) : spec.trajectory.kind
    if want_tk !== nothing && tk !== want_tk
        push!(
            errs,
            SpecError(
                "trajectory.kind",
                "materialized trajectory kind does not match the spec",
                string(tk),
                [string(want_tk)],
            ),
        )
    end
    if spec.pulse !== nothing && hasproperty(qtraj, :pulse)
        pk = _pulse_kind(qtraj.pulse)
        if pk !== spec.pulse.kind
            push!(
                errs,
                SpecError(
                    "pulse.kind",
                    "materialized pulse kind does not match the spec",
                    string(pk),
                    [string(spec.pulse.kind)],
                ),
            )
        end
    end

    # --- params ---------------------------------------------------------------
    # Only meaningful once the tag matches: params structs are per-template, so a
    # tag mismatch already explains everything and comparing across two different
    # params types would be noise.
    if tag_ok
        want_params = try
            Control.ProblemTemplates._template_params(
                decl.params_type,
                Tuple(decl.passthrough),
                p.template,
                pairs(_template_kwargs(spec)),
            )
        catch
            nothing
        end
        got = Control.template_params(base)
        # Compare only the fields the spec can actually carry: per-template
        # extras (`integrator_type`, `du_bounds`, …) live outside `TemplateBlock`
        # by design ("best-effort, flagged as non-canonical" — see
        # `_params_to_template_kwargs`), and `materialize` may legitimately set
        # them beyond the spec (e.g. the SplinePulseProblem `:pwc` pulse-type
        # guard). A retained-params disagreement on such a field is not a spec
        # violation.
        if want_params !== nothing &&
           typeof(got) === typeof(want_params) &&
           !isequal(
               _spec_carried_params(got),
               _spec_carried_params(want_params),
           )
            push!(
                errs,
                SpecError(
                    "problem",
                    "retained template params disagree with the spec on: " *
                    join(_params_diff(got, want_params), ", "),
                    string(got),
                ),
            )
        end
    end

    return errs
end

"""
    _params_diff(got, want) -> Vector{String}

Which params fields disagree. A type mismatch is reported as a single `::type`
entry rather than a field list, since the two structs need not share fields.
"""
function _params_diff(got, want)
    typeof(got) === typeof(want) && return String[
        string(f) for
        f in fieldnames(typeof(want)) if !isequal(getfield(got, f), getfield(want, f))
    ]
    return String["::type ($(nameof(typeof(got))) vs $(nameof(typeof(want))))"]
end

# The subset of a params struct's fields the spec can carry (the TemplateBlock
# field set). Fields outside it are per-template extras the spec deliberately
# does not round-trip.
function _spec_carried_params(params::Control.AbstractTemplateParams)
    return NamedTuple{filter(
        f -> f in _TEMPLATE_BLOCK_FIELDS,
        fieldnames(typeof(params)),
    )}(
        getfield(params, f) for
        f in fieldnames(typeof(params)) if f in _TEMPLATE_BLOCK_FIELDS
    )
end

# ---------------------------------------------------------------------------
# Best-effort emission for a problem with no retained spec
# ---------------------------------------------------------------------------

_wire_complex(z) = Any[real(z), imag(z)]
_wire_matrix(M::AbstractMatrix) =
    Any[Any[_wire_complex(M[i, j]) for j in axes(M, 2)] for i in axes(M, 1)]

function _best_effort_spec(problem::AbstractQuantumControlProblem)
    base = Control.base_problem(problem)
    tag = Control.template_tag(base)
    if tag isa Control.NoTemplate
        throw(
            SpecValidationError([
                SpecError(
                    "problem.template",
                    "problem carries no template tag (it was assembled from a " *
                    "DirectTrajOptProblem by hand), so it has no spec form; " *
                    "materialize from a spec, or build it through a problem template",
                    "NoTemplate",
                ),
            ]),
        )
    end
    name = _declared_name(tag)
    name === nothing && throw(
        SpecValidationError([
            SpecError(
                "problem.template",
                "template tag is not registered in TEMPLATE_DECLARATIONS",
                string(nameof(typeof(tag))),
            ),
        ]),
    )

    sys = get_system(base)
    qtraj = Control.quantum_trajectory(base)
    traj = get_trajectory(base)

    # system: `raw` matrices — explicitly non-canonical (a built QuantumSystem does
    # not retain the registry template name it may have come from).
    system = SystemSpec(;
        kind = :raw,
        H_drift = _wire_matrix(Matrix(_drift_matrix(sys))),
        H_drives = Any[_wire_matrix(Matrix(H)) for H in _drive_matrices(sys)],
        params = Dict{Symbol,Any}(:drive_bounds => _wire_drive_bounds(sys)),
    )

    # goal: inline matrix / ket — the GATES name is not retained either.
    goal = _best_effort_goal(qtraj)

    pulse =
        hasproperty(qtraj, :pulse) ?
        PulseSpec(;
            kind = something(_pulse_kind(qtraj.pulse), :zero_order),
            T = _traj_duration(traj),
        ) : nothing

    params = Control.template_params(base)
    problem_block =
        TemplateBlock(; template = name, N = traj.N, _params_to_template_kwargs(params)...)

    return ProblemSpec(;
        schema_version = 1,
        kind = :control,
        system = system,
        goal = goal,
        pulse = pulse,
        trajectory = TrajectorySpec(; kind = _trajectory_kind(qtraj)),
        problem = problem_block,
        wrappers = WrapperSpec[WrapperSpec(; kind = k) for k in _wrapper_chain(problem)],
        solver = SolverSpec(),
    )
end

function _declared_name(tag)
    for (name, d) in Control.ProblemTemplates.TEMPLATE_DECLARATIONS
        d.tag_type === typeof(tag) && return name
    end
    return nothing
end

_drift_matrix(sys) =
    hasproperty(sys, :H_drift) ? sys.H_drift : zeros(ComplexF64, sys.levels, sys.levels)
# `sys.H_drives` holds `LinearDrive` wrappers, not bare matrices.
_drive_matrix(d) = hasproperty(d, :H) ? d.H : d
_drive_matrices(sys) =
    hasproperty(sys, :H_drives) ? [_drive_matrix(d) for d in sys.H_drives] : []
function _wire_drive_bounds(sys)
    hasproperty(sys, :drive_bounds) || return Any[]
    return Any[
        b isa Tuple ? Any[Float64(b[1]), Float64(b[2])] : Float64(b) for
        b in sys.drive_bounds
    ]
end

_traj_duration(traj) = haskey(traj.components, :Δt) ? sum(vec(traj.Δt)) : Float64(traj.N)

function _best_effort_goal(qtraj)
    hasproperty(qtraj, :goal) || return nothing
    g = qtraj.goal
    if g isa AbstractMatrix
        return GoalSpec(; kind = :unitary, matrix = _wire_matrix(g))
    elseif g isa EmbeddedOperator
        return GoalSpec(;
            kind = :unitary,
            matrix = _wire_matrix(Matrix(unembed(g))),
            subsystem_levels = collect(Int, g.subsystem_levels),
        )
    elseif g isa AbstractVector
        return GoalSpec(; kind = :ket, target = string("ComplexF64", collect(g)))
    end
    return nothing
end

# The subset of the params struct that `TemplateBlock` can carry. Fields with no
# `TemplateBlock` counterpart (per-template extras like `du_bounds`,
# `integrator_type`) are dropped — best-effort, and flagged as non-canonical.
const _TEMPLATE_BLOCK_FIELDS = Set{Symbol}(fieldnames(TemplateBlock))

function _params_to_template_kwargs(params::Control.AbstractTemplateParams)
    kw = Dict{Symbol,Any}()
    for f in fieldnames(typeof(params))
        f in _TEMPLATE_BLOCK_FIELDS || continue
        v = getfield(params, f)
        v === nothing && continue
        if f === :R_u || f === :R_du || f === :R_ddu
            v isa Real || continue          # TemplateBlock takes scalars only
            kw[f] = Float64(v)
        elseif f === :global_bounds
            kw[f] = Dict{Symbol,Any}(k => v[k] for k in keys(v))
        else
            kw[f] = v
        end
    end
    return kw
end

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@testitem "extract_spec: spec → materialize → extract_spec is the identity" begin
    using Piccolo, Piccolo.Specs

    TOML_SRC = """
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
    Q = 250.0
    """
    spec = Specs.parse_spec(TOML_SRC; format = :toml)
    qcp = Specs.materialize(spec)

    # the spec is retained on the problem …
    @test retained_spec(qcp) === spec
    # … and comes back out, canonical-JSON identical
    got = Specs.extract_spec(qcp)
    @test Specs.canonical_json(Specs.full_dict(got)) ==
          Specs.canonical_json(Specs.full_dict(spec))
    @test Specs.problem_hash(got) == Specs.problem_hash(spec)
    @test Specs.structure_hash(got) == Specs.structure_hash(spec)

    # the retained params really came from the spec
    @test template_params(qcp).Q == 250.0
end

@testitem "extract_spec: wrapper-nested round-trip + verification has teeth" begin
    using Piccolo, Piccolo.Specs

    SAMPLING_TOML = """
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
    [[wrappers]]
    kind = "sampling"
    variants = [ { "levels" = 3 }, { "levels" = 3 } ]
    """
    spec = Specs.parse_spec(SAMPLING_TOML; format = :toml)
    qcp = Specs.materialize(spec)

    @test qcp isa SamplingProblem
    @test retained_spec(qcp) === spec
    got = Specs.extract_spec(qcp)
    @test Specs.canonical_json(Specs.full_dict(got)) ==
          Specs.canonical_json(Specs.full_dict(spec))

    # the wrapper nesting is verified, not assumed: a spec claiming no wrapper
    # cannot be extracted from a wrapped problem
    unwrapped = Specs.ProblemSpec(;
        schema_version = spec.schema_version,
        kind = :control,
        system = spec.system,
        goal = spec.goal,
        pulse = spec.pulse,
        trajectory = spec.trajectory,
        problem = spec.problem,
        integrator = spec.integrator,
        wrappers = Specs.WrapperSpec[],
        solver = spec.solver,
    )
    retain_spec!(qcp, unwrapped)
    err = try
        Specs.extract_spec(qcp)
        nothing
    catch e
        e
    end
    @test err isa Specs.SpecValidationError
    @test any(e -> e.path == "wrappers", err.errors)

    # and so is the template tag: a spec naming a different template is rejected
    retain_spec!(qcp, unwrapped)
    wrong_template = Specs.ProblemSpec(;
        schema_version = 1,
        kind = :control,
        system = spec.system,
        goal = spec.goal,
        pulse = spec.pulse,
        trajectory = spec.trajectory,
        problem = Specs.TemplateBlock(; template = :SmoothPulseProblem, N = 11),
        wrappers = spec.wrappers,
        solver = spec.solver,
    )
    retain_spec!(qcp, wrong_template)
    err2 = try
        Specs.extract_spec(qcp)
        nothing
    catch e
        e
    end
    @test err2 isa Specs.SpecValidationError
    @test any(e -> e.path == "problem.template", err2.errors)
end

@testitem "extract_spec: best-effort for a template-built problem with no spec" begin
    using Piccolo, Piccolo.Specs
    using DirectTrajOpt
    using NamedTrajectories

    N = 8
    T = 2.0
    sys = QuantumSystem(0.1 * GATES[:Z], [GATES[:X]], [1.0])
    times = collect(range(0.0, T, length = N))
    qtraj = UnitaryTrajectory(sys, ZeroOrderPulse(zeros(1, N), times), GATES[:X])
    qcp = SmoothPulseProblem(
        qtraj,
        N;
        Q = 77.0,
        piccolo_options = PiccoloOptions(display = :silent),
    )
    @test retained_spec(qcp) === nothing

    spec = Specs.extract_spec(qcp)
    @test spec isa Specs.ProblemSpec
    @test spec.system.kind === :raw               # non-canonical by construction
    @test spec.system.H_drift !== nothing
    @test spec.problem.template === :SmoothPulseProblem
    @test spec.problem.N == N
    @test spec.problem.Q == 77.0
    @test spec.pulse.kind === :zero_order
    @test spec.goal.kind === :unitary
    # it is a real spec: it canonicalises and hashes
    @test !isempty(Specs.canonical_json(Specs.full_dict(spec)))
    @test length(Specs.problem_hash(spec)) == 64

    # a problem assembled by hand carries no template tag and therefore no spec form
    traj = NamedTrajectory(qtraj, N)
    bare = QuantumControlProblem(
        qtraj,
        DirectTrajOptProblem(
            traj,
            QuadraticRegularizer(:u, traj, 1.0),
            BilinearIntegrator(qtraj, N),
        ),
    )
    @test template_tag(bare) === NoTemplate()
    err = try
        Specs.extract_spec(bare)
        nothing
    catch e
        e
    end
    @test err isa Specs.SpecValidationError
    @test any(e -> e.path == "problem.template", err.errors)
end
