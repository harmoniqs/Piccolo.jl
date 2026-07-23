export materialize

using Random: MersenneTwister
using DirectTrajOpt:
    DirectTrajOptProblem,
    CompositeObjective,
    TerminalObjective,
    NullObjective,
    QuadraticRegularizer,
    MinimumTimeObjective,
    AbstractObjective,
    AbstractConstraint

# ===========================================================================
# materialize — the function barrier between untyped spec data and the typed
# Piccolo world. For `control` specs it looks up registry factories, builds
# system → goal → pulse → trajectory, calls the existing template function, then
# applies the composition axes (Task 8) and wrappers (Task 9). All validation is
# a pre-pass returning structured `SpecError`s — never a materialization crash.
# ===========================================================================

"""
    materialize(spec::ProblemSpec) -> QuantumControlProblem

Materialize a `control` [`ProblemSpec`](@ref) into a live `QuantumControlProblem`
by looking up the registered system/template factories and calling the existing
Piccolo template function with mapped keyword arguments.

Phase-1 scope (wire-format-first): `system.kind ∈ {:template, :raw}` (`:composite`
deferred); solver backend execution is `ipopt`-only (backend dispatch is the
runner's concern, Task 12); the `:robust` wrapper is schema-only and returns a
structured "deferred" [`SpecError`](@ref). All trait/compatibility violations are
collected and thrown as a [`SpecValidationError`](@ref) *before* any Piccolo
object is constructed.
"""
function materialize(spec::ProblemSpec)
    errs = SpecError[]
    _validate!(spec, errs)
    isempty(errs) || throw(SpecValidationError(errs))

    sys = _build_system(spec.system)
    goal = _build_goal(spec.goal, sys)
    qtraj = _build_pulse_trajectory(spec, sys, goal)
    base_qcp = _call_template(spec, qtraj)
    qcp = _apply_composition(spec, base_qcp)
    qcp = _apply_wrappers(spec, qcp)
    return qcp
end

# ---------------------------------------------------------------------------
# System
# ---------------------------------------------------------------------------

# Narrow a Vector{Any} of homogeneous numbers to a concrete Vector so template
# kwargs like `drive_bounds::Vector{Float64}` typecheck. TOML/JSON parsing boxes
# arrays as Vector{Any}; leave anything non-numeric untouched.
function _concretize(x)
    if x isa AbstractVector && !isempty(x)
        if all(v -> v isa Integer && !(v isa Bool), x)
            return Int[v for v in x]
        elseif all(v -> v isa Real && !(v isa Bool), x)
            return Float64[v for v in x]
        end
    end
    return x
end

_concretize_params(params::AbstractDict) =
    Dict{Symbol,Any}(k => _concretize(v) for (k, v) in params)

function _build_system(s::SystemSpec)
    if s.kind === :template
        entry = lookup_system(s.template)
        return entry.factory(; _concretize_params(s.params)...)
    elseif s.kind === :raw
        # Minimal Phase-1 raw path: build a QuantumSystem from explicit matrices.
        H_drift = _to_matrix(s.H_drift)
        H_drives = [_to_matrix(h) for h in s.H_drives]
        drive_bounds = _concretize(get(s.params, :drive_bounds, fill(1.0, length(H_drives))))
        return Quantum.QuantumSystem(H_drift, H_drives, drive_bounds)
    else
        # :composite deferred in Phase 1 (validation should have flagged it).
        throw(SpecValidationError([SpecError("system.kind",
            "composite systems are deferred in Phase 1 (template + raw only)",
            string(s.kind))]))
    end
end

# nested [re, im] pairs → ComplexF64 matrix (used by :raw H_drift/H_drives and
# inline unitary goal matrices).
function _to_matrix(x)
    x isa AbstractMatrix && return ComplexF64.(x)
    rows = collect(x)
    n = length(rows)
    M = Matrix{ComplexF64}(undef, n, length(rows[1]))
    for (i, row) in enumerate(rows), (j, entry) in enumerate(row)
        M[i, j] = entry isa AbstractVector ? ComplexF64(entry[1], entry[2]) : ComplexF64(entry)
    end
    return M
end

# ---------------------------------------------------------------------------
# Goal
# ---------------------------------------------------------------------------

function _build_goal(goal::GoalSpec, sys)
    if goal.kind === :unitary
        op = goal.gate !== nothing ? GATES[goal.gate] : _to_matrix(goal.matrix)
        if sys isa Quantum.CompositeQuantumSystem
            return EmbeddedOperator(op, sys)
        elseif goal.subspace !== nothing && length(goal.subspace) == 1
            return EmbeddedOperator(op, sys; subspace=goal.subspace[1], levels=sys.levels)
        else
            return EmbeddedOperator(op, sys)
        end
    elseif goal.kind === :ket
        # Phase-1 ket goals: `target`/`initial` are Julia-parseable complex vectors.
        return _parse_ket(goal.target)
    else
        throw(SpecValidationError([SpecError("goal.kind", "unknown goal kind",
            string(goal.kind), ["unitary", "ket"])]))
    end
end

_parse_ket(::Nothing) = throw(SpecValidationError([SpecError("goal.target",
    "ket goal requires a `target`")]))
_parse_ket(s::AbstractString) = ComplexF64.(eval(Meta.parse(s)))

# ---------------------------------------------------------------------------
# Pulse + trajectory (review correction 2: the pulse TYPE is fixed here, not by
# the template — build the concrete pulse ourselves, then wrap it).
# ---------------------------------------------------------------------------

function _controls(p::PulseSpec, n_drives::Int)
    N_seed = 11
    if p.init === :random
        rng = MersenneTwister(p.seed)
        return 0.1 .* randn(rng, n_drives, N_seed)
    else
        return zeros(n_drives, N_seed)
    end
end

function _build_pulse_trajectory(spec::ProblemSpec, sys, goal)
    p = spec.pulse
    N = spec.problem.N
    n_drives = sys.n_drives
    times = collect(range(0.0, p.T, N))
    controls = p.init === :random ? 0.1 .* randn(MersenneTwister(p.seed), n_drives, N) :
               zeros(n_drives, N)
    pulse = if p.kind === :zero_order
        ZeroOrderPulse(controls, times)
    elseif p.kind === :linear_spline
        LinearSplinePulse(controls, times)
    elseif p.kind === :cubic_spline
        CubicSplinePulse(controls, zeros(n_drives, N), times)
    else
        throw(SpecValidationError([SpecError("pulse.kind", "unknown pulse kind",
            string(p.kind), ["zero_order", "linear_spline", "cubic_spline"])]))
    end

    gk = spec.goal.kind
    if gk === :unitary
        return UnitaryTrajectory(sys, pulse, goal)
    elseif gk === :ket
        ψ0 = spec.goal.initial === nothing ?
             ComplexF64[i == 1 ? 1.0 : 0.0 for i in 1:sys.levels] :
             _parse_ket(spec.goal.initial)
        return KetTrajectory(sys, pulse, ψ0, goal)
    else
        throw(SpecValidationError([SpecError("goal.kind", "unknown goal kind",
            string(gk), ["unitary", "ket"])]))
    end
end

# ---------------------------------------------------------------------------
# Integrator block (Task 9): bilinear/omitted → nothing (template default path);
# exponential/spline → Piccolissimo factory (not loaded in this worktree, so
# validation blocks those before we get here).
# ---------------------------------------------------------------------------

function _build_integrator(spec::ProblemSpec, qtraj)
    spec.integrator === nothing && return nothing
    ik = spec.integrator.kind
    ik === :bilinear && return nothing
    entry = lookup_integrator(ik)
    return entry.factory(qtraj, spec.problem.N; alg=spec.integrator.alg)
end

# ---------------------------------------------------------------------------
# Template call (Task 7 happy path). Coercions per review correction 11.
# ---------------------------------------------------------------------------

# global_bounds::Dict{Symbol,Any} → the template's typed bounds dict.
function _coerce_global_bounds(gb::AbstractDict)
    out = Dict{Symbol,Union{Float64,Tuple{Float64,Float64}}}()
    for (k, v) in gb
        if v isa AbstractVector && length(v) == 2
            out[k] = (Float64(v[1]), Float64(v[2]))
        elseif v isa Tuple
            out[k] = (Float64(v[1]), Float64(v[2]))
        else
            out[k] = Float64(v)
        end
    end
    return out
end

function _call_template(spec::ProblemSpec, qtraj)
    p = spec.problem
    fac = lookup_template(p.template).factory
    kwargs = Dict{Symbol,Any}(:Q => p.Q, :R => p.R, :free_phase => p.free_phase)
    p.R_u === nothing || (kwargs[:R_u] = p.R_u)
    p.R_du === nothing || (kwargs[:R_du] = p.R_du)
    p.R_ddu === nothing || (kwargs[:R_ddu] = p.R_ddu)
    isfinite(p.du_bound) && (kwargs[:du_bound] = p.du_bound)
    p.ddu_bound === nothing || (kwargs[:ddu_bound] = Float64(p.ddu_bound))
    p.initial_phases === nothing || (kwargs[:initial_phases] = p.initial_phases)
    isempty(p.calibration_targets) || (kwargs[:calibration_targets] = p.calibration_targets)
    isempty(p.global_names) || (kwargs[:global_names] = p.global_names)
    isempty(p.global_bounds) || (kwargs[:global_bounds] = _coerce_global_bounds(p.global_bounds))
    p.free_dt isa Free && (kwargs[:Δt_bounds] = (p.free_dt.lo, p.free_dt.hi))
    intg = _build_integrator(spec, qtraj)
    intg === nothing || (kwargs[:integrator] = intg)
    return fac(qtraj, p.N; kwargs...)
end

# ---------------------------------------------------------------------------
# Composition (Task 8) — filled in below in the composition section.
# Wrappers (Task 9) — filled in below in the wrappers section.
# These stubs keep the Task-7 happy path a pass-through.
# ---------------------------------------------------------------------------

_apply_composition(spec::ProblemSpec, qcp) = qcp
_apply_wrappers(spec::ProblemSpec, qcp) = qcp
_validate!(spec::ProblemSpec, errs::Vector{SpecError}) = errs

@testitem "materialize: cubic-spline X gate matches hand-built problem" begin
    using Piccolo, Piccolo.Specs
    MAT_CONTROL_TOML = """
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
    spec = Specs.parse_spec(MAT_CONTROL_TOML; format=:toml)
    qcp = Specs.materialize(spec)
    @test qcp isa QuantumControlProblem

    # hand-built equivalent — the pulse TYPE is fixed at trajectory construction,
    # NOT by the template (see plan review correction 2).
    sys = TransmonSystem(; levels=3, drive_bounds=[0.02, 0.02])
    goal = EmbeddedOperator(GATES[:X], sys)
    times = collect(range(0.0, 10.0, 11))
    pulse = CubicSplinePulse(zeros(2, 11), zeros(2, 11), times)
    traj = UnitaryTrajectory(sys, pulse, goal)
    ref = SplinePulseProblem(traj, 11)

    # A single-term objective is a bare objective (no `.objectives`); normalize.
    nterms(J) = J isa DirectTrajOpt.CompositeObjective ? length(J.objectives) : 1
    @test nterms(qcp.prob.objective) == nterms(ref.prob.objective)
    @test get_trajectory(qcp).dims == get_trajectory(ref).dims
end
