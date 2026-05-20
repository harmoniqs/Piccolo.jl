# ============================================================================ #
# Problem inspection — gather all data needed for `display(qcp)`
# ============================================================================ #
#
# `inspect(qcp)` returns a `ProblemInspection`, which is a pre-rendered bundle
# of strings + a few numbers. Rendering (in `show.jl`) reads only from this
# struct, never directly from the QCP — keeping the render layer pure.

"""
    ComponentInfo

One row in the Trajectory section of the rich display.
"""
struct ComponentInfo
    name::Symbol
    dim::Int
    role::Symbol       # :state | :control | :timestep
    bound_repr::String # "±0.0195" / "[0.01, 0.5]" / "—"
    bounded::Bool      # whether actually bounded (vs (-Inf, Inf))
end

"""
    GlobalInfo

One row for a global (free phase, calibration constant, etc.) in the rich display.
"""
struct GlobalInfo
    name::Symbol
    dim::Int
    bound_repr::String
    status::Symbol  # :free | :calibration | :pinned | :free_phase
end

"""
    ObjectiveTermInfo

One row in the Objective section.
"""
struct ObjectiveTermInfo
    name::String        # "KetInfidelity" / "QuadraticRegularizer(u)" / etc.
    weight::Float64     # multiplier
    value::Float64      # current J_i(x₀)
end

"""
    ConstraintInfo

One row in the Constraints section.
"""
struct ConstraintInfo
    name::String        # "FinalKetFidelity F≥0.99" / "BilinearIntegrator(ψ̃)"
    kind::Symbol        # :eq | :ineq | :bnd | :dynamics
    dim::Int            # constraint dimension
    violation::Float64  # ‖g(x₀)‖∞ or 0 for feasible bounds
    feasible::Bool
end

"""
    ProblemInspection

Everything `Base.show` needs to render `display(qcp)`. Built by `inspect(qcp)`.
"""
struct ProblemInspection
    # Identity
    traj_typename::String
    pulse_typename::String
    integrator_summary::Vector{String}

    # System
    sys_dim::Int
    n_drives::Int
    subsystem_levels::Union{Vector{Int},Nothing}
    sys_globals::Vector{Pair{Symbol,Float64}}

    # Trajectory
    N::Int
    T::Float64
    Δt_range::Union{Nothing,Tuple{Float64,Float64}}
    components::Vector{ComponentInfo}
    globals::Vector{GlobalInfo}

    # Goal
    goal_summary::String

    # Objective
    objective_terms::Vector{ObjectiveTermInfo}
    objective_total::Float64

    # Constraints
    constraints::Vector{ConstraintInfo}

    # Fidelity at the current trajectory point. `F_current` is the raw fidelity
    # against `qtraj.goal`. `F_with_phase` is the fidelity with the trajectory's
    # current φ_* globals applied as Z-phase rotations (matches what the
    # free-phase constraint enforces).
    F_current::Union{Nothing,Float64}
    F_with_phase::Union{Nothing,Float64}

    # NLP stats
    n_vars::Int
    n_bounded_vars::Int
    n_eq::Int
    n_ineq::Int
    jac_nnz::Union{Nothing,Int}
    hess_nnz::Union{Nothing,Int}
end

# ---------------------------------------------------------------------------- #
# Inspection entry point
# ---------------------------------------------------------------------------- #

# Strip parameterized types: KetTrajectory{CubicSplinePulse{…}} → "KetTrajectory"
_typename(T::Type) = string(nameof(T))
_typename(x) = _typename(typeof(x))

"""
    inspect(qcp::QuantumControlProblem) -> ProblemInspection

Gather everything needed for the rich `display(qcp)` view. Evaluates the
objective and constraints once at the current trajectory point.
"""
function inspect(qcp::QuantumControlProblem)
    prob = qcp.prob
    traj = prob.trajectory
    qtraj = qcp.qtraj
    sys = QuantumControlProblems.get_system(qcp)

    # Identity
    traj_typename = _typename(qtraj)
    pulse_typename = hasproperty(qtraj, :pulse) ? _typename(qtraj.pulse) : "—"

    integrator_summary = String[_render_integrator(intg) for intg in prob.integrators]

    # System
    sys_dim = _system_dim(sys)
    n_drives =
        hasproperty(sys, :n_drives) ? sys.n_drives :
        length(prob.trajectory.control_names) - 1
    subsystem_levels = _subsystem_levels(sys)
    sys_globals = _system_globals(sys, traj)

    # Trajectory
    N = traj.N
    T = _safe_duration(traj)
    Δt_range = _Δt_range(traj)
    components = _component_infos(traj, qtraj)
    globals = _global_infos(traj, prob.constraints)

    # Goal
    goal_summary = _goal_summary(qtraj)

    # Objective decomposition
    objective_terms, objective_total = _objective_terms(prob.objective, traj)

    # Constraint evaluation (includes dynamics from integrators)
    constraints = _constraint_infos(prob, traj)

    # Fidelity at the current trajectory point (raw + phase-adjusted)
    F_current, F_with_phase = _initial_fidelities(qtraj, traj)

    # NLP stats — variables + bounded count from trajectory; eq/ineq via constraint sums
    n_vars = traj.dim * traj.N + traj.global_dim
    n_bounded_vars = _count_bounded_vars(traj)
    n_eq = sum(c.dim for c in constraints if c.kind in (:eq, :dynamics); init = 0)
    n_ineq = sum(c.dim for c in constraints if c.kind === :ineq; init = 0)

    # Sparsity — best-effort, skip if anything errors (would need Evaluator construction)
    jac_nnz, hess_nnz = nothing, nothing

    return ProblemInspection(
        traj_typename,
        pulse_typename,
        integrator_summary,
        sys_dim,
        n_drives,
        subsystem_levels,
        sys_globals,
        N,
        T,
        Δt_range,
        components,
        globals,
        goal_summary,
        objective_terms,
        objective_total,
        constraints,
        F_current,
        F_with_phase,
        n_vars,
        n_bounded_vars,
        n_eq,
        n_ineq,
        jac_nnz,
        hess_nnz,
    )
end

# ---------------------------------------------------------------------------- #
# Helpers — system
# ---------------------------------------------------------------------------- #

function _system_dim(sys)
    # QuantumSystem types vary; try a few accessors
    if hasproperty(sys, :levels)
        return prod(sys.levels)
    elseif hasproperty(sys, :H_drift)
        return size(sys.H_drift, 1)
    elseif hasmethod(size, Tuple{typeof(sys)})
        return size(sys, 1)
    else
        return 0
    end
end

function _subsystem_levels(sys)
    hasproperty(sys, :levels) || return nothing
    levels = sys.levels
    levels isa AbstractVector{<:Integer} && return collect(Int, levels)
    return nothing
end

# Pull system-level globals (χ, K_q, ...) — those that are also present in
# `traj.global_components` are knobs; print their current value.
function _system_globals(sys, traj)
    isempty(traj.global_components) && return Pair{Symbol,Float64}[]
    if hasproperty(sys, :global_params)
        out = Pair{Symbol,Float64}[]
        for (name, val) in pairs(sys.global_params)
            haskey(traj.global_components, name) || continue
            try
                push!(out, name => Float64(val))
            catch
                # Skip non-scalar globals
            end
        end
        return out
    end
    return Pair{Symbol,Float64}[]
end

# ---------------------------------------------------------------------------- #
# Helpers — trajectory
# ---------------------------------------------------------------------------- #

_safe_duration(traj) =
    try
        # `get_duration(traj) = get_times(traj)[end]` — the correct cumulative end
        # time. `sum(get_timesteps)` overcounts by one Δt and is not the same thing.
        Float64(get_duration(traj))
    catch
        NaN
    end

function _Δt_range(traj)
    haskey(traj.bounds, :Δt) || return nothing
    lo, hi = traj.bounds[:Δt]
    lo isa AbstractVector || return nothing
    isempty(lo) && return nothing
    return (Float64(first(lo)), Float64(first(hi)))
end

function _component_infos(traj, qtraj)
    out = ComponentInfo[]
    state_set = Set(traj.state_names)
    control_set = Set(traj.control_names)
    for (name, idx) in pairs(traj.components)
        dim = length(idx)
        role = if name === traj.timestep
            :timestep
        elseif name in state_set
            :state
        elseif name in control_set
            :control
        else
            :other
        end
        bound_repr, bounded = _bound_repr(traj.bounds, name)
        push!(out, ComponentInfo(name, dim, role, bound_repr, bounded))
    end
    return out
end

function _bound_repr(bounds, name)
    haskey(bounds, name) || return ("—", false)
    lo, hi = bounds[name]
    if all(isinf, lo) && all(isinf, hi)
        return ("—", false)
    end
    # Symmetric vector?
    if lo isa AbstractVector && hi isa AbstractVector
        if length(lo) == length(hi) && all(lo .== .-hi)
            n = length(hi)
            if n <= 4
                return ("±" * string(round.(hi; sigdigits = 3)), true)
            else
                return ("±[" * _summarize_vec(hi) * "]", true)
            end
        else
            return ("[" * _summarize_vec(lo) * ", " * _summarize_vec(hi) * "]", true)
        end
    end
    return (string(lo, "..", hi), true)
end

function _summarize_vec(v)
    if length(v) == 0
        return ""
    elseif length(v) <= 4
        return join((string(round(x; sigdigits = 3)) for x in v), ", ")
    else
        first3 = join((string(round(x; sigdigits = 3)) for x in v[1:3]), ", ")
        return "$first3, … ($(length(v)) total)"
    end
end

function _global_infos(traj, constraints)
    out = GlobalInfo[]
    isempty(traj.global_components) && return out

    # Build a map: global name → constraint involvement
    pinned_names = Set{Symbol}()
    bounded_names = Dict{Symbol,Any}()  # name → bounds tuple
    for c in constraints
        # GlobalEqualityConstraint pins a value
        if _typename(c) == "GlobalEqualityConstraint" && hasproperty(c, :name)
            push!(pinned_names, c.name)
        end
        # GlobalBoundsConstraint adds bounds
        if _typename(c) == "GlobalBoundsConstraint" && hasproperty(c, :name)
            if hasproperty(c, :bounds)
                bounded_names[c.name] = c.bounds
            end
        end
    end

    for name in keys(traj.global_components)
        idx = traj.global_components[name]
        dim = length(idx)

        # Status detection
        status = if name in pinned_names
            :pinned
        elseif startswith(string(name), "φ_")
            :free_phase
        elseif haskey(bounded_names, name)
            :free
        else
            :free
        end

        # Bound rendering
        bound_repr = if haskey(bounded_names, name)
            b = bounded_names[name]
            _repr_global_bound(b)
        elseif status === :pinned
            "pinned"
        else
            "—"
        end

        push!(out, GlobalInfo(name, dim, bound_repr, status))
    end
    return out
end

function _repr_global_bound(b)
    if b isa Tuple && length(b) == 2
        lo, hi = b
        if lo isa AbstractVector && hi isa AbstractVector
            if length(lo) == 1 && lo[1] == -hi[1]
                v = hi[1]
                # Recognize 2π
                if isapprox(v, 2π; rtol = 1e-6)
                    return "±2π"
                end
                return "±$(round(v; sigdigits = 4))"
            end
            return "[$(round.(lo; sigdigits=4)), $(round.(hi; sigdigits=4))]"
        elseif lo isa Number && hi isa Number
            return lo == -hi ? "±$(round(hi; sigdigits=4))" :
                   "[$(round(lo; sigdigits=4)), $(round(hi; sigdigits=4))]"
        end
    end
    return string(b)
end

function _count_bounded_vars(traj)
    # Per-knot bounded components contribute dim*N variables each
    n = 0
    for (name, idx) in pairs(traj.components)
        if haskey(traj.bounds, name)
            lo, hi = traj.bounds[name]
            if lo isa AbstractVector && !all(isinf, lo) && !all(isinf, hi)
                n += length(idx) * traj.N
            end
        end
    end
    # Globals: we have no per-component bounds map readily; treat all globals as bounded
    # if any GlobalBoundsConstraint exists. Conservative: count all globals.
    n += traj.global_dim
    return n
end

# ---------------------------------------------------------------------------- #
# Helpers — goal description
# ---------------------------------------------------------------------------- #

_goal_summary(qtraj::UnitaryTrajectory) = _goal_summary_unitary(qtraj.goal)
_goal_summary(qtraj::KetTrajectory) = "|ψ_init⟩ → |ψ_goal⟩  (dim=$(length(qtraj.goal)))"
_goal_summary(
    qtraj::MultiKetTrajectory,
) = "$(length(qtraj.goals)) state transfers  (dim=$(length(first(qtraj.goals))))"
_goal_summary(qtraj::DensityTrajectory) = "ρ_init → ρ_goal"
_goal_summary(qtraj::SamplingTrajectory) = "sampled ensemble (n=$(length(qtraj.systems)))"
_goal_summary(qtraj::AbstractQuantumTrajectory) = string(_typename(qtraj), " goal")

function _goal_summary_unitary(goal)
    if goal isa EmbeddedOperator
        sublevels = goal.subsystem_levels
        return "EmbeddedOperator on $(sublevels), subspace dim $(length(goal.subspace))"
    end
    return "U_goal  ($(size(goal, 1))×$(size(goal, 2)))"
end

# ---------------------------------------------------------------------------- #
# Helpers — integrator rendering
# ---------------------------------------------------------------------------- #

function _render_integrator(intg)
    tn = _typename(intg)
    # Try to identify the state symbol the integrator acts on
    sym = if hasproperty(intg, :state_name)
        ":$(intg.state_name)"
    elseif hasproperty(intg, :variable)
        ":$(intg.variable)"
    elseif hasproperty(intg, :name)
        ":$(intg.name)"
    else
        ""
    end
    return isempty(sym) ? tn : "$tn($sym)"
end

# ---------------------------------------------------------------------------- #
# Helpers — objective decomposition
# ---------------------------------------------------------------------------- #

function _objective_terms(obj, traj)
    terms = ObjectiveTermInfo[]
    total = 0.0

    # Flat list of (sub-objective, weight) pairs
    pairs_iter = if hasproperty(obj, :objectives) && hasproperty(obj, :weights)
        zip(obj.objectives, obj.weights)
    else
        [(obj, 1.0)]
    end

    for (subobj, w) in pairs_iter
        name = _objective_label(subobj)
        val = _try_eval_objective(subobj, traj)
        weighted = w * val
        total += weighted
        push!(terms, ObjectiveTermInfo(name, w, weighted))
    end

    return terms, total
end

function _objective_label(obj)
    tn = _typename(obj)
    # Append :name if it's specific (e.g. QuadraticRegularizer has obj.name)
    extra = if hasproperty(obj, :name) && obj.name isa Symbol
        "(:" * string(obj.name) * ")"
    elseif hasproperty(obj, :state_name) && obj.state_name isa Symbol
        "(:" * string(obj.state_name) * ")"
    else
        ""
    end
    return tn * extra
end

function _try_eval_objective(obj, traj)
    try
        v = obj(traj)
        return v isa Number ? Float64(v) : 0.0
    catch
        return NaN
    end
end

# ---------------------------------------------------------------------------- #
# Helpers — constraint classification + evaluation
# ---------------------------------------------------------------------------- #

function _constraint_infos(prob, traj)
    out = ConstraintInfo[]

    # Dynamics from integrators (implicit equality constraints)
    for intg in prob.integrators
        viol = _eval_integrator_violation(intg, traj)
        dim = _integrator_dim(intg, traj)
        push!(
            out,
            ConstraintInfo(
                _render_integrator(intg),
                :dynamics,
                dim,
                viol,
                isnan(viol) ? true : viol < 1e-8,
            ),
        )
    end

    # User-facing constraints
    for c in prob.constraints
        push!(out, _classify_and_eval(c, traj))
    end
    return out
end

function _classify_and_eval(c, traj)
    tn = _typename(c)
    name = _constraint_name(c, tn, traj)
    kind = _constraint_kind(c, tn)
    dim = _constraint_dim(c)
    viol = _try_eval_constraint(c, traj, kind)
    feasible = isnan(viol) || viol < 1e-6
    return ConstraintInfo(name, kind, dim, viol, feasible)
end

function _constraint_name(c, tn, traj)
    # Direct properties from Piccolo-specific subclasses
    if hasproperty(c, :final_fidelity)
        return "FinalFidelity F ≥ $(round(c.final_fidelity; sigdigits = 4))"
    elseif hasproperty(c, :fidelity)
        return "FinalFidelity F ≥ $(round(c.fidelity; sigdigits = 4))"
    end

    # FinalKet/Unitary/Density fidelity constraints are built as
    # NonlinearKnotPointConstraint or NonlinearGlobalKnotPointConstraint at
    # the terminal time with equality=false. Detect the pattern.
    var_label = _constraint_var_label(c)
    if occursin("NonlinearKnotPointConstraint", tn) ||
       occursin("NonlinearGlobalKnotPointConstraint", tn)
        is_terminal = hasproperty(c, :times) && !isempty(c.times) && c.times[end] == traj.N
        is_ineq = hasproperty(c, :equality) && !c.equality
        if is_terminal && is_ineq
            return "Terminal fidelity bound$(isempty(var_label) ? "" : " on $var_label")"
        end
        kind_word = is_ineq ? "ineq" : "eq"
        return "Knot-point $kind_word$(isempty(var_label) ? "" : " on $var_label")"
    end

    if !isempty(var_label)
        return "$tn$var_label"
    end
    return tn
end

# Read a printable variable label from a constraint — supports `name`, `var_names`,
# `state_name`, `variable`. Returns "" when none applies.
function _constraint_var_label(c)
    if hasproperty(c, :name)
        n = c.name
        n isa Symbol && return ":$n"
    end
    if hasproperty(c, :var_names)
        ns = c.var_names
        if ns isa AbstractVector && !isempty(ns)
            return ":" * join(string.(ns), ",:")
        end
    end
    if hasproperty(c, :state_name)
        sn = c.state_name
        sn isa Symbol && return ":$sn"
    end
    return ""
end

function _constraint_kind(c, tn)
    if occursin("Bounds", tn)
        return :bnd
    elseif occursin("Equality", tn) ||
           occursin("TimeSteps", tn) ||
           occursin("TimeConsistency", tn)
        return :eq
    elseif hasproperty(c, :equality)
        return c.equality ? :eq : :ineq
    elseif occursin("Fidelity", tn) || occursin("Leakage", tn)
        return :ineq
    end
    return :ineq
end

function _constraint_dim(c)
    for prop in (:g_dim, :dim, :n_constraints)
        hasproperty(c, prop) || continue
        v = getproperty(c, prop)
        v isa Integer && return Int(v)
    end
    return 1
end

function _try_eval_constraint(c, traj, kind)
    kind === :bnd && return 0.0  # bounds: trust feasibility unless we detect otherwise
    try
        n = _constraint_dim(c)
        g = zeros(n)
        DirectTrajOpt.CommonInterface.evaluate!(g, c, traj)
        return Float64(maximum(abs, g))
    catch
        return NaN
    end
end

function _eval_integrator_violation(intg, traj)
    try
        # Integrators expose dynamics defects via the same evaluate! interface
        n = _integrator_dim(intg, traj)
        g = zeros(n)
        DirectTrajOpt.CommonInterface.evaluate!(g, intg, traj)
        return Float64(maximum(abs, g))
    catch
        return NaN
    end
end

function _integrator_dim(intg, traj)
    # Try to read state dimension off the integrator; fall back to traj
    for prop in (:state_dim, :dim, :n)
        hasproperty(intg, prop) || continue
        v = getproperty(intg, prop)
        v isa Integer && return Int(v) * (traj.N - 1)
    end
    return traj.N - 1
end

# ---------------------------------------------------------------------------- #
# Helpers — initial fidelity
# ---------------------------------------------------------------------------- #

function _initial_fidelities(qtraj, traj)
    F = try
        Float64(fidelity(qtraj))
    catch
        nothing
    end
    F_phase = _fidelity_with_stored_phases(qtraj, traj)
    return F, F_phase
end

# Compute F using the φ globals stored in the trajectory. Before solve, those
# are usually zero (so F_phase ≈ F_raw). After solve, the optimizer has driven
# them to their optimum — so this matches what the fidelity-constraint enforces.
function _fidelity_with_stored_phases(qtraj, traj)
    θ_names =
        sort!([n for n in keys(traj.global_components) if startswith(string(n), "φ_")])
    isempty(θ_names) && return nothing

    # Read φ values from the trajectory's global data
    θ = try
        vcat([Vector{Float64}(traj.global_data[traj.global_components[n]]) for n in θ_names]...)
    catch
        return nothing
    end

    # Pull final state + base goal off the quantum trajectory and compute F.
    try
        return _fidelity_at(qtraj, traj, θ)
    catch
        return nothing
    end
end

function _fidelity_at(qtraj::UnitaryTrajectory, traj, θ)
    goal = qtraj.goal
    goal isa EmbeddedOperator || return nothing
    Ũ⃗ = traj[isomorphism_state_name(qtraj)][:, end]
    # `unitary_fidelity_loss` is named misleadingly: it returns F, not 1-F.
    return Float64(QuantumObjectives.unitary_fidelity_loss(Ũ⃗, _phased_goal(goal, θ)))
end

function _fidelity_at(qtraj::KetTrajectory, traj, θ)
    # subsystem_levels for the phase rotation aren't stored on the qtraj; we
    # need them to build the rotation. Best-effort: try to infer from system.
    sys = QuantumControlProblems.get_system(qtraj)
    levels = hasproperty(sys, :levels) ? sys.levels : nothing
    levels === nothing && return nothing
    ψ̃ = traj[isomorphism_state_name(qtraj)][:, end]
    ψ_goal_phased = _phased_ket_goal(qtraj.goal, θ, levels)
    return Float64(QuantumObjectives.ket_fidelity_loss(ψ̃, ψ_goal_phased))
end

_fidelity_at(::AbstractQuantumTrajectory, _, _) = nothing

# Apply Z(θ) phase rotation to an EmbeddedOperator goal. This is the same
# transformation that `_make_free_phase_goal` (in ProblemTemplates) applies;
# we re-derive it here to keep this submodule independent.
function _phased_goal(op::EmbeddedOperator, θ)
    U_base = unembed(op)
    levels = op.subsystem_levels
    n_qubits = length(levels)
    n_sub = size(U_base, 1)
    phase_diag = map(1:n_sub) do i
        bits = i - 1
        phase = sum(
            θ[j] for j = 1:n_qubits if (bits >> (n_qubits - j)) & 1 == 1;
            init = zero(eltype(θ)),
        )
        return exp(im * phase)
    end
    phased = Diagonal(phase_diag) * U_base
    return EmbeddedOperator(Matrix(phased), op.subspace, levels)
end

function _phased_ket_goal(ψ_goal::AbstractVector, θ, levels::AbstractVector{<:Integer})
    n_qubits = length(levels)
    length(θ) == n_qubits || return ψ_goal
    # Build a per-subsystem diagonal phase: dim = prod(levels). For each basis
    # index, compute Σ_j θ_j * <j-th qubit excited indicator>.
    d = prod(levels)
    phased = similar(ψ_goal)
    for i = 1:d
        bits = i - 1
        phase = sum(
            θ[j] for j = 1:n_qubits if (bits >> (n_qubits - j)) & 1 == 1;
            init = zero(eltype(θ)),
        )
        phased[i] = exp(im * phase) * ψ_goal[i]
    end
    return phased
end
