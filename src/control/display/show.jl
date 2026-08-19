# ============================================================================ #
# Rich display rendering for QuantumControlProblem
# ============================================================================ #
#
# `Base.show(io, ::MIME"text/plain", qcp)` → :standard level (tree view)
# `Base.show(io, qcp)`                     → :compact one-line summary
# `show_problem(io, qcp; detail)`          → explicit entry point

using TestItems

# ---------------------------------------------------------------------------- #
# Tree-rendering glyphs
# ---------------------------------------------------------------------------- #

const _TREE_BRANCH = "├─ "
const _TREE_LAST = "└─ "
const _TREE_PIPE = "│  "
const _TREE_BLANK = "   "

# ---------------------------------------------------------------------------- #
# One-line compact display (Base.show with no MIME)
# ---------------------------------------------------------------------------- #

function Base.show(io::IO, qcp::QuantumControlProblem)
    ins = inspect(qcp)
    # QuantumControlProblem · KetTrajectory · CubicSplinePulse · 1793 vars · 1519 eq · F₀=0.378
    parts = String["QuantumControlProblem", ins.traj_typename]
    ins.pulse_typename == "—" || push!(parts, ins.pulse_typename)
    push!(parts, "$(ins.n_vars) vars")
    ins.n_eq > 0 && push!(parts, "$(ins.n_eq) eq")
    ins.n_ineq > 0 && push!(parts, "$(ins.n_ineq) ineq")
    # Prefer phase-adjusted F if available (matches what the free-phase
    # constraint actually enforces); otherwise show raw F.
    F_show = something(ins.F_with_phase, ins.F_current, nothing)
    isnothing(F_show) || push!(parts, @sprintf("F=%.3f", F_show))
    print(io, join(parts, " · "))
end

# ---------------------------------------------------------------------------- #
# Multi-line rich display (MIME"text/plain") — :standard level
# ---------------------------------------------------------------------------- #

function Base.show(io::IO, ::MIME"text/plain", qcp::QuantumControlProblem)
    show_problem(io, qcp; detail = :standard)
end

# ---------------------------------------------------------------------------- #
# `show_problem(io, qcp; detail=:standard|:full)` — explicit entry point
# ---------------------------------------------------------------------------- #

"""
    show_problem(io::IO, qcp::QuantumControlProblem; detail::Symbol = :standard)

Render the rich problem view to `io`. `detail` is one of:

- `:standard` (default) — the tree view: system / trajectory / goal / objective
  / constraints / initial fidelity / NLP stats.
- `:full` — everything in `:standard` plus integrator-by-integrator detail,
  sparsity, and a terminal pulse plot.
"""
function show_problem(io::IO, qcp::QuantumControlProblem; detail::Symbol = :standard)
    detail in (:standard, :full) ||
        throw(ArgumentError("detail must be :standard or :full; got :$detail"))
    ins = inspect(qcp)
    _print_tree(io, ins, qcp, detail)
end

show_problem(qcp::QuantumControlProblem; detail::Symbol = :standard) =
    show_problem(stdout, qcp; detail = detail)

# ---------------------------------------------------------------------------- #
# Tree rendering
# ---------------------------------------------------------------------------- #

function _print_tree(io::IO, ins::ProblemInspection, qcp, detail::Symbol)
    # Header
    printstyled(io, "QuantumControlProblem\n"; bold = true)

    # Top stem
    pulse_part = ins.pulse_typename == "—" ? "" : "  ·  $(ins.pulse_typename)"
    integ_part =
        isempty(ins.integrator_summary) ? "" : "  ·  " * join(ins.integrator_summary, ", ")
    print(io, _TREE_BRANCH, ins.traj_typename, pulse_part, integ_part, "\n")
    print(io, _TREE_PIPE, "\n")

    _print_system(io, ins)
    _print_trajectory(io, ins)
    _print_goal(io, ins)
    _print_objective(io, ins)
    _print_constraints(io, ins)
    _print_status(io, ins, last = (detail !== :full))

    if detail === :full
        _print_full_extras(io, ins, qcp)
    else
        # Hint about how to drill down
        print(io, "\n")
        printstyled(io, "Hint: ", color = :light_black)
        printstyled(
            io,
            "show_problem(qcp; detail=:full) ",
            color = :light_black,
            bold = true,
        )
        printstyled(io, "for pulse plot + sparsity\n", color = :light_black)
    end
end

# ---------------------------------------------------------------------------- #
# Sections
# ---------------------------------------------------------------------------- #

function _print_system(io::IO, ins::ProblemInspection)
    print(io, _TREE_BRANCH)
    printstyled(io, "System"; bold = true)
    print(io, "\n")
    sub_lines = String[]
    parts = String["dim=$(ins.sys_dim)", "drives=$(ins.n_drives)"]
    if !isnothing(ins.subsystem_levels)
        push!(parts, "subsystems=$(ins.subsystem_levels)")
    end
    push!(sub_lines, join(parts, "   "))
    if !isempty(ins.sys_globals)
        gs = join(
            ("$name=$(round(val; sigdigits=4))" for (name, val) in ins.sys_globals),
            "   ",
        )
        push!(sub_lines, "globals: $gs")
    end
    for line in sub_lines
        print(io, _TREE_PIPE, "  ", line, "\n")
    end
    print(io, _TREE_PIPE, "\n")
end

function _print_trajectory(io::IO, ins::ProblemInspection)
    print(io, _TREE_BRANCH)
    printstyled(io, "Trajectory"; bold = true)
    print(io, "\n")
    head = "N=$(ins.N)   T=$(_fmt_T(ins.T))"
    if !isnothing(ins.Δt_range)
        head *= "   Δt∈[$(_fmt_T(ins.Δt_range[1])), $(_fmt_T(ins.Δt_range[2]))]"
    end
    print(io, _TREE_PIPE, "  ", head, "\n")

    # Components
    if !isempty(ins.components)
        max_name = maximum(length(string(c.name)) for c in ins.components)
        max_dim = maximum(length(string(c.dim)) for c in ins.components)
        # Pick a bound column width that fits the longest displayed value (capped)
        bound_w = min(40, maximum(length(c.bound_repr) for c in ins.components; init = 1))
        bound_w = max(bound_w, 6)
        for c in ins.components
            namepad = rpad(string(c.name), max_name)
            dimpad = lpad(string(c.dim), max_dim)
            role_str =
                c.role === :state ? "state" :
                c.role === :control ? "control" : c.role === :timestep ? "timestep" : ""
            bound_text = c.bound_repr == "—" ? "" : c.bound_repr
            bound_part = rpad(bound_text, bound_w)
            tick = c.bounded ? "✓" : "·"
            tick_color = c.bounded ? :green : :light_black
            print(io, _TREE_PIPE, "  ", namepad, "  (", dimpad, ")  ", bound_part, "  ")
            printstyled(io, tick; color = tick_color)
            print(io, "  ", role_str, "\n")
        end
    end

    # Globals
    if !isempty(ins.globals)
        print(io, _TREE_PIPE, "  ", "globals:\n")
        max_name = maximum(length(string(g.name)) for g in ins.globals)
        for g in ins.globals
            namepad = rpad(string(g.name), max_name)
            status_str = string(g.status)
            print(
                io,
                _TREE_PIPE,
                "    ",
                namepad,
                "  ",
                rpad(g.bound_repr, 12),
                " ",
                status_str,
                "\n",
            )
        end
    end
    print(io, _TREE_PIPE, "\n")
end

function _print_goal(io::IO, ins::ProblemInspection)
    print(io, _TREE_BRANCH)
    printstyled(io, "Goal"; bold = true)
    print(io, "\n")
    print(io, _TREE_PIPE, "  ", ins.goal_summary, "\n")
    print(io, _TREE_PIPE, "\n")
end

function _print_objective(io::IO, ins::ProblemInspection)
    print(io, _TREE_BRANCH)
    printstyled(io, "Objective"; bold = true)
    print(io, "   ")
    printstyled(
        io,
        "total = $(_fmt_val(ins.objective_total))  @ current x";
        color = :light_black,
    )
    print(io, "\n")

    if !isempty(ins.objective_terms)
        max_name = maximum(length(t.name) for t in ins.objective_terms)
        for t in ins.objective_terms
            namepad = rpad(t.name, max_name)
            w_part = isfinite(t.weight) ? "w=$(_fmt_val(t.weight))" : ""
            print(
                io,
                _TREE_PIPE,
                "  ",
                namepad,
                "   ",
                rpad(w_part, 16),
                _fmt_val(t.value),
                "\n",
            )
        end
    end
    print(io, _TREE_PIPE, "\n")
end

function _print_constraints(io::IO, ins::ProblemInspection)
    print(io, _TREE_BRANCH)
    printstyled(io, "Constraints"; bold = true)
    feas_count = count(c -> c.feasible, ins.constraints)
    n = length(ins.constraints)
    print(io, "   ")
    if feas_count == n
        printstyled(io, "all feasible at x₀"; color = :green)
    else
        printstyled(io, "$(n - feas_count)/$n violated at x₀"; color = :red)
    end
    print(io, "\n")

    if !isempty(ins.constraints)
        max_name = maximum(length(c.name) for c in ins.constraints)
        for c in ins.constraints
            kind_str =
                c.kind === :dynamics ? "[dyn]" :
                c.kind === :eq ? "[eq]" : c.kind === :ineq ? "[ineq]" : "[bnd]"
            kind_part = rpad(kind_str, 7)
            namepad = rpad(c.name, max_name)
            tick_color = c.feasible ? :green : :red
            tick = c.feasible ? "✓" : "✗"
            v_part = if c.kind === :bnd
                ""
            elseif isnan(c.violation)
                "(no eval)"
            elseif c.kind === :ineq
                "(viol = $(_fmt_val(c.violation)))"
            else
                "(‖c‖∞ = $(_fmt_val(c.violation)))"
            end
            print(io, _TREE_PIPE, "  ", kind_part, namepad, "   ")
            printstyled(io, tick; color = tick_color)
            print(io, "  ", v_part, "\n")
        end
    end
    print(io, _TREE_PIPE, "\n")
end

function _print_status(io::IO, ins::ProblemInspection; last::Bool = true)
    glyph = last ? _TREE_LAST : _TREE_BRANCH
    print(io, glyph)
    printstyled(io, "Status"; bold = true)
    print(io, "\n")
    indent = last ? _TREE_BLANK : _TREE_PIPE
    print(io, indent, "  variables: $(ins.n_vars)   ($(ins.n_bounded_vars) bounded)\n")
    print(io, indent, "  equality:  $(ins.n_eq)\n")
    print(io, indent, "  inequality: $(ins.n_ineq)\n")
    if !isnothing(ins.F_current)
        print(io, indent, "  ")
        printstyled(io, "F (raw)       = "; color = :light_black)
        printstyled(io, @sprintf("%.6f", ins.F_current); bold = true)
        print(io, "\n")
    end
    if !isnothing(ins.F_with_phase)
        print(io, indent, "  ")
        printstyled(io, "F (with φ)    = "; color = :light_black)
        # Bold + green when essentially 1
        F = ins.F_with_phase
        is_target = F >= 0.99
        printstyled(
            io,
            @sprintf("%.6f", F);
            bold = true,
            color = is_target ? :green : :default,
        )
        print(io, "\n")
    end
end

function _print_full_extras(io::IO, ins::ProblemInspection, qcp)
    print(io, "\n")
    printstyled(io, "Pulse (current)\n"; bold = true)
    pulse_plot_str = try
        pulse_lineplot(qcp; height = 12)
    catch e
        "(plot unavailable: $(typeof(e).name.name))"
    end
    println(io, pulse_plot_str)
end

# ---------------------------------------------------------------------------- #
# Small formatters
# ---------------------------------------------------------------------------- #

function _fmt_T(t::Real)
    if isnan(t)
        return "—"
    elseif abs(t) >= 1000
        return @sprintf("%.1f", t)
    elseif abs(t) >= 1
        return @sprintf("%.3f", t)
    else
        return @sprintf("%.3g", t)
    end
end

function _fmt_val(v::Real)
    isnan(v) && return "NaN"
    av = abs(v)
    if av == 0
        return "0"
    elseif av >= 1e4 || (av != 0 && av < 1e-2)
        return @sprintf("%.3e", v)
    else
        return @sprintf("%.4g", v)
    end
end

# ============================================================================ #
# Tests
# ============================================================================ #

@testitem "show.jl small formatters" begin
    using Piccolo

    _fmt_T = Piccolo.Control.ProblemDisplay._fmt_T
    _fmt_val = Piccolo.Control.ProblemDisplay._fmt_val

    # _fmt_T: NaN, large, order-1, sub-1
    @test _fmt_T(NaN) == "—"
    @test _fmt_T(1234.5) == "1234.5"
    @test _fmt_T(12.3456) == "12.346"
    @test _fmt_T(0.0123) == "0.0123"

    # _fmt_val: NaN, exact zero, huge, tiny, ordinary
    @test _fmt_val(NaN) == "NaN"
    @test _fmt_val(0.0) == "0"
    @test _fmt_val(2.5e5) == "2.500e+05"
    @test _fmt_val(3.7e-3) == "3.700e-03"
    @test _fmt_val(0.5) == "0.5"
end

@testitem "inspect.jl bound-rendering helpers" begin
    using Piccolo

    PD = Piccolo.Control.ProblemDisplay

    # _bound_repr: missing key → unbounded dash
    @test PD._bound_repr((x = (-fill(1.0, 2), fill(1.0, 2)),), :u) == ("—", false)

    # all-inf bounds → unbounded dash
    bounds = (u = (fill(-Inf, 2), fill(Inf, 2)),)
    @test PD._bound_repr(bounds, :u) == ("—", false)

    # symmetric vector bounds, short (n ≤ 4): "±[values]"
    r, b = PD._bound_repr((u = ([-1.0, -0.5], [1.0, 0.5]),), :u)
    @test b
    @test occursin("±", r)
    @test occursin("1.0", r)

    # symmetric vector bounds, long (n > 4): summarized with "…"
    lo = -collect(1.0:6.0)
    hi = collect(1.0:6.0)
    r_long, b_long = PD._bound_repr((u = (lo, hi),), :u)
    @test b_long
    @test occursin("…", r_long)
    @test occursin("6 total", r_long)

    # asymmetric vector bounds: "[lo, hi]" both summarized
    r_asym, b_asym = PD._bound_repr((u = (fill(0.0, 6), collect(1.0:6.0)),), :u)
    @test b_asym
    @test occursin("[", r_asym)

    # scalar-ish bounds fall through to "lo..hi"
    r_scalar, b_scalar = PD._bound_repr((Δt = (0.01, 0.5),), :Δt)
    @test b_scalar
    @test occursin("..", r_scalar)

    # _summarize_vec: empty, short, long
    @test PD._summarize_vec(Float64[]) == ""
    @test PD._summarize_vec([1.0, 2.0]) == "1.0, 2.0"
    s = PD._summarize_vec(collect(1.0:7.0))
    @test occursin("…", s) && occursin("7 total", s)

    # _repr_global_bound: 2π recognition (vector form), symmetric scalar,
    # asymmetric scalar, symmetric/asymmetric vectors, non-tuple fallback
    @test PD._repr_global_bound(([-2π], [2π])) == "±2π"
    @test PD._repr_global_bound((-1.5, 1.5)) == "±1.5"
    @test PD._repr_global_bound((-0.5, 1.5)) == "[-0.5, 1.5]"
    @test PD._repr_global_bound(([-1.0], [1.0])) == "±1.0"
    @test occursin("[", PD._repr_global_bound(([-0.5, 0.0], [1.0, 0.5])))
    @test PD._repr_global_bound([0.1, 0.2]) == "[0.1, 0.2]"
end

@testitem "inspect.jl system/trajectory helper branches" begin
    using Piccolo
    using NamedTrajectories

    PD = Piccolo.Control.ProblemDisplay

    # -- _system_dim: levels scalar, H_drift, size-only fallback, nothing --
    sys = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])
    @test PD._system_dim(sys) == 2          # via :levels (a scalar Int)
    struct DriftOnlySys
        H_drift::Matrix{ComplexF64}
    end
    @test PD._system_dim(DriftOnlySys(PAULIS.Z)) == 2
    struct SizedSys end
    Base.size(::SizedSys) = (3, 3)
    Base.size(::SizedSys, d::Integer) = 1 <= d <= 2 ? 3 : 1
    @test PD._system_dim(SizedSys()) == 3
    struct OpaqueSys end
    @test PD._system_dim(OpaqueSys()) == 0

    # -- _subsystem_levels: scalar → nothing, vector → collected, other → nothing --
    @test PD._subsystem_levels(sys) === nothing
    struct LevelsVecSys
        levels::Vector{Int}
    end
    @test PD._subsystem_levels(LevelsVecSys([2, 3])) == [2, 3]
    struct LevelsAnySys
        levels
    end
    @test PD._subsystem_levels(LevelsAnySys("not-levels")) === nothing

    # -- _system_globals: no global components, matching globals, non-scalars --
    traj_plain =
        NamedTrajectory((u = randn(1, 5), Δt = fill(0.1, 5)); timestep = :Δt, controls = :u)
    @test isempty(PD._system_globals(sys, traj_plain))

    traj_glob = NamedTrajectory(
        (u = randn(1, 5), Δt = fill(0.1, 5));
        timestep = :Δt,
        controls = :u,
        global_data = [0.5, 1.0],
        global_components = (δ = 1:1, junk = 2:2),
    )
    struct GlobalsSys
        global_params::NamedTuple
    end
    gsys = GlobalsSys((δ = 0.5, χ = 0.2, M = [1.0 2.0; 3.0 4.0]))
    pairs_out = PD._system_globals(gsys, traj_glob)
    @test pairs_out == [:δ => 0.5]                     # only names shared with traj
    @test isempty(PD._system_globals(sys, traj_glob))  # no global_params on QuantumSystem

    # -- _Δt_range: missing, present, empty --
    @test PD._Δt_range(traj_plain) === nothing
    traj_bnd = NamedTrajectory(
        (u = randn(1, 5), Δt = fill(0.1, 5));
        timestep = :Δt,
        controls = :u,
        bounds = (Δt = ([0.01], [0.5]),),
    )
    @test PD._Δt_range(traj_bnd) == (0.01, 0.5)

    # -- _safe_duration: get_duration works on ordinary trajectories --
    # (N = 5 knots × Δt = 0.1 → 4 intervals)
    @test PD._safe_duration(traj_plain) ≈ 0.4
    struct NoDurationTraj end
    @test isnan(PD._safe_duration(NoDurationTraj()))

    # -- _count_bounded_vars: per-knot bounded components + globals --
    traj_count = NamedTrajectory(
        (x = randn(2, 4), u = randn(1, 4), Δt = fill(0.1, 4));
        timestep = :Δt,
        controls = :u,
        bounds = (x = (-ones(2), ones(2)), u = 1.0),
        global_data = [0.7],
        global_components = (a = 1:1,),
    )
    # x contributes dim 2 × N 4 = 8; u contributes 1 × 4 = 4; globals 1
    @test PD._count_bounded_vars(traj_count) == 13
end

@testitem "inspect.jl goal summaries and integrator rendering" begin
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra

    PD = Piccolo.Control.ProblemDisplay

    sys = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])
    times = collect(range(0.0, 1.0, length = 6))
    pulse = ZeroOrderPulse(0.1 * randn(1, 6), times)

    # Unitary goal summaries: plain matrix vs EmbeddedOperator
    qtraj_u = UnitaryTrajectory(sys, pulse, GATES[:X])
    @test occursin("2×2", PD._goal_summary(qtraj_u))
    op = EmbeddedOperator(GATES[:X], 1:2, [2])
    qtraj_e = UnitaryTrajectory(sys, pulse, op)
    @test occursin("EmbeddedOperator", PD._goal_summary(qtraj_e))
    @test occursin("[2]", PD._goal_summary(qtraj_e))

    # Ket / MultiKet / Density / Sampling summaries
    ψ0, ψ1 = ComplexF64[1.0, 0.0], ComplexF64[0.0, 1.0]
    qtraj_k = KetTrajectory(sys, pulse, ψ0, ψ1)
    @test occursin("ψ_init", PD._goal_summary(qtraj_k))

    qtraj_mk = MultiKetTrajectory(sys, pulse, [ψ0, ψ1], [ψ1, ψ0])
    @test occursin("2 state transfers", PD._goal_summary(qtraj_mk))

    osys = OpenQuantumSystem(sys)
    qtraj_d = DensityTrajectory(osys, pulse, ψ0 * ψ0', ψ1 * ψ1')
    @test PD._goal_summary(qtraj_d) == "ρ_init → ρ_goal"

    qtraj_s = SamplingTrajectory(qtraj_u, [sys, sys])
    @test occursin("sampled ensemble (n=2)", PD._goal_summary(qtraj_s))

    # _render_integrator: state_name / variable / name / bare typename
    struct IntA
        state_name::Symbol
    end
    struct IntB
        variable::Symbol
    end
    struct IntC
        name::Symbol
    end
    struct IntD end
    @test PD._render_integrator(IntA(:ψ̃)) == "IntA(:ψ̃)"
    @test PD._render_integrator(IntB(:u)) == "IntB(:u)"
    @test PD._render_integrator(IntC(:ρ)) == "IntC(:ρ)"
    @test PD._render_integrator(IntD()) == "IntD"
end

@testitem "inspect.jl objective decomposition helpers" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt

    PD = Piccolo.Control.ProblemDisplay

    N = 5
    traj = NamedTrajectory(
        (x = randn(2, N), u = randn(1, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
    )

    # Flat objective: single (obj, 1.0) pair
    J_flat = QuadraticRegularizer(:u, traj, 0.5)
    terms, total = PD._objective_terms(J_flat, traj)
    @test length(terms) == 1
    @test terms[1].weight == 1.0
    @test terms[1].value ≈ objective_value(J_flat, traj)
    @test total ≈ terms[1].value

    # Composite objective: zip(objectives, weights) branch
    J_comp = CompositeObjective([J_flat, QuadraticRegularizer(:x, traj, 0.25)], [2.0, 0.5])
    terms_c, total_c = PD._objective_terms(J_comp, traj)
    @test length(terms_c) == 2
    @test terms_c[1].weight == 2.0
    @test terms_c[2].weight == 0.5
    @test total_c ≈
          2 * terms_c[1].value / 2 + sum(terms_c[t].value for t = 1:2) - terms_c[1].value

    # _objective_label: name symbol, state_name symbol, bare
    @test occursin("(:u)", PD._objective_label(J_flat))
    struct NamedStateObj
        state_name::Symbol
    end
    @test PD._objective_label(NamedStateObj(:ψ̃)) == "NamedStateObj(:ψ̃)"
    struct BareObj end
    @test PD._objective_label(BareObj()) == "BareObj"

    # _try_eval_objective: value, non-number, error
    struct ValObj end
    (::ValObj)(traj) = 3.25
    @test PD._try_eval_objective(ValObj(), traj) == 3.25
    struct NonNumberObj end
    (::NonNumberObj)(traj) = "not a number"
    @test PD._try_eval_objective(NonNumberObj(), traj) == 0.0
    struct ThrowingObj end
    (::ThrowingObj)(traj) = error("boom")
    @test isnan(PD._try_eval_objective(ThrowingObj(), traj))
end

@testitem "inspect.jl constraint classification helpers" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    PD = Piccolo.Control.ProblemDisplay

    N = 6
    traj = NamedTrajectory(
        (ψ̃ = randn(4, N), u = randn(1, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
    )

    # -- _constraint_name: final_fidelity / fidelity property branches --
    struct FinalFidCon
        final_fidelity::Float64
    end
    @test PD._constraint_name(FinalFidCon(0.99), "FinalFidCon", traj) ==
          "FinalFidelity F ≥ 0.99"
    struct FidCon
        fidelity::Float64
    end
    @test PD._constraint_name(FidCon(0.95), "FidCon", traj) == "FinalFidelity F ≥ 0.95"

    # -- NonlinearKnotPointConstraint patterns --
    g_ineq = x -> [1.0 - x[1]]
    terminal_ineq =
        NonlinearKnotPointConstraint(g_ineq, :ψ̃, traj; equality = false, times = [N])
    @test occursin(
        "Terminal fidelity bound",
        PD._constraint_name(terminal_ineq, "NonlinearKnotPointConstraint", traj),
    )

    eq_con = NonlinearKnotPointConstraint(x -> [0.0], [:ψ̃, :u], traj; equality = true)
    @test occursin(
        "Knot-point eq",
        PD._constraint_name(eq_con, "NonlinearKnotPointConstraint", traj),
    )

    ineq_nonterminal = NonlinearKnotPointConstraint(
        x -> [1.0 - x[1]],
        :ψ̃,
        traj;
        equality = false,
        times = 1:(N-2),
    )
    @test occursin(
        "Knot-point ineq",
        PD._constraint_name(ineq_nonterminal, "NonlinearKnotPointConstraint", traj),
    )

    # -- var label extraction: name / var_names / state_name --
    struct NamedCon
        name::Symbol
    end
    @test PD._constraint_name(NamedCon(:my_thing), "NamedCon", traj) == "NamedCon:my_thing"
    struct StateNameCon
        state_name::Symbol
    end
    @test PD._constraint_name(StateNameCon(:ψ̃), "StateNameCon", traj) == "StateNameCon:ψ̃"
    struct BareCon end
    @test PD._constraint_name(BareCon(), "BareCon", traj) == "BareCon"

    # _constraint_var_label through the var_names branch
    @test occursin(":ψ̃", PD._constraint_name(eq_con, "NonlinearKnotPointConstraint", traj))

    # -- _constraint_kind --
    bnd = BoundsConstraint(:u, collect(1:N), 1.0)
    @test PD._constraint_kind(bnd, "BoundsConstraint") == :bnd
    @test PD._constraint_kind(1, "GlobalEqualityConstraint") == :eq
    @test PD._constraint_kind(1, "TimeStepsAllEqualConstraint") == :eq
    @test PD._constraint_kind(1, "TimeConsistencyConstraint") == :eq
    @test PD._constraint_kind(terminal_ineq, "NonlinearKnotPointConstraint") == :ineq
    @test PD._constraint_kind(eq_con, "NonlinearKnotPointConstraint") == :eq
    @test PD._constraint_kind(1, "SomethingFidelity") == :ineq
    @test PD._constraint_kind(1, "SomethingLeakage") == :ineq
    @test PD._constraint_kind(1, "MysteryConstraint") == :ineq

    # -- _constraint_dim: g_dim / dim / n_constraints / fallback --
    @test PD._constraint_dim(terminal_ineq) == 1
    struct DimCon
        dim::Int
    end
    @test PD._constraint_dim(DimCon(7)) == 7
    struct NCon
        n_constraints::Int
    end
    @test PD._constraint_dim(NCon(4)) == 4
    @test PD._constraint_dim(BareCon()) == 1

    # -- _try_eval_constraint: bounds → 0.0, real eval → max|g|, error → NaN --
    @test PD._try_eval_constraint(bnd, traj, :bnd) == 0.0
    v = PD._try_eval_constraint(terminal_ineq, traj, :ineq)
    @test v isa Float64
    @test PD._try_eval_constraint(BareCon(), traj, :ineq) |> isnan

    # -- _classify_and_eval end-to-end on a real constraint --
    ci = PD._classify_and_eval(terminal_ineq, traj)
    @test ci.kind == :ineq
    @test ci.dim == 1
    @test ci.violation == v
    @test ci.feasible == (v < 1e-6)

    # -- _integrator_dim: state_dim / dim / n / fallback --
    struct SDimInt
        state_dim::Int
    end
    struct NInt
        n::Int
    end
    @test PD._integrator_dim(SDimInt(4), traj) == 4 * (N - 1)
    @test PD._integrator_dim(DimCon(3), traj) == 3 * (N - 1)
    @test PD._integrator_dim(NInt(2), traj) == 2 * (N - 1)
    @test PD._integrator_dim(BareCon(), traj) == N - 1
end

@testitem "inspect.jl free-phase fidelity helpers" begin
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra

    PD = Piccolo.Control.ProblemDisplay

    sys = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])
    times = collect(range(0.0, 1.0, length = 5))
    pulse = ZeroOrderPulse(0.1 * randn(1, 5), times)

    # _fidelity_at on a UnitaryTrajectory with a plain-matrix goal: nothing
    qtraj_u = UnitaryTrajectory(sys, pulse, GATES[:X])
    traj = NamedTrajectory(
        (Ũ⃗ = randn(8, 5), u = randn(1, 5), Δt = fill(0.25, 5));
        controls = :u,
        timestep = :Δt,
        global_data = [0.3],
        global_components = (φ_1 = 1:1,),
    )
    @test PD._fidelity_at(qtraj_u, traj, [0.3]) === nothing

    # _fidelity_at on a KetTrajectory: QuantumSystem.levels is a scalar Int, so
    # _phased_ket_goal has no matching method — the MethodError is caught by
    # _fidelity_with_stored_phases (which then reports F_with_phase = nothing).
    qtraj_k = KetTrajectory(sys, pulse, ComplexF64[1.0, 0.0], ComplexF64[0.0, 1.0])
    traj_ket = NamedTrajectory(
        (ψ̃ = randn(4, 5), u = randn(1, 5), Δt = fill(0.25, 5));
        controls = :u,
        timestep = :Δt,
        global_data = [0.3],
        global_components = (φ_1 = 1:1,),
    )
    @test_throws MethodError PD._fidelity_at(qtraj_k, traj_ket, [0.3])
    @test PD._fidelity_with_stored_phases(qtraj_k, traj_ket) === nothing

    # _fidelity_with_stored_phases: no φ_ globals → nothing
    traj_nophi = NamedTrajectory(
        (ψ̃ = randn(4, 5), u = randn(1, 5), Δt = fill(0.25, 5));
        controls = :u,
        timestep = :Δt,
        global_data = [0.7],
        global_components = (δ = 1:1,),
    )
    @test PD._fidelity_with_stored_phases(qtraj_k, traj_nophi) === nothing

    # _phased_ket_goal: matching θ applies e^{iθ} to |1⟩
    ψ = ComplexF64[1.0, 1.0] / √2
    phased = PD._phased_ket_goal(ψ, [π], [2])
    @test phased ≈ ComplexF64[1.0, -1.0] / √2
    # Length mismatch: goal returned unchanged
    @test PD._phased_ket_goal(ψ, [π, π], [2]) === ψ
    # Multi-level number-operator rotation: level s gets e^{isθ}
    ψ3 = ComplexF64[0.0, 1.0, 0.0]
    phased3 = PD._phased_ket_goal(ψ3, [π / 2], [3])
    @test phased3 ≈ ComplexF64[0.0, im, 0.0]

    # _initial_fidelities: F computed, F_phase nothing without φ_ globals
    F, F_phase = PD._initial_fidelities(qtraj_k, traj_nophi)
    @test F isa Float64
    @test 0.0 <= F <= 1.0
    @test F_phase === nothing

    # _typename strips type parameters
    @test PD._typename([1.0, 2.0]) == "Array"
    @test PD._typename(GATES[:X]) == "Array"
end

@testitem "inspect + show_problem on a solved-shape SmoothPulseProblem" begin
    using Piccolo
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra
    using Random

    Random.seed!(11)
    T, N = 5.0, 20
    sys = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    pulse = ZeroOrderPulse(0.1 * randn(2, N), collect(range(0.0, T, length = N)))
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

    qcp = SmoothPulseProblem(
        qtraj,
        N;
        Q = 100.0,
        R = 1e-2,
        piccolo_options = PiccoloOptions(display = :silent),
    )

    ins = inspect(qcp)
    @test ins isa ProblemInspection
    @test ins.traj_typename == "KetTrajectory"
    @test ins.pulse_typename == "ZeroOrderPulse"
    @test ins.sys_dim == 2
    @test ins.n_drives == 2
    @test ins.subsystem_levels === nothing
    @test ins.N == N
    @test ins.T ≈ T
    @test ins.goal_summary isa String
    @test !isempty(ins.components)
    @test !isempty(ins.integrator_summary)
    @test length(ins.constraints) >= 3  # dynamics + derivative integrators
    @test ins.F_current isa Float64
    @test ins.F_with_phase === nothing   # no φ_ globals on this problem
    @test ins.n_vars > 0
    @test ins.n_eq > 0
    @test ins.n_bounded_vars > 0

    # Compact one-line show: Base.show(io, qcp)
    s = sprint(show, qcp)
    @test occursin("QuantumControlProblem", s)
    @test occursin("KetTrajectory", s)
    @test occursin("ZeroOrderPulse", s)
    @test occursin("vars", s)
    @test occursin("eq", s)

    # Rich text/plain show → :standard tree
    s_std = sprint(show, MIME"text/plain"(), qcp)
    @test occursin("System", s_std)
    @test occursin("Trajectory", s_std)
    @test occursin("Goal", s_std)
    @test occursin("Objective", s_std)
    @test occursin("Constraints", s_std)
    @test occursin("Status", s_std)
    @test occursin("Hint:", s_std)      # drill-down hint at :standard level

    # Explicit entry points: :standard, :full, and the ArgumentError guard
    io = IOBuffer()
    show_problem(io, qcp; detail = :standard)
    @test occursin("QuantumControlProblem", String(take!(io)))
    io = IOBuffer()
    show_problem(io, qcp; detail = :full)
    s_full = String(take!(io))
    @test occursin("Pulse (current)", s_full)  # _print_full_extras ran
    @test occursin("time", s_full)             # ASCII pulse plot rendered
    @test !occursin("Hint:", s_full)           # hint replaced by extras at :full
    @test_throws ArgumentError show_problem(IOBuffer(), qcp; detail = :bogus)

    # stdout convenience method runs without error
    show_problem(qcp; detail = :standard)
end

@testitem "inspect reports phase-adjusted fidelity for free-phase problems" begin
    using Piccolo
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    # Free-phase unitary problem: the trajectory carries φ_1, so the inspector
    # must compute F_with_phase via the phased EmbeddedOperator goal.
    σx = ComplexF64[0 1; 1 0]
    H_drift_3 = ComplexF64[0 0 0; 0 1 0; 0 0 2]
    H_drive_3 = ComplexF64[0 1 0; 1 0 1; 0 1 0] / √2
    sys = QuantumSystem(H_drift_3, [H_drive_3], [1.0])

    T, N = 10.0, 21
    times = collect(range(0.0, T, length = N))
    pulse = LinearSplinePulse(0.1 * ones(1, N), times)
    U_goal = EmbeddedOperator(σx, [1, 2], [3])
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)

    qcp = SplinePulseProblem(
        qtraj,
        N;
        Q = 100.0,
        R = 1e-2,
        free_phase = true,
        piccolo_options = PiccoloOptions(display = :silent),
    )

    ins = inspect(qcp)
    @test haskey(get_trajectory(qcp).global_components, :φ_1)
    @test ins.F_current isa Float64
    @test ins.F_with_phase isa Float64
    @test 0.0 <= ins.F_with_phase <= 1.0
    # Before any solve, φ_1 = 0, so the phase-adjusted F equals the raw F
    @test ins.F_with_phase ≈ ins.F_current

    # The Globals section renders the φ_1 row with :free_phase status
    s = sprint(show, MIME"text/plain"(), qcp)
    @test occursin("φ_1", s)
    @test occursin("free_phase", s)

    # Global bounds on φ_1 surface as a ±2π row in the inspection
    @test any(g -> g.name === :φ_1, ins.globals)
end

@testitem "_maybe_display respects the display tier" begin
    using Piccolo
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra
    using Random

    Random.seed!(5)
    T, N = 5.0, 10
    sys = QuantumSystem(GATES[:Z], [GATES[:X]], [1.0])
    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    qtraj = KetTrajectory(sys, pulse, ComplexF64[1.0, 0.0], ComplexF64[0.0, 1.0])

    # :detailed runs show_problem(detail=:full) — the pulse plot + sparsity view
    qcp_detailed =
        SmoothPulseProblem(qtraj, N; piccolo_options = PiccoloOptions(display = :detailed))
    @test qcp_detailed isa QuantumControlProblem

    # :silent constructs with no output at all
    qcp_silent =
        SmoothPulseProblem(qtraj, N; piccolo_options = PiccoloOptions(display = :silent))
    @test qcp_silent isa QuantumControlProblem
end
