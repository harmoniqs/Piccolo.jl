# ============================================================================ #
# Rich display rendering for QuantumControlProblem
# ============================================================================ #
#
# `Base.show(io, ::MIME"text/plain", qcp)` → :standard level (tree view)
# `Base.show(io, qcp)`                     → :compact one-line summary
# `show_problem(io, qcp; detail)`          → explicit entry point

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

# ---------------------------------------------------------------------------- #
# Mutating helpers that keep the displayed view fresh
# ---------------------------------------------------------------------------- #
#
# The auto-display fires once at the end of each outer-template constructor.
# If you mutate the problem afterwards (add an objective, push a constraint),
# the previously-printed view goes stale. These helpers do the mutation AND
# redisplay in one call so the log shows the state the solver actually sees.

"""
    add_objective!(qcp, obj; redisplay=true, detail=:standard) -> qcp

Append `obj` to the inner `prob.objective` (via `+`) and, by default, redisplay
the full problem view so the augmented objective appears in the log.

Pass `redisplay=false` for a silent mutation (useful inside batch updates),
or `detail=:full` to also print the pulse plot after augmenting.
"""
function add_objective!(
    qcp::QuantumControlProblem,
    obj;
    redisplay::Bool = true,
    detail::Symbol = :standard,
)
    qcp.prob.objective = qcp.prob.objective + obj
    if redisplay
        show_problem(stdout, qcp; detail = detail)
        println()
    end
    return qcp
end

"""
    add_constraint!(qcp, c; redisplay=true, detail=:standard) -> qcp

Push `c` onto `prob.constraints` and, by default, redisplay the full problem
view so the new constraint appears in the log. See [`add_objective!`](@ref) for
the `redisplay` and `detail` knobs.
"""
function add_constraint!(
    qcp::QuantumControlProblem,
    c;
    redisplay::Bool = true,
    detail::Symbol = :standard,
)
    push!(qcp.prob.constraints, c)
    if redisplay
        show_problem(stdout, qcp; detail = detail)
        println()
    end
    return qcp
end
