# ============================================================================ #
# Terminal pulse plotting via UnicodePlots
# ============================================================================ #
#
# Rendered as a string so it can be embedded in other display contexts (e.g.
# `show_problem(io, qcp; detail=:full)`).

"""
    pulse_lineplot(qcp::QuantumControlProblem; height=12, width=72) -> String

Render the current control pulse `u(t)` from a `QuantumControlProblem` as an
ASCII line plot via `UnicodePlots`. Returns the plot as a string suitable for
embedding in `show` output.

Each drive channel is one labeled line. For a `CubicSplinePulse`, the plot uses
the knot values directly; tangent (`du`) information is not separately rendered.
"""
function pulse_lineplot(qcp::QuantumControlProblem; height::Int = 8, width::Int = 64)
    traj = qcp.prob.trajectory
    qtraj = qcp.qtraj
    control_sym = drive_name(qtraj)

    u = traj[control_sym]                # (n_drives, N)
    times = collect(get_times(traj))     # length N

    n_drives = size(u, 1)
    n_drives == 0 && return "(no drives)"

    # Build first line, then mutate
    p = lineplot(
        times,
        view(u, 1, :);
        name = _drive_label(control_sym, 1, n_drives),
        height = height,
        width = width,
        xlabel = "time",
        ylabel = "u",
    )
    for k = 2:n_drives
        lineplot!(p, times, view(u, k, :); name = _drive_label(control_sym, k, n_drives))
    end
    return string(p)
end

function _drive_label(control_sym::Symbol, k::Int, n::Int)
    n == 1 && return string(control_sym)
    return string(control_sym, "[", k, "]")
end
