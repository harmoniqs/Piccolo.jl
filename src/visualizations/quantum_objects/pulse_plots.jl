export plot_pulse, plot_pulse!, plot_pulse_IQ, plot_pulse_phases

using LaTeXStrings
using Makie
using NamedTrajectories: NamedTrajectory, get_times
using Piccolo:
    AbstractPulse,
    AbstractSplinePulse,
    AbstractQuantumSystem,
    AbstractQuantumTrajectory,
    QuantumControlProblem,
    ZeroOrderPulse,
    LinearSplinePulse,
    CubicSplinePulse,
    GaussianPulse,
    ErfPulse,
    CompositePulse,
    FunctionPulse,
    duration,
    n_drives,
    sample,
    get_knot_times,
    get_knot_values,
    get_knot_derivatives,
    drive_name,
    get_pulse,
    get_system,
    get_trajectory,
    extract_pulse

using TestItems

# ============================================================================ #
# Theme-aware color helpers
# ============================================================================ #
#
# We don't hardcode dark or light visual choices. Instead the plot reads the
# active Makie theme (via `Makie.current_default_theme()`), so calls like
# `set_theme!(theme_dark())` automatically get a contrasting neutral for
# strokes / zero-lines and a palette-driven color cycle for the drive lines.

const _DEFAULT_PULSE_COLORS = Makie.wong_colors()

function _theme_palette()
    try
        theme = Makie.current_default_theme()
        if haskey(theme, :palette)
            pal = theme[:palette]
            if pal isa Makie.Attributes && haskey(pal, :color)
                colors = Makie.to_value(pal[:color])
                if colors isa AbstractVector && !isempty(colors)
                    return colors
                end
            end
        end
    catch
    end
    return _DEFAULT_PULSE_COLORS
end

function _drive_color(i::Int)
    pal = _theme_palette()
    return pal[mod1(i, length(pal))]
end

# Neutral color that contrasts with the current theme background (text color).
function _theme_neutral()
    try
        theme = Makie.current_default_theme()
        if haskey(theme, :textcolor)
            return Makie.to_value(theme[:textcolor])
        end
    catch
    end
    return :black
end

# ============================================================================ #
# LaTeX label helpers — render labels in math font for publication quality.
# Mirrors NamedTrajectories.jl's `latexstring(b, "_{", i, "}")` pattern: italic
# math letters with proper subscripts. Allows users to pass `LaTeXString`s
# without double-`$`-wrapping, and keeps plain strings working.
# ============================================================================ #

# Pass LaTeXString through unchanged; wrap any other string as math.
_to_latex_label(s::LaTeXString) = s
_to_latex_label(s::AbstractString) = latexstring(s)

# Render a NamedTrajectory component name as a math label, recognizing the
# `:du`, `:ddu`, `:dddu` derivative naming convention used across Piccolo
# problem templates (`smooth_pulse_problem.jl`, etc.) so derivatives render as
# `\dot{u}`, `\ddot{u}`, `\dddot{u}` instead of literally `du`.
function _component_latex(name::Symbol)
    s = String(name)
    if startswith(s, "ddd") && length(s) > 3
        return latexstring("\\dddot{", s[4:end], "}")
    elseif startswith(s, "dd") && length(s) > 2
        return latexstring("\\ddot{", s[3:end], "}")
    elseif startswith(s, "d") && length(s) > 1 && !isdigit(s[2])
        return latexstring("\\dot{", s[2:end], "}")
    end
    return latexstring(s)
end

# ============================================================================ #
# Main entry point: plot_pulse
# ============================================================================ #

"""
    plot_pulse(
        pulse::AbstractPulse;
        n_samples=500,
        labels=nothing,
        title="",
        layout=:stacked,
        show_knots=true,
        show_tangents=false,
        tangent_scale=0.1,
        bounds=nothing,
        figsize=nothing,
        kwargs...
    )

Plot a pulse with type-appropriate rendering.

Each pulse type is rendered to be visually faithful to its actual interpolation:
- `ZeroOrderPulse`: Step function via `stairs!`
- `LinearSplinePulse`: Line segments connecting knot values
- `CubicSplinePulse`: Dense-sampled smooth curve with knot markers
- `GaussianPulse`, `ErfPulse`, `CompositePulse`, `FunctionPulse`: Dense-sampled smooth curve

# Theming

`plot_pulse` honors the active Makie theme. To render every pulse plot in dark
mode (e.g. for a dark documentation build), set the theme once at the top of
your script:

```julia
using CairoMakie
set_theme!(theme_dark())
```

Drive line colors default to the theme's `:palette[:color]` cycle (falling back
to the colorblind-safe Wong palette). Knot strokes, zero-lines, and bounds use
the theme's `:textcolor`, so they remain visible on either background.

# Keyword Arguments
- `n_samples::Int=500`: Number of time samples for smooth curves (cubic/analytic).
  Ignored for ZeroOrderPulse and LinearSplinePulse since the exact shape is known.
- `labels::Union{Nothing,Vector{String}}=nothing`: Labels for each drive channel.
  Defaults to `"Drive 1"`, `"Drive 2"`, etc.
- `title::AbstractString=""`: Plot super-title.
- `layout::Symbol=:stacked`: Plot layout.
  - `:stacked` — one vertically stacked subplot per drive (recommended).
  - `:overlay` — all drives on a single shared axis with legend.
- `show_knots::Bool=true`: Show knot points as markers for spline pulses.
- `show_tangents::Bool=false`: Show tangent whiskers at knot points (CubicSplinePulse only).
- `tangent_scale::Float64=0.1`: Scale factor for tangent whisker length (fraction of duration).
- `bounds::Union{Nothing,AbstractVector}=nothing`: Per-drive hardware bounds as
  `(lower, upper)` tuples. Draws shaded band and dashed limits (stacked layout only).
- `figsize::Union{Nothing,Tuple}=nothing`: Figure size. Defaults to `(800, 400)` for overlay
  or `(800, 100 + 140 * n_drives)` for stacked.
- `kwargs...`: Passed to `Figure(; kwargs...)`.

# Returns
A Makie `Figure`.
"""
function plot_pulse(
    pulse::AbstractPulse;
    n_samples::Int = 500,
    labels::Union{Nothing,AbstractVector{<:AbstractString}} = nothing,
    title::AbstractString = "",
    layout::Symbol = :stacked,
    show_knots::Bool = true,
    show_tangents::Bool = false,
    tangent_scale::Float64 = 0.1,
    bounds::Union{Nothing,AbstractVector} = nothing,
    figsize::Union{Nothing,Tuple} = nothing,
    kwargs...,
)
    @assert layout in (:stacked, :overlay) "layout must be :stacked or :overlay"

    nd = n_drives(pulse)

    # Always end up with a Vector{LaTeXString} so axis ylabels and legend
    # entries render in math font with proper italics + subscripts.
    if isnothing(labels)
        labels = _default_drive_labels(pulse)
    else
        labels = [_to_latex_label(s) for s in labels]
    end

    if layout == :stacked
        sz = isnothing(figsize) ? (800, 100 + 140 * nd) : figsize
        fig = Figure(; size = sz, kwargs...)

        for i = 1:nd
            ax = Axis(
                fig[i, 1];
                titlealign = :left,
                titlesize = 18,
                titlefont = :bold,
                title = (i == 1 && !isempty(title)) ? title : "",
                xticklabelsvisible = i == nd,
                xtickalign = 1,
                xlabel = i == nd ? "Time" : "",
                xlabelsize = 16,
                ylabel = labels[min(i, length(labels))],
                ylabelsize = 16,
            )

            bound_pair = _drive_bound_pair(bounds, i)
            _draw_pulse_subplot!(
                ax,
                pulse,
                i;
                bound_pair,
                n_samples,
                show_knots,
                show_tangents,
                tangent_scale,
            )
        end

        # Tight rowgap between stacked panels — same pattern as
        # NamedTrajectories.jl's `plot(::NamedTrajectory)`.
        for i = 1:(nd-1)
            rowgap!(fig.layout, i, Relative(0.015))
        end

        return fig

    else  # :overlay
        sz = isnothing(figsize) ? (800, 400) : figsize
        fig = Figure(; size = sz, kwargs...)
        ax = Axis(
            fig[1, 1];
            xlabel = "Time",
            xlabelsize = 16,
            ylabel = "Amplitude",
            ylabelsize = 16,
            titlealign = :left,
            titlesize = 18,
            title = isempty(title) ? "Pulse Controls" : title,
        )

        colors = [_drive_color(i) for i = 1:nd]
        plot_pulse!(ax, pulse; n_samples, show_knots, show_tangents, tangent_scale, colors)

        # Legend
        entries = [LineElement(; color = colors[i], linewidth = 2) for i = 1:nd]
        Legend(fig[1, 2], entries, labels; tellheight = false)

        return fig
    end
end

# ============================================================================ #
# plot_pulse! — draw onto an existing Axis
# ============================================================================ #

"""
    plot_pulse!(
        ax, pulse::AbstractPulse;
        drive_indices=1:n_drives(pulse),
        n_samples=500,
        show_knots=true,
        show_tangents=false,
        tangent_scale=0.1,
        colors=nothing,
        t_max=nothing,
        kwargs...
    )

Plot a pulse's control channels onto an existing Makie `Axis`.

Dispatches rendering based on pulse type for visually accurate results.

# Keyword Arguments
- `drive_indices`: Which drives to plot (default: all).
- `n_samples`: Sampling density for smooth curves.
- `show_knots`: Overlay knot markers for interpolated pulses.
- `show_tangents`: Show derivative tangent whiskers (CubicSplinePulse only).
- `tangent_scale`: Length scale for tangent whiskers (fraction of duration).
- `colors`: Vector of colors, one per drive index. Defaults to the active theme palette.
- `t_max`: Optional `Observable{<:Real}` playhead time. When `nothing` (default),
  the full static pulse is drawn. When set, each channel is drawn as an
  observable-driven curve revealed up to `t_max[]`, using the same per-type
  primitive (stairs for `ZeroOrderPulse`, knot-to-knot lines for
  `LinearSplinePulse`, dense curve for cubic/analytic). This is what
  `animate_pulse` drives; static knot/tangent overlays are suppressed while
  revealing. Backward-compatible: omit it for unchanged behavior.
- `kwargs...`: Forwarded to the underlying Makie plot calls.
"""
function plot_pulse!(
    ax,
    pulse::AbstractPulse;
    drive_indices::AbstractVector{Int} = collect(1:n_drives(pulse)),
    n_samples::Int = 500,
    show_knots::Bool = true,
    show_tangents::Bool = false,
    tangent_scale::Float64 = 0.1,
    colors::Union{Nothing,AbstractVector} = nothing,
    t_max::Union{Nothing,Observable} = nothing,
    kwargs...,
)
    if isnothing(colors)
        colors = [_drive_color(i) for i in drive_indices]
    end

    _plot_pulse_type!(
        ax,
        pulse;
        drive_indices,
        n_samples,
        show_knots,
        show_tangents,
        tangent_scale,
        colors,
        t_max,
        kwargs...,
    )
end

# ============================================================================ #
# Type-specific rendering dispatches
# ============================================================================ #

# Animation reveal helper. Returns the (xs, ys) for the portion of one drive of
# `pulse` visible at playhead time `t_max`. A single helper works for every pulse
# type because `pulse(t)` already dispatches to the correct per-type evaluation
# (held value for ZOH, linear/cubic interpolation, analytic form for Gaussian/Erf).
# `render_times` is the type-appropriate base grid (knot times for step/linear,
# dense grid for cubic/analytic); the synthetic endpoint at `t_max` makes the
# revealed curve terminate exactly at the playhead instead of at the last grid
# point behind it. `render_times` must be sorted ascending.
function _reveal_xy(pulse, render_times, drive_index::Int, t_max::Real)
    tmax = float(t_max)
    k = searchsortedlast(render_times, tmax)
    xs = k == 0 ? [tmax] : vcat(render_times[1:k], tmax)
    ys = [pulse(t)[drive_index] for t in xs]
    return xs, ys
end

# ---------------------------------------------------------------------------- #
# ZeroOrderPulse — stairs (step function)
# ---------------------------------------------------------------------------- #

function _plot_pulse_type!(
    ax,
    pulse::ZeroOrderPulse;
    drive_indices,
    n_samples,  # unused — exact shape from knots
    show_knots,
    show_tangents,  # unused
    tangent_scale,  # unused
    colors,
    t_max = nothing,
    kwargs...,
)
    knot_times = collect(get_knot_times(pulse))
    knot_vals = sample(pulse, knot_times)

    for (ci, i) in enumerate(drive_indices)
        if isnothing(t_max)
            stairs!(
                ax,
                knot_times,
                knot_vals[i, :];
                color = colors[ci],
                linewidth = 2,
                step = :post,
                kwargs...,
            )
        else
            xs = @lift(_reveal_xy(pulse, knot_times, i, $t_max)[1])
            ys = @lift(_reveal_xy(pulse, knot_times, i, $t_max)[2])
            stairs!(ax, xs, ys; color = colors[ci], linewidth = 2, step = :post, kwargs...)
        end

        # Knot markers are static structure; suppress during a reveal so future
        # knots don't appear ahead of the playhead.
        if show_knots && isnothing(t_max)
            scatter!(
                ax,
                knot_times,
                knot_vals[i, :];
                color = colors[ci],
                markersize = 6,
                strokewidth = 1,
                strokecolor = _theme_neutral(),
            )
        end
    end
end

# ---------------------------------------------------------------------------- #
# LinearSplinePulse — line segments through knots (exact shape)
# ---------------------------------------------------------------------------- #

function _plot_pulse_type!(
    ax,
    pulse::LinearSplinePulse;
    drive_indices,
    n_samples,  # unused — knot-to-knot lines ARE the exact shape
    show_knots,
    show_tangents,  # unused
    tangent_scale,  # unused
    colors,
    t_max = nothing,
    kwargs...,
)
    knot_times = collect(get_knot_times(pulse))
    knot_vals = get_knot_values(pulse)

    for (ci, i) in enumerate(drive_indices)
        if isnothing(t_max)
            lines!(
                ax,
                knot_times,
                collect(knot_vals[i, :]);
                color = colors[ci],
                linewidth = 2,
                kwargs...,
            )
        else
            xs = @lift(_reveal_xy(pulse, knot_times, i, $t_max)[1])
            ys = @lift(_reveal_xy(pulse, knot_times, i, $t_max)[2])
            lines!(ax, xs, ys; color = colors[ci], linewidth = 2, kwargs...)
        end

        if show_knots && isnothing(t_max)
            scatter!(
                ax,
                knot_times,
                collect(knot_vals[i, :]);
                color = colors[ci],
                markersize = 6,
                strokewidth = 1,
                strokecolor = _theme_neutral(),
            )
        end
    end
end

# ---------------------------------------------------------------------------- #
# CubicSplinePulse — dense-sampled smooth curve + knots + optional tangents
# ---------------------------------------------------------------------------- #

function _plot_pulse_type!(
    ax,
    pulse::CubicSplinePulse;
    drive_indices,
    n_samples,
    show_knots,
    show_tangents,
    tangent_scale,
    colors,
    t_max = nothing,
    kwargs...,
)
    # Dense sampling for smooth curve
    controls, times = sample(pulse, n_samples)

    for (ci, i) in enumerate(drive_indices)
        if isnothing(t_max)
            lines!(ax, times, controls[i, :]; color = colors[ci], linewidth = 2, kwargs...)
        else
            xs = @lift(_reveal_xy(pulse, times, i, $t_max)[1])
            ys = @lift(_reveal_xy(pulse, times, i, $t_max)[2])
            lines!(ax, xs, ys; color = colors[ci], linewidth = 2, kwargs...)
        end
    end

    # Static knot markers and tangent whiskers are suppressed during a reveal
    # so future structure doesn't appear ahead of the playhead.
    if !isnothing(t_max)
        return
    end

    # Knot markers
    if show_knots
        knot_times = collect(get_knot_times(pulse))
        knot_vals = get_knot_values(pulse)

        for (ci, i) in enumerate(drive_indices)
            scatter!(
                ax,
                knot_times,
                collect(knot_vals[i, :]);
                color = colors[ci],
                markersize = 6,
                strokewidth = 1,
                strokecolor = _theme_neutral(),
            )
        end
    end

    # Tangent whiskers
    if show_tangents
        knot_times = collect(get_knot_times(pulse))
        knot_vals = get_knot_values(pulse)
        knot_derivs = get_knot_derivatives(pulse)
        dt = tangent_scale * duration(pulse)

        for (ci, i) in enumerate(drive_indices)
            for k in eachindex(knot_times)
                t_k = knot_times[k]
                v_k = knot_vals[i, k]
                dv_k = knot_derivs[i, k]

                # Draw short line segment centered on knot
                t_lo = t_k - dt / 2
                t_hi = t_k + dt / 2
                v_lo = v_k - dv_k * dt / 2
                v_hi = v_k + dv_k * dt / 2

                lines!(
                    ax,
                    [t_lo, t_hi],
                    [v_lo, v_hi];
                    color = (colors[ci], 0.5),
                    linewidth = 1.5,
                    linestyle = :dash,
                )
            end
        end
    end
end

# ---------------------------------------------------------------------------- #
# Fallback for analytic pulses (Gaussian, Erf, Composite, FunctionPulse, …)
# ---------------------------------------------------------------------------- #

function _plot_pulse_type!(
    ax,
    pulse::AbstractPulse;
    drive_indices,
    n_samples,
    show_knots,  # no knots for analytic pulses
    show_tangents,  # unused
    tangent_scale,  # unused
    colors,
    t_max = nothing,
    kwargs...,
)
    controls, times = sample(pulse, n_samples)

    for (ci, i) in enumerate(drive_indices)
        if isnothing(t_max)
            lines!(ax, times, controls[i, :]; color = colors[ci], linewidth = 2, kwargs...)
        else
            xs = @lift(_reveal_xy(pulse, times, i, $t_max)[1])
            ys = @lift(_reveal_xy(pulse, times, i, $t_max)[2])
            lines!(ax, xs, ys; color = colors[ci], linewidth = 2, kwargs...)
        end
    end
end

# ============================================================================ #
# Subplot helpers — used by both `plot_pulse(::AbstractPulse)`'s stacked layout
# and the high-level `plot_pulse(::QuantumControlProblem)` overload that mixes
# pulse rows with NamedTrajectory component rows. Centralizing the per-axis
# rendering keeps the two paths in lockstep.
# ============================================================================ #

# Pull the (lo, hi) tuple for drive `i` out of a `bounds` argument, or `nothing`
# if not supplied / out of range.
function _drive_bound_pair(bounds, i::Int)
    if isnothing(bounds) || i > length(bounds) || isnothing(bounds[i])
        return nothing
    end
    return bounds[i]
end

# Render a single pulse drive into `ax`: hardware bound band (if any),
# type-aware pulse line via `plot_pulse!`, and a faint zero reference.
function _draw_pulse_subplot!(
    ax,
    pulse::AbstractPulse,
    drive_index::Int;
    bound_pair = nothing,
    n_samples::Int = 500,
    show_knots::Bool = true,
    show_tangents::Bool = false,
    tangent_scale::Float64 = 0.1,
)
    if !isnothing(bound_pair)
        _draw_bounds!(ax, pulse, bound_pair)
    end
    plot_pulse!(
        ax,
        pulse;
        drive_indices = [drive_index],
        n_samples,
        show_knots,
        show_tangents,
        tangent_scale,
        colors = [_drive_color(drive_index)],
    )
    hlines!(ax, [0.0]; color = (_theme_neutral(), 0.4), linestyle = :dash, linewidth = 0.5)
end

# Render a NamedTrajectory component into `ax`: optional shared (lo, hi) band,
# all dimensions overlaid as lines, and a faint zero reference.
function _draw_component_subplot!(
    ax,
    traj::NamedTrajectory,
    name::Symbol;
    bound_pair = nothing,
)
    times = collect(get_times(traj))

    if !isnothing(bound_pair)
        lo, hi = bound_pair
        t0, t1 = times[1], times[end]
        neutral = _theme_neutral()
        band!(
            ax,
            [t0, t1],
            [Float64(lo), Float64(lo)],
            [Float64(hi), Float64(hi)];
            color = (neutral, 0.08),
        )
        hlines!(ax, [lo, hi]; color = (neutral, 0.5), linestyle = :dash, linewidth = 0.8)
    end

    data = traj[name]
    for i = 1:size(data, 1)
        lines!(ax, times, collect(data[i, :]); color = _drive_color(i), linewidth = 2)
    end
    hlines!(ax, [0.0]; color = (_theme_neutral(), 0.4), linestyle = :dash, linewidth = 0.5)
end

# ============================================================================ #
# Helper: draw hardware bounds
# ============================================================================ #

function _draw_bounds!(ax, pulse::AbstractPulse, bound_pair)
    lo, hi = bound_pair
    t_start = 0.0
    t_end = duration(pulse)
    neutral = _theme_neutral()
    band!(
        ax,
        [t_start, t_end],
        [Float64(lo), Float64(lo)],
        [Float64(hi), Float64(hi)];
        color = (neutral, 0.08),
    )
    hlines!(ax, [lo, hi]; color = (neutral, 0.5), linestyle = :dash, linewidth = 0.8)
end

# ============================================================================ #
# IQ-pair helpers (4-drive: Ω_I, Ω_Q, α_I, α_Q)
# ============================================================================ #

"""
    plot_pulse_IQ(pulse::AbstractPulse; n_samples=500, drive_labels=nothing, title=nothing, kwargs...)

Plot a 4-drive pulse as two IQ pairs: drives (Ω_I, Ω_Q) and displacement (α_I, α_Q),
with magnitude envelopes.

Returns a Makie `Figure` with 2 rows: drive IQ and displacement IQ.
"""
function plot_pulse_IQ(
    pulse::AbstractPulse;
    n_samples::Int = 500,
    title::Union{Nothing,String} = nothing,
    figsize::Tuple = (800, 500),
    show_knots::Bool = true,
)
    @assert n_drives(pulse) == 4 "plot_pulse_IQ requires exactly 4 drives (Ω_I, Ω_Q, α_I, α_Q)"

    controls, times = sample(pulse, n_samples)

    Ω_I, Ω_Q = controls[1, :], controls[2, :]
    α_I, α_Q = controls[3, :], controls[4, :]
    Ω_mag = sqrt.(Ω_I .^ 2 .+ Ω_Q .^ 2)
    α_mag = sqrt.(α_I .^ 2 .+ α_Q .^ 2)

    c_I = _drive_color(1)
    c_Q = _drive_color(2)
    c_mag = _theme_neutral()

    fig = Figure(size = figsize)
    if !isnothing(title)
        # Single-column layout: span the one column the axes actually occupy.
        Label(fig[0, 1], title; fontsize = 18, font = :bold, tellwidth = false)
    end

    # Drive IQ
    ax1 = Axis(
        fig[1, 1];
        xlabel = "Time (ns)",
        xlabelsize = 16,
        ylabel = "Amplitude (rad⋅GHz)",
        ylabelsize = 16,
        title = L"\mathrm{Drive}\;(\Omega)",
        titlesize = 18,
        titlealign = :left,
    )
    lines!(ax1, times, Ω_I; label = L"\Omega_I", color = c_I, linewidth = 2)
    lines!(ax1, times, Ω_Q; label = L"\Omega_Q", color = c_Q, linewidth = 2)
    lines!(
        ax1,
        times,
        Ω_mag;
        label = L"|\Omega|",
        color = c_mag,
        linestyle = :dash,
        linewidth = 2,
    )
    axislegend(ax1; position = :rt)

    # Displacement IQ
    ax2 = Axis(
        fig[2, 1];
        xlabel = "Time (ns)",
        xlabelsize = 16,
        ylabel = "Amplitude",
        ylabelsize = 16,
        title = L"\mathrm{Displacement}\;(\alpha)",
        titlesize = 18,
        titlealign = :left,
    )
    lines!(ax2, times, α_I; label = L"\alpha_I", color = c_I, linewidth = 2)
    lines!(ax2, times, α_Q; label = L"\alpha_Q", color = c_Q, linewidth = 2)
    lines!(
        ax2,
        times,
        α_mag;
        label = L"|\alpha|",
        color = c_mag,
        linestyle = :dash,
        linewidth = 2,
    )
    axislegend(ax2; position = :rt)

    # Overlay knot points
    if show_knots && pulse isa AbstractSplinePulse
        knot_times = get_knot_times(pulse)
        knot_controls = sample(pulse, knot_times)
        for (ax_row, idxs) in [(ax1, 1:2), (ax2, 3:4)]
            for i in idxs
                scatter!(
                    ax_row,
                    knot_times,
                    knot_controls[i, :];
                    markersize = 5,
                    color = c_mag,
                )
            end
        end
    end

    return fig
end

"""
    plot_pulse_phases(pulse::AbstractPulse; n_samples=500, title=nothing, figsize=(900, 500), phase_threshold=0.01)

Plot a 4-drive pulse in polar form: magnitude and phase for drive and displacement.

Phase is masked to NaN wherever the corresponding amplitude is below
`phase_threshold * max(amplitude)`, since `atan(Q, I)` is numerically meaningless
in the near-zero region and otherwise produces high-frequency noise that swamps
the meaningful phase signal. Phase axes are pinned to `[-1, 1]` (units of π).

Returns a Makie `Figure` with 4 subplots: |Ω|, φ_Ω, |α|, φ_α.
"""
function plot_pulse_phases(
    pulse::AbstractPulse;
    n_samples::Int = 500,
    title::Union{Nothing,String} = nothing,
    figsize::Tuple = (900, 500),
    phase_threshold::Float64 = 0.01,
)
    @assert n_drives(pulse) == 4 "plot_pulse_phases requires exactly 4 drives"

    controls, times = sample(pulse, n_samples)

    Ω_I, Ω_Q = controls[1, :], controls[2, :]
    α_I, α_Q = controls[3, :], controls[4, :]
    Ω_mag = sqrt.(Ω_I .^ 2 .+ Ω_Q .^ 2)
    α_mag = sqrt.(α_I .^ 2 .+ α_Q .^ 2)

    # Phase is undefined where amplitude → 0; suppress noise where |·| is below
    # `phase_threshold` × max amplitude by NaN-masking. Makie skips NaN segments
    # rather than drawing wild oscillations that swamp the actual phase signal.
    Ω_peak = maximum(Ω_mag)
    α_peak = maximum(α_mag)
    Ω_cut = phase_threshold * (Ω_peak > 0 ? Ω_peak : 1.0)
    α_cut = phase_threshold * (α_peak > 0 ? α_peak : 1.0)
    Ω_phase = [m > Ω_cut ? atan(q, i) / π : NaN for (m, q, i) in zip(Ω_mag, Ω_Q, Ω_I)]
    α_phase = [m > α_cut ? atan(q, i) / π : NaN for (m, q, i) in zip(α_mag, α_Q, α_I)]

    c_drive = _drive_color(1)
    c_disp = _drive_color(2)

    fig = Figure(size = figsize)
    if !isnothing(title)
        Label(fig[0, 1:2], title; fontsize = 18, font = :bold, tellwidth = false)
    end

    ax1 = Axis(
        fig[1, 1];
        ylabel = L"|\Omega|\;(\mathrm{rad}\!\cdot\!\mathrm{GHz})",
        ylabelsize = 16,
        title = "Drive magnitude",
        titlesize = 18,
        titlealign = :left,
    )
    lines!(ax1, times, Ω_mag; color = c_drive, linewidth = 2)

    # Pin phase axes to [-1, 1] (units of π) so the visible range stays
    # consistent across runs and isn't blown out by a stray near-zero point.
    ax2 = Axis(
        fig[1, 2];
        ylabel = L"\varphi_\Omega / \pi",
        ylabelsize = 16,
        title = "Drive phase",
        titlesize = 18,
        titlealign = :left,
        yticks = -1.0:0.5:1.0,
    )
    ylims!(ax2, -1.05, 1.05)
    lines!(ax2, times, Ω_phase; color = c_drive, linewidth = 2)

    ax3 = Axis(
        fig[2, 1];
        xlabel = "Time (ns)",
        xlabelsize = 16,
        ylabel = L"|\alpha|",
        ylabelsize = 16,
        title = "Displacement magnitude",
        titlesize = 18,
        titlealign = :left,
    )
    lines!(ax3, times, α_mag; color = c_disp, linewidth = 2)

    ax4 = Axis(
        fig[2, 2];
        xlabel = "Time (ns)",
        xlabelsize = 16,
        ylabel = L"\varphi_\alpha / \pi",
        ylabelsize = 16,
        title = "Displacement phase",
        titlesize = 18,
        titlealign = :left,
        yticks = -1.0:0.5:1.0,
    )
    ylims!(ax4, -1.05, 1.05)
    lines!(ax4, times, α_phase; color = c_disp, linewidth = 2)

    return fig
end

# ============================================================================ #
# High-level overloads — plot from qtraj / qcp directly
# ============================================================================ #

# Default labels derived from the pulse's drive_name. `:u` → [L"u_{1}", L"u_{2}", …].
# Uses `latexstring` so labels render in math font with italic letters and
# proper subscripts — same convention as NamedTrajectories.jl's plot recipe.
function _default_drive_labels(pulse::AbstractPulse)
    nd = n_drives(pulse)
    base = hasfield(typeof(pulse), :drive_name) ? String(drive_name(pulse)) : "u"
    return [latexstring(base, "_{", i, "}") for i = 1:nd]
end

# Convert system.drive_bounds (`Vector{Tuple{Float64,Float64}}`) into the form
# `plot_pulse` expects, with a length sanity check. Returns `nothing` if the
# pulse drive count doesn't match (e.g. CompositePulse).
function _system_bounds(sys::AbstractQuantumSystem, nd::Int)
    if length(sys.drive_bounds) != nd
        @warn "System has $(length(sys.drive_bounds)) bounds but pulse has $nd drives — bounds skipped"
        return nothing
    end
    return [(Float64(b[1]), Float64(b[2])) for b in sys.drive_bounds]
end

# NamedTrajectory stores `bounds[name] = (lower::Vector, upper::Vector)`.
# Return a per-dimension list of (lo, hi) tuples, or `nothing` if not bounded.
function _component_bounds(traj::NamedTrajectory, name::Symbol)
    if !haskey(traj.bounds, name)
        return nothing
    end
    lo, hi = traj.bounds[name]
    return [(Float64(lo[i]), Float64(hi[i])) for i in eachindex(lo)]
end

"""
    plot_pulse(qtraj::AbstractQuantumTrajectory; bounds=false, labels=nothing, kwargs...)

Plot the pulse stored on `qtraj`. Drive labels default to
`["<drive_name>_1", "<drive_name>_2", …]`. Pass `bounds=true` to shade
hardware bounds derived from `get_system(qtraj).drive_bounds`.

All other `plot_pulse(::AbstractPulse)` keyword arguments are forwarded.
"""
function plot_pulse(
    qtraj::AbstractQuantumTrajectory;
    bounds::Bool = false,
    labels::Union{Nothing,AbstractVector{<:AbstractString}} = nothing,
    kwargs...,
)
    pulse = get_pulse(qtraj)
    actual_bounds = bounds ? _system_bounds(get_system(qtraj), n_drives(pulse)) : nothing
    return plot_pulse(pulse; bounds = actual_bounds, labels, kwargs...)
end

"""
    plot_pulse(qcp::QuantumControlProblem; bounds=false, components=Symbol[], component_bounds=false, labels=nothing, kwargs...)

Plot the (possibly optimized) pulse from a `QuantumControlProblem`.

# Keyword Arguments
- `bounds::Bool=false`: When `true`, derive per-drive bounds from
  `get_system(qcp).drive_bounds` and shade them on the pulse subplots.
- `components::Vector{Symbol}=Symbol[]`: Additional `NamedTrajectory` components
  to stack under the pulse, e.g. `[:du, :ddu]`. Each renders as one row showing
  all dimensions overlaid.
- `component_bounds::Bool=false`: When `true`, read per-component bounds from
  `get_trajectory(qcp).bounds` and shade them under each component panel. The
  first dimension's bounds are used as a representative band.
- `labels::Union{Nothing,Vector{String}}=nothing`: Pulse drive labels (default
  derived from `drive_name`).
- All other `plot_pulse(::AbstractPulse)` kwargs are forwarded.
"""
function plot_pulse(
    qcp::QuantumControlProblem;
    bounds::Bool = false,
    components::Vector{Symbol} = Symbol[],
    component_bounds::Bool = false,
    labels::Union{Nothing,AbstractVector{<:AbstractString}} = nothing,
    figsize::Union{Nothing,Tuple} = nothing,
    title::AbstractString = "",
    kwargs...,
)
    traj = get_trajectory(qcp)
    pulse = extract_pulse(qcp.qtraj, traj)
    nd = n_drives(pulse)

    if isnothing(labels)
        labels = _default_drive_labels(pulse)
    else
        labels = [_to_latex_label(s) for s in labels]
    end

    actual_bounds = bounds ? _system_bounds(get_system(qcp), nd) : nothing

    if isempty(components)
        return plot_pulse(pulse; bounds = actual_bounds, labels, title, figsize, kwargs...)
    end

    return _plot_pulse_qcp_with_components(
        pulse,
        traj,
        components;
        bounds = actual_bounds,
        component_bounds,
        labels,
        title,
        figsize,
        kwargs...,
    )
end

function _plot_pulse_qcp_with_components(
    pulse::AbstractPulse,
    traj::NamedTrajectory,
    components::Vector{Symbol};
    bounds::Union{Nothing,AbstractVector} = nothing,
    component_bounds::Bool = false,
    labels::AbstractVector{<:AbstractString},
    title::AbstractString = "",
    figsize::Union{Nothing,Tuple} = nothing,
    n_samples::Int = 500,
    show_knots::Bool = true,
    show_tangents::Bool = false,
    tangent_scale::Float64 = 0.1,
    kwargs...,
)
    nd = n_drives(pulse)
    n_total_rows = nd + length(components)
    sz = isnothing(figsize) ? (800, 100 + 140 * n_total_rows) : figsize
    fig = Figure(; size = sz, kwargs...)

    # Pulse rows: same per-drive rendering as plot_pulse(::AbstractPulse, layout=:stacked),
    # but with x-tick suppression (components rows below own the x-axis label).
    for i = 1:nd
        ax = Axis(
            fig[i, 1];
            titlealign = :left,
            titlesize = 18,
            titlefont = :bold,
            title = (i == 1 && !isempty(title)) ? title : "",
            xticklabelsvisible = false,
            xtickalign = 1,
            ylabel = labels[min(i, length(labels))],
            ylabelsize = 16,
        )
        bound_pair = _drive_bound_pair(bounds, i)
        _draw_pulse_subplot!(
            ax,
            pulse,
            i;
            bound_pair,
            n_samples,
            show_knots,
            show_tangents,
            tangent_scale,
        )
    end

    # Component rows: one row per component, all dimensions overlaid.
    # Component ylabels are LaTeX-rendered: `:du` → \dot{u}, `:ddu` → \ddot{u}.
    for (k, comp) in enumerate(components)
        row = nd + k
        is_last = row == n_total_rows
        ax = Axis(
            fig[row, 1];
            xticklabelsvisible = is_last,
            xtickalign = 1,
            xlabel = is_last ? "Time" : "",
            xlabelsize = 16,
            ylabel = _component_latex(comp),
            ylabelsize = 16,
        )

        # Use dim-1 bounds as a representative band — components are typically
        # bounded uniformly across drives in our problem templates.
        comp_bound_pair = nothing
        if component_bounds
            cb = _component_bounds(traj, comp)
            if !isnothing(cb) && !isempty(cb)
                comp_bound_pair = cb[1]
            end
        end

        _draw_component_subplot!(ax, traj, comp; bound_pair = comp_bound_pair)
    end

    # Tight rowgap to match plot_pulse(::AbstractPulse) stacked layout.
    for i = 1:(n_total_rows-1)
        rowgap!(fig.layout, i, Relative(0.015))
    end

    return fig
end

# ============================================================================ #
# Tests
# ============================================================================ #

@testitem "plot_pulse ZeroOrderPulse" begin
    using CairoMakie
    using Piccolo

    controls = [0.0 1.0 0.5 0.0; 0.0 -1.0 -0.5 0.0]
    times = [0.0, 0.25, 0.5, 1.0]
    pulse = ZeroOrderPulse(controls, times)

    # Stacked layout
    fig = plot_pulse(pulse)
    @test fig isa Figure

    # Overlay layout
    fig2 = plot_pulse(pulse; layout = :overlay, labels = ["Ωx", "Ωy"])
    @test fig2 isa Figure

    # With bounds
    fig3 = plot_pulse(pulse; bounds = [(-1.5, 1.5), (-1.5, 1.5)])
    @test fig3 isa Figure
end

@testitem "plot_pulse LinearSplinePulse" begin
    using CairoMakie
    using Piccolo

    controls = [0.0 1.0 0.0; 0.0 -1.0 0.0]
    times = [0.0, 0.5, 1.0]
    pulse = LinearSplinePulse(controls, times)

    fig = plot_pulse(pulse)
    @test fig isa Figure

    fig2 = plot_pulse(pulse; layout = :overlay, labels = ["I", "Q"])
    @test fig2 isa Figure
end

@testitem "plot_pulse CubicSplinePulse" begin
    using CairoMakie
    using Piccolo

    controls = [0.0 0.5 1.0 0.5 0.0; 0.0 -0.5 -1.0 -0.5 0.0]
    derivatives = [2.0 2.0 0.0 -2.0 -2.0; -2.0 -2.0 0.0 2.0 2.0]
    times = [0.0, 0.25, 0.5, 0.75, 1.0]
    pulse = CubicSplinePulse(controls, derivatives, times)

    # Without tangents
    fig = plot_pulse(pulse)
    @test fig isa Figure

    # With tangents
    fig2 = plot_pulse(pulse; show_tangents = true, tangent_scale = 0.05)
    @test fig2 isa Figure

    # Overlay
    fig3 = plot_pulse(pulse; layout = :overlay)
    @test fig3 isa Figure
end

@testitem "plot_pulse GaussianPulse" begin
    using CairoMakie
    using Piccolo

    pulse = GaussianPulse([1.0, 2.0], 0.1, 1.0)

    fig = plot_pulse(pulse)
    @test fig isa Figure

    fig2 = plot_pulse(pulse; layout = :overlay, title = "Gaussian")
    @test fig2 isa Figure
end

@testitem "plot_pulse ErfPulse" begin
    using CairoMakie
    using Piccolo

    pulse = ErfPulse([1.0], 0.2, 1.0)

    fig = plot_pulse(pulse)
    @test fig isa Figure
end

@testitem "plot_pulse CompositePulse" begin
    using CairoMakie
    using Piccolo

    amp = GaussianPulse([1.0, 2.0], 1.0, 10.0)
    phase = ErfPulse([0.5, 0.8], 1.0, 10.0)
    pulse = CompositePulse([amp, phase], :interleave)

    fig = plot_pulse(pulse; labels = ["Ω₁", "φ₁", "Ω₂", "φ₂"])
    @test fig isa Figure
end

@testitem "plot_pulse FunctionPulse" begin
    using CairoMakie
    using Piccolo

    T = 1.0
    pulse = FunctionPulse(t -> [sin(π * t / T)^2, cos(π * t / T)^2], T, 2)

    fig = plot_pulse(pulse; title = "FunctionPulse")
    @test fig isa Figure

    fig2 = plot_pulse(pulse; layout = :overlay, labels = ["sin²", "cos²"])
    @test fig2 isa Figure
end

@testitem "plot_pulse! on existing axis" begin
    using CairoMakie
    using Piccolo

    pulse = GaussianPulse([1.0], 0.1, 1.0)

    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_pulse!(ax, pulse)
    @test fig isa Figure
end

@testitem "plot_pulse honors theme_dark" begin
    using CairoMakie
    using Piccolo

    # Render the same pulse under both themes; both should produce a Figure
    # without throwing on missing palette entries / hardcoded colors.
    pulse = GaussianPulse([1.0, 0.5], 0.2, 1.0)

    set_theme!(theme_light())
    fig_light = plot_pulse(pulse; show_knots = true)
    @test fig_light isa Figure

    set_theme!(theme_dark())
    fig_dark = plot_pulse(pulse; bounds = [(-1.5, 1.5), (-1.5, 1.5)])
    @test fig_dark isa Figure

    # Restore default
    set_theme!()
end

@testitem "plot_pulse_IQ 4-drive spline" begin
    using CairoMakie
    using Piccolo
    using Random

    Random.seed!(0)
    T = 50.0
    n_knots = 16
    times = collect(range(0, T, length = n_knots))
    controls = 0.5 * randn(4, n_knots)
    pulse = LinearSplinePulse(controls, times)

    # Default: knot overlay enabled (spline branch)
    fig = plot_pulse_IQ(pulse)
    @test fig isa Figure

    # With title + custom figsize, knots disabled
    fig2 = plot_pulse_IQ(pulse; title = "MS gate", figsize = (700, 500), show_knots = false)
    @test fig2 isa Figure
end

@testitem "plot_pulse_IQ 4-drive analytic" begin
    using CairoMakie
    using Piccolo

    # Analytic 4-drive pulse: spline-knot overlay branch is skipped.
    pulse = GaussianPulse([1.0, 0.7, 0.5, 0.3], 0.2, 1.0)
    fig = plot_pulse_IQ(pulse)
    @test fig isa Figure
end

@testitem "plot_pulse_IQ rejects non-4-drive" begin
    using CairoMakie
    using Piccolo

    pulse = GaussianPulse([1.0, 2.0], 0.1, 1.0)  # 2 drives
    @test_throws AssertionError plot_pulse_IQ(pulse)
end

@testitem "plot_pulse_phases 4-drive" begin
    using CairoMakie
    using Piccolo

    pulse = GaussianPulse([1.0, 0.7, 0.5, 0.3], 0.2, 1.0)

    fig = plot_pulse_phases(pulse)
    @test fig isa Figure

    fig2 = plot_pulse_phases(pulse; title = "Polar view", figsize = (800, 500))
    @test fig2 isa Figure
end

@testitem "plot_pulse_phases rejects non-4-drive" begin
    using CairoMakie
    using Piccolo

    pulse = GaussianPulse([1.0], 0.1, 1.0)  # 1 drive
    @test_throws AssertionError plot_pulse_phases(pulse)
end

@testitem "plot_pulse(qtraj::UnitaryTrajectory) basic + bounds" begin
    using CairoMakie
    using Piccolo
    using Random

    Random.seed!(0)
    sys = QuantumSystem(0.5 * PAULIS[:Z], [PAULIS[:X], PAULIS[:Y]], [1.0, 0.8])
    times = collect(range(0, 1.0, length = 20))
    pulse = ZeroOrderPulse(0.1 * randn(2, 20), times)
    qtraj = UnitaryTrajectory(sys, pulse, GATES[:X])

    fig1 = plot_pulse(qtraj)
    @test fig1 isa Figure

    fig2 = plot_pulse(qtraj; bounds = true)
    @test fig2 isa Figure
end

@testitem "plot_pulse(qtraj::KetTrajectory) smoke" begin
    using CairoMakie
    using Piccolo
    using Random

    Random.seed!(1)
    sys = QuantumSystem(0.5 * PAULIS[:Z], [PAULIS[:X]], [1.5])
    times = collect(range(0, 1.0, length = 20))
    pulse = LinearSplinePulse(0.1 * randn(1, 20), times)
    qtraj = KetTrajectory(sys, pulse, ComplexF64[1, 0], ComplexF64[0, 1])

    fig = plot_pulse(qtraj; bounds = true, labels = ["Ωx"])
    @test fig isa Figure
end

@testitem "plot_pulse(qcp) basic, bounds, components" begin
    using CairoMakie
    using Piccolo
    using Random

    Random.seed!(2)
    N = 30
    T = 1.0
    sys = QuantumSystem(0.5 * PAULIS[:Z], [PAULIS[:X], PAULIS[:Y]], [1.0, 1.0])
    times = collect(range(0, T, length = N))
    pulse = ZeroOrderPulse(0.1 * randn(2, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, GATES[:X])
    qcp = SmoothPulseProblem(qtraj, N; Q = 100.0, R = 1e-2, ddu_bound = 1.0)

    # Don't actually solve — just verify the trajectory plumbing works at the
    # initial-guess state.

    fig1 = plot_pulse(qcp)
    @test fig1 isa Figure

    fig2 = plot_pulse(qcp; bounds = true)
    @test fig2 isa Figure

    fig3 = plot_pulse(qcp; components = [:du, :ddu])
    @test fig3 isa Figure

    fig4 = plot_pulse(qcp; bounds = true, components = [:du, :ddu], component_bounds = true)
    @test fig4 isa Figure
end
