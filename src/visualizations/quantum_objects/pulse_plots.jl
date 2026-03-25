export plot_pulse, plot_pulse!

using Makie
using Piccolo:
    AbstractPulse,
    AbstractSplinePulse,
    ZeroOrderPulse,
    LinearSplinePulse,
    CubicSplinePulse,
    GaussianPulse,
    ErfPulse,
    CompositePulse,
    duration,
    n_drives,
    sample,
    get_knot_times,
    get_knot_values,
    get_knot_derivatives,
    drive_name

using TestItems

# ============================================================================ #
# Color palette
# ============================================================================ #

# Wong colors for accessibility (same palette used in gate_populations, etc.)
const PULSE_COLORS = [
    Makie.wong_colors()[1],  # blue
    Makie.wong_colors()[2],  # orange
    Makie.wong_colors()[3],  # green
    Makie.wong_colors()[4],  # red
    Makie.wong_colors()[5],  # purple
    Makie.wong_colors()[6],  # brown
    Makie.wong_colors()[7],  # pink
]

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
- `GaussianPulse`, `ErfPulse`, `CompositePulse`: Dense-sampled smooth curve

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
    labels::Union{Nothing,Vector{String}} = nothing,
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

    if isnothing(labels)
        labels = ["Drive $i" for i = 1:nd]
    end

    if layout == :stacked
        sz = isnothing(figsize) ? (800, 100 + 140 * nd) : figsize
        fig = Figure(; size = sz, kwargs...)

        for i = 1:nd
            ax = Axis(
                fig[i, 1];
                titlealign = :left,
                titlesize = 16,
                titlefont = :bold,
                title = (i == 1 && !isempty(title)) ? title : "",
                xticklabelsvisible = i == nd,
                xtickalign = 1,
                xlabel = i == nd ? "Time" : "",
                ylabel = labels[min(i, length(labels))],
            )

            # Hardware bounds as shaded band
            if !isnothing(bounds) && i <= length(bounds) && !isnothing(bounds[i])
                _draw_bounds!(ax, pulse, bounds[i])
            end

            plot_pulse!(
                ax,
                pulse;
                drive_indices = [i],
                n_samples,
                show_knots,
                show_tangents,
                tangent_scale,
                colors = [PULSE_COLORS[mod1(i, length(PULSE_COLORS))]],
            )

            hlines!(ax, [0.0]; color = :gray, linestyle = :dash, linewidth = 0.5)
        end

        return fig

    else  # :overlay
        sz = isnothing(figsize) ? (800, 400) : figsize
        fig = Figure(; size = sz, kwargs...)
        ax = Axis(
            fig[1, 1];
            xlabel = "Time",
            ylabel = "Amplitude",
            title = isempty(title) ? "Pulse Controls" : title,
        )

        colors = [PULSE_COLORS[mod1(i, length(PULSE_COLORS))] for i = 1:nd]
        plot_pulse!(ax, pulse; n_samples, show_knots, show_tangents, tangent_scale, colors)

        # Legend
        entries = [LineElement(; color = colors[i], linewidth = 2) for i = 1:nd]
        Legend(fig[1, 2], entries, labels)

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
- `colors`: Vector of colors, one per drive index.
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
    kwargs...,
)
    if isnothing(colors)
        colors = [PULSE_COLORS[mod1(i, length(PULSE_COLORS))] for i in drive_indices]
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
        kwargs...,
    )
end

# ============================================================================ #
# Type-specific rendering dispatches
# ============================================================================ #

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
    kwargs...,
)
    knot_times = collect(get_knot_times(pulse))
    knot_vals = sample(pulse, knot_times)

    for (ci, i) in enumerate(drive_indices)
        stairs!(
            ax,
            knot_times,
            knot_vals[i, :];
            color = colors[ci],
            linewidth = 2,
            step = :pre,
            kwargs...,
        )

        if show_knots
            scatter!(
                ax,
                knot_times,
                knot_vals[i, :];
                color = colors[ci],
                markersize = 6,
                strokewidth = 1,
                strokecolor = :black,
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
    kwargs...,
)
    knot_times = collect(get_knot_times(pulse))
    knot_vals = get_knot_values(pulse)

    for (ci, i) in enumerate(drive_indices)
        lines!(
            ax,
            knot_times,
            collect(knot_vals[i, :]);
            color = colors[ci],
            linewidth = 2,
            kwargs...,
        )

        if show_knots
            scatter!(
                ax,
                knot_times,
                collect(knot_vals[i, :]);
                color = colors[ci],
                markersize = 6,
                strokewidth = 1,
                strokecolor = :black,
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
    kwargs...,
)
    # Dense sampling for smooth curve
    controls, times = sample(pulse; n_samples)

    for (ci, i) in enumerate(drive_indices)
        lines!(ax, times, controls[i, :]; color = colors[ci], linewidth = 2, kwargs...)
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
                strokecolor = :black,
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
# Fallback for analytic pulses (Gaussian, Erf, Composite, any future type)
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
    kwargs...,
)
    controls, times = sample(pulse; n_samples)

    for (ci, i) in enumerate(drive_indices)
        lines!(ax, times, controls[i, :]; color = colors[ci], linewidth = 2, kwargs...)
    end
end

# ============================================================================ #
# Helper: draw hardware bounds
# ============================================================================ #

function _draw_bounds!(ax, pulse::AbstractPulse, bound_pair)
    lo, hi = bound_pair
    t_start = 0.0
    t_end = duration(pulse)
    band!(
        ax,
        [t_start, t_end],
        [Float64(lo), Float64(lo)],
        [Float64(hi), Float64(hi)];
        color = (:steelblue, 0.1),
    )
    hlines!(
        ax,
        [lo, hi];
        color = :steelblue,
        linestyle = :dash,
        linewidth = 0.8,
        alpha = 0.5,
    )
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

@testitem "plot_pulse! on existing axis" begin
    using CairoMakie
    using Piccolo

    pulse = GaussianPulse([1.0], 0.1, 1.0)

    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_pulse!(ax, pulse)
    @test fig isa Figure
end
