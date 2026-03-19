export plot_pulse, plot_pulse!

using Makie
using Piccolo: AbstractPulse, AbstractSplinePulse, CubicSplinePulse,
    duration, n_drives, sample, get_knot_times, get_knot_values, drive_name

"""
    plot_pulse(pulse::AbstractPulse; n_samples=500, labels=nothing, title=nothing, kwargs...)

Plot a pulse's control channels over time.

Returns a Makie `Figure` with one axis showing all drive channels as lines.
For spline pulses, knot points are overlaid as scatter markers.

# Arguments
- `pulse::AbstractPulse`: The pulse to plot

# Keyword Arguments
- `n_samples::Int=500`: Number of time samples for smooth curves
- `labels::Union{Nothing, Vector{String}}=nothing`: Labels for each drive channel
- `title::Union{Nothing, String}=nothing`: Plot title
- `figsize::Tuple=(800, 400)`: Figure size
- `show_knots::Bool=true`: Show knot points for spline pulses
- `kwargs...`: Additional keyword arguments passed to `lines!`
"""
function plot_pulse(
    pulse::AbstractPulse;
    n_samples::Int = 500,
    labels::Union{Nothing, Vector{String}} = nothing,
    title::Union{Nothing, String} = nothing,
    figsize::Tuple = (800, 400),
    show_knots::Bool = true,
    kwargs...,
)
    fig = Figure(size = figsize)
    ax = Axis(fig[1, 1];
        xlabel = "Time (ns)",
        ylabel = "Amplitude",
        title = isnothing(title) ? "Pulse Controls" : title,
    )
    plot_pulse!(ax, pulse; n_samples, labels, show_knots, kwargs...)
    if !isnothing(labels)
        axislegend(ax; position = :rt)
    end
    return fig
end

"""
    plot_pulse!(ax, pulse::AbstractPulse; n_samples=500, labels=nothing, show_knots=true, kwargs...)

Plot a pulse's control channels onto an existing Makie axis.
"""
function plot_pulse!(
    ax,
    pulse::AbstractPulse;
    n_samples::Int = 500,
    labels::Union{Nothing, Vector{String}} = nothing,
    show_knots::Bool = true,
    kwargs...,
)
    controls, times = sample(pulse; n_samples)
    nd = n_drives(pulse)

    for i in 1:nd
        label = isnothing(labels) ? nothing : labels[i]
        lines!(ax, times, controls[i, :]; label, kwargs...)
    end

    # Overlay knot points for spline pulses
    if show_knots && pulse isa AbstractSplinePulse
        knot_times = get_knot_times(pulse)
        knot_controls = sample(pulse, knot_times)
        for i in 1:nd
            scatter!(ax, knot_times, knot_controls[i, :];
                markersize = 6, color = :black)
        end
    end
end

"""
    plot_pulse_IQ(pulse::AbstractPulse; n_samples=500, drive_labels=nothing, title=nothing, kwargs...)

Plot a 4-drive pulse as two IQ pairs: drives (Ω_I, Ω_Q) and displacement (α_I, α_Q),
with magnitude envelopes.

Returns a Makie `Figure` with 2 rows: drive IQ and displacement IQ.
"""
function plot_pulse_IQ(
    pulse::AbstractPulse;
    n_samples::Int = 500,
    title::Union{Nothing, String} = nothing,
    figsize::Tuple = (900, 600),
    show_knots::Bool = true,
)
    @assert n_drives(pulse) == 4 "plot_pulse_IQ requires exactly 4 drives (Ω_I, Ω_Q, α_I, α_Q)"

    controls, times = sample(pulse; n_samples)

    Ω_I, Ω_Q = controls[1, :], controls[2, :]
    α_I, α_Q = controls[3, :], controls[4, :]
    Ω_mag = sqrt.(Ω_I.^2 .+ Ω_Q.^2)
    α_mag = sqrt.(α_I.^2 .+ α_Q.^2)

    fig = Figure(size = figsize)
    if !isnothing(title)
        Label(fig[0, 1:2], title; fontsize = 16, tellwidth = false)
    end

    # Drive IQ
    ax1 = Axis(fig[1, 1]; xlabel = "Time (ns)", ylabel = "Amplitude (rad⋅GHz)",
        title = "Drive (Ω)")
    lines!(ax1, times, Ω_I; label = "Ω_I", color = :blue)
    lines!(ax1, times, Ω_Q; label = "Ω_Q", color = :red)
    lines!(ax1, times, Ω_mag; label = "|Ω|", color = :black, linestyle = :dash)
    axislegend(ax1; position = :rt)

    # Displacement IQ
    ax2 = Axis(fig[2, 1]; xlabel = "Time (ns)", ylabel = "Amplitude",
        title = "Displacement (α)")
    lines!(ax2, times, α_I; label = "α_I", color = :blue)
    lines!(ax2, times, α_Q; label = "α_Q", color = :red)
    lines!(ax2, times, α_mag; label = "|α|", color = :black, linestyle = :dash)
    axislegend(ax2; position = :rt)

    # Overlay knot points
    if show_knots && pulse isa AbstractSplinePulse
        knot_times = get_knot_times(pulse)
        knot_controls = sample(pulse, knot_times)
        for (ax_row, idxs) in [(ax1, 1:2), (ax2, 3:4)]
            for i in idxs
                scatter!(ax_row, knot_times, knot_controls[i, :];
                    markersize = 5, color = :black)
            end
        end
    end

    return fig
end

"""
    plot_pulse_phases(pulse::AbstractPulse; n_samples=500, title=nothing, kwargs...)

Plot a 4-drive pulse in polar form: magnitude and phase for drive and displacement.

Returns a Makie `Figure` with 4 subplots: |Ω|, φ_Ω, |α|, φ_α.
"""
function plot_pulse_phases(
    pulse::AbstractPulse;
    n_samples::Int = 500,
    title::Union{Nothing, String} = nothing,
    figsize::Tuple = (900, 600),
)
    @assert n_drives(pulse) == 4 "plot_pulse_phases requires exactly 4 drives"

    controls, times = sample(pulse; n_samples)

    Ω_I, Ω_Q = controls[1, :], controls[2, :]
    α_I, α_Q = controls[3, :], controls[4, :]
    Ω_mag = sqrt.(Ω_I.^2 .+ Ω_Q.^2)
    α_mag = sqrt.(α_I.^2 .+ α_Q.^2)
    Ω_phase = atan.(Ω_Q, Ω_I) ./ π
    α_phase = atan.(α_Q, α_I) ./ π

    fig = Figure(size = figsize)
    if !isnothing(title)
        Label(fig[0, 1:2], title; fontsize = 16, tellwidth = false)
    end

    ax1 = Axis(fig[1, 1]; ylabel = "|Ω| (rad⋅GHz)", title = "Drive magnitude")
    lines!(ax1, times, Ω_mag; color = :blue)

    ax2 = Axis(fig[1, 2]; ylabel = "φ_Ω / π", title = "Drive phase")
    lines!(ax2, times, Ω_phase; color = :blue)

    ax3 = Axis(fig[2, 1]; xlabel = "Time (ns)", ylabel = "|α|", title = "Displacement magnitude")
    lines!(ax3, times, α_mag; color = :red)

    ax4 = Axis(fig[2, 2]; xlabel = "Time (ns)", ylabel = "φ_α / π", title = "Displacement phase")
    lines!(ax4, times, α_phase; color = :red)

    return fig
end
