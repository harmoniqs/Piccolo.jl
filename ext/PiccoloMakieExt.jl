module PiccoloMakieExt

using Piccolo
using Makie
using NamedTrajectories
using TestItems

# Animation implementations - extend Piccolo stubs
"""
    animate_figure(
        fig::Figure,
        frames::AbstractVector{Int},
        update_frame!::Function;
        mode::Symbol = :inline,
        fps::Int = 30,
        filename::String = "animation.mp4"
    )

Animate a Makie figure by updating frames. Works best with `GLMakie`.
"""
function Piccolo.animate_figure(
    fig::Figure,
    frames::AbstractVector{Int},
    update_frame!::Function;
    mode::Symbol = :inline,
    fps::Int = 24,
    filename::String = "animation.mp4",
)
    if mode == :inline
        display(fig) # open the scene
        if !isopen(fig.scene)
            @warn "Unable to open :inline animation for the current backend. This is a known limitation of CairoMakie. Consider setting mode = :record."
        end
        @async begin # don't block
            while isopen(fig.scene)
                for i in frames
                    update_frame!(i)
                    sleep(1 / fps)
                end
                if !isopen(fig.scene)
                    break # exit after close
                else
                    sleep(10 / fps)
                end
            end
        end
    elseif mode == :record
        record(fig, filename, frames; framerate = fps) do i
            update_frame!(i)
        end
    else
        throw(ArgumentError("Unsupported mode: $mode. Use :inline or :record."))
    end

    return fig
end


"""
    animate_name(
        traj::NamedTrajectory,
        name::Symbol;
        fps::Int = 24,
        mode::Symbol = :inline,
        fname::String = "name_animation.mp4",
        kwargs...
    )

Animate the evolution of a variable in a `NamedTrajectory`.
"""
function Piccolo.animate_name(
    traj::NamedTrajectory,
    name::Symbol;
    fps::Int = 24,
    mode::Symbol = :inline,
    filename = "name_animation.mp4",
    kwargs...,
)
    fig = Figure()
    ax = Axis(fig[1, 1])

    # Set limits
    times = get_times(traj)
    xlims = [times[1], times[end]]
    ylims = collect(extrema(traj[name]))

    scatter!(ax, xlims, ylims, markersize = 0.0)
    QuantumObjectPlots.plot_name!(ax, traj, name, indices = 1:1, kwargs...)

    # TODO: Unclear how to set observables via indices, so redraw
    # WARNING: Known bug with CairoMakie, which cannot re-render
    function update_frame!(i)
        empty!(ax.scene)
        scatter!(ax, xlims, ylims, markersize = 0.0)
        QuantumObjectPlots.plot_name!(ax, traj, name, indices = 1:i)
    end

    return Piccolo.animate_figure(
        fig,
        1:traj.N,
        update_frame!,
        mode = mode,
        fps = fps,
        filename = filename,
    )
end


"""
    animate_pulse(
        pulse::AbstractPulse;
        fps::Int = 24,
        mode::Symbol = :inline,
        filename::String = "pulse_animation.mp4",
        n_samples::Int = 240,
        labels = nothing,
        populations = nothing,
        population_times = nothing,
        population_labels = nothing,
        title::AbstractString = "Pulse evolution"
    )

Animate a pulse being drawn over time, with optional population/state traces in
a second panel. The pulse is sampled uniformly so this works with every
`AbstractPulse` implementation.
"""
function Piccolo.animate_pulse(
    pulse::Piccolo.AbstractPulse;
    fps::Int = 24,
    mode::Symbol = :inline,
    filename::String = "pulse_animation.mp4",
    n_samples::Int = 240,
    labels = nothing,
    populations = nothing,
    population_times = nothing,
    population_labels = nothing,
    title::AbstractString = "Pulse evolution",
)
    n_samples >= 2 || throw(ArgumentError("n_samples must be at least 2."))

    pulse_times = collect(range(0.0, Piccolo.duration(pulse), length = n_samples))
    pulse_values = Piccolo.sample(pulse, pulse_times)
    n_drive = Piccolo.n_drives(pulse)

    size(pulse_values, 1) == n_drive ||
        throw(ArgumentError("sample(pulse, times) must return one row per drive."))
    size(pulse_values, 2) == n_samples ||
        throw(ArgumentError("sample(pulse, times) must return one column per sample."))

    drive_labels = if isnothing(labels)
        _pulse_drive_labels(pulse, n_drive)
    else
        length(labels) == n_drive ||
            throw(ArgumentError("labels must contain one entry per pulse drive."))
        collect(labels)
    end

    has_populations = !isnothing(populations)
    population_values = nothing
    pop_times = nothing
    pop_labels = nothing
    n_pop = 0

    if has_populations
        population_values = Matrix(populations)
        size(population_values, 2) >= 2 ||
            throw(ArgumentError("populations must have at least two time samples."))

        n_pop = size(population_values, 1)
        pop_times =
            isnothing(population_times) ?
            collect(
                range(
                    pulse_times[1],
                    pulse_times[end],
                    length = size(population_values, 2),
                ),
            ) : collect(population_times)

        length(pop_times) == size(population_values, 2) || throw(
            ArgumentError("population_times must match the number of population samples."),
        )

        pop_labels = if isnothing(population_labels)
            ["Population $i" for i = 1:n_pop]
        else
            length(population_labels) == n_pop || throw(
                ArgumentError(
                    "population_labels must contain one entry per population row.",
                ),
            )
            collect(population_labels)
        end
    end

    fig = Figure(size = has_populations ? (900, 620) : (900, 380))
    pulse_axis = Axis(
        fig[1, 1],
        xlabel = has_populations ? "" : "Time",
        xticklabelsvisible = !has_populations,
        ylabel = "Amplitude",
        title = title,
        titlealign = :left,
    )

    pulse_y_min, pulse_y_max = _padded_extrema(pulse_values)
    xlims!(pulse_axis, pulse_times[1], pulse_times[end])
    ylims!(pulse_axis, pulse_y_min, pulse_y_max)
    hlines!(pulse_axis, [0.0], color = (:gray, 0.45), linestyle = :dash, linewidth = 1)

    colors = Makie.wong_colors()
    pulse_xs = [Observable(pulse_times[1:1]) for _ = 1:n_drive]
    pulse_ys = [Observable(vec(pulse_values[i, 1:1])) for i = 1:n_drive]
    marker_xs = [Observable([pulse_times[1]]) for _ = 1:n_drive]
    marker_ys = [Observable([pulse_values[i, 1]]) for i = 1:n_drive]

    for i = 1:n_drive
        color = colors[mod1(i, length(colors))]
        lines!(
            pulse_axis,
            pulse_times,
            vec(pulse_values[i, :]);
            color = (color, 0.18),
            linewidth = 1,
        )
        lines!(
            pulse_axis,
            pulse_xs[i],
            pulse_ys[i];
            color,
            linewidth = 3,
            label = string(drive_labels[i]),
        )
        scatter!(
            pulse_axis,
            marker_xs[i],
            marker_ys[i];
            color,
            markersize = 10,
            strokewidth = 1,
            strokecolor = :white,
        )
    end
    axislegend(pulse_axis, position = :rt)

    population_xs = Observable[]
    population_ys = Observable[]
    pop_marker_xs = Observable[]
    pop_marker_ys = Observable[]
    population_axis = nothing

    if has_populations
        population_axis = Axis(
            fig[2, 1],
            xlabel = "Time",
            ylabel = "Population",
            title = "State populations",
            titlealign = :left,
        )
        pop_y_min, pop_y_max = _padded_extrema(population_values)
        xlims!(population_axis, pop_times[1], pop_times[end])
        ylims!(population_axis, pop_y_min, pop_y_max)

        population_xs = [Observable(pop_times[1:1]) for _ = 1:n_pop]
        population_ys = [Observable(vec(population_values[i, 1:1])) for i = 1:n_pop]
        pop_marker_xs = [Observable([pop_times[1]]) for _ = 1:n_pop]
        pop_marker_ys = [Observable([population_values[i, 1]]) for i = 1:n_pop]

        for i = 1:n_pop
            color = colors[mod1(i + n_drive, length(colors))]
            lines!(
                population_axis,
                pop_times,
                vec(population_values[i, :]);
                color = (color, 0.18),
                linewidth = 1,
            )
            lines!(
                population_axis,
                population_xs[i],
                population_ys[i];
                color,
                linewidth = 3,
                label = string(pop_labels[i]),
            )
            scatter!(
                population_axis,
                pop_marker_xs[i],
                pop_marker_ys[i];
                color,
                markersize = 9,
                strokewidth = 1,
                strokecolor = :white,
            )
        end
        axislegend(population_axis, position = :rt)
        rowgap!(fig.layout, 1, Relative(0.06))
    end

    function update_frame!(i)
        for drive_i = 1:n_drive
            pulse_xs[drive_i][] = pulse_times[1:i]
            pulse_ys[drive_i][] = vec(pulse_values[drive_i, 1:i])
            marker_xs[drive_i][] = [pulse_times[i]]
            marker_ys[drive_i][] = [pulse_values[drive_i, i]]
        end

        if has_populations
            pop_i = clamp(
                round(Int, 1 + (i - 1) * (length(pop_times) - 1) / (n_samples - 1)),
                1,
                length(pop_times),
            )
            for state_i = 1:n_pop
                population_xs[state_i][] = pop_times[1:pop_i]
                population_ys[state_i][] = vec(population_values[state_i, 1:pop_i])
                pop_marker_xs[state_i][] = [pop_times[pop_i]]
                pop_marker_ys[state_i][] = [population_values[state_i, pop_i]]
            end
        end
    end

    return Piccolo.animate_figure(
        fig,
        1:n_samples,
        update_frame!,
        mode = mode,
        fps = fps,
        filename = filename,
    )
end


function _padded_extrema(values)
    lo, hi = extrema(values)
    if lo == hi
        pad = max(abs(lo), 1.0) * 0.1
    else
        pad = 0.08 * (hi - lo)
    end
    return lo - pad, hi + pad
end


function _pulse_drive_labels(pulse, n_drive::Int)
    base = try
        string(Piccolo.drive_name(pulse))
    catch
        "Drive"
    end
    return [base == "Drive" ? "Drive $i" : "$(base)_$i" for i = 1:n_drive]
end


@testitem "Test animate_name for NamedTrajectory" begin
    using QuantumToolbox
    using NamedTrajectories
    using CairoMakie

    x = ComplexF64[1.0; 0.0]
    y = ComplexF64[0.0, 1.0]

    comps = (ψ̃ = hcat(ket_to_iso(x), ket_to_iso(y)), Δt = [1.0; 1.0])
    traj = NamedTrajectory(comps)

    fig = animate_name(traj, :ψ̃, mode = :inline, fps = 10)
    @test fig isa Figure
end


@testitem "Test animate_pulse for AbstractPulse" begin
    using CairoMakie

    pulse = GaussianPulse([1.0, 0.5], 0.2, 1.0)
    populations = [
        range(1.0, 0.0, length = 12)';
        range(0.0, 1.0, length = 12)';
    ]

    fig = animate_pulse(
        pulse;
        mode = :inline,
        fps = 10,
        n_samples = 12,
        populations,
        population_labels = ["|0>", "|1>"],
    )
    @test fig isa Figure
end


end
