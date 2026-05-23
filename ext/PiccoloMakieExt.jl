module PiccoloMakieExt

using Piccolo
using Makie
using NamedTrajectories
using TestItems

using DirectTrajOpt: AbstractIntermediateCallback
using NamedTrajectories: update!
using Piccolo: AbstractQuantumTrajectory, extract_pulse, plot_pulse

# Animation implementations - extend Piccolo stubs.
# Docstrings live on the stubs in `src/visualizations/animations.jl`.

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


function _pulse_sample_data(pulse::Piccolo.AbstractPulse, n_samples::Int)
    controls, times = Piccolo.sample(pulse, n_samples)
    return Matrix(controls), collect(times)
end


function _pulse_y_limits(values)
    finite_values = Float64[x for x in vec(values) if isfinite(Float64(x))]
    if isempty(finite_values)
        return (-1.0, 1.0)
    end

    lo, hi = extrema(finite_values)
    pad = 0.08 * (hi - lo)
    if pad == 0
        pad = 0.1 * max(abs(lo), 1.0)
    end
    return (lo - pad, hi + pad)
end


function _pulse_drive_labels(pulse::Piccolo.AbstractPulse, labels)
    nd = Piccolo.n_drives(pulse)
    if isnothing(labels)
        return ["Drive $i" for i = 1:nd]
    end

    @assert length(labels) >= nd "labels must include one label per drive"
    return collect(labels)
end


const _PULSE_ANIMATION_COLORS = Makie.wong_colors()


function _pulse_animation_colors(nd::Int)
    return [_PULSE_ANIMATION_COLORS[mod1(i, length(_PULSE_ANIMATION_COLORS))] for i = 1:nd]
end


function _animation_neutral()
    try
        theme = Makie.current_default_theme()
        if haskey(theme, :textcolor)
            return Makie.to_value(theme[:textcolor])
        end
    catch
    end
    return :black
end


function _frame_label(title::AbstractString, frame_text::AbstractString)
    return isempty(frame_text) ? title : "$title\n$frame_text"
end


function Piccolo.animate_pulse(
    pulses::AbstractVector{<:Piccolo.AbstractPulse};
    fps::Int = 12,
    mode::Symbol = :inline,
    filename = "pulse_sweep_animation.mp4",
    n_samples::Int = 240,
    layout::Symbol = :stacked,
    title::AbstractString = "Pulse parameter sweep",
    labels = nothing,
    parameter_values = nothing,
    parameter_label::AbstractString = "Frame",
)
    @assert !isempty(pulses) "pulses cannot be empty"
    @assert n_samples >= 2 "n_samples must be at least 2"
    @assert layout in (:stacked, :overlay) "layout must be :stacked or :overlay"

    nd = Piccolo.n_drives(first(pulses))
    @assert all(Piccolo.n_drives(p) == nd for p in pulses) "all pulses must have the same number of drives"
    if !isnothing(parameter_values)
        @assert length(parameter_values) == length(pulses) "parameter_values must match the number of pulses"
    end

    sampled = [_pulse_sample_data(pulse, n_samples) for pulse in pulses]
    controls_by_frame = [item[1] for item in sampled]
    times_by_frame = [item[2] for item in sampled]
    drive_labels = _pulse_drive_labels(first(pulses), labels)
    colors = _pulse_animation_colors(nd)
    x_limits = (minimum(first.(times_by_frame)), maximum(last.(times_by_frame)))
    fig = Figure(size = layout == :stacked ? (800, 120 + 150 * nd) : (850, 460))

    frame_text = isnothing(parameter_values) ? "$parameter_label 1" :
                 "$parameter_label = $(parameter_values[1])"
    header = Label(
        fig[0, 1],
        _frame_label(title, frame_text);
        fontsize = 18,
        font = :bold,
        tellwidth = false,
    )

    function drive_values(drive)
        return vcat([vec(controls[drive, :]) for controls in controls_by_frame]...)
    end

    y_limits_by_drive = [_pulse_y_limits(drive_values(drive)) for drive = 1:nd]
    y_limits_overlay = _pulse_y_limits(vcat([vec(controls) for controls in controls_by_frame]...))

    if layout == :stacked
        axes = Axis[]
        for drive = 1:nd
            ax = Axis(
                fig[drive, 1];
                xlabel = drive == nd ? "Time" : "",
                xticklabelsvisible = drive == nd,
                ylabel = drive_labels[drive],
                titlealign = :left,
            )
            push!(axes, ax)
        end
    else
        ax = Axis(
            fig[1, 1];
            xlabel = "Time",
            ylabel = "Amplitude",
            titlealign = :left,
        )
        Legend(
            fig[1, 2],
            [LineElement(; color = colors[i], linewidth = 2) for i = 1:nd],
            drive_labels;
            tellheight = false,
        )
    end

    function update_frame!(frame)
        idx = clamp(frame, 1, length(pulses))
        frame_text = isnothing(parameter_values) ? "$parameter_label $idx" :
                     "$parameter_label = $(parameter_values[idx])"
        header.text[] = _frame_label(title, frame_text)

        if layout == :stacked
            for drive = 1:nd
                empty!(axes[drive])
                hlines!(
                    axes[drive],
                    [0.0];
                    color = (_animation_neutral(), 0.35),
                    linestyle = :dash,
                    linewidth = 0.5,
                )
                Piccolo.plot_pulse!(
                    axes[drive],
                    pulses[idx];
                    drive_indices = [drive],
                    n_samples = n_samples,
                    show_knots = false,
                    colors = [colors[drive]],
                )
                xlims!(axes[drive], x_limits...)
                ylims!(axes[drive], y_limits_by_drive[drive]...)
            end
        else
            empty!(ax)
            hlines!(
                ax,
                [0.0];
                color = (_animation_neutral(), 0.35),
                linestyle = :dash,
                linewidth = 0.5,
            )
            Piccolo.plot_pulse!(
                ax,
                pulses[idx];
                n_samples = n_samples,
                show_knots = false,
                colors = colors,
            )
            xlims!(ax, x_limits...)
            ylims!(ax, y_limits_overlay...)
        end
    end

    update_frame!(1)

    return Piccolo.animate_figure(
        fig,
        collect(1:length(pulses)),
        update_frame!,
        mode = mode,
        fps = fps,
        filename = filename,
    )
end


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
        1:(traj.N),
        update_frame!,
        mode = mode,
        fps = fps,
        filename = filename,
    )
end


# ---------------------------------------------------------------------------- #
# LivePulsePlotCallback — implementation of stub in src/visualizations/live_callbacks.jl
# ---------------------------------------------------------------------------- #

struct _LivePulsePlotCallback{QT<:AbstractQuantumTrajectory} <: AbstractIntermediateCallback
    qtraj::QT
    trajectory::NamedTrajectory
    every::Int
    save_dir::Union{Nothing,String}
    title_prefix::String
end

function Piccolo.LivePulsePlotCallback(
    qtraj::AbstractQuantumTrajectory,
    traj::NamedTrajectory;
    every::Int = 1,
    save_dir::Union{Nothing,String} = nothing,
    title_prefix::String = "iter",
)
    every >= 1 || throw(ArgumentError("every must be >= 1, got $every"))
    if save_dir !== nothing
        mkpath(save_dir)
    end
    return _LivePulsePlotCallback(qtraj, traj, every, save_dir, title_prefix)
end

function (cb::_LivePulsePlotCallback)(primal::AbstractVector, iter::Integer)
    # Gate skipped iters before any allocation.
    if (iter % cb.every) != 0
        return true
    end

    traj = cb.trajectory
    data_dim = traj.dim * traj.N
    expected = data_dim + traj.global_dim
    if length(primal) != expected
        @warn "primal-vector length mismatch ($(length(primal)) vs $expected); " *
              "trajectory reconstruction skipped. (For MadNLP, this usually " *
              "means `fixed_variable_treatment` was set to `MakeParameter` " *
              "explicitly, overriding DTO's default RelaxBound coupling.)" maxlog = 1
        return true
    end

    if traj.global_dim > 0
        update!(traj, collect(view(primal, 1:expected)); type = :both)
    else
        update!(traj, collect(view(primal, 1:data_dim)); type = :data)
    end

    pulse = extract_pulse(cb.qtraj, traj)
    fig = plot_pulse(pulse; title = "$(cb.title_prefix) $iter")

    if cb.save_dir !== nothing
        Makie.save(joinpath(cb.save_dir, "iter_$(lpad(iter, 5, '0')).png"), fig)
    end
    return true
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


@testitem "LivePulsePlotCallback subtype + direct invocation" begin
    using DirectTrajOpt
    using NamedTrajectories
    using CairoMakie
    using Random

    Random.seed!(42)
    sys = QuantumSystem(0.5 * PAULIS[:Z], [PAULIS[:X], PAULIS[:Y]], [1.0, 1.0])
    T, N = 10.0, 30
    times = collect(range(0, T, length = N))
    pulse = ZeroOrderPulse(0.1 * randn(2, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, GATES[:X])

    opts = PiccoloOptions(timesteps_all_equal = true, verbose = false)
    qcp = SmoothPulseProblem(
        qtraj,
        N;
        Q = 100.0,
        R = 1e-2,
        ddu_bound = 1.0,
        piccolo_options = opts,
    )

    save_dir = mktempdir()
    cb = LivePulsePlotCallback(qtraj, qcp.prob.trajectory; save_dir = save_dir)

    # Type contract: the callback subtypes DTO's solver-agnostic abstract type,
    # so any DTO backend can install it via its intermediate_callback option.
    @test cb isa DirectTrajOpt.AbstractIntermediateCallback

    # Direct invocation with a synthesized primal vector of the right length.
    traj = qcp.prob.trajectory
    expected = traj.dim * traj.N + traj.global_dim
    @test cb(randn(expected), 0) === true
    @test cb(randn(expected), 1) === true

    pngs = filter(f -> endswith(f, ".png"), readdir(save_dir))
    @test length(pngs) == 2

    # Each call produces a distinct frame.
    @test stat(joinpath(save_dir, pngs[1])).size != stat(joinpath(save_dir, pngs[2])).size
end

@testitem "LivePulsePlotCallback every > 1 skips renders" begin
    using DirectTrajOpt
    using NamedTrajectories
    using CairoMakie
    using Random

    Random.seed!(7)
    sys = QuantumSystem(0.5 * PAULIS[:Z], [PAULIS[:X]], [1.0])
    times = collect(range(0, 5.0, length = 20))
    qtraj = UnitaryTrajectory(sys, ZeroOrderPulse(0.1 * randn(1, 20), times), GATES[:X])
    qcp = SmoothPulseProblem(
        qtraj,
        20;
        piccolo_options = PiccoloOptions(timesteps_all_equal = true, verbose = false),
    )

    save_dir = mktempdir()
    cb = LivePulsePlotCallback(qtraj, qcp.prob.trajectory; every = 3, save_dir = save_dir)

    traj = qcp.prob.trajectory
    expected = traj.dim * traj.N + traj.global_dim
    for iter = 0:10
        cb(randn(expected), iter)
    end

    # iter % 3 == 0 fires on iters 0, 3, 6, 9 → 4 PNGs out of 11 invocations.
    pngs = filter(f -> endswith(f, ".png"), readdir(save_dir))
    @test length(pngs) == 4
end

@testitem "LivePulsePlotCallback rejects every <= 0" begin
    using NamedTrajectories
    using CairoMakie

    sys = QuantumSystem(0.5 * PAULIS[:Z], [PAULIS[:X]], [1.0])
    times = collect(range(0, 1.0, length = 5))
    qtraj = UnitaryTrajectory(sys, ZeroOrderPulse(zeros(1, 5), times), GATES[:X])
    qcp = SmoothPulseProblem(
        qtraj,
        5;
        piccolo_options = PiccoloOptions(timesteps_all_equal = true, verbose = false),
    )

    @test_throws ArgumentError LivePulsePlotCallback(qtraj, qcp.prob.trajectory; every = 0)
    @test_throws ArgumentError LivePulsePlotCallback(qtraj, qcp.prob.trajectory; every = -1)
end

@testitem "Test animate_pulse for pulse sweeps" begin
    using CairoMakie
    using Piccolo

    pulses = [GaussianPulse([amp, 0.5 * amp], 0.2, 1.0) for amp in range(0.2, 0.6, length = 3)]
    fig_sweep = animate_pulse(
        pulses;
        mode = :inline,
        fps = 10,
        n_samples = 8,
        parameter_values = collect(range(0.2, 0.6, length = 3)),
        parameter_label = "amplitude",
    )
    @test fig_sweep isa Figure
end


end
