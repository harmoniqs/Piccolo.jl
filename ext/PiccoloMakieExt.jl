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


end
