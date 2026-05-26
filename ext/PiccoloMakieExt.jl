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
    if save_dir !== nothing
        mkpath(save_dir)
    end
    return _LivePulsePlotCallback(qtraj, traj, every, save_dir, title_prefix)
end

function (cb::_LivePulsePlotCallback)(primal::AbstractVector, iter::Integer)
    traj = cb.trajectory
    data_dim = traj.dim * traj.N
    expected = data_dim + traj.global_dim
    if length(primal) != expected
        @warn "primal-vector length mismatch ($(length(primal)) vs $expected); " *
              "pass `fixed_variable_treatment = MadNLP.RelaxBound` to MadNLPOptions"
        return true
    end

    if traj.global_dim > 0
        update!(traj, collect(view(primal, 1:expected)); type = :both)
    else
        update!(traj, collect(view(primal, 1:data_dim)); type = :data)
    end

    if (iter % cb.every) != 0
        return true
    end

    pulse = extract_pulse(cb.qtraj, traj)
    fig = plot_pulse(pulse; title = "$(cb.title_prefix) $iter")

    if cb.save_dir !== nothing
        Makie.save(joinpath(cb.save_dir, "iter_$(lpad(iter, 3, '0')).png"), fig)
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


@testitem "LivePulsePlotCallback fires through MadNLP via abstract adapter" begin
    using DirectTrajOpt
    using NamedTrajectories
    using CairoMakie
    import MadNLP
    using Random

    Random.seed!(42)
    H_drift = 0.5 * PAULIS[:Z]
    H_drives = [PAULIS[:X], PAULIS[:Y]]
    sys = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

    T, N = 10.0, 30
    times = collect(range(0, T, length = N))
    pulse = ZeroOrderPulse(0.1 * randn(2, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, GATES[:X])

    opts = PiccoloOptions(timesteps_all_equal = true, verbose = false)
    qcp = SmoothPulseProblem(qtraj, N; Q = 100.0, R = 1e-2, ddu_bound = 1.0,
                             piccolo_options = opts)

    save_dir = mktempdir()
    cb = LivePulsePlotCallback(qtraj, qcp.prob.trajectory;
                               every = 1, save_dir = save_dir)
    @test cb isa DirectTrajOpt.AbstractIntermediateCallback

    DirectTrajOpt.solve!(qcp; options = DirectTrajOpt.MadNLPOptions(
        max_iter = 3,
        intermediate_callback = cb,
        fixed_variable_treatment = MadNLP.RelaxBound,
    ), verbose = false)

    pngs = filter(f -> endswith(f, ".png"), readdir(save_dir))
    @test !isempty(pngs)
end


end
