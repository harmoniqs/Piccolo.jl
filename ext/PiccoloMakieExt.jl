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
    filename::String = "animation.mp4"
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
        record(fig, filename, frames; framerate=fps) do i
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
    fps::Int=24,
    mode::Symbol=:inline,
    filename="name_animation.mp4",
    kwargs...
)
    fig = Figure()
    ax = Axis(fig[1, 1])

    # Set limits
    times = get_times(traj)
    xlims = [times[1], times[end]]
    ylims = collect(extrema(traj[name]))

    scatter!(ax, xlims, ylims, markersize=0.0)
    QuantumObjectPlots.plot_name!(ax, traj, name, indices=1:1, kwargs...)

    # TODO: Unclear how to set observables via indices, so redraw
    # WARNING: Known bug with CairoMakie, which cannot re-render
    function update_frame!(i)
        empty!(ax.scene)
        scatter!(ax, xlims, ylims, markersize=0.0)
        QuantumObjectPlots.plot_name!(ax, traj, name, indices=1:i)
    end

    return Piccolo.animate_figure(
        fig,
        1:traj.N,
        update_frame!,
        mode=mode,
        fps=fps,
        filename=filename
    )
end


@testitem "Test animate_name for NamedTrajectory" begin
    using QuantumToolbox
    using NamedTrajectories
    using CairoMakie

    x = ComplexF64[1.0; 0.0]
    y = ComplexF64[0.0, 1.0]
    
    comps = (ψ̃ = hcat(ket_to_iso(x), ket_to_iso(y)), Δt = [1.0; 1.0],)
    traj = NamedTrajectory(comps)

    fig = animate_name(traj, :ψ̃, mode=:inline, fps=10)
    @test fig isa Figure
end


end
