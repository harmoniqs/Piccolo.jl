module RydbergChainPlots

export plot_rydberg_chain
export animate_rydberg_chain

using Makie
import Makie: plot
using LinearAlgebra
using NamedTrajectories
using Piccolo
using TestItems

"""
    rydberg_populations(ψ::AbstractVector{<:Number}, N::Int)

Calculate the population in the Rydberg state (|1⟩) for each atom in a 1D chain of N atoms.
"""
function rydberg_populations(ψ::AbstractVector{<:Number}, N::Int)
    pops = zeros(Float64, N)
    for idx = 0:(2^N-1)
        prob = abs2(ψ[idx+1])
        for i = 1:N
            b_i = (idx ÷ 2^(N - i)) % 2
            if b_i == 1
                pops[i] += prob
            end
        end
    end
    return pops
end

"""
    plot_rydberg_chain(
        N::Int, 
        distance::Real; 
        C::Real = 862690 * 2π, 
        cutoff_order::Int = 1,
        populations::Union{Nothing, Observable, AbstractVector} = nothing,
        kwargs...
    )

Plots a 1D chain of `N` neutral atoms with spacing `distance` (μm).
Connects interacting atoms up to `cutoff_order` neighbors with thickness proportional to the Rydberg interaction strength `C / r^6`.
`populations` can be an Observable or a static Vector of values in [0, 1], representing the Rydberg state (excited) population of each atom. 
If `populations` are provided, the atom colors map from ground state to Rydberg state.
"""
function plot_rydberg_chain(
    N::Int,
    distance::Real;
    C::Real = 862690 * 2π,
    cutoff_order::Int = 1,
    populations::Union{Nothing,Observable,AbstractVector} = nothing,
    kwargs...,
)
    fig = Figure(; kwargs...)
    ax = Axis(
        fig[1, 1],
        xlabel = "Position (μm)",
        ylabel = "",
        yticksvisible = false,
        yticklabelsvisible = false,
    )
    hidespines!(ax, :t, :r, :l, :b)
    hidedecorations!(ax, label = false)

    positions = [distance * (i - 1) for i = 1:N]

    # Draw interactions
    for gap = 1:cutoff_order
        for i = 1:(N-gap)
            j = i + gap
            r = distance * gap
            interaction = C / r^6
            # Normalize thickness relative to nearest neighbor
            lw = (interaction / (C / distance^6)) * 5.0

            lines!(
                ax,
                [positions[i], positions[j]],
                [0.0, 0.0],
                linewidth = lw,
                color = :black,
                alpha = 0.3,
            )

            # Draw arcs for next-nearest neighbors to make them visible
            if gap > 1
                arc_x = range(positions[i], positions[j], length = 50)
                # Height proportional to gap length
                h = 0.5 * distance * gap
                arc_y = [
                    h * (1 - (2 * (x - (positions[i] + positions[j]) / 2) / (r))^2) for
                    x in arc_x
                ]
                lines!(
                    ax,
                    collect(arc_x),
                    arc_y,
                    linewidth = lw,
                    color = :black,
                    alpha = 0.3,
                )
            end
        end
    end

    # Handle populations
    pops_obs = if populations isa Observable
        populations
    elseif populations !== nothing
        Observable(populations)
    else
        Observable(zeros(Float64, N))
    end

    scatter!(
        ax,
        positions,
        zeros(N),
        markersize = 30,
        color = pops_obs,
        colormap = :coolwarm,
        colorrange = (0.0, 1.0),
    )

    ylims!(ax, -distance, distance * max(1, cutoff_order) * 0.6)
    xlims!(ax, -distance * 0.5, distance * (N - 0.5))

    return fig
end

"""
    animate_rydberg_chain(
        trajectory::NamedTrajectory, 
        N::Int, 
        distance::Real,
        filename::String = "rydberg_chain.gif"; 
        C::Real = 862690 * 2π, 
        cutoff_order::Int = 1,
        state_name::Symbol = :ψ̃,
        framerate::Int = 30,
        kwargs...
    )

Animate the excitation pattern of a neutral atom chain based on a gate trajectory.
Extracts the population of the Rydberg state at each timestep and updates the visualization.
Saves the animation to `filename`.
"""
function animate_rydberg_chain(
    trajectory::NamedTrajectory,
    N::Int,
    distance::Real,
    filename::String = "rydberg_chain.gif";
    C::Real = 862690 * 2π,
    cutoff_order::Int = 1,
    state_name::Symbol = :ψ̃,
    framerate::Int = 30,
    kwargs...,
)
    # Find the corresponding state vector
    state_prefix = string(state_name)
    matching_states =
        [name for name in trajectory.names if startswith(string(name), state_prefix)]

    trajectory_state_name = isempty(matching_states) ? state_name : first(matching_states)

    if !(trajectory_state_name in trajectory.names)
        error("State `\$(trajectory_state_name)` not found in trajectory!")
    end

    state_data = trajectory[trajectory_state_name]
    num_timesteps = size(state_data, 2)
    pops_obs = Observable(zeros(Float64, N))

    fig = plot_rydberg_chain(
        N,
        distance;
        C = C,
        cutoff_order = cutoff_order,
        populations = pops_obs,
        kwargs...,
    )

    record(fig, filename, 1:num_timesteps; framerate = framerate) do t
        # Convert iso to ket
        ket = iso_to_ket(state_data[:, t])
        pops_obs[] = rydberg_populations(ket, N)
    end

    return fig
end


@testitem "Rydberg chain plots" begin
    using CairoMakie
    using NamedTrajectories
    using Piccolo

    # Test static plot
    fig = plot_rydberg_chain(3, 5.0)
    @test fig isa Figure

    # Setup for animation test
    N = 2
    T = 10
    iso_size = 2 * (2^N)

    # Mock some data
    ψ_data = randn(iso_size, T)
    traj = NamedTrajectory(
        (ψ̃ = ψ_data, u = randn(1, T), Δt = fill(0.1, T));
        controls = :u,
        timestep = :Δt,
    )

    # Test animation creation (use .gif since CairoMakie supports it without FFMPEG)
    fig_anim = animate_rydberg_chain(traj, N, 5.0, "test_rydberg.gif")
    @test fig_anim isa Figure
    @test isfile("test_rydberg.gif")

    # Clean up
    rm("test_rydberg.gif", force = true)
end

end
