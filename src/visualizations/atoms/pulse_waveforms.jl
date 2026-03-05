module PulseWaveformPlots

export plot_pulse_waveforms

using Makie
using NamedTrajectories
using TestItems

"""
    plot_pulse_waveforms(
        traj::NamedTrajectory;
        control_name::Symbol = :u,
        labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
        bounds::Union{Nothing, AbstractVector} = nothing,
        title::AbstractString = "",
        kwargs...
    )

Plot pulse amplitudes vs time with each drive in a vertically stacked subplot.

When `bounds` are provided, each subplot shows a shaded band and dashed lines
marking the hardware-feasible region.

# Keyword Arguments
- `control_name`: Symbol for the control component. Default `:u`.
- `labels`: Y-axis labels per drive (e.g., `["Ωx", "Ωy", "Δ"]`).
  Defaults to `"Drive 1"`, `"Drive 2"`, etc.
- `bounds`: Per-drive hardware bounds as `(lower, upper)` tuples.
  Pass `nothing` for a drive to skip its band.
- `title`: Optional super-title.
- `kwargs...`: Passed to `Figure(; kwargs...)`.

# Returns
A Makie `Figure`.
"""
function plot_pulse_waveforms(
    traj::NamedTrajectory;
    control_name::Symbol = :u,
    labels::Union{Nothing,AbstractVector{<:AbstractString}} = nothing,
    bounds::Union{Nothing,AbstractVector} = nothing,
    title::AbstractString = "",
    kwargs...,
)
    u = traj[control_name]
    n_drives, T = size(u)
    times = get_times(traj)

    if isnothing(labels)
        labels = ["Drive $i" for i in 1:n_drives]
    end

    fig = Figure(; size = (800, 100 + 140 * n_drives), kwargs...)

    for i in 1:n_drives
        ax = Axis(
            fig[i, 1],
            ylabel = labels[min(i, length(labels))] * " (MHz)",
            xlabel = i == n_drives ? "Time (μs)" : "",
        )

        # Hardware bounds as shaded band
        if !isnothing(bounds) && i <= length(bounds) && !isnothing(bounds[i])
            lo, hi = bounds[i]
            band!(ax, times, fill(Float64(lo), T), fill(Float64(hi), T);
                color = (:steelblue, 0.1))
            hlines!(ax, [lo, hi]; color = :steelblue, linestyle = :dash,
                linewidth = 0.8, alpha = 0.5)
        end

        lines!(ax, times, u[i, :]; linewidth = 2)
        hlines!(ax, [0.0]; color = :gray, linestyle = :dash, linewidth = 0.5)

        if i < n_drives
            hidexdecorations!(ax; grid = false)
        end
    end

    if !isempty(title)
        Label(fig[0, 1], title; fontsize = 16, font = :bold)
    end

    return fig
end

# ============================================================================ #

@testitem "Pulse waveform plot basic" begin
    using CairoMakie
    using NamedTrajectories

    T = 20
    traj = NamedTrajectory(
        (u = randn(3, T), Δt = fill(0.1, T));
        controls = :u,
        timestep = :Δt,
    )

    fig = plot_pulse_waveforms(traj)
    @test fig isa Figure
end

@testitem "Pulse waveform plot with labels and bounds" begin
    using CairoMakie
    using NamedTrajectories

    T = 20
    traj = NamedTrajectory(
        (u = randn(3, T), Δt = fill(0.1, T));
        controls = :u,
        timestep = :Δt,
    )

    fig = plot_pulse_waveforms(traj;
        labels = ["Ωx", "Ωy", "Δ"],
        bounds = [(-15.8, 15.8), (-15.8, 15.8), (-124.0, 124.0)],
        title = "Test Pulses",
    )
    @test fig isa Figure
end

end
