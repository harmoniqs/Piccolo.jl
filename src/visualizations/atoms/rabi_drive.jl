module RabiDrivePlots

export plot_rabi_drive

using Makie
using NamedTrajectories
using TestItems

"""
    plot_rabi_drive(
        traj::NamedTrajectory;
        control_name::Symbol = :u,
        omega_x_index::Int = 1,
        omega_y_index::Int = 2,
        title::AbstractString = "",
        kwargs...
    )

Plot the Rabi drive in polar coordinates: magnitude `|Ω(t)| = √(Ωx² + Ωy²)`
and phase `φ(t) = atan(Ωy, Ωx)`.

# Keyword Arguments
- `control_name`: Symbol for the control component. Default `:u`.
- `omega_x_index`: Row index of Ωx in the control matrix. Default `1`.
- `omega_y_index`: Row index of Ωy in the control matrix. Default `2`.
- `title`: Optional super-title.
- `kwargs...`: Passed to `Figure(; kwargs...)`.

# Returns
A Makie `Figure` with two subplots: magnitude (MHz) and phase (rad).
"""
function plot_rabi_drive(
    traj::NamedTrajectory;
    control_name::Symbol = :u,
    omega_x_index::Int = 1,
    omega_y_index::Int = 2,
    title::AbstractString = "",
    kwargs...,
)
    u = traj[control_name]
    times = get_times(traj)

    Ωx = u[omega_x_index, :]
    Ωy = u[omega_y_index, :]

    magnitude = sqrt.(Ωx .^ 2 .+ Ωy .^ 2)
    phase = atan.(Ωy, Ωx)

    fig = Figure(; size = (800, 450), kwargs...)

    # Magnitude subplot
    ax_mag = Axis(fig[1, 1], ylabel = "|Ω| (MHz)")
    lines!(ax_mag, times, magnitude; linewidth = 2, color = :dodgerblue)
    hlines!(ax_mag, [0.0]; color = :gray, linestyle = :dash, linewidth = 0.5)
    hidexdecorations!(ax_mag; grid = false)

    # Phase subplot
    ax_phase = Axis(
        fig[2, 1],
        ylabel = "φ (rad)",
        xlabel = "Time (μs)",
        yticks = ([-π, -π / 2, 0, π / 2, π], ["-π", "-π/2", "0", "π/2", "π"]),
    )
    lines!(ax_phase, times, phase; linewidth = 2, color = :coral)
    hlines!(ax_phase, [0.0]; color = :gray, linestyle = :dash, linewidth = 0.5)

    if !isempty(title)
        Label(fig[0, 1], title; fontsize = 16, font = :bold)
    end

    return fig
end

# ============================================================================ #

@testitem "Rabi drive plot basic" begin
    using CairoMakie
    using NamedTrajectories

    T = 30
    traj = NamedTrajectory(
        (u = randn(3, T), Δt = fill(0.1, T));
        controls = :u,
        timestep = :Δt,
    )

    fig = plot_rabi_drive(traj)
    @test fig isa Figure
end

@testitem "Rabi drive plot with custom indices" begin
    using CairoMakie
    using NamedTrajectories

    T = 20
    traj = NamedTrajectory(
        (u = randn(4, T), Δt = fill(0.05, T));
        controls = :u,
        timestep = :Δt,
    )

    fig = plot_rabi_drive(traj;
        omega_x_index = 2, omega_y_index = 3,
        title = "Custom Rabi Drive",
    )
    @test fig isa Figure
end

end
