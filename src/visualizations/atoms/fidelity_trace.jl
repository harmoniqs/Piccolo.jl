module FidelityTracePlots

export plot_fidelity_trace

using Makie
using LinearAlgebra
using NamedTrajectories
using Piccolo
using TestItems

"""
    plot_fidelity_trace(
        traj::NamedTrajectory,
        U_goal::AbstractMatrix{<:Number};
        unitary_name::Symbol = :Ũ⃗,
        subspace::Union{Nothing, AbstractVector{Int}} = nothing,
        title::AbstractString = "",
        kwargs...
    )

Plot instantaneous gate fidelity ``F(t) = |\\text{tr}(U_{\\text{goal}}^\\dagger U(t))|^2 / n^2`` vs time.

# Arguments
- `traj`: Trajectory containing a unitary in isomorphism representation.
- `U_goal`: Target unitary matrix.

# Keyword Arguments
- `unitary_name`: Symbol for the unitary iso-vec. Default `:Ũ⃗`.
- `subspace`: Compute fidelity only on this subspace. Default `nothing` (full space).
- `title`: Optional super-title.
- `kwargs...`: Passed to `Figure(; kwargs...)`.

# Returns
A Makie `Figure` with a single fidelity-vs-time line plot.
"""
function plot_fidelity_trace(
    traj::NamedTrajectory,
    U_goal::AbstractMatrix{<:Number};
    unitary_name::Symbol = :Ũ⃗,
    subspace::Union{Nothing,AbstractVector{Int}} = nothing,
    title::AbstractString = "",
    kwargs...,
)
    Ũ⃗_data = traj[unitary_name]
    T = size(Ũ⃗_data, 2)
    times = get_times(traj)

    fidelities = zeros(Float64, T)
    for t in 1:T
        U_t = iso_vec_to_operator(Ũ⃗_data[:, t])
        if isnothing(subspace)
            fidelities[t] = unitary_fidelity(U_t, U_goal)
        else
            fidelities[t] = unitary_fidelity(U_t, U_goal; subspace = subspace)
        end
    end

    fig = Figure(; size = (800, 350), kwargs...)
    ax = Axis(fig[1, 1], xlabel = "Time (μs)", ylabel = "Fidelity")

    lines!(ax, times, fidelities; linewidth = 2, color = :black)
    hlines!(ax, [1.0]; color = :green, linestyle = :dash, linewidth = 1.0, alpha = 0.7)

    ylims!(ax, max(0.0, minimum(fidelities) - 0.05), 1.02)

    if !isempty(title)
        Label(fig[0, 1], title; fontsize = 16, font = :bold)
    end

    return fig
end

# ============================================================================ #

@testitem "Fidelity trace with identity" begin
    using CairoMakie
    using NamedTrajectories
    using Piccolo
    using LinearAlgebra

    T = 20
    U = Matrix{ComplexF64}(I, 2, 2)
    Ũ⃗_data = hcat([operator_to_iso_vec(U) for _ in 1:T]...)

    traj = NamedTrajectory(
        (Ũ⃗ = Ũ⃗_data, u = randn(2, T), Δt = fill(0.1, T));
        controls = :u,
        timestep = :Δt,
    )

    fig = plot_fidelity_trace(traj, U)
    @test fig isa Figure
end

@testitem "Fidelity trace with X gate target" begin
    using CairoMakie
    using NamedTrajectories
    using Piccolo
    using LinearAlgebra

    T = 15
    X = ComplexF64[0 1; 1 0]
    Us = [exp(-im * (π / 2) * t / T * X) for t in 1:T]
    Ũ⃗_data = hcat(operator_to_iso_vec.(Us)...)

    traj = NamedTrajectory(
        (Ũ⃗ = Ũ⃗_data, u = randn(1, T), Δt = fill(0.1, T));
        controls = :u,
        timestep = :Δt,
    )

    fig = plot_fidelity_trace(traj, X; title = "X Gate Fidelity")
    @test fig isa Figure
end

end
