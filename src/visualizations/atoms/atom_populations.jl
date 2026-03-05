module AtomPopulationPlots

export plot_atom_populations

using Makie
using LinearAlgebra
using NamedTrajectories
using Piccolo
using TestItems

"""
    rydberg_populations(ψ::AbstractVector{<:Number}, N::Int)

Calculate the population in the Rydberg state (|1⟩) for each atom in a chain of N atoms.
"""
function rydberg_populations(ψ::AbstractVector{<:Number}, N::Int)
    pops = zeros(Float64, N)
    for idx in 0:(2^N - 1)
        prob = abs2(ψ[idx + 1])
        for i in 1:N
            b_i = (idx ÷ 2^(N - i)) % 2
            if b_i == 1
                pops[i] += prob
            end
        end
    end
    return pops
end

"""
    plot_atom_populations(
        traj::NamedTrajectory,
        N_atoms::Int;
        unitary_name::Symbol = :Ũ⃗,
        initial_states::Union{Nothing, AbstractVector{<:AbstractVector}} = nothing,
        title::AbstractString = "",
        kwargs...
    )

Plot per-atom Rydberg state population vs time.

Propagates the unitary trajectory ``U(t)`` on initial state(s) and computes the
single-atom excitation probability ``\\langle n_i \\rangle(t)`` at each timestep.

# Arguments
- `traj`: Trajectory containing a unitary in isomorphism representation.
- `N_atoms`: Number of atoms in the chain.

# Keyword Arguments
- `unitary_name`: Symbol for the unitary iso-vec. Default `:Ũ⃗`.
- `initial_states`: Vector of initial kets to propagate. Defaults to `|0…0⟩`.
- `title`: Optional super-title.
- `kwargs...`: Passed to `Figure(; kwargs...)`.

# Returns
A Makie `Figure` with one subplot per initial state, one line per atom.
"""
function plot_atom_populations(
    traj::NamedTrajectory,
    N_atoms::Int;
    unitary_name::Symbol = :Ũ⃗,
    initial_states::Union{Nothing,AbstractVector{<:AbstractVector}} = nothing,
    title::AbstractString = "",
    kwargs...,
)
    Ũ⃗_data = traj[unitary_name]
    T = size(Ũ⃗_data, 2)
    times = get_times(traj)
    dim = 2^N_atoms

    if isnothing(initial_states)
        ψ0 = zeros(ComplexF64, dim)
        ψ0[1] = 1.0
        initial_states = [ψ0]
    end

    n_states = length(initial_states)
    fig = Figure(; size = (800, 300 * n_states), kwargs...)

    for (s, ψ0) in enumerate(initial_states)
        ax = Axis(
            fig[s, 1],
            ylabel = "Rydberg population",
            xlabel = s == n_states ? "Time (μs)" : "",
        )

        pop_traces = zeros(Float64, N_atoms, T)
        for t in 1:T
            U = iso_vec_to_operator(Ũ⃗_data[:, t])
            ψ_t = U * ψ0
            pop_traces[:, t] = rydberg_populations(ψ_t, N_atoms)
        end

        for i in 1:N_atoms
            lines!(ax, times, pop_traces[i, :]; linewidth = 2, label = "Atom $i")
        end

        ylims!(ax, -0.05, 1.05)

        if N_atoms <= 6
            axislegend(ax; position = :rt)
        end

        if s < n_states
            hidexdecorations!(ax; grid = false)
        end
    end

    if !isempty(title)
        Label(fig[0, 1], title; fontsize = 16, font = :bold)
    end

    return fig
end

# ============================================================================ #

@testitem "Atom populations basic" begin
    using CairoMakie
    using NamedTrajectories
    using Piccolo
    using LinearAlgebra

    N = 2
    T = 15
    dim = 2^N
    U = Matrix{ComplexF64}(I, dim, dim)
    Ũ⃗_data = hcat([operator_to_iso_vec(U) for _ in 1:T]...)

    traj = NamedTrajectory(
        (Ũ⃗ = Ũ⃗_data, u = randn(3, T), Δt = fill(0.1, T));
        controls = :u,
        timestep = :Δt,
    )

    fig = plot_atom_populations(traj, N)
    @test fig isa Figure
end

@testitem "Atom populations with custom initial state" begin
    using CairoMakie
    using NamedTrajectories
    using Piccolo
    using LinearAlgebra

    N = 2
    T = 10
    dim = 2^N
    U = Matrix{ComplexF64}(I, dim, dim)
    Ũ⃗_data = hcat([operator_to_iso_vec(U) for _ in 1:T]...)

    traj = NamedTrajectory(
        (Ũ⃗ = Ũ⃗_data, u = randn(2, T), Δt = fill(0.1, T));
        controls = :u,
        timestep = :Δt,
    )

    ψ0 = ComplexF64[0, 0, 0, 1]
    fig = plot_atom_populations(traj, N; initial_states = [ψ0])
    @test fig isa Figure
end

end
