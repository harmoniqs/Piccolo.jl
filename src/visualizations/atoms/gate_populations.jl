module GatePopulationPlots

export plot_gate_populations

using Makie
using LinearAlgebra
using NamedTrajectories
using Piccolo
using TestItems

"""
    basis_labels(n_qubits::Int)

Generate computational basis state labels: |0⟩, |1⟩, |00⟩, |01⟩, ...
"""
function basis_labels(n_qubits::Int)
    n = 2^n_qubits
    return ["|" * string(i - 1, base = 2, pad = n_qubits) * "⟩" for i in 1:n]
end

"""
    plot_gate_populations(
        traj::NamedTrajectory;
        unitary_name::Symbol = :Ũ⃗,
        columns::Union{Nothing, AbstractVector{Int}} = nothing,
        title::AbstractString = "",
        kwargs...
    )

Plot population dynamics ``|U_{ij}(t)|^2`` in a grid layout with computational basis labels.

Each subplot shows the evolution of all output state populations for a given input state
(column of the unitary).

# Keyword Arguments
- `unitary_name`: Symbol for the unitary iso-vec. Default `:Ũ⃗`.
- `columns`: Which columns of ``U`` to plot. Default `nothing` (all columns).
- `title`: Optional super-title.
- `kwargs...`: Passed to `Figure(; kwargs...)`.

# Returns
A Makie `Figure`.
"""
function plot_gate_populations(
    traj::NamedTrajectory;
    unitary_name::Symbol = :Ũ⃗,
    columns::Union{Nothing,AbstractVector{Int}} = nothing,
    title::AbstractString = "",
    kwargs...,
)
    Ũ⃗_data = traj[unitary_name]
    T = size(Ũ⃗_data, 2)
    times = get_times(traj)

    n_iso = size(Ũ⃗_data, 1)
    n = Int(sqrt(n_iso ÷ 2))
    n_qubits = Int(log2(n))

    cols = isnothing(columns) ? collect(1:n) : columns
    n_panels = length(cols)

    # Grid layout
    if n_panels <= 2
        grid_cols = n_panels
        grid_rows = 1
    elseif n_panels <= 4
        grid_cols = 2
        grid_rows = 2
    else
        grid_cols = 4
        grid_rows = cld(n_panels, 4)
    end

    labels = basis_labels(n_qubits)

    fig_w = 320 * grid_cols + 80
    fig_h = 240 * grid_rows + (isempty(title) ? 20 : 50)
    fig = Figure(; size = (fig_w, fig_h), kwargs...)

    # Precompute populations: pops[i, j, t] = |U_{ij}(t)|²
    pops = zeros(Float64, n, n, T)
    for t in 1:T
        U = iso_vec_to_operator(Ũ⃗_data[:, t])
        for j in 1:n, i in 1:n
            pops[i, j, t] = abs2(U[i, j])
        end
    end

    colors = Makie.wong_colors()

    for (idx, j) in enumerate(cols)
        row = cld(idx, grid_cols)
        col = mod1(idx, grid_cols)

        ax = Axis(
            fig[row, col];
            title = "Input: $(labels[j])",
            titlesize = 13,
            ylabel = col == 1 ? "Population" : "",
            xlabel = row == grid_rows ? "Time (μs)" : "",
        )
        ylims!(ax, -0.05, 1.05)

        for i in 1:n
            c = colors[mod1(i, length(colors))]
            lines!(ax, times, pops[i, j, :]; linewidth = 1.8, color = c, label = labels[i])
        end

        if row < grid_rows
            hidexdecorations!(ax; grid = false)
        end
        if col > 1
            hideydecorations!(ax; grid = false)
        end
    end

    # Shared legend to the right of the grid
    legend_entries = [
        LineElement(; color = colors[mod1(i, length(colors))], linewidth = 2) for i in 1:n
    ]
    Legend(fig[:, grid_cols + 1], legend_entries, labels; framevisible = false, labelsize = 12)

    if !isempty(title)
        Label(fig[0, 1:grid_cols], title; fontsize = 16, font = :bold)
    end

    return fig
end

# ============================================================================ #

@testitem "Gate populations 1Q" begin
    using CairoMakie
    using NamedTrajectories
    using Piccolo
    using LinearAlgebra

    T = 15
    dim = 2
    U = Matrix{ComplexF64}(I, dim, dim)
    Ũ⃗_data = hcat([operator_to_iso_vec(U) for _ in 1:T]...)

    traj = NamedTrajectory(
        (Ũ⃗ = Ũ⃗_data, u = randn(2, T), Δt = fill(0.1, T));
        controls = :u,
        timestep = :Δt,
    )

    fig = plot_gate_populations(traj)
    @test fig isa Figure
end

@testitem "Gate populations 2Q" begin
    using CairoMakie
    using NamedTrajectories
    using Piccolo
    using LinearAlgebra

    T = 20
    dim = 4
    U = Matrix{ComplexF64}(I, dim, dim)
    Ũ⃗_data = hcat([operator_to_iso_vec(U) for _ in 1:T]...)

    traj = NamedTrajectory(
        (Ũ⃗ = Ũ⃗_data, u = randn(3, T), Δt = fill(0.1, T));
        controls = :u,
        timestep = :Δt,
    )

    fig = plot_gate_populations(traj; title = "2Q Identity")
    @test fig isa Figure
end

@testitem "Gate populations selected columns" begin
    using CairoMakie
    using NamedTrajectories
    using Piccolo
    using LinearAlgebra

    T = 15
    dim = 4
    U = Matrix{ComplexF64}(I, dim, dim)
    Ũ⃗_data = hcat([operator_to_iso_vec(U) for _ in 1:T]...)

    traj = NamedTrajectory(
        (Ũ⃗ = Ũ⃗_data, u = randn(3, T), Δt = fill(0.1, T));
        controls = :u,
        timestep = :Δt,
    )

    fig = plot_gate_populations(traj; columns = [1, 4])
    @test fig isa Figure
end

end
