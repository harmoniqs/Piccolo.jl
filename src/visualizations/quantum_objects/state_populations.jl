"""
    plot_state_populations(
        traj::NamedTrajectory;
        state_name::Symbol=:ψ̃,
        state_indices::Union{Nothing, AbstractVector{Int}}=nothing,
        control_name::Symbol=:u,
        subspace::Union{Nothing, AbstractVector{Int}}=nothing,
        kwargs...
    )

Plot populations for multiple quantum states stored in a trajectory.

This function visualizes the time evolution of quantum state populations for trajectories
containing multiple state trajectories (e.g., from sampling problems where multiple initial
states are evolved). States are identified by a common prefix and numeric suffix pattern
(e.g., `:ψ̃1_system_1`, `:ψ̃2_system_1`, etc.).

# Mathematical Background

For a quantum state ``|\\psi(t)\\rangle \\in \\mathcal{H}`` evolving under the Schrödinger equation,
the population in computational basis state ``|i\\rangle`` is given by

```math
P_i(t) = |\\langle i|\\psi(t)\\rangle|^2 = |\\psi_i(t)|^2
```

where ``\\psi_i(t)`` is the ``i``-th component of the state vector in the computational basis.

For a normalized state, we have the conservation property:

```math
\\sum_{i=1}^{n} P_i(t) = \\langle\\psi(t)|\\psi(t)\\rangle = 1
```

When multiple states are being optimized simultaneously (as in sampling problems), this function
plots the populations for each state, allowing comparison of how different initial conditions evolve
under the same control protocol.

**Key Properties:**
- Populations are real and bounded: ``P_i(t) \\in [0,1]``
- Total probability is conserved: ``\\sum_i P_i(t) = 1``
- Can optionally restrict to a subspace (e.g., computational subspace excluding leakage states)

# Arguments
- `traj::NamedTrajectory`: A trajectory containing one or more quantum states in isomorphism representation.

# Keyword Arguments
- `state_name::Symbol`: The base name for state components. The function will find all trajectory
  components matching this prefix (e.g., `:ψ̃` matches `:ψ̃1_system_1`, `:ψ̃2_system_1`, etc.).
  Default is `:ψ̃`.
- `state_indices::Union{Nothing, AbstractVector{Int}}`: If provided, only plot states with these
  indices (e.g., `[1, 3]` plots only the 1st and 3rd states). If `nothing`, plots all states
  matching the prefix. Default is `nothing`.
- `control_name::Symbol`: The name of the control signal component to include in the plot.
  Default is `:u`.
- `subspace::Union{Nothing, AbstractVector{Int}}`: If provided, only plot populations for these
  basis states (e.g., `1:2` for a qubit subspace). Useful for excluding leakage levels.
  Default is `nothing` (plot all levels).
- `kwargs...`: Additional keyword arguments passed to [`NamedTrajectories.plot`](https://docs.harmoniqs.co/NamedTrajectories/dev/generated/plotting/).

# Returns
- A Makie `Figure` object containing the population plots for all selected states.

# Example
```julia
using NamedTrajectories
using Piccolo
using Piccolo

# Example: Two initial states evolving under the same controls
N = 100
Δt = 0.1

# Initial states
ψ1_init = ComplexF64[1, 0, 0]  # |0⟩
ψ2_init = ComplexF64[0, 1, 0]  # |1⟩

# Create trajectory with multiple states
traj = NamedTrajectory(
    (
        ψ̃1_system_1 = hcat([ket_to_iso(ψ1_init) for _ in 1:N]...),
        ψ̃2_system_1 = hcat([ket_to_iso(ψ2_init) for _ in 1:N]...),
        u = randn(2, N),
        Δt = fill(Δt, N),
    );
    controls = :u,
    timestep = :Δt,
)

# Plot populations for all states
plot_state_populations(traj)

# Plot only computational subspace (excluding 3rd level)
plot_state_populations(traj; subspace=1:2)

# Plot only first state
plot_state_populations(traj; state_indices=[1])
```

See also: [`plot_unitary_populations`](@ref), [`NamedTrajectories.plot`](https://docs.harmoniqs.co/NamedTrajectories/dev/generated/plotting/)
"""
function plot_state_populations(
    traj::NamedTrajectory;
    state_name::Symbol = :ψ̃,
    state_indices::Union{Nothing,AbstractVector{Int}} = nothing,
    control_name::Symbol = :u,
    subspace::Union{Nothing,AbstractVector{Int}} = nothing,
    kwargs...,
)
    # Find all state components matching the prefix
    state_prefix = string(state_name)
    matching_states = [
        name for name in traj.names if
        startswith(string(name), state_prefix) && string(name) != state_prefix
    ]

    if isempty(matching_states)
        # Fallback: maybe there's just a single state with exact name
        if state_name in traj.names
            matching_states = [state_name]
        else
            error("No state components found matching prefix :$state_name")
        end
    end

    # Sort states to ensure consistent ordering
    sort!(matching_states)

    # Filter by indices if specified
    if !isnothing(state_indices)
        matching_states = matching_states[state_indices]
    end

    # Create transformations for each state
    transformations = Pair{Symbol,Function}[]
    transformation_labels = String[]
    transformation_titles = String[]

    for (idx, sname) in enumerate(matching_states)
        # Define population extraction function
        pop_fn = if isnothing(subspace)
            x -> abs2.(iso_to_ket(x))
        else
            x -> abs2.(iso_to_ket(x)[subspace])
        end

        push!(transformations, (sname => pop_fn))
        push!(transformation_labels, L"P")

        # Create title with state name
        state_str = string(sname)
        push!(transformation_titles, L"Populations: $|%$(state_str)(t)|^2$")
    end

    # Plot using NamedTrajectories.plot
    plot(
        traj,
        [control_name];
        transformations = transformations,
        transformation_titles = transformation_titles,
        transformation_labels = transformation_labels,
        kwargs...,
    )
end

# ============================================================================ #

@testitem "Plot state populations with multiple states" begin
    using CairoMakie
    using NamedTrajectories
    using Piccolo

    # Setup: Two states evolving (sampling problem style)
    N = 30
    Δt = 0.1

    ψ1_init = ComplexF64[1, 0]
    ψ2_init = ComplexF64[0, 1]

    # Create dummy evolution (just rotate slightly each step)
    θs = range(0, π/4, length = N)
    ψ1s = [ComplexF64[cos(θ), sin(θ)] for θ in θs]
    ψ2s = [ComplexF64[sin(θ), cos(θ)] for θ in θs]

    traj = NamedTrajectory(
        (
            ψ̃1_system_1 = hcat(ket_to_iso.(ψ1s)...),
            ψ̃2_system_1 = hcat(ket_to_iso.(ψ2s)...),
            u = randn(2, N),
            Δt = fill(Δt, N),
        );
        controls = :u,
        timestep = :Δt,
    )

    # Test: Default behavior plots all states
    fig = plot_state_populations(traj)
    @test fig isa Figure
    @test length(fig.content) > 0
end

@testitem "Plot state populations with state index filtering" begin
    using CairoMakie
    using NamedTrajectories
    using Piccolo

    # Setup: Three states
    N = 20
    Δt = 0.1

    ψ_data = [ket_to_iso(ComplexF64[1, 0]) for _ = 1:N]

    traj = NamedTrajectory(
        (
            ψ̃1_system_1 = hcat(ψ_data...),
            ψ̃2_system_1 = hcat(ψ_data...),
            ψ̃3_system_1 = hcat(ψ_data...),
            u = zeros(1, N),
            Δt = fill(Δt, N),
        );
        controls = :u,
        timestep = :Δt,
    )

    # Test: Plot only selected states
    fig1 = plot_state_populations(traj; state_indices = [1])
    @test fig1 isa Figure

    fig2 = plot_state_populations(traj; state_indices = [1, 3])
    @test fig2 isa Figure
end

@testitem "Plot state populations with subspace restriction" begin
    using CairoMakie
    using NamedTrajectories
    using Piccolo

    # Setup: 3-level system (qubit + leakage)
    N = 15
    Δt = 0.05

    ψ_init = ComplexF64[1, 0, 0]  # Start in |0⟩
    ψs = [ψ_init for _ = 1:N]

    traj = NamedTrajectory(
        (ψ̃1_system_1 = hcat(ket_to_iso.(ψs)...), u = zeros(1, N), Δt = fill(Δt, N));
        controls = :u,
        timestep = :Δt,
    )

    # Test: Plot only computational subspace (exclude leakage)
    fig = plot_state_populations(traj; subspace = 1:2)
    @test fig isa Figure
end

@testitem "Plot state populations verifies normalization" begin
    using CairoMakie
    using NamedTrajectories
    using Piccolo
    using LinearAlgebra

    # Setup: Verify populations sum to 1
    N = 25
    Δt = 0.1

    # Create normalized states
    θs = range(0, π, length = N)
    ψs = [normalize(ComplexF64[cos(θ), sin(θ), 0.0]) for θ in θs]

    traj = NamedTrajectory(
        (ψ̃1_system_1 = hcat(ket_to_iso.(ψs)...), u = zeros(1, N), Δt = fill(Δt, N));
        controls = :u,
        timestep = :Δt,
    )

    # Verify normalization at each knot point
    for k = 1:N
        ψ = iso_to_ket(traj[:ψ̃1_system_1][:, k])
        @test norm(ψ) ≈ 1.0 atol=1e-10

        pops = abs2.(ψ)
        @test sum(pops) ≈ 1.0 atol=1e-10
    end

    # Plot still works
    fig = plot_state_populations(traj)
    @test fig isa Figure
end

@testitem "Plot single state with exact name match" begin
    using CairoMakie
    using NamedTrajectories
    using Piccolo

    # Setup: Single state without numeric suffix
    N = 20
    Δt = 0.1

    ψs = [ComplexF64[1, 0] for _ = 1:N]

    traj = NamedTrajectory(
        (ψ̃ = hcat(ket_to_iso.(ψs)...), u = zeros(1, N), Δt = fill(Δt, N));
        controls = :u,
        timestep = :Δt,
    )

    # Test: Should work with exact name match
    fig = plot_state_populations(traj; state_name = :ψ̃)
    @test fig isa Figure
end
