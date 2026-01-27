"""
    plot_unitary_populations(
        traj::NamedTrajectory;
        unitary_columns::AbstractVector{Int}=1:2,
        unitary_name::Symbol=:Ũ⃗,
        control_name::Symbol=:u,
        kwargs...
    )

Plot the state populations for specified columns of a unitary operator trajectory.

This function visualizes how the populations (squared magnitudes) of quantum states evolve
over time for selected columns of a unitary matrix stored in a `NamedTrajectory`.

# Mathematical Background

For a unitary operator ``U(t) \\in \\mathcal{U}(n)`` evolving under a time-dependent Hamiltonian,
this function plots the populations

```math
P_{i,j}(t) = |U_{i,j}(t)|^2
```

where ``U_{i,j}(t)`` is the ``(i,j)``-th element of the unitary matrix at time ``t``.

For a quantum system evolving according to the Schrödinger equation

```math
\\frac{d}{dt}U(t) = -iH(t)U(t), \\quad U(0) = I
```

each column ``j`` of ``U(t)`` represents the time evolution of the initial basis state ``|j\\rangle``:

```math
|\\psi_j(t)\\rangle = U(t)|j\\rangle = \\sum_{i=1}^{n} U_{i,j}(t)|i\\rangle
```

The population ``P_{i,j}(t) = |U_{i,j}(t)|^2`` gives the probability of finding the system in state ``|i\\rangle``
at time ``t`` given that it started in state ``|j\\rangle``.

**Key Properties:**
- Unitarity ensures ``\\sum_{i=1}^{n} P_{i,j}(t) = 1`` for all ``j`` and ``t`` (probability conservation)
- ``P_{i,j}(0) = \\delta_{ij}`` (initially in definite state)
- Populations are real and bounded: ``P_{i,j}(t) \\in [0,1]``

The trajectory stores ``U(t)`` in isomorphism representation ``\\tilde{U}(t)``, a vectorized form that
preserves the operator structure while enabling efficient optimization algorithms.

# Arguments
- `traj::NamedTrajectory`: A trajectory containing a unitary operator in isomorphism representation.

# Keyword Arguments
- `unitary_columns::AbstractVector{Int}`: Indices of unitary matrix columns to plot. Each column ``j``
  corresponds to the evolution of basis state ``|j\\rangle``. Default is `1:2`.
- `unitary_name::Symbol`: The name of the unitary operator component in the trajectory, stored
  as an isomorphism vector (``\\tilde{U}``). Default is `:Ũ⃗`.
- `control_name::Symbol`: The name of the control signal component to include in the plot,
  typically the time-dependent control parameters ``u(t)`` in ``H(t) = H_0 + \\sum_k u_k(t) H_k``.
  Default is `:u`.
- `kwargs...`: Additional keyword arguments passed to [`NamedTrajectories.plot`](https://docs.harmoniqs.co/NamedTrajectories/dev/generated/plotting/),
  such as `xlims`, `ylims`, or Makie-specific plotting options.

# Returns
- A Makie `Figure` object containing the population plots.

# Example
```julia
using NamedTrajectories
using Piccolo
using Piccolo

# Define Hamiltonian: H = X + a₁(t)Z + a₂(t)Y
H_drift = PAULIS[:X]
H_drives = [PAULIS[:Z], PAULIS[:Y]]

# Generate control trajectory
N = 100
Δt = 0.1
ts = collect(0:Δt:Δt*(N-1))
u = 0.1 * randn(2, length(ts))

# Generate unitaries
Us = exp.(-im * [(H_drift + sum(u[:, k] .* H_drives)) * ts[k] for k = 1:N])

# Create trajectory
traj = NamedTrajectory(
    (
        Ũ⃗ = hcat(operator_to_iso_vec.(Us)...),
        u = u,
        Δt = ts,
    );
    controls = :u,
    timestep = :Δt,
)

# Plot populations for first two columns
plot_unitary_populations(traj)

# Plot only the first column
plot_unitary_populations(traj; unitary_columns=[1])
```

See also: [`NamedTrajectories.plot`](https://docs.harmoniqs.co/NamedTrajectories/dev/generated/plotting/)
"""
function plot_unitary_populations(
    traj::NamedTrajectory;
    unitary_columns::AbstractVector{Int} = 1:2,
    unitary_name::Symbol = :Ũ⃗,
    control_name::Symbol = :u,
    kwargs...,
)

    transformations =
        [(unitary_name => x -> abs2.(iso_vec_to_operator(x)[:, i])) for i ∈ unitary_columns]

    transformation_labels = [L"P" for i ∈ unitary_columns]

    transformation_titles =
        [L"Populations: $\left| U_{:, %$(i)}(t) \right|^2$" for i ∈ unitary_columns]

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

@testitem "Plot unitary populations with default parameters" begin
    using CairoMakie
    using LinearAlgebra

    # Setup: Create a simple 2-qubit system with random controls
    H_drift = PAULIS[:X]
    H_drives = [PAULIS[:Z], PAULIS[:Y]]

    N = 50
    Δt = 0.1
    ts = collect(0:Δt:(Δt*(N-1)))
    u = 0.1 * randn(length(H_drives), length(ts))

    # Generate unitaries via time evolution
    Us = exp.(-im * [(H_drift + sum(u[:, k] .* H_drives)) * ts[k] for k = 1:N])

    traj = NamedTrajectory(
        (Ũ⃗ = hcat(operator_to_iso_vec.(Us)...), u = u, Δt = ts);
        controls = :u,
        timestep = :Δt,
    )

    # Test: Default behavior plots first two columns
    fig = plot_unitary_populations(traj)
    save("../../../assets/unitary_populations.png", fig)

    @test fig isa Figure
    @test length(fig.content) > 0  # Figure has content
end

@testitem "Plot unitary populations with custom columns" begin
    using CairoMakie
    using LinearAlgebra

    # Setup: Create identity evolution (trivial case)
    N = 20
    Δt = 0.1
    U = Matrix{ComplexF64}(I, 2, 2)
    Us = [U for _ = 1:N]

    traj = NamedTrajectory(
        (Ũ⃗ = hcat(operator_to_iso_vec.(Us)...), u = zeros(1, N), Δt = fill(Δt, N));
        controls = :u,
        timestep = :Δt,
    )

    # Test: Plot only first column
    fig1 = plot_unitary_populations(traj; unitary_columns = [1])
    @test fig1 isa Figure

    # Test: Plot both columns explicitly
    fig2 = plot_unitary_populations(traj; unitary_columns = [1, 2])
    @test fig2 isa Figure
end

@testitem "Plot unitary populations with custom names" begin
    using CairoMakie
    using LinearAlgebra

    # Setup: Custom component names
    N = 15
    Δt = 0.05
    U = PAULIS[:X]  # Simple X gate
    Us = [U for _ = 1:N]

    traj = NamedTrajectory(
        (
            U_iso = hcat(operator_to_iso_vec.(Us)...),
            control = zeros(1, N),
            Δt = fill(Δt, N),
        );
        controls = :control,
        timestep = :Δt,
    )

    # Test: Use custom unitary and control names
    fig = plot_unitary_populations(
        traj;
        unitary_name = :U_iso,
        control_name = :control,
        unitary_columns = [1, 2],
    )

    @test fig isa Figure
end
