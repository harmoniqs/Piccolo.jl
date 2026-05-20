module PiccoloQuantumToolboxExt

using Piccolo
using QuantumToolbox
using Makie
using NamedTrajectories
using LinearAlgebra

using TestItems

# `animate_bloch`, `animate_wigner`, `plot_bloch!`, and `plot_wigner!` extend Piccolo
# stubs. Docstrings for those live in `src/visualizations/quantum_toolbox.jl`. The
# `plot_bloch` and `plot_wigner` docstrings remain here, attached to the QuantumToolbox-
# owned methods we extend.

function iso_to_bloch(ψ̃::AbstractVector{<:Real}, subspace::AbstractVector{Int})
    ψ = iso_to_ket(ψ̃)[subspace]
    return density_to_bloch(ψ * ψ')
end

function iso_vec_to_bloch(ρ̃⃗::AbstractVector{<:Real}, subspace::AbstractVector{Int})
    ρ = iso_vec_to_density(ρ̃⃗)[subspace, subspace]
    return density_to_bloch(ρ)
end

function density_to_bloch(ρ::AbstractMatrix{<:Number})
    x = real(ρ[1, 2] + ρ[2, 1])
    y = imag(ρ[2, 1] - ρ[1, 2])
    z = real(ρ[1, 1] - ρ[2, 2])
    return Point3f([x, y, z])
end

# Reduce size for correct arrow head position
bloch_arrow(v, arrow_size) = (1 - arrow_size / norm(v)) * Vec3f(v...)


"""
    plot_bloch(
        traj::NamedTrajectory;
        index::Union{Nothing, Int} = nothing,
        state_name::Symbol = :ψ̃,
        state_type::Symbol = :ket,
        subspace::AbstractVector{Int} = 1:2,
        kwargs...
    ) -> Figure

Plot the trajectory of a quantum state on the Bloch sphere.

Visualizes how a two-level quantum state evolves on the Bloch sphere representation.
The trajectory path is shown as a line connecting Bloch vectors at each timestep.

# Mathematical Background

For a two-level quantum state ``|\\psi\\rangle = \\alpha|0\\rangle + \\beta|1\\rangle``, the Bloch vector is:

```math
\\vec{r} = (x, y, z) = (\\langle \\sigma_x \\rangle, \\langle \\sigma_y \\rangle, \\langle \\sigma_z \\rangle)
```

where the components are:
```math
\\begin{aligned}
x &= \\langle \\sigma_x \\rangle = 2\\text{Re}(\\alpha^*\\beta) \\\\
y &= \\langle \\sigma_y \\rangle = 2\\text{Im}(\\alpha^*\\beta) \\\\
z &= \\langle \\sigma_z \\rangle = |\\alpha|^2 - |\\beta|^2
\\end{aligned}
```

Pure states lie on the surface of the Bloch sphere (``||\\vec{r}|| = 1``), while mixed states
lie in the interior.

# Arguments
- `traj::NamedTrajectory`: The trajectory containing quantum states to plot.

# Keyword Arguments
- `index::Union{Nothing, Int}`: If provided, add a vector arrow at this time index.
  If `nothing`, no arrow is shown. Default is `nothing`.
- `state_name::Symbol`: The name of the quantum state component in the trajectory. Default is `:ψ̃`.
- `state_type::Symbol`:
  - `:ket` - State vector ``|\\psi\\rangle`` in isomorphism form
  - `:density` - Density matrix ``\\rho`` in isomorphism form
  Default is `:ket`.
- `subspace::AbstractVector{Int}`: The qubit subspace indices to extract (for multi-level systems).
  Default is `1:2` (first two levels).
- `kwargs...`: Additional keyword arguments passed to `QuantumToolbox.render`.

# Returns
- A Makie `Figure` object containing the Bloch sphere visualization.

# Examples

```julia
using Piccolo
using QuantumToolbox
using CairoMakie

# After solving a ket trajectory problem:
traj_ket = get_trajectory(qcp_ket)

# Plot trajectory
fig = plot_bloch(traj_ket)

# Show vector at timestep 25
fig = plot_bloch(traj_ket; index=25)

# For a multi-level system, restrict to the qubit subspace
fig = plot_bloch(traj_ket; subspace=1:2)

# For a density-matrix trajectory
fig = plot_bloch(traj_ρ; state_name=:ρ̃⃗, state_type=:density)
```

See also: `animate_bloch`, `plot_wigner`.
"""
function QuantumToolbox.plot_bloch(
    traj::NamedTrajectory;
    index::Union{Nothing,Int} = nothing,
    state_name::Symbol = :ψ̃,
    state_type::Symbol = :ket,
    subspace::AbstractVector{Int} = 1:2,
    kwargs...,
)
    @assert state_name in traj.names "$state_name ∉ traj.names"
    bloch_pts = map(eachcol(traj[state_name])) do s
        if state_type == :ket
            iso_to_bloch(s, subspace)
        elseif state_type == :density
            iso_vec_to_bloch(s, subspace)
        else
            raise(ArgumentError("State type must be :ket or :density."))
        end
    end

    # Render Bloch sphere
    b = QuantumToolbox.Bloch()
    QuantumToolbox.add_points!(b, stack(bloch_pts))
    fig, lscene = QuantumToolbox.render(b; kwargs...)

    # Add line connecting points
    lines!(lscene, bloch_pts, color = :black)

    if !isnothing(index)
        @assert 1 ≤ index ≤ length(bloch_pts) "Invalid vector index."

        # Save for animation
        fig.attributes[:bloch] = b
        fig.attributes[:state_name] = state_name
        fig.attributes[:subspace] = subspace
        fig.attributes[:vec] = [bloch_arrow(bloch_pts[index], b.vector_tiplength)]

        # Draw the saved vec observable
        arrows3d!(
            lscene,
            [Point3f(0, 0, 0)],
            fig.attributes[:vec],
            shaftradius = b.vector_width,
            tiplength = b.vector_tiplength,
            tipradius = b.vector_tipradius,
            rasterize = 3,
        )
    end

    return fig
end


function Piccolo.plot_bloch!(fig::Figure, traj::NamedTrajectory, idx::Int; kwargs...)
    @assert 1 ≤ idx ≤ traj.N "Invalid knot point index."

    if haskey(fig.attributes, :vec)
        b = fig.attributes[:bloch][]
        state_name = fig.attributes[:state_name][]
        subspace = fig.attributes[:subspace][]
        v = iso_to_bloch(traj[idx][state_name], subspace)
        fig.attributes[:vec] = [bloch_arrow(v, b.b.vector_tiplength)]
    end

    return fig
end


function Piccolo.animate_bloch(
    traj::NamedTrajectory;
    fps::Int = 24,
    mode::Symbol = :inline,
    filename = "bloch_animation.mp4",
    kwargs...,
)
    fig = QuantumToolbox.plot_bloch(traj; index = 1, kwargs...)

    return Piccolo.animate_figure(
        fig,
        1:traj.N,
        i -> plot_bloch!(fig, traj, i; kwargs...),
        mode = mode,
        fps = fps,
        filename = filename,
    )
end


"""
    plot_wigner(
        traj::NamedTrajectory,
        idx::Int;
        state_name::Symbol = :ψ̃,
        state_type::Symbol = :ket,
        xvec = -5:0.1:5,
        yvec = -5:0.1:5,
        projection::Val = Val(:two_dim),
        colorbar::Bool = true,
        kwargs...
    ) -> Figure

Plot the Wigner function of a quantum state at a specific time index.

The Wigner function is a quasi-probability distribution that represents quantum states
in phase space.

# Mathematical Background

For a quantum state ``\\rho``, the Wigner function is

```math
W(x, p) = \\frac{1}{\\pi\\hbar} \\int_{-\\infty}^{\\infty} \\langle x - y | \\rho | x + y \\rangle e^{2ipy/\\hbar} dy
```

- **Classical states** (e.g., coherent): non-negative everywhere, Gaussian peak.
- **Non-classical states** (e.g., Fock, cat): exhibit negative regions.

# Arguments
- `traj::NamedTrajectory`: The trajectory containing quantum states.
- `idx::Int`: Time index of the state to plot (`1 ≤ idx ≤ traj.N`).

# Keyword Arguments
- `state_name::Symbol`: Component name for the state. Default `:ψ̃`.
- `state_type::Symbol`: `:ket` or `:density`. Default `:ket`.
- `xvec`: Range for x-axis (position / real part). Default `-5:0.1:5`.
- `yvec`: Range for y-axis (momentum / imaginary part). Default `-5:0.1:5`.
- `projection::Val`: Projection. Default `Val(:two_dim)`.
- `colorbar::Bool`: Show colorbar. Default `true`.
- `kwargs...`: Forwarded to `QuantumToolbox.plot_wigner`.

# Returns
- A Makie `Figure` containing the Wigner function visualization.

# Examples

```julia
using Piccolo
using QuantumToolbox
using CairoMakie

# Plot Wigner function at final time
fig = plot_wigner(traj_cavity, traj_cavity.N)

# Plot at a specific timestep with a custom grid
fig = plot_wigner(traj_cavity, 50; xvec=-4:0.05:4, yvec=-4:0.05:4)

# For a density-matrix trajectory
fig = plot_wigner(traj_ρ, 1; state_name=:ρ̃⃗, state_type=:density)
```

See also: `animate_wigner`, `plot_bloch`.
"""
function QuantumToolbox.plot_wigner(
    traj::NamedTrajectory,
    idx::Int;
    state_name::Symbol = :ψ̃,
    state_type::Symbol = :ket,
    kwargs...,
)
    @assert 1 ≤ idx ≤ traj.N "Invalid knot point index."

    state = if state_type == :ket
        QuantumObject(iso_to_ket(traj[idx][state_name]))
    elseif state_type == :density
        QuantumObject(iso_vec_to_density(traj[idx][state_name]))
    else
        raise(ArgumentError("State type must be :ket or :density."))
    end

    fig, ax, hm = QuantumToolbox.plot_wigner(state; library = Val(:Makie), kwargs...)
    colsize!(fig.layout, 1, Aspect(1, 1.0))

    # Save attributes for animations
    fig.attributes[:state_name] = state_name
    fig.attributes[:state_type] = state_type
    fig.attributes[:ax] = ax
    fig.attributes[:hm] = hm
    fig.attributes[:label] = Label(fig[0, 1], "Timestep $idx", tellwidth = false)

    return fig
end

function Piccolo.plot_wigner!(fig::Figure, traj::NamedTrajectory, idx::Int)
    @assert 1 ≤ idx ≤ traj.N "Invalid knot point index."

    # Extract attributes from the figure
    state_name = fig.attributes[:state_name][]
    state_type = fig.attributes[:state_type][]
    hm = fig.attributes[:hm][]
    label = fig.attributes[:label][]

    # Update heatmap with new Wigner function
    state = if state_type == :ket
        QuantumObject(iso_to_ket(traj[idx][state_name]))
    elseif state_type == :density
        QuantumObject(iso_vec_to_density(traj[idx][state_name]))
    else
        raise(ArgumentError("State type must be :ket or :density."))
    end
    W = transpose(wigner(state, hm[1][], hm[2][]))
    hm[3][] = W
    label.text[] = "Timestep $idx"
    return fig
end

function Piccolo.animate_wigner(
    traj::NamedTrajectory;
    mode = :inline,
    fps::Int = 24,
    filename = "wigner_animation.mp4",
    kwargs...,
)
    # Setup: plot first frame and capture observables
    fig = QuantumToolbox.plot_wigner(traj, 1; kwargs...)

    Piccolo.animate_figure(
        fig,
        1:traj.N,
        i -> plot_wigner!(fig, traj, i),
        mode = mode,
        fps = fps,
        filename = filename,
    )

    return fig
end


# ============================================================================ #


@testitem "Test plot_bloch for Bloch sphere ket trajectory" begin
    using QuantumToolbox
    using NamedTrajectories
    using Piccolo
    using CairoMakie

    x = ComplexF64[1.0; 0.0]
    y = ComplexF64[0.0, 1.0]

    comps = (ψ̃ = hcat(ket_to_iso(x), ket_to_iso(y)), Δt = [1.0; 1.0])
    traj = NamedTrajectory(comps)

    fig = QuantumToolbox.plot_bloch(traj)
    @test fig isa Figure
end

@testitem "Test plot_bloch for Bloch sphere density trajectory" begin
    using QuantumToolbox
    using NamedTrajectories
    using Piccolo
    using CairoMakie

    x = ComplexF64[1.0; 0.0]
    y = ComplexF64[0.0, 1.0]
    ρ̃⃗ = hcat(density_to_iso_vec(x * x'), density_to_iso_vec(y * y'))
    traj = NamedTrajectory((ρ̃⃗ = ρ̃⃗, Δt = [1.0; 1.0]))

    fig = QuantumToolbox.plot_bloch(traj, state_name = :ρ̃⃗, state_type = :density)
    @test fig isa Figure
end

@testitem "Test plot_bloch for Bloch sphere trajectory with one vector arrow shown" begin
    using QuantumToolbox
    using NamedTrajectories
    using Piccolo
    using CairoMakie

    x = ComplexF64[1.0; 0.0]
    y = ComplexF64[0.0, 1.0]
    ψ̃ = hcat(ket_to_iso(x), ket_to_iso(y))
    traj = NamedTrajectory((ψ̃ = ψ̃, Δt = [1.0; 1.0]))

    fig = QuantumToolbox.plot_bloch(traj, index = 1)
    @test fig isa Figure
end

@testitem "plot_bloch shows expected curved Bloch path" begin
    using QuantumToolbox
    using NamedTrajectories
    using Piccolo
    using CairoMakie

    T = 20
    ts = range(0, π/2; length = T)

    kets = [
        QuantumObject(cos(θ) * [1.0 + 0im, 0.0 + 0im] + sin(θ) * [0.0 + 0im, 1.0 + 0im]) for θ in ts
    ]

    iso_kets = ket_to_iso.(ψ.data for ψ in kets)

    ψ̃ = hcat(iso_kets...)
    Δt = fill(1.0, T)

    comps = (ψ̃ = ψ̃, Δt = Δt)
    traj = NamedTrajectory(comps)
    fig = QuantumToolbox.plot_bloch(traj, state_name = :ψ̃)
    @test fig isa Figure
end

@testitem "Plot Wigner function of coherent state" begin
    using QuantumToolbox
    using NamedTrajectories
    using Piccolo
    using CairoMakie

    N = 20
    α = 1.5 + 0.5im
    ψ = coherent(N, α)

    traj = NamedTrajectory((ψ̃ = hcat(ket_to_iso(ψ.data)), Δt = [1.0]))
    fig = QuantumToolbox.plot_wigner(traj, 1, state_name = :ψ̃)

    @test fig isa Figure
end

@testitem "Plot Wigner function of density" begin
    using QuantumToolbox
    using NamedTrajectories
    using Piccolo
    using CairoMakie

    N = 20
    α = 1.5 + 0.5im
    ψ = coherent(N, α)
    ρ = ψ * ψ'

    traj = NamedTrajectory((ρ̃⃗ = hcat(density_to_iso_vec(ρ.data)), Δt = [1.0]))
    fig = QuantumToolbox.plot_wigner(traj, 1, state_name = :ρ̃⃗, state_type = :density)

    @test fig isa Figure
end

end
