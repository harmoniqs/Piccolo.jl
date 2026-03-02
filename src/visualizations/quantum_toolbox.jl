module QuantumToolboxPlots

export animate_bloch
export animate_wigner
export plot_bloch!
export plot_wigner!

using TestItems

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
- `index::Union{Nothing, Int}`: If provided, add a vector arrow at this time index. Useful for
  animations. If `nothing`, no arrow is shown. Default is `nothing`.
- `state_name::Symbol`: The name of the quantum state component in the trajectory. Default is `:ψ̃`.
- `state_type::Symbol`: The type of quantum state representation:
  - `:ket` - State vector ``|\\psi\\rangle`` in isomorphism form
  - `:density` - Density matrix ``\\rho`` in isomorphism form
  Default is `:ket`.
- `subspace::AbstractVector{Int}`: The qubit subspace indices to extract (for multi-level systems).
  Default is `1:2` (first two levels).
- `kwargs...`: Additional keyword arguments passed to `QuantumToolbox.render` (e.g., `xlims`, `ylims`).

# Returns
- A Makie `Figure` object containing the Bloch sphere visualization.

# Examples

## Plot qubit trajectory on Bloch sphere
```julia
using Piccolo
using QuantumToolbox
using NamedTrajectories
using CairoMakie

# Create a trajectory of states evolving on Bloch sphere
T = 50
θs = range(0, π, length=T)
ϕs = range(0, 2π, length=T)

# States: |ψ(t)⟩ = cos(θ/2)|0⟩ + e^(iϕ)sin(θ/2)|1⟩
kets = [cos(θ/2) * [1.0+0im, 0.0] + exp(im*ϕ) * sin(θ/2) * [0.0, 1.0+0im] 
        for (θ, ϕ) in zip(θs, ϕs)]

traj = NamedTrajectory(
    (ψ̃ = hcat(ket_to_iso.(kets)...), Δt = fill(0.1, T)),
    timestep = :Δt
)

# Plot trajectory
fig = plot_bloch(traj)
```

## Plot with state vector arrow
```julia
# Show vector at timestep 25
fig = plot_bloch(traj; index=25)
```

## Plot multi-level system (extract qubit subspace)
```julia
# For a qutrit (3-level) system, plot only computational qubit subspace
T = 30
kets = [normalize([1.0+0im, 0.5, 0.2]) for _ in 1:T]  # 3-level states

traj = NamedTrajectory(
    (ψ̃ = hcat(ket_to_iso.(kets)...), Δt = fill(0.1, T)),
    timestep = :Δt
)

# Extract and plot first two levels
fig = plot_bloch(traj; subspace=1:2)
```

## Plot density matrix trajectory
```julia
# Create mixed states
kets = [normalize([cos(t), sin(t)]) for t in range(0, π/2, length=20)]
ρs = [ket * ket' for ket in kets]

traj = NamedTrajectory(
    (ρ̃⃗ = hcat(density_to_iso_vec.(ρs)...), Δt = fill(0.1, 20)),
    timestep = :Δt
)

fig = plot_bloch(traj; state_name=:ρ̃⃗, state_type=:density)
```

See also: [`animate_bloch`](@ref), [`plot_wigner`](@ref)
"""
function animate_bloch end

function plot_bloch! end

"""
    animate_bloch(
        traj::NamedTrajectory;
        fps::Int = 24,
        mode::Symbol = :inline,
        filename::String = "bloch_animation.mp4",
        state_name::Symbol = :ψ̃,
        state_type::Symbol = :ket,
        subspace::AbstractVector{Int} = 1:2,
        kwargs...
    ) -> Figure

Animate the evolution of a quantum state on the Bloch sphere.

Creates an animation showing a quantum state's trajectory on the Bloch sphere,
with a moving vector arrow tracking the current state.

# Arguments
- `traj::NamedTrajectory`: The trajectory containing quantum states to animate.

# Keyword Arguments
- `fps::Int`: Frames per second for animation playback. Default is `24`.
- `mode::Symbol`: Animation mode:
  - `:inline` - Display in interactive window (requires `GLMakie`).
  - `:record` - Save to file.
  Default is `:inline`.
- `filename::String`: Output filename when `mode = :record`. Default is `"bloch_animation.mp4"`.
- `state_name::Symbol`: The name of the quantum state component. Default is `:ψ̃`.
- `state_type::Symbol`: State representation (`:ket` or `:density`). Default is `:ket`.
- `subspace::AbstractVector{Int}`: Subspace indices for multi-level systems. Default is `1:2`.
- `kwargs...`: Additional arguments passed to `plot_bloch` and `QuantumToolbox.render`.

# Returns
- A Makie `Figure` object containing the animation.

# Examples

## Animate Rabi oscillation
```julia
using Piccolo
using QuantumToolbox
using GLMakie

# Rabi oscillation: |ψ(t)⟩ = cos(Ωt/2)|0⟩ + sin(Ωt/2)|1⟩
T = 100
Ω = 2π  # Rabi frequency
times = range(0, 2π/Ω, length=T)

kets = [cos(Ω*t/2) * [1.0+0im, 0.0] + sin(Ω*t/2) * [0.0, 1.0+0im] for t in times]

traj = NamedTrajectory(
    (ψ̃ = hcat(ket_to_iso.(kets)...), Δt = fill(times[2]-times[1], T)),
    timestep = :Δt
)

# Animate evolution
animate_bloch(traj; fps=30)
```

## Save animation to file
```julia
using CairoMakie

# Same trajectory as above
animate_bloch(traj; mode=:record, filename="rabi.mp4", fps=24)
```

## Animate solved quantum control problem
```julia
using GLMakie

# Assume qcp is a solved QuantumControlProblem with KetTrajectory
traj = get_trajectory(qcp)

# Animate the state evolution
animate_bloch(traj; state_name=:ψ̃, fps=20)
```

## Animate qutrit in qubit subspace
```julia
# For 3-level system, animate only computational subspace
T = 50
kets = [normalize([1.0, 0.5*sin(t), 0.2*cos(t)]) for t in range(0, 2π, length=T)]

traj = NamedTrajectory(
    (ψ̃ = hcat(ket_to_iso.(kets)...), Δt = fill(0.1, T)),
    timestep = :Δt
)

animate_bloch(traj; subspace=1:2, fps=25)
```

See also: [`plot_bloch`](@ref), [`animate_figure`](@ref), [`animate_wigner`](@ref)
"""
function animate_wigner end

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
in phase space. It provides insight into the quantum-classical boundary and can reveal
non-classical features like negativity.

# Mathematical Background

For a quantum state ``\\rho``, the Wigner function is defined as:

```math
W(x, p) = \\frac{1}{\\pi\\hbar} \\int_{-\\infty}^{\\infty} \\langle x - y | \\rho | x + y \\rangle e^{2ipy/\\hbar} dy
```

For dimensionless phase space coordinates ``(\\alpha, \\alpha^*)`` or ``(q, p)``:
- **Classical states**: Wigner function is non-negative everywhere
- **Non-classical states**: Exhibit negative regions (e.g., Fock states, Schrödinger cats)
- **Coherent states**: Appear as Gaussian peaks

# Arguments
- `traj::NamedTrajectory`: The trajectory containing quantum states.
- `idx::Int`: The time index of the state to plot (must satisfy `1 ≤ idx ≤ traj.N`).

# Keyword Arguments
- `state_name::Symbol`: The name of the quantum state component in the trajectory. Default is `:ψ̃`.
- `state_type::Symbol`: The state representation:
  - `:ket` - State vector in isomorphism form
  - `:density` - Density matrix in isomorphism form
  Default is `:ket`.
- `xvec`: Range for x-axis (position/real part). Default is `-5:0.1:5`.
- `yvec`: Range for y-axis (momentum/imaginary part). Default is `-5:0.1:5`.
- `projection::Val`: Projection type. Default is `Val(:two_dim)` for 2D heatmap.
- `colorbar::Bool`: Whether to display colorbar. Default is `true`.
- `kwargs...`: Additional keyword arguments passed to `QuantumToolbox.plot_wigner`.

# Returns
- A Makie `Figure` object containing the Wigner function visualization.

# Examples

## Plot Wigner function of coherent state
```julia
using Piccolo
using QuantumToolbox
using CairoMakie

# Create coherent state |α⟩ with α = 2 + i
N = 20  # Fock space cutoff
α = 2.0 + 1.0im
ψ = coherent(N, α)

traj = NamedTrajectory(
    (ψ̃ = hcat(ket_to_iso(ψ.data)), Δt = [1.0]),
    timestep = :Δt
)

# Plot Wigner function - should show Gaussian peak at α
fig = plot_wigner(traj, 1; state_name=:ψ̃)
```

## Plot Fock state (non-classical)
```julia
# Fock state |n=3⟩ - exhibits negative regions
N = 20
ψ = fock(N, 3)

traj = NamedTrajectory(
    (ψ̃ = hcat(ket_to_iso(ψ.data)), Δt = [1.0]),
    timestep = :Δt
)

# Adjust grid for better resolution
fig = plot_wigner(traj, 1; xvec=-4:0.05:4, yvec=-4:0.05:4)
```

## Plot Wigner of density matrix
```julia
# Thermal state
N = 20
n̄ = 2.0  # Average photon number
ρ = thermal_dm(N, n̄)

traj = NamedTrajectory(
    (ρ̃⃗ = hcat(density_to_iso_vec(ρ.data)), Δt = [1.0]),
    timestep = :Δt
)

fig = plot_wigner(traj, 1; state_name=:ρ̃⃗, state_type=:density)
```

## Compare Wigner functions at different times
```julia
# Evolving state
T = 50
kets = [coherent(20, 2*exp(im*2π*t/T)).data for t in 1:T]

traj = NamedTrajectory(
    (ψ̃ = hcat(ket_to_iso.(kets)...), Δt = fill(0.1, T)),
    timestep = :Δt
)

# Initial state
fig1 = plot_wigner(traj, 1)

# Final state
fig2 = plot_wigner(traj, T)
```

See also: [`animate_wigner`](@ref), [`plot_bloch`](@ref)
"""
function plot_wigner! end

"""
    animate_wigner(
        traj::NamedTrajectory;
        mode::Symbol = :inline,
        fps::Int = 24,
        filename::String = "wigner_animation.mp4",
        state_name::Symbol = :ψ̃,
        state_type::Symbol = :ket,
        kwargs...
    ) -> Figure

Animate the time evolution of the Wigner function for a quantum state trajectory.

Creates an animation showing how the Wigner quasi-probability distribution evolves
in phase space over time.

# Arguments
- `traj::NamedTrajectory`: The trajectory containing quantum states to animate.

# Keyword Arguments
- `mode::Symbol`: Animation mode:
  - `:inline` - Display in interactive window (requires `GLMakie`).
  - `:record` - Save to file.
  Default is `:inline`.
- `fps::Int`: Frames per second for animation playback. Default is `24`.
- `filename::String`: Output filename when `mode = :record`. Default is `"wigner_animation.mp4"`.
- `state_name::Symbol`: The name of the quantum state component. Default is `:ψ̃`.
- `state_type::Symbol`: State representation (`:ket` or `:density`). Default is `:ket`.
- `kwargs...`: Additional arguments passed to `plot_wigner` (e.g., `xvec`, `yvec`, `colorbar`).

# Returns
- A Makie `Figure` object containing the animation.

# Examples

## Animate coherent state rotation
```julia
using Piccolo
using QuantumToolbox
using GLMakie

# Coherent state rotating in phase space: |α(t)⟩ with α(t) = 2e^(iωt)
T = 100
ω = 2π
times = range(0, 2π/ω, length=T)
N = 20

kets = [coherent(N, 2*exp(im*ω*t)).data for t in times]

traj = NamedTrajectory(
    (ψ̃ = hcat(ket_to_iso.(kets)...), Δt = fill(times[2]-times[1], T)),
    timestep = :Δt
)

# Animate Wigner function - Gaussian peak circles origin
animate_wigner(traj; fps=30, xvec=-3:0.1:3, yvec=-3:0.1:3)
```

## Save animation to file
```julia
using CairoMakie

# Same trajectory as above
animate_wigner(traj; mode=:record, filename="coherent_rotation.mp4", fps=24)
```

## Animate Fock state evolution
```julia
# Superposition evolving: |ψ(t)⟩ = cos(t)|0⟩ + sin(t)|1⟩
T = 60
times = range(0, π/2, length=T)
N = 20

kets = [normalize(cos(t) * fock(N, 0).data + sin(t) * fock(N, 1).data) for t in times]

traj = NamedTrajectory(
    (ψ̃ = hcat(ket_to_iso.(kets)...), Δt = fill(times[2]-times[1], T)),
    timestep = :Δt
)

# Animate - watch negative regions appear
animate_wigner(traj; fps=20, xvec=-4:0.1:4, yvec=-4:0.1:4)
```

## Animate quantum control result
```julia
using GLMakie

# Assume qcp is a solved QuantumControlProblem for cavity/oscillator
traj = get_trajectory(qcp)

# Animate Wigner function evolution
animate_wigner(traj; state_name=:ψ̃, fps=15)
```

## Animate with custom grid
```julia
# High-resolution animation for detailed features
animate_wigner(
    traj;
    xvec = -6:0.05:6,
    yvec = -6:0.05:6,
    fps = 30,
    colorbar = true
)
```

See also: [`plot_wigner`](@ref), [`animate_bloch`](@ref), [`animate_figure`](@ref)
"""


end
