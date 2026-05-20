# Piccolo.jl Visualization API - AI Agent Context

This document provides a comprehensive guide to Piccolo.jl's visualization capabilities, specifically structured for AI coding assistants. It includes common patterns, function signatures, and complete working examples.

## Overview

Piccolo.jl provides visualization tools for quantum optimal control problems, built on top of the Makie.jl ecosystem. The visualization API has three main components:

1. **Quantum Object Plots** - Plot populations for unitaries and states
2. **Animation Functions** - Animate trajectory evolution
3. **Quantum Toolbox Plots** - Bloch sphere and Wigner function visualizations

## Backend Selection

```julia
# For static plots (publication-quality, fast)
using CairoMakie

# For interactive plots and :inline animations (recommended for development)
using GLMakie

# For web-based interactive plots
using WGLMakie
```

**Important**: CairoMakie cannot handle `:inline` animations. Use `:record` mode or switch to GLMakie.

## Core Visualization Functions

### 1. plot_unitary_populations

**Purpose**: Visualize how unitary operator populations evolve during gate synthesis.

**Function Signature**:
```julia
plot_unitary_populations(
    traj::NamedTrajectory;
    unitary_columns::AbstractVector{Int} = 1:2,
    unitary_name::Symbol = :Ũ⃗,
    control_name::Symbol = :u,
    kwargs...
) -> Figure
```

**What it plots**: For each column `j` of the unitary matrix `U(t)`, plots `|U_{i,j}(t)|²` for all rows `i`. This shows the probability distribution of basis states when the system starts in state `|j⟩`.

**Common Usage Patterns**:

```julia
# After solving a gate synthesis problem
using Piccolo
using CairoMakie

qcp = # ... your solved quantum control problem
traj = get_trajectory(qcp)

# Plot populations for first two columns (default)
fig = plot_unitary_populations(traj)

# Plot only first column
fig = plot_unitary_populations(traj; unitary_columns=[1])

# Plot with custom unitary name
fig = plot_unitary_populations(traj; unitary_name=:U_iso, unitary_columns=[1,2,3])
```

**Key Points for AI Agents**:
- Always use `get_trajectory(qcp)` to extract trajectory from a problem
- `unitary_name` must match the actual component name in the trajectory
- Returns a Makie Figure that can be saved with `save("filename.png", fig)`
- Automatically includes control signals in the plot

### 2. plot_state_populations

**Purpose**: Visualize populations for multiple quantum states (e.g., from sampling problems).

**Function Signature**:
```julia
plot_state_populations(
    traj::NamedTrajectory;
    state_name::Symbol = :ψ̃,
    state_indices::Union{Nothing, AbstractVector{Int}} = nothing,
    control_name::Symbol = :u,
    subspace::Union{Nothing, AbstractVector{Int}} = nothing,
    kwargs...
) -> Figure
```

**What it plots**: For state `|ψ⟩`, plots `|ψᵢ|²` for each basis component `i`.

**Common Usage Patterns**:

```julia
# For a ket trajectory problem
traj = get_trajectory(qcp_ket)

# Plot all states matching :ψ̃ prefix
fig = plot_state_populations(traj)

# Plot only computational subspace (exclude leakage)
fig = plot_state_populations(traj; subspace=1:2)

# Plot specific states by index
fig = plot_state_populations(traj; state_indices=[1, 3, 5])
```

**Key Points for AI Agents**:
- Automatically finds states matching the prefix (e.g., `:ψ̃1_system_1`, `:ψ̃2_system_1`)
- Use `subspace` to filter out leakage levels in multi-level systems
- `state_indices` allows selecting specific states from the full set

### 3. plot_bloch

**Purpose**: Visualize qubit state trajectory on the Bloch sphere.

**Function Signature**:
```julia
plot_bloch(
    traj::NamedTrajectory;
    index::Union{Nothing, Int} = nothing,
    state_name::Symbol = :ψ̃,
    state_type::Symbol = :ket,
    subspace::AbstractVector{Int} = 1:2,
    kwargs...
) -> Figure
```

**What it plots**: Trajectory of quantum state on Bloch sphere, with optional vector arrow at specific time.

**Common Usage Patterns**:

```julia
using QuantumToolbox
using CairoMakie

# Plot trajectory path
fig = plot_bloch(traj)

# Add vector arrow at timestep 50
fig = plot_bloch(traj; index=50)

# For multi-level system, extract qubit subspace
fig = plot_bloch(traj; subspace=1:2)

# Plot density matrix trajectory
fig = plot_bloch(traj; state_name=:ρ̃⃗, state_type=:density)
```

**Requirements**:
- Requires `QuantumToolbox.jl` extension
- State must be 2-level (or extract 2-level subspace)
- Works with both ket and density matrix representations

### 4. animate_bloch

**Purpose**: Animate Bloch sphere trajectory with moving vector.

**Function Signature**:
```julia
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
```

**Common Usage Patterns**:

```julia
using GLMakie  # Required for :inline mode

# Interactive animation
fig = animate_bloch(traj; fps=30)

# Save to file
using CairoMakie
fig = animate_bloch(traj; mode=:record, filename="evolution.mp4", fps=24)
```

### 5. plot_wigner

**Purpose**: Plot Wigner quasi-probability distribution in phase space.

**Function Signature**:
```julia
plot_wigner(
    traj::NamedTrajectory,
    idx::Int;
    state_name::Symbol = :ψ̃,
    state_type::Symbol = :ket,
    xvec = -5:0.1:5,
    yvec = -5:0.1:5,
    kwargs...
) -> Figure
```

**Common Usage Patterns**:

```julia
# Plot Wigner function at final time
fig = plot_wigner(traj, traj.N)

# Plot at specific timestep with custom grid
fig = plot_wigner(traj, 50; xvec=-4:0.05:4, yvec=-4:0.05:4)

# For density matrix
fig = plot_wigner(traj, 1; state_name=:ρ̃⃗, state_type=:density)
```

### 6. animate_wigner

**Purpose**: Animate Wigner function evolution.

**Function Signature**:
```julia
animate_wigner(
    traj::NamedTrajectory;
    mode::Symbol = :inline,
    fps::Int = 24,
    filename::String = "wigner_animation.mp4",
    kwargs...
) -> Figure
```

### 7. animate_figure (Low-level)

**Purpose**: Generic animation function for custom visualizations.

**Function Signature**:
```julia
animate_figure(
    fig::Figure,
    frames::AbstractVector{Int},
    update_frame!::Function;
    mode::Symbol = :inline,
    fps::Int = 24,
    filename::String = "animation.mp4"
) -> Figure
```

**Common Usage Patterns**:

```julia
using GLMakie

# Custom animation
fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, Float64[], Float64[])

function update_frame!(k)
    empty!(ax)
    lines!(ax, times[1:k], data[1:k])
end

animate_figure(fig, 1:100, update_frame!; fps=30)
```

### 8. animate_name

**Purpose**: Animate trajectory component appearing progressively.

**Function Signature**:
```julia
animate_name(
    traj::NamedTrajectory,
    name::Symbol;
    fps::Int = 24,
    mode::Symbol = :inline,
    filename::String = "name_animation.mp4",
    kwargs...
) -> Figure
```

**Common Usage Patterns**:

```julia
# Animate controls appearing
animate_name(traj, :u; fps=30)

# Save to file
animate_name(traj, :u; mode=:record, filename="controls.mp4")
```

## Complete Workflow Examples

### Example 1: Basic Gate Synthesis Visualization

```julia
using Piccolo
using CairoMakie

# Define system
H_drift = 0.5 * PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
sys = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

# Setup problem
T = 10.0
N = 100
pulse = ZeroOrderPulse(0.1 * randn(2, N), T, N)
qtraj = UnitaryTrajectory(sys, pulse, GATES[:X])

# Solve
qcp = SmoothPulseProblem(qtraj, N; Q=100.0, R=1e-2)
solve!(qcp; max_iter=100)

# Visualize
traj = get_trajectory(qcp)
fig = plot_unitary_populations(traj)
save("populations.png", fig)
```

### Example 2: Ket State Visualization with Bloch Sphere

```julia
using Piccolo
using QuantumToolbox
using GLMakie

# Define system and problem for state transfer
ψ_init = ComplexF64[1, 0]
ψ_goal = ComplexF64[0, 1]

pulse = ZeroOrderPulse(0.1 * randn(2, N), T, N)
qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

qcp = SmoothPulseProblem(qtraj, N; Q=100.0, R=1e-2)
solve!(qcp; max_iter=100)

# Extract trajectory
traj = get_trajectory(qcp)

# Plot populations
fig1 = plot_state_populations(traj)

# Animate Bloch sphere
fig2 = animate_bloch(traj; fps=30)
```

### Example 3: Multi-Level System with Leakage Suppression

```julia
using Piccolo
using CairoMakie

# 3-level system (qutrit)
H_drift = create(3) * destroy(3)
H_drives = [PAULIS[:X] ⊗ I(3), PAULIS[:Y] ⊗ I(3)]
sys = QuantumSystem(H_drift, H_drives)

# ... setup and solve problem ...

traj = get_trajectory(qcp)

# Plot all populations
fig1 = plot_unitary_populations(traj; unitary_columns=[1,2,3])

# Plot only computational subspace
fig2 = plot_state_populations(traj; subspace=1:2)

# Bloch sphere for qubit subspace
fig3 = plot_bloch(traj; subspace=1:2)
```

### Example 4: Saving Multiple Visualizations

```julia
using Piccolo
using CairoMakie

traj = get_trajectory(qcp)

# Static plots
fig1 = plot_unitary_populations(traj)
save("populations.png", fig1)

# If QuantumToolbox available
fig2 = plot_bloch(traj)
save("bloch.png", fig2)

# Animations (switch to record mode for CairoMakie)
fig3 = animate_bloch(traj; mode=:record, filename="evolution.mp4")
```

## Common Pitfalls and Solutions

### Issue: CairoMakie animation error
**Problem**: `:inline` animation doesn't work with CairoMakie
**Solution**: Either use `mode=:record` or switch to `using GLMakie`

### Issue: State name not found
**Problem**: `state_name=:ψ̃` error
**Solution**: Check actual names with `traj.names`. Use correct symbol (e.g., `:ψ̃`, `:Ũ⃗`)

### Issue: Bloch plot error for 3-level system
**Problem**: Can't plot 3-level state on Bloch sphere
**Solution**: Use `subspace=1:2` to extract qubit subspace

### Issue: Wrong populations displayed
**Problem**: `unitary_columns` doesn't match expected output
**Solution**: Verify unitary matrix size matches gate dimension. Use `length(sys.levels)` to check

## Type Reference

Key types used in visualization:
- `NamedTrajectory` - Container for trajectory data
- `Figure` - Makie figure object (return type of all plot functions)
- `Symbol` - Component names (e.g., `:u`, `:ψ̃`, `:Ũ⃗`, `:Δt`)
- Animation modes: `:inline` or `:record`
- State types: `:ket` or `:density`

## Quick Reference Card

```julia
# Import visualization backend
using CairoMakie  # or GLMakie

# Basic plots (no animation)
plot_unitary_populations(traj)
plot_state_populations(traj)
plot_bloch(traj)  # requires QuantumToolbox
plot_wigner(traj, idx)  # requires QuantumToolbox

# Animations (need GLMakie for :inline)
animate_bloch(traj)
animate_wigner(traj)
animate_name(traj, :u)

# Save figures
save("output.png", fig)  # Static plot
# Animations auto-save with mode=:record
```

## Additional Resources

- [Piccolo.jl Documentation](https://docs.harmoniqs.co/Piccolo/dev/)
- [NamedTrajectories.jl Plotting Guide](https://docs.harmoniqs.co/NamedTrajectories/dev/generated/plotting/)
- [Makie.jl Documentation](https://docs.makie.org/stable/)
- [QuantumToolbox.jl](https://qutip.github.io/QuantumToolbox.jl/stable/)
