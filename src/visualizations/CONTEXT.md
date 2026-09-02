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

## AI-agent quick start

Use this section when generating Piccolo visualization code for users. Prefer the high-level helpers below before writing manual Makie code.

```julia
using Piccolo
using CairoMakie          # static figures, docs, CI-friendly scripts
# using GLMakie           # interactive local windows / :inline animations
# using WGLMakie          # notebooks / browser output
```

Load `QuantumToolbox` only for Bloch-sphere and Wigner plotting helpers:

```julia
using QuantumToolbox
```

Common solved-problem workflow:

```julia
qcp = SmoothPulseProblem(qtraj, N; Q=100.0, R=1e-2, ddu_bound=1.0)
solve!(qcp)

traj = get_trajectory(qcp)
fig = plot_pulse(qcp; bounds=true, components=[:du, :ddu])
fig = plot_unitary_populations(traj)
```

For ket-state transfer problems, use state population helpers:

```julia
traj = get_trajectory(qcp)
fig = plot_state_populations(traj)
```

For Bloch and Wigner visualizations:

```julia
using QuantumToolbox
fig = plot_bloch(traj; index=traj.N)
fig = plot_wigner(traj, traj.N; xvec=-3:0.1:3, yvec=-3:0.1:3)
```

When plotting a physical pulse, do **not** use `plot(traj, [:u])` as the default. That plots optimization variables at knot points. Use `plot_pulse(qcp)` to reconstruct and draw the actual pulse waveform for the pulse type.

For manual Makie plots from a `NamedTrajectory`, derive control times from trajectory timesteps:

```julia
traj = get_trajectory(qcp)
plot_times = cumsum([0; get_timesteps(traj)])[1:end-1]
lines!(ax, plot_times, traj[:u][1, :])
```

If a uniform visual time grid is important, construct the problem with equal timesteps:

```julia
opts = PiccoloOptions(timesteps_all_equal=true, display=:silent)
qcp = SmoothPulseProblem(qtraj, N; piccolo_options=opts)
```

Export figures with Makie `save`:

```julia
save("controls.png", fig)  # raster
save("controls.pdf", fig)  # vector
save("controls.svg", fig)  # vector
```

## Core Visualization Functions

### 0. plot_pulse / plot_pulse_IQ / plot_pulse_phases

**Purpose**: Render an `AbstractPulse` or solved `QuantumControlProblem` with type-appropriate visuals — step functions for `ZeroOrderPulse`, line segments for `LinearSplinePulse`, smooth curves with knots for `CubicSplinePulse`, dense samples for analytic / `FunctionPulse` types.

**Function Signatures** (high-level overloads):

```julia
plot_pulse(qcp::QuantumControlProblem; title="", bounds=false,
           components::Vector{Symbol}=Symbol[],
           component_bounds::Bool=false,
           kwargs...) -> Figure
plot_pulse(pulse::AbstractPulse; title="", kwargs...) -> Figure

# For 4-drive pulses interpreted as two IQ pairs (Ω_I, Ω_Q, α_I, α_Q)
plot_pulse_IQ(pulse_or_qcp; title="", kwargs...) -> Figure

# Same 4-drive structure, polar form: magnitude + unwrapped phase
plot_pulse_phases(pulse_or_qcp; title="",
                  phase_threshold::Real=0.01,
                  kwargs...) -> Figure
```

**Common Usage Patterns**:

```julia
# After solving:
fig = plot_pulse(qcp; title="Optimized Pulse")

# With hardware bounds shaded + smoothness derivatives stacked beneath
fig = plot_pulse(qcp; bounds=true, components=[:du, :ddu], component_bounds=true)

# Render a standalone pulse object (e.g. resampled, reconstructed, or analytic)
fig = plot_pulse(my_pulse; title="My Pulse")

# IQ view for 4-drive cat-qubit / oscillator problems
fig = plot_pulse_IQ(qcp; title="IQ (Ω, α)")

# Magnitude + phase polar view (phase masked where |amp| < 1% of peak)
fig = plot_pulse_phases(qcp; title="(|·|, ∠·)")
```

**Key Points for AI Agents**:
- This is the canonical pulse visualization — prefer over manually extracting `traj[:u]` and calling `lines!`.
- The `qcp` overload pulls the optimized pulse out for you; the `pulse::AbstractPulse` form is for when you've reconstructed/resampled.
- See the [Pulses concept page](@ref pulses-concept) for per-pulse-type rendering examples and the `:stacked` vs `:overlay` layouts.

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
drive_bounds = [1.0, 1.0]
sys = QuantumSystem(H_drift, H_drives, drive_bounds)

# Setup problem
T = 10.0
N = 100
times = collect(range(0.0, T, length=N))
pulse = ZeroOrderPulse(0.1 * randn(2, N), times)
qtraj = UnitaryTrajectory(sys, pulse, GATES[:X])

# Solve
qcp = SmoothPulseProblem(qtraj, N; Q=100.0, R=1e-2)
solve!(qcp; options=IpoptOptions(max_iter=100))
# Inside docs/literate guides, prefer `cached_solve!(qcp, "cache_name"; max_iter=100)`
# which transparently caches successful solves.

# Visualize
traj = get_trajectory(qcp)
fig = plot_unitary_populations(traj)
save("populations.png", fig)
```

### Example 2: Ket State Visualization with Bloch Sphere

```julia
using Piccolo
using QuantumToolbox
using CairoMakie  # GLMakie also works (required for :inline animation)

# Define system and problem for state transfer
ψ_init = ComplexF64[1, 0]
ψ_goal = ComplexF64[0, 1]

times = collect(range(0.0, T, length=N))
pulse = ZeroOrderPulse(0.1 * randn(2, N), times)
qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

qcp = SmoothPulseProblem(qtraj, N; Q=100.0, R=1e-2)
solve!(qcp; options=IpoptOptions(max_iter=100))

# Extract trajectory
traj = get_trajectory(qcp)

# Plot populations
fig1 = plot_state_populations(traj)

# Bloch sphere (static)
fig2 = plot_bloch(traj)

# Animate Bloch sphere — :inline needs GLMakie, :record works with CairoMakie
fig3 = animate_bloch(traj; mode=:record, filename="bloch.mp4", fps=30)
```

### Example 3: Multi-Level System with Leakage Suppression

For a 3-level transmon, use the `TransmonSystem` template (which constructs a properly
sized drift Hamiltonian + drive operators) and restrict plots to the computational
subspace via `subspace=1:2`:

```julia
using Piccolo
using CairoMakie

# 3-level transmon: drift = anharmonic ladder, drives = (a + a†), -i(a - a†)
sys = TransmonSystem(ω_q=5.0, α=-0.2, n_levels=3; drive_bound=0.1)

# Setup and solve (target = X gate embedded in the 2-level subspace)
T = 50.0
N = 100
times = collect(range(0.0, T, length=N))
pulse = ZeroOrderPulse(0.01 * randn(sys.n_drives, N), times)
qtraj = UnitaryTrajectory(sys, pulse, EmbeddedOperator(GATES[:X], sys; subspace=1:2))

qcp = SmoothPulseProblem(qtraj, N; Q=100.0, R=1e-2, ddu_bound=1.0)
solve!(qcp; options=IpoptOptions(max_iter=200))

traj = get_trajectory(qcp)

# Plot all 3 populations (includes leakage level)
fig1 = plot_unitary_populations(traj; unitary_columns=[1,2,3])

# Plot only the computational subspace (no leakage row)
fig2 = plot_state_populations(traj; subspace=1:2)

# Bloch sphere — restrict to qubit subspace
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

# Pulses (canonical)
plot_pulse(qcp)                                    # optimized pulse straight from problem
plot_pulse(qcp; bounds=true, components=[:du, :ddu], component_bounds=true)
plot_pulse(my_pulse)                               # any AbstractPulse
plot_pulse_IQ(qcp)                                 # 4-drive IQ pair view
plot_pulse_phases(qcp)                             # 4-drive magnitude + phase view

# Populations
plot_unitary_populations(traj)                     # |U_{i,j}(t)|²
plot_state_populations(traj)                       # |ψᵢ(t)|²
plot_state_populations(traj; subspace=1:2)         # restrict to computational subspace

# QuantumToolbox plots (require `using QuantumToolbox`)
plot_bloch(traj)                                   # qubit trajectory on Bloch sphere
plot_bloch(traj; index=k)                          # add vector arrow at frame k
plot_wigner(traj, idx)                             # Wigner function at timestep idx

# Animations (`:inline` needs GLMakie; `:record` works with CairoMakie)
animate_bloch(traj;   mode=:record, filename="bloch.mp4",   fps=30)
animate_wigner(traj;  mode=:record, filename="wigner.mp4",  fps=24)
animate_name(traj, :u; mode=:record, filename="ctrls.mp4",  fps=30)

# Save static figures
save("output.png", fig)                            # PNG / PDF / SVG by extension
```

## Additional Resources

- [Piccolo.jl Documentation](https://docs.harmoniqs.co/Piccolo/dev/)
- [NamedTrajectories.jl Plotting Guide](https://docs.harmoniqs.co/NamedTrajectories/dev/generated/plotting/)
- [Makie.jl Documentation](https://docs.makie.org/stable/)
- [QuantumToolbox.jl](https://qutip.github.io/QuantumToolbox.jl/stable/)
