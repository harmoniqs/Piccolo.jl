# ```@copybutton
# literate/problem-templates/spline_pulse.jl
# ```
#
# # [SplinePulseProblem](@id spline-pulse)
#
# `SplinePulseProblem` sets up trajectory optimization with spline-based pulses where derivative variables represent spline slopes or tangents. This is ideal for inherently smooth control pulses and warm-starting from previous solutions.
#
# ## When to Use
#
# Use `SplinePulseProblem` when:
# - You need inherently smooth control pulses
# - You're warm-starting from a previously optimized solution
# - You want cubic spline smoothness without derivative regularization
# - You're working with hardware that expects smooth pulse shapes
#
# ## Pulse Requirements
#
# `SplinePulseProblem` works with spline pulse types:
#
# | Pulse Type | Derivative Meaning |
# |------------|-------------------|
# | `LinearSplinePulse` | `:du` represents constrained slope between knots |
# | `CubicSplinePulse` | `:du` represents independent Hermite tangents |
#
# ```julia
# # Linear spline
# pulse = LinearSplinePulse(controls, times)
# qtraj = UnitaryTrajectory(sys, pulse, U_goal)
# qcp = SplinePulseProblem(qtraj)  # spline-faithful default (native SplineIntegrator)
#
# # Cubic spline — the default is the correct pairing since #334
# pulse = CubicSplinePulse(controls, tangents, times)
# qtraj = UnitaryTrajectory(sys, pulse, U_goal)
# qcp = SplinePulseProblem(qtraj)  # spline-faithful by default
# # PWC escape hatch, when explicitly wanted: qcp = SplinePulseProblem(qtraj; integrator_type = :pwc)
# ```
#
# ## Constructor Variants
#
# ### Use Native Knot Times (Recommended for Warm-Starting)
#
# ```julia
# SplinePulseProblem(qtraj::AbstractQuantumTrajectory{<:AbstractSplinePulse}; kwargs...)
# ```
#
# Uses the pulse's native knot times without resampling. Best for warm-starting from a previous solution.
#
# ### Resample to N Timesteps
#
# ```julia
# SplinePulseProblem(qtraj::AbstractQuantumTrajectory{<:AbstractSplinePulse}, N::Int; kwargs...)
# ```
#
# Resamples the pulse to `N` uniformly spaced timesteps.
#
# ### Use Specific Times
#
# ```julia
# SplinePulseProblem(qtraj::AbstractQuantumTrajectory{<:AbstractSplinePulse}, times::AbstractVector; kwargs...)
# ```
#
# Resamples the pulse to the specified time points.
#
# ## Parameter Reference
#
# ### Objective Weights
#
# | Parameter | Type | Default | Description |
# |-----------|------|---------|-------------|
# | `Q` | `Float64` | `100.0` | Weight on infidelity objective |
# | `R` | `Float64` | `1e-2` | Base regularization weight |
# | `R_u` | `Float64` or `Vector{Float64}` | `R` | Regularization on control values |
# | `R_du` | `Float64` or `Vector{Float64}` | `R` | Regularization on derivatives/tangents |
#
# ### Bounds
#
# | Parameter | Type | Default | Description |
# |-----------|------|---------|-------------|
# | `du_bound` | `Float64` | `Inf` | Uniform maximum derivative/slope bound for all drives. |
# | `du_bounds` | `Vector{Float64}` | `nothing` | Per-drive maximum derivative/slope bounds (takes precedence over `du_bound`). |
# | `Δt_bounds` | `Tuple{Float64, Float64}` | `nothing` | Time-step bounds for free-time optimization |
#
# ### Advanced Options
#
# | Parameter | Type | Default | Description |
# |-----------|------|---------|-------------|
# | `integrator` | `AbstractIntegrator` | `nothing` | Custom integrator. If `nothing`, uses the native spline-faithful `SplineIntegrator` (the #334 default — globals auto-detected from `sys.global_params`). `integrator_type = :pwc` requests the piecewise-constant escape hatch; see below. |
# | `global_names` | `Vector{Symbol}` | `nothing` | Global parameters to optimize |
# | `global_bounds` | `Dict{Symbol, ...}` | `nothing` | Bounds on global variables |
# | `constraints` | `Vector{AbstractConstraint}` | `[]` | Additional constraints |
# | `piccolo_options` | `PiccoloOptions` | `PiccoloOptions()` | Solver options |
# | `free_phase` | `Bool` | `false` | Optimize a per-subsystem frame phase alongside the pulse. |
# | `initial_phases` | `Vector{Float64}` | `nothing` | Initial values for the per-subsystem phase variables when `free_phase=true`. |
#
# ## Examples
#
# ### Basic Spline Optimization

using Piccolo

## Define system
H_drift = PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
sys = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

## Create cubic spline pulse
T, N = 10.0, 50
times = collect(range(0, T, length = N))
controls = 0.1 * randn(2, N)
tangents = zeros(2, N)  ## Initial tangents
pulse = CubicSplinePulse(controls, tangents, times)

qtraj = UnitaryTrajectory(sys, pulse, GATES[:X])

## Solve using native knot times. Since #334 the default integrator is the
## native spline-faithful SplineIntegrator — the optimized waveform is the
## spline the pulse object describes.
qcp = SplinePulseProblem(qtraj; Q = 100.0, du_bound = 10.0)
cached_solve!(qcp, "spline_pulse_basic"; max_iter = 100)

# !!! note "The default integrator is spline-faithful — the PWC path is the escape hatch"
#     Since #334 the default backend is Piccolo's native `SplineIntegrator`: the
#     collocation integrates the same interpolated waveform the pulse object
#     describes. The piecewise-constant `BilinearIntegrator` remains available as
#     the **explicit** `integrator_type = :pwc` choice — and it *warns* for
#     `CubicSplinePulse`, because it never reads the Hermite tangents (`:du`):
#     the optimizer would minimize against a different waveform than the one your
#     pulse actually produces. Historically that backend was the silent default
#     (measured: optimizer-reported ~1e-8 infidelity for pulses actually achieving
#     ~1e-3, before #275 turned it into an error and #334 replaced it with the
#     correct default). Declare it only when you mean it, and then check the
#     divergence:
#
#     This is easy to *see* rather than take on faith. `rollout_divergence` compares the
#     optimizer's collocation solution against an ODE re-rollout of the same pulse at the final
#     time, and `sync_trajectory!` warns automatically when the two part company:

rollout_divergence(qcp)

# Reference points, so you can judge the number above: a *correctly* paired
# spline-pulse-plus-spline-integrator problem measures `3.5e-7 … 6.4e-4`, while a spline pulse
# under a piecewise-constant integrator measures `5.9e-2 … 1.9e-1` — roughly a 90× gap, which is
# why `ROLLOUT_DIVERGENCE_RTOL[]` defaults to `5e-3`. Refining the grid does **not** help; the
# divergence *grows* with `N`, because a finer collocation grid gives the optimizer more room to
# exploit the model.
#
# `verify` says the same thing in fidelity terms — what the optimizer believes, versus what the
# pulse achieves:

verify(qcp)

# To see the trap the old default set, declare the PWC backend explicitly and
# compare — this is the acknowledged escape hatch in action:

qcp_pwc = SplinePulseProblem(qtraj; Q = 100.0, du_bound = 10.0, integrator_type = :pwc)
cached_solve!(qcp_pwc, "spline_pulse_pwc_escape_hatch"; max_iter = 100)
verify(qcp_pwc)

# Under the PWC backend the optimizer is confident to six digits
# (`F_optimizer ≈ 0.999998`) about a pulse whose rollout misses by roughly
# `0.23` in absolute fidelity — the `F_optimizer` figure is the one the old
# default would have shown you. Prefer `F_rollout`, and keep the
# spline-faithful default unless you specifically need PWC.

# !!! note "Spline-integrator variants"
#     The default `SplineIntegrator` ships in Piccolo (analytic knot-aware
#     sensitivities, adaptive Tsit5 default). Parallel and Magnus-scheme variants
#     (`SplineIntegrator(qtraj, N; alg = MagnusGL4Alg(...))`, ...) live in the
#     Piccolissimo stack — pass those explicitly via `integrator = ...`.

# ### Warm-Starting from Previous Solution
#
# ```julia
# # Load previously optimized pulse
# using JLD2
# @load "optimized_pulse.jld2" saved_pulse
#
# # Create new trajectory with saved pulse
# qtraj = UnitaryTrajectory(sys, saved_pulse, U_goal)
#
# # Use native knot times (no resampling) — spline-faithful default
# qcp = SplinePulseProblem(qtraj)
# solve!(qcp; max_iter=50)  # Converges quickly from good initial guess
# ```
#
# ### Resampling to Different Resolution
#
# The original pulse above has 50 knots. We can resample to 100 for finer control:

qcp_resampled = SplinePulseProblem(qtraj, 100; Q = 100.0)
cached_solve!(qcp_resampled, "spline_pulse_resampled"; max_iter = 100)

# ### Linear vs Cubic Splines
#
# **Linear Splines**: The derivative `:du` represents the slope between knots. A `DerivativeIntegrator` constraint enforces `du[k] = (u[k+1] - u[k]) / Δt`.

pulse_linear = LinearSplinePulse(controls, times)
qtraj_linear = UnitaryTrajectory(sys, pulse_linear, GATES[:X])
qcp_linear = SplinePulseProblem(qtraj_linear)
## du is constrained to be consistent with u

# **Cubic Splines (Hermite)**: The derivative `:du` represents independent Hermite tangents at each knot. These are free optimization variables with no inter-knot constraint.

pulse_cubic = CubicSplinePulse(controls, tangents, times)
qtraj_cubic = UnitaryTrajectory(sys, pulse_cubic, GATES[:X])
qcp_cubic = SplinePulseProblem(qtraj_cubic)  # spline-faithful default — the correct cubic pairing
## du (tangents) are independent variables

# ### Per-Drive Derivative Bounds
#
# Set independent slope limits for each drive, e.g. when hardware channels
# have different slew-rate constraints:

qcp_per_drive = SplinePulseProblem(
    qtraj;
    Q = 100.0,
    du_bounds = [5.0, 2.0],  ## drive 1 allows faster slopes than drive 2
)

# ## Trajectory Structure
#
# Unlike `SmoothPulseProblem` which has three derivative levels (`:u`, `:du`, `:ddu`), `SplinePulseProblem` only has one:
#
# | Variable | Description |
# |----------|-------------|
# | `:u` | Control values at knot points |
# | `:du` | Derivatives/tangents (meaning depends on spline type) |
#
# The reduced number of variables can lead to faster optimization while maintaining smooth pulses through the spline interpolation.
#
# ## See Also
#
# - [SmoothPulseProblem](@ref smooth-pulse) - For piecewise constant controls
# - [BangBangPulseProblem](@ref bang-bang-pulse) - For L1-regularized bang-bang controls
# - [MinimumTimeProblem](@ref minimum-time) - Time-optimal control
# - [Pulses](@ref pulses-concept) - Detailed pulse type documentation
