# # [SmoothPulseProblem](@id smooth-pulse)
#
# `SmoothPulseProblem` is the most commonly used problem template in Piccolo.jl. It sets up trajectory optimization with piecewise constant controls (`ZeroOrderPulse`) where control smoothness is enforced through discrete derivative variables.
#
# ## When to Use
#
# Use `SmoothPulseProblem` when:
# - You want piecewise constant control pulses
# - You need smooth transitions between control values
# - You're starting a new optimization (not warm-starting from a previous solution)
# - You want fine control over derivative bounds
#
# ## Pulse Requirement
#
# `SmoothPulseProblem` requires a trajectory with a `ZeroOrderPulse`:
#
# ```julia
# pulse = ZeroOrderPulse(controls, times)
# qtraj = UnitaryTrajectory(sys, pulse, U_goal)
# qcp = SmoothPulseProblem(qtraj, N)  # Works
# ```
#
# Using a spline pulse will result in an error directing you to `SplinePulseProblem`.
#
# ## Constructor
#
# ```julia
# SmoothPulseProblem(
#     qtraj::AbstractQuantumTrajectory{<:ZeroOrderPulse},
#     N::Int;
#     kwargs...
# )
# ```
#
# ## Parameter Reference
#
# ### Required Parameters
#
# | Parameter | Type | Description |
# |-----------|------|-------------|
# | `qtraj` | `AbstractQuantumTrajectory{ZeroOrderPulse}` | Quantum trajectory containing system, pulse, and goal |
# | `N` | `Int` | Number of discretization timesteps |
#
# ### Objective Weights
#
# | Parameter | Type | Default | Description |
# |-----------|------|---------|-------------|
# | `Q` | `Float64` | `100.0` | Weight on infidelity objective. Higher values prioritize achieving target fidelity. |
# | `R` | `Float64` | `1e-2` | Base regularization weight applied to all derivative terms. |
# | `R_u` | `Float64` or `Vector{Float64}` | `R` | Regularization on control values. Can be per-drive. |
# | `R_du` | `Float64` or `Vector{Float64}` | `R` | Regularization on first derivative (control jumps). |
# | `R_ddu` | `Float64` or `Vector{Float64}` | `R` | Regularization on second derivative (control acceleration). |
#
# ### Bounds
#
# | Parameter | Type | Default | Description |
# |-----------|------|---------|-------------|
# | `du_bound` | `Float64` | `Inf` | Maximum allowed control jump between timesteps. |
# | `ddu_bound` | `Float64` | `1.0` | Maximum allowed control acceleration. |
# | `Δt_bounds` | `Tuple{Float64, Float64}` | `nothing` | Time-step bounds `(min, max)` for free-time optimization. Required for `MinimumTimeProblem`. |
#
# ### Advanced Options
#
# | Parameter | Type | Default | Description |
# |-----------|------|---------|-------------|
# | `integrator` | `AbstractIntegrator` | `nothing` | Custom integrator. If `nothing`, uses `BilinearIntegrator`. |
# | `global_names` | `Vector{Symbol}` | `nothing` | Names of global parameters to optimize (requires custom integrator). |
# | `global_bounds` | `Dict{Symbol, ...}` | `nothing` | Bounds on global variables. Values can be `Float64` (symmetric ±) or `Tuple{Float64, Float64}`. |
# | `constraints` | `Vector{AbstractConstraint}` | `[]` | Additional constraints to add to the problem. |
# | `piccolo_options` | `PiccoloOptions` | `PiccoloOptions()` | Solver and leakage options. |
#
# ## Examples
#
# ### Basic Gate Synthesis

using Piccolo

## Define system
H_drift = PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
sys = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

## Create trajectory
T, N = 10.0, 100
times = collect(range(0, T, length = N))
pulse = ZeroOrderPulse(0.1 * randn(2, N), times)
qtraj = UnitaryTrajectory(sys, pulse, GATES[:X])

## Solve
qcp = SmoothPulseProblem(qtraj, N; Q = 100.0, R = 1e-2)
solve!(qcp; max_iter = 100)

println("Fidelity: ", fidelity(qcp))

# ### With Derivative Bounds
#
# Constrain how fast controls can change:

qcp = SmoothPulseProblem(
    qtraj, N;
    Q=100.0,
    du_bound=0.5,    ## Limit control jumps
    ddu_bound=0.1    ## Limit control acceleration
)

# ### Enabling Free Time for MinimumTimeProblem
#
# To later use `MinimumTimeProblem`, give bounds on variable timesteps:
#
qcp = SmoothPulseProblem(
    qtraj, N;
    Q=100.0,
    Δt_bounds=(0.01, 0.5)  ## Timesteps can vary from 0.01 to 0.5
)
solve!(qcp; max_iter=100)

## Now can minimize time
qcp_mintime = MinimumTimeProblem(qcp; final_fidelity=0.99)
solve!(qcp_mintime; max_iter=100)

# ### Per-Drive Regularization
#
# Apply different regularization to different control channels:

qcp = SmoothPulseProblem(
    qtraj, N;
    Q=100.0,
    R_u=[1e-3, 1e-2],     ## Less regularization on drive 1
    R_du=[1e-2, 1e-1],    ## Different smoothness weights
    R_ddu=[1e-2, 1e-1]
)

# ### With Leakage Suppression
#
# For multilevel systems, use `PiccoloOptions` to enable leakage handling, in this example a 3-level transmon system:
sys = TransmonSystem(levels=3, δ=0.2, drive_bounds=[0.2, 0.2])

## Embedded X gate (only affects computational subspace)
U_goal = EmbeddedOperator(:X, sys)

pulse = ZeroOrderPulse(0.1 * randn(2, N), times)
qtraj = UnitaryTrajectory(sys, pulse, U_goal)

## Enable leakage suppression
opts = PiccoloOptions(
    leakage_constraint=true,
    leakage_constraint_value=1e-3
)

qcp = SmoothPulseProblem(qtraj, N; Q=100.0, piccolo_options=opts)
solve!(qcp; max_iter=100)

# ### State Transfer
#
# Works with `KetTrajectory` for state preparation. Here we go back to the 2-level system:
#
sys = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

ψ_init = ComplexF64[1, 0]  ## |0⟩
ψ_goal = ComplexF64[0, 1]  ## |1⟩

pulse = ZeroOrderPulse(0.1 * randn(2, N), times)
qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

qcp = SmoothPulseProblem(qtraj, N; Q=100.0)
solve!(qcp; max_iter=100)

# ### Multiple State Transfers
#
# Use `MultiKetTrajectory` for gates defined by state mappings:

ψ0, ψ1 = ComplexF64[1, 0], ComplexF64[0, 1]

## X gate: |0⟩ → |1⟩ and |1⟩ → |0⟩
qtraj = MultiKetTrajectory(sys, pulse, [ψ0, ψ1], [ψ1, ψ0])

qcp = SmoothPulseProblem(qtraj, N; Q=100.0)
solve!(qcp; max_iter=100)

# ## How It Works
#
# `SmoothPulseProblem` internally:
#
# 1. **Adds derivative variables**: Creates `:du` (first derivative) and `:ddu` (second derivative) variables alongside controls `:u`
#
# 2. **Sets up integrators**: Uses `DerivativeIntegrator` to enforce consistency:
#    - `du[k] = (u[k+1] - u[k]) / Δt`
#    - `ddu[k] = (du[k+1] - du[k]) / Δt`
#
# 3. **Configures objectives**:
#    - Infidelity objective with weight `Q`
#    - Quadratic regularization on `u`, `du`, `ddu` with weights `R_u`, `R_du`, `R_ddu`
#
# 4. **Applies bounds**: Enforces `du_bound` and `ddu_bound` as hard constraints
#
# The derivative variables act as auxiliary optimization variables that encourage smooth control transitions without requiring explicit smoothness constraints on the original controls.
#
# ## See Also
#
# - [SplinePulseProblem](@ref spline-pulse) - For inherently smooth spline-based controls
# - [MinimumTimeProblem](@ref minimum-time) - Time-optimal control
# - [SamplingProblem](@ref sampling) - Robust optimization
# - [Pulses](@ref pulses-concept) - Pulse type documentation
