# # [BangBangPulseProblem](@id bang-bang-pulse)
#
# `BangBangPulseProblem` promotes **bang-bang** (piecewise-constant, few-switch) controls by penalizing ``\|du\|_1`` via an **exact slack reformulation**. Unlike `SmoothPulseProblem` (which uses 2 derivative levels with L2 regularization), this stores only 1 derivative (`du`) and uses slack variables to impose an exact L1 penalty, promoting sparsity in `du` and thus fewer switches.
#
# ## When to Use
#
# Use `BangBangPulseProblem` when:
# - You want piecewise constant control pulses with minimal switching
# - You want to promote bang-bang style controls (long constant segments)
# - You prefer exact L1 regularization over smooth L2 approximations
#
# ## Comparison with SmoothPulseProblem
#
# | | SmoothPulseProblem | BangBangPulseProblem |
# |---|---|---|
# | Derivatives stored | `du`, `ddu` | `du` only |
# | Regularization on `du` | L2 (`QuadraticRegularizer`) | L1 (`LinearRegularizer` on slacks) |
# | Regularization on `u` | L2 | L2 (same) |
# | Extra variables | — | slack `s_du ≥ 0` |
# | Extra constraints | — | `L1SlackConstraint`: ``\|du\| \leq s`` |
# | Bound params | `du_bound`, `ddu_bound` | `du_bound` only |
#
# ## How the L1 Penalty Works
#
# Instead of a smooth approximation, the L1 penalty uses an exact slack reformulation. Slack variables ``s \geq 0`` (same dimension as `du`) are introduced with the constraint:
# ```math
# |du_{k,i}| \leq s_{k,i}
# ```
# Then the linear cost ``R_{du} \sum_k \sum_i s_{k,i} \Delta t_k`` is minimized. At optimality, ``s = |du|``, giving the exact L1 norm.
#
# ## Pulse Requirement
#
# `BangBangPulseProblem` requires a trajectory with a `ZeroOrderPulse`:
#
# ```julia
# pulse = ZeroOrderPulse(controls, times)
# qtraj = UnitaryTrajectory(sys, pulse, U_goal)
# qcp = BangBangPulseProblem(qtraj, N)  # Works
# ```
#
# Using a spline pulse will result in an error directing you to `SplinePulseProblem`.
#
# ## Constructor
#
# ```julia
# BangBangPulseProblem(
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
# | `R` | `Float64` | `1e-2` | Base regularization weight applied to all terms. |
# | `R_u` | `Float64` or `Vector{Float64}` | `R` | L2 regularization on control values. Can be per-drive. |
# | `R_du` | `Float64` or `Vector{Float64}` | `R` | L1 weight on first derivative (applied to slacks). Higher values produce fewer switches. |
#
# ### Bounds
#
# | Parameter | Type | Default | Description |
# |-----------|------|---------|-------------|
# | `du_bound` | `Float64` | `Inf` | Maximum allowed control jump between timesteps. |
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
using CairoMakie

## Define system
H_drift = PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
sys = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

## Create trajectory
T, N = 10.0, 100
times = collect(range(0, T, length = N))
pulse = ZeroOrderPulse(0.1 * randn(2, N), times)
qtraj = UnitaryTrajectory(sys, pulse, GATES[:X])

## Solve with L1 regularization on du
qcp = BangBangPulseProblem(qtraj, N; Q = 100.0, R_du = 1e-1)
cached_solve!(qcp, "bang_bang_basic"; max_iter = 100)

fidelity(qcp)

# ## Visualize the optimized trajectory
traj = get_trajectory(qcp)
fig = plot(traj)

# ### Tuning the L1 Weight
#
# The `R_du` parameter controls how aggressively switches are penalized. Higher values produce fewer switches (more bang-bang):

qcp = BangBangPulseProblem(
    qtraj,
    N;
    Q = 100.0,
    R_du = 1.0,    ## Strong L1 penalty → fewer switches
)

# ### With Derivative Bounds
#
# Constrain the maximum control jump size:

qcp = BangBangPulseProblem(
    qtraj,
    N;
    Q = 100.0,
    R_du = 1e-1,
    du_bound = 0.5,    ## Limit control jumps
)

# ### Per-Drive Regularization
#
# Apply different L1 weights to different control channels:

qcp = BangBangPulseProblem(
    qtraj,
    N;
    Q = 100.0,
    R_u = [1e-3, 1e-2],     ## Less L2 regularization on drive 1
    R_du = [1e-1, 1.0],     ## Stronger L1 on drive 2
)

# ### State Transfer

ψ_init = ComplexF64[1, 0]  ## |0⟩
ψ_goal = ComplexF64[0, 1]  ## |1⟩

pulse = ZeroOrderPulse(0.1 * randn(2, N), times)
qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

qcp = BangBangPulseProblem(qtraj, N; Q = 100.0, R_du = 1e-1)
cached_solve!(qcp, "bang_bang_state_transfer"; max_iter = 100)

# ### Multiple State Transfers
#
# Use `MultiKetTrajectory` for gates defined by state mappings:

ψ0, ψ1 = ComplexF64[1, 0], ComplexF64[0, 1]

## X gate: |0⟩ → |1⟩ and |1⟩ → |0⟩
pulse = ZeroOrderPulse(0.1 * randn(2, N), times)
qtraj = MultiKetTrajectory(sys, pulse, [ψ0, ψ1], [ψ1, ψ0])

qcp = BangBangPulseProblem(qtraj, N; Q = 100.0, R_du = 1e-1)
cached_solve!(qcp, "bang_bang_multi_ket"; max_iter = 100)

# ## How It Works
#
# `BangBangPulseProblem` internally:
#
# 1. **Adds 1 derivative variable**: Creates `:du` (first derivative) alongside controls `:u` (compared to 2 in `SmoothPulseProblem`)
#
# 2. **Adds slack variables**: Creates `:s_du` (non-negative, same dimension as `:du`), initialized at `|du|`
#
# 3. **Sets up integrators**: Uses `DerivativeIntegrator` to enforce:
#    - `du[k] = (u[k+1] - u[k]) / Δt`
#
# 4. **Configures objectives**:
#    - Infidelity objective with weight `Q`
#    - Quadratic (L2) regularization on `u` with weight `R_u`
#    - Linear (L1) regularization on `s_du` with weight `R_du`
#
# 5. **Applies constraints**: `L1SlackConstraint` enforces `|du| ≤ s_du`, and `du_bound` as hard constraints
#
# At optimality the slacks satisfy `s = |du|`, so the linear cost on `s` equals the exact L1 norm of `du`.
#
# ## See Also
#
# - [SmoothPulseProblem](@ref smooth-pulse) - For smooth (L2-regularized) controls
# - [SplinePulseProblem](@ref spline-pulse) - For inherently smooth spline-based controls
# - [MinimumTimeProblem](@ref minimum-time) - Time-optimal control
