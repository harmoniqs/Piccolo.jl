# # [Constraints](@id constraints-concept)
#
# Constraints are hard requirements enforced exactly by the NLP solver (Ipopt).
# Unlike objectives (which are minimized), constraints must be satisfied at
# every feasible solution.
#
# ## Overview
#
# The full NLP (see [Concepts Overview](@ref concepts-overview)) has three
# classes of constraints:
#
# | Class | Mathematical Form | Piccolo Source |
# |-------|-------------------|----------------|
# | **Dynamics** (equality) | ``x_{k+1} = \exp(\Delta t_k\, G(\boldsymbol{u}_k))\, x_k`` | `BilinearIntegrator` |
# | **Box bounds** (inequality) | ``\boldsymbol{u}_{\min} \leq \boldsymbol{u}_k \leq \boldsymbol{u}_{\max}`` | `QuantumSystem.drive_bounds` |
# | **Custom** (equality / inequality) | fidelity floors, leakage ceilings, … | `PiccoloOptions`, manual |
#
# ## Bound Constraints
#
# ### Control Bounds
#
# Control bounds ``u_i^{\min} \leq u_{i,k} \leq u_i^{\max}`` are inherited
# from the `QuantumSystem` and enforced at every knot point:

using Piccolo

## Bounds specified at system creation
H_drift = PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
drive_bounds = [1.0, 0.5]  # Drive 1: ±1.0, Drive 2: ±0.5
sys = QuantumSystem(H_drift, H_drives, drive_bounds)

# ### Derivative Bounds
#
# Derivative bounds limit how fast controls can change.  For
# `SmoothPulseProblem` the discrete differences
# ``\Delta\boldsymbol{u}_k = \boldsymbol{u}_k - \boldsymbol{u}_{k-1}`` and
# ``\Delta^2\boldsymbol{u}_k`` are bounded element-wise:
#
# ```math
# |\Delta u_{i,k}| \leq b_{du}, \qquad |\Delta^2 u_{i,k}| \leq b_{ddu}
# ```

T, N = 10.0, 100
times = collect(range(0, T, length = N))
pulse = ZeroOrderPulse(0.1 * randn(2, N), times)
U_goal = GATES[:X]
qtraj = UnitaryTrajectory(sys, pulse, U_goal)

qcp = SmoothPulseProblem(
    qtraj,
    N;
    du_bound = 0.5,    # Max control jump per timestep
    ddu_bound = 0.1,   # Max control acceleration
)
cached_solve!(qcp, "constraints_bounds"; max_iter = 50)
fidelity(qcp)

# ### Timestep Bounds
#
# For free-time optimization the timestep ``\Delta t_k`` becomes a decision
# variable with box constraints:
#
# ```math
# \Delta t_{\min} \leq \Delta t_k \leq \Delta t_{\max}
# ```
#
# ```julia
# qcp = SmoothPulseProblem(qtraj, N; Δt_bounds=(0.01, 0.5))
# ```
#
# ## Fidelity Constraints
#
# Used with `MinimumTimeProblem` to enforce a minimum gate quality while
# minimizing total time ``T = \sum_k \Delta t_k``:
#
# ```math
# \min_{\boldsymbol{u},\,\Delta t} \;\sum_{k=1}^{N} \Delta t_k
# \quad \text{s.t.} \quad F(x_N) \geq F_{\min}
# ```
#
# ### Automatic Setup in MinimumTimeProblem

## First solve a base problem with variable timesteps
qcp_base = SmoothPulseProblem(qtraj, N; Δt_bounds = (0.01, 0.5))
cached_solve!(qcp_base, "constraints_base_freetime"; max_iter = 100)
fidelity(qcp_base)

## Automatically adds FinalUnitaryFidelityConstraint
qcp_mintime = MinimumTimeProblem(qcp_base; final_fidelity = 0.99)
cached_solve!(qcp_mintime, "constraints_mintime"; max_iter = 100)
fidelity(qcp_mintime)

# ## Leakage Constraints
#
# For multilevel systems where ``d > d_{\text{comp}}``, leakage to
# non-computational states can be bounded:
#
# ```math
# 1 - \sum_{j \in \mathcal{S}_{\text{comp}}} \langle j | \rho_N | j \rangle \;\leq\; \epsilon_{\text{leak}}
# ```

opts = PiccoloOptions(leakage_constraint = true, leakage_constraint_value = 1e-3)

## Example with a transmon system
sys_transmon = TransmonSystem(levels = 3, δ = 0.2, drive_bounds = [0.2, 0.2])
pulse_t = ZeroOrderPulse(0.01 * randn(2, N), times)
U_X = EmbeddedOperator(:X, sys_transmon)
qtraj_t = UnitaryTrajectory(sys_transmon, pulse_t, U_X)

qcp_leak = SmoothPulseProblem(qtraj_t, N; piccolo_options = opts)
cached_solve!(qcp_leak, "constraints_leakage"; max_iter = 100)
fidelity(qcp_leak)

# ## PiccoloOptions
#
# `PiccoloOptions` provides a convenient way to configure common constraint
# and objective settings:
#
# | Option | Type | Default | Description |
# |--------|------|---------|-------------|
# | `leakage_constraint` | `Bool` | `false` | Enable leakage constraint |
# | `leakage_constraint_value` | `Float64` | `1e-3` | Maximum leakage ``\epsilon_{\text{leak}}`` |
# | `leakage_cost` | `Float64` | `0.0` | Leakage objective weight |
# | `timesteps_all_equal` | `Bool` | `false` | Force uniform timesteps |
# | `verbose` | `Bool` | `false` | Print solver progress |
#
# ## Constraints vs Objectives
#
# | Use Case | Constraint | Objective |
# |----------|------------|-----------|
# | Must achieve ``F \geq 0.99`` | Fidelity constraint | |
# | Prefer higher fidelity | | Infidelity objective |
# | Control must be ``\leq 1.0`` | Bound constraint | |
# | Prefer smaller controls | | Regularization |
# | Leakage must be ``< 10^{-3}`` | Leakage constraint | |
# | Prefer less leakage | | Leakage objective |
#
# **Trade-off**: Constraints guarantee satisfaction (if feasible) but can make
# the problem harder to solve.  Objectives are softer but offer no guarantees.
# A practical strategy is to start with objectives alone, then add constraints
# for hard requirements.
#
# ## Constraint Feasibility
#
# Common causes of infeasibility:
# 1. Target fidelity too high for the given gate time
# 2. Insufficient time for the gate (quantum speed limit)
# 3. Control bounds too tight
# 4. Conflicting constraints (e.g., low leakage + high fidelity + short time)
#
# ## Best Practices
#
# ### 1. Start Without Constraints

## First, find a good solution with just objectives
qcp_simple = SmoothPulseProblem(UnitaryTrajectory(sys, pulse, U_goal), N; Q = 100.0)
cached_solve!(qcp_simple, "constraints_simple"; max_iter = 100)
fidelity(qcp_simple)

# ### 2. Add Constraints Gradually
#
# ```julia
# if fidelity(qcp_simple) > 0.99
#     qcp_constrained = MinimumTimeProblem(qcp_simple; final_fidelity=0.99)
#     solve!(qcp_constrained)
# end
# ```
#
# ### 3. Use Margin for Robustness
#
# ```julia
# # Target slightly higher than needed
# qcp = MinimumTimeProblem(qcp_base; final_fidelity=0.995)  # Want 0.99
# ```
#
# ## See Also
#
# - [Objectives](@ref objectives-concept) - Soft optimization targets
# - [Problem Templates](@ref problem-templates-overview) - Using constraints in practice
# - [MinimumTimeProblem](@ref minimum-time) - Fidelity-constrained time optimization
