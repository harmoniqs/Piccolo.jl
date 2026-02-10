# # [Constraints](@id constraints-concept)
#
# Constraints define hard requirements that optimized solutions must satisfy.
# Unlike objectives (which are minimized), constraints are enforced exactly.
#
# ## Overview
#
# Piccolo.jl supports several constraint types:
#
# | Constraint Type | Description |
# |-----------------|-------------|
# | Bound constraints | Limits on variable values |
# | Fidelity constraints | Minimum fidelity requirements |
# | Leakage constraints | Maximum allowed leakage |
# | Equality constraints | Custom equality requirements |
#
# ## Bound Constraints
#
# ### Control Bounds
#
# Control bounds are specified in the `QuantumSystem` and automatically enforced:

using Piccolo

## Bounds specified at system creation
H_drift = PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
drive_bounds = [1.0, 0.5]  # Drive 1: ±1.0, Drive 2: ±0.5
sys = QuantumSystem(H_drift, H_drives, drive_bounds)

# ### Derivative Bounds
#
# Derivative bounds limit how fast controls can change:

T, N = 10.0, 100
times = collect(range(0, T, length = N))
pulse = ZeroOrderPulse(0.1 * randn(2, N), times)
U_goal = GATES[:X]
qtraj = UnitaryTrajectory(sys, pulse, U_goal)

qcp = SmoothPulseProblem(
    qtraj, N;
    du_bound = 0.5,    # Max control jump per timestep
    ddu_bound = 0.1,   # Max control acceleration
)
cached_solve!(qcp, "constraints_bounds"; max_iter = 50)
fidelity(qcp)

# ### Timestep Bounds
#
# For free-time optimization:
#
# ```julia
# qcp = SmoothPulseProblem(
#     qtraj, N;
#     Δt_bounds=(0.01, 0.5)  # Min and max timestep
# )
# ```
#
# ## Fidelity Constraints
#
# Used with `MinimumTimeProblem` to enforce minimum gate quality.
#
# ### FinalUnitaryFidelityConstraint
#
# ```julia
# constraint = FinalUnitaryFidelityConstraint(
#     :Ũ⃗,        # State variable name
#     U_goal,    # Target unitary
#     0.99       # Minimum fidelity
# )
# ```
#
# ### Automatic Setup in MinimumTimeProblem
#
# You typically don't create fidelity constraints manually.
# `MinimumTimeProblem` adds them automatically:

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
# ### Via PiccoloOptions
#
# The easiest approach to handle leakage:

opts = PiccoloOptions(
    leakage_constraint = true,
    leakage_constraint_value = 1e-3,
)

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
# settings:
#
# ```julia
# opts = PiccoloOptions(
#     # Leakage handling
#     leakage_constraint=true,
#     leakage_constraint_value=1e-3,
#     leakage_cost=10.0,  # Also adds objective
#
#     # Timestep handling
#     timesteps_all_equal=false,
#
#     # Verbosity
#     verbose=true
# )
# ```
#
# ### Key Options
#
# | Option | Type | Default | Description |
# |--------|------|---------|-------------|
# | `leakage_constraint` | `Bool` | `false` | Enable leakage constraint |
# | `leakage_constraint_value` | `Float64` | `1e-3` | Maximum leakage |
# | `leakage_cost` | `Float64` | `0.0` | Leakage objective weight |
# | `timesteps_all_equal` | `Bool` | `false` | Force uniform timesteps |
# | `verbose` | `Bool` | `false` | Print solver progress |
#
# ## Constraints vs Objectives
#
# Understanding when to use each:
#
# | Use Case | Constraint | Objective |
# |----------|------------|-----------|
# | Must achieve F ≥ 0.99 | ✓ Fidelity constraint | |
# | Prefer higher fidelity | | ✓ Infidelity objective |
# | Control must be ≤ 1.0 | ✓ Bound constraint | |
# | Prefer smaller controls | | ✓ Regularization |
# | Leakage must be < 1e-3 | ✓ Leakage constraint | |
# | Prefer less leakage | | ✓ Leakage objective |
#
# ### Trade-offs
#
# **Constraints:**
# - Guarantee satisfaction (if feasible)
# - Can make problem harder to solve
# - May be infeasible
#
# **Objectives:**
# - More flexible
# - Easier optimization
# - No guarantees
#
# ### Recommendation
#
# Start with objectives, add constraints for hard requirements.
#
# ## Constraint Feasibility
#
# ### Common Causes of Infeasibility
#
# 1. **Fidelity too high**: Target fidelity may not be achievable
# 2. **Time too short**: Insufficient time for the gate
# 3. **Bounds too tight**: Controls can't reach required values
# 4. **Conflicting constraints**: e.g., low leakage with high fidelity in short time
#
# ### Solutions
#
# 1. **Relax constraints**: Lower fidelity target, increase time
# 2. **Better initialization**: Start from a working solution
# 3. **Adjust bounds**: Allow larger controls or longer time
# 4. **Use objectives first**: Find a good solution, then add constraints
#
# ## Best Practices
#
# ### 1. Start Without Constraints

## First, find a good solution with just objectives
qcp_simple = SmoothPulseProblem(
    UnitaryTrajectory(sys, pulse, U_goal), N;
    Q = 100.0,
)
cached_solve!(qcp_simple, "constraints_simple"; max_iter = 100)
fidelity(qcp_simple)

# ### 2. Add Constraints Gradually
#
# ```julia
# # Then add constraints if needed
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
