# # [Global Variables](@id global-variables)
#
# Global variables allow you to optimize system parameters alongside control pulses.
# This is useful for calibrating system properties or finding optimal operating points.

# ## When to Use Global Variables
#
# - **System calibration**: Finding optimal qubit frequencies or coupling strengths
# - **Gate design**: Optimizing hardware parameters for better gates
# - **Sensitivity analysis**: Understanding parameter effects

# ## Basic Concept
#
# Standard optimization only varies control pulses ``u(t)``. With global variables,
# you can also optimize system parameters that appear in the Hamiltonian.
#
# ```math
# H(u, t; \theta) = H_{\text{drift}}(\theta) + \sum_i u_i(t) H_{\text{drive},i}(\theta)
# ```
#
# Where ``\theta`` are global parameters to optimize.

# ## Setup Requirements
#
# Global variable optimization requires:
# 1. A `QuantumSystem` with `global_params`
# 2. A **custom integrator** that supports globals (e.g., from Piccolissimo)
# 3. **Global bounds** specification
#
# Piccolo's built-in `BilinearIntegrator` does **not** support global variables —
# it has no global-aware Jacobian columns or Hessian blocks. A custom integrator
# from Piccolissimo handles the extended control vector `[controls..., globals...]`
# and provides the correct derivative information to the optimizer.
#
# This guide demonstrates the global variable API using Piccolo's built-in
# integrator. The globals are stored in the trajectory and bounded, but are not
# coupled to the dynamics. For full global optimization, use a Piccolissimo
# integrator.

# ## Defining a System with Global Parameters
#
# Global parameters are stored on the `QuantumSystem` via the `global_params`
# keyword argument. With a custom integrator from Piccolissimo, the Hamiltonian
# function receives `u = [controls..., globals...]`:
#
# ```julia
# ## Function-based system (requires Piccolissimo integrator for dynamics)
# H = (u, t) -> u[3] * PAULIS[:Z] + u[1] * PAULIS[:X] + u[2] * PAULIS[:Y]
# sys = QuantumSystem(H, [1.0, 1.0]; time_dependent=true, global_params=(δ=0.5,))
# ```
#
# For this guide we use a matrix-based system, which works with the built-in
# `BilinearIntegrator`:

using Piccolo

sys = QuantumSystem(
    PAULIS[:Z],
    [PAULIS[:X], PAULIS[:Y]],
    [1.0, 1.0];
    global_params = (δ = 0.5,),
)

sys.global_params

# ## Setting Up a Problem with Globals
#
# When `global_params` is set on the system, `SmoothPulseProblem` automatically
# includes the global variables in the trajectory as optimization variables.
# The `global_bounds` keyword constrains their range:

T = 10.0
N = 50
U_goal = GATES[:X]

qtraj = UnitaryTrajectory(sys, U_goal, T)

qcp = SmoothPulseProblem(
    qtraj,
    N;
    Q = 100.0,
    R = 1e-2,
    global_bounds = Dict(:δ => 1.0),  ## symmetric bounds: δ ∈ [δ₀ - 1.0, δ₀ + 1.0]
)

cached_solve!(
    qcp,
    "global_variables_single";
    max_iter = 100,
    verbose = false,
    print_level = 1,
)

# ## Accessing Global Variables
#
# After solving, global values are stored in the trajectory's `global_data`
# vector, indexed by `global_components`:

traj = get_trajectory(qcp)
optimized_δ = traj.global_data[traj.global_components[:δ]][1]
optimized_δ

# !!! note
#     With the built-in `BilinearIntegrator`, the global variable is not coupled
#     to the dynamics, so its value will remain near the initial value. To
#     actually optimize globals through the Hamiltonian, use a Piccolissimo
#     integrator that provides global-aware Jacobians and Hessians.

# ## Global Bounds Format
#
# Global bounds can be specified in two ways:
#
# **Symmetric bounds** (scalar): applied as ±value around the initial value.
# ```julia
# global_bounds = Dict(:δ => 0.05)  # δ ∈ [δ₀ - 0.05, δ₀ + 0.05]
# ```
#
# **Asymmetric bounds** (tuple): explicit (lower, upper) bounds.
# ```julia
# global_bounds = Dict(:δ => (0.05, 0.2))  # δ ∈ [0.05, 0.2]
# ```
#
# **Multiple globals** with mixed bound types:
# ```julia
# global_bounds = Dict(
#     :δ => (0.05, 0.2),
#     :J => 0.01,  # J ∈ [J₀ - 0.01, J₀ + 0.01]
# )
# ```

# ## Multiple Global Variables
#
# You can define several system parameters simultaneously:

sys_multi =
    QuantumSystem(PAULIS[:Z], [PAULIS[:X]], [1.0]; global_params = (ω = 1.0, J = 0.05))

qtraj_multi = UnitaryTrajectory(sys_multi, U_goal, T)

qcp_multi = SmoothPulseProblem(
    qtraj_multi,
    N;
    Q = 100.0,
    R = 1e-2,
    global_bounds = Dict(:ω => (0.9, 1.1), :J => (0.02, 0.1)),
)

cached_solve!(
    qcp_multi,
    "global_variables_multi";
    max_iter = 100,
    verbose = false,
    print_level = 1,
)

## Access optimized global values
traj_multi = get_trajectory(qcp_multi)
optimized_global_data = traj_multi.global_data
(optimized_global_data)

# ## Best Practices

# ### 1. Start with Good Initial Values
#
# Use physically reasonable initial parameters. The initial values come from
# `global_params` in the `QuantumSystem` constructor:
#
# ```julia
# sys = QuantumSystem(H_drift, H_drives, bounds; global_params = (δ = 0.15,))
# ```

# ### 2. Use Reasonable Bounds
#
# Don't make bounds too wide — ~50% variation is usually sufficient:
#
# ```julia
# global_bounds = Dict(:δ => (0.1, 0.2))  # Not 10x variation
# ```

# ### 3. Consider Sensitivity
#
# Global parameters affect all timesteps, so they have strong influence on the
# solution. Start with tighter bounds and relax if needed.

# ### 4. Verify Physical Meaning
#
# After optimization, check that global values are physically reasonable:

# ## Comparison with SamplingProblem
#
# | Feature | Global Variables | SamplingProblem |
# |---------|-----------------|-----------------|
# | Parameter role | Optimized | Fixed (sampled) |
# | Goal | Find optimal parameter | Robust to parameter variation |
# | Use case | Calibration, design | Uncertainty handling |
#
# For robustness to uncertainty, use `SamplingProblem`. For finding optimal
# parameters, use global variables.

# ## Limitations
#
# - Full global optimization requires a custom integrator from Piccolissimo that
#   provides global-aware Jacobians and Hessians
# - Piccolo's built-in `BilinearIntegrator` stores globals in the trajectory but
#   does not couple them to the dynamics
# - More complex optimization landscape
# - Convergence can be slower

# ## See Also
#
# - [SamplingProblem](@ref sampling) - For parameter robustness
# - [Quantum Systems](@ref quantum-systems) - System construction
# - [Problem Templates](@ref problem-templates-overview) - Main optimization API
