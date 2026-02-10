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
# 1. A **time-dependent system** with parameters
# 2. A **custom integrator** that supports globals (e.g., from Piccolissimo)
# 3. **Global bounds** specification

# ## Example: Defining a Parametric System
#
# The first step is to create a time-dependent `QuantumSystem` with named global parameters:

using Piccolo

## Hamiltonian with global parameter δ (detuning)
## The global parameter is passed as an extra element of u
function H(u, t)
    δ = u[end]
    return δ * PAULIS[:Z] + u[1] * PAULIS[:X] + u[2] * PAULIS[:Y]
end

## Create system with named global parameter and initial value
sys = QuantumSystem(H, [1.0, 1.0]; time_dependent = true, global_params = (δ = 0.1,))

sys.global_params

# ## Setting Up a Problem with Globals
#
# Global variable optimization requires a custom integrator that handles global
# parameters. The `HermitianExponentialIntegrator` from Piccolissimo supports this:

using Piccolissimo

T, N = 10.0, 100
times = collect(range(0, T, length = N))
pulse = ZeroOrderPulse(0.1 * randn(2, N), times)
qtraj = UnitaryTrajectory(sys, pulse, GATES[:X])

## Create integrator with global support
integrator = HermitianExponentialIntegrator(qtraj, N; global_names = [:δ])

## Set up problem with global bounds
qcp = SmoothPulseProblem(
    qtraj,
    N;
    integrator = integrator,
    global_names = [:δ],
    global_bounds = Dict(:δ => (0.05, 0.2)),  # δ ∈ [0.05, 0.2]
)

solve!(qcp; max_iter = 100)

## Access optimized global value
traj = get_trajectory(qcp)
optimized_δ = traj[:δ][1]  # Global variables are constant over time
optimized_δ

# ## Global Bounds Format
#
# Global bounds can be specified in two ways:

## Symmetric bounds: δ ∈ [δ₀ - 0.05, δ₀ + 0.05]
global_bounds = Dict(:δ => 0.05)

## Asymmetric bounds: δ ∈ [0.05, 0.2]
global_bounds = Dict(:δ => (0.05, 0.2))

## Multiple globals
global_bounds = Dict(
    :δ => (0.05, 0.2),
    :J => 0.01,  # J ∈ [J₀ - 0.01, J₀ + 0.01]
)

# ## Multiple Global Variables
#
# You can optimize several system parameters simultaneously:

## Hamiltonian with multiple parameters
function H_multi(u, t)
    ω = u[end-1]
    J = u[end]
    return ω * PAULIS[:Z] + J * PAULIS[:X] + u[1] * PAULIS[:X]
end

sys_multi = QuantumSystem(
    H_multi,
    [1.0];
    time_dependent = true,
    global_params = (ω = 1.0, J = 0.05),
)

## Optimize ω and J
qcp = SmoothPulseProblem(
    qtraj,
    N;
    integrator = integrator,
    global_names = [:ω, :J],
    global_bounds = Dict(:ω => (0.9, 1.1), :J => 0.02),
)

# ## Best Practices

# ### 1. Start with Good Initial Values
#
# Use physically reasonable initial parameters:

sys = QuantumSystem(H, [1.0, 1.0]; time_dependent = true, global_params = (δ = 0.15,))

# ### 2. Use Reasonable Bounds
#
# Don't make bounds too wide — ~50% variation is usually sufficient:

global_bounds = Dict(:δ => (0.1, 0.2))  # Not 10x variation

# ### 3. Consider Sensitivity
#
# Global parameters affect all timesteps, so they have strong influence on the
# solution. Start with tighter bounds and relax if needed.

# ### 4. Verify Physical Meaning
#
# After optimization, check that global values are physically reasonable:

optimized_δ = traj[:δ][1]
@assert 0.05 < optimized_δ < 0.25 "Detuning outside expected range"

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
# - Requires custom integrator support (e.g., Piccolissimo)
# - More complex optimization landscape
# - Convergence can be slower

# ## See Also
#
# - [SamplingProblem](@ref sampling) - For parameter robustness
# - [Quantum Systems](@ref quantum-systems) - System construction
# - [Problem Templates](@ref problem-templates-overview) - Main optimization API
