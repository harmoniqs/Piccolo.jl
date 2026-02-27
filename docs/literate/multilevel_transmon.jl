# # [Multilevel Transmon](@id multilevel-transmon-tutorial)

# In this example we will look at a multilevel transmon qubit with a Hamiltonian given by
#
# ```math
# \hat{H}(t) = -\frac{\delta}{2} \hat{n}(\hat{n} - 1) + u_1(t) (\hat{a} + \hat{a}^\dagger) + u_2(t) i (\hat{a} - \hat{a}^\dagger)
# ```
# where $\hat{n} = \hat{a}^\dagger \hat{a}$ is the number operator, $\hat{a}$ is the annihilation operator, $\delta$ is the anharmonicity, and $u_1(t)$ and $u_2(t)$ are control fields.
#
# We will use the following parameter values:
#
# ```math
# \begin{aligned}
# \delta &= 0.2 \text{ GHz}\\
# \abs{u_i(t)} &\leq 0.2 \text{ GHz}\\
# T_0 &= 10 \text{ ns}\\
# \end{aligned}
# ```
#
# For convenience, we have defined the `TransmonSystem` function in the `QuantumSystemTemplates` module, which returns a `QuantumSystem` object for a transmon qubit. We will use this function to define the system.

# ## Setting up the problem

# To begin, let's load the necessary packages, define the system parameters, and create a `QuantumSystem` object using the `TransmonSystem` function.

using Piccolo
using SparseArrays
using Random;
Random.seed!(123);

using CairoMakie

## define the time parameters

T₀ = 10     # total time in ns
N = 50      # number of time steps
Δt = T₀ / N # time step
times = collect(range(0, T₀, length = N))

## define the system parameters
levels = 5
δ = 0.2

## add a bound to the controls
u_bound = [0.2, 0.2]
ddu_bound = 1.0

## create the system
sys = TransmonSystem(drive_bounds = u_bound, levels = levels, δ = δ)

## let's look at the drives of the system
get_drives(sys)[1] |> sparse


# Since this is a multilevel transmon and we want to implement an, let's say, $X$ gate on the qubit subspace, i.e., the first two levels we can utilize the `EmbeddedOperator` type to define the target operator.

## define the target operator
op = EmbeddedOperator(:X, sys)

## show the full operator
op.operator |> sparse

# We can then create a pulse, trajectory, and optimization problem using the new API:

## create a random initial pulse
initial_controls = 0.1 * randn(2, N)
pulse = ZeroOrderPulse(initial_controls, times)

## create a unitary trajectory with the embedded operator as goal
qtraj = UnitaryTrajectory(sys, pulse, op)

## create the optimization problem
qcp = SmoothPulseProblem(qtraj, N; ddu_bound = ddu_bound, Q = 100.0, R = 1e-2)

# ## Solving the problem

## We solve the problem using `cached_solve!`, which transparently caches the
## optimized trajectory and solver output for docs purposes. In practice, you can use `solve!` directly.

cached_solve!(qcp, "multilevel_transmon"; max_iter = 50)

# After optimization, we can check the fidelity in the subspace:

fid = fidelity(qcp)
println("Fidelity: ", fid)

# ## Leakage suppression

# As can be seen from multilevel systems, there can be substantial leakage into higher energy levels
# during the evolution. To mitigate this, Piccolo provides options to add leakage suppression
# via the `PiccoloOptions` type.

# To implement leakage suppression, pass `leakage_constraint=true` and configure the leakage
# parameters:

## create a leakage suppression problem
qcp_leakage = SmoothPulseProblem(
    qtraj,
    N;
    ddu_bound = ddu_bound,
    Q = 100.0,
    R = 1e-2,
    piccolo_options = PiccoloOptions(
        leakage_constraint = true,
        leakage_constraint_value = 1e-2,
        leakage_cost = 1e-2,
    ),
)

## solve the problem
cached_solve!(qcp_leakage, "multilevel_transmon_leakage"; max_iter = 50)

# The leakage suppression adds:
# - An L1-norm cost on populating leakage levels (drives populations toward zero)
# - A constraint that keeps leakage below the specified threshold

fid_leakage = fidelity(qcp_leakage)
println("Fidelity with leakage suppression: ", fid_leakage)

# ## Visualizing Results

# Piccolo provides plotting utilities to visualize the unitary evolution.
# First, without leakage suppression:

plot_unitary_populations(get_trajectory(qcp); fig_size = (900, 700))

# And with leakage suppression — you should see the population staying mostly
# within the computational subspace (first two levels):

plot_unitary_populations(get_trajectory(qcp_leakage); fig_size = (900, 700))

# ## Summary

# This tutorial demonstrated:
# 1. Using `TransmonSystem` to create a realistic multilevel transmon system
# 2. Using `EmbeddedOperator` to define gates on a subspace
# 3. Creating `UnitaryTrajectory` with `ZeroOrderPulse`
# 4. Using `SmoothPulseProblem` to set up the optimization
# 5. Adding leakage suppression via `PiccoloOptions`

# For more details on:
# - Problem templates: See [SmoothPulseProblem](@ref smooth-pulse)
# - Leakage suppression: See [Leakage Suppression](@ref leakage-suppression)
# - Quantum systems: See [Transmon Qubits](@ref transmon-systems)
