# # Multilevel Transmon

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
times = collect(range(0, T₀, length=N))

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
qcp = SmoothPulseProblem(qtraj, N; ddu_bound=ddu_bound, Q=100.0, R=1e-2)

#=
To solve this problem, we would run:

```julia
solve!(qcp; max_iter=50)
```

The output would look something like:

```
    initializing optimizer...
    applying constraint: timesteps all equal constraint
    applying constraint: initial value of Ũ⃗
    applying constraint: initial value of u
    applying constraint: final value of u
    applying constraint: bounds on u
    applying constraint: bounds on du
    applying constraint: bounds on ddu
    applying constraint: bounds on Δt
This is Ipopt version 3.14.19, running with linear solver MUMPS 5.8.1.

Number of nonzeros in equality constraint Jacobian...:   130578
Number of nonzeros in inequality constraint Jacobian.:        0
Number of nonzeros in Lagrangian Hessian.............:    11223

Total number of variables............................:     2796

Number of Iterations....: 50

EXIT: Maximum Number of Iterations Exceeded.
```
=#

# For this documentation, we skip the actual solve and show pre-computed results below.
# In practice, you would run `solve!(qcp; max_iter=50)` to optimize the pulse.

# After optimization, you can check the fidelity in the subspace:

#=
```julia
fid = fidelity(qcp)
println("Fidelity: ", fid)
```

With sufficient iterations, you should see fidelity > 0.99.
=#

# ## Leakage suppression

# As can be seen from multilevel systems, there can be substantial leakage into higher energy levels
# during the evolution. To mitigate this, Piccolo provides options to add leakage suppression
# via the `PiccoloOptions` type.

# To implement leakage suppression, pass `leakage_constraint=true` and configure the leakage
# parameters:

#=
```julia
## create a leakage suppression problem
qcp_leakage = SmoothPulseProblem(
    qtraj, N;
    ddu_bound=ddu_bound,
    Q=100.0,
    R=1e-2,
    piccolo_options=PiccoloOptions(
        leakage_constraint=true,
        leakage_constraint_value=1e-2,
        leakage_cost=1e-2,
    ),
)

## solve the problem
solve!(qcp_leakage; max_iter=50)
```
=#

# The leakage suppression adds:
# - An L1-norm cost on populating leakage levels (drives populations toward zero)
# - A constraint that keeps leakage below the specified threshold

# After optimization with leakage suppression, you should see substantially reduced
# population in the higher energy levels during the gate evolution.

# ## Visualizing Results

# Piccolo provides plotting utilities to visualize the unitary evolution:

#=
```julia
## plot the state populations over time
plot_unitary_populations(get_trajectory(qcp); fig_size=(900, 700))
```
=#

# This shows how the population flows between energy levels during the gate.
# With leakage suppression, you should see the population staying mostly
# within the computational subspace (first two levels).

# ## Summary

# This tutorial demonstrated:
# 1. Using `TransmonSystem` to create a realistic multilevel transmon system
# 2. Using `EmbeddedOperator` to define gates on a subspace
# 3. Creating `UnitaryTrajectory` with `ZeroOrderPulse`
# 4. Using `SmoothPulseProblem` to set up the optimization
# 5. Adding leakage suppression via `PiccoloOptions`

# For more details on:
# - Problem templates: See [SmoothPulseProblem](@ref)
# - Leakage suppression: See [Leakage Suppression Guide](@ref)
# - System templates: See [System Templates](@ref)
