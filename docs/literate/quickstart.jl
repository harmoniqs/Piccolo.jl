pretty_print(X::AbstractMatrix) = Base.show(stdout, "text/plain", X); # Helper function

```
* compute a rollout with zero controls
* compute a rollout with random controls
* compute a fidelity show its bad
* compute good controls
    * trajectory initialization
    * objectives
    * dynamics (technically constraints)
    * constraints
    * note: integrators, [link to integrators]
    * run solve!
    * get controls out
* compute a fidelity show its good
** note: I use grape for this! [link to more on grape vs Piccolo section]

Build single qubit system from paper
* skip intermediate steps and use problem template
* run min-time
```

#=

# Quickstart Guide

To set up and solve a quantum optimal control problems we provide high level problem templates to quickly get started. For unitary gate problems, where we want to realize a gate $U_{\text{goal}}$, with a system Hamiltonian of the form,
```math
H(t) = H_0 + \sum_i a^i(t) H_i
```
there is the `UnitarySmoothPulseProblem` constructor which only requires
- the drift Hamiltonian, $H_0$
- the drive Hamiltonians, $\qty{H_i}$
- the target unitary, $U_{\text{goal}}$
- the number of timesteps, $T$
- the (initial) time step size, $\Delta t$

## Smooth Pulse Problems

For example, to create a problem for a single qubit $X$ gate (with a bound on the drive of $|a^i| < a_{\text{bound}}$), i.e., with system hamiltonian
```math
H(t) = \frac{\omega}{2} \sigma_z + a^1(t) \sigma_x + a^2(t) \sigma_y
```
we can do the following:

=#

using CairoMakie
using Piccolo

## set time parameters
T = 50
Δt = 0.2

## define drift and drive Hamiltonians
H_drift = 0.2 * PAULIS.Z
H_drives = [PAULIS.X, PAULIS.Y]

## create a QuantumSystem from the Hamiltonians
system = QuantumSystem(H_drift, H_drives)

## define target unitary
U_goal = GATES.X

## set bounds on the drive
a_bound = 0.2
dda_bound = 5.0

## build the problem
prob = UnitarySmoothPulseProblem(
    system,
    U_goal,
    T,
    Δt;
    a_bound = a_bound,
    dda_bound = dda_bound,
)

## solve the problem
solve!(prob; max_iter = 50)

#=
The above output comes from the Ipopt.jl solver. The problem's trajectory has been updated with the solution. 
We can see the final control amplitudes and the final unitary by accessing the `a` and `Ũ⃗` fields of the trajectory.
=#
size(prob.trajectory.a) |> println
size(prob.trajectory.Ũ⃗) |> println

#=
The `Ũ⃗` field is a vectorized representation of the unitary, which we can convert back to a 
matrix using the `iso_vec_to_operator` function exported by PiccoloQuantumObjects.jl.
=#
iso_vec_to_operator(prob.trajectory.Ũ⃗[:, end]) |> pretty_print

# To see the final fidelity we can use the `unitary_rollout_fidelity` function exported by QuantumCollocation.jl.
println("Final fidelity: ", unitary_rollout_fidelity(prob.trajectory, system))

# We can also easily plot the solutions using the `plot` function exported by NamedTrajectories.jl.
plot(prob.trajectory, [:Ũ⃗, :a])

# ## Minimum Time Problems

# We can also easily set up and solve a minimum time problem, where we enforce a constraint on the final fidelity:
# ```math
# \mathcal{F}(U_T, U_{\text{goal}}) \geq \mathcal{F}_{\text{min}}
# ```
# Using the problem we just solved we can do the following:

## final fidelity constraint
final_fidelity = 0.99

min_time_prob = UnitaryMinimumTimeProblem(prob, U_goal; final_fidelity = final_fidelity)

solve!(min_time_prob; max_iter = 50)

# We can see that the final fidelity is indeed greater than the minimum fidelity we set.

println("Final fidelity:    ", unitary_rollout_fidelity(min_time_prob.trajectory, system))

# and that the duration of the pulse has decreased.

initial_duration = get_times(prob.trajectory)[end]
min_time_duration = get_times(min_time_prob.trajectory)[end]

println("Initial duration:  ", initial_duration)
println("Minimum duration:  ", min_time_duration)
println("Duration decrease: ", initial_duration - min_time_duration)

## We can also plot the solutions for the minimum time problem, and see that the control amplitudes saturate the bound.
plot(min_time_prob.trajectory, [:Ũ⃗, :a])
