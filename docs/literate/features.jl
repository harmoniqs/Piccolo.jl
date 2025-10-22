# ```@meta
# CollapsedDocStrings = true
# ```

# Feature Listing

# ## Quantum Systems

#=
General overview of the Piccolo.jl workflow:
1. Define a Hamiltonian (drift and drive)
2. Create a `QuantumSystem` from the Hamiltonian
3. Define a target unitary (or state)
4. Do one of the following:
    a. Build an optimization problem using initial trajectory, Objectives, Constraints+Dynamics, Initialization
    b. Use a problem template to quickly set up a problem

5. Finally, solve and get trajectory/pulse
=#

# ### Some Basic Systems
# Building a hamiltonian
using Piccolo
using SparseArrays

# We can build a hamiltonian from operators, e.g. the Pauli operators, creation/annihilation operators, etc.
#=
```@docs; canonical = false
PAULIS
create
annihilate
```
=#
PAULIS.X

# Number operator (on 3 level qubit) is just:
number_op =  create(3) * annihilate(3)
number_op

# Build single qubit system from paper
# ## TODO

# Build two qubit system from paper
levels = [2, 2]

a = sparse(lift_operator(create(levels[1]), 1, levels))
b = sparse(lift_operator(create(levels[2]), 2, levels))

g = 0.1 * 2 * π # 100 MHz

function two_qubit_2_2_level_system()
    H_drift = g * (a'a) * (b'b)
    H_drives = [
        a + a',
        im * (a - a'),
        b + b',
        im * (b - b')
    ]
    return [H_drift, H_drives]
end

H_drift, H_drives = two_qubit_2_2_level_system()
system = QuantumSystem(H_drift, H_drives)

# Build multilevel_transmon example from QC docs
levels = 5

a = sparse(create(levels)) # no lift operator required here as there is only one qubit

ω = 4.0 # 4 GHz
δ = 0.2 # 0.2 GHz

function single_qubit_5_level_transmon()
    H_drift = -δ / 2 * a' * a' * a * a  # rotating frame
    H_drives = [
        a + a',
        im * (a - a')
    ]
    return [H_drift, H_drives]
end

H_drift, H_drives = single_qubit_5_level_transmon()
system = QuantumSystem(H_drift, H_drives)

# or you could use one of our built-in systems like [`TransmonSystem`](@extref QuantumCollocation.QuantumSystemTemplates.TransmonSystem)

#=
```@docs; canonical = false
TransmonSystem
```
=#

system = TransmonSystem(levels=levels, δ=δ)
system

# Build three qubit system from paper
# Build Cat system (TODO: (aaron) write about this)

# ### Time-dependent Systems
# Example of defining hamiltonian as lambda,
# shortcuts for carrier waves and offets -> special constructor
# for time-dependent systems

# ### Open Systems
# From PQO can rip example

# ### Composite systems
# From PQO can rip example - fix up wording? (TODO (andy) review wording)

# ### Variational Quantum Systems
# TODO (gennadi)

# ## Basic Quantum System operations
# ### Internal representations: isomorphisms (state, operator, density matrices, dynamics, etc.
# #### Isomorphisms
# Explain its effect on the interface for users
# ** DETAILS: for those who want to know more, we can link to the PQO docs
# #### Gates
# * PQO Gates consts
# #### Operators
# * Annihilation, creation, number
# * Ladder operators
# * Displacement, squeezing, etc.
# * Pauli operators

# ### Drift vs. drive

# ### Rollouts & fidelities
# ### Reachability tests
# ### Direct sums
# ### GRAPE

# ## Trajectory Initialization
#=
```@docs; canonical = false
initialize_trajectory(
        U_goal::AbstractPiccoloOperator,
        T::Int,
        Δt::Union{Real, AbstractVecOrMat{<:Real}}
    )

QuantumCollocation.TrajectoryInitialization.unitary_geodesic
QuantumCollocation.TrajectoryInitialization.unitary_linear_interpolation
```
=#

# more on interpolations [here](@extref NamedTrajectories.trajectory_interpolation)

# ## Hamiltonian construction + limitations

# ## Optimization Parameters
# ### Custom Objectives
# ### Custom Constraints

# ## Optimization solvers
# ### TODO: IPOPT vs. MadNLP
# ### TODO: Custom callbacks, trajectory checkpointing, etc.

# ## Problem Templates usage
# ### Min time
# ### smooth
# ### …
# ### State vs Unitary problems

# ## Global Parameters
# ## Callbacks
# ## Rollouts

# ## Fidelities
# ## Subspaces/leakage
# ## Integrators
# ## Plotting
