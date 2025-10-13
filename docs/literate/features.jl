# Feature Listing

# ## Quantum Systems

#=
General overview of the workflow to use this package:
1. Define a Hamiltonian (drift and drive)
2. Create a `QuantumSystem` from the Hamiltonian
3. Define a target unitary (or state)
4. Do one of the following:

    a. Build an optimization problem using initial trajectory, Objectives, Constraints+Dynamics, Initialization
    
    b. Use a problem template to quickly set up a problem

5. Solve and get trajectory
=#

# ### Some Basic Systems
# Building a hamiltonian
using Piccolo

# We can build a hamiltonian from operators, e.g. the Pauli operators, creation/annihilation operators, etc.
# They are defined in the `PiccoloQuantumObjects.Gates.PAULIS` NamedTuple and `PiccoloQuantumObjects.QuantumObjectUtils` module.
#md # ```@docs
#md # Piccolo.Gates.PAULIS
#md # Piccolo.QuantumObjectUtils
#md # ```

# Annihilation, creation, number operators
# [Link to operators]

# Build single qubit system from paper
# Build two qubit system from paper
# Build multilevel_transmon example from QC docs
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
