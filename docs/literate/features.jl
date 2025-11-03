# ```@meta
# CollapsedDocStrings = true
# ```

# # Features
# Please use the sidebar and the table of contents below to navigate this long single page
# This page's goal is to provide a repository for commonly used features and usage patterns
# that users frequently use.
#=
```@contents
Pages = ["features.md"]
Depth = 2:3
```
=#

# ## General Overview
#=
General overview of the Piccolo.jl workflow:
1. Define a hamiltonian (drift and drive)
2. Create a `QuantumSystem` from the hamiltonian
3. Define a target unitary (or state)
4. Do one of the following:
    a. Build an optimization problem using initial trajectory, Objectives, Constraints+Dynamics, Initialization
    b. Use a problem template to quickly set up a problem

5. Finally, solve and get trajectory/pulse
=#

# ## Quantum Systems
# ### Systems Basics
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

# Number operator (on a 3-level qutrit) is just:
number_op =  create(3) * annihilate(3)
number_op

# Pauli X, Y, Z, I are also defined there
PAULIS.Y

# Operators often act on certain subsystems, specific qubits, or on specific subspaces in the system.
# Placing them correctly can be done via the [`lift_operator`](@extref PiccoloQuantumObjects.QuantumSystems.lift_operator) function

# Example: lifting the Pauli X to act on the second qubit of a 2, 2 level quantum system. Equivalent to the below
#=
```math
I \otimes X = \begin{bmatrix} 0 & 1 & 0 & 0 \\ 1 & 0 & 0 & 0 \\ 0 & 0 & 0 & 1 \\ 0 & 0 & 1 & 0 \end{bmatrix}
```
=#
levels = [2, 2]
sparse(lift_operator(PAULIS.X, 2, levels))

# For systems with more levels, its as simple as defining the number of levels for each subsystem, and then placing the operator on
# the right subsystem.
levels = [2, 3]
sparse(lift_operator(annihilate(levels[2]), 2, levels))

# ### Example multi-qudit system
# Here are 2 qudits with 5 levels.
# We use [`lift_operator`](@extref PiccoloQuantumObjects.QuantumSystems.lift_operator) to place the operator for each respective subsystem into a sparse matrix.
levels = [5, 5]

a = sparse(lift_operator(annihilate(levels[1]), 1, levels))
b = sparse(lift_operator(annihilate(levels[2]), 2, levels))

g = 0.1 * 2 * π # 100 MHz

function my_system()
    H_drift = g * (a'a) * (b'b)
    H_drives = [
        a + a',
        im * (a - a'),
        b + b',
        im * (b - b')
    ]
    return [H_drift, H_drives]
end

H_drift, H_drives = my_system()
system = QuantumSystem(H_drift, H_drives)

# ### Example multi-level transmon system
# As in the [multilevel_transmon example from the QC docs](@extref generated/examples/multilevel_transmon)
levels = 5

a = sparse(annihilate(levels)) # no lift operator required here as there is only one qudit

ω = 4.0 # 4 GHz
δ = 0.2 # 0.2 GHz

function single_5_level_transmon()
    H_drift = -δ / 2 * a' * a' * a * a  # rotating frame
    H_drives = [
        a + a',
        im * (a - a')
    ]
    return [H_drift, H_drives]
end

H_drift, H_drives = single_5_level_transmon()
system = QuantumSystem(H_drift, H_drives)

# or you could use one of our built-in systems like [`TransmonSystem`](@extref QuantumCollocation.QuantumSystemTemplates.TransmonSystem)

#=
```@docs; canonical = false
TransmonSystem
```
=#
levels = 5 # hide
system = TransmonSystem(levels=levels, δ=δ)

# ### (TODO) Example: Build three qubit system from paper

# ### (TODO) Example: Build Cat system (TODO: (aaron) write about this)

# ### Example: Time-dependent Systems
# You can define hamiltonians that are dependent on time as well and pass them into the [`TimeDependentQuantumSystem`](@extref PiccoloQuantumObjects.QuantumSystems.TimeDependentQuantumSystem) constructor.
# For example, this hamiltonian:
#=
```math
H(a, t) = Z + a_1 \cos \left( t \right)\, X
```
=#
H(a, t) = PAULIS.Z + a[1] * cos(t) * PAULIS.X
n_drives = 1  # our `a` vector has one element
TimeDependentQuantumSystem(H, n_drives)

# See more on PiccoloQuantumObjects: [Time-dependent quantum systems examples](@extref Time-Dependent-Quantum-Systems)

# ### (TODO: Gennadi) shortcuts for carrier waves and offets -> special constructor
# for time-dependent systems

# ### Open Systems
# You can define hamiltonians and dissapation operators for open quantum systems defined via Lindblad dynamics
#=
```math
\rho
```
=#

H_drives = [PAULIS[:X]]
a = annihilate(2)
dissipation_operators = [a'a, a]
system = OpenQuantumSystem(H_drives, dissipation_operators=dissipation_operators)

# See more on PiccoloQuantumObjects: [Open Quantum Systems](@extref Open-quantum-systems)

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

# more on interpolations [here](@extref NamedTrajectories.MethodsNamedTrajectory.trajectory_interpolation)

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
