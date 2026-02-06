# # [Quantum Module API](@id quantum-api)
#
# The Quantum module provides all quantum mechanical building blocks: systems, trajectories, pulses, operators, and isomorphisms.
# This guide demonstrates the API with runnable examples.

# ## Setup
#
# First, load the required packages:

using Piccolo
using CairoMakie
using SparseArrays # hide

# # Quantum Systems
#
# ## Basic QuantumSystem
#
# A `QuantumSystem` represents a quantum system with drift and drive Hamiltonians.
# Let's create a simple qubit system:

## Define the Hamiltonian components
H_drift = 0.5 * PAULIS[:Z]  # Drift: ω/2 σ_z
H_drives = [PAULIS[:X], PAULIS[:Y]]  # Drives: σ_x and σ_y
drive_bounds = [1.0, 1.0]  # Maximum amplitudes

## Create the system
sys = QuantumSystem(H_drift, H_drives, drive_bounds)

println("System has $(sys.levels) energy levels")
println("System has $(sys.n_drives) control drives")

# ## Accessing System Properties
#
# You can extract the components of a quantum system:

H_d = get_drift(sys)
H_c = get_drives(sys) # nothing

# The drift Hamiltonian:
H_d

# The drive Hamiltonians (first drive shown):
H_c[1]

# ## System Templates
#
# Piccolo provides pre-built constructors for common hardware platforms.
#
# ### TransmonSystem
#
# Create a superconducting transmon qubit with anharmonicity:

transmon_sys = TransmonSystem(
    levels = 3,           # Include first 3 levels
    δ = -0.2,             # Anharmonicity (relative to ω)
    drive_bounds = [0.2, 0.2]  # Control amplitudes
)

# This creates a 3-level transmon system:
transmon_sys.levels

# # Pulses
#
# ## ZeroOrderPulse
#
# Piecewise constant controls (use with `SmoothPulseProblem`):

T = 10.0
N = 50
times = collect(range(0, T, length = N))
controls = 0.1 * randn(2, N)

pulse = ZeroOrderPulse(controls, times)

# Check the pulse properties:

println("Pulse duration: ", duration(pulse))
println("Number of drives: ", n_drives(pulse))
println("Control at t=5: ", pulse(5.0))

# ## GaussianPulse
#
# Analytic Gaussian envelope for smooth pulses:

gaussian_pulse = GaussianPulse(
    [0.5, 0.3],  # Amplitudes for each drive
    2.0,         # Sigma (width)
    10.0         # Duration
)

duration(gaussian_pulse)

# Let's visualize the Gaussian pulse:

fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1], xlabel = "Time", ylabel = "Amplitude", title = "Gaussian Pulse")

t_plot = range(0, duration(gaussian_pulse), length = 200)
u_plot = hcat([gaussian_pulse(t) for t in t_plot]...)'

lines!(ax, t_plot, u_plot[:, 1], label = "Drive 1", linewidth = 2)
lines!(ax, t_plot, u_plot[:, 2], label = "Drive 2", linewidth = 2)
axislegend(ax, position = :rt)

fig

# # Quantum Trajectories
#
# ## UnitaryTrajectory
#
# For unitary gate synthesis, we specify a target gate:

U_goal = GATES[:X]  # Target: X gate
qtraj = UnitaryTrajectory(sys, pulse, U_goal)

# The goal operator:
U_goal |> sparse

# ## KetTrajectory
#
# For state-to-state transfer:

ψ_init = [1.0 + 0im, 0.0 + 0im]  # |0⟩
ψ_goal = [0.0 + 0im, 1.0 + 0im]  # |1⟩

ket_traj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

println("Initial state: ")
println(sparse(ψ_init))
println("Goal state: ")
println(sparse(ψ_goal))

# ## MultiKetTrajectory
#
# For implementing gates via multiple state mappings:

## Define computational basis states
ψ0 = [1.0 + 0im, 0.0 + 0im]
ψ1 = [0.0 + 0im, 1.0 + 0im]

## X gate should map |0⟩ → |1⟩ and |1⟩ → |0⟩
initial_states = [ψ0, ψ1]
goal_states = [ψ1, ψ0]

multi_traj = MultiKetTrajectory(sys, pulse, initial_states, goal_states)

# # Operators
#
# ## Pauli Matrices
#
# Access standard Pauli operators. The X matrix:

PAULIS[:X] |> sparse

# The Y matrix:
PAULIS[:Y] |> sparse

# The Z matrix:
PAULIS[:Z] |> sparse

# ## Common Gates
#
# The `GATES` dictionary provides standard quantum gates:

keys(GATES)

# The Hadamard gate:
GATES[:H] |> sparse

# The CNOT gate:
GATES[:CX] |> sparse

# ## Bosonic Operators
#
# For oscillator/cavity systems, use creation and annihilation operators:

n_levels = 4
a = annihilate(n_levels)
a_dag = create(n_levels)

# The annihilation operator (4 levels):
a |> sparse

# Verify the commutator [a, a†] = I:
a * a_dag - a_dag * a |> sparse

# ## EmbeddedOperator
#
# For gates on a subspace of a larger Hilbert space.
# This is useful when you have a transmon with leakage levels
# but want to target a gate on the qubit subspace:

## Create a 3-level transmon
transmon_3level = TransmonSystem(levels = 3, δ = -0.2, drive_bounds = [0.2, 0.2])

## X gate embedded in the qubit subspace (levels 0,1)
embedded_X = EmbeddedOperator(:X, transmon_3level)

## Get the full-space operator by embedding the X gate into the 3-level system
X_2level = GATES[:X]
X_embedded_full = embed(X_2level, embedded_X) 
X_embedded_full |> sparse

# # Isomorphisms
#
# For optimization, we convert complex quantum objects to real vectors.
#
# ## Ket States

ψ = [1.0 + 0im, 0.0 + 0im]  # Complex ket
ψ̃ = ket_to_iso(ψ)           # Real isomorphic vector

# The isomorphic representation:
(ψ, ψ̃)

# Convert back to complex ket:
iso_to_ket(ψ̃) 

# ## Operators/Unitaries

U = GATES[:H]                # Complex operator
Ũ⃗ = operator_to_iso_vec(U)   # Real vector

# The Hadamard operator:
U

# Convert back from the isomorphic vector:
iso_vec_to_operator(Ũ⃗)

# # Trajectory Functions
#
# Once you have a trajectory, you can access various properties:

println("Initial fidelity: ", fidelity(qtraj))
println("Pulse duration: ", duration(get_pulse(qtraj)))
println("System levels: ", get_system(qtraj).levels)
println("State variable name: ", state_name(qtraj))
println("Control variable name: ", drive_name(qtraj))

# # Summary
#
# The Quantum module provides:
#
# - **Systems**: `QuantumSystem`, `OpenQuantumSystem`, `CompositeQuantumSystem`, `VariationalQuantumSystem`
# - **System Templates**: `TransmonSystem`, `MultiTransmonSystem`, `IonChainSystem`, `RydbergChainSystem`, `CatSystem`
# - **Trajectories**: `UnitaryTrajectory`, `KetTrajectory`, `MultiKetTrajectory`, `DensityTrajectory`, `SamplingTrajectory`
# - **Pulses**: `ZeroOrderPulse`, `LinearSplinePulse`, `CubicSplinePulse`, `GaussianPulse`, `CompositePulse`
# - **Operators**: `PAULIS`, `GATES`, `annihilate`, `create`, `EmbeddedOperator`
# - **Isomorphisms**: `ket_to_iso`, `operator_to_iso_vec`, `density_to_iso_vec` (and their inverses)
#
# For complete function signatures and advanced usage, see the [API Reference](@ref) overview.
