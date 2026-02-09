# # [Trajectories](@id trajectories-concept)
#
# Trajectories in Piccolo.jl represent a complete quantum control scenario: a
# system, a control pulse, and a goal. They are the central object passed to
# problem templates.
#
# ## Overview
#
# A trajectory encapsulates:
# - **Quantum system**: The hardware being controlled
# - **Pulse**: How controls vary over time
# - **Goal**: The desired outcome (gate, state, etc.)
# - **State evolution**: The quantum state over time
#
# ## Trajectory Types
#
# | Type | Use Case | State Variable |
# |------|----------|----------------|
# | `UnitaryTrajectory` | Gate synthesis | Unitary matrix `U(t)` |
# | `KetTrajectory` | State preparation | State vector `\|ψ(t)⟩` |
# | `DensityTrajectory` | Open systems | Density matrix `ρ(t)` |
# | `MultiKetTrajectory` | Multiple state transfers | Multiple `\|ψᵢ(t)⟩` |
# | `SamplingTrajectory` | Robust optimization | States for multiple systems |
#
# ## UnitaryTrajectory
#
# For synthesizing quantum gates.
#
# ### Construction

using Piccolo

## Define system
H_drift = PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
sys = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

## Create pulse
T, N = 10.0, 100
times = collect(range(0, T, length = N))
pulse = ZeroOrderPulse(0.1 * randn(2, N), times)

## Create trajectory with goal
U_goal = GATES[:X]
qtraj = UnitaryTrajectory(sys, pulse, U_goal)

# ### Solve and Analyze

qcp = SmoothPulseProblem(qtraj, N; Q = 100.0)
solve!(qcp; max_iter = 50)
fidelity(qcp)

# ### Extracting the Pulse

optimized_pulse = get_pulse(qcp.qtraj)
duration(optimized_pulse)

# ## KetTrajectory
#
# For state preparation (state-to-state transfer).
#
# ### Construction

## Initial and goal states
ψ_init = ComplexF64[1, 0]  # |0⟩
ψ_goal = ComplexF64[0, 1]  # |1⟩

qtraj_ket = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

# ### When to Use
#
# Use `KetTrajectory` when:
# - Preparing a specific quantum state
# - Global phase doesn't matter
# - Single state transfer

qcp_ket = SmoothPulseProblem(qtraj_ket, N; Q = 100.0)
solve!(qcp_ket; max_iter = 50)
fidelity(qcp_ket)

# ## MultiKetTrajectory
#
# For gates defined by multiple state mappings with coherent phases.
#
# ### Construction

## Define state pairs: X gate maps |0⟩ → |1⟩ and |1⟩ → |0⟩
ψ0 = ComplexF64[1, 0]
ψ1 = ComplexF64[0, 1]

initial_states = [ψ0, ψ1]
goal_states = [ψ1, ψ0]

qtraj_multi = MultiKetTrajectory(sys, pulse, initial_states, goal_states)

# ### Coherent Fidelity
#
# `MultiKetTrajectory` uses `CoherentKetInfidelityObjective`, which ensures:
# - Each state reaches its target
# - Relative phases between states are preserved
#
# This is important for gates where phase relationships matter (unlike
# `UnitaryTrajectory` which tracks the full unitary).

qcp_multi = SmoothPulseProblem(qtraj_multi, N; Q = 100.0)
solve!(qcp_multi; max_iter = 50)
fidelity(qcp_multi)

# ## DensityTrajectory
#
# For open quantum systems with dissipation.
#
# ```julia
# # Open system with collapse operators
# open_sys = OpenQuantumSystem(H_drift, H_drives, bounds, c_ops)
#
# # Initial density matrix
# ρ_init = [1.0 0.0; 0.0 0.0]  # Pure |0⟩
# ρ_goal = [0.0 0.0; 0.0 1.0]  # Pure |1⟩
#
# qtraj = DensityTrajectory(open_sys, pulse, ρ_init, ρ_goal)
# ```
#
# !!! note
#     `DensityTrajectory` support for fidelity objectives is still in
#     development. For most open-system problems, consider using `KetTrajectory`
#     with an effective Hamiltonian.
#
# ## SamplingTrajectory
#
# For robust optimization over parameter variations. This is created internally
# by `SamplingProblem`.
#
# ### How It Works
#
# `SamplingTrajectory` wraps multiple system variants:
#
# ```julia
# # Created internally by SamplingProblem
# systems = [sys_nominal, sys_high, sys_low]
# qcp_robust = SamplingProblem(qcp_base, systems)
# # This creates a SamplingTrajectory internally
# ```
#
# The resulting trajectory has:
# - Single shared control pulse
# - Multiple state trajectories (one per system)
#
# ## Common Operations
#
# ### Named Trajectory Integration
#
# After solving, the `NamedTrajectory` stores discrete time points:

traj = get_trajectory(qcp)

## Controls at timestep k
u_1 = traj[1][:u]
u_1

## All timesteps
Δts = get_timesteps(traj)
length(Δts)
# ### Internal Representation
#
# Piccolo.jl uses real isomorphisms internally for optimization. Complex states
# are converted to real vectors:
#
# ```julia
# # Complex ket |ψ⟩ → real vector ψ̃
# ψ = ComplexF64[1, im] / √2
# ψ̃ = ket_to_iso(ψ)  # [Re(ψ); Im(ψ)]
#
# # Unitary U → real vector Ũ⃗
# U = GATES[:H]
# Ũ = operator_to_iso_vec(U)
# ```
#
# See [Isomorphisms](@ref isomorphisms-concept) for details.
#
# ## Best Practices
#
# ### 1. Match Pulse Type to Problem
#
# ```julia
# # For SmoothPulseProblem
# pulse = ZeroOrderPulse(controls, times)
# qtraj = UnitaryTrajectory(sys, pulse, U_goal)
# qcp = SmoothPulseProblem(qtraj, N)  # ✓
#
# # For SplinePulseProblem
# pulse = CubicSplinePulse(controls, tangents, times)
# qtraj = UnitaryTrajectory(sys, pulse, U_goal)
# qcp = SplinePulseProblem(qtraj)  # ✓
# ```
#
# ### 2. Initialize with Reasonable Controls

## Scale by drive bounds
max_amp = 0.1 * 1.0
initial_controls = max_amp * randn(2, N)
extrema(initial_controls)

# ## See Also
#
# - [Quantum Systems](@ref quantum-systems) - System definitions
# - [Pulses](@ref pulses-concept) - Control parameterizations
# - [Problem Templates](@ref problem-templates-overview) - Using trajectories in optimization
