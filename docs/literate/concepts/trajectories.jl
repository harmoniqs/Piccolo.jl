# # [Trajectories](@id trajectories-concept)
#
# A trajectory bundles a quantum system, a control pulse, and a goal into the
# central object that problem templates consume.
#
# ## Overview
#
# Every trajectory type defines:
# - The **state** ``x_k`` and its dynamics ``x_{k+1} = \exp(\Delta t_k\, G(\boldsymbol{u}_k))\, x_k``
# - The **initial condition** ``x_1 = x_{\text{init}}``
# - The **goal** ``x_{\text{goal}}`` and associated fidelity metric
#
# ## Trajectory Types
#
# | Type | Dynamics | State ``x_k`` | Fidelity ``F`` | ``\dim(x_k)`` |
# |------|----------|---------------|----------------|---------------|
# | `UnitaryTrajectory` | ``\dot{U} = -iH\,U`` | ``\tilde{U} \in \mathbb{R}^{2d^2}`` | ``\lvert\operatorname{tr}(U_g^\dagger U)\rvert^2 / d^2`` | ``2d^2`` |
# | `KetTrajectory` | ``\dot{\psi} = -iH\,\psi`` | ``\tilde{\psi} \in \mathbb{R}^{2d}`` | ``\lvert\langle\psi_g\mid\psi\rangle\rvert^2`` | ``2d`` |
# | `DensityTrajectory` | ``\dot{\rho} = \mathcal{G}\,\text{vec}(\rho)`` | ``\tilde{\rho} \in \mathbb{R}^{d^2}`` (compact) | ``\operatorname{tr}(\rho_g\,\rho)`` | ``d^2`` |
# | `MultiKetTrajectory` | ``\dot{\psi}_j = -iH\,\psi_j`` | multiple ``\tilde{\psi}_j`` | coherent (see below) | ``2d \times n_{\text{states}}`` |
# | `SamplingTrajectory` | per-system dynamics | per-system states | average fidelity | varies |
#
# The isomorphic state ``x_k`` is always a **real** vector; see
# [Isomorphisms](@ref isomorphisms-concept) for the conversions.
#
# ## UnitaryTrajectory
#
# For synthesizing quantum gates.  The propagator satisfies
# ``\dot{U} = G(\boldsymbol{u})\,U`` with ``U(0) = I``.
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
cached_solve!(qcp, "trajectories_unitary"; max_iter = 50)
fidelity(qcp)

# ### Extracting the Pulse

optimized_pulse = get_pulse(qcp.qtraj)
duration(optimized_pulse)

# ## KetTrajectory
#
# For state preparation: find controls that map
# ``|\psi_{\text{init}}\rangle \to |\psi_{\text{goal}}\rangle`` (up to
# global phase).  The fidelity is
# ``F = |\langle\psi_{\text{goal}}|\psi(T)\rangle|^2``.
#
# ### Construction

## Initial and goal states
ψ_init = ComplexF64[1, 0]  # |0⟩
ψ_goal = ComplexF64[0, 1]  # |1⟩

qtraj_ket = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

# ### Solve

qcp_ket = SmoothPulseProblem(qtraj_ket, N; Q = 100.0)
cached_solve!(qcp_ket, "trajectories_ket"; max_iter = 50)
fidelity(qcp_ket)

# ## MultiKetTrajectory
#
# For gates defined by multiple state mappings with coherent phases.  The
# fidelity enforces phase alignment across all state pairs:
#
# ```math
# F = \left| \frac{1}{n} \sum_{j=1}^{n} \langle \psi_{\text{goal},j} | \psi_j(T) \rangle \right|^2
# ```
#
# This is strictly harder than per-state fidelity because relative phases must
# be correct for the mapping to represent a valid gate.
#
# ### Construction

## Define state pairs: X gate maps |0⟩ → |1⟩ and |1⟩ → |0⟩
ψ0 = ComplexF64[1, 0]
ψ1 = ComplexF64[0, 1]

initial_states = [ψ0, ψ1]
goal_states = [ψ1, ψ0]

qtraj_multi = MultiKetTrajectory(sys, pulse, initial_states, goal_states)

# ### Solve

qcp_multi = SmoothPulseProblem(qtraj_multi, N; Q = 100.0)
cached_solve!(qcp_multi, "trajectories_multi"; max_iter = 50)
fidelity(qcp_multi)

# ## DensityTrajectory
#
# For open quantum systems governed by the Lindblad master equation.  The
# state ``\rho(t)`` evolves as ``\dot{\vec\rho} = \mathcal{G}(\boldsymbol{u})\,\vec\rho``
# where ``\mathcal{G}`` is the Lindbladian superoperator (see [Systems](@ref systems-overview)).
#
# Internally the state is stored in the **compact isomorphism**: a real vector
# of dimension ``d^2`` (not ``2d^2``) that exploits Hermiticity.
# The Lindbladian generators are also compacted via
# ``\mathcal{G}_c = P\,\mathcal{G}\,L`` (size ``d^2 \times d^2`` instead of
# ``2d^2 \times 2d^2``), giving roughly a 4× speedup in integration.
# See [Isomorphisms](@ref isomorphisms-concept) for details.
#
# ### Construction

## Open system with a weak dissipation operator
L = ComplexF64[0.1 0.0; 0.0 0.0]
open_sys = OpenQuantumSystem(
    PAULIS[:Z],
    [PAULIS[:X], PAULIS[:Y]],
    [1.0, 1.0];
    dissipation_operators = [L],
)

## Initial and goal density matrices
ρ_init = ComplexF64[1.0 0.0; 0.0 0.0]  # |0⟩⟨0|
ρ_goal = ComplexF64[0.0 0.0; 0.0 1.0]  # |1⟩⟨1|

T_density, N_density = 10.0, 50
times_density = collect(range(0, T_density, length = N_density))
pulse_density = ZeroOrderPulse(0.1 * randn(2, N_density), times_density)
qtraj_density = DensityTrajectory(open_sys, pulse_density, ρ_init, ρ_goal)

# ### Solve and Analyze
#
# The fidelity for density matrices is ``F = \operatorname{tr}(\rho_{\text{goal}}\,\rho(T))``.

qcp_density = SmoothPulseProblem(qtraj_density, N_density; Q = 100.0, R = 1e-2)
cached_solve!(qcp_density, "trajectories_density"; max_iter = 150)
fidelity(qcp_density)
# ## SamplingTrajectory
#
# For robust optimization over parameter variations.  Created internally by
# `SamplingProblem`, this trajectory type stores multiple system variants
# sharing a single control pulse.  The objective averages fidelity across all
# samples:
#
# ```math
# \bar{\ell} = \frac{1}{S}\sum_{s=1}^{S} \ell\!\bigl(x_N^{(s)},\, x_{\text{goal}}\bigr)
# ```
#
# ```julia
# systems = [sys_nominal, sys_high, sys_low]
# qcp_robust = SamplingProblem(qcp_base, systems)
# ```
#
# ## Common Operations
#
# ### Named Trajectory Integration
#
# After solving, the `NamedTrajectory` stores the NLP decision vector at
# each of the ``N`` knot points:

traj = get_trajectory(qcp)

## Controls at timestep k
u_1 = traj[1][:u]
u_1

## All timesteps
Δts = get_timesteps(traj)
length(Δts)

# ### Internal Representation
#
# States in the trajectory are real isomorphic vectors.  Convert back to
# complex form with:
#
# ```julia
# U = iso_vec_to_operator(traj[:Ũ⃗][:, end])     # unitary
# ψ = iso_to_ket(traj[:ψ̃][:, end])               # ket
# ρ = compact_iso_to_density(traj[:ρ⃗̃][:, end])   # density matrix
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
# - [Quantum Systems](@ref systems-overview) - System definitions
# - [Pulses](@ref pulses-concept) - Control parameterizations
# - [Problem Templates](@ref problem-templates-overview) - Using trajectories in optimization
