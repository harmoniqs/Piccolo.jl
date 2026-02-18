# # [Objectives](@id objectives-concept)
#
# Objectives define what the optimization minimizes.  The total cost function
# evaluated at the ``N``-point trajectory is a weighted sum of a terminal
# infidelity and running regularization penalties:
#
# ```math
# J(\boldsymbol{z}) = Q \cdot \ell(x_N,\, x_{\text{goal}})
# \;+\; \sum_{k=1}^{N}\!\left(
#     R_u \lVert \boldsymbol{u}_k \rVert^2
#   + R_{du} \lVert \Delta\boldsymbol{u}_k \rVert^2
#   + R_{ddu} \lVert \Delta^2\boldsymbol{u}_k \rVert^2
# \right)
# ```
#
# where ``\boldsymbol{z}`` is the full NLP decision vector (see
# [Concepts Overview](@ref concepts-overview)).
#
# ## Fidelity Objectives
#
# The terminal cost ``\ell = 1 - F`` for a trajectory-dependent fidelity
# ``F``.  Gradients are computed via `ForwardDiff` automatic differentiation
# through the isomorphic state representation.
#
# ### UnitaryInfidelityObjective
#
# Gate synthesis fidelity for a ``d``-level system:
#
# ```math
# F_U = \frac{1}{d^2} \left| \operatorname{tr}(U_{\text{goal}}^\dagger\, U_N) \right|^2
# ```
#
# ```julia
# obj = UnitaryInfidelityObjective(:Ũ⃗, U_goal; Q=100.0)
# ```
#
# ### KetInfidelityObjective
#
# State-transfer fidelity:
#
# ```math
# F_\psi = \left| \langle \psi_{\text{goal}} | \psi_N \rangle \right|^2
# ```
#
# ```julia
# obj = KetInfidelityObjective(:ψ̃, ψ_goal; Q=100.0)
# ```
#
# ### CoherentKetInfidelityObjective
#
# For ``n`` state pairs with phase coherence (used by `MultiKetTrajectory`):
#
# ```math
# F_{\text{coh}} = \left| \frac{1}{n} \sum_{j=1}^{n} \langle \psi_{\text{goal},j} | \psi_{j,N} \rangle \right|^2
# ```
#
# This is strictly harder than per-state fidelity because relative phases
# must be correct.
#
# ```julia
# obj = CoherentKetInfidelityObjective([:ψ̃1, :ψ̃2], [ψ_goal1, ψ_goal2]; Q=100.0)
# ```
#
# ### DensityMatrixInfidelityObjective
#
# For open-system optimization with the compact density isomorphism:
#
# ```math
# F_\rho = \operatorname{tr}(\rho_{\text{goal}}\, \rho_N)
# ```
#
# The state ``\tilde{\rho}_N \in \mathbb{R}^{d^2}`` is converted back to a
# Hermitian matrix via `compact_iso_to_density` before computing the trace.
#
# ```julia
# obj = DensityMatrixInfidelityObjective(:ρ⃗̃, ρ_goal, traj; Q=100.0)
# ```
#
# ### UnitaryFreePhaseInfidelityObjective
#
# When global phase doesn't matter, optimizes over ``\phi``:
#
# ```math
# F = \max_{\phi} \frac{1}{d^2} \left| \operatorname{tr}(e^{i\phi} U_{\text{goal}}^\dagger\, U_N) \right|^2
# ```
#
# ```julia
# obj = UnitaryFreePhaseInfidelityObjective(:Ũ⃗, U_goal; Q=100.0)
# ```
#
# ## Regularization Objectives
#
# Regularization penalizes large or rapidly-varying controls via quadratic
# running costs:
#
# ```math
# J_u = \sum_{k=1}^{N} \lVert \boldsymbol{u}_k \rVert^2, \qquad
# J_{du} = \sum_{k=1}^{N} \lVert \Delta\boldsymbol{u}_k \rVert^2, \qquad
# J_{ddu} = \sum_{k=1}^{N} \lVert \Delta^2\boldsymbol{u}_k \rVert^2
# ```
#
# where ``\Delta\boldsymbol{u}_k = \boldsymbol{u}_k - \boldsymbol{u}_{k-1}``
# are discrete differences.
#
# ```julia
# reg_u   = QuadraticRegularizer(:u, traj, R)
# reg_du  = QuadraticRegularizer(:du, traj, R)
# reg_ddu = QuadraticRegularizer(:ddu, traj, R)
# ```
#
# ### Why Regularize?
#
# 1. **Smoothness**: Derivative regularization encourages smooth pulses
# 2. **Robustness**: Prevents exploiting numerical precision
# 3. **Hardware-friendliness**: Bounded, smooth controls are easier to implement
# 4. **Convergence**: Regularization improves the optimization landscape
#
# ## Leakage Objectives
#
# For multilevel systems, leakage to non-computational states can be penalized.
#
# ### LeakageObjective
#
# ```julia
# op = EmbeddedOperator(:X, sys)  # X gate in computational subspace
# obj = LeakageObjective(:Ũ⃗, op; Q=10.0)
# ```
#
# Penalizes population outside the computational subspace at the final time.
#
# ### Via PiccoloOptions
#
# ```julia
# opts = PiccoloOptions(
#     leakage_constraint=true,
#     leakage_constraint_value=1e-3,
#     leakage_cost=10.0
# )
#
# qcp = SmoothPulseProblem(qtraj, N; piccolo_options=opts)
# ```
#
# ## Using Objectives in Problem Templates
#
# Problem templates automatically set up objectives.  The fidelity objective
# is selected based on the trajectory type:
#
# | Trajectory Type | Default Objective |
# |-----------------|-------------------|
# | `UnitaryTrajectory` | `UnitaryInfidelityObjective` |
# | `KetTrajectory` | `KetInfidelityObjective` |
# | `MultiKetTrajectory` | `CoherentKetInfidelityObjective` |
# | `DensityTrajectory` | `DensityMatrixInfidelityObjective` |
#
# ### Example

using Piccolo

## Set up a system
H_drift = PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
sys = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

T, N = 10.0, 100
times = collect(range(0, T, length = N))
pulse = ZeroOrderPulse(0.1 * randn(2, N), times)
U_goal = GATES[:X]
qtraj = UnitaryTrajectory(sys, pulse, U_goal)

## SmoothPulseProblem automatically creates:
## - UnitaryInfidelityObjective (for UnitaryTrajectory)
## - QuadraticRegularizer for :u, :du, :ddu
qcp = SmoothPulseProblem(
    qtraj,
    N;
    Q = 100.0,      # Fidelity weight
    R_u = 1e-3,     # Control regularization
    R_du = 1e-2,    # First derivative regularization
    R_ddu = 1e-2,   # Second derivative regularization
)
cached_solve!(qcp, "objectives_example"; max_iter = 50)
fidelity(qcp)

# ## Objective Weights
#
# ### The ``Q`` Parameter
#
# ``Q`` scales the infidelity term ``Q \cdot (1 - F)``:
# - **Higher Q** (e.g., 1000): Prioritize fidelity over smoothness
# - **Lower Q** (e.g., 10): Allow more flexibility in controls
#
# ### The ``R`` Parameters
#
# ``R_u``, ``R_{du}``, ``R_{ddu}`` scale the regularization terms:
# - **Higher R**: Smoother, smaller controls
# - **Lower R**: More aggressive controls allowed
#
# ### Balancing Trade-offs

## High fidelity weight
qcp_high_Q =
    SmoothPulseProblem(UnitaryTrajectory(sys, pulse, U_goal), N; Q = 1000.0, R = 1e-2)
cached_solve!(qcp_high_Q, "objectives_high_Q"; max_iter = 100)
fidelity(qcp_high_Q)

# High regularization
qcp_high_R =
    SmoothPulseProblem(UnitaryTrajectory(sys, pulse, U_goal), N; Q = 100.0, R = 0.1)
cached_solve!(qcp_high_R, "objectives_high_R"; max_iter = 100)
fidelity(qcp_high_R)

# ### Typical Starting Values
#
# | Parameter | Typical Range | Starting Point |
# |-----------|---------------|----------------|
# | `Q` | 10 - 10000 | 100 |
# | `R` | 1e-6 - 1.0 | 1e-2 |
# | `R_u` | same as R | R |
# | `R_du` | same as R | R |
# | `R_ddu` | same as R | R |
#
# ## Best Practices
#
# ### 1. Start with Defaults
#
# Problem templates have sensible defaults. Start there:
#
# ```julia
# qcp = SmoothPulseProblem(qtraj, N)  # Uses Q=100, R=1e-2
# ```
#
# ### 2. Tune Q First
#
# If fidelity is too low, increase Q.
#
# ### 3. Tune R if Controls are Problematic
#
# If controls are too noisy or large, increase R.
#
# ### 4. Use Per-Derivative Tuning for Fine Control
#
# ```julia
# qcp = SmoothPulseProblem(
#     qtraj, N;
#     R_u=1e-4,    # Allow larger control values
#     R_du=1e-2,   # Penalize jumps moderately
#     R_ddu=0.1    # Strongly penalize acceleration
# )
# ```
#
# ## See Also
#
# - [Constraints](@ref constraints-concept) - Hard constraints on solutions
# - [Problem Templates](@ref problem-templates-overview) - How objectives are used
# - [SmoothPulseProblem](@ref smooth-pulse) - Parameter reference
