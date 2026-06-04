# ```@copybutton
# literate/two_qubit_gate_validation.jl
# ```
#
# # [Two-Qubit Gate Validation](@id two-qubit-gate-validation)
#
# This tutorial walks through synthesizing a CNOT gate (Controlled NOT or Controlled X gate) on
# two coupled transmon qubits using three Pulse types with varying degrees of smoothness (numbers
# continuous derivatives).
#
# Furthermore, we validate the fidelities reported by Piccolo.jl when using
# those pulses by rolling out the same pulses in
# [QuantumToolbox.jl](https://github.com/qutip/QuantumToolbox.jl).
#

# CLAUDE GENERATED DELETE
# We solve the same gate three times:
#
# | Pulse type          | Problem template     |
# |:--------------------|:---------------------|
# | `ZeroOrderPulse`    | `SmoothPulseProblem` |
# | `LinearSplinePulse` | `SplinePulseProblem` |
# | `CubicSplinePulse`  | `SplinePulseProblem` |
#
# Three ideas run through the page:
#
# 1. **Independent validation.** A high *reported* fidelity only means the
#    optimizer satisfied its own discretized dynamics. Rolling the optimized
#    pulse out through QuantumToolbox's adaptive Schrödinger solver (`sesolve`)
#    and recomputing the gate fidelity is an independent check that the pulse
#    really implements the gate — and a regression canary for integrator parity.
# 2. **Integration accuracy depends on the pulse type.** Piccolo's collocation
#    integrator is *exact* for a piecewise-constant (`ZeroOrderPulse`) control,
#    but only approximate for a control that varies *within* a timestep, so the
#    continuous spline pulses need a finer time grid.
# 3. **Bootstrapping.** We solve the zero-order pulse first, use it to seed the
#    linear spline, and use that to seed the cubic spline. Each step starts close
#    to the answer, so the more expensive fine-grid spline solves converge in far
#    fewer iterations.

# ## Setup

using Piccolo
using QuantumToolbox
using CairoMakie
using LinearAlgebra
using Printf
using Random
Random.seed!(1234)  # For reproducibility

# !!! note "A naming clash"
#     Both Piccolo and QuantumToolbox export `fidelity`. With both packages in
#     scope we write `Piccolo.fidelity` for Piccolo's trajectory fidelity. The
#     gate fidelity from the QuantumToolbox rollout we compute by hand below, and
#     for the honest continuous-time Piccolo fidelity we use
#     `unitary_rollout_fidelity` (whose name does not clash).

# ## Step 1: Define the Two-Qubit System
#
# The system Hamiltonian for ``N`` transmon qubits with pairwise dipole coupling in
# the rotating frame is
# ```math
# H_{\textrm{drift}} = \sum_{i=1}^{N} \left[
#     (\omega_i - \omega_{\text{frame}})\, a_i^\dagger a_i
#     - \frac{\delta_i}{2}\, a_i^{\dagger 2} a_i^2
# \right]
# + \sum_{i < j} g_{ij} \left( a_i\, a_j^\dagger + a_i^\dagger\, a_j \right),
# ```
# and the control Hamiltonians are
# ```math
# H_{\textrm{control}}(t) = \sum_{i=1}^{N} \left[ u_{x,i}(t)(a_i + a_i^\dagger) + u_y(t)i(a_i^\dagger - a_i) \right]
# ```
# where ``g_{ij}`` determines the coupling strength between qubits ``i`` and ``j``.
#
# We could build the drift and control Hamiltonians ourselves and
# contruct the quantum system directly, but Piccolo.jl
# provides the `MultiTransmonSystem` to easily construct the system.

ωs = [4.0, 4.1]            # transmon frequencies (GHz)
δs = [0.2, 0.2]            # anharmonicities (GHz) — unused at 2 levels, kept for realism
g = 0.1                    # exchange coupling (GHz) — artificially large (see note)
gs = [0.0 g; g 0.0]

sys = MultiTransmonSystem(ωs, δs, gs; levels_per_transmon = 2, drive_bounds = 0.1)

# 
# !!! note "Accuracy of model"
#     We ignore two [best practices](@ref transmon-best-practices) to reduce the
#     computational cost of the gate synthesis in this tutorial.
#
#     1. We only model two levels of each qubit. This reduces the size of
#        the Hilbert space, which reduces the number of 
#        [decision variables](@ref decision-variables) in the NLP optimization.
#     2. The coupling between the qubits is artificially high, which
#        greatly reduces the pulse duration needed to synthesize the gate with high
#        fidelity.
#
#     When synthesizing a gate for real hardware, at least 3 levels should be
#     modeled for each qubit in order to suppress leakege to excited states
#     outside of the computational basis. Additionally, ``g`` will by in the
#     1-10 MHz range, which will require a much longer gate duration.




# ## Step 2: Define the Gate
#
# Our target is the CNOT gate. Piccolo already defines this gate in `GATES[:CX]`.
# Because we only model two levels per qubit, we could use `GATES[:CX]` as our
# gate directly, but we will use an `EmbeddedOperator` to place the gate in
# the system's Hilbert space and record the indices of the computational
# subspace, so that the tutorial is still applicable when `levels_per_transmon` 
# ``\geq 2``.

U_goal = EmbeddedOperator(GATES[:CX], sys)

# The indices of computational subspace basis states (``|00\rangle, |01\rangle, |10\rangle, |11\rangle``) are:

U_goal.subspace

# We set the gate duration and some options **TODO ARE THESE OPTIONS NECESSARY?** 

T = 10.0    # gate duration (ns)

opts = PiccoloOptions(timesteps_all_equal = true, display = :silent)
