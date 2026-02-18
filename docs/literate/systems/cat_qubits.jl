# # [Cat Qubits](@id cat-qubit-systems)
#
# Cat qubits encode quantum information in superpositions of coherent states of
# a bosonic (cavity) mode, stabilized by engineered dissipation.  Because bit
# flips are exponentially suppressed in the coherent state amplitude
# ``|\alpha|^2``, cat qubits are a promising route to hardware-efficient quantum
# error correction.
#
# Piccolo's `CatSystem` models a **two-mode cat qubit** (cat cavity + buffer
# mode) as an `OpenQuantumSystem` with Lindblad dissipation.
#
# ## Hamiltonian
#
# The drift Hamiltonian includes self-Kerr, cross-Kerr, and two-photon exchange:
#
# ```math
# H = -\frac{\chi_{aa}}{2}\, a^{\dagger 2} a^2
#     \;-\; \frac{\chi_{bb}}{2}\, b^{\dagger 2} b^2
#     \;-\; \chi_{ab}\, a^\dagger a\, b^\dagger b
#     \;+\; g_2\, a^{\dagger 2} b
#     \;+\; g_2^*\, a^2 b^\dagger
# ```
#
# where ``a`` (``b``) annihilates a cat (buffer) photon, ``\chi_{aa},
# \chi_{bb}`` are self-Kerr nonlinearities, ``\chi_{ab}`` is the cross-Kerr,
# and ``g_2`` is the two-photon exchange coupling that stabilizes the cat
# manifold ``\{|\alpha\rangle, |-\alpha\rangle\}``.
#
# The two drives are:
# 1. Buffer displacement: ``b + b^\dagger``
# 2. Kerr correction: ``a^\dagger a``
#
# ## Dissipation
#
# Lindblad jump operators model photon loss in both modes:
#
# ```math
# L_a = \sqrt{\kappa_a}\, a, \qquad L_b = \sqrt{\kappa_b}\, b
# ```
#
# The buffer decay rate ``\kappa_b`` is typically much larger than the cat decay
# rate ``\kappa_a``, ensuring the buffer mode is rapidly reset.  The full
# dynamics follow the Lindblad master equation:
#
# ```math
# \dot{\rho} = -i[H, \rho] + \sum_{k \in \{a,b\}} \left(
#     L_k \rho L_k^\dagger - \tfrac{1}{2}\{L_k^\dagger L_k, \rho\}
# \right)
# ```
#
# This returns an `OpenQuantumSystem`, which is used with `DensityTrajectory`
# for optimization.
#
# ## Construction

using Piccolo

sys_cat = CatSystem(
    cat_levels = 13,       # Cat mode Fock space truncation
    buffer_levels = 3,     # Buffer mode Fock space truncation
    g2 = 0.36,             # Two-photon exchange (MHz·2π)
    χ_aa = -7e-3,          # Cat self-Kerr (MHz·2π)
    χ_bb = -32.0,          # Buffer self-Kerr (MHz·2π)
    χ_ab = 0.79,           # Cross-Kerr (MHz·2π)
    κa = 53e-3,            # Cat decay rate (MHz·2π)
    κb = 13.0,             # Buffer decay rate (MHz·2π)
    drive_bounds = [1.0, 1.0],
)
sys_cat.levels, sys_cat.n_drives

# The Hilbert space dimension is `cat_levels × buffer_levels`.
#
# ## Parameters
#
# | Parameter | Default | Description |
# |-----------|---------|-------------|
# | `cat_levels` | `13` | Cat mode Fock space truncation |
# | `buffer_levels` | `3` | Buffer mode Fock space truncation |
# | `g2` | `0.36` | Two-photon exchange (MHz·``2\pi``) |
# | `χ_aa` | `-7e-3` | Cat self-Kerr (MHz·``2\pi``) |
# | `χ_bb` | `-32` | Buffer self-Kerr (MHz·``2\pi``) |
# | `χ_ab` | `0.79` | Cross-Kerr (MHz·``2\pi``) |
# | `κa` | `53e-3` | Cat photon loss rate (MHz·``2\pi``) |
# | `κb` | `13` | Buffer photon loss rate (MHz·``2\pi``) |
# | `prefactor` | `1` | Global scaling of all couplings and rates |
#
# ## Helper Functions
#
# ### Coherent States

## Coherent state |α⟩ in the Fock basis
α = 2.0
ψ_coherent = coherent_ket(α, 13)
length(ψ_coherent)

# The coherent state ``|\alpha\rangle = e^{-|\alpha|^2/2} \sum_{n=0}^{d-1}
# \frac{\alpha^n}{\sqrt{n!}} |n\rangle`` is the starting point for cat state
# preparation.
#
# ### Steady-State Controls

## Controls that maintain |α⟩ in the cat mode
N_steps = 100
u_steady = get_cat_controls(sys_cat, α, N_steps)
u_steady[:, 1]

# `get_cat_controls(sys, α, N)` returns the buffer drive and Kerr correction
# that balance the Hamiltonian at coherent amplitude ``\alpha``, reading the
# coupling constants from `sys.global_params`:
# - Buffer drive: ``g_2 |\alpha|^2``
# - Kerr correction: ``\chi_{aa}(2|\alpha|^2 + 1)``
#
# ## Example: Density Matrix Optimization

using LinearAlgebra

## Small system for demonstration
cat_levels, buffer_levels = 3, 2
sys = CatSystem(cat_levels = cat_levels, buffer_levels = buffer_levels)
n = sys.levels

## Initial state: vacuum |0⟩ ⊗ |0⟩
ρ0 = zeros(ComplexF64, n, n)
ρ0[1, 1] = 1.0

## Goal: coherent state |α⟩ ⊗ |0⟩
α = 0.5
ψ_cat = coherent_ket(α, cat_levels)
ψ_buf = ComplexF64[1, 0]
ψ_goal = kron(ψ_cat, ψ_buf)
ρg = ψ_goal * ψ_goal'

T, N_steps = 1.0, 11
times = collect(range(0, T, length = N_steps))

## Initialize with steady-state controls + small perturbation
u_init = get_cat_controls(sys, α, N_steps) + 0.01 * randn(2, N_steps)
pulse = ZeroOrderPulse(u_init, times)

## DensityTrajectory for open system optimization
qtraj = DensityTrajectory(sys, pulse, ρ0, ρg)

qcp = SmoothPulseProblem(qtraj, N_steps; Q = 100.0, R = 1e-2)
cached_solve!(qcp, "systems_cat_density"; max_iter = 50, print_level = 1)
fidelity(qcp)

# ## Compact Representation
#
# Because ``\rho`` is Hermitian, `DensityTrajectory` uses the **compact
# isomorphism** internally: a real vector of dimension ``d^2`` instead of
# ``2d^2``.  The compact Lindbladian generators
# ``\mathcal{G}_c = P\,\mathcal{G}\,L`` are ``d^2 \times d^2`` (instead of
# ``2d^2 \times 2d^2``), giving roughly a ``4\times`` speedup.
# See [Isomorphisms](@ref isomorphisms-concept) for details.
#
# ## Typical Parameters
#
# | Parameter | Typical Value | Unit |
# |-----------|---------------|------|
# | ``g_2 / 2\pi`` | 0.1–1.0 | MHz |
# | ``\chi_{aa} / 2\pi`` | 1–10 | kHz |
# | ``\chi_{bb} / 2\pi`` | 10–50 | MHz |
# | ``\chi_{ab} / 2\pi`` | 0.5–2.0 | MHz |
# | ``\kappa_a / 2\pi`` | 10–100 | kHz |
# | ``\kappa_b / 2\pi`` | 5–20 | MHz |
# | Coherent amplitude |``\alpha``| | 1–3 | — |
# | Cat levels | 10–20 | — |
# | Buffer levels | 2–4 | — |
#
# ## Best Practices
#
# ### 1. Initialize with Steady-State Controls
#
# Use `get_cat_controls` to warm-start the optimization.  Adding small random
# perturbations helps the optimizer escape local minima.
#
# ### 2. Truncate Carefully
#
# The cat mode Fock space must be large enough to support the coherent state:
# ``n_{\max} \gtrsim 2|\alpha|^2 + 5``.  For ``|\alpha| = 2``, use at least
# ``n_{\max} = 13``.
#
# ### 3. Use DensityTrajectory
#
# Since `CatSystem` returns an `OpenQuantumSystem`, it must be used with
# `DensityTrajectory` (not `UnitaryTrajectory`).  The fidelity is
# ``F = \operatorname{tr}(\rho_{\text{goal}}\, \rho(T))``.
#
# ## References
#
# - Mirrahimi et al., "Dynamically protected cat-qubits: a new paradigm for
#   universal quantum computation," *New J. Phys.* **16**, 045014 (2014)
# - Grimm et al., "Stabilization and operation of a Kerr-cat qubit,"
#   *Nature* **584**, 205 (2020)
# - Lescanne et al., "Exponential suppression of bit-flips in a qubit encoded
#   in an oscillator," *Nature Phys.* **16**, 509 (2020)
#
# ## See Also
#
# - [Quantum Systems Overview](@ref systems-overview) — General system API
# - [Isomorphisms](@ref isomorphisms-concept) — Compact density representation
# - [Trajectories](@ref trajectories-concept) — `DensityTrajectory` for open systems
