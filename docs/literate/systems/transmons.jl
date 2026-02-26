# # [Transmon Qubits](@id transmon-systems)
#
# Superconducting transmon qubits are charge-insensitive artificial atoms formed
# by a Josephson junction shunted by a large capacitance.  The transmon spectrum
# is weakly anharmonic, so at least one level beyond the computational subspace
# must be included to capture leakage.
#
# ## Single Transmon
#
# ### Hamiltonian
#
# In the rotating frame at frequency ``\omega``, the transmon Hamiltonian
# truncated to ``d`` levels is:
#
# ```math
# H = (\omega - \omega_{\text{frame}})\, a^\dagger a
#     \;-\; \frac{\delta}{2}\, a^{\dagger 2} a^2
#     \;+\; u_x(t)\,(a + a^\dagger)
#     \;+\; u_y(t)\, i(a^\dagger - a)
# ```
#
# where ``\delta`` is the anharmonicity (``\approx 200`` MHz for typical
# transmons) and ``a`` is the annihilation operator truncated to ``d`` levels.
# In the default rotating frame ``\omega_{\text{frame}} = \omega``, the first
# term vanishes.
#
# The energy spectrum is
# ``E_n = n\omega + \tfrac{n(n-1)}{2}\,\delta``, giving transition frequencies
# ``\omega_{01} = \omega`` and ``\omega_{12} = \omega + \delta``.
#
# ### Construction

using Piccolo

## 3-level transmon in the rotating frame (default)
sys = TransmonSystem(
    levels = 3,           # Truncation level (≥2)
    δ = 0.2,              # Anharmonicity (GHz)
    ω = 4.0,              # Qubit frequency (GHz) — cancels in rotating frame
    drive_bounds = [0.2, 0.2],  # Max amplitude for X, Y drives (GHz)
)
sys.levels, sys.n_drives

# ### Parameters
#
# | Parameter | Default | Description |
# |-----------|---------|-------------|
# | `levels` | `3` | Number of transmon levels |
# | `δ` | `0.2` | Anharmonicity (GHz) |
# | `ω` | `4.0` | Qubit frequency (GHz) |
# | `drive_bounds` | `[1.0, 1.0]` | Bounds for ``u_x, u_y`` |
# | `lab_frame` | `false` | Use lab frame instead of rotating frame |
# | `lab_frame_type` | `:duffing` | Lab frame model (`:duffing`, `:quartic`, `:cosine`) |
# | `multiply_by_2π` | `true` | Multiply by ``2\pi`` (GHz → rad/ns) |
#
# ### Lab Frame Variants
#
# The `lab_frame_type` keyword selects the lab-frame Hamiltonian model:
#
# | Type | Hamiltonian |
# |------|-------------|
# | `:duffing` | ``H = \omega\, a^\dagger a - \frac{\delta}{2}\, a^{\dagger 2} a^2`` |
# | `:quartic` | ``H = \omega_0\, a^\dagger a - \frac{\delta}{12}\, (a + a^\dagger)^4`` |
# | `:cosine` | ``H = 4 E_C\, \hat{n}^2 - E_J\, \cos\hat{\varphi}`` |
#
# ### Example: X Gate on a 3-Level Transmon

sys = TransmonSystem(levels = 3, δ = 0.2, drive_bounds = [0.2, 0.2])

## Goal: X gate in the computational subspace {|0⟩, |1⟩}
U_goal = EmbeddedOperator(:X, sys)

T, N = 20.0, 100
times = collect(range(0, T, length = N))
pulse = ZeroOrderPulse(0.05 * randn(2, N), times)
qtraj = UnitaryTrajectory(sys, pulse, U_goal)

qcp = SmoothPulseProblem(qtraj, N; Q = 100.0)
cached_solve!(qcp, "systems_transmon_X"; max_iter = 100, print_level = 1)
fidelity(qcp)

# !!! tip "Leakage"
#     With `levels ≥ 3`, use `EmbeddedOperator` to define the gate in the
#     computational subspace, and add leakage penalties via
#     `PiccoloOptions(leakage_cost=10.0)` or
#     `PiccoloOptions(leakage_constraint=true)`.
#     See [Leakage Suppression](@ref leakage-suppression).
#
# ## Multi-Transmon Systems
#
# ### Hamiltonian
#
# ``N`` transmon qubits with pairwise dipole coupling in the rotating frame:
#
# ```math
# H = \sum_{i=1}^{N} \left[
#     (\omega_i - \omega_{\text{frame}})\, a_i^\dagger a_i
#     - \frac{\delta_i}{2}\, a_i^{\dagger 2} a_i^2
# \right]
# + \sum_{i < j} g_{ij} \left( a_i\, a_j^\dagger + a_i^\dagger\, a_j \right)
# ```
#
# The coupling ``g_{ij}(a_i a_j^\dagger + a_i^\dagger a_j)`` is the
# rotating-wave approximation of the transmon dipole interaction.  In the lab
# frame the full form ``g_{ij}(a_i + a_i^\dagger)(a_j + a_j^\dagger)`` is used
# instead.
#
# ### Construction

sys_multi = MultiTransmonSystem(
    [4.0, 4.1],               # Qubit frequencies ωᵢ (GHz)
    [0.2, 0.22],              # Anharmonicities δᵢ (GHz)
    [0.0 0.01; 0.01 0.0];     # Coupling matrix gᵢⱼ (GHz)
    drive_bounds = 0.2,
    levels_per_transmon = 3,
)

# Returns a `CompositeQuantumSystem` with tensor-product Hilbert space
# ``\mathcal{H} = \mathcal{H}_1 \otimes \mathcal{H}_2``.
#
# ### Parameters
#
# | Parameter | Description |
# |-----------|-------------|
# | `ωs` | Qubit frequencies (vector, GHz) |
# | `δs` | Anharmonicities (vector, GHz) |
# | `gs` | Symmetric coupling matrix (GHz) |
# | `levels_per_transmon` | Levels per qubit (default 3) |
# | `drive_bounds` | Scalar or vector of bounds |
# | `subsystems` | Subset of indices to include |
# | `subsystem_drive_indices` | Which qubits get drives |
#
# ## Transmon-Cavity (Circuit QED)
#
# ### Hamiltonian
#
# A transmon dispersively coupled to a microwave cavity:
#
# ```math
# H = \tilde{\Delta}\, \hat{b}^\dagger \hat{b}
#     - \chi\, \hat{a}^\dagger \hat{a}\, \hat{b}^\dagger \hat{b}
#     - \chi'\, \hat{b}^{\dagger 2} \hat{b}^2\, \hat{a}^\dagger \hat{a}
#     - K_q\, \hat{a}^{\dagger 2} \hat{a}^2
#     - K_c\, \hat{b}^{\dagger 2} \hat{b}^2
# ```
#
# where ``\hat{a}`` (``\hat{b}``) annihilates a transmon (cavity) excitation,
# ``\chi`` is the dispersive shift, ``K_q, K_c`` are self-Kerr nonlinearities,
# and ``\tilde{\Delta} = \chi/2``.
#
# The four drives are:
# 1. ``\hat{a} + \hat{a}^\dagger`` — real transmon drive
# 2. ``i(\hat{a}^\dagger - \hat{a})`` — imaginary transmon drive
# 3. ``\hat{b} + \hat{b}^\dagger`` — real cavity drive
# 4. ``i(\hat{b}^\dagger - \hat{b})`` — imaginary cavity drive
#
# ### Construction

sys_cqed = TransmonCavitySystem(
    qubit_levels = 4,
    cavity_levels = 12,
    χ = 2π * 32.8e-6,       # Dispersive shift (GHz)
    K_q = 2π * 193e-3 / 2,  # Qubit self-Kerr (GHz)
    K_c = 2π * 1e-9 / 2,    # Cavity self-Kerr (GHz)
    drive_bounds = [0.5, 0.5, 1.0, 1.0],
)
sys_cqed.levels, sys_cqed.n_drives

# ### Parameters
#
# | Parameter | Default | Description |
# |-----------|---------|-------------|
# | `qubit_levels` | `4` | Transmon Fock states |
# | `cavity_levels` | `12` | Cavity Fock states |
# | `χ` | ``2\pi \times 32.8`` kHz | Dispersive shift |
# | `χ′` | ``2\pi \times 1.5`` Hz | Higher-order correction |
# | `K_q` | ``2\pi \times 96.5`` MHz | Qubit self-Kerr |
# | `K_c` | ``2\pi \times 0.5`` Hz | Cavity self-Kerr |
#
# ## Best Practices
#
# ### 1. Include Enough Levels
#
# For single-qubit gates: ``d \geq 3`` (captures leakage to ``|2\rangle``).
# For high-fidelity (``F > 0.9999``): ``d \geq 4``.

sys_3 = TransmonSystem(levels = 3, δ = 0.2)
sys_4 = TransmonSystem(levels = 4, δ = 0.2)
sys_3.levels, sys_4.levels

# ### 2. Use Realistic Parameters
#
# | Parameter | Typical Value | Unit |
# |-----------|---------------|------|
# | ``\omega / 2\pi`` | 4–5 | GHz |
# | ``\delta / 2\pi`` | 150–250 | MHz |
# | Drive amplitude | 10–50 | MHz |
# | ``g_{ij} / 2\pi`` | 1–10 | MHz |
#
# ### 3. Use EmbeddedOperator for Gates
#
# When ``d > 2``, define the gate in the computational subspace:
#
# ```julia
# U_goal = EmbeddedOperator(:X, sys)  # X on {|0⟩, |1⟩}, identity on |2⟩, ...
# ```
#
# ## References
#
# - Koch et al., "Charge-insensitive qubit design derived from the Cooper pair
#   box," *Phys. Rev. A* **76**, 042319 (2007)
# - Blais et al., "Circuit quantum electrodynamics," *Rev. Mod. Phys.* **93**,
#   025005 (2021)
#
# ## See Also
#
# - [Quantum Systems Overview](@ref systems-overview) — General system API
# - [Leakage Suppression](@ref leakage-suppression) — Handling higher levels
# - [Operators](@ref operators-concept) — `EmbeddedOperator` details
