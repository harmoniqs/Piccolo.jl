# # [Trapped Ions](@id trapped-ion-systems)
#
# Trapped-ion quantum computers encode qubits in internal states of individual
# ions confined in electromagnetic traps.  Entangling operations are mediated by
# shared motional (phonon) modes of the ion chain.
#
# ## IonChainSystem
#
# ### Hamiltonian
#
# A chain of ``N`` ions, each with internal levels coupled to ``M`` shared
# motional modes.  In the rotating frame at frequency ``\omega_{\text{frame}}``:
#
# ```math
# H = \sum_{i=1}^{N} (\omega_{q,i} - \omega_{\text{frame}})\, \sigma_i^+ \sigma_i^-
#     \;+\; \sum_{m=1}^{M} \omega_m\, a_m^\dagger a_m
#     \;+\; \sum_{i,m} \eta_{i,m}\, \sigma_i^x\,(a_m + a_m^\dagger)
#     \;+\; \sum_i \bigl[ \Omega_{x,i}(t)\, \sigma_i^x + \Omega_{y,i}(t)\, \sigma_i^y \bigr]
# ```
#
# where:
# - ``\sigma_i^{\pm}`` are raising/lowering operators for ion ``i``
# - ``a_m, a_m^\dagger`` are phonon operators for motional mode ``m``
# - ``\omega_{q,i}`` is the qubit transition frequency
# - ``\omega_m`` is the motional mode frequency
# - ``\eta_{i,m}`` is the **Lamb-Dicke parameter** coupling ion ``i`` to mode ``m``
# - ``\Omega_{x,i}(t), \Omega_{y,i}(t)`` are the laser drive amplitudes
#
# The Hilbert space is ``\mathcal{H} = \bigotimes_i \mathbb{C}^{d_{\text{ion}}}
# \otimes \bigotimes_m \mathbb{C}^{n_{\text{max}}}``, with total dimension
# ``d_{\text{ion}}^N \times n_{\text{max}}^M``.
#
# ### Construction

using Piccolo

## Two ions, one motional mode
sys_ion = IonChainSystem(
    N_ions = 2,
    N_modes = 1,
    mode_levels = 5,        # Fock space truncation
    ωq = 1.0,               # Qubit frequency (GHz)
    ωm = 0.1,               # Motional mode frequency (GHz)
    η = 0.1,                # Lamb-Dicke parameter
    drive_bounds = fill(0.5, 4),  # [Ωx₁, Ωy₁, Ωx₂, Ωy₂]
)
sys_ion.levels, sys_ion.n_drives

# ### Parameters
#
# | Parameter | Default | Description |
# |-----------|---------|-------------|
# | `N_ions` | `2` | Number of ions |
# | `ion_levels` | `2` | Internal levels per ion |
# | `N_modes` | `1` | Number of motional modes |
# | `mode_levels` | `10` | Fock space truncation |
# | `ωq` | `1.0` | Qubit frequency (scalar or vector, GHz) |
# | `ωm` | `0.1` | Mode frequency (scalar or vector, GHz) |
# | `η` | `0.1` | Lamb-Dicke parameter (scalar or ``N \times M`` matrix) |
# | `lab_frame` | `false` | Use lab frame |
# | `multiply_by_2π` | `true` | Multiply Hamiltonian by ``2\pi`` |
#
# ### Example: Two-Ion Entangling Gate

T, N = 50.0, 100
times = collect(range(0, T, length = N))
pulse = ZeroOrderPulse(0.01 * randn(4, N), times)

## MS-type entangling gate: exp(-iπ/4 σx⊗σx) on ions ⊗ |0⟩ on mode
ψ_init = zeros(ComplexF64, sys_ion.levels)
ψ_init[1] = 1.0  # |00⟩ ⊗ |0⟩_phonon

## Target: Bell state in qubit subspace ⊗ ground motional state
MS(θ) = exp(-im * θ / 2 * kron(ComplexF64[0 1; 1 0], ComplexF64[0 1; 1 0]))
ψ_goal_qubit = MS(π / 2) * ComplexF64[1, 0, 0, 0]
ψ_goal = zeros(ComplexF64, sys_ion.levels)
ψ_goal[1:4] = ψ_goal_qubit

qtraj_ion = KetTrajectory(sys_ion, pulse, ψ_init, ψ_goal)

qcp_ion = SmoothPulseProblem(qtraj_ion, N; Q = 100.0)
cached_solve!(qcp_ion, "systems_ion_ms"; max_iter = 50, print_level = 1)
fidelity(qcp_ion)

# ## RadialMSGateSystem
#
# ### Hamiltonian
#
# Specialized system for the **radial-mode Mølmer-Sørensen gate**.  In the
# interaction picture, the time-dependent Hamiltonian is:
#
# ```math
# H(t) = -\frac{i}{2} \sum_{i=1}^{N} \sum_{k=1}^{2N}
#     \eta_{k,i}\, \Omega_i(t)\, \sigma_{x,i}
#     \left( a_k\, e^{-i \delta_k t} - a_k^\dagger\, e^{i \delta_k t} \right)
# ```
#
# where ``k`` indexes the ``2N`` **radial modes** (``N`` modes along each
# transverse axis), ``\Omega_i(t)`` is the Rabi frequency for ion ``i``, and
# ``\delta_k`` is the detuning from mode ``k``.
#
# This produces a time-dependent `QuantumSystem` where `sys.H(u, t)` evaluates
# the Hamiltonian at controls `u` and time `t`.
#
# ### Construction

sys_ms = RadialMSGateSystem(
    N_ions = 2,
    mode_levels = 3,
    ωm_radial = [5.0, 5.0, 5.1, 5.1],  # 2N radial mode frequencies (GHz)
    δ = 0.2,                             # Detuning (GHz)
    η = 0.1,                             # Lamb-Dicke parameter
    drive_bounds = [1.0, 1.0],           # Per-ion amplitude bounds
)
sys_ms.levels, sys_ms.n_drives

# ### Parameters
#
# | Parameter | Default | Description |
# |-----------|---------|-------------|
# | `N_ions` | `2` | Number of ions |
# | `mode_levels` | `5` | Fock states per radial mode |
# | `ωm_radial` | `[5, 5, 5.1, 5.1]` | Radial mode frequencies (vector of length ``2N``) |
# | `δ` | `0.2` | Detuning from sideband (scalar or vector) |
# | `η` | `0.1` | Lamb-Dicke parameter (scalar or ``N \times 2N`` matrix) |
#
# ## RadialMSGateSystemWithPhase
#
# ### Hamiltonian
#
# Extends `RadialMSGateSystem` with **phase controls** ``\phi_i(t)`` for AC
# Stark shift compensation:
#
# ```math
# H(t) = \frac{1}{2} \sum_{i,k} \eta_{k,i}\, \Omega_i(t)
#     \left( \sigma_i^+ e^{i\phi_i(t)} + \sigma_i^- e^{-i\phi_i(t)} \right)
#     \left( a_k\, e^{-i \delta_k t} + a_k^\dagger\, e^{i \delta_k t} \right)
# ```
#
# Off-resonant coupling to spectator modes creates AC Stark shifts
# ``\Delta E_{\text{Stark}} \sim \eta^2 \Omega^2 / \delta_{\text{spectator}}``.
# Modulating ``\phi_i(t)`` cancels this time-varying phase accumulation,
# enabling high-fidelity gates.
#
# Controls: ``[\Omega_1, \phi_1, \Omega_2, \phi_2, \ldots]`` (interleaved
# amplitudes and phases).
#
# ### Construction

sys_ms_phase = RadialMSGateSystemWithPhase(
    N_ions = 2,
    mode_levels = 3,
    ωm_radial = [5.0, 5.0, 5.1, 5.1],
    δ = 0.2,
    η = 0.1,
    amplitude_bounds = [1.0, 1.0],
    phase_bounds = [(-Float64(π), Float64(π)), (-Float64(π), Float64(π))],
)
sys_ms_phase.n_drives  # 4: [Ω₁, φ₁, Ω₂, φ₂]

# ## Typical Parameters
#
# Typical values for ``{}^{171}``Yb``^+`` ions (e.g., Q-SCOUT platform):
#
# | Parameter | Typical Value | Unit |
# |-----------|---------------|------|
# | ``\omega_q / 2\pi`` | 12.6 | GHz |
# | ``\omega_m / 2\pi`` (axial) | 1–3 | MHz |
# | ``\omega_m / 2\pi`` (radial) | 3–6 | MHz |
# | Lamb-Dicke ``\eta`` | 0.05–0.15 | — |
# | Gate time | 50–200 | μs |
# | ``n_{\text{max}}`` | 3–5 | — |
#
# ## Best Practices
#
# ### 1. Truncate Motional Modes Carefully
#
# The phonon Fock space must be large enough to capture displaced states during
# the gate.  Start with `mode_levels = 5` and increase if dynamics
# require larger excursions.
#
# ### 2. Ensure Motional Closure
#
# For MS gates, the motional state must return to vacuum at the end of the gate:
# ``|\alpha_k(\tau)| \approx 0`` for all modes ``k``.  This is a necessary
# condition for separating the qubit and motional degrees of freedom.
#
# ### 3. Use Individual Ion Addressing
#
# Each ion has independent ``\Omega_{x,i}, \Omega_{y,i}`` controls, enabling
# individually optimized Rabi frequencies.
#
# ## References
#
# - Sørensen & Mølmer, "Quantum computation with ions in thermal motion,"
#   *Phys. Rev. Lett.* **82**, 1971 (1999)
# - Sørensen & Mølmer, "Entanglement and quantum computation with ions in
#   thermal motion," *Phys. Rev. A* **62**, 022311 (2000)
# - Mizrahi et al., "Realization and Calibration of Continuously Parameterized
#   Two-Qubit Gates on a Trapped-Ion Quantum Processor," *IEEE TQE* (2024)
#
# ## See Also
#
# - [Quantum Systems Overview](@ref systems-overview) — General system API
# - [Pulses](@ref pulses-concept) — `GaussianPulse` for analytical pulse shaping
