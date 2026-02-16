# # [Rydberg Atoms](@id rydberg-atom-systems)
#
# Rydberg atom arrays use highly excited atomic states with strong van der Waals
# interactions to implement quantum gates and many-body simulations.  The atoms
# are individually trapped (e.g., in optical tweezers) with global laser drives.
#
# ## Hamiltonian
#
# For ``N`` atoms in a 1D chain with spacing ``d``, the Rydberg Hamiltonian in
# the spin basis ``|g\rangle = |0\rangle``, ``|r\rangle = |1\rangle`` is:
#
# ```math
# H = \frac{\Omega(t)}{2} \sum_i \sigma_i^x
#     \;-\; \Delta(t) \sum_i n_i
#     \;+\; \sum_{i < j} \frac{C_6}{|r_i - r_j|^6}\, n_i\, n_j
# ```
#
# where:
# - ``\Omega(t)`` is the global Rabi frequency (laser drive)
# - ``\Delta(t)`` is the global detuning
# - ``n_i = |r\rangle\langle r|_i`` is the Rydberg population operator
# - ``C_6`` is the van der Waals coefficient (``\approx 862\,690 \times 2\pi``
#   MHz·μm``^6`` for rubidium)
#
# The interaction ``C_6 / |r_i - r_j|^6`` gives rise to the **Rydberg blockade**:
# within the blockade radius ``r_b = (C_6 / \Omega)^{1/6}``, double excitation
# is energetically forbidden.
#
# If the Y drive is included (`ignore_Y_drive = false`), an additional term
# ``\frac{\Omega_y(t)}{2} \sum_i \sigma_i^y`` appears, enabling full
# phase control of the drive.
#
# ## Construction

using Piccolo

sys_rydberg = RydbergChainSystem(
    N = 3,                     # Number of atoms
    C = 862690 * 2π,           # van der Waals coefficient (MHz·μm⁶)
    distance = 8.7,            # Atom spacing (μm)
    cutoff_order = 1,          # 1 = nearest-neighbor only
    drive_bounds = [1.0, 1.0, 1.0],  # [Ωx, Ωy, Δ]
)
sys_rydberg.levels, sys_rydberg.n_drives

# The Hilbert space is ``(\mathbb{C}^2)^{\otimes N}`` with dimension ``2^N``.
#
# ## Parameters
#
# | Parameter | Default | Description |
# |-----------|---------|-------------|
# | `N` | `3` | Number of atoms |
# | `C` | ``862\,690 \times 2\pi`` | Van der Waals coefficient (MHz·μm``^6``) |
# | `distance` | `8.7` | Atom spacing (μm) |
# | `cutoff_order` | `1` | Interaction range (1 = nearest neighbor, 2 = next-nearest) |
# | `all2all` | `true` | Include all-to-all interactions |
# | `ignore_Y_drive` | `false` | Drop the ``\sigma^y`` drive |
# | `local_detune` | `false` | Include local detuning pattern |
#
# ## Interaction Range
#
# The `cutoff_order` parameter controls which interaction terms are included:
#
# | `cutoff_order` | Pairs included | Scaling |
# |----------------|----------------|---------|
# | 1 | nearest-neighbor ``(i, i+1)`` | ``C_6 / d^6`` |
# | 2 | + next-nearest ``(i, i+2)`` | ``C_6 / (2d)^6`` |
# | `all2all = true` | all pairs ``(i, j)`` | ``C_6 / ((j-i) \cdot d)^6`` |
#
# For ``C_6 / d^6 \gg \Omega``, the nearest-neighbor approximation is often
# sufficient.  The ratio ``V_{\text{nn}} / V_{\text{nnn}} = 2^6 = 64``.

## All-to-all interactions (default)
sys_all = RydbergChainSystem(N = 4, all2all = true, drive_bounds = [1.0, 1.0, 1.0])

## Nearest-neighbor only
sys_nn = RydbergChainSystem(N = 4, all2all = false, cutoff_order = 1, drive_bounds = [1.0, 1.0, 1.0])

sys_all.levels, sys_nn.levels

# ## Example: GHZ State Preparation

N_atoms = 3
sys = RydbergChainSystem(N = N_atoms, drive_bounds = [1.0, 1.0, 1.0])

T, N_steps = 5.0, 80
times = collect(range(0, T, length = N_steps))
pulse = ZeroOrderPulse(0.05 * randn(3, N_steps), times)

## GHZ state: (|000⟩ + |111⟩)/√2
ψ_init = zeros(ComplexF64, 2^N_atoms)
ψ_init[1] = 1.0  # |000⟩

ψ_ghz = zeros(ComplexF64, 2^N_atoms)
ψ_ghz[1] = 1 / √2       # |000⟩
ψ_ghz[end] = 1 / √2     # |111⟩

qtraj_rydberg = KetTrajectory(sys, pulse, ψ_init, ψ_ghz)

qcp_rydberg = SmoothPulseProblem(qtraj_rydberg, N_steps; Q = 100.0)
cached_solve!(qcp_rydberg, "systems_rydberg_ghz"; max_iter = 100, print_level = 1)
fidelity(qcp_rydberg)

# ## Typical Parameters
#
# | Parameter | Typical Value | Unit |
# |-----------|---------------|------|
# | ``C_6 / 2\pi`` (Rb) | ``862\,690`` | MHz·μm``^6`` |
# | Atom spacing | 4–10 | μm |
# | ``\Omega / 2\pi`` | 1–10 | MHz |
# | ``\Delta / 2\pi`` | 0–30 | MHz |
# | Blockade radius | 5–10 | μm |
#
# ## Best Practices
#
# ### 1. Check the Blockade Regime
#
# Ensure ``V_{\text{nn}} = C_6 / d^6 \gg \Omega`` for the blockade mechanism
# to be effective.  If ``V \sim \Omega``, the system operates in the
# intermediate (anti-blockade) regime.
#
# ### 2. Start with Nearest-Neighbor Interactions
#
# Use `cutoff_order = 1` or `all2all = false` for initial explorations.
# Add longer-range interactions (`all2all = true`) when needed for accuracy.
#
# ### 3. Scaling
#
# The Hilbert space grows as ``2^N``, so direct Piccolo optimization is
# practical for ``N \leq 8-10`` atoms.
#
# ## References
#
# - Saffman, Walker & Mølmer, "Quantum information with Rydberg atoms,"
#   *Rev. Mod. Phys.* **82**, 2313 (2010)
# - Bernien et al., "Probing many-body dynamics on a 51-atom quantum
#   simulator," *Nature* **551**, 579 (2017)
#
# ## See Also
#
# - [Quantum Systems Overview](@ref systems-overview) — General system API
# - [Trajectories](@ref trajectories-concept) — `KetTrajectory` for state preparation
