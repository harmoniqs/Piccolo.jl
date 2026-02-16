# # [Silicon Spins](@id silicon-spin-systems)
#
# Silicon spin qubits use electron or hole spins in semiconductor quantum dots
# for quantum information processing.  **Exchange-only (EO) qubits** encode a
# logical qubit in three electron spins, using only exchange interactions
# (controlled by gate voltages) for universal single-qubit control.
#
# !!! note "Coming Soon"
#     The `SiliconSpinSystem` template is under development.  The Hamiltonian
#     and API below describe the planned implementation.
#
# ## Exchange-Only Qubit
#
# ### Hamiltonian
#
# Three electron spins in a triple quantum dot with pairwise Heisenberg
# exchange:
#
# ```math
# H = J_{12}(\epsilon)\, \boldsymbol{S}_1 \cdot \boldsymbol{S}_2
#     \;+\; J_{23}(\epsilon)\, \boldsymbol{S}_2 \cdot \boldsymbol{S}_3
# ```
#
# where ``J_{ij}(\epsilon)`` are the exchange couplings controlled by detuning
# voltages ``\epsilon``.  In the ``S = 1/2`` subspace (total spin ``S_{\text{tot}}
# = 1/2``), this gives a logical qubit with two states:
#
# ```math
# |0_L\rangle = |S_{12}\rangle|\!\uparrow\rangle, \qquad
# |1_L\rangle = \sqrt{\tfrac{1}{3}}|T_{12}\rangle|\!\downarrow\rangle
#     - \sqrt{\tfrac{2}{3}}|T_{12}^0\rangle|\!\uparrow\rangle
# ```
#
# where ``|S_{12}\rangle`` and ``|T_{12}\rangle`` are singlet and triplet states
# of the first two spins.
#
# ### Projected Hamiltonian
#
# Projecting into the logical subspace, the effective Hamiltonian is:
#
# ```math
# H_{\text{eff}} = \frac{J_{12} + J_{23}}{4}\, I
#     + \frac{J_{12} - J_{23}}{4}\, \sigma_z
#     + \frac{\sqrt{3}}{4}\, J_{23}\, \sigma_x
# ```
#
# The exchange couplings ``J_{12}`` and ``J_{23}`` serve as control
# parameters, providing full single-qubit control via ``\sigma_x`` and
# ``\sigma_z`` rotations.
#
# ## Two Exchange-Only Qubits
#
# For two EO qubits (six electrons in six dots), the coupling between logical
# qubits arises from exchange between the adjacent dots of each triple:
#
# ```math
# H = H_{\text{EO},1} + H_{\text{EO},2} + J_{36}\, \boldsymbol{S}_3 \cdot \boldsymbol{S}_4
# ```
#
# where ``J_{36}`` is the inter-qubit exchange coupling (between dot 3 of
# qubit 1 and dot 4 of qubit 2).
#
# ## Planned API
#
# ```julia
# # Single EO qubit
# sys = ExchangeOnlyQubit(
#     J12_bounds = (0.0, 10.0),   # Exchange coupling range (GHz)
#     J23_bounds = (0.0, 10.0),
# )
#
# # Two coupled EO qubits
# sys = TwoExchangeOnlyQubits(
#     J12_bounds = (0.0, 10.0),
#     J23_bounds = (0.0, 10.0),
#     J45_bounds = (0.0, 10.0),
#     J56_bounds = (0.0, 10.0),
#     J34_bounds = (0.0, 5.0),    # Inter-qubit coupling
# )
# ```
#
# ## Typical Parameters
#
# | Parameter | Typical Value | Unit |
# |-----------|---------------|------|
# | Exchange ``J`` | 0.1–10 | GHz |
# | Gate time | 1–100 | ns |
# | ``T_2^*`` (dephasing) | 1–10 | μs |
# | ``T_1`` (relaxation) | 0.1–10 | ms |
#
# ## References
#
# - DiVincenzo et al., "Universal quantum computation with the exchange
#   interaction," *Nature* **408**, 339 (2000)
# - Weinstein et al., "Universal logic with encoded spin qubits in silicon,"
#   *Nature* **615**, 817 (2023)
#
# ## See Also
#
# - [Quantum Systems Overview](@ref systems-overview) — General system API
