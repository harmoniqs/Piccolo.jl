# # [SamplingProblem](@id sampling)
#
# `SamplingProblem` enables robust optimization over multiple system variants with shared controls. This is essential for designing pulses that perform well despite parameter uncertainty or variation.
#
# ## When to Use
#
# Use `SamplingProblem` when:
# - Your system parameters have uncertainty (e.g., qubit frequency drift)
# - You need pulses robust to fabrication variations
# - You want to optimize across an ensemble of similar systems
# - You're designing calibration-free gates
#
# ## Key Design: Composition Pattern
#
# `SamplingProblem` **wraps** an existing `QuantumControlProblem` and extends it with multiple system variants:
#
# ```julia
# # Create base problem
# qcp_base = SmoothPulseProblem(qtraj, N; Q=100.0)
# solve!(qcp_base; max_iter=100)
#
# # Create perturbed systems
# sys_nominal = get_system(qcp_base)
# sys_high = QuantumSystem(1.05 * H_drift, H_drives, drive_bounds)
# sys_low = QuantumSystem(0.95 * H_drift, H_drives, drive_bounds)
#
# # Robust optimization
# qcp_robust = SamplingProblem(qcp_base, [sys_nominal, sys_high, sys_low])
# solve!(qcp_robust; max_iter=100)
# ```
#
# ## Constructor
#
# ```julia
# SamplingProblem(
#     qcp::QuantumControlProblem,
#     systems::Vector{<:AbstractQuantumSystem};
#     weights = fill(1.0, length(systems)),
#     Q = 100.0,
#     piccolo_options = PiccoloOptions()
# )
# ```
#
# ## Parameter Reference
#
# | Parameter | Type | Default | Description |
# |-----------|------|---------|-------------|
# | `qcp` | `QuantumControlProblem` | required | Base problem providing trajectory structure |
# | `systems` | `Vector{AbstractQuantumSystem}` | required | System variants to optimize over |
# | `weights` | `Vector{Float64}` | `fill(1.0, length(systems))` | Relative importance of each system |
# | `Q` | `Float64` | `100.0` | Infidelity weight (applied to all systems) |
# | `piccolo_options` | `PiccoloOptions` | `PiccoloOptions()` | Solver options |
#
# ## Optimization Objective
#
# `SamplingProblem` creates a weighted sum of objectives:
#
# ```
# minimize: Σᵢ wᵢ × Qᵢ × (1 - Fᵢ) + regularization
# ```
#
# Where:
# - `wᵢ` is the weight for system `i`
# - `Qᵢ` is the infidelity weight
# - `Fᵢ` is the fidelity for system `i`
#
# All systems share the same control pulse, but each has its own state trajectory.
#
# ## Examples
#
# ### Robust to Frequency Drift

using Piccolo

## Nominal system
H_drift = 0.5 * PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
sys_nominal = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

## Create trajectory and base problem
T, N = 10.0, 100
times = collect(range(0, T, length = N))
pulse = ZeroOrderPulse(0.1 * randn(2, N), times)
qtraj = UnitaryTrajectory(sys_nominal, pulse, GATES[:X])

qcp_base = SmoothPulseProblem(qtraj, N; Q = 100.0)
cached_solve!(qcp_base, "sampling_base"; max_iter = 100)

## ±5% frequency variation
sys_high = QuantumSystem(1.05 * H_drift, H_drives, [1.0, 1.0])
sys_low = QuantumSystem(0.95 * H_drift, H_drives, [1.0, 1.0])

## Robust optimization
qcp_robust = SamplingProblem(
    qcp_base,
    [sys_nominal, sys_high, sys_low];
    weights = [1.0, 1.0, 1.0],
)
cached_solve!(qcp_robust, "sampling_robust"; max_iter = 100)

# ### Weighted Sampling
#
# Prioritize certain parameter values:

qcp_weighted = SamplingProblem(
    qcp_base,
    [sys_nominal, sys_high, sys_low];
    weights = [2.0, 1.0, 1.0],  ## Nominal weighted 2x
)
cached_solve!(qcp_weighted, "sampling_weighted"; max_iter = 100)

# ### Dense Parameter Sampling
#
# For smooth performance across a parameter range:

scales = range(0.9, 1.1, length = 3)
systems = [QuantumSystem(s * H_drift, H_drives, [1.0, 1.0]) for s in scales]

qcp_dense = SamplingProblem(qcp_base, systems)
cached_solve!(qcp_dense, "sampling_dense"; max_iter = 100)
#
# ### With Time Optimization
#
# Chain with `MinimumTimeProblem`:
#
# ```julia
# # Base problem with free time
# qcp_base = SmoothPulseProblem(qtraj, N; Q=100.0, Δt_bounds=(0.05, 0.5))
# solve!(qcp_base; max_iter=100)
#
# # Add robustness
# qcp_robust = SamplingProblem(qcp_base, [sys_nominal, sys_high, sys_low])
# solve!(qcp_robust; max_iter=100)
#
# # Minimize time
# qcp_mintime = MinimumTimeProblem(qcp_robust; final_fidelity=0.95)
# solve!(qcp_mintime; max_iter=100)
# ```
#
# ## Trajectory Structure
#
# `SamplingProblem` creates a `SamplingTrajectory` internally with:
#
# | Variable | Description |
# |----------|-------------|
# | `:u` | Shared control values |
# | `:Ũ⃗1`, `:Ũ⃗2`, ... | State for each system (unitary case) |
# | or `:ψ̃1`, `:ψ̃2`, ... | State for each system (ket case) |
#
# Note: Derivative variables (`:du`, `:ddu`) from the base problem are **not** carried over. The robustness is achieved through multiple dynamics integrators, one per system.
#
# ## Difference from MultiKetTrajectory
#
# These serve different purposes:
#
# | Feature | `SamplingProblem` | `MultiKetTrajectory` |
# |---------|-------------------|---------------------|
# | Systems | Multiple different systems | Single system |
# | States | Same goal across systems | Different initial/goal pairs |
# | Use case | Parameter uncertainty | Multi-state gates |
#
# **`SamplingProblem`**: "Same gate, different systems"
# **`MultiKetTrajectory`**: "Same system, different state transfers"
#
# ## Evaluating Robustness
#
# After solving, check fidelity across the parameter range:
#
# ```julia
# # Sample more densely for evaluation
# eval_scales = range(0.8, 1.2, length=21)
# eval_systems = [QuantumSystem(s * H_drift, H_drives, drive_bounds) for s in eval_scales]
#
# # Get optimized pulse
# optimized_pulse = get_pulse(qcp_robust.qtraj)
#
# # Evaluate
# for (s, sys) in zip(eval_scales, eval_systems)
#     qtraj_eval = UnitaryTrajectory(sys, optimized_pulse, GATES[:X])
#     fid = fidelity(qtraj_eval)
#     println("Scale $s: Fidelity = $fid")
# end
# ```
#
# ## See Also
#
# - [SmoothPulseProblem](@ref smooth-pulse) - Base problem for piecewise constant controls
# - [MinimumTimeProblem](@ref minimum-time) - Minimize duration of robust pulses
# - [Composing Templates](@ref composition) - Advanced composition patterns
