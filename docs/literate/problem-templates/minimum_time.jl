# # [MinimumTimeProblem](@id minimum-time)
#
# `MinimumTimeProblem` converts an existing optimization problem into a time-optimal control problem. It minimizes the total gate duration while maintaining a minimum fidelity constraint.
#
# ## When to Use
#
# Use `MinimumTimeProblem` when:
# - You have a working solution and want to minimize its duration
# - You need the fastest possible gate that achieves a target fidelity
# - You're exploring the fidelity-time trade-off
#
# ## Key Design: Composition Pattern
#
# `MinimumTimeProblem` **wraps** an existing `QuantumControlProblem`. It does not create a problem directly from a trajectory:
#
# ```julia
# # This does NOT work
# qcp_mintime = MinimumTimeProblem(qtraj, N)  # Error!
#
# # This works
# qcp_base = SmoothPulseProblem(qtraj, N; Δt_bounds=(0.01, 0.5))
# solve!(qcp_base)
# qcp_mintime = MinimumTimeProblem(qcp_base; final_fidelity=0.99)
# ```
#
# ## Prerequisites
#
# The base problem **must** have `Δt_bounds` set to enable variable timesteps:
#
# ```julia
# # Enable free-time optimization in the base problem
# qcp_base = SmoothPulseProblem(qtraj, N; Δt_bounds=(0.01, 0.5))
# ```
#
# ## Constructor
#
# ```julia
# MinimumTimeProblem(
#     qcp::QuantumControlProblem;
#     goal = nothing,
#     final_fidelity = 0.99,
#     D = 100.0,
#     piccolo_options = PiccoloOptions()
# )
# ```
#
# ## Parameter Reference
#
# | Parameter | Type | Default | Description |
# |-----------|------|---------|-------------|
# | `qcp` | `QuantumControlProblem` | required | Base problem to convert (must have `Δt_bounds`) |
# | `goal` | `AbstractPiccoloOperator` or `AbstractVector` | `nothing` | Optional new goal (uses base problem's goal if `nothing`) |
# | `final_fidelity` | `Float64` | `0.99` | Minimum fidelity constraint |
# | `D` | `Float64` | `100.0` | Weight on total time objective |
# | `piccolo_options` | `PiccoloOptions` | `PiccoloOptions()` | Solver options |
#
# ## Optimization Problem
#
# `MinimumTimeProblem` sets up:
#
# ```
# minimize:    J_original + D × Σ Δtₖ
# subject to:  Original dynamics and constraints
#              F_final ≥ final_fidelity
#              Δt_min ≤ Δtₖ ≤ Δt_max
# ```
#
# The parameter `D` controls the trade-off between the original objective and time minimization. Higher `D` values prioritize shorter duration.
#
# ## Examples
#
# ### Basic Time-Optimal Gate

using Piccolo

## Setup
H_drift = PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
sys = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

T, N = 20.0, 100
times = collect(range(0, T, length = N))
pulse = ZeroOrderPulse(0.1 * randn(2, N), times)
qtraj = UnitaryTrajectory(sys, pulse, GATES[:X])

## Step 1: Solve base problem with free time enabled
qcp_base = SmoothPulseProblem(qtraj, N; Q = 100.0, Δt_bounds = (0.05, 0.5))
solve!(qcp_base; max_iter = 100)

println("Initial duration: ", sum(get_timesteps(get_trajectory(qcp_base))))
println("Initial fidelity: ", fidelity(qcp_base))

## Step 2: Minimize time
qcp_mintime = MinimumTimeProblem(qcp_base; final_fidelity = 0.99, D = 100.0)
solve!(qcp_mintime; max_iter = 100)

println("Final duration: ", sum(get_timesteps(get_trajectory(qcp_mintime))))
println("Final fidelity: ", fidelity(qcp_mintime))

# ### Exploring Fidelity-Time Trade-off

for target_fidelity in [0.999, 0.99, 0.95, 0.90]
    qcp_mt = MinimumTimeProblem(qcp_base; final_fidelity = target_fidelity)
    solve!(qcp_mt; max_iter = 100)

    dur = sum(get_timesteps(get_trajectory(qcp_mt)))
    actual_fid = fidelity(qcp_mt)
    println("Target: $target_fidelity, Duration: $(round(dur, digits=3)), Achieved: $(round(actual_fid, digits=4))")
end
#
# ### With Robust Optimization
#
# Chain with `SamplingProblem` for robust time-optimal control:
#
# ```julia
# # Base problem
# qcp_base = SmoothPulseProblem(qtraj, N; Q=100.0, Δt_bounds=(0.05, 0.5))
# solve!(qcp_base; max_iter=100)
#
# # Add robustness
# sys_perturbed = QuantumSystem(1.1 * H_drift, H_drives, [1.0, 1.0])
# qcp_robust = SamplingProblem(qcp_base, [sys, sys_perturbed])
# solve!(qcp_robust; max_iter=100)
#
# # Minimize time while maintaining robustness
# qcp_mintime = MinimumTimeProblem(qcp_robust; final_fidelity=0.95)
# solve!(qcp_mintime; max_iter=100)
# ```
#
# ## Fidelity Constraints by Trajectory Type
#
# `MinimumTimeProblem` automatically adds the appropriate fidelity constraint based on the trajectory type:
#
# | Trajectory Type | Constraint Type |
# |-----------------|-----------------|
# | `UnitaryTrajectory` | `FinalUnitaryFidelityConstraint` |
# | `KetTrajectory` | `FinalKetFidelityConstraint` |
# | `MultiKetTrajectory` | `FinalCoherentKetFidelityConstraint` |
# | `DensityTrajectory` | Not yet implemented |
#
# ## Tips
#
# ### Setting D
#
# The `D` parameter controls the weight on total time:
# - **Higher D** (e.g., 1000): Aggressively minimize time, may sacrifice fidelity margin
# - **Lower D** (e.g., 10): More conservative, maintains fidelity buffer
#
# ### Initial Solution Quality
#
# `MinimumTimeProblem` works best when starting from a good solution. If your base problem has low fidelity, solve it first:
#
# ```julia
# # Ensure good initial solution
# qcp_base = SmoothPulseProblem(qtraj, N; Q=1000.0, Δt_bounds=(0.05, 0.5))
# solve!(qcp_base; max_iter=200)
#
# # Only then minimize time
# if fidelity(qcp_base) > 0.99
#     qcp_mintime = MinimumTimeProblem(qcp_base; final_fidelity=0.99)
#     solve!(qcp_mintime; max_iter=100)
# end
# ```
#
# ### Changing the Goal
#
# You can optimize for a different goal without recreating the base problem:

qcp_y = MinimumTimeProblem(qcp_base; goal = GATES[:Y], final_fidelity = 0.99)
solve!(qcp_y; max_iter = 100)
#
# ## See Also
#
# - [SmoothPulseProblem](@ref smooth-pulse) - Base problem for piecewise constant controls
# - [SplinePulseProblem](@ref spline-pulse) - Base problem for spline controls
# - [SamplingProblem](@ref sampling) - Add robustness before minimizing time
# - [Composing Templates](@ref composition) - Advanced composition patterns
