# # [Saving and Loading Pulses](@id saving-loading)
#
# Optimized pulses are valuable — they can take minutes or hours to compute.
# This guide shows how to save pulses to disk and reload them for warm-starting,
# analysis, or deployment to hardware.
#
# ## Why Save Pulses?
#
# - **Warm-starting**: Load a solved pulse as the initial guess for a harder
#   problem (e.g., minimum-time optimization, higher fidelity targets, or
#   robustness sweeps).
# - **Reproducibility**: Archive optimized pulses alongside your scripts so
#   results can be inspected or extended later.
# - **Hardware deployment**: Export the final pulse for upload to an AWG or
#   other control electronics.
# - **Sharing**: Send a `.jld2` file to a collaborator — they can rebuild the
#   trajectory and continue optimizing without re-running the original solve.
#
# ## Saving a Pulse
#
# After solving, extract the optimized pulse with `get_pulse` and save it
# with JLD2:

using Piccolo
using JLD2

## Set up and solve a simple problem
H_drift = PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
sys = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

T, N = 10.0, 100
times = collect(range(0, T, length = N))
pulse = ZeroOrderPulse(0.1 * randn(2, N), times)
qtraj = UnitaryTrajectory(sys, pulse, GATES[:X])

qcp = SmoothPulseProblem(qtraj, N; Q = 100.0, R = 1e-2, ddu_bound = 1.0)
cached_solve!(qcp, "saving_loading_base"; max_iter = 50)

## Extract the optimized pulse
optimized_pulse = get_pulse(qcp.qtraj)

# Save to disk:

jldsave("my_pulse.jld2"; pulse = optimized_pulse)

# You can also save metadata alongside the pulse — fidelity, duration,
# system parameters, or anything else you want to remember:

jldsave("my_pulse_with_metadata.jld2";
    pulse = optimized_pulse,
    fidelity = fidelity(qcp),
    duration = sum(get_timesteps(get_trajectory(qcp))),
    gate = "X",
    system_drift = H_drift,
)

# ## Loading a Pulse
#
# Load the pulse back and use it to construct a new trajectory:

saved_data = load("my_pulse_with_metadata.jld2")
loaded_pulse = saved_data["pulse"]

loaded_pulse

# The loaded pulse is a regular pulse object — you can evaluate it, check its
# duration, and pass it to any trajectory constructor:

duration(loaded_pulse)

#-

n_drives(loaded_pulse)

#-

loaded_pulse(5.0)

# ## Warm-Starting from a Saved Pulse
#
# The most common use case: load a previously optimized pulse and use it as
# the starting point for a new optimization. This works because the pulse
# already carries good control values — the optimizer just needs to refine
# them.

## Load the saved pulse
loaded = load("my_pulse.jld2", "pulse")

## Build a new trajectory using the loaded pulse
qtraj_warm = UnitaryTrajectory(sys, loaded, GATES[:X])

## Optimize — starting from a good initial guess converges faster
qcp_warm = SmoothPulseProblem(qtraj_warm, N; Q = 1000.0, R = 1e-3, ddu_bound = 1.0)
cached_solve!(qcp_warm, "saving_loading_warm"; max_iter = 20)
fidelity(qcp_warm)

# ### Warm-Starting for Minimum-Time Optimization
#
# A particularly powerful pattern: solve a fixed-duration problem first, save
# the pulse, then load it for time-optimal compression.
#
# ```julia
# # Step 1: Solve and save (in one script)
# qcp = SmoothPulseProblem(qtraj, N; Q=100.0, Δt_bounds=(0.01, 0.5))
# solve!(qcp; max_iter=200)
# jldsave("gate_pulse.jld2"; pulse=get_pulse(qcp.qtraj))
#
# # Step 2: Load and minimize time (in another script)
# saved_pulse = load("gate_pulse.jld2", "pulse")
# qtraj = UnitaryTrajectory(sys, saved_pulse, U_goal)
# qcp = SmoothPulseProblem(qtraj, N; Q=100.0, Δt_bounds=(0.01, 0.5))
# qcp_mintime = MinimumTimeProblem(qcp; final_fidelity=0.999)
# solve!(qcp_mintime; max_iter=100)
# ```

# ## Saving and Loading Trajectories
#
# If you need to save the full optimization state (controls, quantum states
# at every timestep, and timesteps), save the `NamedTrajectory` instead:

traj = get_trajectory(qcp)
save("my_trajectory.jld2", traj)

# Load it back:

loaded_traj = load_traj("my_trajectory.jld2")
typeof(loaded_traj)

# The trajectory contains all NLP decision variables. You can inspect
# controls, states, and timesteps:

loaded_traj[:u][:, 1:3]  # First 3 timesteps of controls

# ## Organizing Saved Pulses
#
# For projects with many gates, use a consistent directory structure:
#
# ```
# my_project/
# ├── scripts/
# │   ├── optimize_X.jl
# │   ├── optimize_CZ.jl
# │   └── mintime_CZ.jl      # loads from data/CZ.jld2
# └── data/
#     ├── X.jld2
#     ├── CZ.jld2
#     └── CZ_mintime.jld2
# ```
#
# Each optimization script saves its result to `data/`, and downstream
# scripts (like minimum-time solves) load from there. This creates a clear
# dependency chain and makes results easy to find.

# ## Tips
#
# ### Always Save After a Successful Solve
#
# Optimization can be expensive. Get in the habit of saving immediately after
# `solve!` returns a good result:
#
# ```julia
# solve!(qcp; max_iter=300)
# if fidelity(qcp) > 0.999
#     jldsave("data/my_gate.jld2"; pulse=get_pulse(qcp.qtraj))
# end
# ```
#
# ### Include Metadata
#
# Future-you will thank present-you for saving context:
#
# ```julia
# jldsave("data/CZ.jld2";
#     pulse = get_pulse(qcp.qtraj),
#     fidelity = fidelity(qcp),
#     duration = sum(get_timesteps(get_trajectory(qcp))),
#     gate_name = "CZ",
#     system_config = "2-atom Rydberg, 8.7 μm spacing",
#     date = string(Dates.now()),
# )
# ```
#
# ### Converting Pulse Types Before Saving
#
# If you optimized with `ZeroOrderPulse` but need a smooth pulse for
# hardware, convert before saving:
#
# ```julia
# zop = get_pulse(qcp.qtraj)  # ZeroOrderPulse from SmoothPulseProblem
#
# # Convert to cubic spline
# ctrl = hcat([zop(t) for t in times]...)
# tgts = zeros(size(ctrl))
# for k in 1:(N-1)
#     tgts[:, k] = (ctrl[:, k+1] - ctrl[:, k]) / (times[k+1] - times[k])
# end
# tgts[:, N] = tgts[:, N-1]
# smooth_pulse = CubicSplinePulse(ctrl, tgts, times)
#
# jldsave("data/CZ_smooth.jld2"; pulse=smooth_pulse)
# ```

## Clean up temporary files
rm("my_pulse.jld2"; force = true)
rm("my_pulse_with_metadata.jld2"; force = true)
rm("my_trajectory.jld2"; force = true)

# ## See Also
#
# - [Pulses](@ref pulses-concept) — Pulse types and conversions
# - [Trajectories](@ref trajectories-concept) — Full trajectory structure
# - [SplinePulseProblem](@ref spline-pulse) — Warm-starting with spline pulses
# - [MinimumTimeProblem](@ref minimum-time) — Time-optimal control from saved solutions
