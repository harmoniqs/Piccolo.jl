# ============================================================================ #
# Abstract Type
# ============================================================================ #

"""
    AbstractQuantumTrajectory{P<:AbstractPulse}

Abstract type for quantum trajectories that wrap physics (system, pulse, solution, goal).
Parametric on pulse type `P` to enable dispatch in problem templates.

All concrete subtypes should implement:
- `state_name(traj)` — user-facing symbol addressing rollout states
- `isomorphism_state_name(traj)` — internal symbol used in NamedTrajectory
- `drive_name(traj)` — drive/control variable symbol (from pulse)
- `time_name(traj)` — time variable symbol (fixed `:t`)
- `timestep_name(traj)` — timestep variable symbol (fixed `:Δt`)
- `duration(traj)` — duration (from pulse)

Subtypes inherit lifted accessors `traj.t`, `traj.sol`, and
`traj[state_name(traj)]` for reading rollout states without touching the
underlying ODE solution directly.
"""
abstract type AbstractQuantumTrajectory{P<:AbstractPulse} end

export AbstractQuantumTrajectory
export UnitaryTrajectory, KetTrajectory, MultiKetTrajectory, DensityTrajectory
export MultiDensityTrajectory, SamplingTrajectory
export state_name, state_names, isomorphism_state_name, isomorphism_state_names
export drive_name, time_name, timestep_name
export get_system, get_pulse, get_initial, get_goal, get_solution
