# ============================================================================ #
# Abstract Type
# ============================================================================ #

"""
    AbstractQuantumTrajectory{P<:AbstractPulse}

Abstract type for quantum trajectories that wrap physics (system, pulse, solution, goal).
Parametric on pulse type `P` to enable dispatch in problem templates.

All concrete subtypes should implement:
- `state_name(traj)` - Get the state variable symbol (fixed per type)
- `drive_name(traj)` - Get the drive variable symbol (from pulse)
- `time_name(traj)` - Get the time variable symbol (fixed `:t`)
- `timestep_name(traj)` - Get the timestep variable symbol (fixed `:Δt`)
- `duration(traj)` - Get the duration (from pulse)
"""
abstract type AbstractQuantumTrajectory{P<:AbstractPulse} end

export AbstractQuantumTrajectory
export UnitaryTrajectory, KetTrajectory, MultiKetTrajectory, DensityTrajectory
export SamplingTrajectory
export state_name, state_names, drive_name, time_name, timestep_name
export get_system, get_pulse, get_initial, get_goal, get_solution
