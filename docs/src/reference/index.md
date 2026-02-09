# API Reference

This section provides complete documentation of all public types and functions in Piccolo.jl.

## Module Organization

Piccolo.jl is organized into three main modules:

| Module | Contents |
|--------|----------|
| Quantum | Systems, trajectories, pulses, operators, isomorphisms |
| Control | Problems, objectives, constraints, integrators, templates |
| Visualizations | Plotting and animation utilities |

See the [Library Reference](@ref) page for full docstrings of all exported types and functions.

## Reexported Packages

Piccolo.jl reexports several foundation packages. Their functionality is available when you `using Piccolo`:

| Package | Purpose | Key Types/Functions |
|---------|---------|---------------------|
| `DirectTrajOpt` | Trajectory optimization | `solve!`, `add_objective!`, `add_constraint!` |
| `NamedTrajectories` | Trajectory containers | `NamedTrajectory`, `get_timesteps`, `plot` |
| `TrajectoryIndexingUtils` | Indexing utilities | Trajectory slicing and indexing |

## Quick Reference

### Creating Systems

```julia
# Basic system
sys = QuantumSystem(H_drift, H_drives, drive_bounds)

# From templates
sys = TransmonSystem(levels=3, δ=0.2, drive_bounds=[0.2, 0.2])
```

### Creating Trajectories

```julia
# Unitary gate synthesis
qtraj = UnitaryTrajectory(sys, pulse, U_goal)

# State preparation
qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)
```

### Creating Problems

```julia
# Smooth pulse optimization
qcp = SmoothPulseProblem(qtraj, N; Q=100.0)

# Time-optimal control
qcp_mintime = MinimumTimeProblem(qcp; final_fidelity=0.99)

# Robust optimization
qcp_robust = SamplingProblem(qcp, systems)
```

### Solving and Analysis

```julia
# Solve
solve!(qcp; max_iter=100)

# Analyze
fid = fidelity(qcp)
traj = get_trajectory(qcp)
```
