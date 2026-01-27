# Concepts Overview

Piccolo.jl is organized into three main modules that work together to enable quantum optimal control:

## Architecture

```
Piccolo.jl
├── Quantum          # Quantum mechanical building blocks
│   ├── Systems      # Hamiltonian representations
│   ├── Trajectories # Time evolution containers
│   ├── Pulses       # Control parameterizations
│   ├── Operators    # Embedded and lifted operators
│   └── Isomorphisms # Real vector representations
│
├── Control          # Optimal control framework
│   ├── Problems     # QuantumControlProblem wrapper
│   ├── Objectives   # Fidelity and regularization
│   ├── Constraints  # Bounds and equality constraints
│   └── Templates    # High-level problem constructors
│
└── Visualizations   # Plotting and analysis
    ├── Trajectories # State and control plots
    └── Populations  # Population dynamics
```

## Workflow Overview

A typical Piccolo.jl workflow follows these steps:

```julia
# 1. Define the quantum system
sys = QuantumSystem(H_drift, H_drives, drive_bounds)

# 2. Create an initial control pulse
pulse = ZeroOrderPulse(initial_controls, times)

# 3. Define the optimization goal via a trajectory
qtraj = UnitaryTrajectory(sys, pulse, U_goal)

# 4. Set up the optimization problem
qcp = SmoothPulseProblem(qtraj, N; Q=100.0, R=1e-2)

# 5. Solve
solve!(qcp; max_iter=100)

# 6. Analyze and use results
fidelity(qcp)
optimized_pulse = get_pulse(qcp.qtraj)
```

## Core Concepts

### [Quantum Systems](@ref)

Quantum systems represent the physical hardware you're controlling. They encapsulate:
- **Drift Hamiltonian** `H₀`: Always-on system dynamics
- **Drive Hamiltonians** `Hᵢ`: Controllable interactions
- **Drive bounds**: Hardware limits on control amplitudes

```julia
# The Hamiltonian: H(u,t) = H_drift + Σᵢ uᵢ(t) H_drives[i]
sys = QuantumSystem(H_drift, H_drives, drive_bounds)
```

### [Trajectories](@ref)

Trajectories combine a system, pulse, and goal to represent a complete optimization task:

| Type | Use Case |
|------|----------|
| `UnitaryTrajectory` | Gate synthesis |
| `KetTrajectory` | State preparation |
| `DensityTrajectory` | Open system evolution |
| `MultiKetTrajectory` | Multiple state transfers |
| `SamplingTrajectory` | Robust optimization |

### [Pulses](@ref)

Pulses parameterize how controls vary in time:

| Type | Description |
|------|-------------|
| `ZeroOrderPulse` | Piecewise constant (for `SmoothPulseProblem`) |
| `LinearSplinePulse` | Linear interpolation between knots |
| `CubicSplinePulse` | Smooth cubic Hermite splines |
| `GaussianPulse` | Parametric Gaussian envelope |

### [Objectives](@ref)

Objectives define what the optimization minimizes:
- **Infidelity objectives**: `1 - F` where `F` is fidelity
- **Regularization**: Penalize large or rapidly changing controls
- **Leakage objectives**: Penalize population outside computational subspace

### [Constraints](@ref)

Constraints define what solutions must satisfy:
- **Bound constraints**: Limits on control values
- **Fidelity constraints**: Minimum fidelity requirements
- **Leakage constraints**: Maximum allowed leakage

## Reexported Packages

Piccolo.jl reexports several foundation packages:

| Package | Purpose |
|---------|---------|
| `DirectTrajOpt` | Trajectory optimization solver |
| `NamedTrajectories` | Named trajectory data structures |
| `TrajectoryIndexingUtils` | Trajectory slicing and indexing |

These are available when you `using Piccolo` without additional imports.

## Next Steps

- **New to Piccolo?** Start with the [Getting Started](@ref) guide
- **Ready to optimize?** See [Problem Templates](@ref) for the main API
- **Need details?** Explore the individual concept pages below

## Concept Pages

- [Quantum Systems](@ref) - Hamiltonian representations
- [Trajectories](@ref) - Time evolution containers
- [Pulses](@ref) - Control parameterizations
- [Objectives](@ref) - Optimization objectives
- [Constraints](@ref) - Problem constraints
- [Operators](@ref) - Embedded and lifted operators
- [Isomorphisms](@ref) - Real vector representations
