# Piccolo Control Module Context

> AI-friendly context for maintaining consistency. Update this when making significant changes.

## Module Purpose

The `Control` submodule provides **problem templates** for quantum optimal control using direct collocation. It builds on:
- **DirectTrajOpt.jl** - Core optimization infrastructure (objectives, constraints, integrators, NLP solver)
- **Piccolo.Quantum** - Quantum systems, trajectories, pulses, and utilities
- **NamedTrajectories.jl** - Named variable trajectory data structure

The main abstractions are:
1. **QuantumControlProblem** - Wrapper combining quantum trajectory info with optimization problem
2. **Problem Templates** - Constructors that build complete optimization problems

## Canonical Workflow

### Basic Gate Synthesis

```julia
using Piccolo

# 1. Define quantum system (no T_max - duration comes from trajectory)
H_drift = GATES[:Z]
H_drives = [GATES[:X], GATES[:Y]]
drive_bounds = [1.0, 1.0]  # symmetric bounds
sys = QuantumSystem(H_drift, H_drives, drive_bounds)

# 2. Create quantum trajectory (defines the control problem)
U_goal = GATES[:H]  # target: Hadamard gate
T = 10.0  # duration
qtraj = UnitaryTrajectory(sys, U_goal, T)  # creates zero pulse internally

# 3. Build optimization problem
N = 51  # number of timesteps
qcp = SmoothPulseProblem(qtraj, N; Q=100.0, R=1e-2)

# 4. Solve (use 'options' keyword for IpoptOptions)
solve!(qcp; options=IpoptOptions(max_iter=200))

# 5. Extract results
traj = get_trajectory(qcp)
fid = fidelity(qcp)
```

### State Transfer (Ket)

```julia
sys = QuantumSystem(H_drift, H_drives, drive_bounds)

ψ_init = ComplexF64[1.0, 0.0]  # |0⟩
ψ_goal = ComplexF64[0.0, 1.0]  # |1⟩
T = 10.0
qtraj = KetTrajectory(sys, ψ_init, ψ_goal, T)

qcp = SmoothPulseProblem(qtraj, N; Q=50.0, R=1e-3)
solve!(qcp; max_iter=100)
```

### Ensemble Control (Multiple Initial States)

```julia
sys = QuantumSystem(H_drift, H_drives, drive_bounds)

# X gate: |0⟩→|1⟩ and |1⟩→|0⟩
initials = [ComplexF64[1,0], ComplexF64[0,1]]
goals = [ComplexF64[0,1], ComplexF64[1,0]]
T = 10.0
qtraj = MultiKetTrajectory(sys, initials, goals, T)

qcp = SmoothPulseProblem(qtraj, N; Q=100.0, R=1e-2)
solve!(qcp; max_iter=200)
```

### With Explicit Pulse (Spline-Based)

```julia
sys = QuantumSystem(H_drift, H_drives, drive_bounds)

# Create explicit pulse with initial guess
T = 10.0
times = collect(range(0.0, T, length=N))
controls = 0.1 * randn(sys.n_drives, N)
pulse = LinearSplinePulse(controls, times)

qtraj = UnitaryTrajectory(sys, pulse, U_goal)

# Use SplinePulseProblem for spline pulses
qcp = SplinePulseProblem(qtraj, N; Q=100.0, R=1e-2, du_bound=10.0)
solve!(qcp; max_iter=200)
```

### Minimum Time Optimization

```julia
# First solve with fixed time
qcp = SmoothPulseProblem(qtraj, N; Q=100.0, Δt_bounds=(0.01, 0.5))
solve!(qcp; max_iter=100)

# Convert to minimum-time problem
qcp_mintime = MinimumTimeProblem(qcp; final_fidelity=0.99, D=50.0)
solve!(qcp_mintime; max_iter=100)
```

### Bootstrapping Across System Sizes

**Best Practice:** When optimizing similar systems of different sizes, **pass the pulse object** between solves to preserve both controls and derivatives:

```julia
# Optimize for N=3 atoms
sys3 = build_system(N=3)
qtraj3 = UnitaryTrajectory(sys3, pulse_init, U_goal)
qcp3 = SplinePulseProblem(qtraj3, N_timesteps; Q=100.0, R_du=1e-4)
solve!(qcp3; max_iter=200)

# Bootstrap to N=4 using optimized pulse directly
optimized_pulse = qcp3.qtraj.pulse  # Contains both u and du
sys4 = build_system(N=4)
qtraj4 = UnitaryTrajectory(sys4, optimized_pulse, U_goal)  # Reuse pulse
qcp4 = SplinePulseProblem(qtraj4, N_timesteps; Q=100.0, R_du=1e-4)
solve!(qcp4; max_iter=200)

# Save pulse for next bootstrap
using JLD2
save("optimized_N4.jld2", "pulse", qcp4.qtraj.pulse)
```

**Why this works:**
- `Pulse` objects are system-independent (just time → control mappings)
- Preserves full optimization state: `u` (controls) + `du` (derivatives/tangents)
- The new system's dynamics are enforced during `solve!` via rollout constraints

### Robust Control (Parameter Sampling)

```julia
# Create base problem
qtraj = UnitaryTrajectory(sys_nominal, pulse, U_goal)
qcp = SmoothPulseProblem(qtraj, N; Q=100.0)

# Add robustness over parameter variations
sampling_prob = SamplingProblem(qcp, [sys_nominal, sys_perturbed]; Q=100.0)
solve!(sampling_prob; max_iter=200)
```

## Key Abstractions

### QuantumControlProblem

The central type wrapping optimization:

```julia
struct QuantumControlProblem{Q<:AbstractQuantumTrajectory}
    qtraj::Q                    # Quantum trajectory (physics info)
    prob::DirectTrajOptProblem  # Optimization problem (objectives, constraints)
end
```

**Key functions:**
- `solve!(qcp; options=...)` - Solve the optimization
- `get_trajectory(qcp)` - Get optimized NamedTrajectory
- `fidelity(qcp)` - Compute final fidelity (uses updated qtraj)
- `qcp.qtraj.pulse` - Access optimized pulse for bootstrapping

### Problem Templates

| Template | Use Case | Pulse Type |
|----------|----------|------------|
| `SmoothPulseProblem` | Standard optimization with smooth pulses | `ZeroOrderPulse` (piecewise constant) |
| `SplinePulseProblem` | Spline-based pulse optimization | `LinearSplinePulse`, `CubicSplinePulse` |
| `MinimumTimeProblem` | Time-optimal control | Any (converts existing problem) |
| `SamplingProblem` | Robust control over parameters | Any (wraps existing problem) |

**SmoothPulseProblem** (for ZeroOrderPulse):
- Adds derivative variables `:du`, `:ddu` for smoothness
- Creates `DerivativeIntegrator` constraints enforcing `u[k+1] - u[k] = Δt * du[k]`
- Applies quadratic regularization on `u`, `du`, `ddu`
- Supports `global_bounds` for time-invariant parameters (requires custom integrator)

**SplinePulseProblem** (for spline pulses):
- For `LinearSplinePulse`: `:du` represents slopes (added automatically)
- For `CubicSplinePulse`: `:du` represents Hermite tangents (built into pulse)
- Uses `DerivativeIntegrator` with spline semantics
- Supports `global_bounds` for time-invariant parameters (requires custom integrator)

### Global Variable Bounds

Both `SmoothPulseProblem` and `SplinePulseProblem` support bounds on global (time-invariant) optimization variables:

```julia
using Piccolissimo  # For HermitianExponentialIntegrator or SplineIntegrator

# System with global parameter
H = (u, t) -> u[2] * GATES.Z + u[1] * GATES.X  # δ = u[2] is global
sys = QuantumSystem(H, [1.0]; time_dependent=true, global_params=(δ=0.1,))
qtraj = UnitaryTrajectory(sys, pulse, U_goal)

# Integrator that supports globals
integrator = HermitianExponentialIntegrator(qtraj, N)

# Symmetric bounds: -0.5 ≤ δ ≤ 0.5
qcp = SmoothPulseProblem(qtraj, N; 
    integrator=integrator, 
    global_bounds=Dict(:δ => 0.5))

# Asymmetric bounds: 0.1 ≤ δ ≤ 0.8
qcp = SmoothPulseProblem(qtraj, N; 
    integrator=integrator,
    global_bounds=Dict(:δ => (0.1, 0.8)))
```

### Quantum Integrators (`quantum_integrators.jl`)

`BilinearIntegrator` dispatches on trajectory type:
```julia
BilinearIntegrator(qtraj::UnitaryTrajectory, N)      # → single integrator
BilinearIntegrator(qtraj::KetTrajectory, N)          # → single integrator
BilinearIntegrator(qtraj::MultiKetTrajectory, N)     # → Vector of integrators
BilinearIntegrator(qtraj::SamplingTrajectory, N)     # → Vector of integrators
```

**Signature pattern:** `(qtraj, N::Int)` - creates NamedTrajectory internally

### Quantum Constraints (`quantum_constraints.jl`)

Fidelity constraints for minimum-time problems:
- `FinalUnitaryFidelityConstraint(qtraj, traj, fidelity)`
- `FinalKetFidelityConstraint(qtraj, traj, fidelity)`
- `FinalCoherentKetFidelityConstraint(qtraj, traj, fidelity)` - For ensembles, preserves phase coherence

## Component Naming Conventions

| Component | Symbol | Description |
|-----------|--------|-------------|
| Unitary state | `:Ũ⃗` | Isomorphism-vectorized unitary |
| Ket state | `:ψ̃` | Isomorphism-vectorized ket |
| Ensemble kets | `:ψ̃1`, `:ψ̃2`, ... | Multiple states |
| Sampling states | `:Ũ⃗1`, `:Ũ⃗2`, ... | States for each system sample |
| Controls | `:u` | Pulse amplitudes (canonical name) |
| Control derivative | `:du` | First derivative / slope |
| Control 2nd derivative | `:ddu` | Second derivative / acceleration |
| Timestep | `:Δt` | Duration of each timestep |
| Time | `:t` | Accumulated time (always present) |

## File Structure

```
control/
├── _control.jl                    # Submodule definition
├── quantum_control_problem.jl     # QuantumControlProblem wrapper
├── quantum_integrators.jl         # BilinearIntegrator dispatch
├── quantum_constraints.jl         # Fidelity constraints
├── quantum_objectives.jl          # Fidelity objectives
└── problem_templates/
    ├── _problem_templates.jl      # Submodule definition
    ├── smooth_pulse_problem.jl    # ZeroOrderPulse optimization
    ├── spline_pulse_problem.jl    # Spline pulse optimization
    ├── minimum_time_problem.jl    # Time-optimal control
    └── sampling_problem.jl        # Robust control
```

## Common Gotchas

1. **Pulse type determines problem template**: Use `SmoothPulseProblem` for `ZeroOrderPulse`, `SplinePulseProblem` for spline pulses
2. **N is timesteps, not knot points**: For `N=51`, you get 50 intervals
3. **Trajectories always have `:t`**: Time is accumulated automatically, no need for `time_dependent` flag
4. **MultiKetTrajectory vs SamplingTrajectory**: Ensemble = same system, different initial states; Sampling = different systems, same initial state
5. **Fidelity vs Infidelity**: Objectives minimize infidelity, constraints bound infidelity (e.g., `1 - fidelity ≤ 1 - 0.99`)
6. **Bootstrapping between system sizes**: Pass the **pulse object** (`qcp.qtraj.pulse`), not the `NamedTrajectory` or extracted controls. This preserves both `u` and `du` optimization state.
7. **Don't recreate DirectTrajOptProblem**: When modifying a `SplinePulseProblem`, add constraints to the existing problem rather than recreating it, as this breaks the dynamical timestep mechanism.

## Removed Patterns (Don't Reintroduce!)

- ❌ `adapt_trajectory` / `unadapt_trajectory` - Control scaling was removed
- ❌ `@test sys.time_dependent` in tests - Not needed since `:t` is always present
- ❌ `QuantumSystem(...; time_dependent=true)` in tests - Redundant
- ❌ Control name `:a` - Canonical control name is `:u`
