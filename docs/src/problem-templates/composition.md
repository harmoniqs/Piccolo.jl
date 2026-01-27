# [Composing Templates](@id composition)

Problem templates in Piccolo.jl are designed to be composable. You can chain them together to build sophisticated optimization pipelines that combine multiple capabilities.

## Composition Overview

The templates form a hierarchy:

```
Base Problems (create from trajectory):
├── SmoothPulseProblem
└── SplinePulseProblem
        │
        ▼
Wrapper Problems (wrap existing problem):
├── SamplingProblem (adds robustness)
└── MinimumTimeProblem (adds time optimization)
```

Any wrapper can wrap another wrapper, enabling combinations like:
- `MinimumTimeProblem(SamplingProblem(SmoothPulseProblem(...)))`
- `MinimumTimeProblem(SplinePulseProblem(...))`

## Full Pipeline Example

Here's a complete pipeline: base optimization → robust optimization → time-optimal control.

```julia
using Piccolo

# ============================================================
# Step 1: Setup
# ============================================================

# Define nominal system
H_drift = 0.5 * PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
sys_nominal = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

# Create initial trajectory
T, N = 20.0, 100
times = collect(range(0, T, length=N))
pulse = ZeroOrderPulse(0.1 * randn(2, N), times)
qtraj = UnitaryTrajectory(sys_nominal, pulse, GATES[:X])

# ============================================================
# Step 2: Base Problem (with free time enabled)
# ============================================================

qcp_base = SmoothPulseProblem(
    qtraj, N;
    Q=100.0,
    R=1e-2,
    Δt_bounds=(0.05, 0.5)  # Required for MinimumTimeProblem
)
solve!(qcp_base; max_iter=100)

println("Step 2 - Base:")
println("  Fidelity: ", fidelity(qcp_base))
println("  Duration: ", sum(get_timesteps(get_trajectory(qcp_base))))

# ============================================================
# Step 3: Add Robustness
# ============================================================

# Create perturbed systems (±5% drift variation)
sys_high = QuantumSystem(1.05 * H_drift, H_drives, [1.0, 1.0])
sys_low = QuantumSystem(0.95 * H_drift, H_drives, [1.0, 1.0])

qcp_robust = SamplingProblem(
    qcp_base,
    [sys_nominal, sys_high, sys_low];
    Q=100.0
)
solve!(qcp_robust; max_iter=100)

println("\nStep 3 - Robust:")
println("  Fidelity (nominal): ", fidelity(qcp_robust))

# ============================================================
# Step 4: Minimize Time
# ============================================================

qcp_mintime = MinimumTimeProblem(
    qcp_robust;
    final_fidelity=0.95,
    D=100.0
)
solve!(qcp_mintime; max_iter=100)

println("\nStep 4 - Time-Optimal:")
println("  Fidelity: ", fidelity(qcp_mintime))
println("  Duration: ", sum(get_timesteps(get_trajectory(qcp_mintime))))
```

## Common Composition Patterns

### Pattern 1: Robust Gate

Optimize for parameter uncertainty without time constraints.

```julia
qcp_base = SmoothPulseProblem(qtraj, N; Q=100.0)
solve!(qcp_base)

qcp_robust = SamplingProblem(qcp_base, systems)
solve!(qcp_robust)
```

### Pattern 2: Fast Gate

Minimize time without robustness requirements.

```julia
qcp_base = SmoothPulseProblem(qtraj, N; Q=100.0, Δt_bounds=(0.01, 0.5))
solve!(qcp_base)

qcp_fast = MinimumTimeProblem(qcp_base; final_fidelity=0.99)
solve!(qcp_fast)
```

### Pattern 3: Fast + Robust Gate

The full pipeline for production-quality gates.

```julia
qcp_base = SmoothPulseProblem(qtraj, N; Q=100.0, Δt_bounds=(0.01, 0.5))
solve!(qcp_base)

qcp_robust = SamplingProblem(qcp_base, systems)
solve!(qcp_robust)

qcp_final = MinimumTimeProblem(qcp_robust; final_fidelity=0.95)
solve!(qcp_final)
```

### Pattern 4: Spline Warm-Start Pipeline

Start with smooth problem, refine with splines.

```julia
# Initial optimization with piecewise constant
qcp_smooth = SmoothPulseProblem(qtraj_smooth, N; Q=100.0)
solve!(qcp_smooth)

# Extract optimized pulse and convert to spline
optimized_pulse = get_pulse(qcp_smooth.qtraj)
spline_pulse = CubicSplinePulse(optimized_pulse)
qtraj_spline = UnitaryTrajectory(sys, spline_pulse, U_goal)

# Refine with spline problem
qcp_spline = SplinePulseProblem(qtraj_spline; Q=100.0)
solve!(qcp_spline; max_iter=50)  # Quick refinement
```

## Iteration and Refinement

You can iteratively refine solutions:

```julia
# First pass: coarse optimization
qcp = SmoothPulseProblem(qtraj, 50; Q=10.0)
solve!(qcp; max_iter=50)

# Second pass: increase resolution
# ... resample to higher N ...

# Third pass: tighten tolerances
solve!(qcp; max_iter=100, tol=1e-8)
```

## Accessing Results Through the Chain

Each wrapper preserves access to the underlying trajectory:

```julia
# Access trajectory at any level
traj = get_trajectory(qcp_mintime)
sys = get_system(qcp_mintime)

# Fidelity evaluation uses the innermost trajectory type
fid = fidelity(qcp_mintime)

# Get optimized pulse
pulse = get_pulse(qcp_mintime.qtraj)
```

## Order Matters

The order of composition affects the optimization:

**`MinimumTimeProblem(SamplingProblem(base))`**:
- First achieves robustness, then minimizes time
- Time minimization respects the robust solution
- Generally preferred for production gates

**`SamplingProblem(MinimumTimeProblem(base))`**:
- First minimizes time, then adds robustness
- May require re-solving if time-optimal solution isn't robust
- Less common, but useful for exploring trade-offs

## Tips for Complex Pipelines

### 1. Solve Each Stage

Always `solve!` after each composition step:

```julia
qcp_base = SmoothPulseProblem(qtraj, N)
solve!(qcp_base)  # Important!

qcp_robust = SamplingProblem(qcp_base, systems)
solve!(qcp_robust)  # Important!

qcp_mintime = MinimumTimeProblem(qcp_robust)
solve!(qcp_mintime)  # Important!
```

### 2. Monitor Progress

Check fidelity at each stage to catch issues early:

```julia
for (name, qcp) in [("Base", qcp_base), ("Robust", qcp_robust), ("MinTime", qcp_mintime)]
    println("$name: Fidelity = $(fidelity(qcp))")
end
```

### 3. Adjust Parameters at Each Stage

Different stages may need different settings:

```julia
# Base: prioritize finding a solution
qcp_base = SmoothPulseProblem(qtraj, N; Q=100.0, R=1e-2)
solve!(qcp_base; max_iter=200)

# Robust: may need more iterations
qcp_robust = SamplingProblem(qcp_base, systems; Q=100.0)
solve!(qcp_robust; max_iter=300)

# MinTime: typically faster since starting from good solution
qcp_mintime = MinimumTimeProblem(qcp_robust; final_fidelity=0.95, D=100.0)
solve!(qcp_mintime; max_iter=100)
```

## See Also

- [SmoothPulseProblem](@ref smooth-pulse) - Base problem template
- [SplinePulseProblem](@ref spline-pulse) - Spline-based base problem
- [MinimumTimeProblem](@ref minimum-time) - Time optimization wrapper
- [SamplingProblem](@ref sampling) - Robustness wrapper
