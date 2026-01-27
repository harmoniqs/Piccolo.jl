# [Objectives](@id objectives-concept)

Objectives define what the optimization minimizes. Piccolo.jl provides objectives for fidelity, regularization, and leakage suppression.

## Overview

The total objective is a weighted sum:

```math
J = Q \cdot J_{\text{fidelity}} + R_u \cdot J_u + R_{du} \cdot J_{du} + R_{ddu} \cdot J_{ddu} + \cdots
```

Where:
- `J_fidelity`: Infidelity (1 - F)
- `J_u`, `J_du`, `J_ddu`: Control regularization terms

## Fidelity Objectives

### UnitaryInfidelityObjective

For unitary gate synthesis:

```julia
obj = UnitaryInfidelityObjective(:Ũ⃗, U_goal; Q=100.0)
```

**Fidelity metric:**
```math
F = \frac{1}{d^2} \left| \text{Tr}(U_{\text{goal}}^\dagger U) \right|^2
```

Where `d` is the Hilbert space dimension.

### KetInfidelityObjective

For single state transfer:

```julia
obj = KetInfidelityObjective(:ψ̃, ψ_goal; Q=100.0)
```

**Fidelity metric:**
```math
F = \left| \langle \psi_{\text{goal}} | \psi \rangle \right|^2
```

### CoherentKetInfidelityObjective

For multiple state transfers with phase coherence (used by `MultiKetTrajectory`):

```julia
obj = CoherentKetInfidelityObjective(
    [:ψ̃1, :ψ̃2],
    [ψ_goal1, ψ_goal2];
    Q=100.0
)
```

This ensures both:
1. Each state reaches its target
2. Relative phases between states are preserved

### UnitaryFreePhaseInfidelityObjective

When global phase doesn't matter:

```julia
obj = UnitaryFreePhaseInfidelityObjective(:Ũ⃗, U_goal; Q=100.0)
```

Optimizes over the global phase to find the best match.

## Regularization Objectives

Regularization penalizes large or rapidly-varying controls.

### QuadraticRegularizer

The standard regularization form:

```julia
# Penalize control magnitude
reg_u = QuadraticRegularizer(:u, traj; R=1e-2)

# Penalize control derivatives
reg_du = QuadraticRegularizer(:du, traj; R=1e-2)
reg_ddu = QuadraticRegularizer(:ddu, traj; R=1e-2)
```

**Mathematical form:**
```math
J_x = \sum_k \| x_k \|^2
```

### Per-Variable Weights

Apply different regularization to different drives:

```julia
# Vector of weights (one per drive)
R_u = [1e-3, 1e-2]  # Less regularization on drive 1
```

### Why Regularize?

1. **Smoothness**: Derivative regularization encourages smooth pulses
2. **Robustness**: Prevents exploiting numerical precision
3. **Hardware-friendliness**: Bounded, smooth controls are easier to implement
4. **Convergence**: Regularization improves optimization landscape

## Leakage Objectives

For multilevel systems, leakage to non-computational states can be penalized.

### LeakageObjective

```julia
# For EmbeddedOperator with defined subspace
op = EmbeddedOperator(:X, sys)  # X gate in computational subspace
obj = LeakageObjective(:Ũ⃗, op; Q=10.0)
```

**Mathematical form:**
Penalizes population outside the computational subspace at the final time.

### Via PiccoloOptions

The easier way is through `PiccoloOptions`:

```julia
opts = PiccoloOptions(
    leakage_constraint=true,
    leakage_constraint_value=1e-3,
    leakage_cost=10.0
)

qcp = SmoothPulseProblem(qtraj, N; piccolo_options=opts)
```

## Sensitivity Objectives

For analyzing robustness to parameter variations.

### UnitarySensitivityObjective

Penalizes sensitivity of fidelity to system parameters:

```julia
obj = UnitarySensitivityObjective(:Ũ⃗, U_goal, parameter_variations; Q=1.0)
```

This is related to but different from `SamplingProblem`:
- `SamplingProblem`: Optimize for multiple systems simultaneously
- `UnitarySensitivityObjective`: Penalize gradient of fidelity w.r.t. parameters

## Using Objectives in Problem Templates

Problem templates automatically set up objectives. You typically don't create them manually.

### Automatic Setup

```julia
# SmoothPulseProblem automatically creates:
# - UnitaryInfidelityObjective (for UnitaryTrajectory)
# - QuadraticRegularizer for :u, :du, :ddu
# - LeakageObjective if piccolo_options.leakage_cost > 0

qcp = SmoothPulseProblem(
    qtraj, N;
    Q=100.0,      # Fidelity weight
    R=1e-2,       # Base regularization
    R_u=1e-3,     # Control regularization (overrides R)
    R_du=1e-2,    # First derivative regularization
    R_ddu=1e-2    # Second derivative regularization
)
```

### Trajectory-Dependent Objectives

The objective type is chosen based on the trajectory:

| Trajectory Type | Default Objective |
|-----------------|-------------------|
| `UnitaryTrajectory` | `UnitaryInfidelityObjective` |
| `KetTrajectory` | `KetInfidelityObjective` |
| `MultiKetTrajectory` | `CoherentKetInfidelityObjective` |
| `DensityTrajectory` | None (manual setup) |

### Custom Objectives

For advanced use, you can add custom objectives:

```julia
# Access the underlying DirectTrajOpt problem
prob = qcp.prob

# Add custom objective
custom_obj = MyCustomObjective(...)
add_objective!(prob, custom_obj)
```

## Objective Weights

### The Q Parameter

`Q` weights the fidelity objective:
- **Higher Q** (e.g., 1000): Prioritize fidelity over smoothness
- **Lower Q** (e.g., 10): Allow more flexibility in controls

### The R Parameters

`R`, `R_u`, `R_du`, `R_ddu` weight regularization:
- **Higher R**: Smoother, smaller controls
- **Lower R**: More aggressive controls allowed

### Balancing Trade-offs

```julia
# High fidelity, some smoothness
qcp = SmoothPulseProblem(qtraj, N; Q=1000.0, R=1e-2)

# Very smooth, may sacrifice fidelity
qcp = SmoothPulseProblem(qtraj, N; Q=100.0, R=1.0)

# Minimal regularization (may overfit)
qcp = SmoothPulseProblem(qtraj, N; Q=100.0, R=1e-6)
```

### Typical Starting Values

| Parameter | Typical Range | Starting Point |
|-----------|---------------|----------------|
| `Q` | 10 - 10000 | 100 |
| `R` | 1e-6 - 1.0 | 1e-2 |
| `R_u` | same as R | R |
| `R_du` | same as R | R |
| `R_ddu` | same as R | R |

## Best Practices

### 1. Start with Defaults

Problem templates have sensible defaults. Start there:

```julia
qcp = SmoothPulseProblem(qtraj, N)  # Uses Q=100, R=1e-2
```

### 2. Tune Q First

If fidelity is too low, increase Q:

```julia
qcp = SmoothPulseProblem(qtraj, N; Q=1000.0)
```

### 3. Tune R if Controls are Problematic

If controls are too noisy or large:

```julia
qcp = SmoothPulseProblem(qtraj, N; R=0.1)
```

### 4. Use Per-Derivative Tuning for Fine Control

```julia
qcp = SmoothPulseProblem(
    qtraj, N;
    R_u=1e-4,    # Allow larger control values
    R_du=1e-2,   # Penalize jumps moderately
    R_ddu=0.1    # Strongly penalize acceleration
)
```

## See Also

- [Constraints](@ref) - Hard constraints on solutions
- [Problem Templates](@ref) - How objectives are used
- [SmoothPulseProblem](@ref smooth-pulse) - Parameter reference
