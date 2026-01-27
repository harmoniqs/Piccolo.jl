# [Tutorials](@id tutorials-overview)

Step-by-step tutorials for learning Piccolo.jl. Each tutorial builds on previous ones, introducing new concepts progressively.

## Learning Path

### For Beginners

If you're new to quantum optimal control or Piccolo.jl:

1. **[Your First Gate](@ref first-gate-tutorial)** - Complete walkthrough of X gate synthesis
2. **[State Transfer](@ref state-transfer-tutorial)** - Learn about KetTrajectory for state preparation
3. **[Multilevel Transmon](@ref multilevel-transmon-tutorial)** - Work with realistic multilevel systems

### For Experienced Users

If you're familiar with quantum control concepts:

1. **[Multilevel Transmon](@ref multilevel-transmon-tutorial)** - EmbeddedOperator and leakage suppression
2. **[Robust Control](@ref robust-control-tutorial)** - SamplingProblem for parameter robustness

## Tutorial Overview

| Tutorial | Duration | Topics |
|----------|----------|--------|
| [Your First Gate](@ref first-gate-tutorial) | 15 min | Systems, pulses, trajectories, solving, analysis |
| [State Transfer](@ref state-transfer-tutorial) | 15 min | KetTrajectory, state preparation, coherent gates |
| [Multilevel Transmon](@ref multilevel-transmon-tutorial) | 20 min | TransmonSystem, EmbeddedOperator, leakage |
| [Robust Control](@ref robust-control-tutorial) | 20 min | SamplingProblem, parameter uncertainty |

## Prerequisites

Before starting the tutorials, ensure you have:

1. **Julia 1.9+** installed
2. **Piccolo.jl** installed: `using Pkg; Pkg.add("Piccolo")`
3. **CairoMakie** for plotting: `using Pkg; Pkg.add("CairoMakie")`

Test your setup:

```julia
using Piccolo
using CairoMakie

println("Setup complete!")
```

## What You'll Learn

### Your First Gate

- Defining quantum systems with `QuantumSystem`
- Creating control pulses with `ZeroOrderPulse`
- Setting optimization goals with `UnitaryTrajectory`
- Solving with `SmoothPulseProblem`
- Analyzing results with `fidelity()` and plotting

### State Transfer

- Using `KetTrajectory` for state-to-state transfer
- Working with `MultiKetTrajectory` for multiple state mappings
- Understanding coherent fidelity for phase-sensitive gates

### Multilevel Transmon

- Using `TransmonSystem` for realistic superconducting qubits
- Defining subspace gates with `EmbeddedOperator`
- Suppressing leakage to higher energy levels
- Using `PiccoloOptions` for leakage handling

### Robust Control

- Creating system variants for parameter uncertainty
- Using `SamplingProblem` for robust optimization
- Evaluating robustness across parameter ranges
- Combining robustness with time optimization

## Running the Tutorials

Each tutorial is a complete, runnable Julia script. You can:

1. **Read online**: Follow along in the documentation
2. **Run locally**: Download the `.jl` file from the `docs/literate/` folder
3. **Copy-paste**: Copy code blocks into your Julia REPL

## Next Steps

After completing the tutorials:

- Explore [Problem Templates](@ref problem-templates-overview) for the full API
- Read [Concepts](@ref concepts-overview) for detailed documentation
- Check [How-To Guides](@ref guides-overview) for specific tasks
