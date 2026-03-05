# [How-To Guides](@id guides-overview)

This section contains task-oriented guides for specific problems in quantum optimal control.

## Available Guides

| Guide | Description |
|-------|-------------|
| [Saving and Loading Pulses](@ref saving-loading) | Saving optimized pulses, warm-starting, and organizing results |
| [Leakage Suppression](@ref leakage-suppression) | Controlling population leakage in multilevel systems |
| [Global Variables](@ref global-variables) | Optimizing system parameters alongside controls |
| [Visualization](@ref visualization) | Plotting trajectories and analyzing results |
| [Custom Objectives](@ref custom-objectives) | Creating custom optimization objectives |

## Guide Format

Each guide follows a consistent format:
1. **Problem statement**: What are we trying to accomplish?
2. **Prerequisites**: What you need to know first
3. **Step-by-step solution**: How to solve the problem
4. **Complete example**: Full working code
5. **Variations**: Common modifications

## Quick Links by Task

### Saving and Reusing Results

- [Saving optimized pulses to disk](@ref saving-loading)
- [Warm-starting from a saved pulse](@ref saving-loading)
- [Organizing pulse files across a project](@ref saving-loading)

### Handling Multilevel Systems

- [Defining gates in subspaces with EmbeddedOperator](@ref leakage-suppression)
- [Adding leakage penalties and constraints](@ref leakage-suppression)

### Advanced Optimization

- [Optimizing system parameters](@ref global-variables)
- [Creating custom fidelity measures](@ref custom-objectives)

### Analysis and Visualization

- [Plotting control pulses](@ref visualization)
- [Visualizing state populations](@ref visualization)
- [Creating animations](@ref visualization)

## See Also

- [Quantum Systems](@ref systems-overview) - Platform-specific system templates
- [Problem Templates](@ref problem-templates-overview) - Main optimization API
- [Concepts](@ref concepts-overview) - Detailed type documentation
- [Tutorials](@ref tutorials-overview) - Step-by-step learning examples
