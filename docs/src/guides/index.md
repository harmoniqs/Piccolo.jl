# How-To Guides

This section contains task-oriented guides for specific problems in quantum optimal control.

## Available Guides

| Guide | Description |
|-------|-------------|
| [System Templates](@ref) | Using pre-built systems for transmons, ions, and atoms |
| [Leakage Suppression](@ref) | Controlling population leakage in multilevel systems |
| [Global Variables](@ref) | Optimizing system parameters alongside controls |
| [Visualization](@ref) | Plotting trajectories and analyzing results |
| [Custom Objectives](@ref) | Creating custom optimization objectives |

## Guide Format

Each guide follows a consistent format:
1. **Problem statement**: What are we trying to accomplish?
2. **Prerequisites**: What you need to know first
3. **Step-by-step solution**: How to solve the problem
4. **Complete example**: Full working code
5. **Variations**: Common modifications

## Quick Links by Task

### Setting Up Systems

- [Using TransmonSystem for superconducting qubits](@ref system-templates)
- [Using IonChainSystem for trapped ions](@ref system-templates)
- [Using RydbergChainSystem for neutral atoms](@ref system-templates)

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

- [Problem Templates](@ref) - Main optimization API
- [Concepts](@ref) - Detailed type documentation
- [Tutorials](@ref) - Step-by-step learning examples
