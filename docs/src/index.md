# Piccolo.jl

<p align="center">
<img src="assets/piccolo_logo.svg" alt="logo" width="50%"/>
</p>

<div align="center">
  
| **Documentation** | **Build Status** | **Support** | **License** |
|:-----------------:|:----------------:|:-----------:|:---------:|
| [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kestrelquantum.github.io/QuantumCollocation.jl/dev/) | [![Build Status](https://github.com/aarontrowbridge/QuantumCollocation.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/kestrelquantum/QuantumCollocation.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Coverage](https://codecov.io/gh/kestrelquantum/QuantumCollocation.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/kestrelquantum/QuantumCollocation.jl)| [![Unitary Fund](https://img.shields.io/badge/Supported%20By-Unitary%20Fund-FFFF00.svg)](https://unitary.fund) | [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

</div>

> *"Technologies are ways of commandeering nature: the sky belongs to those who know how to fly; the sea belongs to those who know how to swim and navigate." &ndash; Simone de Beauvoir*

## Description
Piccolo.jl is a meta-package for quantum optimal control using the Pade Integrator COllocation (PICO) method.

For documentation please see the individual packages.

### JuliaCon 2023 Talk
To see an overview of the PICO method and a demo of how to use this package, check out our recorded talk at JuliaCon 2023 [here](https://www.youtube.com/watch?v=NBdck6UX0Tc).

### Direct Collocation for Quantum Optimal Control
To see a detailed description of the PICO method, check out our paper [here](https://arxiv.org/abs/2305.03261). It won 2nd best paper in the QTEM category at IEEE QCE 2023! 


## Usage

Just run
```julia
using Piccolo
```

and this package reexports the following packages

- [QuantumCollocation.jl](https://github.com/aarontrowbridge/QuantumCollocation.jl)
- [NamedTrajectories.jl](https://github.com/aarontrowbridge/NamedTrajectories.jl)
- [TrajectoryIndexingUtils.jl](https://github.com/aarontrowbridge/TrajectoryIndexingUtils.jl)

## Installation
This package is registered! To install enter the Julia REPL, type `]` to enter pkg mode, and then run 
```julia
pkg> add Piccolo
```

## Local Development

To develop locally, clone this repo and then instantiate the environment in the REPL by running first
```julia
pkg> activate .
```
and then, 
```julia
pkg> instantiate 
```
Both commands should be run in pkg mode, which is activated by typing `]` in the REPL.

To start julia with the current environment, run
```bash
julia --project
```
from the Piccolo.jl directory.

To run the scripts use, e.g.,
```bash
julia -t <num_threads> --project examples/three_qubit_swap/swap.jl
```
where `<num_threads>` is the number of threads you want to use as QuantumCollocation.jl takes advantage of multithreading.  The `--project` flag is necessary to make sure the correct environment is loaded.


## Docstrings

```@autodocs
Modules = [Piccolo]
```
