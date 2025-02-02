<!--```@raw html-->
<div align="center">

<a href="https://github.com/kestrelquantum/Piccolo.jl">
  <img src="assets/piccolo_logo.svg" alt="Piccolo.jl" width="25%"/>
</a> 

<table>
  <tr>
    <td align="center">
      <b>Documentation</b>
      <br>
      <a href="https://kestrelquantum.github.io/Piccolo.jl/stable/">
        <img src="https://img.shields.io/badge/docs-stable-blue.svg" alt="Stable"/>
      </a>
      <a href="https://kestrelquantum.github.io/Piccolo.jl/dev/">
        <img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="Dev"/>
      </a>
    </td>
    <td align="center">
      <b>Build Status</b>
      <br>
      <a href="https://github.com/kestrelquantum/Piccolo.jl/actions/workflows/CI.yml?query=branch%3Amain">
        <img src="https://github.com/kestrelquantum/Piccolo.jl/actions/workflows/CI.yml/badge.svg?branch=main" alt="Build Status"/>
      </a>
      <a href="https://codecov.io/gh/kestrelquantum/Piccolo.jl">
        <img src="https://codecov.io/gh/kestrelquantum/Piccolo.jl/branch/main/graph/badge.svg" alt="Coverage"/>
      </a>
    </td>
    <td align="center">
      <b>License</b>
      <br>
      <a href="https://opensource.org/licenses/MIT">
        <img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="MIT License"/>
      </a>
    </td>
    <td align="center">
      <b>Support</b>
      <br>
      <a href="https://unitary.fund">
        <img src="https://img.shields.io/badge/Supported%20By-Unitary%20Fund-FFFF00.svg" alt="Unitary Fund"/>
      </a>
    </td>
  </tr>
</table>

</div>
<!--```-->

## Description
**Piccolo.jl** is a meta-package for quantum optimal control using the Pade Integrator Collocation (Piccolo) method. This package reexports the following packages

- [QuantumCollocation.jl](https://github.com/kestrelquantum/QuantumCollocation.jl)
- [NamedTrajectories.jl](https://github.com/kestrelquantum/NamedTrajectories.jl)
- [TrajectoryIndexingUtils.jl](https://github.com/kestrelquantum/TrajectoryIndexingUtils.jl)
- [PiccoloQuantumObjects.jl](https://github.com/kestrelquantum/PiccoloQuantumObjects.jl)
- [PiccoloPlots.jl](https://github.com/kestrelquantum/PiccoloPlots.jl)

For documentation please see the individual packages.

## Usage

Just run
```Julia
using Piccolo
```

## Installation
This package is registered! To install enter the Julia REPL, type `]` to enter pkg mode, activate your environment `activate`, and then run 
```Julia
pkg> add Piccolo
```

-----

*"Technologies are ways of commandeering nature: the sky belongs to those who know how to fly; the sea belongs to those who know how to swim and navigate." – Simone de Beauvoir*
