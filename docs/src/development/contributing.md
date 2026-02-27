# Contributing

We welcome contributions to Piccolo.jl! This document outlines the guidelines for contributing to the project. If you know what you want to see but are unsure of the best way to achieve it, [open an issue](https://github.com/harmoniqs/Piccolo.jl/issues) and start a discussion with the community.

## Development Setup

### Install Julia

[Juliaup](https://github.com/JuliaLang/juliaup) is an installer and version manager for Julia. After installing, run `julia` to obtain the Julia REPL.

### Julia Environments

Your project's environment is stored in `Project.toml`. You can interactively manage packages using the Julia REPL and package manager:

1. Start Julia in the project folder
2. Type `]` to enter the package manager
3. Type `activate .` to activate the environment
4. Run `instantiate` to install dependencies

### Development Installation

For development, clone the repository and use `dev` mode:

```julia
using Pkg
Pkg.dev("path/to/Piccolo.jl")
```

### Using Revise

[Revise.jl](https://timholy.github.io/Revise.jl/stable/) automatically reloads code changes during development:

```julia
using Revise
using Piccolo
```

Load `Revise` before any packages you intend to edit.

## Tips for Visual Studio Code

### Julia Extension

The [Julia extension](https://code.visualstudio.com/docs/languages/julia) provides excellent Julia support including notebooks, REPL integration, and debugging.

### Fonts

VS Code may not display all Julia characters correctly. Change the editor font family to `'JuliaMono'` for full Unicode support. You can create a VS Code settings profile for Julia at *File > Preferences > Profile*.

### Tests

Tests automatically populate in VS Code when working with Piccolo packages. Click the *Testing* sidebar icon to see and run tests. Sometimes you may need to restart the Julia kernel to see changes reflected in tests.

## Writing Tests

Tests are implemented using [TestItems.jl](https://www.julia-vscode.org/docs/stable/userguide/testitems/):

```julia
@testitem "X gate synthesis" begin
    H_drift = PAULIS[:Z]
    H_drives = [PAULIS[:X], PAULIS[:Y]]
    sys = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

    T, N = 10.0, 100
    times = collect(range(0, T, length=N))
    pulse = ZeroOrderPulse(0.1 * randn(2, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, GATES[:X])

    qcp = SmoothPulseProblem(qtraj, N)
    solve!(qcp; max_iter=100)

    @test fidelity(qcp) > 0.99
end
```

Tests should be included in the same file as the code they test. Individual tests populate in the Testing panel in VS Code and are run in CI on each PR.

## Building Documentation

Documentation is built using [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) with [Literate.jl](https://github.com/fredrikekre/Literate.jl) for executable examples.

### Build Locally

```bash
julia --project=docs docs/make.jl
```

### Live Development Server

For interactive documentation development:

```julia
julia --project=docs
```

```julia
using Revise, LiveServer, Piccolo
servedocs(
    literate_dir="docs/literate",
    skip_dirs=["docs/src/generated"],
    skip_files=["docs/src/index.md"]
)
```

Changes to documentation files are automatically reflected. To see source code changes (e.g., docstrings), restart the live server.

## Reporting Issues

Use the [GitHub issue templates](https://github.com/harmoniqs/Piccolo.jl/issues) for bug reports and feature requests.
