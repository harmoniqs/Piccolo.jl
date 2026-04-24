# Contributing to Piccolo.jl

We welcome contributions! Much of our developer practices are based on and reference the [SciML ecosystem standards](https://github.com/SciML/SciMLDocs/blob/main/docs/src/highlevels/developer_documentation.md?plain=1). We use tooling and well-defined practices to support users and developers alike.

> **Do not be deterred from contributing if you think you do not know everything. No one knows everything. These rules and styles are designed for iterative contributions. Open pull requests and contribute what you can with what you know, and the maintainers will help you learn and do the rest!**
> — [SciML](https://sciml.ai)

If you need help contributing, open a PR or issue and we can invite you to our *#community* channel in Slack.

## Development Setup

### Install Julia

[Juliaup](https://github.com/JuliaLang/juliaup) is the recommended installer and version manager for Julia. After installing, run `julia` to open the REPL.

### Clone and Activate

Clone the repository and activate the project environment:

```
] activate .
] instantiate
```

Or equivalently:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

For development mode (tracking local changes):

```julia
using Pkg
Pkg.dev("path/to/Piccolo.jl")
```

### Using Revise

[Revise.jl](https://timholy.github.io/Revise.jl/stable/) automatically reloads code changes during development. Load it *before* any packages you intend to edit:

```julia
using Revise
using Piccolo
```

## Code Style & Formatting

We currently use the Julia default style. In the future we will adopt the [SciML Style Guide](https://github.com/SciML/SciMLStyle):

[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

### Formatting with JuliaFormatter

[JuliaFormatter.jl](https://github.com/domluna/JuliaFormatter.jl) enforces consistent formatting. The CI formatter check runs on every pull request and will **fail** if any files are not properly formatted.

**Run locally before pushing:**

```julia
using JuliaFormatter
format(".")
```

Or from the command line:

```bash
julia -e 'using JuliaFormatter; format(".", verbose=true)'
```

**Maintainers** can also trigger the formatter check on any branch via [workflow dispatch](https://github.com/harmoniqs/Piccolo.jl/actions/workflows/CI.yml) in GitHub Actions.

## COLPRAC

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

We follow [COLPRAC](https://github.com/SciML/ColPrac) for collaborative practices — it defines rules for when PRs should be merged and whether to tag a major, minor, or patch release.

## Writing Tests

Tests use [TestItems.jl](https://www.julia-vscode.org/docs/stable/userguide/testitems/) and should be included in the same file as the code they test:

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

Individual test items populate in the VS Code Testing panel and run in CI on each pull request.

## Building Documentation

Documentation is built using [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) with [Literate.jl](https://github.com/fredrikekre/Literate.jl) for executable examples. The build configuration is managed by [PiccoloDocsTemplate.jl](https://github.com/harmoniqs/PiccoloDocsTemplate.jl), an internal tool that wraps Documenter's `makedocs` and `deploydocs` calls.

### First-Time Setup

Fetch the documentation template:

```bash
./docs/get_docs_utils.sh
```

### Build

```bash
julia --project=docs docs/make.jl
```

### Live Development Server

For interactive documentation editing with hot-reload:

```julia
julia --project=docs
```

```julia
using Revise, LiveServer, Piccolo
servedocs(
    literate_dir="docs/literate",
    skip_dirs=["docs/src/generated", "docs/src/assets/"],
    skip_files=["docs/src/index.md"]
)
```

> **Note:** `servedocs` watches files in the docs folder. Generated files (from Literate.jl) must be excluded via `skip_dirs` and `skip_files` to prevent continuous rebuild loops.

## VS Code Tips

The [Julia VS Code extension](https://code.visualstudio.com/docs/languages/julia) provides REPL integration, notebooks, and debugging support.

- **Fonts:** Set the editor font to `JuliaMono` for full Unicode support (*File > Preferences > Profile* for a Julia-specific settings profile).
- **Tests:** Test items populate in the Testing sidebar. Restart the Julia kernel to pick up new changes.

## Reporting Issues

Use [GitHub Issues](https://github.com/harmoniqs/Piccolo.jl/issues) for bug reports and feature requests.

## CI Workflows

The repository uses [Julia Actions](https://github.com/julia-actions) for GitHub Actions CI:

| Workflow | Purpose |
|----------|---------|
| `CI.yml` | Runs tests across Julia versions and checks formatting |
| `docs.yml` | Builds and deploys documentation |
| `CompatHelper.yml` | Opens PRs when dependency upper bounds can be raised |
| `TagBot.yml` | Creates GitHub releases when new versions are registered |

## Additional Resources

- [Developing Julia Packages](https://julialang.org/contribute/developing_package/)
- [Julia Contributing Guide](https://julialang.org/contribute/)
