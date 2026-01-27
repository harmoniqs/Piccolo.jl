# [Installation](@id installation)

## Requirements

- Julia 1.9 or later
- A working internet connection (for package installation)

## Installing Piccolo.jl

Piccolo.jl is registered in the Julia General registry. Install it using the Julia package manager:

```julia
using Pkg
Pkg.add("Piccolo")
```

Or from the Julia REPL, press `]` to enter package mode and run:

```
pkg> add Piccolo
```

## Verifying Installation

Test your installation by running:

```julia
using Piccolo

# Check that core types are available
H_drift = PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
sys = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

println("Piccolo.jl installed successfully!")
println("System has $(sys.levels) levels and $(sys.n_drives) drives")
```

## Optional: Visualization Support

For plotting, you'll need a Makie backend. CairoMakie is recommended:

```julia
using Pkg
Pkg.add("CairoMakie")
```

Then in your code:

```julia
using Piccolo
using CairoMakie
```

## Development Installation

For contributing to Piccolo.jl or using the latest development version:

```julia
using Pkg
Pkg.develop(url="https://github.com/harmoniqs/Piccolo.jl")
```

Or clone the repository and use `dev` mode:

```bash
git clone https://github.com/harmoniqs/Piccolo.jl.git
cd Piccolo.jl
julia --project=.
```

```julia
using Pkg
Pkg.instantiate()
```

## Troubleshooting

### Precompilation Takes Long

First-time precompilation can take several minutes due to dependencies. This is normal and only happens once.

### Missing Dependencies

If you encounter missing dependency errors, try:

```julia
using Pkg
Pkg.instantiate()
Pkg.precompile()
```

### Version Conflicts

If you have version conflicts with other packages:

```julia
using Pkg
Pkg.update()
```

Or create a fresh environment:

```julia
using Pkg
Pkg.activate("my_project")
Pkg.add("Piccolo")
```

## Next Steps

Once installed, continue to the [Quickstart](@ref quickstart) guide to run your first optimization.
