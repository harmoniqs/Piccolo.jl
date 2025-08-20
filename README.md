<!--```@raw html-->
<div align="center">
  <a href="https://github.com/harmoniqs/Piccolo.jl">
    <img src="assets/piccolo_logo_high_contrast.svg" alt="Piccolo.jl" width="45%"/>
  </a> 
</div>

<div align="center">
  <table>
    <tr>
      <td align="center">
        <b>Documentation</b>
        <br>
        <a href="https://docs.harmoniqs.co/Piccolo/dev/">
          <img src="https://img.shields.io/badge/docs-stable-blue.svg" alt="Stable"/>
        </a>
        <a href="https://docs.harmoniqs.co/Piccolo/dev/">
          <img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="Dev"/>
        </a>
      </td>
      <td align="center">
        <b>Build Status</b>
        <br>
        <a href="https://github.com/harmoniqs/Piccolo.jl/actions/workflows/CI.yml?query=branch%3Amain">
          <img src="https://github.com/harmoniqs/Piccolo.jl/actions/workflows/CI.yml/badge.svg?branch=main" alt="Build Status"/>
        </a>
        <a href="https://codecov.io/gh/harmoniqs/Piccolo.jl">
          <img src="https://codecov.io/gh/harmoniqs/Piccolo.jl/branch/main/graph/badge.svg" alt="Coverage"/>
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

### Description
**Piccolo.jl** is a meta-package for quantum optimal control using the Pade Integrator Collocation (Piccolo) method. This package reexports the following packages

- [DirectTrajOpt.jl](https://github.com/harmoniqs/DirectTrajOpt.jl)
- [NamedTrajectories.jl](https://github.com/harmoniqs/NamedTrajectories.jl)
- [TrajectoryIndexingUtils.jl](https://github.com/harmoniqs/TrajectoryIndexingUtils.jl)

For documentation please see the individual packages.

### Usage

Just run
```Julia
using Piccolo
```

### Installation
This package is registered! To install enter the Julia REPL, type `]` to enter pkg mode, activate your environment `activate`, and then run 
```Julia
pkg> add Piccolo
```

### Building Documentation
This package uses a Documenter config that is shared with many of our other repositories. To build the docs, you will need to run the docs setup script to clone and pull down the utility. 
```
# first time only
./docs/get_docs_utils.sh   # or ./get_docs_utils.sh if cwd is in ./docs/
```

To build the docs pages:
```
julia --project=docs docs/make.jl
```

or editing the docs live:
```
julia --project=docs
> using LiveServer, Piccolo, Revise
> servedocs(literate_dir="docs/literate", skip_dirs=["docs/src/generated", "docs/src/assets/"], skip_files=["docs/src/index.md"])
```

> **Note:** `servedocs` needs to watch a subset of the files in the `docs/` folder. If it watches files that are generated on a docs build/re-build, `servedocs` will continuously try to re-serve the pages.
> 
> To prevent this, ensure all generated files are included in the skip dirs or skip files args for `servedocs`.

For example, if we forget index.md like so:
```
julia --project=docs
> using LiveServer, Piccolo, Revise
> servedocs(literate_dir="docs/literate", skip_dirs=["docs/src/generated", "docs/src/assets/"])
```
it will not build and serve.

-----

<!--## Star History-->

<!--[![Star History Chart](https://api.star-history.com/svg?repos=harmoniqs/piccolo.jl,harmoniqs/namedtrajectories.jl,harmoniqs/quantumcollocation.jl,harmoniqs/directtrajopt.jl&type=Date)](https://www.star-history.com/#harmoniqs/piccolo.jl&harmoniqs/namedtrajectories.jl&harmoniqs/quantumcollocation.jl&harmoniqs/directtrajopt.jl&Date)-->

*"Technologies are ways of commandeering nature: the sky belongs to those who know how to fly; the sea belongs to those who know how to swim and navigate." – Simone de Beauvoir*
