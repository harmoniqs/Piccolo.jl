# Developer Documentation

Much of these developer docs are based on and reference the SciML ecosystem and their standards.
As much as possible, we want to encourage an open and collaborative atmosphere while also
maintaining standards of quality and support through our tooling. We use this tooling and a set
of well-defined standards and practicies to support users and developers alike. As stated by SciML:

> For uniformity and clarity, the SciML Open-Source Software Organization has many
> well-defined rules and practices for its development. However, we stress one
> important principle:
> 
> **Do not be deterred from contributing if you think you do not know everything. No
> one knows everything. These rules and styles are designed for iterative contributions.
> Open pull requests and contribute what you can with what you know, and the maintainers
> will help you learn and do the rest!**
> 
> If you need any help contributing, please feel welcome by opening and PR/issue and we can invite you
> to our *#community* channel in slack!
> 
> We welcome everybody.

## Getting Started With Contributing to Piccolo

The following are a list of conventions that we use for development. There great resources online for Julia
development - if your are new. Give these a read/watch:

- [Developing Julia Packages](https://julialang.org/contribute/developing_package/)
- [More Tooling and Information](https://julialang.org/contribute/)

## Style:

We are currently using the Julia default style and formatting for the repository. 
In the future we will encourage use of the SciML Style Guide for Julia

[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

As described by SciML's contribution guidelines:

> This is a style guide for how to program in Julia for SciML contributions. It describes
> everything one needs to know, from preferred naming schemes of functions to fundamental
> dogmas for designing traits. We stress that this style guide is meant to be comprehensive
> for the sake of designing automatic formatters and teaching desired rules, but complete
> knowledge and adherence to the style guide is not required for contributions!

## COLPRAC: Contributor's Guide on Collaborative Practices for Community Packages

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

What are the rules for when PRs should be merged? What are the rules for whether to tag
a major, minor, or patch release? All of these development rules are defined in COLPRAC.

---

# Building Documentation

## Setup
This package uses a Documenter config that is shared with many of our other repositories. To build the docs, you will need to run the docs setup script to clone and pull down the utility. 
```
# first time only
./docs/get_docs_utils.sh   # or ./get_docs_utils.sh if cwd is in ./docs/
```

## Build
To build the docs pages:
```
julia --project=docs docs/make.jl
```

## Live-editing!
or editing the docs live:
```
julia --project=docs
> using LiveServer, Piccolo, Revise
> servedocs(literate_dir="docs/literate", skip_dirs=["docs/src/generated", "docs/src/assets/"], skip_files=["docs/src/index.md"])
```

> **Note:** `servedocs` needs to watch a subset of the files in the `docs/` folder. If it watches files that are generated on a docs build/re-build (by Literate.jl for example), `servedocs` will continuously try to re-serve the pages.
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

# Third-Party Libraries to Note

## PiccoloDocsTemplate.jl

Configuration for our docs building process is in `docs/make.jl` and relies upon an internal tool that wraps
[Documenter.jl](#documenterjl)'s makedocs and deploydocs calls: [PiccoloDocsTemplate.jl](https://github.com/harmoniqs/PiccoloDocsTemplate.jl). To bring this into the repository, before first build always run `./docs/get_docs_utils.sh` or `get_docs_utils.sh`.

## Literate.jl
[Literate.jl](https://fredrikekre.github.io/Literate.jl/v2/) is a tool that allows us to write
nearly all of our docs using [Literate Programming](https://en.wikipedia.org/wiki/Literate_programming).
Its source are julia scripts, which are parsed and then converted into markdown files with the relevant
example blocks, which are then executed when building the docs to provide example output by [Documenter.jl](#documenterjl)

## Documenter.jl

[Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) is the documentation generation
library that our organization uses, and thus its documentation is the documentation
of the documentation.

## JuliaFormatter.jl

[JuliaFormatter.jl](https://github.com/domluna/JuliaFormatter.jl) is the formatter used by the
our organization to enforce the a consistent style.

To run JuliaFormatter, do:

```julia
import JuliaFormatter, DevedPackage
JuliaFormatter.format(pkgdir(DevedPackage))
```

which will reformat the code according to the SciML Style.

## GitHub Actions Continuous Integrations

The Harmoniqs Organization uses continuous integration testing to always ensure tests are passing when merging
pull requests. The organization uses the GitHub Actions supplied by [Julia Actions](https://github.com/julia-actions)
to accomplish this. Common continuous integration scripts are:

  - CI.yml, the standard CI script
  - docs.yml, the full documentation build and deployment script - generates documentation via [Literate](#literatejl) and builds and runs examples with [Documenter.jl](#documenterjl).
  - CompatHelper.yml, checks for updated versions of dependencies and makes pull requests automatically when new versions
  are available.
  - TagBot.yml, triggers after new versions of package are registered on the [General registry](https://github.com/JuliaRegistries/General) to tag, release, (and also build the docs with the new tag!)

## CompatHelper

[CompatHelper](https://github.com/JuliaRegistries/CompatHelper.jl) is used to automatically create pull requests whenever
a dependent package is upper bounded. The results of CompatHelper PRs should be checked to ensure that the latest version
of the dependencies are grabbed for the test process. After successful CompatHelper PRs, i.e. if the increase of the upper
bound did not cause a break to the tests, a new version tag should follow. It is set up by adding the CompatHelper.yml GitHub action.

## TagBot

[TagBot](https://github.com/JuliaRegistries/TagBot) automatically creates tags in the GitHub repository whenever a package
is registered to the Julia General repository. It is set up by adding the TagBot.yml GitHub action.
