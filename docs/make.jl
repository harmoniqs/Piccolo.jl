using Documenter
using Literate
using Piccolo

# Process Literate files
literate_dir = joinpath(@__DIR__, "literate")
output_dir = joinpath(@__DIR__, "src", "generated")

# Ensure output directory exists
mkpath(output_dir)

# Files that are fast enough to execute during doc builds
executable_files = Set(["multilevel_transmon.jl"])

# Process all .jl files in literate/
for file in readdir(literate_dir)
    if endswith(file, ".jl")
        # Use DocumenterFlavor for fast files, CommonMarkFlavor for slow ones
        if file in executable_files
            Literate.markdown(
                joinpath(literate_dir, file),
                output_dir;
                documenter = true,
                credit = false
            )
        else
            # CommonMark flavor creates fenced code blocks that don't execute
            Literate.markdown(
                joinpath(literate_dir, file),
                output_dir;
                flavor = Literate.CommonMarkFlavor(),
                credit = false
            )
        end
    end
end

# Build documentation
makedocs(
    sitename = "Piccolo.jl",
    modules = [Piccolo],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://docs.harmoniqs.co/Piccolo.jl",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => [
            "Installation" => "getting-started/installation.md",
            "Quickstart" => "generated/quickstart.md",
            "Core Concepts" => "getting-started/concepts.md",
        ],
        "Tutorials" => [
            "Overview" => "tutorials/index.md",
            "Your First Gate" => "generated/first_gate.md",
            "Multilevel Transmon" => "generated/multilevel_transmon.md",
            "State Transfer" => "generated/state_transfer.md",
            "Robust Control" => "generated/robust_control.md",
        ],
        "Problem Templates" => [
            "Overview" => "problem-templates/index.md",
            "SmoothPulseProblem" => "problem-templates/smooth-pulse.md",
            "SplinePulseProblem" => "problem-templates/spline-pulse.md",
            "MinimumTimeProblem" => "problem-templates/minimum-time.md",
            "SamplingProblem" => "problem-templates/sampling.md",
            "Composing Templates" => "problem-templates/composition.md",
        ],
        "Concepts" => [
            "Overview" => "concepts/index.md",
            "Quantum Systems" => "concepts/systems.md",
            "Trajectories" => "concepts/trajectories.md",
            "Pulses" => "concepts/pulses.md",
            "Objectives" => "concepts/objectives.md",
            "Constraints" => "concepts/constraints.md",
            "Operators" => "concepts/operators.md",
            "Isomorphisms" => "concepts/isomorphisms.md",
        ],
        "How-To Guides" => [
            "Overview" => "guides/index.md",
            "System Templates" => "guides/system-templates.md",
            "Leakage Suppression" => "guides/leakage-suppression.md",
            "Global Variables" => "guides/global-variables.md",
            "Visualization" => "guides/visualization.md",
            "Custom Objectives" => "guides/custom-objectives.md",
        ],
        "API Reference" => [
            "Overview" => "reference/index.md",
            "Quantum Module" => "reference/quantum.md",
            "Control Module" => "reference/control.md",
            "Visualizations" => "reference/visualizations.md",
        ],
        "Development" => [
            "Contributing" => "development/contributing.md",
            "Release Notes" => "development/release-notes.md",
        ],
    ],
    checkdocs = :exports,
    warnonly = [:missing_docs, :cross_references],
)

deploydocs(
    repo = "github.com/harmoniqs/Piccolo.jl.git",
    devbranch = "main",
    push_preview = true,
)
