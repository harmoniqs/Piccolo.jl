using Piccolo
using PiccoloDocsTemplate

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
]

generate_docs(
    @__DIR__,
    "Piccolo",
    Piccolo,
    pages;
    make_index = true,
    make_literate = true,
    make_assets = true,
    format_kwargs = (canonical = "https://docs.harmoniqs.co/Piccolo.jl",),
)
