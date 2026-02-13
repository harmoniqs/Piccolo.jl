using Piccolo
using PiccoloDocsTemplate

draft_mode_pages = [
# "first_gate.jl",
# "multilevel_transmon.jl",
# "quickstart.jl",
# "robust_control.jl",
# "state_transfer.jl",
# # "quantum.jl",
# "visualizations.jl",
# "control.jl",
# "custom_objectives.jl",
# "global_variables.jl",
# "leakage_suppression.jl",
# "system_templates.jl",
# "visualization.jl",
# "smooth_pulse.jl",
# "spline_pulse.jl",
# "minimum_time.jl",
# "sampling.jl",
# "composition.jl",
# "concepts.jl",
# "installation.jl",
# # "systems.jl",
# "trajectories.jl",
# "pulses.jl",
# "objectives.jl",
# "constraints.jl",
# "operators.jl",
# "isomorphisms.jl",
]

pages = [
    "Home" => "index.md",
    "Getting Started" => [
        "Installation" => "generated/getting-started/installation.md",
        "Quickstart" => "generated/quickstart.md",
        "Core Concepts" => "generated/getting-started/concepts.md",
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
        "SmoothPulseProblem" => "generated/problem-templates/smooth_pulse.md",
        "SplinePulseProblem" => "generated/problem-templates/spline_pulse.md",
        "MinimumTimeProblem" => "generated/problem-templates/minimum_time.md",
        "SamplingProblem" => "generated/problem-templates/sampling.md",
        "Composing Templates" => "generated/problem-templates/composition.md",
    ],
    "Concepts" => [
        "Overview" => "concepts/index.md",
        "Quantum Systems" => "generated/concepts/systems.md",
        "Trajectories" => "generated/concepts/trajectories.md",
        "Pulses" => "generated/concepts/pulses.md",
        "Objectives" => "generated/concepts/objectives.md",
        "Constraints" => "generated/concepts/constraints.md",
        "Operators" => "generated/concepts/operators.md",
        "Isomorphisms" => "generated/concepts/isomorphisms.md",
    ],
    "How-To Guides" => [
        "Overview" => "guides/index.md",
        "System Templates" => "generated/guides/system_templates.md",
        "Leakage Suppression" => "generated/guides/leakage_suppression.md",
        "Global Variables" => "generated/guides/global_variables.md",
        "Visualization" => "generated/guides/visualization.md",
        "Custom Objectives" => "generated/guides/custom_objectives.md",
    ],
    "API Reference" => ["Overview" => "reference/index.md", "Library" => "lib.md"],
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
    literate_draft_pages = draft_mode_pages,
    literate_kwargs = (execute = false,),
    format_kwargs = (
        canonical = "https://docs.harmoniqs.co/Piccolo.jl",
        size_threshold = 400 * 2^10,  # 400 KiB for lib.md
    ),
)
