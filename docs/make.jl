using Piccolo
using PiccoloDocsTemplate
using DocumenterCopyButton: CopyButton

draft_mode_pages = [
# "first_gate.jl",
# "multilevel_transmon.jl",
# "quickstart.jl",
# "robust_control.jl",
# "state_transfer.jl",
# "quantum.jl",
# "visualizations.jl",
# "control.jl",
# "custom_objectives.jl",
# "global_variables.jl",
# "leakage_suppression.jl",
# "system_templates.jl",
# "visualization.jl",
# "bang_bang_pulse.jl",
# "smooth_pulse.jl",
# "spline_pulse.jl",
# "minimum_time.jl",
# "sampling.jl",
# "composition.jl",
# "concepts.jl",
# "installation.jl",
# "systems.jl",
# "trajectories.jl",
# "pulses.jl",
# "objectives.jl",
# "constraints.jl",
# "operators.jl",
# "isomorphisms.jl",
# "transmons.jl",
# "trapped_ions.jl",
# "rydberg_atoms.jl",
# "cat_qubits.jl",
# "silicon_spins.jl",
]

pages = [
    "Home" => "index.md",
    "Getting Started" => [
        "Installation" => "generated/getting-started/installation.md",
        "Quickstart" => "generated/quickstart.md",
    ],
    "Concepts" => [
        "Overview" => "concepts/index.md",
        "Trajectories" => "generated/concepts/trajectories.md",
        "Pulses" => "generated/concepts/pulses.md",
        "Objectives" => "generated/concepts/objectives.md",
        "Constraints" => "generated/concepts/constraints.md",
        "Operators" => "generated/concepts/operators.md",
        "Isomorphisms" => "generated/concepts/isomorphisms.md",
    ],
    "Tutorials" => [
        "Overview" => "tutorials/index.md",
        "Your First Gate" => "generated/first_gate.md",
        "Multilevel Transmon" => "generated/multilevel_transmon.md",
        "State Transfer" => "generated/state_transfer.md",
        "Robust Control" => "generated/robust_control.md",
        "Two-Qubit Gate Validation" => "generated/two_qubit_gate_validation.md",
    ],
    "Problem Templates" => [
        "Overview" => "problem-templates/index.md",
        "SmoothPulseProblem" => "generated/problem-templates/smooth_pulse.md",
        "BangBangPulseProblem" => "generated/problem-templates/bang_bang_pulse.md",
        "SplinePulseProblem" => "generated/problem-templates/spline_pulse.md",
        "MinimumTimeProblem" => "generated/problem-templates/minimum_time.md",
        "SamplingProblem" => "generated/problem-templates/sampling.md",
        "Composing Templates" => "generated/problem-templates/composition.md",
    ],
    "Quantum Systems" => [
        "Overview" => "systems/index.md",
        "Transmon Qubits" => "generated/systems/transmons.md",
        "Trapped Ions" => "generated/systems/trapped_ions.md",
        "Rydberg Atoms" => "generated/systems/rydberg_atoms.md",
        "Cat Qubits" => "generated/systems/cat_qubits.md",
        "Silicon Spins" => "generated/systems/silicon_spins.md",
    ],
    "How-To Guides" => [
        "Overview" => "guides/index.md",
        "Saving and Loading Pulses" => "generated/guides/saving_loading.md",
        "Leakage Suppression" => "generated/guides/leakage_suppression.md",
        "Global Variables" => "generated/guides/global_variables.md",
        "Visualization" => "generated/guides/visualization.md",
        "Wigner Bosonic Qubits" => "generated/guides/wigner_bosonic_qubits.md",
        "Custom Objectives" => "generated/guides/custom_objectives.md",
    ],
    "API Reference" => ["Overview" => "reference/index.md", "Library" => "lib.md"],
    "Development" => [
        "Contributing" => "development/contributing.md",
        "Release Notes" => "development/release-notes.md",
    ],
]

draft = get(ENV, "DOCS_DRAFT", "false") == "true"

generate_docs(
    @__DIR__,
    "Piccolo",
    Piccolo,
    pages;
    make_index = true,
    make_literate = true,
    make_assets = true,
    make_contributing = true,
    literate_draft_pages = draft_mode_pages,
    literate_kwargs = (execute = false,),
    format_kwargs = (
        canonical = "https://docs.harmoniqs.co/Piccolo.jl",
        # Headroom for lib.md, which is the whole auto-generated API reference and is therefore
        # inherently the largest page. It grows monotonically with every docstring added anywhere
        # in the package. Phase 1b's parametric-typing surface (the tag types, constrained
        # aliases, trait methods and `@problem_template` itself) pushed the generated HTML past
        # the previous 500 KiB ceiling — #260 and #262 each raised this independently to 700 KiB
        # and 1 MiB. Resolved to 1 MiB with headroom; `size_threshold_warn` still flags growth.
        #
        # The durable fix is to split lib.md per-module rather than keep raising this.
        size_threshold = 1024 * 2^10,      # 1 MiB hard limit for lib.md
        size_threshold_warn = 500 * 2^10,  # warn at the old ceiling
    ),
    mask_cached_solve = true,
    makedocs_kwargs = (draft = draft, plugins = [CopyButton()]),
)
