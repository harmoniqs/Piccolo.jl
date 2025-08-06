using Piccolo

pages = [
    "Home" => "index.md",
    # "Quickstart" => "generated/quickstart.md",
    # "Examples" => ["generated/multilevel_transmon.md"],
    "Features" => "generated/features.md",
    "Contribution Guide" => "contribution_guide.md",
    "Release Notes" => "release_notes.md",
]

include("utils.jl")

generate_docs(
    @__DIR__,
    "Piccolo",
    Piccolo,
    pages;
    format_kwargs = (
        canonical = "https://docs.harmoniqs.co/Piccolo.jl",
    ),
)
