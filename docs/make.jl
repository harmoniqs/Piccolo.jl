using Piccolo
using PiccoloDocsTemplate

pages = [
    "Home" => "index.md",
    # "Quickstart" => "generated/quickstart.md",
    "Features" => "generated/features.md",
    "Other Considerations" => "generated/other.md",
    # "Examples" => ["generated/multilevel_transmon.md"],
    "Contribution Guide" => "contribution_guide.md",
    "Release Notes" => "release_notes.md",
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
