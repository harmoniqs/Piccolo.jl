using Piccolo
using DocumenterInterLinks
using PiccoloQuantumObjects
using QuantumCollocation
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

links = InterLinks(
    "PiccoloQuantumObjects" => "https://docs.harmoniqs.co/PiccoloQuantumObjects/dev/",
    "QuantumCollocation" => "https://docs.harmoniqs.co/QuantumCollocation/dev/",
    "NamedTrajectories" => "https://docs.harmoniqs.co/NamedTrajectories/dev/",
)

makedocs_kwargs = (
    plugins=[links],
)

generate_docs(
    @__DIR__,
    "Piccolo",
    [Piccolo, QuantumCollocation, PiccoloQuantumObjects],
    pages;
    make_index = true,
    make_literate = true,
    make_assets = true,
    format_kwargs = (canonical = "https://docs.harmoniqs.co/Piccolo.jl",),
    makedocs_kwargs=makedocs_kwargs
)
