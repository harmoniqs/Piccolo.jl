using Piccolo

pages = [
    "Home" => "index.md",
    "Quickstart" => "generated/quickstart.md",
    "Examples" => ["generated/multilevel_transmon.md"],
    "Contribution Guide" => "contribution_guide.md",
    "Release Notes" => "release_notes.md",
]

# Check if utils.jl exists and warn if not found
utils_path = joinpath(@__DIR__, "utils.jl")
if !isfile(utils_path)
    error("docs/utils.jl is required but not found. Please run ./get_docs_utils.sh")
end

include("utils.jl")

generate_docs(
    @__DIR__,
    "Piccolo",
    Piccolo,
    pages;
    format_kwargs = (
        canonical = "https://docs.harmoniqs.co/Piccolo.jl",
    ),
    makedocs_kwargs = (
        draft = true,
    ),
)
