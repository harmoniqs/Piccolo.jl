module Specs

using Reexport
using ..Quantum
using ..Control
using TOML
using JSON3
using SHA
using TestItems

include("spec_structs.jl")
@reexport using .SpecStructs

include("registries.jl")
@reexport using .SpecRegistries

# includes added in later tasks:
# include("errors.jl");       @reexport using .SpecErrors
# include("hashes.jl")
# include("parse.jl")
# include("materialize.jl")
# include("rollout_kind.jl")
# include("run.jl")
# include("schema.jl")

# Piccolo self-registers its systems/templates/etc. at load time (idempotent).
register_all!()

end
