module Specs

using Reexport
using ..Quantum
using ..Control
using TOML
using JSON3
using SHA
using Printf
using TestItems

include("spec_structs.jl")
@reexport using .SpecStructs

include("registries.jl")
@reexport using .SpecRegistries

include("errors.jl")
@reexport using .SpecErrors

include("parse.jl")
include("hashes.jl")
include("materialize.jl")
include("rollout_kind.jl")
include("run.jl")
include("schema.jl")

# Piccolo self-registers its systems/templates/etc. at load time (idempotent).
register_all!()

end
