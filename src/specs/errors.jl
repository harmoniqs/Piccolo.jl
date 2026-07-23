module SpecErrors

using JSON3

export SpecError, SpecValidationError, emit_error_json

"""
    SpecError(path, msg, [got, [allowed]])

A single structured validation error. `path` is a dotted field path
(e.g. `"problem.Q"`), `msg` is a human-readable message, `got` is the offending
value (optional), and `allowed` is the set of acceptable values (optional, for
enum-style errors). Convenience constructors default the trailing fields to
`nothing`.
"""
struct SpecError
    path::String
    msg::String
    got::Any
    allowed::Union{Nothing,Vector}
end
SpecError(path, msg) = SpecError(path, msg, nothing, nothing)
SpecError(path, msg, got) = SpecError(path, msg, got, nothing)

"""
    SpecValidationError(errors::Vector{SpecError}) <: Exception

Raised when a spec fails validation. Carries **all** collected errors (the
parser does not fail on the first), and renders as the machine-readable
`emit_error_json` payload. The CLI entrypoint maps this to exit code 2.
"""
struct SpecValidationError <: Exception
    errors::Vector{SpecError}
end
Base.showerror(io::IO, e::SpecValidationError) =
    print(io, "SpecValidationError:\n", emit_error_json(e.errors))

"""
    emit_error_json(errs::Vector{SpecError}) -> String

Serialize a list of [`SpecError`](@ref)s to the `{"ok":false,"errors":[…]}`
JSON contract. `got`/`allowed` are omitted when unset.
"""
function emit_error_json(errs::Vector{SpecError})
    payload = Dict("ok" => false, "errors" => map(errs) do e
        d = Dict{String,Any}("path" => e.path, "msg" => e.msg)
        e.got === nothing || (d["got"] = e.got)
        e.allowed === nothing || (d["allowed"] = e.allowed)
        d
    end)
    JSON3.write(payload)
end

end

@testitem "structured errors carry field paths and emit JSON" begin
    using Piccolo.Specs
    errs = [Specs.SpecError("problem.Q", "must be > 0", -5.0, nothing),
        Specs.SpecError("system.template", "unknown system template", "TransmonSytem",
            ["TransmonSystem", "MultiTransmonSystem"])]
    json = Specs.emit_error_json(errs)
    @test occursin("\"ok\":false", replace(json, " " => ""))
    @test occursin("problem.Q", json)
    @test_throws Specs.SpecValidationError throw(Specs.SpecValidationError(errs))
end
