# Regenerate the checked-in OSS ProblemSpec JSON Schema variant from the registries.
#
# OSS variant only: Piccolo is loaded but Piccolissimo is NOT, so `emit_schema()`
# reflects only the public Piccolo registrations (the `bilinear` integrator
# sentinel, the `sampling` wrapper, the public objective terms, the seven systems,
# and the three problem templates). The FULL variant (Piccolo + private
# Piccolissimo) is emitted by a sibling script that lives in the Piccolissimo repo
# — never here — so that private capability names (exponential/spline integrators,
# hermite_* objective terms, the robust wrapper, the altissimo solver, the
# continuation/staged strategies) are never written into the public Piccolo repo
# or exercised in public Piccolo CI. See plan review correction (Blocking #1),
# 2026-07-22.
#
# Run:   cd packages/Piccolo.jl && julia --project=. src/specs/schema/regenerate.jl
# Drift: `git diff --exit-code src/specs/schema/problemspec.oss.schema.json`
#        (the schema-drift @testitem in drift.jl asserts the same invariant).

using Piccolo

Piccolo.Specs.register_all!()   # idempotent; Piccolo already self-registers at load

here = @__DIR__
out = joinpath(here, "problemspec.oss.schema.json")
open(out, "w") do io
    write(io, Piccolo.Specs.emit_schema())
end
println("wrote ", out)
