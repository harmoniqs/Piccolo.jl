# OSS schema drift gate (Plan 2 Task 1).
#
# Asserts the checked-in OSS variant `problemspec.oss.schema.json` is byte-for-byte
# what `emit_schema()` produces today (both normalized through a JSON3 round-trip so
# an editor's trailing newline is not spurious drift). If a registry registration or
# a parser allowed-key set changes, `emit_schema()` moves and this test fails —
# forcing `src/specs/schema/regenerate.jl` to be re-run and the file re-committed.
#
# The FULL variant's drift gate lives in the private Piccolissimo repo (it needs
# Piccolissimo loaded); this gate covers ONLY the public OSS variant so Piccolo's
# public CI never loads Piccolissimo. See plan review correction (Blocking #1).

@testitem "OSS schema drift: emit_schema() == checked-in problemspec.oss.schema.json" begin
    using Piccolo.Specs, JSON3
    Specs.register_all!()
    emitted = JSON3.write(JSON3.read(Specs.emit_schema()))
    file =
        joinpath(pkgdir(Piccolo), "src", "specs", "schema", "problemspec.oss.schema.json")
    @test isfile(file)
    checked = JSON3.write(JSON3.read(read(file, String)))
    @test emitted == checked
end

@testitem "OSS schema excludes private Piccolissimo enum values (scoped)" begin
    using Piccolo.Specs, JSON3
    Specs.register_all!()
    sch = JSON3.read(Specs.emit_schema())
    ctrl = sch.oneOf[1].properties.kind.const == "control" ? sch.oneOf[1] : sch.oneOf[2]
    p = ctrl.properties
    _enum(x) = String[string(v) for v in x]
    # Scope the exclusion to the specific registry-reflected enums. A full-text grep
    # would false-positive on PUBLIC names: `SplinePulseProblem` (problem.template),
    # `cubic_spline`/`linear_spline` (pulse.kind), and `magnus_*` (integrator.alg is a
    # hardcoded ODE-alg enum, present even in OSS) are all legitimately public.
    intk = _enum(p.integrator.properties.kind.enum)
    objk = _enum(p.problem.properties.objectives.items.properties.kind.enum)
    wrapk = _enum(p.wrappers.items.properties.kind.enum)
    solb = _enum(p.solver.properties.backend.enum)
    strat = _enum(p.solver.properties.strategy.enum)
    @test isempty(intersect(intk, ["exponential", "spline", "magnus_gl4", "magnus_adapt4"]))
    @test isempty(filter(x -> startswith(x, "hermite"), objk))
    @test !("robust" in wrapk)
    @test !("altissimo" in solb)
    @test isempty(intersect(strat, ["continuation", "staged"]))
end
