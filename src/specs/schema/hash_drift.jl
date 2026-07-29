# Cross-language hash-sidecar drift gate.
#
# `src/specs/schema/hash_fixtures/*.hashes.json` are the AUTHORITATIVE structure/
# problem hashes for the shared spec fixtures — the values amicode's TypeScript
# `hashing.ts` mirror and Prova's vendored fixtures are checked against. They are
# produced by `emit_hash_fixtures.jl`, which a human has to remember to re-run.
#
# Without this gate, a change to `hashes.jl` (the canonical-JSON rule, a
# structure-field carve-out, a numeric formatting detail) leaves the checked-in
# sidecars stale, and every downstream repo keeps validating against the OLD
# digests — so both sides agree with each other while both have drifted from the
# Julia implementation they are supposed to mirror. That is the exact failure the
# schema's `drift.jl` gate exists to prevent; this is the same gate for hashes.
#
# The test recomputes each digest through the SAME path as the emitter
# (`parse_spec(toml; format = :toml)` → `structure_hash`/`problem_hash`) and
# compares against the sidecar on disk. If it fails, re-run:
#
#   julia --project=. src/specs/schema/emit_hash_fixtures.jl
#
# and re-vendor the sidecars into amicode (packages/schema/test/fixtures/hashes/)
# and Prova (test/fixtures/hashes/) — a failure here is a CROSS-REPO change.

@testitem "hash sidecar drift: checked-in *.hashes.json == recomputed digests" begin
    using Piccolo, Piccolo.Specs, JSON3
    Specs.register_all!()

    fixdir = joinpath(pkgdir(Piccolo), "src", "specs", "schema", "hash_fixtures")
    @test isdir(fixdir)

    tomls = sort(filter(f -> endswith(f, ".toml"), readdir(fixdir)))
    @test !isempty(tomls)

    for f in tomls
        name = f[1:(end-5)]  # strip ".toml"
        sidecar = joinpath(fixdir, "$(name).hashes.json")
        # Every fixture must ship a sidecar: adding a .toml without re-running the
        # emitter is itself drift (amicode/Prova would have nothing to pin to).
        @test isfile(sidecar)
        isfile(sidecar) || continue

        spec = Specs.parse_spec(read(joinpath(fixdir, f), String); format = :toml)
        checked = JSON3.read(read(sidecar, String))
        @test String(checked.structure_hash) == Specs.structure_hash(spec)
        @test String(checked.problem_hash) == Specs.problem_hash(spec)
    end

    # ...and no orphan sidecars whose fixture was deleted or renamed.
    sidecars = filter(f -> endswith(f, ".hashes.json"), readdir(fixdir))
    for s in sidecars
        @test isfile(joinpath(fixdir, s[1:(end-length(".hashes.json"))] * ".toml"))
    end
end

# The digests are only meaningful if the canonical JSON behind them is stable, so
# pin that too: `structure_hash` must be computed from `structure_fields` and
# `problem_hash` from `full_dict`, each through `canonical_json`. This is what the
# TypeScript mirror reimplements, so a silent change of *which dict* is hashed
# would break amicode without moving a single digest in an obvious way.
@testitem "hash sidecar drift: digests are sha256 over the canonical JSON of the documented dicts" begin
    using Piccolo, Piccolo.Specs, SHA

    fixdir = joinpath(pkgdir(Piccolo), "src", "specs", "schema", "hash_fixtures")
    tomls = sort(filter(f -> endswith(f, ".toml"), readdir(fixdir)))

    for f in tomls
        spec = Specs.parse_spec(read(joinpath(fixdir, f), String); format = :toml)
        sh = bytes2hex(sha256(Specs.canonical_json(Specs.structure_fields(spec))))
        ph = bytes2hex(sha256(Specs.canonical_json(Specs.full_dict(spec))))
        @test sh == Specs.structure_hash(spec)
        @test ph == Specs.problem_hash(spec)
    end
end
