#!/usr/bin/env julia
# Emit the AUTHORITATIVE cross-language hash sidecars (Plan 2 Task 5).
#
# For every shared spec fixture in ./hash_fixtures/*.toml, parse it with
# `Piccolo.Specs.parse_spec` and write a `<name>.hashes.json` sidecar carrying the
# `structure_hash` and `problem_hash`. These sidecars are the SINGLE SOURCE OF
# TRUTH the TypeScript `hashing.ts` mirror is checked against (packages/schema/
# test/hashing.test.ts) — the TS side must reproduce them byte-for-byte from the
# same TOML. Copy the fixtures + sidecars into amicode's
# packages/schema/test/fixtures/hashes/ after running this.
#
# NOTE: this is distinct from julia/emit_fixture.jl in amicode, which is a solver
# STUB for the producer round-trip lane — NOT a hash emitter.
#
# Run:  cd packages/Piccolo.jl && julia --project=. src/specs/schema/emit_hash_fixtures.jl
using Piccolo
using Piccolo.Specs

const HERE = @__DIR__
const FIXDIR = joinpath(HERE, "..", "src", "specs", "schema", "hash_fixtures")

# Minimal JSON string escaper for the two hex digests (always ASCII hex, but keep
# it honest). We hand-write the object so no JSON dep is needed and the sidecar is
# a stable, human-diffable two-key object.
_jstr(s) = '"' * replace(String(s), '\\' => "\\\\", '"' => "\\\"") * '"'

function main()
    files = sort(filter(f -> endswith(f, ".toml"), readdir(FIXDIR)))
    isempty(files) && error("no *.toml fixtures found in $FIXDIR")
    for f in files
        name = f[1:(end-5)]  # strip ".toml"
        toml = read(joinpath(FIXDIR, f), String)
        spec = Specs.parse_spec(toml; format = :toml)
        sh = Specs.structure_hash(spec)
        ph = Specs.problem_hash(spec)
        sidecar = joinpath(FIXDIR, "$(name).hashes.json")
        open(sidecar, "w") do io
            print(
                io,
                "{\n  ",
                _jstr("structure_hash"),
                ": ",
                _jstr(sh),
                ",\n  ",
                _jstr("problem_hash"),
                ": ",
                _jstr(ph),
                "\n}\n",
            )
        end
        # Diagnostics: the canonical JSON strings behind each digest, so a
        # cross-language mismatch can be compared side-by-side.
        println("== $f ==")
        println("  structure_json = ", Specs.canonical_json(Specs.structure_fields(spec)))
        println("  problem_json   = ", Specs.canonical_json(Specs.full_dict(spec)))
        println("  structure_hash = ", sh)
        println("  problem_hash   = ", ph)
    end
    println("wrote $(length(files)) *.hashes.json sidecars to $FIXDIR")
end

main()
