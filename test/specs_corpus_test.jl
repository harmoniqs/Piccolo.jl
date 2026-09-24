@testitem "specs corpus: canonical wire emission is byte-stable" tags=[:specs_corpus] begin
    using Piccolo: Specs
    using TOML

    corpus_dir = joinpath(@__DIR__, "specs_corpus")
    include(joinpath(corpus_dir, "generate.jl"))
    using .PiccoloSpecsCorpus: corpus_specs, wire

    recorded = TOML.parsefile(joinpath(corpus_dir, "hashes.toml"))
    @test !isempty(recorded)

    n_fixtures = 0
    for (slug, spec) in corpus_specs()
        path = joinpath(corpus_dir, "$slug.json")
        @test isfile(path)
        expected = read(path, String)
        emitted = wire(spec)
        @test emitted == expected
        # the recorded hashes must match what the live code computes
        @test recorded[slug]["structure_hash"] == Specs.structure_hash(spec)
        @test recorded[slug]["problem_hash"] == Specs.problem_hash(spec)
        # the round-trip identity the wire format promises
        @test Specs.parse_spec(emitted; format = :json) isa Specs.ProblemSpec
        global n_fixtures += 1
    end
    @test n_fixtures == length(recorded)  # no orphan files, no missing fixtures
end
