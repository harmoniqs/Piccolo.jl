@testitem "play reproduces strumento's compiled int16 envelopes and gain codes exactly" tags =
    [:playback] begin
    using JSON3
    using SHA

    path = joinpath(@__DIR__, "fixtures", "playback_parity", "loopback_demo.json")
    # pinned in fixtures/playback_parity/README.md: a refreshed corpus must update the README
    @test bytes2hex(open(sha256, path)) ==
          "f3fbafc9061aca41f333658ab72700f42748a1e02dfdd17b5b7f40e89512d241"
    corpus = JSON3.read(read(path, String), Dict{String,Any})
    @test corpus["schema"] == "strumento.parity-corpus/1"
    playback = corpus["playback"]
    generator = playback["generator"]
    @test generator["gain_rounding"] == "trunc"

    channels(d) = [Float64.(d["i"]); isnothing(d["q"]) ? Float64[] : Float64.(d["q"])]
    rows(d) = isnothing(d["q"]) ? 1 : 2
    function corpus_pulse(p)
        kind = p["kind"]
        if kind in ("gaussian", "erf")
            P = p["params"]
            a, σ, μ = Float64.(P["amplitude"]), Float64.(P["sigma"]), Float64.(P["center"])
            T = p["duration_us"]
            return kind == "gaussian" ? GaussianPulse(a, σ, μ, T) : ErfPulse(a, σ, μ, T)
        end
        knots = Float64.(p["times_us"]) .- p["times_us"][1]
        u = permutedims(reshape(channels(p["values"]), length(knots), rows(p["values"])))
        kind == "hold" && return ZeroOrderPulse(u, knots)
        kind == "linear" && return LinearSplinePulse(u, knots)
        du = permutedims(reshape(channels(p["tangents"]), length(knots), rows(p["values"])))
        kind == "cubic" && return CubicSplinePulse(u, du, knots)
        error("unknown corpus pulse kind $kind")
    end

    @test length(corpus["pulses"]) == 13
    matched = String[]
    for p in corpus["pulses"]
        # the recipe: play the RAW pulse with full_scale = peak, the DDS gain post-quantization
        d = merge(
            playback["spec"],
            Dict{String,Any}(
                "full_scale" => p["peak"],
                "post_gain" => Dict{String,Any}(
                    "fraction" => p["gain_frac"],
                    "max_code" => generator["gain_code_max"],
                    "rounding" => generator["gain_rounding"],
                ),
            ),
        )
        spec = PlaybackSpec(d)
        _, info = play(spec, corpus_pulse(p))
        env = p["envelope"]
        want = permutedims(reshape(Int.(channels(env)), p["n_padded"], rows(env)))
        ok = info.codes == want && info.post_gain_code == p["gain_code"]
        @test info.codes == want
        @test size(info.codes, 2) == p["n_padded"]
        @test info.post_gain_code == p["gain_code"]
        @test info.post_gain_code / generator["gain_code_max"] ≈ p["gain_frac_compiled"]
        # the played amplitude is the code times the as-compiled gain, not the requested one
        @test info.samples ≈
              info.codes .* (p["peak"] / spec.max_code) .* p["gain_frac_compiled"]
        ok || @info "playback parity mismatch" p["name"]
        ok && push!(matched, p["name"])
    end
    @test length(matched) == 13
end
