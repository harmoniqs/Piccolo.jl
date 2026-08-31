# ============================================================================ #
# #323 — `duration` binding under a NamedTrajectories ≥ 0.9.3 co-resolve
#
# NT 0.9.3 added the TimeWarp module, which exports `duration`. Piccolo
# reexports both NamedTrajectories and Quantum (→ Quantum.Pulses, the defining
# module of `duration`), so with NT ≥ 0.9.3 in the environment the name used to
# resolve as ambiguous in Piccolo's own namespace — UndefVarError at top level,
# at `Piccolo.duration`, and at `Piccolo.Quantum.duration`, instead of binding
# to the Pulses function. Regression contract: `duration` binds to
# `Quantum.Pulses.duration` — the exact binding a non-colliding environment
# (NT < 0.9.3) has always provided — and is NOT the TimeWarp function.
#
# The test item reproduces the issue's namespace faithfully: `using Piccolo`
# only. `import NamedTrajectories as NT` deliberately binds just the module
# name — a plain `using NamedTrajectories` would re-create the two-package
# export collision in this module, where bare `duration` is inherently
# ambiguous (both bindings are exported; consumers bind from the defining
# module, as harmoniqs/Strumento.jl#19 does).
# ============================================================================ #

@testitem "duration binds to Quantum.Pulses with NamedTrajectories >= 0.9.3 co-resolved (#323)" begin
    using Piccolo
    import NamedTrajectories as NT

    # Collision precondition — the suite env must actually co-resolve NT ≥ 0.9.3
    # whose TimeWarp exports `duration`. If a resolve ever pins NT below 0.9.3
    # these fail loudly, so this test cannot pass vacuously.
    nt_version = pkgversion(NT)
    @test nt_version >= v"0.9.3"
    @test :duration in names(NT)

    # Top-level binding (the issue's repro surface: `using Piccolo; duration`).
    @test @isdefined(duration)
    @test duration === Piccolo.Quantum.Pulses.duration

    # The package's own namespaces must resolve the name too.
    @test isdefined(Piccolo, :duration)
    @test Piccolo.duration === Piccolo.Quantum.Pulses.duration
    @test isdefined(Piccolo.Quantum, :duration)
    @test Piccolo.Quantum.duration === Piccolo.Quantum.Pulses.duration

    # Explicit non-collision contract: the binding is the Pulses function
    # object, never the TimeWarp one.
    @test Piccolo.duration !== NT.TimeWarp.duration
end
