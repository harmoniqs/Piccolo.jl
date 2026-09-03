# JET integrates tightly with the Julia compiler and is gated to v1.12.x only.
# Older Julia ships older JET (0.9.x) which can't analyze the current Piccolo
# source. On 1.13+ (nightly / pre) the compiler shifts faster than JET can
# track, so findings are noisy and unactionable until JET catches up.
# Don't widen this gate.
#
# Performance analysis (type instabilities / runtime dispatch) — run manually:
#     julia --project=. -e 'using Piccolo, JET; JET.@report_opt rollout(...)'

@testitem "JET correctness analysis" tags=[:jet] begin
    if v"1.12" <= VERSION < v"1.13"
        using JET, Piccolo
        # STILL BROKEN, but no longer for the reason this comment used to give.
        #
        # The old note claimed the *single* "remaining finding" was a call to
        # `EnsembleSplineIntegrator`, defined nowhere in Piccolo or its deps. That call was
        # deleted on 2026-07-25 (`integrator_type` now offers only `:pwc`, with `:spline`
        # and `:ensemble` raising informative errors). Removing it did clear that finding —
        # and revealed that it was never the only one. JET reports **7** findings, none of
        # them `EnsembleSplineIntegrator`:
        #
        #   1-3. `ZeroOrderPulse` / `LinearSplinePulse` / `CubicSplinePulse` inner
        #        constructors: no matching method on the 1/2 union split where the
        #        interpolant is an `InterpolationsExt._ConstantExtrapolationFix{I}`
        #        (`quantum/primitives/pulses.jl:228, 338, 430`).
        #   4-5. `fidelity(::UnitaryTrajectory)`: `unembed(::AbstractMatrix{<:Number})` has
        #        no method — the goal is typed `AbstractPiccoloOperator`, so the non-
        #        `EmbeddedOperator` branch of the union is unhandled
        #        (`quantum/trajectories/rollouts_extensions.jl:900, 903`).
        #   6.   `_constraint_var_label`: `*(::String, ::Nothing)` on the 1/3 union split
        #        where `join(...)` returns `Nothing` (`control/display/inspect.jl:571`).
        #   7.   `_fidelity_at(::KetTrajectory, …)`: no matching
        #        `_phased_ket_goal(::Vector{ComplexF64}, ::Any, ::Int64)`
        #        (`control/display/inspect.jl:691`) — the free-phase display path.
        #
        # So "flip `broken = true`" is NOT a valid acceptance test for the integrator-type
        # cleanup; it needs those 7 fixed. Items 4-5 and 7 touch fidelity reporting and are
        # the ones worth doing first. Update this list whenever the count changes.
        JET.test_package(Piccolo; target_modules = (Piccolo,), broken = true)
    else
        @info "Skipping JET correctness analysis on Julia $VERSION (gated to 1.12.x)"
        @test true
    end
end
