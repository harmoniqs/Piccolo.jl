# JET integrates tightly with the Julia compiler and is gated to v1.12 only.
# Older Julia ships older JET (0.9.x) which can't analyze the current Piccolo
# source. Don't relax this gate.
#
# Performance analysis (type instabilities / runtime dispatch) — run manually:
#     julia --project=. -e 'using Piccolo, JET; JET.@report_opt rollout(...)'

@testitem "JET correctness analysis" tags=[:jet] begin
    if VERSION < v"1.12"
        @info "Skipping JET correctness analysis on Julia $VERSION (requires >= 1.12)"
        return
    end
    using JET, Piccolo
    # Remaining finding: `EnsembleSplineIntegrator` is referenced inside
    # `SplinePulseProblem(...)` when `integrator_type = :ensemble`, but no method
    # is defined in Piccolo or any of its declared deps. The constructor errors
    # the user toward Piccolissimo's `SplineIntegrator`. Left as broken until the
    # ensemble-integrator path is either deleted or reified as an extension.
    JET.test_package(Piccolo; target_modules = (Piccolo,), broken = true)
end
