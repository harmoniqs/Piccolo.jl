# Performance analysis (type instabilities / runtime dispatch) — run manually:
#     julia --project=. -e 'using Piccolo, JET; JET.@report_opt rollout(...)'

@testitem "JET correctness analysis" tags=[:jet] begin
    using JET, Piccolo
    # Remaining finding: `EnsembleSplineIntegrator` is referenced inside
    # `SplinePulseProblem(...)` when `integrator_type = :ensemble`, but no method
    # is defined in Piccolo or any of its declared deps. The constructor errors
    # the user toward Piccolissimo's `SplineIntegrator`. Left as broken until the
    # ensemble-integrator path is either deleted or reified as an extension.
    JET.test_package(Piccolo; target_modules = (Piccolo,), broken = true)
end
