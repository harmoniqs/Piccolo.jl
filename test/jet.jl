# Performance analysis (type instabilities / runtime dispatch) — run manually:
#     julia --project=. -e 'using Piccolo, JET; JET.@report_opt rollout(...)'

@testitem "JET correctness analysis" tags=[:jet] begin
    using JET, Piccolo
    JET.test_package(Piccolo; target_modules = (Piccolo,), broken = true)  # TODO: 31 findings
end
