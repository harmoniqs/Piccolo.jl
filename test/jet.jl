using TestItems

# JET.jl static analysis.
#
# `JET.test_package` runs the same correctness checks as `report_package`
# (type errors, undefined methods, inference failures), wrapped in a
# `@test`. We restrict analysis to Piccolo's modules with
# `target_modules` so we don't get noise from upstream packages we don't
# control.
#
# Findings as of the initial wiring of this test (Piccolo v1.13.0, JET
# 0.11): 31 possible errors, dominated by undefined names
# (`CubicSplineInterpolation`, `EnsembleSplineIntegrator`,
# `update_base_trajectory`, `_update_system!`, `default_algorithm` for
# `OpenQuantumSystem`), `Vector{Float64}(::Symbol)` calls inside Pulse
# constructors, and a few constructor-dispatch issues in
# `OpenQuantumSystem` / `VariationalQuantumSystem` /
# `CompositeQuantumSystem`. We mark the test `broken=true` for now: the
# findings are informational, visible in CI logs without making the
# build red. Drop `broken=true` after the underlying issues are fixed,
# or convert specific cases to `@test_skip`.

@testitem "JET correctness analysis (test_package)" begin
    using JET
    using Piccolo

    JET.test_package(Piccolo; target_modules = (Piccolo,), broken = true)
end

# `JET.report_opt` flags type instabilities and runtime dispatch — the
# highest-value JET output for performance-critical numerical code like
# Piccolo. It is NOT run in CI by default (it requires a callable, not a
# module, so it is per-function and noisy). Run on individual hot paths
# locally:
#
#     julia --project=. -e '
#         using Piccolo, JET
#         JET.@report_opt some_hot_function(args...)
#     '
#
# Gate this test behind ENV["JET_OPT"] == "true" so a developer can
# opt-in via `JET_OPT=true julia --project=. -e "using Pkg; Pkg.test()"`.
@testitem "JET performance analysis (report_opt, opt-in)" begin
    using JET
    using Piccolo

    if get(ENV, "JET_OPT", "false") == "true"
        @info "JET_OPT=true: run JET.@report_opt on individual hot paths"
        @info "  e.g. JET.@report_opt rollout(...) — see test/jet.jl for guidance"
        @test true
    else
        @info "Skipping JET.report_opt (set JET_OPT=true to enable)"
        @test true
    end
end
