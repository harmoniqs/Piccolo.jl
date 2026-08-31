# ============================================================================ #
# #326 — the SplineIntegrators re-export seam must be complete
#
# The 3b extraction (Piccolissimo #430) moved the spline family into the
# `SplineIntegrators` submodule; `src/control/integrators.jl` re-exports the
# headline names at top level. The second export block (rollout helpers +
# sensitivity-kick family, incl. `unitary_rollout_trajectory`) was missed, so
# downstream `using Piccolo` callers — Piccolissimo's benchmark/correctness
# suite — hit UndefVarError on names the package ships.
#
# Regression contract: EVERY name exported by the submodule resolves at (is
# re-exported through) `Piccolo`, so the seam cannot silently drift again.
# ============================================================================ #

@testitem "SplineIntegrators re-export seam is complete (#326)" begin
    using Piccolo
    using Test

    sub = Piccolo.Control.QuantumIntegrators.SplineIntegrators
    exported_from_submodule = filter(n -> Base.isexported(sub, n), names(sub; all = false))
    @test !isempty(exported_from_submodule)

    # The CommonInterface-owned methods the submodule re-exports (eval_constraint_*)
    # are deliberately NOT re-exported at Piccolo top level — same convention as
    # the constraint types: interface functions come from CommonInterface.
    interface_owned = filter(n -> !isdefined(sub, n), exported_from_submodule)
    spline_owned = setdiff(exported_from_submodule, interface_owned)

    missing_at_top =
        filter(n -> !isdefined(Piccolo, n) || !Base.isexported(Piccolo, n), spline_owned)
    @test isempty(missing_at_top)

    # The name that started this: Piccolissimo's correctness benchmark calls it
    # bare after `using Piccolo`.
    @test isdefined(Piccolo, :unitary_rollout_trajectory)
end
