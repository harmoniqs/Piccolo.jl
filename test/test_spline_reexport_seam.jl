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

# ============================================================================ #
# #431 — the SplineConstraints re-export seam must be complete
#
# The 3c extraction (Piccolissimo #431) moved the spline-shape constraint
# family into the `SplineConstraints` submodule; `src/control/constraints.jl`
# re-exports the family's names at top level. Same regression contract as the
# SplineIntegrators guard above: EVERY name the submodule exports resolves at
# (is re-exported through) `Piccolo`, so the seam cannot silently drift — the
# exact failure mode #326 recorded for the 3b module's second export block.
# ============================================================================ #

@testitem "SplineConstraints re-export seam is complete (#431)" begin
    using Piccolo
    using Test

    sub = Piccolo.Control.QuantumConstraints.SplineConstraints
    exported_from_submodule = filter(n -> Base.isexported(sub, n), names(sub; all = false))
    @test !isempty(exported_from_submodule)

    # The CommonInterface-owned methods the submodule imports but does not
    # export (evaluate!, jacobian!, ...) never enter this list — the submodule
    # imports them without export, per the seam convention: interface functions
    # come from CommonInterface, which `using Piccolo` already reaches through
    # the DirectTrajOpt re-export.
    interface_owned = filter(n -> !isdefined(sub, n), exported_from_submodule)
    spline_owned = setdiff(exported_from_submodule, interface_owned)

    missing_at_top =
        filter(n -> !isdefined(Piccolo, n) || !Base.isexported(Piccolo, n), spline_owned)
    @test isempty(missing_at_top)

    # The names that started this slice: the template binds interior spline
    # bounds through the top-level seam (AC 2), and Piccolissimo's re-export
    # surface reads the same names.
    for n in (
        :CubicSplineBoundConstraint,
        :HermiteSmoothAccelerationConstraint,
        :NonlinearSegmentConstraint,
        :OptimizedNonlinearKnotPointConstraint,
        :ConstraintStencilTable,
    )
        @test isdefined(Piccolo, n) && Base.isexported(Piccolo, n)
    end
end
