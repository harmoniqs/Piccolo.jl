# ============================================================================ #
# #432 — the AdjointRobustness re-export seam must be complete
#
# The 3d extraction (Piccolissimo #432) moved the adjoint-robustness objective
# family (AdjointRobustnessObjective / KetAdjointRobustnessObjective and
# siblings), the propagation methods, and the RobustControlProblem wrapper
# into the `AdjointRobustness` submodule; `src/control/_control.jl` re-exports
# the family's names at top level. Same regression contract as the
# SplineIntegrators (#326) and SplineConstraints (#431) guards above: EVERY
# name the submodule exports resolves at (is re-exported through) `Piccolo`,
# so the seam cannot silently drift.
# ============================================================================ #

@testitem "AdjointRobustness re-export seam is complete (#432)" begin
    using Piccolo
    using Test

    sub = Piccolo.Control.AdjointRobustness
    exported_from_submodule = filter(n -> Base.isexported(sub, n), names(sub; all = false))
    @test !isempty(exported_from_submodule)

    # Interface functions the submodule extends on DirectTrajOpt generics
    # (objective_value, gradient!, hessian_structure, get_full_hessian) are
    # deliberately NOT re-exported here — same seam convention as the constraint
    # modules: interface functions come from DirectTrajOpt, which `using Piccolo`
    # already re-exports.
    interface_owned = filter(
        n ->
            !isdefined(sub, n) || n in (:objective_value, :gradient!, :get_full_hessian),
        exported_from_submodule,
    )
    family_owned = setdiff(exported_from_submodule, interface_owned)

    missing_at_top =
        filter(n -> !isdefined(Piccolo, n) || !Base.isexported(Piccolo, n), family_owned)
    @test isempty(missing_at_top)

    # The names that started this slice: Piccolissimo's re-export surface and
    # its Altissimo backend glue read exactly these.
    for n in (
        :AdjointRobustnessObjective,
        :KetAdjointRobustnessObjective,
        :RobustControlProblem,
        :ZKickPropagation,
        :test_robustness_objective,
        :objective_hvp!,
    )
        @test isdefined(Piccolo, n) && Base.isexported(Piccolo, n)
    end

    # The submodule name itself rides the seam (the SplineConstraints convention):
    # Piccolissimo's re-export seam reaches the module through this path.
    @test Piccolo.Control.AdjointRobustness isa Module
end
