@testitem "Aqua quality assurance" tags=[:aqua] begin
    using Aqua, Piccolo

    Aqua.test_all(
        Piccolo;
        # `hessian_structure` is exported by three different DirectTrajOpt submodules
        # (CommonInterface, Constraints, Integrators); the conflict makes it appear
        # undefined at Piccolo's surface even though all three sub-definitions exist.
        # This is a DirectTrajOpt-side issue, not Piccolo's to fix.
        undefined_exports = (broken = true,),
    )
end
