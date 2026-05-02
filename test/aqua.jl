@testitem "Aqua quality assurance" tags=[:aqua] begin
    using Aqua, Piccolo

    Aqua.test_all(
        Piccolo;
        ambiguities = false,                                   # TODO: triage 4 ambiguities
        deps_compat = (check_extras = false, ignore = [:Libdl]),
        undefined_exports = (broken = true,),                  # TODO: 8 stale exports
    )
end
