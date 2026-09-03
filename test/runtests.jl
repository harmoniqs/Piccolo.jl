using TestItemRunner

# @testsnippet: inlined into every @testitem that lists it in setup=[...].
# Mirrors DirectTrajOpt's own DTOTestHelpers — the template testitems
# (src/control/templates/*.jl) use the same bare DTO exports.
@testsnippet PiccoloTemplateHelpers begin
    using Piccolo
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra
    using Test
end

@run_package_tests
