module PiccoloSpecsCorpus

#
# Frozen wire-format corpus for the Piccolo 2.0 program (Piccolissimo #427 / #423).
#
# Each fixture is one registered spec family rendered to the canonical wire form
# `canonical_json(full_dict(spec))`, with `structure_hash` / `problem_hash`
# recorded in `hashes.toml`. The corpus is the byte-stability baseline for the
# `spec_hash_stability` CI job: NOTHING in this directory may change bytes
# without a deliberate, reviewed re-freeze.
#
# Regeneration (must reproduce every byte):
#   julia --project=<env with this checkout dev'd> test/specs_corpus/generate.jl
#
# Disclosure note: fixtures reference public Piccolo names only (feature-level
# content) — safe for the open repo by construction.

using Piccolo: Specs
using Piccolo.Specs:
    ProblemSpec,
    SystemSpec,
    GoalSpec,
    PulseSpec,
    TemplateBlock,
    IntegratorSpec,
    WrapperSpec,
    ObjectiveTermSpec,
    Free,
    Fixed
using TOML

const TRANSMON_PARAMS =
    Dict{Symbol,Any}(:n_qubits => 1, :levels => 3, :ω => 4.0e0, :δ => -2.0e-1)
const RYDBERG_PARAMS = Dict{Symbol,Any}(:n_atoms => 2, :levels => 3)

"""
    corpus_specs() -> Vector{Pair{String,ProblemSpec}}

The frozen fixture set: slug => spec. Deterministic — no RNG, no timestamps.
Any edit here is a deliberate re-freeze and MUST be reviewed against the
wire-format freeze invariant before merging.
"""
function corpus_specs()
    specs = Vector{Pair{String,ProblemSpec}}()

    transmon(args...; kwargs...) = SystemSpec(;
        kind = :template,
        template = :TransmonSystem,
        params = TRANSMON_PARAMS,
        args...,
        kwargs...,
    )
    push!(
        specs,
        "smooth-unitary-zo" => ProblemSpec(;
            system = transmon(),
            goal = GoalSpec(kind = :unitary, gate = :X),
            pulse = PulseSpec(kind = :zero_order, T = 10.0, init = :zeros, seed = 42),
            problem = TemplateBlock(template = :SmoothPulseProblem, N = 50),
            integrator = IntegratorSpec(kind = :bilinear, alg = :tsit5),
        ),
    )
    push!(
        specs,
        "smooth-ket-zo" => ProblemSpec(;
            system = transmon(),
            goal = GoalSpec(kind = :ket, target = "e0"; initial = "g0"),
            pulse = PulseSpec(kind = :zero_order, T = 20.0, init = :random, seed = 7),
            problem = TemplateBlock(template = :SmoothPulseProblem, N = 100, Q = 100.0),
        ),
    )
    push!(
        specs,
        "spline-unitary-linear-freephase" => ProblemSpec(;
            system = transmon(),
            goal = GoalSpec(kind = :unitary, gate = :X),
            pulse = PulseSpec(kind = :linear_spline, T = 30.0, init = :default, seed = 0),
            problem = TemplateBlock(;
                template = :SplinePulseProblem,
                N = 200,
                free_phase = true,
                R_du = 1.0e-2,
                du_bound = 0.3,
            ),
            integrator = IntegratorSpec(kind = :spline, alg = :tsit5),
        ),
    )
    push!(
        specs,
        "spline-unitary-cubic-freephase" => ProblemSpec(;
            system = transmon(),
            goal = GoalSpec(kind = :unitary, gate = :CZ, subsystem_levels = [3, 3]),
            pulse = PulseSpec(kind = :cubic_spline, T = 40.0, init = :default, seed = 0),
            problem = TemplateBlock(;
                template = :SplinePulseProblem,
                N = 200,
                free_phase = true,
                objectives = [ObjectiveTermSpec(kind = :reg_ddu, weight = 1.0e-3)],
            ),
            integrator = IntegratorSpec(kind = :spline, alg = :magnus_gl4),
        ),
    )
    push!(
        specs,
        "spline-ket-cubic" => ProblemSpec(;
            system = transmon(),
            goal = GoalSpec(kind = :ket, target = "(g0+e0)/√2"),
            pulse = PulseSpec(kind = :cubic_spline, T = 15.0, init = :default, seed = 0),
            problem = TemplateBlock(
                template = :SplinePulseProblem,
                N = 75,
                free_phase = true,
            ),
            integrator = IntegratorSpec(kind = :spline, alg = :tsit5),
        ),
    )
    push!(
        specs,
        "bangbang-unitary-zo" => ProblemSpec(;
            system = transmon(),
            goal = GoalSpec(kind = :unitary, gate = :X),
            pulse = PulseSpec(kind = :zero_order, T = 8.0, init = :random, seed = 11),
            problem = TemplateBlock(template = :BangBangPulseProblem, N = 40, R = 1.0e-2),
        ),
    )
    push!(
        specs,
        "bangbang-ket-zo-freephase" => ProblemSpec(;
            system = transmon(),
            goal = GoalSpec(kind = :ket, target = "e0"),
            pulse = PulseSpec(kind = :zero_order, T = 12.0, init = :random, seed = 3),
            problem = TemplateBlock(;
                template = :BangBangPulseProblem,
                N = 60,
                free_phase = true,
                initial_phases = [0.0, 0.0],
            ),
        ),
    )
    push!(
        specs,
        "smooth-mintime-recipe" => ProblemSpec(;
            system = transmon(),
            goal = GoalSpec(kind = :unitary, gate = :X),
            pulse = PulseSpec(kind = :zero_order, T = 20.0, init = :zeros, seed = 0),
            problem = TemplateBlock(;
                template = :SmoothPulseProblem,
                N = 100,
                goal_treatment = :both,
                final_fidelity = 0.999,
                free_dt = Free(5.0, 25.0),
                objectives = [ObjectiveTermSpec(kind = :time, weight = 1.0)],
            ),
        ),
    )
    push!(
        specs,
        "sampling-wrapped" => ProblemSpec(;
            system = SystemSpec(;
                kind = :template,
                template = :TransmonSystem,
                params = TRANSMON_PARAMS,
            ),
            goal = GoalSpec(kind = :unitary, gate = :X),
            pulse = PulseSpec(kind = :zero_order, T = 10.0, init = :zeros, seed = 0),
            problem = TemplateBlock(template = :SmoothPulseProblem, N = 50),
            wrappers = [
                WrapperSpec(;
                    kind = :sampling,
                    variants = [
                        Dict{Symbol,Any}(:ω => 4.0e0),
                        Dict{Symbol,Any}(:ω => 4.0e-2 + 4.0e0),
                    ],
                    weights = [0.5, 0.5],
                ),
            ],
        ),
    )
    push!(
        specs,
        "transmon-globals-unpinned" => ProblemSpec(;
            system = SystemSpec(;
                kind = :template,
                template = :TransmonSystem,
                params = TRANSMON_PARAMS,
                global_params = Dict{Symbol,Any}(:ω => 4.0e0),
            ),
            goal = GoalSpec(kind = :unitary, gate = :X),
            pulse = PulseSpec(kind = :zero_order, T = 10.0, init = :zeros, seed = 0),
            problem = TemplateBlock(;
                template = :SmoothPulseProblem,
                N = 50,
                global_names = [:ω],
                global_bounds = Dict{Symbol,Any}(:ω => [3.9e0, 4.1e0]),
            ),
            integrator = IntegratorSpec(kind = :exponential, alg = :magnus_adapt4),
        ),
    )
    push!(
        specs,
        "rydberg-chain-cz-3level" => ProblemSpec(;
            system = SystemSpec(;
                kind = :template,
                template = :RydbergChainSystem,
                params = RYDBERG_PARAMS,
            ),
            goal = GoalSpec(kind = :unitary, gate = :CZ, subsystem_levels = [3, 3]),
            pulse = PulseSpec(kind = :linear_spline, T = 0.5, init = :default, seed = 0),
            problem = TemplateBlock(;
                template = :SplinePulseProblem,
                N = 100,
                free_phase = true,
                objectives = [ObjectiveTermSpec(kind = :leakage, weight = 1.0e-2)],
            ),
            integrator = IntegratorSpec(kind = :spline, alg = :magnus_adapt4),
        ),
    )
    # arbitrary inline unitary (nested [re,im] pairs) — the no-gate-name branch
    # 1-qubit Z rotation by pi/2 as nested [re,im]: diag(1, i)
    matrix = Any[Any[Any[1.0, 0.0], Any[0.0, 0.0]], Any[Any[0.0, 0.0], Any[0.0, 1.0]]]
    push!(
        specs,
        "smooth-unitary-inline-matrix" => ProblemSpec(;
            system = transmon(),
            goal = GoalSpec(kind = :unitary, matrix = matrix),
            pulse = PulseSpec(kind = :zero_order, T = 10.0, init = :zeros, seed = 0),
            problem = TemplateBlock(template = :SmoothPulseProblem, N = 50),
        ),
    )
    return specs
end

wire(spec::ProblemSpec) = Specs.canonical_json(Specs.full_dict(spec))

"""    main(dir = @__DIR__) — write every fixture + hashes.toml into `dir`."""
function main(dir = @__DIR__)
    hashes = Dict{String,Any}()
    for (slug, spec) in corpus_specs()
        json = wire(spec)
        open(joinpath(dir, "$slug.json"), "w") do io
            write(io, json)
        end
        hashes[slug] = Dict{String,Any}(
            "structure_hash" => Specs.structure_hash(spec),
            "problem_hash" => Specs.problem_hash(spec),
        )
    end
    open(joinpath(dir, "hashes.toml"), "w") do io
        TOML.print(io, hashes)
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end # module
