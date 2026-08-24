module Piccolo

using Reexport

# Foundation packages (stay separate - they're generic)
@reexport using TrajectoryIndexingUtils
@reexport using NamedTrajectories
@reexport using DirectTrajOpt

# Quantum objects: systems, gates, pulses, trajectories
include("quantum/_quantum.jl")
@reexport using .Quantum

# Optimal control: objectives, constraints, problem templates
include("control/_control.jl")
@reexport using .Control

# Declarative problem specs: wire-format ProblemSpec, parser, registries, hashes
include("specs/_specs.jl")
@reexport using .Specs

# Visualizations
include("visualizations/_visualizations.jl")
@reexport using .Visualizations

# Documentation caching utilities
include("docs_cache.jl")
@reexport using .DocsCache

# ============================================================================= #
# Precompile workload
# ============================================================================= #
#
# `Specs.precompile_workload()` sweeps the **registry-enumerated type universe**:
# one tiny problem per `(template, pulse_kind, trajectory_kind)` triple the
# registered templates admit, plus the `sampling` wrapper, plus the declarative
# parse → materialize → one-Ipopt-step path.
#
# This is the point of the parametric-type rewrite. A template is a constrained
# alias of `QuantumControlProblem{Tag, QT}`, so the set of concrete problem types a
# spec can name is computable from the registry — and `structure_hash` covers
# exactly the type-determining spec fields. Therefore a cloud worker routed by
# `structure_hash` is guaranteed to have the concrete types already compiled: same
# `structure_hash` ⇒ same concrete types ⇒ no JIT.
#
# The workload body runs during precompilation, i.e. *before* any `__init__`, so it
# calls `Specs.register_all!()` itself (it is idempotent). Every combination is
# individually guarded, so a failure degrades the workload rather than breaking
# precompilation.
using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    @compile_workload begin
        Specs.precompile_workload(; tier1_only = true)
    end
end

end
