# Tasks

## Todo

### Add plot recipes for pulse types
- **Labels:** feature, visualization

Add plotting support for pulse objects (`SplinePulse`, `ZeroOrderHoldPulse`, `LinearPulse`),
akin to the existing unitary population plots. Should integrate with the Makie extension.

### Improved testing for re-solving problems with global variables
- **Labels:** testing

Expand and improve test coverage for re-solving optimal control problems that use
global variables, including edge cases and regression tests.

### Saving and loading pulse objects
- **Labels:** feature

Now that typed pulse objects (`SplinePulse`, `ZeroOrderHoldPulse`, `LinearPulse`) exist
and are an important work product, add serialization support (save/load) so pulses can
be persisted and resampled across sessions (likely via JLD2).

### Sort out integrator usage with SplinePulseProblem and pulse type constraints
- **Labels:** refactor

Clarify and enforce which integrators are valid for each problem type:
- Private Piccolissimo integrators should work with `SplinePulseProblem`
- `SmoothPulseProblem` and `BangBangProblem` should only accept zero-order hold pulses

## In Progress

## Done

### Implement tasksmd-sync for Piccolo.jl
- **Labels:** meta, tooling

Add `TASKS.md` and `.github/workflows/tasksmd-sync.yml` to sync tasks to the
harmoniqs GitHub Project board. This PR.

### Typed drive system (LinearDrive, NonlinearDrive) with auto-Jacobian

Merged in PR #79. Added `LinearDrive` and `NonlinearDrive` types with automatic
Jacobian generation via ForwardDiff.
