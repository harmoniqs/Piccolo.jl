# Release Notes

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v2.0.0] - 2026-08-18

The multiple-shooting notation release: **N stays N** — the number of knot
points (nodes), aligned with the direct multiple-shooting literature — and
the API stops pretending silent defaults are safe.

### Removed (breaking)

- **`PiccoloOptions` inert fields removed** — `rollout_integrator`, `geodesic`,
  `zero_initial_and_final_derivative` were stored but never consulted
  (removal publicly promised since v1.x, #255/#283). Constructing with these
  kwargs now errors.
- **`SplinePulseProblem` no longer silently defaults spline pulses to
  piecewise-constant dynamics** (#275, #282): `BilinearIntegrator` never reads
  `:du`, so a `CubicSplinePulse` problem was optimizing a different waveform
  than its name promised (measured: optimizer-reported ~1e-8 infidelity for
  pulses that actually achieve ~1e-3). Cubic-spline problems must now either
  pass a spline-faithful integrator (`Piccolissimo.SplineIntegrator`) or
  explicitly request the PWC backend with `integrator_type = :pwc` (which
  warns). Linear-spline problems warn on the default. The deprecated
  `integrator_type = :spline` (which silently returned the PWC integrator)
  and `:ensemble` (never implemented) now raise informative errors.

### Added

- **Cubic-spline interior bounds** (#276): the between-knots spline envelope
  is now constrained (knot-only bounds miss interior extrema).
- **`solve!` returns `SolveStats`** — termination status, raw status string,
  NLP objective, IPM iterations, solve wall time, solver symbol (forwarded
  through `QuantumControlProblem.solve!`; live once DirectTrajOpt's next
  release lands).
- **`Piccolo.Specs` frozen wire-format corpus** (`test/specs_corpus/`) — the
  canonical ProblemSpec projection with recorded structure/problem hashes and
  a byte-stability test; the wire format is now frozen by CI.
- **Canonical Notation section** in the concepts guide: **N** = number of knot
  points (nodes); **"timestep"** exclusively the per-knot Δt_k; **T** =
  duration; PWC = piecewise-constant (zero-order hold).

### Fixed

- `RolloutStates`/sampling coverage: divergent pairings surface loudly
  (`rollout_divergence`); the constructor-side guards above are the fix.
- `workflow_dispatch` triggers on CI and Documentation workflows (retrigger
  without a push).

## [v1.0.0] - Upcoming

### Changed

- Complete documentation overhaul with concept-driven organization
- Migrated from PiccoloDocsTemplate to standard Documenter.jl

## [v0.3.1] - 2024-10-17

### Fixed

- Fixed and added tests to `RydbergChainSystem`

## [v0.3.0] - 2024-10-10

### Added

- `PiccoloOptions` to handle custom problem settings

### Changed

- Refactored trajectory initialization functions
- Improved documentation
- Typo fixes

## [v0.2.0] - 2024-02-22

### Added

- `EmbeddedOperator` to handle subspace gate optimization and leakage suppression
- Plotting methods for unitary populations

### Changed

- New quantum systems interface
  - Transmon system template
- Restructured the code base for easier quantum system and problem template development

### Removed

- Stale examples

### Fixed

- Robustness improvements objective test fixes
