# AGENTS.md — Piccolo.jl

Piccolo.jl is the quantum optimal control package at the center of the
Harmoniqs stack: trajectory types (`UnitaryTrajectory`, `KetTrajectory`),
quantum systems, the `EmbeddedOperator` gate-embedding layer, and the
`SmoothPulseProblem` / `SplinePulseProblem` / `BangBangPulseProblem`
family, solved through DirectTrajOpt backends.

## Conventions

- Standard Julia package: `Pkg.test()` (the test suite includes Aqua and
  JET runs — do not skip them locally and expect CI to agree).
- Docs: docs.harmoniqs.co/Piccolo/dev.
- Siblings: NamedTrajectories.jl (trajectory datastructures),
  DirectTrajOpt.jl (solver layer), Piccolissimo.jl (extensions — has its
  own ADR corpus, read it before changing shared interfaces).

## Toolchain lessons (load-bearing, learned the hard way)

1. **`EmbeddedOperator` takes GATES-dict symbols, never bare matrices.**
   Passing a raw unitary matrix (`EmbeddedOperator(CNOT_matrix, sys)`)
   fails at *runtime* with `ArgumentError: Operator must be a valid
   gate`. Use a recognized Symbol or constructor path, e.g.
   `EmbeddedOperator(:CNOT, sys)`. Evidence: failed CNOT solve
   r20260713-212934Z-0d81.

2. **`BilinearIntegrator` is the default-reliable integrator.** For
   standard transmon gate optimization, the Piccolissimo `MagnusGL4`
   path matched baseline fidelity (F = 0.999981) only with substantially
   more complexity (variable time, D tuning); a minimal variant reached
   F = 0.9775. Reserve Magnus-series paths for strong-drive /
   variable-time regimes that specifically need them.

3. **Smooth parameterization reliably reaches F > 0.9999** for
   transmon X-class gates. When a smooth solve will not converge,
   suspect the problem setup (levels, drive bounds, timestep
   resolution) before the parameterization.

4. **Evaluating a `CubicSplinePulse` AT its last knot is a
   known-broken path** in pulses.jl — clamp evaluation strictly inside
   the knot range. Lab-frame truth evaluations hit this silently.

5. **Warm-start idiom:** `traj = load_traj("path/to/pulse.jld2")` passed
   as the initial guess to the problem constructor. `load_traj` is the
   correct loader for saved pulses in this codebase.

## Provenance

Seeded 2026-08-22 from amicode memory cards (piccolo-gate-operator-api,
piccolissimo-underperforms, smooth-pulses, pulse-shape-rwa) under the
knowledge-routing law: lessons live in the repo where agents use them.
Decisions with history go to ADRs in `docs/adr/`; this file carries the
standing lessons an agent needs before touching code.
