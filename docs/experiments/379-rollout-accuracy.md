# Rollout accuracy — exact ZOH, spline default (#379)

Found by Piccolo #375 / Legatissimo #24 (stack-harmony W1): the closed-system default
rollout (`MagnusAdapt4`) did not meet its requested tolerance. Measured on
`379-rollout` (branched from `fcbe0f15`), Julia 1.12, `OPENBLAS_NUM_THREADS=1`.

## What changed

- **ZOH pulses roll out exactly.** With `algorithm = nothing` (the default), a
  `ZeroOrderPulse` on a time-independent system is propagated by
  `PiecewiseExponential`: one exponential per segment between consecutive knot/save
  times (`eigen(Hermitian(H))` for Hermitian closed systems, `exp(-iHΔt)` for
  non-Hermitian ones, `exp(𝓛Δt)` of the Lindbladian superoperator on `vec(ρ)` for open
  systems). Save times inside a knot interval are reached by the exact
  partial-interval exponential; `sol.t == save_times` and `n_save` are unchanged. This
  covers `UnitaryTrajectory`, `KetTrajectory`, `MultiKetTrajectory`,
  `DensityTrajectory`, `MultiDensityTrajectory` (and so `SamplingTrajectory`), their
  `rollout` / `rollout!`, `fidelity(qtraj)`, `sync_trajectory!`, and the
  `NamedTrajectory` helpers (`rollout_fidelity`, `ket_rollout`, `unitary_rollout`, …)
  when the interpolation resolves to `:constant`. "Time-independent" means
  `time_dependent == false` on a system assembled from drive operators; a system built
  from a bare `H(u, t)` function keeps the ODE path, because its `time_dependent` flag
  is an unchecked declaration (the time-dependent `MinimumTimeProblem` template test
  builds `H = Z + u₁cos(ωt)X + …` with the flag left `false`).
- **Linear / cubic spline pulses on a time-independent system default to `Tsit5()`**
  (was `MagnusAdapt4()` for Hermitian closed systems). Everything else — analytic
  pulses, time-dependent systems — keeps `default_algorithm(system)` (`MagnusAdapt4()`
  for Hermitian closed systems, `Tsit5()` otherwise), which is unchanged.
- In-place `rollout!` with `algorithm = nothing` reuses the algorithm the stored
  solution was produced with (the field is typed by the first solve; this also fixes
  the old `MethodError` when a trajectory built with e.g. `Tsit5()` went through
  `solve!` → `sync_trajectory!`).
- An explicit `algorithm` is always honoured; `MagnusAdapt4()` stays available.

## Why `MagnusAdapt4`'s tolerance is inert

OrdinaryDiffEqLinear's `MagnusAdapt4` embedded error estimate is
`u − exp(y₆)·uprev`. When the generator does not depend on the state, its stage
evaluations at `t + dt/2` coincide (`Q₃ = Q₅ = 0`) and those at `t + dt` coincide
(`Q₄ = Q₆`), so `v₄ − y₆ = −[Q₁, Q₄]/12` with `Q₄ = dt·(A(t+dt) − 2A(t+dt/2) + A(t))`
— a second difference, **identically zero whenever `H(t)` is affine in `t`** (every
linear-spline interval, every ZOH interval). The step size is then set by the `tstops`
alone (101 accepted steps, 0 rejected, at every tolerance from 1e-6 to 1e-12 on the
9-level fixture), so the error has a tolerance-independent floor. On ZOH pulses the
`t + dt` stage additionally samples the *next* knot's value at each interval end.
(The `MagnusAdapt4Alg` "tol is inert on KetTrajectory" note in
`integrators_spline/algorithms.jl` is a different mechanism — the Ket spline cell is a
fixed-step Chebyshev cell — but the `UnitaryTrajectory` spline cell's adaptive
`MagnusAdapt4` solve has this same degenerate estimator on linear splines.)

## Tolerance semantics

`abstol` / `reltol` are SciML's **local, per-step** tolerances. The global terminal
error tracks them up to a factor that grows with `‖H‖·T`: ≤ tol on the qubit/qutrit
fixtures, ~150× tol on the 9-level 200 ns block. Request roughly
`tol ≈ ε / (‖H‖·T)` for a global error ε. ZOH rollouts ignore the tolerances (exact).

## Measurements

Fixtures: `N = 11` knots, controls `0.4·(sin, cos)` over a half period, cubic tangents
`0.3·(cos, −sin)`. `qubit_2`: `H = Z/2 + u₁X + u₂Y`, T = 2. `qutrit_3`: anharmonic
transmon, `α = −0.2`, T = 2. `block_9`: two coupled qutrits (`2π·0.01` exchange,
`2π·0.02` drives), T = 200. Error = ‖U(T) − U_ref(T)‖_F against
`Tsit5(abstol = reltol = 1e-14)` (whose own ZOH error is the ~1e-12 floor in the first
column).

| system | basis | tol | default (new) | MagnusAdapt4 (old default) |
|---|---|---|---|---|
| qubit_2 | zoh | 1e-08 | 9.1e-13 | 1.8e-05 |
| qubit_2 | zoh | 1e-10 | 9.1e-13 | 1.4e-06 |
| qubit_2 | zoh | 1e-12 | 9.1e-13 | 3.0e-07 |
| qubit_2 | linear | 1e-08 | 2.8e-14 | 1.5e-09 |
| qubit_2 | linear | 1e-10 | 4.7e-15 | 1.5e-09 |
| qubit_2 | linear | 1e-12 | 1.1e-15 | 1.5e-09 |
| qubit_2 | cubic | 1e-08 | 9.3e-13 | 6.3e-10 |
| qubit_2 | cubic | 1e-10 | 2.7e-13 | 3.4e-11 |
| qubit_2 | cubic | 1e-12 | 5.0e-15 | 8.7e-13 |
| qutrit_3 | zoh | 1e-08 | 5.4e-12 | 3.1e-04 |
| qutrit_3 | zoh | 1e-10 | 5.4e-12 | 2.9e-05 |
| qutrit_3 | zoh | 1e-12 | 5.4e-12 | 1.8e-06 |
| qutrit_3 | linear | 1e-08 | 3.0e-13 | 1.6e-09 |
| qutrit_3 | linear | 1e-10 | 4.5e-14 | 1.6e-09 |
| qutrit_3 | linear | 1e-12 | 2.0e-15 | 1.6e-09 |
| qutrit_3 | cubic | 1e-08 | 3.0e-12 | 5.2e-09 |
| qutrit_3 | cubic | 1e-10 | 8.6e-13 | 1.9e-10 |
| qutrit_3 | cubic | 1e-12 | 1.7e-14 | 4.2e-12 |
| block_9 | zoh | 1e-08 | 3.5e-12 | 7.8e-06 |
| block_9 | zoh | 1e-10 | 3.5e-12 | 1.7e-06 |
| block_9 | zoh | 1e-12 | 3.5e-12 | 2.9e-07 |
| block_9 | linear | 1e-08 | 1.6e-06 | 1.8e-04 |
| block_9 | linear | 1e-10 | 1.5e-08 | 2.0e-04 |
| block_9 | linear | 1e-12 | 1.5e-10 | 1.6e-04 |
| block_9 | cubic | 1e-08 | 1.6e-06 | 2.3e-07 |
| block_9 | cubic | 1e-10 | 1.5e-08 | 2.4e-09 |
| block_9 | cubic | 1e-12 | 1.5e-10 | 2.6e-11 |

Tsit5 tracks the tolerance everywhere; MagnusAdapt4 does not on ZOH or linear splines
(the `block_9` linear floor is 1.6e-4 at every tolerance). MagnusAdapt4 is better on
`block_9` cubic (its estimator is non-degenerate there): pass it explicitly if that is
your regime.

ZOH rollout vs the exponential integrators' collocation propagation — the largest
per-interval defect `x_{k+1} − exp(G(u_k)Δt_k)x_k` of `HermitianExponentialIntegrator` /
`NonHermitianExponentialIntegrator` evaluated on `NamedTrajectory(qtraj, N)` (default
`n_save = 101`):

| ZOH fixture | new default | old default |
|---|---|---|
| unitary 2 | 6.2e-15 | 1.0e-05 (MagnusAdapt4) |
| unitary 9 (200 ns) | 6.7e-14 | 3.3e-06 (MagnusAdapt4) |
| ket 3 | 4.2e-15 | 9.5e-06 (MagnusAdapt4) |
| multi-ket 3 | 4.2e-15 | 3.0e-05 (MagnusAdapt4) |
| density 2 | 2.9e-15 | 7.9e-08 (Tsit5) |

Regression tests: `test/test_rollout_exact.jl`.
