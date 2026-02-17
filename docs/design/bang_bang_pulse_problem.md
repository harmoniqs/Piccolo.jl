# BangBangPulseProblem — Design Document

## Overview

A new problem template that promotes **bang-bang** (piecewise-constant, few-switch) controls by penalizing ‖du‖₁ via **exact slack reformulation**. Unlike `SmoothPulseProblem` (2 derivative levels + L2 regularization), this stores only 1 derivative (`du`) and uses slack variables to impose an exact L1 penalty on it.

### Comparison with SmoothPulseProblem

| | SmoothPulseProblem | BangBangPulseProblem |
|---|---|---|
| Derivatives stored | `du`, `ddu` | `du` only |
| `add_control_derivatives` | `n=2` | `n=1` |
| `DerivativeIntegrator`s | 2 (u→du, du→ddu) | 1 (u→du) |
| Regularization on `du` | `QuadraticRegularizer` (L2) | `LinearRegularizer` on slacks (L1) |
| Regularization on `u` | `QuadraticRegularizer` (L2) | same |
| Extra variables | — | slack `s_du ≥ 0` |
| Extra constraints | — | `L1SlackConstraint`: \|du\| ≤ s |
| Bound params | `du_bound`, `ddu_bound` | `du_bound` only |

---

## L1 Penalty via Slack Variables

Instead of a smooth approximation √(x² + ε²), we use an **exact** slack reformulation:

**Introduce** slack variables `s ≥ 0` (same dimension as `du`, at every timestep) and enforce:

$$
du_{k,i} \leq s_{k,i}, \quad -du_{k,i} \leq s_{k,i} \quad \Longleftrightarrow \quad |du_{k,i}| \leq s_{k,i}
$$

**Then minimize** the linear cost on slacks:

$$
J_{\text{L1}} = R \sum_{k} \sum_i s_{k,i} \cdot \Delta t_k
$$

At optimality, `s = |du|`, so `J_L1 = R · ‖du‖₁` (weighted by timestep).

### Why slacks over smooth approximation?

| | Smooth approx √(x² + ε²) | Slack variables |
|---|---|---|
| **Exactness** | Approximate (bias ~ε) | Exact L1 |
| **Extra variables** | None | `n_drives × N` slacks |
| **Extra constraints** | None | `2 × n_drives × N` linear ineq |
| **IPOPT behavior** | Smooth Hessian everywhere | Linear obj + constraints (IPOPT loves this) |
| **Tuning** | Need to pick ε | Nothing to tune |
| **Complexity** | 1 new objective | 1 new constraint type + 1 objective + trajectory setup |

---

## Implementation Plan

### Phase 1: DirectTrajOpt.jl — New Objective & Constraint

#### 1A. `L1SlackConstraint` — new linear constraint

**File**: `DirectTrajOpt.jl/src/constraints/linear/l1_slack_constraint.jl` (NEW)

```julia
struct L1SlackConstraint <: AbstractLinearConstraint
    var_name::Symbol     # variable to penalize (e.g. :du)
    slack_name::Symbol   # slack variable name (e.g. :s_du)
    label::String
end
```

The `constrain!` method (in `DirectTrajOpt.jl/src/solvers/ipopt_solver/constraints.jl`) adds for each timestep k and each component i:

$$
du_{k,i} - s_{k,i} \leq 0 \quad \text{and} \quad -du_{k,i} - s_{k,i} \leq 0
$$

(i.e., `|du| ≤ s`)

The bounds `s ≥ 0` come from the trajectory's bounds on the slack component (set in Phase 2).

This is purely linear — no Jacobian/Hessian of Lagrangian needed. Follows the same pattern as `EqualityConstraint`, `BoundsConstraint`, etc.

**Wire-up**:
- Include in `_constraints.jl`: `include("linear/l1_slack_constraint.jl")`
- Add `export L1SlackConstraint`
- Add the `constrain!` dispatch in `ipopt_solver/constraints.jl` (replacing the commented-out code)

#### 1B. `LinearRegularizer` — new objective

**File**: Append to `DirectTrajOpt.jl/src/objectives/regularizers.jl`

```julia
struct LinearRegularizer <: AbstractObjective
    name::Symbol          # variable to penalize (e.g. :s_du)
    R::Vector{Float64}    # per-component weights
    times::Vector{Int}    # timesteps
end
```

Computes:

$$
J = \sum_{k \in \text{times}} \sum_i R_i \cdot v_{k,i} \cdot \Delta t_k
$$

- **Gradient**: ∂J/∂v_{k,i} = R_i · Δt_k, ∂J/∂Δt_k = Σᵢ R_i · v_{k,i}
- **Hessian**: only cross-terms ∂²J/∂v_{k,i}∂Δt_k = R_i (no diagonal)
- Much simpler than `QuadraticRegularizer`

**Wire-up**:
- Add `export LinearRegularizer` to `_objectives.jl`
- Add `@testitem` in the same file

#### 1C. Update `remove_slack_variables!`

In `DirectTrajOpt.jl/src/solvers/ipopt_solver/solver.jl`, the existing function already looks for `L1SlackConstraint` — it just needs the struct to exist. Verify the `slack_names` field name matches or update accordingly (the old code uses `con.slack_names`, we'll use `con.slack_name` as a single Symbol → wrap in vector).

---

### Phase 2: Piccolo.jl — `BangBangPulseProblem`

#### 2A. New problem template file

**File**: `Piccolo.jl/src/control/templates/bang_bang_pulse_problem.jl` (NEW)

Three methods mirroring `SmoothPulseProblem`:

**Method 1**: `BangBangPulseProblem(qtraj::AbstractQuantumTrajectory{<:ZeroOrderPulse}, N::Int; ...)`

Handles `UnitaryTrajectory`, `KetTrajectory`, `DensityTrajectory`.

| Keyword | Default | Purpose |
|---|---|---|
| `integrator` | `nothing` | Custom integrator(s) |
| `global_names` | `nothing` | Global variables |
| `global_bounds` | `nothing` | Global variable bounds |
| `du_bound` | `Inf` | Bound on first derivative |
| `Δt_bounds` | `nothing` | Timestep bounds |
| `Q` | `100.0` | Infidelity weight |
| `R` | `1e-2` | Default regularization weight |
| `R_u` | `R` | L2 weight on control amplitude |
| `R_du` | `R` | L1 weight on `du` (applied to slacks) |
| `constraints` | `[]` | Additional constraints |
| `piccolo_options` | `PiccoloOptions()` | Options |

**Steps inside the constructor** (paralleling `SmoothPulseProblem`):

1. Extract `sys`, `state_sym`, `control_sym` from `qtraj`
2. Build `global_data` from `sys.global_params` if present
3. Convert `qtraj` → `base_traj` via `NamedTrajectory(qtraj, N; ...)`
4. **Add 1 control derivative**: `add_control_derivatives(base_traj, 1; control_name=control_sym, derivative_bounds=(du_bounds,))`
5. **Add slack component** `:s_du` to the trajectory:
   - Data: `abs.(traj_bb[:du])` (initialize at |du| to be feasible)
   - Dimension: same as `du` (`n_drives × N`)
   - Bounds: `(zeros(n_drives), fill(Inf, n_drives))` (non-negative)
   - Registered as a control (so it's optimized)
6. Build integrators:
   - Dynamics: `BilinearIntegrator(qtraj, N)` (or custom)
   - 1 `DerivativeIntegrator(:u, :du, traj)` (instead of 2)
7. Build objective:
   - `_state_objective(qtraj, traj, state_sym, Q)` — same helper as `SmoothPulseProblem`
   - `QuadraticRegularizer(:u, traj, R_u)` — L2 on control amplitude
   - `LinearRegularizer(:s_du, traj, R_du)` — **L1 on switching** (via slacks)
   - `_apply_piccolo_options(...)` — same helper
8. Build constraints:
   - `L1SlackConstraint(:du, :s_du)` — ties slacks to |du|
   - Plus user-provided constraints + global bounds
9. Return `QuantumControlProblem(qtraj, DirectTrajOptProblem(...))`

**Method 2**: `BangBangPulseProblem(qtraj::MultiKetTrajectory{<:ZeroOrderPulse}, N::Int; ...)`

Same pattern but uses `_ensemble_ket_objective` and `BilinearIntegrator` for ensemble. Identical structural differences as `SmoothPulseProblem`'s ensemble method.

**Method 3**: Fallback error for non-`ZeroOrderPulse` types.

#### 2B. Wire into module

In `_problem_templates.jl`, add:
```julia
include("bang_bang_pulse_problem.jl")
```

No changes needed to `_control.jl` or `Piccolo.jl` — exports propagate via `@reexport`.

---

### Phase 3: Helper for adding slack components

The slack variable needs to be added to the `NamedTrajectory` after `add_control_derivatives`. Two approaches:

**Option A** — Inline in `BangBangPulseProblem`: Reconstruct the trajectory with an added component using the same `NamedTrajectory(comps_data, global_data; ...)` pattern that `add_control_derivatives` uses internally. Extract all existing components, append `:s_du`, rebuild.

**Option B** — New utility `add_slack_variable(traj, :du, :s_du)` in `NamedTrajectories.jl` that generalizes the pattern.

**Recommendation**: Option A first (keep it self-contained), refactor to Option B later if reused.

---

### Phase 4: Tests

#### 4A. Unit tests in DirectTrajOpt.jl

- **`L1SlackConstraint`**: `@testitem` in the constraint file — build a trajectory with `:du` and `:s_du`, create the constraint, solve a simple problem, verify `s ≥ |du|` at solution
- **`LinearRegularizer`**: `@testitem` using `test_objective` (finite-diff gradient/Hessian validation) — same pattern as `QuadraticRegularizer` tests

#### 4B. Integration tests in Piccolo.jl

- `@testitem` in `bang_bang_pulse_problem.jl` — same structure as existing `SmoothPulseProblem` tests:
  - Construct a simple 2-level system
  - Create `ZeroOrderPulse` + `UnitaryTrajectory`
  - Solve `BangBangPulseProblem` for ~50 iterations
  - Verify infidelity decreases
  - Verify `du` is sparse (many near-zero entries)
  - Optionally test `MultiKetTrajectory` variant

---

## File Change Summary

| File | Action | LOC est. |
|---|---|---|
| `DirectTrajOpt.jl/src/constraints/linear/l1_slack_constraint.jl` | **New** — `L1SlackConstraint` struct | ~40 |
| `DirectTrajOpt.jl/src/constraints/_constraints.jl` | Add `include` + `export` | ~2 |
| `DirectTrajOpt.jl/src/solvers/ipopt_solver/constraints.jl` | Add `constrain!` dispatch for `L1SlackConstraint` (replace commented code) | ~25 |
| `DirectTrajOpt.jl/src/solvers/ipopt_solver/solver.jl` | Update `remove_slack_variables!` field name | ~3 |
| `DirectTrajOpt.jl/src/objectives/regularizers.jl` | **Append** — `LinearRegularizer` struct + interface | ~120 |
| `DirectTrajOpt.jl/src/objectives/_objectives.jl` | Add `export LinearRegularizer` | ~1 |
| `Piccolo.jl/src/control/templates/bang_bang_pulse_problem.jl` | **New** — 3 methods + helpers reuse + tests | ~350 |
| `Piccolo.jl/src/control/templates/_problem_templates.jl` | Add `include(...)` | ~1 |

**Total**: ~540 lines across 8 files (3 new, 5 modified).

---

## Dependency Order

```
1. LinearRegularizer (objective, no deps)
2. L1SlackConstraint (constraint struct, no deps)
   ↓
3. IPOPT constrain! dispatch (needs L1SlackConstraint)
4. remove_slack_variables! update (needs L1SlackConstraint)
   ↓
5. BangBangPulseProblem (needs LinearRegularizer + L1SlackConstraint)
   ↓
6. Tests
```

Steps 1–2 are independent and can be done in parallel. Steps 3–4 depend on 2. Step 5 depends on 1–4. Step 6 depends on 5.
