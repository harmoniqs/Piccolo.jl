# Modulated Drives Documentation — Design Spec

**Date:** 2026-04-03
**Status:** Approved
**Scope:** Add `## Modulated Drives` section to `docs/literate/concepts/systems.jl`

---

## Background

PR #132 (`c0dbe9b`) introduced `DriftTerm`, `ModulatedDrive`, and a Pair-based
`QuantumSystem` constructor for structured time modulation of drift and drive terms.
No documentation page was added. The feature is only incidentally mentioned in the
trapped ions system page.

---

## Decision: Extend Existing Page

Add a `## Modulated Drives` section to `docs/literate/concepts/systems.jl`, inserted
**after** the existing `## Nonlinear Drives` section (currently ending around line 209).

**Why not a new page?**
`ModulatedDrive` / `DriftTerm` are constructor-level features of `QuantumSystem`, not a
separate problem template or system type. All drive construction knowledge lives in the
systems concept page. A standalone page would be thin and disconnected.

---

## Section Content

### Order within the file

```
## The Hamiltonian Model
## QuantumSystem
## Drive Bounds
## Accessing System Properties
## Nonlinear Drives          ← existing
## Modulated Drives          ← NEW (inserted here)
## OpenQuantumSystem
## CompositeQuantumSystem
...
```

### Subsections

#### 1. Motivation (prose)

Explain the distinction from `NonlinearDrive`:
- `NonlinearDrive`: coefficient is an *unknown function of the controls* — to be optimized
- `ModulatedDrive` / `DriftTerm`: coefficient is a *known function of time* — not optimized, just evaluated

Use cases: rotating-frame oscillation `cos(ωt)`, Floquet modulation, sideband envelopes,
AC Stark shifts, any term whose time dependence is analytically known.

#### 2. Drive types table extension

Extend the drive type table to include:

| Type | Coefficient | Use case |
|------|-------------|----------|
| `ModulatedDrive(base, t->f(t))` | `f(t) · c_base(u)` | Rotating frame, sideband drives |
| `DriftTerm(H, t->f(t))` | `f(t)` (drift only) | Time-varying drift (AC Stark, Floquet) |

#### 3. Pair-based constructor syntax

Show the sugar syntax:

```julia
# Modulated drive channel: H_x oscillates at ω
sys = QuantumSystem(H_drift, [H_x => t -> cos(ω * t), H_y], [1.0, 1.0])

# Modulated drift: AC Stark shift on the drift term
sys = QuantumSystem(H_z => t -> cos(ω * t), [H_x, H_y], [1.0, 1.0])

# Mixed: static + modulated drift terms
sys = QuantumSystem([H_static, H_ac => t -> cos(ω * t)], [H_drive], [1.0])
```

Also show the explicit `ModulatedDrive` constructor for when you need to wrap a
`NonlinearDrive`:

```julia
nd = NonlinearDrive(PAULIS[:Z], u -> u[1]^2 + u[2]^2)
sys = QuantumSystem(H_drift, [nd => t -> cos(ω * t)], [1.0, 1.0])
```

#### 4. Auto-detection of `time_dependent`

Note that Piccolo automatically detects the presence of modulation and switches to
`TimeDependentBilinearIntegrator`. The user does not configure this manually.

```julia
QuantumSystem(H_drift, [H_x], [1.0]).time_dependent          # false
QuantumSystem(H_drift, [H_x => t -> cos(t)], [1.0]).time_dependent  # true
```

#### 5. Physical example

A transmon in the lab frame where the drive is explicitly `cos(ωt) σ_x` (before
making the rotating-wave approximation), showing how to set it up as a modulated drive.
This grounds the feature in a concrete, recognizable use case.

---

## Out of Scope

- A separate How-To guide page (can be added later if needed)
- Documentation of `drive_coeff_dt` as a public API (internal; not called by users)
- Derivation of rotating-wave approximation

---

## Files Changed

| File | Change |
|------|--------|
| `docs/literate/concepts/systems.jl` | Add `## Modulated Drives` section |

No changes to `make.jl` (the page is already in the nav).
