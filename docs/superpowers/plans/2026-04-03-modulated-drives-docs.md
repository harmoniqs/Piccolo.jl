# Modulated Drives Docs Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `## Modulated Drives` section to `docs/literate/concepts/systems.jl` documenting `DriftTerm`, `ModulatedDrive`, and the Pair-based `QuantumSystem` constructor introduced in PR #132.

**Architecture:** Single-file edit — insert a new section between `## Nonlinear Drives` and `## OpenQuantumSystem` in the existing literate concepts page. No new files, no `make.jl` changes.

**Tech Stack:** Julia literate docs (Literate.jl), Documenter.jl. Code blocks in the literate file are executable Julia; comments become Markdown prose.

---

## File Map

| File | Change |
|------|--------|
| `docs/literate/concepts/systems.jl` | Insert `## Modulated Drives` section after line 209 |

---

### Task 1: Insert the `## Modulated Drives` section

**Files:**
- Modify: `docs/literate/concepts/systems.jl` (after line 209, before `# ## OpenQuantumSystem`)

The section covers: motivation, drive types table, Pair syntax, `time_dependent` auto-detection, typed-drive wrapping, physical example, and cross-links.

- [ ] **Step 1: Insert the section**

In `docs/literate/concepts/systems.jl`, find the line:

```julia
sys_linear.drives  ## auto-populated LinearDrives
```

Add the following block **immediately after** that line (before the blank line and `# ## OpenQuantumSystem`):

```julia
# ## Modulated Drives
#
# When the coefficient of a drive or drift term is a **known function of time** —
# rather than an unknown control amplitude to be optimized — use `ModulatedDrive`
# or a modulated `DriftTerm`.
#
# Common use cases:
# - Rotating-frame oscillation: `cos(ωt)` carrier on a drive channel
# - Floquet modulation: periodic envelope on drift or drive terms
# - AC Stark shifts: time-varying correction to the drift Hamiltonian
# - Sideband envelopes: analytically known time-dependent coupling profile
#
# This is distinct from [`NonlinearDrive`](@ref): a `NonlinearDrive` coefficient
# is an *unknown function of the controls* that the optimizer finds. Here the
# time dependence is specified analytically — Piccolo evaluates it at each
# timestep, not optimizes over it.
#
# ### Drive Types
#
# | Type | Coefficient | Use case |
# |------|-------------|----------|
# | `LinearDrive(H, i)` | ``u_i`` | Standard bilinear control |
# | `NonlinearDrive(H, f)` | ``f(\boldsymbol{u})`` | Displaced frames, cross-Kerr |
# | `ModulatedDrive(base, t->b(t))` | ``b(t) \cdot c_{\text{base}}(\boldsymbol{u})`` | Rotating frame, sideband drives |
# | `DriftTerm(H, t->b(t))` | ``b(t)`` (drift only) | Time-varying drift (AC Stark, Floquet) |
#
# ### Pair-Based Constructor Syntax
#
# The simplest way to build a modulated system is `H => t -> b(t)` in the drives
# vector, or the same for the drift argument:

omega = 2π * 0.1  # modulation frequency (GHz)

## Modulated drive channel: H_x oscillates at ω
sys_mod = QuantumSystem(PAULIS[:Z], [PAULIS[:X] => t -> cos(omega * t), PAULIS[:Y]], [1.0, 1.0])
sys_mod.time_dependent

#-

## Modulated drift: time-varying drift term
sys_mod_drift = QuantumSystem(PAULIS[:Z] => t -> cos(omega * t), [PAULIS[:X], PAULIS[:Y]], [1.0, 1.0])
sys_mod_drift.time_dependent

#-

## Mixed: static + modulated drift terms as a vector
sys_mixed = QuantumSystem([PAULIS[:Z], 0.05 * PAULIS[:X] => t -> cos(omega * t)], [PAULIS[:Y]], [1.0])
length(sys_mixed.drift_terms)

# When any modulation is present, Piccolo automatically sets `time_dependent = true`
# and dispatches `TimeDependentBilinearIntegrator` — no manual configuration needed.
#
# ### Wrapping Typed Drives
#
# `NonlinearDrive` can also be modulated using the same `Pair` syntax:

nd_mod = NonlinearDrive(PAULIS[:Z], u -> u[1]^2 + u[2]^2)
sys_nl_mod = QuantumSystem(PAULIS[:Z], [nd_mod => t -> cos(omega * t)], [1.0, 1.0])
sys_nl_mod.H_drives[1]

# Or wrap explicitly: `ModulatedDrive(nd_mod, t -> cos(omega * t))`.
#
# ### Physical Example: Lab-Frame Drive
#
# In the lab frame a microwave drive on a transmon takes the form
# ``\Omega(t)\cos(\omega_d t)\,\sigma_x``. Before applying the rotating-wave
# approximation you can model it directly:
#
# ```math
# H(t) = \tfrac{\omega_q}{2}\,\sigma_z
#        + u(t)\,\cos(\omega_d t)\,\sigma_x
# ```
#
# where ``u(t)`` is the envelope amplitude the optimizer finds and
# ``\cos(\omega_d t)`` is the known carrier:

omega_q = 5.0   # qubit frequency (GHz)
omega_d = 5.0   # drive frequency, on resonance

sys_lab = QuantumSystem(
    (omega_q / 2) * PAULIS[:Z],
    [PAULIS[:X] => t -> cos(omega_d * t)],
    [0.5],
)
sys_lab.time_dependent

# The optimizer finds the envelope ``u(t)`` while the carrier is evaluated
# exactly at each timestep, enabling sub-RWA modeling for short high-power gates.
#
# ## See Also
#
# - [Nonlinear Drives](@ref) — coefficient is an unknown function of the controls
# - [Trapped Ions](@ref trapped-ion-systems) — `RadialMSGateSystem` uses time-dependent
#   Hamiltonians as a canonical physical example

```

- [ ] **Step 2: Verify the section looks right in context**

Open `docs/literate/concepts/systems.jl` and confirm:
- The new section appears between `sys_linear.drives  ## auto-populated LinearDrives` and `# ## OpenQuantumSystem`
- All variable names (`omega`, `sys_mod`, `sys_mod_drift`, `sys_mixed`, `nd_mod`, `sys_nl_mod`, `sys_lab`) are defined before use and do not shadow earlier locals
- The `omega_q` / `omega_d` variables defined here do not conflict with any earlier definitions in the file (search for `omega_q` above the insertion point — it should not exist)

- [ ] **Step 3: Smoke-test the new code blocks**

Run in a Julia REPL (from the repo root) to confirm no errors:

```julia
using Piccolo

omega = 2π * 0.1

sys_mod = QuantumSystem(PAULIS[:Z], [PAULIS[:X] => t -> cos(omega * t), PAULIS[:Y]], [1.0, 1.0])
@assert sys_mod.time_dependent

sys_mod_drift = QuantumSystem(PAULIS[:Z] => t -> cos(omega * t), [PAULIS[:X], PAULIS[:Y]], [1.0, 1.0])
@assert sys_mod_drift.time_dependent

sys_mixed = QuantumSystem([PAULIS[:Z], 0.05 * PAULIS[:X] => t -> cos(omega * t)], [PAULIS[:Y]], [1.0])
@assert length(sys_mixed.drift_terms) == 2

nd_mod = NonlinearDrive(PAULIS[:Z], u -> u[1]^2 + u[2]^2)
sys_nl_mod = QuantumSystem(PAULIS[:Z], [nd_mod => t -> cos(omega * t)], [1.0, 1.0])
@assert sys_nl_mod.H_drives[1] isa ModulatedDrive

sys_lab = QuantumSystem((5.0 / 2) * PAULIS[:Z], [PAULIS[:X] => t -> cos(5.0 * t)], [0.5])
@assert sys_lab.time_dependent

println("All assertions passed.")
```

Expected output: `All assertions passed.`

- [ ] **Step 4: Commit**

```bash
git add docs/literate/concepts/systems.jl
git commit -m "docs: add Modulated Drives section to systems concept page"
```
