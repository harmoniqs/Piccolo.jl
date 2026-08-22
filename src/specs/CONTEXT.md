# `Piccolo.Specs` — submodule context

`Piccolo.Specs` is the declarative problem-spec layer: it parses a versioned,
platform-neutral `ProblemSpec` (TOML on disk / JSON on the wire), **materializes**
it into a real `QuantumControlProblem` by *calling the existing template
functions*, solves it through one vetted generic runner, and emits both a
structured validation-error contract and a JSON Schema reflected from the
registries. It exists so that a CZ-on-transmons optimization can run end-to-end
from a single artifact with zero hand-authored Julia, and so that every run has a
durable, hashable identity.

This file orients a contributor to the **Phase 1 (wire-format-first)** cut. Read
it before touching anything under `src/specs/`.

## The wire-format-first boundary

Phase 1 delivers the *declarative spec + interpreter + error contract + schema
emission* over **today's function-templates**. `materialize` looks up a registry
factory and calls `SmoothPulseProblem(...)` / `SplinePulseProblem(...)` /
`SamplingProblem(...)` etc. directly. There is **no new problem type** — the
result is the same `QuantumControlProblem` you would build by hand, so the
existing template call surface is untouched (call-surface neutrality; see the
full-suite gate in Task 15).

What this buys us now: reliability (strict validation), iteration speed (one
runner, no per-problem scripts), durable artifacts (hashes + `result.toml`), and
a UI-facing schema. What it defers: the parametric-type rewrite (below).

## The registry entry is the single source of truth

Everything the interpreter and the schema need about a template/system/
integrator/wrapper/objective-term/solver/strategy/kind lives in one
`RegistryEntry` (`registries.jl`):

- `factory`  — the callable that materializes the object (an existing function).
- `params`   — hand-declared parameter schema (`name => (; type, default, range, doc)`).
- `compat`   — compatibility metadata (`:pulse_kinds`, `:trajectory_kinds`,
  `:ket_free_phase`, …).

`materialize` validates against `compat`; `emit_schema` (`schema.jl`) **reflects**
`params` + `compat` into the JSON Schema (enums, per-template param schemas,
conditional compatibility branches). Because the schema is emitted *from* the
registry, any registration change moves the schema — that is the wire-format-first
anti-drift property. (Plan 2 checks the emitted schema into the repo and adds the
CI diff gate + ajv / JSONSchema.jl fixture runs.)

Registration is idempotent with conflict detection: re-registering a
value-identical entry is a no-op; a conflicting entry under an existing name
errors. Piccolo self-registers at load (`register_all!()`); Piccolissimo augments
the same registries at load time.

## Deferred to Phase 1b (tracked, NOT here)

Phase 1b is the `@problem_template` parametric-type rewrite, gated on the
cloud-latency work actually starting. Doing it now would incur the
`isa QuantumControlProblem` migration blast radius before it buys anything. When
1b lands the macro, the macro *generates* the same `RegistryEntry`, so this
layer's anti-drift property is preserved unchanged. Explicitly out of Phase 1:

- `@problem_template` macro, **tagged aliases**, parametric **wrapper types**.
- `AbstractQuantumControlProblem` (there is no such supertype today; do not add
  one — the migration is 1b's).
- `extract_spec(problem)` round-tripped **against a materialized object** (needs a
  `spec` field retained on a typed problem). Phase 1 substitutes the
  **spec → canonical → parse → canonical** identity (see `run.jl`) plus the
  materialize golden tests; the runner retains the canonical spec in scope and
  writes it under `result.toml`'s `params.spec`.
- `PrecompileTools` `@compile_workload` wiring + sysimage (Phase 4 cloud).
  `schema.jl` ships a **callable** `precompile_workload(; tier1_only=true)` that
  exercises the parse → materialize → one-Ipopt-step path on tiny instances, but
  it is **not** wired to `@compile_workload` and PrecompileTools is not a dep.

Also reduced in Phase 1 (spec v1 lists them; not implemented here): `composite`
systems (`template` + `raw` only), the `robust` wrapper (register-for-schema-only
once Piccolissimo loads; `materialize` returns a structured "deferred" error),
the `tuning`/`compile` kinds, warm-start catalog resolution, and GPU ensemble
backends. The `altissimo` solver is **register-for-schema-only** — `solve_spec`
executes the `ipopt` backend only; `solver.backend = "altissimo"` raises a
structured `SpecValidationError`.

## Closed subset + the script-tier escape hatch

The spec covers a **closed subset** of problems: the `control` and `rollout`
kinds, the three base templates, the composition algebra (`goal_treatment`,
`[[objectives]]`, `free_dt`, the min-time recipe), the `sampling` wrapper, and the
`ipopt` solver. Anything outside the subset is a *strict* error — the parser
rejects unknown fields (with dotted field paths) and `materialize` rejects
incompatible combinations via `compat` queries, rather than silently doing
something surprising. The escape hatch for out-of-subset work is the ordinary
**script tier**: call the template functions directly in Julia. The spec layer is
additive; it never removes the hand-authored path.

## The two-hash carve-out

`hashes.jl` computes two identity keys from spec data (both `sha256` over
`canonical_json`):

- `structure_hash` — the **type-determining** subset only (system kind/template/
  levels, trajectory/pulse kinds, problem template + `goal_treatment` +
  `free_dt` *arm only* + `free_phase` + sorted objective kinds +
  `options.leakage_constraint`, integrator kind/alg, wrapper kinds, solver
  backend/device/precision/strategy). Deliberately **excludes** $N$, $T$, $Q$, and
  regularizer weights, so a resize/reweight preserves `structure_hash`.
- `problem_hash` — the whole canonical wire form; changes when any wire field
  changes.

`canonical_json` is hand-rolled to be **byte-identical** across Julia and
TypeScript (Plan 2's `hashing.ts`): sorted keys, JSON-escaped strings, and an
int/float-agnostic numeric rule (integer-valued numbers render as bare integers;
non-integers use the ECMAScript `Number::toString` algorithm). Do **not** swap in
Julia's native `repr`/Ryu float formatting — it does not match JS.

## House rules baked in

- **Save the pulse, never the trajectory.** The runner writes `extract_pulse(...)`
  to `pulse-<hash>.jld2`, not the `NamedTrajectory`.
- **Globals are pinned by default.** `TemplateBlock.global_names` defaults to
  `Symbol[]`; globals are only co-optimized when explicitly listed.
- **`T` is scalar duration, `N` is knot count** — the field names follow this.
- **Never construct/modify `BilinearIntegrator`.** The `bilinear` integrator entry
  is a *sentinel* meaning "pass `integrator = nothing`" so the template builds its
  own default; we never build one ourselves.

## OSS vs full schema variant

`emit_schema` reflects the *live* registries, so its output depends on what is
loaded:

- **OSS variant** (Piccolo only, Piccolissimo NOT loaded): Piccolo-only names —
  integrator `bilinear`; wrapper `sampling`; solver `ipopt`; strategy `direct`;
  the six Piccolo objective terms; seven systems; three templates.
- **Full variant** (Piccolissimo loaded): the same registries gain the
  `exponential`/`spline` integrators, the `hermite_*` objective terms, the
  `robust` wrapper, and the `altissimo` solver — with **no code change** in
  `schema.jl`. Emitting the full schema therefore *requires Piccolissimo to be
  loaded* in the emitting process.

Note the `free_phase ∨ globals → integrator.kind ∈ {exponential, spline}`
conditional is emitted with the literal non-bilinear integrator families even in
the OSS variant (matching the `materialize` rule and its `SpecError` message), so
the OSS schema correctly states that free-phase requires an integrator that only
becomes available once Piccolissimo is loaded.

`emit_schema` uses **per-branch / per-object** `additionalProperties:false`, never
a single top-level one. This is deliberate: in draft-07 a top-level
`additionalProperties:false` does not see properties introduced by a sibling
`if`/`then`/`oneOf`, so combining strict-unknown-field rejection with the
conditional-compat branches would reject *valid* specs. The conditional branches
only *narrow enums* on already-declared properties (`pulse.kind`,
`integrator.kind`) — they never introduce a property.

## Flagged for Plan 2 (schema owner + amicode)

Two assumptions here should be confirmed on the amicode/Plan-2 side:

1. **`result/v1` `additionalProperties` on `[params]`.** The runner writes
   `result.toml` (the amicode `result/v1` artifact) with `structure_hash`,
   `problem_hash`, `converged`, `spec`, and the `Pkg` `versions` map all under the
   `[params]` sub-table — assumed to be the only `additionalProperties:true` table
   in `result/v1`, so this is schema-valid today with no cross-repo change. Plan 3's
   learning-loops ledger reads `params.versions`. Plan 2 may later promote versions
   to a top-level `[versions]` block via a `result` schema v2 bump; doing it
   top-level now would fail the existing `amico-validate` CI gate.
2. **`AMICODE_PULSE_META` reconciliation.** The runner emits, around the solve,
   `AMICODE_PULSE_META n_drives=… knots=… pulse=…` for the extension's run-dir
   pulse-preview reader. But `amico-run`'s `telemetry.ts classifyLine` treats it as
   a plain log line, and the *real* extension's pulse-preview line is prefix-less.
   The `AMICODE_ITER …` progress lines and the terminal `DONE …` line **are** parsed
   by `classifyLine` and are the load-bearing contract. The `AMICODE_PULSE_META`
   line is forward-looking; reconciling its shape with the real extension is a
   Plan 2 / amicode task. **Do not change the runner's sentinels to chase this** —
   flag it, do not fix it here.

## File map

| File | Responsibility |
|---|---|
| `_specs.jl` | `module Specs`; includes + reexports; `register_all!()` at load. |
| `spec_structs.jl` | Plain-data structs (`ProblemSpec`, …) + the `FreeDt` sum type. |
| `registries.jl` | `SpecRegistries`: the eight registry globals, register/lookup API, `RegistryEntry`, Piccolo's `register_all!`. |
| `errors.jl` | `SpecError`, `SpecValidationError`, `emit_error_json` (exit-2 contract). |
| `parse.jl` | Strict TOML/JSON → structs; unknown-field rejection with field paths; kind dispatch. |
| `hashes.jl` | `canonical_json`, `full_dict`, `structure_fields`, `structure_hash`, `problem_hash`. |
| `materialize.jl` | `materialize(::ProblemSpec)` → registry lookup → existing template call + composition axes + wrappers + trait validations. |
| `rollout_kind.jl` | `rollout` materialize/run; `referee_rollout`; independence checks; `Verdict`. |
| `run.jl` | `solve_spec`/`run_spec` generic runner; `result.toml`; `AMICODE_ITER`/`DONE`/`AMICODE_PULSE_META` sentinels; JLD2 pulse save. |
| `schema.jl` | `emit_schema()` (draft-07 discriminated union from registries) + `precompile_workload()`. |
