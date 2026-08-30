# ============================================================================ #
# Algorithm-Specific Data Types
#
# Pre-allocated buffers and ODE problems for each integration algorithm.
# These are stored in the `alg_data` field of SplineIntegrator.
# ============================================================================ #

export Tsit5Data, MagnusGL4Data, MagnusAdapt4Data, Rodas5PData
export ChebyshevData, expl_discretization_error, suggest_n_sub

"""
    Tsit5Data{R<:Number}

Pre-built ODE problems and sparsity structure for Tsit5-based forward propagation.

# Fields
- `Φ_probs::Vector{ODEProblem}`: One ODE problem per interval, remade with new parameters
- `Φ_structure::SparseMatrixCSC{Float64,Int}`: Sparsity structure of propagator
- `jvp_probs::Union{Nothing,Vector{ODEProblem}}`: Phase 2 ket-level JVP ODE
  problems — one per interval, state `[ψ; Y] ∈ ℂ^{2d}`. Populated only for
  `KetTrajectory` integrators built via the standard constructor; `nothing`
  for trajectory variants that don't yet use the ket-level JVP path
  (e.g., `MultiKetTrajectory`, `MultiDensityTrajectory`).
- `vjp_probs::Union{Nothing,Vector{ODEProblem}}`: Phase 3 ket-level adjoint
  VJP ODE problems — one per interval, state `[ψ; λ; g] ∈ ℂ^{2d + n_p}`,
  integrated backward (`tspan = (1.0, 0.0)`). Populated only for
  `KetTrajectory` integrators built via the standard constructor; `nothing`
  otherwise.
- `hvp_fwd_probs::Union{Nothing,Vector{ODEProblem}}`: Phase 4 ket-level HVP
  forward-pass ODE problems — one per interval, state `[ψ; δψ] ∈ ℂ^{2d}`.
  Structurally identical to `jvp_probs` (both use `build_ket_jvp_ode`-style
  dynamics) but kept as a separate slot to allow future specialization
  (e.g., dense-output retention) without disturbing Phase 2's JVP path.
- `hvp_bwd_probs::Union{Nothing,Vector{ODEProblem}}`: Phase 4 ket-level
  second-order adjoint ODE problems — one per interval, state
  `[ψ'; δψ'; λ; ν; g] ∈ ℂ^{4d + n_p}`, integrated backward
  (`tspan = (1.0, 0.0)`). Includes the second-derivative-of-coefficient term
  ∂²H/∂u_i∂u_j that Gauss-Newton drops.
- `hvp_solvers::Union{Nothing,Vector{Any}}`: **#342** — one PREALLOCATED
  integrator bundle per interval for the HVP's forward and backward solves, built
  lazily on the interval's first `step_hvp!` and reused across every subsequent
  product. `nothing` for trajectory variants with no ket-level HVP; otherwise a
  `Vector` of length `N-1` whose slots start as `nothing`. Element type is `Any`
  because a bundle carries two `ODEIntegrator`s whose parameterisation depends on
  the per-interval RHS closure. **#356**: the bundle is therefore UNPACKED at the
  call site into concretely typed buffers, and the knot body is split into three
  concrete calls with only the two solves left `Any`-dispatched — a single call
  taking the bundle as an `Any` made that call dynamic and boxed the driver's
  per-knot `KnotColumns` and views at 144 B/knot (see `step_hvp!`). The bundle also
  caches the interval's initial ODE step, determined once per iterate. Slot `k` is
  touched only by the task that owns knot `k`, so this is per-INTERVAL state — never
  `threadid()`-keyed.
"""
struct Tsit5Data{R<:Number}
    Φ_probs::Vector{ODEProblem}
    Φ_structure::SparseMatrixCSC{Float64,Int}
    jvp_probs::Union{Nothing,Vector{ODEProblem}}
    vjp_probs::Union{Nothing,Vector{ODEProblem}}
    hvp_fwd_probs::Union{Nothing,Vector{ODEProblem}}
    hvp_bwd_probs::Union{Nothing,Vector{ODEProblem}}
    hvp_solvers::Union{Nothing,Vector{Any}}
end

# #342: the pre-#342 six-field form. Slot vector is derived from
# `hvp_bwd_probs`, so every existing caller keeps working and no construction
# site has to know about the solver cache.
Tsit5Data{R}(
    Φ_probs::AbstractVector{<:ODEProblem},
    Φ_structure::SparseMatrixCSC{Float64,Int},
    jvp_probs::Union{Nothing,AbstractVector{<:ODEProblem}},
    vjp_probs::Union{Nothing,AbstractVector{<:ODEProblem}},
    hvp_fwd_probs::Union{Nothing,AbstractVector{<:ODEProblem}},
    hvp_bwd_probs::Union{Nothing,AbstractVector{<:ODEProblem}},
) where {R} = Tsit5Data{R}(
    Vector{ODEProblem}(Φ_probs),
    Φ_structure,
    isnothing(jvp_probs) ? nothing : Vector{ODEProblem}(jvp_probs),
    isnothing(vjp_probs) ? nothing : Vector{ODEProblem}(vjp_probs),
    isnothing(hvp_fwd_probs) ? nothing : Vector{ODEProblem}(hvp_fwd_probs),
    isnothing(hvp_bwd_probs) ? nothing : Vector{ODEProblem}(hvp_bwd_probs),
    (isnothing(hvp_fwd_probs) || isnothing(hvp_bwd_probs)) ? nothing :
    Vector{Any}(nothing, length(hvp_bwd_probs)),
)

# Backwards-compatible constructors:
# - existing callers without Phase 2/3/4 problems can omit them
# - Phase 2 callers may pass jvp_probs only
# - Phase 3 callers may pass jvp_probs + vjp_probs
# - Phase 4 callers pass all four matrix-free problem vectors
# Accept any concrete `Vector{<:ODEProblem}` shape — the inner ODEProblem
# parameterization differs per trajectory variant (Matrix-state for unitary,
# vector-state for ket, etc.).
Tsit5Data{R}(
    Φ_probs::AbstractVector{<:ODEProblem},
    Φ_structure::SparseMatrixCSC{Float64,Int},
) where {R} = Tsit5Data{R}(
    Vector{ODEProblem}(Φ_probs),
    Φ_structure,
    nothing,
    nothing,
    nothing,
    nothing,
)

Tsit5Data{R}(
    Φ_probs::AbstractVector{<:ODEProblem},
    Φ_structure::SparseMatrixCSC{Float64,Int},
    jvp_probs::Union{Nothing,AbstractVector{<:ODEProblem}},
) where {R} = Tsit5Data{R}(
    Vector{ODEProblem}(Φ_probs),
    Φ_structure,
    isnothing(jvp_probs) ? nothing : Vector{ODEProblem}(jvp_probs),
    nothing,
    nothing,
    nothing,
)

Tsit5Data{R}(
    Φ_probs::AbstractVector{<:ODEProblem},
    Φ_structure::SparseMatrixCSC{Float64,Int},
    jvp_probs::Union{Nothing,AbstractVector{<:ODEProblem}},
    vjp_probs::Union{Nothing,AbstractVector{<:ODEProblem}},
) where {R} = Tsit5Data{R}(
    Vector{ODEProblem}(Φ_probs),
    Φ_structure,
    isnothing(jvp_probs) ? nothing : Vector{ODEProblem}(jvp_probs),
    isnothing(vjp_probs) ? nothing : Vector{ODEProblem}(vjp_probs),
    nothing,
    nothing,
)

"""
    MagnusGL4Data

Pre-allocated dense buffers for MagnusGL4 forward propagation.
Thread-safe: one set of buffers per interval.

# Fields
- `G_drift::SparseMatrixCSC{Float64,Int}`: Reference drift generator (isomorphic real form)
- `G_drives::Vector{SparseMatrixCSC{Float64,Int}}`: Reference drive generators
- `G_drift_copies::Vector{SparseMatrixCSC{Float64,Int}}`: Per-interval drift copies (threading)
- `G_drives_copies::Vector{Vector{SparseMatrixCSC{Float64,Int}}}`: Per-interval drive copies
- `G_buffers::Vector{Matrix{Float64}}`: Dense buffers for exp(Ω) computation
- `drives::Vector{AbstractDrive}`: Drive terms for computing coefficients (supports nonlinear drives)
"""
struct MagnusGL4Data
    G_drift::SparseMatrixCSC{Float64,Int}
    G_drives::Vector{SparseMatrixCSC{Float64,Int}}
    G_drift_copies::Vector{SparseMatrixCSC{Float64,Int}}
    G_drives_copies::Vector{Vector{SparseMatrixCSC{Float64,Int}}}
    G_buffers::Vector{Matrix{Float64}}
    drives::Vector{AbstractDrive}
end

"""
    MagnusAdapt4Data

Pre-allocated dense buffers for MagnusAdapt4 forward propagation.
Same structure as MagnusGL4Data — the adaptive solver uses the same operator construction.
"""
const MagnusAdapt4Data = MagnusGL4Data

"""
    Rodas5PData{R<:Number}

Pre-built ODE problems and sparsity structure for Rodas5P-based forward propagation.
Mirrors `Tsit5Data` shape so it can be a drop-in peer at construction sites.

# Fields
- `Φ_probs::Vector{ODEProblem}`: One ODE problem per interval, remade with new parameters
- `Φ_structure::SparseMatrixCSC{Float64,Int}`: Sparsity structure of propagator (dense placeholder for stiff path)
"""
struct Rodas5PData{R<:Number}
    Φ_probs::Vector{ODEProblem}
    Φ_structure::SparseMatrixCSC{Float64,Int}
end

# ============================================================================ #

# ─────────────────────────────────────────────────────────────────────────── #
# ChebyshevData — the matrix-free cell's workspace CONTAINER (slice 3b: the
# struct moves so the shared dispatch gates can check it; construction and
# consumers stay in Piccolissimo).
# ─────────────────────────────────────────────────────────────────────────── #

# ============================================================================ #

"""
    ChebyshevData

Per-interval workspaces for [`ChebyshevAlg`](@ref) forward propagation on
`KetTrajectory` (#132). One `MagnusProductCore` + coefficient evaluator +
checkpoint buffer per knot interval (intervals may have different `n_sub`;
per-interval workspaces also make the threaded `evaluate!` safe). The
propagator `Φ` is never materialized — only `ketdim`-sized states are stored.

# Fields
- `cores`: per-interval `MagnusProductCore` (fixed `n_sub`, shared bracket)
- `coeff_fns`: per-interval callable coefficient evaluators (packed-parameter →
  drive coefficients via `SplineIntervalCoeffs`)
- `checkpoints`: per-interval substep states, length `n_sub[k] + 1`
- `n_sub`: per-interval substep counts, **frozen at construction**
- `suggested_n_sub`: the `:auto` suggestion (equal to `n_sub` under `:auto`)
- `bracket`: the spectral bracket `(a, b)` all exp-actions claim
- `phase_budget`: per-panel phase budget `φ_max` used by sizing and the
  per-iterate re-check
- `spec_drift`: frozen Gershgorin interval `(glo, ghi)` of `H_drift` (precomputed
  once), the base of the current-controls bracket re-check (#223)
- `spec_radius_drives`: frozen per-drive Gershgorin spectral-radius bounds of the
  `H_l` (precomputed once); the re-check widens `spec_drift` by the CURRENT
  `Σ|c_l|·ρ(H_l)` and checks containment in `bracket` (#223)
- `u_dim`: controls + globals dimension (packed-parameter layout)
- `phase_warned`: warn-once latch for the free-time-dilation re-check
- `vjp_λ`: per-interval costate buffer for the VJP adjoint sweep (#133)
- `vjp_scatters`: per-interval packed-gradient scatter callables (#133); each
  wraps the SAME per-interval `SplineIntervalCoeffs` as `coeff_fns`
- `vjp_lin_du`: preallocated per-interval linear-spline `du` fold-back buffer
  (length `2·u_dim+2`), reused (zeroed in place) by
  `_fold_packed_grad(::LinearSpline)` so the linear VJP hot path
  allocates nothing (per-interval ⇒ race-free) (#223)
"""
struct ChebyshevData{C,F,V,Sc,CD,HS,MS}
    cores::Vector{C}
    coeff_fns::Vector{F}
    checkpoints::Vector{Vector{V}}
    n_sub::Vector{Int}
    suggested_n_sub::Vector{Int}
    bracket::Tuple{Float64,Float64}
    phase_budget::Float64
    # Frozen Gershgorin data for the current-controls bracket re-check (#223): the
    # H_drift Gershgorin interval + per-drive |H_l| spectral-radius bounds,
    # precomputed once. The forward re-check widens the drift interval by the
    # CURRENT Σ|c_l|·ρ(H_l) (O(n_drives), allocation-free) and warns if the
    # resulting spectrum interval escapes the frozen bracket [a, b]. (Interval
    # containment, not an origin-radius vs half-width test — the latter would
    # false-positive on asymmetric-spectrum drifts, e.g. a number operator.)
    spec_drift::Tuple{Float64,Float64}
    spec_radius_drives::Vector{Float64}
    u_dim::Int
    phase_warned::Threads.Atomic{Bool}
    vjp_λ::Vector{V}
    vjp_scatters::Vector{Sc}
    # Preallocated per-interval linear-spline `du` fold-back buffer (#223, length
    # 2·u_dim+2): reused (zeroed in place) by `_fold_packed_grad(::LinearSpline)` so
    # the linear VJP hot path allocates nothing (per-interval ⇒ race-free).
    vjp_lin_du::Vector{Vector{Float64}}
    # Per-interval HVP scratch (#134): tangent-node checkpoints, the forward-over-
    # reverse working buffers, the packed direction/gradient, and the three
    # per-interval callables (coeff-direction, δP scatter, coeff-Hessian scatter),
    # all sharing the matching `sics[k]`. Distinct per interval ⇒ chunk-parallel
    # per-knot accumulation is race-free (same discipline as the VJP scratch).
    hvp_dcheckpoints::Vector{Vector{V}}
    # (slice 3b: `MagnusHvpScratch{V}` — the MF HVP scratch type — stays
    # proprietary in Piccolissimo; the field is typed by a dedicated parameter
    # `MS`, whose concrete binding is exactly `MagnusHvpScratch{V}` at every
    # construction site, so the layout is unchanged.)
    hvp_scratch::Vector{MS}
    hvp_g::Vector{Vector{Float64}}
    hvp_δp::Vector{Vector{Float64}}
    hvp_coeff_dirs::Vector{CD}
    hvp_vjp_scatters::Vector{Sc}
    hvp_hvp_scatters::Vector{HS}
end
