# ============================================================================ #
# Integration Algorithm Types
#
# Abstract hierarchy for forward propagator solvers.
# The algorithm determines HOW the propagator ODE dΦ/dτ = A(τ)Φ is solved,
# and is also used for the corresponding sensitivity ODE.
# ============================================================================ #

export IntegrationAlgorithm, Tsit5Alg, MagnusGL4Alg, MagnusAdapt4Alg, Rodas5PAlg
export ChebyshevAlg

"""
    IntegrationAlgorithm

Abstract type for forward propagation algorithms used by `SplineIntegrator`.

Concrete subtypes:
- [`Tsit5Alg`](@ref): Adaptive Runge-Kutta (Tsit5) — the default workhorse
- [`MagnusGL4Alg`](@ref): 4th-order Magnus via Gauss-Legendre — Lie group structure preserving
- [`Rodas5PAlg`](@ref): 5th-order L-stable Rosenbrock — for stiff regimes
"""
abstract type IntegrationAlgorithm end

"""
    Tsit5Alg(; adaptive=true, tol=1e-6, ode_h=0.1)

Adaptive Tsit5 Runge-Kutta for the forward propagator ODE.

Robust and efficient for smooth Hamiltonians at moderate dimensions.
Supports both adaptive and fixed-step modes: with `adaptive=false`, the forward
propagator solve for `UnitaryTrajectory`, `KetTrajectory`, and
`MultiKetTrajectory` steps with fixed size `ode_h` in normalized interval time
τ ∈ [0, 1]. Sensitivity ODEs are always solved adaptively regardless of this
setting.

# Keyword Arguments
- `adaptive::Bool=true`: Use adaptive step size control
- `tol::Float64=1e-6`: Absolute and relative tolerance for adaptive mode
- `ode_h::Float64=0.1`: Fixed step size when `adaptive=false`, in normalized
  interval time (e.g. `ode_h=0.1` takes 10 steps per knot interval)
"""
struct Tsit5Alg <: IntegrationAlgorithm
    adaptive::Bool
    tol::Float64
    ode_h::Float64
end
Tsit5Alg(; adaptive = true, tol = 1e-6, ode_h = 0.1) = Tsit5Alg(adaptive, tol, ode_h)

"""
    _solve_forward_tsit5(prob, alg, tol)

Solve a forward propagator ODE problem with Tsit5, honoring the stepping mode
of `alg::Tsit5Alg`: adaptive stepping at tolerance `tol`, or fixed steps of
size `alg.ode_h` when `alg.adaptive == false` (#180). The fallback method
keeps non-`Tsit5Alg` algorithms that reach a Tsit5 forward solve (e.g.
`Rodas5PAlg` on the ket-sensitivity path) on the adaptive solve they used
before. Sensitivity ODEs are always solved adaptively and must not route
through the fixed-step branch.
"""
_solve_forward_tsit5(prob, alg::Tsit5Alg, tol::Float64) =
    alg.adaptive ? solve(prob, Tsit5(); abstol = tol, reltol = tol, saveat = 1.0) :
    solve(prob, Tsit5(); dt = alg.ode_h, adaptive = false, saveat = 1.0)
_solve_forward_tsit5(prob, ::IntegrationAlgorithm, tol::Float64) =
    solve(prob, Tsit5(); abstol = tol, reltol = tol, saveat = 1.0)

"""
    MagnusGL4Alg(; n_steps=10, tol=1e-6)

4th-order Magnus integrator via Gauss-Legendre quadrature.

Preserves Lie group structure (Φ stays unitary to machine precision).
Requires explicit matrix Hamiltonians (H_drift, H_drives).

!!! note "KetTrajectory uses a different cell"
    On a `UnitaryTrajectory` this is the 4th-order fixed-step Gauss–Legendre
    Magnus path described here (materializing `Φ` via OrdinaryDiffEq). On a
    `KetTrajectory`, `SplineIntegrator` instead routes this alg through the
    **matrix-free 2nd-order-per-substep Magnus + Chebyshev exp-action cell**
    (#223): each knot interval is sub-stepped by the midpoint-Magnus product and
    `exp(Ω)` is applied through the Chebyshev expansion — `Φ` is never
    materialized. The substep count is phase-budget sized with `n_steps` as a
    floor (raise `n_steps` for accuracy, as below); the analytic Duhamel VJP is
    the exact adjoint of that discrete forward. This is a distinct discretization
    from the Unitary 4th-order path.

# Keyword Arguments
- `n_steps::Int=10`: Number of fixed integration steps per knot interval. This is
  the **sole accuracy knob** for this algorithm: the substep is `Δtₖ/n_steps`, so
  accuracy improves (and cost grows) with `n_steps`. Raise it for stiff/large-‖H‖
  regimes — e.g. `n_steps≈50` for a deep Rydberg blockade (V/Ω ≫ 1) — where the
  default 10 under-resolves and the optimizer's fidelity diverges from a fine
  re-rollout.
- `tol::Float64=1e-6`: **Effectively inert** for `MagnusGL4Alg` — the integration
  is fixed-step (governed by `n_steps`), so this tolerance does not tighten it.
  Retained for interface uniformity with the adaptive algorithms.
"""
struct MagnusGL4Alg <: IntegrationAlgorithm
    n_steps::Int
    tol::Float64
end
MagnusGL4Alg(; n_steps = 10, tol = 1e-6) = MagnusGL4Alg(n_steps, tol)

"""
    MagnusAdapt4Alg(; tol=1e-6)

4th-order adaptive Magnus integrator with automatic step size control.

Preserves Lie group structure (Φ stays unitary to machine precision) like `MagnusGL4Alg`,
but uses embedded error estimation to adapt step sizes for guaranteed accuracy.

!!! note "KetTrajectory uses a different cell — `tol` is inert there"
    The adaptive, embedded-error-controlled path described here applies on a
    `UnitaryTrajectory`. On a `KetTrajectory`, `SplineIntegrator` routes this alg
    through the **matrix-free 2nd-order-per-substep Magnus + Chebyshev exp-action
    cell** (#223) shared with `MagnusGL4Alg` — a FIXED-step, phase-budget-sized
    discretization with no embedded error control. There `tol` is **inert** (a
    construction-time `@warn` fires): control accuracy with
    `MagnusGL4Alg(n_steps=…)` instead of this adaptive alg.

# Keyword Arguments
- `tol::Float64=1e-6`: Absolute and relative tolerance for adaptive stepping
  (`UnitaryTrajectory` only; **inert on `KetTrajectory`** — see the note above)
"""
struct MagnusAdapt4Alg <: IntegrationAlgorithm
    tol::Float64
end
MagnusAdapt4Alg(; tol = 1e-6) = MagnusAdapt4Alg(tol)

"""
    ChebyshevAlg(; n_sub=:auto, sub_dt=nothing, bracket=nothing, kwargs...)

Matrix-free Chebyshev exp-action forward propagation for `KetTrajectory`
(#132, ADR-0003). Each knot interval is discretized into `n_sub` midpoint-Magnus
substeps; every substep applies `exp(−iθH(sⱼ))` to the state through the fixed
Chebyshev polynomial expansion (`chebyshev_expv!`) — the propagator `Φ` is never
materialized. PWC/constant-`H` intervals are the `n_sub = 1` corner of the same
scaffold.

# Sub-step sizing
- `n_sub = :auto` (default): each interval is sized at construction from the
  per-panel phase budget `θ·r ≤ phase_budget` (`θ = Δtₖ/n_sub`, `r` the spectral
  bracket radius), then verified by Richardson rollout and probe-direction
  gradient checks (`n_sub` vs `2·n_sub`, 2nd-order contraction) against
  `dyn_tol` / `grad_tol`, doubling until both pass (see [`suggest_n_sub`](@ref)).
- `n_sub::Int`: uniform manual sub-step count, honored with a loud warning when
  it falls below the suggestion. Takes precedence over `sub_dt` if both are given.
- `sub_dt::Real`: manual target substep duration; per-interval
  `n_subₖ = max(1, ceil(Δtₖ/sub_dt))`, honored with the same warning.

Sizing is **frozen after construction**: a cheap per-iterate re-check warns once
if free-time dilation pushes the per-panel phase past `phase_budget`, but never
re-sizes (no discretization jumps along an outer optimizer's search direction).

# Spectral bracket
- `bracket = nothing` (default): a control-independent bracket is derived from
  Gershgorin bounds over the control box — `H_drift`'s Gershgorin interval
  widened by `Σ_l max|c_l|·ρ(H_l)` with per-drive coefficient extremes taken over
  the knot-box hull (`LinearSpline`) or the Hermite-overshoot-inflated hull
  (`CubicSpline`) — then safety-inflated ×1.1 to tolerate transient
  augmented-Lagrangian box violations. Construction **errors** when a
  dynamics-entering control is unbounded or a drive coefficient cannot be
  bounded over the box (e.g. `NonlinearDrive`) — supply `bracket` explicitly.
- `bracket = (a, b)`: used verbatim (no inflation).

!!! warning "A too-narrow bracket degrades silently"
    The Bessel-tail guard inside `chebyshev_expv!` sees only the *claimed*
    phase `τ·r`: it throws loudly when the bracket/phase is too **wide** for the
    workspace degree, but a bracket that is too **narrow** (spectrum outside
    `[a, b]`) shows up only as slow/failed coefficient decay in lucky cases —
    otherwise the result quietly loses accuracy. If you override `bracket`,
    make sure it truly contains the spectrum over every admissible control.

# Keyword Arguments
- `n_sub::Union{Symbol,Int} = :auto`: sub-step count (`:auto` or uniform manual)
- `sub_dt::Union{Nothing,Real} = nothing`: manual target substep duration
- `bracket::Union{Nothing,Tuple} = nothing`: explicit spectral bracket `(a, b)`
- `n_quad::Int = 12`: Gauss–Legendre nodes per substep (Duhamel probe/sensitivity)
- `cheb_tol::Float64 = 1e-13`: Chebyshev coefficient-tail tolerance
- `maxdeg::Int = 512`: Chebyshev workspace degree cap (Bessel-tail guard trips
  loudly past it)
- `dyn_tol::Float64 = 1e-8`: `:auto` Richardson rollout tolerance
- `grad_tol::Float64 = 1e-6`: `:auto` probe-direction gradient tolerance
- `phase_budget::Float64 = 2.0`: per-panel phase budget `φ_max`
- `n_sub_max::Int = 4096`: `:auto` sizing cap (construction errors past it)
- `tol::Float64 = 1e-6`: interface-uniform tolerance knob (the forward action is
  governed by `cheb_tol`/`n_sub`; retained for constructor-seam uniformity)

Sensitivities (JVP/VJP/HVP) are **not yet available** for `ChebyshevAlg`: per
ADR-0003 Decision 2, requesting one errors instead of silently falling back to
the Tsit5 augmented sensitivity ODE. The analytic Chebyshev/Magnus Duhamel
sensitivity lands in #133.
"""
struct ChebyshevAlg <: IntegrationAlgorithm
    n_sub::Union{Symbol,Int}
    sub_dt::Union{Nothing,Float64}
    bracket::Union{Nothing,Tuple{Float64,Float64}}
    n_quad::Int
    cheb_tol::Float64
    maxdeg::Int
    dyn_tol::Float64
    grad_tol::Float64
    phase_budget::Float64
    n_sub_max::Int
    tol::Float64
end

function ChebyshevAlg(;
    n_sub::Union{Symbol,Int} = :auto,
    sub_dt::Union{Nothing,Real} = nothing,
    bracket::Union{Nothing,Tuple{<:Real,<:Real}} = nothing,
    n_quad::Int = 12,
    cheb_tol::Real = 1e-13,
    maxdeg::Int = 512,
    dyn_tol::Real = 1e-8,
    grad_tol::Real = 1e-6,
    phase_budget::Real = 2.0,
    n_sub_max::Int = 4096,
    tol::Real = 1e-6,
)
    if n_sub isa Symbol && n_sub !== :auto
        throw(
            ArgumentError(
                "ChebyshevAlg: n_sub must be :auto or a positive Int (got :$n_sub)",
            ),
        )
    end
    if n_sub isa Int && n_sub < 1
        throw(ArgumentError("ChebyshevAlg: n_sub must be ≥ 1 (got $n_sub)"))
    end
    if !isnothing(sub_dt) && sub_dt <= 0
        throw(ArgumentError("ChebyshevAlg: sub_dt must be positive (got $sub_dt)"))
    end
    if !isnothing(bracket) && !(bracket[2] > bracket[1])
        throw(ArgumentError("ChebyshevAlg: bracket must satisfy b > a (got $bracket)"))
    end
    phase_budget > 0 || throw(ArgumentError("ChebyshevAlg: phase_budget must be positive"))
    return ChebyshevAlg(
        n_sub,
        isnothing(sub_dt) ? nothing : float(sub_dt),
        isnothing(bracket) ? nothing : (float(bracket[1]), float(bracket[2])),
        n_quad,
        float(cheb_tol),
        maxdeg,
        float(dyn_tol),
        float(grad_tol),
        float(phase_budget),
        n_sub_max,
        float(tol),
    )
end

"""
    Rodas5PAlg(; tol=1e-6)

Stiff-aware 5th-order L-stable Rosenbrock-Wanner integrator (Rodas5P from OrdinaryDiffEq).

Use when standard non-stiff solvers (`Tsit5Alg`) hit max-iterations on stiff Hamiltonians,
or when an outer optimizer's line-search trial points probe stiff regimes that Tsit5
cannot handle without `maxiters` overrides.

Used by both forward propagation and sensitivity equations.

# Keyword Arguments
- `tol::Float64=1e-6`: Absolute and relative tolerance
"""
struct Rodas5PAlg <: IntegrationAlgorithm
    tol::Float64
end
Rodas5PAlg(; tol = 1e-6) = Rodas5PAlg(tol)
