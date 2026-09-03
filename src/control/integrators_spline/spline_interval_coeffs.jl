# =============================================================================
# Cubic-Hermite interval coefficient layer — ALTISSIMO 1.6 slice 3a (#198)
#
# Maps the packed per-interval parameter vector
#
#     p = [uₖ (u_dim); uₖ₊₁ (u_dim); duₖ (u_dim); duₖ₊₁ (u_dim); Δtₖ; tₖ]
#
# (the `build_ode_params` layout — u-blocks INCLUDE appended globals, whose
# du-entries are zero) to the drive-coefficient closures the Magnus product core
# consumes, replicating the forward RHS EXACTLY (PiccolissimoCUDAExt `rhs!`):
#
#     h00 = 2τ³−3τ²+1        h01 = −2τ³+3τ²
#     h10 = (τ³−2τ²+τ)·Δt    h11 = (τ³−τ²)·Δt          ← derivative basis Δt-SCALED
#     u(τ) = h00·uₖ + h01·uₖ₊₁ + h10·duₖ + h11·duₖ₊₁
#     c_l(τ) = drive_coeff(drives[l], u(τ))
#
# Globals ride the u-blocks as constants automatically (h00+h01 ≡ 1, du = 0).
#
# Parameter-class derivatives (all exact chain rules through the basis):
#   ∂c_l/∂uₖ[i]   = J_l[i]·h00      ∂c_l/∂uₖ₊₁[i] = J_l[i]·h01
#   ∂c_l/∂duₖ[i]  = J_l[i]·h10      ∂c_l/∂duₖ₊₁[i]= J_l[i]·h11
#   ∂c_l/∂Δt      = J_l · ∂u/∂Δt,   ∂u/∂Δt = (τ³−2τ²+τ)·duₖ + (τ³−τ²)·duₖ₊₁
# where J_l[i] = drive_coeff_jac(drives[l], u(τ), i). Δt ALSO enters the
# generator scale θ = Δt/n_sub — that exact (quadrature-free) part lives in the
# Magnus core; this layer supplies only the coefficient part. `tₖ` derivatives
# are zero for the in-scope systems (drive coefficients ignore t; time-modulated
# coefficients are out of #198's scope).
# =============================================================================

using Piccolo: drive_coeff, drive_coeff_jac, drive_coeff_hess
using TestItemRunner: @testitem

export SplineIntervalCoeffs, interval_coeff!, interval_coeff_dir!, interval_vjp_scatter!
export interval_hvp_scatter!

"""
    SplineIntervalCoeffs(drives, u_dim, p)

Per-interval coefficient evaluator over the packed cubic-Hermite parameter
vector `p` (length `4·u_dim + 2`). Reusable buffers; update `p` in place (or
`copyto!(sic.p, ...)`) between intervals.
"""
struct SplineIntervalCoeffs{DR}
    drives::DR
    u_dim::Int
    p::Vector{Float64}
    u::Vector{Float64}          # u(τ) buffer
    dudΔt::Vector{Float64}      # ∂u/∂Δt buffer
    δu::Vector{Float64}         # directional δu(τ) buffer
end

function SplineIntervalCoeffs(drives, u_dim::Int, p::Vector{Float64})
    @assert length(p) == 4u_dim + 2 "packed interval params must have length 4·u_dim+2"
    return SplineIntervalCoeffs(drives, u_dim, p, zeros(u_dim), zeros(u_dim), zeros(u_dim))
end

@inline function _hermite(τ::Float64, Δt::Float64)
    τ2 = τ * τ
    τ3 = τ2 * τ
    h00 = 2τ3 - 3τ2 + 1
    h01 = -2τ3 + 3τ2
    g10 = τ3 - 2τ2 + τ            # UNSCALED derivative bases
    g11 = τ3 - τ2
    return h00, h01, g10 * Δt, g11 * Δt, g10, g11
end

@inline _Δt(sic::SplineIntervalCoeffs) = sic.p[4*sic.u_dim+1]

function _u_at!(sic::SplineIntervalCoeffs, τ::Float64)
    u_dim = sic.u_dim
    p = sic.p
    h00, h01, h10, h11, _, _ = _hermite(τ, _Δt(sic))
    @inbounds for i = 1:u_dim
        sic.u[i] = h00 * p[i] + h01 * p[u_dim+i] + h10 * p[2u_dim+i] + h11 * p[3u_dim+i]
    end
    return sic.u
end

"""
    interval_coeff!(c, sic, τ)

Fill `c[l] = drive_coeff(drives[l], u(τ))` — the Magnus core's `coeff!`.
"""
function interval_coeff!(c, sic::SplineIntervalCoeffs, τ::Real)
    u = _u_at!(sic, float(τ))
    @inbounds for (l, d) in enumerate(sic.drives)
        c[l] = drive_coeff(d, u)
    end
    return c
end

"""
    interval_coeff_dir!(dc, sic, τ, δp)

Fill the DIRECTIONAL coefficient derivative `dc[l] = (∂c_l/∂p)·δp` at `τ`, for a
packed direction `δp` (length `4·u_dim+2`; the `tₖ` slot is ignored). Includes
the Δt-through-the-basis coefficient term (`δp`'s Δt slot); the exact θ-scale
part of Δt is handled by the Magnus core, not here.
"""
function interval_coeff_dir!(dc, sic::SplineIntervalCoeffs, τ::Real, δp::AbstractVector)
    u_dim = sic.u_dim
    p = sic.p
    τf = float(τ)
    Δt = _Δt(sic)
    h00, h01, h10, h11, g10, g11 = _hermite(τf, Δt)
    u = _u_at!(sic, τf)
    δΔt = δp[4u_dim+1]
    @inbounds for i = 1:u_dim
        # δu from the four control blocks + Δt-scaling of the derivative bases
        sic.δu[i] =
            h00 * δp[i] +
            h01 * δp[u_dim+i] +
            h10 * δp[2u_dim+i] +
            h11 * δp[3u_dim+i] +
            δΔt * (g10 * p[2u_dim+i] + g11 * p[3u_dim+i])
    end
    @inbounds for (l, d) in enumerate(sic.drives)
        acc = 0.0
        for i = 1:u_dim
            δui = sic.δu[i]
            if δui != 0.0
                acc += drive_coeff_jac(d, u, i) * δui
            end
        end
        dc[l] = acc
    end
    return dc
end

"""
    interval_vjp_scatter!(g, sic, l, τ, val)

Scatter one adjoint pairing `val = −i·θ·w_q·⟨μ_q, H_l·φ_q⟩` (from
`magnus_vjp!`'s accumulation callback) into the packed-parameter gradient `g`
(length `4·u_dim+2`, real): all four control blocks via the Hermite bases, plus
the Δt slot's coefficient term. The caller adds the Magnus core's exact `g_Δt`
(θ-scale term) to `g[4·u_dim+1]` separately. Real-part convention: the iso-real
pairing used by the equality-constraint Jacobian upstream.
"""
function interval_vjp_scatter!(
    g::AbstractVector,
    sic::SplineIntervalCoeffs,
    l::Int,
    τ::Real,
    val::Complex,
)
    u_dim = sic.u_dim
    p = sic.p
    τf = float(τ)
    Δt = _Δt(sic)
    h00, h01, h10, h11, g10, g11 = _hermite(τf, Δt)
    u = _u_at!(sic, τf)
    d = sic.drives[l]
    rv = real(val)
    @inbounds for i = 1:u_dim
        Ji = drive_coeff_jac(d, u, i)
        if Ji != 0.0
            g[i] += rv * Ji * h00
            g[u_dim+i] += rv * Ji * h01
            g[2u_dim+i] += rv * Ji * h10
            g[3u_dim+i] += rv * Ji * h11
            g[4u_dim+1] += rv * Ji * (g10 * p[2u_dim+i] + g11 * p[3u_dim+i])
        end
    end
    return g
end

"""
    interval_hvp_scatter!(g, sic, l, τ, val, δp)

Second-order (HVP) scatter: accumulate `Re(val)·δ(∂c_l/∂p)` into the packed
gradient `g` (length `4·u_dim+2`, real), where `δ(∂c_l/∂p)` is the directional
derivative (in the packed direction `δp`) of the chain-rule map `∂c_l/∂p` that
[`interval_vjp_scatter!`](@ref) contracts. Two sources (#134):

- **coefficient Hessian** (`drive_coeff_hess`): `δ(J_l[i]) = Σ_m ∂²c_l/∂u_i∂u_m ·
  δu(τ)_m` — the `∂²H/∂c²` same-time term Gauss–Newton drops (nonzero for
  `NonlinearDrive`); scattered through the same four Hermite bases + the Δt-slot
  coefficient term as the first-order scatter;
- **basis second derivatives**: the cubic `du` bases carry an explicit `Δt`
  factor (`h10 = g10·Δt`, `h11 = g11·Δt`) and the Δt-slot's coefficient term
  (`g10·duₖ + g11·duₖ₊₁`) depends on the `du` params — so `δ(∂u_i/∂p)` couples the
  `du` blocks with the Δt slot. `val = ⟨costate, ∂_{c_l}E·state⟩` is the same
  first-order pairing `magnus_vjp!` produces; the caller passes it verbatim.
"""
function interval_hvp_scatter!(
    g::AbstractVector,
    sic::SplineIntervalCoeffs,
    l::Int,
    τ::Real,
    val::Complex,
    δp::AbstractVector,
)
    u_dim = sic.u_dim
    p = sic.p
    τf = float(τ)
    Δt = _Δt(sic)
    h00, h01, h10, h11, g10, g11 = _hermite(τf, Δt)
    u = _u_at!(sic, τf)
    d = sic.drives[l]
    δΔt = δp[4u_dim+1]
    rv = real(val)
    # directional δu(τ) (same as interval_coeff_dir!'s δu; includes Δt-through-basis)
    @inbounds for i = 1:u_dim
        sic.δu[i] =
            h00 * δp[i] +
            h01 * δp[u_dim+i] +
            h10 * δp[2u_dim+i] +
            h11 * δp[3u_dim+i] +
            δΔt * (g10 * p[2u_dim+i] + g11 * p[3u_dim+i])
    end
    # (a) coefficient-Hessian part: δ(J_l[i]) = Σ_m coeff_hess(d, u, i, m)·δu_m
    @inbounds for i = 1:u_dim
        δJi = 0.0
        for m = 1:u_dim
            δum = sic.δu[m]
            if δum != 0.0
                δJi += drive_coeff_hess(d, u, i, m) * δum
            end
        end
        if δJi != 0.0
            g[i] += rv * δJi * h00
            g[u_dim+i] += rv * δJi * h01
            g[2u_dim+i] += rv * δJi * h10
            g[3u_dim+i] += rv * δJi * h11
            g[4u_dim+1] += rv * δJi * (g10 * p[2u_dim+i] + g11 * p[3u_dim+i])
        end
    end
    # (b) basis-derivative part: δ(∂u_i/∂p) — du bases carry Δt (→ δΔt), and the
    # Δt-slot coefficient term depends on du params (→ δduₖ, δduₖ₊₁).
    @inbounds for i = 1:u_dim
        Ji = drive_coeff_jac(d, u, i)
        if Ji != 0.0
            g[2u_dim+i] += rv * Ji * g10 * δΔt
            g[3u_dim+i] += rv * Ji * g11 * δΔt
            g[4u_dim+1] += rv * Ji * (g10 * δp[2u_dim+i] + g11 * δp[3u_dim+i])
        end
    end
    return g
end
