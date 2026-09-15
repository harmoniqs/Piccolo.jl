# ============================================================================= #
#   Cubic Hermite spline primitives — the ONE definition                        #
# ============================================================================= #
#
# These were duplicated across two constraint files: `spline_bound_constraints.jl`
# defined the basis functions and the critical-point solvers, and
# `cubic_spline_bound_constraint.jl` re-derived the same four basis polynomials
# inline inside `eval_cubic_hermite` / `eval_cubic_hermite_jacobian`. Both now call
# these — and as of #331 those two adapters are gone entirely, the value gradient having
# moved here as `hermite_value_gradient`. Everything here is generic in the scalar type
# (ForwardDiff.Dual included) and allocation-free except the two root finders, which
# return a small vector.
#
# Part of #330 (ADR-0010): a representation change needs one source of truth for
# the stencil math, and two independently-maintained copies of the Hermite basis
# is where a wrong coefficient hides.

# ----------------------------------------------------------------------------- #
# Basis
# ----------------------------------------------------------------------------- #

"""
    hermite_basis_functions(τ)

Evaluate cubic Hermite basis functions at normalized time τ ∈ [0,1].

Returns (h00, h10, h01, h11) where:
- h00(τ): basis for u_k (position at start)
- h10(τ): basis for du_k * Δt (tangent at start)
- h01(τ): basis for u_{k+1} (position at end)
- h11(τ): basis for du_{k+1} * Δt (tangent at end)
"""
function hermite_basis_functions(τ::T) where {T}
    τ² = τ * τ
    τ³ = τ² * τ

    h00 = 2τ³ - 3τ² + 1
    h10 = τ³ - 2τ² + τ
    h01 = -2τ³ + 3τ²
    h11 = τ³ - τ²

    return (h00, h10, h01, h11)
end

"""
    hermite_derivative_basis(τ)

Evaluate derivative of cubic Hermite basis functions: dh/dτ at τ ∈ [0,1].
"""
function hermite_derivative_basis(τ::T) where {T}
    τ² = τ * τ

    dh00 = 6τ² - 6τ
    dh10 = 3τ² - 4τ + 1
    dh01 = -6τ² + 6τ
    dh11 = 3τ² - 2τ

    return (dh00, dh10, dh01, dh11)
end

"""
    evaluate_hermite_spline(τ, u_k, u_kp1, du_k, du_kp1, Δt)

Evaluate cubic Hermite spline at normalized time τ ∈ [0,1].

# Arguments
- `τ`: Normalized time in [0, 1]
- `u_k`: Control value at knot k
- `u_kp1`: Control value at knot k+1
- `du_k`: Control derivative at knot k
- `du_kp1`: Control derivative at knot k+1
- `Δt`: Time step size
"""
function evaluate_hermite_spline(τ, u_k, u_kp1, du_k, du_kp1, Δt)
    h00, h10, h01, h11 = hermite_basis_functions(τ)
    return h00 * u_k + h10 * Δt * du_k + h01 * u_kp1 + h11 * Δt * du_kp1
end

"""
    hermite_value_gradient(du_k, du_kp1, Δt, τ)

Gradient of [`evaluate_hermite_spline`](@ref) at a FIXED τ with respect to the segment's
decision variables, returned as the tuple

    (∂s/∂u_k, ∂s/∂u_{k+1}, ∂s/∂du_k, ∂s/∂du_{k+1}, ∂s/∂Δt)
    = (h00, h01, Δt·h10, Δt·h11, h10·du_k + h11·du_{k+1})

The tuple order is the column order the interior-point bound constraint declares to its
stencil table, so a coefficient refresh can splat it straight into `coeffs`. Returned as
a tuple, not a vector: `cubic_spline_bound_constraint.jl` previously allocated a
five-element heap vector per sample point per Jacobian (#331 / ADR-0010).

`s` is linear in `u`/`du`, so this gradient depends on `u` not at all and on `du` only
through the `Δt` partial.

Generic in the scalar type (ForwardDiff.Dual included); allocation-free.
"""
@inline function hermite_value_gradient(du_k, du_kp1, Δt, τ)
    h00, h10, h01, h11 = hermite_basis_functions(τ)
    return (h00, h01, Δt * h10, Δt * h11, h10 * du_k + h11 * du_kp1)
end

"""
    evaluate_hermite_derivative(τ, u_k, u_kp1, du_k, du_kp1, Δt)

Evaluate derivative of cubic Hermite spline at normalized time τ ∈ [0,1].

Returns du/dt at parameter τ.
"""
function evaluate_hermite_derivative(τ, u_k, u_kp1, du_k, du_kp1, Δt)
    dh00, dh10, dh01, dh11 = hermite_derivative_basis(τ)
    # du/dt = (du/dτ) / Δt
    du_dτ = dh00 * u_k + dh10 * Δt * du_k + dh01 * u_kp1 + dh11 * Δt * du_kp1
    return du_dτ / Δt
end

# ----------------------------------------------------------------------------- #
# Second derivative (acceleration) at the two segment endpoints
# ----------------------------------------------------------------------------- #

"""
    hermite_accel_start(u_k, u_kp1, du_k, du_kp1, Δt)

Acceleration `d²u/dt²` of the cubic Hermite segment `k → k+1` at its START (τ=0):

    a_start = (6/Δt²)(u_{k+1} - u_k) - (2/Δt)(2 du_k + du_{k+1})

Generic in the scalar type; allocation-free.
"""
@inline function hermite_accel_start(u_k, u_kp1, du_k, du_kp1, Δt)
    # The association `6 * (inv²) * Δu` is deliberate, not incidental: it reproduces
    # the pre-port expression BIT-FOR-BIT, which is what lets the residual-parity
    # witness assert `==` rather than a tolerance. `(6 * inv) * inv * Δu` differs in
    # the last digit.
    inv_Δt = inv(Δt)
    inv_Δt² = inv_Δt * inv_Δt
    return 6 * inv_Δt² * (u_kp1 - u_k) - 2 * inv_Δt * (2 * du_k + du_kp1)
end

"""
    hermite_accel_end(u_k, u_kp1, du_k, du_kp1, Δt)

Acceleration `d²u/dt²` of the cubic Hermite segment `k → k+1` at its END (τ=1):

    a_end = -(6/Δt²)(u_{k+1} - u_k) + (2/Δt)(du_k + 2 du_{k+1})

Generic in the scalar type; allocation-free.
"""
@inline function hermite_accel_end(u_k, u_kp1, du_k, du_kp1, Δt)
    # See `hermite_accel_start` on the association.
    inv_Δt = inv(Δt)
    inv_Δt² = inv_Δt * inv_Δt
    return -6 * inv_Δt² * (u_kp1 - u_k) + 2 * inv_Δt * (du_k + 2 * du_kp1)
end

# ----------------------------------------------------------------------------- #
# Acceleration gradients as tuples (#333)
# ----------------------------------------------------------------------------- #
#
# WHY THESE EXIST AND WHY THEY ARE NOT A REFACTOR. `hermite_value_gradient` above already
# returns the cubic-bound constraint's stencil row as a tuple, so #333's device coefficient
# refresh can BROADCAST the very function the host refresh calls — one source of truth, no
# divergence possible. The acceleration constraint had no such function: its eighteen
# partials live inline inside `_hsa_refresh_coefficients!`. A device kernel that re-derived
# them would be the second copy ADR-0010 warns about ("divergence here is the main
# long-term risk").
#
# These three functions are that missing tuple form, ADDED rather than refactored in:
# `_hsa_refresh_coefficients!` is the host reference and stays byte-for-byte as #330 left
# it (#333 AC6). The expressions here reproduce it EXACTLY — same inverse powers, same
# association, same operand order — and a testitem asserts that bit-for-bit over randomized
# inputs for all three kinds, so the duplication is pinned rather than trusted.
#
# Generic in the scalar type, allocation-free, and GPU-safe (tuple return, no branches).

"""
    hermite_accel_start_gradient(u_k, u_kp1, du_k, du_kp1, Δt) -> NTuple{5}

Gradient of [`hermite_accel_start`](@ref) with respect to
`(u_k, u_{k+1}, du_k, du_{k+1}, Δt)` — the column order the acceleration constraint
declares to its stencil table for an `HSA_A_START` functional.
"""
@inline function hermite_accel_start_gradient(u_k, u_kp1, du_k, du_kp1, Δt)
    i = 1.0 / Δt
    i2 = i * i
    i3 = i2 * i
    return (
        -6 * i2,
        6 * i2,
        -4 * i,
        -2 * i,
        -12 * i3 * (u_kp1 - u_k) + 2 * i2 * (2 * du_k + du_kp1),
    )
end

"""
    hermite_accel_end_gradient(u_k, u_kp1, du_k, du_kp1, Δt) -> NTuple{5}

Gradient of [`hermite_accel_end`](@ref) with respect to
`(u_k, u_{k+1}, du_k, du_{k+1}, Δt)` — the `HSA_A_END` column order.
"""
@inline function hermite_accel_end_gradient(u_k, u_kp1, du_k, du_kp1, Δt)
    i = 1.0 / Δt
    i2 = i * i
    i3 = i2 * i
    return (
        6 * i2,
        -6 * i2,
        2 * i,
        4 * i,
        12 * i3 * (u_kp1 - u_k) - 2 * i2 * (du_k + 2 * du_kp1),
    )
end

"""
    hermite_accel_jump_gradient(u_km1, u_k, u_kp1, du_km1, du_k, du_kp1, Δt_km1, Δt_k)

Gradient of the acceleration JUMP
`hermite_accel_end(u_{k-1}, u_k, du_{k-1}, du_k, Δt_{k-1}) −
hermite_accel_start(u_k, u_{k+1}, du_k, du_{k+1}, Δt_k)`
with respect to `(u_{k-1}, u_k, u_{k+1}, du_{k-1}, du_k, du_{k+1}, Δt_{k-1}, Δt_k)` — the
`HSA_JUMP` column order. Returned as an `NTuple{8}`.
"""
@inline function hermite_accel_jump_gradient(
    u_km1,
    u_k,
    u_kp1,
    du_km1,
    du_k,
    du_kp1,
    Δt_km1,
    Δt_k,
)
    i_km1 = 1.0 / Δt_km1
    i_k = 1.0 / Δt_k
    i_km1_2 = i_km1 * i_km1
    i_km1_3 = i_km1_2 * i_km1
    i_k_2 = i_k * i_k
    i_k_3 = i_k_2 * i_k
    return (
        6 * i_km1_2,
        -6 * i_km1_2 + 6 * i_k_2,
        -6 * i_k_2,
        2 * i_km1,
        4 * i_km1 + 4 * i_k,
        2 * i_k,
        12 * i_km1_3 * (u_k - u_km1) - 2 * i_km1_2 * (du_km1 + 2 * du_k),
        12 * i_k_3 * (u_kp1 - u_k) - 2 * i_k_2 * (2 * du_k + du_kp1),
    )
end

# ----------------------------------------------------------------------------- #
# Critical points
# ----------------------------------------------------------------------------- #

"""
    find_cubic_critical_points(u_k, u_kp1, du_k, du_kp1, Δt)

Find critical points of cubic Hermite spline in [0, 1] where du/dτ = 0.

Returns a vector of τ values in [0, 1] where the spline has extrema.

# Algorithm
The derivative du/dτ is a quadratic: aτ² + bτ + c = 0
Solve using quadratic formula and keep roots in (0, 1).

# Note
Generic implementation compatible with ForwardDiff.Dual numbers.
"""
function find_cubic_critical_points(u_k, u_kp1, du_k, du_kp1, Δt)
    # Coefficients of du/dτ = aτ² + bτ + c, from `hermite_derivative_basis`:
    #   du/dτ = dh00·u_k + dh10·Δt·du_k + dh01·u_kp1 + dh11·Δt·du_kp1

    # Collecting terms in τ²:
    a = 6*u_k + 3*Δt*du_k - 6*u_kp1 + 3*Δt*du_kp1

    # Collecting terms in τ:
    b = -6*u_k - 4*Δt*du_k + 6*u_kp1 - 2*Δt*du_kp1

    # Constant term:
    c = Δt * du_k

    # Use generic type for array (compatible with Dual numbers)
    T = promote_type(typeof(u_k), typeof(u_kp1), typeof(du_k), typeof(du_kp1))
    critical_points = T[]

    if abs(a) < 1e-12
        # Linear case: bτ + c = 0
        if abs(b) > 1e-12
            τ = -c / b
            if 0 < τ < 1
                push!(critical_points, τ)
            end
        end
    else
        # Quadratic case
        discriminant = b^2 - 4*a*c

        if discriminant >= 0
            sqrt_disc = sqrt(discriminant)
            τ1 = (-b + sqrt_disc) / (2*a)
            τ2 = (-b - sqrt_disc) / (2*a)

            # Only keep roots in (0, 1) - exclude endpoints since we already check them
            if 0 < τ1 < 1
                push!(critical_points, τ1)
            end
            if 0 < τ2 < 1
                push!(critical_points, τ2)
            end
        end
    end

    return critical_points
end

"""
    find_hermite_extrema(u_k, u_kp1, du_k, du_kp1, Δt)

Find all extrema (critical points + endpoints) and return their values.

Returns a vector of u values at all critical points and endpoints.
"""
function find_hermite_extrema(u_k, u_kp1, du_k, du_kp1, Δt)
    extrema_values = Float64[u_k, u_kp1]  # Always check endpoints

    τ_crits = find_cubic_critical_points(u_k, u_kp1, du_k, du_kp1, Δt)

    for τ in τ_crits
        u_τ = evaluate_hermite_spline(τ, u_k, u_kp1, du_k, du_kp1, Δt)
        push!(extrema_values, u_τ)
    end

    return extrema_values
end

"""
    find_slope_critical_points(u_k, u_kp1, du_k, du_kp1, Δt)

Find critical points of the derivative (where d²u/dt² = 0) in [0, 1].

Returns a vector of τ values where the slope has local extrema.

# Algorithm
The second derivative d²u/dτ² is linear: 2aτ + b = 0
where a and b are the quadratic coefficients from du/dτ.
"""
function find_slope_critical_points(u_k, u_kp1, du_k, du_kp1, Δt)
    # From find_cubic_critical_points, we have du/dτ = aτ² + bτ + c
    # So d²u/dτ² = 2aτ + b

    a = 6*u_k + 3*Δt*du_k - 6*u_kp1 + 3*Δt*du_kp1
    b = -6*u_k - 4*Δt*du_k + 6*u_kp1 - 2*Δt*du_kp1

    T = promote_type(typeof(u_k), typeof(u_kp1), typeof(du_k), typeof(du_kp1))
    critical_points = T[]

    if abs(a) > 1e-12
        # Linear d²u/dτ² = 0: 2aτ + b = 0
        τ = -b / (2*a)
        if 0 < τ < 1
            push!(critical_points, τ)
        end
    end
    # If a ≈ 0, derivative is constant (no extrema in interior)

    return critical_points
end

@testitem "Hermite primitives: one definition, endpoints, and the two accelerations" begin

    using Piccolo:
        evaluate_hermite_derivative,
        evaluate_hermite_spline,
        hermite_accel_end,
        hermite_accel_start,
        hermite_basis_functions,
        hermite_derivative_basis
    using ForwardDiff
    using Random

    Random.seed!(0x330ABC1)

    # Basis partition-of-unity / endpoint interpolation
    h0 = hermite_basis_functions(0.0)
    @test h0 == (1.0, 0.0, 0.0, 0.0)
    h1 = hermite_basis_functions(1.0)
    @test h1 == (0.0, 0.0, 1.0, 0.0)

    u_k, u_kp1, du_k, du_kp1, Δt = randn(), randn(), randn(), randn(), 0.17
    @test evaluate_hermite_spline(0.0, u_k, u_kp1, du_k, du_kp1, Δt) ≈ u_k
    @test evaluate_hermite_spline(1.0, u_k, u_kp1, du_k, du_kp1, Δt) ≈ u_kp1
    @test evaluate_hermite_derivative(0.0, u_k, u_kp1, du_k, du_kp1, Δt) ≈ du_k
    @test evaluate_hermite_derivative(1.0, u_k, u_kp1, du_k, du_kp1, Δt) ≈ du_kp1

    # The accelerations ARE the second derivative of the spline at τ = 0, 1 —
    # this is the identity the smooth-acceleration constraint's stencil rests on.
    s = τ -> evaluate_hermite_spline(τ, u_k, u_kp1, du_k, du_kp1, Δt)
    d2 = τ -> ForwardDiff.derivative(t -> ForwardDiff.derivative(s, t), τ) / Δt^2
    @test hermite_accel_start(u_k, u_kp1, du_k, du_kp1, Δt) ≈ d2(0.0)
    @test hermite_accel_end(u_k, u_kp1, du_k, du_kp1, Δt) ≈ d2(1.0)

    # dh basis is the τ-derivative of the h basis
    for τ in (0.0, 0.25, 0.5, 1.0)
        dh = hermite_derivative_basis(τ)
        dh_ad = ForwardDiff.derivative(t -> collect(hermite_basis_functions(t)), τ)
        @test all(isapprox.(collect(dh), dh_ad; atol = 1e-12))
    end
end
