# B-spline basis for control-point-native pulse parameterization.
# Port of Drake's math::BsplineBasis<T> (see drake/math/bspline_basis.h).
#
# Convention: `order = degree + 1` (Drake). A cubic B-spline has order=4, degree=3.
# Lives on normalized parameter τ ∈ [initial_parameter_value, final_parameter_value]
# — by convention, [0, 1] when consumed by SplineIntegrator's per-segment ODE.

module BsplineBases

using SparseArrays
using TestItems: @testitem

export BsplineBasis
export evaluate_curve, evaluate_linear_in_control_points, as_linear_in_control_points
export num_basis_functions, degree

"""
    BsplineBasis{T<:Real}

B-spline basis with order `k` (= degree + 1) and a non-descending knot vector.

# Constructor preconditions
- `order >= 2` (no degree-0 splines; ZeroOrderPulse handles that)
- `length(knots) >= 2 * order` (Drake's DRAKE_DEMAND at bspline_basis.h:148)
- `knots` non-descending
- `num_basis_functions = length(knots) - order` must be `>= order`

The valid parameter range is `[knots[order], knots[end - order + 1]]`.
For a clamped basis (knot multiplicity `order` at endpoints), this is `[t_min, t_max]`
and `evaluate_curve(b, C, t_min) == C[:, 1]`, `evaluate_curve(b, C, t_max) == C[:, end]`.
"""
struct BsplineBasis{T<:Real}
    order::Int                  # k = degree + 1
    knots::Vector{T}            # non-descending, length = num_basis_functions + order
    # Caches (computed lazily; depend only on basis, not on control points):
    _as_linear_cache::Dict{Int,SparseMatrixCSC{T,Int}}

    function BsplineBasis{T}(order::Int, knots::Vector{T}) where {T<:Real}
        order >= 2 || throw(ArgumentError(
            "BsplineBasis requires order >= 2 (got order=$order); for order=1 use LinearSplinePulse"))
        length(knots) >= 2 * order || throw(ArgumentError(
            "BsplineBasis requires length(knots) >= 2 * order (got length=$(length(knots)), order=$order)"))
        issorted(knots) || throw(ArgumentError(
            "BsplineBasis knots must be non-descending"))
        nbf = length(knots) - order
        nbf >= order || throw(ArgumentError(
            "BsplineBasis num_basis_functions ($nbf) must be >= order ($order)"))
        new{T}(order, knots, Dict{Int,SparseMatrixCSC{T,Int}}())
    end
end

BsplineBasis(order::Int, knots::Vector{T}) where {T<:Real} =
    BsplineBasis{T}(order, knots)

"""
    BsplineBasis(order, num_basis_functions; knot_vector_type, initial_parameter_value, final_parameter_value)

Convenience constructor that auto-generates a knot vector. Default `:clamped_uniform`
gives endpoint interpolation. `:uniform` gives a non-clamped uniform vector (rare;
endpoints not interpolated).
"""
function BsplineBasis(
    order::Int,
    num_basis_functions::Int;
    knot_vector_type::Symbol = :clamped_uniform,
    initial_parameter_value::Real = 0.0,
    final_parameter_value::Real = 1.0,
)
    num_basis_functions >= order || throw(ArgumentError(
        "num_basis_functions ($num_basis_functions) must be >= order ($order)"))
    initial_parameter_value < final_parameter_value || throw(ArgumentError(
        "initial_parameter_value ($initial_parameter_value) must be < final ($final_parameter_value)"))

    T = promote_type(typeof(float(initial_parameter_value)), typeof(float(final_parameter_value)))
    a = T(initial_parameter_value)
    b = T(final_parameter_value)
    nbf = num_basis_functions
    k = order

    knots = if knot_vector_type == :clamped_uniform
        # Multiplicity k at each endpoint, uniform interior.
        n_interior = nbf - k + 1  # number of *unique* interior parameter values
        if n_interior == 1
            # Single Bézier segment: no interior knots
            vcat(fill(a, k), fill(b, k))
        else
            # Interior knots evenly spaced strictly between a and b.
            # range(a, b, length = n_interior + 1)[2:end-1] yields (n_interior - 1) = nbf - k interior knots.
            interior = collect(range(a, b; length = n_interior + 1))[2:end-1]
            # Total knots = k (left clamp) + (nbf - k) interior + k (right clamp) = nbf + k ✓
            vcat(fill(a, k), interior, fill(b, k))
        end
    elseif knot_vector_type == :uniform
        # Non-clamped uniform: evenly spaced across [a, b]
        collect(range(a, b; length = nbf + k))
    else
        throw(ArgumentError(
            "Unsupported knot_vector_type=$knot_vector_type; use :clamped_uniform or :uniform"))
    end
    return BsplineBasis{T}(k, knots)
end

@testitem "BsplineBasis: constructor preconditions" begin
    using Piccolo
    # Reject degree 0 / order 1
    @test_throws ArgumentError BsplineBasis(1, 10)
    # Reject too-few knots
    @test_throws ArgumentError BsplineBasis(4, collect(0.0:0.1:0.5))  # length < 2*order
    # Reject non-descending knots
    @test_throws ArgumentError BsplineBasis(4, [0.0, 0.1, 0.05, 0.2, 0.3, 0.4, 0.5, 0.6])
    # Reject n_basis_functions < order
    @test_throws ArgumentError BsplineBasis(4, 3)
    # Convenience constructor: bad arg order
    @test_throws ArgumentError BsplineBasis(4, 8; initial_parameter_value = 1.0, final_parameter_value = 0.0)
    # Unknown knot type
    @test_throws ArgumentError BsplineBasis(4, 8; knot_vector_type = :foo)

    # Sanity: clamped cubic with 8 CPs
    b = BsplineBasis(4, 8)
    @test b.order == 4
    @test length(b.knots) == 8 + 4
    @test b.knots[1:4] == zeros(4)
    @test b.knots[end-3:end] == ones(4)
end

end # module BsplineBases
