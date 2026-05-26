module Pulses

"""
Pulse types for quantum control.

Provides interpolated and analytic pulse types:
- `ZeroOrderPulse`: Piecewise constant (zero-order hold)
- `LinearSplinePulse`: Linear interpolation between samples
- `CubicSplinePulse`: Cubic spline interpolation
- `GaussianPulse`: Analytic Gaussian envelope
- `ErfPulse`: Analytic error function (erf) envelope for phase compensation
- `MSPhaseCompensationPulse`: Scaled erf pulse [0,1] for MS gate phase compensation (matches IEEE TQE 2024 paper)
- `CompositePulse`: Combine multiple pulses by interleaving their drives

All pulses are callable: `pulse(t)` returns the control vector at time `t`.
"""

export AbstractPulse, AbstractSplinePulse
export ZeroOrderPulse,
    LinearSplinePulse,
    CubicSplinePulse,
    BSplinePulse,
    GaussianPulse,
    ErfPulse,
    CompositePulse,
    FunctionPulse
export duration, n_drives, sample, drive_name
export load_pulse
export get_control_points, get_knots

using DataInterpolations
using DataInterpolations:
    ConstantInterpolation, LinearInterpolation, CubicHermiteSpline, ExtrapolationType
using ForwardDiff
using JLD2
using LinearAlgebra: norm
using SpecialFunctions: erf
using TestItems

# Sibling submodule for B-spline basis support (BSplinePulse)
using ..BsplineBases:
    BsplineBasis,
    evaluate_curve,
    evaluate_linear_in_control_points,
    num_basis_functions

const DEFAULT_SAMPLES::Int = 100

# ============================================================================ #
# Abstract type
# ============================================================================ #

"""
    AbstractPulse

Abstract type for all pulse types. All pulses are callable: `pulse(t)` returns
the control vector at time `t`.
"""
abstract type AbstractPulse end

"""
    AbstractSplinePulse <: AbstractPulse

Abstract type for spline-based pulses (linear and cubic interpolation).
These pulses use the spline coefficients as optimization variables.
"""
abstract type AbstractSplinePulse <: AbstractPulse end

# Make all pulses callable
(pulse::AbstractPulse)(t::Real) = evaluate(pulse, t)

"""
    duration(pulse::AbstractPulse)

Return the duration of the pulse.
"""
duration(pulse::AbstractPulse) = pulse.duration

"""
    n_drives(pulse::AbstractPulse)

Return the number of control drives in the pulse.
"""
n_drives(pulse::AbstractPulse) = pulse.n_drives

"""
    drive_name(pulse::AbstractPulse)

Return the name of the drive variable for this pulse.
"""
drive_name(pulse::AbstractPulse) = pulse.drive_name

"""
    sample(pulse::AbstractPulse, times::AbstractVector)

Sample the pulse at the given times. Returns a matrix of size `(n_drives, length(times))`.
"""
function sample(pulse::AbstractPulse, times::AbstractVector)
    return hcat([pulse(t) for t in times]...)
end

"""
    sample(pulse::AbstractPulse, n_samples::Int)

Sample the pulse uniformly with `n_samples` points. Returns `(controls, times)`.
"""
function sample(pulse::AbstractPulse, n_samples::Int)
    times = collect(range(0.0, duration(pulse), length = n_samples))
    return sample(pulse, times), times
end

"""
    sample(bounds::AbstractVector{Tuple}, n_samples::Int)

Sample the bounds uniformly with `n_samples` points.
"""
function sample(bounds::AbstractVector{Tuple{Float64,Float64}}, n_samples::Int)
    return [(high - low) * rand() + low for (low, high) in bounds, _ = 1:n_samples]
end

"""
    derivative(pulse::AbstractPulse, time::Real)

Return the time derivative of the pulse at `time`.

Falls back to ForwardDiff for generic pulses. Specialized pulse types may benefit from custom derivative methods (e.g., DataInterpolations.derivative) for better handling of knot points and boundaries.
"""
function derivative(pulse::AbstractPulse, time::Real)
    return ForwardDiff.derivative(pulse, time)
end

Base.summary(io::IO, p::AbstractPulse) = print(
    io,
    "$(nameof(typeof(p)))(Number of drives = ",
    n_drives(p),
    ", ",
    "T = ",
    duration(p),
    ")",
)

Base.show(io::IO, p::AbstractPulse) = summary(io, p)

function Base.show(io::IO, ::MIME"text/plain", p::AbstractPulse)
    println(io, nameof(typeof(p)))
    println(io, "  drives: ", n_drives(p))
    print(io, "  duration: ", duration(p))
end

# ============================================================================ #
# ZeroOrderPulse (piecewise constant / zero-order hold)
# ============================================================================ #

"""
    ZeroOrderPulse{I<:ConstantInterpolation} <: AbstractPulse

Piecewise constant pulse (zero-order hold). The control value at time `t` is
the value at the most recent sample point.

# Fields
- `controls::I`: ConstantInterpolation from DataInterpolations
- `duration::Float64`: Total pulse duration
- `n_drives::Int`: Number of control drives
- `drive_name::Symbol`: Name of the drive variable (default `:u`)
- `initial_value::Vector{Float64}`: Initial boundary condition (default: zeros)
- `final_value::Vector{Float64}`: Final boundary condition (default: zeros)
"""
struct ZeroOrderPulse{I<:ConstantInterpolation} <: AbstractPulse
    controls::I
    duration::Float64
    n_drives::Int
    drive_name::Symbol
    initial_value::Vector{Float64}
    final_value::Vector{Float64}
    snap_to_knots::Bool
end

"""
    ZeroOrderPulse(controls::AbstractMatrix, times::AbstractVector; drive_name=:u, initial_value=nothing, final_value=nothing, snap_to_knots=true)

Create a zero-order hold pulse from control samples and times.

# Arguments
- `controls`: Matrix of size `(n_drives, n_times)` with control values
- `times`: Vector of sample times (must start at 0)

# Keyword Arguments
- `drive_name`: Name of the drive variable (default `:u`)
- `initial_value`: Initial boundary condition (default: zeros(n_drives))
- `final_value`: Final boundary condition (default: zeros(n_drives))
- `snap_to_knots`: When `true` (default), `evaluate` snaps query times within
  `1e-12` of a stored knot time to that knot's value. This avoids the
  off-by-one that arises when recomputed times (e.g. from `range()` or a
  `cumsum` round-trip through `NamedTrajectory`) differ from stored times by
  float roundoff at `ConstantInterpolation`'s left-continuous discontinuities,
  silently shifting the entire sampled control vector by one knot.
  Pass `snap_to_knots=false` to opt back into the raw left-continuous
  `ConstantInterpolation` semantics.
"""
function ZeroOrderPulse(
    controls::AbstractMatrix,
    times::AbstractVector;
    drive_name::Symbol = :u,
    initial_value::Union{Nothing,Symbol,Vector{<:Real}} = nothing,
    final_value::Union{Nothing,Symbol,Vector{<:Real}} = nothing,
    snap_to_knots::Bool = true,
)
    n_drives = size(controls, 1)
    # :free means "no boundary constraint" (stored as NaN sentinel)
    # nothing defaults to zeros for backward compatibility
    init_val = if initial_value === :free
        fill(NaN, n_drives)
    elseif initial_value isa Symbol
        error("Unknown initial_value symbol: $initial_value (only :free is supported)")
    elseif isnothing(initial_value)
        zeros(n_drives)
    else
        Vector{Float64}(initial_value)
    end
    final_val = if final_value === :free
        fill(NaN, n_drives)
    elseif final_value isa Symbol
        error("Unknown final_value symbol: $final_value (only :free is supported)")
    elseif isnothing(final_value)
        zeros(n_drives)
    else
        Vector{Float64}(final_value)
    end
    # Materialize to Matrix/Vector to ensure consistent type parameters.
    # Constant extrapolation on both sides matches the zero-order-hold semantics
    # and lets ODE integrators query slightly outside [t₁, tₙ] (e.g. at
    # `t = tₙ + ε` due to step rounding) without throwing.
    interp = ConstantInterpolation(
        Matrix(controls),
        collect(times);
        extrapolation = ExtrapolationType.Constant,
    )
    return ZeroOrderPulse(
        interp,
        Float64(times[end]),
        n_drives,
        drive_name,
        init_val,
        final_val,
        snap_to_knots,
    )
end

derivative(p::ZeroOrderPulse, t::Real) = DataInterpolations.derivative(p.controls, t)

function evaluate(p::ZeroOrderPulse, t)
    if p.snap_to_knots
        knots = p.controls.t
        idx = searchsortedfirst(knots, t)
        if idx <= length(knots) && abs(t - knots[idx]) < 1e-12
            return p.controls.u[:, idx]
        end
    end
    return p.controls(t)
end

function sample(pulse::ZeroOrderPulse, times::AbstractVector)
    knot_times = collect(pulse.controls.t)
    is_native =
        length(times) == length(knot_times) &&
        all(isapprox.(times, knot_times; atol = 1e-12))
    if is_native
        return Matrix(pulse.controls.u)
    else
        return hcat([pulse(t) for t in times]...)
    end
end

# ============================================================================ #
# LinearSplinePulse
# ============================================================================ #

"""
    LinearSplinePulse{I<:LinearInterpolation} <: AbstractSplinePulse

Pulse with linear interpolation between sample points.

# Fields
- `controls::I`: LinearInterpolation from DataInterpolations
- `duration::Float64`: Total pulse duration
- `n_drives::Int`: Number of control drives
- `drive_name::Symbol`: Name of the drive variable (default `:u`)
- `initial_value::Vector{Float64}`: Initial boundary condition (default: zeros)
- `final_value::Vector{Float64}`: Final boundary condition (default: zeros)
"""
struct LinearSplinePulse{I<:LinearInterpolation} <: AbstractSplinePulse
    controls::I
    duration::Float64
    n_drives::Int
    drive_name::Symbol
    initial_value::Vector{Float64}
    final_value::Vector{Float64}
end

"""
    LinearSplinePulse(controls::AbstractMatrix, times::AbstractVector; drive_name=:u, initial_value=nothing, final_value=nothing)

Create a linearly interpolated pulse from control samples and times.

# Arguments
- `controls`: Matrix of size `(n_drives, n_times)` with control values
- `times`: Vector of sample times (must start at 0)

# Keyword Arguments
- `drive_name`: Name of the drive variable (default `:u`)
- `initial_value`: Initial boundary condition (default: zeros(n_drives))
- `final_value`: Final boundary condition (default: zeros(n_drives))
"""
function LinearSplinePulse(
    controls::AbstractMatrix,
    times::AbstractVector;
    drive_name::Symbol = :u,
    initial_value::Union{Nothing,Symbol,Vector{<:Real}} = nothing,
    final_value::Union{Nothing,Symbol,Vector{<:Real}} = nothing,
)
    n_drives = size(controls, 1)
    init_val = if initial_value === :free
        fill(NaN, n_drives)
    elseif initial_value isa Symbol
        error("Unknown initial_value symbol: $initial_value (only :free is supported)")
    elseif isnothing(initial_value)
        zeros(n_drives)
    else
        Vector{Float64}(initial_value)
    end
    final_val = if final_value === :free
        fill(NaN, n_drives)
    elseif final_value isa Symbol
        error("Unknown final_value symbol: $final_value (only :free is supported)")
    elseif isnothing(final_value)
        zeros(n_drives)
    else
        Vector{Float64}(final_value)
    end
    # Materialize to Matrix/Vector to ensure consistent type parameters.
    # Constant extrapolation prevents DataInterpolations from throwing when
    # ODE integrators query slightly outside [t₁, tₙ] due to step rounding.
    interp = LinearInterpolation(
        Matrix(controls),
        collect(times);
        extrapolation = ExtrapolationType.Constant,
    )
    return LinearSplinePulse(
        interp,
        Float64(times[end]),
        n_drives,
        drive_name,
        init_val,
        final_val,
    )
end

derivative(p::LinearSplinePulse, t::Real) = DataInterpolations.derivative(p.controls, t)

evaluate(p::LinearSplinePulse, t) = p.controls(t)

# ============================================================================ #
# CubicSplinePulse (Hermite spline with explicit derivatives)
# ============================================================================ #

"""
    CubicSplinePulse{I<:CubicHermiteSpline} <: AbstractPulse

Pulse with cubic Hermite spline interpolation. Uses both control values AND 
derivatives for exact reconstruction after optimization.

# Fields
- `controls::I`: CubicHermiteSpline from DataInterpolations
- `duration::Float64`: Total pulse duration
- `n_drives::Int`: Number of control drives
- `drive_name::Symbol`: Name of the drive variable (default `:u`)
- `initial_value::Vector{Float64}`: Initial boundary condition (default: zeros)
- `final_value::Vector{Float64}`: Final boundary condition (default: zeros)
"""
struct CubicSplinePulse{I<:CubicHermiteSpline} <: AbstractSplinePulse
    controls::I
    duration::Float64
    n_drives::Int
    drive_name::Symbol
    initial_value::Vector{Float64}
    final_value::Vector{Float64}
end

"""
    CubicSplinePulse(controls::AbstractMatrix, derivatives::AbstractMatrix, times::AbstractVector; drive_name=:u, initial_value=nothing, final_value=nothing)

Create a cubic Hermite spline pulse from control values, derivatives, and times.

# Arguments
- `controls`: Matrix of size `(n_drives, n_times)` with control values
- `derivatives`: Matrix of size `(n_drives, n_times)` with control derivatives
- `times`: Vector of sample times (must start at 0)

# Keyword Arguments
- `drive_name`: Name of the drive variable (default `:u`)
- `initial_value`: Initial boundary condition (default: zeros(n_drives))
- `final_value`: Final boundary condition (default: zeros(n_drives))
"""
function CubicSplinePulse(
    controls::AbstractMatrix,
    derivatives::AbstractMatrix,
    times::AbstractVector;
    drive_name::Symbol = :u,
    initial_value::Union{Nothing,Symbol,Vector{<:Real}} = nothing,
    final_value::Union{Nothing,Symbol,Vector{<:Real}} = nothing,
)
    n_drives = size(controls, 1)
    init_val = if initial_value === :free
        fill(NaN, n_drives)
    elseif initial_value isa Symbol
        error("Unknown initial_value symbol: $initial_value (only :free is supported)")
    elseif isnothing(initial_value)
        zeros(n_drives)
    else
        Vector{Float64}(initial_value)
    end
    final_val = if final_value === :free
        fill(NaN, n_drives)
    elseif final_value isa Symbol
        error("Unknown final_value symbol: $final_value (only :free is supported)")
    elseif isnothing(final_value)
        zeros(n_drives)
    else
        Vector{Float64}(final_value)
    end
    # Materialize to Matrix to ensure consistent type parameters across construction methods.
    # Constant extrapolation prevents DataInterpolations from throwing when
    # ODE integrators query slightly outside [t₁, tₙ] due to step rounding.
    interp = CubicHermiteSpline(
        Matrix(derivatives),
        Matrix(controls),
        collect(times);
        extrapolation = ExtrapolationType.Constant,
    )
    return CubicSplinePulse(
        interp,
        Float64(times[end]),
        n_drives,
        drive_name,
        init_val,
        final_val,
    )
end

"""
    CubicSplinePulse(controls::AbstractMatrix, times::AbstractVector; drive_name=:u, initial_value=nothing, final_value=nothing)

Create a cubic Hermite spline pulse with zero derivatives at all knot points.
Useful for initial guesses where smoothness constraints will be enforced by optimizer.

# Arguments
- `controls`: Matrix of size `(n_drives, n_times)` with control values
- `times`: Vector of sample times (must start at 0)

# Keyword Arguments
- `drive_name`: Name of the drive variable (default `:u`)
- `initial_value`: Initial boundary condition (default: zeros(n_drives))
- `final_value`: Final boundary condition (default: zeros(n_drives))
"""
function CubicSplinePulse(
    controls::AbstractMatrix,
    times::AbstractVector;
    drive_name::Symbol = :u,
    initial_value::Union{Nothing,Vector{<:Real}} = nothing,
    final_value::Union{Nothing,Vector{<:Real}} = nothing,
)
    derivatives = zeros(size(controls))
    return CubicSplinePulse(
        controls,
        derivatives,
        times;
        drive_name,
        initial_value,
        final_value,
    )
end

# TODO: Unimplemented by DataInterpolations
# derivative(p::CubicSplinePulse, t::Real) = DataInterpolations.derivative(p.controls, t)

evaluate(p::CubicSplinePulse, t) = p.controls(t)

# ============================================================================ #
# Spline pulse knot accessors
# ============================================================================ #

"""
    get_knot_times(pulse::AbstractPulse)

Return the knot/discontinuity times for the pulse.

For piecewise pulses (ZeroOrderPulse, spline pulses), these are the sample times
where the control value or its derivative may be discontinuous.
For smooth analytic pulses (GaussianPulse, ErfPulse), returns just the endpoints.
For CompositePulse, returns the sorted union of all sub-pulse knot times.
"""
get_knot_times(p::ZeroOrderPulse) = p.controls.t
get_knot_times(p::LinearSplinePulse) = p.controls.t
get_knot_times(p::CubicSplinePulse) = p.controls.t
# get_knot_times for BSplinePulse is defined below, after the type declaration.

"""
    get_knot_count(pulse::AbstractSplinePulse)

Return the number of knots in the spline pulse.
"""
get_knot_count(p::AbstractSplinePulse) = length(get_knot_times(p))

"""
    get_knot_values(pulse::CubicSplinePulse)

Return the control values at knot points (the `u` matrix).
"""
get_knot_values(p::LinearSplinePulse) = p.controls.u
get_knot_values(p::CubicSplinePulse) = p.controls.u

"""
    get_knot_derivatives(pulse::CubicSplinePulse)

Return the Hermite tangents at knot points (the `du` matrix).
Only available for CubicSplinePulse.
"""
get_knot_derivatives(p::CubicSplinePulse) = p.controls.du

export get_knot_times, get_knot_count, get_knot_values, get_knot_derivatives

# ============================================================================ #
# Conversion methods (analytic → spline pulses)
# ============================================================================ #

"""
    LinearSplinePulse(pulse::AbstractPulse, n_samples::Int; kwargs...)
    LinearSplinePulse(pulse::AbstractPulse, times::AbstractVector; kwargs...)

Convert any pulse to a LinearSplinePulse by sampling at specified times.

Useful for initializing optimization problems with analytic pulse shapes.

# Arguments
- `pulse`: Source pulse (GaussianPulse, ErfPulse, CompositePulse, etc.)
- `n_samples`: Number of uniformly spaced samples (alternative to `times`)
- `times`: Specific sample times (alternative to `n_samples`)

# Keyword Arguments
- `drive_name`: Name for the drive variable (default: `:u`)
- `initial_value`: Initial boundary condition (default: pulse(0.0))
- `final_value`: Final boundary condition (default: pulse(duration))

# Example
```julia
gaussian = GaussianPulse([1.0, 2.0], 0.1, 1.0)
linear = LinearSplinePulse(gaussian, 50)  # 50 samples
```
"""
function LinearSplinePulse(
    pulse::AbstractPulse,
    n_samples::Int;
    drive_name::Symbol = :u,
    initial_value::Union{Nothing,Vector{<:Real}} = nothing,
    final_value::Union{Nothing,Vector{<:Real}} = nothing,
)
    times = collect(range(0.0, duration(pulse), length = n_samples))
    return LinearSplinePulse(pulse, times; drive_name, initial_value, final_value)
end

function LinearSplinePulse(
    pulse::AbstractPulse,
    times::AbstractVector;
    drive_name::Symbol = :u,
    initial_value::Union{Nothing,Vector{<:Real}} = nothing,
    final_value::Union{Nothing,Vector{<:Real}} = nothing,
)
    controls = sample(pulse, times)
    init_val = isnothing(initial_value) ? pulse(times[1]) : initial_value
    final_val = isnothing(final_value) ? pulse(times[end]) : final_value
    return LinearSplinePulse(
        controls,
        times;
        drive_name,
        initial_value = init_val,
        final_value = final_val,
    )
end

"""
    CubicSplinePulse(pulse::AbstractPulse, n_samples::Int; kwargs...)
    CubicSplinePulse(pulse::AbstractPulse, times::AbstractVector; kwargs...)

Convert any pulse to a CubicSplinePulse by sampling at specified times.
Derivatives are computed using ForwardDiff for automatic differentiation.

Useful for initializing optimization problems with smooth analytic pulse shapes.

# Arguments
- `pulse`: Source pulse (GaussianPulse, ErfPulse, CompositePulse, etc.)
- `n_samples`: Number of uniformly spaced samples (alternative to `times`)
- `times`: Specific sample times (alternative to `n_samples`)

# Keyword Arguments
- `drive_name`: Name for the drive variable (default: `:du`)
- `initial_value`: Initial boundary condition (default: pulse(0.0))
- `final_value`: Final boundary condition (default: pulse(duration))

# Example
```julia
gaussian = GaussianPulse([1.0, 2.0], 0.1, 1.0)
cubic = CubicSplinePulse(gaussian, 50)  # 50 samples with ForwardDiff derivatives
```
"""
function CubicSplinePulse(
    pulse::AbstractPulse,
    n_samples::Int;
    drive_name::Symbol = :du,
    initial_value::Union{Nothing,Vector{<:Real}} = nothing,
    final_value::Union{Nothing,Vector{<:Real}} = nothing,
)
    times = collect(range(0.0, duration(pulse), length = n_samples))
    return CubicSplinePulse(pulse, times; drive_name, initial_value, final_value)
end

function CubicSplinePulse(
    pulse::AbstractPulse,
    times::AbstractVector;
    drive_name::Symbol = :du,
    initial_value::Union{Nothing,Vector{<:Real}} = nothing,
    final_value::Union{Nothing,Vector{<:Real}} = nothing,
)
    controls = sample(pulse, times)
    derivatives = stack(map(t -> derivative(pulse, t), times))

    init_val = isnothing(initial_value) ? pulse(times[1]) : initial_value
    final_val = isnothing(final_value) ? pulse(times[end]) : final_value

    return CubicSplinePulse(
        controls,
        derivatives,
        times;
        drive_name,
        initial_value = init_val,
        final_value = final_val,
    )
end

# ============================================================================ #
# BSplinePulse
# ============================================================================ #

"""
    BSplinePulse{T, B<:BsplineBasis{T}} <: AbstractSplinePulse

Control-point-native B-spline pulse parameterization. Decision variables are
control points; the curve `u(t) = control_points * M_0(τ)` lies in the convex
hull of its control points. C^{degree-1} continuity is structural (no
regularizer needed).

# Fields
- `basis::B`: defines order, knots, normalized parameter τ ∈ [0, 1]
- `control_points::Matrix{T}`: `n_drives × num_basis_functions(basis)`
- `duration::Float64`
- `n_drives::Int`
- `drive_name::Symbol`
- `initial_value::Union{Vector{Float64},Symbol}`: clamped to `control_points[:, 1]` or `:free`
- `final_value::Union{Vector{Float64},Symbol}`: clamped to `control_points[:, end]` or `:free`

# Convention
The basis lives on normalized parameter τ ∈ [0, 1]; see the spec's Conventions
section for free-time (per-segment Δt) compatibility.
"""
struct BSplinePulse{T<:Real,B<:BsplineBasis{T}} <: AbstractSplinePulse
    basis::B
    control_points::Matrix{T}
    duration::Float64
    n_drives::Int
    drive_name::Symbol
    initial_value::Union{Vector{Float64},Symbol}
    final_value::Union{Vector{Float64},Symbol}
end

# Basis knots live on normalized τ ∈ [0,1]; map to physical time
# and dedupe the clamped boundary repetitions.
get_knot_times(p::BSplinePulse) = unique(p.basis.knots) .* p.duration

function _resolve_bspline_boundary(spec, cps::AbstractMatrix, which::Symbol)
    if spec === nothing
        col = which == :first ? cps[:, 1] : cps[:, end]
        return Vector{Float64}(col)
    elseif spec === :free
        return :free
    elseif spec isa Symbol
        throw(ArgumentError("Invalid boundary symbol: $spec (only :free is supported)"))
    elseif spec isa AbstractVector
        return Vector{Float64}(spec)
    else
        throw(ArgumentError("Invalid boundary spec: $spec"))
    end
end

"""
    BSplinePulse(control_points::AbstractMatrix, times::AbstractVector;
                 degree=3, drive_name=:u, initial_value=nothing, final_value=nothing)

Primary constructor. The `times` vector is used only for its duration
(`times[end] - times[1]`) — the actual time parameterization lives in
the basis on normalized τ ∈ [0, 1].

`initial_value` / `final_value`:
- `nothing` (default): set to `control_points[:, 1]` / `control_points[:, end]`
- A vector: enforces invariant `control_points[:, 1] ≈ initial_value` (atol 1e-12)
- `:free`: no pinning; the boundary CP is a free decision variable
"""
function BSplinePulse(
    control_points::AbstractMatrix{T},
    times::AbstractVector{<:Real};
    degree::Int = 3,
    drive_name::Symbol = :u,
    initial_value::Union{Nothing,AbstractVector,Symbol} = nothing,
    final_value::Union{Nothing,AbstractVector,Symbol} = nothing,
) where {T<:Real}
    n_drv = size(control_points, 1)
    nbf = size(control_points, 2)
    nbf >= degree + 1 || throw(ArgumentError(
        "BSplinePulse requires n_control_points ($nbf) >= degree + 1 ($(degree+1))"))
    length(times) >= 2 || throw(ArgumentError(
        "BSplinePulse requires at least 2 time entries"))

    dur = Float64(times[end] - times[1])
    dur > 0 || throw(ArgumentError("BSplinePulse duration must be positive"))

    iv = _resolve_bspline_boundary(initial_value, control_points, :first)
    fv = _resolve_bspline_boundary(final_value, control_points, :last)

    # Invariant check (only when iv/fv are concrete vectors and user-supplied)
    if initial_value isa AbstractVector
        isapprox(control_points[:, 1], iv; atol = 1e-12) || throw(ArgumentError(
            "BSplinePulse: initial_value $iv must match control_points[:, 1] = $(control_points[:, 1])"))
    end
    if final_value isa AbstractVector
        isapprox(control_points[:, end], fv; atol = 1e-12) || throw(ArgumentError(
            "BSplinePulse: final_value $fv must match control_points[:, end] = $(control_points[:, end])"))
    end

    Tp = promote_type(T, Float64)
    cps_p = Matrix{Tp}(control_points)
    # Build a basis with matching T parameter for compatibility with cps storage type
    basis_p = BsplineBasis(degree + 1, nbf;
        knot_vector_type = :clamped_uniform,
        initial_parameter_value = zero(Tp),
        final_parameter_value = one(Tp),
    )
    return BSplinePulse{Tp,typeof(basis_p)}(
        basis_p, cps_p, dur, n_drv, drive_name, iv, fv,
    )
end

"""
    BSplinePulse(pulse::AbstractPulse, n_control_points::Int; degree=3)

Warm-start constructor. Samples `pulse` densely, then uses
`DataInterpolations.BSplineApprox` to compute a least-squares-fit set of
control points. Inherits duration / n_drives / drive_name from `pulse`.
"""
function BSplinePulse(
    pulse::AbstractPulse,
    n_control_points::Int;
    degree::Int = 3,
)
    n_control_points >= degree + 1 || throw(ArgumentError(
        "n_control_points must be >= degree + 1"))

    T_dur = duration(pulse)
    n_drv = n_drives(pulse)
    n_samples = max(4 * n_control_points, 200)
    sample_ts = collect(range(0.0, T_dur, length = n_samples))

    samples = Matrix{Float64}(undef, n_drv, n_samples)
    for (k, t) in enumerate(sample_ts)
        samples[:, k] = evaluate(pulse, t)
    end

    C = Matrix{Float64}(undef, n_drv, n_control_points)
    for j in 1:n_drv
        approx = DataInterpolations.BSplineApprox(
            samples[j, :], sample_ts, degree, n_control_points, :Uniform, :Uniform,
        )
        C[j, :] = approx.c
    end

    times = collect(range(0.0, T_dur, length = n_control_points))
    return BSplinePulse(C, times; degree = degree, drive_name = drive_name(pulse))
end

# --- evaluate / derivative / accessors ---

function evaluate(p::BSplinePulse, t::Real)
    τ = t / p.duration
    return evaluate_curve(p.basis, p.control_points, τ)
end

function derivative(p::BSplinePulse, t::Real)
    τ = t / p.duration
    M1 = evaluate_linear_in_control_points(p.basis, τ; derivative_order = 1)
    # Chain rule: du/dt = (du/dτ) * (dτ/dt) = (C * M1) / duration
    return (p.control_points * M1) ./ p.duration
end

get_control_points(p::BSplinePulse) = p.control_points
get_knots(p::BSplinePulse) = p.basis.knots

# ============================================================================ #
# GaussianPulse (analytic)
# ============================================================================ #

"""
    GaussianPulse{F<:Function} <: AbstractPulse

Analytic Gaussian pulse. Each drive has its own amplitude, width (sigma), and center.

    u_i(t) = amplitudes[i] * exp(-(t - centers[i])² / (2 * sigmas[i]²))

# Fields
- `f::F`: Function that evaluates the pulse
- `amplitudes::Vector{Float64}`: Peak amplitude for each drive
- `sigmas::Vector{Float64}`: Gaussian width for each drive
- `centers::Vector{Float64}`: Center time for each drive
- `duration::Float64`: Total pulse duration
- `n_drives::Int`: Number of control drives
"""
struct GaussianPulse{F<:Function} <: AbstractPulse
    f::F
    amplitudes::Vector{Float64}
    sigmas::Vector{Float64}
    centers::Vector{Float64}
    duration::Float64
    n_drives::Int
end

"""
    GaussianPulse(amplitudes, sigmas, centers, duration)

Create a Gaussian pulse with per-drive parameters.

# Arguments
- `amplitudes`: Peak amplitude for each drive
- `sigmas`: Gaussian width (standard deviation) for each drive
- `centers`: Center time for each drive
- `duration`: Total pulse duration
"""
function GaussianPulse(
    amplitudes::AbstractVector{<:Real},
    sigmas::AbstractVector{<:Real},
    centers::AbstractVector{<:Real},
    duration::Real,
)
    n = length(amplitudes)
    @assert length(sigmas) == n "sigmas must have same length as amplitudes"
    @assert length(centers) == n "centers must have same length as amplitudes"

    amps = Vector{Float64}(amplitudes)
    sigs = Vector{Float64}(sigmas)
    ctrs = Vector{Float64}(centers)
    dur = Float64(duration)

    f = t -> [amps[i] * exp(-(t - ctrs[i])^2 / (2 * sigs[i]^2)) for i = 1:n]
    return GaussianPulse(f, amps, sigs, ctrs, dur, n)
end

"""
    GaussianPulse(amplitudes, sigma, duration; center=duration/2)

Create a Gaussian pulse with shared sigma and center across all drives.

# Arguments
- `amplitudes`: Peak amplitude for each drive
- `sigma`: Shared Gaussian width for all drives
- `duration`: Total pulse duration

# Keyword Arguments
- `center`: Shared center time (default: `duration/2`)
"""
function GaussianPulse(
    amplitudes::AbstractVector{<:Real},
    sigma::Real,
    duration::Real;
    center::Real = duration / 2,
)
    n = length(amplitudes)
    return GaussianPulse(
        amplitudes,
        fill(Float64(sigma), n),
        fill(Float64(center), n),
        duration,
    )
end

evaluate(p::GaussianPulse, t) = p.f(t)

# ============================================================================ #
# ErfPulse (analytic error function)
# ============================================================================ #

"""
    ErfPulse{F<:Function} <: AbstractPulse

Analytic error function pulse for phase compensation in trapped ion gates.

The error function profile is commonly used to compensate AC Stark shifts in 
Mølmer-Sørensen gates, where φ(t) ∝ erf(√2 (t - t₀)/σ) cancels time-varying
phases from off-resonant spectator modes.

    u_i(t) = amplitudes[i] * erf(√2 * (t - centers[i]) / sigmas[i])

Typically scaled to range [0, 1] or [-1, 1] by adjusting amplitude.

# Fields
- `f::F`: Function that evaluates the pulse
- `amplitudes::Vector{Float64}`: Peak amplitude for each drive
- `sigmas::Vector{Float64}`: Width parameter for each drive
- `centers::Vector{Float64}`: Center time for each drive
- `duration::Float64`: Total pulse duration
- `n_drives::Int`: Number of control drives

# References
- Mizrahi et al., "Realization and Calibration of Continuously Parameterized 
  Two-Qubit Gates...", IEEE TQE (2024), Figure 7b
"""
struct ErfPulse{F<:Function} <: AbstractPulse
    f::F
    amplitudes::Vector{Float64}
    sigmas::Vector{Float64}
    centers::Vector{Float64}
    duration::Float64
    n_drives::Int
end

"""
    ErfPulse(amplitudes, sigmas, centers, duration)

Create an error function pulse with per-drive parameters.

# Arguments
- `amplitudes`: Peak amplitude for each drive
- `sigmas`: Width parameter for each drive (controls steepness)
- `centers`: Center time for each drive (inflection point)
- `duration`: Total pulse duration

# Example
```julia
using SpecialFunctions: erf

# Phase compensation for MS gate
φ_max = π/4  # Maximum phase shift
T = 50.0     # Gate duration
σ = T/4      # Width parameter

pulse = ErfPulse([φ_max], [σ], [T/2], T)
```
"""
function ErfPulse(
    amplitudes::AbstractVector{<:Real},
    sigmas::AbstractVector{<:Real},
    centers::AbstractVector{<:Real},
    duration::Real,
)
    n = length(amplitudes)
    @assert length(sigmas) == n "sigmas must have same length as amplitudes"
    @assert length(centers) == n "centers must have same length as amplitudes"

    amps = Vector{Float64}(amplitudes)
    sigs = Vector{Float64}(sigmas)
    ctrs = Vector{Float64}(centers)
    dur = Float64(duration)

    f = t -> [amps[i] * erf(√2 * (t - ctrs[i]) / sigs[i]) for i = 1:n]
    return ErfPulse(f, amps, sigs, ctrs, dur, n)
end

"""
    ErfPulse(amplitudes, sigma, duration; center=duration/2)

Create an error function pulse with shared sigma and center across all drives.

# Arguments
- `amplitudes`: Peak amplitude for each drive
- `sigma`: Shared width parameter for all drives
- `duration`: Total pulse duration

# Keyword Arguments
- `center`: Shared center time (default: `duration/2`)
"""
function ErfPulse(
    amplitudes::AbstractVector{<:Real},
    sigma::Real,
    duration::Real;
    center::Real = duration / 2,
)
    n = length(amplitudes)
    return ErfPulse(amplitudes, fill(Float64(sigma), n), fill(Float64(center), n), duration)
end

evaluate(p::ErfPulse, t) = p.f(t)

# ============================================================================ #
# CompositePulse (combine multiple pulses)
# ============================================================================ #

"""
    CompositePulse{F<:Function} <: AbstractPulse

Composite pulse that combines multiple pulse objects by interleaving their drives.

Useful for creating pulses with different shapes for different control types,
such as Gaussian amplitude + erf phase for trapped ion gates.

# Fields
- `f::F`: Function that evaluates the composite pulse
- `pulses::Vector{<:AbstractPulse}`: Component pulses
- `drive_mapping::Vector{Vector{Int}}`: Maps pulse i, drive j to composite drive index
- `duration::Float64`: Total pulse duration (must match for all components)
- `n_drives::Int`: Total number of drives across all pulses

# Example
```julia
# Amplitude: Gaussian (2 drives for 2 ions)
Ω_pulse = GaussianPulse([Ω_max, Ω_max], σ, T)

# Phase: Error function (2 drives for 2 ions)
φ_pulse = ErfPulse([φ_max, φ_max], σ, T)

# Composite: [Ω₁, φ₁, Ω₂, φ₂] - interleaved
pulse = CompositePulse([Ω_pulse, φ_pulse], :interleave)
```
"""
struct CompositePulse{F<:Function} <: AbstractPulse
    f::F
    pulses::Vector{<:AbstractPulse}
    drive_mapping::Vector{Vector{Int}}
    duration::Float64
    n_drives::Int
end

"""
    CompositePulse(pulses::Vector{<:AbstractPulse}, mode::Symbol=:interleave)

Create a composite pulse from multiple component pulses.

# Arguments
- `pulses`: Vector of pulse objects to combine
- `mode`: How to combine the drives
  - `:interleave` - Interleave drives: [p1_d1, p2_d1, p1_d2, p2_d2, ...]
  - `:concatenate` - Concatenate drives: [p1_d1, p1_d2, ..., p2_d1, p2_d2, ...]

# Example
```julia
# For MS gate with 2 ions: [Ω₁, φ₁, Ω₂, φ₂]
Ω_pulse = GaussianPulse([Ω₁, Ω₂], σ, T)  # 2 drives
φ_pulse = ErfPulse([φ₁, φ₂], σ, T)        # 2 drives
pulse = CompositePulse([Ω_pulse, φ_pulse], :interleave)
# Result: pulse(t) = [Ω₁(t), φ₁(t), Ω₂(t), φ₂(t)]
```
"""
function CompositePulse(pulses::Vector{<:AbstractPulse}, mode::Symbol = :interleave)
    @assert !isempty(pulses) "Must provide at least one pulse"
    @assert mode in [:interleave, :concatenate] "mode must be :interleave or :concatenate"

    # Check all pulses have same duration
    dur = duration(pulses[1])
    for p in pulses
        @assert abs(duration(p) - dur) < 1e-10 "All pulses must have same duration"
    end

    # Build drive mapping based on mode
    n_pulses = length(pulses)
    n_drives_per_pulse = [n_drives(p) for p in pulses]
    total_drives = sum(n_drives_per_pulse)

    drive_mapping = [Vector{Int}() for _ = 1:n_pulses]

    if mode == :interleave
        # Check all pulses have same number of drives for interleaving
        @assert allequal(n_drives_per_pulse) "All pulses must have same n_drives for :interleave mode"

        n_per_pulse = n_drives_per_pulse[1]
        composite_idx = 1

        # Interleave: [p1_d1, p2_d1, ..., pN_d1, p1_d2, p2_d2, ..., pN_d2, ...]
        for drive_idx = 1:n_per_pulse
            for pulse_idx = 1:n_pulses
                push!(drive_mapping[pulse_idx], composite_idx)
                composite_idx += 1
            end
        end
    else  # :concatenate
        # Concatenate: [p1_d1, p1_d2, ..., p2_d1, p2_d2, ...]
        composite_idx = 1
        for pulse_idx = 1:n_pulses
            for _ = 1:n_drives_per_pulse[pulse_idx]
                push!(drive_mapping[pulse_idx], composite_idx)
                composite_idx += 1
            end
        end
    end

    # Create evaluation function (type-generic for ForwardDiff compatibility)
    f = function (t)
        # Get first pulse values to determine element type
        first_vals = pulses[1](t)
        T = eltype(first_vals)
        result = zeros(T, total_drives)

        # Fill in values from first pulse
        for (local_idx, composite_idx) in enumerate(drive_mapping[1])
            result[composite_idx] = first_vals[local_idx]
        end

        # Fill in values from remaining pulses
        for pulse_idx = 2:n_pulses
            pulse_vals = pulses[pulse_idx](t)
            for (local_idx, composite_idx) in enumerate(drive_mapping[pulse_idx])
                result[composite_idx] = pulse_vals[local_idx]
            end
        end
        return result
    end

    return CompositePulse(f, pulses, drive_mapping, dur, total_drives)
end

evaluate(p::CompositePulse, t) = p.f(t)

# ============================================================================ #
# FunctionPulse (arbitrary user-defined function)
# ============================================================================ #

"""
    FunctionPulse{F<:Function} <: AbstractPulse

Pulse defined by an arbitrary function `f(t) -> Vector{Float64}`.

Useful for testing analytic pulse shapes (e.g. sin² envelopes) with the
`rollout` / `fidelity` interface without discretizing into spline knots.

# Fields
- `f::F`: Function mapping time to control vector
- `duration::Float64`: Total pulse duration
- `n_drives::Int`: Number of control drives
- `drive_name::Symbol`: Name of the drive variable (default `:u`)

# Example
```julia
T = 1000.0
pulse = FunctionPulse(t -> [0.0, 0.0, 1.5 * sin(π*t/T)^2, 0.0], T, 4)
qtraj = MultiKetTrajectory(sys, pulse, initials, goals)
qtraj_out = rollout(qtraj)
fid = fidelity(qtraj_out)
```
"""
struct FunctionPulse{F<:Function} <: AbstractPulse
    f::F
    duration::Float64
    n_drives::Int
    drive_name::Symbol
end

function FunctionPulse(f::Function, duration::Real, n_drives::Int; drive_name::Symbol = :u)
    return FunctionPulse(f, Float64(duration), n_drives, drive_name)
end

evaluate(p::FunctionPulse, t) = p.f(t)

# Knot time accessors for analytic and composite pulses
# (defined here because GaussianPulse, ErfPulse, CompositePulse are defined above)
get_knot_times(p::FunctionPulse) = [0.0, p.duration]
get_knot_times(p::GaussianPulse) = [0.0, p.duration]
get_knot_times(p::ErfPulse) = [0.0, p.duration]
get_knot_times(p::CompositePulse) =
    sort(unique(vcat([get_knot_times(sub) for sub in p.pulses]...)))

# ============================================================================ #
# Persistence (JLD2)
# ============================================================================ #

"""
    save(filename::String, pulse::AbstractPulse)

Save `pulse` to `filename` (must end in `.jld2`) under the JLD2 key `"pulse"`.
Reload with [`load_pulse`](@ref).

To bundle a pulse with metadata (fidelity, gate name, etc.), use `JLD2.jldsave`
directly:

```julia
jldsave("my_gate.jld2"; pulse=optimized_pulse, fidelity=fidelity(qcp))
```
"""
function JLD2.save(filename::String, pulse::AbstractPulse)
    @assert split(filename, ".")[end] == "jld2"
    save(filename, "pulse", pulse)
end

"""
    load_pulse(filename::String) -> AbstractPulse

Load a pulse previously written with [`save`](@ref) (or with
`jldsave(filename; pulse=...)`) from `filename` (must end in `.jld2`).
"""
function load_pulse(filename::String)
    @assert split(filename, ".")[end] == "jld2"
    return load(filename, "pulse")
end

# ============================================================================ #
# NamedTrajectory Constructors
# ============================================================================ #

using NamedTrajectories: NamedTrajectory, get_times

"""
    ZeroOrderPulse(traj::NamedTrajectory; drive_name=:u, snap_to_knots=true)

Construct a ZeroOrderPulse from a NamedTrajectory.

# Arguments
- `traj`: NamedTrajectory with control data

# Keyword Arguments
- `drive_name`: Name of the drive component (default: `:u`)
- `snap_to_knots`: See primary constructor docstring (default: `true`).
  Round-tripping through `NamedTrajectory` produces `cumsum`-based knot
  times that differ from `range()`-based resample times by float
  roundoff, which would otherwise trigger the off-by-one at every knot.
"""
function ZeroOrderPulse(
    traj::NamedTrajectory;
    drive_name::Symbol = :u,
    snap_to_knots::Bool = true,
)
    controls = traj[drive_name]
    times = get_times(traj)
    return ZeroOrderPulse(controls, times; drive_name, snap_to_knots)
end

"""
    LinearSplinePulse(traj::NamedTrajectory; drive_name=:u)

Construct a LinearSplinePulse from a NamedTrajectory.

# Arguments
- `traj`: NamedTrajectory with control data

# Keyword Arguments
- `drive_name`: Name of the drive component (default: `:u`)
"""
function LinearSplinePulse(traj::NamedTrajectory; drive_name::Symbol = :u)
    controls = traj[drive_name]
    times = get_times(traj)
    return LinearSplinePulse(controls, times; drive_name)
end

"""
    CubicSplinePulse(traj::NamedTrajectory; drive_name=:u, derivative_name=:du)

Construct a CubicSplinePulse (Hermite) from a NamedTrajectory using both 
control values and derivatives.

# Arguments
- `traj`: NamedTrajectory with control and derivative data

# Keyword Arguments
- `drive_name`: Name of the drive component (default: `:u`)
- `derivative_name`: Name of the derivative component (default: `:du`)
"""
function CubicSplinePulse(
    traj::NamedTrajectory;
    drive_name::Symbol = :u,
    derivative_name::Symbol = :du,
)
    controls = traj[drive_name]
    derivatives = traj[derivative_name]
    times = get_times(traj)
    return CubicSplinePulse(controls, derivatives, times; drive_name)
end

# ============================================================================ #
# _rebuild_like: NamedTuple-keyed warm-restart reconstruction
# ============================================================================ #

"""
    _rebuild_like(p::AbstractPulse, vars::NamedTuple) -> AbstractPulse

Reconstruct a pulse of the same concrete type with new decision-variable values.
Used by the optimizer's warm-restart utility. The `vars` NamedTuple keys must
match the decision-variable component symbols owned by the pulse type:

- `LinearSplinePulse`: `(u = new_values,)`
- `CubicSplinePulse`:  `(u = new_values, du = new_tangents)`
- `ZeroOrderPulse`:    `(u = new_values,)`
- `BSplinePulse`:      `(u = new_control_points,)` (added in BSplinePulse PR)

Additional pulse-type implementations may be added later under the same
NamedTuple contract. The pulse's basis / time grid / duration / boundary
metadata is preserved.
"""
function _rebuild_like end

function _rebuild_like(p::LinearSplinePulse, vars::NamedTuple)
    haskey(vars, :u) || throw(ArgumentError(
        "LinearSplinePulse._rebuild_like requires :u in NamedTuple; got keys $(keys(vars))"))
    return LinearSplinePulse(
        vars.u,
        get_knot_times(p);
        drive_name = p.drive_name,
        initial_value = p.initial_value,
        final_value = p.final_value,
    )
end

function _rebuild_like(p::CubicSplinePulse, vars::NamedTuple)
    (haskey(vars, :u) && haskey(vars, :du)) || throw(ArgumentError(
        "CubicSplinePulse._rebuild_like requires :u and :du in NamedTuple; got keys $(keys(vars))"))
    return CubicSplinePulse(
        vars.u,
        vars.du,
        get_knot_times(p);
        drive_name = p.drive_name,
        initial_value = p.initial_value,
        final_value = p.final_value,
    )
end

function _rebuild_like(p::ZeroOrderPulse, vars::NamedTuple)
    haskey(vars, :u) || throw(ArgumentError(
        "ZeroOrderPulse._rebuild_like requires :u in NamedTuple; got keys $(keys(vars))"))
    return ZeroOrderPulse(
        vars.u,
        get_knot_times(p);
        drive_name = p.drive_name,
        initial_value = p.initial_value,
        final_value = p.final_value,
        snap_to_knots = p.snap_to_knots,
    )
end

function _rebuild_like(p::BSplinePulse, vars::NamedTuple)
    haskey(vars, :u) || throw(ArgumentError(
        "BSplinePulse._rebuild_like requires :u in NamedTuple; got keys $(keys(vars))"))
    new_cps = vars.u
    size(new_cps) == size(p.control_points) || throw(DimensionMismatch(
        "new control_points size $(size(new_cps)) != original $(size(p.control_points))"))
    Tnew = promote_type(eltype(new_cps), Float64)
    return BSplinePulse{Tnew,typeof(p.basis)}(
        p.basis, Matrix{Tnew}(new_cps), p.duration, p.n_drives,
        p.drive_name, p.initial_value, p.final_value,
    )
end

# ============================================================================ #
# Tests
# ============================================================================ #

@testitem "ZeroOrderPulse" begin

    # Create simple pulse
    controls = [0.0 1.0 0.5 0.0; 0.0 -1.0 -0.5 0.0]
    times = [0.0, 0.25, 0.5, 1.0]
    pulse = ZeroOrderPulse(controls, times)

    @test duration(pulse) == 1.0
    @test n_drives(pulse) == 2
    @test drive_name(pulse) == :u  # Default

    # Test evaluation (zero-order hold)
    @test pulse(0.0) ≈ [0.0, 0.0]
    @test pulse(0.1) ≈ [0.0, 0.0]  # Before first transition
    @test pulse(0.3) ≈ [1.0, -1.0] # After first transition
    @test pulse(1.0) ≈ [0.0, 0.0]

    # Test sampling
    n_samples = 5
    sampled, ts = sample(pulse, n_samples)
    @test size(sampled) == (2, 5)
    @test length(ts) == 5

    # Test custom drive_name
    pulse_custom = ZeroOrderPulse(controls, times; drive_name = :Ω)
    @test drive_name(pulse_custom) == :Ω
end

@testitem "LinearSplinePulse" begin

    # Create simple pulse
    controls = [0.0 1.0 0.0; 0.0 -1.0 0.0]
    times = [0.0, 0.5, 1.0]
    pulse = LinearSplinePulse(controls, times)

    @test duration(pulse) == 1.0
    @test n_drives(pulse) == 2
    @test drive_name(pulse) == :u  # Default

    # Test linear interpolation
    @test pulse(0.0) ≈ [0.0, 0.0]
    @test pulse(0.25) ≈ [0.5, -0.5]  # Midpoint of first segment
    @test pulse(0.5) ≈ [1.0, -1.0]
    @test pulse(0.75) ≈ [0.5, -0.5]  # Midpoint of second segment
    @test pulse(1.0) ≈ [0.0, 0.0]

    # Test sampling
    n_samples = 5
    sampled, ts = sample(pulse, n_samples)
    @test size(sampled) == (2, 5)

    # Test custom drive_name
    pulse_custom = LinearSplinePulse(controls, times; drive_name = :amplitude)
    @test drive_name(pulse_custom) == :amplitude
end

@testitem "CubicSplinePulse" begin

    # Create pulse with explicit derivatives (Hermite spline)
    controls = [0.0 0.5 1.0 0.5 0.0; 0.0 -0.5 -1.0 -0.5 0.0]
    derivatives = [2.0 2.0 0.0 -2.0 -2.0; -2.0 -2.0 0.0 2.0 2.0]  # Symmetric derivatives
    times = [0.0, 0.25, 0.5, 0.75, 1.0]
    pulse = CubicSplinePulse(controls, derivatives, times)

    @test duration(pulse) == 1.0
    @test n_drives(pulse) == 2
    @test drive_name(pulse) == :u  # Default

    # Test that endpoints are correct
    @test pulse(0.0) ≈ [0.0, 0.0]
    @test pulse(1.0) ≈ [0.0, 0.0]

    # Test that it passes through sample points
    @test pulse(0.5) ≈ [1.0, -1.0]

    # Cubic spline should be smooth (no discontinuities)
    @test pulse(0.4) isa Vector{Float64}
    @test pulse(0.6) isa Vector{Float64}

    # Test custom drive_name
    pulse_custom = CubicSplinePulse(controls, derivatives, times; drive_name = :a)
    @test drive_name(pulse_custom) == :a

    # Test zero-derivative constructor
    pulse_zero_deriv = CubicSplinePulse(controls, times)
    @test pulse_zero_deriv(0.5) ≈ [1.0, -1.0]
end

@testitem "sample interface" begin
    using .Pulses: sample, duration, n_drives

    controls = [0.0 1.0 0.0; 0.0 -1.0 0.0]
    times_ctrl = [0.0, 0.5, 1.0]

    for P in (ZeroOrderPulse, LinearSplinePulse, CubicSplinePulse)
        pulse = P(controls, times_ctrl)

        # sample(pulse, n_samples::Int) -> (values, times)
        n = 7
        vals, ts = sample(pulse, n)
        @test vals isa AbstractMatrix
        @test size(vals) == (n_drives(pulse), n)
        @test length(ts) == n
        @test first(ts) == 0.0
        @test last(ts) ≈ duration(pulse)

        # sample(pulse, times::AbstractVector) -> values only
        query_ts = collect(LinRange(0.0, duration(pulse), 4))
        vals_at_times = sample(pulse, query_ts)
        @test vals_at_times isa AbstractMatrix
        @test size(vals_at_times) == (n_drives(pulse), length(query_ts))

        # Locked-in interface: no keyword-argument form, no zero-arg default
        @test_throws MethodError sample(pulse; n_samples = 5)
        @test_throws MethodError sample(pulse)
    end

    # sample(bounds, n_samples) -> random matrix within bounds
    bounds = [(-1.0, 1.0), (-0.5, 0.5)]
    rand_ctrls = sample(bounds, 10)
    @test size(rand_ctrls) == (length(bounds), 10)
    for (i, (lo, hi)) in enumerate(bounds)
        @test all(lo .<= rand_ctrls[i, :] .<= hi)
    end
end

@testitem "DataInterpolations derivative at boundaries" begin
    using .Pulses: derivative, sample

    n_samples = 100
    bounds = [(-1e-3, 1e-3), (-1e3, 1e3)]
    inits = sample(bounds, n_samples)
    times = LinRange(0, 10.0, n_samples)

    for P in [ZeroOrderPulse, LinearSplinePulse]

        pulse = P(inits, times)

        # The pulse itself should sample cleanly on the closed interval.
        vals = sample(pulse, times)
        @test size(vals) == (n_drives(pulse), length(times))

        # The pulse-level derivative API should also work on the closed interval,
        # including both boundaries.
        derivs = [derivative(pulse, t) for t in times]
        @test length(derivs) == length(times)
        @test all(d -> length(d) == n_drives(pulse), derivs)

        # Explicitly check the endpoints, since these are where ForwardDiff on the
        # raw interpolation used to fail due to extrapolation.
        @test derivative(pulse, first(times)) isa AbstractVector
        @test derivative(pulse, last(times)) isa AbstractVector
        @test length(derivative(pulse, first(times))) == n_drives(pulse)
        @test length(derivative(pulse, last(times))) == n_drives(pulse)
    end

    # CubicHermiteSpline known failure mode
    pulse = CubicSplinePulse(inits, times)
    @test_broken derivative(pulse, last(times)) isa AbstractVector
end

@testitem "GaussianPulse" begin

    # Test with per-drive parameters
    amplitudes = [1.0, 2.0]
    sigmas = [0.1, 0.2]
    centers = [0.5, 0.6]
    dur = 1.0

    pulse = GaussianPulse(amplitudes, sigmas, centers, dur)

    @test duration(pulse) == 1.0
    @test n_drives(pulse) == 2
    @test pulse.amplitudes == amplitudes
    @test pulse.sigmas == sigmas
    @test pulse.centers == centers

    # Test that peak is at center
    @test pulse(0.5)[1] ≈ 1.0  # First drive peaks at t=0.5
    @test pulse(0.6)[2] ≈ 2.0  # Second drive peaks at t=0.6

    # Test Gaussian decay
    @test pulse(0.0)[1] < 0.01  # Should be near zero far from center
    @test pulse(1.0)[1] < 0.01

    # Test convenience constructor with shared parameters
    pulse2 = GaussianPulse([1.0, 2.0], 0.1, 1.0)
    @test duration(pulse2) == 1.0
    @test n_drives(pulse2) == 2
    @test pulse2.centers == [0.5, 0.5]  # Default center is duration/2
    @test pulse2.sigmas == [0.1, 0.1]

    # Test with custom center
    pulse3 = GaussianPulse([1.0], 0.1, 1.0; center = 0.3)
    @test pulse3.centers == [0.3]
    @test pulse3(0.3)[1] ≈ 1.0  # Peak at specified center
end

@testitem "ErfPulse" begin
    using SpecialFunctions: erf

    # Test with per-drive parameters
    amplitudes = [1.0, 2.0]
    sigmas = [0.2, 0.3]
    centers = [0.5, 0.6]
    dur = 1.0

    pulse = ErfPulse(amplitudes, sigmas, centers, dur)

    @test duration(pulse) == 1.0
    @test n_drives(pulse) == 2
    @test pulse.amplitudes == amplitudes
    @test pulse.sigmas == sigmas
    @test pulse.centers == centers

    # Test that center is inflection point (erf(0) = 0)
    @test pulse(0.5)[1] ≈ 0.0  # First drive at center
    @test pulse(0.6)[2] ≈ 0.0  # Second drive at center

    # Test asymptotic behavior
    @test pulse(0.0)[1] < -0.9 * amplitudes[1]  # Near -1 far before center
    @test pulse(1.0)[1] > 0.9 * amplitudes[1]   # Near +1 far after center

    # Test monotonic increase
    @test pulse(0.4)[1] < pulse(0.5)[1]
    @test pulse(0.5)[1] < pulse(0.6)[1]

    # Test convenience constructor with shared parameters
    pulse2 = ErfPulse([1.0, 2.0], 0.2, 1.0)
    @test duration(pulse2) == 1.0
    @test n_drives(pulse2) == 2
    @test pulse2.centers == [0.5, 0.5]  # Default center is duration/2
    @test pulse2.sigmas == [0.2, 0.2]

    # Test with custom center
    pulse3 = ErfPulse([1.0], 0.2, 1.0; center = 0.3)
    @test pulse3.centers == [0.3]
    @test pulse3(0.3)[1] ≈ 0.0  # Zero at center

    # Test phase compensation use case
    φ_max = π / 4
    T = 50.0
    σ = T / 4
    phase_pulse = ErfPulse([φ_max], [σ], [T / 2], T)
    @test phase_pulse(T / 2)[1] ≈ 0.0
    @test phase_pulse(0.0)[1] < -0.8 * φ_max
    @test phase_pulse(T)[1] > 0.8 * φ_max
end

@testitem "CompositePulse" begin

    # Create component pulses
    T = 10.0
    amp_pulse = GaussianPulse([1.0, 2.0], 1.0, T)     # 2 drives
    phase_pulse = ErfPulse([0.5, 0.8], 1.0, T)        # 2 drives

    # Test interleave mode (default)
    composite = CompositePulse([amp_pulse, phase_pulse], :interleave)

    @test duration(composite) == T
    @test n_drives(composite) == 4  # 2 + 2

    # Test evaluation at center (t=5.0)
    # Gaussian peaks at center, erf is zero at center
    vals = composite(5.0)
    @test length(vals) == 4

    # Interleaved order: [amp_d1, phase_d1, amp_d2, phase_d2]
    @test vals[1] ≈ amp_pulse(5.0)[1]      # Amplitude drive 1
    @test vals[2] ≈ phase_pulse(5.0)[1]    # Phase drive 1
    @test vals[3] ≈ amp_pulse(5.0)[2]      # Amplitude drive 2
    @test vals[4] ≈ phase_pulse(5.0)[2]    # Phase drive 2

    @test vals[1] ≈ 1.0  # Gaussian peak
    @test vals[2] ≈ 0.0  # Erf zero
    @test vals[3] ≈ 2.0  # Gaussian peak
    @test vals[4] ≈ 0.0  # Erf zero

    # Test concatenate mode
    composite2 = CompositePulse([amp_pulse, phase_pulse], :concatenate)
    @test n_drives(composite2) == 4

    vals2 = composite2(5.0)
    # Concatenated order: [amp_d1, amp_d2, phase_d1, phase_d2]
    @test vals2[1] ≈ amp_pulse(5.0)[1]
    @test vals2[2] ≈ amp_pulse(5.0)[2]
    @test vals2[3] ≈ phase_pulse(5.0)[1]
    @test vals2[4] ≈ phase_pulse(5.0)[2]

    # Test with different n_drives (concatenate only)
    pulse_3drive = GaussianPulse([1.0, 2.0, 3.0], 1.0, T)
    pulse_2drive = ErfPulse([0.5, 0.8], 1.0, T)
    composite3 = CompositePulse([pulse_3drive, pulse_2drive], :concatenate)

    @test n_drives(composite3) == 5  # 3 + 2
    vals3 = composite3(5.0)
    @test length(vals3) == 5
    @test vals3[1:3] ≈ pulse_3drive(5.0)
    @test vals3[4:5] ≈ pulse_2drive(5.0)

    # Test that interleave with different n_drives fails
    @test_throws AssertionError CompositePulse([pulse_3drive, pulse_2drive], :interleave)
end

@testitem "Pulse conversion to splines" begin

    # Create analytic pulses
    T = 10.0
    gaussian = GaussianPulse([1.0, 2.0], 1.0, T)
    erf_pulse = ErfPulse([0.5, 0.8], 1.0, T)

    # Test LinearSplinePulse conversion
    linear = LinearSplinePulse(gaussian, 20)
    @test linear isa LinearSplinePulse
    @test duration(linear) == T
    @test n_drives(linear) == 2

    # Check that sampled values match original at knot points
    times = range(0, T, length = 20)
    for t in times
        @test linear(t) ≈ gaussian(t) atol = 1e-10
    end

    # Test with custom times
    custom_times = [0.0, 2.5, 5.0, 7.5, 10.0]
    linear2 = LinearSplinePulse(gaussian, custom_times)
    @test duration(linear2) == T
    for t in custom_times
        @test linear2(t) ≈ gaussian(t) atol = 1e-10
    end

    # Test CubicSplinePulse conversion with derivatives
    cubic = CubicSplinePulse(gaussian, 20)
    @test cubic isa CubicSplinePulse
    @test duration(cubic) == T
    @test n_drives(cubic) == 2

    # Cubic spline should be very close to original at sample points
    for t in times
        @test cubic(t) ≈ gaussian(t) atol = 1e-10
    end

    # Test with CompositePulse
    composite = CompositePulse([gaussian, erf_pulse], :interleave)
    composite_linear = LinearSplinePulse(composite, 30)
    @test n_drives(composite_linear) == 4  # 2 + 2 interleaved
    @test duration(composite_linear) == T

    composite_cubic = CubicSplinePulse(composite, 30)
    @test n_drives(composite_cubic) == 4
    @test duration(composite_cubic) == T

    # Check interleaving is preserved
    test_times = [0.0, 5.0, 10.0]
    for t in test_times
        orig = composite(t)
        linear_val = composite_linear(t)
        cubic_val = composite_cubic(t)
        @test linear_val ≈ orig atol = 1e-1
        @test cubic_val ≈ orig atol = 1e-1
    end
end

@testitem "Pulse callability" begin

    # All pulse types should be callable
    controls = [0.0 1.0 0.0]
    times = [0.0, 0.5, 1.0]

    zop = ZeroOrderPulse(controls, times)
    lp = LinearSplinePulse(controls, times)
    gp = GaussianPulse([1.0], 0.1, 1.0)

    # Test callable interface
    for pulse in [zop, lp, gp]
        @test pulse(0.5) isa Vector{Float64}
        @test length(pulse(0.5)) == 1
    end
end

@testitem "Pulse from NamedTrajectory" begin
    using NamedTrajectories: NamedTrajectory

    # Create a simple NamedTrajectory with controls and derivatives
    T = 5
    times = collect(range(0, 1, length = T))
    Δt = diff(times)
    Δt = [Δt; Δt[end]]

    u = [sin(π * t) for t in times]
    du = [π * cos(π * t) for t in times]

    traj = NamedTrajectory(
        (u = reshape(u, 1, :), du = reshape(du, 1, :), Δt = Δt);
        timestep = :Δt,
        controls = (:u,),
    )

    # Test ZeroOrderPulse from NamedTrajectory
    zop = ZeroOrderPulse(traj)
    @test duration(zop) ≈ 1.0
    @test n_drives(zop) == 1
    @test zop(0.0) ≈ [0.0] atol = 1e-10

    # Test LinearSplinePulse from NamedTrajectory
    lp = LinearSplinePulse(traj)
    @test duration(lp) ≈ 1.0
    @test n_drives(lp) == 1
    @test lp(0.5) ≈ [1.0] atol = 1e-10

    # Test CubicSplinePulse from NamedTrajectory (uses both u and du)
    csp = CubicSplinePulse(traj)
    @test duration(csp) ≈ 1.0
    @test n_drives(csp) == 1
    @test csp(0.5) ≈ [1.0] atol = 1e-10
    @test csp(0.0) ≈ [0.0] atol = 1e-10
end

@testitem "Pulse Boundary Conditions" begin
    using LinearAlgebra

    # Test ZeroOrderPulse boundary conditions
    n_drives = 2
    T = 1.0
    N = 10
    controls = rand(n_drives, N)
    times = range(0, T, length = N)

    # Default: zeros
    pulse_default = ZeroOrderPulse(controls, times)
    @test pulse_default.initial_value == zeros(n_drives)
    @test pulse_default.final_value == zeros(n_drives)

    # Custom values
    init_val = [1.0, 2.0]
    final_val = [3.0, 4.0]
    pulse_custom =
        ZeroOrderPulse(controls, times; initial_value = init_val, final_value = final_val)
    @test pulse_custom.initial_value == init_val
    @test pulse_custom.final_value == final_val

    # Test LinearSplinePulse boundary conditions
    pulse_linear_default = LinearSplinePulse(controls, times)
    @test pulse_linear_default.initial_value == zeros(n_drives)
    @test pulse_linear_default.final_value == zeros(n_drives)

    pulse_linear_custom = LinearSplinePulse(
        controls,
        times;
        initial_value = init_val,
        final_value = final_val,
    )
    @test pulse_linear_custom.initial_value == init_val
    @test pulse_linear_custom.final_value == final_val

    # Test CubicSplinePulse boundary conditions
    derivatives = rand(n_drives, N)

    pulse_cubic_default = CubicSplinePulse(controls, derivatives, times)
    @test pulse_cubic_default.initial_value == zeros(n_drives)
    @test pulse_cubic_default.final_value == zeros(n_drives)

    pulse_cubic_custom = CubicSplinePulse(
        controls,
        derivatives,
        times;
        initial_value = init_val,
        final_value = final_val,
    )
    @test pulse_cubic_custom.initial_value == init_val
    @test pulse_cubic_custom.final_value == final_val

    # Test CubicSplinePulse constructor without derivatives
    pulse_cubic_no_deriv =
        CubicSplinePulse(controls, times; initial_value = init_val, final_value = final_val)
    @test pulse_cubic_no_deriv.initial_value == init_val
    @test pulse_cubic_no_deriv.final_value == final_val
end

@testitem "NamedTrajectory Boundary Condition Extraction" begin
    using LinearAlgebra
    using NamedTrajectories

    # Create a simple quantum system
    H_drift = zeros(ComplexF64, 2, 2)
    H_drives = [ComplexF64[0 1; 1 0]]
    sys = QuantumSystem(H_drift, H_drives, [(-1.0, 1.0)])

    # Create a pulse with custom boundary conditions
    n_drives = 1
    T = 1.0
    N = 10
    controls = rand(n_drives, N)
    derivatives = rand(n_drives, N)
    times = range(0, T, length = N)

    init_val = [0.5]
    final_val = [-0.5]
    pulse = CubicSplinePulse(
        controls,
        derivatives,
        times;
        initial_value = init_val,
        final_value = final_val,
    )

    # Create a quantum trajectory and convert to NamedTrajectory
    U_goal = exp(-im * 0.5 * H_drives[1])
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)
    traj = NamedTrajectory(qtraj, N)

    # Verify boundary conditions are extracted
    @test haskey(traj.initial, :u)
    @test haskey(traj.final, :u)
    @test haskey(traj.initial, :du)
    @test haskey(traj.final, :du)
    @test traj.initial[:u] == init_val
    @test traj.final[:u] == final_val
    @test traj.initial[:du] == zeros(n_drives)  # Derivatives constrained to zero
    @test traj.final[:du] == zeros(n_drives)

    # Test with zero boundaries (default)
    pulse_zero = CubicSplinePulse(controls, derivatives, times)
    qtraj_zero = UnitaryTrajectory(sys, pulse_zero, U_goal)
    traj_zero = NamedTrajectory(qtraj_zero, N)

    @test traj_zero.initial[:u] == zeros(n_drives)
    @test traj_zero.final[:u] == zeros(n_drives)
    @test traj_zero.initial[:du] == zeros(n_drives)
    @test traj_zero.final[:du] == zeros(n_drives)
end

@testitem "ZeroOrderPulse snap_to_knots" begin
    using NamedTrajectories: NamedTrajectory, get_times

    N = 51
    T = 10.0
    controls = Float64.(reshape(1:N, 1, N))

    # cumsum-based times (what get_times(traj) produces for variable-Δt)
    Δt = T / (N - 1)
    times_cumsum = cumsum([0.0; fill(Δt, N - 1)])
    times_range = collect(range(0.0, T, length = N))

    # Confirm time vectors actually differ
    @test times_cumsum != times_range
    @test count(!iszero, times_range .- times_cumsum) > 0

    # Default is snap_to_knots=true: pointwise evaluation at range times is
    # correct even though the stored knots are cumsum-based.
    pulse_default = ZeroOrderPulse(controls, times_cumsum)
    @test pulse_default.snap_to_knots == true
    stored = Matrix(pulse_default.controls.u)
    sampled_default = hcat([pulse_default(t) for t in times_range]...)
    @test sampled_default == stored

    # Explicit snap_to_knots=true matches the default.
    pulse = ZeroOrderPulse(controls, times_cumsum; snap_to_knots = true)
    @test pulse.snap_to_knots == true
    sampled = hcat([pulse(t) for t in times_range]...)
    @test sampled == stored

    # Batch sample() correct (is_native bypass, independent of snap_to_knots)
    sampled_batch, _ = sample(pulse, N)
    @test sampled_batch == stored

    # Full round-trip: NamedTrajectory → ZeroOrderPulse → resample at
    # range-based times. This is the path that triggered the bug in the
    # cryoCMOS spin demo (knot times are reconstructed via get_times(traj),
    # which uses cumsum).
    Δt_vec = fill(Δt, N)
    data = (u = controls, Δt = reshape(Δt_vec, 1, N))
    traj = NamedTrajectory(data; timestep = :Δt, controls = (:Δt, :u))
    pulse_rt = ZeroOrderPulse(traj; drive_name = :u)
    @test pulse_rt.snap_to_knots == true
    rt_range = collect(range(0.0, duration(pulse_rt), length = N))
    rt_sampled = hcat([pulse_rt(t) for t in rt_range]...)
    @test rt_sampled == Matrix(pulse_rt.controls.u)

    # snap_to_knots=false is the legacy left-continuous path, retained as an
    # explicit opt-out. It still reproduces the off-by-one (every mismatch
    # is exactly u[k-1]), which we lock in here so that anyone who flips this
    # flag back on knows what they're getting.
    pulse_raw = ZeroOrderPulse(controls, times_cumsum; snap_to_knots = false)
    @test pulse_raw.snap_to_knots == false
    sampled_raw = hcat([pulse_raw(t) for t in times_range]...)
    n_mismatches = count(k -> sampled_raw[1, k] != stored[1, k], 1:N)
    @test n_mismatches > 0
    for k = 1:N
        if sampled_raw[1, k] != stored[1, k]
            @test k > 1
            @test sampled_raw[1, k] == stored[1, k-1]
        end
    end

    # Batch sample() with is_native bypass already returned correct values
    # regardless of the snap flag — verify that holds for the legacy path too.
    sampled_legacy_batch, _ = sample(pulse_raw, N)
    @test sampled_legacy_batch == stored
end

@testitem "ZeroOrderPulse → NamedTrajectory round-trip preserves knot values" begin
    # Regression test for the off-by-one: when a NamedTrajectory is built via
    # the trajectory→pulse→trajectory path with N matching the pulse, every
    # knot's control should equal the pulse's stored knot value.
    using NamedTrajectories: NamedTrajectory, get_times

    N = 51
    T = 10.0
    # Sharply varying control: u[k] = k makes any shift visible immediately.
    controls = Float64.(reshape(1:N, 1, N))
    Δt = T / (N - 1)

    # The path that triggered the bug in the wild: a saved pulse whose
    # stored knot times come from cumsum is later sampled at range-based
    # times by NamedTrajectory(qtraj, N).
    times_cumsum = cumsum([0.0; fill(Δt, N - 1)])
    pulse = ZeroOrderPulse(controls, times_cumsum)
    stored = Matrix(pulse.controls.u)

    # Pointwise evaluation at range times: every knot's value is the stored
    # knot value, not u[k-1].
    times_range = collect(range(0.0, T, length = N))
    sampled = hcat([pulse(t) for t in times_range]...)
    @test sampled == stored

    # Physical-scale alternating ±600 rad/μs (the cryoCMOS J-channel scale).
    # Without snapping this would sign-flip controls at affected knots.
    phys = reshape([k % 2 == 0 ? 600.0 : -600.0 for k = 1:N], 1, N)
    phys_pulse = ZeroOrderPulse(phys, times_cumsum)
    phys_sampled = hcat([phys_pulse(t) for t in times_range]...)
    @test maximum(abs, phys_sampled .- Matrix(phys_pulse.controls.u)) == 0.0
end

@testitem "save / load_pulse round-trip" begin
    using JLD2

    times = collect(range(0.0, 1.0, length = 8))
    controls = randn(2, 8)

    pulses = (
        ZeroOrderPulse(controls, times),
        LinearSplinePulse(controls, times),
        CubicSplinePulse(controls, randn(2, 8), times),
        GaussianPulse([1.0, 0.5], [0.1, 0.2], [0.3, 0.6], 1.0),
    )

    mktempdir() do dir
        for pulse in pulses
            path = joinpath(dir, "p.jld2")
            save(path, pulse)
            @test isfile(path)
            loaded = load_pulse(path)
            @test typeof(loaded) == typeof(pulse)
            @test duration(loaded) == duration(pulse)
            @test n_drives(loaded) == n_drives(pulse)
            for t in range(0.0, duration(pulse), length = 5)
                @test loaded(t) ≈ pulse(t)
            end
        end
    end
end

@testitem "save / load_pulse filename extension assertion" begin
    pulse = ZeroOrderPulse(randn(1, 4), collect(range(0.0, 1.0, length = 4)))
    @test_throws AssertionError save("bad_extension.txt", pulse)
    @test_throws AssertionError load_pulse("missing.txt")
end

@testitem "jldsave bundle with metadata round-trips pulse" begin
    using JLD2

    times = collect(range(0.0, 1.0, length = 6))
    pulse = ZeroOrderPulse(randn(2, 6), times)

    mktempdir() do dir
        path = joinpath(dir, "bundle.jld2")
        jldsave(path; pulse = pulse, fidelity = 0.9876, gate = "X")
        data = load(path)
        @test data["fidelity"] == 0.9876
        @test data["gate"] == "X"
        @test data["pulse"](0.5) ≈ pulse(0.5)
    end
end

@testitem "_rebuild_like NamedTuple dispatch for existing pulse types" begin
    using Piccolo

    N = 8
    times = collect(range(0.0, 1.0, length = N))
    u0 = randn(2, N)

    # LinearSplinePulse
    p_lin = LinearSplinePulse(u0, times)
    p_lin_new = Piccolo.Pulses._rebuild_like(p_lin, (u = u0 .+ 0.1,))
    @test p_lin_new isa LinearSplinePulse
    @test get_knot_values(p_lin_new) ≈ u0 .+ 0.1
    @test duration(p_lin_new) == duration(p_lin)

    # CubicSplinePulse
    du0 = randn(2, N)
    p_cub = CubicSplinePulse(u0, du0, times)
    p_cub_new = Piccolo.Pulses._rebuild_like(p_cub, (u = u0 .+ 0.1, du = du0 .+ 0.01))
    @test p_cub_new isa CubicSplinePulse
    @test get_knot_values(p_cub_new) ≈ u0 .+ 0.1
    @test get_knot_derivatives(p_cub_new) ≈ du0 .+ 0.01

    # ZeroOrderPulse
    p_zoh = ZeroOrderPulse(u0, times)
    p_zoh_new = Piccolo.Pulses._rebuild_like(p_zoh, (u = u0 .+ 0.1,))
    @test p_zoh_new isa ZeroOrderPulse

    # Missing-key error contracts
    @test_throws ArgumentError Piccolo.Pulses._rebuild_like(p_lin, (v = u0,))
    @test_throws ArgumentError Piccolo.Pulses._rebuild_like(p_cub, (u = u0,))
end

@testitem "BSplinePulse: primary constructor + struct invariants" begin
    using Piccolo
    using Random
    Random.seed!(20260525)

    nbf = 22  # for d=3, with arbitrary >= d+1
    times = collect(range(0.0, 1.0, length = 20))
    C = randn(2, nbf)
    iv = randn(2)
    C[:, 1] = iv
    fv = randn(2)
    C[:, end] = fv

    p = BSplinePulse(C, times; degree = 3,
        initial_value = iv, final_value = fv, drive_name = :u)
    @test p isa BSplinePulse
    @test n_drives(p) == 2
    @test duration(p) == 1.0
    @test p.initial_value ≈ iv
    @test p.final_value ≈ fv
    @test p.basis.order == 4
    @test num_basis_functions(p.basis) == nbf

    # Inconsistent clamp should error
    C_bad = copy(C)
    C_bad[:, 1] .+= 0.1
    @test_throws ArgumentError BSplinePulse(C_bad, times; degree = 3,
        initial_value = iv, final_value = fv)

    # Too-few CPs should error
    @test_throws ArgumentError BSplinePulse(randn(2, 3), times; degree = 3)

    # :free boundary semantics
    p_free = BSplinePulse(C, times; degree = 3, initial_value = :free, final_value = :free)
    @test p_free.initial_value === :free
    @test p_free.final_value === :free
end

@testitem "BSplinePulse: evaluate / derivative / accessors" begin
    using Piccolo

    using Random
    Random.seed!(20260525)

    nbf, T_dur = 18, 2.5
    times = collect(range(0.0, T_dur, length = 16))
    C = randn(2, nbf)
    p = BSplinePulse(C, times; degree = 3, initial_value = :free, final_value = :free)

    # Endpoint exactness via call form
    @test p(0.0) ≈ C[:, 1] atol = 1e-12
    @test p(T_dur) ≈ C[:, end] atol = 1e-12

    # Pointwise consistency with basis
    for _ in 1:30
        t = T_dur * rand()
        τ = t / T_dur
        @test p(t) ≈ evaluate_curve(p.basis, p.control_points, τ) atol = 1e-12
    end

    # Derivative consistency: du/dt = (1/T) du/dτ via FD
    h = 1e-6
    for _ in 1:10
        t = 0.1 * T_dur + 0.8 * T_dur * rand()
        analytic = Piccolo.Pulses.derivative(p, t)
        fd = (p(t + h) .- p(t - h)) ./ (2h)
        @test analytic ≈ fd atol = 1e-4
    end

    # Accessors
    @test get_control_points(p) === p.control_points
    @test get_knots(p) === p.basis.knots
    @test drive_name(p) == :u
    @test duration(p) == T_dur
    @test n_drives(p) == 2
end

@testitem "BSplinePulse: warm-start from AbstractPulse via BSplineApprox" begin
    using Piccolo
    using Random
    Random.seed!(20260525)

    N = 20
    times = collect(range(0.0, 1.0, length = N))
    u_vals = randn(2, N)
    du_vals = randn(2, N) .* 0.5
    p_cub = CubicSplinePulse(u_vals, du_vals, times)

    # Fit with 1.5x oversampling
    n_cp = Int(round(1.5 * N))
    p_bsp = BSplinePulse(p_cub, n_cp; degree = 3)
    @test p_bsp isa BSplinePulse
    @test num_basis_functions(p_bsp.basis) == n_cp

    # Fit error: sample original at 50 points, compare to bspline approx
    diffs = Float64[]
    for t in range(0.05, 0.95, length = 50)
        diff = p_cub(t) .- p_bsp(t)
        push!(diffs, sqrt(sum(abs2, diff)))
    end
    @test maximum(diffs) < 5e-1  # generous: BSplineApprox is L2-fit, not interp
end

@testitem "BSplinePulse: JLD2 round-trip" begin
    using Piccolo
    using JLD2
    using Random
    Random.seed!(20260525)

    nbf = 16
    times = collect(range(0.0, 1.0, length = 14))
    C = randn(2, nbf)
    p = BSplinePulse(C, times; degree = 3, initial_value = :free, final_value = :free)

    mktempdir() do dir
        path = joinpath(dir, "bspline_pulse.jld2")
        JLD2.save(path, p)
        p_loaded = load_pulse(path)
        @test p_loaded isa BSplinePulse
        @test get_control_points(p_loaded) ≈ get_control_points(p)
        @test p_loaded.basis.order == p.basis.order
        @test p_loaded.basis.knots ≈ p.basis.knots
        @test p_loaded.duration == p.duration
        # Evaluate at a few points to verify functional equivalence after load
        for t in (0.0, 0.25, 0.5, 1.0)
            @test p_loaded(t) ≈ p(t) atol = 1e-12
        end
    end
end

@testitem "BSplinePulse: _rebuild_like" begin
    using Piccolo
    using Random
    Random.seed!(20260525)

    nbf = 16
    times = collect(range(0.0, 1.0, length = 14))
    C = randn(2, nbf)
    p = BSplinePulse(C, times; degree = 3, initial_value = :free, final_value = :free)

    # _rebuild_like with :u key
    C_new = C .+ 0.1
    p_new = Piccolo.Pulses._rebuild_like(p, (u = C_new,))
    @test p_new isa BSplinePulse
    @test get_control_points(p_new) ≈ C_new
    @test duration(p_new) == duration(p)
    @test p_new.basis.order == p.basis.order
    @test p_new.basis === p.basis  # basis is shared (immutable structure)
    @test_throws ArgumentError Piccolo.Pulses._rebuild_like(p, (du = C_new,))
end

@testitem "BSplinePulse ↔ CubicSplinePulse round-trip identity (degree 3)" begin
    using Piccolo
    using Random
    Random.seed!(20260525)

    # Construct a known B-spline with random CPs
    N = 20
    n_cp = N + 2  # n_cp = N + d - 1 for d=3
    T_dur = 1.0
    times = collect(range(0.0, T_dur, length = N))
    C = randn(2, n_cp)
    p_bsp = BSplinePulse(C, times; degree = 3, initial_value = :free, final_value = :free)

    # Sample at knot times (the Hermite knots — the N data-point times)
    u_vals = Matrix{Float64}(undef, 2, N)
    du_vals = Matrix{Float64}(undef, 2, N)
    for (k, t) in enumerate(times)
        u_vals[:, k] = p_bsp(t)
        du_vals[:, k] = Piccolo.Pulses.derivative(p_bsp, t)
    end

    # Build CubicSplinePulse from those samples
    p_cub = CubicSplinePulse(u_vals, du_vals, times)

    # Round-trip identity: must agree pointwise on each segment (same cubic poly)
    for _ in 1:300
        t = T_dur * rand()
        v_bsp = p_bsp(t)
        v_cub = p_cub(t)
        @test v_bsp ≈ v_cub atol = 1e-10
    end

    # Knot-time exactness
    for k in 1:N
        @test p_bsp(times[k]) ≈ u_vals[:, k] atol = 1e-14
    end
end

end # module Pulses
