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
export ZeroOrderPulse, LinearSplinePulse, CubicSplinePulse, GaussianPulse, ErfPulse, CompositePulse
export duration, n_drives, sample, drive_name

using DataInterpolations: ConstantInterpolation, LinearInterpolation, CubicHermiteSpline
using ForwardDiff
using SpecialFunctions: erf
using TestItems

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
    sample(pulse::AbstractPulse; n_samples::Int=100)

Sample the pulse uniformly with `n_samples` points. Returns `(controls, times)`.
"""
function sample(pulse::AbstractPulse; n_samples::Int=100)
    times = collect(range(0.0, duration(pulse), length=n_samples))
    return sample(pulse, times), times
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
end

"""
    ZeroOrderPulse(controls::AbstractMatrix, times::AbstractVector; drive_name=:u, initial_value=nothing, final_value=nothing)

Create a zero-order hold pulse from control samples and times.

# Arguments
- `controls`: Matrix of size `(n_drives, n_times)` with control values
- `times`: Vector of sample times (must start at 0)

# Keyword Arguments
- `drive_name`: Name of the drive variable (default `:u`)
- `initial_value`: Initial boundary condition (default: zeros(n_drives))
- `final_value`: Final boundary condition (default: zeros(n_drives))
"""
function ZeroOrderPulse(controls::AbstractMatrix, times::AbstractVector; drive_name::Symbol=:u, initial_value::Union{Nothing, Vector{<:Real}}=nothing, final_value::Union{Nothing, Vector{<:Real}}=nothing)
    n_drives = size(controls, 1)
    init_val = isnothing(initial_value) ? zeros(n_drives) : Vector{Float64}(initial_value)
    final_val = isnothing(final_value) ? zeros(n_drives) : Vector{Float64}(final_value)
    # Materialize to Matrix/Vector to ensure consistent type parameters
    interp = ConstantInterpolation(Matrix(controls), collect(times))
    return ZeroOrderPulse(interp, Float64(times[end]), n_drives, drive_name, init_val, final_val)
end

evaluate(p::ZeroOrderPulse, t) = p.controls(t)

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
function LinearSplinePulse(controls::AbstractMatrix, times::AbstractVector; drive_name::Symbol=:u, initial_value::Union{Nothing, Vector{<:Real}}=nothing, final_value::Union{Nothing, Vector{<:Real}}=nothing)
    n_drives = size(controls, 1)
    init_val = isnothing(initial_value) ? zeros(n_drives) : Vector{Float64}(initial_value)
    final_val = isnothing(final_value) ? zeros(n_drives) : Vector{Float64}(final_value)
    # Materialize to Matrix/Vector to ensure consistent type parameters
    interp = LinearInterpolation(Matrix(controls), collect(times))
    return LinearSplinePulse(interp, Float64(times[end]), n_drives, drive_name, init_val, final_val)
end

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
function CubicSplinePulse(controls::AbstractMatrix, derivatives::AbstractMatrix, times::AbstractVector; drive_name::Symbol=:u, initial_value::Union{Nothing, Vector{<:Real}}=nothing, final_value::Union{Nothing, Vector{<:Real}}=nothing)
    n_drives = size(controls, 1)
    init_val = isnothing(initial_value) ? zeros(n_drives) : Vector{Float64}(initial_value)
    final_val = isnothing(final_value) ? zeros(n_drives) : Vector{Float64}(final_value)
    # Materialize to Matrix to ensure consistent type parameters across construction methods
    interp = CubicHermiteSpline(Matrix(derivatives), Matrix(controls), collect(times))
    return CubicSplinePulse(interp, Float64(times[end]), n_drives, drive_name, init_val, final_val)
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
function CubicSplinePulse(controls::AbstractMatrix, times::AbstractVector; drive_name::Symbol=:u, initial_value::Union{Nothing, Vector{<:Real}}=nothing, final_value::Union{Nothing, Vector{<:Real}}=nothing)
    derivatives = zeros(size(controls))
    return CubicSplinePulse(controls, derivatives, times; drive_name, initial_value, final_value)
end

evaluate(p::CubicSplinePulse, t) = p.controls(t)

# ============================================================================ #
# Spline pulse knot accessors
# ============================================================================ #

"""
    get_knot_times(pulse::AbstractSplinePulse)

Return the knot times stored in the spline pulse interpolant.
"""
get_knot_times(p::LinearSplinePulse) = p.controls.t
get_knot_times(p::CubicSplinePulse) = p.controls.t

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
    drive_name::Symbol=:u,
    initial_value::Union{Nothing, Vector{<:Real}}=nothing,
    final_value::Union{Nothing, Vector{<:Real}}=nothing
)
    times = collect(range(0.0, duration(pulse), length=n_samples))
    return LinearSplinePulse(pulse, times; drive_name, initial_value, final_value)
end

function LinearSplinePulse(
    pulse::AbstractPulse,
    times::AbstractVector;
    drive_name::Symbol=:u,
    initial_value::Union{Nothing, Vector{<:Real}}=nothing,
    final_value::Union{Nothing, Vector{<:Real}}=nothing
)
    controls = sample(pulse, times)
    init_val = isnothing(initial_value) ? pulse(times[1]) : initial_value
    final_val = isnothing(final_value) ? pulse(times[end]) : final_value
    return LinearSplinePulse(controls, times; drive_name, initial_value=init_val, final_value=final_val)
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
    drive_name::Symbol=:du,
    initial_value::Union{Nothing, Vector{<:Real}}=nothing,
    final_value::Union{Nothing, Vector{<:Real}}=nothing
)
    times = collect(range(0.0, duration(pulse), length=n_samples))
    return CubicSplinePulse(pulse, times; drive_name, initial_value, final_value)
end

function CubicSplinePulse(
    pulse::AbstractPulse,
    times::AbstractVector;
    drive_name::Symbol=:du,
    initial_value::Union{Nothing, Vector{<:Real}}=nothing,
    final_value::Union{Nothing, Vector{<:Real}}=nothing
)
    controls = sample(pulse, times)
    
    # Compute derivatives using ForwardDiff
    derivatives = hcat([ForwardDiff.derivative(pulse, t) for t in times]...)
    
    init_val = isnothing(initial_value) ? pulse(times[1]) : initial_value
    final_val = isnothing(final_value) ? pulse(times[end]) : final_value
    
    return CubicSplinePulse(controls, derivatives, times; drive_name, initial_value=init_val, final_value=final_val)
end

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
    duration::Real
)
    n = length(amplitudes)
    @assert length(sigmas) == n "sigmas must have same length as amplitudes"
    @assert length(centers) == n "centers must have same length as amplitudes"
    
    amps = Vector{Float64}(amplitudes)
    sigs = Vector{Float64}(sigmas)
    ctrs = Vector{Float64}(centers)
    dur = Float64(duration)
    
    f = t -> [amps[i] * exp(-(t - ctrs[i])^2 / (2 * sigs[i]^2)) for i in 1:n]
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
    center::Real=duration/2
)
    n = length(amplitudes)
    return GaussianPulse(
        amplitudes, 
        fill(Float64(sigma), n), 
        fill(Float64(center), n), 
        duration
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
    duration::Real
)
    n = length(amplitudes)
    @assert length(sigmas) == n "sigmas must have same length as amplitudes"
    @assert length(centers) == n "centers must have same length as amplitudes"
    
    amps = Vector{Float64}(amplitudes)
    sigs = Vector{Float64}(sigmas)
    ctrs = Vector{Float64}(centers)
    dur = Float64(duration)
    
    f = t -> [amps[i] * erf(√2 * (t - ctrs[i]) / sigs[i]) for i in 1:n]
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
    center::Real=duration/2
)
    n = length(amplitudes)
    return ErfPulse(
        amplitudes, 
        fill(Float64(sigma), n), 
        fill(Float64(center), n), 
        duration
    )
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
function CompositePulse(
    pulses::Vector{<:AbstractPulse},
    mode::Symbol=:interleave
)
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
    
    drive_mapping = [Vector{Int}() for _ in 1:n_pulses]
    
    if mode == :interleave
        # Check all pulses have same number of drives for interleaving
        @assert allequal(n_drives_per_pulse) "All pulses must have same n_drives for :interleave mode"
        
        n_per_pulse = n_drives_per_pulse[1]
        composite_idx = 1
        
        # Interleave: [p1_d1, p2_d1, ..., pN_d1, p1_d2, p2_d2, ..., pN_d2, ...]
        for drive_idx in 1:n_per_pulse
            for pulse_idx in 1:n_pulses
                push!(drive_mapping[pulse_idx], composite_idx)
                composite_idx += 1
            end
        end
    else  # :concatenate
        # Concatenate: [p1_d1, p1_d2, ..., p2_d1, p2_d2, ...]
        composite_idx = 1
        for pulse_idx in 1:n_pulses
            for _ in 1:n_drives_per_pulse[pulse_idx]
                push!(drive_mapping[pulse_idx], composite_idx)
                composite_idx += 1
            end
        end
    end
    
    # Create evaluation function (type-generic for ForwardDiff compatibility)
    f = function(t)
        # Get first pulse values to determine element type
        first_vals = pulses[1](t)
        T = eltype(first_vals)
        result = zeros(T, total_drives)
        
        # Fill in values from first pulse
        for (local_idx, composite_idx) in enumerate(drive_mapping[1])
            result[composite_idx] = first_vals[local_idx]
        end
        
        # Fill in values from remaining pulses
        for pulse_idx in 2:n_pulses
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
# NamedTrajectory Constructors
# ============================================================================ #

using NamedTrajectories: NamedTrajectory, get_times

"""
    ZeroOrderPulse(traj::NamedTrajectory; drive_name=:u)

Construct a ZeroOrderPulse from a NamedTrajectory.

# Arguments
- `traj`: NamedTrajectory with control data

# Keyword Arguments
- `drive_name`: Name of the drive component (default: `:u`)
"""
function ZeroOrderPulse(traj::NamedTrajectory; drive_name::Symbol=:u)
    controls = traj[drive_name]
    times = get_times(traj)
    return ZeroOrderPulse(controls, times; drive_name)
end

"""
    LinearSplinePulse(traj::NamedTrajectory; drive_name=:u)

Construct a LinearSplinePulse from a NamedTrajectory.

# Arguments
- `traj`: NamedTrajectory with control data

# Keyword Arguments
- `drive_name`: Name of the drive component (default: `:u`)
"""
function LinearSplinePulse(traj::NamedTrajectory; drive_name::Symbol=:u)
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
function CubicSplinePulse(traj::NamedTrajectory; drive_name::Symbol=:u, derivative_name::Symbol=:du)
    controls = traj[drive_name]
    derivatives = traj[derivative_name]
    times = get_times(traj)
    return CubicSplinePulse(controls, derivatives, times; drive_name)
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
    sampled, ts = sample(pulse; n_samples=5)
    @test size(sampled) == (2, 5)
    @test length(ts) == 5
    
    # Test custom drive_name
    pulse_custom = ZeroOrderPulse(controls, times; drive_name=:Ω)
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
    sampled, ts = sample(pulse; n_samples=5)
    @test size(sampled) == (2, 5)
    
    # Test custom drive_name
    pulse_custom = LinearSplinePulse(controls, times; drive_name=:amplitude)
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
    pulse_custom = CubicSplinePulse(controls, derivatives, times; drive_name=:a)
    @test drive_name(pulse_custom) == :a
    
    # Test zero-derivative constructor
    pulse_zero_deriv = CubicSplinePulse(controls, times)
    @test pulse_zero_deriv(0.5) ≈ [1.0, -1.0]
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
    pulse3 = GaussianPulse([1.0], 0.1, 1.0; center=0.3)
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
    pulse3 = ErfPulse([1.0], 0.2, 1.0; center=0.3)
    @test pulse3.centers == [0.3]
    @test pulse3(0.3)[1] ≈ 0.0  # Zero at center
    
    # Test phase compensation use case
    φ_max = π/4
    T = 50.0
    σ = T/4
    phase_pulse = ErfPulse([φ_max], [σ], [T/2], T)
    @test phase_pulse(T/2)[1] ≈ 0.0
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
    times = range(0, T, length=20)
    for t in times
        @test linear(t) ≈ gaussian(t) atol=1e-10
    end
    
    # Test with custom times
    custom_times = [0.0, 2.5, 5.0, 7.5, 10.0]
    linear2 = LinearSplinePulse(gaussian, custom_times)
    @test duration(linear2) == T
    for t in custom_times
        @test linear2(t) ≈ gaussian(t) atol=1e-10
    end
    
    # Test CubicSplinePulse conversion with derivatives
    cubic = CubicSplinePulse(gaussian, 20)
    @test cubic isa CubicSplinePulse
    @test duration(cubic) == T
    @test n_drives(cubic) == 2
    
    # Cubic spline should be very close to original at sample points
    for t in times
        @test cubic(t) ≈ gaussian(t) atol=1e-10
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
        @test linear_val ≈ orig atol=1e-1
        @test cubic_val ≈ orig atol=1e-1
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
    times = collect(range(0, 1, length=T))
    Δt = diff(times)
    Δt = [Δt; Δt[end]]
    
    u = [sin(π * t) for t in times]
    du = [π * cos(π * t) for t in times]
    
    traj = NamedTrajectory(
        (u = reshape(u, 1, :), du = reshape(du, 1, :), Δt = Δt);
        timestep=:Δt,
        controls=(:u,)
    )
    
    # Test ZeroOrderPulse from NamedTrajectory
    zop = ZeroOrderPulse(traj)
    @test duration(zop) ≈ 1.0
    @test n_drives(zop) == 1
    @test zop(0.0) ≈ [0.0] atol=1e-10
    
    # Test LinearSplinePulse from NamedTrajectory
    lp = LinearSplinePulse(traj)
    @test duration(lp) ≈ 1.0
    @test n_drives(lp) == 1
    @test lp(0.5) ≈ [1.0] atol=1e-10
    
    # Test CubicSplinePulse from NamedTrajectory (uses both u and du)
    csp = CubicSplinePulse(traj)
    @test duration(csp) ≈ 1.0
    @test n_drives(csp) == 1
    @test csp(0.5) ≈ [1.0] atol=1e-10
    @test csp(0.0) ≈ [0.0] atol=1e-10
end

@testitem "Pulse Boundary Conditions" begin
    using LinearAlgebra
    
    # Test ZeroOrderPulse boundary conditions
    n_drives = 2
    T = 1.0
    N = 10
    controls = rand(n_drives, N)
    times = range(0, T, length=N)
    
    # Default: zeros
    pulse_default = ZeroOrderPulse(controls, times)
    @test pulse_default.initial_value == zeros(n_drives)
    @test pulse_default.final_value == zeros(n_drives)
    
    # Custom values
    init_val = [1.0, 2.0]
    final_val = [3.0, 4.0]
    pulse_custom = ZeroOrderPulse(controls, times; initial_value=init_val, final_value=final_val)
    @test pulse_custom.initial_value == init_val
    @test pulse_custom.final_value == final_val
    
    # Test LinearSplinePulse boundary conditions
    pulse_linear_default = LinearSplinePulse(controls, times)
    @test pulse_linear_default.initial_value == zeros(n_drives)
    @test pulse_linear_default.final_value == zeros(n_drives)
    
    pulse_linear_custom = LinearSplinePulse(controls, times; initial_value=init_val, final_value=final_val)
    @test pulse_linear_custom.initial_value == init_val
    @test pulse_linear_custom.final_value == final_val
    
    # Test CubicSplinePulse boundary conditions
    derivatives = rand(n_drives, N)
    
    pulse_cubic_default = CubicSplinePulse(controls, derivatives, times)
    @test pulse_cubic_default.initial_value == zeros(n_drives)
    @test pulse_cubic_default.final_value == zeros(n_drives)
    
    pulse_cubic_custom = CubicSplinePulse(controls, derivatives, times; initial_value=init_val, final_value=final_val)
    @test pulse_cubic_custom.initial_value == init_val
    @test pulse_cubic_custom.final_value == final_val
    
    # Test CubicSplinePulse constructor without derivatives
    pulse_cubic_no_deriv = CubicSplinePulse(controls, times; initial_value=init_val, final_value=final_val)
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
    times = range(0, T, length=N)
    
    init_val = [0.5]
    final_val = [-0.5]
    pulse = CubicSplinePulse(controls, derivatives, times; initial_value=init_val, final_value=final_val)
    
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

end # module Pulses
