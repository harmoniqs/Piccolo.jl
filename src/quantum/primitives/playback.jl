module Playback

"""
Hardware playback of designed pulses.

A designed pulse is a continuous-time envelope (or a held-sample pulse). Control
hardware plays something else: samples on a DAC clock lattice, held or interpolated
by the DAC, quantized to the DAC's integer codes, clipped at full scale, analog
filtered and delayed. This module models that chain:

- [`PlaybackSpec`](@ref) — the hardware description (sample period, DAC
  interpolation, bits / full scale / clip policy, causal filter, delay, gain);
- [`play`](@ref) — render a pulse into the waveform the hardware plays;
- [`playback_jacobian`](@ref) — the sparse linear map ∂(played samples)/∂(design
  parameters), exact everywhere except at clip and quantization, which are flagged;
- [`design_parameters`](@ref) / [`with_design_parameters`](@ref) — the documented
  design-parameter vector of each pulse type, and its inverse.

Reference hardware facts (the presets a spec is built from):

- **QICK** (tProc v2 / v6 generators): one envelope sample per
  `1/(samps_per_clk · f_fabric)` (= `1/fs`), int16 envelopes with full-scale code
  `maxv = 32766` (`bits = 16, max_code = 32766`); the `sg_int4` generator
  interpolates ×4 with a FIR (`dac_interp = taps, interp_factor = 4`).
- **Pasqal**: a 4 ns-clock zero-order hold followed by an ≈ 8 MHz output modulation
  filter (rise ≈ 85 ns) — `dt_sample = 4` (ns) with a causal low-pass `filter`.
- **Spin PLL lattices**: 1–10 ns holds, `dac_interp = :hold`, plus clipping at full
  scale.
"""

export PlaybackSpec, play, playback_jacobian
export AbstractPlaybackFilter, FIRFilter, IIRFilter, RCLowpass, ButterworthLowpass
export design_parameters, with_design_parameters
export playback_to_json, playback_from_json

using JSON3
using LinearAlgebra
using SparseArrays
using SpecialFunctions: erf
using TestItems

using ..Pulses

# ============================================================================ #
# Filters
# ============================================================================ #

"""
    AbstractPlaybackFilter

A causal, linear time-invariant output filter. Every filter is applied in
**discrete time** to the played sample sequence, discretized at the output sample
period (`dt_sample`, or `dt_sample / interp_factor` behind a FIR interpolator).
"""
abstract type AbstractPlaybackFilter end

"""
    FIRFilter(taps)

Causal FIR filter `y[i] = Σₖ taps[k] · x[i-k+1]`, taps at the output sample period.
"""
struct FIRFilter <: AbstractPlaybackFilter
    taps::Vector{Float64}
end

"""
    IIRFilter(b, a)

Causal IIR filter in direct form, `a[1] y[i] = Σₖ b[k] x[i-k+1] - Σₖ₍ₖ≥₂₎ a[k] y[i-k+1]`,
coefficients at the output sample period.
"""
struct IIRFilter <: AbstractPlaybackFilter
    b::Vector{Float64}
    a::Vector{Float64}
    function IIRFilter(b, a)
        (isempty(a) || a[1] == 0) && throw(ArgumentError("IIRFilter needs a[1] ≠ 0"))
        return new(b, a)
    end
end

"""
    RCLowpass(tau)

First-order RC low-pass `1/(τs + 1)` with time constant `tau` (spec time units),
discretized step-invariantly: `y[k] = α y[k-1] + (1-α) x[k-1]`, `α = exp(-dt/τ)`.
Driven by a held input it is exact at the sample times (a step plays
`1 - exp(-t/τ)`), at the cost of one sample of latency.
"""
struct RCLowpass <: AbstractPlaybackFilter
    tau::Float64
end

"""
    ButterworthLowpass(order, cutoff)

Butterworth low-pass of the given `order` with −3 dB `cutoff` in cycles per spec time
unit, discretized by the bilinear transform with cutoff prewarping (unity DC gain).
`cutoff` must lie below the output Nyquist frequency.
"""
struct ButterworthLowpass <: AbstractPlaybackFilter
    order::Int
    cutoff::Float64
end

# (b, a) at output period `dt` (spec time units)
_discretize(f::FIRFilter, dt) = (f.taps, [1.0])
_discretize(f::IIRFilter, dt) = (f.b, f.a)

function _discretize(f::RCLowpass, dt)
    α = exp(-dt / f.tau)
    return ([0.0, 1 - α], [1.0, -α])
end

function _discretize(f::ButterworthLowpass, dt)
    n = f.order
    n ≥ 1 || throw(ArgumentError("Butterworth order must be ≥ 1"))
    0 < f.cutoff * dt < 0.5 || throw(
        ArgumentError(
            "Butterworth cutoff $(f.cutoff) must lie below the Nyquist frequency $(0.5 / dt)",
        ),
    )
    ωa = 2 / dt * tan(π * f.cutoff * dt)
    poles = [ωa * cis(π * (2k + n - 1) / (2n)) for k = 1:n]
    a = real.(_poly([(1 + p * dt / 2) / (1 - p * dt / 2) for p in poles]))
    b = real.(_poly(fill(-1.0 + 0im, n)))
    return (b .* (sum(a) / sum(b)), a)
end

# coefficients in z⁻¹ of Π (1 - r z⁻¹)
_poly(roots) = foldl((c, r) -> [c; 0] .- r .* [0; c], roots; init = ComplexF64[1])

function _lfilter(b, a, x)
    y = zeros(length(x))
    for i in eachindex(x)
        acc = 0.0
        for k in eachindex(b)
            i - k + 1 ≥ 1 && (acc += b[k] * x[i-k+1])
        end
        for k = 2:length(a)
            i - k + 1 ≥ 1 && (acc -= a[k] * y[i-k+1])
        end
        y[i] = acc / a[1]
    end
    return y
end

# ============================================================================ #
# PlaybackSpec
# ============================================================================ #

"""
    PlaybackSpec(; dt_sample, kwargs...)

How control hardware plays a pulse. All times (`dt_sample`, `t0`, `delay`, filter
time constants / cutoffs) are in **spec time units**; `time_unit` converts them to
the pulse's time units (`t_pulse = time_unit · t_spec`, e.g. a spec in µs played
against a pulse in ns has `time_unit = 1000`).

# Keyword arguments
- `dt_sample`: DAC sample period (one envelope sample per `dt_sample`).
- `time_unit = 1.0`: pulse time units per spec time unit.
- `t0 = 0.0`: lattice origin, `0 ≤ t0 < dt_sample`; samples are taken at
  `t0 + k·dt_sample`.
- `t0_align = 1`: the played sample count is padded with zeros to a multiple of
  `t0_align` (e.g. a generator's samples per fabric clock).
- `dac_interp = :hold`: the DAC's own interpolation — `:hold` (zero-order hold),
  `:linear`, or a vector of FIR interpolation taps applied after zero-stuffing by
  `interp_factor` (QICK `sg_int4`: `interp_factor = 4`). Boxcar taps `ones(L)`
  reproduce the hold.
- `interp_factor = 1`: upsampling factor of a FIR interpolator.
- `bits = nothing`: DAC resolution; `nothing` disables quantization.
- `max_code = 2^(bits-1) - 1`: the integer code of `full_scale` (QICK: 32766).
- `full_scale = Inf`: amplitude (after `gain`) of the largest code; values beyond
  it are clipped.
- `clip = :clip`: `:clip` saturates at `±full_scale` and flags the samples,
  `:error` throws.
- `filter = nothing`: causal output filter ([`AbstractPlaybackFilter`](@ref)).
- `delay = 0.0`: output latency; the played waveform is shifted later by `delay`.
- `gain = nothing`: per-channel gain applied before clipping (`nothing` = unity).
"""
struct PlaybackSpec
    dt_sample::Float64
    time_unit::Float64
    t0::Float64
    t0_align::Int
    dac_interp::Union{Symbol,Vector{Float64}}
    interp_factor::Int
    bits::Union{Nothing,Int}
    max_code::Union{Nothing,Int}
    full_scale::Float64
    clip::Symbol
    filter::Union{Nothing,AbstractPlaybackFilter}
    delay::Float64
    gain::Union{Nothing,Vector{Float64}}
end

function PlaybackSpec(;
    dt_sample::Real,
    time_unit::Real = 1.0,
    t0::Real = 0.0,
    t0_align::Integer = 1,
    dac_interp::Union{Symbol,AbstractVector{<:Real}} = :hold,
    interp_factor::Integer = 1,
    bits::Union{Nothing,Integer} = nothing,
    max_code::Union{Nothing,Integer} = nothing,
    full_scale::Union{Nothing,Real} = Inf,
    clip::Symbol = :clip,
    filter::Union{Nothing,AbstractPlaybackFilter} = nothing,
    delay::Real = 0.0,
    gain::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    fs = isnothing(full_scale) ? Inf : Float64(full_scale)
    dt_sample > 0 || throw(ArgumentError("dt_sample must be positive"))
    time_unit > 0 || throw(ArgumentError("time_unit must be positive"))
    0 ≤ t0 < dt_sample || throw(ArgumentError("t0 must satisfy 0 ≤ t0 < dt_sample"))
    t0_align ≥ 1 || throw(ArgumentError("t0_align must be ≥ 1"))
    dac_interp isa Symbol &&
        dac_interp ∉ (:hold, :linear) &&
        throw(ArgumentError("dac_interp must be :hold, :linear or FIR taps"))
    dac_interp isa AbstractVector &&
        isempty(dac_interp) &&
        throw(ArgumentError("FIR dac_interp needs at least one tap"))
    interp_factor ≥ 1 || throw(ArgumentError("interp_factor must be ≥ 1"))
    interp_factor > 1 &&
        dac_interp isa Symbol &&
        throw(ArgumentError("interp_factor > 1 requires FIR dac_interp taps"))
    fs > 0 || throw(ArgumentError("full_scale must be positive"))
    clip ∈ (:clip, :error) || throw(ArgumentError("clip must be :clip or :error"))
    isnothing(bits) ||
        (bits ≥ 2 && isfinite(fs)) ||
        throw(ArgumentError("quantization needs bits ≥ 2 and a finite full_scale"))
    delay ≥ 0 || throw(ArgumentError("delay must be ≥ 0"))
    code = isnothing(bits) ? nothing : something(max_code, 2^(bits - 1) - 1)
    return PlaybackSpec(
        dt_sample,
        time_unit,
        t0,
        t0_align,
        dac_interp isa Symbol ? dac_interp : Vector{Float64}(dac_interp),
        interp_factor,
        bits,
        code,
        fs,
        clip,
        filter,
        delay,
        isnothing(gain) ? nothing : Vector{Float64}(gain),
    )
end

const _ValueTypes = Union{PlaybackSpec,FIRFilter,IIRFilter,RCLowpass,ButterworthLowpass}

Base.:(==)(a::T, b::T) where {T<:_ValueTypes} =
    all(getfield(a, f) == getfield(b, f) for f in fieldnames(T))
Base.hash(a::T, h::UInt) where {T<:_ValueTypes} =
    foldr(hash, [getfield(a, f) for f in fieldnames(T)]; init = hash(T, h))

# output sample period in spec time units
_dt_out(spec::PlaybackSpec) =
    spec.dac_interp isa Symbol ? spec.dt_sample : spec.dt_sample / spec.interp_factor

# ============================================================================ #
# JSON / Dict
# ============================================================================ #

"""
    Dict(spec::PlaybackSpec) -> Dict{String,Any}

JSON-ready representation (`full_scale = Inf` is stored as `nothing`; filters as
`{"type": "fir"|"iir"|"rc"|"butterworth", ...}`). Inverse: `PlaybackSpec(dict)`.
"""
Base.Dict(spec::PlaybackSpec) = Dict{String,Any}(
    "dt_sample" => spec.dt_sample,
    "time_unit" => spec.time_unit,
    "t0" => spec.t0,
    "t0_align" => spec.t0_align,
    "dac_interp" =>
        spec.dac_interp isa Symbol ? String(spec.dac_interp) : spec.dac_interp,
    "interp_factor" => spec.interp_factor,
    "bits" => spec.bits,
    "max_code" => spec.max_code,
    "full_scale" => isfinite(spec.full_scale) ? spec.full_scale : nothing,
    "clip" => String(spec.clip),
    "filter" => _filter_dict(spec.filter),
    "delay" => spec.delay,
    "gain" => spec.gain,
)

_filter_dict(::Nothing) = nothing
_filter_dict(f::FIRFilter) = Dict{String,Any}("type" => "fir", "taps" => f.taps)
_filter_dict(f::IIRFilter) = Dict{String,Any}("type" => "iir", "b" => f.b, "a" => f.a)
_filter_dict(f::RCLowpass) = Dict{String,Any}("type" => "rc", "tau" => f.tau)
_filter_dict(f::ButterworthLowpass) =
    Dict{String,Any}("type" => "butterworth", "order" => f.order, "cutoff" => f.cutoff)

"""
    PlaybackSpec(d::AbstractDict)

Build a spec from a `Dict` (String or Symbol keys, e.g. a parsed JSON object).
Missing keys take the keyword defaults; `dt_sample` is required.
"""
function PlaybackSpec(d::AbstractDict)
    s = Dict{String,Any}(string(k) => v for (k, v) in d)
    floats(v) = isnothing(v) ? nothing : Vector{Float64}(v)
    interp = get(s, "dac_interp", "hold")
    return PlaybackSpec(;
        dt_sample = s["dt_sample"],
        time_unit = get(s, "time_unit", 1.0),
        t0 = get(s, "t0", 0.0),
        t0_align = get(s, "t0_align", 1),
        dac_interp = interp isa AbstractString ? Symbol(interp) : floats(interp),
        interp_factor = get(s, "interp_factor", 1),
        bits = get(s, "bits", nothing),
        max_code = get(s, "max_code", nothing),
        full_scale = get(s, "full_scale", nothing),
        clip = Symbol(get(s, "clip", "clip")),
        filter = _filter_from(get(s, "filter", nothing)),
        delay = get(s, "delay", 0.0),
        gain = floats(get(s, "gain", nothing)),
    )
end

_filter_from(::Nothing) = nothing

function _filter_from(d::AbstractDict)
    s = Dict{String,Any}(string(k) => v for (k, v) in d)
    type = s["type"]
    type == "fir" && return FIRFilter(Vector{Float64}(s["taps"]))
    type == "iir" && return IIRFilter(Vector{Float64}(s["b"]), Vector{Float64}(s["a"]))
    type == "rc" && return RCLowpass(s["tau"])
    type == "butterworth" && return ButterworthLowpass(s["order"], s["cutoff"])
    throw(ArgumentError("unknown playback filter type \"$type\""))
end

"""
    playback_to_json(spec::PlaybackSpec) -> String

Serialize a spec to JSON (see `Dict(spec)` for the layout).
"""
playback_to_json(spec::PlaybackSpec) = JSON3.write(Dict(spec))

"""
    playback_from_json(json::AbstractString) -> PlaybackSpec

Parse a spec from a JSON string.
"""
playback_from_json(json::AbstractString) = PlaybackSpec(JSON3.read(json, Dict{String,Any}))

# ============================================================================ #
# play
# ============================================================================ #

"""
    play(spec::PlaybackSpec, pulse::AbstractPulse) -> (played, info)

Render `pulse` into the waveform the hardware plays.

1. **Lattice.** Sample `pulse` at `t0 + k·dt` (`dt = time_unit · dt_sample`; never
   `T/(n-1)`). A `:hold` / FIR DAC takes the samples in `[0, T)` (each held for one
   period, played duration in `[T, T + dt)`); a `:linear` DAC takes those in
   `[0, T]` (played duration in `(T - dt, T]`). The count is zero-padded to a multiple
   of `t0_align`.
2. **Digital.** Multiply by the per-channel `gain`; clip to `±full_scale` (flagging
   the samples, or throwing under `clip = :error`); quantize to
   `round(v / full_scale · max_code)` codes when `bits` is set.
3. **DAC interpolation.** FIR taps zero-stuff by `interp_factor` and convolve
   causally (output period `dt / interp_factor`); `:hold` and `:linear` pass the
   sequence through.
4. **Filter.** The causal `filter` is applied to the output sample sequence in
   discrete time, discretized at the output period. The output is truncated at the
   end of the played window (no ring-out tail).
5. **Delay.** The waveform is shifted later by `delay` (zero before it).

`played` is a `ZeroOrderPulse` (`:hold` / FIR) or a `LinearSplinePulse` (`:linear`)
through the output samples. `info` carries:

- `samples`, `times`: output samples (`n_drives × N_out`) and their play times;
- `dt_out`: output sample period (pulse time units);
- `input_samples`, `input_times`: the designed pulse on the lattice (after padding);
- `clipped::BitMatrix`: which lattice samples saturated (`n_drives × N`);
- `codes`: integer DAC codes (`nothing` without quantization);
- `quantized::Bool`.
"""
function play(spec::PlaybackSpec, pulse::AbstractPulse)
    r = _render(spec, pulse)
    n = size(r.y, 1)
    dt_o = spec.time_unit * _dt_out(spec)
    Δ = spec.time_unit * (spec.t0 + spec.delay)
    N = size(r.y, 2)
    knots = collect(Δ .+ dt_o .* (0:(N-1)))
    name = hasproperty(pulse, :drive_name) ? pulse.drive_name : :u
    played = if spec.dac_interp === :linear
        # the DAC ramps up from the (zero) sample before the first one
        pre = Δ ≤ 0 ? Float64[] : (Δ - dt_o > 0 ? [0.0, Δ - dt_o] : [0.0])
        LinearSplinePulse([zeros(n, length(pre)) r.y], [pre; knots]; drive_name = name)
    else
        pre = Δ > 0 ? [0.0] : Float64[]
        ZeroOrderPulse(
            [zeros(n, length(pre)) r.y r.y[:, end]],
            [pre; knots; Δ + N * dt_o];
            drive_name = name,
        )
    end
    info = (
        samples = r.y,
        times = knots,
        dt_out = dt_o,
        input_samples = r.x,
        input_times = r.times,
        clipped = r.clipped,
        codes = r.codes,
        quantized = !isnothing(r.codes),
    )
    return played, info
end

function _render(spec::PlaybackSpec, pulse::AbstractPulse)
    dt = spec.time_unit * spec.dt_sample
    t0 = spec.time_unit * spec.t0
    T = duration(pulse)
    tol = 1e-9 * dt
    n_in = if spec.dac_interp === :linear
        floor(Int, (T - t0 + tol) / dt) + 1
    else
        max(ceil(Int, (T - t0 - tol) / dt), 1)
    end
    n_pad = cld(n_in, spec.t0_align) * spec.t0_align
    times = collect(t0 .+ dt .* (0:(n_pad-1)))
    n = n_drives(pulse)
    x = [sample(pulse, times[1:n_in]) zeros(n, n_pad - n_in)]
    g = something(spec.gain, ones(n))
    length(g) == n ||
        throw(DimensionMismatch("gain has $(length(g)) channels, pulse has $n drives"))
    fs = spec.full_scale
    v = g .* x
    clipped = abs.(v) .> fs
    spec.clip === :error &&
        any(clipped) &&
        throw(
            ArgumentError(
                "$(count(clipped)) played samples exceed full_scale = $fs (clip = :error)",
            ),
        )
    v = clamp.(v, -fs, fs)
    codes = isnothing(spec.bits) ? nothing : round.(Int, v ./ fs .* spec.max_code)
    v = isnothing(codes) ? v : codes .* (fs / spec.max_code)
    y = reduce(vcat, [permutedims(_chain(spec, v[c, :])) for c = 1:n])
    return (; times, x, y, clipped, codes, n_in, n_pad, g)
end

# DAC interpolation + filter, one channel (linear, time-invariant)
function _chain(spec::PlaybackSpec, s::AbstractVector)
    w = if spec.dac_interp isa Symbol
        Vector{Float64}(s)
    else
        L = spec.interp_factor
        z = zeros(length(s) * L)
        z[1:L:end] .= s
        _lfilter(spec.dac_interp, [1.0], z)
    end
    isnothing(spec.filter) && return w
    b, a = _discretize(spec.filter, _dt_out(spec))
    return _lfilter(b, a, w)
end

# ============================================================================ #
# Jacobian
# ============================================================================ #

"""
    playback_jacobian(spec, pulse; quantization = :error) -> (J, flags)

Sparse Jacobian `J = ∂ vec(info.samples) / ∂ design_parameters(pulse)` of
[`play`](@ref), with rows ordered like `vec(info.samples)` (drive fastest, then
output sample) and columns like [`design_parameters`](@ref).

`J = C · D · S`: `S` samples the pulse on the lattice (local: two knots per sample
for linear, four parameters for cubic), `D` is the per-channel gain, and `C` is the
DAC-interpolation + filter convolution (block-Toeplitz per channel; impulse-response
entries below `1e-20` of its peak are dropped).

Clipping and quantization are not differentiable and are never differentiated
through silently:

- clipped samples carry the saturation derivative (zero) and are reported in
  `flags.clipped` (`n_drives × N` lattice samples);
- with `bits` set, the call throws unless `quantization = :straight_through`, which
  treats rounding as the identity (`flags.quantization === :straight_through`).
"""
function playback_jacobian(
    spec::PlaybackSpec,
    pulse::AbstractPulse;
    quantization::Symbol = :error,
)
    quantization ∈ (:error, :straight_through) ||
        throw(ArgumentError("quantization must be :error or :straight_through"))
    !isnothing(spec.bits) &&
        quantization === :error &&
        throw(
            ArgumentError(
                "quantization is not differentiable; pass `quantization = :straight_through` " *
                "to differentiate it as the identity",
            ),
        )
    r = _render(spec, pulse)
    n = n_drives(pulse)
    S0 = _sample_jacobian(pulse, r.times[1:r.n_in])
    S = [S0; spzeros(n * (r.n_pad - r.n_in), size(S0, 2))]
    h = _chain(spec, [1.0; zeros(r.n_pad - 1)])
    peak = maximum(abs, h)
    keep = peak == 0 ? 0 : findlast(x -> abs(x) > 1e-20 * peak, h)
    L = spec.dac_interp isa Symbol ? 1 : spec.interp_factor
    N_out = length(h)
    rows, cols, vals = Int[], Int[], Float64[]
    for j = 1:r.n_pad, c = 1:n
        r.clipped[c, j] && continue
        for m = 1:min(keep, N_out-(j-1)*L)
            iszero(h[m]) && continue
            push!(rows, ((j - 1) * L + m - 1) * n + c)
            push!(cols, (j - 1) * n + c)
            push!(vals, h[m] * r.g[c])
        end
    end
    C = sparse(rows, cols, vals, n * N_out, n * r.n_pad)
    flags = (
        clipped = r.clipped,
        quantization = isnothing(spec.bits) ? :none : :straight_through,
    )
    return C * S, flags
end

# ∂ vec(sample(pulse, times)) / ∂ design_parameters(pulse)
function _sample_jacobian(p::ZeroOrderPulse, times)
    knots = get_knot_times(p)
    n = n_drives(p)
    # left hold, with `evaluate`'s 1e-12 snap onto a knot
    ks = [clamp(searchsortedlast(knots, t + 1e-12), 1, length(knots)) for t in times]
    rows = [(j - 1) * n + c for j in eachindex(times) for c = 1:n]
    cols = [(ks[j] - 1) * n + c for j in eachindex(times) for c = 1:n]
    return sparse(rows, cols, 1.0, n * length(times), n * length(knots))
end

function _sample_jacobian(p::LinearSplinePulse, times)
    knots = get_knot_times(p)
    n = n_drives(p)
    rows, cols, vals = Int[], Int[], Float64[]
    for (j, t) in enumerate(times)
        k, w = _interval(knots, t)
        for c = 1:n
            append!(rows, ((j - 1) * n + c, (j - 1) * n + c))
            append!(cols, ((k - 1) * n + c, k * n + c))
            append!(vals, (1 - w, w))
        end
    end
    return sparse(rows, cols, vals, n * length(times), n * length(knots))
end

function _sample_jacobian(p::CubicSplinePulse, times)
    knots = get_knot_times(p)
    n = n_drives(p)
    rows, cols, vals = Int[], Int[], Float64[]
    for (j, t) in enumerate(times)
        k, s = _interval(knots, t)
        h = knots[k+1] - knots[k]
        # Hermite basis: u_k, du_k, u_{k+1}, du_{k+1}
        basis = (2s^3 - 3s^2 + 1, (s^3 - 2s^2 + s) * h, -2s^3 + 3s^2, (s^3 - s^2) * h)
        for c = 1:n
            row = (j - 1) * n + c
            append!(rows, (row, row, row, row))
            append!(
                cols,
                ((k - 1) * 2n + c, (k - 1) * 2n + n + c, k * 2n + c, k * 2n + n + c),
            )
            append!(vals, basis)
        end
    end
    return sparse(rows, cols, vals, n * length(times), 2n * length(knots))
end

function _interval(knots, t)
    k = clamp(searchsortedlast(knots, t), 1, length(knots) - 1)
    return k, clamp((t - knots[k]) / (knots[k+1] - knots[k]), 0.0, 1.0)
end

function _sample_jacobian(p::GaussianPulse, times)
    n = n_drives(p)
    a, σ, μ = p.amplitudes, p.sigmas, p.centers
    rows, cols, vals = Int[], Int[], Float64[]
    for (j, t) in enumerate(times), i = 1:n
        e = exp(-(t - μ[i])^2 / (2σ[i]^2))
        row = (j - 1) * n + i
        append!(rows, (row, row, row))
        append!(cols, (i, n + i, 2n + i))
        append!(vals, (e, a[i] * e * (t - μ[i])^2 / σ[i]^3, a[i] * e * (t - μ[i]) / σ[i]^2))
    end
    return sparse(rows, cols, vals, n * length(times), 3n)
end

function _sample_jacobian(p::ErfPulse, times)
    n = n_drives(p)
    a, σ, μ = p.amplitudes, p.sigmas, p.centers
    rows, cols, vals = Int[], Int[], Float64[]
    for (j, t) in enumerate(times), i = 1:n
        z = √2 * (t - μ[i]) / σ[i]
        d = 2 / √π * exp(-z^2)
        row = (j - 1) * n + i
        append!(rows, (row, row, row))
        append!(cols, (i, n + i, 2n + i))
        append!(vals, (erf(z), -a[i] * d * z / σ[i], -a[i] * d * √2 / σ[i]))
    end
    return sparse(rows, cols, vals, n * length(times), 3n)
end

function _sample_jacobian(p::CompositePulse, times)
    n = n_drives(p)
    nt = length(times)
    blocks = map(p.pulses, p.drive_mapping) do sub, mapping
        m = length(mapping)
        rows = [(j - 1) * n + mapping[d] for j = 1:nt for d = 1:m]
        sparse(rows, 1:(m*nt), 1.0, n * nt, m * nt) * _sample_jacobian(sub, times)
    end
    return reduce(hcat, blocks)
end

_sample_jacobian(p::FunctionPulse, times) = spzeros(n_drives(p) * length(times), 0)

_sample_jacobian(p::AbstractPulse, _) =
    throw(ArgumentError("no design parameterization is defined for $(typeof(p))"))

# ============================================================================ #
# Design parameters
# ============================================================================ #

"""
    design_parameters(pulse) -> Vector{Float64}

The design-parameter vector that [`playback_jacobian`](@ref) differentiates against,
knot-major like a `NamedTrajectory`'s packed data:

- `ZeroOrderPulse`, `LinearSplinePulse`: `vec(u)` (knot values, drive fastest);
- `CubicSplinePulse`: `vec(vcat(u, du))` — per knot, the values then the Hermite
  tangents (the `[u; du]` component order of a spline trajectory);
- `GaussianPulse`, `ErfPulse`: `[amplitudes; sigmas; centers]`;
- `CompositePulse`: the components' vectors concatenated in order;
- `FunctionPulse`: empty (no design parameters).
"""
design_parameters(p::ZeroOrderPulse) = vec(Matrix(p.controls.u))
design_parameters(p::LinearSplinePulse) = vec(Matrix(get_knot_values(p)))
design_parameters(p::CubicSplinePulse) =
    vec(vcat(get_knot_values(p), get_knot_derivatives(p)))
design_parameters(p::Union{GaussianPulse,ErfPulse}) = [p.amplitudes; p.sigmas; p.centers]
design_parameters(p::CompositePulse) = reduce(vcat, design_parameters.(p.pulses))
design_parameters(::FunctionPulse) = Float64[]

"""
    with_design_parameters(pulse, θ) -> AbstractPulse

A pulse of the same type, knot times and metadata as `pulse` with design parameters
`θ` (the inverse of [`design_parameters`](@ref)).
"""
function with_design_parameters(p::ZeroOrderPulse, θ::AbstractVector)
    knots = collect(get_knot_times(p))
    return ZeroOrderPulse(
        reshape(θ, n_drives(p), length(knots)),
        knots;
        drive_name = p.drive_name,
        initial_value = p.initial_value,
        final_value = p.final_value,
        snap_to_knots = p.snap_to_knots,
    )
end

function with_design_parameters(p::LinearSplinePulse, θ::AbstractVector)
    knots = collect(get_knot_times(p))
    return LinearSplinePulse(
        reshape(θ, n_drives(p), length(knots)),
        knots;
        drive_name = p.drive_name,
        initial_value = p.initial_value,
        final_value = p.final_value,
    )
end

function with_design_parameters(p::CubicSplinePulse, θ::AbstractVector)
    knots = collect(get_knot_times(p))
    n = n_drives(p)
    M = reshape(θ, 2n, length(knots))
    return CubicSplinePulse(
        M[1:n, :],
        M[(n+1):end, :],
        knots;
        drive_name = p.drive_name,
        initial_value = p.initial_value,
        final_value = p.final_value,
    )
end

function with_design_parameters(p::GaussianPulse, θ::AbstractVector)
    n = n_drives(p)
    return GaussianPulse(
        θ[1:n],
        θ[(n+1):2n],
        θ[(2n+1):3n],
        duration(p);
        drive_name = p.drive_name,
    )
end

function with_design_parameters(p::ErfPulse, θ::AbstractVector)
    n = n_drives(p)
    return ErfPulse(θ[1:n], θ[(n+1):2n], θ[(2n+1):3n], duration(p))
end

function with_design_parameters(p::CompositePulse, θ::AbstractVector)
    ends = cumsum(length.(design_parameters.(p.pulses)))
    subs = [
        with_design_parameters(sub, θ[(e-length(design_parameters(sub))+1):e]) for
        (sub, e) in zip(p.pulses, ends)
    ]
    widths = cumsum(n_drives.(subs))
    concatenated = [collect((w-n_drives(s)+1):w) for (s, w) in zip(subs, widths)]
    return CompositePulse(
        subs,
        p.drive_mapping == concatenated ? :concatenate : :interleave,
    )
end

with_design_parameters(p::FunctionPulse, θ::AbstractVector) = p

# ============================================================================= #
# Tests
# ============================================================================= #

@testitem "play: samples on the clock lattice t0 + k·dt for every pulse type" begin
    T = 10.3
    dt = 1.0
    knots = collect(range(0, T, length = 7))
    u = [sin(0.4 * t) for _ = 1:1, t in knots]
    pulses = [
        ZeroOrderPulse(u, knots),
        LinearSplinePulse(u, knots),
        CubicSplinePulse(u, 0.4 .* cos.(0.4 .* u), knots),
        GaussianPulse([0.7], 2.0, T),
        ErfPulse([0.3], 2.0, T),
        FunctionPulse(t -> [sin(t)], T, 1),
        CompositePulse([GaussianPulse([0.7], 2.0, T), ErfPulse([0.3], 2.0, T)]),
    ]
    for p in pulses
        # hold: samples strictly inside [0, T), each held for one period
        spec = PlaybackSpec(dt_sample = dt)
        played, info = play(spec, p)
        @test info.input_times ≈ collect(0.0:dt:10.0)
        @test length(info.input_times) == 11            # not range(0, T, length = n)
        @test info.input_samples ≈ hcat([p(t) for t in info.input_times]...)
        @test abs(duration(played) - duration(p)) ≤ dt
        @test played isa ZeroOrderPulse

        # linear: samples on [0, T], interpolated
        played_l, info_l = play(PlaybackSpec(dt_sample = dt, dac_interp = :linear), p)
        @test played_l isa LinearSplinePulse
        @test info_l.input_times ≈ collect(0.0:dt:10.0)
        @test abs(duration(played_l) - duration(p)) ≤ dt

        # lattice offset t0: samples at t0 + k·dt
        _, info_o = play(PlaybackSpec(dt_sample = dt, t0 = 0.25), p)
        @test info_o.input_times ≈ collect(0.25:dt:10.25)
        @test info_o.input_samples ≈ hcat([p(t) for t in info_o.input_times]...)
    end

    # time_unit: spec times in µs, pulse in ns
    p = GaussianPulse([0.7], 20.0, 100.0)
    _, info = play(PlaybackSpec(dt_sample = 0.004, time_unit = 1000.0), p)
    @test info.input_times ≈ collect(0.0:4.0:96.0)

    # t0_align pads the sample count to a clock multiple with zeros
    _, info =
        play(PlaybackSpec(dt_sample = 1.0, t0_align = 4), GaussianPulse([0.7], 2.0, T))
    @test size(info.samples, 2) == 12
    @test all(iszero, info.samples[:, 12])
end

@testitem "play: hold, linear and FIR DAC interpolation" begin
    T = 8.0
    knots = collect(0.0:1.0:T)
    u = [0.1 * k - 0.02 * k^2 for _ = 1:2, k = 0:8]
    u[2, :] .*= -1
    zop = ZeroOrderPulse(u, knots)

    # hold of a ZOH pulse on its own grid reproduces it exactly
    played, info = play(PlaybackSpec(dt_sample = 1.0), zop)
    @test info.samples == u[:, 1:8]
    @test duration(played) == T
    for t in range(0, T - 1e-9, length = 97)
        @test played(t) == zop(t)
    end

    # linear interpolation of a linear pulse on its own grid reproduces it
    lin = LinearSplinePulse(u, knots)
    played_l, info_l = play(PlaybackSpec(dt_sample = 1.0, dac_interp = :linear), lin)
    @test info_l.samples ≈ u
    for t in range(0, T, length = 97)
        @test played_l(t) ≈ lin(t) atol = 1e-14
    end

    # linear DAC on a coarser grid linearly interpolates between lattice samples
    lin_c, info_c = play(PlaybackSpec(dt_sample = 2.0, dac_interp = :linear), lin)
    @test info_c.input_times ≈ collect(0.0:2.0:8.0)
    @test lin_c(3.0) ≈ (lin(2.0) + lin(4.0)) / 2

    # FIR ×4 interpolation: boxcar taps reproduce the hold at dt/4
    L = 4
    fir, info_f =
        play(PlaybackSpec(dt_sample = 1.0, dac_interp = ones(L), interp_factor = L), zop)
    @test info_f.dt_out ≈ 0.25
    @test size(info_f.samples, 2) == 8L
    @test info_f.samples ≈ repeat(u[:, 1:8], inner = (1, L))
    @test fir isa ZeroOrderPulse

    # FIR triangle taps: causal linear interpolation (one input period of group delay)
    tri = [1:L; (L-1):-1:1] ./ L
    _, info_t =
        play(PlaybackSpec(dt_sample = 1.0, dac_interp = tri, interp_factor = L), zop)
    s = u[1, 1:8]
    expected = [
        let k = fld(i, L), r = mod(i, L)
            (k ≥ 1 ? s[k] : 0.0) * (1 - (r + 1) / L) +
            (k + 1 ≤ 8 ? s[k+1] : 0.0) * ((r + 1) / L)
        end for i = 0:(8L-1)
    ]
    @test info_t.samples[1, :] ≈ expected
end

@testitem "play: quantization, clip policy and per-channel gain" begin
    T = 4.0
    knots = collect(0.0:1.0:T)
    u = [0.10 0.50 1.30 -1.70 0.0; 0.2 0.3 0.4 0.5 0.0]
    zop = ZeroOrderPulse(u, knots)

    # quantization to int16 at QICK full scale (maxv = 32766)
    spec = PlaybackSpec(dt_sample = 1.0, bits = 16, max_code = 32766, full_scale = 2.0)
    _, info = play(spec, zop)
    @test info.codes == round.(Int, u[:, 1:4] ./ 2.0 .* 32766)
    @test info.samples ≈ info.codes .* (2.0 / 32766)
    @test maximum(abs.(info.samples .- u[:, 1:4])) ≤ 2.0 / 32766 / 2 + 1e-15
    @test !any(info.clipped)

    # default max code is 2^(bits-1) - 1
    _, info3 = play(PlaybackSpec(dt_sample = 1.0, bits = 3, full_scale = 1.0), zop)
    @test info3.samples[2, :] ≈ round.(u[2, 1:4] .* 3) ./ 3

    # :clip saturates at ±full_scale and records which samples clipped
    _, info_c = play(PlaybackSpec(dt_sample = 1.0, full_scale = 1.0, clip = :clip), zop)
    @test info_c.clipped == Bool[0 0 1 1; 0 0 0 0]
    @test info_c.samples[1, :] ≈ [0.10, 0.50, 1.0, -1.0]

    # :error throws
    @test_throws ArgumentError play(
        PlaybackSpec(dt_sample = 1.0, full_scale = 1.0, clip = :error),
        zop,
    )

    # per-channel gain is applied before clipping
    _, info_g =
        play(PlaybackSpec(dt_sample = 1.0, gain = [0.5, 3.0], full_scale = 1.0), zop)
    @test info_g.samples[1, :] ≈ 0.5 .* u[1, 1:4]
    @test info_g.clipped[2, :] == Bool[0, 0, 1, 1]
    @test info_g.samples[2, :] ≈ [0.6, 0.9, 1.0, 1.0]
end

@testitem "play: causal low-pass filters and delay" begin
    T = 200.0
    dt = 0.25
    step = FunctionPulse(t -> [1.0], T, 1)

    # first-order RC (step-invariant discretization): exact at the sample times
    τ = 7.0
    _, info = play(PlaybackSpec(dt_sample = dt, filter = RCLowpass(τ)), step)
    t = info.times
    @test info.samples[1, :] ≈ 1 .- exp.(-t ./ τ) atol = 1e-12

    # 2nd-order Butterworth (bilinear, prewarped): matches the analog step response
    fc = 0.01
    _, info_b = play(PlaybackSpec(dt_sample = dt, filter = ButterworthLowpass(2, fc)), step)
    ω = 2π * fc
    ζ = 1 / √2
    ωd = ω * √(1 - ζ^2)
    analog = @. 1 - exp(-ζ * ω * t) * (cos(ωd * t) + ζ / √(1 - ζ^2) * sin(ωd * t))
    @test maximum(abs.(info_b.samples[1, :] .- analog)) < 1e-2
    @test maximum(info_b.samples) ≈ 1 + exp(-π * ζ / √(1 - ζ^2)) atol = 2e-3  # 4.3% overshoot
    @test info_b.samples[1, end] ≈ 1 atol = 1e-3                              # unity DC gain

    # 4th-order Butterworth settles (more slowly) to unity
    long_step = FunctionPulse(t -> [1.0], 1000.0, 1)
    _, info_4 =
        play(PlaybackSpec(dt_sample = dt, filter = ButterworthLowpass(4, fc)), long_step)
    @test info_4.samples[1, end] ≈ 1 atol = 1e-3
    @test all(abs.(info_4.samples[1, 1:4]) .< 1e-3)                            # causal

    # explicit FIR and IIR filters agree with their difference equations
    _, info_f =
        play(PlaybackSpec(dt_sample = dt, filter = FIRFilter([0.5, 0.25, 0.25])), step)
    @test info_f.samples[1, 1:4] ≈ [0.5, 0.75, 1.0, 1.0]
    _, info_i =
        play(PlaybackSpec(dt_sample = dt, filter = IIRFilter([0.5], [1.0, -0.5])), step)
    @test info_i.samples[1, 1:4] ≈ [0.5, 0.75, 0.875, 0.9375]

    # delay shifts the played waveform by the specified amount
    g = GaussianPulse([1.0], 5.0, 40.0)
    p0, _ = play(PlaybackSpec(dt_sample = 1.0, filter = RCLowpass(2.0)), g)
    for d in (3.0, 2.5)
        pd, _ = play(PlaybackSpec(dt_sample = 1.0, filter = RCLowpass(2.0), delay = d), g)
        @test duration(pd) ≈ duration(p0) + d
        @test all(iszero, pd(t)[1] for t in range(0, d - 1e-9, length = 9))
        for t in range(0, duration(p0) - 1e-9, length = 83)
            @test pd(t + d) ≈ p0(t)
        end
    end
end

@testitem "playback_jacobian matches central finite differences" begin
    using SparseArrays
    using LinearAlgebra

    T = 20.0
    knots = collect(range(0, T, length = 9))
    u = [0.3 * sin(0.3 * t + c) for c = 1:2, t in knots]
    du = [0.09 * cos(0.3 * t + c) for c = 1:2, t in knots]
    pulses = [
        ZeroOrderPulse(u, knots),
        LinearSplinePulse(u, knots),
        CubicSplinePulse(u, du, knots),
        GaussianPulse([0.5, -0.4], [3.0, 4.0], [9.0, 11.0], T),
        ErfPulse([0.2, 0.3], [5.0, 4.0], [10.0, 9.0], T),
        CompositePulse([GaussianPulse([0.5], 3.0, T), LinearSplinePulse(u[1:1, :], knots)]),
    ]
    specs = [
        PlaybackSpec(dt_sample = 0.7, t0 = 0.1),
        PlaybackSpec(dt_sample = 0.7, dac_interp = :linear, filter = RCLowpass(1.5)),
        PlaybackSpec(
            dt_sample = 0.7,
            dac_interp = [0.25, 0.5, 0.75, 1.0, 0.75, 0.5, 0.25],
            interp_factor = 4,
            filter = ButterworthLowpass(3, 0.2),
            gain = [0.9, 1.1],
            delay = 0.3,
            full_scale = 10.0,
        ),
    ]
    for p in pulses, spec in specs
        θ = design_parameters(p)
        J, flags = playback_jacobian(spec, p)
        @test J isa SparseMatrixCSC
        @test !any(flags.clipped)
        _, info = play(spec, p)
        @test size(J) == (length(info.samples), length(θ))
        f(θ′) = vec(play(spec, with_design_parameters(p, θ′))[2].samples)
        @test f(θ) == vec(info.samples)
        h = 1e-5
        J_fd = hcat(
            [
                (f(θ .+ h .* e) .- f(θ .- h .* e)) ./ (2h) for
                e in eachcol(Matrix(I, length(θ), length(θ)))
            ]...,
        )
        @test norm(Matrix(J) .- J_fd) ≤ 1e-8 * max(norm(J_fd), 1.0)
    end

    # knot-major ordering: θ = vec(vcat(u, du)) for cubic, vec(u) for ZOH / linear
    @test design_parameters(pulses[3]) == vec(vcat(u, du))
    @test design_parameters(pulses[2]) == vec(u)

    # the map is sparse for a hold DAC with no filter: one nonzero per played sample
    J, _ = playback_jacobian(PlaybackSpec(dt_sample = 0.7), pulses[1])
    @test nnz(J) == size(J, 1)

    # clipped samples are flagged and their rows carry the saturation derivative (zero)
    spec_c = PlaybackSpec(dt_sample = 0.7, full_scale = 0.2)
    Jc, flags = playback_jacobian(spec_c, pulses[2])
    @test any(flags.clipped)
    @test all(iszero, Jc[findall(vec(flags.clipped)), :])

    # quantization is non-differentiable: refused unless straight-through is requested
    spec_q = PlaybackSpec(dt_sample = 0.7, bits = 16, full_scale = 1.0)
    @test_throws ArgumentError playback_jacobian(spec_q, pulses[2])
    Jq, flags_q = playback_jacobian(spec_q, pulses[2]; quantization = :straight_through)
    @test flags_q.quantization === :straight_through
    @test Jq ≈
          playback_jacobian(PlaybackSpec(dt_sample = 0.7, full_scale = 1.0), pulses[2])[1]
end

@testitem "PlaybackSpec round-trips through JSON" begin
    specs = [
        PlaybackSpec(dt_sample = 0.002604, time_unit = 1000.0),
        PlaybackSpec(
            dt_sample = 4.0,
            t0 = 1.0,
            t0_align = 16,
            dac_interp = [0.25, 0.5, 0.75, 1.0, 0.75, 0.5, 0.25],
            interp_factor = 4,
            bits = 16,
            max_code = 32766,
            full_scale = 1.5,
            clip = :error,
            filter = ButterworthLowpass(2, 0.008),
            delay = 12.0,
            gain = [1.0, 0.98],
        ),
        PlaybackSpec(dt_sample = 2.0, dac_interp = :linear, filter = RCLowpass(20.0)),
        PlaybackSpec(dt_sample = 1.0, filter = FIRFilter([0.5, 0.5])),
        PlaybackSpec(dt_sample = 1.0, filter = IIRFilter([0.5], [1.0, -0.5])),
    ]
    for spec in specs
        json = playback_to_json(spec)
        @test json isa String
        @test playback_from_json(json) == spec
        @test PlaybackSpec(Dict(spec)) == spec
    end
end

@testitem "played_fidelity: design-vs-played gap shrinks with the lattice step" begin
    using LinearAlgebra
    sys = QuantumSystem(zeros(ComplexF64, 2, 2), [ComplexF64[0 1; 1 0]], [(-2.0, 2.0)])
    N, T = 11, 5.0
    times = collect(range(0, T, length = N))
    # smooth envelope with area ≈ π/2 · (x drive), so F is far from trivial
    u = [0.5 * π / T * 2 * sin(π * t / T)^2 for _ = 1:1, t in times]
    du = [0.5 * π / T * 2 * π / T * sin(2π * t / T) for _ = 1:1, t in times]
    pulse = CubicSplinePulse(u, du, times)
    qtraj = UnitaryTrajectory(sys, pulse, GATES[:X])
    qcp = SplinePulseProblem(qtraj, N)

    gaps = map([0.5, 0.1, 0.02]) do dt
        r = played_fidelity(qcp, PlaybackSpec(dt_sample = dt))
        @test 0 ≤ r.F_played ≤ 1
        @test r.gap ≈ r.F_design - r.F_played
        @test r.played isa ZeroOrderPulse
        abs(r.gap)
    end
    @test gaps[3] < 1e-3
    @test gaps[1] > gaps[2] > gaps[3]
end

end
