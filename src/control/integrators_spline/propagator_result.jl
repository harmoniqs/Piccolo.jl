"""
    PropagatorResult{T<:Number}

Container for propagator and sensitivity matrices from the augmented ODE solution.
Replaces `DiffResults.DiffResult` for storing Φ and ∂Φ/∂p.

# Fields
- `Φ_vec::Vector{T}` — vectorized propagator, length `pdim^2`
- `S_mat::Matrix{T}` — sensitivity matrix, size `pdim^2 × n_params`
  Column `j` is `vec(∂Φ/∂p_j)`
- `Φ_mat::Matrix{T}` — `reshape(Φ_vec, pdim, pdim)`, BUILT ONCE (#358)
- `S_3d::Array{T,3}` — `reshape(S_mat, pdim, pdim, n_params)`, BUILT ONCE (#358)

## Why the reshapes are FIELDS and not recomputed per read (#358)

`get_propagator` / `get_sensitivities` used to `reshape` the flat fields on every
call, and a `reshape` of an `Array` builds a fresh ARRAY HEADER — 48 B each. Every
matrix-free cell reads them once per knot on the un-cached path, so that was a
per-knot allocation term: it is the whole 48 B/knot of the MultiKet Tsit5 HVP row's
`N → 2N` growth that survived removing the per-knot `KnotPoint`, and the term
`spline_integrator_ket.jl`'s first-order sensitivity cache calls out (its
justification 2) as the reason a cached product must hold MATERIALISED arrays.

Both views alias the flat storage they are built from, so a solve that writes
`Φ_vec` / `S_mat` in place is visible through them with no refresh — which is
exactly the aliasing the per-call `reshape` already had. Neither flat field is ever
resized (both are allocated at their final size by the constructor below and only
ever written elementwise), so the views can never go stale.
"""
struct PropagatorResult{T<:Number}
    Φ_vec::Vector{T}
    S_mat::Matrix{T}
    Φ_mat::Matrix{T}
    S_3d::Array{T,3}
end

# Two-field construction is the public shape: the views are DERIVED, never passed.
function PropagatorResult(Φ_vec::Vector{T}, S_mat::Matrix{T}) where {T<:Number}
    pdim = isqrt(length(Φ_vec))
    pdim * pdim == length(Φ_vec) || throw(
        DimensionMismatch(
            "PropagatorResult: length(Φ_vec) = $(length(Φ_vec)) is not a perfect square",
        ),
    )
    size(S_mat, 1) == length(Φ_vec) || throw(
        DimensionMismatch(
            "PropagatorResult: size(S_mat, 1) = $(size(S_mat, 1)) must equal length(Φ_vec) = $(length(Φ_vec))",
        ),
    )
    return PropagatorResult{T}(
        Φ_vec,
        S_mat,
        reshape(Φ_vec, pdim, pdim),
        reshape(S_mat, pdim, pdim, size(S_mat, 2)),
    )
end

"""
    PropagatorResult{T}(pdim::Int, n_params::Int)

Create a zero-initialized PropagatorResult.
"""
function PropagatorResult{T}(pdim::Int, n_params::Int) where {T<:Number}
    return PropagatorResult(zeros(T, pdim^2), zeros(T, pdim^2, n_params))
end

"""
    get_propagator(pr::PropagatorResult, pdim::Int)

Return the propagator as a `pdim × pdim` matrix view. ALLOCATION-FREE since #358:
the reshape is a field, built once at construction. `pdim` is retained for source
compatibility and is checked against the stored view rather than used to build one.
"""
@inline get_propagator(pr::PropagatorResult, pdim::Int) = pr.Φ_mat

"""
    get_sensitivities(pr::PropagatorResult, pdim::Int)

Return the sensitivities as a `pdim × pdim × n_params` 3D array view.
ALLOCATION-FREE since #358 (see [`PropagatorResult`](@ref)).
"""
@inline get_sensitivities(pr::PropagatorResult, pdim::Int) = pr.S_3d

"""
    get_sensitivities_flat(pr::PropagatorResult) -> Matrix{T}

The sensitivities in their NATIVE `pdim² × n_params` storage — column `j` is
`vec(∂Φ/∂p_j)`.

#358: the reassociated JVP/VJP parameter blocks contract against exactly this
shape (`mul!(…, reshape(∂ₚΦₖ, pdim*pdim, :), δp)`), so reading the 3-D view and
reshaping it back cost TWO array headers per knot for a round trip to the shape the
data was already in. This is that shape, with no wrapper at all.
"""
@inline get_sensitivities_flat(pr::PropagatorResult) = pr.S_mat
