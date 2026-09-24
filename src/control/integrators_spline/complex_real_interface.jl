#= ============================================================================
Complex ↔ Real Interface Utilities

Converts between complex-domain propagators/sensitivities and the real
isomorphic representations used by Ipopt.

The internal ODE runs in complex C^{n×n}. The external Jacobian/Hessian
seen by Ipopt is real. These utilities perform the conversion at the boundary.
============================================================================ =#

"""
    propagator_to_iso(Φ::AbstractMatrix{<:Complex})

Convert a complex n×n propagator to a real 2n×2n isomorphic block matrix.

Returns:
    [Re(Φ) -Im(Φ)]
    [Im(Φ)  Re(Φ)]

This is the real representation that acts on isomorphic ket vectors
ψ̃ = [Re(ψ); Im(ψ)].
"""
function propagator_to_iso(Φ::AbstractMatrix{<:Complex})
    n = size(Φ, 1)
    R = real(Φ)
    I_part = imag(Φ)
    Φ̃ = Matrix{Float64}(undef, 2n, 2n)
    Φ̃[1:n, 1:n] .= R
    Φ̃[1:n, (n+1):2n] .= .-I_part
    Φ̃[(n+1):2n, 1:n] .= I_part
    Φ̃[(n+1):2n, (n+1):2n] .= R
    return Φ̃
end

"""
    complex_ket_to_iso(ψ::AbstractVector{<:Complex})

Convert a complex ket to real isomorphic form: [Re(ψ); Im(ψ)].
"""
function complex_ket_to_iso(ψ::AbstractVector{<:Complex})
    return [real(ψ); imag(ψ)]
end

"""
    complex_ket_to_iso!(ψ̃::AbstractVector{<:Real}, ψ::AbstractVector{<:Complex})

In-place conversion of a complex ket to real isomorphic form. Zero allocations —
the matrix-free MultiKet kernels call this once per (knot, ket) per matvec, where
the allocating form showed up as knot-scaled hot-path allocation (#307).
"""
function complex_ket_to_iso!(ψ̃::AbstractVector{<:Real}, ψ::AbstractVector{<:Complex})
    n = length(ψ)
    @inbounds for i = 1:n
        ψ̃[i] = real(ψ[i])
        ψ̃[n+i] = imag(ψ[i])
    end
    return ψ̃
end

"""
    iso_to_complex_ket(ψ̃::AbstractVector{<:Real})

Convert a real isomorphic ket to complex: ψ̃[1:n] + im * ψ̃[n+1:2n].
"""
function iso_to_complex_ket(ψ̃::AbstractVector{<:Real})
    n = length(ψ̃) ÷ 2
    return ψ̃[1:n] + im * ψ̃[(n+1):2n]
end

"""
    iso_to_complex_ket!(ψ::AbstractVector{<:Complex}, ψ̃::AbstractVector{<:Real})

In-place conversion of real isomorphic ket to complex. Zero allocations.
"""
function iso_to_complex_ket!(ψ::AbstractVector{<:Complex}, ψ̃::AbstractVector{<:Real})
    n = length(ψ)
    @inbounds for i = 1:n
        ψ[i] = ψ̃[i] + im * ψ̃[n+i]
    end
    return ψ
end

"""
    iso_gather_ket!(ψ, src, cols) -> ψ

Gather a blocked-`[Re; Im]` iso ket out of `src` at the slab-relative columns
`cols` (length `2n`), writing the complex result into `ψ` (length `n`).

The allocation-free composition of `view(src, cols)` and
`iso_to_complex_ket!`: the matrix-free inner kernels read a knot's state,
tangent or cotangent block through a precomputed column table, and constructing a
`SubArray` over a `Vector{Int}` index per knot per matvec is exactly the kind of
per-knot term the knot-flat contract forbids (#338).
"""
@inline function iso_gather_ket!(
    ψ::AbstractVector{<:Complex},
    src::AbstractVector{<:Real},
    cols::AbstractVector{Int},
)
    n = length(ψ)
    @inbounds for i = 1:n
        ψ[i] = src[cols[i]] + im * src[cols[n+i]]
    end
    return ψ
end

"""
    iso_vec_to_complex_operator(Ũ⃗::AbstractVector{<:Real}, n::Int)

Convert a real isomorphic unitary vector to a complex n×n matrix.
The iso vec layout is column-major with interleaved [Re(col); Im(col)] blocks.
"""
function iso_vec_to_complex_operator(Ũ⃗::AbstractVector{<:Real}, n::Int)
    U = Matrix{ComplexF64}(undef, n, n)
    for i = 0:(n-1)
        U[:, i+1] .= @view(Ũ⃗[i*2n .+ (1:n)]) .+ im .* @view(Ũ⃗[i*2n .+ ((n+1):2n)])
    end
    return U
end

"""
    complex_operator_to_iso_vec!(Ũ⃗::AbstractVector{<:Real}, U::AbstractMatrix{<:Complex})

Convert a complex n×n operator to a real isomorphic vector, writing in-place.
"""
function complex_operator_to_iso_vec!(
    Ũ⃗::AbstractVector{<:Real},
    U::AbstractMatrix{<:Complex},
)
    n = size(U, 1)
    for i = 0:(n-1)
        Ũ⃗[i*2n .+ (1:n)] .= real.(@view(U[:, i+1]))
        Ũ⃗[i*2n .+ ((n+1):2n)] .= imag.(@view(U[:, i+1]))
    end
    return Ũ⃗
end

# ============================================================================
# Jacobian Assembly Utilities
# ============================================================================

"""
    sensitivity_to_jac_col!(col::AbstractVector{<:Real}, S_j, ψ_complex, tmp_complex)

Compute a real Jacobian column from complex sensitivity and state.

Computes: `col = ket_to_iso(-S_j * ψ_complex)`

# Arguments
- `col`: output real vector of length 2n
- `S_j`: complex n×n sensitivity matrix (∂Φ/∂p_j)
- `ψ_complex`: complex n-vector (state in complex form)
- `tmp_complex`: preallocated complex n-vector for intermediate result
"""
function sensitivity_to_jac_col!(
    col::AbstractVector{<:Real},
    S_j::AbstractMatrix{<:Complex},
    ψ_complex::AbstractVector{<:Complex},
    tmp_complex::AbstractVector{<:Complex},
)
    mul!(tmp_complex, S_j, ψ_complex)
    n = length(tmp_complex)
    @inbounds for i = 1:n
        col[i] = -real(tmp_complex[i])
        col[n+i] = -imag(tmp_complex[i])
    end
    return col
end

# ============================================================================
# Hessian Assembly Utilities
# ============================================================================

"""
    sensitivity_to_hess_col!(col::AbstractVector{<:Real}, S_j, μ_real, tmp1, tmp2)

Compute a real Hessian cross-term column from complex sensitivity and real multiplier.

Computes: `col = ket_to_iso(-S_j' * iso_to_ket(μ_real))`

# Arguments
- `col`: output real vector of length 2n
- `S_j`: complex n×n sensitivity matrix
- `μ_real`: real 2n-vector (Lagrange multiplier in iso form)
- `tmp1`: preallocated complex n-vector (used for μ_C input)
- `tmp2`: preallocated complex n-vector (used for result)
"""
function sensitivity_to_hess_col!(
    col::AbstractVector{<:Real},
    S_j::AbstractMatrix{<:Complex},
    μ_real::AbstractVector{<:Real},
    tmp1::AbstractVector{<:Complex},
    tmp2::AbstractVector{<:Complex},
)
    n = length(tmp1)
    # Convert real multiplier to complex: tmp1 = iso_to_ket(μ_real)
    @inbounds for i = 1:n
        tmp1[i] = μ_real[i] + im * μ_real[n+i]
    end
    # Compute -S_j' * μ_C (using BLAS gemv! to avoid adjoint wrapper allocations)
    BLAS.gemv!('C', one(ComplexF64), S_j, tmp1, zero(ComplexF64), tmp2)
    @inbounds for i = 1:n
        col[i] = -real(tmp2[i])
        col[n+i] = -imag(tmp2[i])
    end
    return col
end

"""
    hessian_pp_contraction(T_ij, μ_real, ψ_real, tmp1, tmp2)

Compute the (p,p) Hessian contraction: `-real(μ_C' * T_{ij} * ψ_C)`.

This is the scalar contribution of the second-order sensitivity matrix T_{ij}
to the Hessian of Lagrangian entry `∂²L/∂pᵢ∂pⱼ`.

# Arguments
- `T_ij`: complex n×n second-order sensitivity matrix (∂²Φ/∂pᵢ∂pⱼ)
- `μ_real`: real 2n-vector (Lagrange multiplier in iso form)
- `ψ_real`: real 2n-vector (state in iso form)
- `tmp1`: preallocated complex n-vector workspace
- `tmp2`: preallocated complex n-vector workspace
"""
function hessian_pp_contraction(
    T_ij::AbstractMatrix{<:Complex},
    μ_real::AbstractVector{<:Real},
    ψ_real::AbstractVector{<:Real},
    tmp1::AbstractVector{<:Complex},
    tmp2::AbstractVector{<:Complex},
)
    n = length(tmp1)
    # Convert state to complex: tmp1 = iso_to_ket(ψ_real)
    @inbounds for i = 1:n
        tmp1[i] = ψ_real[i] + im * ψ_real[n+i]
    end
    # Compute T_{ij} * ψ_C
    mul!(tmp2, T_ij, tmp1)
    # Convert multiplier to complex: tmp1 = iso_to_ket(μ_real)
    @inbounds for i = 1:n
        tmp1[i] = μ_real[i] + im * μ_real[n+i]
    end
    # Return -real(μ_C' * T_{ij} * ψ_C)
    val = zero(Float64)
    @inbounds for i = 1:n
        val += real(conj(tmp1[i]) * tmp2[i])
    end
    return -val
end
