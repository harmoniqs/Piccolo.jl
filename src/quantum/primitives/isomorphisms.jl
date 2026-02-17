module Isomorphisms

# Export state isomorphisms
export mat
export ket_to_iso
export iso_to_ket
export density_to_iso_vec
export iso_vec_to_density
export density_to_compact_iso
export compact_iso_to_density
export density_lift_matrix
export density_projection_matrix
export iso_vec_to_operator
export iso_vec_to_iso_operator
export operator_to_iso_vec
export iso_operator_to_iso_vec
export iso_operator_to_operator
export operator_to_iso_operator
export ket_to_bloch
export bloch_to_ket
export density_to_compact_iso
export compact_iso_to_density
export density_lift_matrix
export density_projection_matrix

# Do not export Hamiltonian isomorphisms

using LinearAlgebra
using SparseArrays
using TestItems
using ..Gates


@doc raw"""
    mat(x::AbstractVector)

Convert a vector `x` into a square matrix. The length of `x` must be a perfect square.
"""
function mat(x::AbstractVector)
    n = isqrt(length(x))
    @assert n^2 ≈ length(x) "Vector length must be a perfect square"
    return reshape(x, n, n)
end


# ----------------------------------------------------------------------------- #
#                                Kets                                           #
# ----------------------------------------------------------------------------- #

@doc raw"""
    ket_to_iso(ψ::AbstractVector{<:Number})

Convert a ket vector `ψ` into a complex vector with real and imaginary parts.
"""
ket_to_iso(ψ::AbstractVector{<:Number}) = [real(ψ); imag(ψ)]

@doc raw"""
    iso_to_ket(ψ̃::AbstractVector{<:Real})

Convert a real isomorphism vector `ψ̃` into a ket vector.
"""
iso_to_ket(ψ̃::AbstractVector{<:Real}) =
    ψ̃[1:div(length(ψ̃), 2)] + im * ψ̃[(div(length(ψ̃), 2)+1):end]

# ----------------------------------------------------------------------------- #
#                             Unitaries                                         #
# ----------------------------------------------------------------------------- #

@doc raw"""
    iso_vec_to_operator(Ũ⃗::AbstractVector{ℝ}) where ℝ <: Real

Convert a real vector `Ũ⃗` into a complex matrix representing an operator.
"""
function iso_vec_to_operator(Ũ⃗::AbstractVector{ℝ}) where {ℝ<:Real}
    Ũ⃗_dim = div(length(Ũ⃗), 2)
    N = Int(sqrt(Ũ⃗_dim))
    U = Matrix{complex(ℝ)}(undef, N, N)
    for i = 0:(N-1)
        U[:, i+1] .= @view(Ũ⃗[i*2N.+(1:N)]) + one(ℝ) * im * @view(Ũ⃗[i*2N.+((N+1):2N)])
    end
    return U
end

@doc raw"""
    iso_vec_to_iso_operator(Ũ⃗::AbstractVector{ℝ}) where ℝ <: Real

Convert a real vector `Ũ⃗` into a real matrix representing an isomorphism operator.
"""
function iso_vec_to_iso_operator(Ũ⃗::AbstractVector{ℝ}) where {ℝ<:Real}
    N = Int(sqrt(length(Ũ⃗) ÷ 2))
    Ũ = Matrix{ℝ}(undef, 2N, 2N)
    U_real = Matrix{ℝ}(undef, N, N)
    U_imag = Matrix{ℝ}(undef, N, N)
    for i = 0:(N-1)
        U_real[:, i+1] .= @view(Ũ⃗[i*2N.+(1:N)])
        U_imag[:, i+1] .= @view(Ũ⃗[i*2N.+((N+1):2N)])
    end
    Ũ[1:N, 1:N] .= U_real
    Ũ[1:N, (N+1):end] .= -U_imag
    Ũ[(N+1):end, 1:N] .= U_imag
    Ũ[(N+1):end, (N+1):end] .= U_real
    return Ũ
end

@doc raw"""
    operator_to_iso_vec(U::AbstractMatrix{ℂ}) where ℂ <: Number

Convert a complex matrix `U` representing an operator into a real vector.
"""
function operator_to_iso_vec(U::AbstractMatrix{ℂ}) where {ℂ<:Number}
    N = size(U, 1)
    Ũ⃗ = Vector{real(ℂ)}(undef, N^2 * 2)
    for i = 0:(N-1)
        Ũ⃗[i*2N.+(1:N)] .= real(@view(U[:, i+1]))
        Ũ⃗[i*2N.+((N+1):2N)] .= imag(@view(U[:, i+1]))
    end
    return Ũ⃗
end

@doc raw"""
    iso_operator_to_iso_vec(Ũ::AbstractMatrix{ℝ}) where ℝ <: Real

Convert a real matrix `Ũ` representing an isomorphism operator into a real vector.
"""
function iso_operator_to_iso_vec(Ũ::AbstractMatrix{ℝ}) where {ℝ<:Real}
    N = size(Ũ, 1) ÷ 2
    Ũ⃗ = Vector{ℝ}(undef, N^2 * 2)
    for i = 0:(N-1)
        Ũ⃗[i*2N.+(1:2N)] .= @view Ũ[:, i+1]
    end
    return Ũ⃗
end

@doc raw"""
    iso_operator_to_operator(Ũ)
"""
iso_operator_to_operator(Ũ) = iso_vec_to_operator(iso_operator_to_iso_vec(Ũ))

@doc raw"""
    operator_to_iso_operator(U)
"""
operator_to_iso_operator(U) = iso_vec_to_iso_operator(operator_to_iso_vec(U))

# ----------------------------------------------------------------------------- #
#                             Density matrix                                    #
# ----------------------------------------------------------------------------- #

@doc raw"""
    density_to_iso_vec(ρ::AbstractMatrix{<:Number})

Returns the isomorphism `ρ⃗̃ = ket_to_iso(vec(ρ))` of a density matrix `ρ`
"""
density_to_iso_vec(ρ::AbstractMatrix{<:Number}) = ket_to_iso(vec(ρ))

@doc raw"""
    iso_vec_to_density(ρ⃗̃::AbstractVector{<:Real})

Returns the density matrix `ρ` from its isomorphism `ρ⃗̃`
"""
iso_vec_to_density(ρ⃗̃::AbstractVector{<:Real}) = mat(iso_to_ket(ρ⃗̃))

@doc raw"""
    density_to_compact_iso(ρ::AbstractMatrix{<:Number})

Returns a compact real isomorphism vector of length ``n^2`` for a Hermitian density matrix
``ρ`` of size ``n \times n``. Exploits Hermiticity (``ρ = ρ^\dagger``) to halve the
representation size compared to [`density_to_iso_vec`](@ref).

The compact vector is ordered as:
1. ``\text{Re}(ρ_{jk})`` for ``j \leq k`` (upper triangle, column-major): ``n(n+1)/2`` entries
2. ``\text{Im}(ρ_{jk})`` for ``j < k`` (strict upper triangle, column-major): ``n(n-1)/2`` entries

See also [`compact_iso_to_density`](@ref), [`density_lift_matrix`](@ref),
[`density_projection_matrix`](@ref).
"""
function density_to_compact_iso(ρ::AbstractMatrix{<:Number})
    n = size(ρ, 1)
    x = Vector{Float64}(undef, n^2)
    idx = 0
    # Re upper triangle (column-major)
    @inbounds for k = 1:n, j = 1:k
        idx += 1
        x[idx] = real(ρ[j, k])
    end
    # Im strict upper triangle (column-major)
    @inbounds for k = 2:n, j = 1:k-1
        idx += 1
        x[idx] = imag(ρ[j, k])
    end
    return x
end

@doc raw"""
    compact_iso_to_density(x::AbstractVector{<:Real})

Reconstructs a Hermitian density matrix ``ρ`` from its compact real isomorphism vector ``x``
of length ``n^2``.

See also [`density_to_compact_iso`](@ref).
"""
function compact_iso_to_density(x::AbstractVector{<:Real})
    n = isqrt(length(x))
    T = complex(eltype(x))
    ρ = Matrix{T}(undef, n, n)
    idx = 0
    # Re upper triangle (column-major)
    @inbounds for k = 1:n, j = 1:k
        idx += 1
        ρ[j, k] = x[idx]
        if j != k
            ρ[k, j] = x[idx]
        end
    end
    # Im strict upper triangle (column-major)
    @inbounds for k = 2:n, j = 1:k-1
        idx += 1
        ρ[j, k] += im * x[idx]
        ρ[k, j] -= im * x[idx]
    end
    return ρ
end

@doc raw"""
    density_lift_matrix(n::Int)

Returns the sparse lift matrix ``L \in \mathbb{R}^{2n^2 \times n^2}`` that maps the compact
density isomorphism to the full ``\text{iso\_vec}`` representation.

Given a compact vector ``x = \text{density\_to\_compact\_iso}(ρ)``, the full representation
is recovered as ``L x = \text{density\_to\_iso\_vec}(ρ)``.

The matrix has ``n(2n-1)`` nonzero entries.

See also [`density_projection_matrix`](@ref), [`density_to_compact_iso`](@ref).
"""
function density_lift_matrix(n::Int)
    n² = n^2
    rows = Int[]
    cols = Int[]
    vals = Float64[]

    col_idx = 0

    # Re upper triangle (column-major)
    @inbounds for k = 1:n, j = 1:k
        col_idx += 1
        # Re(ρ[j,k]) position in iso_vec: (k-1)*n + j
        re_jk = (k - 1) * n + j
        push!(rows, re_jk)
        push!(cols, col_idx)
        push!(vals, 1.0)
        if j != k
            # Re(ρ[k,j]) = Re(ρ[j,k]) — symmetric
            re_kj = (j - 1) * n + k
            push!(rows, re_kj)
            push!(cols, col_idx)
            push!(vals, 1.0)
        end
    end

    # Im strict upper triangle (column-major)
    @inbounds for k = 2:n, j = 1:k-1
        col_idx += 1
        # Im(ρ[j,k]) position in iso_vec: n² + (k-1)*n + j
        im_jk = n² + (k - 1) * n + j
        push!(rows, im_jk)
        push!(cols, col_idx)
        push!(vals, 1.0)
        # Im(ρ[k,j]) = -Im(ρ[j,k]) — antisymmetric
        im_kj = n² + (j - 1) * n + k
        push!(rows, im_kj)
        push!(cols, col_idx)
        push!(vals, -1.0)
    end

    return sparse(rows, cols, vals, 2n², n²)
end

@doc raw"""
    density_projection_matrix(n::Int)

Returns the sparse projection matrix ``P \in \mathbb{R}^{n^2 \times 2n^2}`` that maps the
full ``\text{iso\_vec}`` representation to the compact density isomorphism.

Given a full vector ``v = \text{density\_to\_iso\_vec}(ρ)`` for Hermitian ``ρ``, the compact
representation is ``P v = \text{density\_to\_compact\_iso}(ρ)``.

Satisfies ``PL = I_{n^2}`` where ``L`` = [`density_lift_matrix`](@ref).

See also [`density_lift_matrix`](@ref), [`density_to_compact_iso`](@ref).
"""
function density_projection_matrix(n::Int)
    n² = n^2
    rows = Int[]
    cols = Int[]
    vals = Float64[]

    row_idx = 0

    # Re upper triangle (column-major)
    @inbounds for k = 1:n, j = 1:k
        row_idx += 1
        # Re(ρ[j,k]) position in iso_vec: (k-1)*n + j
        re_jk = (k - 1) * n + j
        push!(rows, row_idx)
        push!(cols, re_jk)
        push!(vals, 1.0)
    end

    # Im strict upper triangle (column-major)
    @inbounds for k = 2:n, j = 1:k-1
        row_idx += 1
        # Im(ρ[j,k]) position in iso_vec: n² + (k-1)*n + j
        im_jk = n² + (k - 1) * n + j
        push!(rows, row_idx)
        push!(cols, im_jk)
        push!(vals, 1.0)
    end

    return sparse(rows, cols, vals, n², 2n²)
end

# ----------------------------------------------------------------------------- #
#                             Hamiltonians                                      #
# ----------------------------------------------------------------------------- #

const Im2 = [
    0 -1
    1 0
]

@doc raw"""
    iso(H::AbstractMatrix{<:Number})

Returns the isomorphism of ``H``:

```math
iso(H) = \widetilde{H} = \mqty(1 & 0 \\ 0 & 1) \otimes \Re(H) + \mqty(0 & -1 \\ 1 & 0) \otimes \Im(H)
```

where ``\Im(H)`` and ``\Re(H)`` are the imaginary and real parts of ``H`` and the tilde 
indicates the standard isomorphism of a complex valued matrix:

```math
\widetilde{H} = \mqty(1 & 0 \\ 0 & 1) \otimes \Re(H) + \mqty(0 & -1 \\ 1 & 0) \otimes \Im(H)
```

See also [`Isomorphisms.G`](@ref), [`Isomorphisms.H`](@ref).
"""
iso(H::AbstractMatrix{<:Number}) = kron(I(2), real(H)) + kron(Im2, imag(H))

@doc raw"""
    G(H::AbstractMatrix)::Matrix{Float64}

Returns the isomorphism of ``-iH``, i.e. ``G(H) = \text{iso}(-iH)``.

See also [`Isomorphisms.iso`](@ref), [`Isomorphisms.H`](@ref).
"""
G(H::AbstractMatrix{<:Number}) = iso(-im * H)

@doc raw"""
    H(G::AbstractMatrix{<:Real})

Returns the inverse of ``G(H) = iso(-iH)``, i.e. returns H.

See also [`Isomorphisms.iso`](@ref), [`Isomorphisms.G`](@ref).
"""
function H(G::AbstractMatrix{<:Real})
    dim = size(G, 1) ÷ 2
    H_imag = G[1:dim, 1:dim]
    H_real = -G[(dim+1):end, 1:dim]
    return H_real + 1.0im * H_imag
end

@doc raw"""
    ad_vec(H::AbstractMatrix{ℂ}; anti::Bool=false) where ℂ <: Number

Returns the vectorized adjoint action of a matrix `H`:

```math
\text{ad_vec}(H) = \mqty(1 & 0 \\ 0 & 1) \otimes H - (-1)^{\text{anti}} \mqty(0 & 1 \\ 1 & 0) \otimes H^*
```
"""
function ad_vec(H::AbstractMatrix{ℂ}; anti::Bool = false) where {ℂ<:Number}
    Id = sparse(ℂ, I, size(H)...)
    return kron(Id, H) - (-1)^anti * kron(conj(H)', Id)
end

@doc raw"""
    iso_D(L::AbstractMatrix{ℂ}) where ℂ <: Number

Returns the isomorphic representation of the Lindblad dissipator `L`.
"""
function iso_D(L::AbstractMatrix{ℂ}) where {ℂ<:Number}
    return iso(kron(conj(L), L) - 1 / 2 * ad_vec(L'L, anti = true))
end

@doc raw"""
    var_G(G::AbstractMatrix{<:Real}, G_vars::AbstractVector{<:AbstractMatrix{<:Real}})

Returns the variational generator of `G` with variational derivatives, `G_vars`.

The variational generator is 
```math
\text{var}_G(G, [G_a, G_b]) = \mqty( G & 0 & 0 \\ G_a & G & 0 \\ G_b & 0 & G )
```
where `G` is the isomorphism of a Hamiltonian and `G_a` and `G_b` are the variational 
derivatives of `G` for parameters `a` and `b`, respectively.
"""
function var_G(
    G::AbstractMatrix{ℝ1},
    G_vars::AbstractVector{<:AbstractMatrix{ℝ2}},
) where {ℝ1<:Real,ℝ2<:Real}
    n, m = size(G)
    v = length(G_vars)
    G_0 = kron(I(v + 1), G)
    G_V = spzeros(ℝ2, (v + 1) * n, (v + 1) * m)
    for i in eachindex(G_vars)
        G_V[(i*n+1):((i+1)*n), 1:m] = G_vars[i]
    end
    return G_0 + G_V
end

# ----------------------------------------------------------------------------- #
#                             Bloch Sphere                                      #
# ----------------------------------------------------------------------------- #

@doc raw"""
    ket_to_bloch(ψ::AbstractVector{<:Number})

Convert a ket to a Bloch vector representation.
"""
function ket_to_bloch(ψ::AbstractVector{<:Number})
    @assert length(ψ) == 2
    ρ = ψ * ψ'
    bloch_vector = [real(tr(ρ * P)) for P in [PAULIS.X, PAULIS.Y, PAULIS.Z]]

    return bloch_vector / norm(bloch_vector)
end

@doc raw"""
    bloch_to_ket(v::AbstractVector{<:Real}; digits=6)


Convert a Bloch vector to a ket (up to global phase).
"""
function bloch_to_ket(bloch::AbstractVector{R}; digits::Integer = 6) where {R<:Real}
    @assert length(bloch) == 3
    x, y, z = bloch

    θ = acos(z)
    φ = atan(y, x)

    return Complex{R}[cos(θ / 2), exp(im * φ)*sin(θ / 2)]

end

function density_to_bloch(ψ::AbstractVector{<:Number})
    @assert length(ψ) == 2
    ρ = ψ * ψ'
    return [real(tr(ρ * P)) for P in [PAULIS.X, PAULIS.Y, PAULIS.Z]]
end

function bloch_to_density(v::AbstractVector{<:Real})
    @assert length(v) == 3
    return 0.5 * (I(2) + v[1] * PAULIS.X + v[2] * PAULIS.Y + v[3] * PAULIS.Z)
end

# *************************************************************************** #

@testitem "Test ket isomorphisms" begin
    @test ket_to_iso([1.0, 2.0]) ≈ [1.0, 2.0, 0.0, 0.0]
    @test ket_to_iso([-im, 2.0 + 3.0im]) ≈ [0.0, 2.0, -1.0, 3.0]
    @test iso_to_ket([1.0, 2.0, 0.0, 0.0]) ≈ [1.0, 2.0]
    @test iso_to_ket([0.0, 2.0, -1.0, 3.0]) ≈ [-im, 2.0 + 3.0im]
end

@testitem "Test operator isomorphisms" begin
    iso_vec_I = [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    @test mat([1.0, 2.0, 3.0, 4.0]) ≈ [
        1.0 3.0
        2.0 4.0
    ]
    @test iso_vec_to_operator(iso_vec_I) ≈ [
        1.0 0.0
        0.0 1.0
    ]
    @test iso_vec_to_iso_operator(iso_vec_I) ≈ [
        1.0 0.0 0.0 0.0
        0.0 1.0 0.0 0.0
        0.0 0.0 1.0 0.0
        0.0 0.0 0.0 1.0
    ]
    @test operator_to_iso_vec(Complex[1.0 0.0; 0.0 1.0]) ≈ iso_vec_I
    @test iso_operator_to_iso_vec(iso_vec_to_iso_operator(iso_vec_I)) ≈ iso_vec_I

    iso_vec_XY = [0, 1, 0, 1, 1, 0, -1, 0]
    @test iso_vec_to_operator(iso_vec_XY) ≈ [
        0 1-im
        1+im 0
    ]
    @test iso_vec_to_iso_operator(iso_vec_XY) ≈ [
        0 1 0 1
        1 0 -1 0
        0 -1 0 1
        1 0 1 0
    ]
    @test operator_to_iso_vec(Complex[0.0 1-im; 1+im 0.0]) ≈ iso_vec_XY
    @test iso_operator_to_iso_vec(iso_vec_to_iso_operator(iso_vec_XY)) ≈ iso_vec_XY
end

@testitem "Test density matrix isomorphisms" begin
    # Totally mixed state
    ρ = [1.0 0.0; 0.0 1.0]
    ρ_iso = [1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
    @test density_to_iso_vec(ρ) ≈ ρ_iso
    @test iso_vec_to_density(ρ_iso) ≈ ρ

    # Density matrix of a Bell state
    ρ = [1.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 1.0 0.0 0.0 1.0] / 2
    @test iso_vec_to_density(density_to_iso_vec(ρ)) ≈ ρ

    # Random
    ρ1 = [1.0 1.0; 1.0 1.0] / 2
    U1 = [
        -0.831976-0.101652im -0.422559-0.344857im
        -0.527557+0.138444im 0.799158+0.252713im
    ]
    ρ2 = [1.0 0.0; 0.0 0.0]
    U2 = [
        -0.784966-0.163279im -0.597246-0.0215881im
        0.597536+0.0109124im -0.792681+0.120364im
    ]
    ρ = (U1 * ρ1 * U1' + U2 * ρ2 * U2') / 2
    @test iso_vec_to_density(density_to_iso_vec(ρ)) ≈ ρ
    @test iso_vec_to_density(density_to_iso_vec(ρ)) ≈ ρ
end

@testitem "Test compact density isomorphisms" begin
    using LinearAlgebra
    using SparseArrays

    # --- Round-trip for 2×2 ---
    ρ2 = ComplexF64[0.7 0.3-0.1im; 0.3+0.1im 0.3]
    x2 = density_to_compact_iso(ρ2)
    @test length(x2) == 4
    @test compact_iso_to_density(x2) ≈ ρ2

    # --- Round-trip for 3×3 ---
    A = randn(ComplexF64, 3, 3)
    ρ3 = (A + A') / 2
    x3 = density_to_compact_iso(ρ3)
    @test length(x3) == 9
    @test compact_iso_to_density(x3) ≈ ρ3

    # --- Round-trip for 4×4 ---
    A = randn(ComplexF64, 4, 4)
    ρ4 = (A + A') / 2
    x4 = density_to_compact_iso(ρ4)
    @test length(x4) == 16
    @test compact_iso_to_density(x4) ≈ ρ4

    # --- L and P dimensions ---
    for n in [2, 3, 4, 7]
        L = density_lift_matrix(n)
        P = density_projection_matrix(n)

        @test size(L) == (2n^2, n^2)
        @test size(P) == (n^2, 2n^2)

        # nnz counts
        @test nnz(L) == n * (2n - 1)
        @test nnz(P) == n^2

        # P * L = I
        @test sparse(P * L) ≈ sparse(I, n^2, n^2)
    end

    # --- L * P projects onto Hermitian subspace ---
    for n in [2, 3, 4]
        L = density_lift_matrix(n)
        P = density_projection_matrix(n)
        LP = L * P

        # Idempotent
        @test LP * LP ≈ LP

        # Preserves valid Hermitian iso_vecs
        A = randn(ComplexF64, n, n)
        ρ = (A + A') / 2
        v = density_to_iso_vec(ρ)
        @test LP * v ≈ v
    end

    # --- Compatibility: iso_vec = L * compact ---
    for n in [2, 3, 5]
        L = density_lift_matrix(n)
        P = density_projection_matrix(n)

        A = randn(ComplexF64, n, n)
        ρ = (A + A') / 2

        v_full = density_to_iso_vec(ρ)
        v_compact = density_to_compact_iso(ρ)

        @test L * v_compact ≈ v_full
        @test P * v_full ≈ v_compact
    end

    # --- Explicit 2×2 ordering check ---
    # ρ = [a  c+di; c-di  b] with a,b,c,d real (Hermitian: ρ[1,2] = c+di, ρ[2,1] = c-di)
    a, b, c, d = 0.6, 0.4, 0.2, 0.1
    ρ = ComplexF64[a c+d*im; c-d*im b]
    x = density_to_compact_iso(ρ)
    # Re upper triangle (col-major): Re(ρ[1,1]), Re(ρ[1,2]), Re(ρ[2,2])
    # Im strict upper triangle (col-major): Im(ρ[1,2])
    @test x ≈ [a, c, b, d]
end

@testitem "Test Hamiltonian isomorphisms" begin
    using .Isomorphisms

    H_real = [1.0 2.0; 3.0 4.0]
    H_imag = [0.0 1.0; 1.0 0.0]
    H_complex = H_real + 1.0im * H_imag
    G_H = Isomorphisms.G(H_complex)

    @test Isomorphisms.H(G_H) ≈ H_complex

    @test G_H ≈ [0 1 1 2; 1 0 3 4; -1 -2 0 1; -3 -4 1 0]

    @test Isomorphisms.iso(H_complex) ≈ [1 2 0 -1; 3 4 -1 0; 0 1 1 2; 1 0 3 4]

    @test Isomorphisms.iso(-im * H_complex) ≈ G_H

    op = [0 1; 1 0]
    ad_H = Isomorphisms.ad_vec(op)
    @test ad_H ≈ [0 1 -1 0; 1 0 0 -1; -1 0 0 1; 0 -1 1 0]

    op = [0 -im; im 0]
    ad_H = Isomorphisms.ad_vec(op)
    @test ad_H ≈ [0 -im -im 0; im 0 0 -im; im 0 0 -im; 0 im im 0]
end

@testitem "Test variational G isomorphism" begin
    using .Isomorphisms

    G = [1.0 2.0; 3.0 4.0]
    G_var1 = [0.0 1.0; 1.0 0.0]
    G_var2 = [0.0 0.0; 1.0 1.0]

    G_vars = [G_var1]
    Ĝ = Isomorphisms.var_G(G, G_vars)
    @test Ĝ ≈ [
        1.0 2.0 0.0 0.0
        3.0 4.0 0.0 0.0
        0.0 1.0 1.0 2.0
        1.0 0.0 3.0 4.0
    ]

    G_vars = [G_var1, G_var2]
    Ĝ = Isomorphisms.var_G(G, G_vars)
    @test Ĝ ≈ [
        1.0 2.0 0.0 0.0 0.0 0.0
        3.0 4.0 0.0 0.0 0.0 0.0
        0.0 1.0 1.0 2.0 0.0 0.0
        1.0 0.0 3.0 4.0 0.0 0.0
        0.0 0.0 0.0 0.0 1.0 2.0
        1.0 1.0 0.0 0.0 3.0 4.0
    ]
end

@testitem "Test Bloch vector to ket and ket to Bloch vector" begin
    using LinearAlgebra: dot

    ψ₁ = [1.0, 0.0]
    ψ₂ = [0.0, 1.0]
    ψ₃ = [1 / √2, 1 / √2]
    ψ₄ = [1 / √2, -1im / √2]

    for ψ in (ψ₁, ψ₂, ψ₃, ψ₄)
        bloch = ket_to_bloch(ψ)
        ψ′ = bloch_to_ket(bloch)
        @test abs2(dot(ψ′, ψ)) ≈ 1.0 atol = 1e-10
    end
end

end
