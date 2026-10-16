"""
    gauss_legendre_01(n)

Return Gauss-Legendre quadrature nodes and weights on [0, 1].
Uses the Golub-Welsch algorithm (eigenvalues of the symmetric tridiagonal matrix).
"""
function gauss_legendre_01(n::Int)
    β = [i / sqrt(4i^2 - 1) for i = 1:(n-1)]
    T = SymTridiagonal(zeros(n), β)
    nodes, V = eigen(T)
    weights = 2 .* V[1, :] .^ 2
    # Transform [-1, 1] → [0, 1]
    return (nodes .+ 1) ./ 2, weights ./ 2
end

"""
    compute_sensitivity_kick_first_order(H_drive, ψ_next, Δt)

First-order sensitivity kick: `-iΔt · H_drive · ψ_{j+1}`.

This is the leading-order approximation to ∂(U_j ψ_j)/∂u_k,
valid when Δt is small.
"""
function compute_sensitivity_kick_first_order(
    H_drive::AbstractMatrix,
    ψ_next::AbstractVector,
    Δt::Real,
)
    return (-im * Δt) .* (H_drive * ψ_next)
end

"""
    compute_sensitivity_kick_exact(H_drive, ψ, Δt, λ, V; n_quad=4)

Exact sensitivity kick via Gauss-Legendre quadrature over the eigendecomposition.

Computes `∂(exp(-iΔtH)ψ)/∂u_k` using:

    -iΔt ∫₀¹ exp(-iΔt(1-s)H) · H_drive · exp(-iΔt·s·H) · ψ  ds

where `H = V · Diagonal(λ) · V'` is the cached eigendecomposition of the total
Hamiltonian at this knot point.
"""
function compute_sensitivity_kick_exact(
    H_drive::AbstractMatrix,
    ψ::AbstractVector,
    Δt::Real,
    λ::AbstractVector,
    V::AbstractMatrix;
    n_quad::Int = 4,
)
    nodes, weights = gauss_legendre_01(n_quad)
    n = length(ψ)
    kick = zeros(ComplexF64, n)

    # Precompute transformed quantities
    VtHkV = V' * H_drive * V
    Vtψ = V' * ψ

    for (s, w) in zip(nodes, weights)
        # exp(-iΔt·s·H) · ψ  in eigenbasis
        Ds_Vtψ = cis.(-Δt .* s .* λ) .* Vtψ
        # H_drive · exp(-iΔt·s·H) · ψ  in eigenbasis
        mid = VtHkV * Ds_Vtψ
        # exp(-iΔt·(1-s)·H) · H_drive · exp(-iΔt·s·H) · ψ  in eigenbasis
        result = cis.(-Δt .* (1 - s) .* λ) .* mid
        kick .+= w .* result
    end

    return (-im * Δt) .* (V * kick)
end
