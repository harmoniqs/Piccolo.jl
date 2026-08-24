# ============================================================================ #
# Daleckii–Krein (divided-difference) first-order derivative kernel — Layer A
#
# Trajectory-agnostic math for the analytic derivative of the Hermitian
# exponential step Φ = exp(-iΔt H) with respect to an affine-in-parameter drive.
# This is the shared kernel every DK slice consumes (Piccolissimo.jl#199); it does
# NOT wire into the integrators (later slices). No finite differences and no
# autodiff live here — analytic only. All apply math is array-type-generic
# (no hardcoded `Vector`/`Matrix{Float64}`) so a device block is a drop-in
# follow-on.
#
# Convention (pinned by #199):
#   ∂Φ/∂aⱼ = V (M ∘ (V† Hⱼ V)) V†
# where H = V diag(λ) V† is the eigendecomposition the forward step already
# holds, Hⱼ = ∂H/∂aⱼ is the CONSTANT drive matrix (NO -iΔt factor on the
# direction), and M is the divided-difference matrix of φ(λ)=exp(-iΔt λ):
#   off-diagonal   M_ab = (φ_a - φ_b)/(λ_a - λ_b)
#   confluent diag M_aa = φ'(λ_a) = -iΔt·exp(-iΔt λ_a)
# The -iΔt lives entirely in φ/M; putting it on the direction too double-counts.
# ============================================================================ #

"""
    DK_DEGENERACY_RTOL

Relative, Δt-aware switchover tolerance for the divided-difference degeneracy
band (pinned in Piccolissimo.jl#199 slice ①). An entry `(a,b)` uses the confluent
(derivative) limit instead of the raw divided difference when the dimensionless
phase gap satisfies

    |Δt · (λ_a - λ_b)| < DK_DEGENERACY_RTOL.

`Δt·(λ_a - λ_b)` is the argument gap of `φ(λ)=exp(-iΔt λ)`, so the band is
naturally relative (dimensionless) and Δt-aware.

# Why this value
The raw divided difference `(φ_a - φ_b)/(λ_a - λ_b)` suffers catastrophic
cancellation for a small phase gap `θ = Δt(λ_a-λ_b)`: the numerator has true
magnitude `≈|θ|` but absolute rounding error `≈ε` (the operands have unit
modulus), giving a relative error `≈ ε/|θ|`. Just OUTSIDE the band (still using
the divided difference) that error is `≤ ε/DK_DEGENERACY_RTOL ≈ 2·10⁻¹⁰`, below
the 1e-9 correctness gate. INSIDE the band we substitute the confluent value at
the midpoint `m=(λ_a+λ_b)/2`, whose approximation error to the true divided
difference is `≈ θ²/24` — at the boundary `≈ 4·10⁻¹⁴`. So both branches stay
well under 1e-9 and the switchover is continuous (no jump).
"""
const DK_DEGENERACY_RTOL = 1e-6

"""
    DaleckiiKreinWorkspace(n; T=ComplexF64)
    DaleckiiKreinWorkspace(M, tmp1, tmp2)

Per-thread scratch for the first-order Daleckii–Krein kernel: the
divided-difference matrix `M` and the two conjugation temporaries used by
[`dk_apply!`](@ref), each sized `(n, n)`. Sized to ketdim and intended to be
allocated once alongside the `exp_eigen!` buffer set (later slices thread one
workspace per OS thread; this slice is pure kernel).

The struct is array-type-generic: pass device matrices to the `(M, tmp1, tmp2)`
constructor and the kernel runs unchanged. The `n` convenience constructor
allocates CPU `Matrix{T}` scratch.
"""
struct DaleckiiKreinWorkspace{T<:Number,MT<:AbstractMatrix{T}}
    M::MT      # divided-difference matrix   (n, n)
    tmp1::MT   # conjugation temporary       (n, n)
    tmp2::MT   # conjugation temporary       (n, n)
end

# The 3-arg `DaleckiiKreinWorkspace(M, tmp1, tmp2)` form is the struct's
# auto-generated constructor — do NOT add an explicit one with the same
# signature: it overwrites the auto-generated method, which errors during
# module precompilation ("Method overwriting is not permitted...").

function DaleckiiKreinWorkspace(n::Integer; T::Type{<:Number} = ComplexF64)
    return DaleckiiKreinWorkspace(zeros(T, n, n), zeros(T, n, n), zeros(T, n, n))
end

"""
    dk_divided_difference!(M, λ, Δt; degeneracy_rtol=DK_DEGENERACY_RTOL)

Fill `M` (an `n×n` complex matrix) in place with the divided-difference matrix
of `φ(λ)=exp(-iΔt λ)` over the eigenvalues `λ`:

- off-diagonal, well-separated: `M[a,b] = (φ_a - φ_b)/(λ_a - λ_b)`
- confluent (diagonal `a==b`, or `|Δt·(λ_a-λ_b)| < degeneracy_rtol`):
  `M[a,b] = -iΔt·exp(-iΔt·m)`, the derivative `φ'` evaluated at the midpoint
  `m = (λ_a + λ_b)/2`. For `a==b` this is exactly `-iΔt·exp(-iΔt λ_a)`.

The midpoint confluent value makes the switchover continuous (error `≈ θ²/24`
to the true divided difference, `θ = Δt(λ_a-λ_b)`), so the sweep `λ_b → λ_a`
has no jump and no `NaN`. Analytic only — no finite differences, no autodiff.
Array-type-generic and zero-allocation after the caller preallocates `M`.
"""
function dk_divided_difference!(
    M::AbstractMatrix,
    λ::AbstractVector,
    Δt::Real;
    degeneracy_rtol::Real = DK_DEGENERACY_RTOL,
)
    n = length(λ)
    @boundscheck size(M) == (n, n) ||
                 throw(DimensionMismatch("M must be $(n)×$(n) for length-$(n) λ"))
    @inbounds for b = 1:n, a = 1:n
        λa = λ[a]
        λb = λ[b]
        d = λa - λb
        if abs(Δt * d) < degeneracy_rtol
            # Confluent (derivative) limit at the midpoint — continuous across
            # the band and equal to φ'(λ_a) on the exact diagonal.
            m = (λa + λb) / 2
            M[a, b] = -im * Δt * cis(-Δt * m)
        else
            M[a, b] = (cis(-Δt * λa) - cis(-Δt * λb)) / d
        end
    end
    return M
end

"""
    dk_apply!(∂Φ, V, M, Hⱼ, tmp1, tmp2)

Apply the first-order Daleckii–Krein action

    ∂Φ = V (M ∘ (V† Hⱼ V)) V†

into `∂Φ`, where `V` are the eigenvectors of `H` (from the forward step), `M` is
the divided-difference matrix (see [`dk_divided_difference!`](@ref)), `Hⱼ` is the
constant drive direction `∂H/∂aⱼ`, and `tmp1`/`tmp2` are `n×n` scratch matrices.
`∘` is the elementwise (Hadamard) product. Zero-allocation after warmup and
array-type-generic (only `mul!` + broadcast on `AbstractMatrix`).
"""
function dk_apply!(
    ∂Φ::AbstractMatrix,
    V::AbstractMatrix,
    M::AbstractMatrix,
    Hⱼ::AbstractMatrix,
    tmp1::AbstractMatrix,
    tmp2::AbstractMatrix,
)
    mul!(tmp1, V', Hⱼ)      # tmp1 = V† Hⱼ
    mul!(tmp2, tmp1, V)     # tmp2 = V† Hⱼ V
    tmp2 .*= M              # tmp2 = M ∘ (V† Hⱼ V)
    mul!(tmp1, V, tmp2)     # tmp1 = V (M ∘ (V† Hⱼ V))
    mul!(∂Φ, tmp1, V')      # ∂Φ   = V (M ∘ (V† Hⱼ V)) V†
    return ∂Φ
end

"""
    dk_first_order_derivative!(∂Φ, ws, V, λ, Hⱼ, Δt; degeneracy_rtol=DK_DEGENERACY_RTOL)

Convenience driver: build the divided-difference matrix into `ws.M` and apply the
first-order action into `∂Φ`, reusing the eigenbasis `(V, λ)` the forward step
already produced. `ws` is a [`DaleckiiKreinWorkspace`](@ref). Zero-allocation
after warmup.
"""
function dk_first_order_derivative!(
    ∂Φ::AbstractMatrix,
    ws::DaleckiiKreinWorkspace,
    V::AbstractMatrix,
    λ::AbstractVector,
    Hⱼ::AbstractMatrix,
    Δt::Real;
    degeneracy_rtol::Real = DK_DEGENERACY_RTOL,
)
    dk_divided_difference!(ws.M, λ, Δt; degeneracy_rtol = degeneracy_rtol)
    dk_apply!(∂Φ, V, ws.M, Hⱼ, ws.tmp1, ws.tmp2)
    return ∂Φ
end

# ============================================================================ #
# Second-order Daleckii–Krein kernel — Layer A (Piccolissimo.jl#201, part of #199)
#
# The second Fréchet derivative of Φ = exp(-iΔt H) is a genuine THREE-index
# second-divided-difference contraction (NOT a Hadamard mask — that is the
# first-order form). With H = V diag(λ) V†, drive directions E₁,E₂ = ∂H/∂aᵢ,
# ∂H/∂aⱼ (constant, no -iΔt factor), and Ẽₖ = V† Eₖ V:
#
#   D²Φ[E₁,E₂] = V C V†,
#   C[a,c] = Σ_b φ⁽²⁾(λ_a, λ_b, λ_c) (Ẽ₁[a,b]Ẽ₂[b,c] + Ẽ₂[a,b]Ẽ₁[b,c])
#
# summed over the intermediate eigen-index b. φ⁽²⁾ is the second divided
# difference of φ(λ)=exp(-iΔt λ) in the STANDARD (1/k!) convention, so on the
# full diagonal φ⁽²⁾(x,x,x) = φ″(x)/2. The symmetrization over (E₁,E₂) makes the
# block D²Φ[Eᵢ,Eⱼ] symmetric in (i,j) by construction; for a pure second
# derivative ∂²Φ/∂a² take E₁=E₂ and D²Φ[E,E] = ∂²Φ/∂a² directly.
#
# The confluent limits (any two — or all three — of λ_a,λ_b,λ_c collide) are the
# delicate part: a naive divided difference cancels catastrophically. We keep it
# finite and accurate by (i) dividing by the LARGEST of the three pairwise gaps
# (the outer pair), which minimises cancellation and stays exact unless all
# three nodes are clustered, and (ii) an all-confluent Taylor branch (a short
# convergent series in the centroid-shifted spread) when even the largest gap is
# inside the degeneracy band. Reuses slice ①'s `DK_DEGENERACY_RTOL` philosophy.
#
# This is the pure Layer-A kernel + tests; no integrator wiring (slice ④), no
# matrix-free HVP (slice ⑥). Analytic only — FD/ForwardDiff live in the tests as
# witnesses. All apply math is array-type-generic (only `mul!` + scalar fills,
# matching slice ①), so a device block is a drop-in follow-on.
# ============================================================================ #

# Stable, cancellation-free first divided difference of φ(λ)=exp(-iΔt λ):
#   φ¹(p,q) = -iΔt·exp(-iΔt·m)·sinc₀(Δt·δ),  m=(p+q)/2, δ=(p-q)/2, sinc₀(t)=sin t/t.
# Uniformly accurate for all p,q (equal to φ′(p) at p==q); used for the inner
# divided differences of the second-order recursion so only the OUTER division
# can lose digits — and that is guarded by choosing the largest gap.
_dk_sinc0(x) = iszero(x) ? one(x) : sin(x) / x

function _dk_phi1(p, q, Δt::Real)
    m = (p + q) / 2
    δ = (p - q) / 2
    return -im * Δt * cis(-Δt * m) * _dk_sinc0(Δt * δ)
end

# All-confluent branch: nodes clustered inside the degeneracy band. Shift to the
# centroid m (so the shifted nodes uᵢ sum to zero) and expand
#   φ⁽²⁾(x,y,z) = e^{-iΔt·m} · Σ_{k≥0} c^{k+2}/(k+2)! · h_k(u),   c = -iΔt,
# where h_k is the complete homogeneous symmetric polynomial of the uᵢ. Because
# Σuᵢ = 0 we have h₁ = 0, and the uᵢ are O(band), so a handful of terms is exact
# to machine precision (term ratio ≈ (Δt·gap)² ≲ 1e-12 at the band edge). Newton's
# identity gives h_k from the power sums pₖ = Σuᵢᵏ. Allocation-free (scalars only).
function _dk_phi2_confluent(λa, λb, λc, Δt::Real)
    c = -im * Δt
    m = (λa + λb + λc) / 3
    u1 = λa - m
    u2 = λb - m
    u3 = λc - m
    p2 = u1^2 + u2^2 + u3^2
    p3 = u1^3 + u2^3 + u3^3
    p4 = u1^4 + u2^4 + u3^4
    h2 = p2 / 2                 # h₁ = p₁ = 0 (centroid shift)
    h3 = p3 / 3
    h4 = (p2^2 / 2 + p4) / 4
    S = c^2 / 2 + c^4 / 24 * h2 + c^5 / 120 * h3 + c^6 / 720 * h4
    return cis(-Δt * m) * S
end

"""
    dk_second_divided_difference(λa, λb, λc, Δt; degeneracy_rtol=DK_DEGENERACY_RTOL)

Second divided difference `φ⁽²⁾(λ_a, λ_b, λ_c)` of `φ(λ)=exp(-iΔt λ)` in the
standard (`1/k!`) convention, so `φ⁽²⁾(x,x,x) = φ″(x)/2 = -Δt²·exp(-iΔt x)/2`.
Symmetric in its three arguments.

Evaluated so it stays finite and accurate through every confluent limit:

- **Well-separated** (largest pairwise gap outside the band): the recursion
  `(φ¹(x,y) - φ¹(y,z))/(x - z)` with the outer pair `(x,z)` chosen as the pair
  with the LARGEST separation (`|x-z|` = max gap). The inner `φ¹` use the stable
  `sinc₀` form (`_dk_phi1`), so the only place digits can cancel is the
  outer division — and dividing by the largest gap makes that as well-conditioned
  as the spectrum allows.
- **All clustered** (even the largest gap `|Δt·gap| < degeneracy_rtol`): a short
  convergent Taylor series in the centroid-shifted node spread — finite,
  continuous, no cancellation.

Analytic only (no finite differences / autodiff); allocation-free.
"""
function dk_second_divided_difference(
    λa,
    λb,
    λc,
    Δt::Real;
    degeneracy_rtol::Real = DK_DEGENERACY_RTOL,
)
    gab = abs(λa - λb)
    gbc = abs(λb - λc)
    gac = abs(λa - λc)
    gmax = max(gab, gbc, gac)
    if abs(Δt) * gmax < degeneracy_rtol
        return _dk_phi2_confluent(λa, λb, λc, Δt)
    end
    # Reorder so the outer division uses the largest gap (minimal cancellation):
    #   φ⁽²⁾(x,y,z) = (φ¹(x,y) - φ¹(y,z)) / (x - z),  |x-z| = gmax, y = middle node.
    if gac ≥ gab && gac ≥ gbc
        x, y, z = λa, λb, λc        # outer pair (a,c)
    elseif gab ≥ gbc
        x, y, z = λa, λc, λb        # outer pair (a,b)
    else
        x, y, z = λb, λa, λc        # outer pair (b,c)
    end
    return (_dk_phi1(x, y, Δt) - _dk_phi1(y, z, Δt)) / (x - z)
end

"""
    DaleckiiKreinSecondOrderWorkspace(n; T=ComplexF64)
    DaleckiiKreinSecondOrderWorkspace(E1t, E2t, core, tmp, D2)

Per-thread scratch for the second-order Daleckii–Krein kernel: the two rotated
directions `Ẽₖ = V† Eₖ V` (`E1t`, `E2t`), the middle contraction tensor `C`
(`core`), a conjugation temporary `tmp`, and an output buffer `D2` used by
`dk_second_order_block!` — each `(n, n)`. Sized to ketdim and intended to
be allocated once alongside the `exp_eigen!` buffer set (later slices thread one
workspace per OS thread; this slice is pure kernel).

Array-type-generic: pass device matrices to the 5-arg constructor and the kernel
runs unchanged; the `n` convenience constructor allocates CPU `Matrix{T}` scratch.
"""
struct DaleckiiKreinSecondOrderWorkspace{T<:Number,MT<:AbstractMatrix{T}}
    E1t::MT    # Ẽ₁ = V† E₁ V              (n, n)
    E2t::MT    # Ẽ₂ = V† E₂ V              (n, n)
    core::MT   # middle contraction tensor C (n, n)
    tmp::MT    # conjugation temporary      (n, n)
    D2::MT     # output buffer for block!   (n, n)
end

# The 5-arg form is the struct's auto-generated constructor — do NOT add an
# explicit one with the same signature (it would overwrite the auto-generated
# method and break module precompilation, cf. slice ①'s hard-won lesson).
function DaleckiiKreinSecondOrderWorkspace(n::Integer; T::Type{<:Number} = ComplexF64)
    return DaleckiiKreinSecondOrderWorkspace(
        zeros(T, n, n),
        zeros(T, n, n),
        zeros(T, n, n),
        zeros(T, n, n),
        zeros(T, n, n),
    )
end

"""
    dk_second_order_apply!(∂²Φ, V, λ, E1, E2, Δt, E1t, E2t, core, tmp; degeneracy_rtol=DK_DEGENERACY_RTOL)

Apply the second Fréchet derivative of `Φ = exp(-iΔt H)`,

    ∂²Φ = D²Φ[E₁,E₂] = V C V†,
    C[a,c] = Σ_b φ⁽²⁾(λ_a,λ_b,λ_c)·(Ẽ₁[a,b]Ẽ₂[b,c] + Ẽ₂[a,b]Ẽ₁[b,c]),

into `∂²Φ`, reusing the eigenbasis `(V, λ)` the forward step already produced.
`E1`,`E2` are the constant drive directions `∂H/∂aᵢ`, `∂H/∂aⱼ`; `E1t`,`E2t`,
`core`,`tmp` are `n×n` scratch (`Ẽₖ = V† Eₖ V`, the contraction tensor, and a
conjugation temporary). For a pure second derivative `∂²Φ/∂a²`, pass `E1 == E2`.

Cost is `O(n³)` (the three-index contraction), materially heavier than the
first-order Hadamard — the scratch and cost are sized accordingly (`n` = ketdim,
tens). Zero-allocation after warmup; array-type-generic (only `mul!` + scalar
fills of `core`, mirroring [`dk_divided_difference!`](@ref)).
"""
function dk_second_order_apply!(
    ∂²Φ::AbstractMatrix,
    V::AbstractMatrix,
    λ::AbstractVector,
    E1::AbstractMatrix,
    E2::AbstractMatrix,
    Δt::Real,
    E1t::AbstractMatrix,
    E2t::AbstractMatrix,
    core::AbstractMatrix,
    tmp::AbstractMatrix;
    degeneracy_rtol::Real = DK_DEGENERACY_RTOL,
)
    n = length(λ)
    # Rotate the drive directions into the eigenbasis: Ẽₖ = V† Eₖ V.
    mul!(tmp, V', E1)
    mul!(E1t, tmp, V)
    mul!(tmp, V', E2)
    mul!(E2t, tmp, V)
    # Three-index contraction over the intermediate eigen-index b.
    @inbounds for c = 1:n, a = 1:n
        λa = λ[a]
        λc = λ[c]
        s = zero(eltype(core))
        for b = 1:n
            φ2 = dk_second_divided_difference(
                λa,
                λ[b],
                λc,
                Δt;
                degeneracy_rtol = degeneracy_rtol,
            )
            s += φ2 * (E1t[a, b] * E2t[b, c] + E2t[a, b] * E1t[b, c])
        end
        core[a, c] = s
    end
    # ∂²Φ = V C V†
    mul!(tmp, V, core)
    mul!(∂²Φ, tmp, V')
    return ∂²Φ
end

"""
    dk_second_order_derivative!(∂²Φ, ws, V, λ, E1, E2, Δt; degeneracy_rtol=DK_DEGENERACY_RTOL)

Convenience driver: apply the second-order action `D²Φ[E₁,E₂]` into `∂²Φ` using
the scratch on `ws` (a [`DaleckiiKreinSecondOrderWorkspace`](@ref)), reusing the
eigenbasis `(V, λ)`. Zero-allocation after warmup.
"""
function dk_second_order_derivative!(
    ∂²Φ::AbstractMatrix,
    ws::DaleckiiKreinSecondOrderWorkspace,
    V::AbstractMatrix,
    λ::AbstractVector,
    E1::AbstractMatrix,
    E2::AbstractMatrix,
    Δt::Real;
    degeneracy_rtol::Real = DK_DEGENERACY_RTOL,
)
    dk_second_order_apply!(
        ∂²Φ,
        V,
        λ,
        E1,
        E2,
        Δt,
        ws.E1t,
        ws.E2t,
        ws.core,
        ws.tmp;
        degeneracy_rtol = degeneracy_rtol,
    )
    return ∂²Φ
end

"""
    dk_second_order_block!(B, ws, V, λ, Hs, μ, Δt; degeneracy_rtol=DK_DEGENERACY_RTOL)

Assemble the μ-contracted parameter-parameter Hessian block

    B[i,j] = ⟨μ, D²Φ[Hᵢ,Hⱼ]⟩ = Re tr(μ† D²Φ[Hᵢ,Hⱼ]),

the `gauss_newton=false` block consumed by the exact-Hessian path (slice ⑥).
`Hs` is an indexable collection of drive directions `∂H/∂aⱼ` and `μ` is the
`n×n` multiplier. The block is symmetric to machine precision by construction
(the second Fréchet derivative is symmetrized over its two directions, and each
`B[i,j]` / `B[j,i]` pair is written from one contraction). Uses `ws.D2` as the
per-pair output buffer; `B` is `length(Hs) × length(Hs)` real.
"""
function dk_second_order_block!(
    B::AbstractMatrix,
    ws::DaleckiiKreinSecondOrderWorkspace,
    V::AbstractMatrix,
    λ::AbstractVector,
    Hs,
    μ::AbstractMatrix,
    Δt::Real;
    degeneracy_rtol::Real = DK_DEGENERACY_RTOL,
)
    m = length(Hs)
    @boundscheck size(B) == (m, m) ||
                 throw(DimensionMismatch("B must be $(m)×$(m) for $(m) directions"))
    for j = 1:m, i = 1:j
        dk_second_order_apply!(
            ws.D2,
            V,
            λ,
            Hs[i],
            Hs[j],
            Δt,
            ws.E1t,
            ws.E2t,
            ws.core,
            ws.tmp;
            degeneracy_rtol = degeneracy_rtol,
        )
        val = real(dot(μ, ws.D2))
        B[i, j] = val
        B[j, i] = val
    end
    return B
end

# ============================================================================ #
# Tests (AC1–AC5 of Piccolissimo.jl#200). Random-Hermitian fixtures mirror the
# `exp_eigen!` testitems above (Random.seed! + randn ComplexF64, Hermitian(A+A')).
# ============================================================================ #

@testitem "DK divided-difference M matches naive elementwise reference (AC1)" begin
    using LinearAlgebra
    using Random
    using ..QuantumIntegrators: dk_divided_difference!

    Random.seed!(42)
    for _ = 1:25
        n = rand(3:8)
        λ = sort(randn(n) .* 3)          # random, well-separated spectrum
        Δt = 0.3 + rand()
        # guard the fixture away from the degeneracy band (naive ref is only
        # valid where the raw divided difference is used)
        @assert minimum(diff(λ)) * Δt > 1e-3

        M = zeros(ComplexF64, n, n)
        dk_divided_difference!(M, λ, Δt)

        φ = cis.(-Δt .* λ)
        Mref = zeros(ComplexF64, n, n)
        for a = 1:n, b = 1:n
            Mref[a, b] = a == b ? -im * Δt * φ[a] : (φ[a] - φ[b]) / (λ[a] - λ[b])
        end

        @test maximum(abs.(M .- Mref)) ≤ 1e-14
        # confluent diagonal is exactly -iΔt·exp(-iΔt λ_a)
        for a = 1:n
            @test M[a, a] == -im * Δt * cis(-Δt * λ[a])
        end
    end
end

@testitem "DK divided-difference is continuous through the degeneracy band (AC2)" begin
    using LinearAlgebra
    using ..QuantumIntegrators: dk_divided_difference!, DK_DEGENERACY_RTOL

    Δt = 0.85
    λa = 1.234

    # Well-conditioned closed form of the exact divided difference (no
    # cancellation): with m=(λa+λb)/2 and δ=(λa-λb)/2,
    #   (φa-φb)/(λa-λb) = -iΔt·exp(-iΔt m)·sinc0(Δt·δ),  sinc0(x)=sin(x)/x.
    sinc0(x) = x == 0 ? one(x) : sin(x) / x
    exact_dd(λ1, λ2) = begin
        m = (λ1 + λ2) / 2
        δ = (λ1 - λ2) / 2
        -im * Δt * cis(-Δt * m) * sinc0(Δt * δ)
    end

    # Sweep λb → λa across many decades, INCLUDING exact degeneracy (ε=0) and
    # straddling the band boundary ε* = DK_DEGENERACY_RTOL/Δt.
    εs = vcat(10.0 .^ range(-1, -13; length = 400), 0.0)
    prev = nothing
    for ε in εs
        M = zeros(ComplexF64, 2, 2)
        dk_divided_difference!(M, [λa, λa + ε], Δt)
        v = M[1, 2]
        @test isfinite(v)                             # no NaN/Inf anywhere
        @test isapprox(v, exact_dd(λa, λa + ε); atol = 1e-10)  # tracks the smooth truth
    end

    # At exact degeneracy the confluent (derivative) limit is used.
    Mdeg = zeros(ComplexF64, 2, 2)
    dk_divided_difference!(Mdeg, [λa, λa], Δt)
    @test Mdeg[1, 2] == -im * Δt * cis(-Δt * λa)

    # No jump AT the switchover: straddle the band boundary ε* by ±0.01% and
    # confirm the confluent side (just inside) and the divided-difference side
    # (just outside) agree to well below the 1e-9 gate — the only smooth drift
    # over this sub-nm interval is O(|φ''|·Δε) ≈ 1e-10.
    ε_star = DK_DEGENERACY_RTOL / Δt
    Mlo = zeros(ComplexF64, 2, 2)   # ε just inside band → confluent branch
    Mhi = zeros(ComplexF64, 2, 2)   # ε just outside band → divided-difference branch
    dk_divided_difference!(Mlo, [λa, λa + ε_star * 0.9999], Δt)
    dk_divided_difference!(Mhi, [λa, λa + ε_star * 1.0001], Δt)
    @test isapprox(Mlo[1, 2], Mhi[1, 2]; atol = 1e-9)
end

@testitem "DK first-order action matches central finite differences (AC3)" begin
    using LinearAlgebra
    using Random
    using ..QuantumIntegrators: dk_divided_difference!, dk_apply!

    Random.seed!(7)
    for _ = 1:12
        n = rand(3:6)
        A0 = randn(ComplexF64, n, n)
        H0 = Matrix(Hermitian(A0 + A0'))
        A1 = randn(ComplexF64, n, n)
        H1 = Matrix(Hermitian(A1 + A1'))   # affine drive direction ∂H/∂a
        Δt = 0.3 + rand()
        a0 = randn()

        H = H0 + a0 * H1
        F = eigen(Hermitian(H))
        V = Matrix(F.vectors)
        λ = F.values

        M = zeros(ComplexF64, n, n)
        dk_divided_difference!(M, λ, Δt)
        ∂Φ = zeros(ComplexF64, n, n)
        dk_apply!(∂Φ, V, M, H1, similar(∂Φ), similar(∂Φ))

        # central FD of exp(-iΔt H(a)) w.r.t. the affine control a
        h = 1e-6
        Φp = exp(-im * Δt * (H0 + (a0 + h) * H1))
        Φm = exp(-im * Δt * (H0 + (a0 - h) * H1))
        fd = (Φp .- Φm) ./ (2h)

        @test norm(∂Φ - fd) / norm(fd) < 1e-9
    end
end

@testitem "DK first-order action matches ForwardDiff-∘-expv witness (AC4)" begin
    using LinearAlgebra
    using Random
    using ForwardDiff
    using ExponentialAction: expv
    using Piccolo: Isomorphisms
    using ..QuantumIntegrators: dk_divided_difference!, dk_apply!

    Random.seed!(11)
    for _ = 1:12
        n = rand(2:5)
        A0 = randn(ComplexF64, n, n)
        H0 = Matrix(Hermitian(A0 + A0'))
        A1 = randn(ComplexF64, n, n)
        H1 = Matrix(Hermitian(A1 + A1'))
        Δt = 0.4 + rand()
        a0 = randn()

        H = H0 + a0 * H1
        F = eigen(Hermitian(H))
        V = Matrix(F.vectors)
        λ = F.values

        M = zeros(ComplexF64, n, n)
        dk_divided_difference!(M, λ, Δt)
        ∂Φ = zeros(ComplexF64, n, n)
        dk_apply!(∂Φ, V, M, H1, similar(∂Φ), similar(∂Φ))

        # Witness: ForwardDiff through `expv` on the real isomorphic generator,
        # mirroring the production ForwardDiff-∘-expv derivative. iso is a ring
        # homomorphism with iso(-iH)=G, so
        #   expv(Δt, iso(-iH0) + a·iso(-iH1), ψ̃) = iso(exp(-iΔt H(a)))·ψ,
        # whose a-derivative is [Re(∂Φ ψ); Im(∂Φ ψ)].
        G0 = Isomorphisms.iso(-im * H0)
        G1 = Isomorphisms.iso(-im * H1)
        ψ = randn(ComplexF64, n)
        ψ ./= norm(ψ)
        ψ̃ = vcat(real(ψ), imag(ψ))

        fdd = ForwardDiff.derivative(a -> expv(Δt, G0 + a * G1, ψ̃), a0)
        dkψ = ∂Φ * ψ
        dkiso = vcat(real(dkψ), imag(dkψ))

        @test norm(dkiso - fdd) / norm(fdd) < 1e-9
    end
end

@testitem "DK first-order derivative is zero-alloc after warmup (AC5)" begin
    using LinearAlgebra
    using Random
    using BenchmarkTools
    using ..QuantumIntegrators:
        dk_divided_difference!,
        dk_apply!,
        dk_first_order_derivative!,
        DaleckiiKreinWorkspace

    Random.seed!(3)
    n = 6
    A = randn(ComplexF64, n, n)
    H = Hermitian(A + A')
    Δt = 0.7
    F = eigen(H)
    V = Matrix(F.vectors)   # Matrix{ComplexF64}
    λ = F.values            # Vector{Float64}
    B = randn(ComplexF64, n, n)
    Hⱼ = Matrix(Hermitian(B + B'))

    ws = DaleckiiKreinWorkspace(n)   # M + 2 conjugation temporaries, sized to ketdim
    ∂Φ = zeros(ComplexF64, n, n)

    # warmup
    dk_first_order_derivative!(∂Φ, ws, V, λ, Hⱼ, Δt)

    a_build = @ballocated dk_divided_difference!($(ws.M), $λ, $Δt)
    a_apply = @ballocated dk_apply!($∂Φ, $V, $(ws.M), $Hⱼ, $(ws.tmp1), $(ws.tmp2))
    a_full = @ballocated dk_first_order_derivative!($∂Φ, $ws, $V, $λ, $Hⱼ, $Δt)

    @test a_build == 0
    @test a_apply == 0
    @test a_full == 0
end

# ============================================================================ #
# Second-order kernel tests (AC1–AC4 of Piccolissimo.jl#201). Reuse slice ①'s
# random-Hermitian fixtures and degeneracy-sweep harness.
# ============================================================================ #

@testitem "DK second-order action matches central-FD of the first-order action (AC1)" begin
    using LinearAlgebra
    using Random
    using ..QuantumIntegrators:
        dk_first_order_derivative!,
        dk_second_order_derivative!,
        DaleckiiKreinWorkspace,
        DaleckiiKreinSecondOrderWorkspace

    Random.seed!(2001)
    for _ = 1:12
        n = rand(3:6)
        A0 = randn(ComplexF64, n, n)
        H0 = Matrix(Hermitian(A0 + A0'))
        A1 = randn(ComplexF64, n, n)
        H1 = Matrix(Hermitian(A1 + A1'))     # affine drive direction ∂H/∂a
        Δt = 0.3 + rand()
        a0 = randn()
        μ = randn(ComplexF64, n, n)          # multiplier (any matrix; pairing is Re tr(μ†·))

        ws1 = DaleckiiKreinWorkspace(n)
        ws2 = DaleckiiKreinSecondOrderWorkspace(n)

        # analytic second derivative ∂²Φ/∂a² = D²Φ[H1,H1] at a0
        H = H0 + a0 * H1
        F = eigen(Hermitian(H))
        V = Matrix(F.vectors)
        λ = F.values
        ∂²Φ = zeros(ComplexF64, n, n)
        dk_second_order_derivative!(∂²Φ, ws2, V, λ, H1, H1, Δt)

        # central FD of the (slice-1) first-order analytic action w.r.t. a
        first_order(a) = begin
            Fa = eigen(Hermitian(H0 + a * H1))
            d = zeros(ComplexF64, n, n)
            dk_first_order_derivative!(d, ws1, Matrix(Fa.vectors), Fa.values, H1, Δt)
            d
        end
        h = 1e-5
        fd2 = (first_order(a0 + h) - first_order(a0 - h)) / (2h)

        # full-matrix action
        @test norm(∂²Φ - fd2) / norm(fd2) < 1e-8
        # μ-contracted action (the object #201 names): scalar ≡ Re tr(μ† ∂²Φ).
        # Measure the pairing error against its natural Cauchy–Schwarz scale
        # norm(μ)·norm(fd2) — the raw scalar can partially cancel to near zero,
        # which would spuriously inflate a value-relative error even though the
        # contracted matrices agree to the full-matrix tolerance above.
        @test abs(real(dot(μ, ∂²Φ)) - real(dot(μ, fd2))) / (norm(μ) * norm(fd2)) < 1e-8
    end
end

@testitem "DK second-order μ-block matches ForwardDiff.hessian witness (AC2)" begin
    using LinearAlgebra
    using Random
    using ForwardDiff
    using ExponentialAction: expv
    using Piccolo: Isomorphisms
    using ..QuantumIntegrators: dk_second_order_block!, DaleckiiKreinSecondOrderWorkspace

    Random.seed!(2002)
    for _ = 1:8
        n = rand(2:5)
        A0 = randn(ComplexF64, n, n)
        H0 = Matrix(Hermitian(A0 + A0'))
        A1 = randn(ComplexF64, n, n)
        H1 = Matrix(Hermitian(A1 + A1'))
        A2 = randn(ComplexF64, n, n)
        H2 = Matrix(Hermitian(A2 + A2'))
        Δt = 0.4 + rand()
        apt = randn(2)

        # μ = χψ† realises the real pairing Re tr(μ† X) = Re(χ† X ψ), matching the
        # ForwardDiff witness scalar s(a) = Re(χ† Φ(a) ψ) built via the iso `expv`
        # path (iso ring-homomorphism ⇒ expv(Δt, iso(-iH(a)), ψ̃) = iso(Φ(a)ψ)).
        χ = randn(ComplexF64, n)
        χ ./= norm(χ)
        ψ = randn(ComplexF64, n)
        ψ ./= norm(ψ)
        μ = χ * ψ'
        χ̃ = vcat(real(χ), imag(χ))
        ψ̃ = vcat(real(ψ), imag(ψ))
        G0 = Isomorphisms.iso(-im * H0)
        G1 = Isomorphisms.iso(-im * H1)
        G2 = Isomorphisms.iso(-im * H2)
        s(a) = dot(χ̃, expv(Δt, G0 + a[1] * G1 + a[2] * G2, ψ̃))
        Hess = ForwardDiff.hessian(s, apt)

        # analytic μ-contracted parameter-parameter block at apt
        H = H0 + apt[1] * H1 + apt[2] * H2
        F = eigen(Hermitian(H))
        V = Matrix(F.vectors)
        λ = F.values
        ws2 = DaleckiiKreinSecondOrderWorkspace(n)
        B = zeros(2, 2)
        dk_second_order_block!(B, ws2, V, λ, [H1, H2], μ, Δt)

        @test norm(B - Hess) / norm(Hess) < 1e-8
    end
end

@testitem "DK second-order handles the near-degenerate confluent limit (AC3)" begin
    using LinearAlgebra
    using Random
    using ForwardDiff
    using ExponentialAction: expv
    using Piccolo: Isomorphisms
    using ..QuantumIntegrators:
        dk_first_order_derivative!,
        dk_second_order_derivative!,
        dk_second_order_block!,
        DaleckiiKreinWorkspace,
        DaleckiiKreinSecondOrderWorkspace

    Random.seed!(2003)
    Δt = 0.6
    # Test both an EXACT 2-fold degeneracy and a tiny (in-band) gap.
    for gap in (0.0, 1e-9)
        n = 5
        Q = Matrix(qr(randn(ComplexF64, n, n)).Q)
        λ0 = sort(randn(n) .* 2.0)
        λ0[2] = λ0[1] + gap                 # (near-)degenerate pair inside the band
        H0m = Q * Diagonal(λ0) * Q'
        H0 = Matrix(Hermitian((H0m + H0m') / 2))
        B1 = randn(ComplexF64, n, n)
        H1 = Matrix(Hermitian(B1 + B1'))
        B2 = randn(ComplexF64, n, n)
        H2 = Matrix(Hermitian(B2 + B2'))

        # kernel second derivative at a0 = 0 (spectrum is (near-)degenerate here,
        # so the confluent triples are exercised)
        F = eigen(Hermitian(H0))
        V = Matrix(F.vectors)
        λ = F.values
        ws2 = DaleckiiKreinSecondOrderWorkspace(n)
        ∂²Φ = zeros(ComplexF64, n, n)
        dk_second_order_derivative!(∂²Φ, ws2, V, λ, H1, H1, Δt)
        @test all(isfinite, ∂²Φ)            # finite — no catastrophic cancellation

        # oracle: central FD of the slice-1 first-order action (stable near
        # degeneracy; the ±h perturbation splits the pair into the recursion band)
        ws1 = DaleckiiKreinWorkspace(n)
        first_order(a) = begin
            Fa = eigen(Hermitian(H0 + a * H1))
            d = zeros(ComplexF64, n, n)
            dk_first_order_derivative!(d, ws1, Matrix(Fa.vectors), Fa.values, H1, Δt)
            d
        end
        h = 1e-5
        fd2 = (first_order(h) - first_order(-h)) / (2h)
        @test norm(∂²Φ - fd2) / norm(fd2) < 1e-8

        # μ-contracted mixed block vs a stable ForwardDiff.hessian witness at the
        # (near-)degenerate point (expv uses no eigenbasis, so it stays stable)
        χ = randn(ComplexF64, n)
        χ ./= norm(χ)
        ψ = randn(ComplexF64, n)
        ψ ./= norm(ψ)
        μ = χ * ψ'
        χ̃ = vcat(real(χ), imag(χ))
        ψ̃ = vcat(real(ψ), imag(ψ))
        G0 = Isomorphisms.iso(-im * H0)
        G1 = Isomorphisms.iso(-im * H1)
        G2 = Isomorphisms.iso(-im * H2)
        s(a) = dot(χ̃, expv(Δt, G0 + a[1] * G1 + a[2] * G2, ψ̃))
        Hess = ForwardDiff.hessian(s, [0.0, 0.0])
        B = zeros(2, 2)
        dk_second_order_block!(B, ws2, V, λ, [H1, H2], μ, Δt)
        @test all(isfinite, B)
        @test norm(B - Hess) / norm(Hess) < 1e-8
    end
end

@testitem "DK second divided difference is finite and continuous through confluence (AC3)" begin
    using LinearAlgebra
    using ..QuantumIntegrators: dk_second_divided_difference, DK_DEGENERACY_RTOL

    Δt = 0.85
    λ0 = 1.234

    # Full triple degeneracy → confluent (derivative) limit φ⁽²⁾(x,x,x)=φ″(x)/2.
    φ2diag = -Δt^2 * cis(-Δt * λ0) / 2
    @test isapprox(dk_second_divided_difference(λ0, λ0, λ0, Δt), φ2diag; rtol = 1e-13)

    # Symmetry in all three arguments (permutation invariance).
    v = dk_second_divided_difference(λ0 + 0.7, λ0 - 0.3, λ0 + 1.1, Δt)
    @test dk_second_divided_difference(λ0 - 0.3, λ0 + 1.1, λ0 + 0.7, Δt) ≈ v rtol = 1e-14
    @test dk_second_divided_difference(λ0 + 1.1, λ0 + 0.7, λ0 - 0.3, Δt) ≈ v rtol = 1e-14

    # Sweep all three nodes colliding, (λ0+ε, λ0, λ0-ε), across many decades
    # (incl. exact confluence ε=0) — nothing may go NaN/Inf.
    for ε in vcat(10.0 .^ range(-1, -13; length = 300), 0.0)
        val = dk_second_divided_difference(λ0 + ε, λ0, λ0 - ε, Δt)
        @test isfinite(val)
    end

    # No jump AT the all-three switchover: straddle the band boundary (gmax = 2ε)
    # by ±0.01% — confluent side (inside) vs recursion side (outside) agree well
    # below the 1e-8 gate.
    ε_star = DK_DEGENERACY_RTOL / (2Δt)
    vlo = dk_second_divided_difference(λ0 + 0.9999ε_star, λ0, λ0 - 0.9999ε_star, Δt)
    vhi = dk_second_divided_difference(λ0 + 1.0001ε_star, λ0, λ0 - 1.0001ε_star, Δt)
    @test isapprox(vlo, vhi; atol = 1e-9)

    # A single degenerate pair with the third node far uses the recursion via the
    # largest gap; the inner sinc-form φ¹ keeps it stable as the pair collides.
    # As ε→0 the value converges CONTINUOUSLY to the exactly-confluent-pair value,
    # differing by O(ε) (slope ≈ |∂φ²/∂node| ≈ 0.09) — a Lipschitz-in-ε bound. A
    # cancelling divided difference would instead deviate like ε_machine/ε, which
    # GROWS as ε→0; the bound below rules that out across 11 decades of collision.
    exact_pair = dk_second_divided_difference(λ0, λ0, λ0 + 3.0, Δt)   # a==b exactly
    for ε in 10.0 .^ range(-2, -13; length = 50)
        val = dk_second_divided_difference(λ0 + ε, λ0, λ0 + 3.0, Δt)
        @test isfinite(val)
        @test abs(val - exact_pair) ≤ ε + 1e-13
    end
end

@testitem "DK second-order block is symmetric to machine precision (AC4)" begin
    using LinearAlgebra
    using Random
    using ..QuantumIntegrators:
        dk_second_order_apply!, dk_second_order_block!, DaleckiiKreinSecondOrderWorkspace

    Random.seed!(2004)
    n = 5
    A0 = randn(ComplexF64, n, n)
    H0 = Matrix(Hermitian(A0 + A0'))
    Hs = map(1:3) do _
        B = randn(ComplexF64, n, n)
        Matrix(Hermitian(B + B'))
    end
    Δt = 0.7
    F = eigen(Hermitian(H0))
    V = Matrix(F.vectors)
    λ = F.values
    μc = randn(ComplexF64, n, n)

    ws = DaleckiiKreinSecondOrderWorkspace(n)

    # (a) structural swap symmetry: D²Φ[Hᵢ,Hⱼ] == D²Φ[Hⱼ,Hᵢ] to machine precision
    D12 = zeros(ComplexF64, n, n)
    D21 = zeros(ComplexF64, n, n)
    dk_second_order_apply!(D12, V, λ, Hs[1], Hs[2], Δt, ws.E1t, ws.E2t, ws.core, ws.tmp)
    dk_second_order_apply!(D21, V, λ, Hs[2], Hs[1], Δt, ws.E1t, ws.E2t, ws.core, ws.tmp)
    @test norm(D12 - D21) ≤ 1e-13 * norm(D12)

    # (b) assembled block is symmetric to machine precision
    B = zeros(3, 3)
    dk_second_order_block!(B, ws, V, λ, Hs, μc, Δt)
    @test norm(B - B') ≤ 1e-12 * max(norm(B), eps())
    @test all(isfinite, B)
end

@testitem "DK second-order derivative is zero-alloc after warmup" begin
    using LinearAlgebra
    using Random
    using BenchmarkTools
    using ..QuantumIntegrators:
        dk_second_order_apply!,
        dk_second_order_derivative!,
        DaleckiiKreinSecondOrderWorkspace

    Random.seed!(2005)
    n = 6
    A = randn(ComplexF64, n, n)
    H = Hermitian(A + A')
    Δt = 0.7
    F = eigen(H)
    V = Matrix(F.vectors)
    λ = F.values
    B1 = randn(ComplexF64, n, n)
    H1 = Matrix(Hermitian(B1 + B1'))
    B2 = randn(ComplexF64, n, n)
    H2 = Matrix(Hermitian(B2 + B2'))

    ws = DaleckiiKreinSecondOrderWorkspace(n)
    ∂²Φ = zeros(ComplexF64, n, n)

    # warmup
    dk_second_order_derivative!(∂²Φ, ws, V, λ, H1, H2, Δt)

    a_apply = @ballocated dk_second_order_apply!(
        $∂²Φ,
        $V,
        $λ,
        $H1,
        $H2,
        $Δt,
        $(ws.E1t),
        $(ws.E2t),
        $(ws.core),
        $(ws.tmp),
    )
    a_full = @ballocated dk_second_order_derivative!($∂²Φ, $ws, $V, $λ, $H1, $H2, $Δt)

    @test a_apply == 0
    @test a_full == 0
end
