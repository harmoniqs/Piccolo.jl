# ----------------------------------------------------------------------------- #
# AbstractDissipator hierarchy
# ----------------------------------------------------------------------------- #
# Dissipators generalize the constant-matrix `L` jump operators to
# `rate(u)·L` jump operators, letting dissipation rates participate in
# QILC / global-param calibration the same way NonlinearDrive coefficients do.
#
# Per-dissipator Lindblad contribution to 𝒟:
#     𝒟_j(u) = rate_coeff(d_j, u) · iso_D(dissipator_matrix(d_j))
#
# The Hamiltonian is unchanged: H(u,t) = H_drift + Σ_d drive_coeff(d,u)·d.H
# Assemblers compose 𝒢_drift + Σ drive_coeff·𝒢_drive + Σ rate_coeff·𝒟.

"""
    AbstractDissipator

Abstract type for dissipators: jump operators with a scalar rate coefficient.
See [`LinearDissipator`](@ref) and [`NonlinearDissipator`](@ref).
"""
abstract type AbstractDissipator end

"""
    LinearDissipator <: AbstractDissipator

Constant-rate dissipator. `rate_coeff(d, u) = d.rate` for all `u`.
Absorb a constant prefactor (e.g. `γ`) into `d.rate`, or keep `L` scaled and
leave `rate = 1.0`.

# Fields
- `L::SparseMatrixCSC{ComplexF64,Int}`: The jump operator
- `rate::Float64`: Constant scalar rate
"""
struct LinearDissipator <: AbstractDissipator
    L::SparseMatrixCSC{ComplexF64,Int}
    rate::Float64
end

function LinearDissipator(L::AbstractMatrix, rate::Real = 1.0)
    return LinearDissipator(sparse(ComplexF64.(L)), Float64(rate))
end

"""
    NonlinearDissipator{F,DF,HF} <: AbstractDissipator

Dissipator with a `u`-dependent rate. `rate_coeff(d, u) = d.rate(u)`.

Analogous to [`NonlinearDrive`](@ref): supports auto-generated Jacobians /
Hessians via ForwardDiff or user-supplied closures, with an `active_controls`
list for sparsity.

# Fields
- `L::SparseMatrixCSC{ComplexF64,Int}`: The jump operator
- `rate::F`: `(u::AbstractVector) -> scalar` — rate coefficient function
- `rate_jac::DF`: `(u::AbstractVector, j::Int) -> scalar` — ∂rate/∂u_j
- `rate_hess::HF`: `(u::AbstractVector, i::Int, j::Int) -> scalar` — ∂²rate/∂u_i∂u_j
- `active_controls::Vector{Int}`: Control indices where ∂rate/∂u_j can be nonzero.
  Empty means "all controls" (no structural sparsity info).
"""
struct NonlinearDissipator{F,DF,HF} <: AbstractDissipator
    L::SparseMatrixCSC{ComplexF64,Int}
    rate::F
    rate_jac::DF
    rate_hess::HF
    active_controls::Vector{Int}
end

"""
    NonlinearDissipator(L, rate; active_controls=Int[])

Construct a `NonlinearDissipator` with auto-generated Jacobian and Hessian via
ForwardDiff. `rate` is a scalar-valued function `u -> Number`.
"""
function NonlinearDissipator(
    L::AbstractMatrix,
    rate::F;
    active_controls::Vector{Int} = Int[],
) where {F}
    rate_jac = (u, j) -> ForwardDiff.gradient(rate, u)[j]
    rate_hess = (u, i, j) -> ForwardDiff.hessian(rate, u)[i, j]
    return NonlinearDissipator(
        sparse(ComplexF64.(L)),
        rate,
        rate_jac,
        rate_hess,
        active_controls,
    )
end

"""
    NonlinearDissipator(L, rate, rate_jac; active_controls=Int[], rate_hess=nothing)

Construct a `NonlinearDissipator` with an explicit Jacobian function; Hessian
is auto-generated via ForwardDiff unless supplied.
"""
function NonlinearDissipator(
    L::AbstractMatrix,
    rate::F,
    rate_jac::DF;
    active_controls::Vector{Int} = Int[],
    rate_hess = nothing,
) where {F,DF}
    hess =
        isnothing(rate_hess) ? (u, i, j) -> ForwardDiff.hessian(rate, u)[i, j] :
        rate_hess
    return NonlinearDissipator(
        sparse(ComplexF64.(L)),
        rate,
        rate_jac,
        hess,
        active_controls,
    )
end

# ----------------------------------------------------------------------------- #
# Interface
# ----------------------------------------------------------------------------- #

"""
    rate_coeff(d::AbstractDissipator, u::AbstractVector) -> Number

Scalar rate of this dissipator at controls `u`.
"""
@inline rate_coeff(d::LinearDissipator, ::AbstractVector) = d.rate
@inline rate_coeff(d::NonlinearDissipator, u::AbstractVector) = d.rate(u)

"""
    rate_coeff_jac(d::AbstractDissipator, u::AbstractVector, j::Int) -> Number

Partial derivative `∂rate/∂u_j` at controls `u`.
"""
@inline rate_coeff_jac(::LinearDissipator, ::AbstractVector, ::Int) = 0.0
@inline rate_coeff_jac(d::NonlinearDissipator, u::AbstractVector, j::Int) =
    d.rate_jac(u, j)

"""
    rate_coeff_hess(d::AbstractDissipator, u::AbstractVector, i::Int, j::Int) -> Number

Second partial `∂²rate/∂u_i∂u_j` at controls `u`.
"""
@inline rate_coeff_hess(::LinearDissipator, ::AbstractVector, ::Int, ::Int) = 0.0
@inline rate_coeff_hess(d::NonlinearDissipator, u::AbstractVector, i::Int, j::Int) =
    d.rate_hess(u, i, j)

"""
    active_controls(d::AbstractDissipator) -> Vector{Int}

Control indices where `rate_coeff_jac(d, u, j)` can be nonzero.
"""
active_controls(::LinearDissipator) = Int[]
active_controls(d::NonlinearDissipator) = d.active_controls

"""
    has_nonlinear_dissipators(dissipators) -> Bool
"""
has_nonlinear_dissipators(ds::AbstractVector{<:AbstractDissipator}) =
    any(d -> d isa NonlinearDissipator, ds)

# ----------------------------------------------------------------------------- #
# Matrix access
# ----------------------------------------------------------------------------- #

"""
    dissipator_matrix(d::AbstractDissipator) -> AbstractMatrix

Materialize the jump-operator `L` as a matrix.
"""
dissipator_matrix(d::AbstractDissipator) = d.L

# ----------------------------------------------------------------------------- #
# Tests
# ----------------------------------------------------------------------------- #

@testitem "AbstractDissipator: LinearDissipator dispatch" begin
    using Piccolo
    using SparseArrays
    L = sqrt(0.3) * PAULIS.Z
    d = LinearDissipator(L)
    @test rate_coeff(d, [0.1, 0.2]) == 1.0
    @test rate_coeff_jac(d, [0.1, 0.2], 1) == 0.0
    @test rate_coeff_hess(d, [0.1, 0.2], 1, 2) == 0.0
    @test active_controls(d) == Int[]
    @test dissipator_matrix(d) == sparse(ComplexF64.(L))

    d2 = LinearDissipator(PAULIS.X, 0.7)
    @test rate_coeff(d2, [0.0]) == 0.7
end

@testitem "AbstractDissipator: NonlinearDissipator with explicit jac" begin
    using Piccolo
    using SparseArrays
    L = PAULIS.Z / sqrt(2)
    rate = u -> u[2]
    rate_jac = (u, j) -> j == 2 ? 1.0 : 0.0
    d = NonlinearDissipator(L, rate, rate_jac; active_controls=[2])
    @test rate_coeff(d, [0.3, 0.7]) ≈ 0.7
    @test rate_coeff_jac(d, [0.3, 0.7], 2) ≈ 1.0
    @test rate_coeff_jac(d, [0.3, 0.7], 1) ≈ 0.0
    @test active_controls(d) == [2]
end

@testitem "AbstractDissipator: NonlinearDissipator ForwardDiff auto-jac matches hand" begin
    using Piccolo
    using SparseArrays
    L = PAULIS.Z / sqrt(2)
    # Quadratic rate: γ(u) = u[1]^2 + 2·u[2]
    rate = u -> u[1]^2 + 2 * u[2]
    d = NonlinearDissipator(L, rate; active_controls=[1, 2])
    u = [0.3, 0.5]
    @test rate_coeff(d, u) ≈ 0.3^2 + 2 * 0.5
    @test rate_coeff_jac(d, u, 1) ≈ 2 * 0.3
    @test rate_coeff_jac(d, u, 2) ≈ 2.0
    @test rate_coeff_hess(d, u, 1, 1) ≈ 2.0
    @test rate_coeff_hess(d, u, 1, 2) ≈ 0.0
    @test rate_coeff_hess(d, u, 2, 2) ≈ 0.0
end

@testitem "AbstractDissipator: has_nonlinear_dissipators" begin
    using Piccolo
    lin = LinearDissipator(PAULIS.Z)
    nonlin = NonlinearDissipator(PAULIS.Z, u -> u[1]; active_controls=[1])
    @test !has_nonlinear_dissipators([lin, lin])
    @test has_nonlinear_dissipators(AbstractDissipator[lin, nonlin])
end
