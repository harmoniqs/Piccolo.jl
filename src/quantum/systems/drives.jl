# ----------------------------------------------------------------------------- #
# Drive Types: Operator + Coefficient for Hamiltonian Terms
#
# A Drive pairs a Hamiltonian matrix with a scalar coefficient function of the
# control vector u. This enables both linear drives (u[i] * H_i) and nonlinear
# drives (f(u) * H) in the same framework.
#
# The Hamiltonian is: H(u, t) = H_drift + Σ_d drive_coeff(d, u) * d.H
#
# For analytical sensitivity equations (e.g., spline integrators), the chain
# rule requires: ∂H/∂u_j = Σ_d drive_coeff_jac(d, u, j) * d.H
# ----------------------------------------------------------------------------- #

"""
    AbstractDrive

Abstract supertype for Hamiltonian drive terms.

A drive pairs a Hermitian matrix `H` with a scalar coefficient that depends on the
control vector `u`. The full Hamiltonian is:

    H(u, t) = H_drift + Σ_d drive_coeff(d, u) * d.H

Subtypes must implement:
- `drive_coeff(d, u)` — compute the scalar coefficient
- `drive_coeff_jac(d, u, j)` — compute ∂coeff/∂u_j

See [`LinearDrive`](@ref) and [`NonlinearDrive`](@ref).
"""
abstract type AbstractDrive end

"""
    LinearDrive <: AbstractDrive

Standard linear drive: coefficient is `u[index]`.

This is the default representation when constructing `QuantumSystem(H_drift, H_drives, bounds)`.

# Fields
- `H::SparseMatrixCSC{ComplexF64,Int}`: The Hermitian drive operator
- `index::Int`: Index into the control vector `u`

# Example
```julia
LinearDrive(sparse(PAULIS.X), 1)  # u[1] * σx
```
"""
struct LinearDrive <: AbstractDrive
    H::SparseMatrixCSC{ComplexF64,Int}
    index::Int
end

"""
    NonlinearDrive{F,DF} <: AbstractDrive

Drive term with a nonlinear scalar coefficient: `f(u) * H`.

The user provides both the coefficient function and its Jacobian (derivative
with respect to each control component). This enables analytical sensitivity
equations without finite differences.

# Fields
- `H::SparseMatrixCSC{ComplexF64,Int}`: The Hermitian drive operator
- `coeff::F`: `(u::AbstractVector) -> scalar` — coefficient function
- `coeff_jac::DF`: `(u::AbstractVector, j::Int) -> scalar` — ∂coeff/∂u_j

# Example: Displaced frame |α|² term
```julia
NonlinearDrive(
    sparse(σz / 2),
    u -> u[3]^2 + u[4]^2,                                    # αI² + αQ²
    (u, j) -> j == 3 ? 2u[3] : j == 4 ? 2u[4] : 0.0        # ∂/∂u_j
)
```

# Example: Product of controls
```julia
NonlinearDrive(
    sparse(H_coupling),
    u -> u[1] * u[2],                                         # u₁ · u₂
    (u, j) -> j == 1 ? u[2] : j == 2 ? u[1] : 0.0           # ∂/∂u_j
)
```
"""
struct NonlinearDrive{F,DF} <: AbstractDrive
    H::SparseMatrixCSC{ComplexF64,Int}
    coeff::F
    coeff_jac::DF
end

# Convenience: accept any AbstractMatrix and sparsify
function NonlinearDrive(H::AbstractMatrix, coeff::F, coeff_jac::DF) where {F,DF}
    return NonlinearDrive(sparse(ComplexF64.(H)), coeff, coeff_jac)
end

# ----------------------------------------------------------------------------- #
# Interface
# ----------------------------------------------------------------------------- #

"""
    drive_coeff(d::AbstractDrive, u::AbstractVector) -> Number

Compute the scalar coefficient of this drive at controls `u`.
"""
@inline drive_coeff(d::LinearDrive, u::AbstractVector) = u[d.index]
@inline drive_coeff(d::NonlinearDrive, u::AbstractVector) = d.coeff(u)

"""
    drive_coeff_jac(d::AbstractDrive, u::AbstractVector, j::Int) -> Number

Compute ∂coeff/∂u_j at controls `u`.

For `LinearDrive`: returns 1.0 if `j == d.index`, 0.0 otherwise (Kronecker delta).
For `NonlinearDrive`: evaluates the user-provided Jacobian function.
"""
@inline drive_coeff_jac(d::LinearDrive, ::AbstractVector, j::Int) = j == d.index ? 1.0 : 0.0
@inline drive_coeff_jac(d::NonlinearDrive, u::AbstractVector, j::Int) = d.coeff_jac(u, j)

"""
    has_nonlinear_drives(drives::Vector{<:AbstractDrive}) -> Bool

Check if any drive terms are nonlinear.
"""
has_nonlinear_drives(drives::AbstractVector{<:AbstractDrive}) = any(d -> d isa NonlinearDrive, drives)

# ----------------------------------------------------------------------------- #
# Tests
# ----------------------------------------------------------------------------- #

@testitem "LinearDrive" begin
    using Piccolo
    using SparseArrays

    H = sparse([0.0+0im 1.0+0im; 1.0+0im 0.0+0im])
    d = LinearDrive(H, 2)

    u = [0.3, 0.7, 0.1]
    @test drive_coeff(d, u) == 0.7
    @test drive_coeff_jac(d, u, 1) == 0.0
    @test drive_coeff_jac(d, u, 2) == 1.0
    @test drive_coeff_jac(d, u, 3) == 0.0
    @test d.H === H
end

@testitem "NonlinearDrive" begin
    using Piccolo
    using SparseArrays

    H = sparse([1.0+0im 0.0+0im; 0.0+0im -1.0+0im])

    # Quadratic: u[1]^2 + u[2]^2
    d = NonlinearDrive(
        H,
        u -> u[1]^2 + u[2]^2,
        (u, j) -> j == 1 ? 2u[1] : j == 2 ? 2u[2] : 0.0
    )

    u = [3.0, 4.0, 0.0]
    @test drive_coeff(d, u) == 25.0
    @test drive_coeff_jac(d, u, 1) == 6.0
    @test drive_coeff_jac(d, u, 2) == 8.0
    @test drive_coeff_jac(d, u, 3) == 0.0
    @test d.H === H

    # Product: u[1] * u[2]
    d2 = NonlinearDrive(
        H,
        u -> u[1] * u[2],
        (u, j) -> j == 1 ? u[2] : j == 2 ? u[1] : 0.0
    )

    @test drive_coeff(d2, u) == 12.0
    @test drive_coeff_jac(d2, u, 1) == 4.0
    @test drive_coeff_jac(d2, u, 2) == 3.0
end

@testitem "NonlinearDrive accepts dense matrix" begin
    using Piccolo
    using SparseArrays

    H_dense = [1.0+0im 0.0+0im; 0.0+0im -1.0+0im]
    d = NonlinearDrive(H_dense, u -> u[1], (u, j) -> j == 1 ? 1.0 : 0.0)

    @test d.H isa SparseMatrixCSC
    @test drive_coeff(d, [0.5]) == 0.5
end
