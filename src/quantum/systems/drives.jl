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
# For second-order sensitivity equations (exact Hessian):
# ∂²H/∂u_i∂u_j = Σ_d drive_coeff_hess(d, u, i, j) * d.H
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
- `drive_coeff_hess(d, u, i, j)` — compute ∂²coeff/∂u_i∂u_j
- `active_controls(d)` — return indices where ∂coeff/∂u_j can be nonzero

See [`LinearDrive`](@ref) and [`NonlinearDrive`](@ref).
"""
abstract type AbstractDrive end

"""
    LinearDrive{H} <: AbstractDrive

Standard linear drive: coefficient is `u[index]`.

This is the default representation when constructing `QuantumSystem(H_drift, H_drives, bounds)`.

# Fields
- `H::H`: The Hermitian drive operator (matrix or `AbstractDynamicsOperator`)
- `index::Int`: Index into the control vector `u`

# Example
```julia
LinearDrive(sparse(PAULIS.X), 1)  # u[1] * σx
```
"""
struct LinearDrive{H} <: AbstractDrive
    H::H
    index::Int
end

"""
    NonlinearDrive{H,F,DF,HF} <: AbstractDrive

Drive term with a nonlinear scalar coefficient: `f(u) * H`.

The user provides either:
1. Both the coefficient function and its Jacobian (3-arg form), or
2. Just the coefficient function (2-arg form) — the Jacobian and Hessian are
   computed automatically via ForwardDiff.

The Hessian (`coeff_hess`) is used by the second-order sensitivity ODE for exact
Hessian computation in SplineIntegrator (`exact_hessian=true`).

# Fields
- `H::H`: The Hermitian drive operator (matrix or `AbstractDynamicsOperator`)
- `coeff::F`: `(u::AbstractVector) -> scalar` — coefficient function
- `coeff_jac::DF`: `(u::AbstractVector, j::Int) -> scalar` — ∂coeff/∂u_j
- `coeff_hess::HF`: `(u::AbstractVector, i::Int, j::Int) -> scalar` — ∂²coeff/∂u_i∂u_j
- `active_controls::Vector{Int}`: Control indices where ∂coeff/∂u_j can be nonzero.
  Empty means "all controls" (no structural sparsity info).

# Example: Auto-Jacobian (recommended)
```julia
NonlinearDrive(σz / 2, u -> u[3]^2 + u[4]^2)
```

# Example: Manual Jacobian with active controls
```julia
NonlinearDrive(
    σz / 2,
    u -> u[3]^2 + u[4]^2,
    (u, j) -> j == 3 ? 2u[3] : j == 4 ? 2u[4] : 0.0;
    active_controls = [3, 4]
)
```
"""
struct NonlinearDrive{H,F,DF,HF} <: AbstractDrive
    H::H
    coeff::F
    coeff_jac::DF
    coeff_hess::HF
    active_controls::Vector{Int}
end

"""
    NonlinearDrive(H::AbstractMatrix, coeff, coeff_jac; active_controls=Int[], coeff_hess=nothing)

Construct a `NonlinearDrive` with an explicit Jacobian function.
`active_controls` lists which control indices have nonzero ∂coeff/∂u_j (empty = all).
If `coeff_hess` is not provided, it is auto-generated via ForwardDiff from `coeff`.
"""
function NonlinearDrive(
    H::AbstractMatrix,
    coeff::F,
    coeff_jac::DF;
    active_controls::Vector{Int} = Int[],
    coeff_hess = nothing,
) where {F,DF}
    hess = isnothing(coeff_hess) ? _forwarddiff_hessian(coeff) : coeff_hess
    return NonlinearDrive(sparse(ComplexF64.(H)), coeff, coeff_jac, hess, active_controls)
end

"""
    NonlinearDrive(H::AbstractMatrix, coeff; active_controls=Int[])

Construct a `NonlinearDrive` with an auto-generated Jacobian via ForwardDiff.

This is the recommended constructor — the Jacobian is computed automatically from
the coefficient function, eliminating the risk of hand-written Jacobian errors.

# Example
```julia
# Displaced-frame Stark shift: coefficient = u₃² + u₄²
drive = NonlinearDrive(σz / 2, u -> u[3]^2 + u[4]^2; active_controls = [3, 4])
```
"""
function NonlinearDrive(
    H::AbstractMatrix,
    coeff::F;
    active_controls::Vector{Int} = Int[],
) where {F}
    coeff_jac = _forwarddiff_jacobian(coeff)
    coeff_hess = _forwarddiff_hessian(coeff)
    return NonlinearDrive(
        sparse(ComplexF64.(H)),
        coeff,
        coeff_jac,
        coeff_hess,
        active_controls,
    )
end

# ── Generic constructors (for structured operators or any non-matrix Hamiltonian) ──

"""
    NonlinearDrive(H, coeff; active_controls=Int[])

Construct a `NonlinearDrive` with any Hamiltonian type (e.g., `AbstractDynamicsOperator`).
The Hamiltonian is stored as-is (not materialized or sparsified).
"""
function NonlinearDrive(H, coeff::F; active_controls::Vector{Int} = Int[]) where {F}
    coeff_jac = _forwarddiff_jacobian(coeff)
    coeff_hess = _forwarddiff_hessian(coeff)
    return NonlinearDrive(H, coeff, coeff_jac, coeff_hess, active_controls)
end

"""
    NonlinearDrive(H, coeff, coeff_jac; active_controls=Int[], coeff_hess=nothing)

Construct a `NonlinearDrive` with any Hamiltonian type and explicit Jacobian.
"""
function NonlinearDrive(
    H,
    coeff::F,
    coeff_jac::DF;
    active_controls::Vector{Int} = Int[],
    coeff_hess = nothing,
) where {F,DF}
    hess = isnothing(coeff_hess) ? _forwarddiff_hessian(coeff) : coeff_hess
    return NonlinearDrive(H, coeff, coeff_jac, hess, active_controls)
end

"""
    _forwarddiff_jacobian(f) -> (u, j) -> ∂f/∂u_j

Create a Jacobian function from a scalar-valued function `f(u)` using ForwardDiff.
"""
function _forwarddiff_jacobian(f::F) where {F}
    return (u, j) -> ForwardDiff.gradient(f, u)[j]
end

"""
    _forwarddiff_hessian(f) -> (u, i, j) -> ∂²f/∂u_i∂u_j

Create a Hessian function from a scalar-valued function `f(u)` using ForwardDiff.
"""
function _forwarddiff_hessian(f::F) where {F}
    return (u, i, j) -> ForwardDiff.hessian(f, u)[i, j]
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
For `NonlinearDrive`: evaluates the Jacobian function (user-provided or auto-generated).
"""
@inline drive_coeff_jac(d::LinearDrive, ::AbstractVector, j::Int) = j == d.index ? 1.0 : 0.0
@inline drive_coeff_jac(d::NonlinearDrive, u::AbstractVector, j::Int) = d.coeff_jac(u, j)

"""
    drive_coeff_hess(d::AbstractDrive, u::AbstractVector, i::Int, j::Int) -> Number

Compute ∂²coeff/∂u_i∂u_j at controls `u`.

For `LinearDrive`: always returns 0.0 (linear coefficient has zero Hessian).
For `NonlinearDrive`: evaluates the Hessian function (user-provided or auto-generated).
"""
@inline drive_coeff_hess(d::LinearDrive, ::AbstractVector, ::Int, ::Int) = 0.0
@inline drive_coeff_hess(d::NonlinearDrive, u::AbstractVector, i::Int, j::Int) =
    d.coeff_hess(u, i, j)

"""
    active_controls(d::AbstractDrive) -> Vector{Int}

Return the control indices where `drive_coeff_jac(d, u, j)` can be nonzero.
An empty vector means "all controls" (no structural sparsity information available).

For `LinearDrive`, returns `[d.index]`.
For `NonlinearDrive`, returns `d.active_controls` (user-specified or empty).
"""
active_controls(d::LinearDrive) = [d.index]
active_controls(d::NonlinearDrive) = d.active_controls

"""
    has_nonlinear_drives(drives::Vector{<:AbstractDrive}) -> Bool

Check if any drive terms are nonlinear.
"""
has_nonlinear_drives(drives::AbstractVector{<:AbstractDrive}) =
    any(d -> d isa NonlinearDrive, drives)

# ----------------------------------------------------------------------------- #
# Matrix access
# ----------------------------------------------------------------------------- #

"""
    _ensure_matrix(H) -> AbstractMatrix

Convert `H` to a matrix if it isn't one already. For `AbstractMatrix` inputs,
returns `H` directly. For operator types (e.g., from Piccolissimo's operators
module), falls back to `Matrix(H)` — operator packages should define
`Base.Matrix(op)` to return a dense matrix representation.
"""
_ensure_matrix(H::AbstractMatrix) = H
_ensure_matrix(H) = Matrix(H)

"""
    drive_matrix(d::AbstractDrive) -> AbstractMatrix

Return the drive operator as a matrix. Use this instead of `d.H` when a matrix
is required (e.g., for `Isomorphisms.G`, Hermiticity checks, or `H(u,t)` closures).
"""
drive_matrix(d::AbstractDrive) = _ensure_matrix(d.H)

"""
    drive_dim(d::AbstractDrive) -> Int

Return the Hilbert space dimension of the drive operator.
"""
drive_dim(d::AbstractDrive) = size(drive_matrix(d), 1)

# ----------------------------------------------------------------------------- #
# Isomorphism dispatch for drives
# ----------------------------------------------------------------------------- #

"""
    Isomorphisms.G(d::AbstractDrive)

Delegate `G` to the underlying Hamiltonian matrix `d.H`, so that broadcasting
`G.(sys.H_drives)` works when `H_drives` contains `AbstractDrive` objects.
"""
Isomorphisms.G(d::AbstractDrive) = Isomorphisms.G(drive_matrix(d))

# ----------------------------------------------------------------------------- #
# Jacobian Validation
# ----------------------------------------------------------------------------- #

"""
    validate_drive_jacobian(d::NonlinearDrive, n_controls::Int; atol=1e-6, n_samples=3)

Spot-check the Jacobian of a `NonlinearDrive` against ForwardDiff at random control vectors.
Throws an `AssertionError` if the user-provided Jacobian disagrees with the AD Jacobian.

This is called automatically during `QuantumSystem` construction for all `NonlinearDrive`
terms, catching sign errors or off-by-one bugs early.
"""
function validate_drive_jacobian(
    d::NonlinearDrive,
    n_controls::Int;
    atol::Float64 = 1e-6,
    n_samples::Int = 3,
)
    for _ = 1:n_samples
        u = randn(n_controls)
        grad_ad = ForwardDiff.gradient(d.coeff, u)
        for j = 1:n_controls
            user_val = d.coeff_jac(u, j)
            @assert abs(user_val - grad_ad[j]) < atol (
                "NonlinearDrive Jacobian mismatch at u=$u, j=$j: " *
                "user=$(user_val), ForwardDiff=$(grad_ad[j])"
            )
        end
    end
end

"""
    validate_drive_hessian(d::NonlinearDrive, n_controls::Int; atol=1e-6, n_samples=3)

Spot-check the Hessian of a `NonlinearDrive` against ForwardDiff at random control vectors.
Throws an `AssertionError` if the user-provided Hessian disagrees with the AD Hessian.
"""
function validate_drive_hessian(
    d::NonlinearDrive,
    n_controls::Int;
    atol::Float64 = 1e-6,
    n_samples::Int = 3,
)
    for _ = 1:n_samples
        u = randn(n_controls)
        hess_ad = ForwardDiff.hessian(d.coeff, u)
        for i = 1:n_controls, j = i:n_controls
            user_val = d.coeff_hess(u, i, j)
            @assert abs(user_val - hess_ad[i, j]) < atol (
                "NonlinearDrive Hessian mismatch at u=$u, (i,j)=($i,$j): " *
                "user=$(user_val), ForwardDiff=$(hess_ad[i, j])"
            )
        end
    end
end

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

    # Hessian is always zero for linear drives
    @test drive_coeff_hess(d, u, 1, 1) == 0.0
    @test drive_coeff_hess(d, u, 2, 2) == 0.0
    @test drive_coeff_hess(d, u, 1, 2) == 0.0

    # active_controls
    @test active_controls(d) == [2]
end

@testitem "NonlinearDrive with explicit Jacobian" begin
    using Piccolo
    using SparseArrays

    H = sparse([1.0+0im 0.0+0im; 0.0+0im -1.0+0im])

    # Quadratic: u[1]^2 + u[2]^2
    d = NonlinearDrive(
        H,
        u -> u[1]^2 + u[2]^2,
        (u, j) -> j == 1 ? 2u[1] : j == 2 ? 2u[2] : 0.0,
    )

    u = [3.0, 4.0, 0.0]
    @test drive_coeff(d, u) == 25.0
    @test drive_coeff_jac(d, u, 1) == 6.0
    @test drive_coeff_jac(d, u, 2) == 8.0
    @test drive_coeff_jac(d, u, 3) == 0.0
    @test d.H == H
    @test isempty(active_controls(d))  # default: empty = all controls

    # Hessian: ∂²(u₁² + u₂²)/∂u_i∂u_j
    @test drive_coeff_hess(d, u, 1, 1) ≈ 2.0  # ∂²/∂u₁² = 2
    @test drive_coeff_hess(d, u, 2, 2) ≈ 2.0  # ∂²/∂u₂² = 2
    @test drive_coeff_hess(d, u, 1, 2) ≈ 0.0  # ∂²/∂u₁∂u₂ = 0
    @test drive_coeff_hess(d, u, 3, 3) ≈ 0.0  # u₃ not involved

    # Product: u[1] * u[2]
    d2 = NonlinearDrive(H, u -> u[1] * u[2], (u, j) -> j == 1 ? u[2] : j == 2 ? u[1] : 0.0)

    @test drive_coeff(d2, u) == 12.0
    @test drive_coeff_jac(d2, u, 1) == 4.0
    @test drive_coeff_jac(d2, u, 2) == 3.0

    # Hessian: ∂²(u₁u₂)/∂u_i∂u_j
    @test drive_coeff_hess(d2, u, 1, 2) ≈ 1.0  # ∂²/∂u₁∂u₂ = 1
    @test drive_coeff_hess(d2, u, 2, 1) ≈ 1.0  # symmetric
    @test drive_coeff_hess(d2, u, 1, 1) ≈ 0.0  # ∂²/∂u₁² = 0
    @test drive_coeff_hess(d2, u, 2, 2) ≈ 0.0  # ∂²/∂u₂² = 0

    # With active_controls
    d3 = NonlinearDrive(
        H,
        u -> u[1]^2 + u[2]^2,
        (u, j) -> j == 1 ? 2u[1] : j == 2 ? 2u[2] : 0.0;
        active_controls = [1, 2],
    )
    @test active_controls(d3) == [1, 2]
end

@testitem "NonlinearDrive auto-Jacobian" begin
    using Piccolo
    using SparseArrays

    H = sparse([1.0+0im 0.0+0im; 0.0+0im -1.0+0im])

    # 2-arg constructor: auto-Jacobian via ForwardDiff
    d = NonlinearDrive(H, u -> u[1]^2 + u[2]^2)

    u = [3.0, 4.0, 0.0]
    @test drive_coeff(d, u) == 25.0
    @test drive_coeff_jac(d, u, 1) ≈ 6.0
    @test drive_coeff_jac(d, u, 2) ≈ 8.0
    @test drive_coeff_jac(d, u, 3) ≈ 0.0

    # Product of controls
    d2 = NonlinearDrive(H, u -> u[1] * u[2])
    @test drive_coeff_jac(d2, [3.0, 4.0], 1) ≈ 4.0
    @test drive_coeff_jac(d2, [3.0, 4.0], 2) ≈ 3.0

    # Cubic term
    d3 = NonlinearDrive(H, u -> u[1]^3)
    @test drive_coeff_jac(d3, [2.0], 1) ≈ 12.0
    @test drive_coeff_hess(d3, [2.0], 1, 1) ≈ 12.0  # ∂²(u³)/∂u² = 6u = 12

    # Hessian of auto-generated quadratic
    @test drive_coeff_hess(d, [3.0, 4.0, 0.0], 1, 1) ≈ 2.0
    @test drive_coeff_hess(d, [3.0, 4.0, 0.0], 1, 2) ≈ 0.0

    # With active_controls
    d4 = NonlinearDrive(H, u -> u[2]^2; active_controls = [2])
    @test active_controls(d4) == [2]
    @test drive_coeff_jac(d4, [1.0, 3.0], 2) ≈ 6.0
    @test drive_coeff_hess(d4, [1.0, 3.0], 2, 2) ≈ 2.0
end

@testitem "NonlinearDrive accepts dense matrix" begin
    using Piccolo
    using SparseArrays

    H_dense = [1.0+0im 0.0+0im; 0.0+0im -1.0+0im]

    # 3-arg with dense matrix
    d = NonlinearDrive(H_dense, u -> u[1], (u, j) -> j == 1 ? 1.0 : 0.0)
    @test d.H isa SparseMatrixCSC
    @test drive_coeff(d, [0.5]) == 0.5

    # 2-arg with dense matrix
    d2 = NonlinearDrive(H_dense, u -> u[1]^2)
    @test d2.H isa SparseMatrixCSC
    @test drive_coeff(d2, [3.0]) == 9.0
    @test drive_coeff_jac(d2, [3.0], 1) ≈ 6.0
    @test drive_coeff_hess(d2, [3.0], 1, 1) ≈ 2.0
end

@testitem "NonlinearDrive with explicit Hessian" begin
    using Piccolo
    using SparseArrays

    H = sparse([1.0+0im 0.0+0im; 0.0+0im -1.0+0im])

    # Provide all three: coeff, jac, hess
    d = NonlinearDrive(
        H,
        u -> u[1]^2 + u[2]^2,
        (u, j) -> j == 1 ? 2u[1] : j == 2 ? 2u[2] : 0.0;
        coeff_hess = (u, i, j) -> (i == j && i <= 2) ? 2.0 : 0.0,
    )

    u = [3.0, 4.0]
    @test drive_coeff_hess(d, u, 1, 1) == 2.0
    @test drive_coeff_hess(d, u, 2, 2) == 2.0
    @test drive_coeff_hess(d, u, 1, 2) == 0.0

    # Cubic with explicit hess: u[1]^3
    d2 = NonlinearDrive(
        H,
        u -> u[1]^3,
        (u, j) -> j == 1 ? 3u[1]^2 : 0.0;
        coeff_hess = (u, i, j) -> (i == 1 && j == 1) ? 6u[1] : 0.0,
    )
    @test drive_coeff_hess(d2, [2.0], 1, 1) == 12.0
    @test drive_coeff_hess(d2, [5.0], 1, 1) == 30.0
end

@testitem "validate_drive_jacobian" begin
    using Piccolo
    using SparseArrays

    H = sparse([1.0+0im 0.0+0im; 0.0+0im -1.0+0im])

    # Correct Jacobian should pass
    d_correct = NonlinearDrive(
        H,
        u -> u[1]^2 + u[2]^2,
        (u, j) -> j == 1 ? 2u[1] : j == 2 ? 2u[2] : 0.0,
    )
    validate_drive_jacobian(d_correct, 2)  # should not throw

    # Wrong Jacobian should fail
    d_wrong = NonlinearDrive(
        H,
        u -> u[1]^2 + u[2]^2,
        (u, j) -> j == 1 ? u[1] : j == 2 ? u[2] : 0.0,  # missing factor of 2
    )
    @test_throws AssertionError validate_drive_jacobian(d_wrong, 2)

    # Auto-Jacobian should always pass validation
    d_auto = NonlinearDrive(H, u -> u[1]^2 + u[2]^2)
    validate_drive_jacobian(d_auto, 2)
end

@testitem "validate_drive_hessian" begin
    using Piccolo
    using SparseArrays

    H = sparse([1.0+0im 0.0+0im; 0.0+0im -1.0+0im])

    # Auto-Hessian should always pass
    d_auto = NonlinearDrive(H, u -> u[1]^2 + u[2]^2)
    validate_drive_hessian(d_auto, 2)

    # Correct explicit Hessian should pass
    d_correct = NonlinearDrive(
        H,
        u -> u[1]^2 + u[2]^2,
        (u, j) -> j == 1 ? 2u[1] : j == 2 ? 2u[2] : 0.0;
        coeff_hess = (u, i, j) -> (i == j && i <= 2) ? 2.0 : 0.0,
    )
    validate_drive_hessian(d_correct, 2)

    # Wrong explicit Hessian should fail
    d_wrong = NonlinearDrive(
        H,
        u -> u[1]^2 + u[2]^2,
        (u, j) -> j == 1 ? 2u[1] : j == 2 ? 2u[2] : 0.0;
        coeff_hess = (u, i, j) -> 0.0,  # wrong: should be 2.0 on diagonal
    )
    @test_throws AssertionError validate_drive_hessian(d_wrong, 2)
end

@testitem "G works on AbstractDrive types" begin
    using Piccolo
    using SparseArrays

    H = sparse(ComplexF64[0 1; 1 0])

    # G on LinearDrive should match G on the matrix
    ld = LinearDrive(H, 1)
    @test Piccolo.Isomorphisms.G(ld) == Piccolo.Isomorphisms.G(H)

    # G on NonlinearDrive
    nd = NonlinearDrive(H, u -> u[1]^2)
    @test Piccolo.Isomorphisms.G(nd) == Piccolo.Isomorphisms.G(H)

    # Broadcasting over Vector{AbstractDrive}
    drives = AbstractDrive[ld, nd]
    G_mats = Piccolo.Isomorphisms.G.(drives)
    @test length(G_mats) == 2
    @test all(G_mats .== Ref(Piccolo.Isomorphisms.G(H)))
end
