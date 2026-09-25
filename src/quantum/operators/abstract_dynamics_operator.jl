# Slice 3b (#430): moved verbatim from Piccolissimo's Operators module — the shared
# operator seam the dense spline cells build on. Exports here mirror the
# original module-level export block for the moved surface.
export AbstractDynamicsOperator
export apply!, state_dim, materialize, _to_operator
export dynamics_operator

"""
    AbstractDynamicsOperator

Abstract supertype for operators in quantum dynamics.

All concrete subtypes must implement:
- `apply!(y, op, x, α, β)` — compute `y = α*(op*x) + β*y` in-place
- `state_dim(op)` — return the input/output dimension
- `Base.eltype(op)` — return the element type

Optional interface:
- `materialize(op)` — return a dense/sparse `AbstractMatrix` representation
"""
abstract type AbstractDynamicsOperator end

"""
    apply!(y, op::AbstractDynamicsOperator, x, α, β)

Compute `y = α * (op * x) + β * y` in-place.

This is the fundamental operation for all dynamics operators, matching
the 5-argument `mul!` BLAS interface.
"""
function apply! end

"""
    state_dim(op::AbstractDynamicsOperator) → Int

Return the input/output dimension of the operator.
"""
function state_dim end

"""
    materialize(op::AbstractDynamicsOperator) → AbstractMatrix

Return an explicit matrix representation of the operator.
Falls back to error for operators that cannot be materialized.
"""
function materialize end
materialize(op::AbstractDynamicsOperator) =
    error("materialize not implemented for $(typeof(op))")

# ── Matrix interop ──────────────────────────────────────────────────────────

"""
    Base.size(op::AbstractDynamicsOperator)

Return `(n, n)` where `n = state_dim(op)`, for compatibility with code that
calls `size(H, 1)` on drive Hamiltonians.
"""
Base.size(op::AbstractDynamicsOperator) = (state_dim(op), state_dim(op))
Base.size(op::AbstractDynamicsOperator, d::Int) = state_dim(op)

"""
    Base.Matrix(op::AbstractDynamicsOperator)

Materialize the operator as a dense `Matrix`. Used by `_ensure_matrix` in
Piccolo.jl when a plain matrix is needed (e.g., for `exp(G)` in Magnus).
"""
Base.Matrix(op::AbstractDynamicsOperator) = Matrix(materialize(op))

# ── Operator→operator conversion ────────────────────────────────────────────

"""
    _to_operator(H)

Convert a drive Hamiltonian to an `AbstractDynamicsOperator` for use in the
sensitivity ODE. If `H` is already an operator, return it directly.
If `H` is an `AbstractMatrix`, wrap it in a `MatrixOperator`.
"""
_to_operator(H::AbstractDynamicsOperator) = H
# _to_operator(H::AbstractMatrix) is defined in matrix_operator.jl after MatrixOperator

# ── LinearAlgebra.mul! adapters ─────────────────────────────────────────────
#
# Make every `AbstractDynamicsOperator` subtype usable wherever `mul!` is
# expected — most importantly in `ExponentialAction.expv` / `KrylovKit`'s
# Krylov methods, which only require `mul!(y, op, x)`. This is the foundation
# for matrix-free `expv` on structured Hamiltonians (e.g., bosonic
# `LadderOperator`, multi-qubit `KronIdentityOperator`).

# Note: `ExponentialAction.expv` defaults to `shift=true` which requires
# `tr(A)` and `A - I·shift` (UniformScaling subtraction). Pass `shift=false`
# when calling `expv` on an `AbstractDynamicsOperator` — the `shift_matrix`
# detour is a no-op for matrix-free operators anyway.

LinearAlgebra.mul!(y::AbstractVector, op::AbstractDynamicsOperator, x::AbstractVector) =
    apply!(y, op, x, true, false)
LinearAlgebra.mul!(
    y::AbstractVector,
    op::AbstractDynamicsOperator,
    x::AbstractVector,
    α,
    β,
) = apply!(y, op, x, α, β)
LinearAlgebra.mul!(Y::AbstractMatrix, op::AbstractDynamicsOperator, X::AbstractMatrix) =
    apply!(Y, op, X, true, false)
LinearAlgebra.mul!(
    Y::AbstractMatrix,
    op::AbstractDynamicsOperator,
    X::AbstractMatrix,
    α,
    β,
) = apply!(Y, op, X, α, β)

# `*` overload for convenience — useful for one-off uses outside hot paths.
function Base.:*(op::AbstractDynamicsOperator, x::AbstractVector)
    y = similar(x, promote_type(eltype(op), eltype(x)), state_dim(op))
    return mul!(y, op, x)
end
function Base.:*(op::AbstractDynamicsOperator, X::AbstractMatrix)
    Y = similar(X, promote_type(eltype(op), eltype(X)), state_dim(op), size(X, 2))
    return mul!(Y, op, X)
end

# ── Thread-safe copying ─────────────────────────────────────────────────────

"""
    Base.copy(op::AbstractDynamicsOperator)

Return a copy of the operator with independent mutable state. Operators with
internal buffers (e.g., `KroneckerOperator`) need this for thread-safe use in
parallel sensitivity ODE closures.
"""
Base.copy(op::AbstractDynamicsOperator) = op  # default: immutable operators return self

# AbstractDrive ↔ AbstractDynamicsOperator bridge
# ---------------------------------------------------------------------------- #

using ..Quantum.QuantumSystems: AbstractDrive, ModulatedDrive

"""
    dynamics_operator(d::AbstractDrive) → AbstractDynamicsOperator

Unwrap the operator payload of an `AbstractDrive` into an `AbstractDynamicsOperator`
suitable for use with the matrix-free `apply!` interface.

- If `d.H isa AbstractDynamicsOperator`, return it unchanged (preserves matrix-free
  structure — e.g. `DiagonalOperator`, `SumOperator`, `KroneckerOperator`).
- If `d.H isa AbstractMatrix`, wrap it with `MatrixOperator(d.H)`.
- Any other payload raises an error: integrators that need a richer treatment
  (e.g. `NonlinearDrive` with non-affine coefficients) must dispatch to a
  dedicated path.

This is the canonical way for spline / exponential integrators to consume
`sys.H_drives::Vector{AbstractDrive}` without materializing dense matrices.
"""
dynamics_operator(d::AbstractDrive) = _dyn_op(d.H)
# ModulatedDrive carries no `.H` of its own — everything delegates to the base
# drive (the drive_matrix/drive_dim/G convention), and the operator bridge does too.
dynamics_operator(d::ModulatedDrive) = dynamics_operator(d.base)

_dyn_op(H::AbstractDynamicsOperator) = H
_dyn_op(H::AbstractMatrix) = MatrixOperator(H)
_dyn_op(::Any) = error("""
                       dynamics_operator only supports drives whose `.H` field is an AbstractMatrix
                       or an AbstractDynamicsOperator. NonlinearDrive or other drive variants with
                       non-affine coefficient structure must be handled by a dedicated integrator
                       path (not yet implemented).
                       """)

# ── Tests ───────────────────────────────────────────────────────────────────

@testitem "AbstractDynamicsOperator: LinearAlgebra interop via MatrixOperator" begin
    using Piccolo
    using LinearAlgebra

    H = ComplexF64[0 1; 1 0]
    op = MatrixOperator(H)
    x = ComplexF64[1.0, 0.0]

    # 5-arg mul! delegates to apply! (the BLAS-style interface)
    y = zeros(ComplexF64, 2)
    mul!(y, op, x, -im, false)
    @test y ≈ -im * (H * x)

    # 3-arg mul! overwrites: y = op * x
    mul!(y, op, x)
    @test y ≈ H * x

    # α/β accumulation: y = α * (op * x) + β * y   (y holds H*x going in)
    mul!(y, op, x, 2.0, 3.0)
    @test y ≈ 5 * (H * x)

    # Matrix operands go through the same adapter
    X = ComplexF64[1.0 0.0; 0.0 0.0]
    Y = zeros(ComplexF64, 2, 2)
    mul!(Y, op, X)
    @test Y ≈ H * X

    # Base surface: size (both forms), Matrix, *, copy, eltype, state_dim, materialize
    @test size(op) == (2, 2)
    @test size(op, 1) == 2
    @test Matrix(op) == H
    @test op * x ≈ H * x
    @test op * X ≈ H * X
    @test copy(op) === op      # immutable operators return self by contract
    @test eltype(op) == ComplexF64
    @test state_dim(op) == 2
    @test Piccolo.Quantum.materialize(op) == H
    @test Piccolo.Quantum._to_operator(op) === op
    @test Piccolo.Quantum._to_operator(H) isa MatrixOperator
end

@testitem "dynamics_operator bridges AbstractDrive payloads to operators" begin
    using Piccolo
    using LinearAlgebra
    using SparseArrays

    H = sparse(ComplexF64[0 1; 1 0])

    # AbstractDrive with a matrix payload wraps into a MatrixOperator
    ld = LinearDrive(H, 1)
    op = dynamics_operator(ld)
    @test op isa MatrixOperator
    @test Piccolo.Quantum.materialize(op) == H
    y = zeros(ComplexF64, 2)
    mul!(y, op, ComplexF64[1.0, 0.0])
    @test y ≈ H * ComplexF64[1.0, 0.0]

    # ModulatedDrive unwraps to its base drive's payload
    md = ModulatedDrive(ld, t -> cos(t))
    @test Piccolo.Quantum.materialize(dynamics_operator(md)) == H
end
