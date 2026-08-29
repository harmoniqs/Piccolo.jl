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

# ── Tests ───────────────────────────────────────────────────────────────────

@testitem "AbstractDynamicsOperator mul!/* adapters" begin
    using Piccolissimo: apply!, state_dim
    using LinearAlgebra

    # Build a small DiagonalOperator and verify mul!/* delegate through apply!.
    # These adapters make `AbstractDynamicsOperator` usable wherever `mul!` is
    # the abstraction (KrylovKit, IterativeSolvers, our own AugmentedAction
    # struct, etc.). Note: `ExponentialAction.expv` additionally needs `tr`,
    # `opnorm`, and `-(::Op, ::UniformScaling)` for its `shift=true` Krylov
    # normalization — those aren't in scope here.
    d = 4
    diag_vals = ComplexF64[1.0, 2.0, 3.0, 4.0]
    op = Piccolissimo.Operators.DiagonalOperator(diag_vals)

    @test state_dim(op) == d
    @test size(op) == (d, d)

    x = ComplexF64.(randn(d) .+ im .* randn(d))
    y_ref = diag_vals .* x

    # 3-arg mul!
    y = similar(x)
    mul!(y, op, x)
    @test isapprox(y, y_ref; atol = 1e-12)

    # 5-arg mul! (α, β)
    y2 = ComplexF64.(randn(d) .+ im .* randn(d))
    y2_ref = 2.0 * y_ref + 0.5 * y2
    mul!(y2, op, x, 2.0, 0.5)
    @test isapprox(y2, y2_ref; atol = 1e-12)

    # `*` overload
    @test isapprox(op * x, y_ref; atol = 1e-12)

    # Matrix-form mul! via apply!
    n = 3
    X = ComplexF64.(randn(d, n) .+ im .* randn(d, n))
    Y = similar(X)
    mul!(Y, op, X)
    @test isapprox(Y, diag_vals .* X; atol = 1e-12)

    println("✓ AbstractDynamicsOperator mul!/* adapter test passed")
end

# ---------------------------------------------------------------------------- #
# AbstractDrive ↔ AbstractDynamicsOperator bridge
# ---------------------------------------------------------------------------- #

using ..Quantum.QuantumSystems: AbstractDrive

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

_dyn_op(H::AbstractDynamicsOperator) = H
_dyn_op(H::AbstractMatrix) = MatrixOperator(H)
_dyn_op(::Any) = error("""
                       dynamics_operator only supports drives whose `.H` field is an AbstractMatrix
                       or an AbstractDynamicsOperator. NonlinearDrive or other drive variants with
                       non-affine coefficient structure must be handled by a dedicated integrator
                       path (not yet implemented).
                       """)
