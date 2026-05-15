module AbstractOperators

export AbstractDynamicsOperator

"""
    AbstractDynamicsOperator

Abstract supertype for structured operator types that participate in quantum
dynamics. Concrete subtypes (defined in Piccolissimo: `DiagonalOperator`,
`LadderOperator`, `KronIdentityOperator`, `KroneckerOperator`, `FockBanded`,
`SumOperator`, `ScaledOperator`, `MatrixOperator`, ...) implement
matrix-free `apply!(y, op, x, α, β)` kernels that the spline integrator
dispatches on to avoid generic sparse `mul!` in the per-substep hot loop.

This file in Piccolo declares the abstract supertype only, so
`QuantumSystem` constructor signatures can accept structured payloads as
`H_drift` (and drive `.H` fields can hold them too) without Piccolo needing
to depend on Piccolissimo. The concrete implementations and their algebra
(`+`, `*`, `materialize`, `apply!`) live in Piccolissimo where the
integrator that consumes them also lives.

The contract for concrete subtypes is:
- `apply!(y, op, x, α, β)` — in-place 5-arg multiply
- `Base.size(op)` — `(n, n)` where `n` is the state dimension
- `Base.Matrix(op)` — fallback dense materialization (used by
  `_ensure_matrix(H)` and Hermiticity checks)
"""
abstract type AbstractDynamicsOperator end

end
