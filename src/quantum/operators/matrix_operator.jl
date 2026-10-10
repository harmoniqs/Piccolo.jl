# Slice 3b (#430): moved verbatim from Piccolissimo's Operators module.
export MatrixOperator

"""
    MatrixOperator{T<:AbstractMatrix} <: AbstractDynamicsOperator

Wraps an existing matrix as an operator. `apply!` dispatches to `mul!`.

This is the migration bridge: every `H_drift`, `H_drive` matrix from
`QuantumSystem` becomes a `MatrixOperator` with zero behavior change.

# Example
```julia
H = [0.0+0im 1.0+0im; 1.0+0im 0.0+0im]
op = MatrixOperator(H)
y = zeros(ComplexF64, 2)
x = [1.0+0im, 0.0+0im]
apply!(y, op, x, -im, false)  # y = -i * H * x
```
"""
struct MatrixOperator{T<:AbstractMatrix} <: AbstractDynamicsOperator
    M::T
end

apply!(y, op::MatrixOperator, x, α, β) = mul!(y, op.M, x, α, β)

state_dim(op::MatrixOperator) = size(op.M, 1)

Base.eltype(op::MatrixOperator) = eltype(op.M)

materialize(op::MatrixOperator) = op.M

_to_operator(H::AbstractMatrix) = MatrixOperator(H)
