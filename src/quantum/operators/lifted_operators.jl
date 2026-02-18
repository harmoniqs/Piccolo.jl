module LiftedOperators

export lift_operator

using LinearAlgebra
using TestItems


@doc raw"""
    lift_operator(operator::AbstractMatrix{<:Number}, i::Int, subsystem_levels::Vector{Int})
    lift_operator(operator::AbstractMatrix{<:Number}, i::Int, n_qubits::Int; kwargs...)
    lift_operator(operators::AbstractVector{<:AbstractMatrix{T}}, indices::AbstractVector{Int}, subsystem_levels::Vector{Int})
    lift_operator(operators::AbstractVector{<:AbstractMatrix{T}}, indices::AbstractVector{Int}, n_qubits::Int; kwargs...)
    lift_operator(operator::AbstractMatrix{T}, indices::AbstractVector{Int}, subsystem_levels::AbstractVector{Int})
    lift_operator(operator::AbstractMatrix{T}, indices::AbstractVector{Int}, n_qubits::Int; kwargs...)

Lift an `operator` acting on the `i`-th subsystem within `subsystem_levels` to an operator
acting on the entire system spanning `subsystem_levels`.
"""
function lift_operator end

function lift_operator(
    operator::AbstractMatrix{T},
    i::Int,
    subsystem_levels::AbstractVector{Int},
) where {T<:Number}
    @assert size(operator, 1) == subsystem_levels[i] "Operator must match subsystem level."
    Is = [Matrix{T}(I(l)) for l ∈ subsystem_levels]
    Is[i] = operator
    return reduce(kron, Is)
end

function lift_operator(
    operator::AbstractMatrix{T},
    i::Int,
    n_qubits::Int;
    levels::Int = size(operator, 1),
) where {T<:Number}
    return lift_operator(operator, i, fill(levels, n_qubits))
end

function lift_operator(
    operators::AbstractVector{<:AbstractMatrix{T}},
    indices::AbstractVector{Int},
    subsystem_levels::AbstractVector{Int},
) where {T<:Number}
    @assert length(operators) == length(indices)
    return prod([
        lift_operator(op, i, subsystem_levels) for (op, i) ∈ zip(operators, indices)
    ])
end

function lift_operator(
    operators::AbstractVector{<:AbstractMatrix{T}},
    indices::AbstractVector{Int},
    n_qubits::Int;
    levels::Int = size(operators[1], 1),
) where {T<:Number}
    return prod([
        lift_operator(op, i, n_qubits, levels = levels) for
        (op, i) ∈ zip(operators, indices)
    ])
end

function lift_operator(
    operator::AbstractMatrix{T},
    indices::AbstractVector{Int},
    subsystem_levels::AbstractVector{Int},
) where {T<:Number}
    N = length(subsystem_levels)
    L = [subsystem_levels[i] for i ∈ indices]
    Lᶜ = [subsystem_levels[i] for i ∈ setdiff(1:N, indices)]
    @assert prod(L) == size(operator, 1)

    # Start with operator at the leading position
    shape = vcat(L, Lᶜ)
    array_shape = reverse(vcat(shape, shape))
    full_operator = kron(operator, [Matrix{T}(I(l)) for l ∈ Lᶜ]...)

    # Permute the array to match the actual subsystem order
    order = vcat(indices, [i for i ∈ setdiff(1:N, indices)])
    perm = sortperm(order)
    array_perm = reverse(2length(perm) + 1 .- vcat(perm, perm .+ length(perm)))

    return reshape(
        PermutedDimsArray(reshape(full_operator, array_shape...), array_perm),
        size(full_operator),
    )
end

function lift_operator(
    operator::AbstractMatrix{T},
    indices::AbstractVector{Int},
    n_qubits::Int;
    levels::Int = 2,
) where {T<:Number}
    return lift_operator(operator, indices, fill(levels, n_qubits))
end

# ****************************************************************************** #

@testitem "lift_operator subsystems" begin
    using LinearAlgebra
    @test lift_operator(PAULIS.X, 1, [2, 3]) ≈ kron(PAULIS.X, I(3))
    @test lift_operator(PAULIS.Y, 2, [4, 2]) ≈ kron(I(4), PAULIS.Y)
    @test lift_operator(PAULIS.X, 2, [3, 2, 4]) ≈ reduce(kron, [I(3), PAULIS.X, I(4)])
end

@testitem "lift_operator qubits" begin
    using LinearAlgebra
    @test lift_operator(PAULIS.X, 1, 2) ≈ kron(PAULIS.X, I(2))
    @test lift_operator(PAULIS.Y, 2, 2) ≈ kron(I(2), PAULIS.Y)
    @test lift_operator(PAULIS.X, 2, 3) ≈ reduce(kron, [I(2), PAULIS.X, I(2)])
end

@testitem "lift_operator multiple operators" begin
    using LinearAlgebra
    pair = [PAULIS.X, PAULIS.Y]
    @test lift_operator(pair, [1, 2], [2, 2]) ≈ kron(PAULIS.X, PAULIS.Y)
    @test lift_operator(pair, [2, 1], [2, 2]) ≈ kron(PAULIS.Y, PAULIS.X)
    @test lift_operator(pair, [1, 2], [2, 2, 3]) ≈ kron(PAULIS.X, PAULIS.Y, I(3))
    @test lift_operator(pair, [2, 3], [4, 2, 2]) ≈ kron(I(4), PAULIS.X, PAULIS.Y)

    # Pass number of qubits
    @test lift_operator(pair, [1, 2], 3) ≈ kron(PAULIS.X, PAULIS.Y, I(2))
end

@testitem "lift_operator single operator into disjoint levels" begin
    U = haar_random(2)
    UU = kron(U, U)
    UUU = kron(UU, U)
    I2 = [1 0; 0 1]
    I3 = [1 0 0; 0 1 0; 0 0 1]

    @test lift_operator(UU, [1, 2], [2, 2, 2]) ≈ kron(U, U, I2)
    @test lift_operator(UU, [2, 3], [2, 2, 2]) ≈ kron(I2, U, U)
    @test lift_operator(UU, [1, 3], [2, 2, 2]) ≈ kron(U, I2, U)

    @test lift_operator(UU, [1], [4, 2, 3]) ≈ kron(UU, I2, I3)
    @test lift_operator(UU, [2], [2, 4, 3]) ≈ kron(I2, UU, I3)
    @test lift_operator(UU, [3], [2, 3, 4]) ≈ kron(I2, I3, UU)

    @test lift_operator(UUU, [1, 3, 4], [2, 3, 2, 2]) ≈ kron(U, I3, UU)
    @test lift_operator(UUU, [1, 2, 4], [2, 2, 3, 2]) ≈ kron(UU, I3, U)

    # Test qubit interface
    @test lift_operator(U, [1], 3) ≈ kron(U, I2, I2)
    @test lift_operator(UU, [2, 3], 3) ≈ kron(I2, U, U)
    @test lift_operator(UU, [1, 3], 3) ≈ kron(U, I2, U)
end


end
