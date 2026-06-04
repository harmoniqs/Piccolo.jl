module DualRailEncodings

export DualRailEncoding
export logical_basis_states
export reduce_to_subspace
export subspace_transform
export target_states

using ..Gates

using LinearAlgebra
using SparseArrays
using TestItems

@doc raw"""
    DualRailEncoding(; n_qubits, levels_per_rail, conservation=:exact_N, N=n_qubits)

Describe a dual-rail encoding of `n_qubits` logical qubits in `2n_qubits`
physical modes. Logical zero and one are represented by occupations `|1,0>` and
`|0,1>`, respectively.

Use [`subspace_transform`](@ref) to reduce full-space Hamiltonians, then construct
a system from the reduced Hamiltonians and pass `EmbeddedOperator(:CX, enc)` as
the goal of a `UnitaryTrajectory`.
"""
struct DualRailEncoding
    n_qubits::Int
    levels_per_rail::Int
    n_rails::Int
    subspace_levels::Vector{Int}
    conservation::Symbol
    N::Int

    function DualRailEncoding(;
        n_qubits::Int,
        levels_per_rail::Int,
        conservation::Symbol = :exact_N,
        N::Int = n_qubits,
    )
        n_qubits > 0 || throw(ArgumentError("n_qubits must be positive."))
        levels_per_rail >= 2 || throw(ArgumentError("levels_per_rail must be at least 2."))
        conservation in (:exact_N, :upto_N) ||
            throw(ArgumentError("conservation must be :exact_N or :upto_N."))
        N >= 0 || throw(ArgumentError("N must be non-negative."))
        N <= 2n_qubits * (levels_per_rail - 1) ||
            throw(ArgumentError("N exceeds the maximum excitation number."))

        n_rails = 2n_qubits
        return new(
            n_qubits,
            levels_per_rail,
            n_rails,
            fill(levels_per_rail, n_rails),
            conservation,
            N,
        )
    end
end

function _occupation_index(occupations::AbstractVector{Int}, levels::AbstractVector{Int})
    index = 1
    stride = 1
    for rail in Iterators.reverse(eachindex(levels))
        index += occupations[rail] * stride
        stride *= levels[rail]
    end
    return index
end

"""
    subspace_transform(enc::DualRailEncoding)

Return the sparse isometry from the selected excitation sector to the full
tensor-product space and the corresponding full-space basis indices.
"""
function subspace_transform(enc::DualRailEncoding)
    full_dimension = prod(enc.subspace_levels)
    good_indices = Int[]

    for linear_index = 0:(full_dimension-1)
        remainder = linear_index
        excitation_number = 0
        for level in Iterators.reverse(enc.subspace_levels)
            excitation_number += remainder % level
            remainder = div(remainder, level)
        end
        keep =
            enc.conservation === :exact_N ? excitation_number == enc.N :
            excitation_number <= enc.N
        keep && push!(good_indices, linear_index + 1)
    end

    T = sparse(
        good_indices,
        eachindex(good_indices),
        ones(ComplexF64, length(good_indices)),
        full_dimension,
        length(good_indices),
    )
    return T, good_indices
end

function reduce_to_subspace(operator::AbstractMatrix, enc::DualRailEncoding)
    T, _ = subspace_transform(enc)
    return T' * operator * T
end

function reduce_to_subspace(state::AbstractVector, enc::DualRailEncoding)
    T, _ = subspace_transform(enc)
    return T' * state
end

function _logical_basis_indices(enc::DualRailEncoding)
    indices = Int[]
    occupations = zeros(Int, enc.n_rails)
    for logical_index = 0:(2^enc.n_qubits-1)
        fill!(occupations, 0)
        for qubit = 1:enc.n_qubits
            bit = (logical_index >> (enc.n_qubits - qubit)) & 1
            occupations[2qubit-1+bit] = 1
        end
        push!(indices, _occupation_index(occupations, enc.subspace_levels))
    end
    return indices
end

"""
    logical_basis_states(enc::DualRailEncoding)

Return full physical-space kets for the computational basis of the logical
register, ordered from `|0...0>` through `|1...1>`.
"""
function logical_basis_states(enc::DualRailEncoding)
    full_dimension = prod(enc.subspace_levels)
    return [
        sparsevec([index], ComplexF64[1], full_dimension) for
        index in _logical_basis_indices(enc)
    ]
end

"""
    target_states(gate::Symbol, enc::DualRailEncoding)

Apply a logical gate to every logical computational-basis input and encode the
result in the full physical space.
"""
function target_states(gate::Symbol, enc::DualRailEncoding)
    gate in keys(GATES) || throw(ArgumentError("Unknown gate: $gate."))
    logical_gate = GATES[gate]
    logical_dimension = 2^enc.n_qubits
    size(logical_gate) == (logical_dimension, logical_dimension) || throw(
        DimensionMismatch(
            "Gate $gate acts on $(size(logical_gate, 1)) states, " *
            "but the encoding has $logical_dimension logical states.",
        ),
    )

    encoded_basis = logical_basis_states(enc)
    return [
        sum(logical_gate[row, column] * encoded_basis[row] for row = 1:logical_dimension)
        for column = 1:logical_dimension
    ]
end

@testitem "Dual-rail encoding" begin
    using LinearAlgebra
    using SparseArrays

    enc = DualRailEncoding(n_qubits = 2, levels_per_rail = 2)
    T, indices = subspace_transform(enc)
    @test T isa SparseMatrixCSC{ComplexF64}
    @test T' * T == I(size(T, 2))
    @test nnz(T) == 6

    full_operator = spdiagm(0 => ComplexF64.(1:size(T, 1)))
    @test reduce_to_subspace(full_operator, enc) == T' * full_operator * T
    full_state = sparsevec([indices[1]], ComplexF64[1], size(T, 1))
    @test reduce_to_subspace(full_state, enc) == T' * full_state

    basis = logical_basis_states(enc)
    targets = target_states(:CX, enc)
    @test targets == basis[[1, 2, 4, 3]]

    open_enc = DualRailEncoding(n_qubits = 2, levels_per_rail = 2, conservation = :upto_N)
    _, open_indices = subspace_transform(open_enc)
    @test issubset(indices, open_indices)
    @test length(open_indices) > length(indices)
end

@testitem "Three-qubit dual-rail targets" begin
    using SparseArrays

    single_qubit_enc = DualRailEncoding(n_qubits = 1, levels_per_rail = 2)
    for gate in (:I, :X, :Y, :Z, :H)
        @test length(target_states(gate, single_qubit_enc)) == 2
    end

    two_qubit_enc = DualRailEncoding(n_qubits = 2, levels_per_rail = 2)
    for gate in (:CX, :CZ)
        @test length(target_states(gate, two_qubit_enc)) == 4
    end

    enc = DualRailEncoding(n_qubits = 3, levels_per_rail = 2)
    basis = logical_basis_states(enc)
    T, _ = subspace_transform(enc)
    @test nnz(T) == binomial(6, 3)
    @test target_states(:CCX, enc) == basis[[1, 2, 3, 4, 5, 6, 8, 7]]
    @test target_states(:CCZ, enc) == [basis[1:7]..., -basis[8]]
end

end
