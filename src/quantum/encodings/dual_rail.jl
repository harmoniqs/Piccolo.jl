module DualRailEncodings

export DualRailEncoding
export subspace_transform
export reduce_to_subspace
export logical_basis_states
export target_states

using ..Gates

using LinearAlgebra
using SparseArrays

@doc raw"""
    DualRailEncoding(; n_qubits, levels_per_rail, conservation=:exact_N, N=n_qubits)

Describe a dual-rail logical encoding where each logical qubit is represented by
two physical rails. Logical `0` occupies the first rail in a pair and logical `1`
occupies the second rail.

# Example

```julia
enc = DualRailEncoding(n_qubits=2, levels_per_rail=2, conservation=:exact_N, N=2)
T, idxs = subspace_transform(enc)
goal = EmbeddedOperator(:CX, enc)
```
"""
struct DualRailEncoding
    n_qubits::Int
    levels_per_rail::Int
    n_rails::Int
    subspace_levels::Vector{Int}
    conservation::Symbol
    N::Int
end

function DualRailEncoding(;
    n_qubits::Int,
    levels_per_rail::Int,
    conservation::Symbol = :exact_N,
    N::Int = n_qubits,
)
    n_qubits > 0 || throw(ArgumentError("n_qubits must be positive."))
    levels_per_rail > 1 || throw(ArgumentError("levels_per_rail must be at least 2."))
    conservation in (:exact_N, :upto_N) ||
        throw(ArgumentError("conservation must be :exact_N or :upto_N."))
    N ≥ 0 || throw(ArgumentError("N must be non-negative."))

    n_rails = 2n_qubits
    return DualRailEncoding(
        n_qubits,
        levels_per_rail,
        n_rails,
        fill(levels_per_rail, n_rails),
        conservation,
        N,
    )
end

function _occupations(index::Int, levels::AbstractVector{Int})
    offset = index - 1
    occupations = Vector{Int}(undef, length(levels))

    for i in eachindex(levels)
        stride = prod(@view levels[(i+1):end]; init = 1)
        occupations[i] = offset ÷ stride
        offset %= stride
    end

    return occupations
end

function _basis_index(occupations::AbstractVector{Int}, levels::AbstractVector{Int})
    length(occupations) == length(levels) ||
        throw(ArgumentError("occupations and levels must have the same length."))

    index = 1
    for i in eachindex(levels)
        0 ≤ occupations[i] < levels[i] ||
            throw(ArgumentError("occupation $(occupations[i]) is outside rail level $(levels[i])."))
        stride = prod(@view levels[(i+1):end]; init = 1)
        index += occupations[i] * stride
    end
    return index
end

function _subspace_indices(enc::DualRailEncoding)
    levels = enc.subspace_levels
    indices = Int[]

    for index in 1:prod(levels)
        excitation_count = sum(_occupations(index, levels))
        if enc.conservation == :exact_N
            excitation_count == enc.N && push!(indices, index)
        else
            excitation_count ≤ enc.N && push!(indices, index)
        end
    end

    return indices
end

"""
    subspace_transform(enc) -> (T, idxs)

Return the sparse transform `T` from encoded subspace coordinates to full-space
coordinates plus the full-space basis indices in the selected excitation sector.
"""
function subspace_transform(enc::DualRailEncoding)
    indices = _subspace_indices(enc)
    rows = indices
    cols = collect(1:length(indices))
    vals = ones(ComplexF64, length(indices))
    transform = sparse(rows, cols, vals, prod(enc.subspace_levels), length(indices))
    return transform, indices
end

function reduce_to_subspace(operator::AbstractMatrix, enc::DualRailEncoding)
    transform, _ = subspace_transform(enc)
    return transform' * operator * transform
end

function reduce_to_subspace(state::AbstractVector, enc::DualRailEncoding)
    transform, _ = subspace_transform(enc)
    return transform' * state
end

function _logical_bits(logical_index::Int, n_qubits::Int)
    value = logical_index - 1
    bits = Vector{Int}(undef, n_qubits)

    for i in 1:n_qubits
        stride = 2^(n_qubits - i)
        bits[i] = value ÷ stride
        value %= stride
    end

    return bits
end

function logical_subspace_indices(enc::DualRailEncoding)
    indices = Int[]

    for logical_index in 1:(2^enc.n_qubits)
        bits = _logical_bits(logical_index, enc.n_qubits)
        occupations = zeros(Int, enc.n_rails)

        for (qubit, bit) in enumerate(bits)
            first_rail = 2qubit - 1
            occupations[first_rail + bit] = 1
        end

        push!(indices, _basis_index(occupations, enc.subspace_levels))
    end

    return indices
end

function logical_basis_states(enc::DualRailEncoding)
    states = Vector{Vector{ComplexF64}}()
    full_dimension = prod(enc.subspace_levels)

    for index in logical_subspace_indices(enc)
        state = zeros(ComplexF64, full_dimension)
        state[index] = 1
        push!(states, state)
    end

    return states
end

function target_states(gate::Symbol, enc::DualRailEncoding)
    haskey(GATES, gate) ||
        throw(ArgumentError("Operator must be a valid gate. See Piccolo.GATES."))
    return target_states(GATES[gate], enc)
end

function target_states(gate::AbstractMatrix{<:Number}, enc::DualRailEncoding)
    expected_dimension = 2^enc.n_qubits
    size(gate) == (expected_dimension, expected_dimension) ||
        throw(DimensionMismatch("gate must act on all $(enc.n_qubits) logical qubits."))

    logical_states = logical_basis_states(enc)
    targets = Vector{Vector{ComplexF64}}()

    for col in 1:expected_dimension
        target = zeros(ComplexF64, length(logical_states[1]))
        for row in 1:expected_dimension
            coefficient = gate[row, col]
            !iszero(coefficient) && (target .+= coefficient .* logical_states[row])
        end
        push!(targets, target)
    end

    return targets
end

end
