module DualRailEncodings

export DualRailEncoding
export subspace_transform
export reduce_to_subspace
export logical_basis_states
export target_states

using ..Gates

using LinearAlgebra
using SparseArrays

# ----------------------------------------------------------------------------- #
#                              Dual-rail encoding                                #
# ----------------------------------------------------------------------------- #

@doc raw"""
    DualRailEncoding(;
        n_qubits::Int,
        levels_per_rail::Int = 2,
        conservation::Symbol = :exact_N,
        N::Int = n_qubits,
    )

Dual-rail encoding of `n_qubits` logical qubits into `n_rails = 2 * n_qubits` physical
modes ("rails"), where logical qubit `q` lives on the rail pair `(2q - 1, 2q)`. With
`m = N ÷ n_qubits` excitations per logical qubit, the logical basis states are

    |0⟩_q ↔ |m, 0⟩,    |1⟩_q ↔ |0, m⟩

on that rail pair (`m = 1` for the typical single-excitation encoding). The pattern
applies to any excitation-conserving system — bosonic, spin, or fermionic modes.

The encoded register lives in an excitation-number subspace of the full
tensor-product space spanned by `subspace_levels`:
- `conservation = :exact_N` keeps the closed-system sector `Σᵢ nᵢ = N`.
- `conservation = :upto_N` keeps `Σᵢ nᵢ ≤ N`, which is required for open / lossy
  systems where loss channels couple the conserved sector to lower-excitation sectors.

# Fields
- `n_qubits::Int`: number of logical qubits.
- `levels_per_rail::Int`: local Hilbert space dimension of each rail.
- `n_rails::Int`: number of physical modes, `2 * n_qubits`.
- `subspace_levels::Vector{Int}`: full tensor-product dimensions, length `n_rails`.
- `conservation::Symbol`: `:exact_N` (closed) or `:upto_N` (open / lossy).
- `N::Int`: target total excitation number (typically `n_qubits`).

# Example

The encoding intentionally does not carry a Hamiltonian — build a `QuantumSystem`
from your own `H_drift` / `H_drives`, then use the encoding to construct goals:

```julia
using Piccolo

enc = DualRailEncoding(n_qubits = 2, levels_per_rail = 2, conservation = :exact_N, N = 2)

T, idxs = subspace_transform(enc)        # sparse projector + "good" basis indices
H_sub = reduce_to_subspace(H_full, enc)  # user-built Hamiltonian, reduced

sys = QuantumSystem(H_drift, H_drives, drive_bounds)
goal = EmbeddedOperator(:CX, enc)        # encoding-aware goal
qtraj = UnitaryTrajectory(sys, pulse, goal)
```
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
        levels_per_rail::Int = 2,
        conservation::Symbol = :exact_N,
        N::Int = n_qubits,
    )
        n_qubits ≥ 1 || throw(ArgumentError("n_qubits must be at least 1, got $n_qubits"))
        levels_per_rail ≥ 2 ||
            throw(ArgumentError("levels_per_rail must be at least 2, got $levels_per_rail"))
        conservation ∈ (:exact_N, :upto_N) || throw(
            ArgumentError(
                "conservation must be :exact_N (closed system) or :upto_N " *
                "(open / lossy system), got :$conservation",
            ),
        )
        N % n_qubits == 0 || throw(
            ArgumentError(
                "N = $N must be divisible by n_qubits = $n_qubits so that each " *
                "logical qubit carries an integer number of excitations",
            ),
        )
        m = N ÷ n_qubits
        1 ≤ m ≤ levels_per_rail - 1 || throw(
            ArgumentError(
                "each logical qubit carries m = N ÷ n_qubits = $m excitations, " *
                "which must satisfy 1 ≤ m ≤ levels_per_rail - 1 = $(levels_per_rail - 1)",
            ),
        )
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

# ----------------------------------------------------------------------------- #
#                          Occupation-number bookkeeping                         #
# ----------------------------------------------------------------------------- #

# Map a 0-based occupation vector to its 1-based index in the kron basis, where the
# first subsystem is the most significant digit (Julia's kron convention, matching
# `basis_labels` in EmbeddedOperators).
function _occupations_to_index(
    occupations::AbstractVector{Int},
    levels::AbstractVector{Int},
)
    index = 0
    for (n, l) ∈ zip(occupations, levels)
        index = index * l + n
    end
    return index + 1
end

# Inverse of `_occupations_to_index`: decode a 1-based kron-basis index into the
# 0-based occupation of each subsystem.
function _index_to_occupations(index::Int, levels::AbstractVector{Int})
    k = index - 1
    occupations = Vector{Int}(undef, length(levels))
    for r ∈ reverse(eachindex(levels))
        occupations[r] = k % levels[r]
        k ÷= levels[r]
    end
    return occupations
end

# Ascending full-space indices of the encoding's excitation sector.
function _subspace_indices(enc::DualRailEncoding)
    levels = enc.subspace_levels
    d_full = prod(levels)
    if enc.conservation == :exact_N
        return [k for k ∈ 1:d_full if sum(_index_to_occupations(k, levels)) == enc.N]
    else # :upto_N (enforced by the constructor)
        return [k for k ∈ 1:d_full if sum(_index_to_occupations(k, levels)) ≤ enc.N]
    end
end

# ----------------------------------------------------------------------------- #
#                              Subspace projection                               #
# ----------------------------------------------------------------------------- #

@doc raw"""
    subspace_transform(enc::DualRailEncoding) -> (T, idxs)

Sparse projector `T::SparseMatrixCSC{ComplexF64}` mapping subspace coordinates to
full-space coordinates, `|ψ_full⟩ = T * |ψ_sub⟩`, together with the ascending list
`idxs::Vector{Int}` of "good" full-space basis indices (the excitation sector
`Σᵢ nᵢ = N` for `:exact_N`, or `Σᵢ nᵢ ≤ N` for `:upto_N`).

`T` is a 0/1 selection isometry: `T' * T = I` on the subspace, and `T * T'` is the
orthogonal projector onto the sector in the full space.
"""
function subspace_transform(enc::DualRailEncoding)
    idxs = _subspace_indices(enc)
    d_full = prod(enc.subspace_levels)
    T = sparse(idxs, 1:length(idxs), ones(ComplexF64, length(idxs)), d_full, length(idxs))
    return T, idxs
end

@doc raw"""
    reduce_to_subspace(O::AbstractMatrix, enc::DualRailEncoding)
    reduce_to_subspace(ψ::AbstractVector, enc::DualRailEncoding)

Reduce a full-space operator `O` (or state `ψ`) to the encoded subspace of `enc`:
`T' * O * T` (or `T' * ψ`) with `T` from [`subspace_transform`](@ref).

Because `T` is a 0/1 selection isometry this equals `O[idxs, idxs]` (or `ψ[idxs]`),
which is how it is computed — the indexing is sparse-aware, so a
`SparseMatrixCSC` input yields a `SparseMatrixCSC` output and the full-space product
is never materialized.
"""
function reduce_to_subspace(O::AbstractMatrix{<:Number}, enc::DualRailEncoding)
    d_full = prod(enc.subspace_levels)
    size(O) == (d_full, d_full) || throw(
        DimensionMismatch(
            "operator size $(size(O)) does not match the full encoding " *
            "dimension ($d_full, $d_full)",
        ),
    )
    idxs = _subspace_indices(enc)
    return O[idxs, idxs]
end

function reduce_to_subspace(ψ::AbstractVector{<:Number}, enc::DualRailEncoding)
    d_full = prod(enc.subspace_levels)
    length(ψ) == d_full || throw(
        DimensionMismatch(
            "state length $(length(ψ)) does not match the full encoding " *
            "dimension $d_full",
        ),
    )
    idxs = _subspace_indices(enc)
    return ψ[idxs]
end

# ----------------------------------------------------------------------------- #
#                          Logical ↔ physical state maps                         #
# ----------------------------------------------------------------------------- #

# Full-space index of the encoded logical basis state with the given bits
# (`bits[q]` is the state of logical qubit `q`, qubit 1 most significant).
function _logical_state_index(bits::AbstractVector{Int}, enc::DualRailEncoding)
    m = enc.N ÷ enc.n_qubits
    occupations = Vector{Int}(undef, enc.n_rails)
    for (q, b) ∈ enumerate(bits)
        occupations[2q-1] = b == 0 ? m : 0
        occupations[2q] = b == 0 ? 0 : m
    end
    return _occupations_to_index(occupations, enc.subspace_levels)
end

@doc raw"""
    logical_state_indices(enc::DualRailEncoding) -> Vector{Int}

Full-space basis indices of the `2^n_qubits` encoded logical basis states, ordered
by logical basis state `|0…0⟩, …, |1…1⟩` (qubit 1 most significant). Note that this
logical ordering is *not* ascending in the full space. Used by the
encoding-aware `EmbeddedOperator` constructors; see also
[`logical_basis_states`](@ref).
"""
function logical_state_indices(enc::DualRailEncoding)
    n = enc.n_qubits
    return [_logical_state_index([(ℓ>>(n-q))&1 for q ∈ 1:n], enc) for ℓ ∈ 0:(2^n-1)]
end

@doc raw"""
    logical_basis_states(enc::DualRailEncoding) -> Vector{Vector{ComplexF64}}

The `2^n_qubits` physical (full-space) kets corresponding to the logical basis
states `|0…0⟩, …, |1…1⟩` (qubit 1 most significant). Use
[`reduce_to_subspace`](@ref) to express them in subspace coordinates.
"""
function logical_basis_states(enc::DualRailEncoding)
    d_full = prod(enc.subspace_levels)
    states = Vector{Vector{ComplexF64}}()
    for index ∈ logical_state_indices(enc)
        ψ = zeros(ComplexF64, d_full)
        ψ[index] = 1.0
        push!(states, ψ)
    end
    return states
end

@doc raw"""
    target_states(gate::Symbol, enc::DualRailEncoding) -> Vector{Vector{ComplexF64}}
    target_states(U::AbstractMatrix{<:Number}, enc::DualRailEncoding) -> Vector{Vector{ComplexF64}}

For each logical basis input `|0…0⟩, …, |1…1⟩`, the encoded physical output ket
under the logical unitary `U` (or `GATES[gate]` for a `Symbol`, e.g. `:I, :X, :Y,
:Z, :H` for one qubit, `:CX, :CZ` for two, `:CCX, :CCZ` for three). The unitary
must act on all `n_qubits` logical qubits; to act on a subset of qubits, see the
`EmbeddedOperator(gate, enc; qubit_indices = ...)` constructor.
"""
function target_states(gate::Symbol, enc::DualRailEncoding)
    gate ∈ keys(GATES) || throw(
        ArgumentError(
            "gate :$gate is not a valid gate. " *
            "See the Piccolo.GATES dict for available gates.",
        ),
    )
    return target_states(GATES[gate], enc)
end

function target_states(U::AbstractMatrix{<:Number}, enc::DualRailEncoding)
    d_logical = 2^enc.n_qubits
    size(U) == (d_logical, d_logical) || throw(
        ArgumentError(
            "the logical unitary has size $(size(U)), but the encoding has " *
            "$(enc.n_qubits) logical qubits (dimension $d_logical); to embed a " *
            "smaller gate on specific qubits use " *
            "EmbeddedOperator(gate, enc; qubit_indices = ...)",
        ),
    )
    ψs = logical_basis_states(enc)
    return [sum(U[j, l] * ψs[j] for j ∈ 1:d_logical) for l ∈ 1:d_logical]
end

end
