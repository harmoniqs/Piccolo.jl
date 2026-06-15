module Encodings

export DualRailEncoding

export subspace_transform
export reduce_to_subspace
export logical_basis_states
export logical_state_indices
export target_states

using ..Gates

using SparseArrays
using TestItems

# ----------------------------------------------------------------------------- #
#                              DualRailEncoding                                 #
# ----------------------------------------------------------------------------- #

@doc raw"""
    DualRailEncoding

Encoding of `n_qubits` logical qubits into `2 * n_qubits` physical modes (rails),
each truncated to `levels_per_rail` levels. Each logical qubit occupies a pair of
rails, with the convention

```math
|0\rangle_L \leftrightarrow |n, 0\rangle, \qquad |1\rangle_L \leftrightarrow |0, n\rangle
```

where ``n`` is the per-qubit excitation number `N ÷ n_qubits` (typically `1`, the
single-excitation dual-rail code). The pattern works for any excitation-conserving
system (bosonic, spin, fermionic), not just photonics.

# Fields
- `n_qubits::Int`: number of logical qubits.
- `levels_per_rail::Int`: local Hilbert dimension per rail.
- `n_rails::Int`: `2 * n_qubits`.
- `subspace_levels::Vector{Int}`: full tensor-product dims, `fill(levels_per_rail, n_rails)`.
- `conservation::Symbol`: `:exact_N` (closed, ``\sum_i n_i = N``) or `:upto_N`
  (open / lossy, ``\sum_i n_i \le N``).
- `N::Int`: target excitation number (typically `n_qubits`).

# Keyword constructor
    DualRailEncoding(; n_qubits, levels_per_rail=2, conservation=:exact_N, N=n_qubits)

# Notes
The encoding intentionally does *not* carry a Hamiltonian — that is the user's
responsibility. Build a `QuantumSystem` from your own `H_drift` / `H_drives`, then
pass an encoding-aware goal into a trajectory:

```julia
using Piccolo

enc = DualRailEncoding(n_qubits = 2, levels_per_rail = 2, conservation = :exact_N, N = 2)

# user-supplied Hamiltonian on the full rail space
system = QuantumSystem(H_drift, H_drives, drive_bounds)

# encoding-aware goal for a logical CX
U_goal = EmbeddedOperator(:CX, enc)
traj = UnitaryTrajectory(system, pulse, U_goal)
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
    levels_per_rail::Int = 2,
    conservation::Symbol = :exact_N,
    N::Int = n_qubits,
)
    n_qubits ≥ 1 ||
        throw(ArgumentError("n_qubits must be ≥ 1, got $n_qubits"))
    levels_per_rail ≥ 2 ||
        throw(ArgumentError("levels_per_rail must be ≥ 2, got $levels_per_rail"))
    conservation ∈ (:exact_N, :upto_N) ||
        throw(ArgumentError("conservation must be :exact_N or :upto_N, got $conservation"))
    N % n_qubits == 0 ||
        throw(ArgumentError("N ($N) must be divisible by n_qubits ($n_qubits)"))
    n_per_qubit = N ÷ n_qubits
    n_per_qubit ≤ levels_per_rail - 1 || throw(
        ArgumentError(
            "per-qubit excitation ($n_per_qubit) exceeds levels_per_rail - 1 ($(levels_per_rail - 1))",
        ),
    )
    n_rails = 2 * n_qubits
    subspace_levels = fill(levels_per_rail, n_rails)
    return DualRailEncoding(
        n_qubits,
        levels_per_rail,
        n_rails,
        subspace_levels,
        conservation,
        N,
    )
end

"""
    excitation_per_qubit(enc::DualRailEncoding)

Per-qubit excitation number, `N ÷ n_qubits`.
"""
excitation_per_qubit(enc::DualRailEncoding) = enc.N ÷ enc.n_qubits

# ----------------------------------------------------------------------------- #
#                            Index <-> occupation                              #
# ----------------------------------------------------------------------------- #
# Convention matches `EmbeddedOperators.basis_labels(...; baseline=0)`: the first
# rail is the most significant digit (changes slowest), indices are 1-based.

function index_to_occupation(idx::Int, levels::AbstractVector{Int})
    r = idx - 1
    n = length(levels)
    occ = zeros(Int, n)
    for k = n:-1:1
        occ[k] = r % levels[k]
        r ÷= levels[k]
    end
    return occ
end

function occupation_to_index(occ::AbstractVector{Int}, levels::AbstractVector{Int})
    idx = 0
    for k in eachindex(levels)
        idx = idx * levels[k] + occ[k]
    end
    return idx + 1
end

# ----------------------------------------------------------------------------- #
#                            Subspace transform                                #
# ----------------------------------------------------------------------------- #

"""
    good_indices(enc::DualRailEncoding) -> Vector{Int}

The "good" full-space basis indices of the excitation-number subspace: occupation
tuples with ``\\sum_i n_i = N`` for `:exact_N`, or ``\\le N`` for `:upto_N`.
"""
function good_indices(enc::DualRailEncoding)
    levels = enc.subspace_levels
    full_dim = prod(levels)
    keep = enc.conservation == :exact_N ? (s -> s == enc.N) : (s -> s ≤ enc.N)
    idxs = Int[]
    for idx = 1:full_dim
        occ = index_to_occupation(idx, levels)
        if keep(sum(occ))
            push!(idxs, idx)
        end
    end
    return idxs
end

@doc raw"""
    subspace_transform(enc::DualRailEncoding) -> (T::SparseMatrixCSC{ComplexF64}, idxs::Vector{Int})

Projector mapping subspace coordinates to full coordinates, ``|\psi_\text{full}\rangle = T |\psi_\text{sub}\rangle``,
together with the list of "good" basis indices `idxs`. `T` is `full_dim × length(idxs)`
with a single unit entry per column. Supports both `:exact_N` and `:upto_N`.
"""
function subspace_transform(enc::DualRailEncoding)
    idxs = good_indices(enc)
    full_dim = prod(enc.subspace_levels)
    k = length(idxs)
    T = sparse(idxs, collect(1:k), ones(ComplexF64, k), full_dim, k)
    return T, idxs
end

"""
    reduce_to_subspace(O::AbstractMatrix, enc::DualRailEncoding)
    reduce_to_subspace(ψ::AbstractVector, enc::DualRailEncoding)

Reduce a full-space operator (`T' * O * T`) or state (`T' * ψ`) to the encoded
subspace. Sparse-aware.
"""
function reduce_to_subspace(O::AbstractMatrix, enc::DualRailEncoding)
    T, _ = subspace_transform(enc)
    return T' * O * T
end

function reduce_to_subspace(ψ::AbstractVector, enc::DualRailEncoding)
    T, _ = subspace_transform(enc)
    return T' * ψ
end

# ----------------------------------------------------------------------------- #
#                            Logical states                                    #
# ----------------------------------------------------------------------------- #

"""
    logical_state_indices(enc::DualRailEncoding) -> Vector{Int}

Full-space indices of the `2^n_qubits` logical basis states, ordered
`|0…0⟩, …, |1…1⟩` (first qubit most significant). For qubit `i` (rails `2i-1`, `2i`),
bit `0` places `n` excitations in the first rail (`|n,0⟩`) and bit `1` in the second
(`|0,n⟩`), where `n = N ÷ n_qubits`.
"""
function logical_state_indices(enc::DualRailEncoding)
    n = enc.n_qubits
    n_per = excitation_per_qubit(enc)
    levels = enc.subspace_levels
    indices = Vector{Int}(undef, 2^n)
    for j = 0:(2^n-1)
        occ = zeros(Int, enc.n_rails)
        for i = 1:n
            bit = (j >> (n - i)) & 1
            if bit == 0
                occ[2i-1] = n_per
            else
                occ[2i] = n_per
            end
        end
        indices[j+1] = occupation_to_index(occ, levels)
    end
    return indices
end

"""
    logical_basis_states(enc::DualRailEncoding) -> Vector{Vector{ComplexF64}}

The `2^n_qubits` physical kets (full-space dimension) corresponding to the logical
basis states `|0…0⟩, …, |1…1⟩`.
"""
function logical_basis_states(enc::DualRailEncoding)
    full_dim = prod(enc.subspace_levels)
    idxs = logical_state_indices(enc)
    states = Vector{Vector{ComplexF64}}(undef, length(idxs))
    for (j, idx) in enumerate(idxs)
        ψ = zeros(ComplexF64, full_dim)
        ψ[idx] = 1
        states[j] = ψ
    end
    return states
end

"""
    target_states(gate::Symbol, enc::DualRailEncoding) -> Vector{Vector{ComplexF64}}

For each logical input basis state, the encoded physical output ket under the logical
`gate`. The logical unitary is looked up in `Piccolo.GATES` and must act on the full
register (`size == 2^n_qubits`). Supports `:I, :X, :Y, :Z, :H, :CX, :CZ, :CCX, :CCZ`
(any gate present in `GATES` of the right dimension).
"""
function target_states(gate::Symbol, enc::DualRailEncoding)
    if gate ∉ keys(GATES)
        throw(
            ArgumentError(
                "Gate must be a valid gate. See Piccolo.GATES for available gates.",
            ),
        )
    end
    U = GATES[gate]
    n = enc.n_qubits
    size(U, 1) == 2^n || throw(
        ArgumentError(
            "Gate dimension $(size(U, 1)) does not match 2^n_qubits = $(2^n).",
        ),
    )
    ψ_logical = logical_basis_states(enc)
    full_dim = length(ψ_logical[1])
    targets = Vector{Vector{ComplexF64}}(undef, 2^n)
    for j = 1:(2^n)
        out = zeros(ComplexF64, full_dim)
        for i = 1:(2^n)
            out .+= U[i, j] .* ψ_logical[i]
        end
        targets[j] = out
    end
    return targets
end

end
