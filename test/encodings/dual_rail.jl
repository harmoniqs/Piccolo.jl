@testitem "DualRailEncoding construction and validation" begin
    enc =
        DualRailEncoding(n_qubits = 3, levels_per_rail = 2, conservation = :exact_N, N = 3)
    @test enc.n_rails == 6
    @test enc.subspace_levels == fill(2, 6)
    @test enc.conservation == :exact_N
    @test enc.N == 3

    # defaults: qubit rails, closed system, one excitation per logical qubit
    enc_d = DualRailEncoding(n_qubits = 2)
    @test enc_d.levels_per_rail == 2
    @test enc_d.conservation == :exact_N
    @test enc_d.N == 2

    @test_throws ArgumentError DualRailEncoding(n_qubits = 0)
    @test_throws ArgumentError DualRailEncoding(n_qubits = 2, levels_per_rail = 1)
    @test_throws ArgumentError DualRailEncoding(n_qubits = 2, conservation = :sometimes)
    # N must give an integer number of excitations per logical qubit
    @test_throws ArgumentError DualRailEncoding(n_qubits = 2, N = 3)
    @test_throws ArgumentError DualRailEncoding(n_qubits = 2, N = 0)
    # m = N ÷ n_qubits excitations must fit in a rail: m ≤ levels_per_rail - 1
    @test_throws ArgumentError DualRailEncoding(n_qubits = 1, levels_per_rail = 2, N = 2)
end

@testitem "Dual-rail subspace transform round-trip" begin
    using LinearAlgebra: I
    using SparseArrays: SparseMatrixCSC, nnz

    enc =
        DualRailEncoding(n_qubits = 2, levels_per_rail = 2, conservation = :exact_N, N = 2)
    T, idxs = subspace_transform(enc)
    @test T isa SparseMatrixCSC{ComplexF64}
    @test size(T) == (16, 6)
    @test nnz(T) == 6
    # hand-enumerated Σᵢnᵢ = 2 sector of 4 qubit rails (kron order, rail 1 most
    # significant): occupations 0011, 0101, 0110, 1001, 1010, 1100
    @test idxs == [4, 6, 7, 10, 11, 13]
    @test issorted(idxs)
    # round-trip: T' * T = I on the subspace
    @test T' * T == I

    enc_u = DualRailEncoding(n_qubits = 2, conservation = :upto_N)
    T_u, idxs_u = subspace_transform(enc_u)
    @test T_u' * T_u == I
    # hand-enumerated Σᵢnᵢ ≤ 2 sector: sectors 0, 1, and 2
    @test idxs_u == [1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 13]
    # :upto_N strictly contains the :exact_N index set
    @test issubset(idxs, idxs_u)
    @test length(idxs_u) > length(idxs)
end

@testitem "Dual-rail subspace transform sparsity counts" begin
    using SparseArrays: nnz

    # n_qubits = 3, levels_per_rail = 2: rails are qubits, so the exact-N sector
    # has binomial(6, 3) basis states and T has one entry per kept state
    enc = DualRailEncoding(n_qubits = 3, levels_per_rail = 2, N = 3)
    T, idxs = subspace_transform(enc)
    @test size(T) == (64, binomial(6, 3))
    @test nnz(T) == binomial(6, 3)
    @test length(idxs) == binomial(6, 3)

    enc_u = DualRailEncoding(n_qubits = 3, conservation = :upto_N)
    T_u, idxs_u = subspace_transform(enc_u)
    @test nnz(T_u) == sum(binomial(6, k) for k ∈ 0:3)
    @test issubset(idxs, idxs_u)
    @test length(idxs_u) > length(idxs)
end

@testitem "Dual-rail reduce_to_subspace" begin
    using SparseArrays: spdiagm, issparse

    enc = DualRailEncoding(n_qubits = 2)
    T, idxs = subspace_transform(enc)

    # agrees with the naive dense projection T' * O * T
    O = reshape(collect(1.0:256.0), 16, 16)
    O_sub = reduce_to_subspace(O, enc)
    @test size(O_sub) == (6, 6)
    @test O_sub ≈ Matrix(T)' * O * Matrix(T)

    # sparse-aware: sparse in, sparse out, same values
    O_sparse = spdiagm(0 => collect(1.0:16.0), 1 => fill(0.5, 15))
    O_sparse_sub = reduce_to_subspace(O_sparse, enc)
    @test issparse(O_sparse_sub)
    @test Matrix(O_sparse_sub) ≈ Matrix(T)' * Matrix(O_sparse) * Matrix(T)

    # states reduce consistently, and subspace states re-expand exactly
    ψ = collect(1.0:16.0)
    @test reduce_to_subspace(ψ, enc) == ψ[idxs]
    for ψ_l ∈ logical_basis_states(enc)
        @test T * reduce_to_subspace(ψ_l, enc) == ψ_l
    end

    @test_throws DimensionMismatch reduce_to_subspace(zeros(3, 3), enc)
    @test_throws DimensionMismatch reduce_to_subspace(zeros(3), enc)
end

@testitem "Dual-rail logical basis states" begin
    using LinearAlgebra: I, norm

    # single qubit: |0⟩ ↦ |1,0⟩ (index 3), |1⟩ ↦ |0,1⟩ (index 2)
    enc1 = DualRailEncoding(n_qubits = 1)
    ψs1 = logical_basis_states(enc1)
    @test length(ψs1) == 2
    @test ψs1[1][3] == 1.0 && norm(ψs1[1]) == 1.0
    @test ψs1[2][2] == 1.0 && norm(ψs1[2]) == 1.0

    # two qubits: |00⟩, |01⟩, |10⟩, |11⟩ ↦ occupations 1010, 1001, 0110, 0101
    enc = DualRailEncoding(n_qubits = 2)
    ψs = logical_basis_states(enc)
    @test length(ψs) == 4
    @test [findfirst(x -> x ≠ 0, ψ) for ψ ∈ ψs] == [11, 10, 7, 6]

    # orthonormal, and inside the :exact_N sector (projector acts as identity)
    Ψ = reduce(hcat, ψs)
    @test Ψ' * Ψ == I
    T, _ = subspace_transform(enc)
    for ψ ∈ ψs
        @test T * T' * ψ == ψ
    end

    # generalized rails: m = 2 excitations per qubit on 3-level rails,
    # |0⟩ ↦ |2,0⟩ (index 7), |1⟩ ↦ |0,2⟩ (index 3)
    enc_g = DualRailEncoding(n_qubits = 1, levels_per_rail = 3, N = 2)
    ψs_g = logical_basis_states(enc_g)
    @test [findfirst(x -> x ≠ 0, ψ) for ψ ∈ ψs_g] == [7, 3]
    _, idxs_g = subspace_transform(enc_g)
    @test idxs_g == [3, 5, 7]
end

@testitem "Dual-rail target states CX" begin
    enc =
        DualRailEncoding(n_qubits = 2, levels_per_rail = 2, conservation = :exact_N, N = 2)
    ψs = logical_basis_states(enc)
    ts = target_states(:CX, enc)
    @test length(ts) == 4
    # CX: |00⟩ → |00⟩, |01⟩ → |01⟩, |10⟩ → |11⟩, |11⟩ → |10⟩
    @test ts[1] == ψs[1]
    @test ts[2] == ψs[2]
    @test ts[3] == ψs[4]
    @test ts[4] == ψs[3]
    # explicit physical indices of the outputs (kron order, rail 1 most significant)
    @test [findfirst(x -> x ≠ 0, ψ) for ψ ∈ ts] == [11, 10, 6, 7]
end

@testitem "Dual-rail target states CCX" begin
    enc = DualRailEncoding(n_qubits = 3)
    ψs = logical_basis_states(enc)
    ts = target_states(:CCX, enc)
    @test length(ts) == 8
    # CCX swaps |110⟩ ↔ |111⟩ and is the identity elsewhere
    for l ∈ 1:6
        @test ts[l] == ψs[l]
    end
    @test ts[7] == ψs[8]
    @test ts[8] == ψs[7]
end

@testitem "Dual-rail target states superpositions and errors" begin
    # Hadamard exercises non-permutation (superposition) re-encoding
    enc1 = DualRailEncoding(n_qubits = 1)
    ψs = logical_basis_states(enc1)
    ts = target_states(:H, enc1)
    @test ts[1] ≈ (ψs[1] + ψs[2]) / √2
    @test ts[2] ≈ (ψs[1] - ψs[2]) / √2

    # matrix and Symbol methods agree
    @test target_states(GATES[:H], enc1) == ts

    enc = DualRailEncoding(n_qubits = 2)
    @test_throws ArgumentError target_states(:X, enc)        # 1-qubit gate, 2 qubits
    @test_throws ArgumentError target_states(:NOTAGATE, enc)
    @test_throws ArgumentError target_states(zeros(ComplexF64, 3, 3), enc)
end

@testitem "CCX and CCZ gates" begin
    using LinearAlgebra: I, Diagonal

    @test GATES[:CCX] isa Matrix{ComplexF64}
    @test GATES[:CCZ] isa Matrix{ComplexF64}
    @test size(GATES[:CCX]) == (8, 8)
    @test size(GATES[:CCZ]) == (8, 8)
    # CCX swaps |110⟩ ↔ |111⟩ (indices 7, 8), identity elsewhere
    @test GATES[:CCX][1:6, 1:6] == I
    @test GATES[:CCX][7, 8] == 1 && GATES[:CCX][8, 7] == 1 && GATES[:CCX][7, 7] == 0
    @test GATES[:CCX]' * GATES[:CCX] ≈ I
    # CCZ flips the phase of |111⟩
    @test GATES[:CCZ] == Diagonal([1, 1, 1, 1, 1, 1, 1, -1])
end

@testitem "Dual-rail EmbeddedOperator constructors" begin
    using LinearAlgebra: I

    enc =
        DualRailEncoding(n_qubits = 2, levels_per_rail = 2, conservation = :exact_N, N = 2)
    op = EmbeddedOperator(:CX, enc)
    @test op.subsystem_levels == [2, 2, 2, 2]
    @test size(op.operator) == (16, 16)
    # subspace is ordered by logical basis state, so unembed returns the gate
    @test op.subspace == [11, 10, 7, 6]
    U_sub = unembed(op)
    @test U_sub == GATES[:CX]
    # unitary on the embedded space
    @test U_sub' * U_sub ≈ I

    # Symbol and matrix constructors agree
    @test EmbeddedOperator(GATES[:CX], enc).operator == op.operator

    # consistency: the embedded operator maps encoded inputs to target states
    ψs = logical_basis_states(enc)
    ts = target_states(:CX, enc)
    for (ψ, ψ_target) ∈ zip(ψs, ts)
        @test op.operator * ψ ≈ ψ_target
    end

    # embedding on a specific qubit pair within a 3-qubit register
    enc3 = DualRailEncoding(n_qubits = 3)
    op12 = EmbeddedOperator(:CX, enc3; qubit_indices = [1, 2])
    @test unembed(op12) ≈ kron(GATES[:CX], PAULIS.I)
    op23 = EmbeddedOperator(:CX, enc3; qubit_indices = [2, 3])
    @test unembed(op23) ≈ kron(PAULIS.I, GATES[:CX])
    @test op12.subspace == [findfirst(x -> x ≠ 0, ψ) for ψ ∈ logical_basis_states(enc3)]

    @test_throws ArgumentError EmbeddedOperator(:CX, enc3)  # 2-qubit gate, 3 qubits
    @test_throws ArgumentError EmbeddedOperator(:CX, enc3; qubit_indices = [1, 1])
    @test_throws ArgumentError EmbeddedOperator(:CX, enc3; qubit_indices = [0, 1])
    @test_throws ArgumentError EmbeddedOperator(:NOTAGATE, enc)
end
