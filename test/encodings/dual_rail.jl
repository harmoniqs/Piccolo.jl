@testitem "DualRailEncoding: construction and validation" begin
    enc =
        DualRailEncoding(n_qubits = 3, levels_per_rail = 2, conservation = :exact_N, N = 3)
    @test enc.n_qubits == 3
    @test enc.levels_per_rail == 2
    @test enc.n_rails == 6
    @test enc.subspace_levels == fill(2, 6)
    @test enc.conservation == :exact_N
    @test enc.N == 3

    # Defaults: levels_per_rail = 2, conservation = :exact_N, N = n_qubits
    enc_default = DualRailEncoding(n_qubits = 2)
    @test enc_default.levels_per_rail == 2
    @test enc_default.conservation == :exact_N
    @test enc_default.N == 2

    # Invalid conservation symbol
    @test_throws AssertionError DualRailEncoding(n_qubits = 2, conservation = :bogus)
    # N not divisible by n_qubits
    @test_throws AssertionError DualRailEncoding(n_qubits = 2, N = 3)
    # per-qubit excitation exceeds levels_per_rail - 1
    @test_throws AssertionError DualRailEncoding(n_qubits = 1, levels_per_rail = 2, N = 2)
end

@testitem "DualRailEncoding: subspace_transform round-trip T'T = I" begin
    using LinearAlgebra: I
    using SparseArrays

    for conservation in (:exact_N, :upto_N)
        enc = DualRailEncoding(
            n_qubits = 3,
            levels_per_rail = 2,
            conservation = conservation,
            N = 3,
        )
        T, idxs = subspace_transform(enc)
        @test T isa SparseMatrixCSC{ComplexF64}
        @test size(T) == (prod(enc.subspace_levels), length(idxs))
        @test Matrix(T' * T) ≈ I(length(idxs))
    end
end

@testitem "DualRailEncoding: reduce_to_subspace matches naive T'OT" begin
    using LinearAlgebra
    using SparseArrays

    enc =
        DualRailEncoding(n_qubits = 2, levels_per_rail = 2, conservation = :exact_N, N = 2)
    T, idxs = subspace_transform(enc)
    full_dim = prod(enc.subspace_levels)

    O = rand(ComplexF64, full_dim, full_dim)
    @test Matrix(reduce_to_subspace(O, enc)) ≈ Matrix(T' * O * T)

    ψ = rand(ComplexF64, full_dim)
    @test Vector(reduce_to_subspace(ψ, enc)) ≈ Vector(T' * ψ)

    # Reduced operator picks out exactly the good-index block
    @test Matrix(reduce_to_subspace(O, enc)) ≈ O[idxs, idxs]
end

@testitem "DualRailEncoding: target_states(:CX) for n=2" begin
    using LinearAlgebra: norm

    enc =
        DualRailEncoding(n_qubits = 2, levels_per_rail = 2, conservation = :exact_N, N = 2)
    ψ = logical_basis_states(enc)
    @test length(ψ) == 4
    @test all(norm(s) ≈ 1 for s in ψ)

    tg = target_states(:CX, enc)
    @test length(tg) == 4
    # CX: |00>->|00>, |01>->|01>, |10>->|11>, |11>->|10>
    @test tg[1] == ψ[1]
    @test tg[2] == ψ[2]
    @test tg[3] == ψ[4]
    @test tg[4] == ψ[3]
end

@testitem "DualRailEncoding: target_states(:CCX) for n=3" begin
    enc =
        DualRailEncoding(n_qubits = 3, levels_per_rail = 2, conservation = :exact_N, N = 3)
    ψ = logical_basis_states(enc)
    @test length(ψ) == 8

    tg = target_states(:CCX, enc)
    @test length(tg) == 8
    # CCX swaps |110> (idx 7) and |111> (idx 8), identity elsewhere
    for j = 1:6
        @test tg[j] == ψ[j]
    end
    @test tg[7] == ψ[8]
    @test tg[8] == ψ[7]
end

@testitem "DualRailEncoding: EmbeddedOperator(:CX) is unitary on subspace" begin
    using LinearAlgebra: I

    enc =
        DualRailEncoding(n_qubits = 2, levels_per_rail = 2, conservation = :exact_N, N = 2)
    op = EmbeddedOperator(:CX, enc)
    @test op.subsystem_levels == enc.subspace_levels
    @test length(op.subspace) == 4

    U = unembed(op)
    @test U ≈ GATES[:CX]
    @test U'U ≈ I(4)

    # Embedding a 2-qubit gate on a specific pair of a 3-qubit register
    enc3 =
        DualRailEncoding(n_qubits = 3, levels_per_rail = 2, conservation = :exact_N, N = 3)
    op_pair = EmbeddedOperator(GATES[:CX], enc3; qubit_indices = [1, 2])
    Up = unembed(op_pair)
    @test Up ≈ kron(GATES[:CX], GATES[:I])
    @test Up'Up ≈ I(8)

    op_pair23 = EmbeddedOperator(GATES[:CX], enc3; qubit_indices = [2, 3])
    @test unembed(op_pair23) ≈ kron(GATES[:I], GATES[:CX])
end

@testitem "DualRailEncoding: sparsity nnz(T) matches combinatorial count" begin
    using SparseArrays

    # n_qubits=3, levels_per_rail=2, exact_N, N=3:
    # weight-3 occupation tuples over 6 binary rails = C(6,3) = 20
    enc =
        DualRailEncoding(n_qubits = 3, levels_per_rail = 2, conservation = :exact_N, N = 3)
    T, idxs = subspace_transform(enc)
    @test nnz(T) == 20
    @test length(idxs) == 20
end

@testitem "DualRailEncoding: upto_N strictly contains exact_N" begin
    enc_exact =
        DualRailEncoding(n_qubits = 2, levels_per_rail = 2, conservation = :exact_N, N = 2)
    enc_upto =
        DualRailEncoding(n_qubits = 2, levels_per_rail = 2, conservation = :upto_N, N = 2)

    _, idx_exact = subspace_transform(enc_exact)
    _, idx_upto = subspace_transform(enc_upto)

    @test issubset(Set(idx_exact), Set(idx_upto))
    @test length(idx_upto) > length(idx_exact)
end
