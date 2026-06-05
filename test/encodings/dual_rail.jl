using TestItems

@testitem "DualRailEncoding subspace transforms" begin
    using LinearAlgebra
    using Piccolo
    using SparseArrays

    enc = DualRailEncoding(n_qubits = 2, levels_per_rail = 2, conservation = :exact_N, N = 2)
    T, idxs = subspace_transform(enc)

    @test issparse(T)
    @test size(T) == (16, 6)
    @test nnz(T) == 6
    @test T' * T == I(6)
    @test idxs == [4, 6, 7, 10, 11, 13]

    upto = DualRailEncoding(n_qubits = 2, levels_per_rail = 2, conservation = :upto_N, N = 2)
    _, upto_idxs = subspace_transform(upto)
    @test issubset(idxs, upto_idxs)
    @test length(upto_idxs) > length(idxs)
end

@testitem "reduce_to_subspace matches explicit projection" begin
    using LinearAlgebra
    using Piccolo
    using SparseArrays

    enc = DualRailEncoding(n_qubits = 1, levels_per_rail = 2, conservation = :exact_N, N = 1)
    T, _ = subspace_transform(enc)
    full_op = sparse(Diagonal(ComplexF64[1, 2, 3, 4]))
    full_state = ComplexF64[1, 2, 3, 4]

    @test reduce_to_subspace(full_op, enc) == T' * full_op * T
    @test reduce_to_subspace(full_state, enc) == T' * full_state
end

@testitem "target_states returns encoded CX and CCX outputs" begin
    using Piccolo

    enc2 = DualRailEncoding(n_qubits = 2, levels_per_rail = 2, conservation = :exact_N, N = 2)
    logical2 = logical_basis_states(enc2)
    cx_targets = target_states(:CX, enc2)

    @test cx_targets[1] == logical2[1] # 00 -> 00
    @test cx_targets[2] == logical2[2] # 01 -> 01
    @test cx_targets[3] == logical2[4] # 10 -> 11
    @test cx_targets[4] == logical2[3] # 11 -> 10

    enc3 = DualRailEncoding(n_qubits = 3, levels_per_rail = 2, conservation = :exact_N, N = 3)
    logical3 = logical_basis_states(enc3)
    ccx_targets = target_states(:CCX, enc3)

    @test ccx_targets[7] == logical3[8] # 110 -> 111
    @test ccx_targets[8] == logical3[7] # 111 -> 110
end

@testitem "EmbeddedOperator accepts dual rail encodings" begin
    using LinearAlgebra
    using Piccolo
    using SparseArrays

    enc = DualRailEncoding(n_qubits = 2, levels_per_rail = 2, conservation = :exact_N, N = 2)
    op = EmbeddedOperator(:CX, enc)

    @test issparse(op.operator)
    @test op.subsystem_levels == [2, 2, 2, 2]

    sub = unembed(op)
    @test sub' * sub ≈ I(size(sub, 1))
end
