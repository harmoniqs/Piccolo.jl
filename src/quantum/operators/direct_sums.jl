module DirectSums

export direct_sum

using SparseArrays
using SparseArrays: blockdiag
using TestItems

using ..QuantumSystems
using ..Isomorphisms: operator_to_iso_vec, iso_vec_to_operator


"""
    direct_sum(A::AbstractMatrix, B::AbstractMatrix)

Returns the direct sum of two matrices.
"""
function direct_sum(A::AbstractMatrix, B::AbstractMatrix)
    return [A spzeros((size(A, 1), size(B, 2))); spzeros((size(B, 1), size(A, 2))) B]
end

"""
    direct_sum(A::SparseMatrixCSC, B::SparseMatrixCSC)

Returns the direct sum of two sparse matrices.
"""
function direct_sum(A::SparseMatrixCSC, B::SparseMatrixCSC)
    return blockdiag(A, B)
end

"""
    direct_sum(Ã⃗::AbstractVector, B̃⃗::AbstractVector)

Returns the direct sum of two iso_vec operators.
"""
function direct_sum(Ã⃗::AbstractVector, B̃⃗::AbstractVector)
    return operator_to_iso_vec(
        direct_sum(iso_vec_to_operator(Ã⃗), iso_vec_to_operator(B̃⃗)),
    )
end

"""
    direct_sum(sys1::QuantumSystem, sys2::QuantumSystem)

Returns the direct sum of two `QuantumSystem` objects.

Constructs a new system where the Hilbert space is the direct sum of the two input systems:
H = H₁ ⊕ H₂ = [H₁  0 ]
               [0   H₂]

Both systems must have the same number of drives. The resulting system uses sys1's drive_bounds.

# Example
```julia
sys1 = QuantumSystem([PAULIS[:X]], [(-1.0, 1.0)])
sys2 = QuantumSystem([PAULIS[:Y]], [(-1.0, 1.0)])
sys_combined = direct_sum(sys1, sys2)
```
"""
function direct_sum(sys1::QuantumSystem, sys2::QuantumSystem)
    @assert sys1.n_drives == sys2.n_drives "System 1 drives ($(sys1.n_drives)) must equal System 2 drives ($(sys2.n_drives))"
    n_drives = sys1.n_drives
    H = (u, t) -> direct_sum(sys1.H(u, t), sys2.H(u, t))

    # Get combined drive bounds from both systems - but they should be the same length since n_drives must match
    drive_bounds = sys1.drive_bounds  # They should be the same as sys2.drive_bounds

    # Create new QuantumSystem with the Hamiltonian function
    return QuantumSystem(H, drive_bounds)
end

direct_sum(systems::AbstractVector{<:QuantumSystem}) = reduce(direct_sum, systems)


# *************************************************************************** #

@testitem "Test matrix direct sum" begin
    using SparseArrays
    using LinearAlgebra

    # Test with Hermitian matrices (valid Hamiltonians)
    A = ComplexF64[1 1+im; 1-im 2]
    B = ComplexF64[3 2-im; 2+im 4]
    result = direct_sum(A, B)
    @test result == ComplexF64[1 1+im 0 0; 1-im 2 0 0; 0 0 3 2-im; 0 0 2+im 4]
    @test ishermitian(result)

    # Test with sparse Hermitian matrices
    A = sparse(ComplexF64[1 1+im; 1-im 2])
    B = sparse(ComplexF64[3 2-im; 2+im 4])
    result = direct_sum(A, B)
    @test result == sparse(ComplexF64[1 1+im 0 0; 1-im 2 0 0; 0 0 3 2-im; 0 0 2+im 4])
    @test ishermitian(result)
end



@testitem "Test quantum system direct sum" begin
    using LinearAlgebra

    # Test with no drives - use Hermitian matrices
    H1 = ComplexF64[1 1+im; 1-im 2]
    H2 = ComplexF64[3 2-im; 2+im 4]
    sys1 = QuantumSystem(H1)
    sys2 = QuantumSystem(H2)
    sys = direct_sum(sys1, sys2)
    result = sys.H(Float64[], 0.0)
    @test result == ComplexF64[1 1+im 0 0; 1-im 2 0 0; 0 0 3 2-im; 0 0 2+im 4]
    @test ishermitian(result)
    @test sys.n_drives == 0

    # Test with drives
    sys1 = QuantumSystem([PAULIS[:X]], [(-1.0, 1.0)])
    sys2 = QuantumSystem([PAULIS[:Y]], [(-1.0, 1.0)])
    sys = direct_sum(sys1, sys2)
    @test sys.n_drives == 1
    @test sys.drive_bounds == [(-1.0, 1.0)]

    # Check Hamiltonian structure
    u = [0.5]
    H = sys.H(u, 0.0)
    @test size(H) == (4, 4)
    @test ishermitian(H)

    # Test with multiple drives
    sys1 = QuantumSystem([PAULIS[:X], PAULIS[:Y]], [1.0, 1.0])
    sys2 = QuantumSystem([PAULIS[:Z], PAULIS[:X]], [1.0, 1.0])
    sys = direct_sum(sys1, sys2)
    @test sys.n_drives == 2

    # Test with vector of systems
    sys1 = QuantumSystem(PAULIS[:Z])
    sys2 = QuantumSystem(PAULIS[:X])
    sys3 = QuantumSystem(PAULIS[:Y])
    sys = direct_sum([sys1, sys2, sys3])
    @test sys.n_drives == 0
    H = sys.H(Float64[], 0.0)
    @test size(H) == (6, 6)
    @test ishermitian(H)
end

@testitem "Test direct sum error handling" begin

    # Test mismatched n_drives
    sys1 = QuantumSystem([PAULIS[:X]], [(-1.0, 1.0)])
    sys2 = QuantumSystem([PAULIS[:Y], PAULIS[:Z]], [(-1.0, 1.0), (-1.0, 1.0)])

    @test_throws AssertionError direct_sum(sys1, sys2)
end

end
