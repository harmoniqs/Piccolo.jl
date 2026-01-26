module Gates

export PAULIS
export GATES

using TestItems

@doc raw"""
The 2×2 Pauli matrics and identity.
"""
const PAULIS = (
    I = ComplexF64[1 0;
                   0 1],
    X = ComplexF64[0 1;
                   1 0],
    Y = ComplexF64[0 -im;
                   im 0],
    Z = ComplexF64[1 0;
                   0 -1],
)

@doc raw"""
A constant dictionary `GATES` containing common quantum gate matrices as complex-valued matrices. Each gate is represented by its unitary matrix.

- `GATES[:I]` - Identity: Leaves the state unchanged.
- `GATES[:X]` - Pauli-X (NOT): Flips the qubit state.
- `GATES[:Y]` - Pauli-Y: Rotates the qubit state around the Y-axis of the Bloch sphere.
- `GATES[:Z]` - Pauli-Z: Flips the phase of the qubit state.
- `GATES[:H]` - Hadamard: Creates superposition by transforming basis states.
- `GATES[:CX]` - Controlled-X (CNOT): Flips the 2nd qubit (target) if the first qubit (control) is |1⟩.
- `GATES[:CZ]` - Controlled-Z (CZ): Flips the phase of the 2nd qubit (target) if the 1st qubit (control) is |1⟩.
- `GATES[:XI]` - Complex: A gate for complex operations.
- `GATES[:sqrtiSWAP]` - Square root of iSWAP: Partially swaps two qubits with a phase.
"""
const GATES = (
    I = PAULIS.I,

    X = PAULIS.X,

    Y = PAULIS.Y,

    Z = PAULIS.Z,

    H = ComplexF64[1 1;
                    1 -1]/√2,

    CX = ComplexF64[1 0 0 0;
                    0 1 0 0;
                    0 0 0 1;
                    0 0 1 0],

    CZ = ComplexF64[1 0 0 0;
                    0 1 0 0;
                    0 0 1 0;
                    0 0 0 -1],

    XI = ComplexF64[0 0 -im 0;
                    0 0 0 -im;
                    -im 0 0 0;
                    0 -im 0 0],

    sqrtiSWAP = ComplexF64[1 0 0 0;
                            0 1/√2 1im/√2 0;
                            0 1im/√2 1/√2 0;
                            0 0 0 1],
)

# ******************************************************************************* #

@testitem "Gates are complex matrices" begin
    # test call
    @test GATES.X == [0 1; 1 0]
    @test GATES[:X] == [0 1; 1 0]

    for k in keys(GATES)
        @test typeof(GATES[k]) == Matrix{ComplexF64}
    end
end

@testitem "Paulis are complex matrices" begin
    # test call
    @test PAULIS.X == [0 1; 1 0]
    @test PAULIS[:X] == [0 1; 1 0]
    
    for k in keys(PAULIS)
        @test typeof(PAULIS[k]) == Matrix{ComplexF64}
    end
end

end
