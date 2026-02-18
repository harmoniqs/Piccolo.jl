module QuantumObjectUtils

export operator_from_string
export ket_from_string
export ket_from_bitstring

export haar_random
export haar_identity

export create
export annihilate

using ..Gates

using LinearAlgebra
using TestItems

@doc raw"""
    operator_from_string(operator::String; lookup=PAULIS)

Reduce the string (each character is one key) via operators from a dictionary.
"""
function operator_from_string(
    operator::String;
    lookup::NamedTuple{names,<:Tuple{Vararg{AbstractMatrix{<:Number}}}} where {names} = PAULIS,
)::Matrix{ComplexF64}
    # split string into keys and replace with operators
    characters = [lookup[Symbol(c)] for c ∈ operator]
    return foldr(kron, characters)
end

@doc raw"""
    ket_from_string(
        ket::String,
        levels::Vector{Int};
        level_dict=Dict(:g => 0, :e => 1, :f => 2, :h => 3, :i => 4, :j => 5, :k => 6, :l => 7),
        return_states=false
    )

Construct a quantum state from a string ket representation.
"""
function ket_from_string(
    ket::String,
    levels::Vector{Int};
    level_dict = Dict(
        :g => 0,
        :e => 1,
        :f => 2,
        :h => 3,
        :i => 4,
        :j => 5,
        :k => 6,
        :l => 7,
    ),
    return_states = false,
)::Vector{ComplexF64}
    kets = []

    for x ∈ split(ket, ['(', ')'])
        if x == ""
            continue
        elseif all(Symbol(xᵢ) ∈ keys(level_dict) for xᵢ ∈ x)
            append!(kets, x)
        elseif occursin("+", x)
            superposition = split(x, '+')
            @assert all(all(Symbol(xᵢ) ∈ keys(level_dict) for xᵢ ∈ x) for x ∈ superposition) "Invalid ket: $x"
            @assert length(superposition) == 2 "Only two states can be superposed for now"
            push!(kets, x)
        else
            error("Invalid ket: $x")
        end
    end

    states = []

    for (ψᵢ, l) ∈ zip(kets, levels)
        if ψᵢ isa AbstractString && occursin("+", ψᵢ)
            superposition = split(ψᵢ, '+')
            superposition_states = [level_dict[Symbol(x)] for x ∈ superposition]
            @assert all(state ≤ l - 1 for state ∈ superposition_states) "Level $ψᵢ is not allowed for $l levels"
            superposition_state =
                sum([ComplexF64.(I[1:l, state+1]) for state ∈ superposition_states])
            normalize!(superposition_state)
            push!(states, superposition_state)
        else
            state = level_dict[Symbol(ψᵢ)]
            @assert state ≤ l - 1 "Level $ψᵢ is not allowed for $l levels"
            push!(states, ComplexF64.(I[1:l, state+1]))
        end
    end

    if return_states
        return states
    else
        return kron([1.0], states...)
    end
end

@doc raw"""
    ket_from_bitstring(ket::String)

Get the state vector for a qubit system given a ket string `ket` of 0s and 1s.
"""
function ket_from_bitstring(ket::String)::Vector{ComplexF64}
    cs = [c for c ∈ ket]
    @assert all(c ∈ "01" for c ∈ cs)
    states = [c == '0' ? [1, 0] : [0, 1] for c ∈ cs]
    return foldr(kron, states)
end

# --------------------------------------------------------------------------- #
# Random operators
# --------------------------------------------------------------------------- #

@doc raw"""
    haar_random(n::Int)

Generate a random unitary matrix using the Haar measure for an `n`-dimensional system.
"""
function haar_random(n::Int)
    # Ginibre matrix
    Z = (randn(n, n) + im * randn(n, n)) / √2
    F = qr(Z)
    # QR correction (R main diagonal is real, strictly positive)
    Λ = diagm(diag(F.R) ./ abs.(diag(F.R)))
    return F.Q * Λ
end

@doc raw"""
    haar_identity(n::Int, radius::Number)

Generate a random unitary matrix close to the identity matrix using the Haar measure for
an `n`-dimensional system with a given `radius`. The smaller the radius, the closer the
matrix will be to the identity.
"""
function haar_identity(n::Int, radius::Number)
    # Ginibre matrix
    Z = (I + radius * (randn(n, n) + im * randn(n, n)) / √2) / (1 + radius)
    F = qr(Z)
    # QR correction (R main diagonal is real, strictly positive)
    Λ = diagm(diag(F.R) ./ abs.(diag(F.R)))
    return F.Q * Λ
end

# --------------------------------------------------------------------------- #
# Oscillator operators
# --------------------------------------------------------------------------- #

@doc raw"""
    annihilate(levels::Int)

Get the annihilation operator for a system with `levels`.
"""
annihilate(levels::Int)::Matrix{ComplexF64} = diagm(1 => map(sqrt, 1:(levels-1)))

@doc raw"""
    create(levels::Int)

Get the creation operator for a system with `levels`.
"""
create(levels::Int) = collect(annihilate(levels)')

# ****************************************************************************** #

@testitem "Test ket_from_bitstring function" begin
    @test ket_from_bitstring("0") == [1, 0]
    @test ket_from_bitstring("1") == [0, 1]
    @test ket_from_bitstring("00") == [1, 0, 0, 0]
    @test ket_from_bitstring("01") == [0, 1, 0, 0]
    @test ket_from_bitstring("10") == [0, 0, 1, 0]
    @test ket_from_bitstring("11") == [0, 0, 0, 1]
end

@testitem "Test annihilate" begin
    @test annihilate(2) ≈ [0 1; 0 0]
    @test annihilate(3) ≈ [0 1 0; 0 0 √2; 0 0 0]
end

@testitem "Test create" begin
    @test create(2) ≈ [0 0; 1 0]
    @test create(3) ≈ [0 0 0; 1 0 0; 0 √2 0]
end

@testitem "Test operator_from_string" begin
    @test operator_from_string("X") ≈ [0 1; 1 0]
    @test operator_from_string("XZ") ≈ [0 0 1 0; 0 0 0 -1; 1 0 0 0; 0 -1 0 0]
end

@testitem "Test ket_from_string" begin
    @test ket_from_string("g", [2]) ≈ [1, 0]
    @test ket_from_string("gg", [2, 2]) ≈ [1, 0, 0, 0]
    @test ket_from_string("(g+e)g", [2, 2]) ≈ [1 / √2, 0, 1 / √2, 0]
end

@testitem "Test Haar random" begin
    using LinearAlgebra: I
    U = haar_random(2)
    @test size(U) == (2, 2)
    @test U'U ≈ I atol = 1e-10

    U = haar_random(3)
    @test size(U) == (3, 3)
    @test U'U ≈ I atol = 1e-10

    U = haar_identity(2, 0.1)
    @test size(U) == (2, 2)
    @test U'U ≈ I atol = 1e-10

    Id = haar_identity(3, 0.0)
    @test size(Id) == (3, 3)
    @test Id ≈ I atol = 1e-10
end

end
