# ----------------------------------------------------------------------------- #
# Composite Quantum Systems
# ----------------------------------------------------------------------------- #

"""
    CompositeQuantumSystem <: AbstractQuantumSystem

A composite quantum system consisting of multiple `subsystems` with optional coupling terms.

Composite systems represent multiple quantum subsystems (e.g., multiple qubits or oscillators)
that may be coupled together. Each subsystem's Hamiltonians are automatically lifted to the 
full tensor product space, and subsystem drives are appended to any coupling drives.

# Fields
- `H::Function`: The total Hamiltonian function: (u, t) -> H(u, t)
- `G::Function`: The isomorphic generator function: (u, t) -> G(u, t)
- `H_drift::SparseMatrixCSC{ComplexF64, Int}`: The total drift Hamiltonian including subsystem drifts and couplings
- `H_drives::Vector{SparseMatrixCSC{ComplexF64, Int}}`: All drive Hamiltonians (coupling drives + subsystem drives)
- `drive_bounds::Vector{Tuple{Float64, Float64}}`: Drive amplitude bounds for each control
- `n_drives::Int`: Total number of control drives
- `levels::Int`: Total dimension of the composite system (product of subsystem dimensions)
- `subsystem_levels::Vector{Int}`: Dimensions of each subsystem
- `subsystems::Vector{QuantumSystem}`: The individual quantum subsystems

See also [`LiftedOperators`](@ref), [`lift_operator`](@ref).

# Example
```julia
# Two qubits with ZZ coupling
sys1 = QuantumSystem([PAULIS[:X]], [(-1.0, 1.0)])
sys2 = QuantumSystem([PAULIS[:Y]], [(-1.0, 1.0)])
H_coupling = 0.1 * kron(PAULIS[:Z], PAULIS[:Z])
csys = CompositeQuantumSystem(H_coupling, [sys1, sys2], Float64[])
```
"""
struct CompositeQuantumSystem{F1<:Function, F2<:Function} <: AbstractQuantumSystem
    H::F1
    G::F2
    H_drift::SparseMatrixCSC{ComplexF64, Int}
    H_drives::Vector{SparseMatrixCSC{ComplexF64, Int}}
    drive_bounds::Vector{Tuple{Float64, Float64}}
    n_drives::Int
    levels::Int
    subsystem_levels::Vector{Int}
    subsystems::Vector{QuantumSystem}
end

"""
    CompositeQuantumSystem(
        H_drift::AbstractMatrix,
        H_drives::AbstractVector{<:AbstractMatrix},
        subsystems::AbstractVector{<:QuantumSystem},
        drive_bounds::DriveBounds
    )

Construct a CompositeQuantumSystem with coupling drift and drive terms.

# Arguments
- `H_drift::AbstractMatrix`: Coupling drift Hamiltonian (in full tensor product space)
- `H_drives::AbstractVector{<:AbstractMatrix}`: Coupling drive Hamiltonians
- `subsystems::AbstractVector{<:QuantumSystem}`: Vector of subsystems to compose
- `drive_bounds::DriveBounds`: Drive bounds for the coupling drives (subsystem bounds are inherited). Can be:
  - Tuples `(lower, upper)` for asymmetric bounds
  - Scalars which are interpreted as symmetric bounds `(-value, value)`

The total drift includes both the coupling drift and all subsystem drifts (automatically lifted).
The total drives include coupling drives followed by all subsystem drives (automatically lifted).

# Example
```julia
sys1 = QuantumSystem(PAULIS[:Z], [PAULIS[:X]], [1.0])
sys2 = QuantumSystem([PAULIS[:Y]], [1.0])
g12 = 0.1 * kron(PAULIS[:X], PAULIS[:X])  # coupling drift
csys = CompositeQuantumSystem(g12, Matrix{ComplexF64}[], [sys1, sys2], Float64[])
```
"""
function CompositeQuantumSystem(
    H_drift::AbstractMatrix{<:Number},
    H_drives::AbstractVector{<:AbstractMatrix{<:Number}},
    subsystems::AbstractVector{<:QuantumSystem},
    drive_bounds::DriveBounds
)
    # Normalize drive bounds to tuples
    drive_bounds = normalize_drive_bounds(drive_bounds)
    
    subsystem_levels = [sys.levels for sys ∈ subsystems]
    levels = prod(subsystem_levels)

    H_drift = sparse(H_drift)
    for (i, sys) ∈ enumerate(subsystems)
        H_drift += lift_operator(get_drift(sys), i, subsystem_levels)
    end

    H_drives = sparse.(H_drives)
    for (i, sys) ∈ enumerate(subsystems)
        for H_drive ∈ get_drives(sys)
            push!(H_drives, lift_operator(H_drive, i, subsystem_levels))
        end
    end

    n_drives = length(H_drives)
    H_drives = sparse.(H_drives)
    G_drift = sparse(Isomorphisms.G(H_drift))
    G_drives = sparse.(Isomorphisms.G.(H_drives))

    if n_drives == 0
        H = (u, t) -> H_drift
        G = (u, t) -> G_drift
    else
        H = (u, t) -> H_drift + sum(u .* H_drives)
        G = (u, t) -> G_drift + sum(u .* G_drives)
    end

    return CompositeQuantumSystem{typeof(H), typeof(G)}(
        H,
        G,
        H_drift,
        H_drives,
        drive_bounds,
        n_drives,
        levels,
        subsystem_levels,
        subsystems
    )
end

"""
    CompositeQuantumSystem(
        H_drives::AbstractVector{<:AbstractMatrix},
        subsystems::AbstractVector{<:QuantumSystem},
        drive_bounds::DriveBounds
    )

Convenience constructor for a composite system with coupling drives but no coupling drift.

# Arguments
- `H_drives::AbstractVector{<:AbstractMatrix}`: Coupling drive Hamiltonians
- `subsystems::AbstractVector{<:QuantumSystem}`: Vector of subsystems to compose
- `drive_bounds::DriveBounds`: Drive bounds for the coupling drives. Can be:
  - Tuples `(lower, upper)` for asymmetric bounds
  - Scalars which are interpreted as symmetric bounds `(-value, value)`

# Example
```julia
sys1 = QuantumSystem([PAULIS[:X]], [1.0])
sys2 = QuantumSystem([PAULIS[:Y]], [1.0])
g12 = 0.1 * kron(PAULIS[:X], PAULIS[:X])  # coupling drive
csys = CompositeQuantumSystem([g12], [sys1, sys2], [1.0])  # symmetric bound
```
"""
function CompositeQuantumSystem(
    H_drives::AbstractVector{<:AbstractMatrix{T}},
    subsystems::AbstractVector{<:QuantumSystem},
    drive_bounds::DriveBounds
) where T <: Number
    @assert !isempty(H_drives) "At least one drive is required"
    return CompositeQuantumSystem(
        spzeros(T, size(H_drives[1])),
        H_drives,
        subsystems,
        drive_bounds
    )
end

"""
    CompositeQuantumSystem(
        H_drift::AbstractMatrix,
        subsystems::AbstractVector{<:QuantumSystem},
        drive_bounds::DriveBounds
    )

Convenience constructor for a composite system with coupling drift but no coupling drives.

# Arguments
- `H_drift::AbstractMatrix`: Coupling drift Hamiltonian
- `subsystems::AbstractVector{<:QuantumSystem}`: Vector of subsystems to compose
- `drive_bounds::DriveBounds`: Drive bounds for the coupling drives (typically empty). Can be:
  - Tuples `(lower, upper)` for asymmetric bounds
  - Scalars which are interpreted as symmetric bounds `(-value, value)`

# Example
```julia
sys1 = QuantumSystem([PAULIS[:X]], [1.0])
sys2 = QuantumSystem([PAULIS[:Y]], [1.0])
H_coupling = 0.1 * kron(PAULIS[:Z], PAULIS[:Z])  # coupling drift
csys = CompositeQuantumSystem(H_coupling, [sys1, sys2], Float64[])
```
"""
function CompositeQuantumSystem(
    H_drift::AbstractMatrix{T},
    subsystems::AbstractVector{<:QuantumSystem},
    drive_bounds::DriveBounds
) where T <: Number
    return CompositeQuantumSystem(
        H_drift, 
        Matrix{T}[], 
        subsystems,
        drive_bounds
    )
end

"""
    CompositeQuantumSystem(
        subsystems::AbstractVector{<:QuantumSystem},
        drive_bounds::DriveBounds
    )

Convenience constructor for a composite system with no coupling terms (neither drift nor drives).

Use this when you have independent subsystems that you want to represent in a single
composite space, but without any direct coupling between them.

# Arguments
- `subsystems::AbstractVector{<:QuantumSystem}`: Vector of subsystems to compose
- `drive_bounds::DriveBounds`: Drive bounds for the coupling drives (typically empty). Can be:
  - Tuples `(lower, upper)` for asymmetric bounds
  - Scalars which are interpreted as symmetric bounds `(-value, value)`

# Example
```julia
sys1 = QuantumSystem([PAULIS[:X]], [1.0])
sys2 = QuantumSystem([PAULIS[:Y]], [1.0])
csys = CompositeQuantumSystem([sys1, sys2], Float64[])
```
"""
function CompositeQuantumSystem(
    subsystems::AbstractVector{<:QuantumSystem},
    drive_bounds::DriveBounds
)
    @assert !isempty(subsystems) "At least one subsystem is required"
    T = eltype(get_drift(subsystems[1]))
    levels = prod([sys.levels for sys ∈ subsystems])
    return CompositeQuantumSystem(
        spzeros(T, (levels, levels)), 
        Matrix{T}[], 
        subsystems,
        drive_bounds
    )
end

# ****************************************************************************** #

@testitem "Composite system" begin
    subsystem_levels = [4, 2, 2]
    sys1 = QuantumSystem(kron(PAULIS[:Z], PAULIS[:Z]), [kron(PAULIS[:X], PAULIS[:Y])], [(-1.0, 1.0)])
    sys2 = QuantumSystem([PAULIS[:Y], PAULIS[:Z]], [(-1.0, 1.0), (-1.0, 1.0)])
    sys3 = QuantumSystem(zeros(ComplexF64, 2, 2))
    subsystems = [sys1, sys2, sys3]
    g12 = 0.1 * lift_operator([kron(PAULIS[:X], PAULIS[:X]), PAULIS[:X]], [1, 2], subsystem_levels)
    g23 = 0.2 * lift_operator([PAULIS[:Y], PAULIS[:Y]], [2, 3], subsystem_levels)

    # Construct composite system
    csys = CompositeQuantumSystem(g12, [g23], [sys1, sys2, sys3], [(-1.0, 1.0)])
    @test csys.levels == prod(subsystem_levels)
    @test csys.n_drives == 1 + sum([sys.n_drives for sys ∈ subsystems])
    @test csys.subsystems == subsystems
    @test csys.subsystem_levels == subsystem_levels
    @test get_drift(csys) ≈ g12 + lift_operator(kron(PAULIS[:Z], PAULIS[:Z]), 1, subsystem_levels)
end

@testitem "Composite system from drift" begin
    using LinearAlgebra

    subsystem_levels = [2, 2]
    sys1 = QuantumSystem([PAULIS[:X], PAULIS[:Y]], [(-1.0, 1.0), (-1.0, 1.0)])
    sys2 = QuantumSystem([PAULIS[:Y], PAULIS[:Z]], [(-1.0, 1.0), (-1.0, 1.0)])
    subsystems = [sys1, sys2]
    g12 = 0.1 * kron(PAULIS[:X], PAULIS[:X])

    # Construct composite system from drift
    csys = CompositeQuantumSystem(g12, [sys1, sys2], Float64[])
    @test csys.levels == prod(subsystem_levels)
    @test csys.n_drives == sum([sys.n_drives for sys ∈ subsystems])
    @test csys.subsystems == subsystems
    @test csys.subsystem_levels == subsystem_levels
    @test get_drift(csys) ≈ g12
end

@testitem "Composite system from drives" begin
    subsystem_levels = [2, 2, 2]
    sys1 = QuantumSystem(PAULIS[:Z], [PAULIS[:X], PAULIS[:Y]], [(-1.0, 1.0), (-1.0, 1.0)])
    sys2 = QuantumSystem([PAULIS[:Y], PAULIS[:Z]], [(-1.0, 1.0), (-1.0, 1.0)])
    sys3 = QuantumSystem(zeros(ComplexF64, 2, 2))
    subsystems = [sys1, sys2, sys3]
    g12 = 0.1 * lift_operator([PAULIS[:X], PAULIS[:X]], [1, 2], subsystem_levels)
    g23 = 0.2 * lift_operator([PAULIS[:Y], PAULIS[:Y]], [2, 3], subsystem_levels)

    csys = CompositeQuantumSystem([g12, g23], [sys1, sys2, sys3], [(-1.0, 1.0), (-1.0, 1.0)])
    @test csys.levels == prod(subsystem_levels)
    @test csys.n_drives == 2 + sum([sys.n_drives for sys ∈ subsystems])
    @test csys.subsystems == subsystems
    @test csys.subsystem_levels == subsystem_levels
    @test get_drift(csys) ≈ lift_operator(PAULIS[:Z], 1, subsystem_levels)
end

@testitem "CompositeQuantumSystem drive_bounds conversion" begin
    using LinearAlgebra

    # Test scalar bounds are converted to symmetric tuples
    sys1 = QuantumSystem([PAULIS[:X]], [1.0])
    sys2 = QuantumSystem([PAULIS[:Y]], [1.0])
    subsystems = [sys1, sys2]
    g12 = 0.1 * kron(PAULIS[:X], PAULIS[:X])
    
    # Test with scalar bounds for coupling drives
    csys_scalar = CompositeQuantumSystem([g12], subsystems, [0.5])
    # The system should have subsystem bounds appended automatically
    # First drive is the coupling drive with bound (-0.5, 0.5)
    @test csys_scalar.drive_bounds[1] == (-0.5, 0.5)
    
    # Test with tuple bounds for coupling drives
    csys_tuple = CompositeQuantumSystem([g12], subsystems, [(-0.3, 0.7)])
    @test csys_tuple.drive_bounds[1] == (-0.3, 0.7)
    
    # Test with mixed bounds (scalars and tuples) - requires explicit type annotation
    g23 = 0.2 * kron(PAULIS[:Y], PAULIS[:Y])
    mixed_bounds = Union{Float64, Tuple{Float64,Float64}}[0.5, (-0.2, 0.8)]
    csys_mixed = CompositeQuantumSystem([g12, g23], subsystems, mixed_bounds)
    @test csys_mixed.drive_bounds[1] == (-0.5, 0.5)
    @test csys_mixed.drive_bounds[2] == (-0.2, 0.8)
end
