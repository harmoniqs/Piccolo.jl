# ----------------------------------------------------------------------------- #
# OpenQuantumSystem
# ----------------------------------------------------------------------------- #

"""
    OpenQuantumSystem <: AbstractQuantumSystem

A struct for storing open quantum dynamics.

# Fields
- `H::Function`: The Hamiltonian function: (u, t) -> H(u, t)
- `𝒢::Function`: The Lindbladian generator function: u -> 𝒢(u)
- `H_drift::SparseMatrixCSC{ComplexF64, Int}`: The drift Hamiltonian
- `H_drives::Vector{SparseMatrixCSC{ComplexF64, Int}}`: The drive Hamiltonians
- `drive_bounds::Vector{Tuple{Float64, Float64}}`: Drive amplitude bounds
- `n_drives::Int`: The number of control drives
- `levels::Int`: The number of levels in the system
- `dissipation_operators::Vector{SparseMatrixCSC{ComplexF64, Int}}`: The dissipation operators
- `time_dependent::Bool`: Whether the Hamiltonian has explicit time dependence

See also [`QuantumSystem`](@ref Piccolo.QuantumSystem).
"""
struct OpenQuantumSystem{F1<:Function,F2<:Function} <: AbstractQuantumSystem
    H::F1
    𝒢::F2
    H_drift::SparseMatrixCSC{ComplexF64,Int}
    H_drives::Vector{SparseMatrixCSC{ComplexF64,Int}}
    drive_bounds::Vector{Tuple{Float64,Float64}}
    n_drives::Int
    levels::Int
    dissipation_operators::Vector{SparseMatrixCSC{ComplexF64,Int}}
    time_dependent::Bool
end

"""
    OpenQuantumSystem(
        H_drift::AbstractMatrix{<:Number},
        H_drives::AbstractVector{<:AbstractMatrix{<:Number}},
        drive_bounds::DriveBounds;
        dissipation_operators::AbstractVector{<:AbstractMatrix{<:Number}}=Matrix{ComplexF64}[],
        time_dependent::Bool=false
    )
    OpenQuantumSystem(
        H_drift::AbstractMatrix{<:Number};
        dissipation_operators::AbstractVector{<:AbstractMatrix{<:Number}}=Matrix{ComplexF64}[],
        time_dependent::Bool=false
    )
    OpenQuantumSystem(
        H_drives::Vector{<:AbstractMatrix{<:Number}},
        drive_bounds::DriveBounds;
        dissipation_operators::AbstractVector{<:AbstractMatrix{<:Number}}=Matrix{ComplexF64}[],
        time_dependent::Bool=false
    )
    OpenQuantumSystem(
        H::Function, 
        drive_bounds::DriveBounds;
        dissipation_operators::Vector{<:AbstractMatrix{<:Number}}=Matrix{ComplexF64}[],
        time_dependent::Bool=false
    )
    OpenQuantumSystem(
        system::QuantumSystem; 
        dissipation_operators::Vector{<:AbstractMatrix{<:Number}}=Matrix{ComplexF64}[]
    )

Constructs an OpenQuantumSystem object from the drift and drive Hamiltonian terms and
dissipation operators.

# Drive Bounds
The `drive_bounds` parameter can be:
- Tuples `(lower, upper)` for asymmetric bounds
- Scalars which are interpreted as symmetric bounds `(-value, value)`
"""
function OpenQuantumSystem(
    H_drift::AbstractMatrix{<:Number},
    H_drives::Vector{<:AbstractMatrix{<:Number}},
    drive_bounds::Vector{<:Union{Tuple{Float64,Float64},Float64}};
    dissipation_operators::Vector{<:AbstractMatrix{<:Number}} = Matrix{ComplexF64}[],
    time_dependent::Bool = false,
)
    drive_bounds = normalize_drive_bounds(drive_bounds)

    H_drift_sparse = sparse(H_drift)
    𝒢_drift = Isomorphisms.G(Isomorphisms.ad_vec(H_drift_sparse))

    n_drives = length(H_drives)
    H_drives_sparse = sparse.(H_drives)
    𝒢_drives = [Isomorphisms.G(Isomorphisms.ad_vec(H_drive)) for H_drive in H_drives_sparse]

    # Build dissipator
    if isempty(dissipation_operators)
        𝒟 = spzeros(size(𝒢_drift))
    else
        𝒟 = sum(Isomorphisms.iso_D(sparse(L)) for L in dissipation_operators)
    end

    if n_drives == 0
        H = (u, t) -> H_drift_sparse
        𝒢 = u -> 𝒢_drift + 𝒟
    else
        H = (u, t) -> H_drift_sparse + sum(u .* H_drives_sparse)
        𝒢 = u -> 𝒢_drift + sum(u .* 𝒢_drives) + 𝒟
    end

    levels = size(H_drift, 1)

    return OpenQuantumSystem(
        H,
        𝒢,
        H_drift_sparse,
        H_drives_sparse,
        drive_bounds,
        n_drives,
        levels,
        sparse.(dissipation_operators),
        time_dependent,
    )
end

# Convenience constructors
"""
    OpenQuantumSystem(H_drives, drive_bounds; dissipation_operators=[], time_dependent=false)

Construct an OpenQuantumSystem with no drift.
"""
function OpenQuantumSystem(
    H_drives::Vector{<:AbstractMatrix{ℂ}},
    drive_bounds::Vector{<:Union{Tuple{Float64,Float64},Float64}};
    dissipation_operators::Vector{<:AbstractMatrix{<:Number}} = Matrix{ComplexF64}[],
    time_dependent::Bool = false,
) where {ℂ<:Number}
    @assert !isempty(H_drives) "At least one drive is required"
    return OpenQuantumSystem(
        spzeros(ℂ, size(H_drives[1])),
        H_drives,
        drive_bounds;
        dissipation_operators = dissipation_operators,
        time_dependent = time_dependent,
    )
end

"""
    OpenQuantumSystem(H_drift; dissipation_operators=[], time_dependent=false)

Construct an OpenQuantumSystem with only drift (no drives).
"""
function OpenQuantumSystem(
    H_drift::AbstractMatrix{ℂ};
    dissipation_operators::Vector{<:AbstractMatrix{<:Number}} = Matrix{ComplexF64}[],
    time_dependent::Bool = false,
) where {ℂ<:Number}
    return OpenQuantumSystem(
        H_drift,
        Matrix{ℂ}[],
        Float64[];
        dissipation_operators = dissipation_operators,
        time_dependent = time_dependent,
    )
end

"""
    OpenQuantumSystem(H::Function, drive_bounds; dissipation_operators=[], time_dependent=false)

Construct an OpenQuantumSystem from a Hamiltonian function.
"""
function OpenQuantumSystem(
    H::F,
    drive_bounds::Vector{<:Union{Tuple{Float64,Float64},Float64}};
    dissipation_operators::Vector{<:AbstractMatrix{ℂ}} = Matrix{ComplexF64}[],
    time_dependent::Bool = false,
) where {F<:Function,ℂ<:Number}

    drive_bounds = normalize_drive_bounds(drive_bounds)

    n_drives = length(drive_bounds)

    # Extract drift by evaluating with zero controls
    H_drift = H(zeros(n_drives), 0.0)
    levels = size(H_drift, 1)

    # Build dissipator
    if isempty(dissipation_operators)
        𝒟 = spzeros(ComplexF64, levels^2, levels^2)
    else
        𝒟 = sum(Isomorphisms.iso_D(sparse(L)) for L in dissipation_operators)
    end

    return OpenQuantumSystem(
        H,
        u -> Isomorphisms.G(Isomorphisms.ad_vec(sparse(H(u, 0.0)))) + 𝒟,
        sparse(H_drift),
        Vector{SparseMatrixCSC{ComplexF64,Int}}(),  # Empty drives vector for function-based systems
        drive_bounds,
        n_drives,
        levels,
        sparse.(dissipation_operators),
        time_dependent,
    )
end

"""
    OpenQuantumSystem(system::QuantumSystem; dissipation_operators=[])

Construct an OpenQuantumSystem from a QuantumSystem by adding dissipation operators.
"""
function OpenQuantumSystem(
    system::QuantumSystem;
    dissipation_operators::Vector{<:AbstractMatrix{<:Number}} = Matrix{ComplexF64}[],
)
    return OpenQuantumSystem(
        system.H_drift,
        system.H_drives,
        system.drive_bounds;
        dissipation_operators = dissipation_operators,
        time_dependent = system.time_dependent,
    )
end

# ******************************************************************************* #

@testitem "Open system creation" begin

    H_drift = PAULIS.Z
    # don't want drives == levels
    H_drives = [PAULIS.X]
    dissipation_operators = [PAULIS.Z, PAULIS.X]
    drive_bounds = [1.0]

    system = OpenQuantumSystem(
        H_drift,
        H_drives,
        drive_bounds,
        dissipation_operators = dissipation_operators,
    )
    @test system isa OpenQuantumSystem
    @test get_drift(system) == H_drift
    @test get_drives(system) == H_drives
    @test system.dissipation_operators == dissipation_operators

    # test dissipation
    𝒢_drift = Isomorphisms.G(Isomorphisms.ad_vec(H_drift))
    @test system.𝒢(zeros(system.n_drives)) != 𝒢_drift
end

@testitem "Open system alternate constructors" begin

    H_drift = PAULIS.Z
    # don't want drives == levels
    H_drives = [PAULIS.X]
    dissipation_operators = [PAULIS.Z, PAULIS.X]
    drive_bounds = [1.0]

    system = OpenQuantumSystem(
        H_drift,
        H_drives,
        drive_bounds,
        dissipation_operators = dissipation_operators,
    )
    @test system isa OpenQuantumSystem
    @test get_drift(system) == H_drift
    @test get_drives(system) == H_drives
    @test system.dissipation_operators == dissipation_operators

    # no drift
    system = OpenQuantumSystem(
        H_drives,
        drive_bounds,
        dissipation_operators = dissipation_operators,
    )
    @test system isa OpenQuantumSystem
    @test get_drift(system) == zeros(size(H_drift))
    @test get_drives(system) == H_drives
    @test system.dissipation_operators == dissipation_operators

    # no drives
    system = OpenQuantumSystem(H_drift, dissipation_operators = dissipation_operators)
    @test system isa OpenQuantumSystem
    @test system isa OpenQuantumSystem
    @test get_drift(system) == H_drift
    @test get_drives(system) == []
    @test system.dissipation_operators == dissipation_operators

    # function
    H = (u, t) -> PAULIS.Z + u[1] * PAULIS.X
    system =
        OpenQuantumSystem(H, drive_bounds, dissipation_operators = dissipation_operators)
    @test system isa OpenQuantumSystem
    @test get_drift(system) == H_drift
    @test get_drives(system) == H_drives
    @test system.dissipation_operators == dissipation_operators

    # from QuantumSystem
    qsys = QuantumSystem(H_drift, H_drives, drive_bounds)
    system = OpenQuantumSystem(qsys, dissipation_operators = dissipation_operators)
    @test system isa OpenQuantumSystem
    @test get_drift(system) == H_drift
    @test get_drives(system) == H_drives
    @test system.dissipation_operators == dissipation_operators

end

@testitem "OpenQuantumSystem drive_bounds conversion" begin

    # Test scalar bounds are converted to symmetric tuples
    H_drift = PAULIS.Z
    H_drives = [PAULIS.X, PAULIS.Y]
    dissipation_operators = [PAULIS.Z]

    # Test with scalar bounds
    sys_scalar = OpenQuantumSystem(
        H_drift,
        H_drives,
        [1.0, 1.5],
        dissipation_operators = dissipation_operators,
    )
    @test sys_scalar.drive_bounds == [(-1.0, 1.0), (-1.5, 1.5)]

    # Test with tuple bounds
    sys_tuple = OpenQuantumSystem(
        H_drift,
        H_drives,
        [(-0.5, 1.0), (-1.5, 0.5)],
        dissipation_operators = dissipation_operators,
    )
    @test sys_tuple.drive_bounds == [(-0.5, 1.0), (-1.5, 0.5)]

    # Test with mixed bounds (scalars and tuples) - requires explicit type annotation
    mixed_bounds = Union{Float64,Tuple{Float64,Float64}}[1.0, (-0.5, 1.5)]
    sys_mixed = OpenQuantumSystem(
        H_drift,
        H_drives,
        mixed_bounds,
        dissipation_operators = dissipation_operators,
    )
    @test sys_mixed.drive_bounds == [(-1.0, 1.0), (-0.5, 1.5)]

    # Test with function-based Hamiltonian
    H = (u, t) -> H_drift + sum(u .* H_drives)
    sys_func =
        OpenQuantumSystem(H, [0.8, 1.2], dissipation_operators = dissipation_operators)
    @test sys_func.drive_bounds == [(-0.8, 0.8), (-1.2, 1.2)]
end
