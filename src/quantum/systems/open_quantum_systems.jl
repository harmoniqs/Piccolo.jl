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
- `drives::Vector{AbstractDrive}`: Typed drive terms pairing operators with coefficient functions. For matrix-based constructors, auto-populated as `LinearDrive` objects. For function-based systems, empty.
- `H_drives::Vector{SparseMatrixCSC{ComplexF64, Int}}`: The drive Hamiltonians (backward compat). Populated for linear-only systems; empty when nonlinear drives are present or for function-based systems.
- `drive_bounds::Vector{Tuple{Float64, Float64}}`: Drive amplitude bounds
- `n_drives::Int`: The number of control drives
- `levels::Int`: The number of levels in the system
- `dissipation_operators::Vector{SparseMatrixCSC{ComplexF64, Int}}`: The dissipation operators
- `time_dependent::Bool`: Whether the Hamiltonian has explicit time dependence
- `global_params::NamedTuple`: Global parameters stored with the system (e.g., physical constants)

See also [`QuantumSystem`](@ref Piccolo.Quantum.QuantumSystems.QuantumSystem).
"""
struct OpenQuantumSystem{F1<:Function,F2<:Function,PT<:NamedTuple} <: AbstractQuantumSystem
    H::F1
    𝒢::F2
    H_drift::SparseMatrixCSC{ComplexF64,Int}
    drives::Vector{AbstractDrive}
    H_drives::Vector{SparseMatrixCSC{ComplexF64,Int}}
    drive_bounds::Vector{Tuple{Float64,Float64}}
    n_drives::Int
    levels::Int
    dissipation_operators::Vector{SparseMatrixCSC{ComplexF64,Int}}
    time_dependent::Bool
    global_params::PT
end

"""
    OpenQuantumSystem(
        H_drift::AbstractMatrix{<:Number},
        H_drives::AbstractVector{<:AbstractMatrix{<:Number}},
        drive_bounds::DriveBounds;
        dissipation_operators=Matrix{ComplexF64}[],
        time_dependent::Bool=false,
        global_params::NamedTuple=NamedTuple()
    )
    OpenQuantumSystem(
        H_drift::AbstractMatrix{<:Number};
        dissipation_operators=Matrix{ComplexF64}[],
        time_dependent::Bool=false,
        global_params::NamedTuple=NamedTuple()
    )
    OpenQuantumSystem(
        H_drives::Vector{<:AbstractMatrix{<:Number}},
        drive_bounds::DriveBounds;
        dissipation_operators=Matrix{ComplexF64}[],
        time_dependent::Bool=false,
        global_params::NamedTuple=NamedTuple()
    )
    OpenQuantumSystem(
        H::Function,
        drive_bounds::DriveBounds;
        dissipation_operators=Matrix{ComplexF64}[],
        time_dependent::Bool=false,
        global_params::NamedTuple=NamedTuple()
    )
    OpenQuantumSystem(
        system::QuantumSystem;
        dissipation_operators=Matrix{ComplexF64}[]
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
    global_params::NamedTuple = NamedTuple(),
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

    # Build LinearDrive objects from H_drives
    drives = AbstractDrive[LinearDrive(H_drives_sparse[i], i) for i in 1:n_drives]

    return OpenQuantumSystem(
        H,
        𝒢,
        H_drift_sparse,
        drives,
        H_drives_sparse,
        drive_bounds,
        n_drives,
        levels,
        sparse.(dissipation_operators),
        time_dependent,
        _float_params(global_params),
    )
end

"""
    OpenQuantumSystem(
        H_drift::AbstractMatrix,
        drives::Vector{<:AbstractDrive},
        drive_bounds::DriveBounds;
        dissipation_operators=Matrix{ComplexF64}[],
        time_dependent::Bool=false,
        global_params::NamedTuple=NamedTuple()
    )

Construct an OpenQuantumSystem from a drift Hamiltonian and typed drive terms.

This constructor supports both linear and nonlinear drives. The control dimension
is determined by `length(drive_bounds)`, which may differ from `length(drives)`
when nonlinear drives combine multiple controls.

# Example
```julia
L_ops = [sqrt(γ) * annihilate(levels)]
drives = [
    LinearDrive(sparse(σx), 1),
    NonlinearDrive(σz, u -> u[1]^2 + u[2]^2),
]
sys = OpenQuantumSystem(H_drift, drives, [1.0, 1.0]; dissipation_operators=L_ops)
```
"""
function OpenQuantumSystem(
    H_drift::AbstractMatrix{<:Number},
    drives::Vector{<:AbstractDrive},
    drive_bounds::Vector{<:Union{Tuple{Float64,Float64},Float64}};
    dissipation_operators::Vector{<:AbstractMatrix{<:Number}} = Matrix{ComplexF64}[],
    time_dependent::Bool = false,
    global_params::NamedTuple = NamedTuple(),
)
    drive_bounds = normalize_drive_bounds(drive_bounds)

    # Check that H_drift is Hermitian
    @assert is_hermitian(H_drift) "Drift Hamiltonian H_drift is not Hermitian"

    # Check that all drive operators are Hermitian
    for (i, d) in enumerate(drives)
        @assert is_hermitian(d.H) "Drive operator drives[$i].H is not Hermitian"
    end

    H_drift_sparse = sparse(ComplexF64.(H_drift))
    n_drives = length(drive_bounds)
    levels = size(H_drift, 1)

    # Validate LinearDrive indices
    for (i, d) in enumerate(drives)
        if d isa LinearDrive
            @assert 1 <= d.index <= n_drives "LinearDrive at drives[$i] has index $(d.index) but control dimension is $n_drives (length of drive_bounds)"
        end
    end

    # Validate NonlinearDrive Jacobians
    for d in drives
        if d isa NonlinearDrive
            validate_drive_jacobian(d, n_drives)
        end
    end

    # Build H(u,t) and 𝒢(u) from drives
    𝒢_drift = Isomorphisms.G(Isomorphisms.ad_vec(H_drift_sparse))
    𝒢_drive_mats = [Isomorphisms.G(Isomorphisms.ad_vec(d.H)) for d in drives]

    # Build dissipator
    if isempty(dissipation_operators)
        𝒟 = spzeros(size(𝒢_drift))
    else
        𝒟 = sum(Isomorphisms.iso_D(sparse(L)) for L in dissipation_operators)
    end

    if isempty(drives)
        H_fn = (u, t) -> H_drift_sparse
        𝒢_fn = u -> 𝒢_drift + 𝒟
    else
        H_fn = (u, t) -> H_drift_sparse + sum(drive_coeff(d, u) * d.H for d in drives)
        𝒢_fn = u -> 𝒢_drift + sum(drive_coeff(d, u) * 𝒢_d for (d, 𝒢_d) in zip(drives, 𝒢_drive_mats)) + 𝒟
    end

    # H_drives: populated only for purely linear systems (backward compat)
    if has_nonlinear_drives(drives)
        H_drives_compat = Vector{SparseMatrixCSC{ComplexF64,Int}}()
    else
        H_drives_compat = SparseMatrixCSC{ComplexF64,Int}[d.H for d in drives]
    end

    return OpenQuantumSystem(
        H_fn,
        𝒢_fn,
        H_drift_sparse,
        collect(AbstractDrive, drives),
        H_drives_compat,
        drive_bounds,
        n_drives,
        levels,
        sparse.(dissipation_operators),
        time_dependent,
        _float_params(global_params),
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
    global_params::NamedTuple = NamedTuple(),
) where {ℂ<:Number}
    @assert !isempty(H_drives) "At least one drive is required"
    return OpenQuantumSystem(
        spzeros(ℂ, size(H_drives[1])),
        H_drives,
        drive_bounds;
        dissipation_operators = dissipation_operators,
        time_dependent = time_dependent,
        global_params = global_params,
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
    global_params::NamedTuple = NamedTuple(),
) where {ℂ<:Number}
    return OpenQuantumSystem(
        H_drift,
        Matrix{ℂ}[],
        Float64[];
        dissipation_operators = dissipation_operators,
        time_dependent = time_dependent,
        global_params = global_params,
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
    global_params::NamedTuple = NamedTuple(),
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
        AbstractDrive[],                                  # No drives for function-based systems
        Vector{SparseMatrixCSC{ComplexF64,Int}}(),        # No H_drives for function-based systems
        drive_bounds,
        n_drives,
        levels,
        sparse.(dissipation_operators),
        time_dependent,
        _float_params(global_params),
    )
end

"""
    OpenQuantumSystem(system::QuantumSystem; dissipation_operators=[])

Construct an OpenQuantumSystem from a QuantumSystem by adding dissipation operators.

When the `QuantumSystem` has typed drives (including nonlinear), they are preserved
in the resulting `OpenQuantumSystem`.
"""
function OpenQuantumSystem(
    system::QuantumSystem;
    dissipation_operators::Vector{<:AbstractMatrix{<:Number}} = Matrix{ComplexF64}[],
)
    if !isempty(system.drives)
        # Use drives-based constructor (handles both linear and nonlinear)
        return OpenQuantumSystem(
            system.H_drift,
            system.drives,
            system.drive_bounds;
            dissipation_operators = dissipation_operators,
            time_dependent = system.time_dependent,
            global_params = system.global_params,
        )
    elseif !isempty(system.H_drives)
        return OpenQuantumSystem(
            system.H_drift,
            system.H_drives,
            system.drive_bounds;
            dissipation_operators = dissipation_operators,
            time_dependent = system.time_dependent,
            global_params = system.global_params,
        )
    else
        # Function-based system
        return OpenQuantumSystem(
            system.H,
            system.drive_bounds;
            dissipation_operators = dissipation_operators,
            time_dependent = system.time_dependent,
            global_params = system.global_params,
        )
    end
end

# ----------------------------------------------------------------------------- #
# Compact Lindbladian generators
# ----------------------------------------------------------------------------- #

"""
    compact_lindbladian_generators(sys::OpenQuantumSystem)

Compute the compact Lindbladian generators for use with the compact density
isomorphism. Returns `(𝒢c_drift, 𝒢c_drives)` where:

- `𝒢c_drift = P * (𝒢_drift + 𝒟) * L` — compact drift + dissipation generator (n² × n²)
- `𝒢c_drives[i]` — compact generator for each drive term (each n² × n²)

For linear-only systems, `𝒢c_drives[i]` corresponds to control `u[i]`.
For systems with nonlinear drives, `𝒢c_drives[i]` corresponds to drive term `i`,
and must be weighted by `drive_coeff(sys.drives[i], u)` instead of `u[i]`.
"""
function compact_lindbladian_generators(sys::OpenQuantumSystem)
    n = sys.levels

    # Reconstruct full Lindbladian components from stored fields
    𝒢_drift = Isomorphisms.G(Isomorphisms.ad_vec(sys.H_drift))

    # Use drives when available, fall back to H_drives
    if !isempty(sys.drives)
        𝒢_drive_terms = [Isomorphisms.G(Isomorphisms.ad_vec(d.H)) for d in sys.drives]
    else
        𝒢_drive_terms = [Isomorphisms.G(Isomorphisms.ad_vec(H_d)) for H_d in sys.H_drives]
    end

    if isempty(sys.dissipation_operators)
        𝒟 = spzeros(size(𝒢_drift))
    else
        𝒟 = sum(Isomorphisms.iso_D(L) for L in sys.dissipation_operators)
    end

    L = Isomorphisms.density_lift_matrix(n)
    P = Isomorphisms.density_projection_matrix(n)

    𝒢c_drift = P * (𝒢_drift + 𝒟) * L
    𝒢c_drives = [P * 𝒢d * L for 𝒢d in 𝒢_drive_terms]

    return 𝒢c_drift, 𝒢c_drives
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

    # drives field should be auto-populated with LinearDrives
    @test length(system.drives) == 1
    @test system.drives[1] isa LinearDrive
    @test system.drives[1].index == 1

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
    @test length(system.drives) == 1
    @test system.drives[1] isa LinearDrive

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

@testitem "OpenQuantumSystem with nonlinear drives" begin
    using SparseArrays
    using LinearAlgebra

    H_drift = PAULIS.Z
    dissipation_operators = [0.1 * PAULIS[:X]]  # simple dephasing

    # Mixed drives: linear + nonlinear
    drives = AbstractDrive[
        LinearDrive(sparse(ComplexF64.(PAULIS.X)), 1),
        LinearDrive(sparse(ComplexF64.(PAULIS.Y)), 2),
        NonlinearDrive(PAULIS.Z, u -> u[1]^2 + u[2]^2),  # auto-Jacobian
    ]

    sys = OpenQuantumSystem(
        H_drift, drives, [1.0, 1.0];
        dissipation_operators = dissipation_operators,
    )

    @test sys isa OpenQuantumSystem
    @test sys.n_drives == 2
    @test length(sys.drives) == 3
    @test has_nonlinear_drives(sys.drives)
    @test isempty(sys.H_drives)  # nonlinear → H_drives empty

    # H function should work correctly
    u_test = [0.3, 0.7]
    H_result = sys.H(u_test, 0.0)
    coeff_nl = 0.3^2 + 0.7^2
    H_expected = PAULIS.Z + 0.3 * PAULIS.X + 0.7 * PAULIS.Y + coeff_nl * PAULIS.Z
    @test norm(H_result - H_expected) < 1e-10

    # 𝒢 function should include dissipation
    𝒢_result = sys.𝒢(u_test)
    @test size(𝒢_result) == (8, 8)  # 2² × 2² × 2 (iso)

    # From QuantumSystem with nonlinear drives
    qsys = QuantumSystem(H_drift, drives, [1.0, 1.0])
    osys = OpenQuantumSystem(qsys; dissipation_operators = dissipation_operators)
    @test osys isa OpenQuantumSystem
    @test length(osys.drives) == 3
    @test has_nonlinear_drives(osys.drives)
end
