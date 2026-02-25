
# ----------------------------------------------------------------------------- #
# QuantumSystem
# ----------------------------------------------------------------------------- #

"""
    is_hermitian(H::AbstractMatrix; tol=1e-10)

Check if a matrix is Hermitian within a tolerance.
"""
function is_hermitian(H::AbstractMatrix; tol = 1e-10)
    return norm(H - H') < tol
end

"""
    QuantumSystem <: AbstractQuantumSystem

A struct for storing quantum dynamics.

# Fields
- `H::Function`: The Hamiltonian function: (u, t) -> H(u, t), where u is the control vector and t is time
- `G::Function`: The isomorphic generator function: (u, t) -> G(u, t), including the Hamiltonian mapped to superoperator space
- `H_drift::SparseMatrixCSC{ComplexF64, Int}`: The drift Hamiltonian (time-independent component)
- `H_drives::Vector{SparseMatrixCSC{ComplexF64, Int}}`: The drive Hamiltonians (control-dependent components)
- `drive_bounds::Vector{Tuple{Float64, Float64}}`: Drive amplitude bounds for each control (lower, upper)
- `n_drives::Int`: The number of control drives in the system
- `levels::Int`: The number of levels (dimension) in the system
- `time_dependent::Bool`: Whether the Hamiltonian has explicit time dependence beyond control modulation
- `global_params::NamedTuple`: Global parameters that the Hamiltonian may depend on (e.g., (δ=0.5, Ω=1.0))

See also [`OpenQuantumSystem`](@ref), [`VariationalQuantumSystem`](@ref).
"""
struct QuantumSystem{F1<:Function,F2<:Function,PT<:NamedTuple} <: AbstractQuantumSystem
    H::F1
    G::F2
    H_drift::SparseMatrixCSC{ComplexF64,Int}
    drives::Vector{AbstractDrive}
    H_drives::Vector{SparseMatrixCSC{ComplexF64,Int}}
    drive_bounds::Vector{Tuple{Float64,Float64}}
    n_drives::Int
    levels::Int
    time_dependent::Bool
    global_params::PT
end

"""
    QuantumSystem(H::Function, drive_bounds::Vector; time_dependent::Bool=false)

Construct a QuantumSystem from a Hamiltonian function.

# Arguments
- `H::Function`: Hamiltonian function with signature (u, t) -> H(u, t) where:
  - `u` is a vector containing `[controls..., globals...]` (if system has global parameters)
  - For matrix-based systems, only the first n_drives elements are used for controls
  - For function-based systems, handle globals via closure or by accessing u beyond control indices
  - `t` is time
- `drive_bounds::DriveBounds`: Drive amplitude bounds for each control. Can be:
  - Tuples `(lower, upper)` for asymmetric bounds
  - Scalars which are interpreted as symmetric bounds `(-value, value)`

# Keyword Arguments
- `time_dependent::Bool=false`: Set to `true` if the Hamiltonian has explicit time dependence (e.g., cos(ωt) modulation)
- `global_params::NamedTuple=NamedTuple()`: Global parameters stored with the system for bookkeeping

# Example
```julia
# Define a time-dependent Hamiltonian
H = (u, t) -> PAULIS[:Z] + u[1] * cos(ω * t) * PAULIS[:X]
sys = QuantumSystem(H, [(-1.0, 1.0)]; time_dependent=true)
```
"""
function QuantumSystem(
    H::Function,
    drive_bounds::Vector{<:Union{Tuple{Float64,Float64},Float64}};
    time_dependent::Bool = false,
    global_params::NamedTuple = NamedTuple(),
)
    drive_bounds = normalize_drive_bounds(drive_bounds)

    n_drives = length(drive_bounds)
    n_globals = length(global_params)

    # Build test vector u = [controls..., globals...]
    u_zeros =
        n_globals > 0 ? vcat(zeros(n_drives), collect(values(global_params))) :
        zeros(n_drives)

    # Extract drift by evaluating with zero controls (and initial globals if present)
    H_drift = H(u_zeros, 0.0)
    levels = size(H_drift, 1)

    # Check that H_drift is Hermitian
    @assert is_hermitian(H_drift) "Drift Hamiltonian H(u=0, t=0) is not Hermitian"

    # Check that Hamiltonian is Hermitian for sample control values
    u_test_controls = [b isa Tuple ? (b[1] + b[2]) / 2 : 0.0 for b in drive_bounds]
    u_test =
        n_globals > 0 ? vcat(u_test_controls, collect(values(global_params))) :
        u_test_controls
    H_test = H(u_test, 0.0)
    @assert is_hermitian(H_test) "Hamiltonian H(u, t=0) is not Hermitian for test control values u=$u_test"

    return QuantumSystem(
        (u, t) -> H(u, t),
        (u, t) -> Isomorphisms.G(H(u, t)),
        sparse(H_drift),
        AbstractDrive[],                                  # No drives for function-based systems
        Vector{SparseMatrixCSC{ComplexF64,Int}}(),        # No H_drives for function-based systems
        drive_bounds,
        n_drives,
        levels,
        time_dependent,
        _float_params(global_params),
    )
end

"""
    QuantumSystem(
        H_drift::AbstractMatrix{<:Number},
        H_drives::Vector{<:AbstractMatrix{<:Number}},
        drive_bounds::Vector{<:Union{Tuple{Float64, Float64}, Float64}};
        time_dependent::Bool=false
    )

Construct a QuantumSystem from drift and drive Hamiltonian terms.

# Arguments
- `H_drift::AbstractMatrix`: The drift (time-independent) Hamiltonian
- `H_drives::Vector{<:AbstractMatrix}`: Vector of drive Hamiltonians, one for each control
- `drive_bounds::DriveBounds`: Drive amplitude bounds for each control. Can be:
  - Tuples `(lower, upper)` for asymmetric bounds
  - Scalars which are interpreted as symmetric bounds `(-value, value)`

# Keyword Arguments
- `time_dependent::Bool=false`: Set to `true` if using time-dependent modulation (typically handled at a higher level)
- `global_params::NamedTuple=NamedTuple()`: Global parameters stored with the system. Note: for matrix-based systems,
  matrices are fixed at construction, so global_params are mainly for storage/bookkeeping and later updates via `update_global_params!`

The resulting Hamiltonian is: H(u, t) = H_drift + Σᵢ uᵢ * H_drives[i]

# Example
```julia
sys = QuantumSystem(
    PAULIS[:Z],                    # drift
    [PAULIS[:X], PAULIS[:Y]],      # drives
    [1.0, 1.0]                     # symmetric bounds: [(-1.0, 1.0), (-1.0, 1.0)]
)
```
"""
function QuantumSystem(
    H_drift::AbstractMatrix{<:Number},
    H_drives::Vector{<:AbstractMatrix{<:Number}},
    drive_bounds::Vector{<:Union{Tuple{Float64,Float64},Float64}};
    time_dependent::Bool = false,
    global_params::NamedTuple = NamedTuple(),
)
    drive_bounds = [b isa Tuple ? b : (-b, b) for b in drive_bounds]

    # Check that H_drift is Hermitian
    @assert is_hermitian(H_drift) "Drift Hamiltonian H_drift is not Hermitian"

    # Check that all drive Hamiltonians are Hermitian
    for (i, H_drive) in enumerate(H_drives)
        @assert is_hermitian(H_drive) "Drive Hamiltonian H_drives[$i] is not Hermitian"
    end

    H_drift = sparse(H_drift)
    G_drift = sparse(Isomorphisms.G(H_drift))

    n_drives = length(H_drives)
    H_drives = sparse.(H_drives)
    G_drives = sparse.(Isomorphisms.G.(H_drives))

    # Note: u may contain [controls..., globals...] where globals are extra elements beyond n_drives
    # The integrator handles splitting u appropriately
    if n_drives == 0
        H = (u, t) -> H_drift
        G = (u, t) -> G_drift
    else
        H = (u, t) -> H_drift + sum(view(u, 1:n_drives) .* H_drives)
        G = (u, t) -> G_drift + sum(view(u, 1:n_drives) .* G_drives)
    end

    levels = size(H_drift, 1)

    # Build LinearDrive objects from H_drives
    drives = AbstractDrive[LinearDrive(H_drives[i], i) for i in 1:n_drives]

    return QuantumSystem(
        H,
        G,
        H_drift,
        drives,
        H_drives,
        drive_bounds,
        n_drives,
        levels,
        time_dependent,
        _float_params(global_params),
    )
end

# Convenience constructors
"""
    QuantumSystem(H_drives::Vector{<:AbstractMatrix}, drive_bounds::Vector; time_dependent::Bool=false)

Convenience constructor for a system with no drift Hamiltonian (H_drift = 0).

# Arguments
- `H_drives::Vector{<:AbstractMatrix}`: Vector of drive Hamiltonians
- `drive_bounds::DriveBounds`: Drive amplitude bounds for each control. Can be:
  - Tuples `(lower, upper)` for asymmetric bounds
  - Scalars which are interpreted as symmetric bounds `(-value, value)`

# Example
```julia
# Using scalars for symmetric bounds
sys = QuantumSystem([PAULIS[:X], PAULIS[:Y]], [1.0, 1.0])
# Equivalent to: drive_bounds = [(-1.0, 1.0), (-1.0, 1.0)]
```
"""
function QuantumSystem(
    H_drives::Vector{<:AbstractMatrix{ℂ}},
    drive_bounds::Vector{<:Union{Tuple{Float64,Float64},Float64}};
    time_dependent::Bool = false,
    global_params::NamedTuple = NamedTuple(),
) where {ℂ<:Number}
    @assert !isempty(H_drives) "At least one drive is required"
    return QuantumSystem(
        spzeros(ℂ, size(H_drives[1])),
        H_drives,
        drive_bounds;
        time_dependent = time_dependent,
        global_params = global_params,
    )
end

"""
    QuantumSystem(H_drift::AbstractMatrix; time_dependent::Bool=false)

Convenience constructor for a system with only a drift Hamiltonian (no drives).

# Example
```julia
sys = QuantumSystem(PAULIS[:Z])
```
"""
function QuantumSystem(
    H_drift::AbstractMatrix{ℂ};
    time_dependent::Bool = false,
    global_params::NamedTuple = NamedTuple(),
) where {ℂ<:Number}
    QuantumSystem(
        H_drift,
        Matrix{ℂ}[],
        Float64[];
        time_dependent = time_dependent,
        global_params = global_params,
    )
end

"""
    QuantumSystem(
        H_drift::AbstractMatrix,
        drives::Vector{<:AbstractDrive},
        drive_bounds::Vector;
        time_dependent::Bool=false,
        global_params::NamedTuple=NamedTuple()
    )

Construct a QuantumSystem from a drift Hamiltonian and typed drive terms.

This constructor supports both linear and nonlinear drives. The control dimension
is determined by `length(drive_bounds)`, which may differ from `length(drives)`
when nonlinear drives combine multiple controls.

The resulting Hamiltonian is: H(u, t) = H_drift + Σ_d drive_coeff(d, u) * d.H

# Arguments
- `H_drift::AbstractMatrix`: The drift (time-independent) Hamiltonian
- `drives::Vector{<:AbstractDrive}`: Vector of drive terms (LinearDrive or NonlinearDrive)
- `drive_bounds::DriveBounds`: Bounds for each physical control. Length = control dimension.

# Example: Displaced frame with nonlinear |α|² term
```julia
sys = QuantumSystem(
    H_drift,
    [
        LinearDrive(sparse(σx), 1),                    # u[1] * σx (qubit I)
        LinearDrive(sparse(σy), 2),                    # u[2] * σy (qubit Q)
        LinearDrive(sparse(χ * Xa * σz), 3),           # u[3] * χ·Xa·σz (displacement I)
        LinearDrive(sparse(χ * Pa * σz), 4),           # u[4] * χ·Pa·σz (displacement Q)
        NonlinearDrive(                                 # (u[3]²+u[4]²) * χ·σz/2
            sparse(χ * σz / 2),
            u -> u[3]^2 + u[4]^2,
            (u, j) -> j == 3 ? 2u[3] : j == 4 ? 2u[4] : 0.0
        ),
    ],
    [Ω_bound, Ω_bound, α_bound, α_bound]  # 4 physical controls
)
```
"""
function QuantumSystem(
    H_drift::AbstractMatrix{<:Number},
    drives::Vector{<:AbstractDrive},
    drive_bounds::Vector{<:Union{Tuple{Float64,Float64},Float64}};
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

    H_drift = sparse(ComplexF64.(H_drift))
    n_drives = length(drive_bounds)
    levels = size(H_drift, 1)

    # Build H(u,t) and G(u,t) from drives
    G_drift = sparse(Isomorphisms.G(H_drift))
    G_drive_mats = [sparse(Isomorphisms.G(d.H)) for d in drives]

    if isempty(drives)
        H_fn = (u, t) -> H_drift
        G_fn = (u, t) -> G_drift
    else
        H_fn = (u, t) -> H_drift + sum(drive_coeff(d, u) * d.H for d in drives)
        G_fn = (u, t) -> G_drift + sum(drive_coeff(d, u) * G_d for (d, G_d) in zip(drives, G_drive_mats))
    end

    # H_drives: populated only for purely linear systems (backward compat)
    if has_nonlinear_drives(drives)
        H_drives_compat = Vector{SparseMatrixCSC{ComplexF64,Int}}()
    else
        H_drives_compat = SparseMatrixCSC{ComplexF64,Int}[d.H for d in drives]
    end

    drives_stored = AbstractDrive[d for d in drives]

    return QuantumSystem(
        H_fn,
        G_fn,
        H_drift,
        drives_stored,
        H_drives_compat,
        drive_bounds,
        n_drives,
        levels,
        time_dependent,
        _float_params(global_params),
    )
end

# ******************************************************************************* #

@testitem "System creation" begin
    using SparseArrays: sparse

    H_drift = PAULIS.Z
    H_drives = [PAULIS.X, PAULIS.Y]
    n_drives = length(H_drives)
    u_bounds = ones(n_drives)

    system = QuantumSystem(H_drift, H_drives, u_bounds)
    @test system isa QuantumSystem
    @test get_drift(system) == H_drift
    @test get_drives(system) == H_drives

    # repeat with a bigger system
    H_drift = kron(PAULIS.Z, PAULIS.Z)
    H_drives = [
        kron(PAULIS.X, PAULIS.I),
        kron(PAULIS.I, PAULIS.X),
        kron(PAULIS.Y, PAULIS.I),
        kron(PAULIS.I, PAULIS.Y),
    ]
    n_drives = length(H_drives)
    u_bounds = ones(n_drives)

    system = QuantumSystem(H_drift, H_drives, u_bounds)
    @test system isa QuantumSystem
    @test get_drift(system) == H_drift
    @test get_drives(system) == H_drives
end

@testitem "No drift system creation" begin
    using SparseArrays: spzeros

    H_drift = zeros(ComplexF64, 2, 2)
    H_drives = [PAULIS.X, PAULIS.Y]
    u_bounds = [1.0, 1.0]

    sys1 = QuantumSystem(H_drift, H_drives, u_bounds)
    sys2 = QuantumSystem(H_drives, u_bounds)

    @test get_drift(sys1) == get_drift(sys2) == H_drift
    @test get_drives(sys1) == get_drives(sys2) == H_drives
end

@testitem "No drive system creation" begin

    H_drift = PAULIS.Z
    H_drives = Matrix{ComplexF64}[]
    u_bounds = Float64[]

    sys1 = QuantumSystem(H_drift, H_drives, u_bounds)
    sys2 = QuantumSystem(H_drift)

    @test get_drift(sys1) == get_drift(sys2) == H_drift
    @test get_drives(sys1) == get_drives(sys2) == H_drives
end

@testitem "System creation with Hamiltonian function" begin

    # test one drive

    H_drift = PAULIS.Z
    H_drives = [PAULIS.X]

    system = QuantumSystem((a, t) -> H_drift + sum(a .* H_drives), [1.0])
    @test system isa QuantumSystem
    @test get_drift(system) == H_drift
    @test get_drives(system) == H_drives

    # test no drift + three drives

    H_drives = [PAULIS.X, PAULIS.Y, PAULIS.Z]
    system =
        QuantumSystem((a, t) -> sum(a .* H_drives), [1.0, 1.0, 1.0], time_dependent = false)
    @test system isa QuantumSystem
    @test get_drift(system) == zeros(2, 2)
    @test get_drives(system) == H_drives
end

@testitem "Hermiticity check" begin
    using LinearAlgebra: I

    # Non-Hermitian drift should fail
    H_drift_bad = [1.0 1.0im; 0.0 1.0]  # Not Hermitian
    @test_throws AssertionError QuantumSystem(H_drift_bad, [PAULIS.X], [1.0])

    # Non-Hermitian drive should fail  
    H_drive_bad = [1.0 1.0im; 0.0 1.0]  # Not Hermitian
    @test_throws AssertionError QuantumSystem(PAULIS.Z, [H_drive_bad], [1.0])

    # Hermitian matrices should succeed
    H_drift = PAULIS.Z
    H_drives = [PAULIS.X, PAULIS.Y]
    sys = QuantumSystem(H_drift, H_drives, [1.0, 1.0])
    @test sys isa QuantumSystem

    # Function-based: non-Hermitian should fail
    H_bad = (u, t) -> [1.0 1.0im; 0.0 1.0]
    @test_throws AssertionError QuantumSystem(H_bad, [1.0])

    # Function-based: Hermitian should succeed
    H_good = (u, t) -> PAULIS.Z + u[1] * PAULIS.X
    sys2 = QuantumSystem(H_good, [1.0])
    @test sys2 isa QuantumSystem
end

@testitem "System creation variants" begin

    # Test with drift, drives, and bounds
    H_drift = PAULIS.Z
    H_drives = [PAULIS.X, PAULIS.Y]
    u_bounds = [1.0, 1.0]

    sys = QuantumSystem(H_drift, H_drives, u_bounds)
    @test sys isa QuantumSystem
    @test get_drift(sys) == H_drift
    @test get_drives(sys) == H_drives
    @test sys.n_drives == 2

    # drives field should be auto-populated with LinearDrives
    @test length(sys.drives) == 2
    @test all(d -> d isa LinearDrive, sys.drives)
    @test sys.drives[1].index == 1
    @test sys.drives[2].index == 2

    # Test with drives only (no drift)
    sys2 = QuantumSystem(H_drives, u_bounds)
    @test sys2 isa QuantumSystem
    @test get_drift(sys2) == zeros(ComplexF64, 2, 2)
    @test get_drives(sys2) == H_drives

    # Test with drift only (no drives)
    sys3 = QuantumSystem(H_drift)
    @test sys3 isa QuantumSystem
    @test get_drift(sys3) == H_drift
    @test isempty(get_drives(sys3))
    @test sys3.n_drives == 0
    @test isempty(sys3.drives)
end

@testitem "Drives-based system creation" begin
    using SparseArrays
    using LinearAlgebra

    # Purely linear drives — should be equivalent to matrix-based constructor
    H_drift = PAULIS.Z
    drives = [
        LinearDrive(sparse(ComplexF64.(PAULIS.X)), 1),
        LinearDrive(sparse(ComplexF64.(PAULIS.Y)), 2),
    ]
    sys = QuantumSystem(H_drift, drives, [1.0, 1.0])
    @test sys isa QuantumSystem
    @test sys.n_drives == 2
    @test length(sys.drives) == 2
    @test !isempty(sys.H_drives)  # linear-only → H_drives populated

    u_test = [0.3, 0.7]
    H_expected = PAULIS.Z + 0.3 * PAULIS.X + 0.7 * PAULIS.Y
    @test norm(sys.H(u_test, 0.0) - H_expected) < 1e-10

    # Mixed drives (linear + nonlinear) — H_drives should be empty
    nonlinear_drive = NonlinearDrive(
        sparse(ComplexF64.([1.0 0.0; 0.0 -1.0])),
        u -> u[1]^2 + u[2]^2,
        (u, j) -> j == 1 ? 2u[1] : j == 2 ? 2u[2] : 0.0
    )
    mixed_drives = [drives..., nonlinear_drive]
    sys2 = QuantumSystem(H_drift, mixed_drives, [1.0, 1.0])

    @test sys2 isa QuantumSystem
    @test sys2.n_drives == 2  # control dimension = 2
    @test length(sys2.drives) == 3  # 3 drive terms
    @test isempty(sys2.H_drives)  # nonlinear → H_drives empty

    # H function should work correctly
    H_result = sys2.H(u_test, 0.0)
    coeff_nonlinear = 0.3^2 + 0.7^2  # = 0.58
    H_expected2 = PAULIS.Z + 0.3 * PAULIS.X + 0.7 * PAULIS.Y + coeff_nonlinear * [1.0 0.0; 0.0 -1.0]
    @test norm(H_result - H_expected2) < 1e-10

    # has_nonlinear_drives
    @test !has_nonlinear_drives(sys.drives)
    @test has_nonlinear_drives(sys2.drives)
end

@testitem "Global parameters" begin
    using LinearAlgebra

    # Test default empty global_params
    H_drives = [PAULIS[:X], PAULIS[:Y]]
    sys1 = QuantumSystem(H_drives, [1.0, 1.0])
    @test isempty(sys1.global_params)
    @test sys1.global_params isa NamedTuple

    # Test with global parameters
    global_params = (δ = 0.5, Ω = 1.0, α = -0.2)
    sys2 = QuantumSystem(H_drives, [1.0, 1.0]; global_params = global_params)
    @test sys2.global_params === global_params
    @test sys2.global_params.δ == 0.5
    @test sys2.global_params.Ω == 1.0
    @test sys2.global_params.α == -0.2

    # Test with function-based constructor
    H(u, t) = u[1] * PAULIS[:X] + u[2] * PAULIS[:Y]
    sys3 = QuantumSystem(H, [1.0, 1.0]; global_params = (β = 2.5,))
    @test sys3.global_params.β == 2.5

    # Test that function-based system can use global params via closure
    # Users should capture global_params in their H function definition
    gp = (scale = 2.0,)
    H_with_global(u, t) = gp.scale * (u[1] * PAULIS[:X] + u[2] * PAULIS[:Y])
    sys4 = QuantumSystem(H_with_global, [1.0, 1.0]; global_params = gp)
    @test sys4.global_params.scale == 2.0
    # Verify H function uses the global parameter via closure
    u_test = [0.5, 0.5]
    H_result = sys4.H(u_test, 0.0)
    H_expected = 2.0 * (0.5 * PAULIS[:X] + 0.5 * PAULIS[:Y])
    @test norm(H_result - H_expected) < 1e-10
end
