module QuantumSystems

export AbstractQuantumSystem
export QuantumSystem
export OpenQuantumSystem
export VariationalQuantumSystem
export CompositeQuantumSystem

export AbstractDrive, LinearDrive, NonlinearDrive
export drive_coeff, drive_coeff_jac
export active_controls
export has_nonlinear_drives
export validate_drive_jacobian

export get_drift
export get_drives
export get_drive_terms
export get_c_ops
export compact_lindbladian_generators

using ..Isomorphisms
using ..QuantumObjectUtils
using ..LiftedOperators

using LinearAlgebra
using SparseArrays
using TestItems
using ForwardDiff

# ----------------------------------------------------------------------------- #
# Drive Bounds Types and Utilities
# ----------------------------------------------------------------------------- #

"""
    DriveBounds

Type alias for drive amplitude bounds input. Bounds can be specified as:
- A tuple `(lower, upper)` for asymmetric bounds
- A scalar value which is interpreted as symmetric bounds `(-value, value)`

# Examples
```julia
drive_bounds = [(-1.0, 1.0), 0.5, (-0.3, 0.7)]
# Interpreted as: [(-1.0, 1.0), (-0.5, 0.5), (-0.3, 0.7)]
```
"""
const DriveBounds = Vector{<:Union{Tuple{Float64,Float64},Float64}}

"""
    normalize_drive_bounds(bounds::DriveBounds)

Convert drive bounds to a consistent tuple format. Scalar values are converted to 
symmetric bounds around zero: `b` becomes `(-b, b)`.

# Arguments
- `bounds::DriveBounds`: Input bounds, can be tuples or scalars

# Returns
- `Vector{Tuple{Float64, Float64}}`: Normalized bounds as tuples

# Examples
```julia
# All scalars (symmetric bounds)
normalize_drive_bounds([1.0, 1.5, 0.5])
# Returns: [(-1.0, 1.0), (-1.5, 1.5), (-0.5, 0.5)]

# All tuples (asymmetric bounds)
normalize_drive_bounds([(-2.0, 3.0), (-1.0, 1.0)])
# Returns: [(-2.0, 3.0), (-1.0, 1.0)]

# Mixed types (requires explicit type annotation)
normalize_drive_bounds(Union{Float64, Tuple{Float64,Float64}}[1.0, (-2.0, 3.0)])
# Returns: [(-1.0, 1.0), (-2.0, 3.0)]
```
"""
function normalize_drive_bounds(bounds::DriveBounds)
    return [b isa Tuple ? b : (-b, b) for b in bounds]
end

# ----------------------------------------------------------------------------- #
# AbstractQuantumSystem
# ----------------------------------------------------------------------------- #

"""
    AbstractQuantumSystem

Abstract type for defining systems.
"""
abstract type AbstractQuantumSystem end

# ----------------------------------------------------------------------------- #
# AbstractQuantumSystem methods
# ----------------------------------------------------------------------------- #

"""
    get_drift(sys::AbstractQuantumSystem)

Returns the drift Hamiltonian of the system.
"""
get_drift(sys::AbstractQuantumSystem) = sys.H(zeros(sys.n_drives), 0.0)

"""
    get_drives(sys::AbstractQuantumSystem)

Returns the drive Hamiltonian matrices of the system.

For systems with typed drives (`sys.drives`), returns the `H` matrix from each
drive term. For function-based systems without drives, extracts operators via
basis vector evaluation.

!!! note
    When nonlinear drives are present, the number of returned matrices may differ
    from `sys.n_drives` (the control dimension). For example, a system with 2
    controls and 3 drive terms (2 linear + 1 nonlinear) returns 3 matrices.
    Use [`get_drive_terms`](@ref) to access the full `AbstractDrive` objects with
    coefficient functions, Jacobians, and active control indices.
"""
function get_drives(sys::AbstractQuantumSystem)
    if hasproperty(sys, :drives) && !isempty(sys.drives)
        return [d.H for d in sys.drives]
    end
    H_drift = get_drift(sys)
    # Basis vectors for controls will extract drive operators (linear systems only)
    return [sys.H(I[1:sys.n_drives, i], 0.0) - H_drift for i ∈ 1:sys.n_drives]
end

"""
    get_drive_terms(sys::AbstractQuantumSystem) -> Vector{AbstractDrive}

Return the typed drive terms from the system. Each `AbstractDrive` pairs a
Hamiltonian operator with a scalar coefficient function and Jacobian.

Returns `sys.drives` directly when available, or an empty vector for
function-based systems that don't use typed drives.

See also [`get_drives`](@ref) for just the Hamiltonian matrices.
"""
function get_drive_terms(sys::AbstractQuantumSystem)
    if hasproperty(sys, :drives)
        return sys.drives
    end
    return AbstractDrive[]
end

function Base.show(io::IO, sys::AbstractQuantumSystem)
    print(io, "$(nameof(typeof(sys))): levels = $(sys.levels), n_drives = $(sys.n_drives)")
end

# ----------------------------------------------------------------------------- #
# Quantum Toolbox ext
# ----------------------------------------------------------------------------- #

function get_c_ops end

# ----------------------------------------------------------------------------- #
# Global parameter utilities
# ----------------------------------------------------------------------------- #

"""
    _float_params(nt::NamedTuple)

Convert all values in a NamedTuple to their floating-point equivalents.
Ensures type-stable ODE solutions when global parameters are updated during optimization.
"""
_float_params(nt::NamedTuple{K}) where {K} = NamedTuple{K}(float.(values(nt)))

# ----------------------------------------------------------------------------- #
# Quantum System Types
# ----------------------------------------------------------------------------- #

include("drives.jl")
include("quantum_systems.jl")
include("open_quantum_systems.jl")
include("variational_quantum_systems.jl")
include("composite_quantum_systems.jl")

end
