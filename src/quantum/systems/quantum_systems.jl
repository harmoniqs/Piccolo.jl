
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
- `H::Function`: The Hamiltonian function: `(u, t) -> H(u, t)`
- `G::Function`: The isomorphic generator function: `(u, t) -> G(u, t)`
- `H_drift`: The drift Hamiltonian. Typically `SparseMatrixCSC{ComplexF64,Int}` for matrix-based systems, but can be any type (e.g., an `AbstractDynamicsOperator` from Piccolissimo) for structured operator acceleration.
- `drift_terms::Vector{DriftTerm}`: Source of truth for drift terms with optional time modulation
- `H_drives::Vector{AbstractDrive}`: Drive terms (`LinearDrive`, `NonlinearDrive`, or `ModulatedDrive`)
- `drive_bounds::Vector{Tuple{Float64, Float64}}`: Drive amplitude bounds per control
- `n_drives::Int`: Number of control channels
- `levels::Int`: System dimension
- `time_dependent::Bool`: Set automatically when modulation is present
- `global_params::NamedTuple`: Global parameters
- `hermitian::Bool`: Whether the Hamiltonian is Hermitian (governs default ODE algorithm)

# Time Modulation (Pair syntax)

Use `=>` to attach time-dependent modulation to drift or drive terms:

```julia
# Modulated drive: u[1] * cos(ωt) * H_x
sys = QuantumSystem(H_z, [H_x => t -> cos(ω*t)], [1.0])

# Modulated drift: cos(ωt) * H_z
sys = QuantumSystem(H_z => t -> cos(ω*t), [H_x], [1.0])

# Multiple drift terms
sys = QuantumSystem([H_z, H_x => t -> cos(ω*t)], [H_y], [1.0])
```

See also [`OpenQuantumSystem`](@ref), [`VariationalQuantumSystem`](@ref).
"""
struct QuantumSystem{F1<:Function,F2<:Function,PT<:NamedTuple,DT,HD,DD<:AbstractVector{<:AbstractDrive}} <:
       AbstractQuantumSystem
    H::F1
    G::F2
    H_drift::HD
    drift_terms::DT
    H_drives::DD   # parametric — preserves caller's concrete eltype (e.g. Vector{LinearDrive}
                   # or Vector{Union{LinearDrive,NonlinearDrive}}), avoiding dynamic dispatch
                   # on drive_coeff / drive_coeff_jac in per-substep RHS hot loops.
    drive_bounds::Vector{Tuple{Float64,Float64}}
    n_drives::Int
    levels::Int
    time_dependent::Bool
    global_params::PT
    hermitian::Bool
end

"""
    default_algorithm(sys::AbstractQuantumSystem)

Return the default ODE algorithm for trajectory rollouts.
Uses `Tsit5()` for non-Hermitian systems (where Magnus adaptive error
control fails), `MagnusAdapt4()` for Hermitian systems.
"""
function default_algorithm end

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
    hermitian::Bool = true,
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
    if hermitian
        @assert is_hermitian(H_drift) "Drift Hamiltonian H(u=0, t=0) is not Hermitian"
    end

    # Check that Hamiltonian is Hermitian for sample control values
    u_test_controls = [b isa Tuple ? (b[1] + b[2]) / 2 : 0.0 for b in drive_bounds]
    u_test =
        n_globals > 0 ? vcat(u_test_controls, collect(values(global_params))) :
        u_test_controls
    H_test = H(u_test, 0.0)
    if hermitian
        @assert is_hermitian(H_test) "Hamiltonian H(u, t=0) is not Hermitian for test control values u=$u_test"
    end

    return QuantumSystem(
        (u, t) -> H(u, t),
        (u, t) -> Isomorphisms.G(H(u, t)),
        sparse(H_drift),
        [DriftTerm(sparse(H_drift))],
        AbstractDrive[],
        drive_bounds,
        n_drives,
        levels,
        time_dependent,
        _float_params(global_params),
        hermitian,
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
    hermitian::Bool = true,
)
    drive_bounds = [b isa Tuple ? b : (-b, b) for b in drive_bounds]

    # Check that H_drift is Hermitian
    if hermitian
        @assert is_hermitian(H_drift) "Drift Hamiltonian H_drift is not Hermitian"
    end

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

    # Build LinearDrive objects from H_drives matrices.
    # Use concrete eltype LinearDrive (not AbstractDrive) for type-stable
    # iteration in per-substep RHS loops.
    linear_drives = [LinearDrive(H_drives[i], i) for i = 1:n_drives]

    return QuantumSystem(
        H,
        G,
        H_drift,
        [DriftTerm(H_drift)],
        linear_drives,
        drive_bounds,
        n_drives,
        levels,
        time_dependent,
        _float_params(global_params),
        hermitian,
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
    hermitian::Bool = true,
) where {ℂ<:Number}
    @assert !isempty(H_drives) "At least one drive is required"
    return QuantumSystem(
        spzeros(ℂ, size(H_drives[1])),
        H_drives,
        drive_bounds;
        time_dependent = time_dependent,
        global_params = global_params,
        hermitian = hermitian,
    )
end

"""
    QuantumSystem(drives::Vector{<:AbstractDrive}, drive_bounds::Vector; time_dependent::Bool=false)

Convenience constructor for a typed-drive system with no drift Hamiltonian (H_drift = 0).

# Example
```julia
sys = QuantumSystem([LinearDrive(sparse(PAULIS[:X]), 1)], [1.0])
```
"""
function QuantumSystem(
    drives::Vector{<:AbstractDrive},
    drive_bounds::Vector{<:Union{Tuple{Float64,Float64},Float64}};
    time_dependent::Bool = false,
    global_params::NamedTuple = NamedTuple(),
    hermitian::Bool = true,
)
    @assert !isempty(drives) "At least one drive is required"
    levels = drive_dim(first(drives))
    return QuantumSystem(
        spzeros(ComplexF64, levels, levels),
        drives,
        drive_bounds;
        time_dependent = time_dependent,
        global_params = global_params,
        hermitian = hermitian,
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
    hermitian::Bool = true,
) where {ℂ<:Number}
    QuantumSystem(
        H_drift,
        Matrix{ℂ}[],
        Float64[];
        time_dependent = time_dependent,
        global_params = global_params,
        hermitian = hermitian,
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
    hermitian::Bool = true,
)
    drive_bounds = normalize_drive_bounds(drive_bounds)

    # Check that H_drift is Hermitian
    if hermitian
        @assert is_hermitian(H_drift) "Drift Hamiltonian H_drift is not Hermitian"
    end

    # Check that all drive operators are Hermitian
    for (i, d) in enumerate(drives)
        @assert is_hermitian(drive_matrix(d)) "Drive operator drives[$i].H is not Hermitian"
    end

    H_drift = sparse(ComplexF64.(H_drift))
    n_drives = length(drive_bounds)
    levels = size(H_drift, 1)

    # Validate LinearDrive indices are within the control dimension (unwrap ModulatedDrive)
    for (i, d) in enumerate(drives)
        base = d isa ModulatedDrive ? d.base : d
        if base isa LinearDrive
            @assert 1 <= base.index <= n_drives "LinearDrive at drives[$i] has index $(base.index) but control dimension is $n_drives (length of drive_bounds)"
        end
    end

    # Validate NonlinearDrive Jacobians against ForwardDiff (unwrap ModulatedDrive)
    for d in drives
        base = d isa ModulatedDrive ? d.base : d
        if base isa NonlinearDrive
            validate_drive_jacobian(base, n_drives + length(global_params))
        end
    end

    # Build H(u,t) and G(u,t) from drives
    G_drift = sparse(Isomorphisms.G(H_drift))
    H_drive_mats = [sparse(ComplexF64.(drive_matrix(d))) for d in drives]
    G_drive_mats = [sparse(Isomorphisms.G(H_d)) for H_d in H_drive_mats]

    if isempty(drives)
        H_fn = (u, t) -> H_drift
        G_fn = (u, t) -> G_drift
    else
        H_fn =
            (u, t) ->
                H_drift +
                sum(drive_coeff(d, u, t) * H_d for (d, H_d) in zip(drives, H_drive_mats))
        G_fn =
            (u, t) ->
                G_drift +
                sum(drive_coeff(d, u, t) * G_d for (d, G_d) in zip(drives, G_drive_mats))
    end

    # Type-stability fix: preserve the caller's concrete eltype on `drives`
    # instead of widening to AbstractDrive. The struct's `DD` type parameter
    # captures it. Callers can pass Vector{LinearDrive}, Vector{NonlinearDrive},
    # or Vector{Union{LinearDrive,NonlinearDrive}} for fully type-stable iteration.
    return QuantumSystem(
        H_fn,
        G_fn,
        H_drift,
        [DriftTerm(H_drift)],
        drives,
        drive_bounds,
        n_drives,
        levels,
        time_dependent,
        _float_params(global_params),
        hermitian,
    )
end

# ----------------------------------------------------------------------------- #
# Modulation normalization helpers
# ----------------------------------------------------------------------------- #

"""Normalize a drift argument into a Vector{DriftTerm}."""
_normalize_drift(H::AbstractMatrix) = [DriftTerm(sparse(ComplexF64.(H)))]
_normalize_drift(p::Pair{<:AbstractMatrix,<:Function}) =
    [DriftTerm(sparse(ComplexF64.(p.first)), p.second)]
function _normalize_drift(terms::AbstractVector)
    return [_normalize_drift_entry(t) for t in terms]
end
_normalize_drift_entry(H::AbstractMatrix) = DriftTerm(sparse(ComplexF64.(H)))
_normalize_drift_entry(p::Pair{<:AbstractMatrix,<:Function}) =
    DriftTerm(sparse(ComplexF64.(p.first)), p.second)

"""Normalize a single drive input into an AbstractDrive."""
_normalize_drive(H::AbstractMatrix, index::Int) = LinearDrive(sparse(ComplexF64.(H)), index)
function _normalize_drive(p::Pair{<:AbstractMatrix,<:Function}, index::Int)
    return ModulatedDrive(LinearDrive(sparse(ComplexF64.(p.first)), index), p.second)
end
_normalize_drive(d::AbstractDrive, ::Int) = d
_normalize_drive(p::Pair{<:AbstractDrive,<:Function}, ::Int) =
    ModulatedDrive(p.first, p.second)

"""Check if any drift_terms or drives have non-identity modulation."""
function _has_any_modulation(drift_terms, drives)
    any(has_modulation, drift_terms) || any(has_modulation, drives)
end

const _DriftInput =
    Union{AbstractMatrix{<:Number},Pair{<:AbstractMatrix{<:Number},<:Function}}
const _DriftInputs = Union{_DriftInput,AbstractVector}

"""
    QuantumSystem(drift, H_drives_input, drive_bounds; ...)

Pair-based constructor supporting modulated drift terms and modulated drive terms.

Drift may be:
- A plain matrix: `H_z`
- A modulated matrix pair: `H_z => t -> cos(ω*t)`
- A vector of the above: `[H_z, H_x => t -> cos(ω*t)]`

Each element of `H_drives_input` may be:
- A plain matrix: `H_x`
- A modulated matrix pair: `H_x => t -> cos(ω*t)`
- An `AbstractDrive` (e.g. `NonlinearDrive`)
- A modulated drive pair: `nd => t -> cos(ω*t)`

`time_dependent` is auto-detected from the presence of modulation and need not be
set manually.
"""
function QuantumSystem(
    drift::_DriftInputs,
    H_drives_input::AbstractVector,
    drive_bounds::Vector{<:Union{Tuple{Float64,Float64},Float64}};
    time_dependent::Bool = false,
    global_params::NamedTuple = NamedTuple(),
    hermitian::Bool = true,
)
    # Normalize drift into Vector{DriftTerm}
    drift_terms = _normalize_drift(drift)
    H_drift_sum =
        isempty(drift_terms) ? spzeros(ComplexF64, 0, 0) : sum(dt.H for dt in drift_terms)

    # Check drift is Hermitian
    if hermitian
        @assert is_hermitian(H_drift_sum) "Drift Hamiltonian is not Hermitian"
    end

    # Normalize drives — track linear index separately for LinearDrive(index).
    # Accumulator is AbstractDrive[] during construction (push! needs flexibility),
    # then narrowed to the most-specific Union eltype via `identity.()` before
    # passing to QuantumSystem. This preserves type stability in the per-substep
    # RHS hot loops.
    drives = AbstractDrive[]
    linear_idx = 1
    for d_input in H_drives_input
        if d_input isa Pair{<:AbstractDrive,<:Function}
            push!(drives, ModulatedDrive(d_input.first, d_input.second))
        elseif d_input isa AbstractDrive
            push!(drives, d_input)
        elseif d_input isa Pair{<:AbstractMatrix,<:Function}
            push!(
                drives,
                ModulatedDrive(
                    LinearDrive(sparse(ComplexF64.(d_input.first)), linear_idx),
                    d_input.second,
                ),
            )
            linear_idx += 1
        else
            # Plain matrix
            push!(drives, LinearDrive(sparse(ComplexF64.(d_input)), linear_idx))
            linear_idx += 1
        end
    end
    # Narrow eltype to the most specific Union after dynamic construction.
    drives = identity.(drives)

    # Check drive operators are Hermitian
    if hermitian
        for (i, d) in enumerate(drives)
            @assert is_hermitian(drive_matrix(d)) "Drive operator drives[$i].H is not Hermitian"
        end
    end

    # Detect time dependence
    td = time_dependent || _has_any_modulation(drift_terms, drives)

    n_drives = length(drive_bounds)
    drive_bounds_vec = normalize_drive_bounds(drive_bounds)

    H_drive_mats = [drive_matrix(d) for d in drives]
    G_drift_terms = [sparse(Isomorphisms.G(dt.H)) for dt in drift_terms]
    G_drive_mats = [sparse(Isomorphisms.G(d)) for d in drives]

    levels = size(H_drift_sum, 1)

    if isempty(drives)
        H_fn = (u, t) -> sum(dt.modulation(t) * dt.H for dt in drift_terms)
        G_fn =
            (u, t) -> sum(
                dt.modulation(t) * G_dt for (dt, G_dt) in zip(drift_terms, G_drift_terms)
            )
    else
        H_fn =
            (u, t) ->
                sum(dt.modulation(t) * dt.H for dt in drift_terms) +
                sum(drive_coeff(d, u, t) * H_d for (d, H_d) in zip(drives, H_drive_mats))
        G_fn =
            (u, t) ->
                sum(
                    dt.modulation(t) * G_dt for
                    (dt, G_dt) in zip(drift_terms, G_drift_terms)
                ) +
                sum(drive_coeff(d, u, t) * G_d for (d, G_d) in zip(drives, G_drive_mats))
    end

    return QuantumSystem(
        H_fn,
        G_fn,
        H_drift_sum,
        drift_terms,
        collect(AbstractDrive, drives),
        drive_bounds_vec,
        n_drives,
        levels,
        td,
        _float_params(global_params),
        hermitian,
    )
end

# Disambiguator: a `(drift::_DriftInputs, drives::Vector{<:AbstractDrive}, ...)` call
# matches both the pair-based method above and the structured-drift method below.
# Both methods share the same intent for this input shape; route through the
# pair-based path (it handles `Vector{<:AbstractDrive}` natively in its drive loop)
# so dispatch is unambiguous.
function QuantumSystem(
    drift::_DriftInputs,
    drives::Vector{<:AbstractDrive},
    drive_bounds::Vector{<:Union{Tuple{Float64,Float64},Float64}};
    time_dependent::Bool = false,
    global_params::NamedTuple = NamedTuple(),
    hermitian::Bool = true,
)
    return invoke(
        QuantumSystem,
        Tuple{_DriftInputs,AbstractVector,Vector{<:Union{Tuple{Float64,Float64},Float64}}},
        drift,
        drives,
        drive_bounds;
        time_dependent = time_dependent,
        global_params = global_params,
        hermitian = hermitian,
    )
end

"""
    QuantumSystem(
        H_drift,
        drives::Vector{<:AbstractDrive},
        drive_bounds::Vector;
        time_dependent::Bool=false,
        global_params::NamedTuple=NamedTuple()
    )

Construct a QuantumSystem from a non-matrix drift Hamiltonian (e.g., a structured
operator from Piccolissimo) and typed drive terms.

The drift is stored directly without sparse conversion, enabling structured operator
acceleration in the spline integrator. The `H(u,t)` and `G(u,t)` closures are built
from materialized matrices for backward compatibility with function-based paths.

# Example
```julia
using Piccolissimo: DiagonalOperator
H_drift_op = DiagonalOperator(ComplexF64[0.0, 1.0, 2.0])
sys = QuantumSystem(H_drift_op, [LinearDrive(H_x_op, 1)], [1.0])
```
"""
function QuantumSystem(
    H_drift_op,
    drives::Vector{<:AbstractDrive},
    drive_bounds::Vector{<:Union{Tuple{Float64,Float64},Float64}};
    time_dependent::Bool = false,
    global_params::NamedTuple = NamedTuple(),
    hermitian::Bool = true,
)
    # Only dispatch here for non-AbstractMatrix types; matrices use the method above
    H_drift_op isa AbstractMatrix && return QuantumSystem(
        convert(AbstractMatrix{ComplexF64}, H_drift_op),
        drives,
        drive_bounds;
        time_dependent = time_dependent,
        global_params = global_params,
        hermitian = hermitian,
    )

    drive_bounds = normalize_drive_bounds(drive_bounds)

    # Check that all drive operators are Hermitian (drive_matrix materializes)
    for (i, d) in enumerate(drives)
        @assert is_hermitian(drive_matrix(d)) "Drive operator drives[$i].H is not Hermitian"
    end

    n_drives = length(drive_bounds)
    levels = size(H_drift_op, 1)

    # Validate LinearDrive indices
    for (i, d) in enumerate(drives)
        if d isa LinearDrive
            @assert 1 <= d.index <= n_drives "LinearDrive at drives[$i] has index $(d.index) but control dimension is $n_drives"
        end
    end

    # Validate NonlinearDrive Jacobians
    for d in drives
        if d isa NonlinearDrive
            validate_drive_jacobian(d, n_drives + length(global_params))
        end
    end

    # Build H(u,t) and G(u,t) from materialized matrices (for backward-compat closures)
    H_drift_mat = sparse(ComplexF64.(_ensure_matrix(H_drift_op)))
    G_drift = sparse(Isomorphisms.G(H_drift_mat))
    H_drive_mats = [sparse(ComplexF64.(drive_matrix(d))) for d in drives]
    G_drive_mats = [sparse(Isomorphisms.G(H_d)) for H_d in H_drive_mats]

    if isempty(drives)
        H_fn = (u, t) -> H_drift_mat
        G_fn = (u, t) -> G_drift
    else
        H_fn =
            (u, t) ->
                H_drift_mat +
                sum(drive_coeff(d, u, t) * H_d for (d, H_d) in zip(drives, H_drive_mats))
        G_fn =
            (u, t) ->
                G_drift +
                sum(drive_coeff(d, u, t) * G_d for (d, G_d) in zip(drives, G_drive_mats))
    end

    return QuantumSystem(
        H_fn,
        G_fn,
        H_drift_op,    # store the structured operator, not the sparse matrix
        [DriftTerm(H_drift_mat)],
        collect(AbstractDrive, drives),
        drive_bounds,
        n_drives,
        levels,
        time_dependent,
        _float_params(global_params),
        hermitian,
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

    # Non-Hermitian drift with hermitian=false should succeed
    H_drift_nh = ComplexF64[0 1; 1 0] + ComplexF64[-0.5im 0; 0 0]
    H_drives_nh = [ComplexF64[0 1; 1 0]]
    bounds_nh = [(-1.0, 1.0)]
    sys_nh = QuantumSystem(H_drift_nh, H_drives_nh, bounds_nh; hermitian = false)
    @test sys_nh isa QuantumSystem
    @test sys_nh.n_drives == 1
    @test sys_nh.levels == 2
    @test sys_nh.hermitian == false

    # Hermitian system should have hermitian=true by default
    sys_h = QuantumSystem(ComplexF64[1 0; 0 -1], [ComplexF64[0 1; 1 0]], [(-1.0, 1.0)])
    @test sys_h.hermitian == true

    # Non-Hermitian drift with hermitian=true (default) should still fail
    @test_throws AssertionError QuantumSystem(H_drift_nh, H_drives_nh, bounds_nh)
    @test_throws AssertionError QuantumSystem(
        H_drift_nh,
        H_drives_nh,
        bounds_nh;
        hermitian = true,
    )
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

    # H_drives field should be auto-populated with LinearDrives
    @test length(sys.H_drives) == 2
    @test all(d -> d isa LinearDrive, sys.H_drives)
    @test sys.H_drives[1].index == 1
    @test sys.H_drives[2].index == 2

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
    @test isempty(sys3.H_drives)
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
    @test length(sys.H_drives) == 2

    u_test = [0.3, 0.7]
    H_expected = PAULIS.Z + 0.3 * PAULIS.X + 0.7 * PAULIS.Y
    @test norm(sys.H(u_test, 0.0) - H_expected) < 1e-10

    # Mixed drives (linear + nonlinear)
    nonlinear_drive = NonlinearDrive(
        sparse(ComplexF64.([1.0 0.0; 0.0 -1.0])),
        u -> u[1]^2 + u[2]^2,
        (u, j) -> j == 1 ? 2u[1] : j == 2 ? 2u[2] : 0.0,
    )
    mixed_drives = [drives..., nonlinear_drive]
    sys2 = QuantumSystem(H_drift, mixed_drives, [1.0, 1.0])

    @test sys2 isa QuantumSystem
    @test sys2.n_drives == 2  # control dimension = 2
    @test length(sys2.H_drives) == 3  # 3 drive terms

    # H function should work correctly
    H_result = sys2.H(u_test, 0.0)
    coeff_nonlinear = 0.3^2 + 0.7^2  # = 0.58
    H_expected2 =
        PAULIS.Z + 0.3 * PAULIS.X + 0.7 * PAULIS.Y + coeff_nonlinear * [1.0 0.0; 0.0 -1.0]
    @test norm(H_result - H_expected2) < 1e-10

    # has_nonlinear_drives
    @test !has_nonlinear_drives(sys.H_drives)
    @test has_nonlinear_drives(sys2.H_drives)
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

@testitem "QuantumSystem: NonlinearDrive reading appended global params" begin
    # Regression: validator must accept drives whose coefficient reads into
    # the appended global-parameter slots. The typed-drives constructor passes
    # n_drives + length(global_params) into validate_drive_jacobian.
    using Piccolo
    using SparseArrays

    H = sparse([0.0+0im 1.0+0im; 1.0+0im 0.0+0im])
    # Coefficient reads u[3] = the single appended global parameter.
    d = NonlinearDrive(H, u -> u[1] * u[2] * u[3])

    sys = QuantumSystem(H, [d], [(0.0, 1.0), (0.0, 1.0)]; global_params = (g1 = 0.5,))
    @test sys.n_drives == 2
    @test length(sys.global_params) == 1
end

@testitem "QuantumSystem Pair-based modulation constructors" begin
    using Piccolo
    using SparseArrays

    H_z = sparse(ComplexF64[1 0; 0 -1])
    H_x = sparse(ComplexF64[0 1; 1 0])
    H_y = sparse(ComplexF64[0 -im; im 0])
    omega = 2.0

    # Modulated drive: H_x => t -> cos(omega*t)
    sys1 = QuantumSystem(H_z, [H_x => t -> cos(omega * t), H_y], [1.0, 1.0])
    @test sys1.n_drives == 2
    @test sys1.H_drives[1] isa ModulatedDrive
    @test sys1.H_drives[2] isa LinearDrive
    @test length(sys1.drift_terms) == 1
    @test !has_modulation(sys1.drift_terms[1])

    # Modulated drift: H_z => t -> cos(omega*t)
    sys2 = QuantumSystem(H_z => t -> cos(omega * t), [H_x, H_y], [1.0, 1.0])
    @test sys2.n_drives == 2
    @test length(sys2.drift_terms) == 1
    @test has_modulation(sys2.drift_terms[1])
    @test sys2.H_drift ≈ H_z  # static sum

    # Multiple drift terms
    sys3 = QuantumSystem([H_z, H_x => t -> cos(omega * t)], [H_y], [1.0])
    @test length(sys3.drift_terms) == 2
    @test !has_modulation(sys3.drift_terms[1])
    @test has_modulation(sys3.drift_terms[2])
    @test sys3.H_drift ≈ H_z + H_x  # sum of all drift matrices

    # Typed drive with modulation: NonlinearDrive => modulation
    nd = NonlinearDrive(H_x, u -> u[1]^2; active_controls = [1])
    sys4 = QuantumSystem(H_z, [nd => t -> cos(omega * t)], [1.0])
    @test sys4.H_drives[1] isa ModulatedDrive
    @test sys4.H_drives[1].base === nd

    # time_dependent auto-detection
    @test !QuantumSystem(H_z, [H_x], [1.0]).time_dependent
    @test sys1.time_dependent
    @test sys2.time_dependent
    @test sys3.time_dependent

    # H(u, t) closure evaluates modulation
    u = [1.0, 0.0]
    @test sys2.H(u, 0.0) ≈ H_z * cos(0.0) + 1.0 * H_x
    t_test = 0.3
    @test sys2.H(u, t_test) ≈ H_z * cos(omega * t_test) + 1.0 * H_x
end
