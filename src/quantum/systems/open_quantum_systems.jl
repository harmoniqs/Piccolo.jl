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
- `H_drives::Vector{AbstractDrive}`: The canonical drive terms, each pairing an operator with a coefficient function. Matrix-based constructors auto-populate this with `LinearDrive` objects; function-based systems leave it empty.
- `drive_bounds::Vector{Tuple{Float64, Float64}}`: Drive amplitude bounds
- `n_drives::Int`: The number of control drives
- `levels::Int`: The number of levels in the system
- `dissipators::Vector{AbstractDissipator}`: The typed dissipator terms, each pairing a jump operator with a scalar rate coefficient function. Matrix-based constructors auto-wrap as `LinearDissipator` with unit rate.
- `time_dependent::Bool`: Whether the Hamiltonian has explicit time dependence
- `global_params::NamedTuple`: Global parameters stored with the system (e.g., physical constants)

!!! note
    `dissipation_operators` is available as a derived read-only property via
    `Base.getproperty` shim (materializing a `Vector{SparseMatrixCSC{ComplexF64,Int}}`
    from `dissipators` on demand), not a stored field. It exists for backward
    compatibility; new code should consume `sys.dissipators` directly.

See also [`QuantumSystem`](@ref Piccolo.Quantum.QuantumSystems.QuantumSystem).
"""
struct OpenQuantumSystem{F1<:Function,F2<:Function,PT<:NamedTuple} <: AbstractQuantumSystem
    H::F1
    𝒢::F2
    H_drift::SparseMatrixCSC{ComplexF64,Int}
    H_drives::Vector{AbstractDrive}
    drive_bounds::Vector{Tuple{Float64,Float64}}
    n_drives::Int
    levels::Int
    dissipators::Vector{AbstractDissipator}
    time_dependent::Bool
    global_params::PT
end

# Backward-compatible property access: `sys.dissipation_operators` returns a
# fresh Vector{SparseMatrixCSC{ComplexF64,Int}} materialized from the typed
# dissipators. Read-only — mutations do not propagate back.
function Base.getproperty(sys::OpenQuantumSystem, s::Symbol)
    if s === :dissipation_operators
        return SparseMatrixCSC{ComplexF64,Int}[
            sparse(ComplexF64.(dissipator_matrix(d))) for d in getfield(sys, :dissipators)
        ]
    end
    return getfield(sys, s)
end

# Propertynames: expose both the real field and the derived view for discoverability.
# Overload the 2-arg form so that `propertynames(sys, true)` and
# `propertynames(sys; private=true)` route through and still surface the
# `:dissipation_operators` shim instead of falling back to Base's default.
Base.propertynames(::OpenQuantumSystem, private::Bool = false) = (
    :H,
    :𝒢,
    :H_drift,
    :H_drives,
    :drive_bounds,
    :n_drives,
    :levels,
    :dissipators,
    :dissipation_operators,
    :time_dependent,
    :global_params,
)

"""
    _normalize_dissipators(arg) -> Vector{AbstractDissipator}

Coerce a dissipator-list argument (matrices OR typed dissipators) to a
canonical `Vector{AbstractDissipator}`. Matrix inputs are auto-wrapped as
`LinearDissipator` with unit rate.
"""
function _normalize_dissipators(arg::AbstractVector)
    out = AbstractDissipator[]
    for x in arg
        if x isa AbstractDissipator
            push!(out, x)
        elseif x isa AbstractMatrix
            push!(out, LinearDissipator(sparse(ComplexF64.(x))))
        else
            error(
                "_normalize_dissipators: expected AbstractDissipator or AbstractMatrix, got $(typeof(x))",
            )
        end
    end
    return out
end
_normalize_dissipators(::Nothing) = AbstractDissipator[]

"""
    _select_dissipators(dissipation_operators, dissipators) -> Vector{AbstractDissipator}

Internal helper: enforce "pass at most one of `dissipation_operators=` / `dissipators=`",
normalize the chosen input to `Vector{AbstractDissipator}`.
"""
function _select_dissipators(dissipation_operators, dissipators)
    if !isnothing(dissipation_operators) && !isnothing(dissipators)
        error("Pass either `dissipation_operators=` or `dissipators=`, not both.")
    end
    if !isnothing(dissipators)
        return _normalize_dissipators(dissipators)
    elseif !isnothing(dissipation_operators)
        return _normalize_dissipators(dissipation_operators)
    else
        return AbstractDissipator[]
    end
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
    dissipation_operators = nothing,
    dissipators = nothing,
    time_dependent::Bool = false,
    global_params::NamedTuple = NamedTuple(),
)
    drive_bounds = normalize_drive_bounds(drive_bounds)

    diss_list = _select_dissipators(dissipation_operators, dissipators)
    diss_mats = SparseMatrixCSC{ComplexF64,Int}[
        sparse(ComplexF64.(dissipator_matrix(d))) for d in diss_list
    ]

    H_drift_sparse = sparse(H_drift)
    𝒢_drift = Isomorphisms.G(Isomorphisms.ad_vec(H_drift_sparse))

    n_drives = length(H_drives)
    H_drives_sparse = sparse.(H_drives)
    𝒢_drives = [Isomorphisms.G(Isomorphisms.ad_vec(H_drive)) for H_drive in H_drives_sparse]

    levels = size(H_drift, 1)

    # Build LinearDrive objects from H_drives matrices
    linear_drives = AbstractDrive[LinearDrive(H_drives_sparse[i], i) for i = 1:n_drives]

    # Precompute iso_D factors for each dissipator (matched to diss_list)
    diss_iso_mats = [Isomorphisms.iso_D(L) for L in diss_mats]

    # Build H(u,t) and 𝒢(u) — assembled at call time from drive_coeff + rate_coeff
    if n_drives == 0
        H = (u, t) -> H_drift_sparse
    else
        H = (u, t) -> begin
            out = H_drift_sparse
            for (d, H_d) in zip(linear_drives, H_drives_sparse)
                out = out + drive_coeff(d, u, t) * H_d
            end
            return out
        end
    end

    𝒢 =
        let 𝒢d = 𝒢_drift,
            drs = linear_drives,
            𝒢drs = 𝒢_drives,
            dss = diss_list,
            𝒟s = diss_iso_mats

            u -> begin
                out = 𝒢d + zero(𝒢d)
                for (i, d) in enumerate(drs)
                    c = drive_coeff(d, u, 0.0)
                    out = out + c * 𝒢drs[i]
                end
                for (j, diss) in enumerate(dss)
                    γ = rate_coeff(diss, u)
                    out = out + γ * 𝒟s[j]
                end
                return out
            end
        end

    return OpenQuantumSystem(
        H,
        𝒢,
        H_drift_sparse,
        linear_drives,
        drive_bounds,
        n_drives,
        levels,
        diss_list,
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
    dissipation_operators = nothing,
    dissipators = nothing,
    time_dependent::Bool = false,
    global_params::NamedTuple = NamedTuple(),
)
    drive_bounds = normalize_drive_bounds(drive_bounds)

    diss_list = _select_dissipators(dissipation_operators, dissipators)
    diss_mats = SparseMatrixCSC{ComplexF64,Int}[
        sparse(ComplexF64.(dissipator_matrix(d))) for d in diss_list
    ]

    # Check that H_drift is Hermitian
    @assert is_hermitian(H_drift) "Drift Hamiltonian H_drift is not Hermitian"

    # Check that all drive operators are Hermitian
    for (i, d) in enumerate(drives)
        @assert is_hermitian(drive_matrix(d)) "Drive operator drives[$i].H is not Hermitian"
    end

    H_drift_sparse = sparse(ComplexF64.(H_drift))
    n_drives = length(drive_bounds)
    levels = size(H_drift, 1)

    # Validate LinearDrive indices (unwrap ModulatedDrive)
    for (i, d) in enumerate(drives)
        base = d isa ModulatedDrive ? d.base : d
        if base isa LinearDrive
            @assert 1 <= base.index <= n_drives "LinearDrive at drives[$i] has index $(base.index) but control dimension is $n_drives (length of drive_bounds)"
        end
    end

    # Validate NonlinearDrive Jacobians (unwrap ModulatedDrive)
    for d in drives
        base = d isa ModulatedDrive ? d.base : d
        if base isa NonlinearDrive
            validate_drive_jacobian(base, n_drives + length(global_params))
        end
    end

    # Build H(u,t) and 𝒢(u) from drives + dissipators. Use drive_matrix(d)
    # so ModulatedDrive / NonlinearDrive unwrap to their underlying operator.
    H_drive_mats = [sparse(ComplexF64.(drive_matrix(d))) for d in drives]
    𝒢_drift_ham = Isomorphisms.G(Isomorphisms.ad_vec(H_drift_sparse))
    𝒢_drive_mats = [Isomorphisms.G(Isomorphisms.ad_vec(H_d)) for H_d in H_drive_mats]
    diss_iso_mats = [Isomorphisms.iso_D(L) for L in diss_mats]

    if isempty(drives)
        H_fn = (u, t) -> H_drift_sparse
    else
        H_fn = (u, t) -> begin
            out = H_drift_sparse
            for (d, H_d) in zip(drives, H_drive_mats)
                out = out + drive_coeff(d, u, t) * H_d
            end
            return out
        end
    end

    𝒢_fn =
        let 𝒢d = 𝒢_drift_ham,
            drs = drives,
            𝒢drs = 𝒢_drive_mats,
            dss = diss_list,
            𝒟s = diss_iso_mats

            u -> begin
                out = 𝒢d + zero(𝒢d)
                for (i, d) in enumerate(drs)
                    c = drive_coeff(d, u, 0.0)
                    out = out + c * 𝒢drs[i]
                end
                for (j, diss) in enumerate(dss)
                    γ = rate_coeff(diss, u)
                    out = out + γ * 𝒟s[j]
                end
                return out
            end
        end

    return OpenQuantumSystem(
        H_fn,
        𝒢_fn,
        H_drift_sparse,
        collect(AbstractDrive, drives),
        drive_bounds,
        n_drives,
        levels,
        diss_list,
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
    dissipation_operators = nothing,
    dissipators = nothing,
    time_dependent::Bool = false,
    global_params::NamedTuple = NamedTuple(),
) where {ℂ<:Number}
    @assert !isempty(H_drives) "At least one drive is required"
    # Build the zero drift as a 2-D sparse matrix; splatting `size(...)` (a Tuple)
    # without binding the rank lets `spzeros` resolve to its `SparseVector`
    # method in JET's union split, hence the explicit (m, n) form.
    m, n = size(H_drives[1])
    return OpenQuantumSystem(
        spzeros(ℂ, m, n),
        H_drives,
        drive_bounds;
        dissipation_operators = dissipation_operators,
        dissipators = dissipators,
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
    dissipation_operators = nothing,
    dissipators = nothing,
    time_dependent::Bool = false,
    global_params::NamedTuple = NamedTuple(),
) where {ℂ<:Number}
    return OpenQuantumSystem(
        H_drift,
        Matrix{ℂ}[],
        Float64[];
        dissipation_operators = dissipation_operators,
        dissipators = dissipators,
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
    dissipation_operators = nothing,
    dissipators = nothing,
    time_dependent::Bool = false,
    global_params::NamedTuple = NamedTuple(),
) where {F<:Function}

    drive_bounds = normalize_drive_bounds(drive_bounds)

    diss_list = _select_dissipators(dissipation_operators, dissipators)
    diss_mats = SparseMatrixCSC{ComplexF64,Int}[
        sparse(ComplexF64.(dissipator_matrix(d))) for d in diss_list
    ]

    n_drives = length(drive_bounds)
    n_globals = length(global_params)

    # Build test vector u = [controls..., globals...] — matches QuantumSystem(H::Function, ...)
    u_zeros =
        n_globals > 0 ? vcat(zeros(n_drives), collect(values(global_params))) :
        zeros(n_drives)

    # Extract drift by evaluating with zero controls (and initial globals if present)
    H_drift = H(u_zeros, 0.0)
    levels = size(H_drift, 1)

    diss_iso_mats = [Isomorphisms.iso_D(L) for L in diss_mats]

    𝒢_fn = let 𝒟s = diss_iso_mats, dss = diss_list, H_user = H
        u -> begin
            out = Isomorphisms.G(Isomorphisms.ad_vec(sparse(H_user(u, 0.0))))
            for (j, diss) in enumerate(dss)
                γ = rate_coeff(diss, u)
                out = out + γ * 𝒟s[j]
            end
            return out
        end
    end

    return OpenQuantumSystem(
        H,
        𝒢_fn,
        sparse(H_drift),
        AbstractDrive[],
        drive_bounds,
        n_drives,
        levels,
        diss_list,
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
    dissipation_operators = nothing,
    dissipators = nothing,
)
    if !isempty(system.H_drives)
        # Use drives-based constructor (handles both linear and nonlinear)
        return OpenQuantumSystem(
            system.H_drift,
            system.H_drives,
            system.drive_bounds;
            dissipation_operators = dissipation_operators,
            dissipators = dissipators,
            time_dependent = system.time_dependent,
            global_params = system.global_params,
        )
    else
        # Function-based system
        return OpenQuantumSystem(
            system.H,
            system.drive_bounds;
            dissipation_operators = dissipation_operators,
            dissipators = dissipators,
            time_dependent = system.time_dependent,
            global_params = system.global_params,
        )
    end
end

# ----------------------------------------------------------------------------- #
# Compact Lindbladian generators
# ----------------------------------------------------------------------------- #

"""
    compact_lindbladian_parts(sys::OpenQuantumSystem)

Compute compact-iso matrix factors for Lindbladian assembly. Returns

    (𝒢c_drift_ham, 𝒢c_drives, 𝒢c_dissipators)

where each is projected into the compact density iso via `P · M · L`:

- `𝒢c_drift_ham` — Hamiltonian drift (no dissipator contribution baked in).
- `𝒢c_drives[i]` — generator for drive term `i` in `sys.H_drives`. Weight by
  `drive_coeff(sys.H_drives[i], u)`.
- `𝒢c_dissipators[j]` — generator for dissipator `j` in `sys.dissipators`.
  Weight by `rate_coeff(sys.dissipators[j], u)`.

Consumers assemble the full generator:

    G(u) = 𝒢c_drift_ham + Σ_i drive_coeff(d_i, u)·𝒢c_drives[i]
                        + Σ_j rate_coeff(diss_j, u)·𝒢c_dissipators[j]

See [`compact_generator_closure`](@ref) for a helper that builds this closure.

This is the canonical 3-tuple API. The legacy 2-tuple
[`compact_lindbladian_generators`](@ref) is a thin wrapper that fuses dissipators
back into the drift at `u=0` for back-compat with pre-1.12 consumers.
"""
function compact_lindbladian_parts(sys::OpenQuantumSystem)
    n = sys.levels
    P = Isomorphisms.density_projection_matrix(n)
    L = Isomorphisms.density_lift_matrix(n)

    𝒢_drift_ham = Isomorphisms.G(Isomorphisms.ad_vec(sys.H_drift))
    𝒢c_drift_ham = P * 𝒢_drift_ham * L

    # Eltype is left to Julia's inference (concrete `SparseMatrixCSC{ComplexF64,Int}`
    # in practice). When the system has no drives/dissipators the comprehension
    # produces `Vector{Any}` — that's fine for the helper's `AbstractVector`
    # signature, since each element is indexed and used as a matrix.
    𝒢c_drives =
        [P * Isomorphisms.G(Isomorphisms.ad_vec(drive_matrix(d))) * L for d in sys.H_drives]

    𝒢c_dissipators = [
        P * Isomorphisms.iso_D(sparse(ComplexF64.(dissipator_matrix(diss)))) * L for
        diss in getfield(sys, :dissipators)
    ]

    return 𝒢c_drift_ham, 𝒢c_drives, 𝒢c_dissipators
end

"""
    compact_lindbladian_generators(sys::OpenQuantumSystem)

Legacy 2-tuple API: returns `(𝒢c_drift, 𝒢c_drives)` where `𝒢c_drift` fuses the
Hamiltonian drift with the dissipator contribution evaluated at `u = zeros(n_drives)`.

Preserved as a back-compat wrapper for pre-1.12 consumers; new code should call
[`compact_lindbladian_parts`](@ref) and weight per-dissipator factors via
[`rate_coeff`](@ref). Exact for `LinearDissipator`s; for systems with
`NonlinearDissipator`s the wrapper returns the rate at `u=0`, which generally
differs from the runtime rate — those callers must migrate to the 3-tuple API.
"""
function compact_lindbladian_generators(sys::OpenQuantumSystem)
    drift_ham, drives, dissipators = compact_lindbladian_parts(sys)
    sys_dissipators = getfield(sys, :dissipators)
    if isempty(sys_dissipators)
        return drift_ham, drives
    end
    u0 = zeros(sys.n_drives)
    𝒟_at_u0 = sum(
        rate_coeff(d, u0) * M for (d, M) in zip(sys_dissipators, dissipators);
        init = zero(drift_ham),
    )
    return drift_ham + 𝒟_at_u0, drives
end

"""
    compact_generator_closure(sys, 𝒢c_drift_ham, 𝒢c_drives, 𝒢c_dissipators) -> Function

Return a closure `u_ -> G(u_)` assembling the compact Lindbladian generator as

    G(u_) = 𝒢c_drift_ham
          + Σ_i drive_coeff(sys.H_drives[i], u_)·𝒢c_drives[i]
          + Σ_j rate_coeff(sys.dissipators[j], u_)·𝒢c_dissipators[j]

Used by `BilinearIntegrator{DensityTrajectory}` (Piccolo) and
`NonHermitianExponentialIntegrator{DensityTrajectory}` + its multi-density / spline
cousins (Piccolissimo). Typed drives / dissipators whose coefficients read
global-parameter slots in the extended control vector propagate automatically.

Pre-materialized `Matrix` inputs keep the closure body allocation-light
(single-alloc per call from the additive rebuild — small relative to `expv`).
"""
function compact_generator_closure(
    sys::OpenQuantumSystem,
    𝒢c_drift_ham::AbstractMatrix,
    𝒢c_drives::AbstractVector,
    𝒢c_dissipators::AbstractVector,
)
    drives = sys.H_drives
    dissipators = getfield(sys, :dissipators)
    return let 𝒢d = 𝒢c_drift_ham,
        drs = drives,
        𝒢drs = 𝒢c_drives,
        dss = dissipators,
        𝒟s = 𝒢c_dissipators

        u_ -> begin
            out = copy(𝒢d)
            @inbounds for i in eachindex(drs)
                c = drive_coeff(drs[i], u_, 0.0)
                c == 0 && continue
                out = out + c * 𝒢drs[i]
            end
            @inbounds for j in eachindex(dss)
                γ = rate_coeff(dss[j], u_)
                γ == 0 && continue
                out = out + γ * 𝒟s[j]
            end
            return out
        end
    end
end

# ******************************************************************************* #

@testitem "compact_lindbladian_parts: 3-tuple (drift_ham, drives, dissipators)" begin
    using Piccolo
    using SparseArrays
    using LinearAlgebra

    γ = 0.3
    diss = LinearDissipator(sqrt(γ) * PAULIS.Z)
    sys = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipators = [diss])

    out = compact_lindbladian_parts(sys)
    @test length(out) == 3
    𝒢c_drift_ham, 𝒢c_drives, 𝒢c_dissipators = out
    @test length(𝒢c_drives) == 1        # one drive
    @test length(𝒢c_dissipators) == 1   # one dissipator

    # Drift is Hamiltonian only — does not contain dissipator iso_D contribution
    n = sys.levels
    P = Piccolo.Quantum.Isomorphisms.density_projection_matrix(n)
    L = Piccolo.Quantum.Isomorphisms.density_lift_matrix(n)
    expected_drift_ham =
        P *
        Piccolo.Quantum.Isomorphisms.G(Piccolo.Quantum.Isomorphisms.ad_vec(sys.H_drift)) *
        L
    @test isapprox(𝒢c_drift_ham, expected_drift_ham; atol = 1e-10)

    # Dissipator factor is P · iso_D(sqrt(γ)·Z) · L
    expected_diss =
        P * Piccolo.Quantum.Isomorphisms.iso_D(sparse(ComplexF64.(sqrt(γ) * PAULIS.Z))) * L
    @test isapprox(𝒢c_dissipators[1], expected_diss; atol = 1e-10)
end

@testitem "compact_lindbladian_generators: legacy 2-tuple bridge" begin
    using Piccolo
    using SparseArrays

    # No-dissipator system: 2-tuple drift equals 3-tuple drift_ham.
    sys_closed = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])
    drift, drives = compact_lindbladian_generators(sys_closed)
    drift_ham, drives_parts, _ = compact_lindbladian_parts(sys_closed)
    @test isapprox(drift, drift_ham; atol = 1e-12)
    @test length(drives) == length(drives_parts) == 1

    # LinearDissipator system: 2-tuple drift fuses iso_D contribution at u=0.
    γ = 0.3
    sys_open = OpenQuantumSystem(
        PAULIS.Z,
        [PAULIS.X],
        [1.0];
        dissipators = [LinearDissipator(sqrt(γ) * PAULIS.Z)],
    )
    drift_legacy, _ = compact_lindbladian_generators(sys_open)
    drift_ham, _, dissipators = compact_lindbladian_parts(sys_open)
    @test isapprox(drift_legacy, drift_ham + dissipators[1]; atol = 1e-12)
end

@testitem "compact_generator_closure: assembles from per-drive and per-dissipator factors" begin
    using Piccolo
    using SparseArrays

    # Use 3-drive bounds so validate_drive_jacobian can sample len-3 u vectors
    drive_global = NonlinearDrive(PAULIS.X, u -> u[1] * u[2]; active_controls = [1, 2])
    diss_global = NonlinearDissipator(PAULIS.Z / sqrt(2), u -> u[3]; active_controls = [3])
    sys = OpenQuantumSystem(
        PAULIS.Z,
        AbstractDrive[drive_global],
        [1.0, 1.0, 1.0];
        dissipators = [diss_global],
        global_params = (θ = 1.0, γ = 0.1),
    )

    𝒢c_drift_ham, 𝒢c_drives, 𝒢c_dissipators = compact_lindbladian_parts(sys)
    G_fn = compact_generator_closure(sys, 𝒢c_drift_ham, 𝒢c_drives, 𝒢c_dissipators)

    u_base = [0.3, 1.0, 0.1]
    u_θ2 = [0.3, 2.0, 0.1]
    u_γ2 = [0.3, 1.0, 0.5]

    @test G_fn(u_θ2) ≠ G_fn(u_base)
    @test G_fn(u_γ2) ≠ G_fn(u_base)
    expected =
        Matrix(𝒢c_drift_ham) +
        (0.3 * 1.0) * Matrix(𝒢c_drives[1]) +
        0.1 * Matrix(𝒢c_dissipators[1])
    @test isapprox(Matrix(G_fn(u_base)), expected; atol = 1e-10)
end

@testitem "OpenQuantumSystem.𝒢: respects NonlinearDrive with global-reading coeff" begin
    using Piccolo
    using SparseArrays

    # NonlinearDrive with coefficient u[1]*u[2] — the extended control vector
    # from the integrator will place a "global" read into u[2].
    # drive_bounds has two entries so validate_drive_jacobian samples u of len 2.
    drive_global = NonlinearDrive(PAULIS.X, u -> u[1] * u[2]; active_controls = [1, 2])
    sys = OpenQuantumSystem(
        PAULIS.Z,
        AbstractDrive[drive_global],
        [1.0, 1.0];
        global_params = (θ = 1.0,),
    )

    G1 = sys.𝒢([0.3, 1.0])
    G2 = sys.𝒢([0.3, 2.0])
    @test !isapprox(G1, G2; atol = 1e-10)
    # Δ = (0.3 * 2.0 - 0.3 * 1.0) · 𝒢_drive(X) = 0.3 · 𝒢_drive(X)
    𝒢X = Isomorphisms.G(Isomorphisms.ad_vec(sparse(ComplexF64.(PAULIS.X))))
    @test isapprox(G2 - G1, (2.0 - 1.0) * 0.3 * 𝒢X; atol = 1e-10)
end

@testitem "OpenQuantumSystem.𝒢: respects NonlinearDissipator with global-reading rate" begin
    using Piccolo
    using SparseArrays

    L = PAULIS.Z / sqrt(2)
    diss_global = NonlinearDissipator(L, u -> u[2]; active_controls = [2])
    sys = OpenQuantumSystem(
        PAULIS.X,
        [PAULIS.Y],
        [1.0];
        dissipators = [diss_global],
        global_params = (γ = 0.1,),
    )

    G1 = sys.𝒢([0.0, 0.1])
    G2 = sys.𝒢([0.0, 0.5])
    @test !isapprox(G1, G2; atol = 1e-10)
    # Δ = (0.5 - 0.1) · iso_D(L)
    expected_delta = (0.5 - 0.1) * Isomorphisms.iso_D(sparse(ComplexF64.(L)))
    @test isapprox(G2 - G1, expected_delta; atol = 1e-10)
end

@testitem "OpenQuantumSystem.𝒢: regression — closed-system LinearDrives only" begin
    using Piccolo
    using SparseArrays
    # Before-the-rewrite semantics: linear drives + static dissipators still
    # yield the correct 𝒢(u) for arbitrary u.
    γ = 0.3
    sys = OpenQuantumSystem(
        PAULIS.Z,
        [PAULIS.X, PAULIS.Y],
        [1.0, 1.0];
        dissipation_operators = [sqrt(γ) * PAULIS.Z],
    )
    u = [0.4, -0.2]
    # Manually assemble: 𝒢_drift + u[1]·𝒢_X + u[2]·𝒢_Y + iso_D(sqrt(γ)·Z)
    𝒢d = Isomorphisms.G(Isomorphisms.ad_vec(sparse(ComplexF64.(PAULIS.Z))))
    𝒢x = Isomorphisms.G(Isomorphisms.ad_vec(sparse(ComplexF64.(PAULIS.X))))
    𝒢y = Isomorphisms.G(Isomorphisms.ad_vec(sparse(ComplexF64.(PAULIS.Y))))
    𝒟 = Isomorphisms.iso_D(sparse(ComplexF64.(sqrt(γ) * PAULIS.Z)))
    expected = 𝒢d + u[1] * 𝒢x + u[2] * 𝒢y + 𝒟
    @test isapprox(sys.𝒢(u), expected; atol = 1e-10)
end

@testitem "OpenQuantumSystem: accepts typed dissipators" begin
    using Piccolo

    # Matrix input → auto-wrapped to LinearDissipator
    sys_mat =
        OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [PAULIS.Z])
    @test length(sys_mat.dissipators) == 1
    @test sys_mat.dissipators[1] isa LinearDissipator
    @test sys_mat.dissipators[1].rate == 1.0
    @test length(sys_mat.dissipation_operators) == 1  # backward-compat view

    # Typed input → preserved as-is
    diss = NonlinearDissipator(PAULIS.Z, u -> u[2]; active_controls = [2])
    sys_typed = OpenQuantumSystem(
        PAULIS.Z,
        [PAULIS.X],
        [1.0];
        dissipators = [diss],
        global_params = (γ = 0.1,),
    )
    @test sys_typed.dissipators[1] === diss
    @test length(sys_typed.dissipation_operators) == 1

    # Passing both kwargs should error
    @test_throws ErrorException OpenQuantumSystem(
        PAULIS.Z,
        [PAULIS.X],
        [1.0];
        dissipation_operators = [PAULIS.Z],
        dissipators = [LinearDissipator(PAULIS.Z)],
    )
end

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

    # H_drives field should be auto-populated with LinearDrives
    @test length(system.H_drives) == 1
    @test system.H_drives[1] isa LinearDrive
    @test system.H_drives[1].index == 1

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
    @test length(system.H_drives) == 1
    @test system.H_drives[1] isa LinearDrive

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
        H_drift,
        drives,
        [1.0, 1.0];
        dissipation_operators = dissipation_operators,
    )

    @test sys isa OpenQuantumSystem
    @test sys.n_drives == 2
    @test length(sys.H_drives) == 3
    @test has_nonlinear_drives(sys.H_drives)

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
    @test length(osys.H_drives) == 3
    @test has_nonlinear_drives(osys.H_drives)
end

@testitem "OpenQuantumSystem with modulated drives" begin
    using SparseArrays
    using LinearAlgebra

    H_drift = PAULIS.Z
    dissipation_operators = [0.1 * PAULIS[:X]]

    # Modulated linear drive via Pair syntax on QuantumSystem
    omega = 2.0
    qsys = QuantumSystem(H_drift, [PAULIS.X => t -> cos(omega * t)], [1.0])
    osys = OpenQuantumSystem(qsys; dissipation_operators = dissipation_operators)

    @test osys isa OpenQuantumSystem
    @test osys.H_drives[1] isa ModulatedDrive

    # H(u, t) should reflect time modulation
    u_test = [0.5]
    t_test = 0.3
    H_result = osys.H(u_test, t_test)
    H_expected = PAULIS.Z + 0.5 * cos(omega * t_test) * PAULIS.X
    @test norm(H_result - H_expected) < 1e-10

    # Different time should give different result
    H_result2 = osys.H(u_test, 0.0)
    H_expected2 = PAULIS.Z + 0.5 * PAULIS.X
    @test norm(H_result2 - H_expected2) < 1e-10

    # 𝒢 function should include dissipation
    𝒢_result = osys.𝒢(u_test)
    @test size(𝒢_result) == (8, 8)

    # Modulated nonlinear drive
    nd = NonlinearDrive(PAULIS.Y, u -> u[1]^2 + u[2]^2)
    drives = AbstractDrive[
        LinearDrive(sparse(ComplexF64.(PAULIS.X)), 1),
        ModulatedDrive(nd, t -> sin(omega * t)),
    ]
    osys2 = OpenQuantumSystem(
        H_drift,
        drives,
        [1.0, 1.0];
        dissipation_operators = dissipation_operators,
    )

    @test osys2 isa OpenQuantumSystem
    u2 = [0.3, 0.7]
    t2 = 0.5
    coeff_nl = (0.3^2 + 0.7^2) * sin(omega * t2)
    H_expected3 = PAULIS.Z + 0.3 * PAULIS.X + coeff_nl * PAULIS.Y
    @test norm(osys2.H(u2, t2) - H_expected3) < 1e-10
end
