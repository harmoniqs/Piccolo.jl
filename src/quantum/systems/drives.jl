# ----------------------------------------------------------------------------- #
# Drive Types: Operator + Coefficient for Hamiltonian Terms
#
# A Drive pairs a Hamiltonian matrix with a scalar coefficient function of the
# control vector u and time t. This enables linear drives (u[i] * H_i),
# nonlinear drives (f(u) * H), and time-modulated drives (f(u) * b(t) * H).
#
# The Hamiltonian is: H(u, t) = Σ_dt a(t)*H_dt + Σ_d drive_coeff(d, u, t) * H_d
#
# For analytical sensitivity equations (e.g., spline integrators), the chain
# rule requires: ∂H/∂u_j = Σ_d drive_coeff_jac(d, u, t, j) * H_d
# For second-order sensitivity equations (exact Hessian):
# ∂²H/∂u_i∂u_j = Σ_d drive_coeff_hess(d, u, i, j) * d.H
# ----------------------------------------------------------------------------- #

"""
    AbstractDrive

Abstract supertype for Hamiltonian drive terms.

A drive pairs a Hermitian matrix `H` with a scalar coefficient that depends on the
control vector `u` and (optionally) time `t`. The full Hamiltonian is:

    H(u, t) = Σ_k a_k(t) * H_drift_k + Σ_d drive_coeff(d, u, t) * H_d

Subtypes must implement:
- `drive_coeff(d, u, t)` — compute the scalar coefficient
- `drive_coeff_jac(d, u, t, j)` — compute ∂coeff/∂u_j
- `drive_coeff_hess(d, u, i, j)` — compute ∂²coeff/∂u_i∂u_j
- `active_controls(d)` — return indices where ∂coeff/∂u_j can be nonzero

See [`LinearDrive`](@ref), [`NonlinearDrive`](@ref), [`ModulatedDrive`](@ref).
"""
abstract type AbstractDrive end

"""
    LinearDrive{H} <: AbstractDrive

Standard linear drive: coefficient is `u[index]`.

This is the default representation when constructing `QuantumSystem(H_drift, H_drives, bounds)`.

# Fields
- `H::H`: The Hermitian drive operator (matrix or `AbstractDynamicsOperator`)
- `index::Int`: Index into the control vector `u`

# Example
```julia
LinearDrive(sparse(PAULIS.X), 1)  # u[1] * σx
```
"""
struct LinearDrive{H} <: AbstractDrive
    H::H
    index::Int
end

"""
    NonlinearDrive{H,F,DF,HF} <: AbstractDrive

Drive term with a nonlinear scalar coefficient: `f(u) * H`.

The user provides either:
1. Both the coefficient function and its Jacobian (3-arg form), or
2. Just the coefficient function (2-arg form) — the Jacobian and Hessian are
   computed automatically via ForwardDiff.

The Hessian (`coeff_hess`) is used by the second-order sensitivity ODE for exact
Hessian computation in SplineIntegrator (`exact_hessian=true`).

# Fields
- `H::H`: The Hermitian drive operator (matrix or `AbstractDynamicsOperator`)
- `coeff::F`: `(u::AbstractVector) -> scalar` — coefficient function
- `coeff_jac::DF`: `(u::AbstractVector, j::Int) -> scalar` — ∂coeff/∂u_j
- `coeff_hess::HF`: `(u::AbstractVector, i::Int, j::Int) -> scalar` — ∂²coeff/∂u_i∂u_j
- `active_controls::Vector{Int}`: Control indices where ∂coeff/∂u_j can be nonzero.
  Empty means "all controls" (no structural sparsity info).

# Example: Auto-Jacobian (recommended)
```julia
NonlinearDrive(σz / 2, u -> u[3]^2 + u[4]^2)
```

# Example: Manual Jacobian with active controls
```julia
NonlinearDrive(
    σz / 2,
    u -> u[3]^2 + u[4]^2,
    (u, j) -> j == 3 ? 2u[3] : j == 4 ? 2u[4] : 0.0;
    active_controls = [3, 4]
)
```
"""
struct NonlinearDrive{H,F,DF,HF} <: AbstractDrive
    H::H
    coeff::F
    coeff_jac::DF
    coeff_hess::HF
    active_controls::Vector{Int}
end

"""
    NonlinearDrive(H::AbstractMatrix, coeff, coeff_jac; active_controls=Int[], coeff_hess=nothing)

Construct a `NonlinearDrive` with an explicit Jacobian function.
`active_controls` lists which control indices have nonzero ∂coeff/∂u_j (empty = all).
If `coeff_hess` is not provided, it is auto-generated via ForwardDiff from `coeff`.
"""
function NonlinearDrive(
    H::AbstractMatrix,
    coeff::F,
    coeff_jac::DF;
    active_controls::Vector{Int} = Int[],
    coeff_hess = nothing,
) where {F,DF}
    hess = isnothing(coeff_hess) ? _forwarddiff_hessian(coeff) : coeff_hess
    return NonlinearDrive(sparse(ComplexF64.(H)), coeff, coeff_jac, hess, active_controls)
end

"""
    NonlinearDrive(H::AbstractMatrix, coeff; active_controls=Int[])

Construct a `NonlinearDrive` with an auto-generated Jacobian via ForwardDiff.

This is the recommended constructor — the Jacobian is computed automatically from
the coefficient function, eliminating the risk of hand-written Jacobian errors.

# Example
```julia
# Displaced-frame Stark shift: coefficient = u₃² + u₄²
drive = NonlinearDrive(σz / 2, u -> u[3]^2 + u[4]^2; active_controls = [3, 4])
```
"""
function NonlinearDrive(
    H::AbstractMatrix,
    coeff::F;
    active_controls::Vector{Int} = Int[],
) where {F}
    coeff_jac = _forwarddiff_jacobian(coeff)
    coeff_hess = _forwarddiff_hessian(coeff)
    return NonlinearDrive(
        sparse(ComplexF64.(H)),
        coeff,
        coeff_jac,
        coeff_hess,
        active_controls,
    )
end

# ── Generic constructors (for structured operators or any non-matrix Hamiltonian) ──

"""
    NonlinearDrive(H, coeff; active_controls=Int[])

Construct a `NonlinearDrive` with any Hamiltonian type (e.g., `AbstractDynamicsOperator`).
The Hamiltonian is stored as-is (not materialized or sparsified).
"""
function NonlinearDrive(H, coeff::F; active_controls::Vector{Int} = Int[]) where {F}
    coeff_jac = _forwarddiff_jacobian(coeff)
    coeff_hess = _forwarddiff_hessian(coeff)
    return NonlinearDrive(H, coeff, coeff_jac, coeff_hess, active_controls)
end

"""
    NonlinearDrive(H, coeff, coeff_jac; active_controls=Int[], coeff_hess=nothing)

Construct a `NonlinearDrive` with any Hamiltonian type and explicit Jacobian.
"""
function NonlinearDrive(
    H,
    coeff::F,
    coeff_jac::DF;
    active_controls::Vector{Int} = Int[],
    coeff_hess = nothing,
) where {F,DF}
    hess = isnothing(coeff_hess) ? _forwarddiff_hessian(coeff) : coeff_hess
    return NonlinearDrive(H, coeff, coeff_jac, hess, active_controls)
end

"""
    BilinearCouplerCoeff

Coefficient functor for [`coupling_drive`](@ref): `a · u[i] · u[j]`.

The drive term it parameterizes is *nonlinear* in the control vector (a
product of two controls — hence it lives inside a [`NonlinearDrive`](@ref))
but *bilinear* with respect to the two coupler controls `u[i]` and `u[j]`,
with a fixed (baked) coupling strength `a`.

A named functor — rather than an anonymous closure — keeps the drive
introspectable (`d.coeff.i`, `d.coeff.j`, `d.coeff.a`), serializable
(plain-data fields, no closure types), and concretely typed (all couplings
share one `NonlinearDrive{H,BilinearCouplerCoeff,...}` type).
"""
struct BilinearCouplerCoeff
    i::Int
    j::Int
    a::Float64
end
(c::BilinearCouplerCoeff)(u) = c.a * u[c.i] * u[c.j]

"""
    BilinearCouplerCoeffJac

Analytic Jacobian functor companion to [`BilinearCouplerCoeff`](@ref):
`∂(a·u[i]·u[j])/∂u_p`. Avoids the per-query ForwardDiff gradient fallback.
"""
struct BilinearCouplerCoeffJac
    i::Int
    j::Int
    a::Float64
end
(c::BilinearCouplerCoeffJac)(u, p) = p == c.i ? c.a * u[c.j] : p == c.j ? c.a * u[c.i] : 0.0

"""
    BilinearCouplerCoeffHess

Analytic Hessian functor companion to [`BilinearCouplerCoeff`](@ref):
`∂²(a·u[i]·u[j])/∂u_p∂u_q` — constant `a` on the `(i,j)`/`(j,i)` entries,
zero elsewhere.
"""
struct BilinearCouplerCoeffHess
    i::Int
    j::Int
    a::Float64
end
(c::BilinearCouplerCoeffHess)(u, p, q) =
    ((p == c.i && q == c.j) || (p == c.j && q == c.i)) ? c.a : 0.0

"""
    coupling_drive(H, i::Int, j::Int; strength=1.0)

Construct a [`NonlinearDrive`](@ref) for a pairwise coupling term with a
**fixed (baked) coupling strength**: `strength · u[i] · u[j] · H`.

The returned drive *is* a `NonlinearDrive` — parameterized by the named
functors [`BilinearCouplerCoeff`](@ref) / [`BilinearCouplerCoeffJac`](@ref) /
[`BilinearCouplerCoeffHess`](@ref) instead of anonymous closures, so coupled
systems stay introspectable, serializable, and concretely typed.

This is the recommended way to express coupler-mediated interactions of the
form `A_ij · u_i(t) · u_j(t) · H_ij` (e.g. photonic networks where a mixing
circuit couples modes `i` and `j` with controllable rates, as in
recirculating quantum photonic networks) when the coupling matrix entry
`A_ij` is a calibration constant rather than an optimization variable.

Baking the strength into the drive — instead of appending it to the control
vector as a global parameter — removes it from the NLP entirely: no extra
sensitivity column in every interval's Jacobian ODE, and no gauge-flat
Hessian directions from the rescaling redundancy
`u_i → α·u_i, A_ij → A_ij/α`. Empirically (rqpnn C-matrix campaign,
2026-06-10) this improves both convergence quality and wall time relative
to co-optimizing the coupling strengths.

The analytic Jacobian and Hessian are provided (no ForwardDiff fallback),
and `active_controls = [i, j]` gives the integrator exact structural
sparsity.

# Example
```julia
# A_12 = 0.85 coupling between modes 1 and 2:
d = coupling_drive(ad_2 * a_1 + ad_1 * a_2, 1, 2; strength = 0.85)
d.coeff.a  # 0.85 — readable back out of any saved system
```
"""
function coupling_drive(H, i::Int, j::Int; strength::Real = 1.0)
    i == j && throw(ArgumentError("coupling_drive requires i ≠ j (got i = j = $i)"))
    a = Float64(strength)
    return NonlinearDrive(
        H,
        BilinearCouplerCoeff(i, j, a),
        BilinearCouplerCoeffJac(i, j, a);
        coeff_hess = BilinearCouplerCoeffHess(i, j, a),
        active_controls = [i, j],
    )
end

"""
    _forwarddiff_jacobian(f) -> (u, j) -> ∂f/∂u_j

Create a Jacobian function from a scalar-valued function `f(u)` using ForwardDiff.
"""
function _forwarddiff_jacobian(f::F) where {F}
    return (u, j) -> ForwardDiff.gradient(f, u)[j]
end

"""
    _forwarddiff_hessian(f) -> (u, i, j) -> ∂²f/∂u_i∂u_j

Create a Hessian function from a scalar-valued function `f(u)` using ForwardDiff.
"""
function _forwarddiff_hessian(f::F) where {F}
    return (u, i, j) -> ForwardDiff.hessian(f, u)[i, j]
end

# ----------------------------------------------------------------------------- #
# Modulation Types
# ----------------------------------------------------------------------------- #

const _identity_modulation = Returns(1.0)
const _zero_modulation_deriv = Returns(0.0)

"""
    DriftTerm{H, F, DF}

A drift Hamiltonian term with an optional time-dependent modulation: `a(t) * H`.

# Fields
- `H::H`: The Hermitian drift operator (matrix or `AbstractDynamicsOperator`)
- `modulation::F`: `t -> scalar` modulation function (default: `t -> 1`)
- `modulation_deriv::DF`: `t -> da/dt`, pre-computed via ForwardDiff

# Constructors
```julia
DriftTerm(H)                    # identity modulation
DriftTerm(H, t -> cos(omega*t)) # auto-computes derivative
```
"""
struct DriftTerm{H,F,DF}
    H::H
    modulation::F
    modulation_deriv::DF
end

function DriftTerm(H, modulation::F) where {F}
    modulation_deriv = t -> ForwardDiff.derivative(modulation, t)
    # Validate at construction
    @assert modulation(0.0) isa Real "Modulation must return a real scalar"
    modulation_deriv(0.0)  # trigger ForwardDiff to catch errors early
    return DriftTerm(H, modulation, modulation_deriv)
end

DriftTerm(H) = DriftTerm(H, _identity_modulation, _zero_modulation_deriv)

"""
    has_modulation(dt::DriftTerm) -> Bool

Return `true` if this drift term has non-identity modulation.
"""
has_modulation(dt::DriftTerm) = dt.modulation !== _identity_modulation

"""
    ModulatedDrive{D<:AbstractDrive, F, DF} <: AbstractDrive

A drive term with time-dependent modulation: `drive_coeff(base, u, t) * b(t) * H`.

Composes with any `AbstractDrive` subtype. The modulation `b(t)` is independent
of the control vector `u`.

# Fields
- `base::D`: The underlying drive (e.g., `LinearDrive`, `NonlinearDrive`)
- `modulation::F`: `t -> scalar` modulation function
- `modulation_deriv::DF`: `t -> db/dt`, pre-computed via ForwardDiff

# Example
```julia
# u[1] * cos(omega*t) * H_x
ModulatedDrive(LinearDrive(H_x, 1), t -> cos(omega * t))
```
"""
struct ModulatedDrive{D<:AbstractDrive,F,DF} <: AbstractDrive
    base::D
    modulation::F
    modulation_deriv::DF
end

function ModulatedDrive(base::D, modulation::F) where {D<:AbstractDrive,F}
    modulation_deriv = t -> ForwardDiff.derivative(modulation, t)
    @assert modulation(0.0) isa Real "Modulation must return a real scalar"
    modulation_deriv(0.0)
    return ModulatedDrive(base, modulation, modulation_deriv)
end

has_modulation(d::ModulatedDrive) = true
has_modulation(::AbstractDrive) = false

# ----------------------------------------------------------------------------- #
# Interface (all signatures require t for time-modulation support)
# ----------------------------------------------------------------------------- #

"""
    drive_coeff(d::AbstractDrive, u::AbstractVector, t) -> Number

Compute the scalar coefficient of this drive at controls `u` and time `t`.
"""
@inline drive_coeff(d::LinearDrive, u::AbstractVector, _t) = u[d.index]
@inline drive_coeff(d::NonlinearDrive, u::AbstractVector, _t) = d.coeff(u)
@inline drive_coeff(d::ModulatedDrive, u::AbstractVector, t) =
    drive_coeff(d.base, u, t) * d.modulation(t)

"""
    drive_coeff_jac(d::AbstractDrive, u::AbstractVector, t, j::Int) -> Number

Compute ∂coeff/∂u_j at controls `u` and time `t`.

For `LinearDrive`: returns 1.0 if `j == d.index`, 0.0 otherwise (Kronecker delta).
For `NonlinearDrive`: evaluates the Jacobian function (user-provided or auto-generated).
For `ModulatedDrive`: scales the base Jacobian by the modulation at time `t`.
"""
@inline drive_coeff_jac(d::LinearDrive, ::AbstractVector, _t, j::Int) =
    j == d.index ? 1.0 : 0.0
@inline drive_coeff_jac(d::NonlinearDrive, u::AbstractVector, _t, j::Int) =
    d.coeff_jac(u, j)
@inline drive_coeff_jac(d::ModulatedDrive, u::AbstractVector, t, j::Int) =
    drive_coeff_jac(d.base, u, t, j) * d.modulation(t)

"""
    drive_coeff_dt(d::AbstractDrive, u, t) -> Number

Time derivative of the drive coefficient: d/dt[drive_coeff(d, u, t)].
Returns zero for unmodulated drives.
"""
@inline drive_coeff_dt(::LinearDrive, u::AbstractVector, t) = 0.0
@inline drive_coeff_dt(::NonlinearDrive, u::AbstractVector, t) = 0.0
@inline drive_coeff_dt(d::ModulatedDrive, u::AbstractVector, t) =
    drive_coeff(d.base, u, t) * d.modulation_deriv(t)

# Backward-compatible 2-arg shims (remove in PR 2 when Piccolissimo gains t)
@inline drive_coeff(d::AbstractDrive, u::AbstractVector) = drive_coeff(d, u, 0.0)
@inline drive_coeff_jac(d::AbstractDrive, u::AbstractVector, j::Int) =
    drive_coeff_jac(d, u, 0.0, j)

"""
    drive_coeff_hess(d::AbstractDrive, u::AbstractVector, i::Int, j::Int) -> Number

Compute ∂²coeff/∂u_i∂u_j at controls `u`.

For `LinearDrive`: always returns 0.0 (linear coefficient has zero Hessian).
For `NonlinearDrive`: evaluates the Hessian function (user-provided or auto-generated).
"""
@inline drive_coeff_hess(d::LinearDrive, ::AbstractVector, ::Int, ::Int) = 0.0
@inline drive_coeff_hess(d::NonlinearDrive, u::AbstractVector, i::Int, j::Int) =
    d.coeff_hess(u, i, j)
@inline drive_coeff_hess(d::ModulatedDrive, u::AbstractVector, i::Int, j::Int) =
    drive_coeff_hess(d.base, u, i, j) * d.modulation(0.0)

"""
    active_controls(d::AbstractDrive) -> Vector{Int}

Return the control indices where `drive_coeff_jac(d, u, j)` can be nonzero.
An empty vector means "all controls" (no structural sparsity information available).

For `LinearDrive`, returns `[d.index]`.
For `NonlinearDrive`, returns `d.active_controls` (user-specified or empty).
"""
active_controls(d::LinearDrive) = [d.index]
active_controls(d::NonlinearDrive) = d.active_controls

"""
    drive_matrix(d::AbstractDrive) -> AbstractMatrix

Return the Hermitian operator matrix for this drive term.
"""
drive_matrix(d::LinearDrive) = d.H
drive_matrix(d::NonlinearDrive) = d.H

"""
    drive_dim(d::AbstractDrive) -> Int

Return the size (number of rows) of the drive operator.
"""
drive_dim(d::LinearDrive) = size(d.H, 1)
drive_dim(d::NonlinearDrive) = size(d.H, 1)

# ModulatedDrive interface delegation
drive_matrix(d::ModulatedDrive) = drive_matrix(d.base)
drive_dim(d::ModulatedDrive) = drive_dim(d.base)
active_controls(d::ModulatedDrive) = active_controls(d.base)
Isomorphisms.G(d::ModulatedDrive) = Isomorphisms.G(d.base)

# Update has_nonlinear_drives to unwrap ModulatedDrive
_is_nonlinear(::AbstractDrive) = false
_is_nonlinear(::NonlinearDrive) = true
_is_nonlinear(d::ModulatedDrive) = _is_nonlinear(d.base)

"""
    has_nonlinear_drives(drives::Vector{<:AbstractDrive}) -> Bool

Check if any drive terms are nonlinear.
"""
has_nonlinear_drives(drives::AbstractVector{<:AbstractDrive}) = any(_is_nonlinear, drives)

# ----------------------------------------------------------------------------- #
# Matrix access
# ----------------------------------------------------------------------------- #

"""
    _ensure_matrix(H) -> AbstractMatrix

Convert `H` to a matrix if it isn't one already. For `AbstractMatrix` inputs,
returns `H` directly. For operator types (e.g., from Piccolissimo's operators
module), falls back to `Matrix(H)` — operator packages should define
`Base.Matrix(op)` to return a dense matrix representation.
"""
_ensure_matrix(H::AbstractMatrix) = H
_ensure_matrix(H) = Matrix(H)

"""
    drive_matrix(d::AbstractDrive) -> AbstractMatrix

Return the drive operator as a matrix. Use this instead of `d.H` when a matrix
is required (e.g., for `Isomorphisms.G`, Hermiticity checks, or `H(u,t)` closures).
"""
drive_matrix(d::AbstractDrive) = _ensure_matrix(d.H)

"""
    drive_dim(d::AbstractDrive) -> Int

Return the Hilbert space dimension of the drive operator.
"""
drive_dim(d::AbstractDrive) = size(drive_matrix(d), 1)

# ----------------------------------------------------------------------------- #
# Isomorphism dispatch for drives
# ----------------------------------------------------------------------------- #

"""
    Isomorphisms.G(d::AbstractDrive)

Delegate `G` to the underlying Hamiltonian matrix `d.H`, so that broadcasting
`G.(sys.H_drives)` works when `H_drives` contains `AbstractDrive` objects.
"""
Isomorphisms.G(d::AbstractDrive) = Isomorphisms.G(drive_matrix(d))

# ----------------------------------------------------------------------------- #
# Jacobian Validation
# ----------------------------------------------------------------------------- #

"""
    validate_drive_jacobian(d::NonlinearDrive, u_dim::Int; atol=1e-6, n_samples=3)

Spot-check the Jacobian of a `NonlinearDrive` against ForwardDiff at random control vectors.
Throws an `AssertionError` if the user-provided Jacobian disagrees with the AD Jacobian.

This is called automatically during `QuantumSystem` construction for all `NonlinearDrive`
terms, catching sign errors or off-by-one bugs early. Sampling uses a local RNG so
construction does not advance the global `Random` stream as a side effect.
"""
function validate_drive_jacobian(
    d::NonlinearDrive,
    u_dim::Int;
    atol::Float64 = 1e-6,
    n_samples::Int = 3,
)
    rng = MersenneTwister(0)
    for _ = 1:n_samples
        u = randn(rng, u_dim)
        grad_ad = ForwardDiff.gradient(d.coeff, u)
        for j = 1:u_dim
            user_val = d.coeff_jac(u, j)
            @assert abs(user_val - grad_ad[j]) < atol (
                "NonlinearDrive Jacobian mismatch at u=$u, j=$j: " *
                "user=$(user_val), ForwardDiff=$(grad_ad[j])"
            )
        end
    end
end

"""
    validate_drive_hessian(d::NonlinearDrive, u_dim::Int; atol=1e-6, n_samples=3)

Spot-check the Hessian of a `NonlinearDrive` against ForwardDiff at random control vectors.
Throws an `AssertionError` if the user-provided Hessian disagrees with the AD Hessian.

Sampling uses a local RNG so construction does not advance the global `Random` stream
as a side effect.
"""
function validate_drive_hessian(
    d::NonlinearDrive,
    u_dim::Int;
    atol::Float64 = 1e-6,
    n_samples::Int = 3,
)
    rng = MersenneTwister(0)
    for _ = 1:n_samples
        u = randn(rng, u_dim)
        hess_ad = ForwardDiff.hessian(d.coeff, u)
        for i = 1:u_dim, j = i:u_dim
            user_val = d.coeff_hess(u, i, j)
            @assert abs(user_val - hess_ad[i, j]) < atol (
                "NonlinearDrive Hessian mismatch at u=$u, (i,j)=($i,$j): " *
                "user=$(user_val), ForwardDiff=$(hess_ad[i, j])"
            )
        end
    end
end

# ----------------------------------------------------------------------------- #
# Tests
# ----------------------------------------------------------------------------- #

@testitem "LinearDrive" begin
    using Piccolo
    using SparseArrays

    H = sparse([0.0+0im 1.0+0im; 1.0+0im 0.0+0im])
    d = LinearDrive(H, 2)

    u = [0.3, 0.7, 0.1]
    @test drive_coeff(d, u) == 0.7
    @test drive_coeff_jac(d, u, 1) == 0.0
    @test drive_coeff_jac(d, u, 2) == 1.0
    @test drive_coeff_jac(d, u, 3) == 0.0
    @test d.H === H

    # Hessian is always zero for linear drives
    @test drive_coeff_hess(d, u, 1, 1) == 0.0
    @test drive_coeff_hess(d, u, 2, 2) == 0.0
    @test drive_coeff_hess(d, u, 1, 2) == 0.0

    # active_controls
    @test active_controls(d) == [2]
end

@testitem "NonlinearDrive with explicit Jacobian" begin
    using Piccolo
    using SparseArrays

    H = sparse([1.0+0im 0.0+0im; 0.0+0im -1.0+0im])

    # Quadratic: u[1]^2 + u[2]^2
    d = NonlinearDrive(
        H,
        u -> u[1]^2 + u[2]^2,
        (u, j) -> j == 1 ? 2u[1] : j == 2 ? 2u[2] : 0.0,
    )

    u = [3.0, 4.0, 0.0]
    @test drive_coeff(d, u) == 25.0
    @test drive_coeff_jac(d, u, 1) == 6.0
    @test drive_coeff_jac(d, u, 2) == 8.0
    @test drive_coeff_jac(d, u, 3) == 0.0
    @test d.H == H
    @test isempty(active_controls(d))  # default: empty = all controls

    # Hessian: ∂²(u₁² + u₂²)/∂u_i∂u_j
    @test drive_coeff_hess(d, u, 1, 1) ≈ 2.0  # ∂²/∂u₁² = 2
    @test drive_coeff_hess(d, u, 2, 2) ≈ 2.0  # ∂²/∂u₂² = 2
    @test drive_coeff_hess(d, u, 1, 2) ≈ 0.0  # ∂²/∂u₁∂u₂ = 0
    @test drive_coeff_hess(d, u, 3, 3) ≈ 0.0  # u₃ not involved

    # Product: u[1] * u[2]
    d2 = NonlinearDrive(H, u -> u[1] * u[2], (u, j) -> j == 1 ? u[2] : j == 2 ? u[1] : 0.0)

    @test drive_coeff(d2, u) == 12.0
    @test drive_coeff_jac(d2, u, 1) == 4.0
    @test drive_coeff_jac(d2, u, 2) == 3.0

    # Hessian: ∂²(u₁u₂)/∂u_i∂u_j
    @test drive_coeff_hess(d2, u, 1, 2) ≈ 1.0  # ∂²/∂u₁∂u₂ = 1
    @test drive_coeff_hess(d2, u, 2, 1) ≈ 1.0  # symmetric
    @test drive_coeff_hess(d2, u, 1, 1) ≈ 0.0  # ∂²/∂u₁² = 0
    @test drive_coeff_hess(d2, u, 2, 2) ≈ 0.0  # ∂²/∂u₂² = 0

    # With active_controls
    d3 = NonlinearDrive(
        H,
        u -> u[1]^2 + u[2]^2,
        (u, j) -> j == 1 ? 2u[1] : j == 2 ? 2u[2] : 0.0;
        active_controls = [1, 2],
    )
    @test active_controls(d3) == [1, 2]
end

@testitem "NonlinearDrive auto-Jacobian" begin
    using Piccolo
    using SparseArrays

    H = sparse([1.0+0im 0.0+0im; 0.0+0im -1.0+0im])

    # 2-arg constructor: auto-Jacobian via ForwardDiff
    d = NonlinearDrive(H, u -> u[1]^2 + u[2]^2)

    u = [3.0, 4.0, 0.0]
    @test drive_coeff(d, u) == 25.0
    @test drive_coeff_jac(d, u, 1) ≈ 6.0
    @test drive_coeff_jac(d, u, 2) ≈ 8.0
    @test drive_coeff_jac(d, u, 3) ≈ 0.0

    # Product of controls
    d2 = NonlinearDrive(H, u -> u[1] * u[2])
    @test drive_coeff_jac(d2, [3.0, 4.0], 1) ≈ 4.0
    @test drive_coeff_jac(d2, [3.0, 4.0], 2) ≈ 3.0

    # Cubic term
    d3 = NonlinearDrive(H, u -> u[1]^3)
    @test drive_coeff_jac(d3, [2.0], 1) ≈ 12.0
    @test drive_coeff_hess(d3, [2.0], 1, 1) ≈ 12.0  # ∂²(u³)/∂u² = 6u = 12

    # Hessian of auto-generated quadratic
    @test drive_coeff_hess(d, [3.0, 4.0, 0.0], 1, 1) ≈ 2.0
    @test drive_coeff_hess(d, [3.0, 4.0, 0.0], 1, 2) ≈ 0.0

    # With active_controls
    d4 = NonlinearDrive(H, u -> u[2]^2; active_controls = [2])
    @test active_controls(d4) == [2]
    @test drive_coeff_jac(d4, [1.0, 3.0], 2) ≈ 6.0
    @test drive_coeff_hess(d4, [1.0, 3.0], 2, 2) ≈ 2.0
end

@testitem "NonlinearDrive accepts dense matrix" begin
    using Piccolo
    using SparseArrays

    H_dense = [1.0+0im 0.0+0im; 0.0+0im -1.0+0im]

    # 3-arg with dense matrix
    d = NonlinearDrive(H_dense, u -> u[1], (u, j) -> j == 1 ? 1.0 : 0.0)
    @test d.H isa SparseMatrixCSC
    @test drive_coeff(d, [0.5]) == 0.5

    # 2-arg with dense matrix
    d2 = NonlinearDrive(H_dense, u -> u[1]^2)
    @test d2.H isa SparseMatrixCSC
    @test drive_coeff(d2, [3.0]) == 9.0
    @test drive_coeff_jac(d2, [3.0], 1) ≈ 6.0
    @test drive_coeff_hess(d2, [3.0], 1, 1) ≈ 2.0
end

@testitem "NonlinearDrive with explicit Hessian" begin
    using Piccolo
    using SparseArrays

    H = sparse([1.0+0im 0.0+0im; 0.0+0im -1.0+0im])

    # Provide all three: coeff, jac, hess
    d = NonlinearDrive(
        H,
        u -> u[1]^2 + u[2]^2,
        (u, j) -> j == 1 ? 2u[1] : j == 2 ? 2u[2] : 0.0;
        coeff_hess = (u, i, j) -> (i == j && i <= 2) ? 2.0 : 0.0,
    )

    u = [3.0, 4.0]
    @test drive_coeff_hess(d, u, 1, 1) == 2.0
    @test drive_coeff_hess(d, u, 2, 2) == 2.0
    @test drive_coeff_hess(d, u, 1, 2) == 0.0

    # Cubic with explicit hess: u[1]^3
    d2 = NonlinearDrive(
        H,
        u -> u[1]^3,
        (u, j) -> j == 1 ? 3u[1]^2 : 0.0;
        coeff_hess = (u, i, j) -> (i == 1 && j == 1) ? 6u[1] : 0.0,
    )
    @test drive_coeff_hess(d2, [2.0], 1, 1) == 12.0
    @test drive_coeff_hess(d2, [5.0], 1, 1) == 30.0
end

@testitem "validate_drive_jacobian" begin
    using Piccolo
    using SparseArrays

    H = sparse([1.0+0im 0.0+0im; 0.0+0im -1.0+0im])

    # Correct Jacobian should pass
    d_correct = NonlinearDrive(
        H,
        u -> u[1]^2 + u[2]^2,
        (u, j) -> j == 1 ? 2u[1] : j == 2 ? 2u[2] : 0.0,
    )
    validate_drive_jacobian(d_correct, 2)  # should not throw

    # Wrong Jacobian should fail
    d_wrong = NonlinearDrive(
        H,
        u -> u[1]^2 + u[2]^2,
        (u, j) -> j == 1 ? u[1] : j == 2 ? u[2] : 0.0,  # missing factor of 2
    )
    @test_throws AssertionError validate_drive_jacobian(d_wrong, 2)

    # Auto-Jacobian should always pass validation
    d_auto = NonlinearDrive(H, u -> u[1]^2 + u[2]^2)
    validate_drive_jacobian(d_auto, 2)
end

@testitem "validate_drive_hessian" begin
    using Piccolo
    using SparseArrays

    H = sparse([1.0+0im 0.0+0im; 0.0+0im -1.0+0im])

    # Auto-Hessian should always pass
    d_auto = NonlinearDrive(H, u -> u[1]^2 + u[2]^2)
    validate_drive_hessian(d_auto, 2)

    # Correct explicit Hessian should pass
    d_correct = NonlinearDrive(
        H,
        u -> u[1]^2 + u[2]^2,
        (u, j) -> j == 1 ? 2u[1] : j == 2 ? 2u[2] : 0.0;
        coeff_hess = (u, i, j) -> (i == j && i <= 2) ? 2.0 : 0.0,
    )
    validate_drive_hessian(d_correct, 2)

    # Wrong explicit Hessian should fail
    d_wrong = NonlinearDrive(
        H,
        u -> u[1]^2 + u[2]^2,
        (u, j) -> j == 1 ? 2u[1] : j == 2 ? 2u[2] : 0.0;
        coeff_hess = (u, i, j) -> 0.0,  # wrong: should be 2.0 on diagonal
    )
    @test_throws AssertionError validate_drive_hessian(d_wrong, 2)
end

@testitem "validate_drive: u_dim extending past n_drives (global params)" begin
    # Regression: when global_params are appended, callers pass
    # u_dim = n_drives + length(global_params) so coefficients that read into
    # the global-parameter slots can be validated without BoundsError.
    using Piccolo
    using SparseArrays

    H = sparse([0.0+0im 1.0+0im; 1.0+0im 0.0+0im])

    # Coefficient reads u[3] — i.e. the appended global slot when n_drives = 2
    # and one global parameter is present.
    d_auto = NonlinearDrive(H, u -> u[1] * u[2] * u[3])
    validate_drive_jacobian(d_auto, 3)  # must not BoundsError
    validate_drive_hessian(d_auto, 3)

    # Explicit jacobian + hessian agreeing with ForwardDiff across the full
    # u_dim = 3 must pass.
    d_correct = NonlinearDrive(
        H,
        u -> u[1] * u[2] * u[3],
        (u, j) ->
            j == 1 ? u[2] * u[3] : j == 2 ? u[1] * u[3] : j == 3 ? u[1] * u[2] : 0.0;
        coeff_hess = (u, i, j) ->
            (i == 1 && j == 2) || (i == 2 && j == 1) ? u[3] :
            (i == 1 && j == 3) || (i == 3 && j == 1) ? u[2] :
            (i == 2 && j == 3) || (i == 3 && j == 2) ? u[1] : 0.0,
    )
    validate_drive_jacobian(d_correct, 3)
    validate_drive_hessian(d_correct, 3)

    # An explicit jacobian wrong at j = 3 (the global slot) is still caught —
    # the validator must check the extended range, not just 1:n_drives.
    d_wrong_global = NonlinearDrive(
        H,
        u -> u[1] * u[2] * u[3],
        (u, j) -> j == 1 ? u[2] * u[3] : j == 2 ? u[1] * u[3] : 0.0,  # wrong at j=3
    )
    @test_throws AssertionError validate_drive_jacobian(d_wrong_global, 3)
end

@testitem "QuantumSystem construction does not advance the global RNG" begin
    # Regression: validators previously sampled with `randn(u_dim)` on the
    # global stream during `QuantumSystem` construction, silently shifting
    # downstream `rand`/`randn` results in user scripts that seed once at
    # the top.

    using Piccolo
    using Random
    using SparseArrays

    H_drift = sparse(ComplexF64[1.0 0.0; 0.0 -1.0])
    H1 = sparse(ComplexF64[0.0 1.0; 1.0 0.0])
    H2 = sparse(ComplexF64[0.0 -1.0im; 1.0im 0.0])

    # 1. Direct validator entry points.
    d_simple = NonlinearDrive(H1, u -> u[1] * u[2] * u[3])
    Random.seed!(42)
    baseline = rand(5)
    Random.seed!(42)
    validate_drive_jacobian(d_simple, 3)
    validate_drive_hessian(d_simple, 3)
    @test rand(5) == baseline

    # 2. Full QuantumSystem path with NonlinearDrives + globals — this is
    #    the surface that hit the user. Multiple nonlinear drives × a
    #    larger u_dim is the worst case for randn consumption.
    drives = AbstractDrive[
        LinearDrive(H1, 1),
        NonlinearDrive(H2, u -> u[1]^2 + u[2] * u[3]),
        NonlinearDrive(H1, u -> u[2] * u[3]),
    ]
    bounds = [(-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0)]
    Random.seed!(42)
    baseline2 = rand(5)
    Random.seed!(42)
    QuantumSystem(H_drift, drives, bounds; global_params = (g1 = 0.5,))
    @test rand(5) == baseline2

    # 3. Same again but with no globals — verifies the no-globals path is
    #    also RNG-clean.
    Random.seed!(42)
    baseline3 = rand(5)
    Random.seed!(42)
    QuantumSystem(H_drift, drives, bounds)
    @test rand(5) == baseline3
end

@testitem "G works on AbstractDrive types" begin
    using Piccolo
    using SparseArrays

    H = sparse(ComplexF64[0 1; 1 0])

    # G on LinearDrive should match G on the matrix
    ld = LinearDrive(H, 1)
    @test Piccolo.Isomorphisms.G(ld) == Piccolo.Isomorphisms.G(H)

    # G on NonlinearDrive
    nd = NonlinearDrive(H, u -> u[1]^2)
    @test Piccolo.Isomorphisms.G(nd) == Piccolo.Isomorphisms.G(H)

    # Broadcasting over Vector{AbstractDrive}
    drives = AbstractDrive[ld, nd]
    G_mats = Piccolo.Isomorphisms.G.(drives)
    @test length(G_mats) == 2
    @test all(G_mats .== Ref(Piccolo.Isomorphisms.G(H)))
end

@testitem "DriftTerm" begin
    using Piccolo
    using SparseArrays
    using ForwardDiff

    H = sparse(ComplexF64[1 0; 0 -1])

    # Default (identity modulation)
    dt = DriftTerm(H)
    @test dt.H === H
    @test dt.modulation(0.0) == 1.0
    @test dt.modulation(5.0) == 1.0
    @test dt.modulation_deriv(0.0) == 0.0
    @test dt.modulation_deriv(5.0) == 0.0

    # With modulation
    omega = 2.0
    dt2 = DriftTerm(H, t -> cos(omega * t))
    @test dt2.modulation(0.0) == 1.0
    @test dt2.modulation(pi / (2 * omega)) ≈ 0.0 atol = 1e-14
    @test dt2.modulation_deriv(0.0) ≈ 0.0 atol = 1e-14
    @test dt2.modulation_deriv(pi / (2 * omega)) ≈ -omega atol = 1e-10

    # Derivative matches ForwardDiff
    for t in [0.0, 0.5, 1.0, 2.3]
        @test dt2.modulation_deriv(t) ≈ ForwardDiff.derivative(dt2.modulation, t) atol =
            1e-10
    end

    # Non-real modulation should fail at construction
    @test_throws AssertionError DriftTerm(H, t -> im * t)
end

@testitem "ModulatedDrive" begin
    using Piccolo
    using SparseArrays
    using ForwardDiff

    H = sparse(ComplexF64[0 1; 1 0])
    omega = 3.0

    # Wrap a LinearDrive
    ld = LinearDrive(H, 2)
    md = ModulatedDrive(ld, t -> cos(omega * t))

    @test md.base === ld
    @test md.modulation(0.0) == 1.0
    @test md.modulation(pi / (2 * omega)) ≈ 0.0 atol = 1e-14

    # Derivative pre-computed correctly
    for t in [0.0, 0.5, 1.0]
        @test md.modulation_deriv(t) ≈ ForwardDiff.derivative(md.modulation, t) atol = 1e-10
    end

    # drive_matrix delegates to base
    @test drive_matrix(md) == drive_matrix(ld)
    @test drive_dim(md) == drive_dim(ld)
    @test active_controls(md) == active_controls(ld)

    # Isomorphisms.G delegates
    @test Piccolo.Isomorphisms.G(md) == Piccolo.Isomorphisms.G(ld)

    # Wrap a NonlinearDrive
    nd = NonlinearDrive(H, u -> u[1]^2; active_controls = [1])
    mnd = ModulatedDrive(nd, t -> sin(omega * t))

    @test drive_matrix(mnd) == drive_matrix(nd)
    @test active_controls(mnd) == [1]

    # has_nonlinear_drives detects wrapped NonlinearDrive
    drives = AbstractDrive[md, mnd]
    @test has_nonlinear_drives(drives)

    drives2 = AbstractDrive[md]
    @test !has_nonlinear_drives(drives2)
end

@testitem "drive_coeff with time argument" begin
    using Piccolo
    using SparseArrays

    H = sparse(ComplexF64[0 1; 1 0])
    omega = 2.0
    u = [0.3, 0.7, 0.1]

    # LinearDrive ignores t
    ld = LinearDrive(H, 2)
    @test drive_coeff(ld, u, 0.0) == 0.7
    @test drive_coeff(ld, u, 99.0) == 0.7
    @test drive_coeff_jac(ld, u, 0.0, 2) == 1.0
    @test drive_coeff_jac(ld, u, 0.0, 1) == 0.0

    # NonlinearDrive ignores t
    nd = NonlinearDrive(H, u -> u[1]^2 + u[2]^2)
    @test drive_coeff(nd, u, 0.0) ≈ 0.09 + 0.49
    @test drive_coeff(nd, u, 99.0) ≈ 0.09 + 0.49
    @test drive_coeff_jac(nd, u, 0.0, 1) ≈ 0.6

    # ModulatedDrive uses t
    md = ModulatedDrive(ld, t -> cos(omega * t))
    @test drive_coeff(md, u, 0.0) ≈ 0.7 * 1.0
    @test drive_coeff(md, u, pi / (2 * omega)) ≈ 0.0 atol = 1e-14
    @test drive_coeff_jac(md, u, 0.0, 2) ≈ 1.0 * 1.0
    @test drive_coeff_jac(md, u, pi / (2 * omega), 2) ≈ 0.0 atol = 1e-14

    # ModulatedDrive wrapping NonlinearDrive
    mnd = ModulatedDrive(nd, t -> cos(omega * t))
    @test drive_coeff(mnd, u, 0.0) ≈ 0.58
    t_test = 0.3
    @test drive_coeff(mnd, u, t_test) ≈ (u[1]^2 + u[2]^2) * cos(omega * t_test)
    @test drive_coeff_jac(mnd, u, t_test, 1) ≈ 2 * u[1] * cos(omega * t_test)

    # drive_coeff_dt
    @test drive_coeff_dt(ld, u, 0.0) == 0.0
    @test drive_coeff_dt(nd, u, 0.0) == 0.0
    @test drive_coeff_dt(md, u, 0.0) ≈ 0.0 atol = 1e-14  # cos'(0) = 0
    @test drive_coeff_dt(md, u, 0.5) ≈ u[2] * (-omega * sin(omega * 0.5)) atol = 1e-10

    # 2-arg backward-compatible shims
    @test drive_coeff(ld, u) == drive_coeff(ld, u, 0.0)
    @test drive_coeff_jac(ld, u, 2) == drive_coeff_jac(ld, u, 0.0, 2)
end

@testitem "coupling_drive: values, Jacobian, Hessian, sparsity" begin
    using Piccolo
    using SparseArrays

    H = sparse([0.0+0im 1.0+0im; 1.0+0im 0.0+0im])
    d = coupling_drive(H, 1, 3; strength = 0.85)

    u = [2.0, 5.0, 3.0, 7.0]
    @test drive_coeff(d, u) ≈ 0.85 * 2.0 * 3.0
    @test drive_coeff_jac(d, u, 1) ≈ 0.85 * 3.0
    @test drive_coeff_jac(d, u, 3) ≈ 0.85 * 2.0
    @test drive_coeff_jac(d, u, 2) == 0.0
    @test drive_coeff_jac(d, u, 4) == 0.0
    @test drive_coeff_hess(d, u, 1, 3) ≈ 0.85
    @test drive_coeff_hess(d, u, 3, 1) ≈ 0.85
    @test drive_coeff_hess(d, u, 1, 1) == 0.0
    @test drive_coeff_hess(d, u, 2, 3) == 0.0
    @test active_controls(d) == [1, 3]

    # default strength = 1.0
    d1 = coupling_drive(H, 2, 4)
    @test drive_coeff(d1, u) ≈ 5.0 * 7.0

    # i == j is rejected
    @test_throws ArgumentError coupling_drive(H, 2, 2)

    # functor parameterization: introspectable, concrete, serializable
    @test d.coeff isa BilinearCouplerCoeff
    @test (d.coeff.i, d.coeff.j, d.coeff.a) == (1, 3, 0.85)
    @test typeof(d) == typeof(d1)  # one concrete type for all couplings

    using Serialization
    buf = IOBuffer()
    serialize(buf, d)
    seekstart(buf)
    d2 = deserialize(buf)
    @test drive_coeff(d2, u) ≈ drive_coeff(d, u)
    @test d2.coeff.a == 0.85
end

@testitem "coupling_drive ≡ trilinear NonlinearDrive with frozen global" begin
    using Piccolo
    using SparseArrays

    # The baked form a·u[i]·u[j] must match the trilinear form
    # u[i]·u[j]·u[g] evaluated at a control vector with u[g] frozen to a.
    H = sparse([1.0+0im 0.0+0im; 0.0+0im -1.0+0im])
    i, j, g, a = 1, 2, 5, 0.7321

    baked = coupling_drive(H, i, j; strength = a)
    trilinear = NonlinearDrive(
        H,
        u -> u[i] * u[j] * u[g],
        (u, p) ->
            p == i ? u[j] * u[g] : p == j ? u[i] * u[g] : p == g ? u[i] * u[j] : 0.0;
        coeff_hess = (u, p, q) -> begin
            if (p == i && q == j) || (p == j && q == i)
                u[g]
            elseif (p == i && q == g) || (p == g && q == i)
                u[j]
            elseif (p == j && q == g) || (p == g && q == j)
                u[i]
            else
                0.0
            end
        end,
        active_controls = [i, j, g],
    )

    u = [1.3, -0.4, 9.9, 2.2, a]  # u[g] frozen at the baked strength
    @test drive_coeff(baked, u) ≈ drive_coeff(trilinear, u)
    for p in (i, j)
        @test drive_coeff_jac(baked, u, p) ≈ drive_coeff_jac(trilinear, u, p)
    end
    for p in (i, j), q in (i, j)
        @test drive_coeff_hess(baked, u, p, q) ≈ drive_coeff_hess(trilinear, u, p, q)
    end
end

@testitem "coupling_drive: AD cross-check + QuantumSystem integration" begin
    using Piccolo
    using SparseArrays

    H = sparse([0.0+0im 1.0+0im; 1.0+0im 0.0+0im])

    # Independent ForwardDiff cross-check of the analytic functors
    # (validate_* sample random u and compare against AD; throws on mismatch)
    for (i, j, a) in ((1, 3, 0.85), (2, 4, -1.3), (1, 2, 2))  # incl. negative + Int strength
        d = coupling_drive(H, i, j; strength = a)
        validate_drive_jacobian(d, 4)
        validate_drive_hessian(d, 4)
    end

    # End-to-end: coupling_drive inside a QuantumSystem
    σz = sparse([1.0+0im 0.0+0im; 0.0+0im -1.0+0im])
    drives = AbstractDrive[LinearDrive(σz, 1), coupling_drive(H, 1, 2; strength = 0.5)]
    sys = QuantumSystem(zeros(ComplexF64, 2, 2), drives, [(-1.0, 1.0), (-1.0, 1.0)])
    u = [0.3, -0.8]
    @test sys.H(u, 0.0) ≈ 0.3 * σz + 0.5 * 0.3 * (-0.8) * H
end
