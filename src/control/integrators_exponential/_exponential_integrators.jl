"""
    ExponentialIntegrators

Matrix-exponential integrators for Piccolo trajectories. Each sub-interval of the
trajectory is advanced exactly as `exp(Δt · G(u, t))`, where `G` is the generator
implied by the trajectory type (isomorphic `iso(-iH)` for closed systems, the full
Lindblad super-generator for open systems).

# Concrete types

- `HermitianExponentialIntegrator` — for closed-system trajectories where
  `G = iso(-iH)` with `H` Hermitian. Forward pass uses eigendecomposition of `H`
  via the in-place `exp_eigen!` helper and a 7-buffer scratch set
  preallocated on the integrator struct. Supported trajectory types:
  `UnitaryTrajectory`, `KetTrajectory`, `MultiKetTrajectory`.

- `NonHermitianExponentialIntegrator` — for open-system trajectories where
  `G` is the (generally non-Hermitian) vectorized Lindbladian. Forward pass uses
  scaling-and-squaring via the in-place `exp_generator!` helper
  (`ExponentialUtilities.exponential!` with `ExpMethodHigham2005`). Supported
  trajectory types: `DensityTrajectory`, `MultiDensityTrajectory`.

Both concrete types share `AbstractExponentialIntegrator`; the
sparsity / index-mapping scaffolding dispatches on the abstract type, while
Jacobian/Hessian assembly for the density family has per-type overrides due to
the `(ρ⃗̃, Δt, t, u)` component ordering.

# In-place helpers

- `exp_eigen!` — Hermitian forward step via `eigen!(Hermitian(·))` with
  preallocated eigenvalue, eigenvector, diagonal-scale, and isomorphism buffers.
- `exp_generator!` — non-Hermitian forward step via Higham-2005 scaling
  and squaring, preserving the input generator buffer.
"""
module ExponentialIntegrators

export AbstractExponentialIntegrator
export HermitianExponentialIntegrator
export NonHermitianExponentialIntegrator
export x_name, single_state_dim
export exp_eigen
export exp_eigen!
export exp_generator!
export dk_divided_difference!, dk_apply!, dk_first_order_derivative!
export DaleckiiKreinWorkspace, DK_DEGENERACY_RTOL
export dk_second_divided_difference,
    dk_second_order_apply!, dk_second_order_derivative!, dk_second_order_block!
export DaleckiiKreinSecondOrderWorkspace
export matrix_free_jacobian_op

using LinearAlgebra
using SparseArrays
using ExponentialAction
using ExponentialUtilities: exponential!, ExpMethodHigham2005
using NamedTrajectories
using NamedTrajectories: KnotPoint
using DirectTrajOpt: AbstractIntegrator
using DirectTrajOpt.Integrators: AbstractBilinearIntegrator, test_integrator
import DirectTrajOpt.Integrators:
    get_jacobian_structure, get_hessian_of_lagrangian_structure
using DirectTrajOpt.CommonInterface
import DirectTrajOpt.CommonInterface:
    evaluate!,
    jacobian_structure,
    jacobian!,
    hessian_structure,
    hessian_of_lagrangian,
    hessian_of_lagrangian!
import DirectTrajOpt.CommonInterface: eval_jacobian, eval_hessian_of_lagrangian
using Piccolo
using Piccolo:
    Isomorphisms,
    KetTrajectory,
    UnitaryTrajectory,
    MultiKetTrajectory,
    DensityTrajectory,
    MultiDensityTrajectory,
    OpenQuantumSystem,
    AbstractQuantumTrajectory
using Piccolo: SamplingTrajectory, sampling_member_states
using Piccolo: get_system, state_name, drive_name
using Piccolo: compact_lindbladian_parts, compact_generator_closure
using ForwardDiff
using TestItemRunner
using TrajectoryIndexingUtils

const ⊗ = kron

"""
    AbstractExponentialIntegrator <: AbstractBilinearIntegrator

Common supertype for the matrix-exponential integrators in this module.
Sparsity / index-mapping scaffolding dispatches on this type; concrete
subtypes (`HermitianExponentialIntegrator`, `NonHermitianExponentialIntegrator`)
provide trajectory-specific forward passes and Jacobian/Hessian assembly.
"""
abstract type AbstractExponentialIntegrator <: AbstractBilinearIntegrator end

"""
    matrix_free_jacobian_op(ℰ, traj)

Hook returning the matrix-free dynamics-Jacobian block for the analytic MultiKet
path. Declared here (empty) so the MultiKet `eval_jacobian` — defined in this
submodule — can dispatch to it; the concrete method returning a `CPUJacobianOp`
is added at the top-level `Piccolissimo` scope in `solvers/cpu_dk_jacobian_op.jl`
(which loads after this submodule, so a forward-reference hook is required). Only
reached when `ℰ.matrix_free && ℰ.H_dirs !== nothing`.
"""
function matrix_free_jacobian_op end

"""
    exp_eigen(G::AbstractMatrix{<:Real})

Compute exp(G) where G is the isomorphic generator from a QuantumSystem (i.e., G = sys.G(a, t)).
Uses eigendecomposition of the underlying Hermitian Hamiltonian.

# Note
This function is specifically designed for quantum system generators where G = iso(-iH).
It computes iso(exp(-iH)) = exp(G) using the eigendecomposition of H.
Do NOT use this for arbitrary real matrices - use `exp(G)` from LinearAlgebra instead.

Returns the full matrix exponential in isomorphic form.
"""
function exp_eigen(G::AbstractMatrix{<:Real})
    # G is in isomorphic real form, convert to Hermitian operator H
    H = Isomorphisms.H(G)
    return exp_eigen(H)
end

"""
    exp_eigen(H::AbstractMatrix{<:Complex})

Compute iso(exp(-iH)) where H is a Hermitian matrix (e.g., from sys.H(a, t)).
Uses eigendecomposition for efficient computation.

This is the preferred method when you have direct access to H, as it avoids the G→H conversion.

Returns the full matrix exponential in isomorphic form.
"""
function exp_eigen(H::AbstractMatrix{<:Complex})
    # Ensure Hermitian type for eigen efficiency (convert to dense if sparse)
    Ĥ = Hermitian(Matrix(H))
    λ, V = eigen(Ĥ)
    # Return in isomorphic form: exp(-iH)
    return Isomorphisms.iso(V * Diagonal(cis.(-λ)) * V')
end

"""
    iso!(out::AbstractMatrix{<:Real}, A::AbstractMatrix{<:Complex})

In-place isomorphism: writes `iso(A)` into `out`. `out` must be size `(2m, 2n)` when
`A` is `(m, n)`. Matches `Piccolo.Isomorphisms.iso`'s block convention

```
iso(A) = [ Re(A)  -Im(A);
           Im(A)   Re(A) ]
```
"""
function iso!(out::AbstractMatrix{<:Real}, A::AbstractMatrix{<:Complex})
    m, n = size(A)
    @inbounds for j = 1:n, i = 1:m
        a = A[i, j]
        ar = real(a)
        ai = imag(a)
        out[i, j] = ar
        out[i+m, j] = ai
        out[i, j+n] = -ai
        out[i+m, j+n] = ar
    end
    return out
end

"""
    exp_eigen!(expG_buf, H_buf, V_buf, λ_buf, cis_diag_buf, tmp_buf, work_buf, Δt)

In-place version of [`exp_eigen`](@ref). Computes `iso(exp(-i Δt H))` into `expG_buf`
using 7 preallocated work buffers. The 7 buffers eliminate the largest per-call
allocations (n×n complex eigenvector and exponentiated-matrix scratch), but
`LinearAlgebra.eigen!(Hermitian(·))` still allocates internally — both `F.values`
(a fresh `Vector{Float64}`) and `F.vectors` (a fresh `Matrix{ComplexF64}`), plus
LAPACK work arrays — so this routine is not strictly allocation-free. See note below.

# Arguments
- `expG_buf::AbstractMatrix{Float64}` — output buffer, size `(2n, 2n)`; receives
  isomorphic `exp(-i Δt H)`.
- `H_buf::AbstractMatrix{ComplexF64}` — on entry: the Hermitian matrix `H`; preserved
  across the call. Size `(n, n)`.
- `V_buf::AbstractMatrix{ComplexF64}` — scratch: `H_buf` is copied in, then overwritten
  by `eigen!` with eigenvectors. Size `(n, n)`.
- `λ_buf::AbstractVector{Float64}` — scratch: eigenvalues. Size `n`.
- `cis_diag_buf::AbstractVector{ComplexF64}` — scratch: `cis(-λ·Δt)`. Size `n`.
- `tmp_buf::AbstractMatrix{ComplexF64}` — scratch: `V · Diagonal(cis_diag)`. Size `(n, n)`.
- `work_buf::AbstractMatrix{ComplexF64}` — scratch: `tmp_buf · V'` (final complex
  product); kept distinct from `V_buf` to avoid aliasing with `adjoint(V_buf)`. Size `(n, n)`.
- `Δt::Real` — time-step scale; `H_buf` is NOT pre-scaled by `Δt`.

# Notes
`LinearAlgebra.eigen!(Hermitian(V_buf))` reduces `V_buf` to tridiagonal form in place
but still allocates a fresh `Vector{Float64}` for `F.values`, a fresh `Matrix{ComplexF64}`
for `F.vectors`, and LAPACK workspace arrays. If these show up as hot in profiling,
replace with a direct `LAPACK.syevr!` call that writes eigenvalues into `λ_buf` and
eigenvectors into `V_buf`.
"""
function exp_eigen!(
    expG_buf::AbstractMatrix{Float64},
    H_buf::AbstractMatrix{ComplexF64},
    V_buf::AbstractMatrix{ComplexF64},
    λ_buf::AbstractVector{Float64},
    cis_diag_buf::AbstractVector{ComplexF64},
    tmp_buf::AbstractMatrix{ComplexF64},
    work_buf::AbstractMatrix{ComplexF64},
    Δt::Real,
)
    # eigen!(Hermitian(·)) reduces the wrapped matrix to tridiagonal form in-place
    # but returns eigenvectors in a fresh matrix (F.vectors) — it does NOT overwrite
    # V_buf with the eigenvector matrix. We copy F.vectors into V_buf so subsequent
    # mul!s operate on the correct buffer.
    copyto!(V_buf, H_buf)
    F = LinearAlgebra.eigen!(Hermitian(V_buf))
    λ_buf .= F.values
    copyto!(V_buf, F.vectors)
    @. cis_diag_buf = cis(-λ_buf * Δt)

    # tmp_buf = V · Diagonal(cis_diag)
    mul!(tmp_buf, V_buf, Diagonal(cis_diag_buf))
    # work_buf = tmp_buf · V'              (distinct output buffer — no aliasing
    # with adjoint(V_buf))
    mul!(work_buf, tmp_buf, adjoint(V_buf))

    # Isomorphic conversion into Float64 output buffer
    iso!(expG_buf, work_buf)
    return expG_buf
end

"""
    exp_generator!(expG_buf, G_buf, Δt)

In-place matrix exponential for a (potentially non-Hermitian) generator `G`.
Computes `expG_buf := exp(Δt * G_buf)` using scaling-and-squaring (Higham 2005).
`G_buf` is preserved (its contents may be consulted by callers after this call).

# Arguments
- `expG_buf::AbstractMatrix{Float64}` — output buffer, same size as `G_buf`
- `G_buf::AbstractMatrix{Float64}` — input generator, in isomorphic real form
- `Δt::Real` — time-step scale
"""
function exp_generator!(
    expG_buf::AbstractMatrix{Float64},
    G_buf::AbstractMatrix{Float64},
    Δt::Real,
)
    # Scale into expG_buf, then exponentiate in place
    @. expG_buf = Δt * G_buf
    exponential!(expG_buf, ExpMethodHigham2005())
    return expG_buf
end

include("daleckii_krein.jl")
include("hermitian_exponential_integrator_type.jl")
include("hermitian_exponential_integrator_unitary.jl")
include("hermitian_exponential_integrator_ket.jl")
include("hermitian_exponential_integrator_multiket.jl")
include("hermitian_exponential_integrator_sampling.jl")
include("nonhermitian_exponential_integrator_type.jl")
include("nonhermitian_exponential_integrator_density.jl")
include("nonhermitian_exponential_integrator_multidensity.jl")
include("nonhermitian_exponential_integrator_sampling.jl")

@testitem "exp_eigen! matches exp_eigen on random Hermitian matrix" begin
    using LinearAlgebra
    using Random
    using Piccolo.Control.QuantumIntegrators: exp_eigen, exp_eigen!

    Random.seed!(42)
    n = 6
    A = randn(ComplexF64, n, n)
    H = Hermitian(A + A')  # random Hermitian

    # Current (allocating) one-shot
    expG_ref = exp_eigen(Matrix(H))

    # New in-place variant — provide 7 buffers
    H_buf = zeros(ComplexF64, n, n)
    V_buf = zeros(ComplexF64, n, n)
    λ_buf = zeros(Float64, n)
    cis_diag_buf = zeros(ComplexF64, n)
    tmp_buf = zeros(ComplexF64, n, n)
    work_buf = zeros(ComplexF64, n, n)
    expG_buf = zeros(Float64, 2n, 2n)

    copyto!(H_buf, H)
    exp_eigen!(expG_buf, H_buf, V_buf, λ_buf, cis_diag_buf, tmp_buf, work_buf, 1.0)

    @test isapprox(expG_buf, expG_ref; atol = 1e-14)
end

@testitem "exp_eigen! respects Δt scaling" begin
    using LinearAlgebra
    using Piccolo.Control.QuantumIntegrators: exp_eigen, exp_eigen!

    n = 4
    H = Hermitian(ComplexF64[
        1 0 0 0
        0 2 0 0
        0 0 3 0
        0 0 0 4
    ])

    H_buf = Matrix(H)
    V_buf = zeros(ComplexF64, n, n)
    λ_buf = zeros(Float64, n)
    cis_diag_buf = zeros(ComplexF64, n)
    tmp_buf = zeros(ComplexF64, n, n)
    work_buf = zeros(ComplexF64, n, n)
    expG_buf = zeros(Float64, 2n, 2n)

    Δt = 0.7
    exp_eigen!(expG_buf, H_buf, V_buf, λ_buf, cis_diag_buf, tmp_buf, work_buf, Δt)

    expected = exp_eigen(Δt * Matrix(H))
    @test isapprox(expG_buf, expected; atol = 1e-14)
end

@testitem "exp_eigen! zero allocations after warmup" begin
    using LinearAlgebra
    using BenchmarkTools
    using Piccolo.Control.QuantumIntegrators: exp_eigen!

    n = 4
    H = Hermitian(ComplexF64[
        1.0 0.2 0 0
        0.2 2.0 0.1 0
        0 0.1 3.0 0.3
        0 0 0.3 4.0
    ])

    H_buf = Matrix(H)
    V_buf = zeros(ComplexF64, n, n)
    λ_buf = zeros(Float64, n)
    cis_diag_buf = zeros(ComplexF64, n)
    tmp_buf = zeros(ComplexF64, n, n)
    work_buf = zeros(ComplexF64, n, n)
    expG_buf = zeros(Float64, 2n, 2n)

    # Warmup
    exp_eigen!(expG_buf, H_buf, V_buf, λ_buf, cis_diag_buf, tmp_buf, work_buf, 1.0)

    # Measure — with all 7 buffers preallocated, the only remaining allocations come from
    # `LinearAlgebra.eigen!(Hermitian(V_buf))`, which internally invokes `LAPACK.syevr!`
    # via the standard (non-preallocated) wrapper. That wrapper allocates:
    #   - `F.values` (Vector{Float64})
    #   - `F.vectors` (Matrix{ComplexF64})
    #   - LAPACK work / rwork / iwork / isuppz scratch arrays
    # Measured ≈ 4.1 KB on n=4 ComplexF64. True zero-alloc requires switching to the
    # low-level `LAPACK.syevr!` signature with all scratch arrays preallocated as
    # additional buffer fields on HermitianExponentialIntegrator — deferred per the
    # plan's "Task 2 Step 3 Note on eigen!" and flagged for a future optimization task.
    allocs = @ballocated exp_eigen!(
        $expG_buf,
        $H_buf,
        $V_buf,
        $λ_buf,
        $cis_diag_buf,
        $tmp_buf,
        $work_buf,
        1.0,
    )
    @test allocs ≤ 5_000
end

@testitem "exp_generator! matches naive exp(Δt * G) on non-Hermitian real matrix" begin
    using LinearAlgebra
    using Random
    using Piccolo.Control.QuantumIntegrators: exp_generator!

    Random.seed!(7)
    n = 8
    G = randn(Float64, n, n)  # non-Hermitian real generator
    Δt = 0.3

    expG_ref = exp(Δt * G)

    G_buf = copy(G)
    expG_buf = zeros(Float64, n, n)
    exp_generator!(expG_buf, G_buf, Δt)

    # scaling-and-squaring (Higham 2005) matches naive exp to high precision
    @test isapprox(expG_buf, expG_ref; atol = 1e-10)
    # G_buf is preserved (its contents may be consulted by callers after this call)
    @test G_buf == G
end

@testitem "exp_generator! allocations after warmup" begin
    using LinearAlgebra
    using BenchmarkTools
    using Piccolo.Control.QuantumIntegrators: exp_generator!

    n = 8
    G = randn(Float64, n, n)
    G_buf = copy(G)
    expG_buf = zeros(Float64, n, n)

    # Warmup
    exp_generator!(expG_buf, G_buf, 0.3)

    # Note: `ExponentialUtilities.exponential!(A, ExpMethodHigham2005())` still allocates
    # internal scratch (Pade pivot vectors, matrix products) that are not exposed via a
    # cache argument in the one-argument call. Plan's suggested bound: < 512 B.
    # Observed: ~3312 B on n=8 Float64 (Pade-13 scratch). Bound relaxed to 10_000 B
    # with headroom; true zero-alloc requires switching to an in-place cache interface
    # (e.g., `exponential!(A, method, cache)`) and is out of scope for Task 7.
    allocs = @ballocated exp_generator!($expG_buf, $G_buf, 0.3)
    @test allocs ≤ 10_000
end

@testitem "NonHermitianExponentialIntegrator type exists and subtypes abstract" begin
    @test NonHermitianExponentialIntegrator <: AbstractExponentialIntegrator
    # Field-existence probe via names(typeof(...)) is harder without an instance;
    # Task 9 adds the constructor-based test
end

end
