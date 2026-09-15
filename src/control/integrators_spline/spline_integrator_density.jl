# ============================================================================ #
# SplineIntegrator for DensityTrajectory
#
# Uses the compact density isomorphism: state x ∈ ℝ^{n²} (Hermiticity halves 2n²)
# via density_to_compact_iso / compact_iso_to_density
#
# For constraint evaluation, integrates the state forward directly (like Ket):
#   ẋ = ℒ(u) · x   where x = compact_iso(ρ), via column-wise complex matrix ops
#   achieving O(n³K) per knot point (K = number of Lindblad jump operators,
#   the `Ks`, not the knot count).
#
# For Jacobian/Hessian, uses compact Lindbladian generators 𝒢c in the
# sensitivity-augmented ODE to compute the full n²×n² propagator and its
# parameter derivatives simultaneously.
#
# The Lindblad dynamics are:
#   ρ̇ = -i[H(u),ρ] + Σⱼ (Lⱼ ρ Lⱼ† - ½{Lⱼ†Lⱼ, ρ})
#
# In compact isomorphism: ẋ = P 𝒢(u) L x = 𝒢c(u) x
# where L = density_lift_matrix, P = density_projection_matrix, PL = I_{n²}
# ============================================================================ #

# ============================================================================ #
# Lindbladian RHS: Apply ℒ[M] via complex n×n matrix products — O(n³)
# ============================================================================ #

"""
    lindblad_apply!(dM, M, H_eff, Δt, Ls, Ks, tmp)

In-place Lindbladian: dM = Δt * (-i[H_eff, M] + Σⱼ(Lⱼ M Lⱼ† - ½{Kⱼ, M}))

# Arguments
- `dM::Matrix{ComplexF64}`: Output (n×n)
- `M::Matrix{ComplexF64}`: Input density/perturbation (n×n)
- `H_eff::Matrix{ComplexF64}`: Effective Hamiltonian H_drift + Σ uᵢ H_drives[i]
- `Δt::Real`: Time scaling factor
- `Ls::Vector`: Lindblad operators
- `Ks::Vector`: Precomputed Lⱼ†Lⱼ
- `tmp::Matrix{ComplexF64}`: Scratch buffer (n×n)
"""
function lindblad_apply!(
    dM::AbstractMatrix{<:Complex},
    M::AbstractMatrix{<:Complex},
    H_eff::AbstractMatrix{<:Complex},
    Δt::Real,
    Ls::Vector,
    Ks::Vector,
    tmp::AbstractMatrix{<:Complex},
)
    # dM = Δt * (-i * H_eff * M + i * M * H_eff)
    mul!(dM, H_eff, M, -im * Δt, false)   # dM = -iΔt H M
    mul!(dM, M, H_eff, im * Δt, true)     # dM += iΔt M H

    # dM += Δt * Σⱼ (Lⱼ M Lⱼ† - ½ Kⱼ M - ½ M Kⱼ)
    _dissipator_accumulate!(dM, M, Δt, Ls, Ks, tmp)

    return nothing
end

# ---------------------------------------------------------------------------
# Dissipator-only kernels.
#
# `lindblad_apply!` / `lindblad_adjoint_apply!` delegate their dissipator loop
# here, so there is exactly ONE implementation of it. The rate-carrying
# overloads add a per-jump multiplier `rates[j]` — the SAME convention as the
# rate-absorbed form (`Lⱼ = √γⱼ Aⱼ`, `Kⱼ = Lⱼ†Lⱼ`), just factored out so a
# control-dependent rate can vary per substep without rebuilding the operator
# products. `rates = fill(1, nL)` reproduces the rate-absorbed path exactly.
#
# Split out for the Lindblad product core (Piccolissimo#313), whose Strang
# substep needs `exp(θ𝒟)` on its own, with the Hamiltonian handled by the
# Chebyshev coherent factor.
# ---------------------------------------------------------------------------

@inline function _dissipator_accumulate!(
    dM::AbstractMatrix{<:Complex},
    M::AbstractMatrix{<:Complex},
    Δt::Real,
    Ls::Vector,
    Ks::Vector,
    tmp::AbstractMatrix{<:Complex},
)
    @inbounds for (L, K) in zip(Ls, Ks)
        mul!(tmp, L, M)                     # tmp = Lⱼ M
        mul!(dM, tmp, L', Δt, true)         # dM += Δt * tmp * Lⱼ†
        mul!(dM, K, M, -Δt / 2, true)       # dM -= Δt/2 * Kⱼ M
        mul!(dM, M, K, -Δt / 2, true)       # dM -= Δt/2 * M Kⱼ
    end
    return nothing
end

@inline function _dissipator_accumulate!(
    dM::AbstractMatrix{<:Complex},
    M::AbstractMatrix{<:Complex},
    Δt::Real,
    Ls::Vector,
    Ks::Vector,
    rates::AbstractVector{<:Real},
    tmp::AbstractMatrix{<:Complex},
)
    @inbounds for j in eachindex(Ls, Ks, rates)
        L, K, s = Ls[j], Ks[j], Δt * rates[j]
        mul!(tmp, L, M)
        mul!(dM, tmp, L', s, true)
        mul!(dM, K, M, -s / 2, true)
        mul!(dM, M, K, -s / 2, true)
    end
    return nothing
end

@inline function _dissipator_adjoint_accumulate!(
    dA::AbstractMatrix{<:Complex},
    A::AbstractMatrix{<:Complex},
    Δt::Real,
    Ls::Vector,
    Ks::Vector,
    tmp::AbstractMatrix{<:Complex},
)
    @inbounds for (L, K) in zip(Ls, Ks)
        mul!(tmp, L', A)                    # tmp = Lⱼ† A
        mul!(dA, tmp, L, Δt, true)          # dA += Δt * tmp * Lⱼ
        mul!(dA, K', A, -Δt / 2, true)      # dA -= Δt/2 * Kⱼ† A
        mul!(dA, A, K', -Δt / 2, true)      # dA -= Δt/2 * A Kⱼ†
    end
    return nothing
end

@inline function _dissipator_adjoint_accumulate!(
    dA::AbstractMatrix{<:Complex},
    A::AbstractMatrix{<:Complex},
    Δt::Real,
    Ls::Vector,
    Ks::Vector,
    rates::AbstractVector{<:Real},
    tmp::AbstractMatrix{<:Complex},
)
    @inbounds for j in eachindex(Ls, Ks, rates)
        L, K, s = Ls[j], Ks[j], Δt * rates[j]
        mul!(tmp, L', A)
        mul!(dA, tmp, L, s, true)
        mul!(dA, K', A, -s / 2, true)
        mul!(dA, A, K', -s / 2, true)
    end
    return nothing
end

"""
    dissipator_apply!(dM, M, Δt, Ls, Ks, rates, tmp)

`dM = Δt · Σⱼ rates[j] · (Lⱼ M Lⱼ† − ½{Kⱼ, M})` — the dissipative half of
`lindblad_apply!`, with no Hamiltonian term and a per-jump rate
multiplier. First argument mutated, caller-owned scratch, zero allocation.
`dM` must not alias `M`.
"""
function dissipator_apply!(
    dM::AbstractMatrix{<:Complex},
    M::AbstractMatrix{<:Complex},
    Δt::Real,
    Ls::Vector,
    Ks::Vector,
    rates::AbstractVector{<:Real},
    tmp::AbstractMatrix{<:Complex},
)
    fill!(dM, 0)
    _dissipator_accumulate!(dM, M, Δt, Ls, Ks, rates, tmp)
    return nothing
end

"""
    dissipator_adjoint_apply!(dA, A, Δt, Ls, Ks, rates, tmp)

Hilbert–Schmidt adjoint of `dissipator_apply!`:
`dA = Δt · Σⱼ rates[j] · (Lⱼ† A Lⱼ − ½{Kⱼ†, A})`. Same contract; `dA` must not
alias `A`.
"""
function dissipator_adjoint_apply!(
    dA::AbstractMatrix{<:Complex},
    A::AbstractMatrix{<:Complex},
    Δt::Real,
    Ls::Vector,
    Ks::Vector,
    rates::AbstractVector{<:Real},
    tmp::AbstractMatrix{<:Complex},
)
    fill!(dA, 0)
    _dissipator_adjoint_accumulate!(dA, A, Δt, Ls, Ks, rates, tmp)
    return nothing
end

"""
    lindblad_adjoint_apply!(dA, A, H_eff, Δt, Ls, Ks, tmp)

In-place Hilbert–Schmidt adjoint of `lindblad_apply!` — the Heisenberg
generator:

    dA = Δt * (+i[H_eff, A] + Σⱼ(Lⱼ† A Lⱼ - ½{Kⱼ, A}))

Defined by `⟨A, ℒ[M]⟩ = ⟨ℒ†[A], M⟩` under the Hilbert–Schmidt inner product
`⟨X, Y⟩ = tr(X' * Y)`, where `ℒ` is exactly what `lindblad_apply!` computes for
the *same* `(H_eff, Δt, Ls, Ks)`.

Same calling contract as the forward action: the dissipator products `Ks` are
run-invariant and supplied by the caller, the scratch buffer is caller-owned,
the function holds no state, and the first argument is mutated in place. Same
rate convention too — rates are absorbed into `Ls` (`Lⱼ = √γⱼ Aⱼ`) with
`Kⱼ = Lⱼ†Lⱼ`; no separate rate vector is introduced.

`H_eff'` and `Kⱼ'` (rather than `H_eff` and `Kⱼ`) appear below deliberately:
they make the identity exact to machine precision for *any* inputs, including
an `H_eff` or a floating-point `Kⱼ = Lⱼ'Lⱼ` that is not bitwise Hermitian. For
the Hermitian inputs the integrators actually pass, the two forms coincide.

# Arguments
- `dA::AbstractMatrix{<:Complex}`: Output (n×n), overwritten
- `A::AbstractMatrix{<:Complex}`: Input observable / costate (n×n)
- `H_eff::AbstractMatrix{<:Complex}`: Effective Hamiltonian H_drift + Σ uᵢ H_drives[i]
- `Δt::Real`: Time scaling factor
- `Ls::Vector`: Lindblad operators
- `Ks::Vector`: Precomputed Lⱼ†Lⱼ
- `tmp::AbstractMatrix{<:Complex}`: Scratch buffer (n×n)
"""
function lindblad_adjoint_apply!(
    dA::AbstractMatrix{<:Complex},
    A::AbstractMatrix{<:Complex},
    H_eff::AbstractMatrix{<:Complex},
    Δt::Real,
    Ls::Vector,
    Ks::Vector,
    tmp::AbstractMatrix{<:Complex},
)
    # dA = Δt * (i * H† * A - i * A * H†)
    mul!(dA, H_eff', A, im * Δt, false)   # dA = iΔt H† A
    mul!(dA, A, H_eff', -im * Δt, true)   # dA -= iΔt A H†

    # dA += Δt * Σⱼ (Lⱼ† A Lⱼ - ½ Kⱼ† A - ½ A Kⱼ†)
    _dissipator_adjoint_accumulate!(dA, A, Δt, Ls, Ks, tmp)

    return nothing
end

# ============================================================================ #
# In-place compact iso ↔ density conversions (allocation-free)
#
# Layout: [Re upper triangle col-major; Im strict upper triangle col-major]
# Matches Isomorphisms.compact_iso_to_density / density_to_compact_iso exactly.
# ============================================================================ #

"""
    compact_iso_to_density!(M, x, n)

In-place: fill n×n Hermitian complex matrix `M` from compact iso vector `x` (length n²).
"""
function compact_iso_to_density!(
    M::AbstractMatrix{<:Complex},
    x::AbstractVector{<:Real},
    n::Int,
)
    idx = 0
    @inbounds for k = 1:n, j = 1:k
        idx += 1
        M[j, k] = x[idx]
        if j != k
            M[k, j] = x[idx]
        end
    end
    @inbounds for k = 2:n, j = 1:(k-1)
        idx += 1
        M[j, k] += im * x[idx]
        M[k, j] -= im * x[idx]
    end
    return nothing
end

"""
    density_to_compact_iso!(x, M, n)

In-place: fill compact iso vector `x` (length n²) from n×n Hermitian complex matrix `M`.
"""
function density_to_compact_iso!(
    x::AbstractVector{<:Real},
    M::AbstractMatrix{<:Complex},
    n::Int,
)
    idx = 0
    @inbounds for k = 1:n, j = 1:k
        idx += 1
        x[idx] = real(M[j, k])
    end
    @inbounds for k = 2:n, j = 1:(k-1)
        idx += 1
        x[idx] = imag(M[j, k])
    end
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────── #
# The Hilbert–Schmidt metric of the packed-triangular packing
# ─────────────────────────────────────────────────────────────────────────── #

"""
    compact_iso_hs_weights(n) -> Vector{Float64}

The diagonal Hilbert–Schmidt metric `W` of the compact (packed-triangular)
density isomorphism: `⟨A,B⟩ = tr(A'B) = aᵀ W b` for Hermitian `A`, `B` with
`a = pack(A)`, `b = pack(B)`.

`1.0` on the `n` real-diagonal slots, `2.0` on every real and imaginary
strict-upper slot — because a Hermitian matrix's off-diagonal entry is stored
once and counted twice. Ordered to match `density_to_compact_iso!`
exactly.

This is the factor that makes the compact-iso adjoint differ from a plain
transpose; see this file's header and `adjoint_compact_lindblad!`.
"""
function compact_iso_hs_weights(n::Int)
    w = Vector{Float64}(undef, n^2)
    idx = 0
    @inbounds for k = 1:n, j = 1:k
        idx += 1
        w[idx] = (j == k) ? 1.0 : 2.0
    end
    @inbounds for k = 2:n, j = 1:(k-1)
        idx += 1
        w[idx] = 2.0
    end
    @assert idx == n^2
    return w
end

# ─────────────────────────────────────────────────────────────────────────── #
# The cell's construction-time tape + alg_data container (slice 3b: moved to
# Piccolo with the density cell — the tape is plain construction-time data over
# Piccolo types; its matrix-free CONSUMERS stay in Piccolissimo).
# ─────────────────────────────────────────────────────────────────────────── #

# ─────────────────────────────────────────────────────────────────────────── #
# The cell's construction-time tape
# ─────────────────────────────────────────────────────────────────────────── #

"""
    LindbladDuhamelTape

Everything the density inner kernels need about the generator, baked once at
construction. Read-only and therefore shareable across the driver's tasks; every
mutable buffer lives in [`DensityDuhamelScratch`](@ref) instead.

`nb` is the number of density blocks the cell carries — `1` for
`DensityTrajectory`, `K` for `MultiDensityTrajectory`. It is the COMPONENT count,
never the core width: the core width is `nb·n` (that many `n`-dim columns the
two-sided generator acts on at once), and the two are declared separately per
ADR-0009.

Dissipators carry no `u` dependence on this path (see the file header), so `Ls`
and `Ks` are constants and there is no rate hook.
"""
struct LindbladDuhamelTape{DR,DA,HT,LT,KT}
    n::Int                              # Hilbert-space dimension
    n²::Int                             # compact-iso length per density block
    nb::Int                             # density blocks (1 for Density, K for MultiDensity)
    u_dim::Int                          # controls + globals
    order::Int                          # spline order: 1 or 3
    n_sub::Int                          # fixed RK4 substeps per knot interval
    n_terms::Int                        # drive terms (≥ n_drives with NonlinearDrive)
    H_drift::Matrix{ComplexF64}
    H_drives::Vector{Matrix{ComplexF64}}
    H_drives_t::HT                      # the same matrices as a TUPLE (see below)
    drives::DR                          # a TUPLE, never a Vector{AbstractDrive}
    drive_active::DA                    # active_controls per term as a TUPLE ([] ⇒ all)
    Ls::LT
    Ks::KT
    hs_weight::Vector{Float64}          # the W of `compact_iso_hs_weights`
end

"""
    duhamel_n_sub(tol) -> Int

Substep count for the fixed-step RK4 Duhamel quadrature, derived from the
integrator tolerance rather than hard-coded: RK4 is globally 4th order, so
`h = 1/n_sub` with `n_sub ≈ tol^(-1/4)` puts the sensitivity's own
discretization error at the same order as the residual solver's tolerance.
Clamped to `[8, 64]`.
"""
duhamel_n_sub(tol::Real) = clamp(ceil(Int, tol^(-1 / 4)), 8, 64)

function LindbladDuhamelTape(
    n::Int,
    nb::Int,
    u_dim::Int,
    order::Int,
    tol::Real,
    H_drift,
    H_drives,
    drives,
    Ls,
    Ks,
)
    return LindbladDuhamelTape(
        n,
        n^2,
        nb,
        u_dim,
        order,
        duhamel_n_sub(tol),
        length(drives),
        Matrix{ComplexF64}(H_drift),
        [Matrix{ComplexF64}(Hd) for Hd in H_drives],
        Tuple(Matrix{ComplexF64}(Hd) for Hd in H_drives),
        # THE DRIVES ARE A TUPLE, NOT A `Vector{AbstractDrive}`, AND THAT IS AN
        # ALLOCATION CONTRACT, NOT A STYLE CHOICE. `sys.H_drives` has abstract
        # element type, so `drive_coeff(drives[t], u)` is a DYNAMIC dispatch whose
        # `Float64` return is heap-boxed: measured at 16 B per call, 4 calls per RK4
        # stage, 128 stages per knot = 8 192 B PER KNOT — a linear-in-`N` term, and
        # exactly the kind of thing a kernel-scope gate cannot see. A tuple makes
        # `map(d -> drive_coeff(d, u), tape.drives)` fully unrolled and STATICALLY
        # dispatched, which is why every drive-coefficient call below goes through
        # `map` over a tuple rather than an indexed loop.
        Tuple(drives),
        Tuple(collect(active_controls(d)) for d in drives),
        Ls,
        Ks,
        compact_iso_hs_weights(n),
    )
end

"""
    DensityLindbladData{R}

Algorithm data for the density cells: `Tsit5Data`'s fields plus the cell's
[`LindbladDuhamelTape`](@ref). A distinct `D` type parameter rather than a new
field on `SplineIntegrator`, so nothing shared changes and the tape travels with
the integrator instead of living in identity-keyed global state (the hazard class
#307 had to fix and ADR-0009 rules out).
"""
struct DensityLindbladData{R<:Number,T<:LindbladDuhamelTape}
    Φ_probs::Vector{ODEProblem}
    Φ_structure::SparseMatrixCSC{Float64,Int}
    jvp_probs::Nothing
    vjp_probs::Nothing
    hvp_fwd_probs::Nothing
    hvp_bwd_probs::Nothing
    tape::T
end

DensityLindbladData{R}(d::Tsit5Data, tape::T) where {R<:Number,T<:LindbladDuhamelTape} =
    DensityLindbladData{R,T}(
        d.Φ_probs,
        d.Φ_structure,
        nothing,
        nothing,
        nothing,
        nothing,
        tape,
    )

"""
    duhamel_tape(𝒮) -> LindbladDuhamelTape

The density cell's generator tape. Errors with a pointer to the constructor for a
density integrator built before #346 (i.e. one whose `alg_data` is a bare
`Tsit5Data`).
"""
duhamel_tape(𝒮::SplineIntegrator) =
    𝒮.alg_data isa DensityLindbladData ? 𝒮.alg_data.tape :
    error(
        "no LindbladDuhamelTape on this integrator ($(typeof(𝒮.alg_data))). The " *
        "matrix-free density products need the construction-time generator tape — " *
        "rebuild the integrator with the current `SplineIntegrator(::DensityTrajectory, …)`.",
    )


# ============================================================================ #
# Constructor for DensityTrajectory
# ============================================================================ #

"""
    SplineIntegrator(qtraj::DensityTrajectory, N::Int; kwargs...)

Construct a SplineIntegrator for open quantum system density matrix dynamics.

Uses the hybrid Lindblad approach:
- Column-wise complex matrix operations for O(n³) per-column RHS
- Compact n²-dimensional state (exploiting Hermiticity)
- Per-knot buffer allocation for thread-safe parallel evaluation
- State-forward integration for efficient constraint evaluation

# Arguments
- `qtraj::DensityTrajectory`: The density matrix quantum trajectory
- `N::Int`: Number of knot points for discretization

# Keyword Arguments
- `spline_order::Union{Int,Nothing}=nothing`: Spline interpolation order (1=linear, 3=cubic);
  ignored when the pulse type already fixes the order (`LinearSplinePulse`/`CubicSplinePulse`).
- `alg::IntegrationAlgorithm=Tsit5Alg()`: Integration algorithm. `MagnusGL4Alg` is not
  supported for density trajectories (the constructor errors).
- `global_names::Union{Vector{Symbol},Nothing}=nothing`: Names of global optimization variables
"""
function SplineIntegrator(
    qtraj::DensityTrajectory,
    N::Int;
    spline_order::Union{Int,Nothing} = nothing,
    alg::IntegrationAlgorithm = Tsit5Alg(),
    global_names::Union{Vector{Symbol},Nothing} = nothing,
)
    if alg isa MagnusGL4Alg
        error("MagnusGL4Alg currently only supports UnitaryTrajectory")
    end
    if alg isa ChebyshevAlg
        error(
            "ChebyshevAlg is Hermitian-generator only (real spectrum): Lindbladian " *
            "density dynamics stay on the Tsit5 augmented path (ADR-0003); an analytic " *
            "non-Hermitian exp-action is a separate future track.",
        )
    end
    traj = NamedTrajectory(qtraj, N)
    return SplineIntegrator(
        qtraj,
        traj;
        spline_order = spline_order,
        alg = alg,
        global_names = global_names,
    )
end

function SplineIntegrator(
    qtraj::DensityTrajectory,
    traj::NamedTrajectory;
    spline_order::Union{Int,Nothing} = nothing,
    alg::IntegrationAlgorithm = Tsit5Alg(),
    global_names::Union{Vector{Symbol},Nothing} = nothing,
)
    if alg isa MagnusGL4Alg
        error("MagnusGL4Alg currently only supports UnitaryTrajectory")
    end
    if alg isa ChebyshevAlg
        error(
            "ChebyshevAlg is Hermitian-generator only (real spectrum): Lindbladian " *
            "density dynamics stay on the Tsit5 augmented path (ADR-0003); an analytic " *
            "non-Hermitian exp-action is a separate future track.",
        )
    end
    # Infer spline type from pulse if possible
    pulse = get_pulse(qtraj)
    if pulse isa LinearSplinePulse
        S = LinearSpline
        order = 1
    elseif pulse isa CubicSplinePulse
        S = CubicSpline
        order = 3
    else
        order = isnothing(spline_order) ? 1 : spline_order
        S = order == 1 ? LinearSpline : CubicSpline
    end

    sys = get_system(qtraj)::OpenQuantumSystem
    x = state_name(qtraj)  # :ρ⃗̃
    u = drive_name(qtraj)

    n = sys.levels  # Hilbert space dimension
    n² = n^2        # compact state dimension
    tol = alg.tol

    # Auto-detect globals if not specified
    if isnothing(global_names)
        global_names = Symbol[]
    end

    @assert traj.N > 1 "Trajectory must have at least two timesteps."
    @assert x in traj.names "State name $x must be in trajectory"

    control_dim = traj.dims[u]

    # Determine global dimensions
    if !isempty(global_names)
        global_comps = vcat([traj.global_components[name] for name in global_names]...)
        global_dim = length(global_comps)
    else
        global_dim = 0
    end

    # u_dim includes controls and globals
    u_dim = control_dim + global_dim

    # ================================================================
    # Extract Hamiltonian and dissipation components
    # ================================================================
    H_drift = collect(sys.H_drift)  # dense for mul!
    drives = sys.H_drives  # AbstractDrive objects (may include NonlinearDrive)
    n_terms = length(drives)  # number of drive terms (≥ n_drives for NonlinearDrive)
    H_drives = [collect(drive_matrix(d)) for d in drives]
    Ls = [collect(L) for L in sys.dissipation_operators]
    Ks = [sparse(L' * L) for L in sys.dissipation_operators]

    # Build control-to-drives mapping for NonlinearDrive support
    control_to_drives = [Int[] for _ = 1:u_dim]
    for (t_idx, d) in enumerate(drives)
        ac = active_controls(d)
        if isempty(ac)
            for j = 1:u_dim
                push!(control_to_drives[j], t_idx)
            end
        else
            for j in ac
                push!(control_to_drives[j], t_idx)
            end
        end
    end

    # Compact Lindbladian generators (for sensitivity ODE).
    # `compact_lindbladian_parts` returns 3-tuple: Hamiltonian-only drift,
    # per-drive factors, per-dissipator factors. This integrator's sensitivity
    # ODE treats `𝒢c_drift` as the "full" drift (Hamiltonian + dissipators),
    # so reassemble by summing dissipator factors into the drift. Dissipator
    # rates are *not* u-dependent on this path (only LinearDissipator is fully
    # supported); a NonlinearDissipator whose rate reads from u would need to
    # move into the per-drive sensitivity loop below — tracked as follow-up
    # in Task 11 of plan-20260419-110000-open-system-globals.
    𝒢c_drift_ham, 𝒢c_drives, 𝒢c_dissipators = compact_lindbladian_parts(sys)
    𝒢c_drift = isempty(𝒢c_dissipators) ? 𝒢c_drift_ham : 𝒢c_drift_ham + sum(𝒢c_dissipators)

    # ================================================================
    # Dimensions
    # ================================================================
    x_dim = traj.dims[x]  # should be n² (compact iso)
    @assert x_dim == n² "State dimension ($x_dim) must equal n² = $n²"

    N_traj = traj.N
    dim = x_dim * (N_traj - 1)

    if order == 1
        p_dim = 2 * u_dim
    elseif order == 3
        p_dim = 4 * u_dim
    end

    # ================================================================
    # Initial ODE parameter vector
    # ================================================================
    Φc_init = vec(Matrix{Float64}(I, n², n²))  # n⁴ entries (propagator)

    u_control_init = traj.bounds[u][2]

    if global_dim > 0
        g₀ = vcat(
            [traj.global_data[traj.global_components[name]] for name in global_names]...,
        )
        u_init = [u_control_init; g₀]
    else
        u_init = u_control_init
    end

    if order == 1
        p₀ = [u_init; u_init]
    elseif order == 3
        du = Symbol("d", u)
        if haskey(traj.bounds, du)
            du_control_init = traj.bounds[du][2]
        else
            du_control_init = zeros(control_dim)
        end
        du_init = [du_control_init; zeros(global_dim)]
        p₀ = [u_init; u_init; du_init; du_init]
    end
    Δt₀ = 1.0
    t₀ = 0.0

    # ================================================================
    # Per-knot ODE closures with dedicated buffers (thread-safe)
    #
    # Each Φ_probs[k] gets its own f! closure capturing its own set of
    # n×n complex buffers. This ensures that Threads.@threads in
    # evaluate! and eval_jacobian cannot cause data races.
    #
    # The f! handles variable-length state vectors:
    # - Propagator (sparsity calc): length n⁴ → n² columns
    # - State-forward (constraint eval): length n² → 1 column
    # ================================================================
    Φ_probs = Vector{ODEProblem}(undef, N_traj - 1)

    for kk = 1:(N_traj-1)
        # Per-knot buffers
        H_eff_kk = Matrix{ComplexF64}(undef, n, n)
        M_kk = Matrix{ComplexF64}(undef, n, n)
        dM_kk = Matrix{ComplexF64}(undef, n, n)
        tmp_kk = Matrix{ComplexF64}(undef, n, n)
        u_interp_kk = zeros(u_dim)  # interpolated controls for drive_coeff

        f_kk! = if !isempty(H_drives) && order == 1

            (dx, xₖ, p, τ) -> begin
                uₖ = @view p[1:u_dim]
                uₖ₊₁ = @view p[(u_dim+1):2u_dim]
                Δtₖ = p[end-1]

                # Interpolate controls
                @inbounds for i = 1:u_dim
                    u_interp_kk[i] = (1 - τ) * uₖ[i] + τ * uₖ₊₁[i]
                end

                # Build effective Hamiltonian via drive coefficients
                @. H_eff_kk = H_drift
                @inbounds for t_idx = 1:n_terms
                    c = drive_coeff(drives[t_idx], u_interp_kk)
                    @. H_eff_kk += c * H_drives[t_idx]
                end

                # Apply Lindbladian column-by-column
                n_cols = length(xₖ) ÷ n²
                @inbounds for j = 1:n_cols
                    col = @view xₖ[((j-1)*n²+1):(j*n²)]
                    d_col = @view dx[((j-1)*n²+1):(j*n²)]
                    compact_iso_to_density!(M_kk, col, n)
                    lindblad_apply!(dM_kk, M_kk, H_eff_kk, Δtₖ, Ls, Ks, tmp_kk)
                    density_to_compact_iso!(d_col, dM_kk, n)
                end
                return nothing
            end

        elseif !isempty(H_drives) && order == 3

            (dx, xₖ, p, τ) -> begin
                uₖ = @view p[1:u_dim]
                uₖ₊₁ = @view p[(u_dim+1):2u_dim]
                duₖ = @view p[(2u_dim+1):3u_dim]
                duₖ₊₁ = @view p[(3u_dim+1):4u_dim]
                Δtₖ = p[end-1]

                # Cubic Hermite basis functions
                τ2 = τ * τ
                τ3 = τ2 * τ
                h00 = 2τ3 - 3τ2 + 1
                h10 = (τ3 - 2τ2 + τ) * Δtₖ
                h01 = -2τ3 + 3τ2
                h11 = (τ3 - τ2) * Δtₖ

                # Interpolate controls
                @inbounds for i = 1:u_dim
                    u_interp_kk[i] =
                        h00 * uₖ[i] + h10 * duₖ[i] + h01 * uₖ₊₁[i] + h11 * duₖ₊₁[i]
                end

                # Build effective Hamiltonian via drive coefficients
                @. H_eff_kk = H_drift
                @inbounds for t_idx = 1:n_terms
                    c = drive_coeff(drives[t_idx], u_interp_kk)
                    @. H_eff_kk += c * H_drives[t_idx]
                end

                # Apply Lindbladian column-by-column
                n_cols = length(xₖ) ÷ n²
                @inbounds for j = 1:n_cols
                    col = @view xₖ[((j-1)*n²+1):(j*n²)]
                    d_col = @view dx[((j-1)*n²+1):(j*n²)]
                    compact_iso_to_density!(M_kk, col, n)
                    lindblad_apply!(dM_kk, M_kk, H_eff_kk, Δtₖ, Ls, Ks, tmp_kk)
                    density_to_compact_iso!(d_col, dM_kk, n)
                end
                return nothing
            end

        else
            # No drives: just drift Lindbladian
            (dx, xₖ, p, τ) -> begin
                Δtₖ = p[end-1]

                H_eff_kk .= H_drift

                n_cols = length(xₖ) ÷ n²
                @inbounds for j = 1:n_cols
                    col = @view xₖ[((j-1)*n²+1):(j*n²)]
                    d_col = @view dx[((j-1)*n²+1):(j*n²)]
                    compact_iso_to_density!(M_kk, col, n)
                    lindblad_apply!(dM_kk, M_kk, H_eff_kk, Δtₖ, Ls, Ks, tmp_kk)
                    density_to_compact_iso!(d_col, dM_kk, n)
                end
                return nothing
            end
        end

        Φ_probs[kk] = ODEProblem(f_kk!, Φc_init, (0.0, 1.0), [p₀; Δt₀; t₀])
    end

    # ================================================================
    # Sparsity structure and PropagatorResults
    # ================================================================
    calculate_sparsity = alg isa Tsit5Alg && !alg.adaptive
    if calculate_sparsity
        Φ_sol = solve(Φ_probs[1], Tsit5(); abstol = tol, reltol = tol, saveat = 1.0).u[end]
        Φ_mat = reshape(Φ_sol, n², n²)
        Φ_structure = sparse(Φ_mat)
    else
        Φ_structure = sparse(ones(n², n²))
    end

    # The generator tape for the matrix-free analytic Lindbladian Duhamel
    # sensitivity (#346). Carried on a density-specific `alg_data` type rather than
    # a new `SplineIntegrator` field, so nothing shared changes and the tape travels
    # with the integrator instead of living in identity-keyed global state.
    # `nb = 1`: one density block.
    alg_data = DensityLindbladData{Float64}(
        Tsit5Data{Float64}(Φ_probs, Φ_structure),
        LindbladDuhamelTape(n, 1, u_dim, order, tol, H_drift, H_drives, drives, Ls, Ks),
    )

    ode_param_count = p_dim + 2

    # Real-domain PropagatorResults (density stays real via compact 𝒢c generators)
    prop_results = [PropagatorResult{Float64}(n², ode_param_count) for _ = 1:(N_traj-1)]

    # ================================================================
    # Sensitivity ODE (analytical Jacobian via compact 𝒢c generators)
    # ================================================================
    make_f!_sens, n_params = build_density_sensitivity_ode(
        𝒢c_drift,
        𝒢c_drives,
        drives,
        control_to_drives,
        u_dim,
        n,
        order,
    )

    if !isnothing(make_f!_sens)
        sens_probs, sens_state = build_density_sensitivity_problems(
            make_f!_sens,
            n_params,
            n,
            [p₀; Δt₀; t₀],
            N_traj,
        )
    else
        sens_probs = nothing
        sens_state = nothing
    end

    return SplineIntegrator{DensityTrajectory,S,Float64,typeof(alg),typeof(alg_data)}(
        [x],       # x_names: single state
        u,         # u_name
        x_dim,     # x_dim: n² (compact)
        u_dim,     # u_dim
        dim,       # total constraint dim
        tol,
        prop_results,
        n,          # ketdim field stores n (Hilbert space dimension)
        global_names,
        global_dim,
        sens_probs,
        sens_state,
        false,
        nothing,
        nothing,
        nothing,
        0,  # exact_hessian
        false,
        nothing,  # ket-level sensitivity (MultiKetTrajectory only)
        alg,
        alg_data,
    )
end

# ============================================================================ #
# Density-specific sensitivity ODE builder
#
# Uses compact Lindbladian generators 𝒢c ∈ ℝ^{n²×n²} directly.
# Same bilinear structure as the ket/unitary case: dΦ/dτ = Δt·𝒢c(u)·Φ
# Augmented state: [Φc_vec; S₁_vec; ...; Sₙ_vec]
# ============================================================================ #

"""
    build_density_sensitivity_ode(𝒢c_drift, 𝒢c_drives, drives, control_to_drives, u_dim, n, order)

Build the sensitivity-augmented ODE for the compact Lindbladian propagator.

Returns `(make_f!_sens, n_params)` analogous to `build_sensitivity_ode` for
ket/unitary — a CLOSURE FACTORY, not a closure: `build_density_sensitivity_problems`
calls it once per knot problem so each knot's RHS owns its own `u_interp` scratch.
Sharing one closure across the `Threads.@threads` knot loop in `eval_jacobian` was
issue #354. The key difference from the ket/unitary builder is the state dimension:
`(n²)²` for the propagator instead of `(2n)²`.

Supports NonlinearDrive via `drive_coeff`/`drive_coeff_jac` dispatch and `control_to_drives`
mapping (same pattern as the ket/unitary sensitivity ODE in spline_integrator_type.jl).
"""
function build_density_sensitivity_ode(
    𝒢c_drift::AbstractMatrix,
    𝒢c_drives::AbstractVector,
    drives::AbstractVector,
    control_to_drives::Vector{Vector{Int}},
    u_dim::Int,
    n::Int,    # Hilbert space dimension
    order::Int,
)
    n² = n^2
    n_terms = length(𝒢c_drives)
    Φ_dim = n² * n²  # n⁴

    if order == 1
        n_params = 2 * u_dim + 2  # [uₖ, uₖ₊₁, Δt, t]
    elseif order == 3
        n_params = 4 * u_dim + 2  # [uₖ, uₖ₊₁, duₖ, duₖ₊₁, Δt, t]
    else
        error("Unsupported spline order: $order")
    end

    if isempty(𝒢c_drives)
        return nothing, n_params
    end

    # Pre-compute active controls per drive term (for Δt sensitivity chain rule)
    drives_active = [active_controls(d) for d in drives]

    make_f!_sens = if order == 1
        # CLOSURE FACTORY — one call per knot problem, each with its OWN `u_interp`.
        #
        # This MUST be a factory (`() -> begin u_interp = zeros(…); (dx,x,p,τ) -> …`)
        # and not a `let`-bound buffer. `eval_jacobian` solves the `N-1` knot
        # sensitivity problems under `Threads.@threads`, so a single closure shared
        # across the problems is a WRITE RACE on `u_interp`: every thread overwrites
        # the interpolated control of every other, and each `𝒢c(u(τ))` is then built
        # from a mix of knots. That was issue #354 — see
        # `build_density_sensitivity_problems` for the measured signature.
        () -> begin
            u_interp = zeros(u_dim)
            (dx, x, p, τ) -> begin
                uₖ = @view p[1:u_dim]
                uₖ₊₁ = @view p[(u_dim+1):2u_dim]
                Δtₖ = p[end-1]

                onemτ = 1 - τ

                # Interpolate full control vector
                @inbounds for i = 1:u_dim
                    u_interp[i] = onemτ * uₖ[i] + τ * uₖ₊₁[i]
                end

                @inline function apply_𝒢c!(dM, M)
                    mul!(dM, 𝒢c_drift, M, Δtₖ, false)
                    @inbounds for t_idx = 1:n_terms
                        c = drive_coeff(drives[t_idx], u_interp)
                        mul!(dM, 𝒢c_drives[t_idx], M, Δtₖ * c, true)
                    end
                end

                # dΦ/dτ = Δt·𝒢c·Φ
                Φ_vec = @view x[1:Φ_dim]
                Φ_mat = reshape(Φ_vec, n², n²)
                dΦ_vec = @view dx[1:Φ_dim]
                dΦ_mat = reshape(dΦ_vec, n², n²)
                apply_𝒢c!(dΦ_mat, Φ_mat)

                # Sensitivities for uₖ[1..u_dim]: spline basis = (1-τ)
                @inbounds for j = 1:u_dim
                    offset = Φ_dim + (j - 1) * Φ_dim
                    Sⱼ_mat = reshape(@view(x[(offset+1):(offset+Φ_dim)]), n², n²)
                    dSⱼ_mat = reshape(@view(dx[(offset+1):(offset+Φ_dim)]), n², n²)
                    apply_𝒢c!(dSⱼ_mat, Sⱼ_mat)
                    for t_idx in control_to_drives[j]
                        dc = drive_coeff_jac(drives[t_idx], u_interp, j)
                        if dc != 0.0
                            mul!(dSⱼ_mat, 𝒢c_drives[t_idx], Φ_mat, Δtₖ * onemτ * dc, true)
                        end
                    end
                end

                # Sensitivities for uₖ₊₁[1..u_dim]: spline basis = τ
                @inbounds for j = 1:u_dim
                    offset = Φ_dim + (u_dim + j - 1) * Φ_dim
                    Sⱼ_mat = reshape(@view(x[(offset+1):(offset+Φ_dim)]), n², n²)
                    dSⱼ_mat = reshape(@view(dx[(offset+1):(offset+Φ_dim)]), n², n²)
                    apply_𝒢c!(dSⱼ_mat, Sⱼ_mat)
                    for t_idx in control_to_drives[j]
                        dc = drive_coeff_jac(drives[t_idx], u_interp, j)
                        if dc != 0.0
                            mul!(dSⱼ_mat, 𝒢c_drives[t_idx], Φ_mat, Δtₖ * τ * dc, true)
                        end
                    end
                end

                # Sensitivity for Δt
                let j = 2 * u_dim + 1
                    offset = Φ_dim + (j - 1) * Φ_dim
                    Sⱼ_mat = reshape(@view(x[(offset+1):(offset+Φ_dim)]), n², n²)
                    dSⱼ_mat = reshape(@view(dx[(offset+1):(offset+Φ_dim)]), n², n²)
                    apply_𝒢c!(dSⱼ_mat, Sⱼ_mat)
                    # Additional forcing: (1/Δt) * dΦ/dτ = 𝒢c(u) · Φ (without Δt scaling)
                    mul!(dSⱼ_mat, 𝒢c_drift, Φ_mat, 1.0, true)
                    @inbounds for t_idx = 1:n_terms
                        c = drive_coeff(drives[t_idx], u_interp)
                        mul!(dSⱼ_mat, 𝒢c_drives[t_idx], Φ_mat, c, true)
                    end
                end

                # Sensitivity for t (zero forcing for time-independent systems)
                let j = 2 * u_dim + 2
                    offset = Φ_dim + (j - 1) * Φ_dim
                    Sⱼ_mat = reshape(@view(x[(offset+1):(offset+Φ_dim)]), n², n²)
                    dSⱼ_mat = reshape(@view(dx[(offset+1):(offset+Φ_dim)]), n², n²)
                    apply_𝒢c!(dSⱼ_mat, Sⱼ_mat)
                end

                return nothing
            end
        end

    else  # order == 3
        # Same factory contract as the linear branch above (#354).
        () -> begin
            u_interp = zeros(u_dim)
            (dx, x, p, τ) -> begin
                uₖ = @view p[1:u_dim]
                uₖ₊₁ = @view p[(u_dim+1):2u_dim]
                duₖ = @view p[(2u_dim+1):3u_dim]
                duₖ₊₁ = @view p[(3u_dim+1):4u_dim]
                Δtₖ = p[end-1]

                τ2 = τ * τ
                τ3 = τ2 * τ
                h00 = 2τ3 - 3τ2 + 1
                h10 = (τ3 - 2τ2 + τ) * Δtₖ
                h01 = -2τ3 + 3τ2
                h11 = (τ3 - τ2) * Δtₖ
                dh10_dΔt = τ3 - 2τ2 + τ
                dh11_dΔt = τ3 - τ2

                # Interpolate full control vector
                @inbounds for i = 1:u_dim
                    u_interp[i] =
                        h00 * uₖ[i] + h10 * duₖ[i] + h01 * uₖ₊₁[i] + h11 * duₖ₊₁[i]
                end

                # dΦ/dτ
                Φ_vec = @view x[1:Φ_dim]
                Φ_mat = reshape(Φ_vec, n², n²)
                dΦ_vec = @view dx[1:Φ_dim]
                dΦ_mat = reshape(dΦ_vec, n², n²)

                mul!(dΦ_mat, 𝒢c_drift, Φ_mat, Δtₖ, false)
                @inbounds for t_idx = 1:n_terms
                    c = drive_coeff(drives[t_idx], u_interp)
                    mul!(dΦ_mat, 𝒢c_drives[t_idx], Φ_mat, Δtₖ * c, true)
                end

                function compute_𝒢c_times_S!(dS_mat, S_mat)
                    mul!(dS_mat, 𝒢c_drift, S_mat, Δtₖ, false)
                    @inbounds for t_idx = 1:n_terms
                        c = drive_coeff(drives[t_idx], u_interp)
                        mul!(dS_mat, 𝒢c_drives[t_idx], S_mat, Δtₖ * c, true)
                    end
                end

                # uₖ sensitivities: spline basis = h00
                @inbounds for j = 1:u_dim
                    offset = Φ_dim + (j - 1) * Φ_dim
                    Sⱼ_mat = reshape(@view(x[(offset+1):(offset+Φ_dim)]), n², n²)
                    dSⱼ_mat = reshape(@view(dx[(offset+1):(offset+Φ_dim)]), n², n²)
                    compute_𝒢c_times_S!(dSⱼ_mat, Sⱼ_mat)
                    for t_idx in control_to_drives[j]
                        dc = drive_coeff_jac(drives[t_idx], u_interp, j)
                        if dc != 0.0
                            mul!(dSⱼ_mat, 𝒢c_drives[t_idx], Φ_mat, Δtₖ * h00 * dc, true)
                        end
                    end
                end

                # uₖ₊₁ sensitivities: spline basis = h01
                @inbounds for j = 1:u_dim
                    offset = Φ_dim + (u_dim + j - 1) * Φ_dim
                    Sⱼ_mat = reshape(@view(x[(offset+1):(offset+Φ_dim)]), n², n²)
                    dSⱼ_mat = reshape(@view(dx[(offset+1):(offset+Φ_dim)]), n², n²)
                    compute_𝒢c_times_S!(dSⱼ_mat, Sⱼ_mat)
                    for t_idx in control_to_drives[j]
                        dc = drive_coeff_jac(drives[t_idx], u_interp, j)
                        if dc != 0.0
                            mul!(dSⱼ_mat, 𝒢c_drives[t_idx], Φ_mat, Δtₖ * h01 * dc, true)
                        end
                    end
                end

                # duₖ sensitivities: spline basis = h10
                @inbounds for j = 1:u_dim
                    offset = Φ_dim + (2u_dim + j - 1) * Φ_dim
                    Sⱼ_mat = reshape(@view(x[(offset+1):(offset+Φ_dim)]), n², n²)
                    dSⱼ_mat = reshape(@view(dx[(offset+1):(offset+Φ_dim)]), n², n²)
                    compute_𝒢c_times_S!(dSⱼ_mat, Sⱼ_mat)
                    for t_idx in control_to_drives[j]
                        dc = drive_coeff_jac(drives[t_idx], u_interp, j)
                        if dc != 0.0
                            mul!(dSⱼ_mat, 𝒢c_drives[t_idx], Φ_mat, Δtₖ * h10 * dc, true)
                        end
                    end
                end

                # duₖ₊₁ sensitivities: spline basis = h11
                @inbounds for j = 1:u_dim
                    offset = Φ_dim + (3u_dim + j - 1) * Φ_dim
                    Sⱼ_mat = reshape(@view(x[(offset+1):(offset+Φ_dim)]), n², n²)
                    dSⱼ_mat = reshape(@view(dx[(offset+1):(offset+Φ_dim)]), n², n²)
                    compute_𝒢c_times_S!(dSⱼ_mat, Sⱼ_mat)
                    for t_idx in control_to_drives[j]
                        dc = drive_coeff_jac(drives[t_idx], u_interp, j)
                        if dc != 0.0
                            mul!(dSⱼ_mat, 𝒢c_drives[t_idx], Φ_mat, Δtₖ * h11 * dc, true)
                        end
                    end
                end

                # Δt sensitivity
                let j = 4 * u_dim + 1
                    offset = Φ_dim + (j - 1) * Φ_dim
                    Sⱼ_mat = reshape(@view(x[(offset+1):(offset+Φ_dim)]), n², n²)
                    dSⱼ_mat = reshape(@view(dx[(offset+1):(offset+Φ_dim)]), n², n²)
                    compute_𝒢c_times_S!(dSⱼ_mat, Sⱼ_mat)
                    # 𝒢c(u) · Φ (without Δt scaling)
                    mul!(dSⱼ_mat, 𝒢c_drift, Φ_mat, 1.0, true)
                    @inbounds for t_idx = 1:n_terms
                        c = drive_coeff(drives[t_idx], u_interp)
                        mul!(dSⱼ_mat, 𝒢c_drives[t_idx], Φ_mat, c, true)
                    end
                    # Chain rule: Δt · Σ_drives (∂coeff/∂u_j) · (∂u_j/∂Δt) · 𝒢c_drive · Φ
                    @inbounds for t_idx = 1:n_terms
                        ac = drives_active[t_idx]
                        iter = isempty(ac) ? (1:u_dim) : ac
                        total_du_dΔt = 0.0
                        for i in iter
                            dc = drive_coeff_jac(drives[t_idx], u_interp, i)
                            if dc != 0.0
                                du_i_dΔt = dh10_dΔt * duₖ[i] + dh11_dΔt * duₖ₊₁[i]
                                total_du_dΔt += dc * du_i_dΔt
                            end
                        end
                        if total_du_dΔt != 0.0
                            mul!(dSⱼ_mat, 𝒢c_drives[t_idx], Φ_mat, Δtₖ * total_du_dΔt, true)
                        end
                    end
                end

                # t sensitivity (zero forcing for time-independent)
                let j = 4 * u_dim + 2
                    offset = Φ_dim + (j - 1) * Φ_dim
                    Sⱼ_mat = reshape(@view(x[(offset+1):(offset+Φ_dim)]), n², n²)
                    dSⱼ_mat = reshape(@view(dx[(offset+1):(offset+Φ_dim)]), n², n²)
                    compute_𝒢c_times_S!(dSⱼ_mat, Sⱼ_mat)
                end

                return nothing
            end
        end
    end

    return make_f!_sens, n_params
end

"""
    build_density_sensitivity_problems(make_f!_sens, n_params, n, ode_p₀, N)

Create ODE problems and preallocated state buffer for the density sensitivity ODE.

`make_f!_sens` is a CLOSURE FACTORY and is invoked ONCE PER KNOT PROBLEM, exactly
as [`build_sensitivity_problems`](@ref) does for the ket/unitary cells. Each knot's
RHS therefore owns its own `u_interp` scratch.

**This is issue #354.** The density builder previously bound `u_interp` in a single
`let` and handed the SAME closure to all `N-1` problems, while `eval_jacobian`
solves them under `Threads.@threads` — a write race on the interpolated control.
The corruption was not uniform across the Jacobian, which is why it survived the
`test_integrator(…; atol = 1e-4)` gates: the propagator and the control forcings
carry a `Δt` prefactor (and the control forcings a spline-basis weight on top),
whereas the **`Δt` sensitivity's own forcing `𝒢c(u)·Φ` carries no prefactor at
all**, so the same corrupted `u_interp` lands hardest exactly there. Measured
against a cold central difference of the integrator's own residual, threaded, at
`ketdim = 2`, `N = 5`, linear spline:

    column      max |J_asm − J_fd|   (before)   (after)
    Δt                    2.74e-1              4.7e-7
    ρ⃗̃ (state)             1.04e-2              2.3e-7
    u  (control)          6.07e-3              3.0e-8

i.e. the `Δt` column was wrong by the FULL MAGNITUDE of its own entries, which is
what #346 measured as `2.3e-1` and filed rather than absorbed.

# Arguments
- `n::Int`: Hilbert space dimension (propagator is n²×n²)
"""
function build_density_sensitivity_problems(
    make_f!_sens,
    n_params::Int,
    n::Int,    # Hilbert space dim
    ode_p₀::AbstractVector,
    N::Int,
)
    n² = n^2
    Φ_dim = n² * n²  # n⁴
    sens_state_dim = Φ_dim * (1 + n_params)

    # Initial condition: Φ₀ = I_{n²}, all sensitivities = 0
    sens_x₀ = zeros(sens_state_dim)
    sens_x₀[1:Φ_dim] = vec(Matrix{Float64}(I, n², n²))

    # `make_f!_sens()` PER PROBLEM — a fresh `u_interp` per knot. Sharing one
    # closure here is the #354 race; see this function's docstring.
    sens_probs = [ODEProblem(make_f!_sens(), sens_x₀, (0.0, 1.0), ode_p₀) for _ = 1:(N-1)]
    sens_state = zeros(sens_state_dim)

    return sens_probs, sens_state
end


# ============================================================================ #
# Call operator for DensityTrajectory
#
# Uses state-forward integration (like KetTrajectory) for efficiency:
# integrates the n²-dim state directly instead of the full n⁴-dim propagator.
# The f! adapts via n_cols = length(x) ÷ n².
# ============================================================================ #

"""
    (𝒮::SplineIntegrator{DensityTrajectory})(δₖ, zₖ, zₖ₊₁, k, globals)

Evaluate the constraint: δₖ = xₖ₊₁ - f(xₖ, p)

where f integrates the Lindbladian forward from xₖ using the ODE.
"""
@views function (𝒮::SplineIntegrator{DensityTrajectory})(
    δₖ::AbstractVector,
    zₖ::KnotPoint,
    zₖ₊₁::KnotPoint,
    k::Int,
    globals::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    xₖ = zₖ[x_name(𝒮)]
    xₖ₊₁ = zₖ₊₁[x_name(𝒮)]

    # Build parameter vector
    pₖ = build_ode_params(𝒮, zₖ, zₖ₊₁, globals)

    # State-forward integration: remake with u0=xₖ (length n², not n⁴)
    data = 𝒮.alg_data::DensityLindbladData
    Φₖ_prob = remake(data.Φ_probs[k], u0 = xₖ, p = pₖ)
    fₖ = solve(Φₖ_prob, Tsit5(); abstol = 𝒮.tol, reltol = 𝒮.tol, saveat = 1.0).u[end]

    # Constraint
    δₖ[:] = xₖ₊₁ - fₖ

    return nothing
end

# ============================================================================ #
# Dispatch helpers for DensityTrajectory
# ============================================================================ #

"""
    compute_ode_jacobian!(𝒮::SplineIntegrator{DensityTrajectory}, ...)

Compute ODE Jacobian for DensityTrajectory using analytical sensitivity equations.
"""
@views function compute_ode_jacobian!(
    𝒮::SplineIntegrator{DensityTrajectory},
    zₖ::KnotPoint,
    zₖ₊₁::KnotPoint,
    k::Int,
    globals::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    isnothing(𝒮.sens_probs) && error(
        "SplineIntegrator{DensityTrajectory} requires explicit H_drift/H_drives matrices " *
        "for Jacobian computation.",
    )
    pₖ = build_ode_params(𝒮, zₖ, zₖ₊₁, globals)
    return _compute_ode_jacobian_analytical!(𝒮, pₖ, k)
end

"""
    get_state_vector(𝒮::SplineIntegrator{DensityTrajectory}, zₖ)

Extract the compact density state vector from a knot point.
"""
function get_state_vector(𝒮::SplineIntegrator{DensityTrajectory}, zₖ::KnotPoint)
    return zₖ[x_name(𝒮)]
end
