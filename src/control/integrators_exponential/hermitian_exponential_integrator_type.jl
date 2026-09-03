export HermitianExponentialIntegrator, x_name, single_state_dim, extract_globals

"""
    HermitianExponentialIntegrator{T<:AbstractQuantumTrajectory} <: AbstractExponentialIntegrator

Unified Hermitian exponential integrator that dispatches on quantum trajectory type.

The type parameter `T` determines the specific integration method:
- `UnitaryTrajectory`: Block-structured integration for unitary evolution
- `KetTrajectory`: Direct integration for ket state evolution
- `MultiKetTrajectory`: Shared propagator for multiple kets

Uses matrix exponential for Hermitian generators: Ũ(Δt) = exp(-iΔt*H).

# Fields
- `G::Function`: Generator function G(u) returning isomorphic anti-Hermitian operator iso(-iH)
- `H::Function`: Hamiltonian function H(u) returning complex Hermitian matrix (used by exp_eigen)
- `x_names::Vector{Symbol}`: Names of state variables (single element for non-ensemble)
- `u_name::Symbol`: Name of control variable
- `x_dim::Int`: Total constraint dimension per knot point
- `u_dim::Int`: Dimension of control (without globals)
- `var_dim::Int`: Total variable dimension (states + controls + time)
- `dim::Int`: Total constraint dimension
- `ketdim::Int`: Dimension of ket space
- `∂ℰs::Vector{SparseMatrixCSC{Float64, Int}}`: Preallocated Jacobian structures per knot point
- `μ∂²ℰs::Vector{SparseMatrixCSC{Float64, Int}}`: Preallocated Hessian structures per knot point
- `z_dim::Int`: Dimension of single knot point
- `var_comps::Vector{Int}`: Component indices for variables
- `Id::Union{Nothing, Matrix{Float64}}`: Preallocated identity matrix (UnitaryTrajectory only)
- `∂u∂Δt_bufs::Union{Nothing, Vector{Vector{Float64}}}`: Per-thread buffers for u-Δt cross term (MultiKetTrajectory only)
- `∂²u_bufs::Union{Nothing, Vector{Matrix{Float64}}}`: Per-thread buffers for u-u Hessian block (MultiKetTrajectory only)
- `global_names::Vector{Symbol}`: Names of global (time-invariant) variables
- `global_dim::Int`: Dimension of global variables

# Per-thread exp_eigen! buffer set
The seven `exp_eigen!` work buffers are stored as `Vector{...}` with one entry per
OS thread (sized to `Threads.maxthreadid()` at construction). Jacobian/Hessian loops
run under `Threads.@threads`; concurrent knot points on different threads get disjoint
buffers by indexing `ℰ.H_bufs[Threads.threadid()]` etc. Sharing a single buffer
across threads would produce `DimensionMismatch` errors in `λ_buf .= F.values` or
silent numerical corruption via aliased eigendecomposition scratch.
"""
struct HermitianExponentialIntegrator{T<:AbstractQuantumTrajectory} <:
       AbstractExponentialIntegrator
    G::Function
    H::Function
    x_names::Vector{Symbol}  # Vector for ensemble support; single element for non-ensemble
    u_name::Symbol
    x_dim::Int
    u_dim::Int
    var_dim::Int
    dim::Int
    ketdim::Int
    ∂ℰs::Vector{SparseMatrixCSC{Float64,Int}}
    μ∂²ℰs::Vector{SparseMatrixCSC{Float64,Int}}
    z_dim::Int
    var_comps::Vector{Int}
    # Type-specific caches
    Id::Union{Nothing,Matrix{Float64}}
    # Per-thread scratch for MultiKetTrajectory Hessian aggregation (one entry per thread).
    # Shared buffers were racing inside `Threads.@threads for k` loops in
    # eval_hessian_of_lagrangian; see PR #52 follow-up on thread safety.
    ∂u∂Δt_bufs::Union{Nothing,Vector{Vector{Float64}}}
    ∂²u_bufs::Union{Nothing,Vector{Matrix{Float64}}}
    # Global variable support
    global_names::Vector{Symbol}
    global_dim::Int
    # Per-thread preallocated work buffers for exp_eigen! (7 per thread, see Task 2 + thread-race fix).
    # Each outer `Vector` has length `Threads.maxthreadid()` and is indexed by `threadid()`
    # inside per-knot methods running under `Threads.@threads`.
    H_bufs::Vector{Matrix{ComplexF64}}         # size (ketdim, ketdim); preserved Hamiltonian
    λ_bufs::Vector{Vector{Float64}}            # size ketdim
    V_bufs::Vector{Matrix{ComplexF64}}         # size (ketdim, ketdim); eigenvectors after eigen!
    cis_diag_bufs::Vector{Vector{ComplexF64}}  # size ketdim
    tmp_bufs::Vector{Matrix{ComplexF64}}       # size (ketdim, ketdim); V · Diagonal(cis)
    work_bufs::Vector{Matrix{ComplexF64}}      # size (ketdim, ketdim); final complex product
    expG_bufs::Vector{Matrix{Float64}}         # size (2*ketdim, 2*ketdim); isomorphic output
    # Hessian approximation
    gauss_newton::Bool
    # ── Matrix-free Jacobian opt-in (Piccolissimo.jl#205, MultiKet CPU path) ──
    # When `true` (and the drive is affine, i.e. `H_dirs !== nothing`), the MultiKet
    # `eval_jacobian` returns a matrix-free, CPU-threaded `CPUJacobianOp` instead of
    # assembling the ~22M-nnz sparse dynamics Jacobian — wiring the dormant threaded
    # JVP/VJP into the Altissimo NewtonCG hot path. Default `false`: ket/unitary
    # variants ignore it (keep sparse), and the Ipopt/MadNLP KKT-factorization path
    # (which needs the assembled sparse Jacobian) stays intact. Modeled exactly on
    # `gauss_newton`: a kwarg threaded through all three per-trajectory constructors.
    matrix_free::Bool
    # ── Analytic Daleckii–Krein derivative path (Piccolissimo.jl#202, slice ③) ──
    # `use_analytical` toggles the affine-drive DENSE derivative assembly between the
    # analytic first-order DK kernel and the retained ForwardDiff witness. It mirrors
    # the analytical-derivative toggle called for in #199; the default is `true`
    # (analytic DK) as of slice ⑧ (#207, Piccolissimo ADR 0003) now every variant is
    # wired — ForwardDiff stays reachable via `use_analytical=false` as the witness and
    # is selected automatically for nonlinear drives (see `H_dirs` below). `H_dirs` holds
    # the constant drive directions ∂H/∂uₚ for the EXTENDED control [controls; globals]
    # when the Hamiltonian is affine-in-parameter; it is `nothing` when H is nonlinear
    # (e.g. a `NonlinearDrive` system), in which case the DK path is unavailable and
    # ForwardDiff is used even with the flag on. The DK scratch below is per-thread
    # (one entry per OS thread, indexed by `threadid()` inside the `Threads.@threads`
    # Jacobian/Hessian loops), matching the `exp_eigen!` buffer discipline.
    use_analytical::Bool
    H_dirs::Union{Nothing,Vector{Matrix{ComplexF64}}}
    dk_ws_bufs::Vector{DaleckiiKreinWorkspace{ComplexF64,Matrix{ComplexF64}}}
    dk_dΦ_bufs::Vector{Matrix{ComplexF64}}     # size (ketdim, ketdim); ∂Φ/∂uₚ output
    dk_vin_bufs::Vector{Vector{ComplexF64}}    # size ketdim; complex ψ / μ input
    dk_vout_bufs::Vector{Vector{ComplexF64}}   # size ketdim; complex applied output
    # ── Second-order Daleckii–Krein scratch (Piccolissimo.jl#203, slice ④) ──
    # Feeds the analytic exact-Hessian (`!gauss_newton`) parameter-parameter blocks
    # (u-u / u-g / g-g via the three-index second divided difference, plus the
    # closed-form u-Δt / Δt-g cross terms). Per-thread (one entry per OS thread,
    # indexed by `threadid()` inside the `Threads.@threads` Hessian loop), matching
    # the first-order DK buffer discipline above. Unused when `gauss_newton=true`
    # or when the DK path is off (allocated for a uniform struct across variants).
    dk_so_ws_bufs::Vector{DaleckiiKreinSecondOrderWorkspace{ComplexF64,Matrix{ComplexF64}}}
    dk_μmat_bufs::Vector{Matrix{ComplexF64}}   # size (ketdim, ketdim); μ_mat = χ ψ† outer product
    dk_psi_bufs::Vector{Vector{ComplexF64}}    # size ketdim; complex state ψ
    dk_y_bufs::Vector{Vector{ComplexF64}}      # size ketdim; y = Φ ψ (parameter-Δt terms)
    dk_Hc_bufs::Vector{Vector{ComplexF64}}     # size ketdim; H χ (parameter-Δt terms)
end

"""
    _alloc_exp_eigen_bufs(ketdim::Int, nthreads::Int)

Allocate the 7 `exp_eigen!` work buffer vectors, each of length `nthreads`. Returned
as a NamedTuple so callers can `splat` into the `HermitianExponentialIntegrator`
constructor positional-argument block.
"""
function _alloc_exp_eigen_bufs(ketdim::Int, nthreads::Int)
    return (
        H_bufs = [zeros(ComplexF64, ketdim, ketdim) for _ = 1:nthreads],
        λ_bufs = [zeros(Float64, ketdim) for _ = 1:nthreads],
        V_bufs = [zeros(ComplexF64, ketdim, ketdim) for _ = 1:nthreads],
        cis_diag_bufs = [zeros(ComplexF64, ketdim) for _ = 1:nthreads],
        tmp_bufs = [zeros(ComplexF64, ketdim, ketdim) for _ = 1:nthreads],
        work_bufs = [zeros(ComplexF64, ketdim, ketdim) for _ = 1:nthreads],
        expG_bufs = [zeros(Float64, 2 * ketdim, 2 * ketdim) for _ = 1:nthreads],
    )
end

# ============================================================================ #
# Analytic Daleckii–Krein setup (Piccolissimo.jl#202, slice ③)
# ============================================================================ #

"""
    _build_affine_directions(Hfun, nparams, ketdim; rtol=1e-8) -> Vector{Matrix{ComplexF64}} | nothing

Build the constant drive directions `Hₚ = ∂H/∂uₚ` (`p = 1..nparams` over the EXTENDED
control `[controls; globals]`) for an affine-in-parameter Hamiltonian, by evaluating the
exact linear coefficients `Hₚ = H(eₚ) - H(0)` once at construction. For an affine map this
is exact — no finite-difference truncation, no derivative approximation — so it does not
violate the "no FD/autodiff in the analytic path" invariant (it is one-time setup, not a
runtime derivative). Affinity is verified at two generic probe points; if the reconstruction
`H(0) + Σ uₚ Hₚ` disagrees with `H(u)` (nonlinear drive) or any evaluation fails, returns
`nothing`, signalling the caller to keep the ForwardDiff path.
"""
function _build_affine_directions(Hfun, nparams::Int, ketdim::Int; rtol::Real = 1e-8)
    nparams == 0 && return nothing
    try
        u0 = zeros(nparams)
        H0 = Matrix{ComplexF64}(Hfun(u0))
        (size(H0, 1) == ketdim && size(H0, 2) == ketdim) || return nothing
        dirs = Vector{Matrix{ComplexF64}}(undef, nparams)
        e = zeros(nparams)
        for j = 1:nparams
            fill!(e, 0.0)
            e[j] = 1.0
            dirs[j] = Matrix{ComplexF64}(Hfun(e)) .- H0
        end
        # Affinity check at two generic probe points (a single probe could miss a
        # nonlinear term that happens to vanish there).
        for probe in (
            Float64[0.37 + 0.11 * j for j = 1:nparams],
            Float64[-0.53 + 0.19 * j for j = 1:nparams],
        )
            Hpred = copy(H0)
            for j = 1:nparams
                @. Hpred += probe[j] * dirs[j]
            end
            Hact = Matrix{ComplexF64}(Hfun(probe))
            norm(Hact - Hpred) > rtol * max(norm(Hact), one(rtol)) && return nothing
        end
        return dirs
    catch
        return nothing
    end
end

"""
    _alloc_dk_bufs(ketdim, nthreads)

Allocate the per-thread first-order Daleckii–Krein scratch consumed by the analytic
Jacobian / Gauss-Newton Hessian cross-term assembly: a [`DaleckiiKreinWorkspace`](@ref)
(divided-difference matrix + two conjugation temporaries), a `∂Φ` output matrix, and
complex in/out vectors — each entry sized to `ketdim` and indexed by `threadid()` inside
the `Threads.@threads` derivative loops. Returned as a NamedTuple to splat into the
`HermitianExponentialIntegrator` constructor's positional block.
"""
function _alloc_dk_bufs(ketdim::Int, nthreads::Int)
    return (
        dk_ws_bufs = [DaleckiiKreinWorkspace(ketdim) for _ = 1:nthreads],
        dk_dΦ_bufs = [zeros(ComplexF64, ketdim, ketdim) for _ = 1:nthreads],
        dk_vin_bufs = [zeros(ComplexF64, ketdim) for _ = 1:nthreads],
        dk_vout_bufs = [zeros(ComplexF64, ketdim) for _ = 1:nthreads],
    )
end

"""
    _alloc_dk_so_bufs(ketdim, nthreads)

Allocate the per-thread SECOND-order Daleckii–Krein scratch consumed by the analytic
exact-Hessian (`!gauss_newton`) parameter-parameter assembly (Piccolissimo.jl#203, slice ④):
a [`DaleckiiKreinSecondOrderWorkspace`](@ref) (the three-index contraction scratch), the
`μ_mat = χ ψ†` outer-product buffer, and the complex vector temporaries `ψ`, `y = Φ ψ`,
and `H χ` used by the closed-form parameter-Δt cross terms — each entry sized to `ketdim`
and indexed by `threadid()` inside the `Threads.@threads` Hessian loop. Returned as a
NamedTuple to splat into the `HermitianExponentialIntegrator` constructor's positional block.
"""
function _alloc_dk_so_bufs(ketdim::Int, nthreads::Int)
    return (
        dk_so_ws_bufs = [DaleckiiKreinSecondOrderWorkspace(ketdim) for _ = 1:nthreads],
        dk_μmat_bufs = [zeros(ComplexF64, ketdim, ketdim) for _ = 1:nthreads],
        dk_psi_bufs = [zeros(ComplexF64, ketdim) for _ = 1:nthreads],
        dk_y_bufs = [zeros(ComplexF64, ketdim) for _ = 1:nthreads],
        dk_Hc_bufs = [zeros(ComplexF64, ketdim) for _ = 1:nthreads],
    )
end

"""
    _iso_vec_to_complex!(vc, ṽ, n)

Write the complex vector encoded by the length-`2n` isomorphic real vector `ṽ`
(`ṽ = [Re; Im]`) into the length-`n` complex buffer `vc`.
"""
@inline function _iso_vec_to_complex!(
    vc::AbstractVector{ComplexF64},
    ṽ::AbstractVector{<:Real},
    n::Int,
)
    @inbounds for r = 1:n
        vc[r] = complex(ṽ[r], ṽ[r+n])
    end
    return vc
end

"""
    _dk_fill_iso_jac_block!(block, V, M, dirs, vin, dΦ, tmp1, tmp2, vout, n; adjoint_op=false)

Fill an isomorphic derivative block (`2n × length(dirs)`) whose column `p` is
`-[Re(y); Im(y)]`, where `y = (∂ₚΦ) vin` (`adjoint_op=false`) or `y = (∂ₚΦ)† vin`
(`adjoint_op=true`), and `∂ₚΦ = V (M ∘ (V† Hₚ V)) V†` is the first-order Daleckii–Krein
derivative of `Φ = exp(-iΔt H)` in direction `Hₚ = dirs[p]`. `M` is prebuilt once per knot
via [`dk_divided_difference!`](@ref) and reused across directions. This is the analytic
replacement for the multiket ForwardDiff-∘-`expv` Jacobian (`adjoint_op=false`, `vin=ψ`)
and Gauss-Newton Hessian cross-term (`adjoint_op=true`, `vin=μ`) call sites (#202).
"""
@views function _dk_fill_iso_jac_block!(
    block::AbstractMatrix{<:Real},
    V::AbstractMatrix{ComplexF64},
    M::AbstractMatrix{ComplexF64},
    dirs,
    vin::AbstractVector{ComplexF64},
    dΦ::AbstractMatrix{ComplexF64},
    tmp1::AbstractMatrix{ComplexF64},
    tmp2::AbstractMatrix{ComplexF64},
    vout::AbstractVector{ComplexF64},
    n::Int;
    adjoint_op::Bool = false,
)
    for (p, Hₚ) in enumerate(dirs)
        dk_apply!(dΦ, V, M, Hₚ, tmp1, tmp2)
        if adjoint_op
            mul!(vout, dΦ', vin)
        else
            mul!(vout, dΦ, vin)
        end
        @inbounds for r = 1:n
            block[r, p] = -real(vout[r])
            block[r+n, p] = -imag(vout[r])
        end
    end
    return block
end

"""
    _dk_fill_iso_jac_block_unitary!(block, V, M, dirs, src, dΦ, tmp1, tmp2, vin, vout, ketdim; adjoint_op=false)

Block-per-column variant of `_dk_fill_iso_jac_block!` for the UNITARY variant,
whose state `Ũ⃗` is a column-stacked isomorphic MATRIX (`ketdim` isomorphic-ket columns,
each of length `2·ketdim`) rather than a single ket. Fills `block`
(`(2·ketdim·ketdim) × length(dirs)`) whose column `p` is `vec(-(∂ₚΦ) U)` in isomorphic
form (`adjoint_op=false`) or `vec(-(∂ₚΦ)† U)` (`adjoint_op=true`), where `U` is the complex
matrix encoded by the length-`2·ketdim·ketdim` isomorphic source vector `src` and
`∂ₚΦ = V (M ∘ (V† Hₚ V)) V†` is the first-order Daleckii–Krein derivative in direction
`Hₚ = dirs[p]`. `∂ₚΦ` is built ONCE per direction and applied to all `ketdim` columns
(reusing the eigenbasis the forward step produced), so the cost is one `dk_apply!` per
direction plus `ketdim` matrix–vector applies — not one `dk_apply!` per (direction, column).
This is the analytic replacement for the unitary ForwardDiff-∘-`expv` Jacobian
(`adjoint_op=false`, `src=Ũ⃗ₖ`) and Gauss-Newton Hessian cross-term (`adjoint_op=true`,
`src=μₖ`) call sites (#204).
"""
@views function _dk_fill_iso_jac_block_unitary!(
    block::AbstractMatrix{<:Real},
    V::AbstractMatrix{ComplexF64},
    M::AbstractMatrix{ComplexF64},
    dirs,
    src::AbstractVector{<:Real},
    dΦ::AbstractMatrix{ComplexF64},
    tmp1::AbstractMatrix{ComplexF64},
    tmp2::AbstractMatrix{ComplexF64},
    vin::AbstractVector{ComplexF64},
    vout::AbstractVector{ComplexF64},
    ketdim::Int;
    adjoint_op::Bool = false,
)
    col_dim = 2 * ketdim               # length of one isomorphic-ket column
    ncols = length(src) ÷ col_dim      # number of unitary columns (= ketdim)
    for (p, Hₚ) in enumerate(dirs)
        dk_apply!(dΦ, V, M, Hₚ, tmp1, tmp2)
        for c = 1:ncols
            off = (c - 1) * col_dim
            _iso_vec_to_complex!(vin, view(src, (off+1):(off+col_dim)), ketdim)
            if adjoint_op
                mul!(vout, dΦ', vin)
            else
                mul!(vout, dΦ, vin)
            end
            @inbounds for r = 1:ketdim
                block[off+r, p] = -real(vout[r])
                block[off+ketdim+r, p] = -imag(vout[r])
            end
        end
    end
    return block
end

# Convenience accessor for single state name (non-ensemble case)
x_name(ℰ::HermitianExponentialIntegrator) = ℰ.x_names[1]
single_state_dim(ℰ::HermitianExponentialIntegrator) = ℰ.x_dim ÷ length(ℰ.x_names)

# Helper to compute canonical hessian knot dimension for exponential integrator
function canonical_hessian_knot_dim(ℰ::AbstractExponentialIntegrator)
    return ℰ.x_dim + ℰ.u_dim + 1  # x, u, Δt
end

"""
    extract_globals(ℰ::AbstractExponentialIntegrator, traj::NamedTrajectory)

Extract global variable values from trajectory for use in integrator evaluation.
Returns `nothing` if no globals are present.

This is called once per API function (evaluate!, eval_jacobian, etc.) and the
result is passed to inner functions to avoid repeated allocation/lookup.
"""
function extract_globals(ℰ::AbstractExponentialIntegrator, traj::NamedTrajectory)
    ℰ.global_dim == 0 && return nothing
    return vcat(
        [traj.global_data[traj.global_components[name]] for name in ℰ.global_names]...,
    )
end

"""
    build_extended_control(aₖ, globals)

Build the extended control vector [controls..., globals...] for dynamics evaluation.
If globals is nothing, returns the original control vector.
"""
@inline function build_extended_control(aₖ::AbstractVector, globals::Nothing)
    return aₖ
end

@inline function build_extended_control(aₖ::AbstractVector, globals::AbstractVector)
    return vcat(aₖ, globals)
end

"""
    _global_full_cols(ℰ, traj) -> full_cols::Vector{Int}

Positions, in the FULL decision vector `Z`, of the integrator's per-knot
global-derivative block — returned in the SAME order that block is laid out
(`ℰ.global_names` / [`extract_globals`](@ref) order), but pointing at each
global's column under the TRAJECTORY's `global_components` order (the order the
solver/Ipopt indexes `Z`'s globals, base `traj.dim*traj.N`).

Each call site keeps its OWN local (source-block) column range — which differs by
context: `2*ℰ.z_dim` for the preallocated `∂ℰ`/`μ∂²ℰ` matrices, `2*traj.dim` for a
fresh structure template, `2*knot_dim` for the canonical Hessian layout — and
pairs it element-wise with `full_cols`, so `∂F[:, full_cols] = block[:, local]`
lands each global's derivative in the right `Z` column.

Fixes a permuted-`∂c/∂θ` bug: the per-knot block is filled in `global_names`
(insertion) order, but globals frequently reach the trajectory through a `Dict`
(`global_bounds`/`global_params`), whose iteration is hash order — so
`global_names` order ≠ `global_components` order and a naive `i→i` mapping swaps
the global columns, handing Ipopt a wrong Jacobian.
"""
function _global_full_cols(ℰ::AbstractExponentialIntegrator, traj::NamedTrajectory)
    base_full = traj.dim * traj.N
    full_cols = Int[]
    for name in ℰ.global_names
        for c in traj.global_components[name]
            push!(full_cols, base_full + c)
        end
    end
    return full_cols
end

# ============================================================================ #
# Hessian structure - Shared for both types
# Uses CANONICAL ordering: [x_k, u_k, Δt_k, x_{k+1}]
# ============================================================================ #

function hessian_structure(
    x_dim::Int,
    u_dim::Int,
    global_dim::Int = 0;
    gauss_newton::Bool = false,
)
    # Canonical layout per knot:
    # - x: 1:x_dim
    # - u: x_dim+1:x_dim+u_dim
    # - Δt: x_dim+u_dim+1
    # - (knot k+1 block)
    # - g: 2*knot_dim+1:2*knot_dim+global_dim
    knot_dim = x_dim + u_dim + 1

    total_dim = 2 * knot_dim + global_dim
    μ∂²ℰ = spzeros(total_dim, total_dim)

    # Canonical indices for knot k
    x_comps = 1:x_dim
    u_comps = (x_dim+1):(x_dim+u_dim)
    Δt_comp = x_dim + u_dim + 1

    # Cross-terms (always kept)
    # μ∂ₓₖ∂u block
    μ∂²ℰ[x_comps, u_comps] .= 1.0

    # μ∂ₓₖ∂Δtₖ block
    μ∂²ℰ[x_comps, Δt_comp] .= 1.0

    # Parameter-parameter blocks (dropped in GN)
    if !gauss_newton
        # μ∂u∂Δtₖ block
        μ∂²ℰ[u_comps, Δt_comp] .= 1.0

        # μ∂²u block
        μ∂²ℰ[u_comps, u_comps] .= 1.0

        # μ∂²Δtₖ block
        μ∂²ℰ[Δt_comp, Δt_comp] = 1.0
    end

    # Global variable blocks
    if global_dim > 0
        g_comps = (2*knot_dim+1):(2*knot_dim+global_dim)
        μ∂²ℰ[x_comps, g_comps] .= 1.0    # x-g (cross-term, always kept)
        if !gauss_newton
            μ∂²ℰ[u_comps, g_comps] .= 1.0    # u-g (p,p block)
            μ∂²ℰ[Δt_comp:Δt_comp, g_comps] .= 1.0  # Δt-g (p,p block)
            μ∂²ℰ[g_comps, g_comps] .= 1.0    # g-g (p,p block)
        end
    end

    # Full symmetric structure needed because canonical→trajectory index mapping
    # can swap upper/lower triangle positions. triu() taken in eval_hessian_of_lagrangian.
    return sparse(Symmetric(μ∂²ℰ))
end

# ============================================================================ #
# API Methods - Shared implementations
# ============================================================================ #

"""
Build canonical-to-trajectory index mapping for exponential integrator hessian.
"""
function build_hessian_index_mapping(
    ℰ::AbstractExponentialIntegrator,
    traj::NamedTrajectory,
)
    x_dim = ℰ.x_dim
    u_dim = ℰ.u_dim
    z_dim = traj.dim
    knot_dim = canonical_hessian_knot_dim(ℰ)

    # Canonical indices
    x_can_k = 1:x_dim
    u_can_k = (x_dim+1):(x_dim+u_dim)
    Δt_can_k = x_dim + u_dim + 1
    x_can_k1 = knot_dim .+ (1:x_dim)

    # Trajectory indices - handle ensemble (multiple x_names) vs single
    x_traj_k = vcat([collect(traj.components[name]) for name in ℰ.x_names]...)
    u_traj_k = collect(traj.components[ℰ.u_name])
    Δt_traj_k = traj.components[traj.timestep][1]
    x_traj_k1 = z_dim .+ x_traj_k

    canonical_comps = [collect(x_can_k); collect(u_can_k); Δt_can_k; collect(x_can_k1)]
    traj_comps = [x_traj_k; u_traj_k; Δt_traj_k; x_traj_k1]

    return canonical_comps, traj_comps
end

function get_hessian_of_lagrangian_structure(
    ℰ::AbstractExponentialIntegrator,
    traj::NamedTrajectory,
)
    # Derived-Δt dynamics (Piccolo.jl#321): packed coordinates + the warp
    # column under a warp; historical structure otherwise. (Hermitian cells
    # only — the warp twin is Hermitian-specific.)
    ℰ isa HermitianExponentialIntegrator &&
        traj.warp !== nothing &&
        return _get_hessian_of_lagrangian_structure_warped(ℰ, traj)

    N = traj.N
    global_dim = traj.global_dim
    z_dim = traj.dim
    Z_dim = z_dim * N + global_dim
    μ∂²F = spzeros(Z_dim, Z_dim)

    # Get structure in canonical ordering (includes global blocks if global_dim > 0)
    μ∂²ℰ_canonical =
        hessian_structure(ℰ.x_dim, ℰ.u_dim, ℰ.global_dim; gauss_newton = ℰ.gauss_newton)

    # Build index mapping: canonical → trajectory (per-knot vars only)
    canonical_comps, traj_comps = build_hessian_index_mapping(ℰ, traj)
    knot_dim = ℰ.x_dim + ℰ.u_dim + 1

    # Map and place in full structure for all knot points
    for k = 1:(N-1)
        μ∂²F[slice(k, traj_comps, z_dim), slice(k, traj_comps, z_dim)] .+=
            μ∂²ℰ_canonical[canonical_comps, canonical_comps]
    end

    # Global structure: cross-terms + g-g block.
    # `g_can` is the canonical (global_names-order) block; `g_traj` maps each
    # canonical global to its Z column in the TRAJECTORY's global_components order
    # (see `_global_full_cols`) so values and structure agree under a reordering.
    if global_dim > 0
        g_can = collect((2*knot_dim+1):(2*knot_dim+global_dim))
        g_traj = _global_full_cols(ℰ, traj)
        for k = 1:(N-1)
            μ∂²F[slice(k, traj_comps, z_dim), g_traj] .+=
                μ∂²ℰ_canonical[canonical_comps, g_can]
            μ∂²F[g_traj, slice(k, traj_comps, z_dim)] .+=
                μ∂²ℰ_canonical[g_can, canonical_comps]
        end
        μ∂²F[g_traj, g_traj] .+= μ∂²ℰ_canonical[g_can, g_can]
    end

    return μ∂²F
end

# NOTE: evaluate! methods are defined per trajectory type in
# hermitian_exponential_integrator_ket.jl and hermitian_exponential_integrator_unitary.jl

@views function eval_jacobian(ℰ::AbstractExponentialIntegrator, traj::NamedTrajectory)
    N = traj.N
    x_dim = ℰ.x_dim
    z_dim = traj.dim
    F_dim = x_dim * (N - 1)
    Z_dim = z_dim * N + traj.global_dim

    # Fill preallocated structures in parallel
    Threads.@threads for k = 1:(N-1)
        jacobian!(ℰ.∂ℰs[k], ℰ, traj[k], traj[k+1], k)
    end

    # Build final Jacobian by slicing out the relevant var_comps from preallocated structures
    # Handle ensemble (multiple x_names) vs single
    x_comps_now = vcat([collect(traj.components[name]) for name in ℰ.x_names]...)
    u_comps_now = traj.components[ℰ.u_name]
    Δt_comp_now = traj.components[traj.timestep][1]
    var_comps_now = [x_comps_now; collect(u_comps_now); Δt_comp_now]

    ∂F = spzeros(F_dim, Z_dim)
    for k = 1:(N-1)
        ∂F[slice(k, x_dim), slice(k, var_comps_now, z_dim)] = ℰ.∂ℰs[k][:, var_comps_now]
        ∂F[slice(k, x_dim), slice(k, z_dim .+ var_comps_now, z_dim)] =
            ℰ.∂ℰs[k][:, ℰ.z_dim .+ var_comps_now]
    end

    return ∂F
end

function eval_hessian_of_lagrangian(
    ℰ::AbstractExponentialIntegrator,
    traj::NamedTrajectory,
    μ::AbstractVector,
)
    N = traj.N
    x_dim = ℰ.x_dim
    z_dim = traj.dim
    Z_dim = z_dim * N + traj.global_dim

    # Fill preallocated Hessian structures in parallel
    Threads.@threads for k = 1:(N-1)
        μₖ = μ[slice(k, x_dim)]
        hessian_of_lagrangian!(ℰ.μ∂²ℰs[k], ℰ, μₖ, traj[k], traj[k+1], k)
    end

    # Build index mapping: canonical → trajectory
    canonical_comps, traj_comps = build_hessian_index_mapping(ℰ, traj)

    # Assemble final Hessian from preallocated structures with index mapping
    # The filled matrices are already symmetric from hessian_of_lagrangian!, just take triu
    μ∂²F = spzeros(Z_dim, Z_dim)
    for k = 1:(N-1)
        μ∂²ℰ_triu = triu(ℰ.μ∂²ℰs[k])
        μ∂²F[slice(k, traj_comps, z_dim), slice(k, traj_comps, z_dim)] .=
            μ∂²ℰ_triu[canonical_comps, canonical_comps]
    end

    return μ∂²F
end

@testitem "HermitianExponentialIntegrator carries preallocated exp buffers (per-thread)" begin
    using LinearAlgebra
    using NamedTrajectories
    using Piccolo

    sys = OpenQuantumSystem(PAULIS.Z, [PAULIS.X, PAULIS.Y], [1.0, 1.0])
    N = 11
    times = collect(range(0, 1.0, length = N))
    pulse = LinearSplinePulse(zeros(2, N), times)
    qtraj = KetTrajectory(sys, pulse, ComplexF64[1, 0], ComplexF64[0, 1])

    ℰ = HermitianExponentialIntegrator(qtraj, N)

    # Per-thread buffer vectors (7 groups total — see Task 2 + thread-race fix).
    # Each is a Vector{...} with length Threads.maxthreadid() so that
    # Threads.@threads loops don't race on shared eigendecomposition scratch.
    @test hasproperty(ℰ, :H_bufs)
    @test hasproperty(ℰ, :λ_bufs)
    @test hasproperty(ℰ, :V_bufs)
    @test hasproperty(ℰ, :cis_diag_bufs)
    @test hasproperty(ℰ, :tmp_bufs)
    @test hasproperty(ℰ, :work_bufs)
    @test hasproperty(ℰ, :expG_bufs)

    # Buffer-vector length matches thread upper bound
    nthreads_alloc = length(ℰ.H_bufs)
    @test nthreads_alloc >= Threads.nthreads()
    @test nthreads_alloc >= Threads.maxthreadid()

    # Per-thread sizes match ketdim (= sys.levels)
    n = sys.levels
    for t = 1:nthreads_alloc
        @test size(ℰ.H_bufs[t]) == (n, n)
        @test size(ℰ.V_bufs[t]) == (n, n)
        @test length(ℰ.λ_bufs[t]) == n
        @test length(ℰ.cis_diag_bufs[t]) == n
        @test size(ℰ.tmp_bufs[t]) == (n, n)
        @test size(ℰ.work_bufs[t]) == (n, n)
        @test size(ℰ.expG_bufs[t]) == (2n, 2n)
    end
end
