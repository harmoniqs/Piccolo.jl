# ============================================================================ #
# Warp plumbing — derived-Δt support for the Hermitian exponential family
# (NamedTrajectories#161 time warp; Piccolo.jl#321, Phase 1b of the TR
# integration, mirroring DirectTrajOpt#150's WarpPlumbing shape).
#
# Under a time warp the trajectory's timestep rows are DERIVED: present as
# components (so `traj[k].timestep` keeps reading) but excluded from the packed
# decision vector. The packed layout of `vec(traj)` is
#
#     [ non-derived rows of knot 1 … knot N ] ++ [ global data ] ++ [ warp params ]
#
# with `n_params(warp)` trailing entries — for `GlobalScale`, the scalar
# duration `T`. Every integrator surface that addresses the packed vector
# (`eval_jacobian`, `eval_hessian_of_lagrangian`, the structure functions) must
# use these packed positions and expose the warp-parameter column with the
# EXACT chain entries `∂Δtₖ/∂θⱼ` (`ddeltats_dparams`, the NT#161 API — exact
# rational lattice weights for `GlobalScale`, generic ForwardDiff through
# `with_params` otherwise). Without a warp every helper below is unused: the
# warp-free code paths are bit-unchanged.
#
# NOTE: these are defined HERE, not imported from `DirectTrajOpt.WarpPlumbing`,
# because the pinned DTO 0.10.1 predates the warp release — the DTO-side
# Evaluator plumbing (post-#150) is only exercised by the cross-repo
# end-to-end smoke (two-repo dev env, see the #321 receipt).
# ============================================================================ #

"""
    _packed_knot_dim(traj)::Int

Per-knot width of the packed vector: all component rows except the derived
timestep. Without a warp this equals `traj.dim` (the helper is then unused).
"""
_packed_knot_dim(traj::NamedTrajectory) =
    traj.warp === nothing ? traj.dim : traj.dim - traj.dims[traj.timestep]

"""
    _packed_length(traj)::Int

Total length of the packed decision vector: `== length(vec(traj))`. Without a
warp this is the historical `traj.dim * traj.N + traj.global_dim`.
"""
_packed_length(traj::NamedTrajectory) = length(traj)

"""
    _packed_row_index(traj, k, r)::Int

Packed position of knot `k`'s datavec row `r` (1-based within the knot). The
derived timestep row has NO packed position — it is never decision data — and
asking for it errors. Without a warp this is the historical
`index(k, r, traj.dim)`.
"""
function _packed_row_index(traj::NamedTrajectory, k::Int, r::Int)
    if traj.warp === nothing
        return TrajectoryIndexingUtils.index(k, r, traj.dim)
    end
    zdim = _packed_knot_dim(traj)
    dt_first = first(traj.components[traj.timestep])
    dt_dim = traj.dims[traj.timestep]
    dt_first ≤ r < dt_first + dt_dim && error(
        "the derived timestep row has no packed index — it is never decision data " *
        "under a time warp (NamedTrajectories#161)",
    )
    return (k - 1) * zdim + (r ≤ dt_first ? r : r - dt_dim)
end

"""
    _packed_slice(traj, k, comps)::Vector{Int}

Packed positions of knot `k`'s component rows `comps` (datavec rows, 1-based
within a knot, none derived). Without a warp this is the historical
`slice(k, comps, traj.dim)` — bit-identical.
"""
_packed_slice(traj::NamedTrajectory, k::Int, comps) =
    [_packed_row_index(traj, k, r) for r in comps]

"""
    _warp_param_indices(traj)::Vector{Int}

Packed positions of the warp parameters (empty without a warp) — the trailing
block of `vec(traj)`. For `GlobalScale` this is the single duration variable
`T`.
"""
function _warp_param_indices(traj::NamedTrajectory)
    warp = traj.warp
    warp === nothing && return Int[]
    base = _packed_knot_dim(traj) * traj.N + traj.global_dim
    return [base + j for j = 1:n_params(warp)]
end

"""
    _warp_chain(traj)

Exact chain rule `∂Δtₖ/∂θⱼ` for the trajectory's warp as an `(N-1) × n_params`
matrix — `NamedTrajectories.ddeltats_dparams`, the NT#161 API: the
EXACT rational lattice weights `1/(N-1)` for `GlobalScale` (its override — no
floating subtraction), the generic ForwardDiff pass through `with_params`
(the Dual rebuild AD seam) for any other warp.
"""
_warp_chain(traj::NamedTrajectory) = ddeltats_dparams(traj.warp, traj.N)

# ---------------------------------------------------------------------------- #
# Shared warp-path assembly for the Hermitian exponential family. All three
# cells (Unitary/Ket/MultiKet) share the same canonical per-knot derivative
# layout `[xₖ, uₖ, Δtₖ, xₖ₊₁]` (+ globals), so the packed re-assembly — with
# the Δt column/row mapped onto the warp-parameter column through the exact
# chain rule — is written once here and dispatched to from each cell's
# `eval_jacobian` / `eval_hessian_of_lagrangian` / structure functions.
#
# Routing stays ENUMERATED: each cell's public API function guards
# `traj.warp !== nothing` and routes to its warp twin; no warp → the
# historical body, bit-unchanged.
# ---------------------------------------------------------------------------- #

"""
    _hermitian_warp_comps(ℰ, traj) -> (x_cols, u_cols, Δt_col, x1_cols, g_cols_full)

Per-knot matrix column indices for the state, control, derived-timestep,
next-state, and global blocks — CONSTRUCTION-time components (`ℰ.var_comps`,
`ℰ.z_dim`): the preallocated per-knot `∂ℰs` were laid out by
`jacobian_structure` on the integrator's construction trajectory, and the
eval-time trajectory may carry extra components (e.g. du/ddu added by the
smooth-pulse template) that shift the eval-time component rows.
"""
function _hermitian_warp_comps(ℰ::HermitianExponentialIntegrator, traj::NamedTrajectory)
    x_cols = collect(ℰ.var_comps[1:ℰ.x_dim])
    u_cols = collect(ℰ.var_comps[(ℰ.x_dim+1):(ℰ.x_dim+ℰ.u_dim)])
    Δt_col = ℰ.var_comps[end]
    x1_cols = ℰ.z_dim .+ x_cols
    g_cols_full = ℰ.global_dim > 0 ? _global_full_cols(ℰ, traj) : Int[]
    return x_cols, u_cols, Δt_col, x1_cols, g_cols_full
end

"""
    _hermitian_eval_comps(ℰ, traj) -> (x_cols, u_cols)

Eval-time trajectory component rows of the state and control blocks — the
TARGET side of the packed re-assembly (the packed positions of the eval-time
trajectory's own data).
"""
function _hermitian_eval_comps(ℰ::HermitianExponentialIntegrator, traj::NamedTrajectory)
    x_cols = vcat([collect(traj.components[name]) for name in ℰ.x_names]...)
    u_cols = collect(traj.components[ℰ.u_name])
    return x_cols, u_cols
end

"""
    _eval_jacobian_warped(ℰ, traj)

Warp-path Jacobian of the defects. Per knot the preallocated `∂ℰs[k]` is filled
by the cell's own `jacobian!` (unchanged — the derived Δtₖ is read from the
trajectory rows, which `sync_timesteps!` re-derives from the warp), and
the per-knot block is scattered into PACKED coordinates with the derived-Δt
column replaced by the EXACT warp column:

    ∂cₖ/∂θⱼ = (∂Δtₖ/∂θⱼ) · ∂cₖ/∂Δtₖ ,   ∂Δtₖ/∂θⱼ = ddeltats_dparams(warp, N)[k, j]

i.e. `wₙ · ∂cₙ/∂Δtₙ` per warp parameter — for `GlobalScale`, `∂cₖ/∂T =
∂cₖ/∂Δtₖ / (N-1)`. On the satisfied defect `∂cₙ/∂Δtₙ = −G(uₙ)·xₙ₊₁` exactly.
"""
function _eval_jacobian_warped(ℰ::HermitianExponentialIntegrator, traj::NamedTrajectory)
    sync_timesteps!(traj)
    N = traj.N
    x_dim = ℰ.x_dim
    F_dim = x_dim * (N - 1)
    Z_dim = _packed_length(traj)
    warp_cols = _warp_param_indices(traj)
    chain = _warp_chain(traj)   # (N-1) × n_params

    globals = extract_globals(ℰ, traj)
    Threads.@threads for k = 1:(N-1)
        jacobian!(ℰ.∂ℰs[k], ℰ, traj[k], traj[k+1], k, globals)
    end

    x_cols, u_cols, Δt_col, x1_cols, g_cols_full = _hermitian_warp_comps(ℰ, traj)
    x_cols_e, u_cols_e = _hermitian_eval_comps(ℰ, traj)

    ∂F = spzeros(F_dim, Z_dim)
    @inbounds for k = 1:(N-1)
        Mₖ = ℰ.∂ℰs[k]
        ∂F[slice(k, x_dim), _packed_slice(traj, k, x_cols_e)] = Mₖ[:, x_cols]
        ∂F[slice(k, x_dim), _packed_slice(traj, k, u_cols_e)] = Mₖ[:, u_cols]
        ∂F[slice(k, x_dim), _packed_slice(traj, k+1, x_cols_e)] = Mₖ[:, x1_cols]
        # the warp column: chain-rule products of the per-knot ∂cₖ/∂Δtₖ entries
        ∂F[slice(k, x_dim), warp_cols] = Mₖ[:, Δt_col] * transpose(chain[k, :])
    end
    if ℰ.global_dim > 0
        g_cols_local = (2*ℰ.z_dim+1):(2*ℰ.z_dim+ℰ.global_dim)
        @inbounds for k = 1:(N-1)
            ∂F[slice(k, x_dim), g_cols_full] = ℰ.∂ℰs[k][:, g_cols_local]
        end
    end
    return ∂F
end

"""
    _eval_hessian_of_lagrangian_warped(ℰ, traj, μ)

Warp-path Hessian of the Lagrangian. Per knot the canonical `μ∂²ℰs[k]` (layout
`[xₖ, uₖ, Δtₖ, xₖ₊₁]` + globals) is filled by the cell's own
`hessian_of_lagrangian!` (unchanged), then re-assembled in PACKED coordinates
with the canonical Δt row/column mapped onto the warp-parameter column through
the exact chain rule:

    ∂²(μ'cₖ)/∂v∂θⱼ = (∂Δtₖ/∂θⱼ) · ∂²(μ'cₖ)/∂v∂Δtₖ ,
    ∂²(μ'cₖ)/∂θⱼ∂θⱼ' = (∂Δtₖ/∂θⱼ)(∂Δtₖ/∂θⱼ') · ∂²(μ'cₖ)/∂Δtₖ² .

Covers the x-Δt / u-Δt / Δt-g cross-terms and the Δt-Δt term exactly as the
canonical structure declares them (u-Δt / Δt-g / Δt-Δt dropped in
Gauss–Newton mode, mirroring the warp-free path).
"""
function _eval_hessian_of_lagrangian_warped(
    ℰ::HermitianExponentialIntegrator,
    traj::NamedTrajectory,
    μ::AbstractVector,
)
    sync_timesteps!(traj)
    N = traj.N
    x_dim, u_dim = ℰ.x_dim, ℰ.u_dim
    knot_dim = x_dim + u_dim + 1
    Z_dim = _packed_length(traj)
    warp_cols = _warp_param_indices(traj)
    nθ = length(warp_cols)
    chain = _warp_chain(traj)
    Δt_can = knot_dim

    globals = extract_globals(ℰ, traj)
    Threads.@threads for k = 1:(N-1)
        μₖ = μ[slice(k, x_dim)]
        hessian_of_lagrangian!(ℰ.μ∂²ℰs[k], ℰ, μₖ, traj[k], traj[k+1], k, globals)
    end

    x_cols, u_cols, _, _, g_cols_full = _hermitian_warp_comps(ℰ, traj)
    x_cols_e, u_cols_e = _hermitian_eval_comps(ℰ, traj)

    # canonical rows/cols of the per-knot Hessian EXCLUDING Δt, in [x, u, x₊₁, g] order
    C_plain = [collect(1:(x_dim+u_dim)); collect((knot_dim+1):(knot_dim+x_dim))]
    C =
        ℰ.global_dim > 0 ? [C_plain; collect((2*knot_dim+1):(2*knot_dim+ℰ.global_dim))] :
        C_plain

    μ∂²F = spzeros(Z_dim, Z_dim)
    @inbounds for k = 1:(N-1)
        Mₖ = ℰ.μ∂²ℰs[k]
        Z = [
            _packed_slice(traj, k, x_cols_e);
            _packed_slice(traj, k, u_cols_e);
            _packed_slice(traj, k+1, x_cols_e);
            g_cols_full;
        ]
        μ∂²F[Z, Z] .+= Mₖ[C, C]
        for j = 1:nθ
            w = chain[k, j]
            iszero(w) && continue
            μ∂²F[Z, warp_cols[j]] .+= w .* Mₖ[C, Δt_can]
            μ∂²F[warp_cols[j], Z] .+= w .* Mₖ[Δt_can, C]
        end
        μ∂²F[warp_cols, warp_cols] .+=
            Mₖ[Δt_can, Δt_can] .* (chain[k, :] * transpose(chain[k, :]))
    end
    return triu(μ∂²F)
end

"""
    _get_jacobian_structure_warped(ℰ, traj)

Warp-path Jacobian STRUCTURE: packed two-knot windows (x, u; x₊₁) plus the
warp-parameter column (declared with ones — a conservative superset of the
exact `chain[k, j] · ∂cₖ/∂Δtₖ` nonzeros, safe for the solver's preallocation),
in the same coordinates the DTO Evaluator's packed path consumes.
"""
function _get_jacobian_structure_warped(
    ℰ::HermitianExponentialIntegrator,
    traj::NamedTrajectory,
)
    N = traj.N
    x_dim = ℰ.x_dim
    ∂F = spzeros(x_dim * (N - 1), _packed_length(traj))
    warp_cols = _warp_param_indices(traj)
    _, _, _, _, g_cols_full = _hermitian_warp_comps(ℰ, traj)
    x_cols_e, u_cols_e = _hermitian_eval_comps(ℰ, traj)

    for k = 1:(N-1)
        cols = [
            _packed_slice(traj, k, x_cols_e);
            _packed_slice(traj, k, u_cols_e);
            warp_cols;
            _packed_slice(traj, k+1, x_cols_e);
        ]
        ∂F[slice(k, x_dim), cols] .= 1.0
    end
    if ℰ.global_dim > 0
        for k = 1:(N-1)
            ∂F[slice(k, x_dim), g_cols_full] .= 1.0
        end
    end
    return ∂F
end

"""
    _get_hessian_of_lagrangian_structure_warped(ℰ, traj)

Warp-path Hessian-of-the-Lagrangian STRUCTURE: the canonical structure
(`hessian_structure` in `[xₖ, uₖ, Δtₖ, xₖ₊₁]` + globals ordering, Gauss–Newton
mode respected) re-assembled in packed coordinates with the canonical Δt
row/column mapped onto the warp-parameter column (conservative ones).
"""
function _get_hessian_of_lagrangian_structure_warped(
    ℰ::HermitianExponentialIntegrator,
    traj::NamedTrajectory,
)
    N = traj.N
    x_dim, u_dim, global_dim = ℰ.x_dim, ℰ.u_dim, ℰ.global_dim
    knot_dim = x_dim + u_dim + 1
    M_can = hessian_structure(x_dim, u_dim, global_dim; gauss_newton = ℰ.gauss_newton)
    warp_cols = _warp_param_indices(traj)
    nθ = length(warp_cols)
    Δt_can = knot_dim

    _, _, _, _, g_cols_full = _hermitian_warp_comps(ℰ, traj)
    x_cols_e, u_cols_e = _hermitian_eval_comps(ℰ, traj)
    C_plain = [collect(1:(x_dim+u_dim)); collect((knot_dim+1):(knot_dim+x_dim))]
    C =
        global_dim > 0 ? [C_plain; collect((2*knot_dim+1):(2*knot_dim+global_dim))] :
        C_plain

    μ∂²F = spzeros(_packed_length(traj), _packed_length(traj))
    for k = 1:(N-1)
        Z = [
            _packed_slice(traj, k, x_cols_e);
            _packed_slice(traj, k, u_cols_e);
            _packed_slice(traj, k+1, x_cols_e);
            g_cols_full;
        ]
        μ∂²F[Z, Z] .+= M_can[C, C]
        for j = 1:nθ
            μ∂²F[Z, warp_cols[j]] .+= 1.0
            μ∂²F[warp_cols[j], Z] .+= 1.0
        end
        μ∂²F[warp_cols, warp_cols] .+= 1.0
    end
    return μ∂²F
end
