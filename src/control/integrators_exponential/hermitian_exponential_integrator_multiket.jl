# ============================================================================ #
# HermitianExponentialIntegrator for MultiKetTrajectory
# 
# This integrator exploits the fact that all kets in an ensemble evolve under
# the same propagator. We compute exp(Δt*G) once per knot point and apply it to all
# kets, significantly reducing computation compared to separate integrators.
# ============================================================================ #

"""
    HermitianExponentialIntegrator(qtraj::MultiKetTrajectory, N::Int; kwargs...)

Construct a Hermitian exponential integrator from an MultiKetTrajectory.
Automatically builds the combined trajectory with separate state variables.

# Arguments
- `qtraj::MultiKetTrajectory`: The ensemble quantum trajectory
- `N::Int`: Number of knot points for discretization

# Keyword Arguments
- `global_names::Union{Vector{Symbol}, Nothing}=nothing`: Names of global (time-invariant) 
  variables to include in dynamics. If `nothing`, auto-detects from `sys.global_params`.
"""
function HermitianExponentialIntegrator(
    qtraj::MultiKetTrajectory,
    N::Int;
    global_names::Union{Vector{Symbol},Nothing} = nothing,
    gauss_newton::Bool = false,
    matrix_free::Bool = false,
    use_analytical::Bool = true,
)
    # Build trajectory WITH globals so structure functions get correct dimensions
    sys = get_system(qtraj)
    resolved_names = if isnothing(global_names)
        !isempty(sys.global_params) ? collect(keys(sys.global_params)) : Symbol[]
    else
        global_names
    end

    traj = if !isempty(resolved_names)
        global_data = Dict{Symbol,Vector{Float64}}()
        for name in resolved_names
            if hasproperty(sys, :global_params) && haskey(sys.global_params, name)
                global_data[name] = [sys.global_params[name]]
            else
                global_data[name] = [0.0]
            end
        end
        NamedTrajectory(qtraj, N; global_data = global_data)
    else
        NamedTrajectory(qtraj, N)
    end

    return HermitianExponentialIntegrator(
        qtraj,
        traj;
        global_names = global_names,
        gauss_newton = gauss_newton,
        matrix_free = matrix_free,
        use_analytical = use_analytical,
    )
end

"""
    HermitianExponentialIntegrator(qtraj::MultiKetTrajectory, traj::NamedTrajectory; kwargs...)

Construct a Hermitian exponential integrator for ensemble ket trajectory evolution.
The propagator is computed ONCE per knot point and applied to all kets.

# Keyword Arguments
- `global_names::Union{Vector{Symbol}, Nothing}=nothing`: Names of global (time-invariant) 
  variables to include in dynamics. If `nothing`, auto-detects from `sys.global_params`.
"""
# Per-member inner constructor (#408's sampling lane): all pieces
# pre-materialized. Delegates to the (qtraj, traj) inner method via the
# member-qtraj shim (which supplies get_system/state_names/drive_name only).
function _hermitian_exp_multiket(
    sys,
    x_names::Vector{Symbol},
    u::Symbol,
    traj::NamedTrajectory,
    global_names::Vector{Symbol};
    kwargs...,
)
    return HermitianExponentialIntegrator(
        _mk_exp_member_qtraj(sys, x_names, u),
        traj;
        global_names = global_names,
        kwargs...,
    )
end

struct _ExpMemberQTraj{S,X<:Vector{Symbol},U<:Symbol}
    sys::S
    x_names::X
    u::U
end
Piccolo.get_system(q::_ExpMemberQTraj) = q.sys
Piccolo.state_names(q::_ExpMemberQTraj) = q.x_names
Piccolo.drive_name(q::_ExpMemberQTraj) = q.u

_mk_exp_member_qtraj(sys, x_names, u) = _ExpMemberQTraj(sys, x_names, u)

function HermitianExponentialIntegrator(
    qtraj::Union{MultiKetTrajectory,_ExpMemberQTraj},
    traj::NamedTrajectory;
    global_names::Union{Vector{Symbol},Nothing} = nothing,
    gauss_newton::Bool = false,
    matrix_free::Bool = false,
    use_analytical::Bool = true,
)
    sys = get_system(qtraj)
    x_names = state_names(qtraj)
    u = drive_name(qtraj)
    n_kets = length(x_names)

    @assert traj.N > 1 "Trajectory must have at least two timesteps."
    @assert all(name in traj.names for name in x_names) "All ensemble state names must be in trajectory"

    # Auto-detect globals if not specified
    if isnothing(global_names)
        if !isempty(sys.global_params)
            global_names = collect(keys(sys.global_params))
        else
            global_names = Symbol[]
        end
    end

    global_dim = traj.global_dim
    ketdim = sys.levels

    # Dimensions for API - total constraint dimension per knot point for all kets
    single_ket_dim = 2 * ketdim  # Single ket dimension
    x_dim = n_kets * single_ket_dim  # Total per-knot constraint dimension
    u_dim = traj.dims[u]
    N = traj.N
    dim = x_dim * (N - 1)

    # Variables: [xₖ (all kets), uₖ, Δtₖ, xₖ₊₁ (all kets)]
    var_dim = 2 * x_dim + u_dim + 1

    # Use ensemble structure functions. Pass the resolved global_dim so the
    # template is sized correctly even when the caller-provided trajectory
    # hasn't had globals attached yet (they may arrive later via
    # SmoothPulseProblem(...; free_phase=true)).
    ∂ℰ_template = jacobian_structure(
        MultiKetTrajectory,
        x_names,
        u,
        ketdim,
        traj;
        global_dim = global_dim,
    )
    # For ensemble: pass individual state dimensions
    x_dims = [traj.dims[name] for name in x_names]
    μ∂²ℰ_template =
        hessian_structure(x_dims, u_dim, global_dim; gauss_newton = gauss_newton)

    # Preallocate one copy per knot point
    ∂ℰs = [copy(∂ℰ_template) for _ = 1:(N-1)]
    μ∂²ℰs = [copy(μ∂²ℰ_template) for _ = 1:(N-1)]

    # Build var_comps - all ket components + controls + time
    x_comps_all = vcat([collect(traj.components[name]) for name in x_names]...)
    u_comps = traj.components[u]
    Δt_comp = traj.components[traj.timestep][1]
    var_comps = [x_comps_all; collect(u_comps); Δt_comp]

    z_dim = traj.dim

    # Per-thread scratch for hessian aggregation (one slot per thread).
    # These are reset+accumulated inside hessian_of_lagrangian! and would race
    # if shared across Threads.@threads-scheduled knot points.
    # When gauss_newton=true, these buffers aren't exercised so we skip alloc.
    nthr = Threads.maxthreadid()
    ∂u∂Δt_bufs = gauss_newton ? nothing : [zeros(u_dim) for _ = 1:nthr]
    ∂²u_bufs = gauss_newton ? nothing : [zeros(u_dim, u_dim) for _ = 1:nthr]

    # Preallocate per-thread exp_eigen! buffers (7 groups, sized to ketdim).
    # See HermitianExponentialIntegrator docstring for thread-safety rationale.
    bufs = _alloc_exp_eigen_bufs(ketdim, nthr)

    # Analytic Daleckii–Krein setup (#202): constant affine drive directions over the
    # extended control [controls; globals], plus per-thread DK scratch. `H_dirs` is
    # `nothing` for nonlinear-drive Hamiltonians (→ ForwardDiff fallback).
    H_dirs = _build_affine_directions(u_ -> sys.H(u_, 0.0), u_dim + global_dim, ketdim)
    dk_bufs = _alloc_dk_bufs(ketdim, nthr)
    dk_so_bufs = _alloc_dk_so_bufs(ketdim, nthr)

    return HermitianExponentialIntegrator{MultiKetTrajectory}(
        u_ -> sys.G(u_, 0.0),
        u_ -> sys.H(u_, 0.0),
        x_names,
        u,
        x_dim,
        u_dim,
        var_dim,
        dim,
        ketdim,
        ∂ℰs,
        μ∂²ℰs,
        z_dim,
        var_comps,
        nothing,     # Id (not used for multiket)
        ∂u∂Δt_bufs,  # Per-thread buffers for u-Δt cross term
        ∂²u_bufs,    # Per-thread buffers for u-u Hessian block
        global_names,
        global_dim,
        bufs.H_bufs,
        bufs.λ_bufs,
        bufs.V_bufs,
        bufs.cis_diag_bufs,
        bufs.tmp_bufs,
        bufs.work_bufs,
        bufs.expG_bufs,
        gauss_newton,
        matrix_free,
        use_analytical,
        H_dirs,
        dk_bufs.dk_ws_bufs,
        dk_bufs.dk_dΦ_bufs,
        dk_bufs.dk_vin_bufs,
        dk_bufs.dk_vout_bufs,
        dk_so_bufs.dk_so_ws_bufs,
        dk_so_bufs.dk_μmat_bufs,
        dk_so_bufs.dk_psi_bufs,
        dk_so_bufs.dk_y_bufs,
        dk_so_bufs.dk_Hc_bufs,
    )
end

# ============================================================================ #
# Call operator for MultiKetTrajectory
# ============================================================================ #

@views function (ℰ::HermitianExponentialIntegrator{MultiKetTrajectory})(
    δₖ::AbstractVector,
    zₖ::KnotPoint,
    zₖ₊₁::KnotPoint,
    k::Int,
    globals::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    aₖ = zₖ[ℰ.u_name]
    Δtₖ = zₖ.timestep

    # Build extended control with globals
    uₖ = build_extended_control(aₖ, globals)

    # Per-thread buffer indexing — see HermitianExponentialIntegrator docstring.
    tid = Threads.threadid()
    H_buf = ℰ.H_bufs[tid]
    V_buf = ℰ.V_bufs[tid]
    λ_buf = ℰ.λ_bufs[tid]
    cis_diag_buf = ℰ.cis_diag_bufs[tid]
    tmp_buf = ℰ.tmp_bufs[tid]
    work_buf = ℰ.work_bufs[tid]
    expG_buf = ℰ.expG_bufs[tid]

    # Compute propagator ONCE for all kets
    Gₖ = ℰ.G(uₖ)
    copyto!(H_buf, ℰ.H(uₖ))
    exp_eigen!(expG_buf, H_buf, V_buf, λ_buf, cis_diag_buf, tmp_buf, work_buf, Δtₖ)
    expGₖ = expG_buf

    # Apply to each ket
    single_dim = single_state_dim(ℰ)
    @inbounds for (i, name) in enumerate(ℰ.x_names)
        ψ̃ₖ = zₖ[name]
        ψ̃ₖ₊₁ = zₖ₊₁[name]
        δ_range = slice(i, single_dim)
        # δₖ = ψ̃ₖ₊₁ - exp(Δt*G)*ψ̃ₖ
        mul!(δₖ[δ_range], expGₖ, ψ̃ₖ, -1.0, 0.0)
        δₖ[δ_range] .+= ψ̃ₖ₊₁
    end

    return nothing
end

# ============================================================================ #
# Jacobian methods for MultiKetTrajectory
# ============================================================================ #

@views function jacobian!(
    ∂ℰ::AbstractMatrix,
    ℰ::HermitianExponentialIntegrator{MultiKetTrajectory},
    zₖ::KnotPoint,
    zₖ₊₁::KnotPoint,
    k::Int,
    globals::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    aₖ = zₖ[ℰ.u_name]
    Δtₖ = zₖ.timestep

    # Build extended control with globals
    uₖ = build_extended_control(aₖ, globals)

    # Per-thread buffer indexing — see HermitianExponentialIntegrator docstring.
    tid = Threads.threadid()
    H_buf = ℰ.H_bufs[tid]
    V_buf = ℰ.V_bufs[tid]
    λ_buf = ℰ.λ_bufs[tid]
    cis_diag_buf = ℰ.cis_diag_bufs[tid]
    tmp_buf = ℰ.tmp_bufs[tid]
    work_buf = ℰ.work_bufs[tid]
    expG_buf = ℰ.expG_bufs[tid]

    Gₖ = ℰ.G(uₖ)
    copyto!(H_buf, ℰ.H(uₖ))
    exp_eigen!(expG_buf, H_buf, V_buf, λ_buf, cis_diag_buf, tmp_buf, work_buf, Δtₖ)
    expGₖ = expG_buf

    # Get component indices
    u_comps = zₖ.components[ℰ.u_name]
    Δt_comp = zₖ.components.Δt[1]
    z_dim_matrix = size(∂ℰ, 2) ÷ 2
    single_dim = single_state_dim(ℰ)
    n_kets = length(ℰ.x_names)
    ketdim = ℰ.ketdim

    # Analytic Daleckii–Krein path for affine drives (#202): build the divided-difference
    # matrix M once per knot from the eigenbasis exp_eigen! just produced, then reuse it
    # across kets/directions. Falls back to ForwardDiff for nonlinear drives (H_dirs===nothing).
    use_dk = ℰ.use_analytical && !isnothing(ℰ.H_dirs)
    if use_dk
        dk_ws = ℰ.dk_ws_bufs[tid]
        dk_dΦ = ℰ.dk_dΦ_bufs[tid]
        dk_vin = ℰ.dk_vin_bufs[tid]
        dk_vout = ℰ.dk_vout_bufs[tid]
        dk_divided_difference!(dk_ws.M, λ_buf, Δtₖ)
        u_dirs = view(ℰ.H_dirs, 1:(ℰ.u_dim))
        g_dirs = view(ℰ.H_dirs, (ℰ.u_dim+1):(ℰ.u_dim+ℰ.global_dim))
    end

    # Process each ket - shared propagator but different state indices
    @inbounds for (i, name) in enumerate(ℰ.x_names)
        ψ̃ₖ = zₖ[name]
        x_comps = zₖ.components[name]
        δ_rows = slice(i, single_dim)
        expGₖψ̃ₖ = expGₖ * ψ̃ₖ  # Cache for ∂Δt term

        # ∂ψ̃ₖℰ: -exp(ΔtG)
        ∂ℰ[δ_rows, x_comps] .= .-expGₖ

        # ∂aₖℰ: control derivatives — analytic DK or ForwardDiff witness
        if use_dk
            _iso_vec_to_complex!(dk_vin, ψ̃ₖ, ketdim)
            _dk_fill_iso_jac_block!(
                ∂ℰ[δ_rows, u_comps],
                V_buf,
                dk_ws.M,
                u_dirs,
                dk_vin,
                dk_dΦ,
                dk_ws.tmp1,
                dk_ws.tmp2,
                dk_vout,
                ketdim,
            )
        else
            ForwardDiff.jacobian!(
                ∂ℰ[δ_rows, u_comps],
                a -> -expv(Δtₖ, ℰ.G(build_extended_control(a, globals)), ψ̃ₖ),
                aₖ,
            )
        end

        # ∂Δtₖℰ: -G * exp(ΔtG) * ψ̃ₖ
        mul!(∂ℰ[δ_rows, Δt_comp], Gₖ, expGₖψ̃ₖ, -1.0, 0.0)

        # ∂gℰ: Global variable derivatives
        if !isnothing(globals) && ℰ.global_dim > 0
            g_cols = (2*ℰ.z_dim+1):(2*ℰ.z_dim+ℰ.global_dim)
            if use_dk
                _iso_vec_to_complex!(dk_vin, ψ̃ₖ, ketdim)
                _dk_fill_iso_jac_block!(
                    ∂ℰ[δ_rows, g_cols],
                    V_buf,
                    dk_ws.M,
                    g_dirs,
                    dk_vin,
                    dk_dΦ,
                    dk_ws.tmp1,
                    dk_ws.tmp2,
                    dk_vout,
                    ketdim,
                )
            else
                ForwardDiff.jacobian!(
                    ∂ℰ[δ_rows, g_cols],
                    g -> -expv(Δtₖ, ℰ.G(vcat(aₖ, g)), ψ̃ₖ),
                    globals,
                )
            end
        end

        # ∂ψ̃ₖ₊₁ℰ: Identity (already in sparsity pattern)
    end

    return nothing
end

function jacobian_structure(
    ::Type{MultiKetTrajectory},
    x_names::Vector{Symbol},
    u_name::Symbol,
    ketdim::Int,
    traj::NamedTrajectory;
    global_dim::Union{Nothing,Int} = nothing,
)
    single_dim = 2 * ketdim  # Single ket dimension in isomorphism
    n_kets = length(x_names)
    x_dim = n_kets * single_dim  # Total constraint rows
    u_dim = traj.dims[u_name]
    z_dim = traj.dim
    # Default to the trajectory's current global_dim, but let the caller
    # override when constructing an integrator before globals are attached
    # to the trajectory (e.g., when passing `global_names` kwarg).
    global_dim = isnothing(global_dim) ? traj.global_dim : global_dim
    u_comps = traj.components[u_name]
    Δt_comp = traj.components[traj.timestep][1]

    ∂ℰ = spzeros(x_dim, 2z_dim + global_dim)

    for (i, name) in enumerate(x_names)
        x_comps = traj.components[name]
        δ_rows = slice(i, single_dim)

        # ∂ψ̃ₖ₊₁ℰ: Identity (per ket)
        ∂ℰ[δ_rows, z_dim .+ x_comps] = I(single_dim)

        # ∂ψ̃ₖℰ: Dense structure for ket
        ∂ℰ[δ_rows, x_comps] .= 1.0

        # ∂aₖℰ: Each ket contributes to control derivatives
        ∂ℰ[δ_rows, u_comps] .= 1.0

        # ∂Δtₖℰ
        ∂ℰ[δ_rows, Δt_comp] .= 1.0

        # ∂gℰ: Global variable derivatives
        if global_dim > 0
            g_cols = (2z_dim+1):(2z_dim+global_dim)
            ∂ℰ[δ_rows, g_cols] .= 1.0
        end
    end

    return ∂ℰ
end

function get_jacobian_structure(
    ℰ::HermitianExponentialIntegrator{MultiKetTrajectory},
    traj::NamedTrajectory,
)
    N = traj.N
    x_dim = ℰ.x_dim  # Total constraint dim per knot point (all kets)
    z_dim = traj.dim
    F_dim = x_dim * (N - 1)
    Z_dim = z_dim * N + traj.global_dim
    ∂F = spzeros(F_dim, Z_dim)

    single_dim = single_state_dim(ℰ)
    u_comps = collect(traj.components[ℰ.u_name])
    Δt_comp = traj.components[traj.timestep][1]

    for k = 1:(N-1)
        for (i, name) in enumerate(ℰ.x_names)
            x_comps = collect(traj.components[name])
            # Constraint rows for ket i at knot point k
            δ_rows = (k - 1) * x_dim .+ slice(i, single_dim)

            # ∂ψ̃ₖℰ: ket i state at knot point k (block-diagonal)
            ∂F[δ_rows, slice(k, x_comps, z_dim)] .= 1.0

            # ∂ψ̃ₖ₊₁ℰ: ket i state at knot point k+1 (block-diagonal)
            ∂F[δ_rows, slice(k, z_dim .+ x_comps, z_dim)] .= 1.0

            # ∂aₖℰ: controls (shared across all kets)
            ∂F[δ_rows, slice(k, u_comps, z_dim)] .= 1.0

            # ∂Δtₖℰ: timestep (shared)
            ∂F[δ_rows, slice(k, [Δt_comp], z_dim)] .= 1.0
        end
    end

    # Global columns (dynamics globals only, at their trajectory-order positions —
    # see `_global_full_cols`; fixes permuted / not-first global columns).
    if ℰ.global_dim > 0
        g_cols_full = _global_full_cols(ℰ, traj)
        for k = 1:(N-1)
            ∂F[slice(k, x_dim), g_cols_full] .= 1.0
        end
    end

    return ∂F
end

# ============================================================================ #
# API methods for MultiKetTrajectory
# ============================================================================ #

@views function evaluate!(
    δ::AbstractVector,
    ℰ::HermitianExponentialIntegrator{MultiKetTrajectory},
    traj::NamedTrajectory,
)
    # Extract globals once (constant across knot points)
    globals = extract_globals(ℰ, traj)

    Threads.@threads for k = 1:(traj.N-1)
        δₖ = δ[slice(k, ℰ.x_dim)]
        ℰ(δₖ, traj[k], traj[k+1], k, globals)
    end
    return nothing
end

@views function eval_jacobian(
    ℰ::HermitianExponentialIntegrator{MultiKetTrajectory},
    traj::NamedTrajectory,
)
    # Matrix-free opt-in (#205): return the CPU-threaded `CPUJacobianOp` instead of
    # assembling the dense-blocked sparse dynamics Jacobian. `_is_matrix_free_block=true`
    # ⇒ the Altissimo backend routes JVP/VJP straight into the threaded `mul!`/`mul!'`.
    # Gated on `matrix_free` AND an affine drive (`H_dirs !== nothing`); nonlinear drives
    # fall through to the retained sparse path (out of scope for the op, #205). The
    # default `matrix_free=false` keeps the assembled sparse Jacobian for the
    # Ipopt/MadNLP KKT-factorization path that shares this method.
    if ℰ.matrix_free && !isnothing(ℰ.H_dirs)
        return matrix_free_jacobian_op(ℰ, traj)
    end

    N = traj.N
    z_dim = traj.dim
    F_dim = ℰ.x_dim * (N - 1)
    Z_dim = z_dim * N + traj.global_dim

    # Extract globals once (constant across knot points)
    globals = extract_globals(ℰ, traj)

    # Fill preallocated structures in parallel
    Threads.@threads for k = 1:(N-1)
        jacobian!(ℰ.∂ℰs[k], ℰ, traj[k], traj[k+1], k, globals)
    end

    # Canonical var_comps (from construction-time trajectory) for indexing into ∂ℰs[k]
    var_comps_canonical = ℰ.var_comps

    # Eval-time var_comps for placing results into ∂F (traj may have extra columns like du/ddu)
    x_comps_all = vcat([collect(traj.components[name]) for name in ℰ.x_names]...)
    u_comps_now = traj.components[ℰ.u_name]
    Δt_comp_now = traj.components[traj.timestep][1]
    var_comps_now = [x_comps_all; collect(u_comps_now); Δt_comp_now]

    ∂F = spzeros(F_dim, Z_dim)
    @inbounds for k = 1:(N-1)
        ∂F[slice(k, ℰ.x_dim), slice(k, var_comps_now, z_dim)] =
            ℰ.∂ℰs[k][:, var_comps_canonical]
        ∂F[slice(k, ℰ.x_dim), slice(k, z_dim .+ var_comps_now, z_dim)] =
            ℰ.∂ℰs[k][:, ℰ.z_dim .+ var_comps_canonical]
    end

    # Global columns: only for globals the integrator was built with (dynamics globals),
    # NOT for globals added later (e.g., free-phase variables from SplinePulseProblem)
    if ℰ.global_dim > 0
        g_cols_local = (2*ℰ.z_dim+1):(2*ℰ.z_dim+ℰ.global_dim)
        g_cols_full = _global_full_cols(ℰ, traj)   # trajectory-order (permuted-∂c/∂θ fix)
        @inbounds for k = 1:(N-1)
            ∂F[slice(k, ℰ.x_dim), g_cols_full] = ℰ.∂ℰs[k][:, g_cols_local]
        end
    end

    return ∂F
end

@views function eval_hessian_of_lagrangian(
    ℰ::HermitianExponentialIntegrator{MultiKetTrajectory},
    traj::NamedTrajectory,
    μ::AbstractVector,
)
    N = traj.N
    x_dim = ℰ.x_dim
    u_dim = ℰ.u_dim
    z_dim = traj.dim
    Z_dim = z_dim * N + traj.global_dim

    # Extract globals once (constant across knot points)
    globals = extract_globals(ℰ, traj)

    # Fill preallocated Hessian structures in parallel
    Threads.@threads for k = 1:(N-1)
        μₖ = μ[slice(k, x_dim)]
        hessian_of_lagrangian!(ℰ.μ∂²ℰs[k], ℰ, μₖ, traj[k], traj[k+1], k, globals)
    end

    # Build index mapping inline for ensemble case
    knot_dim = x_dim + u_dim + 1  # canonical_hessian_knot_dim

    # Canonical indices
    x_can_k = 1:x_dim
    u_can_k = (x_dim+1):(x_dim+u_dim)
    Δt_can_k = x_dim + u_dim + 1
    x_can_k1 = knot_dim .+ (1:x_dim)

    # Trajectory indices - ensemble has multiple state names
    x_traj_k = vcat([collect(traj.components[name]) for name in ℰ.x_names]...)
    u_traj_k = collect(traj.components[ℰ.u_name])
    Δt_traj_k = traj.components[traj.timestep][1]
    x_traj_k1 = z_dim .+ x_traj_k

    canonical_comps = [collect(x_can_k); collect(u_can_k); Δt_can_k; collect(x_can_k1)]
    traj_comps = [x_traj_k; u_traj_k; Δt_traj_k; x_traj_k1]

    # Assemble final Hessian from preallocated structures with index mapping
    # Note: we map the full symmetric matrix first, THEN extract upper triangle
    # This is necessary because the canonical-to-trajectory index mapping can
    # swap upper/lower triangle positions
    μ∂²F = spzeros(Z_dim, Z_dim)
    @inbounds for k = 1:(N-1)
        μ∂²F[slice(k, traj_comps, z_dim), slice(k, traj_comps, z_dim)] .=
            ℰ.μ∂²ℰs[k][canonical_comps, canonical_comps]
    end

    # Global blocks: only for globals the integrator was built with (dynamics globals)
    if ℰ.global_dim > 0
        g_can = (2*knot_dim+1):(2*knot_dim+ℰ.global_dim)
        g_traj = _global_full_cols(ℰ, traj)   # trajectory-order (permuted-∂c/∂θ fix)
        @inbounds for k = 1:(N-1)
            # Per-knot vars × globals (upper triangle since traj indices < global indices)
            μ∂²F[slice(k, traj_comps, z_dim), g_traj] .=
                ℰ.μ∂²ℰs[k][canonical_comps, collect(g_can)]
            # Globals × per-knot vars (lower triangle, for symmetry)
            μ∂²F[g_traj, slice(k, traj_comps, z_dim)] .=
                ℰ.μ∂²ℰs[k][collect(g_can), canonical_comps]
            # g-g block (accumulated across knot points)
            μ∂²F[g_traj, g_traj] .+= ℰ.μ∂²ℰs[k][collect(g_can), collect(g_can)]
        end
    end

    # Return upper triangle (symmetric matrix, only upper needed by Ipopt)
    return triu(μ∂²F)
end

function get_hessian_of_lagrangian_structure(
    ℰ::HermitianExponentialIntegrator{MultiKetTrajectory},
    traj::NamedTrajectory,
)
    N = traj.N
    x_dim = ℰ.x_dim
    u_dim = ℰ.u_dim
    z_dim = traj.dim
    Z_dim = z_dim * N + traj.global_dim
    μ∂²F = spzeros(Z_dim, Z_dim)

    # Get structure in canonical ordering — use ℰ.global_dim (dynamics globals only)
    x_dims = [traj.dims[name] for name in ℰ.x_names]
    μ∂²ℰ_canonical =
        hessian_structure(x_dims, u_dim, ℰ.global_dim; gauss_newton = ℰ.gauss_newton)

    # Build index mapping inline for ensemble case
    knot_dim = x_dim + u_dim + 1

    # Canonical indices
    x_can_k = 1:x_dim
    u_can_k = (x_dim+1):(x_dim+u_dim)
    Δt_can_k = x_dim + u_dim + 1
    x_can_k1 = knot_dim .+ (1:x_dim)

    # Trajectory indices - ensemble has multiple state names
    x_traj_k = vcat([collect(traj.components[name]) for name in ℰ.x_names]...)
    u_traj_k = collect(traj.components[ℰ.u_name])
    Δt_traj_k = traj.components[traj.timestep][1]
    x_traj_k1 = z_dim .+ x_traj_k

    canonical_comps = [collect(x_can_k); collect(u_can_k); Δt_can_k; collect(x_can_k1)]
    traj_comps = [x_traj_k; u_traj_k; Δt_traj_k; x_traj_k1]

    # Map and place in full structure for all knot points
    for k = 1:(N-1)
        μ∂²F[slice(k, traj_comps, z_dim), slice(k, traj_comps, z_dim)] .+=
            μ∂²ℰ_canonical[canonical_comps, canonical_comps]
    end

    # Global structure: only for globals the integrator was built with (dynamics globals)
    if ℰ.global_dim > 0
        g_can = collect((2*knot_dim+1):(2*knot_dim+ℰ.global_dim))
        g_traj = _global_full_cols(ℰ, traj)   # trajectory-order (permuted-∂c/∂θ fix)
        for k = 1:(N-1)
            μ∂²F[slice(k, traj_comps, z_dim), g_traj] .+=
                μ∂²ℰ_canonical[canonical_comps, g_can]
            μ∂²F[g_traj, slice(k, traj_comps, z_dim)] .+=
                μ∂²ℰ_canonical[g_can, canonical_comps]
        end
        # g-g block (appears once, accumulated N-1 times in structure)
        μ∂²F[g_traj, g_traj] .+= μ∂²ℰ_canonical[g_can, g_can]
    end

    return μ∂²F
end

# ============================================================================ #
# Hessian methods for MultiKetTrajectory
# ============================================================================ #

"""
Hessian structure for ensemble with multiple state variables.
x_dims is a vector of individual state dimensions.
"""
function hessian_structure(
    x_dims::Vector{Int},
    u_dim::Int,
    global_dim::Int = 0;
    gauss_newton::Bool = false,
)
    x_dim = sum(x_dims)  # Total state dimension
    knot_dim = x_dim + u_dim + 1

    total_dim = 2 * knot_dim + global_dim
    μ∂²ℰ = spzeros(total_dim, total_dim)

    # Canonical indices for knot k
    x_comps = 1:x_dim
    u_comps = (x_dim+1):(x_dim+u_dim)
    Δt_comp = x_dim + u_dim + 1

    # Cross-terms (always kept)
    μ∂²ℰ[x_comps, u_comps] .= 1.0
    μ∂²ℰ[x_comps, Δt_comp] .= 1.0

    # Parameter-parameter blocks (dropped in GN)
    if !gauss_newton
        μ∂²ℰ[u_comps, Δt_comp] .= 1.0
        μ∂²ℰ[u_comps, u_comps] .= 1.0
        μ∂²ℰ[Δt_comp, Δt_comp] = 1.0
    end

    # Global variable blocks
    if global_dim > 0
        g_comps = (2*knot_dim+1):(2*knot_dim+global_dim)
        μ∂²ℰ[x_comps, g_comps] .= 1.0    # x-g (cross-term, always kept)
        if !gauss_newton
            μ∂²ℰ[u_comps, g_comps] .= 1.0    # u-g
            μ∂²ℰ[Δt_comp:Δt_comp, g_comps] .= 1.0  # Δt-g
            μ∂²ℰ[g_comps, g_comps] .= 1.0    # g-g
        end
    end

    return sparse(Symmetric(μ∂²ℰ))
end

@views function hessian_of_lagrangian!(
    μ∂²ℰ::SparseMatrixCSC,
    ℰ::HermitianExponentialIntegrator{MultiKetTrajectory},
    μₖ::AbstractVector,
    zₖ::KnotPoint,
    zₖ₊₁::KnotPoint,
    k::Int,
    globals::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    x_dim = ℰ.x_dim
    u_dim = ℰ.u_dim
    global_dim = ℰ.global_dim
    single_dim = single_state_dim(ℰ)
    has_globals = !isnothing(globals) && global_dim > 0

    aₖ = zₖ[ℰ.u_name]
    Δtₖ = zₖ.timestep

    # Build extended control with globals
    uₖ = build_extended_control(aₖ, globals)

    # Per-thread buffer indexing — see HermitianExponentialIntegrator docstring.
    tid = Threads.threadid()
    H_buf = ℰ.H_bufs[tid]
    V_buf = ℰ.V_bufs[tid]
    λ_buf = ℰ.λ_bufs[tid]
    cis_diag_buf = ℰ.cis_diag_bufs[tid]
    tmp_buf = ℰ.tmp_bufs[tid]
    work_buf = ℰ.work_bufs[tid]
    expG_buf = ℰ.expG_bufs[tid]

    Gₖ = ℰ.G(uₖ)
    copyto!(H_buf, ℰ.H(uₖ))
    exp_eigen!(expG_buf, H_buf, V_buf, λ_buf, cis_diag_buf, tmp_buf, work_buf, Δtₖ)
    expGₖ = expG_buf
    GₖexpGₖ = Gₖ * expGₖ
    ketdim = ℰ.ketdim

    # Analytic Daleckii–Krein path for the always-computed Gauss-Newton cross-terms
    # (x-u, x-g): these are first derivatives of Φ†, so the slice-① first-order kernel
    # suffices. Build M once per knot from the eigenbasis exp_eigen! produced; the
    # `!gauss_newton` parameter-parameter blocks below stay on ForwardDiff (slice ④).
    # Falls back to ForwardDiff for nonlinear drives (H_dirs===nothing).
    use_dk = ℰ.use_analytical && !isnothing(ℰ.H_dirs)
    if use_dk
        dk_ws = ℰ.dk_ws_bufs[tid]
        dk_dΦ = ℰ.dk_dΦ_bufs[tid]
        dk_vin = ℰ.dk_vin_bufs[tid]
        dk_vout = ℰ.dk_vout_bufs[tid]
        dk_divided_difference!(dk_ws.M, λ_buf, Δtₖ)
        u_dirs = view(ℰ.H_dirs, 1:u_dim)
        g_dirs = view(ℰ.H_dirs, (u_dim+1):(u_dim+global_dim))
    end

    # Canonical component indices for knot k
    knot_dim = x_dim + u_dim + 1
    x_can = 1:x_dim
    u_can = (x_dim+1):(x_dim+u_dim)
    Δt_can = x_dim + u_dim + 1

    # Zero out first
    @inbounds μ∂²ℰ.nzval .= 0.0

    # Use per-thread preallocated buffers for aggregation (only needed for exact Hessian).
    # Per-thread access is required because this function can be called from a
    # Threads.@threads loop over knots, and a shared buffer would race.
    if !ℰ.gauss_newton
        ∂u∂Δt_buf = ℰ.∂u∂Δt_bufs[tid]
        ∂u∂Δt_buf .= 0.0  # Reset for accumulation
        # Analytic exact-Hessian parameter-parameter scratch (#203, slice ④). Per-thread
        # second-order DK buffers + a thread-local m×m real block accumulator (m = full
        # affine-parameter count over [controls; globals]); reused across kets.
        if use_dk
            dk_so_ws = ℰ.dk_so_ws_bufs[tid]
            dk_μmat = ℰ.dk_μmat_bufs[tid]
            dk_psi = ℰ.dk_psi_bufs[tid]
            dk_y = ℰ.dk_y_bufs[tid]
            dk_Hc = ℰ.dk_Hc_bufs[tid]
            dk_B = zeros(u_dim + global_dim, u_dim + global_dim)
        end
    end

    # Global buffers (small, allocated per call)
    if has_globals
        g_can = (2*knot_dim+1):(2*knot_dim+global_dim)
        if !ℰ.gauss_newton
            full_ag_dim = u_dim + global_dim
            ∂²ag_buf = zeros(full_ag_dim, full_ag_dim)
            ∂Δt∂g_buf = zeros(global_dim)
        end
    end

    # Aggregate contributions from all kets
    @inbounds for (i, name) in enumerate(ℰ.x_names)
        ψ̃ₖ = zₖ[name]
        μₖ_i = μₖ[((i-1)*single_dim+1):(i*single_dim)]
        x_can_i = ((i-1)*single_dim+1):(i*single_dim)

        # === Cross-terms (always computed) ===
        # DK cross-terms apply the ADJOINT of ∂ₚΦ to μₖ_i (the derivative routine
        # differentiates -iso(Φ†)·μ), mirroring the ForwardDiff-∘-expv witness.
        if use_dk
            _iso_vec_to_complex!(dk_vin, μₖ_i, ketdim)
        end

        # μₖ∂ψ̃ₖ∂aₖℰ
        if use_dk
            _dk_fill_iso_jac_block!(
                μ∂²ℰ[x_can_i, u_can],
                V_buf,
                dk_ws.M,
                u_dirs,
                dk_vin,
                dk_dΦ,
                dk_ws.tmp1,
                dk_ws.tmp2,
                dk_vout,
                ketdim;
                adjoint_op = true,
            )
        else
            ForwardDiff.jacobian!(
                μ∂²ℰ[x_can_i, u_can],
                a -> -expv(Δtₖ, ℰ.G(build_extended_control(a, globals))', μₖ_i),
                aₖ,
            )
        end

        # μₖ∂ψ̃ₖ∂Δtₖℰ
        mul!(μ∂²ℰ[x_can_i, Δt_can], GₖexpGₖ', μₖ_i, -1.0, 0.0)

        if has_globals
            # μₖ∂ψ̃ₖ∂gℰ - x-g cross-term (per-ket, direct write)
            if use_dk
                _dk_fill_iso_jac_block!(
                    μ∂²ℰ[x_can_i, g_can],
                    V_buf,
                    dk_ws.M,
                    g_dirs,
                    dk_vin,
                    dk_dΦ,
                    dk_ws.tmp1,
                    dk_ws.tmp2,
                    dk_vout,
                    ketdim;
                    adjoint_op = true,
                )
            else
                ForwardDiff.jacobian!(
                    μ∂²ℰ[x_can_i, g_can],
                    g -> -expv(Δtₖ, ℰ.G(vcat(aₖ, g))', μₖ_i),
                    globals,
                )
            end
        end

        # === Parameter-parameter blocks (skipped in GN) ===
        if !ℰ.gauss_newton
            GₖexpGₖψ̃ₖ = GₖexpGₖ * ψ̃ₖ

            if use_dk
                # ── Analytic exact-Hessian parameter-parameter blocks (#203, slice ④) ──
                # Every p-p block is  -∂²/∂p∂q[Re(χ_i† Φ ψ_i)]  with χ_i = iso⁻¹(μₖ_i),
                # ψ_i = iso⁻¹(ψ̃ₖ); μ_mat = χ_i ψ_i† realises the Re tr(μ†·) pairing the
                # second-order kernel expects (matching the ForwardDiff-∘-expv witness).
                # `dk_vin` already holds χ_i (set in the DK cross-term block above).
                _iso_vec_to_complex!(dk_psi, ψ̃ₖ, ketdim)
                dk_μmat .= dk_vin .* dk_psi'      # χ_i ψ_i†
                # y = Φ ψ_i = V (cis(-Δt λ) ∘ (V† ψ_i))   (reuse exp_eigen! eigenbasis)
                mul!(dk_vout, V_buf', dk_psi)
                dk_vout .*= cis_diag_buf
                mul!(dk_y, V_buf, dk_vout)
                mul!(dk_Hc, H_buf, dk_vin)        # H χ_i

                # u-u / u-g / g-g via the three-index second divided difference:
                #   μ∂²ℰ[p,q] += -Re tr(μ_mat† D²Φ[Hₚ,H_q]) = -dk_B[p,q]
                dk_second_order_block!(dk_B, dk_so_ws, V_buf, λ_buf, ℰ.H_dirs, dk_μmat, Δtₖ)
                μ∂²ℰ[u_can, u_can] .-= dk_B[1:u_dim, 1:u_dim]
                if has_globals
                    μ∂²ℰ[u_can, g_can] .-= dk_B[1:u_dim, (u_dim+1):(u_dim+global_dim)]
                    μ∂²ℰ[g_can, g_can] .-=
                        dk_B[(u_dim+1):(u_dim+global_dim), (u_dim+1):(u_dim+global_dim)]
                end

                # Parameter-Δt cross terms (u-Δt, Δt-g). Closed form:
                #   ∂²Φ/∂aₚ∂Δt = -iHₚΦ - iH ∂ₚΦ,  ∂ₚΦ = V(M∘(V†HₚV))V† (first-order DK),
                # so the block value is  -Re(χ_i† ∂²Φ/∂aₚ∂Δt ψ_i)
                #                       = Re(i·[χ_i† Hₚ(Φψ_i) + (Hχ_i)† (∂ₚΦ)ψ_i]).
                for p = 1:(u_dim+global_dim)
                    Hₚ = ℰ.H_dirs[p]
                    dk_apply!(dk_dΦ, V_buf, dk_ws.M, Hₚ, dk_ws.tmp1, dk_ws.tmp2)
                    mul!(dk_vout, Hₚ, dk_y)          # Hₚ (Φ ψ_i)
                    s = dot(dk_vin, dk_vout)         # χ_i† Hₚ Φ ψ_i
                    mul!(dk_vout, dk_dΦ, dk_psi)     # (∂ₚΦ) ψ_i
                    s += dot(dk_Hc, dk_vout)         # (H χ_i)† (∂ₚΦ) ψ_i
                    val = real(im * s)
                    if p <= u_dim
                        μ∂²ℰ[u_can[p], Δt_can] += val
                    else
                        μ∂²ℰ[Δt_can, g_can[p-u_dim]] += val
                    end
                end
            else
                # μₖ∂aₖ∂Δtₖℰ - aggregate from all kets
                ForwardDiff.gradient!(
                    ∂u∂Δt_buf,
                    a ->
                        -μₖ_i' *
                        ℰ.G(build_extended_control(a, globals)) *
                        expv(Δtₖ, ℰ.G(build_extended_control(a, globals)), ψ̃ₖ),
                    aₖ,
                )
                μ∂²ℰ[u_can, Δt_can] .+= ∂u∂Δt_buf

                if has_globals
                    # Full [a,g] Hessian → u-u, u-g, g-g blocks (accumulated)
                    ForwardDiff.hessian!(
                        ∂²ag_buf,
                        ag -> -μₖ_i'expv(Δtₖ, ℰ.G(ag), ψ̃ₖ),
                        vcat(aₖ, globals),
                    )
                    μ∂²ℰ[u_can, u_can] .+= ∂²ag_buf[1:u_dim, 1:u_dim]
                    μ∂²ℰ[u_can, g_can] .+= ∂²ag_buf[1:u_dim, (u_dim+1):end]
                    μ∂²ℰ[g_can, g_can] .+= ∂²ag_buf[(u_dim+1):end, (u_dim+1):end]

                    # μₖ∂Δtₖ∂gℰ - Δt-g block (accumulated)
                    ForwardDiff.gradient!(
                        ∂Δt∂g_buf,
                        g -> -μₖ_i' * ℰ.G(vcat(aₖ, g)) * expv(Δtₖ, ℰ.G(vcat(aₖ, g)), ψ̃ₖ),
                        globals,
                    )
                    for j = 1:global_dim
                        μ∂²ℰ[Δt_can, g_can[j]] += ∂Δt∂g_buf[j]
                    end
                else
                    # μₖ∂²aₖℰ - aggregate from all kets (no globals case); per-thread buf
                    ∂²u_buf_t = ℰ.∂²u_bufs[tid]
                    ForwardDiff.hessian!(
                        ∂²u_buf_t,
                        a -> -μₖ_i'expv(Δtₖ, ℰ.G(build_extended_control(a, globals)), ψ̃ₖ),
                        aₖ,
                    )
                    μ∂²ℰ[u_can, u_can] .+= ∂²u_buf_t
                end
            end

            # μₖ∂²Δtₖℰ - aggregate from all kets (closed form, analytic in both paths)
            μ∂²ℰ[Δt_can, Δt_can] += -μₖ_i' * Gₖ * GₖexpGₖψ̃ₖ
        end
    end

    # Symmetrize: canonical→trajectory index mapping can swap triangle positions
    @inbounds for j = 1:u_dim
        for i = 1:x_dim
            μ∂²ℰ[u_can[j], x_can[i]] = μ∂²ℰ[x_can[i], u_can[j]]
        end
    end
    @inbounds for i = 1:x_dim
        μ∂²ℰ[Δt_can, x_can[i]] = μ∂²ℰ[x_can[i], Δt_can]
    end
    if !ℰ.gauss_newton
        @inbounds for j = 1:u_dim
            μ∂²ℰ[Δt_can, u_can[j]] = μ∂²ℰ[u_can[j], Δt_can]
        end
    end

    # Symmetrize global blocks
    if has_globals
        @inbounds for j = 1:global_dim
            for i = 1:x_dim
                μ∂²ℰ[g_can[j], x_can[i]] = μ∂²ℰ[x_can[i], g_can[j]]
            end
            if !ℰ.gauss_newton
                for i = 1:u_dim
                    μ∂²ℰ[g_can[j], u_can[i]] = μ∂²ℰ[u_can[i], g_can[j]]
                end
                μ∂²ℰ[g_can[j], Δt_can] = μ∂²ℰ[Δt_can, g_can[j]]
            end
        end
    end

    return nothing
end

# ============================================================================ #
# Tests for MultiKetTrajectory
# ============================================================================ #

@testitem "HermitianExponentialIntegrator with MultiKetTrajectory" begin
    using DirectTrajOpt
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra

    T = 10.0
    N = 10
    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])

    # Create initial and goal states for ensemble
    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]

    initials = [ψ0, ψ1]  # |0⟩, |1⟩
    goals = [ψ1, ψ0]     # |1⟩, |0⟩

    # Create ensemble trajectory
    ensemble_qtraj = MultiKetTrajectory(sys, initials, goals, T)

    # Create single HermitianExponentialIntegrator for ensemble
    integrator = HermitianExponentialIntegrator(ensemble_qtraj, N)

    @test integrator isa HermitianExponentialIntegrator{MultiKetTrajectory}
    @test length(integrator.x_names) == 2

    # Verify dimensions
    traj = NamedTrajectory(ensemble_qtraj, N)
    @test integrator.x_dim == 2 * 2 * 2  # 2 kets * 2 ketdim * 2 (isomorphism)
    @test integrator.dim == integrator.x_dim * (traj.N - 1)
end





@testitem "HermitianExponentialIntegrator Jacobian sparsity for MultiKetTrajectory" begin
    using DirectTrajOpt
    using Piccolo
    using NamedTrajectories
    using SparseArrays
    using LinearAlgebra

    T = 1.0
    N = 10
    sys = QuantumSystem(GATES.Z, [GATES.X], [1.0])
    ketdim = 2
    n_kets = 2
    x_dim = 2 * ketdim * n_kets  # 8 total
    u_dim = 1

    # Create states for ensemble
    ψ_inits = [[1.0 + 0.0im, 0.0 + 0.0im], [0.0 + 0.0im, 1.0 + 0.0im]]
    ψ_goals = [[0.0 + 0.0im, 1.0 + 0.0im], [1.0 + 0.0im, 0.0 + 0.0im]]

    qtraj = MultiKetTrajectory(sys, ψ_inits, ψ_goals, T)
    integrator = HermitianExponentialIntegrator(qtraj, N)

    traj = NamedTrajectory(qtraj, N)
    ∂F = Piccolissimo.get_jacobian_structure(integrator, traj)

    # Expected sparsity for MultiKetTrajectory with block-diagonal structure:
    # - Each ket evolves independently → n_kets blocks of (2ketdim × 2ketdim)
    # - Parameters: u_dim + 1 (controls + Δt)
    single_x_dim = 2 * ketdim
    n_params = u_dim + 1  # 2

    # Block-diagonal: each ket has separate (2ketdim × 2ketdim) state blocks
    state_nnz_per_k = 2 * n_kets * single_x_dim^2  # 2 * 2 * 16 = 64
    param_nnz_per_k = x_dim * n_params  # 8 * 2 = 16
    expected_block_diagonal_nnz = (N - 1) * (state_nnz_per_k + param_nnz_per_k)  # 9 * 80 = 720

    # Dense (no structure exploitation)
    dense_state_nnz_per_k = 2 * x_dim * x_dim  # 128
    expected_dense_nnz = (N - 1) * (dense_state_nnz_per_k + param_nnz_per_k)

    actual_nnz = nnz(∂F)
    println("  HermitianExponentialIntegrator MultiKetTrajectory Jacobian nnz: $actual_nnz")
    println("  Block-diagonal expected: $expected_block_diagonal_nnz")
    println("  Dense expected: $expected_dense_nnz")

    @test actual_nnz ≤ expected_block_diagonal_nnz
    @test actual_nnz < expected_dense_nnz
    println("  $(round(actual_nnz / expected_dense_nnz * 100, digits=1))% of dense")
end

@testitem "HermitianExponentialIntegrator{MultiKetTrajectory} forward pass minimal alloc" begin
    using LinearAlgebra
    using NamedTrajectories
    using Piccolo
    using BenchmarkTools

    sys = QuantumSystem(PAULIS.Z, [PAULIS.X, PAULIS.Y], [1.0, 1.0])
    N = 11
    times = collect(range(0, 1.0, length = N))
    pulse = LinearSplinePulse(zeros(2, N), times)
    ψ0 = ComplexF64[1, 0]
    ψ1 = ComplexF64[0, 1]
    qtraj = MultiKetTrajectory(sys, pulse, [ψ0, ψ1], [ψ1, ψ0])
    traj = NamedTrajectory(qtraj, N)

    ℰ = HermitianExponentialIntegrator(qtraj, N)
    δ = zeros(ℰ.x_dim)
    z1 = traj[1]
    z2 = traj[2]

    ℰ(δ, z1, z2, 1)  # warmup

    # Regression guard for Task 6 of plan-20260417-063000-nonhermitian-exp-integrator.md.
    # The forward path calls exp_eigen! with preallocated struct buffers (one shared
    # propagator across all kets). Residual allocation comes from ℰ.H(uₖ) building a
    # fresh matrix each call; same limitation as Tasks 4/5. Threshold is larger than
    # the UnitaryTrajectory bound because the multiket forward path iterates over kets
    # and indexes zₖ[name] per ket (one small slice view per ket).
    allocs = @ballocated $ℰ($δ, $z1, $z2, 1)
    @test allocs < 8192
end

# ============================================================================ #
# NonlinearDrive Tests
# ============================================================================ #

@testitem "HermitianExponentialIntegrator MultiKetTrajectory with NonlinearDrive" begin
    using DirectTrajOpt
    using Piccolissimo
    using Piccolo
    using SparseArrays
    using NamedTrajectories
    using LinearAlgebra

    include("../../../test/test_utils.jl")

    # 1-qubit system: H = σz + u₁·σx + u₁²·σz
    drives = AbstractDrive[
        LinearDrive(sparse(ComplexF64.(PAULIS.X)), 1),
        NonlinearDrive(PAULIS.Z, u -> u[1]^2),
    ]
    drive_bounds = [1.0]

    sys = QuantumSystem(PAULIS.Z, drives, drive_bounds)

    @test has_nonlinear_drives(sys.H_drives)

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]

    initials = [ψ0, ψ1]
    goals = [ψ1, ψ0]

    T = 1.0
    N = 6
    ensemble_qtraj = MultiKetTrajectory(sys, initials, goals, T)

    integrator = HermitianExponentialIntegrator(ensemble_qtraj, N)
    @test integrator isa HermitianExponentialIntegrator{MultiKetTrajectory}
    @test length(integrator.x_names) == 2

    traj = NamedTrajectory(ensemble_qtraj, N)
    test_integrator(integrator, traj; atol = 1e-3, show_hessian_diff = true)
end



@testitem "HermitianExponentialIntegrator MultiKetTrajectory with NonlinearDrive and global variables" begin
    using DirectTrajOpt
    using DirectTrajOpt: BoundsConstraint
    using Piccolissimo
    using Piccolo
    using SparseArrays
    using LinearAlgebra
    using NamedTrajectories

    include("../../../test/test_utils.jl")

    # System with nonlinear drive AND global parameter:
    #   H = u₁·σx + u₁²·σz, with optimizable global δ
    δ_init = 0.01

    drives = AbstractDrive[
        LinearDrive(sparse(ComplexF64.(PAULIS.X)), 1),
        NonlinearDrive(PAULIS.Z, u -> u[1]^2),
    ]
    drive_bounds = [1.0]

    sys = QuantumSystem(PAULIS.Z, drives, drive_bounds; global_params = (δ = δ_init,))

    @test has_nonlinear_drives(sys.H_drives)

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]

    initials = [ψ0, ψ1]
    goals = [ψ1, ψ0]

    T = 10.0
    N = 15
    ensemble_qtraj = MultiKetTrajectory(sys, initials, goals, T)

    integrator = HermitianExponentialIntegrator(ensemble_qtraj, N)

    @test integrator.global_names == [:δ]
    @test integrator.global_dim == 1

    traj = NamedTrajectory(ensemble_qtraj, N)
    test_integrator(integrator, traj; atol = 1e-3, show_hessian_diff = true)

    # Full optimization with global bounds
    δ_bound = 0.5
    qcp = SmoothPulseProblem(
        ensemble_qtraj,
        N;
        Q = 100.0,
        R = 1e-2,
        integrator = integrator,
        global_bounds = Dict(:δ => δ_bound),
    )

    solve!(qcp; max_iter = 200, print_level = 0, verbose = false)

    result_traj = get_trajectory(qcp)
    δ_opt = result_traj.global_data[result_traj.global_components[:δ]][1]
    @test isfinite(δ_opt)
    @test δ_opt >= -δ_bound - 1e-5
    @test δ_opt <= δ_bound + 1e-5
    println("  MultiKet NonlinearDrive + global: δ_init=$δ_init, δ_opt=$δ_opt")
end

# ============================================================================ #
# Gauss-Newton Hessian Tests
# ============================================================================ #

@testitem "HermitianExponentialIntegrator MultiKetTrajectory Gauss-Newton Hessian" begin
    using DirectTrajOpt
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using SparseArrays

    include("../../../test/test_utils.jl")

    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    initials = [ψ0, ψ1]
    goals = [ψ1, ψ0]

    T = 1.0
    N = 6
    ensemble_qtraj = MultiKetTrajectory(sys, initials, goals, T)
    traj = NamedTrajectory(ensemble_qtraj, N)

    # Test 1: GN integrator constructs correctly
    ℰ_gn = HermitianExponentialIntegrator(ensemble_qtraj, N; gauss_newton = true)
    @test ℰ_gn isa HermitianExponentialIntegrator{MultiKetTrajectory}
    @test ℰ_gn.gauss_newton == true
    # Field names are plural `∂u∂Δt_bufs` / `∂²u_bufs` after the per-thread
    # buffer refactor in 307bf6a. For GN the entire vector is `nothing`; for
    # exact-Hessian the entries are allocated per-thread (see merged docstring
    # on HermitianExponentialIntegrator).
    @test ℰ_gn.∂u∂Δt_bufs === nothing
    @test ℰ_gn.∂²u_bufs === nothing

    # Test 2: GN sparsity structure has fewer nonzeros than exact
    ℰ_exact = HermitianExponentialIntegrator(ensemble_qtraj, N; gauss_newton = false)
    H_gn = get_hessian_of_lagrangian_structure(ℰ_gn, traj)
    H_exact = get_hessian_of_lagrangian_structure(ℰ_exact, traj)
    @test nnz(H_gn) < nnz(H_exact)
    println("  GN nnz=$(nnz(H_gn)), Exact nnz=$(nnz(H_exact))")

    # Test 3: GN Hessian cross-terms match FiniteDiff
    test_integrator(ℰ_gn, traj; gauss_newton = true, atol = 1e-3, show_hessian_diff = true)

    # Test 4: GN (p,p) blocks are zero
    μ = rand(ℰ_gn.dim)
    μ∂²f_gn = eval_hessian_of_lagrangian(ℰ_gn, traj, μ)
    μ∂²f_exact = eval_hessian_of_lagrangian(ℰ_exact, traj, μ)
    gn_mask = triu(H_gn) .!= 0
    exact_mask = triu(H_exact) .!= 0
    pp_mask = exact_mask .& .!gn_mask
    @test any(pp_mask)
    @test all(triu(μ∂²f_gn)[pp_mask] .== 0.0)
    @test !all(triu(μ∂²f_exact)[pp_mask] .== 0.0)
    println("  (p,p) block entries: $(count(pp_mask)), all zero in GN: true")
end



# ============================================================================ #
# Analytic Daleckii–Krein dense-wiring tests (Piccolissimo.jl#202, slice ③).
# The multiket integrator is the comoving-OC driver; these guard the analytic
# Jacobian / Gauss-Newton Hessian cross-term path against the retained ForwardDiff
# witness, per-thread safety, the flag toggle, and the nonlinear-drive fallback.
# ============================================================================ #

@testitem "MultiKet DK Jacobian matches ForwardDiff witness — affine drive [#202 AC1]" begin
    using DirectTrajOpt
    using Piccolissimo
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using Random

    Random.seed!(202_001)
    T = 2.0
    N = 8
    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    qtraj = MultiKetTrajectory(sys, [ψ0, ψ1], [ψ1, ψ0], T)

    ℰ_dk = HermitianExponentialIntegrator(qtraj, N; use_analytical = true)
    ℰ_fd = HermitianExponentialIntegrator(qtraj, N; use_analytical = false)

    # Affine drive ⇒ DK path available; the flag records the toggle.
    @test ℰ_dk.H_dirs !== nothing
    @test length(ℰ_dk.H_dirs) == 2          # two linear drives, no globals
    @test ℰ_dk.use_analytical == true
    @test ℰ_fd.use_analytical == false

    traj = NamedTrajectory(qtraj, N)
    traj.datavec .= randn(length(traj.datavec))   # generic states + controls

    Jdk = eval_jacobian(ℰ_dk, traj)
    Jfd = eval_jacobian(ℰ_fd, traj)
    relerr = norm(Jdk - Jfd) / norm(Jfd)
    println("  [#202 AC1] Jacobian DK-vs-ForwardDiff rel err = $relerr")
    @test relerr < 1e-9
end

@testitem "MultiKet DK Gauss-Newton Hessian cross-terms match ForwardDiff [#202 AC2]" begin
    using DirectTrajOpt
    using Piccolissimo
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using Random

    Random.seed!(202_002)
    T = 2.0
    N = 7

    # --- no globals: exercises the x-u cross-term block ---
    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    qtraj = MultiKetTrajectory(sys, [ψ0, ψ1], [ψ1, ψ0], T)

    ℰ_dk =
        HermitianExponentialIntegrator(qtraj, N; use_analytical = true, gauss_newton = true)
    ℰ_fd = HermitianExponentialIntegrator(
        qtraj,
        N;
        use_analytical = false,
        gauss_newton = true,
    )

    traj = NamedTrajectory(qtraj, N)
    traj.datavec .= randn(length(traj.datavec))
    μ = randn(ℰ_dk.dim)

    Hdk = eval_hessian_of_lagrangian(ℰ_dk, traj, μ)
    Hfd = eval_hessian_of_lagrangian(ℰ_fd, traj, μ)
    relerr = norm(Hdk - Hfd) / norm(Hfd)
    println("  [#202 AC2] GN Hessian (x-u) DK-vs-ForwardDiff rel err = $relerr")
    @test relerr < 1e-9

    # --- with globals: additionally exercises the x-g cross-term block ---
    H = (u, t) -> (u[3] + u[4]) * GATES.Z + u[1] * GATES.X + u[2] * GATES.Y
    gsys = QuantumSystem(
        H,
        [1.0, 1.0];
        time_dependent = true,
        global_params = (δ = 0.2, ω = 1.0),
    )
    gqtraj = MultiKetTrajectory(gsys, [ψ0, ψ1], [ψ1, ψ0], T)
    ℰg_dk = HermitianExponentialIntegrator(
        gqtraj,
        N;
        use_analytical = true,
        gauss_newton = true,
    )
    ℰg_fd = HermitianExponentialIntegrator(
        gqtraj,
        N;
        use_analytical = false,
        gauss_newton = true,
    )
    @test ℰg_dk.H_dirs !== nothing
    @test length(ℰg_dk.H_dirs) == 4          # 2 controls + 2 globals
    @test ℰg_dk.global_dim == 2

    gtraj = NamedTrajectory(gqtraj, N)
    gtraj.datavec .= randn(length(gtraj.datavec))
    gtraj.global_data .= randn(length(gtraj.global_data))
    μg = randn(ℰg_dk.dim)

    Hgdk = eval_hessian_of_lagrangian(ℰg_dk, gtraj, μg)
    Hgfd = eval_hessian_of_lagrangian(ℰg_fd, gtraj, μg)
    relerr_g = norm(Hgdk - Hgfd) / norm(Hgfd)
    println("  [#202 AC2] GN Hessian (x-u & x-g) DK-vs-ForwardDiff rel err = $relerr_g")
    @test relerr_g < 1e-9
end

@testitem "MultiKet DK path is thread-safe: parallel == serial [#202 AC3]" begin
    using DirectTrajOpt
    using Piccolissimo
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using Random

    Random.seed!(202_003)
    T = 2.0
    N = 12

    # Affine system with globals so several blocks (x, u, Δt, g) are exercised.
    H = (u, t) -> (u[3] + u[4]) * GATES.Z + u[1] * GATES.X + u[2] * GATES.Y
    sys = QuantumSystem(
        H,
        [1.0, 1.0];
        time_dependent = true,
        global_params = (δ = 0.2, ω = 1.0),
    )
    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    qtraj = MultiKetTrajectory(sys, [ψ0, ψ1], [ψ1, ψ0], T)

    ℰ = HermitianExponentialIntegrator(qtraj, N; use_analytical = true, gauss_newton = true)
    @test ℰ.H_dirs !== nothing

    traj = NamedTrajectory(qtraj, N)
    traj.datavec .= randn(length(traj.datavec))
    traj.global_data .= randn(length(traj.global_data))
    globals = Piccolissimo.extract_globals(ℰ, traj)
    μ = randn(ℰ.dim)

    println("  [#202 AC3] running on $(Threads.nthreads()) thread(s)")

    # Jacobian: threaded fill of the per-knot blocks vs a serial fill — a race on the
    # new per-thread DK scratch would corrupt some knots under Threads.nthreads() > 1.
    Threads.@threads for k = 1:(N-1)
        Piccolissimo.jacobian!(ℰ.∂ℰs[k], ℰ, traj[k], traj[k+1], k, globals)
    end
    J_threaded = [copy(ℰ.∂ℰs[k]) for k = 1:(N-1)]
    for k = 1:(N-1)
        Piccolissimo.jacobian!(ℰ.∂ℰs[k], ℰ, traj[k], traj[k+1], k, globals)
    end
    J_serial = [copy(ℰ.∂ℰs[k]) for k = 1:(N-1)]
    @test all(J_threaded[k] == J_serial[k] for k = 1:(N-1))

    # Hessian: same threaded-vs-serial check on the per-knot Hessian blocks.
    Threads.@threads for k = 1:(N-1)
        μₖ = μ[slice(k, ℰ.x_dim)]
        Piccolissimo.hessian_of_lagrangian!(
            ℰ.μ∂²ℰs[k],
            ℰ,
            μₖ,
            traj[k],
            traj[k+1],
            k,
            globals,
        )
    end
    H_threaded = [copy(ℰ.μ∂²ℰs[k]) for k = 1:(N-1)]
    for k = 1:(N-1)
        μₖ = μ[slice(k, ℰ.x_dim)]
        Piccolissimo.hessian_of_lagrangian!(
            ℰ.μ∂²ℰs[k],
            ℰ,
            μₖ,
            traj[k],
            traj[k+1],
            k,
            globals,
        )
    end
    H_serial = [copy(ℰ.μ∂²ℰs[k]) for k = 1:(N-1)]
    @test all(H_threaded[k] == H_serial[k] for k = 1:(N-1))

    # Determinism of the public threaded assembly across repeated calls.
    @test eval_jacobian(ℰ, traj) == eval_jacobian(ℰ, traj)
    @test eval_hessian_of_lagrangian(ℰ, traj, μ) == eval_hessian_of_lagrangian(ℰ, traj, μ)
end

@testitem "MultiKet DK path is correct vs FiniteDiff oracle [#202 AC4]" begin
    using DirectTrajOpt
    using Piccolissimo
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using Random

    include("../../../test/test_utils.jl")

    Random.seed!(202_004)
    T = 2.0
    N = 6
    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    qtraj = MultiKetTrajectory(sys, [ψ0, ψ1], [ψ1, ψ0], T)
    traj = NamedTrajectory(qtraj, N)

    # The analytic (flag-on) path must agree with an INDEPENDENT FiniteDiff oracle,
    # not merely with the ForwardDiff witness — Jacobian and Gauss-Newton Hessian.
    ℰ_dk =
        HermitianExponentialIntegrator(qtraj, N; use_analytical = true, gauss_newton = true)
    test_integrator(ℰ_dk, traj; gauss_newton = true, atol = 1e-4, show_hessian_diff = true)
end

@testitem "MultiKet nonlinear drive falls back to ForwardDiff [#202 AC5]" begin
    using DirectTrajOpt
    using Piccolissimo
    using Piccolo
    using SparseArrays
    using NamedTrajectories
    using LinearAlgebra
    using Random

    include("../../../test/test_utils.jl")

    Random.seed!(202_005)
    # H = σz + u₁·σx + u₁²·σz — nonlinear in the control.
    drives = AbstractDrive[
        LinearDrive(sparse(ComplexF64.(PAULIS.X)), 1),
        NonlinearDrive(PAULIS.Z, u -> u[1]^2),
    ]
    sys = QuantumSystem(PAULIS.Z, drives, [1.0])
    @test has_nonlinear_drives(sys.H_drives)

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    T = 1.0
    N = 6
    qtraj = MultiKetTrajectory(sys, [ψ0, ψ1], [ψ1, ψ0], T)

    # Even with the analytic flag ON, a nonlinear drive has no constant direction
    # matrix ⇒ H_dirs is nothing ⇒ the DK path is unavailable and ForwardDiff runs.
    ℰ_on = HermitianExponentialIntegrator(qtraj, N; use_analytical = true)
    ℰ_off = HermitianExponentialIntegrator(qtraj, N; use_analytical = false)
    @test ℰ_on.H_dirs === nothing
    @test ℰ_on.use_analytical == true

    traj = NamedTrajectory(qtraj, N)
    traj.datavec .= randn(length(traj.datavec))

    # Flag-on and flag-off produce IDENTICAL results (both take the ForwardDiff path).
    @test eval_jacobian(ℰ_on, traj) == eval_jacobian(ℰ_off, traj)

    # And the fallback path is still correct against the FiniteDiff oracle.
    traj2 = NamedTrajectory(qtraj, N)
    test_integrator(ℰ_on, traj2; atol = 1e-3, show_hessian_diff = true)
end

# ============================================================================ #
# Analytic Daleckii–Krein exact-Hessian tests (Piccolissimo.jl#203, slice ④).
# These extend slice ③'s multiket Hessian surface to the `!gauss_newton`
# parameter-parameter blocks (u-u, u-Δt, Δt-Δt, and the global u-g / g-g / Δt-g
# blocks), now analytic via the second-order (three-index second divided
# difference) kernel + the closed-form parameter-Δt cross terms. The retained
# ForwardDiff-∘-expv path is the witness; an independent FiniteDiff oracle and a
# near-degenerate (confluent) fixture guard correctness.
# ============================================================================ #

@testitem "MultiKet DK exact-Hessian p-p blocks match ForwardDiff witness — affine [#203 AC1]" begin
    using DirectTrajOpt
    using Piccolissimo
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using SparseArrays
    using Random

    Random.seed!(203_001)
    T = 2.0
    N = 7

    # --- no globals: exercises the u-u, u-Δt, Δt-Δt parameter-parameter blocks ---
    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    qtraj = MultiKetTrajectory(sys, [ψ0, ψ1], [ψ1, ψ0], T)

    ℰ_dk = HermitianExponentialIntegrator(
        qtraj,
        N;
        use_analytical = true,
        gauss_newton = false,
    )
    ℰ_fd = HermitianExponentialIntegrator(
        qtraj,
        N;
        use_analytical = false,
        gauss_newton = false,
    )
    @test ℰ_dk.H_dirs !== nothing
    @test ℰ_dk.gauss_newton == false
    @test ℰ_dk.use_analytical == true

    traj = NamedTrajectory(qtraj, N)
    traj.datavec .= randn(length(traj.datavec))
    μ = randn(ℰ_dk.dim)

    Hdk = eval_hessian_of_lagrangian(ℰ_dk, traj, μ)
    Hfd = eval_hessian_of_lagrangian(ℰ_fd, traj, μ)
    relerr = norm(Hdk - Hfd) / norm(Hfd)
    println("  [#203 AC1] exact-Hessian (no globals) DK-vs-ForwardDiff rel err = $relerr")
    @test relerr < 1e-8

    # Isolate the parameter-parameter blocks this slice wires: positions present in
    # the exact structure but ABSENT from the Gauss-Newton structure (slice ③).
    ℰ_gn = HermitianExponentialIntegrator(qtraj, N; gauss_newton = true)
    H_exact_struct = get_hessian_of_lagrangian_structure(ℰ_dk, traj)
    H_gn_struct = get_hessian_of_lagrangian_structure(ℰ_gn, traj)
    pp_mask = (triu(H_exact_struct) .!= 0) .& .!(triu(H_gn_struct) .!= 0)
    @test any(pp_mask)
    ppdk = triu(Hdk)[pp_mask]
    ppfd = triu(Hfd)[pp_mask]
    @test norm(ppdk - ppfd) / norm(ppfd) < 1e-8
    @test !all(ppdk .== 0.0)   # the p-p blocks are actually populated (not silently skipped)

    # --- with globals: adds the u-g, g-g, and Δt-g parameter-parameter blocks ---
    H = (u, t) -> (u[3] + u[4]) * GATES.Z + u[1] * GATES.X + u[2] * GATES.Y
    gsys = QuantumSystem(
        H,
        [1.0, 1.0];
        time_dependent = true,
        global_params = (δ = 0.2, ω = 1.0),
    )
    gqtraj = MultiKetTrajectory(gsys, [ψ0, ψ1], [ψ1, ψ0], T)
    ℰg_dk = HermitianExponentialIntegrator(
        gqtraj,
        N;
        use_analytical = true,
        gauss_newton = false,
    )
    ℰg_fd = HermitianExponentialIntegrator(
        gqtraj,
        N;
        use_analytical = false,
        gauss_newton = false,
    )
    @test ℰg_dk.H_dirs !== nothing
    @test length(ℰg_dk.H_dirs) == 4          # 2 controls + 2 globals
    @test ℰg_dk.global_dim == 2

    gtraj = NamedTrajectory(gqtraj, N)
    gtraj.datavec .= randn(length(gtraj.datavec))
    gtraj.global_data .= randn(length(gtraj.global_data))
    μg = randn(ℰg_dk.dim)

    Hgdk = eval_hessian_of_lagrangian(ℰg_dk, gtraj, μg)
    Hgfd = eval_hessian_of_lagrangian(ℰg_fd, gtraj, μg)
    relerr_g = norm(Hgdk - Hgfd) / norm(Hgfd)
    println("  [#203 AC1] exact-Hessian (globals) DK-vs-ForwardDiff rel err = $relerr_g")
    @test relerr_g < 1e-8
end

@testitem "MultiKet DK exact-Hessian handles the near-degenerate confluent limit [#203 AC1]" begin
    using DirectTrajOpt
    using Piccolissimo
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using Random
    using TrajectoryIndexingUtils: slice

    Random.seed!(203_002)
    T = 2.0
    N = 6

    # 4-level drift with a near-degenerate eigenpair; controls forced to zero at
    # every knot so the spectrum stays (near-)degenerate ⇒ the three-index second
    # divided difference hits its confluent branch. The ForwardDiff-∘-expv witness
    # uses no eigenbasis, so it stays stable and is a valid oracle here.
    n = 4
    Q = Matrix(qr(randn(ComplexF64, n, n)).Q)
    λ0 = [0.3, 0.3 + 1e-9, 1.1, 2.0]
    H0m = Q * Diagonal(λ0) * Q'
    H0 = Matrix(Hermitian((H0m + H0m') / 2))
    B1 = randn(ComplexF64, n, n)
    D1 = Matrix(Hermitian(B1 + B1'))
    B2 = randn(ComplexF64, n, n)
    D2 = Matrix(Hermitian(B2 + B2'))
    sys = QuantumSystem(H0, [D1, D2], [1.0, 1.0])

    kin = [normalize(randn(ComplexF64, n)), normalize(randn(ComplexF64, n))]
    kgo = [normalize(randn(ComplexF64, n)), normalize(randn(ComplexF64, n))]
    qtraj = MultiKetTrajectory(sys, kin, kgo, T)

    ℰ_dk = HermitianExponentialIntegrator(
        qtraj,
        N;
        use_analytical = true,
        gauss_newton = false,
    )
    ℰ_fd = HermitianExponentialIntegrator(
        qtraj,
        N;
        use_analytical = false,
        gauss_newton = false,
    )
    @test ℰ_dk.H_dirs !== nothing

    traj = NamedTrajectory(qtraj, N)
    traj.datavec .= randn(length(traj.datavec))
    uname = ℰ_dk.u_name
    for k = 1:N
        traj.datavec[slice(k, traj.components[uname], traj.dim)] .= 0.0
    end
    μ = randn(ℰ_dk.dim)

    Hdk = eval_hessian_of_lagrangian(ℰ_dk, traj, μ)
    Hfd = eval_hessian_of_lagrangian(ℰ_fd, traj, μ)
    @test all(isfinite, Hdk.nzval)                 # confluent limit: no NaN/Inf
    relerr = norm(Hdk - Hfd) / norm(Hfd)
    println(
        "  [#203 AC1] near-degenerate exact-Hessian DK-vs-ForwardDiff rel err = $relerr",
    )
    @test relerr < 1e-7
end

@testitem "MultiKet DK exact-Hessian assembled block is symmetric [#203 AC2]" begin
    using DirectTrajOpt
    using Piccolissimo
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using Random
    using TrajectoryIndexingUtils: slice

    Random.seed!(203_003)
    T = 2.0
    N = 6

    # affine system with globals so every p-p block (u-u, u-g, g-g, Δt-* ) is present
    H = (u, t) -> (u[3] + u[4]) * GATES.Z + u[1] * GATES.X + u[2] * GATES.Y
    sys = QuantumSystem(
        H,
        [1.0, 1.0];
        time_dependent = true,
        global_params = (δ = 0.2, ω = 1.0),
    )
    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    qtraj = MultiKetTrajectory(sys, [ψ0, ψ1], [ψ1, ψ0], T)

    ℰ = HermitianExponentialIntegrator(
        qtraj,
        N;
        use_analytical = true,
        gauss_newton = false,
    )
    traj = NamedTrajectory(qtraj, N)
    traj.datavec .= randn(length(traj.datavec))
    traj.global_data .= randn(length(traj.global_data))
    globals = Piccolissimo.extract_globals(ℰ, traj)
    μ = randn(ℰ.dim)

    # The per-knot CANONICAL block (before the eval-time triu extraction) is what the
    # canonical→trajectory index mapping expects to be symmetric to machine precision.
    for k = 1:(N-1)
        μₖ = μ[slice(k, ℰ.x_dim)]
        Piccolissimo.hessian_of_lagrangian!(
            ℰ.μ∂²ℰs[k],
            ℰ,
            μₖ,
            traj[k],
            traj[k+1],
            k,
            globals,
        )
        B = ℰ.μ∂²ℰs[k]
        @test norm(B - B') ≤ 1e-10 * max(norm(B), eps())
    end
end

@testitem "MultiKet DK exact-Hessian is thread-safe: parallel == serial [#203 AC3]" begin
    using DirectTrajOpt
    using Piccolissimo
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using Random
    using TrajectoryIndexingUtils: slice

    Random.seed!(203_004)
    T = 2.0
    N = 12

    H = (u, t) -> (u[3] + u[4]) * GATES.Z + u[1] * GATES.X + u[2] * GATES.Y
    sys = QuantumSystem(
        H,
        [1.0, 1.0];
        time_dependent = true,
        global_params = (δ = 0.2, ω = 1.0),
    )
    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    qtraj = MultiKetTrajectory(sys, [ψ0, ψ1], [ψ1, ψ0], T)

    ℰ = HermitianExponentialIntegrator(
        qtraj,
        N;
        use_analytical = true,
        gauss_newton = false,
    )
    @test ℰ.H_dirs !== nothing

    traj = NamedTrajectory(qtraj, N)
    traj.datavec .= randn(length(traj.datavec))
    traj.global_data .= randn(length(traj.global_data))
    globals = Piccolissimo.extract_globals(ℰ, traj)
    μ = randn(ℰ.dim)

    println("  [#203 AC3] running on $(Threads.nthreads()) thread(s)")

    # A race on the new per-thread second-order DK scratch would corrupt some knots
    # under Threads.nthreads() > 1; assert threaded == serial on the exact-Hessian.
    Threads.@threads for k = 1:(N-1)
        μₖ = μ[slice(k, ℰ.x_dim)]
        Piccolissimo.hessian_of_lagrangian!(
            ℰ.μ∂²ℰs[k],
            ℰ,
            μₖ,
            traj[k],
            traj[k+1],
            k,
            globals,
        )
    end
    H_threaded = [copy(ℰ.μ∂²ℰs[k]) for k = 1:(N-1)]
    for k = 1:(N-1)
        μₖ = μ[slice(k, ℰ.x_dim)]
        Piccolissimo.hessian_of_lagrangian!(
            ℰ.μ∂²ℰs[k],
            ℰ,
            μₖ,
            traj[k],
            traj[k+1],
            k,
            globals,
        )
    end
    H_serial = [copy(ℰ.μ∂²ℰs[k]) for k = 1:(N-1)]
    @test all(H_threaded[k] == H_serial[k] for k = 1:(N-1))

    # Public threaded assembly is deterministic across repeated calls.
    @test eval_hessian_of_lagrangian(ℰ, traj, μ) == eval_hessian_of_lagrangian(ℰ, traj, μ)
end

@testitem "MultiKet DK exact-Hessian matches FiniteDiff oracle (no autodiff in path) [#203 AC4]" begin
    using DirectTrajOpt
    using Piccolissimo
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using Random

    include("../../../test/test_utils.jl")

    Random.seed!(203_005)
    T = 2.0
    N = 6
    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    qtraj = MultiKetTrajectory(sys, [ψ0, ψ1], [ψ1, ψ0], T)
    traj = NamedTrajectory(qtraj, N)

    # The analytic exact-Hessian (flag on, gauss_newton=false) must agree with an
    # INDEPENDENT FiniteDiff oracle over the FULL Hessian (all p-p blocks) — this is
    # the "no ForwardDiff/FD in the derivative pipeline" invariant: the analytic path
    # is validated by a witness that is not ForwardDiff.
    ℰ_dk = HermitianExponentialIntegrator(
        qtraj,
        N;
        use_analytical = true,
        gauss_newton = false,
    )
    test_integrator(ℰ_dk, traj; gauss_newton = false, atol = 1e-4, show_hessian_diff = true)
end

@testitem "MultiKet DK exact-Hessian solve reaches the ForwardDiff-path solution [#203 AC4]" begin
    using DirectTrajOpt
    using Piccolissimo
    using Piccolo
    using LinearAlgebra
    using Random

    Random.seed!(203_006)
    T = 8.0
    N = 12
    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    initials = [ψ0, ψ1]
    goals = [ψ1, ψ0]

    function solve_path(use_analytical)
        qtraj = MultiKetTrajectory(sys, initials, goals, T)
        ℰ = HermitianExponentialIntegrator(
            qtraj,
            N;
            use_analytical = use_analytical,
            gauss_newton = false,
        )
        qcp = SmoothPulseProblem(qtraj, N; Q = 100.0, R = 1e-2, integrator = ℰ)
        solve!(qcp; max_iter = 150, print_level = 0, verbose = false)
        traj = get_trajectory(qcp)
        fids = [
            fidelity(iso_to_ket(traj[end][name]), goal) for
            (name, goal) in zip(state_names(qtraj), goals)
        ]
        return sum(fids) / length(fids)
    end

    fid_dk = solve_path(true)     # analytic exact-Hessian path
    fid_fd = solve_path(false)    # retained ForwardDiff exact-Hessian path
    println("  [#203 AC4] avg fidelity: analytic=$fid_dk  ForwardDiff=$fid_fd")
    @test fid_dk > 0.99
    @test isapprox(fid_dk, fid_fd; atol = 1e-2)
end
