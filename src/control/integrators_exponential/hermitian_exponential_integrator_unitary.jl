# ============================================================================ #
# UnitaryTrajectory Methods
# ============================================================================ #

"""
    HermitianExponentialIntegrator(qtraj::UnitaryTrajectory, N::Int; kwargs...)

Construct a Hermitian exponential integrator for unitary trajectory evolution.

# Arguments
- `qtraj::UnitaryTrajectory`: The quantum trajectory with system and pulse
- `N::Int`: Number of knot points for discretization

# Keyword Arguments
- `global_names::Union{Vector{Symbol}, Nothing}=nothing`: Names of global (time-invariant) 
  variables to include in dynamics. If `nothing`, auto-detects from `sys.global_params`.
"""
function HermitianExponentialIntegrator(
    qtraj::UnitaryTrajectory,
    N::Int;
    global_names::Union{Vector{Symbol},Nothing} = nothing,
    gauss_newton::Bool = false,
    matrix_free::Bool = false,
    use_analytical::Bool = true,
)
    sys = get_system(qtraj)
    x = state_name(qtraj)
    u = drive_name(qtraj)

    # Auto-detect globals if not specified
    if isnothing(global_names)
        if !isempty(sys.global_params)
            global_names = collect(keys(sys.global_params))
        else
            global_names = Symbol[]
        end
    end

    # Build trajectory WITH globals so structure functions get correct dimensions
    if !isempty(global_names)
        # Initialize globals from system params if available, else zeros
        global_data = Dict{Symbol,Vector{Float64}}()
        for name in global_names
            if hasproperty(sys, :global_params) && haskey(sys.global_params, name)
                global_data[name] = [sys.global_params[name]]
            else
                global_data[name] = [0.0]
            end
        end
        traj = NamedTrajectory(qtraj, N; global_data = global_data)
    else
        traj = NamedTrajectory(qtraj, N)
    end

    return _hermitian_exp_unitary(
        sys,
        x,
        u,
        traj,
        global_names;
        gauss_newton = gauss_newton,
        matrix_free = matrix_free,
        use_analytical = use_analytical,
    )
end

# Per-member inner constructor (#408's sampling lane): all pieces
# pre-materialized. The (qtraj, N) outer method is this plus trajectory assembly.
function _hermitian_exp_unitary(
    sys,
    x::Symbol,
    u::Symbol,
    traj::NamedTrajectory,
    global_names::Vector{Symbol};
    gauss_newton::Bool = false,
    matrix_free::Bool = false,
    use_analytical::Bool = true,
)
    @assert traj.N > 1 "Trajectory must have at least two timesteps."

    global_dim = traj.global_dim

    # Get ketdim from system - need to sample with extended control if globals present
    sample_controls = traj[1][u]
    if global_dim > 0 && !isempty(sys.global_params)
        # Build sample extended control from sys.global_params
        sample_globals = [sys.global_params[name] for name in global_names]
        sample_u = vcat(sample_controls, sample_globals)
    else
        sample_u = sample_controls
    end
    G_sample = sys.𝒢(sample_u, 0.0)
    ketdim = size(G_sample, 1) ÷ 2

    # Dimensions for API (u_dim is control dimension, not including globals)
    x_dim = traj.dims[x]
    u_dim = traj.dims[u]
    N = traj.N
    dim = x_dim * (N - 1)

    # Variables: [xₖ, uₖ, Δtₖ, xₖ₊₁]
    var_dim = 2*x_dim + u_dim + 1

    # Use structure functions to create templates with canonical ordering
    ∂ℰ_template = jacobian_structure(UnitaryTrajectory, x, u, ketdim, traj)
    μ∂²ℰ_template = hessian_structure(x_dim, u_dim, global_dim)

    # Preallocate one copy per knot point
    ∂ℰs = [copy(∂ℰ_template) for _ = 1:(N-1)]
    μ∂²ℰs = [copy(μ∂²ℰ_template) for _ = 1:(N-1)]

    # Build var_comps
    x_comps = traj.components[x]
    u_comps = traj.components[u]
    Δt_comp = traj.components[traj.timestep][1]
    var_comps = [collect(x_comps); collect(u_comps); Δt_comp]

    z_dim = traj.dim

    # Preallocate identity for Kronecker products
    Id = Matrix{Float64}(I(ketdim))

    # Preallocate per-thread exp_eigen! buffers (7 groups, sized to ketdim).
    # Sized to maxthreadid() so Threads.@threads loops in eval_jacobian /
    # eval_hessian_of_lagrangian can safely index by threadid() without racing.
    nthr = Threads.maxthreadid()
    bufs = _alloc_exp_eigen_bufs(ketdim, nthr)

    # Analytic Daleckii–Krein setup (#202): affine drive directions + per-thread scratch.
    # Populated for all variants at construction; the unitary dense wiring lands in slice ④.
    H_dirs = _build_affine_directions(u_ -> sys.H(u_, 0.0), u_dim + global_dim, ketdim)
    dk_bufs = _alloc_dk_bufs(ketdim, nthr)
    dk_so_bufs = _alloc_dk_so_bufs(ketdim, nthr)

    return HermitianExponentialIntegrator{UnitaryTrajectory}(
        u_ -> sys.𝒢(u_, 0.0),
        u_ -> sys.H(u_, 0.0),
        [x],  # Wrap in vector for unified API
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
        Id,       # Preallocated identity
        nothing,  # ∂u∂Δt_bufs (not used for unitary)
        nothing,  # ∂²u_bufs (not used for unitary)
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
        matrix_free,  # inert for unitary: generic eval_jacobian ignores it → always sparse
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

@views function (ℰ::HermitianExponentialIntegrator{UnitaryTrajectory})(
    δₖ::AbstractVector,
    zₖ::KnotPoint,
    zₖ₊₁::KnotPoint,
    k::Int,
    globals::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    Ũ⃗ₖ = zₖ[x_name(ℰ)]
    Ũ⃗ₖ₊₁ = zₖ₊₁[x_name(ℰ)]
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

    copyto!(H_buf, ℰ.H(uₖ))
    exp_eigen!(expG_buf, H_buf, V_buf, λ_buf, cis_diag_buf, tmp_buf, work_buf, Δtₖ)
    expGₖ = expG_buf

    # Exploit Kronecker structure: (I ⊗ A)·vec(X) = vec(A·X) — the propagator acts on
    # the matrix form of the vectorized unitary. A `mul!` on the reshaped trajectory
    # views cannot be used here: on some builds (Julia 1.10 x64, measured in CI) the
    # ReshapedArray operands fall through to LinearAlgebra's allocating generic path
    # (~21 KB per call), which blew the forward-pass allocation budget. The product is
    # explicitly computed instead: O((2n)²·W) — lower order than the O(n³)
    # eigendecomposition already performed above — and allocation-free by construction
    # on every platform and version.
    ketdim = ℰ.ketdim
    width = length(Ũ⃗ₖ) ÷ (2 * ketdim)
    Ũₖ = reshape(Ũ⃗ₖ, 2 * ketdim, width)
    δₖ_mat = reshape(δₖ, 2 * ketdim, width)
    @inbounds for j = 1:width, i = 1:(2*ketdim)
        s = 0.0
        @simd for l = 1:(2*ketdim)
            s += expGₖ[i, l] * Ũₖ[l, j]
        end
        δₖ_mat[i, j] = -s
    end
    @inbounds δₖ .+= Ũ⃗ₖ₊₁

    return nothing
end

@views function jacobian!(
    ∂ℰ::SparseMatrixCSC,
    ℰ::HermitianExponentialIntegrator{UnitaryTrajectory},
    zₖ::KnotPoint,
    zₖ₊₁::KnotPoint,
    k::Int,
    globals::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    Ũ⃗ₖ = zₖ[x_name(ℰ)]
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
    ketdim = ℰ.ketdim
    x_dim = ℰ.x_dim
    u_dim = length(aₖ)  # Control dimension (not including globals)

    # Reshape for Kronecker-free operations
    Ũₖ = reshape(Ũ⃗ₖ, 2*ketdim, ketdim)

    # Use matrix structure dimension, not stored z_dim
    z_dim_matrix = size(∂ℰ, 2) ÷ 2

    x_comps = ℰ.var_comps[1:x_dim]
    u_comps = ℰ.var_comps[(x_dim+1):(x_dim+u_dim)]
    Δt_comp = ℰ.var_comps[x_dim+u_dim+1]

    # Analytic Daleckii–Krein path for affine drives (#204): build M once per knot from the
    # eigenbasis exp_eigen! just produced, then apply per unitary column (block-per-column,
    # since the unitary state is a matrix). Falls back to ForwardDiff for nonlinear drives
    # (H_dirs===nothing).
    use_dk = ℰ.use_analytical && !isnothing(ℰ.H_dirs)
    if use_dk
        dk_ws = ℰ.dk_ws_bufs[tid]
        dk_dΦ = ℰ.dk_dΦ_bufs[tid]
        dk_vin = ℰ.dk_vin_bufs[tid]
        dk_vout = ℰ.dk_vout_bufs[tid]
        dk_divided_difference!(dk_ws.M, λ_buf, Δtₖ)
        u_dirs = view(ℰ.H_dirs, 1:u_dim)
        g_dirs = view(ℰ.H_dirs, (u_dim+1):(u_dim+ℰ.global_dim))
    end

    # ∂Ũ⃗ₖℰ: -I⊗exp(ΔtG) - still need Kronecker for this block (it's the full matrix)
    @inbounds ∂ℰ[:, x_comps] .= -(ℰ.Id ⊗ expGₖ)

    # ∂aₖℰ - control derivatives — analytic DK (per column) or ForwardDiff witness
    if use_dk
        _dk_fill_iso_jac_block_unitary!(
            ∂ℰ[:, u_comps],
            V_buf,
            dk_ws.M,
            u_dirs,
            Ũ⃗ₖ,
            dk_dΦ,
            dk_ws.tmp1,
            dk_ws.tmp2,
            dk_vin,
            dk_vout,
            ketdim,
        )
    else
        # ForwardDiff needs expv (Krylov), not exp_eigen (eigendecomp doesn't support Duals)
        # Note: derivative is only w.r.t. controls (aₖ), globals are constant
        ForwardDiff.jacobian!(
            ∂ℰ[:, u_comps],
            a -> vec(-expv(Δtₖ, ℰ.G(build_extended_control(a, globals)), Ũₖ)),
            aₖ,
        )
    end

    # ∂Δtₖℰ: -(I ⊗ G*exp(ΔtG)) * Ũ⃗ₖ = vec(-G*exp(ΔtG) * Ũₖ)
    GₖexpGₖŨₖ = Gₖ * expGₖ * Ũₖ  # Matrix form
    @inbounds ∂ℰ[:, Δt_comp] .= -vec(GₖexpGₖŨₖ)

    # ∂gℰ: Global variable derivatives
    # Use vcat directly (not build_extended_control) because we differentiate w.r.t. g
    if !isnothing(globals) && ℰ.global_dim > 0
        g_cols = (2*ℰ.z_dim+1):(2*ℰ.z_dim+ℰ.global_dim)
        if use_dk
            _dk_fill_iso_jac_block_unitary!(
                ∂ℰ[:, g_cols],
                V_buf,
                dk_ws.M,
                g_dirs,
                Ũ⃗ₖ,
                dk_dΦ,
                dk_ws.tmp1,
                dk_ws.tmp2,
                dk_vin,
                dk_vout,
                ketdim,
            )
        else
            ForwardDiff.jacobian!(
                ∂ℰ[:, g_cols],
                g -> vec(-expv(Δtₖ, ℰ.G(vcat(aₖ, g)), Ũₖ)),
                globals,
            )
        end
    end

    # ∂Ũ⃗ₖ₊₁ℰ: Identity (already in sparsity pattern)

    return nothing
end

function jacobian_structure(
    ::Type{UnitaryTrajectory},
    x_name::Symbol,
    u_name::Symbol,
    ketdim::Int,
    traj::NamedTrajectory,
)
    x_dim = traj.dims[x_name]
    u_dim = traj.dims[u_name]
    z_dim = traj.dim
    global_dim = traj.global_dim
    x_comps = traj.components[x_name]
    u_comps = traj.components[u_name]
    Δt_comp = traj.components[traj.timestep][1]

    ∂ℰ = spzeros(x_dim, 2z_dim + global_dim)

    # ∂Ũ⃗ₖ₊₁ℰ: Identity
    ∂ℰ[:, z_dim .+ x_comps] = I(x_dim)

    # ∂Ũ⃗ₖℰ: Block structure for unitary
    ∂ℰ[:, x_comps] = I(ketdim) ⊗ ones(2ketdim, 2ketdim)

    # ∂aₖℰ
    ∂ℰ[:, u_comps] .= 1.0

    # ∂Δtₖℰ
    ∂ℰ[:, Δt_comp] .= 1.0

    # ∂gℰ: Global variable derivatives
    if global_dim > 0
        g_cols = (2z_dim+1):(2z_dim+global_dim)
        ∂ℰ[:, g_cols] .= 1.0
    end

    return ∂ℰ
end

@views function hessian_of_lagrangian!(
    μ∂²ℰ::SparseMatrixCSC,
    ℰ::HermitianExponentialIntegrator{UnitaryTrajectory},
    μₖ::AbstractVector,
    zₖ::KnotPoint,
    zₖ₊₁::KnotPoint,
    k::Int,
    globals::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    x_dim = ℰ.x_dim
    u_dim = ℰ.u_dim
    global_dim = ℰ.global_dim
    has_globals = !isnothing(globals) && global_dim > 0

    Ũ⃗ₖ = zₖ[x_name(ℰ)]
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

    # Reshape for Kronecker-free operations
    Ũₖ = reshape(Ũ⃗ₖ, 2*ketdim, ketdim)
    μₖ_mat = reshape(μₖ, 2*ketdim, ketdim)

    # Canonical component indices for knot k
    knot_dim = x_dim + u_dim + 1
    x_can = 1:x_dim
    u_can = (x_dim+1):(x_dim+u_dim)
    Δt_can = x_dim + u_dim + 1
    g_can = has_globals ? ((2*knot_dim+1):(2*knot_dim+global_dim)) : (1:0)

    # Zero out
    @inbounds μ∂²ℰ.nzval .= 0.0

    # Analytic Daleckii–Krein path for affine drives (#204). The unitary state is a matrix,
    # so the Gauss-Newton cross-terms (x-u, x-g) are filled block-per-column via
    # `_dk_fill_iso_jac_block_unitary!`, and the `!gauss_newton` parameter-parameter blocks
    # use the second-order kernel with the column-summed multiplier M_μ = Σ_c χ_c U_c†.
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

    # === Cross-terms (always computed) ===

    # μₖ∂Ũ⃗ₖ∂aₖℰ - analytic DK (per column) or ForwardDiff (Krylov) witness
    if use_dk
        _dk_fill_iso_jac_block_unitary!(
            μ∂²ℰ[x_can, u_can],
            V_buf,
            dk_ws.M,
            u_dirs,
            μₖ,
            dk_dΦ,
            dk_ws.tmp1,
            dk_ws.tmp2,
            dk_vin,
            dk_vout,
            ketdim;
            adjoint_op = true,
        )
    else
        ForwardDiff.jacobian!(
            μ∂²ℰ[x_can, u_can],
            a -> vec(-expv(Δtₖ, ℰ.G(build_extended_control(a, globals))', μₖ_mat)),
            aₖ,
        )
    end

    # μₖ∂Ũ⃗ₖ∂Δtₖℰ
    @inbounds μ∂²ℰ[x_can, Δt_can] .= -vec(GₖexpGₖ' * μₖ_mat)

    if has_globals
        # μₖ∂Ũ⃗ₖ∂gℰ - x-g block (cross-term, always kept)
        if use_dk
            _dk_fill_iso_jac_block_unitary!(
                μ∂²ℰ[x_can, g_can],
                V_buf,
                dk_ws.M,
                g_dirs,
                μₖ,
                dk_dΦ,
                dk_ws.tmp1,
                dk_ws.tmp2,
                dk_vin,
                dk_vout,
                ketdim;
                adjoint_op = true,
            )
        else
            ForwardDiff.jacobian!(
                μ∂²ℰ[x_can, g_can],
                g -> vec(-expv(Δtₖ, ℰ.G(vcat(aₖ, g))', μₖ_mat)),
                globals,
            )
        end
    end

    # === Parameter-parameter blocks (skipped in GN) ===
    if !ℰ.gauss_newton
        if use_dk
            # ── Analytic exact-Hessian parameter-parameter blocks (#204) ──
            # For the unitary, the objective is s = -Re tr(M_μ† Φ) with the column-summed
            # multiplier M_μ = Σ_c χ_c U_c† = X U† (χ_c/U_c = complex columns of μₖ / Ũ⃗ₖ).
            # Every p-p block is a linear functional of D²Φ against M_μ, matching the
            # ForwardDiff `-dot(μₖ_mat, expv(...))` witness.
            n_ag = u_dim + global_dim
            dk_so_ws = ℰ.dk_so_ws_bufs[tid]
            dk_μmat = ℰ.dk_μmat_bufs[tid]
            dk_B = zeros(n_ag, n_ag)
            Xc = zeros(ComplexF64, ketdim, ketdim)     # multiplier columns χ_c
            Uc = zeros(ComplexF64, ketdim, ketdim)     # unitary columns U_c
            HXc = zeros(ComplexF64, ketdim, ketdim)    # H X
            Yc = zeros(ComplexF64, ketdim, ketdim)     # Φ U
            tmpc = zeros(ComplexF64, ketdim, ketdim)
            for c = 1:ketdim
                off = (c - 1) * 2 * ketdim
                _iso_vec_to_complex!(
                    view(Xc, :, c),
                    view(μₖ, (off+1):(off+2*ketdim)),
                    ketdim,
                )
                _iso_vec_to_complex!(
                    view(Uc, :, c),
                    view(Ũ⃗ₖ, (off+1):(off+2*ketdim)),
                    ketdim,
                )
            end
            mul!(dk_μmat, Xc, Uc')        # M_μ = X U† = Σ_c χ_c U_c†
            mul!(HXc, H_buf, Xc)          # H X
            mul!(tmpc, V_buf', Uc)        # V† U
            tmpc .= cis_diag_buf .* tmpc  # diag(cis(-Δt λ)) V† U
            mul!(Yc, V_buf, tmpc)         # Φ U (reuse exp_eigen! eigenbasis)

            # u-u / u-g / g-g via the three-index second divided difference:
            #   μ∂²ℰ[p,q] += -Re tr(M_μ† D²Φ[Hₚ,H_q]) = -dk_B[p,q]
            dk_second_order_block!(dk_B, dk_so_ws, V_buf, λ_buf, ℰ.H_dirs, dk_μmat, Δtₖ)
            μ∂²ℰ[u_can, u_can] .-= dk_B[1:u_dim, 1:u_dim]
            if has_globals
                μ∂²ℰ[u_can, g_can] .-= dk_B[1:u_dim, (u_dim+1):n_ag]
                μ∂²ℰ[g_can, g_can] .-= dk_B[(u_dim+1):n_ag, (u_dim+1):n_ag]
            end

            # Parameter-Δt cross terms (u-Δt, Δt-g), summed over unitary columns. Closed
            # form (mirrors the multiket per-ket formula, aggregated as traces):
            #   val_p = -Re tr(M_μ† ∂²Φ/∂aₚ∂Δt),  ∂²Φ/∂aₚ∂Δt = -iHₚΦ - iH ∂ₚΦ
            #         = Re(i·[tr(X† Hₚ Y) + tr((HX)† (∂ₚΦ) U)]).
            for p = 1:n_ag
                Hₚ = ℰ.H_dirs[p]
                dk_apply!(dk_dΦ, V_buf, dk_ws.M, Hₚ, dk_ws.tmp1, dk_ws.tmp2)  # ∂ₚΦ
                mul!(tmpc, Hₚ, Yc)
                s = dot(Xc, tmpc)             # tr(X† Hₚ Y)
                mul!(tmpc, dk_dΦ, Uc)
                s += dot(HXc, tmpc)           # tr((HX)† (∂ₚΦ) U)
                val = real(im * s)
                if p <= u_dim
                    μ∂²ℰ[u_can[p], Δt_can] += val
                else
                    μ∂²ℰ[Δt_can, g_can[p-u_dim]] += val
                end
            end
        else
            # μₖ∂aₖ∂Δtₖℰ - use expv on matrix
            ForwardDiff.gradient!(
                μ∂²ℰ[u_can, Δt_can],
                a ->
                    -dot(
                        μₖ_mat,
                        ℰ.G(build_extended_control(a, globals)) *
                        expv(Δtₖ, ℰ.G(build_extended_control(a, globals)), Ũₖ),
                    ),
                aₖ,
            )

            if has_globals
                # Full [a,g] Hessian → u-u, u-g, g-g blocks
                full_ag_dim = u_dim + global_dim
                ∂²ag_buf = zeros(full_ag_dim, full_ag_dim)
                ForwardDiff.hessian!(
                    ∂²ag_buf,
                    ag -> -dot(μₖ_mat, expv(Δtₖ, ℰ.G(ag), Ũₖ)),
                    vcat(aₖ, globals),
                )
                μ∂²ℰ[u_can, u_can] .= ∂²ag_buf[1:u_dim, 1:u_dim]
                μ∂²ℰ[u_can, g_can] .= ∂²ag_buf[1:u_dim, (u_dim+1):end]
                μ∂²ℰ[g_can, g_can] .= ∂²ag_buf[(u_dim+1):end, (u_dim+1):end]

                # μₖ∂Δtₖ∂gℰ - Δt-g block
                ∂Δt∂g_buf = zeros(global_dim)
                ForwardDiff.gradient!(
                    ∂Δt∂g_buf,
                    g -> begin
                        # Use vcat directly (not build_extended_control) because we differentiate w.r.t. g
                        Gₖ_g = ℰ.G(vcat(aₖ, g))
                        -dot(μₖ_mat, Gₖ_g * expv(Δtₖ, Gₖ_g, Ũₖ))
                    end,
                    globals,
                )
                for j = 1:global_dim
                    μ∂²ℰ[Δt_can, g_can[j]] = ∂Δt∂g_buf[j]
                end
            else
                # μₖ∂²aₖℰ (no globals case, p,p block)
                ForwardDiff.hessian!(
                    μ∂²ℰ[u_can, u_can],
                    a ->
                        -dot(
                            μₖ_mat,
                            expv(Δtₖ, ℰ.G(build_extended_control(a, globals)), Ũₖ),
                        ),
                    aₖ,
                )
            end
        end

        # μₖ∂²Δtₖℰ (closed form, analytic in both paths)
        @inbounds μ∂²ℰ[Δt_can, Δt_can] = -dot(μₖ_mat, Gₖ * GₖexpGₖ * Ũₖ)
    end

    # Symmetrize
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
            # x-g is a cross-term, always symmetrized
            for i = 1:x_dim
                μ∂²ℰ[g_can[j], x_can[i]] = μ∂²ℰ[x_can[i], g_can[j]]
            end
            # u-g and Δt-g are (p,p) blocks — skip in GN
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

function get_jacobian_structure(
    ℰ::HermitianExponentialIntegrator{UnitaryTrajectory},
    traj::NamedTrajectory,
)
    # Derived-Δt dynamics (Piccolo.jl#321): packed two-knot windows + the warp
    # column under a warp; historical structure otherwise.
    traj.warp !== nothing && return _get_jacobian_structure_warped(ℰ, traj)

    N = traj.N
    x_dim = traj.dims[x_name(ℰ)]
    z_dim = traj.dim
    F_dim = x_dim * (N - 1)
    global_dim = traj.global_dim
    Z_dim = z_dim * N + global_dim
    ∂F = spzeros(F_dim, Z_dim)

    # Get structure for a single knot point
    ∂ℰ_k = jacobian_structure(UnitaryTrajectory, x_name(ℰ), ℰ.u_name, ℰ.ketdim, traj)

    # Place per-knot structure for all knot points
    for k = 1:(N-1)
        ∂F[slice(k, x_dim), slice(k, 1:2z_dim, z_dim)] = ∂ℰ_k[:, 1:2z_dim]
    end

    # Global columns: all knot points contribute to the same global columns. Map the
    # per-knot block (global_names order) to the trajectory's global_components
    # order via `_global_full_cols` so STRUCTURE matches VALUES (permuted-∂c/∂θ fix).
    if global_dim > 0
        g_cols_local = (2z_dim+1):(2z_dim+global_dim)
        g_cols_full = _global_full_cols(ℰ, traj)
        for k = 1:(N-1)
            ∂F[slice(k, x_dim), g_cols_full] = ∂ℰ_k[:, g_cols_local]
        end
    end

    return ∂F
end

# ============================================================================ #
# API methods for UnitaryTrajectory (explicit dispatch)
# ============================================================================ #

@views function evaluate!(
    δ::AbstractVector,
    ℰ::HermitianExponentialIntegrator{UnitaryTrajectory},
    traj::NamedTrajectory,
)
    # Under a time warp the timestep rows are derived: re-derive them from the
    # warp so Δtₖ = deltats(warp)[k] is what the defect consumes (no-op warp-free).
    traj.warp !== nothing && sync_timesteps!(traj)

    # Extract globals once (constant across knot points)
    globals = extract_globals(ℰ, traj)

    Threads.@threads for k = 1:(traj.N-1)
        δₖ = δ[slice(k, ℰ.x_dim)]
        ℰ(δₖ, traj[k], traj[k+1], k, globals)
    end
    return nothing
end

@views function eval_jacobian(
    ℰ::HermitianExponentialIntegrator{UnitaryTrajectory},
    traj::NamedTrajectory,
)
    # Derived-Δt dynamics (Piccolo.jl#321): with a warp the packed Jacobian
    # drops the derived timestep column and gains the EXACT warp-parameter
    # column (∂cₖ/∂θⱼ = wₖⱼ · ∂cₖ/∂Δtₖ). No warp → historical body, bit-unchanged.
    traj.warp !== nothing && return _eval_jacobian_warped(ℰ, traj)

    N = traj.N
    x_dim = ℰ.x_dim
    z_dim = traj.dim
    F_dim = x_dim * (N - 1)
    Z_dim = z_dim * N + traj.global_dim

    # Extract globals once (constant across knot points)
    globals = extract_globals(ℰ, traj)

    # Fill preallocated structures in parallel
    Threads.@threads for k = 1:(N-1)
        jacobian!(ℰ.∂ℰs[k], ℰ, traj[k], traj[k+1], k, globals)
    end

    # Build var_comps for single state name
    x_comps_now = collect(traj.components[x_name(ℰ)])
    u_comps_now = traj.components[ℰ.u_name]
    Δt_comp_now = traj.components[traj.timestep][1]
    var_comps_now = [x_comps_now; collect(u_comps_now); Δt_comp_now]

    ∂F = spzeros(F_dim, Z_dim)
    @inbounds for k = 1:(N-1)
        ∂F[slice(k, x_dim), slice(k, var_comps_now, z_dim)] = ℰ.∂ℰs[k][:, var_comps_now]
        ∂F[slice(k, x_dim), slice(k, z_dim .+ var_comps_now, z_dim)] =
            ℰ.∂ℰs[k][:, ℰ.z_dim .+ var_comps_now]
    end

    # Global columns: per-knot block is in global_names order; `_global_full_cols`
    # maps it to the trajectory's global_components order (permuted-∂c/∂θ fix).
    if traj.global_dim > 0
        g_cols_local = (2*ℰ.z_dim+1):(2*ℰ.z_dim+traj.global_dim)
        g_cols_full = _global_full_cols(ℰ, traj)
        @inbounds for k = 1:(N-1)
            ∂F[slice(k, x_dim), g_cols_full] = ℰ.∂ℰs[k][:, g_cols_local]
        end
    end

    return ∂F
end

@views function eval_hessian_of_lagrangian(
    ℰ::HermitianExponentialIntegrator{UnitaryTrajectory},
    traj::NamedTrajectory,
    μ::AbstractVector,
)
    # Derived-Δt dynamics (Piccolo.jl#321): with a warp the canonical Δt
    # row/column maps onto the warp-parameter column through the exact chain
    # rule. No warp → historical body, bit-unchanged.
    traj.warp !== nothing && return _eval_hessian_of_lagrangian_warped(ℰ, traj, μ)

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

    # Build index mapping inline
    knot_dim = x_dim + u_dim + 1

    # Canonical indices
    x_can_k = 1:x_dim
    u_can_k = (x_dim+1):(x_dim+u_dim)
    Δt_can_k = x_dim + u_dim + 1
    x_can_k1 = knot_dim .+ (1:x_dim)

    # Trajectory indices - single state name
    x_traj_k = collect(traj.components[x_name(ℰ)])
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

    # Global blocks: cross-terms with per-knot vars + g-g block
    global_dim = traj.global_dim
    if global_dim > 0
        g_can = (2*knot_dim+1):(2*knot_dim+global_dim)
        g_traj = _global_full_cols(ℰ, traj)   # trajectory-order Z columns (permuted-∂c/∂θ fix)
        @inbounds for k = 1:(N-1)
            μ∂²F[slice(k, traj_comps, z_dim), g_traj] .= ℰ.μ∂²ℰs[k][canonical_comps, g_can]
            μ∂²F[g_traj, slice(k, traj_comps, z_dim)] .= ℰ.μ∂²ℰs[k][g_can, canonical_comps]
            μ∂²F[g_traj, g_traj] .+= ℰ.μ∂²ℰs[k][g_can, g_can]
        end
    end

    # Return upper triangle (symmetric matrix, only upper needed by Ipopt)
    return triu(μ∂²F)
end

@testitem "testing HermitianExponentialIntegrator with UnitaryTrajectory" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra

    include("../../../test/test_utils.jl")

    sys = OpenQuantumSystem(
        kron(GATES.Z, GATES.Z),
        [kron(GATES.X, GATES.X), kron(GATES.Y, GATES.Y)],
        [1.0, 1.0],
    )

    Id = Matrix{ComplexF64}(I, sys.levels, sys.levels)
    T = 1.0
    N = 6
    qtraj = UnitaryTrajectory(sys, Id, T)
    traj = NamedTrajectory(qtraj, N)

    ℰ = HermitianExponentialIntegrator(qtraj, N)

    @test ℰ isa HermitianExponentialIntegrator{UnitaryTrajectory}

    test_integrator(ℰ, traj; atol = 1e-3, show_jacobian_diff = false)
end





@testitem "HermitianExponentialIntegrator{Unitary} global Jacobian columns respect trajectory order (permutation regression)" begin
    # Unitary counterpart of the ket permutation regression. `global_params=(b, a)`
    # (non-alphabetical) ⇒ integrator global_names=[:b,:a] but NamedTrajectory sorts
    # globals alphabetically ⇒ trajectory order [:a,:b]; a naive i→i assembly would
    # swap ∂c/∂a and ∂c/∂b. Each global's assembled column, at its trajectory
    # position, must match a finite difference of the forward residual.
    using DirectTrajOpt
    using Piccolo
    using Piccolo.Control.QuantumIntegrators.ExponentialIntegrators
    using NamedTrajectories
    using LinearAlgebra
    using Random

    Random.seed!(90_211)
    H = (u, t) -> u[3] * GATES.Z + u[4] * GATES.Y + u[1] * GATES.X + u[2] * GATES.Y
    sys = OpenQuantumSystem(
        H,
        [1.0, 1.0];
        time_dependent = true,
        global_params = (b = 0.2, a = 0.1),
    )
    U_goal = GATES.X
    T = 2.0
    N = 8
    qtraj = UnitaryTrajectory(sys, U_goal, T)

    ℰ = HermitianExponentialIntegrator(qtraj, N)
    @test ℰ.global_names == [:b, :a]

    traj = NamedTrajectory(qtraj, N)
    @test traj.global_components[:a][1] < traj.global_components[:b][1]  # [:a,:b] ≠ global_names

    traj.datavec .= 0.3 .* randn(length(traj.datavec))
    J = eval_jacobian(ℰ, traj)

    ndof = traj.dim * traj.N
    for name in (:a, :b)
        col = ndof + traj.global_components[name][1]
        base = traj.global_data[traj.global_components[name]][1]
        h = 1e-6
        function resid(v)
            t2 = deepcopy(traj)
            t2.global_data[t2.global_components[name]] .= v
            d = zeros(ℰ.dim)
            DirectTrajOpt.evaluate!(d, ℰ, t2)
            return d
        end
        fd = (resid(base + h) .- resid(base - h)) ./ (2h)
        rel = norm(Vector(J[:, col]) - fd) / max(norm(fd), eps())
        @test rel < 1e-5
    end
end

@testitem "HermitianExponentialIntegrator Jacobian sparsity for UnitaryTrajectory" begin
    using DirectTrajOpt
    using Piccolo
    using NamedTrajectories
    using SparseArrays
    using LinearAlgebra

    T = 1.0
    N = 10
    sys = OpenQuantumSystem(GATES.Z, [GATES.X], [1.0])
    U_goal = GATES.H
    ketdim = 2
    x_dim = 2 * ketdim * ketdim  # 8 (isomorphic unitary)
    u_dim = 1

    qtraj = UnitaryTrajectory(sys, U_goal, T)
    integrator = HermitianExponentialIntegrator(qtraj, N)

    traj = NamedTrajectory(qtraj, N)
    ∂F = Piccolo.Control.QuantumIntegrators.ExponentialIntegrators.get_jacobian_structure(
        integrator,
        traj,
    )

    # Expected sparsity for UnitaryTrajectory with block-diagonal structure:
    # HermitianExponentialIntegrator uses jacobian_structure() which should exploit
    # the block-diagonal structure from independent column evolution.
    # - Parameters per knot point: u_dim + 1 (controls + Δt)
    n_params = u_dim + 1  # 2

    # For block-diagonal state structure:
    # - Each unitary column evolves independently → ketdim blocks of (2ketdim × 2ketdim)
    # - ∂xₖ and ∂xₖ₊₁: ketdim * (2ketdim)² = 2 * 16 = 32 each
    block_state_nnz_per_k = 2 * ketdim * (2*ketdim)^2  # 64
    param_nnz_per_k = x_dim * n_params  # 8 * 2 = 16
    expected_block_diagonal_nnz = (N - 1) * (block_state_nnz_per_k + param_nnz_per_k)

    # Dense would be if we didn't exploit structure
    dense_state_nnz_per_k = 2 * x_dim * x_dim  # 128
    expected_dense_nnz = (N - 1) * (dense_state_nnz_per_k + param_nnz_per_k)

    actual_nnz = nnz(∂F)
    println("  HermitianExponentialIntegrator UnitaryTrajectory Jacobian nnz: $actual_nnz")
    println("  Block-diagonal expected: $expected_block_diagonal_nnz")
    println("  Dense expected: $expected_dense_nnz")

    # Should be at most the block-diagonal count (may be less due to additional sparsity)
    @test actual_nnz ≤ expected_block_diagonal_nnz

    # Verify sparsity is exploited - should be much less than dense
    @test actual_nnz < expected_dense_nnz
    println("  $(round(actual_nnz / expected_dense_nnz * 100, digits=1))% of dense")
end

@testitem "HermitianExponentialIntegrator{UnitaryTrajectory} forward pass minimal alloc" begin
    using LinearAlgebra
    using NamedTrajectories
    using Piccolo
    using BenchmarkTools

    sys = OpenQuantumSystem(PAULIS.Z, [PAULIS.X, PAULIS.Y], [1.0, 1.0])
    N = 11
    times = collect(range(0, 1.0, length = N))
    pulse = LinearSplinePulse(zeros(2, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, GATES[:H])
    traj = NamedTrajectory(qtraj, N)

    ℰ = HermitianExponentialIntegrator(qtraj, N)
    δ = zeros(ℰ.x_dim)
    z1 = traj[1]
    z2 = traj[2]

    ℰ(δ, z1, z2, 1)  # warmup

    # Regression guard for Task 5 of plan-20260417-063000-nonhermitian-exp-integrator.md.
    # The forward path calls exp_eigen! with preallocated struct buffers (no per-call
    # allocation for the eigendecomposition). Residual allocation comes from ℰ.H(uₖ)
    # building a fresh matrix each call; same limitation as Task 4 on the ket forward
    # path. Observed: ~4176 B post-refactor (~6176 B pre-refactor). Threshold mirrors
    # Task 4 (< 8192) pending a follow-on H-buffer refactor.


    allocs = @ballocated $ℰ($δ, $z1, $z2, 1)
    @test allocs < 8192
end

# ============================================================================ #
# NonlinearDrive Tests
# ============================================================================ #

@testitem "HermitianExponentialIntegrator UnitaryTrajectory with NonlinearDrive" begin
    using DirectTrajOpt
    using Piccolo.Control.QuantumIntegrators.ExponentialIntegrators
    using Piccolo
    using SparseArrays
    using LinearAlgebra

    include("../../../test/test_utils.jl")

    # 1-qubit system: H = σz + u₁·σx + u₁²·σz
    drives = AbstractDrive[
        LinearDrive(sparse(ComplexF64.(PAULIS.X)), 1),
        NonlinearDrive(PAULIS.Z, u -> u[1]^2),
    ]
    drive_bounds = [1.0]

    sys = OpenQuantumSystem(PAULIS.Z, drives, drive_bounds)

    @test has_nonlinear_drives(sys.H_drives)

    T = 1.0
    N = 6
    U_goal = GATES.X

    qtraj = UnitaryTrajectory(sys, U_goal, T)
    traj = NamedTrajectory(qtraj, N)

    ℰ = HermitianExponentialIntegrator(qtraj, N)
    @test ℰ isa HermitianExponentialIntegrator{UnitaryTrajectory}

    test_integrator(ℰ, traj; atol = 1e-3, show_jacobian_diff = false)
end



@testitem "HermitianExponentialIntegrator UnitaryTrajectory with NonlinearDrive and global variables" begin
    using DirectTrajOpt
    using DirectTrajOpt: BoundsConstraint
    using Piccolo.Control.QuantumIntegrators.ExponentialIntegrators
    using Piccolo
    using SparseArrays
    using LinearAlgebra
    using NamedTrajectories

    include("../../../test/test_utils.jl")

    # System with nonlinear drive AND global parameter:
    #   H = σz + u₁·σx + u₁²·σz, with optimizable global δ
    δ_init = 0.01

    drives = AbstractDrive[
        LinearDrive(sparse(ComplexF64.(PAULIS.X)), 1),
        NonlinearDrive(PAULIS.Z, u -> u[1]^2),
    ]
    drive_bounds = [1.0]

    sys = OpenQuantumSystem(PAULIS.Z, drives, drive_bounds; global_params = (δ = δ_init,))

    @test has_nonlinear_drives(sys.H_drives)

    T = 10.0
    N = 10
    U_goal = GATES.X

    qtraj = UnitaryTrajectory(sys, U_goal, T)

    integrator = HermitianExponentialIntegrator(qtraj, N)

    @test integrator.global_names == [:δ]
    @test integrator.global_dim == 1

    traj = NamedTrajectory(qtraj, N)
    test_integrator(integrator, traj; atol = 1e-3, show_jacobian_diff = false)

    # Full optimization with global bounds
    δ_bound = 1.0
    qcp = SmoothPulseProblem(
        qtraj,
        N;
        Q = 100.0,
        R = 1e-2,
        integrator = integrator,
        global_bounds = Dict(:δ => δ_bound),
    )

    solve!(qcp; max_iter = 350, print_level = 0)

    result_traj = get_trajectory(qcp)
    δ_opt = result_traj.global_data[result_traj.global_components[:δ]][1]
    @test isfinite(δ_opt)
    @test δ_opt >= -δ_bound - 1e-5
    @test δ_opt <= δ_bound + 1e-5
    println("  Unitary NonlinearDrive + global: δ_init=$δ_init, δ_opt=$δ_opt")
end

# ============================================================================ #
# Gauss-Newton Hessian Tests
# ============================================================================ #

@testitem "HermitianExponentialIntegrator UnitaryTrajectory Gauss-Newton Hessian" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra
    using SparseArrays

    include("../../../test/test_utils.jl")

    sys = OpenQuantumSystem(
        kron(GATES.Z, GATES.Z),
        [kron(GATES.X, GATES.X), kron(GATES.Y, GATES.Y)],
        [1.0, 1.0],
    )

    Id = Matrix{ComplexF64}(I, sys.levels, sys.levels)
    T = 1.0
    N = 6
    qtraj = UnitaryTrajectory(sys, Id, T)
    traj = NamedTrajectory(qtraj, N)

    # Test 1: GN integrator constructs correctly
    ℰ_gn = HermitianExponentialIntegrator(qtraj, N; gauss_newton = true)
    @test ℰ_gn isa HermitianExponentialIntegrator{UnitaryTrajectory}
    @test ℰ_gn.gauss_newton == true

    # Test 2: GN sparsity structure has fewer nonzeros than exact
    ℰ_exact = HermitianExponentialIntegrator(qtraj, N; gauss_newton = false)
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
# Analytic Daleckii–Krein dense-wiring tests (Piccolissimo.jl#204, slice ⑤).
# Mirror the multiket #202/#203 DK tests for the UnitaryTrajectory variant, whose
# state is a MATRIX ⇒ the Jacobian/Hessian cross-terms are filled block-per-column
# and the p-p blocks use the column-summed multiplier M_μ = Σ_c χ_c U_c†. Guards
# the analytic Jacobian, GN Hessian cross-terms, and `!gauss_newton` p-p blocks vs
# the retained ForwardDiff witness + an independent FiniteDiff oracle, per-thread
# safety, the flag toggle, and the nonlinear-drive fallback.
# ============================================================================ #

@testitem "Unitary DK Jacobian matches ForwardDiff witness — affine drive [#204 AC1]" begin
    using DirectTrajOpt
    using Piccolo.Control.QuantumIntegrators.ExponentialIntegrators
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using Random

    Random.seed!(204_201)
    T = 2.0
    N = 8
    sys = OpenQuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    U_goal = GATES.H
    qtraj = UnitaryTrajectory(sys, U_goal, T)

    ℰ_dk = HermitianExponentialIntegrator(qtraj, N; use_analytical = true)
    ℰ_fd = HermitianExponentialIntegrator(qtraj, N; use_analytical = false)

    @test ℰ_dk.H_dirs !== nothing
    @test length(ℰ_dk.H_dirs) == 2          # two linear drives, no globals
    @test ℰ_dk.use_analytical == true

    traj = NamedTrajectory(qtraj, N)
    traj.datavec .= randn(length(traj.datavec))

    Jdk = eval_jacobian(ℰ_dk, traj)
    Jfd = eval_jacobian(ℰ_fd, traj)
    relerr = norm(Jdk - Jfd) / norm(Jfd)
    println("  [#204 AC1] Unitary Jacobian DK-vs-ForwardDiff rel err = $relerr")
    @test relerr < 1e-9
end

@testitem "Unitary DK Gauss-Newton Hessian cross-terms match ForwardDiff [#204 AC1]" begin
    using DirectTrajOpt
    using Piccolo.Control.QuantumIntegrators.ExponentialIntegrators
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using Random

    Random.seed!(204_202)
    T = 2.0
    N = 7

    # --- no globals: x-u cross-term block ---
    sys = OpenQuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    U_goal = GATES.H
    qtraj = UnitaryTrajectory(sys, U_goal, T)

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
    println("  [#204 AC1] Unitary GN Hessian (x-u) DK-vs-ForwardDiff rel err = $relerr")
    @test relerr < 1e-9

    # --- with globals: adds the x-g cross-term block ---
    H = (u, t) -> (u[3] + u[4]) * GATES.Z + u[1] * GATES.X + u[2] * GATES.Y
    gsys = OpenQuantumSystem(
        H,
        [1.0, 1.0];
        time_dependent = true,
        global_params = (δ = 0.2, ω = 1.0),
    )
    gqtraj = UnitaryTrajectory(gsys, U_goal, T)
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
    @test length(ℰg_dk.H_dirs) == 4
    @test ℰg_dk.global_dim == 2

    gtraj = NamedTrajectory(gqtraj, N)
    gtraj.datavec .= randn(length(gtraj.datavec))
    gtraj.global_data .= randn(length(gtraj.global_data))
    μg = randn(ℰg_dk.dim)

    Hgdk = eval_hessian_of_lagrangian(ℰg_dk, gtraj, μg)
    Hgfd = eval_hessian_of_lagrangian(ℰg_fd, gtraj, μg)
    relerr_g = norm(Hgdk - Hgfd) / norm(Hgfd)
    println(
        "  [#204 AC1] Unitary GN Hessian (x-u & x-g) DK-vs-ForwardDiff rel err = $relerr_g",
    )
    @test relerr_g < 1e-9
end

@testitem "Unitary DK exact-Hessian p-p blocks match ForwardDiff witness — affine [#204 AC1]" begin
    using DirectTrajOpt
    using Piccolo.Control.QuantumIntegrators.ExponentialIntegrators
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using SparseArrays
    using Random

    Random.seed!(204_203)
    T = 2.0
    N = 7

    # --- no globals: u-u, u-Δt, Δt-Δt p-p blocks ---
    sys = OpenQuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    U_goal = GATES.H
    qtraj = UnitaryTrajectory(sys, U_goal, T)

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

    traj = NamedTrajectory(qtraj, N)
    traj.datavec .= randn(length(traj.datavec))
    μ = randn(ℰ_dk.dim)

    Hdk = eval_hessian_of_lagrangian(ℰ_dk, traj, μ)
    Hfd = eval_hessian_of_lagrangian(ℰ_fd, traj, μ)
    relerr = norm(Hdk - Hfd) / norm(Hfd)
    println(
        "  [#204 AC1] Unitary exact-Hessian (no globals) DK-vs-ForwardDiff rel err = $relerr",
    )
    @test relerr < 1e-8

    # Isolate the parameter-parameter blocks: present in exact, absent in GN.
    ℰ_gn = HermitianExponentialIntegrator(qtraj, N; gauss_newton = true)
    H_exact_struct = get_hessian_of_lagrangian_structure(ℰ_dk, traj)
    H_gn_struct = get_hessian_of_lagrangian_structure(ℰ_gn, traj)
    pp_mask = (triu(H_exact_struct) .!= 0) .& .!(triu(H_gn_struct) .!= 0)
    @test any(pp_mask)
    ppdk = triu(Hdk)[pp_mask]
    ppfd = triu(Hfd)[pp_mask]
    @test norm(ppdk - ppfd) / norm(ppfd) < 1e-8
    @test !all(ppdk .== 0.0)

    # --- with globals: adds u-g, g-g, Δt-g p-p blocks ---
    H = (u, t) -> (u[3] + u[4]) * GATES.Z + u[1] * GATES.X + u[2] * GATES.Y
    gsys = OpenQuantumSystem(
        H,
        [1.0, 1.0];
        time_dependent = true,
        global_params = (δ = 0.2, ω = 1.0),
    )
    gqtraj = UnitaryTrajectory(gsys, U_goal, T)
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
    @test length(ℰg_dk.H_dirs) == 4
    @test ℰg_dk.global_dim == 2

    gtraj = NamedTrajectory(gqtraj, N)
    gtraj.datavec .= randn(length(gtraj.datavec))
    gtraj.global_data .= randn(length(gtraj.global_data))
    μg = randn(ℰg_dk.dim)

    Hgdk = eval_hessian_of_lagrangian(ℰg_dk, gtraj, μg)
    Hgfd = eval_hessian_of_lagrangian(ℰg_fd, gtraj, μg)
    relerr_g = norm(Hgdk - Hgfd) / norm(Hgfd)
    println(
        "  [#204 AC1] Unitary exact-Hessian (globals) DK-vs-ForwardDiff rel err = $relerr_g",
    )
    @test relerr_g < 1e-8
end

@testitem "Unitary DK path is thread-safe: parallel == serial [#204 AC4]" begin
    using DirectTrajOpt
    using Piccolo.Control.QuantumIntegrators.ExponentialIntegrators
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using Random
    using TrajectoryIndexingUtils: slice

    Random.seed!(204_204)
    T = 2.0
    N = 12

    H = (u, t) -> (u[3] + u[4]) * GATES.Z + u[1] * GATES.X + u[2] * GATES.Y
    sys = OpenQuantumSystem(
        H,
        [1.0, 1.0];
        time_dependent = true,
        global_params = (δ = 0.2, ω = 1.0),
    )
    qtraj = UnitaryTrajectory(sys, GATES.H, T)

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
    globals =
        Piccolo.Control.QuantumIntegrators.ExponentialIntegrators.extract_globals(ℰ, traj)
    μ = randn(ℰ.dim)

    println("  [#204 AC4] Unitary running on $(Threads.nthreads()) thread(s)")

    Threads.@threads for k = 1:(N-1)
        Piccolo.Control.QuantumIntegrators.ExponentialIntegrators.jacobian!(
            ℰ.∂ℰs[k],
            ℰ,
            traj[k],
            traj[k+1],
            k,
            globals,
        )
    end
    J_threaded = [copy(ℰ.∂ℰs[k]) for k = 1:(N-1)]
    for k = 1:(N-1)
        Piccolo.Control.QuantumIntegrators.ExponentialIntegrators.jacobian!(
            ℰ.∂ℰs[k],
            ℰ,
            traj[k],
            traj[k+1],
            k,
            globals,
        )
    end
    J_serial = [copy(ℰ.∂ℰs[k]) for k = 1:(N-1)]
    @test all(J_threaded[k] == J_serial[k] for k = 1:(N-1))

    Threads.@threads for k = 1:(N-1)
        μₖ = μ[slice(k, ℰ.x_dim)]
        Piccolo.Control.QuantumIntegrators.ExponentialIntegrators.hessian_of_lagrangian!(
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
        Piccolo.Control.QuantumIntegrators.ExponentialIntegrators.hessian_of_lagrangian!(
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

    @test eval_jacobian(ℰ, traj) == eval_jacobian(ℰ, traj)
    @test eval_hessian_of_lagrangian(ℰ, traj, μ) == eval_hessian_of_lagrangian(ℰ, traj, μ)
end

@testitem "Unitary DK path matches FiniteDiff oracle [#204 AC1]" begin
    using DirectTrajOpt
    using Piccolo.Control.QuantumIntegrators.ExponentialIntegrators
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using Random

    include("../../../test/test_utils.jl")

    Random.seed!(204_205)
    T = 2.0
    N = 6
    sys = OpenQuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    qtraj = UnitaryTrajectory(sys, GATES.H, T)
    traj = NamedTrajectory(qtraj, N)

    # Gauss-Newton analytic path vs FiniteDiff oracle
    ℰ_gn =
        HermitianExponentialIntegrator(qtraj, N; use_analytical = true, gauss_newton = true)
    test_integrator(
        ℰ_gn,
        traj;
        gauss_newton = true,
        atol = 1e-4,
        show_jacobian_diff = false,
    )

    # Exact-Hessian analytic path vs FiniteDiff oracle (no autodiff in the path)
    ℰ_ex = HermitianExponentialIntegrator(
        qtraj,
        N;
        use_analytical = true,
        gauss_newton = false,
    )
    test_integrator(
        ℰ_ex,
        traj;
        gauss_newton = false,
        atol = 1e-4,
        show_jacobian_diff = false,
    )
end

@testitem "Unitary nonlinear drive falls back to ForwardDiff [#204 AC5]" begin
    using DirectTrajOpt
    using Piccolo.Control.QuantumIntegrators.ExponentialIntegrators
    using Piccolo
    using SparseArrays
    using NamedTrajectories
    using LinearAlgebra
    using Random

    include("../../../test/test_utils.jl")

    Random.seed!(204_206)
    # H = σz + u₁·σx + u₁²·σz — nonlinear in the control.
    drives = AbstractDrive[
        LinearDrive(sparse(ComplexF64.(PAULIS.X)), 1),
        NonlinearDrive(PAULIS.Z, u -> u[1]^2),
    ]
    sys = OpenQuantumSystem(PAULIS.Z, drives, [1.0])
    @test has_nonlinear_drives(sys.H_drives)

    T = 1.0
    N = 6
    qtraj = UnitaryTrajectory(sys, GATES.X, T)

    ℰ_on = HermitianExponentialIntegrator(qtraj, N; use_analytical = true)
    ℰ_off = HermitianExponentialIntegrator(qtraj, N; use_analytical = false)
    @test ℰ_on.H_dirs === nothing
    @test ℰ_on.use_analytical == true

    traj = NamedTrajectory(qtraj, N)
    traj.datavec .= randn(length(traj.datavec))

    # Flag-on and flag-off produce IDENTICAL results (both take the ForwardDiff path).
    @test eval_jacobian(ℰ_on, traj) == eval_jacobian(ℰ_off, traj)

    traj2 = NamedTrajectory(qtraj, N)
    test_integrator(ℰ_on, traj2; atol = 1e-3, show_jacobian_diff = false)
end

# ============================================================================ #
# Phase 1b (Piccolo.jl#321): the defect warp column under a time warp
# (NamedTrajectories#161) — UnitaryTrajectory cell. Exact chain entries
# ∂cₖ/∂θⱼ = (∂Δtₖ/∂θⱼ)·∂cₖ/∂Δtₖ with ∂cₖ/∂Δtₖ = −G(uₖ)·Ũₖ₊₁ on the satisfied
# defect (vec over the stacked iso columns), plus packed FD parity.
# ============================================================================ #

@testitem "HermitianExponentialIntegrator{Unitary} defect warp column: exact entries + FD parity" begin
    using Piccolo
    using Piccolo.Control.QuantumIntegrators.ExponentialIntegrators:
        _hermitian_exp_unitary, exp_eigen
    using DirectTrajOpt
    using NamedTrajectories
    using TrajectoryIndexingUtils: slice
    using Random
    using LinearAlgebra
    using SparseArrays
    using Test

    sys = OpenQuantumSystem(PAULIS.Z, [PAULIS.X, PAULIS.Y], [1.0, 1.0])
    N = 6
    T0 = 1.7
    Random.seed!(20260831)
    Ũ⃗ = randn(8, N)   # stacked iso columns of a 2×2 "unitary" (values arbitrary for FD)
    u = 0.3 .* randn(2, N)

    mk_traj(Ũ⃗) = NamedTrajectory(
        (Ũ⃗ = Ũ⃗, u = u, Δt = fill(T0 / (N - 1), 1, N));  # Δt ignored — derived
        controls = (:u,),
        timestep = :Δt,
        warp = GlobalScale(T0),
    )
    traj = mk_traj(Ũ⃗)
    ℰ = _hermitian_exp_unitary(sys, :Ũ⃗, :u, traj, Symbol[])

    x_dim = ℰ.x_dim                      # 8 (2 iso columns of 4)
    wₖ = 1 / (N - 1)
    T_col = (traj.dim - traj.dims[:Δt]) * traj.N + traj.global_dim + 1

    # evaluate! reads the DERIVED rows (synced from the warp)
    δ = zeros(x_dim * (N - 1))
    evaluate!(δ, ℰ, traj)
    @test !all(iszero, δ)
    for k = 1:(N-1)
        Φₖ = exp_eigen((T0 / (N - 1)) .* sys.H(u[:, k], 0.0))
        @test δ[((k-1)*x_dim) .+ (1:x_dim)] ≈ Ũ⃗[:, k+1] - vec(Φₖ * reshape(Ũ⃗[:, k], 4, 2))
    end

    # Jacobian: packed size + EXACT warp column
    ∂F = eval_jacobian(ℰ, traj)
    @test size(∂F) == (x_dim * (N - 1), length(vec(traj)))
    for k = 1:(N-1)
        Gₖ = sys.𝒢(u[:, k], 0.0)
        Φₖ = exp_eigen((T0 / (N - 1)) .* sys.H(u[:, k], 0.0))
        dcdΔt = -vec(Gₖ * Φₖ * reshape(Ũ⃗[:, k], 4, 2))   # ∂cₖ/∂Δtₖ
        @test ∂F[((k-1)*x_dim) .+ (1:x_dim), T_col] ≈ wₖ .* dcdΔt
    end

    # satisfied defect: ∂cₖ/∂Δtₖ = −G(uₖ)·Ũₖ₊₁ exactly
    Ũ⃗sat = copy(Ũ⃗)
    for k = 1:(N-1)
        Φₖ = exp_eigen((T0 / (N - 1)) .* sys.H(u[:, k], 0.0))
        Ũ⃗sat[:, k+1] = vec(Φₖ * reshape(Ũ⃗sat[:, k], 4, 2))
    end
    traj_sat = mk_traj(Ũ⃗sat)
    ℰ_sat = _hermitian_exp_unitary(sys, :Ũ⃗, :u, traj_sat, Symbol[])
    ∂F_sat = eval_jacobian(ℰ_sat, traj_sat)
    for k = 1:(N-1)
        @test ∂F_sat[((k-1)*x_dim) .+ (1:x_dim), T_col] ≈
              -wₖ .* vec(sys.𝒢(u[:, k], 0.0) * reshape(Ũ⃗sat[:, k+1], 4, 2))
    end

    # FD parity of the full Jacobian over the PACKED vector (perturbs T too)
    function fd_jac(f, z0; h = 1e-6)
        f0 = f(z0)
        J = zeros(length(f0), length(z0))
        for j in eachindex(z0)
            e = zeros(length(z0))
            e[j] = h
            J[:, j] = (f(z0 .+ e) .- f(z0 .- e)) ./ (2h)
        end
        return J
    end
    f̂ = Z⃗ -> begin
        unpack!(traj, Z⃗)
        δ = zeros(x_dim * (N - 1))
        evaluate!(δ, ℰ, traj)
        return δ
    end
    z0 = collect(vec(traj))
    ∂F_fd = fd_jac(f̂, z0)
    @test all(isapprox.(∂F, ∂F_fd; atol = 1e-5, rtol = 1e-5))

    # FD parity of the Hessian of the Lagrangian (exact mode)
    μ = 0.3 .* collect(1.0:(x_dim*(N-1)))
    H = eval_hessian_of_lagrangian(ℰ, traj, μ)
    @test size(H) == (length(vec(traj)), length(vec(traj)))
    function fd_hess(g, z0; h = 1e-4)
        n = length(z0)
        g0 = g(z0)
        Hm = zeros(n, n)
        for i = 1:n
            ei = zeros(n)
            ei[i] = h
            Hm[i, i] = (g(z0 .+ ei) - 2g0 + g(z0 .- ei)) / h^2
            for j = (i+1):n
                ej = zeros(n)
                ej[j] = h
                Hm[i, j] =
                    Hm[j, i] =
                        (
                            g(z0 .+ ei .+ ej) - g(z0 .+ ei .- ej) - g(z0 .- ei .+ ej) +
                            g(z0 .- ei .- ej)
                        ) / (4h^2)
            end
        end
        return Hm
    end
    ĥ = Z⃗ -> begin
        unpack!(traj, Z⃗)
        δ = zeros(x_dim * (N - 1))
        evaluate!(δ, ℰ, traj)
        return μ'δ
    end
    H_fd = fd_hess(ĥ, z0)
    @test all(isapprox.(H, triu(H_fd); atol = 2e-3, rtol = 2e-3))

    # structures declare the warp column in packed coordinates
    S = get_jacobian_structure(ℰ, traj)
    @test size(S) == (x_dim * (N - 1), length(vec(traj)))
    @test !iszero(S[slice(1, x_dim), T_col])
    HS = get_hessian_of_lagrangian_structure(ℰ, traj)
    @test size(HS) == (length(vec(traj)), length(vec(traj)))
    @test any(!iszero, HS[T_col, :])
end
