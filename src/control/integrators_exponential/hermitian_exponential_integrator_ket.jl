# ============================================================================ #
# KetTrajectory Methods
# ============================================================================ #

"""
    HermitianExponentialIntegrator(qtraj::KetTrajectory, N::Int; kwargs...)

Construct a Hermitian exponential integrator for ket state evolution.

# Arguments
- `qtraj::KetTrajectory`: The quantum trajectory with system and pulse
- `N::Int`: Number of knot points for discretization

# Keyword Arguments
- `global_names::Union{Vector{Symbol}, Nothing}=nothing`: Names of global (time-invariant) 
  variables to include in dynamics. If `nothing`, auto-detects from `sys.global_params`.
"""
function HermitianExponentialIntegrator(
    qtraj::KetTrajectory,
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

    return _hermitian_exp_ket(
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

# Per-member inner constructor (#408's sampling lane): all pieces pre-materialized
# (sys/x/u from the member's base trajectory; traj the already-expanded member
# trajectory). The (qtraj, N) outer method is this plus trajectory assembly.
function _hermitian_exp_ket(
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
    ketdim = sys.levels

    # Dimensions for API (u_dim is control dimension, not including globals)
    x_dim = traj.dims[x]
    u_dim = traj.dims[u]
    N = traj.N
    dim = x_dim * (N - 1)

    # Variables: [xₖ, uₖ, Δtₖ, xₖ₊₁]
    var_dim = 2*x_dim + u_dim + 1

    # Use structure functions to create templates with canonical ordering
    ∂ℰ_template = jacobian_structure(KetTrajectory, x, u, ketdim, traj)
    μ∂²ℰ_template = hessian_structure(x_dim, u_dim, global_dim; gauss_newton = gauss_newton)

    # Preallocate one copy per knot point
    ∂ℰs = [copy(∂ℰ_template) for _ = 1:(N-1)]
    μ∂²ℰs = [copy(μ∂²ℰ_template) for _ = 1:(N-1)]

    # Build var_comps
    x_comps = traj.components[x]
    u_comps = traj.components[u]
    Δt_comp = traj.components[traj.timestep][1]
    var_comps = [collect(x_comps); collect(u_comps); Δt_comp]

    z_dim = traj.dim

    # Preallocate per-thread exp_eigen! buffers (7 groups, sized to ketdim).
    # Sized to maxthreadid() so Threads.@threads loops in eval_jacobian /
    # eval_hessian_of_lagrangian can safely index by threadid() without racing.
    nthr = Threads.maxthreadid()
    bufs = _alloc_exp_eigen_bufs(ketdim, nthr)

    # Analytic Daleckii–Krein setup (#202): affine drive directions + per-thread scratch.
    # Populated for all variants at construction; the ket dense wiring lands in slice ④.
    H_dirs = _build_affine_directions(u_ -> sys.H(u_, 0.0), u_dim + global_dim, ketdim)
    dk_bufs = _alloc_dk_bufs(ketdim, nthr)
    dk_so_bufs = _alloc_dk_so_bufs(ketdim, nthr)

    return HermitianExponentialIntegrator{KetTrajectory}(
        u_ -> sys.G(u_, 0.0),
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
        nothing,  # Id (not used for ket)
        nothing,  # ∂u∂Δt_bufs (not used for ket)
        nothing,  # ∂²u_bufs (not used for ket)
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
        matrix_free,  # inert for ket: generic eval_jacobian ignores it → always sparse
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

@views function (ℰ::HermitianExponentialIntegrator{KetTrajectory})(
    δₖ::AbstractVector,
    zₖ::KnotPoint,
    zₖ₊₁::KnotPoint,
    k::Int,
    globals::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    ψ̃ₖ = zₖ[x_name(ℰ)]
    ψ̃ₖ₊₁ = zₖ₊₁[x_name(ℰ)]
    aₖ = zₖ[ℰ.u_name]
    Δtₖ = zₖ.timestep

    # Build extended control with globals
    uₖ = build_extended_control(aₖ, globals)

    # Forward step on the EXACT eigenbasis propagator (exp_eigen!), aligned onto the
    # multiket/unitary variants (#199 "Scope"). Removes the ket's prior forward/derivative
    # method mismatch (forward on Krylov `expv`, derivative on the eigenbasis) — the DK
    # derivative path now reuses the very eigenbasis this forward step computes.
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

    # δₖ = ψ̃ₖ₊₁ - exp(-iΔt H) * ψ̃ₖ
    mul!(δₖ, expG_buf, ψ̃ₖ, -1.0, 0.0)
    @inbounds δₖ .+= ψ̃ₖ₊₁

    return nothing
end

@views function jacobian!(
    ∂ℰ::AbstractMatrix,
    ℰ::HermitianExponentialIntegrator{KetTrajectory},
    zₖ::KnotPoint,
    zₖ₊₁::KnotPoint,
    k::Int,
    globals::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    ψ̃ₖ = zₖ[x_name(ℰ)]
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
    expGₖψ̃ₖ = expGₖ * ψ̃ₖ  # Cache this - used for ∂Δt term

    # Get component indices
    x_comps = zₖ.components[x_name(ℰ)]
    u_comps = zₖ.components[ℰ.u_name]
    Δt_comp = zₖ.components.Δt[1]
    # Use matrix structure dimension, not knot point dimension
    z_dim_matrix = size(∂ℰ, 2) ÷ 2
    x_dim = length(ψ̃ₖ)
    ketdim = ℰ.ketdim

    # Analytic Daleckii–Krein path for affine drives (#204): build the divided-difference
    # matrix M once per knot from the eigenbasis exp_eigen! just produced, then reuse it
    # across directions. `_dk_fill_iso_jac_block!` is reused as-is (ket state is a vector).
    # Falls back to ForwardDiff for nonlinear drives (H_dirs===nothing).
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

    # ∂ψ̃ₖℰ: -exp(ΔtG)
    @inbounds ∂ℰ[:, x_comps] .= .-expGₖ

    # ∂aₖℰ - control derivatives — analytic DK or ForwardDiff witness
    if use_dk
        _iso_vec_to_complex!(dk_vin, ψ̃ₖ, ketdim)
        _dk_fill_iso_jac_block!(
            ∂ℰ[:, u_comps],
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
            ∂ℰ[:, u_comps],
            a -> -expv(Δtₖ, ℰ.G(build_extended_control(a, globals)), ψ̃ₖ),
            aₖ,
        )
    end

    # ∂Δtₖℰ: -G * exp(ΔtG) * ψ̃ₖ
    @inbounds mul!(∂ℰ[:, Δt_comp], Gₖ, expGₖψ̃ₖ, -1.0, 0.0)

    # ∂gℰ: Global variable derivatives
    if !isnothing(globals) && ℰ.global_dim > 0
        g_cols = (2*ℰ.z_dim+1):(2*ℰ.z_dim+ℰ.global_dim)
        if use_dk
            _iso_vec_to_complex!(dk_vin, ψ̃ₖ, ketdim)
            _dk_fill_iso_jac_block!(
                ∂ℰ[:, g_cols],
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
                ∂ℰ[:, g_cols],
                g -> -expv(Δtₖ, ℰ.G(vcat(aₖ, g)), ψ̃ₖ),
                globals,
            )
        end
    end

    # ∂ψ̃ₖ₊₁ℰ: Identity (already in sparsity pattern)

    return nothing
end

function jacobian_structure(
    ::Type{KetTrajectory},
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

    # ∂ψ̃ₖ₊₁ℰ: Identity
    ∂ℰ[:, z_dim .+ x_comps] = I(x_dim)

    # ∂ψ̃ₖℰ: Dense structure for ket
    ∂ℰ[:, x_comps] .= 1.0

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
    ℰ::HermitianExponentialIntegrator{KetTrajectory},
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

    ψ̃ₖ = zₖ[x_name(ℰ)]
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

    # Canonical component indices for knot k
    knot_dim = x_dim + u_dim + 1
    x_can = 1:x_dim
    u_can = (x_dim+1):(x_dim+u_dim)
    Δt_can = x_dim + u_dim + 1
    g_can = has_globals ? ((2*knot_dim+1):(2*knot_dim+global_dim)) : (1:0)

    # Zero out
    @inbounds μ∂²ℰ.nzval .= 0.0

    # Analytic Daleckii–Krein path for affine drives (#204). Build M once per knot from the
    # eigenbasis exp_eigen! produced; `_dk_fill_iso_jac_block!` handles the Gauss-Newton
    # cross-terms (x-u, x-g), and the second-order kernel handles the `!gauss_newton`
    # parameter-parameter blocks. Mirrors the multiket wiring. Falls back to ForwardDiff for
    # nonlinear drives (H_dirs===nothing).
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
    # DK cross-terms apply the ADJOINT of ∂ₚΦ to μₖ (the derivative routine differentiates
    # -iso(Φ†)·μ), mirroring the ForwardDiff-∘-expv witness.
    if use_dk
        _iso_vec_to_complex!(dk_vin, μₖ, ketdim)
    end

    # μₖ∂ψ̃ₖ∂aₖℰ
    if use_dk
        _dk_fill_iso_jac_block!(
            μ∂²ℰ[x_can, u_can],
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
            μ∂²ℰ[x_can, u_can],
            a -> -expv(Δtₖ, ℰ.G(build_extended_control(a, globals))', μₖ),
            aₖ,
        )
    end

    # μₖ∂ψ̃ₖ∂Δtₖℰ
    @inbounds mul!(μ∂²ℰ[x_can, Δt_can], GₖexpGₖ', μₖ, -1.0, 0.0)

    if has_globals
        # μₖ∂ψ̃ₖ∂gℰ - x-g cross-term
        if use_dk
            _dk_fill_iso_jac_block!(
                μ∂²ℰ[x_can, g_can],
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
                μ∂²ℰ[x_can, g_can],
                g -> -expv(Δtₖ, ℰ.G(vcat(aₖ, g))', μₖ),
                globals,
            )
        end
    end

    # === Parameter-parameter blocks (skipped in GN) ===
    if !ℰ.gauss_newton
        GₖexpGₖψ̃ₖ = GₖexpGₖ * ψ̃ₖ  # Cache for ∂²Δt term

        if use_dk
            # ── Analytic exact-Hessian parameter-parameter blocks (#204) ──
            # Every p-p block is -∂²/∂p∂q[Re(χ† Φ ψ)] with χ = iso⁻¹(μₖ), ψ = iso⁻¹(ψ̃ₖ);
            # μ_mat = χ ψ† realises the Re tr(μ†·) pairing the second-order kernel expects
            # (matching the ForwardDiff-∘-expv witness). `dk_vin` already holds χ.
            dk_so_ws = ℰ.dk_so_ws_bufs[tid]
            dk_μmat = ℰ.dk_μmat_bufs[tid]
            dk_psi = ℰ.dk_psi_bufs[tid]
            dk_y = ℰ.dk_y_bufs[tid]
            dk_Hc = ℰ.dk_Hc_bufs[tid]
            dk_B = zeros(u_dim + global_dim, u_dim + global_dim)

            _iso_vec_to_complex!(dk_psi, ψ̃ₖ, ketdim)
            dk_μmat .= dk_vin .* dk_psi'      # χ ψ†
            # y = Φ ψ = V (cis(-Δt λ) ∘ (V† ψ))   (reuse exp_eigen! eigenbasis)
            mul!(dk_vout, V_buf', dk_psi)
            dk_vout .*= cis_diag_buf
            mul!(dk_y, V_buf, dk_vout)
            mul!(dk_Hc, H_buf, dk_vin)        # H χ

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
            # so the block value is  -Re(χ† ∂²Φ/∂aₚ∂Δt ψ)
            #                       = Re(i·[χ† Hₚ(Φψ) + (Hχ)† (∂ₚΦ)ψ]).
            for p = 1:(u_dim+global_dim)
                Hₚ = ℰ.H_dirs[p]
                dk_apply!(dk_dΦ, V_buf, dk_ws.M, Hₚ, dk_ws.tmp1, dk_ws.tmp2)
                mul!(dk_vout, Hₚ, dk_y)          # Hₚ (Φ ψ)
                s = dot(dk_vin, dk_vout)         # χ† Hₚ Φ ψ
                mul!(dk_vout, dk_dΦ, dk_psi)     # (∂ₚΦ) ψ
                s += dot(dk_Hc, dk_vout)         # (H χ)† (∂ₚΦ) ψ
                val = real(im * s)
                if p <= u_dim
                    μ∂²ℰ[u_can[p], Δt_can] += val
                else
                    μ∂²ℰ[Δt_can, g_can[p-u_dim]] += val
                end
            end
        else
            # μₖ∂aₖ∂Δtₖℰ
            ForwardDiff.gradient!(
                μ∂²ℰ[u_can, Δt_can],
                a ->
                    -μₖ' *
                    ℰ.G(build_extended_control(a, globals)) *
                    expv(Δtₖ, ℰ.G(build_extended_control(a, globals)), ψ̃ₖ),
                aₖ,
            )

            if has_globals
                # Full [a,g] Hessian → u-u, u-g, g-g blocks
                full_ag_dim = u_dim + global_dim
                ∂²ag_buf = zeros(full_ag_dim, full_ag_dim)
                ForwardDiff.hessian!(
                    ∂²ag_buf,
                    ag -> -μₖ'expv(Δtₖ, ℰ.G(ag), ψ̃ₖ),
                    vcat(aₖ, globals),
                )
                μ∂²ℰ[u_can, u_can] .= ∂²ag_buf[1:u_dim, 1:u_dim]
                μ∂²ℰ[u_can, g_can] .= ∂²ag_buf[1:u_dim, (u_dim+1):end]
                μ∂²ℰ[g_can, g_can] .= ∂²ag_buf[(u_dim+1):end, (u_dim+1):end]

                # μₖ∂Δtₖ∂gℰ - Δt-g block
                ∂Δt∂g_buf = zeros(global_dim)
                ForwardDiff.gradient!(
                    ∂Δt∂g_buf,
                    g -> -μₖ' * ℰ.G(vcat(aₖ, g)) * expv(Δtₖ, ℰ.G(vcat(aₖ, g)), ψ̃ₖ),
                    globals,
                )
                for j = 1:global_dim
                    μ∂²ℰ[Δt_can, g_can[j]] = ∂Δt∂g_buf[j]
                end
            else
                # μₖ∂²aₖℰ (no globals case)
                ForwardDiff.hessian!(
                    μ∂²ℰ[u_can, u_can],
                    a -> -μₖ'expv(Δtₖ, ℰ.G(build_extended_control(a, globals)), ψ̃ₖ),
                    aₖ,
                )
            end
        end

        # μₖ∂²Δtₖℰ (closed form, analytic in both paths)
        @inbounds μ∂²ℰ[Δt_can, Δt_can] = -μₖ' * Gₖ * GₖexpGₖψ̃ₖ
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
        g_can = (2*knot_dim+1):(2*knot_dim+global_dim)
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

function get_jacobian_structure(
    ℰ::HermitianExponentialIntegrator{KetTrajectory},
    traj::NamedTrajectory,
)
    # Derived-Δt dynamics (Piccolo.jl#321): packed two-knot windows + the warp
    # column under a warp; historical structure otherwise.
    traj.warp !== nothing && return _get_jacobian_structure_warped(ℰ, traj)

    N = traj.N
    x_dim = traj.dims[x_name(ℰ)]
    z_dim = traj.dim
    global_dim = traj.global_dim
    F_dim = x_dim * (N - 1)
    Z_dim = z_dim * N + global_dim
    ∂F = spzeros(F_dim, Z_dim)

    # Get structure for a single knot point (includes global columns)
    ∂ℰ_k = jacobian_structure(KetTrajectory, x_name(ℰ), ℰ.u_name, ℰ.ketdim, traj)

    # Place per-knot structure for all knot points
    for k = 1:(N-1)
        ∂F[slice(k, x_dim), slice(k, 1:2z_dim, z_dim)] = ∂ℰ_k[:, 1:2z_dim]
    end

    # Global columns: all knot points contribute to the same global columns. Map the
    # per-knot block (global_names order) to the trajectory's global_components
    # order via `_global_full_cols` so the STRUCTURE matches the VALUES placement
    # in `eval_jacobian` (fixes the permuted-∂c/∂θ bug).
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
# API methods for KetTrajectory (explicit dispatch)
# ============================================================================ #

@views function evaluate!(
    δ::AbstractVector,
    ℰ::HermitianExponentialIntegrator{KetTrajectory},
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
    ℰ::HermitianExponentialIntegrator{KetTrajectory},
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

    # Global columns: copy from per-knot ∂ℰ to their Z positions. The per-knot
    # block is in global_names order; `_global_full_cols` maps it to the
    # trajectory's global_components order (fixes the permuted-∂c/∂θ bug).
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
    ℰ::HermitianExponentialIntegrator{KetTrajectory},
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

    # Global blocks: cross-terms with per-knot vars + g-g block. `g_can` is the
    # canonical (global_names-order) block; `g_traj` maps it to the trajectory's
    # global_components order in Z (see `_global_full_cols`) — fixes the permuted
    # global Hessian blocks under a reordering.
    global_dim = traj.global_dim
    if global_dim > 0
        g_can = (2*knot_dim+1):(2*knot_dim+global_dim)
        g_traj = _global_full_cols(ℰ, traj)
        @inbounds for k = 1:(N-1)
            μ∂²F[slice(k, traj_comps, z_dim), g_traj] .=
                ℰ.μ∂²ℰs[k][canonical_comps, collect(g_can)]
            μ∂²F[g_traj, slice(k, traj_comps, z_dim)] .=
                ℰ.μ∂²ℰs[k][collect(g_can), canonical_comps]
            μ∂²F[g_traj, g_traj] .+= ℰ.μ∂²ℰs[k][collect(g_can), collect(g_can)]
        end
    end

    # Return upper triangle (symmetric matrix, only upper needed by Ipopt)
    return triu(μ∂²F)
end

@testitem "testing HermitianExponentialIntegrator with KetTrajectory" begin
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
    ψ_init = Id[:, 1]
    ψ_goal = Id[:, 2]

    T = 1.0
    N = 6
    qtraj = KetTrajectory(sys, ψ_init, ψ_goal, T)
    traj = NamedTrajectory(qtraj, N)

    ℰ = HermitianExponentialIntegrator(qtraj, N)
    @test ℰ isa HermitianExponentialIntegrator{KetTrajectory}

    test_integrator(ℰ, traj; atol = 1e-3, show_hessian_diff = true)
end





@testitem "HermitianExponentialIntegrator{Ket} global Jacobian columns respect trajectory order (permutation regression)" begin
    # Regression for the permuted-∂c/∂θ bug. The integrator fills its per-knot
    # global-derivative block in `global_names` order (= the system's
    # `global_params` order, which the forward model requires), but the solver
    # indexes the assembled Jacobian's global columns in the TRAJECTORY's
    # `global_components` order. `NamedTrajectory` sorts globals alphabetically, so
    # a system whose `global_params` are NOT alphabetical produces
    # `global_names ≠ trajectory order` — and a naive i→i placement then SWAPS
    # ∂c/∂a and ∂c/∂b, handing Ipopt a wrong Jacobian (calibration diverges to a
    # bound). Here `global_params=(b, a)` (non-alphabetical) forces exactly that
    # mismatch while keeping the forward model valid; each global's assembled
    # column, at its trajectory position, must match a finite difference of the
    # forward residual.
    using DirectTrajOpt
    using Piccolo
    using Piccolo.Control.QuantumIntegrators.ExponentialIntegrators
    using NamedTrajectories
    using LinearAlgebra
    using Random

    Random.seed!(90_210)
    # globals enter DISTINCT operators so a swap is observable: ∂H/∂b=Z, ∂H/∂a=Y.
    # global_params=(b, a) ⇒ integrator global_names=[:b,:a]; H reads u[3]=b, u[4]=a.
    H = (u, t) -> u[3] * GATES.Z + u[4] * GATES.Y + u[1] * GATES.X + u[2] * GATES.Y
    sys = OpenQuantumSystem(
        H,
        [1.0, 1.0];
        time_dependent = true,
        global_params = (b = 0.2, a = 0.1),
    )
    ψ0 = ComplexF64[1.0, 0.0]
    ψg = ComplexF64[0.0, 1.0]
    T = 2.0
    N = 8
    qtraj = KetTrajectory(sys, ψ0, ψg, T)

    ℰ = HermitianExponentialIntegrator(qtraj, N)
    @test ℰ.global_names == [:b, :a]           # follows global_params order

    traj = NamedTrajectory(qtraj, N)
    # NamedTrajectory sorts globals alphabetically ⇒ [:a,:b] ≠ global_names ⇒ the
    # permutation is genuinely exercised (guard against the test silently passing).
    @test traj.global_components[:a][1] < traj.global_components[:b][1]

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

@testitem "HermitianExponentialIntegrator{KetTrajectory} forward pass minimal alloc" begin
    using LinearAlgebra
    using NamedTrajectories
    using Piccolo
    using BenchmarkTools

    sys = OpenQuantumSystem(PAULIS.Z, [PAULIS.X, PAULIS.Y], [1.0, 1.0])
    N = 11
    times = collect(range(0, 1.0, length = N))
    pulse = LinearSplinePulse(zeros(2, N), times)
    qtraj = KetTrajectory(sys, pulse, ComplexF64[1, 0], ComplexF64[0, 1])
    traj = NamedTrajectory(qtraj, N)

    ℰ = HermitianExponentialIntegrator(qtraj, N)
    δ = zeros(ℰ.x_dim)
    z1 = traj[1]
    z2 = traj[2]

    # Warmup
    ℰ(δ, z1, z2, 1)

    # The ket forward path is aligned onto `exp_eigen!` (#204): it calls exp_eigen! with
    # the preallocated per-thread struct buffers (no per-call allocation for the
    # eigendecomposition scratch), matching the multiket/unitary forward variants.
    # Residual allocation comes from ℰ.H(uₖ) building a fresh matrix each call plus the
    # LAPACK `eigen!` internals (F.values/F.vectors/work arrays); same limitation as the
    # unitary forward path (~4176 B observed there). Threshold mirrors the unitary bound
    # (< 8192) pending a follow-on H-buffer / low-level syevr! refactor.
    allocs = @ballocated $ℰ($δ, $z1, $z2, 1)
    @test allocs < 8192
end

# ============================================================================ #
# NonlinearDrive Tests
# ============================================================================ #

@testitem "HermitianExponentialIntegrator KetTrajectory with NonlinearDrive" begin
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

    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]

    T = 1.0
    N = 6
    qtraj = KetTrajectory(sys, ψ_init, ψ_goal, T)
    traj = NamedTrajectory(qtraj, N)

    ℰ = HermitianExponentialIntegrator(qtraj, N)
    @test ℰ isa HermitianExponentialIntegrator{KetTrajectory}

    test_integrator(ℰ, traj; atol = 1e-3, show_hessian_diff = true)
end



@testitem "HermitianExponentialIntegrator KetTrajectory with NonlinearDrive and global variables" begin
    using DirectTrajOpt
    using DirectTrajOpt: BoundsConstraint
    using Piccolo.Control.QuantumIntegrators.ExponentialIntegrators
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

    sys = OpenQuantumSystem(PAULIS.Z, drives, drive_bounds; global_params = (δ = δ_init,))

    @test has_nonlinear_drives(sys.H_drives)

    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]

    T = 10.0
    N = 15
    qtraj = KetTrajectory(sys, ψ_init, ψ_goal, T)

    integrator = HermitianExponentialIntegrator(qtraj, N)

    @test integrator.global_names == [:δ]
    @test integrator.global_dim == 1

    traj = NamedTrajectory(qtraj, N)
    test_integrator(integrator, traj; atol = 1e-3, show_hessian_diff = true)
end

# ============================================================================ #
# Gauss-Newton Hessian Tests
# ============================================================================ #

@testitem "HermitianExponentialIntegrator KetTrajectory Gauss-Newton Hessian" begin
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
    ψ_init = Id[:, 1]
    ψ_goal = Id[:, 2]

    T = 1.0
    N = 6
    qtraj = KetTrajectory(sys, ψ_init, ψ_goal, T)
    traj = NamedTrajectory(qtraj, N)

    # Test 1: GN integrator constructs correctly
    ℰ_gn = HermitianExponentialIntegrator(qtraj, N; gauss_newton = true)
    @test ℰ_gn isa HermitianExponentialIntegrator{KetTrajectory}
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
    # The exact integrator's (p,p) positions should be nonzero in exact but zero in GN
    μ∂²f_exact = eval_hessian_of_lagrangian(ℰ_exact, traj, μ)
    # Check that GN has zeros where exact has (p,p) values
    gn_mask = triu(H_gn) .!= 0
    exact_mask = triu(H_exact) .!= 0
    pp_mask = exact_mask .& .!gn_mask  # positions in exact but not in GN
    @test any(pp_mask)  # there should be some (p,p) entries
    @test all(triu(μ∂²f_gn)[pp_mask] .== 0.0)  # GN should have zeros there
    @test !all(triu(μ∂²f_exact)[pp_mask] .== 0.0)  # exact should be nonzero
    println("  (p,p) block entries: $(count(pp_mask)), all zero in GN: true")
end



# ============================================================================ #
# Analytic Daleckii–Krein dense-wiring tests (Piccolissimo.jl#204, slice ⑤).
# Mirror the multiket #202/#203 DK tests for the KetTrajectory variant: the
# analytic Jacobian, Gauss-Newton Hessian cross-terms, and `!gauss_newton`
# parameter-parameter blocks vs the retained ForwardDiff witness + an independent
# FiniteDiff oracle, per-thread safety, the flag toggle, and the nonlinear-drive
# fallback. Also the ket forward-step alignment onto exp_eigen! (AC3).
# ============================================================================ #

@testitem "Ket DK Jacobian matches ForwardDiff witness — affine drive [#204 AC2]" begin
    using DirectTrajOpt
    using Piccolo.Control.QuantumIntegrators.ExponentialIntegrators
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using Random

    Random.seed!(204_101)
    T = 2.0
    N = 8
    sys = OpenQuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    qtraj = KetTrajectory(sys, ψ_init, ψ_goal, T)

    ℰ_dk = HermitianExponentialIntegrator(qtraj, N; use_analytical = true)
    ℰ_fd = HermitianExponentialIntegrator(qtraj, N; use_analytical = false)

    @test ℰ_dk.H_dirs !== nothing
    @test length(ℰ_dk.H_dirs) == 2          # two linear drives, no globals
    @test ℰ_dk.use_analytical == true
    @test ℰ_fd.use_analytical == false

    traj = NamedTrajectory(qtraj, N)
    traj.datavec .= randn(length(traj.datavec))

    Jdk = eval_jacobian(ℰ_dk, traj)
    Jfd = eval_jacobian(ℰ_fd, traj)
    relerr = norm(Jdk - Jfd) / norm(Jfd)
    println("  [#204 AC2] Ket Jacobian DK-vs-ForwardDiff rel err = $relerr")
    @test relerr < 1e-9
end

@testitem "Ket DK Gauss-Newton Hessian cross-terms match ForwardDiff [#204 AC2]" begin
    using DirectTrajOpt
    using Piccolo.Control.QuantumIntegrators.ExponentialIntegrators
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using Random

    Random.seed!(204_102)
    T = 2.0
    N = 7

    # --- no globals: exercises the x-u cross-term block ---
    sys = OpenQuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    qtraj = KetTrajectory(sys, ψ_init, ψ_goal, T)

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
    println("  [#204 AC2] Ket GN Hessian (x-u) DK-vs-ForwardDiff rel err = $relerr")
    @test relerr < 1e-9

    # --- with globals: additionally exercises the x-g cross-term block ---
    H = (u, t) -> (u[3] + u[4]) * GATES.Z + u[1] * GATES.X + u[2] * GATES.Y
    gsys = OpenQuantumSystem(
        H,
        [1.0, 1.0];
        time_dependent = true,
        global_params = (δ = 0.2, ω = 1.0),
    )
    gqtraj = KetTrajectory(gsys, ψ_init, ψ_goal, T)
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
    println("  [#204 AC2] Ket GN Hessian (x-u & x-g) DK-vs-ForwardDiff rel err = $relerr_g")
    @test relerr_g < 1e-9
end

@testitem "Ket DK exact-Hessian p-p blocks match ForwardDiff witness — affine [#204 AC2]" begin
    using DirectTrajOpt
    using Piccolo.Control.QuantumIntegrators.ExponentialIntegrators
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using SparseArrays
    using Random

    Random.seed!(204_103)
    T = 2.0
    N = 7

    # --- no globals: u-u, u-Δt, Δt-Δt parameter-parameter blocks ---
    sys = OpenQuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    qtraj = KetTrajectory(sys, ψ_init, ψ_goal, T)

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
        "  [#204 AC2] Ket exact-Hessian (no globals) DK-vs-ForwardDiff rel err = $relerr",
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
    gqtraj = KetTrajectory(gsys, ψ_init, ψ_goal, T)
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
        "  [#204 AC2] Ket exact-Hessian (globals) DK-vs-ForwardDiff rel err = $relerr_g",
    )
    @test relerr_g < 1e-8
end

@testitem "Ket DK path is thread-safe: parallel == serial [#204 AC4]" begin
    using DirectTrajOpt
    using Piccolo.Control.QuantumIntegrators.ExponentialIntegrators
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using Random
    using TrajectoryIndexingUtils: slice

    Random.seed!(204_104)
    T = 2.0
    N = 12

    H = (u, t) -> (u[3] + u[4]) * GATES.Z + u[1] * GATES.X + u[2] * GATES.Y
    sys = OpenQuantumSystem(
        H,
        [1.0, 1.0];
        time_dependent = true,
        global_params = (δ = 0.2, ω = 1.0),
    )
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    qtraj = KetTrajectory(sys, ψ_init, ψ_goal, T)

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

    println("  [#204 AC4] Ket running on $(Threads.nthreads()) thread(s)")

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

@testitem "Ket DK path matches FiniteDiff oracle [#204 AC2]" begin
    using DirectTrajOpt
    using Piccolo.Control.QuantumIntegrators.ExponentialIntegrators
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using Random

    include("../../../test/test_utils.jl")

    Random.seed!(204_105)
    T = 2.0
    N = 6
    sys = OpenQuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    qtraj = KetTrajectory(sys, ψ_init, ψ_goal, T)
    traj = NamedTrajectory(qtraj, N)

    # Gauss-Newton analytic path vs FiniteDiff oracle
    ℰ_gn =
        HermitianExponentialIntegrator(qtraj, N; use_analytical = true, gauss_newton = true)
    test_integrator(ℰ_gn, traj; gauss_newton = true, atol = 1e-4, show_hessian_diff = true)

    # Exact-Hessian analytic path vs FiniteDiff oracle (no autodiff in the path)
    ℰ_ex = HermitianExponentialIntegrator(
        qtraj,
        N;
        use_analytical = true,
        gauss_newton = false,
    )
    test_integrator(ℰ_ex, traj; gauss_newton = false, atol = 1e-4, show_hessian_diff = true)
end

@testitem "Ket nonlinear drive falls back to ForwardDiff [#204 AC5]" begin
    using DirectTrajOpt
    using Piccolo.Control.QuantumIntegrators.ExponentialIntegrators
    using Piccolo
    using SparseArrays
    using NamedTrajectories
    using LinearAlgebra
    using Random

    include("../../../test/test_utils.jl")

    Random.seed!(204_106)
    # H = σz + u₁·σx + u₁²·σz — nonlinear in the control.
    drives = AbstractDrive[
        LinearDrive(sparse(ComplexF64.(PAULIS.X)), 1),
        NonlinearDrive(PAULIS.Z, u -> u[1]^2),
    ]
    sys = OpenQuantumSystem(PAULIS.Z, drives, [1.0])
    @test has_nonlinear_drives(sys.H_drives)

    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    T = 1.0
    N = 6
    qtraj = KetTrajectory(sys, ψ_init, ψ_goal, T)

    # Even with the analytic flag ON, a nonlinear drive has no constant direction matrix
    # ⇒ H_dirs is nothing ⇒ the DK path is unavailable and ForwardDiff runs.
    ℰ_on = HermitianExponentialIntegrator(qtraj, N; use_analytical = true)
    ℰ_off = HermitianExponentialIntegrator(qtraj, N; use_analytical = false)
    @test ℰ_on.H_dirs === nothing
    @test ℰ_on.use_analytical == true

    traj = NamedTrajectory(qtraj, N)
    traj.datavec .= randn(length(traj.datavec))

    # Flag-on and flag-off produce IDENTICAL results (both take the ForwardDiff path).
    @test eval_jacobian(ℰ_on, traj) == eval_jacobian(ℰ_off, traj)

    traj2 = NamedTrajectory(qtraj, N)
    test_integrator(ℰ_on, traj2; atol = 1e-3, show_hessian_diff = true)
end

@testitem "Ket forward step on exp_eigen! matches prior expv forward [#204 AC3]" begin
    using Piccolo.Control.QuantumIntegrators.ExponentialIntegrators
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using ExponentialAction: expv
    using Random

    Random.seed!(204_107)

    # Random Hermitian affine system; the aligned exp_eigen! forward must reproduce the
    # prior Krylov `expv` forward to tight tolerance (no physics change, just method).
    n = 4
    A0 = randn(ComplexF64, n, n)
    H0 = Matrix(Hermitian(A0 + A0'))
    A1 = randn(ComplexF64, n, n)
    H1 = Matrix(Hermitian(A1 + A1'))
    A2 = randn(ComplexF64, n, n)
    H2 = Matrix(Hermitian(A2 + A2'))
    H = (u, t) -> H0 + u[1] * H1 + u[2] * H2
    sys = OpenQuantumSystem(H, [1.0, 1.0]; time_dependent = true)

    ψ_init = normalize(randn(ComplexF64, n))
    ψ_goal = normalize(randn(ComplexF64, n))
    T = 2.0
    N = 8
    qtraj = KetTrajectory(sys, ψ_init, ψ_goal, T)
    traj = NamedTrajectory(qtraj, N)
    traj.datavec .= randn(length(traj.datavec))

    ℰ = HermitianExponentialIntegrator(qtraj, N)

    # `@testitem` runs in a module (hard scope): reassigning an outer variable inside a
    # for-loop makes it loop-local (RHS read → UndefVarError). Accumulate by mutating a
    # vector (push! is mutation, not reassignment), then reduce after the loop.
    relerrs = Float64[]
    for k = 1:(N-1)
        zₖ = traj[k]
        zₖ₊₁ = traj[k+1]
        aₖ = zₖ[ℰ.u_name]
        Δtₖ = zₖ.timestep
        ψ̃ₖ = zₖ[state_name(qtraj)]
        ψ̃ₖ₊₁ = zₖ₊₁[state_name(qtraj)]

        # New (aligned) exp_eigen! forward
        δ_new = zeros(ℰ.x_dim)
        ℰ(δ_new, zₖ, zₖ₊₁, k)

        # Prior Krylov `expv` forward: δ = ψ̃ₖ₊₁ - expv(Δt, G, ψ̃ₖ)
        δ_expv = ψ̃ₖ₊₁ .- expv(Δtₖ, ℰ.G(aₖ), ψ̃ₖ)

        push!(relerrs, norm(δ_new - δ_expv) / max(norm(δ_expv), 1.0))
    end
    max_relerr = maximum(relerrs)
    println("  [#204 AC3] Ket forward exp_eigen!-vs-expv max rel err = $max_relerr")
    @test max_relerr < 1e-10
end

# ============================================================================ #
# Phase 1b (Piccolo.jl#321): the defect warp column under a time warp
# (NamedTrajectories#161). Under a warp the timestep rows are derived
# (Δtₖ = Δtₖ(T)); the packed Jacobian drops the derived Δt column and gains
# the warp-parameter column with the EXACT chain entries
# ∂cₖ/∂θⱼ = (∂Δtₖ/∂θⱼ)·∂cₖ/∂Δtₖ — for GlobalScale, wₖ·∂cₖ/∂Δtₖ with
# wₖ = 1/(N-1), and ∂cₖ/∂Δtₖ = −G(uₖ)·xₖ₊₁ on the satisfied defect.
# ============================================================================ #

@testitem "HermitianExponentialIntegrator{Ket} defect warp column: exact entries + FD parity" begin
    using Piccolo
    using Piccolo.Control.QuantumIntegrators.ExponentialIntegrators:
        _hermitian_exp_ket, exp_eigen
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
    ψ̃ = randn(4, N)
    u = 0.3 .* randn(2, N)

    mk_traj(ψ̃) = NamedTrajectory(
        (ψ̃ = ψ̃, u = u, Δt = fill(T0 / (N - 1), 1, N));  # Δt ignored — derived
        controls = (:u,),
        timestep = :Δt,
        warp = GlobalScale(T0),
    )
    traj = mk_traj(ψ̃)
    ℰ = _hermitian_exp_ket(sys, :ψ̃, :u, traj, Symbol[])

    x_dim = ℰ.x_dim                      # 4 (iso ket)
    wₖ = 1 / (N - 1)                     # GlobalScale lattice weight
    T_col = traj.dim * N - traj.dims[:Δt] * N + traj.global_dim + 1  # packed warp base + 1

    # evaluate! reads the DERIVED rows (synced from the warp)
    δ = zeros(x_dim * (N - 1))
    evaluate!(δ, ℰ, traj)
    @test !all(iszero, δ)
    for k = 1:(N-1)
        @test δ[((k-1)*x_dim) .+ (1:x_dim)] ≈
              ψ̃[:, k+1] - exp_eigen((T0 / (N - 1)) .* sys.H(u[:, k], 0.0)) * ψ̃[:, k]
    end

    # Jacobian: packed size, and the warp column is the EXACT chain column
    ∂F = eval_jacobian(ℰ, traj)
    @test size(∂F) == (x_dim * (N - 1), length(vec(traj)))
    for k = 1:(N-1)
        uₖ = u[:, k]
        Gₖ = sys.G(uₖ, 0.0)
        Φₖ = exp_eigen((T0 / (N - 1)) .* sys.H(uₖ, 0.0))
        dcdΔt = -Gₖ * Φₖ * ψ̃[:, k]                       # ∂cₖ/∂Δtₖ (iso real)
        @test ∂F[((k-1)*x_dim) .+ (1:x_dim), T_col] ≈ wₖ .* dcdΔt
    end
    # the packed width excludes the derived Δt rows and appends the warp param T
    @test length(vec(traj)) ==
          (traj.dim - traj.dims[:Δt]) * traj.N + traj.global_dim + n_params(traj.warp)

    # satisfied defect: ∂cₖ/∂Δtₖ = −G(uₖ)·xₖ₊₁ exactly (the issue's form)
    ψ̃sat = copy(ψ̃)
    for k = 1:(N-1)
        ψ̃sat[:, k+1] = exp_eigen((T0 / (N - 1)) .* sys.H(u[:, k], 0.0)) * ψ̃sat[:, k]
    end
    traj_sat = mk_traj(ψ̃sat)
    ℰ_sat = _hermitian_exp_ket(sys, :ψ̃, :u, traj_sat, Symbol[])
    ∂F_sat = eval_jacobian(ℰ_sat, traj_sat)
    for k = 1:(N-1)
        @test ∂F_sat[((k-1)*x_dim) .+ (1:x_dim), T_col] ≈
              -wₖ .* (sys.G(u[:, k], 0.0) * ψ̃sat[:, k+1])
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
