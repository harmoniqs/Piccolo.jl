# ============================================================================ #
# NonHermitianExponentialIntegrator for MultiDensityTrajectory
#
# This integrator exploits the fact that all densities in an ensemble evolve
# under the same Lindbladian propagator. We compute exp(Δt*G) once per knot point
# and apply it to all densities, significantly reducing computation compared to
# separate integrators.
#
# Structural template: hermitian_exponential_integrator_multiket.jl
# Per-density per-segment transformations: nonhermitian_exponential_integrator_density.jl
# ============================================================================ #

"""
    NonHermitianExponentialIntegrator(qtraj::MultiDensityTrajectory, N::Int; kwargs...)

Construct a non-Hermitian exponential integrator for an ensemble of density
matrices evolving under a shared Lindbladian generator. Returns a
`Vector{NonHermitianExponentialIntegrator{MultiDensityTrajectory}}` — one
integrator per density, each sharing the control variables but owning its own
per-segment Jacobian/Hessian and `exp_generator!` buffers.

The state variable for each density is the compact density-matrix isomorphism
`ρ⃗̃ ∈ ℝ^{n²}` (exploiting Hermiticity), and the generator `G(u)` is the compact
Lindbladian `𝒢c(u) ∈ ℝ^{n² × n²}` assembled at call time by
`compact_generator_closure` from the 3-tuple returned by
`compact_lindbladian_parts(sys)` (drift Hamiltonian, per-drive parts,
per-dissipator parts).

# Arguments
- `qtraj::MultiDensityTrajectory`: the ensemble quantum trajectory
- `N::Int`: number of knot points

# Keyword Arguments
- `global_names::Union{Vector{Symbol},Nothing} = nothing`: auto-detected from
  `sys.global_params` if not given.
- `gauss_newton::Bool = false`
"""
function NonHermitianExponentialIntegrator(
    qtraj::MultiDensityTrajectory,
    N::Int;
    global_names::Union{Vector{Symbol},Nothing} = nothing,
    gauss_newton::Bool = false,
)
    sys = get_system(qtraj)
    x_names = state_names(qtraj)
    u = drive_name(qtraj)
    n_densities = length(x_names)

    # Auto-detect globals
    if isnothing(global_names)
        global_names =
            !isempty(sys.global_params) ? collect(keys(sys.global_params)) : Symbol[]
    end

    # Build NamedTrajectory (with globals if present)
    if !isempty(global_names)
        global_data = Dict{Symbol,Vector{Float64}}()
        for name in global_names
            global_data[name] =
                hasproperty(sys, :global_params) && haskey(sys.global_params, name) ?
                [sys.global_params[name]] : [0.0]
        end
        traj = NamedTrajectory(qtraj, N; global_data = global_data)
    else
        traj = NamedTrajectory(qtraj, N)
    end

    @assert traj.N > 1 "Trajectory must have at least two timesteps."
    @assert all(name in traj.names for name in x_names) "All ensemble state names must be in trajectory"

    # Determine global dim from trajectory
    if traj.global_dim > 0
        global_dim = traj.global_dim
    elseif !isempty(global_names)
        global_dim = length(global_names)
    else
        global_dim = 0
    end

    n = sys.levels
    statedim = n^2  # compact density iso dimension per density

    u_dim = traj.dims[u]
    x_dim = traj.dims[x_names[1]]  # per-density state dim (= statedim)
    dim = x_dim * (N - 1)
    var_dim = 2 * x_dim + u_dim + 1

    # Build the shared compact Lindbladian generator G(u) via the
    # drive_coeff + rate_coeff helper. Mirrors Task 9 in the single-density
    # integrator: `compact_lindbladian_parts` returns a 3-tuple
    # (drift_ham, drives, dissipators) in compact real n²×n² form, and
    # `compact_generator_closure` assembles G(u) at call time so typed
    # NonlinearDrive / NonlinearDissipator with global-reading coefficients
    # propagate through the extended control vector. This is the Lindbladian
    # analogue of multiket's shared `sys.G(u, 0.0)`.
    𝒢c_drift_ham, 𝒢c_drives, 𝒢c_dissipators = compact_lindbladian_parts(sys)
    𝒢c_drift_dense = Matrix(𝒢c_drift_ham)
    𝒢c_drives_dense = Matrix{Float64}[Matrix(M) for M in 𝒢c_drives]
    𝒢c_dissipators_dense = Matrix{Float64}[Matrix(M) for M in 𝒢c_dissipators]
    G_fn = compact_generator_closure(
        sys,
        𝒢c_drift_dense,
        𝒢c_drives_dense,
        𝒢c_dissipators_dense,
    )

    x_comps_all = vcat([collect(traj.components[name]) for name in x_names]...)
    u_comps = traj.components[u]
    Δt_comp = traj.components[traj.timestep][1]
    var_comps = [x_comps_all; collect(u_comps); Δt_comp]
    z_dim = traj.dim

    # Return a vector of per-density integrators (mirrors multiket pattern).
    # Each owns its own ∂ℰs, μ∂²ℰs and per-thread G_bufs / expG_bufs so that
    # forward passes and derivative computations for different densities don't
    # alias scratch, AND so that per-knot loops within one integrator can
    # safely run under Threads.@threads.
    nthr = Threads.maxthreadid()
    integrators = NonHermitianExponentialIntegrator{MultiDensityTrajectory}[]
    for (i, name) in enumerate(x_names)
        # Single-density sparsity templates (per density; index into the full
        # trajectory column space via traj.components[name]).
        ∂ℰ_template = jacobian_structure(MultiDensityTrajectory, name, u, statedim, traj)
        μ∂²ℰ_template = hessian_structure(x_dim, u_dim, global_dim)

        ∂ℰs = [copy(∂ℰ_template) for _ = 1:(N-1)]
        μ∂²ℰs = [copy(μ∂²ℰ_template) for _ = 1:(N-1)]

        # Per-thread buffers sized to the compact density iso dimension
        # (n² × n²). Task 9 caveat: `compact_lindbladian_parts` already
        # returns n²×n² real matrices, NOT 2n²×2n². Do not double it.
        G_bufs = [zeros(Float64, statedim, statedim) for _ = 1:nthr]
        expG_bufs = [zeros(Float64, statedim, statedim) for _ = 1:nthr]

        push!(
            integrators,
            NonHermitianExponentialIntegrator{MultiDensityTrajectory}(
                G_fn,
                [name],       # x_names: single density name for this integrator
                u,
                x_dim,
                u_dim,
                var_dim,
                dim,
                statedim,
                ∂ℰs,
                μ∂²ℰs,
                z_dim,
                var_comps,
                nothing,       # Id (unused)
                nothing,       # ∂u∂Δt_bufs
                nothing,       # ∂²u_bufs
                G_bufs,
                expG_bufs,
                global_names,
                global_dim,
                gauss_newton,
            ),
        )
    end

    return integrators
end

# ============================================================================ #
# Jacobian sparsity structure for MultiDensityTrajectory
# (per-density: one state name, but indexes into the full trajectory column space)
# ============================================================================ #

function jacobian_structure(
    ::Type{MultiDensityTrajectory},
    x_name::Symbol,
    u_name::Symbol,
    statedim::Int,
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

    # ∂ρ⃗ₖ₊₁ℰ: Identity
    ∂ℰ[:, z_dim .+ x_comps] = I(x_dim)

    # ∂ρ⃗ₖℰ: Dense block (expG is dense)
    ∂ℰ[:, x_comps] .= 1.0

    # ∂aₖℰ: Dense dependence on controls via G(u)
    ∂ℰ[:, u_comps] .= 1.0

    # ∂Δtₖℰ: Dense dependence via G · expG
    ∂ℰ[:, Δt_comp] .= 1.0

    # ∂gℰ: Global variable derivatives
    if global_dim > 0
        g_cols = (2z_dim+1):(2z_dim+global_dim)
        ∂ℰ[:, g_cols] .= 1.0
    end

    return ∂ℰ
end

# ============================================================================ #
# Forward evaluation (single-density, per-segment)
# ============================================================================ #

@views function (ℰ::NonHermitianExponentialIntegrator{MultiDensityTrajectory})(
    δₖ::AbstractVector,
    zₖ::KnotPoint,
    zₖ₊₁::KnotPoint,
    k::Int,
    globals::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    ρ⃗ₖ = zₖ[x_name(ℰ)]
    ρ⃗ₖ₊₁ = zₖ₊₁[x_name(ℰ)]
    aₖ = zₖ[ℰ.u_name]
    Δtₖ = zₖ.timestep

    # Build extended control with globals
    uₖ = build_extended_control(aₖ, globals)

    # Per-thread buffer indexing — see NonHermitianExponentialIntegrator docstring.
    tid = Threads.threadid()
    G_buf = ℰ.G_bufs[tid]
    expG_buf = ℰ.expG_bufs[tid]

    # δₖ = ρ⃗ₖ₊₁ − exp(Δt·G(u))·ρ⃗ₖ
    copyto!(G_buf, ℰ.G(uₖ))
    exp_generator!(expG_buf, G_buf, Δtₖ)
    copyto!(δₖ, ρ⃗ₖ₊₁)
    mul!(δₖ, expG_buf, ρ⃗ₖ, -1.0, 1.0)

    return nothing
end

# ============================================================================ #
# Jacobian (per-density)
# ============================================================================ #

@views function jacobian!(
    ∂ℰ::AbstractMatrix,
    ℰ::NonHermitianExponentialIntegrator{MultiDensityTrajectory},
    zₖ::KnotPoint,
    zₖ₊₁::KnotPoint,
    k::Int,
    globals::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    ρ⃗ₖ = zₖ[x_name(ℰ)]
    aₖ = zₖ[ℰ.u_name]
    Δtₖ = zₖ.timestep

    uₖ = build_extended_control(aₖ, globals)

    # Per-thread buffer indexing — see NonHermitianExponentialIntegrator docstring.
    tid = Threads.threadid()
    G_buf = ℰ.G_bufs[tid]
    expG_buf = ℰ.expG_bufs[tid]

    Gₖ = ℰ.G(uₖ)
    copyto!(G_buf, Gₖ)
    exp_generator!(expG_buf, G_buf, Δtₖ)
    expGₖ = expG_buf
    expGₖρ⃗ₖ = expGₖ * ρ⃗ₖ  # cache for ∂Δt term

    x_comps = zₖ.components[x_name(ℰ)]
    u_comps = zₖ.components[ℰ.u_name]
    Δt_comp = zₖ.components.Δt[1]

    # ∂ρ⃗ₖℰ: −exp(Δt·G)
    @inbounds ∂ℰ[:, x_comps] .= .-expGₖ

    # ∂aₖℰ: derivative w.r.t. controls — AD through expv(Δt, G(a), ρ⃗)
    ForwardDiff.jacobian!(
        ∂ℰ[:, u_comps],
        a -> -expv(Δtₖ, ℰ.G(build_extended_control(a, globals)), ρ⃗ₖ),
        aₖ,
    )

    # ∂Δtₖℰ: −G · exp(Δt·G) · ρ⃗ₖ
    @inbounds mul!(∂ℰ[:, Δt_comp], Gₖ, expGₖρ⃗ₖ, -1.0, 0.0)

    # ∂gℰ: Global variable derivatives
    if !isnothing(globals) && ℰ.global_dim > 0
        g_cols = (2*ℰ.z_dim+1):(2*ℰ.z_dim+ℰ.global_dim)
        ForwardDiff.jacobian!(∂ℰ[:, g_cols], g -> -expv(Δtₖ, ℰ.G(vcat(aₖ, g)), ρ⃗ₖ), globals)
    end

    # ∂ρ⃗ₖ₊₁ℰ: Identity (already in sparsity pattern)

    return nothing
end

# ============================================================================ #
# Hessian of the Lagrangian (per-density)
# ============================================================================ #

@views function hessian_of_lagrangian!(
    μ∂²ℰ::SparseMatrixCSC,
    ℰ::NonHermitianExponentialIntegrator{MultiDensityTrajectory},
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

    ρ⃗ₖ = zₖ[x_name(ℰ)]
    aₖ = zₖ[ℰ.u_name]
    Δtₖ = zₖ.timestep

    uₖ = build_extended_control(aₖ, globals)

    # Per-thread buffer indexing — see NonHermitianExponentialIntegrator docstring.
    tid = Threads.threadid()
    G_buf = ℰ.G_bufs[tid]
    expG_buf = ℰ.expG_bufs[tid]

    Gₖ = ℰ.G(uₖ)
    copyto!(G_buf, Gₖ)
    exp_generator!(expG_buf, G_buf, Δtₖ)
    expGₖ = expG_buf
    GₖexpGₖ = Gₖ * expGₖ
    GₖexpGₖρ⃗ₖ = GₖexpGₖ * ρ⃗ₖ  # cache for ∂²Δt term

    # Canonical component indices for knot k
    knot_dim = x_dim + u_dim + 1
    x_can = 1:x_dim
    u_can = (x_dim+1):(x_dim+u_dim)
    Δt_can = x_dim + u_dim + 1

    # μₖ ∂ρ⃗ₖ ∂aₖ ℰ
    ForwardDiff.jacobian!(
        μ∂²ℰ[x_can, u_can],
        a -> -expv(Δtₖ, ℰ.G(build_extended_control(a, globals))', μₖ),
        aₖ,
    )

    # μₖ ∂ρ⃗ₖ ∂Δtₖ ℰ
    @inbounds mul!(μ∂²ℰ[x_can, Δt_can], GₖexpGₖ', μₖ, -1.0, 0.0)

    # μₖ ∂aₖ ∂Δtₖ ℰ
    ForwardDiff.gradient!(
        μ∂²ℰ[u_can, Δt_can],
        a ->
            -μₖ' *
            ℰ.G(build_extended_control(a, globals)) *
            expv(Δtₖ, ℰ.G(build_extended_control(a, globals)), ρ⃗ₖ),
        aₖ,
    )

    if has_globals
        g_can = (2*knot_dim+1):(2*knot_dim+global_dim)

        # μₖ ∂ρ⃗ₖ ∂g ℰ — x-g block
        ForwardDiff.jacobian!(
            μ∂²ℰ[x_can, g_can],
            g -> -expv(Δtₖ, ℰ.G(vcat(aₖ, g))', μₖ),
            globals,
        )

        # Full [a,g] Hessian → u-u, u-g, g-g blocks
        full_ag_dim = u_dim + global_dim
        ∂²ag_buf = zeros(full_ag_dim, full_ag_dim)
        ForwardDiff.hessian!(∂²ag_buf, ag -> -μₖ'expv(Δtₖ, ℰ.G(ag), ρ⃗ₖ), vcat(aₖ, globals))
        μ∂²ℰ[u_can, u_can] .= ∂²ag_buf[1:u_dim, 1:u_dim]
        μ∂²ℰ[u_can, g_can] .= ∂²ag_buf[1:u_dim, (u_dim+1):end]
        μ∂²ℰ[g_can, g_can] .= ∂²ag_buf[(u_dim+1):end, (u_dim+1):end]

        # μₖ ∂Δtₖ ∂g ℰ — Δt-g block
        ∂Δt∂g_buf = zeros(global_dim)
        ForwardDiff.gradient!(
            ∂Δt∂g_buf,
            g -> -μₖ' * ℰ.G(vcat(aₖ, g)) * expv(Δtₖ, ℰ.G(vcat(aₖ, g)), ρ⃗ₖ),
            globals,
        )
        for j = 1:global_dim
            μ∂²ℰ[Δt_can, g_can[j]] = ∂Δt∂g_buf[j]
        end
    else
        # μₖ ∂²aₖ ℰ (no globals case)
        ForwardDiff.hessian!(
            μ∂²ℰ[u_can, u_can],
            a -> -μₖ'expv(Δtₖ, ℰ.G(build_extended_control(a, globals)), ρ⃗ₖ),
            aₖ,
        )
    end

    # μₖ ∂²Δtₖ ℰ
    @inbounds μ∂²ℰ[Δt_can, Δt_can] = -μₖ' * Gₖ * GₖexpGₖρ⃗ₖ

    # Symmetrize: canonical→trajectory index mapping can swap triangle positions,
    # so we fill the full symmetric matrix in canonical coords and let the outer
    # eval_hessian_of_lagrangian take triu AFTER the reorder.
    @inbounds for j = 1:u_dim, i = 1:x_dim
        μ∂²ℰ[u_can[j], x_can[i]] = μ∂²ℰ[x_can[i], u_can[j]]
    end
    @inbounds for i = 1:x_dim
        μ∂²ℰ[Δt_can, x_can[i]] = μ∂²ℰ[x_can[i], Δt_can]
    end
    @inbounds for j = 1:u_dim
        μ∂²ℰ[Δt_can, u_can[j]] = μ∂²ℰ[u_can[j], Δt_can]
    end

    # Symmetrize global blocks
    if has_globals
        @inbounds for j = 1:global_dim
            for i = 1:x_dim
                μ∂²ℰ[g_can[j], x_can[i]] = μ∂²ℰ[x_can[i], g_can[j]]
            end
            for i = 1:u_dim
                μ∂²ℰ[g_can[j], u_can[i]] = μ∂²ℰ[u_can[i], g_can[j]]
            end
            μ∂²ℰ[g_can[j], Δt_can] = μ∂²ℰ[Δt_can, g_can[j]]
        end
    end

    return nothing
end

# ============================================================================ #
# Jacobian structure (full-problem, per-density)
# ============================================================================ #

function get_jacobian_structure(
    ℰ::NonHermitianExponentialIntegrator{MultiDensityTrajectory},
    traj::NamedTrajectory,
)
    N = traj.N
    x_dim = traj.dims[x_name(ℰ)]
    z_dim = traj.dim
    global_dim = traj.global_dim
    F_dim = x_dim * (N - 1)
    Z_dim = z_dim * N + global_dim
    ∂F = spzeros(F_dim, Z_dim)

    ∂ℰ_k = jacobian_structure(MultiDensityTrajectory, x_name(ℰ), ℰ.u_name, ℰ.statedim, traj)

    for k = 1:(N-1)
        ∂F[slice(k, x_dim), slice(k, 1:2z_dim, z_dim)] = ∂ℰ_k[:, 1:2z_dim]
    end

    if global_dim > 0
        g_cols_local = (2z_dim+1):(2z_dim+global_dim)
        g_cols_full = (z_dim*N+1):(z_dim*N+global_dim)
        for k = 1:(N-1)
            ∂F[slice(k, x_dim), g_cols_full] = ∂ℰ_k[:, g_cols_local]
        end
    end

    return ∂F
end

# ============================================================================ #
# API methods for MultiDensityTrajectory (per-type dispatch)
#
# Same motivation as Task 9's density dispatch: the abstract-type generalized
# `eval_jacobian` / `eval_hessian_of_lagrangian` in
# `hermitian_exponential_integrator_type.jl` call `triu()` on the per-segment
# Hessian *before* applying the canonical→trajectory index mapping. For
# density-family trajectories the canonical layout is (ρ⃗, u, Δt) and `u`
# appears BEFORE `Δt` in the trajectory — so entries in the canonical upper
# triangle end up in the trajectory lower triangle after the reorder and get
# zeroed by the early triu. Override per-type: map the full symmetric matrix
# first, triu() at the end.
#
# Also: G_bufs / expG_bufs are per-thread Vector{Matrix} with one entry per
# OS thread, indexed by `threadid()` inside the per-knot methods. That
# makes Threads.@threads safe here (same contract as the multiket Hermitian
# dispatch after the thread-race fix).
# ============================================================================ #

@views function evaluate!(
    δ::AbstractVector,
    ℰ::NonHermitianExponentialIntegrator{MultiDensityTrajectory},
    traj::NamedTrajectory,
)
    globals = extract_globals(ℰ, traj)

    # Per-thread buffers on the integrator make this safe under Threads.@threads.
    Threads.@threads for k = 1:(traj.N-1)
        δₖ = δ[slice(k, ℰ.x_dim)]
        ℰ(δₖ, traj[k], traj[k+1], k, globals)
    end
    return nothing
end

@views function eval_jacobian(
    ℰ::NonHermitianExponentialIntegrator{MultiDensityTrajectory},
    traj::NamedTrajectory,
)
    N = traj.N
    x_dim = ℰ.x_dim
    z_dim = traj.dim
    F_dim = x_dim * (N - 1)
    Z_dim = z_dim * N + traj.global_dim

    globals = extract_globals(ℰ, traj)

    # Parallel fill with per-thread G_bufs / expG_bufs.
    Threads.@threads for k = 1:(N-1)
        jacobian!(ℰ.∂ℰs[k], ℰ, traj[k], traj[k+1], k, globals)
    end

    # Build var_comps for the single per-density state name carried by this
    # integrator (MultiDensityTrajectory ensemble: one integrator per density).
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

    if traj.global_dim > 0
        g_cols_local = (2*ℰ.z_dim+1):(2*ℰ.z_dim+traj.global_dim)
        g_cols_full = (z_dim*N+1):(z_dim*N+traj.global_dim)
        @inbounds for k = 1:(N-1)
            ∂F[slice(k, x_dim), g_cols_full] = ℰ.∂ℰs[k][:, g_cols_local]
        end
    end

    return ∂F
end

@views function eval_hessian_of_lagrangian(
    ℰ::NonHermitianExponentialIntegrator{MultiDensityTrajectory},
    traj::NamedTrajectory,
    μ::AbstractVector,
)
    N = traj.N
    x_dim = ℰ.x_dim
    u_dim = ℰ.u_dim
    z_dim = traj.dim
    Z_dim = z_dim * N + traj.global_dim

    globals = extract_globals(ℰ, traj)

    # Parallel fill with per-thread G_bufs / expG_bufs.
    Threads.@threads for k = 1:(N-1)
        μₖ = μ[slice(k, x_dim)]
        hessian_of_lagrangian!(ℰ.μ∂²ℰs[k], ℰ, μₖ, traj[k], traj[k+1], k, globals)
    end

    knot_dim = x_dim + u_dim + 1

    # Canonical indices
    x_can_k = 1:x_dim
    u_can_k = (x_dim+1):(x_dim+u_dim)
    Δt_can_k = x_dim + u_dim + 1
    x_can_k1 = knot_dim .+ (1:x_dim)

    # Trajectory indices — single per-density state name on this integrator
    x_traj_k = collect(traj.components[x_name(ℰ)])
    u_traj_k = collect(traj.components[ℰ.u_name])
    Δt_traj_k = traj.components[traj.timestep][1]
    x_traj_k1 = z_dim .+ x_traj_k

    canonical_comps = [collect(x_can_k); collect(u_can_k); Δt_can_k; collect(x_can_k1)]
    traj_comps = [x_traj_k; u_traj_k; Δt_traj_k; x_traj_k1]

    # Map full symmetric matrix first, THEN extract upper triangle
    # (canonical→trajectory index mapping can swap triangle positions).
    μ∂²F = spzeros(Z_dim, Z_dim)
    @inbounds for k = 1:(N-1)
        μ∂²F[slice(k, traj_comps, z_dim), slice(k, traj_comps, z_dim)] .=
            ℰ.μ∂²ℰs[k][canonical_comps, canonical_comps]
    end

    global_dim = traj.global_dim
    if global_dim > 0
        g_can = (2*knot_dim+1):(2*knot_dim+global_dim)
        g_traj = (z_dim*N+1):(z_dim*N+global_dim)
        @inbounds for k = 1:(N-1)
            μ∂²F[slice(k, traj_comps, z_dim), g_traj] .=
                ℰ.μ∂²ℰs[k][canonical_comps, collect(g_can)]
            μ∂²F[g_traj, slice(k, traj_comps, z_dim)] .=
                ℰ.μ∂²ℰs[k][collect(g_can), canonical_comps]
            μ∂²F[g_traj, g_traj] .+= ℰ.μ∂²ℰs[k][collect(g_can), collect(g_can)]
        end
    end

    return triu(μ∂²F)
end

# ============================================================================ #
# Tests
# ============================================================================ #

@testitem "NonHermitianExponentialIntegrator{MultiDensityTrajectory} construction" begin
    using LinearAlgebra
    using NamedTrajectories
    using Piccolo
    using Piccolo.Control.QuantumIntegrators: NonHermitianExponentialIntegrator

    L = ComplexF64[0 0.1; 0 0]
    sys = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    ρ0₁ = ComplexF64[1 0; 0 0]
    ρg₁ = ComplexF64[0 0; 0 1]
    ρ0₂ = ComplexF64[0 0; 0 1]
    ρg₂ = ComplexF64[1 0; 0 0]
    N = 11
    times = collect(range(0, 1.0, length = N))
    pulse = ZeroOrderPulse(zeros(1, N), times)
    qtraj = MultiDensityTrajectory(sys, pulse, [ρ0₁, ρ0₂], [ρg₁, ρg₂])

    ℰs = NonHermitianExponentialIntegrator(qtraj, N)

    # Returns a vector, one per density
    @test ℰs isa Vector{<:NonHermitianExponentialIntegrator}
    @test length(ℰs) == length(qtraj.initials)
    @test length(ℰs) == 2

    n = sys.levels
    for ℰ in ℰs
        @test ℰ isa NonHermitianExponentialIntegrator{MultiDensityTrajectory}
        @test ℰ.statedim == n^2
        # Per-thread buffer vectors; each entry sized (n², n²).
        nthreads_alloc = length(ℰ.G_bufs)
        @test nthreads_alloc >= Threads.maxthreadid()
        @test length(ℰ.expG_bufs) == nthreads_alloc
        for t = 1:nthreads_alloc
            @test size(ℰ.G_bufs[t]) == (n^2, n^2)
            @test size(ℰ.G_bufs[t]) == size(ℰ.expG_bufs[t])
        end
        @test length(ℰ.x_names) == 1
    end

    # Integrators carry distinct state names (one per density in the ensemble)
    @test ℰs[1].x_names != ℰs[2].x_names
end

@testitem "NonHermitianExponentialIntegrator{MultiDensityTrajectory} test_integrator" begin
    using LinearAlgebra
    using NamedTrajectories
    using Piccolo
    using DirectTrajOpt
    using Piccolo.Control.QuantumIntegrators: NonHermitianExponentialIntegrator

    # Nontrivial dissipator + NONZERO random controls so the initial trajectory's
    # constraint residual is nonzero (test_integrator's `!all(iszero.(f̂))` check
    # requires this). 2-density ensemble: |0⟩⟨0| → |1⟩⟨1| and |1⟩⟨1| → |0⟩⟨0|.
    L = ComplexF64[0 0.1; 0 0]
    sys = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    ρ0₁ = ComplexF64[1 0; 0 0]
    ρg₁ = ComplexF64[0 0; 0 1]
    ρ0₂ = ComplexF64[0 0; 0 1]
    ρg₂ = ComplexF64[1 0; 0 0]
    N = 11
    times = collect(range(0, 1.0, length = N))
    pulse = ZeroOrderPulse(randn(1, N) .* 0.1, times)
    qtraj = MultiDensityTrajectory(sys, pulse, [ρ0₁, ρ0₂], [ρg₁, ρg₂])
    traj = NamedTrajectory(qtraj, N)

    ℰs = NonHermitianExponentialIntegrator(qtraj, N)
    for ℰ in ℰs
        test_integrator(ℰ, traj; atol = 1e-3)
    end
end

@testitem "NonHermitianExponentialIntegrator{MultiDensityTrajectory} globals propagate" begin
    using LinearAlgebra, NamedTrajectories, Piccolo
    using SparseArrays
    using Piccolo.Control.QuantumIntegrators: NonHermitianExponentialIntegrator

    # Use the typed-drive OpenQuantumSystem constructor so that
    # drive_bounds = length 2 dictates n_drives = 2 (the matrix-form
    # constructor would instead use length(H_drives) = 1, which prevents us
    # from carving a u[2] slot for the dissipator-rate global).
    drive = LinearDrive(sparse(ComplexF64.(PAULIS.X)), 1)
    diss = NonlinearDissipator(PAULIS.Z / sqrt(2), u -> u[2]; active_controls = [2])
    sys = OpenQuantumSystem(
        PAULIS.Z,
        AbstractDrive[drive],
        [1.0, 1.0];
        dissipators = [diss],
        global_params = (γ = 0.1,),
    )

    ρ0_list = [ComplexF64[1 0; 0 0], ComplexF64[0 0; 0 1]]
    ρg_list = [ComplexF64[0 0; 0 1], ComplexF64[1 0; 0 0]]
    N = 11
    times = collect(range(0, 1.0, length = N))
    pulse = LinearSplinePulse(repeat([0.3, 0.1], 1, N), times)
    qtraj = MultiDensityTrajectory(sys, pulse, ρ0_list, ρg_list)

    integrators = NonHermitianExponentialIntegrator(qtraj, N)
    # γ slot (u[2]) changing from 0.1 → 0.5 must change every per-density G(u).
    @test all(
        !isapprox(integ.G([0.3, 0.1]), integ.G([0.3, 0.5]); atol = 1e-10) for
        integ in integrators
    )
end
