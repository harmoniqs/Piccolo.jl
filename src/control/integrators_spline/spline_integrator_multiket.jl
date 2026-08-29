# ============================================================================ #
# SplineIntegrator for MultiKetTrajectory
# 
# This integrator exploits the fact that all kets in an ensemble evolve under
# the same propagator. We compute Φ(p) once per knot point and apply it to all
# kets, significantly reducing computation compared to separate integrators.
# ============================================================================ #

# ============================================================================ #
# Constructor for MultiKetTrajectory
# ============================================================================ #

"""
    SplineIntegrator(qtraj::MultiKetTrajectory, N::Int; kwargs...)

Construct a SplineIntegrator from an MultiKetTrajectory.
Automatically builds the combined trajectory with separate state variables.

# Arguments
- `qtraj::MultiKetTrajectory`: The ensemble quantum trajectory
- `N::Int`: Number of knot points for discretization
"""
function SplineIntegrator(
    qtraj::MultiKetTrajectory,
    N::Int;
    spline_order::Union{Int,Nothing} = nothing,
    alg::IntegrationAlgorithm = Tsit5Alg(),
    global_names::Union{Vector{Symbol},Nothing} = nothing,
    ket_sensitivity::Union{Bool,Nothing} = nothing,
    exact_hessian::Bool = false,
    second_order_shape::Symbol = :auto,
)
    traj = NamedTrajectory(qtraj, N)
    return SplineIntegrator(
        qtraj,
        traj;
        spline_order = spline_order,
        alg = alg,
        global_names = global_names,
        ket_sensitivity = ket_sensitivity,
        exact_hessian = exact_hessian,
        second_order_shape = second_order_shape,
    )
end

# Per-member inner constructor (#408's sampling lane): all pieces
# pre-materialized. The inner (qtraj, traj) method derives pulse/sys/x_names/u
# from qtraj — identical information to what this helper receives.
function _spline_multiket(
    sys,
    pulse,
    x_names::Vector{Symbol},
    u::Symbol,
    traj::NamedTrajectory;
    spline_order::Union{Int,Nothing} = nothing,
    alg::IntegrationAlgorithm = Tsit5Alg(),
    global_names::Union{Vector{Symbol},Nothing} = nothing,
    ket_sensitivity::Union{Bool,Nothing} = nothing,
    exact_hessian::Bool = false,
    second_order_shape::Symbol = :auto,
)
    return SplineIntegrator(
        _mk_multiket_qtraj(sys, pulse, x_names, u, traj),
        traj;
        spline_order = spline_order,
        alg = alg,
        global_names = global_names,
        ket_sensitivity = ket_sensitivity,
        exact_hessian = exact_hessian,
        second_order_shape = second_order_shape,
    )
end

# Minimal MultiKetTrajectory-shaped shim: the inner constructor reads only
# get_pulse/get_system/state_names/drive_name from it.
struct _SamplingMemberQTraj{S,P,X<:Vector{Symbol},U<:Symbol,T}
    sys::S
    pulse::P
    x_names::X
    u::U
    traj::T
end
get_system(q::_SamplingMemberQTraj) = q.sys
get_pulse(q::_SamplingMemberQTraj) = q.pulse
state_names(q::_SamplingMemberQTraj) = q.x_names
drive_name(q::_SamplingMemberQTraj) = q.u
Base.length(q::_SamplingMemberQTraj) = length(q.x_names)

_mk_multiket_qtraj(sys, pulse, x_names, u, traj) =
    _SamplingMemberQTraj(sys, pulse, x_names, u, traj)

function SplineIntegrator(
    qtraj::Union{MultiKetTrajectory,_SamplingMemberQTraj},
    traj::NamedTrajectory;
    spline_order::Union{Int,Nothing} = nothing,
    alg::IntegrationAlgorithm = Tsit5Alg(),
    global_names::Union{Vector{Symbol},Nothing} = nothing,
    ket_sensitivity::Union{Bool,Nothing} = nothing,
    exact_hessian::Bool = false,
    second_order_shape::Symbol = :auto,
)
    # ── #344: which SHAPE the second-order path takes ─────────────────────────
    # `:directional` — forward-over-adjoint HVP, augmented state LINEAR in the
    #   parameter count (2d forward + 4d+n_p backward per ket), one small pair of
    #   solves per (knot, ket) per CG product. Consumed by `eval_constraint_hvp!`.
    # `:pair_indexed` — the legacy `build_second_order_sensitivity_ode` state
    #   `[Φ; S₁..Sₙ; T_active]`, QUADRATIC in the parameter count (`n_active` is the
    #   upper triangle of parameter pairs). This is the shape `eval_hessian_of_lagrangian`
    #   assembles a DENSE exact Hessian from, and the parity oracle for the directional
    #   path. Kept, not deleted.
    # `:auto` — `:directional` where it exists (Tsit5Alg), else `:pair_indexed`.
    #
    # DECISION (#344, parent #335 criterion 20): the shared
    # `build_second_order_sensitivity_ode` is FORKED, not changed — this cell stops
    # calling it under `:directional` and calls the ket-level `build_hvp_forward_ode`
    # / `build_second_order_adjoint_ode` pair instead. Ket's `exact_hessian` path is
    # therefore byte-identical (it still calls the shared builder, untouched), so
    # its parity holds BY CONSTRUCTION rather than by re-assertion.
    second_order_shape ∈ (:auto, :directional, :pair_indexed) || throw(
        ArgumentError(
            "second_order_shape must be :auto, :directional or :pair_indexed; got $second_order_shape",
        ),
    )
    if second_order_shape === :directional && exact_hessian && !(alg isa Tsit5Alg)
        error(
            "second_order_shape = :directional is implemented for Tsit5Alg only " *
            "(#344): the forward-over-adjoint pair is a Tsit5 augmented-ODE cell. " *
            "Use :pair_indexed (or :auto) with $(typeof(alg)).",
        )
    end
    resolved_second_order_shape =
        second_order_shape === :auto ?
        ((exact_hessian && alg isa Tsit5Alg) ? :directional : :pair_indexed) :
        second_order_shape

    if alg isa ChebyshevAlg && exact_hessian
        # ADR-0003 Decision 2: ChebyshevAlg owns its analytic matrix-free HVP; it must
        # not fall back to the Tsit5 second-order sensitivity ODE. The exact
        # Hessian-of-Lagrangian curvature is delivered matrix-free via
        # `eval_constraint_hvp!` (#136), not the dense `eval_hessian_of_lagrangian`.
        error(
            "exact_hessian=true is not supported with ChebyshevAlg on MultiKetTrajectory: " *
            "the analytic second-order curvature is the matrix-free eval_constraint_hvp! " *
            "(#136), and per ADR-0003 Decision 2 a non-Tsit5 algorithm must not fall back " *
            "to the Tsit5 sensitivity path.",
        )
    end
    if alg isa ChebyshevAlg && (spline_order == 1 || get_pulse(qtraj) isa LinearSplinePulse)
        # ADR-0003 / #136: the matrix-free ChebyshevAlg cores are cubic-only (the
        # linear packed-layout du fold-back is a scoped follow-on). Refuse LinearSpline
        # loudly at construction rather than erroring deep in the first VJP/HVP sweep.
        error(
            "ChebyshevAlg on MultiKetTrajectory is implemented for cubic splines only " *
            "(#136 scope, matching #133/#134's Chebyshev VJP/HVP); got a LinearSplinePulse. " *
            "Use a CubicSplinePulse (with bounded :du or an explicit bracket).",
        )
    end
    # Infer spline type from pulse if possible, otherwise use provided order or default
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

    sys = get_system(qtraj)
    x_names = state_names(qtraj)
    u = drive_name(qtraj)

    n_kets = length(x_names)

    # Auto-detect globals if not specified
    if isnothing(global_names)
        if !isempty(sys.global_params)
            global_names = collect(keys(sys.global_params))
        else
            global_names = Symbol[]
        end
    end

    @assert traj.N > 1 "Trajectory must have at least two timesteps."
    @assert all(name in traj.names for name in x_names) "All ensemble state names must be in trajectory"

    control_dim = traj.dims[u]

    # Determine global dimensions
    if !isempty(global_names)
        global_comps = vcat([traj.global_components[name] for name in global_names]...)
        global_dim = length(global_comps)
    else
        global_dim = 0
    end

    u_dim = control_dim + global_dim

    # Read drives from the system (operator-aware: structured ops pass through)
    H_drift_op = _to_operator(sys.H_drift)

    # Pre-compute G matrices for Magnus integrators (real isomorphic form)
    G_drift = Isomorphisms.G(_drift_matrix(sys.H_drift))
    G_drives = [Isomorphisms.G(drive_matrix(d)) for d in sys.H_drives]

    if isempty(sys.H_drives)
        throw(
            ArgumentError(
                "SplineIntegrator requires a QuantumSystem with explicit drive operators.\n" *
                "Function-based systems (QuantumSystem(H::Function, ...)) are not supported because\n" *
                "the spline integrator needs individual drive operators for analytical sensitivity equations.\n" *
                "Use: QuantumSystem(H_drift, H_drives, bounds) or\n" *
                "     QuantumSystem(H_drift, [LinearDrive(...), NonlinearDrive(...)], bounds)",
            ),
        )
    end

    # Operator wrappers for complex-domain forward and sensitivity ODE
    n_drives = sys.n_drives

    # Operators for complex-domain ODE (forward + sensitivity)
    H_drive_ops = [dynamics_operator(d) for d in sys.H_drives]

    # ============================================================================
    # Zero-allocation ODE right-hand side
    # Complex-domain: apply! with H operators (n×n instead of 2n×2n)
    # ============================================================================

    f! = if order == 1
        (dx, x, p, τ) -> begin
            uₖ = @view p[1:u_dim]
            uₖ₊₁ = @view p[(u_dim+1):2u_dim]
            Δtₖ = p[end-1]

            apply!(dx, H_drift_op, x, -im * Δtₖ, false)
            @inbounds for i = 1:n_drives
                u_i = (1 - τ) * uₖ[i] + τ * uₖ₊₁[i]
                apply!(dx, H_drive_ops[i], x, -im * Δtₖ * u_i, true)
            end
            return nothing
        end

    else  # order == 3
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

            apply!(dx, H_drift_op, x, -im * Δtₖ, false)
            @inbounds for i = 1:n_drives
                u_i = h00 * uₖ[i] + h10 * duₖ[i] + h01 * uₖ₊₁[i] + h11 * duₖ₊₁[i]
                apply!(dx, H_drive_ops[i], x, -im * Δtₖ * u_i, true)
            end
            return nothing
        end

    end

    ketdim = sys.levels
    tol = alg.tol

    # Dimensions for API - constraint per ket per knot point
    single_ket_dim = 2 * ketdim
    x_dim = n_kets * single_ket_dim
    N = traj.N
    dim = x_dim * (N - 1)

    if order == 1
        p_dim = 2 * u_dim
    elseif order == 3
        p_dim = 4 * u_dim
    end

    # Initialize ODE parameter vector
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

    # Algorithm-specific data
    if alg isa Tsit5Alg
        # Build forward ODE factory from system drives
        make_f! = build_forward_ode(H_drift_op, sys.H_drives, u_dim, order)
        x₀ = Matrix{ComplexF64}(I, ketdim, ketdim)
        Φ_probs =
            [ODEProblem(make_f!(), x₀, (0.0, 1.0), [p₀; Δt₀; t₀]) for _ = 1:(traj.N-1)]

        # Dense pattern for both adaptive and fixed-step modes. The fixed-step
        # branch used to build the pattern by probing Φ_probs[1], but the
        # complex-valued probe never matched Tsit5Data's Float64-typed
        # Φ_structure (construction crashed — #180), and a single-point probe
        # at the initial parameters is not a sound structural pattern anyway
        # (entries that vanish there may be nonzero elsewhere in control space).
        Φ_structure = sparse(ones(ketdim, ketdim))

        # #344 directional second-order path: per-(knot, ket) forward tangent
        # `[ψ; δψ] ∈ ℂ^{2d}` + backward second-order adjoint
        # `[ψ'; δψ'; λ; ν; g] ∈ ℂ^{4d + n_p}`. Both are LINEAR in n_p and carry NO
        # per-parameter-pair block, which is the whole point (#335 criterion 20).
        # The SAME ket-level builders the Ket Tsit5 cell uses — MultiKet's K kets
        # share one propagator, so μᵀc is additive over kets and each ket's
        # contribution is exactly a Ket-shaped forward-over-adjoint sweep.
        if exact_hessian && resolved_second_order_shape === :directional
            make_f!_hvp_fwd, n_params_hvp_fwd =
                build_hvp_forward_ode(H_drift_op, sys.H_drives, u_dim, ketdim, order)
            make_f!_2nd_adj, n_params_2nd_adj = build_second_order_adjoint_ode(
                H_drift_op,
                sys.H_drives,
                u_dim,
                ketdim,
                order,
            )
            hvp_fwd_probs = if isnothing(make_f!_hvp_fwd)
                nothing
            else
                [
                    ODEProblem(
                        make_f!_hvp_fwd(),
                        zeros(ComplexF64, 2 * ketdim),
                        (0.0, 1.0),
                        [p₀; Δt₀; t₀; zeros(n_params_hvp_fwd)],
                    ) for _ = 1:(traj.N-1)
                ]
            end
            hvp_bwd_probs = if isnothing(make_f!_2nd_adj)
                nothing
            else
                [
                    ODEProblem(
                        make_f!_2nd_adj(),
                        zeros(ComplexF64, 4 * ketdim + n_params_2nd_adj),
                        (1.0, 0.0),
                        [p₀; Δt₀; t₀; zeros(n_params_2nd_adj)],
                    ) for _ = 1:(traj.N-1)
                ]
            end
            alg_data = Tsit5Data{ComplexF64}(
                Φ_probs,
                Φ_structure,
                nothing,
                nothing,
                hvp_fwd_probs,
                hvp_bwd_probs,
            )
        else
            alg_data = Tsit5Data{ComplexF64}(Φ_probs, Φ_structure)
        end
    elseif alg isa ChebyshevAlg
        # #136: matrix-free exp-action cores, one d×K matrix state pushed through the
        # SAME per-interval cores/checkpoints/VJP-HVP scratch as the whole ensemble.
        # `n_cols = n_kets` sizes the state-shaped buffers to d×K (the packed-parameter
        # buffers stay K-independent; the Frobenius `dot` inside the cores auto-sums the
        # parameter-gradient contraction over kets). The n_sub sizing probe inside
        # `_build_alg_data` is a single-ket vector (spectral/bounds-based, K-independent).
        alg_data = _build_alg_data(
            alg,
            sys,
            G_drift,
            G_drives,
            H_drift_op,
            u_dim,
            control_dim,
            ketdim,
            order,
            traj,
            u,
            global_dim,
            global_names,
            p₀,
            Δt₀,
            t₀,
            tol;
            n_cols = n_kets,
        )
    else
        # Magnus path: use shared _build_alg_data dispatch
        alg_data = _build_alg_data(
            alg,
            sys,
            G_drift,
            G_drives,
            H_drift_op,
            u_dim,
            control_dim,
            ketdim,
            order,
            traj,
            u,
            global_dim,
            global_names,
            p₀,
            Δt₀,
            t₀,
            tol,
        )
    end

    ode_param_count = p_dim + 2

    prop_results =
        [PropagatorResult{ComplexF64}(ketdim, ode_param_count) for _ = 1:(traj.N-1)]

    # Decide whether to use ket-level sensitivity ODE.
    # Ket-level propagates K n-dim vectors instead of n×n matrices: O(KnJ) vs O(n²J).
    # Only beneficial when K < n.
    #
    # #136 / ADR-0003 Decision 2: ChebyshevAlg owns its (matrix-free Duhamel)
    # sensitivity — the Tsit5 augmented sensitivity ODE must NOT be built as a silent
    # fallback. Skip the sens-ODE construction entirely (mirrors the KetTrajectory
    # ChebyshevAlg guard); the matrix-free eval_constraint_{jvp,vjp,hvp}! path reads
    # neither sens_probs nor ket_sens_results.
    want_ket_sens = if alg isa ChebyshevAlg
        false
    elseif !isnothing(ket_sensitivity)
        ket_sensitivity
    else
        n_kets < ketdim  # auto: only when there's a speedup
    end

    if want_ket_sens
        make_f!_ket_sens, n_params = build_ket_sensitivity_ode(
            H_drift_op,
            sys.H_drives,
            u_dim,
            ketdim,
            order,
            n_kets,
        )
    else
        make_f!_ket_sens = nothing
        n_params = (order == 1 ? 2 : 4) * u_dim + 2
    end

    if !isnothing(make_f!_ket_sens)
        # Ket-level mode: store ket-level ODE problems in sens_probs
        sens_probs_ket, _ = build_ket_sensitivity_problems(
            make_f!_ket_sens,
            n_params,
            ketdim,
            n_kets,
            [p₀; Δt₀; t₀],
            traj.N - 1,
        )
        ket_sens_results =
            [zeros(ComplexF64, n_kets * ketdim, n_params) for _ = 1:(traj.N-1)]
        sens_probs = sens_probs_ket
        sens_state = nothing
        use_ket_sens = true
    elseif alg isa ChebyshevAlg
        # ChebyshevAlg (ADR-0003 Decision 2): no Tsit5 augmented sensitivity ODE.
        sens_probs = nothing
        sens_state = nothing
        ket_sens_results = nothing
        use_ket_sens = false
    else
        # Standard n×n sensitivity ODE (supports exact Hessian)
        make_f!_sens, _ =
            build_sensitivity_ode(H_drift_op, sys.H_drives, u_dim, ketdim, order)
        if !isnothing(make_f!_sens)
            sens_probs, sens_state = build_sensitivity_problems(
                make_f!_sens,
                n_params,
                ketdim,
                [p₀; Δt₀; t₀],
                traj.N,
            )
        else
            sens_probs = nothing
            sens_state = nothing
        end
        ket_sens_results = nothing
        use_ket_sens = false
    end

    # Build second-order sensitivity ODE for exact Hessian.
    # This uses the standard n×n ODE (same as KetTrajectory) — the second-order
    # sensitivity matrices T_{ij} are n×n and cannot be reduced to ket-level.
    # The per-ket contraction happens at assembly time in eval_hessian_of_lagrangian.
    #
    # #344: ONLY under `:pair_indexed`. Under `:directional` this quadratic state is
    # never built — `hess_active_pairs === nothing` with `exact_hessian == true` is
    # the marker that the directional fwd/bwd pair in `alg_data` is the second-order
    # path (see `_multiket_directional_hvp_probs`).
    if exact_hessian && resolved_second_order_shape === :pair_indexed
        make_f!_hess, hess_n_params, hess_active_pairs, hess_state_dim =
            build_second_order_sensitivity_ode(
                H_drift_op,
                sys.H_drives,
                u_dim,
                ketdim,
                order,
            )

        if !isnothing(make_f!_hess)
            hess_x₀ = zeros(ComplexF64, hess_state_dim)
            # Initialize Φ₀ = I, all sensitivities and Hessians = 0
            hess_x₀[1:(ketdim^2)] = vec(Matrix{ComplexF64}(I, ketdim, ketdim))
            hess_probs = [
                ODEProblem(make_f!_hess(), hess_x₀, (0.0, 1.0), [p₀; Δt₀; t₀]) for
                _ = 1:(traj.N-1)
            ]
            hess_state = zeros(ComplexF64, hess_state_dim)
        else
            hess_probs = nothing
            hess_state = nothing
            hess_active_pairs = nothing
            hess_n_params = 0
        end
    else
        hess_probs = nothing
        hess_state = nothing
        hess_active_pairs = nothing
        hess_n_params = 0
    end

    return SplineIntegrator{MultiKetTrajectory,S,ComplexF64,typeof(alg),typeof(alg_data)}(
        x_names,
        u,
        x_dim,
        u_dim,
        dim,
        tol,
        prop_results,
        ketdim,
        global_names,
        global_dim,
        sens_probs,
        sens_state,
        exact_hessian,
        hess_probs,
        hess_state,
        hess_active_pairs,
        hess_n_params,
        use_ket_sens,
        ket_sens_results,
        alg,
        alg_data,
    )
end

# ============================================================================ #
# Call operator for MultiKetTrajectory
# ============================================================================ #

@views function (𝒮::SplineIntegrator{MultiKetTrajectory})(
    δₖ::AbstractVector,
    zₖ::KnotPoint,
    zₖ₊₁::KnotPoint,
    k::Int,
    globals::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    # Build parameter vector using unified function
    pₖ = build_ode_params(𝒮, zₖ, zₖ₊₁, globals)

    if 𝒮.alg isa ChebyshevAlg
        # #136 matrix-free exp-action forward: pack the K per-ket iso states into a
        # d×K complex matrix, run the SAME per-interval core ONCE, scatter each column
        # back to its ket's constraint block. The propagator Φₖ is never materialized
        # (mirrors the KetTrajectory ChebyshevAlg call operator, one column-generic pass).
        n = 𝒮.ketdim
        K = length(𝒮.x_names)
        single_x = single_state_dim(𝒮)
        Ψ = Matrix{ComplexF64}(undef, n, K)
        for (i, x_name) in enumerate(𝒮.x_names)
            iso_to_complex_ket!(view(Ψ, :, i), zₖ[x_name])
        end
        F = _chebyshev_forward!(𝒮.alg_data::ChebyshevData, _spline_type(𝒮), k, pₖ, Ψ)
        for (i, x_name) in enumerate(𝒮.x_names)
            xₖ₊₁ = zₖ₊₁[x_name]
            i_slice = slice(i, single_x)
            δₖ[i_slice] = xₖ₊₁ - complex_ket_to_iso(F[:, i])
        end
        return nothing
    end

    # Algorithm-dispatched forward solve
    if 𝒮.alg isa Tsit5Alg
        data = 𝒮.alg_data::Tsit5Data
        Φₖ_prob = remake(data.Φ_probs[k], p = pₖ)
        Φₖ_raw = _solve_forward_tsit5(Φₖ_prob, 𝒮.alg, 𝒮.tol).u[end]
        Φₖ = propagator_to_iso(Φₖ_raw)
    else
        Φₖ_raw = _forward_propagate(𝒮, k, pₖ)
        Φₖ = eltype(Φₖ_raw) <: Complex ? propagator_to_iso(Φₖ_raw) : Φₖ_raw
    end

    # Cache complex Φ for ket-level Jacobian (avoids recomputing in compute_ode_jacobian!)
    if 𝒮.use_ket_sensitivity
        if eltype(Φₖ_raw) <: Complex
            copyto!(𝒮.prop_results[k].Φ_vec, vec(Φₖ_raw))
        else
            # Magnus path: convert real isomorphic 2n×2n → complex n×n
            n = 𝒮.ketdim
            Re_block = Φₖ_raw[1:n, 1:n]
            Im_block = Φₖ_raw[(n+1):2n, 1:n]
            Φ_complex = Re_block .+ im .* Im_block
            copyto!(𝒮.prop_results[k].Φ_vec, vec(Φ_complex))
        end
    end

    # Apply propagator to all kets
    for (i, x_name) in enumerate(𝒮.x_names)
        xₖ = zₖ[x_name]
        xₖ₊₁ = zₖ₊₁[x_name]
        i_slice = slice(i, single_state_dim(𝒮))
        δₖ[i_slice] = xₖ₊₁ - Φₖ * xₖ
    end

    return nothing
end

# ============================================================================ #
# Dispatch helpers for MultiKetTrajectory
# ============================================================================ #

"""
    compute_ode_jacobian!(𝒮::SplineIntegrator{MultiKetTrajectory}, zₖ, zₖ₊₁, k, globals)

Compute ODE Jacobian for MultiKetTrajectory.

When `use_ket_sensitivity` is true, solves ket-level sensitivity ODE (K n-dim vectors
per parameter instead of n×n matrices). Φ is already cached in `prop_results[k]`
from `evaluate!`. The ket-level results are stored in `ket_sens_results[k]`.

Otherwise, falls back to the standard n×n augmented sensitivity ODE.
"""
@views function compute_ode_jacobian!(
    𝒮::SplineIntegrator{MultiKetTrajectory},
    zₖ::KnotPoint,
    zₖ₊₁::KnotPoint,
    k::Int,
    globals::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    pₖ = build_ode_params(𝒮, zₖ, zₖ₊₁, globals)

    if 𝒮.use_ket_sensitivity
        # Ket-level sensitivity: Φ is cached from evaluate!, solve for per-ket sensitivities
        _compute_ket_sensitivity!(𝒮, zₖ, pₖ, k)
    else
        isnothing(𝒮.sens_probs) && error(
            "SplineIntegrator requires explicit H_drift/H_drives matrices for Jacobian computation. " *
            "Function-based systems are not supported.",
        )
        _compute_ode_jacobian_analytical!(𝒮, pₖ, k)
    end
    return nothing
end

"""
    _compute_ket_sensitivity!(𝒮, zₖ, pₖ, k)

Solve the ket-level sensitivity ODE for knot point k.

Sets up initial conditions [ψ₁(0), ..., ψ_K(0), 0, ..., 0] from the trajectory,
solves the ODE, and extracts ket-level sensitivity vectors sᵢⱼ into ket_sens_results[k].
"""
function _compute_ket_sensitivity!(
    𝒮::SplineIntegrator{MultiKetTrajectory},
    zₖ::KnotPoint,
    pₖ::AbstractVector,
    k::Int,
)
    n = 𝒮.ketdim
    K = length(𝒮.x_names)
    n_params = size(𝒮.ket_sens_results[k], 2)

    # Build initial condition: [ψ₁, ..., ψ_K, zeros...]
    state_dim = K * n * (1 + n_params)
    u0 = zeros(ComplexF64, state_dim)
    for (i, x_name) in enumerate(𝒮.x_names)
        ψ̃ = zₖ[x_name]
        off = (i - 1) * n
        # Convert iso to complex ket
        iso_to_complex_ket!(@view(u0[(off+1):(off+n)]), ψ̃)
    end

    # Solve the ket-level sensitivity ODE
    prob = remake(𝒮.sens_probs[k], u0 = u0, p = pₖ)
    sol = solve(
        prob,
        Tsit5();
        abstol = 𝒮.tol,
        reltol = 𝒮.tol,
        save_everystep = false,
        maxiters = 1_000_000,
    )
    final_state = sol.u[end]

    # Extract ket-level sensitivities into ket_sens_results[k]
    # Layout: result matrix is (K*n, n_params)
    # ODE state: [ψ₁..ψ_K, s_{1,1}..s_{1,J}, s_{2,1}..s_{2,J}, ..., s_{K,J}]
    ket_block = K * n
    result = 𝒮.ket_sens_results[k]
    @inbounds for i = 1:K
        for j = 1:n_params
            off_s = ket_block + ((i - 1) * n_params + (j - 1)) * n
            row_off = (i - 1) * n
            for r = 1:n
                result[row_off+r, j] = final_state[off_s+r]
            end
        end
    end

    return nothing
end

"""
    get_state_vectors(𝒮::SplineIntegrator{MultiKetTrajectory}, zₖ)

Extract all state vectors from a knot point for MultiKetTrajectory.
Returns a vector of state vectors, one per ket.
"""
function get_state_vectors(𝒮::SplineIntegrator{MultiKetTrajectory}, zₖ::KnotPoint)
    return [zₖ[x_name] for x_name in 𝒮.x_names]
end

# ============================================================================ #
# The matrix-free interface properties `SplineIntegrator{MultiKetTrajectory}`
# declares (#338 AC3).
#
# CORE WIDTH `K` IN `K` COMPONENTS. `K` kets, each its own `2·ketdim` trajectory
# component, so the `d×K` core width and the component count COINCIDE here — which
# is exactly why `GPUJacobianOp`'s `K = length(integ.x_names)` looked correct: it
# is, for this one cell. Unitary is `d`-wide in ONE component, so the two are
# declared separately and neither is derived from the other.
#
# MultiKet's FIRST-ORDER bodies are the reference the inner-kernel contract was
# extracted FROM (`traj.datavec` + precomputed columns, chunk-scoped scratch, the
# Frobenius parameter-block reassociation) and are deliberately NOT retrofitted
# onto the seam by #338 — they are already knot-flat. Its second-order path is
# slice #344.
# ============================================================================ #

block_shape(𝒮::SplineIntegrator{<:MultiKetTrajectory}) =
    MatrixFreeBlockShape(length(𝒮.x_names), length(𝒮.x_names))

generator_action(𝒮::SplineIntegrator{<:MultiKetTrajectory}) = OneSidedAction(𝒮.alg_data)

iso_packing(::SplineIntegrator{<:MultiKetTrajectory}) = BlockedRealImagPacking()

matrix_free_driver(::SplineIntegrator{<:MultiKetTrajectory}) = HostTaskDriver()

# ============================================================================ #
# Jacobian structure for MultiKetTrajectory
# ============================================================================ #

function jacobian_structure(
    ::Type{MultiKetTrajectory},
    x_names::Vector{Symbol},
    u_name::Symbol,
    ketdim::Int,
    Φ_structure::SparseMatrixCSC,
    spline_order::Int,
    traj::NamedTrajectory;
    global_names::Vector{Symbol} = Symbol[],
)
    n_kets = length(x_names)
    x_dim = 2 * ketdim
    total_x_dim = n_kets * x_dim
    u_dim = traj.dims[u_name]
    z_dim = traj.dim
    u_comps = traj.components[u_name]
    Δt_comp = traj.components[traj.timestep][1]
    t_comp = traj.components[:t][1]

    # Determine global dimensions
    if traj.global_dim > 0
        global_dim = traj.global_dim
    elseif !isempty(global_names)
        global_dim = length(global_names)
    else
        global_dim = 0
    end

    # Jacobian columns: [z_k, z_{k+1}, globals]
    ∂𝒮 = spzeros(total_x_dim, 2 * z_dim + global_dim)

    # Build parameter component list
    # Only include du for cubic splines - linear spline dynamics don't depend on du
    if spline_order == 3
        du_name_in_traj = haskey(traj.components, :du) ? :du : :dθ
        if haskey(traj.components, du_name_in_traj)
            du_comps = traj.components[du_name_in_traj]
            pₖ_comps = [u_comps; z_dim .+ u_comps; du_comps; z_dim .+ du_comps]
        else
            pₖ_comps = [u_comps; z_dim .+ u_comps]
        end
    else
        pₖ_comps = [u_comps; z_dim .+ u_comps]
    end
    p_comps = [pₖ_comps; Δt_comp; t_comp]

    # Structure for each ket
    for (i, x_name) in enumerate(x_names)
        x_comps = traj.components[x_name]
        i_slice = slice(i, x_dim)

        # ∂xₖ𝒮: Propagator structure for this ket
        ∂𝒮[i_slice, x_comps] = Φ_structure

        # ∂xₖ₊₁𝒮: Identity for this ket
        ∂𝒮[i_slice, z_dim .+ x_comps] = I(x_dim)

        # ∂pₖ𝒮: Dense columns for parameters
        for j ∈ p_comps
            ∂𝒮[i_slice, j] .= 1.0
        end

        # ∂g𝒮: Global variable dependencies
        if global_dim > 0
            global_offset = 2 * z_dim
            for j = 1:global_dim
                ∂𝒮[i_slice, global_offset+j] .= 1.0
            end
        end
    end

    return ∂𝒮
end

# ============================================================================ #
# API Methods for MultiKetTrajectory
# ============================================================================ #

@views function evaluate!(
    δ::AbstractVector,
    𝒮::SplineIntegrator{MultiKetTrajectory},
    traj::NamedTrajectory,
)
    # Extract globals once
    globals = extract_globals(𝒮, traj)

    Threads.@threads for k = 1:(traj.N-1)
        δₖ = δ[slice(k, 𝒮.x_dim)]
        𝒮(δₖ, traj[k], traj[k+1], k, globals)
    end
    return nothing
end

@views function eval_jacobian(
    𝒮::SplineIntegrator{MultiKetTrajectory},
    traj::NamedTrajectory,
)
    # Tier 4d (2026-05-21): triplet-pattern sparse assembly.
    # The original `spzeros(F_dim, Z_dim) + setindex!` assembly was O(nnz)
    # per scalar write because SparseMatrixCSC's setindex! rebuilds
    # colptr/rowval. For i=2 N=50 K=8, the assembly alone took 14h on CPU
    # (we observed this directly). The triplet construction mirrors the
    # GPU Tier 1 `_assemble_jacobian_sparse` pattern in PiccolissimoCUDAExt
    # and runs in O(nnz) total via a single linear pass + one `sparse(I,J,V)`.
    N = traj.N
    # Layout-derived per-knot index tables, cached across matvecs (#307):
    # rebuilding them per call cost O(N) allocations on a routine that runs
    # once per inner CG iteration.
    layout = matrix_free_layout(𝒮, traj)
    z_dim = traj.dim
    global_dim = traj.global_dim
    F_dim = 𝒮.x_dim * (N - 1)
    Z_dim = z_dim * N + global_dim
    single_x_dim = single_state_dim(𝒮)
    K = length(𝒮.x_names)
    pdim = propagator_side_dim(𝒮)

    # Extract globals once
    globals = extract_globals(𝒮, traj)

    # Compute Jacobians in parallel
    Threads.@threads for k = 1:(N-1)
        compute_ode_jacobian!(𝒮, traj[k], traj[k+1], k, globals)
    end

    # Upper-bound nnz: per interval k, per ket i, we contribute
    #   • dense (2n × 2n) `neg_Φ_iso` block → 4n²
    #   • diagonal I(2n) at xₖ₊₁ columns → 2n
    #   • dense (2n × n_traj_cols) sensitivity block at param columns → 2n·n_traj_cols
    # `get_param_indices` returns the same length for every k. Sample interval 1.
    sample_traj_cols, _ = get_param_indices(𝒮, traj, 1)
    n_traj_cols = length(sample_traj_cols)
    nnz_per_interval = K * (4 * pdim^2 + 2 * pdim + 2 * pdim * n_traj_cols)
    nnz_upper = (N - 1) * nnz_per_interval

    I_buf = Vector{Int}(undef, nnz_upper)
    J_buf = Vector{Int}(undef, nnz_upper)
    V_buf = Vector{Float64}(undef, nnz_upper)
    idx = 0

    temp_vec = Vector{Float64}(undef, single_x_dim)
    tmp_complex = Vector{ComplexF64}(undef, pdim)
    ψ_complex = Vector{ComplexF64}(undef, pdim)

    if 𝒮.use_ket_sensitivity
        # --- Ket-level assembly: use cached Φ and pre-computed sᵢⱼ vectors ---
        @inbounds for k = 1:(N-1)
            F_rows = layout.F_rows[k]
            zₖ = traj[k]
            Φₖ = get_propagator(𝒮.prop_results[k], pdim)
            neg_Φₖ_iso = -propagator_to_iso(Φₖ)
            traj_cols, ode_indices = layout.traj_cols[k], layout.ode_indices[k]
            zk_offset = (k - 1) * z_dim

            for (i, x_name) in enumerate(𝒮.x_names)
                x_comps = layout.x_comps[i]
                i_slice = slice(i, single_x_dim)
                rows = F_rows[i_slice]

                # ∂xₖ𝒮 block: -Φ_iso (single_x_dim × single_x_dim, dense)
                for j = 1:single_x_dim
                    col = zk_offset + x_comps[j]
                    for r_idx = 1:single_x_dim
                        idx += 1
                        I_buf[idx] = rows[r_idx]
                        J_buf[idx] = col
                        V_buf[idx] = neg_Φₖ_iso[r_idx, j]
                    end
                end

                # ∂xₖ₊₁𝒮: identity block
                for r_idx = 1:single_x_dim
                    idx += 1
                    I_buf[idx] = rows[r_idx]
                    J_buf[idx] = zk_offset + z_dim + x_comps[r_idx]
                    V_buf[idx] = 1.0
                end

                # ∂pₖ𝒮: ket-level sensitivity vectors
                row_off = (i - 1) * pdim
                for (col, ode_idx) in zip(traj_cols, ode_indices)
                    sᵢⱼ = 𝒮.ket_sens_results[k][(row_off+1):(row_off+pdim), ode_idx]
                    @inbounds for r = 1:pdim
                        temp_vec[r] = -real(sᵢⱼ[r])
                        temp_vec[pdim+r] = -imag(sᵢⱼ[r])
                    end
                    for r_idx = 1:single_x_dim
                        idx += 1
                        I_buf[idx] = rows[r_idx]
                        J_buf[idx] = col
                        V_buf[idx] = temp_vec[r_idx]
                    end
                end
            end
        end
    else
        # --- Standard assembly: use n×n sensitivity matrices ---
        @inbounds for k = 1:(N-1)
            F_rows = layout.F_rows[k]
            zₖ = traj[k]
            Φₖ = get_propagator(𝒮.prop_results[k], pdim)
            ∂ₚΦₖ = get_sensitivities(𝒮.prop_results[k], pdim)
            neg_Φₖ_iso = -propagator_to_iso(Φₖ)
            traj_cols, ode_indices = layout.traj_cols[k], layout.ode_indices[k]
            zk_offset = (k - 1) * z_dim

            for (i, x_name) in enumerate(𝒮.x_names)
                x_comps = layout.x_comps[i]
                xₖ = zₖ[x_name]
                i_slice = slice(i, single_x_dim)
                rows = F_rows[i_slice]

                # ∂xₖ𝒮 block: -Φ_iso
                for j = 1:single_x_dim
                    col = zk_offset + x_comps[j]
                    for r_idx = 1:single_x_dim
                        idx += 1
                        I_buf[idx] = rows[r_idx]
                        J_buf[idx] = col
                        V_buf[idx] = neg_Φₖ_iso[r_idx, j]
                    end
                end

                # ∂xₖ₊₁𝒮: identity block
                for r_idx = 1:single_x_dim
                    idx += 1
                    I_buf[idx] = rows[r_idx]
                    J_buf[idx] = zk_offset + z_dim + x_comps[r_idx]
                    V_buf[idx] = 1.0
                end

                # ∂pₖ𝒮 via sensitivity matvec
                iso_to_complex_ket!(ψ_complex, xₖ)
                for (col, ode_idx) in zip(traj_cols, ode_indices)
                    Sⱼ = ∂ₚΦₖ[:, :, ode_idx]
                    sensitivity_to_jac_col!(temp_vec, Sⱼ, ψ_complex, tmp_complex)
                    for r_idx = 1:single_x_dim
                        idx += 1
                        I_buf[idx] = rows[r_idx]
                        J_buf[idx] = col
                        V_buf[idx] = temp_vec[r_idx]
                    end
                end
            end
        end
    end

    # Trim and assemble. The upper bound is tight (no `setindex!` collisions),
    # but defend against the rare case where it isn't.
    if idx < nnz_upper
        resize!(I_buf, idx)
        resize!(J_buf, idx)
        resize!(V_buf, idx)
    end
    # `sparse(I, J, V, m, n)` sums duplicates by default — there are no
    # structural duplicates in this assembly so duplicate-sum is a no-op
    # and the output values are bit-identical to the original setindex!
    # path modulo floating-point ordering of identical sums.
    return sparse(I_buf, J_buf, V_buf, F_dim, Z_dim)
end

# ─────────────────────────────────────────────────────────────────────────────
# Matrix-free constraint Jacobian products for MultiKetTrajectory.
# These mirror `eval_jacobian` block-for-block but contract with the input
# vector instead of assembling a sparse `J`. The three per-(knot, ket) blocks
# are: ∂xₖ = -Φ_iso, ∂xₖ₊₁ = I, ∂pₖ = sensitivity columns.
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────── #
# Per-iterate first-order sensitivity cache for the generic MultiKet cell (#310)
#
# `eval_constraint_jvp!` / `eval_constraint_vjp!` re-ran the per-knot sensitivity
# ODE on EVERY invocation. Under Newton-CG that is once per inner CG product
# rather than once per accepted iterate, which is why naively routing this cell
# matrix-free would be SLOWER than the sparse path — the sparse path already skips
# reassembly across inner iterations via the x-version stamp.
#
# The sensitivities are constant for a fixed decision vector, so caching them per
# accepted iterate is exact, not an approximation — the same assumption the
# existing Hessian-vector cache and the sparse Jacobian's reuse path already make.
# ─────────────────────────────────────────────────────────────────────────── #

const _MULTIKET_SENS_REFRESHES = Ref(0)

"""
    multiket_sens_refresh_count() -> Int

Number of per-knot sensitivity REFRESH PASSES run since load: one per
`build_multiket_sens_cache` call, plus one per uncached `eval_constraint_jvp!` /
`eval_constraint_vjp!`. The witness for #310's "two applies at the same decision
vector solve once, not twice" — counted per pass rather than per knot so the
`compute_ode_jacobian!` hot path carries no counter.
"""
multiket_sens_refresh_count() = _MULTIKET_SENS_REFRESHES[]

# One refresh pass over all knots. Bumps the witness.
@views function eval_hessian_of_lagrangian(
    𝒮::SplineIntegrator{MultiKetTrajectory},
    traj::NamedTrajectory,
    μ::AbstractVector,
)
    # #344: `exact_hessian` with NO pair-indexed `hess_active_pairs` means this cell
    # carries the DIRECTIONAL second-order shape, whose curvature is only reachable
    # matrix-free. Silently dropping to Gauss-Newton here is the failure mode that
    # HIDES — every parity test still passes while the (p,p) curvature is zero — so
    # this errors instead of degrading.
    if !isnothing(_multiket_directional_hvp_probs(𝒮))
        error(
            "eval_hessian_of_lagrangian cannot assemble a dense exact Hessian for a " *
            "MultiKet cell with the :directional second-order shape (#344): the " *
            "pair-indexed second-order sensitivities were never built. Use the " *
            "matrix-free eval_constraint_hvp! (which IS the exact second-order " *
            "curvature), or construct with second_order_shape = :pair_indexed to get " *
            "the dense assembly back. Returning a Gauss-Newton Hessian here would " *
            "silently contribute ZERO second-order curvature.",
        )
    end

    if 𝒮.use_ket_sensitivity && !𝒮.exact_hessian
        # Ket-level mode without exact Hessian: GN cross-terms require full n×n sensitivities.
        # Return zero Hessian — caller should use L-BFGS (eval_hessian=false).
        N = traj.N
        z_dim = traj.dim
        global_dim = traj.global_dim
        Z_dim = z_dim * N + global_dim
        return spzeros(Z_dim, Z_dim)
    end

    N = traj.N
    # Layout-derived per-knot index tables, cached across calls (#307): rebuilding
    # them per call cost O(N) allocations.
    layout = matrix_free_layout(𝒮, traj)
    z_dim = traj.dim
    global_dim = traj.global_dim
    Z_dim = z_dim * N + global_dim
    single_x_dim = single_state_dim(𝒮)

    # Build full Hessian
    μ∂²F = spzeros(Z_dim, Z_dim)
    temp_vec = Vector{Float64}(undef, single_x_dim)

    pdim = propagator_side_dim(𝒮)
    tmp1 = Vector{ComplexF64}(undef, pdim)
    tmp2 = Vector{ComplexF64}(undef, pdim)

    @inbounds for k = 1:(N-1)
        μₖ_full = μ[slice(k, 𝒮.x_dim)]
        zₖ = traj[k]
        zk_offset = (k - 1) * z_dim

        # Get parameter mapping for this knot point
        traj_cols, ode_indices = layout.traj_cols[k], layout.ode_indices[k]
        # Lazy zip: `collect` here allocated a fresh pair vector per knot (#307).
        param_pairs = zip(ode_indices, traj_cols)

        if 𝒮.exact_hessian && !isnothing(𝒮.hess_active_pairs)
            # Exact Hessian: solve second-order ODE to get Φ, Sⱼ, T_{ij}
            # (same n×n ODE as KetTrajectory — shared infrastructure)
            T_mats = _compute_ode_hessian_from_traj(𝒮, traj, k)

            # Extract ∂ₚΦₖ from prop_results (populated by _compute_ode_hessian_analytical!)
            ∂ₚΦₖ = get_sensitivities(𝒮.prop_results[k], pdim)

            # Reverse lookup ODE param index → trajectory column, read from the
            # layout table as a preallocated array (#307): a per-knot `Dict` here
            # allocated once per knot and hashed once per (ket, parameter-pair) in
            # the loop below. `0` marks a parameter with no column.
            ode_to_col = layout.ode_to_col[k]

            # Per-ket assembly: (x,p) cross-terms and (p,p) blocks
            for (i, x_name) in enumerate(𝒮.x_names)
                x_comps = zₖ.components[x_name]
                i_slice = slice(i, single_x_dim)
                μₖ = μₖ_full[i_slice]
                x_full_cols = zk_offset .+ x_comps

                # (x,p) cross-terms: -(Sⱼ)' * μ_C for each ket
                for (ode_idx, traj_col) in param_pairs
                    Sⱼ = ∂ₚΦₖ[:, :, ode_idx]
                    sensitivity_to_hess_col!(temp_vec, Sⱼ, μₖ, tmp1, tmp2)
                    for (j, xi) in enumerate(x_full_cols)
                        if xi <= traj_col
                            μ∂²F[xi, traj_col] += temp_vec[j]
                        else
                            μ∂²F[traj_col, xi] += temp_vec[j]
                        end
                    end
                end

                # (p,p) blocks: -real(μ_C' * T_{ij} * ψ_C) for each ket
                ψₖ = zₖ[x_name]  # real iso ket for this ensemble member
                for (pair_k, (ode_i, ode_j)) in enumerate(𝒮.hess_active_pairs)
                    (ode_i <= length(ode_to_col) && ode_j <= length(ode_to_col)) || continue
                    traj_i = ode_to_col[ode_i]
                    traj_j = ode_to_col[ode_j]
                    (traj_i == 0 || traj_j == 0) && continue

                    val = hessian_pp_contraction(T_mats[pair_k], μₖ, ψₖ, tmp1, tmp2)

                    ri, rj = minmax(traj_i, traj_j)
                    μ∂²F[ri, rj] += val
                end
            end
        else
            # Gauss-Newton: reuse prop_results from eval_jacobian (cross-terms only)
            ∂ₚΦₖ = get_sensitivities(𝒮.prop_results[k], pdim)

            for (i, x_name) in enumerate(𝒮.x_names)
                x_comps = zₖ.components[x_name]
                i_slice = slice(i, single_x_dim)
                μₖ = μₖ_full[i_slice]
                x_full_cols = zk_offset .+ x_comps

                for (ode_idx, traj_col) in param_pairs
                    Sⱼ = ∂ₚΦₖ[:, :, ode_idx]
                    sensitivity_to_hess_col!(temp_vec, Sⱼ, μₖ, tmp1, tmp2)
                    for (j, xi) in enumerate(x_full_cols)
                        if xi <= traj_col
                            μ∂²F[xi, traj_col] += temp_vec[j]
                        else
                            μ∂²F[traj_col, xi] += temp_vec[j]
                        end
                    end
                end
            end
        end
    end

    return μ∂²F
end

"""
    build_multiket_hvp_cache(𝒮::SplineIntegrator{MultiKetTrajectory}, traj)

Solve the second-order sensitivity ODE (SOSE) for every knot ONCE and snapshot
the per-knot `T_mats` (second-order) and `∂ₚΦ` (first-order) sensitivities.

The SOSE depends only on `x` (the trajectory), so it is constant across the CG
inner loop of a single Newton solve. Without this cache, `eval_constraint_hvp!`
re-solves the SOSE on every HVP call (~30-60× per outer iter) — the dominant
cost. The cache also snapshots `∂ₚΦ` because the AL Lagrangian-HVP runs JVP/VJP
(which overwrite `prop_results`) *before* the constraint-HVP within the same
call, so reading `prop_results` directly would be stale.
"""
