# ============================================================================ #
# Constructor for KetTrajectory
# ============================================================================ #

function SplineIntegrator(
    qtraj::KetTrajectory,
    N::Int;
    spline_order::Union{Int,Nothing} = nothing,
    alg::IntegrationAlgorithm = Tsit5Alg(),
    global_names::Union{Vector{Symbol},Nothing} = nothing,
    exact_hessian::Bool = false,
    use_ket_sensitivity::Bool = false,
)
    if use_ket_sensitivity && exact_hessian
        error(
            "use_ket_sensitivity=true is incompatible with exact_hessian=true: " *
            "the second-order sensitivity ODE requires the n×n propagator-level " *
            "sensitivity path (T_{ij} matrices cannot be reduced to ket-level).",
        )
    end
    if use_ket_sensitivity && alg isa Rodas5PAlg
        # The matrix-state forward ODE used by use_ket_sensitivity (x₀ = I to
        # materialize Φ for the -Φ Jacobian block) needs a different closed-form
        # Jacobian shape than the build_forward_ode_jacobian helper provides
        # (which targets vector ket states). Disallow the combination for Phase 0;
        # extending the helper to matrix-state ODEs is out of scope here.
        error(
            "use_ket_sensitivity=true with Rodas5PAlg is not yet supported. " *
            "Use Tsit5Alg (the default) for the ket-level sensitivity path.",
        )
    end
    if alg isa ChebyshevAlg && use_ket_sensitivity
        # ADR-0003 Decision 2 transition gate: ChebyshevAlg owns its sensitivity
        # method, and it has not landed yet (#133) — requesting one errors at
        # construction rather than silently differentiating a Tsit5 ODE.
        error(
            "use_ket_sensitivity=true is not supported with ChebyshevAlg: its analytic " *
            "sensitivity method has not landed yet (#133), and per ADR-0003 Decision 2 " *
            "a non-Tsit5 algorithm must not fall back to the Tsit5 sensitivity path.",
        )
    end
    if alg isa ChebyshevAlg && exact_hessian
        error(
            "exact_hessian=true is not supported with ChebyshevAlg: its analytic " *
            "sensitivity method has not landed yet (#133), and per ADR-0003 Decision 2 " *
            "a non-Tsit5 algorithm must not fall back to the Tsit5 sensitivity path.",
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
        # Non-spline pulse: use provided order or default to 1
        order = isnothing(spline_order) ? 1 : spline_order
        S = order == 1 ? LinearSpline : CubicSpline
    end

    sys = get_system(qtraj)
    traj = NamedTrajectory(qtraj, N)
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

    return _spline_ket(
        sys,
        pulse,
        x,
        u,
        traj;
        spline_order = spline_order,
        alg = alg,
        global_names = global_names,
        exact_hessian = exact_hessian,
        use_ket_sensitivity = use_ket_sensitivity,
    )
end

# Per-member inner constructor (#408's sampling lane): all pieces pre-materialized
# (sys/pulse/x/u from the member's base trajectory; traj the already-expanded
# member trajectory). The (qtraj, N) outer method is this plus trajectory
# assembly. S/order inference repeated here — the outer method's validation
# (Rodas5P/Chebyshev guards) has already run by the time this is reached.
function _spline_ket(
    sys,
    pulse,
    x::Symbol,
    u::Symbol,
    traj::NamedTrajectory;
    spline_order::Union{Int,Nothing} = nothing,
    alg::IntegrationAlgorithm = Tsit5Alg(),
    global_names::Union{Vector{Symbol},Nothing} = nothing,
    exact_hessian::Bool = false,
    use_ket_sensitivity::Bool = false,
)
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

    if isnothing(global_names)
        if !isempty(sys.global_params)
            global_names = collect(keys(sys.global_params))
        else
            global_names = Symbol[]
        end
    end

    @assert traj.N > 1 "Trajectory must have at least two timesteps."

    control_dim = traj.dims[u]

    # Determine global dimensions
    if !isempty(global_names)
        global_comps = vcat([traj.global_components[name] for name in global_names]...)
        global_dim = length(global_comps)
    else
        global_dim = 0
    end

    # u_dim includes BOTH controls and globals for unified ODE parameter handling
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
    # Complex-domain: apply! with H operators (n-vec instead of 2n-vec)
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

    # Dimensions for API
    x_dim = traj.dims[x]
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
    if alg isa Tsit5Alg || alg isa Rodas5PAlg
        # Build forward ODE factory from system drives
        make_f! = build_forward_ode(H_drift_op, sys.H_drives, u_dim, order)

        # Initial condition for the forward ODE:
        #   - Default (use_ket_sensitivity=false): x₀ = ψ_init (a ket vector). The
        #     forward path solves dψ/dτ = -iΔt H ψ, returning ψ_{k+1} = Φ_k ψ_k
        #     directly (no propagator ever materialized).
        #   - use_ket_sensitivity=true: x₀ = I (the identity). The forward path
        #     solves dΦ/dτ = -iΔt H Φ, returning Φ_k itself. Required because the
        #     ket-level sensitivity Jacobian assembly needs the full propagator
        #     for the ∂xₖ𝒮 = -Φ_k block.
        x₀ =
            use_ket_sensitivity ? Matrix{ComplexF64}(I, ketdim, ketdim) :
            zeros(ComplexF64, ketdim)

        # For implicit solvers (Rodas5P) on complex states, OrdinaryDiffEq's default
        # Jacobian path goes through ForwardDiff which rejects ComplexF64 — so we
        # provide the closed-form Jacobian J = -iΔt H(τ;p) explicitly. Tsit5 (explicit)
        # never autodiffs the RHS, so it stays on the bare-f! path.
        Φ_probs = if alg isa Rodas5PAlg
            make_jac! = build_forward_ode_jacobian(
                H_drift_op,
                sys.H_drives,
                u_dim,
                ketdim,
                order,
            )
            [
                ODEProblem(
                    ODEFunction(make_f!(); jac = make_jac!()),
                    x₀,
                    (0.0, 1.0),
                    [p₀; Δt₀; t₀],
                ) for _ = 1:(traj.N-1)
            ]
        else
            [ODEProblem(make_f!(), x₀, (0.0, 1.0), [p₀; Δt₀; t₀]) for _ = 1:(traj.N-1)]
        end

        # Dense pattern for both adaptive and fixed-step modes. The fixed-step
        # Tsit5 branch used to build the pattern by probing Φ_probs[1], but the
        # complex-valued probe never matched the Float64-typed Φ_structure
        # field (construction crashed — #180), and a single-point probe at the
        # initial parameters is not a sound structural pattern anyway.
        Φ_structure = sparse(ones(ketdim, ketdim))

        alg_data = if alg isa Rodas5PAlg
            Rodas5PData{ComplexF64}(Φ_probs, Φ_structure)
        else
            # Phase 2: build per-knot ket-level JVP ODE problems (state [ψ; Y] ∈ ℂ^{2d}).
            # Used by the Tsit5Alg-specialized step_jvp! method to avoid materializing
            # propagator-level sensitivities for matrix-free JVP.
            make_f!_jvp, n_params_jvp =
                build_ket_jvp_ode(H_drift_op, sys.H_drives, u_dim, ketdim, order)
            jvp_probs = if isnothing(make_f!_jvp)
                nothing
            else
                jvp_x₀ = zeros(ComplexF64, 2 * ketdim)  # [ψ(0); Y(0)]
                jvp_p₀ = [p₀; Δt₀; t₀; zeros(n_params_jvp)]
                [ODEProblem(make_f!_jvp(), jvp_x₀, (0.0, 1.0), jvp_p₀) for _ = 1:(traj.N-1)]
            end

            # Phase 3: build per-knot ket-level adjoint VJP ODE problems
            # (state [ψ; λ; g] ∈ ℂ^{2d + n_p}, integrated backward via
            # tspan = (1.0, 0.0)). Used by the Tsit5Alg-specialized step_vjp!
            # method to avoid materializing propagator-level sensitivities for
            # matrix-free VJP.
            make_f!_adj, n_params_adj =
                build_adjoint_ode(H_drift_op, sys.H_drives, u_dim, ketdim, order)
            vjp_probs = if isnothing(make_f!_adj)
                nothing
            else
                # Standard params; no v_p tail (VJP doesn't contract over
                # parameters at the parameter-vector level).
                adj_x₀ = zeros(ComplexF64, 2 * ketdim + n_params_adj)
                adj_p₀ = [p₀; Δt₀; t₀]
                [ODEProblem(make_f!_adj(), adj_x₀, (1.0, 0.0), adj_p₀) for _ = 1:(traj.N-1)]
            end

            # Phase 4: build per-knot ket-level HVP ODE problems.
            # Forward path [ψ; δψ] ∈ ℂ^{2d} (structurally identical to JVP) +
            # backward second-order adjoint [ψ'; δψ'; λ; ν; g] ∈ ℂ^{4d + n_p}.
            # Used by the Tsit5Alg-specialized step_hvp! method.
            make_f!_hvp_fwd, n_params_hvp_fwd =
                build_hvp_forward_ode(H_drift_op, sys.H_drives, u_dim, ketdim, order)
            hvp_fwd_probs = if isnothing(make_f!_hvp_fwd)
                nothing
            else
                hvp_fwd_x₀ = zeros(ComplexF64, 2 * ketdim)
                hvp_fwd_p₀ = [p₀; Δt₀; t₀; zeros(n_params_hvp_fwd)]
                [
                    ODEProblem(make_f!_hvp_fwd(), hvp_fwd_x₀, (0.0, 1.0), hvp_fwd_p₀) for _ = 1:(traj.N-1)
                ]
            end

            make_f!_2nd_adj, n_params_2nd_adj = build_second_order_adjoint_ode(
                H_drift_op,
                sys.H_drives,
                u_dim,
                ketdim,
                order,
            )
            hvp_bwd_probs = if isnothing(make_f!_2nd_adj)
                nothing
            else
                # State [ψ'; δψ'; λ; ν; g]: 4d + n_p complex
                # Params [standard; v_p tail]: n_p + n_p
                hvp_bwd_x₀ = zeros(ComplexF64, 4 * ketdim + n_params_2nd_adj)
                hvp_bwd_p₀ = [p₀; Δt₀; t₀; zeros(n_params_2nd_adj)]
                [
                    ODEProblem(make_f!_2nd_adj(), hvp_bwd_x₀, (1.0, 0.0), hvp_bwd_p₀) for _ = 1:(traj.N-1)
                ]
            end

            Tsit5Data{ComplexF64}(
                Φ_probs,
                Φ_structure,
                jvp_probs,
                vjp_probs,
                hvp_fwd_probs,
                hvp_bwd_probs,
            )
        end
    elseif alg isa MagnusGL4Alg || alg isa MagnusAdapt4Alg
        # ADR-0003 Magnus (time-varying-H spline) cell (#223): the interval
        # generator Ω is built by the midpoint-Magnus product and exp(Ω) is applied
        # through the Chebyshev exp-action — matrix-free forward AND its consistent
        # analytic Duhamel VJP (replacing the Tsit5-augmented-ODE fallback). The
        # cores are shared with the Chebyshev cell (#132/#133), so alg_data is a
        # ChebyshevData; _build_magnus_ket_data picks the Magnus bracket/sizing.
        alg_data = _build_magnus_ket_data(
            alg,
            sys,
            u_dim,
            control_dim,
            ketdim,
            order,
            traj,
            u,
            global_dim,
            global_names,
        )
    else
        # ChebyshevAlg path: use shared _build_alg_data dispatch
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

    # Build sensitivity ODE factory from system drives.
    # Two modes:
    #   - use_ket_sensitivity=false (default): n×n augmented sensitivity ODE
    #     (state = [Φ_vec; S_1_vec; ...; S_n_vec], n_blocks * Φ_dim entries)
    #   - use_ket_sensitivity=true: ket-level sensitivity ODE
    #     (state = [ψ; s_1; ...; s_n], n_blocks * ketdim entries — strict reduction)
    if alg isa ChebyshevAlg
        # ADR-0003 Decision 2: ChebyshevAlg owns its sensitivity method — the Tsit5
        # augmented sensitivity ODE must not be built as a silent fallback. Until
        # the analytic Chebyshev Duhamel sensitivity lands (#133), any sensitivity
        # request errors loudly (see compute_ode_jacobian!).
        sens_probs = nothing
        sens_state = nothing
        ket_sens_results = nothing
    elseif use_ket_sensitivity
        make_f!_ket_sens, n_params = build_ket_sensitivity_ode(
            H_drift_op,
            sys.H_drives,
            u_dim,
            ketdim,
            order,
            1,  # n_kets = 1
        )

        if !isnothing(make_f!_ket_sens)
            sens_probs, _ = build_ket_sensitivity_problems(
                make_f!_ket_sens,
                n_params,
                ketdim,
                1,  # n_kets
                [p₀; Δt₀; t₀],
                traj.N - 1,
            )
            ket_sens_results = [zeros(ComplexF64, ketdim, n_params) for _ = 1:(traj.N-1)]
            sens_state = nothing
        else
            sens_probs = nothing
            sens_state = nothing
            ket_sens_results = nothing
        end
    else
        make_f!_sens, n_params =
            build_sensitivity_ode(H_drift_op, sys.H_drives, u_dim, ketdim, order)

        if !isnothing(make_f!_sens)
            sens_probs, sens_state = build_sensitivity_problems(
                make_f!_sens,
                n_params,
                ketdim,
                [p₀; Δt₀; t₀],
                traj.N,
            )

            # For Rodas5P, rebuild the sensitivity problems with the closed-form
            # Jacobian wired in (block-diagonal -iΔt(I⊗H) per the Kronecker structure
            # of vec(H·M)). Same rationale as the forward path: ForwardDiff rejects
            # ComplexF64 ODE states.
            if alg isa Rodas5PAlg
                make_sens_jac! = build_sensitivity_ode_jacobian(
                    H_drift_op,
                    sys.H_drives,
                    u_dim,
                    ketdim,
                    order,
                    n_params,
                )
                sens_x₀ = sens_probs[1].u0  # already initialized to [vec(I); zeros]
                sens_probs = [
                    ODEProblem(
                        ODEFunction(make_f!_sens(); jac = make_sens_jac!()),
                        sens_x₀,
                        (0.0, 1.0),
                        [p₀; Δt₀; t₀],
                    ) for _ = 1:(traj.N-1)
                ]
            end
        else
            sens_probs = nothing
            sens_state = nothing
        end
        ket_sens_results = nothing
    end

    # Build second-order sensitivity ODE for exact Hessian
    if exact_hessian
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

    return SplineIntegrator{KetTrajectory,S,ComplexF64,typeof(alg),typeof(alg_data)}(
        [x],
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
        use_ket_sensitivity,
        ket_sens_results,
        alg,
        alg_data,
    )
end

# ============================================================================ #
# Call operator for KetTrajectory
# ============================================================================ #

@views function (𝒮::SplineIntegrator{KetTrajectory})(
    δₖ::AbstractVector,
    zₖ::KnotPoint,
    zₖ₊₁::KnotPoint,
    k::Int,
    globals::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    xₖ = zₖ[x_name(𝒮)]
    xₖ₊₁ = zₖ₊₁[x_name(𝒮)]

    # Build parameter vector using unified function
    pₖ = build_ode_params(𝒮, zₖ, zₖ₊₁, globals)

    # Algorithm-dispatched forward solve
    if 𝒮.alg isa Tsit5Alg || 𝒮.alg isa Rodas5PAlg
        # State-forward path: convert iso ket → complex, propagate, convert back
        ψₖ = iso_to_complex_ket(xₖ)

        if 𝒮.use_ket_sensitivity
            # use_ket_sensitivity path: Φ_probs are identity-init, so the solve
            # returns Φ_k itself (a matrix). Cache Φ_k for the Jacobian assembly,
            # then apply to ψₖ to get f_k.
            data = 𝒮.alg_data::Tsit5Data
            prob = remake(data.Φ_probs[k], p = pₖ)
            Φₖ_complex = _solve_forward_tsit5(prob, 𝒮.alg, 𝒮.tol).u[end]
            copyto!(𝒮.prop_results[k].Φ_vec, vec(Φₖ_complex))
            fₖ_complex = Φₖ_complex * ψₖ
        elseif 𝒮.alg isa Rodas5PAlg
            # Slice 3b de-scope (director, 2026-08-30): the concrete Rodas5P
            # solve is proprietary — this hook is attached by Piccolissimo
            # (`_stiff_rodas5p_solve`), byte-identical to the pre-split call.
            data = 𝒮.alg_data::Rodas5PData
            prob = remake(data.Φ_probs[k], u0 = ψₖ, p = pₖ)
            fₖ_complex = _stiff_rodas5p_solve(prob, 𝒮.tol; saveat = 1.0).u[end]
        else
            data = 𝒮.alg_data::Tsit5Data
            prob = remake(data.Φ_probs[k], u0 = ψₖ, p = pₖ)
            fₖ_complex = _solve_forward_tsit5(prob, 𝒮.alg, 𝒮.tol).u[end]
        end
        fₖ = complex_ket_to_iso(fₖ_complex)
    elseif 𝒮.alg_data isa ChebyshevData
        # Matrix-free exp-action forward (#132/#223): sub-stepped midpoint-Magnus
        # sweep through the standalone cores (magnus_forward! / chebyshev_expv!) —
        # the propagator Φ is never materialized. Serves both ChebyshevAlg and the
        # Magnus spline cell (MagnusGL4Alg / MagnusAdapt4Alg on KetTrajectory),
        # which share these cores so their VJP is the exact adjoint of THIS forward.
        ψₖ = iso_to_complex_ket(xₖ)
        fₖ_complex =
            _chebyshev_forward!(𝒮.alg_data::ChebyshevData, _spline_type(𝒮), k, pₖ, ψₖ)
        fₖ = complex_ket_to_iso(fₖ_complex)
    else
        # Propagator-forward path for Magnus
        Φₖ_raw = _forward_propagate(𝒮, k, pₖ)
        Φₖ = eltype(Φₖ_raw) <: Complex ? propagator_to_iso(Φₖ_raw) : Φₖ_raw
        fₖ = Φₖ * xₖ
    end

    # Ket: direct constraint
    δₖ[:] = xₖ₊₁ - fₖ

    return nothing
end

# ============================================================================ #
# Dispatch helpers for KetTrajectory
# ============================================================================ #

"""
    compute_ode_jacobian!(𝒮::SplineIntegrator{KetTrajectory}, zₖ, zₖ₊₁, k, globals)

Compute ODE Jacobian for KetTrajectory using analytical sensitivity equations.
Stores Φₖ and ∂ₚΦₖ in `𝒮.ode_jac_results[k]`.
"""
@views function compute_ode_jacobian!(
    𝒮::SplineIntegrator{KetTrajectory},
    zₖ::KnotPoint,
    zₖ₊₁::KnotPoint,
    k::Int,
    globals::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    𝒮.alg isa ChebyshevAlg && error(
        "ChebyshevAlg has no analytic sensitivity method yet (#133): refusing to fall " *
        "back to the Tsit5 augmented sensitivity ODE (ADR-0003 Decision 2). Use " *
        "Tsit5Alg for gradient-based optimization until the Chebyshev sensitivity lands.",
    )
    isnothing(𝒮.sens_probs) && error(
        "SplineIntegrator requires explicit H_drift/H_drives matrices for Jacobian computation. " *
        "Function-based systems are not supported.",
    )
    pₖ = build_ode_params(𝒮, zₖ, zₖ₊₁, globals)
    if 𝒮.use_ket_sensitivity
        # Ket-level sensitivity: Φ is cached from the call operator (evaluate!),
        # solve for per-ket sensitivities only. Mirrors MultiKetTrajectory pattern
        # at K=1.
        _compute_ket_sensitivity!(𝒮, zₖ, pₖ, k)
    else
        _compute_ode_jacobian_analytical!(𝒮, pₖ, k)
    end
    return nothing
end

"""
    _compute_ket_sensitivity!(𝒮::SplineIntegrator{KetTrajectory}, zₖ, pₖ, k)

Solve the ket-level sensitivity ODE for knot point `k` (single-ket case, K=1).

Sets up initial conditions `[ψ(0), 0, ..., 0]` from the trajectory, solves the
augmented ket-level ODE, and extracts ket-level sensitivity vectors `s_j` into
`ket_sens_results[k]`. Mirrors the K>1 implementation in `spline_integrator_multiket.jl`.
"""
function _compute_ket_sensitivity!(
    𝒮::SplineIntegrator{KetTrajectory},
    zₖ::KnotPoint,
    pₖ::AbstractVector,
    k::Int,
)
    n = 𝒮.ketdim
    n_params = size(𝒮.ket_sens_results[k], 2)

    # Solve the propagator-level forward ODE to populate Φ in prop_results[k].
    # This makes eval_jacobian robust to being called without a prior evaluate!,
    # at the cost of one extra n×n forward ODE per knot. The cost is negligible
    # vs the n×n sensitivity ODE the propagator-level path would otherwise run.
    data = 𝒮.alg_data::Tsit5Data
    prob_Φ = remake(data.Φ_probs[k], p = pₖ)
    Φₖ_complex = solve(prob_Φ, Tsit5(); abstol = 𝒮.tol, reltol = 𝒮.tol, saveat = 1.0).u[end]
    copyto!(𝒮.prop_results[k].Φ_vec, vec(Φₖ_complex))

    # Build initial condition for ket-level sensitivity ODE: [ψ; zeros(n * n_params)]
    state_dim = n * (1 + n_params)
    u0 = zeros(ComplexF64, state_dim)
    ψ̃ = zₖ[x_name(𝒮)]
    iso_to_complex_ket!(@view(u0[1:n]), ψ̃)

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
    # ODE state layout (K=1): [ψ, s_1, s_2, ..., s_{n_params}]
    # — first n entries are the ket; subsequent blocks of n entries are sensitivities.
    result = 𝒮.ket_sens_results[k]
    @inbounds for j = 1:n_params
        off_s = n + (j - 1) * n
        for r = 1:n
            result[r, j] = final_state[off_s+r]
        end
    end

    return nothing
end

"""
    compute_ode_hessian!(𝒮::SplineIntegrator{KetTrajectory}, zₖ, zₖ₊₁, k, globals)

Compute exact second-order sensitivities for KetTrajectory.
Stores Φₖ and ∂ₚΦₖ in `𝒮.prop_results[k]`, returns T_{ij} matrices.
"""
@views function compute_ode_hessian!(
    𝒮::SplineIntegrator{KetTrajectory},
    zₖ::KnotPoint,
    zₖ₊₁::KnotPoint,
    k::Int,
    globals::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    isnothing(𝒮.hess_probs) &&
        error("SplineIntegrator exact_hessian=true required for compute_ode_hessian!")
    pₖ = build_ode_params(𝒮, zₖ, zₖ₊₁, globals)
    return _compute_ode_hessian_analytical!(𝒮, pₖ, k)
end

"""
    get_state_vector(𝒮::SplineIntegrator{KetTrajectory}, zₖ)

Extract the state vector from a knot point for KetTrajectory.
"""
function get_state_vector(𝒮::SplineIntegrator{KetTrajectory}, zₖ::KnotPoint)
    return zₖ[x_name(𝒮)]
end

# ============================================================================ #
# eval_jacobian override for use_ket_sensitivity=true on KetTrajectory
#
# When ket-level sensitivity is enabled, the Jacobian assembly differs from the
# shared eval_jacobian path: instead of reading propagator-level ∂ₚΦₖ from
# prop_results, we read ket-level sᵢⱼ vectors from ket_sens_results[k]. The
# −Φₖ block uses Φ cached during the call operator (evaluate!).
#
# Mirrors the K>1 implementation in `spline_integrator_multiket.jl` at K=1.
# ============================================================================ #

@views function eval_jacobian(𝒮::SplineIntegrator{KetTrajectory}, traj::NamedTrajectory)
    if !𝒮.use_ket_sensitivity
        # Default propagator-level path: delegate to the shared method.
        return invoke(eval_jacobian, Tuple{SplineIntegrator,NamedTrajectory}, 𝒮, traj)
    end

    N = traj.N
    z_dim = traj.dim
    global_dim = traj.global_dim
    x_dim = 𝒮.x_dim
    F_dim = x_dim * (N - 1)
    Z_dim = z_dim * N + global_dim

    # Build full Jacobian
    ∂F = spzeros(F_dim, Z_dim)
    temp_vec = Vector{Float64}(undef, x_dim)

    # Extract globals once
    globals = extract_globals(𝒮, traj)

    # Compute ODE Jacobians in parallel (fills ket_sens_results[k] and uses Φ
    # already cached in prop_results[k] from the call operator).
    Threads.@threads for k = 1:(N-1)
        compute_ode_jacobian!(𝒮, traj[k], traj[k+1], k, globals)
    end

    pdim = propagator_side_dim(𝒮)
    n = pdim
    state_name_sym = x_name(𝒮)

    @inbounds for k = 1:(N-1)
        rows = slice(k, x_dim)
        zₖ = traj[k]
        x_comps = zₖ.components[state_name_sym]
        zk_offset = (k - 1) * z_dim

        # ∂xₖ𝒮 = -Φ̃ₖ — cached during evaluate!
        Φₖ = get_propagator(𝒮.prop_results[k], pdim)
        neg_Φₖ_iso = -propagator_to_iso(Φₖ)
        ∂F[rows, zk_offset .+ x_comps] = neg_Φₖ_iso

        # ∂xₖ₊₁𝒮 = I
        ∂F[rows, zk_offset+z_dim .+ x_comps] = I(x_dim)

        # ∂pₖ𝒮 = -ket_to_iso(s_j) — directly from ket-level sensitivity vectors
        traj_cols, ode_indices = get_param_indices(𝒮, traj, k)
        for (col, ode_idx) in zip(traj_cols, ode_indices)
            sⱼ = 𝒮.ket_sens_results[k][1:n, ode_idx]
            for r = 1:n
                temp_vec[r] = -real(sⱼ[r])
                temp_vec[n+r] = -imag(sⱼ[r])
            end
            for (ii, r) in enumerate(rows)
                ∂F[r, col] += temp_vec[ii]
            end
        end
    end

    return ∂F
end
