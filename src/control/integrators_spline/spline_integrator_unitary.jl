# ============================================================================ #
# The matrix-free interface properties `SplineIntegrator{UnitaryTrajectory}`
# declares (#338 AC3).
#
# CORE WIDTH `d` IN **ONE** COMPONENT. The `d×d` propagator IS the state, carried
# in a single trajectory component of dim `2d²`. So the `d×K` core is `K = d`
# while the component count is `1` — exactly the mapping `GPUJacobianOp`'s
# `K = length(integ.x_names)` gets WRONG (it would read `K = 1`), and the reason
# ADR-0009 makes core width and component layout separate inputs.
#
# Reuse is at the CORE, never at the layout: `iso_vec_to_operator` is
# column-interleaved, so a Unitary state is `d` kets of width `2d` concatenated —
# byte-for-byte MultiKet's per-component layout. What a Unitary inner kernel has
# to write is the column mapping (one component, strided), not new numerics.
# Landing that kernel is slice #345; #338 lands only the declared mapping, so
# that nothing derives the width from the component count.
# ============================================================================ #


# ============================================================================ #
# Constructor for UnitaryTrajectory
# ============================================================================ #

function SplineIntegrator(
    qtraj::UnitaryTrajectory,
    N::Int;
    spline_order::Union{Int,Nothing} = nothing,
    alg::IntegrationAlgorithm = Tsit5Alg(),
    global_names::Union{Vector{Symbol},Nothing} = nothing,
)
    return _spline_unitary(
        get_system(qtraj),
        get_pulse(qtraj),
        state_name(qtraj),
        drive_name(qtraj),
        NamedTrajectory(qtraj, N);
        spline_order = spline_order,
        alg = alg,
        global_names = global_names,
    )
end

"""
    _spline_unitary(sys, pulse, x, u, traj; kwargs...)

Core of the UnitaryTrajectory SplineIntegrator constructor, parameterized by the
system, pulse, state/control names, and the FINAL NamedTrajectory. Shared by the
public constructor and the per-member SamplingTrajectory dispatch (#401).
"""
function _spline_unitary(
    sys,
    pulse,
    x::Symbol,
    u::Symbol,
    traj::NamedTrajectory;
    spline_order::Union{Int,Nothing} = nothing,
    alg::IntegrationAlgorithm = Tsit5Alg(),
    global_names::Union{Vector{Symbol},Nothing} = nothing,
)
    if alg isa ChebyshevAlg
        error(
            "ChebyshevAlg supports KetTrajectory only (#132): the unitary state IS the " *
            "d×d operator, so matrix-free Chebyshev buys no memory win (ADR-0003); " *
            "MultiKetTrajectory is the planned fast-follow.",
        )
    end
    # Infer spline type from pulse if possible, otherwise use provided order or default
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

    # Auto-detect globals if not specified
    if isnothing(global_names)
        if !isempty(sys.global_params)
            global_names = collect(keys(sys.global_params))
        else
            global_names = Symbol[]
        end
    end

    @assert traj.N > 1 "Trajectory must have at least two timesteps."

    # Determine global dimensions
    # Use traj.global_dim as the authoritative source when available
    if traj.global_dim > 0
        global_dim = traj.global_dim
    elseif !isempty(global_names)
        # Fallback: compute from global_names (each is scalar)
        global_dim = length(global_names)
    else
        global_dim = 0
    end

    # u_dim includes BOTH controls and globals for unified ODE parameter handling
    control_dim = traj.dims[u]
    u_dim = control_dim + global_dim

    # Pre-compute G matrices for MagnusGL4 (real isomorphic form)
    G_drift = Isomorphisms.G(_drift_matrix(sys.H_drift))
    G_drives = [Isomorphisms.G(drive_matrix(d)) for d in sys.H_drives]

    # Operators for complex-domain sensitivity ODE (operator-aware)
    H_drift_op = _to_operator(sys.H_drift)
    H_drive_ops =
        isempty(sys.H_drives) ? AbstractDynamicsOperator[] :
        [dynamics_operator(d) for d in sys.H_drives]


    ketdim = sys.levels

    # ============================================================================
    # Sensitivity-augmented ODE for Jacobian computation (analytical derivatives)
    # Built via shared builder in spline_integrator_type.jl
    # ============================================================================

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

    # Build sensitivity ODE factory from system drives (supports both linear and nonlinear)
    make_f!_sens, n_params =
        build_sensitivity_ode(H_drift_op, sys.H_drives, u_dim, ketdim, order)

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
    # p₀ format: [uₖ, uₖ₊₁, (duₖ, duₖ₊₁), Δtₖ, tₖ] where uₖ = [controls..., globals...]
    u_control_init = traj.bounds[u][2]

    # Initialize globals part if present
    if global_dim > 0
        g₀ = vcat(
            [traj.global_data[traj.global_components[name]] for name in global_names]...,
        )
        u_init = [u_control_init; g₀]  # Combine controls + globals
    else
        u_init = u_control_init
    end

    if order == 1
        p₀ = [u_init; u_init]  # [uₖ, uₖ₊₁] each containing [controls..., globals...]
    elseif order == 3
        # For cubic splines: [uₖ, uₖ₊₁, duₖ, duₖ₊₁]
        du = Symbol("d", u)
        # Initialize du to zeros if bounds don't exist yet (controls only, not globals)
        if haskey(traj.bounds, du)
            du_control_init = traj.bounds[du][2]
        else
            du_control_init = zeros(control_dim)
        end
        # du has same dimension as u (controls+globals, but globals don't have derivatives)
        du_init = [du_control_init; zeros(global_dim)]

        p₀ = [
            u_init;   # uₖ = [controls..., globals...]
            u_init;   # uₖ₊₁ = [controls..., globals...]
            du_init;  # duₖ = [control_derivs..., zeros for globals]
            du_init   # duₖ₊₁ = [control_derivs..., zeros for globals]
        ]
    end

    Δt₀ = 1.0
    t₀ = 0.0

    # ODE parameter count: [uₖ, uₖ₊₁, (duₖ, duₖ₊₁), Δt, t]
    ode_param_count = p_dim + 2

    # Complex-domain PropagatorResults
    prop_results =
        [PropagatorResult{ComplexF64}(ketdim, ode_param_count) for _ = 1:(traj.N-1)]

    # Build sensitivity ODE for analytical Jacobian computation
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

    # ============================================================================
    # Algorithm-specific data
    # ============================================================================

    tol = alg.tol
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

    return SplineIntegrator{UnitaryTrajectory,S,ComplexF64,typeof(alg),typeof(alg_data)}(
        [x],  # Wrap in vector for unified API
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
        false,
        nothing,
        nothing,
        nothing,
        0,  # exact_hessian (not yet supported for Unitary)
        false,
        nothing,  # ket-level sensitivity (MultiKetTrajectory only)
        alg,
        alg_data,
    )
end

# ============================================================================ #
# Algorithm-specific data builders
# ============================================================================ #

function _build_alg_data(
    alg::Tsit5Alg,
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
    # Build forward ODE factory from system drives (each call produces a thread-safe closure)
    make_f! = build_forward_ode(H_drift_op, sys.H_drives, u_dim, order)

    x₀ = Matrix{ComplexF64}(I, ketdim, ketdim)
    Φ_probs = [ODEProblem(make_f!(), x₀, (0.0, 1.0), [p₀; Δt₀; t₀]) for _ = 1:(traj.N-1)]

    # Dense pattern for both adaptive and fixed-step modes. The fixed-step
    # branch used to build the pattern by probing Φ_probs[1], but the
    # complex-valued probe never matched Tsit5Data's Float64-typed Φ_structure
    # (construction crashed — #180), and a single-point probe at the initial
    # parameters is not a sound structural pattern anyway.
    Φ_structure = sparse(ones(ketdim, ketdim))

    return Tsit5Data{ComplexF64}(Φ_probs, Φ_structure)
end

function _build_alg_data(
    alg::MagnusGL4Alg,
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

    # Pre-allocate per-interval buffers for thread-safe parallel evaluation
    G_drift_copies = [copy(G_drift) for _ = 1:(traj.N-1)]
    G_drives_copies = [[copy(G) for G in G_drives] for _ = 1:(traj.N-1)]
    G_buffers = [Matrix{Float64}(G_drift) for _ = 1:(traj.N-1)]  # Dense for exp(A)

    return MagnusGL4Data(
        G_drift,
        G_drives,
        G_drift_copies,
        G_drives_copies,
        G_buffers,
        sys.H_drives,
    )
end

function _build_alg_data(
    alg::MagnusAdapt4Alg,
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

    # Same buffer structure as MagnusGL4 — the operator construction is identical
    G_drift_copies = [copy(G_drift) for _ = 1:(traj.N-1)]
    G_drives_copies = [[copy(G) for G in G_drives] for _ = 1:(traj.N-1)]
    G_buffers = [Matrix{Float64}(G_drift) for _ = 1:(traj.N-1)]

    return MagnusAdapt4Data(
        G_drift,
        G_drives,
        G_drift_copies,
        G_drives_copies,
        G_buffers,
        sys.H_drives,
    )
end

# ============================================================================ #
# Forward propagation dispatch
# ============================================================================ #

"""
    _forward_propagate(𝒮::SplineIntegrator{UnitaryTrajectory,...,Tsit5Alg}, k, pₖ)

Forward-propagate using Tsit5 Runge-Kutta (adaptive or fixed-step per `𝒮.alg`).
Returns the propagator matrix at τ=1.
"""
function _forward_propagate(
    𝒮::SplineIntegrator{UnitaryTrajectory,S,R,Tsit5Alg},
    k::Int,
    pₖ::AbstractVector,
) where {S,R}
    data = 𝒮.alg_data::Tsit5Data{R}
    prob = remake(data.Φ_probs[k], p = pₖ)
    return _solve_forward_tsit5(prob, 𝒮.alg, 𝒮.tol).u[end]
end

"""
    _forward_propagate(𝒮::SplineIntegrator{UnitaryTrajectory,...,MagnusGL4Alg}, k, pₖ)

Forward-propagate using 4th-order Magnus integrator (MagnusGL4).
Preserves Lie group structure. Returns the propagator as a dense matrix.
"""
function _forward_propagate(
    𝒮::SplineIntegrator{<:Any,S,R,MagnusGL4Alg},
    k::Int,
    pₖ::AbstractVector,
) where {S,R}
    data = 𝒮.alg_data::MagnusGL4Data
    ketdim = 𝒮.ketdim
    u_dim = 𝒮.u_dim
    drives = data.drives
    n_drives = length(data.G_drives)
    order = spline_order(𝒮)
    T = eltype(pₖ)

    # Extract parameters
    uₖ = pₖ[1:u_dim]
    uₖ₊₁ = pₖ[(u_dim+1):2u_dim]
    if order == 3
        duₖ = pₖ[(2u_dim+1):3u_dim]
        duₖ₊₁ = pₖ[(3u_dim+1):4u_dim]
    end
    Δtₖ = pₖ[end-1]

    # Initial propagator (identity)
    Φ₀ = Matrix{T}(I, 2ketdim, 2ketdim)

    # Pre-allocate interpolated control vector (type-matched for ForwardDiff Duals)
    u_interp = Vector{T}(undef, u_dim)

    # Use pre-allocated buffers for Float64, fresh allocation for Dual (ForwardDiff)
    if T === Float64
        G_drift_local = data.G_drift_copies[k]
        G_drives_local = data.G_drives_copies[k]
        G_matrix = data.G_buffers[k]
    else
        G_drift_local = Matrix{T}(data.G_drift)
        G_drives_local = [Matrix{T}(G) for G in data.G_drives]
        G_matrix = Matrix{T}(data.G_drift)
    end

    # Update function for SciMLMatrixOperator
    update_func! =
        (A, _, _, t) -> begin
            τ = t / Δtₖ
            τ = clamp(τ, 0.0, 1.0)

            # Interpolate full control vector
            if order == 1
                onemτ = 1 - τ
                @inbounds for i = 1:u_dim
                    u_interp[i] = onemτ * uₖ[i] + τ * uₖ₊₁[i]
                end
            else  # order == 3
                τ2 = τ * τ
                τ3 = τ2 * τ
                h00 = 2τ3 - 3τ2 + 1
                h10 = (τ3 - 2τ2 + τ) * Δtₖ
                h01 = -2τ3 + 3τ2
                h11 = (τ3 - τ2) * Δtₖ
                @inbounds for i = 1:u_dim
                    u_interp[i] =
                        h00 * uₖ[i] + h10 * duₖ[i] + h01 * uₖ₊₁[i] + h11 * duₖ₊₁[i]
                end
            end

            # Build generator: G = G_drift + Σ drive_coeff(d, u) * G_d
            A .= G_drift_local
            @inbounds for i = 1:n_drives
                c = drive_coeff(drives[i], u_interp)
                A .+= c .* G_drives_local[i]
            end

            return nothing
        end

    A_op = SciMLMatrixOperator(G_matrix; (update_func!) = update_func!)
    prob = ODEProblem(A_op, Φ₀, (0.0, Δtₖ))

    dt = Δtₖ / 𝒮.alg.n_steps
    sol = solve(
        prob,
        MagnusGL4();
        dt = dt,
        save_everystep = false,
        abstol = 𝒮.tol,
        reltol = 𝒮.tol,
    )

    return sol.u[end]
end

"""
    _forward_propagate(𝒮::SplineIntegrator{UnitaryTrajectory,...,MagnusAdapt4Alg}, k, pₖ)

Forward-propagate using 4th-order adaptive Magnus integrator (MagnusAdapt4).
Preserves Lie group structure with automatic step size control.
"""
function _forward_propagate(
    𝒮::SplineIntegrator{<:Any,S,R,MagnusAdapt4Alg},
    k::Int,
    pₖ::AbstractVector,
) where {S,R}
    data = 𝒮.alg_data::MagnusAdapt4Data
    ketdim = 𝒮.ketdim
    u_dim = 𝒮.u_dim
    drives = data.drives
    n_drives = length(data.G_drives)
    order = spline_order(𝒮)
    T = eltype(pₖ)

    # Extract parameters
    uₖ = pₖ[1:u_dim]
    uₖ₊₁ = pₖ[(u_dim+1):2u_dim]
    if order == 3
        duₖ = pₖ[(2u_dim+1):3u_dim]
        duₖ₊₁ = pₖ[(3u_dim+1):4u_dim]
    end
    Δtₖ = pₖ[end-1]

    # Initial propagator (identity)
    Φ₀ = Matrix{T}(I, 2ketdim, 2ketdim)

    # Pre-allocate interpolated control vector
    u_interp = Vector{T}(undef, u_dim)

    # Use pre-allocated buffers for Float64, fresh allocation for Dual (ForwardDiff)
    if T === Float64
        G_drift_local = data.G_drift_copies[k]
        G_drives_local = data.G_drives_copies[k]
        G_matrix = data.G_buffers[k]
    else
        G_drift_local = Matrix{T}(data.G_drift)
        G_drives_local = [Matrix{T}(G) for G in data.G_drives]
        G_matrix = Matrix{T}(data.G_drift)
    end

    # Update function for SciMLMatrixOperator
    update_func! =
        (A, _, _, t) -> begin
            τ = t / Δtₖ
            τ = clamp(τ, 0.0, 1.0)

            # Interpolate full control vector
            if order == 1
                onemτ = 1 - τ
                @inbounds for i = 1:u_dim
                    u_interp[i] = onemτ * uₖ[i] + τ * uₖ₊₁[i]
                end
            else  # order == 3
                τ2 = τ * τ
                τ3 = τ2 * τ
                h00 = 2τ3 - 3τ2 + 1
                h10 = (τ3 - 2τ2 + τ) * Δtₖ
                h01 = -2τ3 + 3τ2
                h11 = (τ3 - τ2) * Δtₖ
                @inbounds for i = 1:u_dim
                    u_interp[i] =
                        h00 * uₖ[i] + h10 * duₖ[i] + h01 * uₖ₊₁[i] + h11 * duₖ₊₁[i]
                end
            end

            # Build generator: G = G_drift + Σ drive_coeff(d, u) * G_d
            A .= G_drift_local
            @inbounds for i = 1:n_drives
                c = drive_coeff(drives[i], u_interp)
                A .+= c .* G_drives_local[i]
            end

            return nothing
        end

    A_op = SciMLMatrixOperator(G_matrix; (update_func!) = update_func!)
    prob = ODEProblem(A_op, Φ₀, (0.0, Δtₖ))

    sol =
        solve(prob, MagnusAdapt4(); save_everystep = false, abstol = 𝒮.tol, reltol = 𝒮.tol)

    return sol.u[end]
end

# ============================================================================ #
# Call operator for UnitaryTrajectory (algorithm-agnostic)
# ============================================================================ #

@views function (𝒮::SplineIntegrator{UnitaryTrajectory})(
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
    Φₖ_raw = _forward_propagate(𝒮, k, pₖ)

    # Convert complex propagator to real iso form (Tsit5 returns complex, Magnus returns real)
    Φₖ = eltype(Φₖ_raw) <: Complex ? propagator_to_iso(Φₖ_raw) : Φₖ_raw

    # Unitary: iterate over blocks
    @inbounds for i = 1:𝒮.ketdim
        i_slice = slice(i, 2𝒮.ketdim)
        δₖ[i_slice] = xₖ₊₁[i_slice] - Φₖ * xₖ[i_slice]
    end

    return nothing
end

# ============================================================================ #
# Dispatch helpers for UnitaryTrajectory
# ============================================================================ #

"""
    compute_ode_jacobian!(𝒮::SplineIntegrator{UnitaryTrajectory}, zₖ, zₖ₊₁, k, globals)

Compute ODE Jacobian for UnitaryTrajectory using analytical sensitivity equations.
Stores complex Φₖ and ∂ₚΦₖ in `𝒮.prop_results[k]`.
"""
@views function compute_ode_jacobian!(
    𝒮::SplineIntegrator{UnitaryTrajectory},
    zₖ::KnotPoint,
    zₖ₊₁::KnotPoint,
    k::Int,
    globals::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    isnothing(𝒮.sens_probs) && error(
        "SplineIntegrator requires explicit H_drift/H_drives matrices for Jacobian computation. " *
        "Function-based systems are not supported.",
    )
    pₖ = build_ode_params(𝒮, zₖ, zₖ₊₁, globals)
    return _compute_ode_jacobian_analytical!(𝒮, pₖ, k)
end

"""
    get_state_vector(𝒮::SplineIntegrator{UnitaryTrajectory}, zₖ)

Extract the state vector from a knot point for UnitaryTrajectory.
"""
function get_state_vector(𝒮::SplineIntegrator{UnitaryTrajectory}, zₖ::KnotPoint)
    return zₖ[x_name(𝒮)]
end

"""
    fill_state_jacobian!(∂F, rows, x_cols, Φₖ, 𝒮::SplineIntegrator{UnitaryTrajectory})

Fill the ∂𝒮/∂xₖ block for UnitaryTrajectory: block-diagonal -propagator_to_iso(Φₖ).
Each column of the unitary evolves independently as a ket.
"""
function fill_state_jacobian!(
    ∂F,
    rows,
    x_cols,
    Φₖ::AbstractMatrix{<:Complex},
    𝒮::SplineIntegrator{UnitaryTrajectory},
)
    neg_Φₖ_iso = -propagator_to_iso(Φₖ)
    @inbounds for i = 1:𝒮.ketdim
        i_slice = slice(i, 2𝒮.ketdim)
        ∂F[rows[i_slice], x_cols[i_slice]] = neg_Φₖ_iso
    end
end

"""
    fill_param_jacobian!(∂F, rows, col, ∂ₚΦₖ, ode_idx, xₖ, temp_vec, tmp_complex, ψ_complex, 𝒮::SplineIntegrator{UnitaryTrajectory})

Fill a single parameter column for UnitaryTrajectory.
Complex sensitivity Sⱼ applied column-wise: each column of U is an iso ket.
"""
function fill_param_jacobian!(
    ∂F,
    rows,
    col,
    ∂ₚΦₖ::AbstractArray{<:Complex,3},
    ode_idx,
    xₖ,
    temp_vec,
    tmp_complex,
    ψ_complex,
    𝒮::SplineIntegrator{UnitaryTrajectory},
)
    n = 𝒮.ketdim
    Sⱼ = @view ∂ₚΦₖ[:, :, ode_idx]
    col_buf = @view temp_vec[1:2n]

    @inbounds for ℓ = 1:n
        iso_to_complex_ket!(ψ_complex, @view(xₖ[slice(ℓ, 2n)]))
        sensitivity_to_jac_col!(col_buf, Sⱼ, ψ_complex, tmp_complex)
        # sensitivity_to_jac_col! computes iso(-Sⱼ·ψ), already includes negative sign
        row_slice = slice(ℓ, 2n)
        for (i, r) in enumerate(rows[row_slice])
            ∂F[r, col] += col_buf[i]
        end
    end
end

"""
    fill_state_param_hessian!(μ∂²F, x_cols, col, ∂ₚΦₖ, ode_idx, μₖ, temp_vec, tmp1, tmp2, 𝒮::SplineIntegrator{UnitaryTrajectory})

Fill the (x_k, p) Hessian block for UnitaryTrajectory.
Complex sensitivity Sⱼ applied column-wise with Lagrange multiplier blocks.
"""
function fill_state_param_hessian!(
    μ∂²F,
    x_cols,
    col,
    ∂ₚΦₖ::AbstractArray{<:Complex,3},
    ode_idx,
    μₖ,
    temp_vec,
    tmp1,
    tmp2,
    𝒮::SplineIntegrator{UnitaryTrajectory},
)
    n = 𝒮.ketdim
    Sⱼ = @view ∂ₚΦₖ[:, :, ode_idx]
    col_buf = @view temp_vec[1:2n]

    @inbounds for ℓ = 1:n
        μₖ_ℓ = @view μₖ[slice(ℓ, 2n)]
        sensitivity_to_hess_col!(col_buf, Sⱼ, μₖ_ℓ, tmp1, tmp2)
        # sensitivity_to_hess_col! computes iso(-Sⱼ'·μ_C), already includes negative sign
        i_slice = slice(ℓ, 2n)
        for (j, xi) in enumerate(@view(x_cols[i_slice]))
            if xi <= col
                μ∂²F[xi, col] += col_buf[j]
            else
                μ∂²F[col, xi] += col_buf[j]
            end
        end
    end
end
