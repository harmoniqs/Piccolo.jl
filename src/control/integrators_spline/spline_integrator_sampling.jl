# ============================================================================ #
# SplineIntegrator dispatch over SamplingTrajectory (closed-system)
#
# Piccolissimo.jl#401 / parent #395. One integrator per ensemble member, built
# from the expanded sampling trajectory — mirroring Piccolo's
# `BilinearIntegrator(qtraj::SamplingTrajectory, N)` precedent and the
# HermitianExponentialIntegrator sampling dispatch alongside it. Per-member
# construction shares NO mutable scratch across members: each member integrator
# allocates its own ODE problem sets, PropagatorResults and alg_data inside the
# per-base cores (the #354 race-fix discipline).
# ============================================================================ #

"""
    SplineIntegrator(qtraj::SamplingTrajectory, N::Int; kwargs...)

Create a vector of SplineIntegrators for each system in a SamplingTrajectory —
one integrator per ensemble member, built from the expanded trajectory. Mirrors
`BilinearIntegrator(qtraj::SamplingTrajectory, N)`; unlike the bilinear path, a
multi-ket member gets ONE integrator over its ket sub-states (shared-propagator),
not one per ket.

Supported closed-system bases: `UnitaryTrajectory`, `KetTrajectory`,
`MultiKetTrajectory`. The spline order follows the base trajectory's pulse
(`LinearSplinePulse` → linear, `CubicSplinePulse` → cubic Hermite).

# Keyword Arguments
- `spline_order::Union{Int, Nothing}=nothing`: explicit order override for
  non-spline pulses (1 = linear, 3 = cubic Hermite).
- `alg::IntegrationAlgorithm=Tsit5Alg()`: propagation algorithm.
- `global_names::Union{Vector{Symbol}, Nothing}=nothing`: global (time-invariant)
  variables in the dynamics. If `nothing`, auto-detects from the NOMINAL system's
  `global_params` (members share the global names).
- `exact_hessian`, `use_ket_sensitivity`/`ket_sensitivity`, `second_order_shape`:
  as in the per-base constructors.
"""
function SplineIntegrator(
    qtraj::SamplingTrajectory,
    N::Int;
    spline_order::Union{Int,Nothing} = nothing,
    alg::IntegrationAlgorithm = Tsit5Alg(),
    global_names::Union{Vector{Symbol},Nothing} = nothing,
    kwargs...,
)
    base = qtraj.base_trajectory
    nominal_sys = get_system(base)

    # Resolve globals from the nominal system (members share the names), matching
    # the per-base constructors' auto-detection.
    resolved_global_names = if !isnothing(global_names)
        global_names
    elseif !isempty(nominal_sys.global_params)
        collect(keys(nominal_sys.global_params))
    else
        Symbol[]
    end

    # Expanded sampling trajectory. The sampling conversion does not attach
    # globals, so attach them here — mirroring the non-sampling conversions (which
    # auto-attach sys.global_params) and SamplingProblem's own propagation.
    traj = NamedTrajectory(qtraj, N)
    if !isempty(resolved_global_names)
        global_data = Dict{Symbol,Vector{Float64}}(
            name => [
                haskey(nominal_sys.global_params, name) ?
                nominal_sys.global_params[name] : 0.0,
            ] for name in resolved_global_names
        )
        traj = NamedTrajectory(
            traj.datavec,
            traj.components,
            traj.N;
            timestep = traj.timestep,
            controls = traj.control_names,
            bounds = traj.bounds,
            initial = traj.initial,
            final = isnothing(traj.final_) ? NamedTuple() : traj.final_,
            goal = traj.goal,
            global_data = global_data,
        )
    end

    pulse = get_pulse(qtraj)
    control_sym = drive_name(qtraj)
    return map(zip(qtraj.systems, sampling_member_states(qtraj))) do (sys, states)
        _sampling_spline_integrator(
            base,
            sys,
            pulse,
            states,
            control_sym,
            traj;
            spline_order = spline_order,
            alg = alg,
            global_names = resolved_global_names,
            kwargs...,
        )
    end
end

# Per-member construction helpers — dispatch on the base trajectory type, one
# integrator per member built from the expanded trajectory (the bilinear
# `_sampling_integrator` precedent).
function _sampling_spline_integrator(
    base::KetTrajectory,
    sys,
    pulse,
    x::Symbol,
    u::Symbol,
    traj::NamedTrajectory;
    kwargs...,
)
    return _spline_ket(sys, pulse, x, u, traj; kwargs...)
end

function _sampling_spline_integrator(
    base::UnitaryTrajectory,
    sys,
    pulse,
    x::Symbol,
    u::Symbol,
    traj::NamedTrajectory;
    kwargs...,
)
    return _spline_unitary(sys, pulse, x, u, traj; kwargs...)
end

function _sampling_spline_integrator(
    base::MultiKetTrajectory,
    sys,
    pulse,
    xs::Vector{Symbol},
    u::Symbol,
    traj::NamedTrajectory;
    kwargs...,
)
    return _spline_multiket(sys, pulse, xs, u, traj; kwargs...)
end

# ============================================================================ #
# Density sampling dispatch (Piccolissimo.jl#402; moved to Piccolo — slice 3b)
# ============================================================================ #

# ============================================================================ #
# Density sampling dispatch (Piccolissimo.jl#402)
# ============================================================================ #

# ── Internal constructors mirroring _spline_ket / _spline_multiket for density ─

function _spline_density(
    sys::OpenQuantumSystem,
    pulse,
    x::Symbol,
    u::Symbol,
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

    # Infer spline type from pulse
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

    n = sys.levels
    n² = n^2
    tol = alg.tol

    if isnothing(global_names)
        global_names = Symbol[]
    end

    @assert traj.N > 1 "Trajectory must have at least two timesteps."
    @assert x in traj.names "State name $x must be in trajectory"

    control_dim = traj.dims[u]

    if !isempty(global_names)
        global_comps = vcat([traj.global_components[name] for name in global_names]...)
        global_dim = length(global_comps)
    else
        global_dim = 0
    end
    u_dim = control_dim + global_dim

    H_drift = collect(sys.H_drift)
    drives = sys.H_drives
    n_terms = length(drives)
    H_drives = [collect(drive_matrix(d)) for d in drives]
    Ls = [collect(L) for L in sys.dissipation_operators]
    Ks = [sparse(L' * L) for L in sys.dissipation_operators]

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

    𝒢c_drift_ham, 𝒢c_drives, 𝒢c_dissipators = compact_lindbladian_parts(sys)
    𝒢c_drift = isempty(𝒢c_dissipators) ? 𝒢c_drift_ham : 𝒢c_drift_ham + sum(𝒢c_dissipators)

    x_dim = traj.dims[x]
    @assert x_dim == n² "State dimension ($x_dim) must equal n² = $n²"

    N_traj = traj.N
    dim = x_dim * (N_traj - 1)

    if order == 1

        p_dim = 2 * u_dim
    elseif order == 3

        p_dim = 4 * u_dim
    end

    Φc_init = vec(Matrix{Float64}(I, n², n²))
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
        du_control_init = haskey(traj.bounds, du) ? traj.bounds[du][2] : zeros(control_dim)
        du_init = [du_control_init; zeros(global_dim)]
        p₀ = [u_init; u_init; du_init; du_init]
    end
    Δt₀ = 1.0
    t₀ = 0.0

    Φ_probs = Vector{ODEProblem}(undef, N_traj - 1)
    for kk = 1:(N_traj-1)
        H_eff_kk = Matrix{ComplexF64}(undef, n, n)
        M_kk = Matrix{ComplexF64}(undef, n, n)
        dM_kk = Matrix{ComplexF64}(undef, n, n)
        tmp_kk = Matrix{ComplexF64}(undef, n, n)
        u_interp_kk = zeros(u_dim)

        f_kk! = if !isempty(H_drives) && order == 1
            (dx, xₖ, p, τ) -> begin
                uₖ = @view p[1:u_dim]
                uₖ₊₁ = @view p[(u_dim+1):2u_dim]
                Δtₖ = p[end-1]
                @inbounds for i = 1:u_dim
                    u_interp_kk[i] = (1 - τ) * uₖ[i] + τ * uₖ₊₁[i]
                end
                @. H_eff_kk = H_drift
                @inbounds for t_idx = 1:n_terms
                    c = drive_coeff(drives[t_idx], u_interp_kk)
                    @. H_eff_kk += c * H_drives[t_idx]
                end
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
                uₖ=@view p[1:u_dim]
                uₖ₊₁=@view p[(u_dim+1):2u_dim]
                duₖ=@view p[(2u_dim+1):3u_dim]
                duₖ₊₁=@view p[(3u_dim+1):4u_dim]
                Δtₖ=p[end-1]
                τ2=τ*τ
                τ3=τ2*τ
                h00=2τ3-3τ2+1
                h10=(τ3-2τ2+τ)*Δtₖ
                h01=-2τ3+3τ2
                h11=(τ3-τ2)*Δtₖ
                @inbounds for i = 1:u_dim
                    u_interp_kk[i]=h00*uₖ[i]+h10*duₖ[i]+h01*uₖ₊₁[i]+h11*duₖ₊₁[i]
                end
                @. H_eff_kk = H_drift
                @inbounds for t_idx = 1:n_terms
                    c=drive_coeff(drives[t_idx], u_interp_kk)
                    @. H_eff_kk+=c*H_drives[t_idx]
                end
                n_cols=length(xₖ)÷n²
                @inbounds for j = 1:n_cols
                    col=@view xₖ[((j-1)*n²+1):(j*n²)]
                    d_col=@view dx[((j-1)*n²+1):(j*n²)]
                    compact_iso_to_density!(M_kk, col, n)
                    lindblad_apply!(dM_kk, M_kk, H_eff_kk, Δtₖ, Ls, Ks, tmp_kk)
                    density_to_compact_iso!(d_col, dM_kk, n)
                end
                return nothing
            end
        else
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

    calculate_sparsity = alg isa Tsit5Alg && !alg.adaptive
    if calculate_sparsity
        Φ_sol = solve(Φ_probs[1], Tsit5(); abstol = tol, reltol = tol, saveat = 1.0).u[end]
        Φ_mat = reshape(Φ_sol, n², n²)
        Φ_structure = sparse(Φ_mat)
    else
        Φ_structure = sparse(ones(n², n²))
    end

    alg_data = DensityLindbladData{Float64}(
        Tsit5Data{Float64}(Φ_probs, Φ_structure),
        LindbladDuhamelTape(n, 1, u_dim, order, tol, H_drift, H_drives, drives, Ls, Ks),
    )
    ode_param_count = p_dim + 2
    prop_results = [PropagatorResult{Float64}(n², ode_param_count) for _ = 1:(N_traj-1)]

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
        [x],
        u,
        x_dim,
        u_dim,
        dim,
        tol,
        prop_results,
        n,
        global_names,
        global_dim,
        sens_probs,
        sens_state,
        false,
        nothing,
        nothing,
        nothing,
        0,
        false,
        nothing,
        alg,
        alg_data,
    )
end

function _spline_multidensity(
    sys::OpenQuantumSystem,
    pulse,
    x_names::Vector{Symbol},
    u::Symbol,
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
            "density dynamics stay on the Tsit5 augmented path (ADR-0003).",
        )
    end

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

    n = sys.levels
    n² = n^2
    tol = alg.tol
    n_densities = length(x_names)

    if isnothing(global_names)

        global_names = Symbol[]
    end

    @assert traj.N > 1 "Trajectory must have at least two timesteps."
    @assert all(name in traj.names for name in x_names) "All state names must be in trajectory"

    control_dim = traj.dims[u]
    if !isempty(global_names)
        global_comps = vcat([traj.global_components[name] for name in global_names]...)
        global_dim = length(global_comps)
    else
        global_dim = 0
    end
    u_dim = control_dim + global_dim

    H_drift = collect(sys.H_drift)
    drives = sys.H_drives
    n_terms = length(drives)
    H_drives = [collect(drive_matrix(d)) for d in drives]
    Ls = [collect(L) for L in sys.dissipation_operators]
    Ks = [sparse(L' * L) for L in sys.dissipation_operators]

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

    𝒢c_drift_ham, 𝒢c_drives, 𝒢c_dissipators = compact_lindbladian_parts(sys)
    𝒢c_drift = isempty(𝒢c_dissipators) ? 𝒢c_drift_ham : 𝒢c_drift_ham + sum(𝒢c_dissipators)

    single_x_dim = n²
    x_dim = n_densities * single_x_dim
    N_traj = traj.N
    dim = x_dim * (N_traj - 1)

    if order == 1

        p_dim = 2 * u_dim
    elseif order == 3

        p_dim = 4 * u_dim
    end

    Φc_init = vec(Matrix{Float64}(I, n², n²))
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
        du_control_init = haskey(traj.bounds, du) ? traj.bounds[du][2] : zeros(control_dim)
        du_init = [du_control_init; zeros(global_dim)]
        p₀ = [u_init; u_init; du_init; du_init]
    end
    Δt₀ = 1.0
    t₀ = 0.0

    Φ_probs = Vector{ODEProblem}(undef, N_traj - 1)
    for kk = 1:(N_traj-1)
        H_eff_kk = Matrix{ComplexF64}(undef, n, n)
        M_kk = Matrix{ComplexF64}(undef, n, n)
        dM_kk = Matrix{ComplexF64}(undef, n, n)
        tmp_kk = Matrix{ComplexF64}(undef, n, n)
        u_interp_kk = zeros(u_dim)

        f_kk! = if !isempty(H_drives) && order == 1
            (dx, xₖ, p, τ) -> begin
                uₖ = @view p[1:u_dim]
                uₖ₊₁ = @view p[(u_dim+1):2u_dim]
                Δtₖ = p[end-1]
                @inbounds for i = 1:u_dim
                    u_interp_kk[i] = (1 - τ) * uₖ[i] + τ * uₖ₊₁[i]
                end
                @. H_eff_kk = H_drift
                @inbounds for t_idx = 1:n_terms
                    c = drive_coeff(drives[t_idx], u_interp_kk)
                    @. H_eff_kk += c * H_drives[t_idx]
                end
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

    calculate_sparsity = alg isa Tsit5Alg && !alg.adaptive
    if calculate_sparsity
        Φ_sol = solve(Φ_probs[1], Tsit5(); abstol = tol, reltol = tol, saveat = 1.0).u[end]
        Φ_mat = reshape(Φ_sol, n², n²)
        Φ_structure = sparse(Φ_mat)
    else
        Φ_structure = sparse(ones(n², n²))
    end

    alg_data = DensityLindbladData{Float64}(
        Tsit5Data{Float64}(Φ_probs, Φ_structure),
        LindbladDuhamelTape(
            n,
            n_densities,
            u_dim,
            order,
            tol,
            H_drift,
            H_drives,
            drives,
            Ls,
            Ks,
        ),
    )
    ode_param_count = p_dim + 2
    prop_results = [PropagatorResult{Float64}(n², ode_param_count) for _ = 1:(N_traj-1)]

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

    return SplineIntegrator{MultiDensityTrajectory,S,Float64,typeof(alg),typeof(alg_data)}(
        x_names,
        u,
        x_dim,
        u_dim,
        dim,
        tol,
        prop_results,
        n,
        global_names,
        global_dim,
        sens_probs,
        sens_state,
        false,
        nothing,
        nothing,
        nothing,
        0,
        false,
        nothing,
        alg,
        alg_data,
    )
end

# ── Per-member sampling dispatch methods ──────────────────────────────────── #

function _sampling_spline_integrator(
    base::DensityTrajectory,
    sys,
    pulse,
    x::Symbol,
    u::Symbol,
    traj::NamedTrajectory;
    kwargs...,
)
    return _spline_density(sys, pulse, x, u, traj; kwargs...)
end

function _sampling_spline_integrator(
    base::MultiDensityTrajectory,
    sys,
    pulse,
    xs::Vector{Symbol},
    u::Symbol,
    traj::NamedTrajectory;
    kwargs...,
)
    return _spline_multidensity(sys, pulse, xs, u, traj; kwargs...)
end
