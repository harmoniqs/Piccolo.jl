# ============================================================================ #
# NonHermitianExponentialIntegrator dispatch over SamplingTrajectory
#
# Piccolissimo.jl#402 / parent #395. One integrator per ensemble member, built
# from the expanded sampling trajectory — mirroring the
# HermitianExponentialIntegrator sampling dispatch (#401) and Piccolo's
# `BilinearIntegrator(qtraj::SamplingTrajectory, N)` precedent.
#
# Supported open-system bases: DensityTrajectory, MultiDensityTrajectory.
# Closed-system bases route to the Hermitian family (out of scope here).
# ============================================================================ #

# ── Per-member construction helpers ────────────────────────────────────────── #
# These mirror `_sampling_hermitian_integrator` from the Hermitian side, building
# one NonHermitianExponentialIntegrator per member from the expanded trajectory.

function _sampling_nonhermitian_density(
    sys,
    x::Symbol,
    u::Symbol,
    traj::NamedTrajectory,
    global_names::Vector{Symbol};
    gauss_newton::Bool = false,
)
    @assert traj.N > 1 "Trajectory must have at least two timesteps."

    # Determine global dim from trajectory
    if !isempty(global_names) && traj.global_dim > 0
        global_dim = traj.global_dim
    elseif !isempty(global_names)
        global_dim = length(global_names)
    else
        global_dim = 0
    end

    n = sys.levels
    statedim = n^2  # compact density iso dimension

    x_dim = traj.dims[x]
    u_dim = traj.dims[u]
    N = traj.N
    dim = x_dim * (N - 1)
    var_dim = 2 * x_dim + u_dim + 1

    # Build compact Lindbladian generator G(u)
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

    # Sparsity templates
    ∂ℰ_template = jacobian_structure(DensityTrajectory, x, u, statedim, traj)
    μ∂²ℰ_template = hessian_structure(x_dim, u_dim, global_dim)
    ∂ℰs = [copy(∂ℰ_template) for _ = 1:(N-1)]
    μ∂²ℰs = [copy(μ∂²ℰ_template) for _ = 1:(N-1)]

    x_comps = traj.components[x]
    u_comps = traj.components[u]
    Δt_comp = traj.components[traj.timestep][1]
    var_comps = [collect(x_comps); collect(u_comps); Δt_comp]
    z_dim = traj.dim

    # Per-thread buffers
    nthr = Threads.maxthreadid()
    G_bufs = [zeros(Float64, statedim, statedim) for _ = 1:nthr]
    expG_bufs = [zeros(Float64, statedim, statedim) for _ = 1:nthr]

    return NonHermitianExponentialIntegrator{DensityTrajectory}(
        G_fn,
        [x],
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
        nothing,       # Id (unused for density)
        nothing,       # ∂u∂Δt_bufs
        nothing,       # ∂²u_bufs
        G_bufs,
        expG_bufs,
        global_names,
        global_dim,
        gauss_newton,
    )
end

function _sampling_nonhermitian_multidensity(
    sys,
    xs::Vector{Symbol},
    u::Symbol,
    traj::NamedTrajectory,
    global_names::Vector{Symbol};
    gauss_newton::Bool = false,
)
    @assert traj.N > 1 "Trajectory must have at least two timesteps."
    @assert all(name in traj.names for name in xs) "All state names must be in trajectory"

    if !isempty(global_names) && traj.global_dim > 0
        global_dim = traj.global_dim
    elseif !isempty(global_names)
        global_dim = length(global_names)
    else
        global_dim = 0
    end

    n = sys.levels
    statedim = n^2
    u_dim = traj.dims[u]
    x_dim = traj.dims[xs[1]]  # per-density state dim (= statedim)
    N = traj.N
    dim = x_dim * (N - 1)
    var_dim = 2 * x_dim + u_dim + 1

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

    nthr = Threads.maxthreadid()
    x_comps_all = vcat([collect(traj.components[name]) for name in xs]...)
    u_comps = traj.components[u]
    Δt_comp = traj.components[traj.timestep][1]
    var_comps = [x_comps_all; collect(u_comps); Δt_comp]
    z_dim = traj.dim

    integrators = NonHermitianExponentialIntegrator{MultiDensityTrajectory}[]
    for name in xs
        ∂ℰ_template = jacobian_structure(MultiDensityTrajectory, name, u, statedim, traj)
        μ∂²ℰ_template = hessian_structure(x_dim, u_dim, global_dim)
        ∂ℰs = [copy(∂ℰ_template) for _ = 1:(N-1)]
        μ∂²ℰs = [copy(μ∂²ℰ_template) for _ = 1:(N-1)]
        G_bufs = [zeros(Float64, statedim, statedim) for _ = 1:nthr]
        expG_bufs = [zeros(Float64, statedim, statedim) for _ = 1:nthr]

        push!(
            integrators,
            NonHermitianExponentialIntegrator{MultiDensityTrajectory}(
                G_fn,
                [name],
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
                nothing,
                nothing,
                nothing,
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

# ── Public sampling constructor ────────────────────────────────────────────── #

"""
    NonHermitianExponentialIntegrator(qtraj::SamplingTrajectory, N::Int; kwargs...)

Create a vector of NonHermitianExponentialIntegrators for each system in a
SamplingTrajectory — one per ensemble member, built from the expanded trajectory.
Mirrors `BilinearIntegrator(qtraj::SamplingTrajectory, N)`.

Supported open-system bases: `DensityTrajectory`, `MultiDensityTrajectory`.
Closed-system bases are not supported (use `HermitianExponentialIntegrator`).

# Keyword Arguments
- `global_names::Union{Vector{Symbol}, Nothing}=nothing`: global (time-invariant)
  variables in the dynamics. If `nothing`, auto-detects from the NOMINAL system's
  `global_params`.
- `gauss_newton::Bool = false`
"""
function NonHermitianExponentialIntegrator(
    qtraj::SamplingTrajectory,
    N::Int;
    global_names::Union{Vector{Symbol},Nothing} = nothing,
    gauss_newton::Bool = false,
)
    base = qtraj.base_trajectory
    nominal_sys = get_system(base)

    resolved_global_names = if !isnothing(global_names)
        global_names
    elseif !isempty(nominal_sys.global_params)
        collect(keys(nominal_sys.global_params))
    else
        Symbol[]
    end

    # Expanded sampling trajectory with globals attached
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

    control_sym = drive_name(qtraj)
    member_results = map(zip(qtraj.systems, sampling_member_states(qtraj))) do (sys, states)
        _nonhermitian_sampling_dispatch(
            base,
            sys,
            traj,
            states,
            control_sym,
            resolved_global_names;
            gauss_newton = gauss_newton,
        )
    end
    # Flatten: MultiDensity base returns one integrator per density per member
    return vcat((r isa AbstractVector ? r : [r] for r in member_results)...)
end

# ── Dispatch on base trajectory type ───────────────────────────────────────── #

function _nonhermitian_sampling_dispatch(
    base::DensityTrajectory,
    sys,
    traj::NamedTrajectory,
    x::Symbol,
    u::Symbol,
    global_names::Vector{Symbol};
    kwargs...,
)
    return _sampling_nonhermitian_density(sys, x, u, traj, global_names; kwargs...)
end

function _nonhermitian_sampling_dispatch(
    base::MultiDensityTrajectory,
    sys,
    traj::NamedTrajectory,
    xs::Vector{Symbol},
    u::Symbol,
    global_names::Vector{Symbol};
    kwargs...,
)
    return _sampling_nonhermitian_multidensity(sys, xs, u, traj, global_names; kwargs...)
end

function _nonhermitian_sampling_dispatch(
    base::AbstractQuantumTrajectory,
    sys,
    traj,
    states,
    control_sym,
    global_names;
    kwargs...,
)
    error(
        "NonHermitianExponentialIntegrator SamplingTrajectory dispatch only supports " *
        "open-system bases (DensityTrajectory, MultiDensityTrajectory). " *
        "For closed-system bases, use HermitianExponentialIntegrator.",
    )
end

# ============================================================================ #
# Tests
# ============================================================================ #

@testitem "NonHermitianExponentialIntegrator dispatch on SamplingTrajectory (Density)" begin
    using DirectTrajOpt, NamedTrajectories, Piccolo

    L = ComplexF64[0 0.1; 0 0]
    sys1 = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    sys2 =
        OpenQuantumSystem(0.95 * PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    ρ0 = ComplexF64[1 0; 0 0]
    ρg = ComplexF64[0 0; 0 1]
    N = 11
    times = collect(range(0, 1.0, length = N))
    pulse = ZeroOrderPulse(zeros(1, N), times)
    base_qtraj = DensityTrajectory(sys1, pulse, ρ0, ρg)
    sampling_qtraj = SamplingTrajectory(base_qtraj, [sys1, sys2])
    expanded_traj = NamedTrajectory(sampling_qtraj, N)

    integrators = NonHermitianExponentialIntegrator(sampling_qtraj, N)
    @test integrators isa Vector{<:NonHermitianExponentialIntegrator}
    @test length(integrators) == 2
    @test [x_name(ℰ) for ℰ in integrators] == [:ρ⃗̃1, :ρ⃗̃2]

    for ℰ in integrators
        test_integrator(ℰ, expanded_traj; atol = 1e-3)
    end
end

@testitem "NonHermitianExponentialIntegrator dispatch on SamplingTrajectory (MultiDensity)" begin
    using DirectTrajOpt, NamedTrajectories, Piccolo

    L = ComplexF64[0 0.1; 0 0]
    sys1 = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    sys2 =
        OpenQuantumSystem(0.95 * PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    ρ0₁ = ComplexF64[1 0; 0 0]
    ρg₁ = ComplexF64[0 0; 0 1]
    ρ0₂ = ComplexF64[0 0; 0 1]
    ρg₂ = ComplexF64[1 0; 0 0]
    N = 11
    times = collect(range(0, 1.0, length = N))
    pulse = ZeroOrderPulse(zeros(1, N), times)
    base_qtraj = MultiDensityTrajectory(sys1, pulse, [ρ0₁, ρ0₂], [ρg₁, ρg₂])
    sampling_qtraj = SamplingTrajectory(base_qtraj, [sys1, sys2])
    expanded_traj = NamedTrajectory(sampling_qtraj, N)

    integrators = NonHermitianExponentialIntegrator(sampling_qtraj, N)
    # 2 systems × 2 densities = 4 integrators
    @test integrators isa Vector{<:NonHermitianExponentialIntegrator}
    @test length(integrators) == 4
    @test [x_name(ℰ) for ℰ in integrators] == [:ρ⃗̃1, :ρ⃗̃2, :ρ⃗̃3, :ρ⃗̃4]

    for ℰ in integrators
        test_integrator(ℰ, expanded_traj; atol = 1e-3)
    end
end

@testitem "NonHermitianExponentialIntegrator sampling Jacobian parity (Density) (AC3)" begin
    using DirectTrajOpt, NamedTrajectories, Piccolo

    L = ComplexF64[0 0.1; 0 0]
    sys = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    ρ0 = ComplexF64[1 0; 0 0]
    ρg = ComplexF64[0 0; 0 1]
    N = 11
    times = range(0, 1.0, length = N)
    pulse = ZeroOrderPulse(randn(1, N) .* 0.1, times)

    qtraj = DensityTrajectory(sys, pulse, ρ0, ρg)
    direct_ℰ = NonHermitianExponentialIntegrator(qtraj, N)
    direct_traj = NamedTrajectory(qtraj, N)

    sampling_qtraj = SamplingTrajectory(qtraj, [sys])
    sampling_ℰ = NonHermitianExponentialIntegrator(sampling_qtraj, N)[1]
    expanded_traj = NamedTrajectory(sampling_qtraj, N)

    test_integrator(direct_ℰ, direct_traj; atol = 1e-3)
    test_integrator(sampling_ℰ, expanded_traj; atol = 1e-3)
end

@testitem "NonHermitianExponentialIntegrator sampling Jacobian parity (MultiDensity) (AC3)" begin
    using DirectTrajOpt, NamedTrajectories, Piccolo

    L = ComplexF64[0 0.1; 0 0]
    sys = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    ρ0₁ = ComplexF64[1 0; 0 0]
    ρg₁ = ComplexF64[0 0; 0 1]
    ρ0₂ = ComplexF64[0 0; 0 1]
    ρg₂ = ComplexF64[1 0; 0 0]
    N = 11
    times = range(0, 1.0, length = N)
    pulse = ZeroOrderPulse(randn(1, N) .* 0.1, times)

    qtraj = MultiDensityTrajectory(sys, pulse, [ρ0₁, ρ0₂], [ρg₁, ρg₂])
    direct_ℰs = NonHermitianExponentialIntegrator(qtraj, N)
    direct_traj = NamedTrajectory(qtraj, N)

    sampling_qtraj = SamplingTrajectory(qtraj, [sys])
    sampling_ℰs = NonHermitianExponentialIntegrator(sampling_qtraj, N)
    expanded_traj = NamedTrajectory(sampling_qtraj, N)

    @test length(direct_ℰs) == length(sampling_ℰs) == 2
    for (dℰ, sℰ) in zip(direct_ℰs, sampling_ℰs)
        test_integrator(dℰ, direct_traj; atol = 1e-3)
        test_integrator(sℰ, expanded_traj; atol = 1e-3)
    end
end

@testitem "NonHermitianExponentialIntegrator sampling per-member buffer isolation" begin
    using DirectTrajOpt, NamedTrajectories, Piccolo

    L = ComplexF64[0 0.1; 0 0]
    sys1 = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    sys2 =
        OpenQuantumSystem(0.95 * PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    ρ0 = ComplexF64[1 0; 0 0]
    ρg = ComplexF64[0 0; 0 1]
    N = 11
    times = range(0, 1.0, length = N)
    pulse = LinearSplinePulse(zeros(1, N), times)
    qtraj = DensityTrajectory(sys1, pulse, ρ0, ρg)
    sampling_qtraj = SamplingTrajectory(qtraj, [sys1, sys2])
    expanded_traj = NamedTrajectory(sampling_qtraj, N)

    ℰs = NonHermitianExponentialIntegrator(sampling_qtraj, N)
    @test length(ℰs) == 2

    F1 = zeros(ℰs[1].x_dim * (N - 1))
    F2 = zeros(ℰs[2].x_dim * (N - 1))
    evaluate!(F1, ℰs[1], expanded_traj)
    evaluate!(F2, ℰs[2], expanded_traj)
    F1_saved = copy(F1)
    F2_saved = copy(F2)
    evaluate!(F1, ℰs[1], expanded_traj)
    evaluate!(F2, ℰs[2], expanded_traj)
    @test F1 ≈ F1_saved
    evaluate!(F2, ℰs[2], expanded_traj)
    @test F1 ≈ F1_saved
end
