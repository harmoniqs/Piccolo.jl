# ============================================================================ #
# Cluster E1 (#347) — SplineIntegrator MultiKetTrajectory cell behavior tests.
#
# Ported/adapted from the pre-split embedded suite (Piccolissimo
# src/integrators/spline/spline_integrator_multiket.jl). The multiket cell sat
# at 20.8% (constructor + call operator only, via the template testitems).
#
# SEAM NOTE (open-core split, slice 3b): the multiket `eval_jacobian` /
# `eval_hessian_of_lagrangian` overrides read their per-knot index tables from
# `matrix_free_layout(𝒮, traj)` — a hook with NO method in Piccolo (the layout
# machinery is proprietary, attached by Piccolissimo's module). So the
# assembled-Jacobian/Hessian conformance items are NOT portable into this
# package's own suite; instead this file covers everything around the seam:
# constructor lanes (all kwargs + error gates), the forward paths (call
# operator, evaluate!, fixed-step vs adaptive), and the per-knot sensitivity
# solves called directly (`compute_ode_jacobian!` fills prop_results /
# ket_sens_results without touching the layout hook).
# ============================================================================ #

@testitem "E1: testing SplineIntegrator{MultiKetTrajectory}" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra

    T = 1.0
    N = 5
    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]

    initials = [ψ0, ψ1]  # |0⟩, |1⟩
    goals = [ψ1, ψ0]     # |1⟩, |0⟩

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.5, 2, N), times)
    ensemble_qtraj = MultiKetTrajectory(sys, pulse, initials, goals)

    𝒮 = SplineIntegrator(ensemble_qtraj, N)
    @test 𝒮 isa SplineIntegrator{MultiKetTrajectory,LinearSpline}
    @test length(𝒮.x_names) == 2
    @test 𝒮.x_names == state_names(ensemble_qtraj)

    combined_traj = NamedTrajectory(ensemble_qtraj, N)

    # Ensemble constraint dim: n_kets * 2*ketdim per knot interval
    @test 𝒮.x_dim == 2 * 2 * sys.levels
    @test 𝒮.dim == 𝒮.x_dim * (combined_traj.N - 1)

    # The forward residual is well-defined and nonzero at the rollout
    δ = zeros(𝒮.dim)
    evaluate!(δ, 𝒮, combined_traj)
    @test !all(iszero, δ)
    # ...and the rollout is dynamically consistent (small residual)
    @test norm(δ, Inf) < 1e-5
end

@testitem "E1: SplineIntegrator{MultiKetTrajectory} with 3 kets" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra

    T = 1.0
    N = 5
    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    ψ_plus = ComplexF64[1.0, 1.0] / sqrt(2)

    initials = [ψ0, ψ1, ψ0]
    goals = [ψ1, ψ0, ψ_plus]

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.5, 2, N), times)
    ensemble_qtraj = MultiKetTrajectory(sys, pulse, initials, goals)

    𝒮 = SplineIntegrator(ensemble_qtraj, N)
    @test 𝒮 isa SplineIntegrator{MultiKetTrajectory,LinearSpline}
    @test length(𝒮.x_names) == 3

    combined_traj = NamedTrajectory(ensemble_qtraj, N)
    δ = zeros(𝒮.dim)
    evaluate!(δ, 𝒮, combined_traj)
    @test norm(δ, Inf) < 1e-5
end

@testitem "E1: SplineIntegrator{MultiKetTrajectory} with CubicSpline" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra

    T = 1.0
    N = 5
    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]

    initials = [ψ0, ψ1]
    goals = [ψ1, ψ0]

    times = collect(range(0.0, T, N))
    pulse = CubicSplinePulse(fill(0.5, 2, N), fill(0.0, 2, N), times)
    ensemble_qtraj = MultiKetTrajectory(sys, pulse, initials, goals)

    𝒮 = SplineIntegrator(ensemble_qtraj, N)
    @test 𝒮 isa SplineIntegrator{MultiKetTrajectory,CubicSpline}
    @test length(𝒮.x_names) == 2

    combined_traj = NamedTrajectory(ensemble_qtraj, N)
    δ = zeros(𝒮.dim)
    evaluate!(δ, 𝒮, combined_traj)
    @test norm(δ, Inf) < 1e-5
end

@testitem "E1: SplineIntegrator{MultiKetTrajectory} efficiency vs separate integrators" begin
    using DirectTrajOpt
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra

    # Shared system
    T = 1.0
    N = 10
    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])

    n_kets = 4
    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]

    initials = [i % 2 == 1 ? ψ0 : ψ1 for i = 1:n_kets]
    goals = [i % 2 == 1 ? ψ1 : ψ0 for i = 1:n_kets]

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.5, 2, N), times)
    ensemble_qtraj = MultiKetTrajectory(sys, pulse, initials, goals)

    qtrajs = [KetTrajectory(sys, pulse, initials[i], goals[i]) for i = 1:n_kets]

    ensemble_integrator = SplineIntegrator(ensemble_qtraj, N)
    separate_integrators = [SplineIntegrator(qtraj, N) for qtraj in qtrajs]

    @test ensemble_integrator isa SplineIntegrator{MultiKetTrajectory,LinearSpline}
    @test length(separate_integrators) == n_kets

    combined_traj = NamedTrajectory(ensemble_qtraj, N)
    # Ensemble should have n_kets * x_dim constraints per knot point
    @test ensemble_integrator.dim == n_kets * 2 * sys.levels * (combined_traj.N - 1)

    # FORWARD PARITY (the shared-propagator contract): the ensemble integrator's
    # residual must equal the STACKED per-ket residuals — all kets ride the same
    # Φₖ, so the two computations must agree to solver precision.
    δ_ens = zeros(ensemble_integrator.dim)
    evaluate!(δ_ens, ensemble_integrator, combined_traj)
    for i = 1:n_kets
        𝒮 = separate_integrators[i]
        t = NamedTrajectory(qtrajs[i], N)
        δ_i = zeros(𝒮.dim)
        evaluate!(δ_i, 𝒮, t)
        # ensemble rows for ket i, knot k = interleaved per-ket blocks
        for k = 1:(N-1)
            rows_ens =
                ((i - 1) * 2 * sys.levels) .+ (1:(2*sys.levels)) .+
                (k - 1) * ensemble_integrator.x_dim
            rows_sep = (2 * sys.levels) .* (k - 1) .+ (1:(2*sys.levels))
            @test δ_ens[rows_ens] ≈ δ_i[rows_sep] atol = 1e-8
        end
    end
end

@testitem "E1: SplineIntegrator{MultiKetTrajectory} per-knot sensitivity solves: ket-level matches propagator-level" begin
    using DirectTrajOpt
    using DirectTrajOpt: evaluate!
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators:
        compute_ode_jacobian!, get_state_vectors

    # 3-level system (ketdim=3) with 2 kets — ket-level mode (K=2 < n=3)
    H_drift = diagm(ComplexF64[0.0, 1.0, 2.0])
    H_x = ComplexF64[0 1 0; 1 0 1; 0 1 0] / sqrt(2)
    sys = QuantumSystem(H_drift, [H_x], [1.0])

    T = 1.0
    N = 5

    ψ0 = ComplexF64[1.0, 0.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0, 0.0]

    initials = [ψ0, ψ1]
    goals = [ψ1, ψ0]

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.5, 1, N), times)
    ensemble_qtraj = MultiKetTrajectory(sys, pulse, initials, goals)
    traj = NamedTrajectory(ensemble_qtraj, N)

    𝒮_prop = SplineIntegrator(ensemble_qtraj, N; ket_sensitivity = false)
    𝒮_ket = SplineIntegrator(ensemble_qtraj, N; ket_sensitivity = true)
    @test 𝒮_ket.use_ket_sensitivity
    @test !isnothing(𝒮_ket.ket_sens_results)
    @test !𝒮_prop.use_ket_sensitivity

    # Direct per-knot sensitivity solves (no assembled Jacobian — that path
    # goes through the proprietary matrix_free_layout seam).
    for k = 1:(N-1)
        compute_ode_jacobian!(𝒮_prop, traj[k], traj[k+1], k, nothing)
        compute_ode_jacobian!(𝒮_ket, traj[k], traj[k+1], k, nothing)
    end

    pdim = 3
    n_params = 2 * 1 + 2  # linear: [u_k; u_k+1; Δt; t]
    for k = 1:(N-1)
        # The propagator-level solve filled Φ and ∂ₚΦ
        S = get_sensitivities(𝒮_prop.prop_results[k], pdim)
        Φ = get_propagator(𝒮_prop.prop_results[k], pdim)
        @test !iszero(Φ)

        # member states at knot k (iso-real → complex for the S_j·ψ contraction)
        states_k = get_state_vectors(𝒮_ket, traj[k])
        @test length(states_k) == 2
        for (i, ψ̃ᵢ) in enumerate(states_k)
            ψᵢ = ψ̃ᵢ[1:pdim] + im * ψ̃ᵢ[(pdim+1):(2pdim)]
            for j = 1:n_params
                Sⱼ = @view S[:, :, j]
                sⱼ_ref = Sⱼ * ψᵢ
                sⱼ = 𝒮_ket.ket_sens_results[k][(pdim*(i-1)+1):(pdim*i), j]
                @test norm(sⱼ - sⱼ_ref) / max(norm(sⱼ_ref), 1e-12) < 1e-6
            end
        end
    end
end

@testitem "E1: SplineIntegrator{MultiKetTrajectory} exact Hessian construction lanes (#344 shapes)" begin
    using DirectTrajOpt
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra

    H_drift = diagm(ComplexF64[0.0, 1.0, 2.0])
    H_x = ComplexF64[0 1 0; 1 0 1; 0 1 0] / sqrt(2)
    sys = QuantumSystem(H_drift, [H_x], [1.0])

    T = 1.0
    N = 5

    ψ0 = ComplexF64[1.0, 0.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0, 0.0]
    initials = [ψ0, ψ1]
    goals = [ψ1, ψ0]

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.5, 1, N), times)
    ensemble_qtraj = MultiKetTrajectory(sys, pulse, initials, goals)

    # :pair_indexed — the legacy dense second-order shape: quadratic state
    # [Φ; S₁..Sₙ; T_active], built at construction.
    𝒮p = SplineIntegrator(
        ensemble_qtraj,
        N;
        exact_hessian = true,
        second_order_shape = :pair_indexed,
    )
    @test 𝒮p.exact_hessian
    @test !isnothing(𝒮p.hess_probs)
    @test !isnothing(𝒮p.hess_active_pairs)
    @test length(𝒮p.hess_probs) == N - 1
    # n_active is the strictly-upper triangle of parameter pairs
    n_p = 2 * 1 + 2
    @test length(𝒮p.hess_active_pairs) == n_p * (n_p - 1) ÷ 2

    # :auto on Tsit5 + exact_hessian resolves to :directional — the
    # forward-over-adjoint shape carries NO per-pair block at all.
    𝒮d = SplineIntegrator(ensemble_qtraj, N; exact_hessian = true)
    @test 𝒮d.exact_hessian
    @test isnothing(𝒮d.hess_probs)
    @test isnothing(𝒮d.hess_active_pairs)

    # Second-order state scaling: pair-indexed is QUADRATIC in the parameter
    # count (grows with n_active), directional is not built here.
    n_ctrl = 2
    H_x2 = ComplexF64[0 1 0; 1 0 1; 0 1 0] / sqrt(2)
    H_y = ComplexF64[0 -im 0; im 0 -im; 0 im 0] / sqrt(2)
    sys2 = QuantumSystem(H_drift, [H_x2, H_y], [1.0, 1.0])
    pulse2 = LinearSplinePulse(fill(0.5, n_ctrl, N), times)
    eq2 = MultiKetTrajectory(sys2, pulse2, initials, goals)
    𝒮p2 = SplineIntegrator(eq2, N; exact_hessian = true, second_order_shape = :pair_indexed)
    n_p2 = 2 * n_ctrl + 2
    @test length(𝒮p2.hess_active_pairs) == n_p2 * (n_p2 - 1) ÷ 2
    @test length(𝒮p2.hess_probs[1].u0) > length(𝒮p.hess_probs[1].u0)

    # ── Error lanes ──────────────────────────────────────────────────────────
    # Invalid shape keyword
    @test_throws ArgumentError SplineIntegrator(
        ensemble_qtraj,
        N;
        exact_hessian = true,
        second_order_shape = :bogus,
    )
    # :directional is a Tsit5 augmented-ODE cell only
    @test_throws ErrorException SplineIntegrator(
        ensemble_qtraj,
        N;
        exact_hessian = true,
        second_order_shape = :directional,
        alg = MagnusGL4Alg(),
    )
    # ChebyshevAlg owns its matrix-free curvature — no silent Tsit5 fallback
    cpulse = CubicSplinePulse(fill(0.5, 1, N), fill(0.0, 1, N), times)
    ceq = MultiKetTrajectory(sys, cpulse, initials, goals)
    @test_throws ErrorException SplineIntegrator(
        ceq,
        N;
        exact_hessian = true,
        alg = ChebyshevAlg(bracket = (-8.0, 8.0), n_sub = 8),
    )
end

@testitem "E1: SplineIntegrator MultiKetTrajectory with NonlinearDrive" begin
    using DirectTrajOpt
    using Piccolo
    using SparseArrays
    using LinearAlgebra

    # 1-qubit system with 1 control and 2 drives:
    #   H = u₁·σx + u₁²·σz
    drives = AbstractDrive[
        LinearDrive(sparse(ComplexF64.(PAULIS.X)), 1),
        NonlinearDrive(PAULIS.Z, u -> u[1]^2),
    ]
    drive_bounds = [1.0]

    sys = QuantumSystem(zeros(ComplexF64, 2, 2), drives, drive_bounds)

    @test has_nonlinear_drives(sys.H_drives)

    T = 1.0
    N = 5

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]

    initials = [ψ0, ψ1]
    goals = [ψ1, ψ0]

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.3, 1, N), times)
    ensemble_qtraj = MultiKetTrajectory(sys, pulse, initials, goals)

    𝒮 = SplineIntegrator(ensemble_qtraj, N)
    @test 𝒮 isa SplineIntegrator{MultiKetTrajectory}
    @test length(𝒮.x_names) == 2

    combined_traj = NamedTrajectory(ensemble_qtraj, N)
    δ = zeros(𝒮.dim)
    evaluate!(δ, 𝒮, combined_traj)
    @test !all(iszero, δ)

    # Per-ket residual parity against single-ket integrators on the same
    # nonlinear system (the shared-propagator contract under a NonlinearDrive)
    for i = 1:2
        q = KetTrajectory(sys, pulse, initials[i], goals[i])
        t = NamedTrajectory(q, N)
        𝒮ᵢ = SplineIntegrator(q, N)
        δᵢ = zeros(𝒮ᵢ.dim)
        evaluate!(δᵢ, 𝒮ᵢ, t)
        for k = 1:(N-1)
            rows_ens = ((i - 1) * 2 * sys.levels) .+ (1:(2*sys.levels)) .+ (k - 1) * 𝒮.x_dim
            rows_sep = (2 * sys.levels) .* (k - 1) .+ (1:(2*sys.levels))
            @test δ[rows_ens] ≈ δᵢ[rows_sep] atol = 1e-8
        end
    end
end

@testitem "E1: SplineIntegrator MultiKetTrajectory NonlinearDrive cubic splines" begin
    using DirectTrajOpt
    using Piccolo
    using SparseArrays
    using LinearAlgebra

    # 2-qubit system: H = H_drift + u₁·σx⊗I + u₂·I⊗σx + (u₁·u₂)·σz⊗σz
    H_drift = kron(PAULIS.Z, PAULIS.I)
    drives = AbstractDrive[
        LinearDrive(sparse(ComplexF64.(kron(PAULIS.X, PAULIS.I))), 1),
        LinearDrive(sparse(ComplexF64.(kron(PAULIS.I, PAULIS.X))), 2),
        NonlinearDrive(kron(PAULIS.Z, PAULIS.Z), u -> u[1] * u[2]),
    ]
    drive_bounds = [1.0, 1.0]

    sys = QuantumSystem(H_drift, drives, drive_bounds)

    T = 1.0
    N = 5

    ψ0 = ComplexF64[1.0, 0.0, 0.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0, 0.0, 0.0]

    initials = [ψ0, ψ1]
    goals = [ψ1, ψ0]

    times = collect(range(0.0, T, N))
    pulse = CubicSplinePulse(fill(0.3, 2, N), fill(0.0, 2, N), times)
    ensemble_qtraj = MultiKetTrajectory(sys, pulse, initials, goals)

    # Force standard n×n sensitivity; cubic forward parity per ket
    𝒮 = SplineIntegrator(ensemble_qtraj, N; ket_sensitivity = false)
    @test 𝒮 isa SplineIntegrator{MultiKetTrajectory,CubicSpline}
    @test !𝒮.use_ket_sensitivity

    combined_traj = NamedTrajectory(ensemble_qtraj, N)
    δ = zeros(𝒮.dim)
    evaluate!(δ, 𝒮, combined_traj)
    @test norm(δ, Inf) < 1e-4
end

@testitem "E1: SplineIntegrator{MultiKetTrajectory} fixed-step Tsit5 (adaptive=false)" begin
    using DirectTrajOpt
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra

    T = 1.0
    N = 5
    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    initials = [ψ0, ψ1]
    goals = [ψ1, ψ0]

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.5, 2, N), times)
    qtraj = MultiKetTrajectory(sys, pulse, initials, goals)

    # Issue #180: adaptive=false crashed at construction (complex Φ_structure
    # vs the Float64-typed Tsit5Data field).
    𝒮 = SplineIntegrator(qtraj, N; alg = Tsit5Alg(adaptive = false, ode_h = 0.01))
    @test 𝒮 isa SplineIntegrator{MultiKetTrajectory,LinearSpline}

    traj = NamedTrajectory(qtraj, N)

    # ode_h is honored (#180 "related gap"): a fine fixed step agrees with the
    # adaptive solve, while a single coarse step over the interval must differ.
    𝒮_adapt = SplineIntegrator(qtraj, N; alg = Tsit5Alg())
    𝒮_coarse = SplineIntegrator(qtraj, N; alg = Tsit5Alg(adaptive = false, ode_h = 1.0))
    δ_fine = zeros(𝒮.x_dim)
    δ_adapt = zeros(𝒮.x_dim)
    δ_coarse = zeros(𝒮.x_dim)
    𝒮(δ_fine, traj[1], traj[2], 1)
    𝒮_adapt(δ_adapt, traj[1], traj[2], 1)
    𝒮_coarse(δ_coarse, traj[1], traj[2], 1)
    @test δ_fine ≈ δ_adapt atol = 1e-5
    @test norm(δ_coarse - δ_adapt) > 1e-9
end

@testitem "E1: SplineIntegrator MultiKetTrajectory with globals: extract_globals and per-knot solves" begin
    using DirectTrajOpt
    using DirectTrajOpt: evaluate!
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators:
        has_global_dependence, compute_ode_jacobian!, extract_globals

    # System with a global parameter: H = δ·σz + u₁·σx
    δ_init = 0.01
    sys = QuantumSystem(GATES.Z, [GATES.X], [1.0]; global_params = (δ = δ_init,))

    T = 5.0
    N = 5

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    initials = [ψ0, ψ1]
    goals = [ψ1, ψ0]

    times = collect(range(0.0, T, N))
    pulse = CubicSplinePulse(fill(0.5, 1, N), fill(0.0, 1, N), times)
    qtraj = MultiKetTrajectory(sys, pulse, initials, goals)

    integrator = SplineIntegrator(qtraj, N; global_names = [:δ])

    @test has_global_dependence(integrator)
    @test integrator.global_dim == 1
    @test integrator.global_names == [:δ]

    traj = NamedTrajectory(qtraj, N)
    @test haskey(traj.global_components, :δ)

    # extract_globals concatenates the global values in global_names order
    g = extract_globals(integrator, traj)
    @test g ≈ [δ_init]

    # The per-knot sensitivity solves run with the globals threaded into the
    # packed ODE parameters (u_dim = control_dim + global_dim = 2)
    @test integrator.u_dim == 2
    for k = 1:(N-1)
        compute_ode_jacobian!(integrator, traj[k], traj[k+1], k, g)
        @test !iszero(get_propagator(integrator.prop_results[k], 2))
    end

    # Forward residual with globals is dynamically consistent
    δ = zeros(integrator.dim)
    evaluate!(δ, integrator, traj)
    @test norm(δ, Inf) < 1e-4
end

# ── E1 resume fill (#347): Magnus lanes, ctor fallback lanes, the standalone
# ── Jacobian structure, and the seam gates reachable from the open core. ──── #

@testitem "E1: MultiKet MagnusGL4: construction + forward parity vs Tsit5" begin
    using DirectTrajOpt
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators: compute_ode_jacobian!

    H_drift = diagm(ComplexF64[0.0, 1.0, 2.0])
    H_x = ComplexF64[0 1 0; 1 0 1; 0 1 0] / sqrt(2)
    H_y = ComplexF64[0 -im 0; im 0 -im; 0 im 0] / sqrt(2)
    sys = QuantumSystem(H_drift, [H_x, H_y], [1.0, 1.0])

    T = 1.0
    N = 5
    ψ0 = ComplexF64[1.0, 0.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0, 0.0]
    initials = [ψ0, ψ1]
    goals = [ψ1, ψ0]

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.4, 2, N), times)
    eq = MultiKetTrajectory(sys, pulse, initials, goals)
    traj = NamedTrajectory(eq, N)

    # MagnusGL4 multiket construction runs the shared _build_alg_data dispatch
    # (the MagnusGL4Data buffer lane). MagnusAdapt4 mirrors it.
    𝒮_m = SplineIntegrator(eq, traj; alg = MagnusGL4Alg(n_steps = 40))
    @test 𝒮_m isa SplineIntegrator{MultiKetTrajectory,LinearSpline}
    @test 𝒮_m.alg isa MagnusGL4Alg

    𝒮_a = SplineIntegrator(eq, traj; alg = MagnusAdapt4Alg())
    @test 𝒮_a.alg isa MagnusAdapt4Alg

    # Forward: the non-Tsit5 lane propagates via the real-isomorphic Magnus
    # propagator and applies it to each ket.
    𝒮_t = SplineIntegrator(eq, traj)
    δ_m = zeros(𝒮_m.dim)
    δ_t = zeros(𝒮_t.dim)
    evaluate!(δ_m, 𝒮_m, traj)
    evaluate!(δ_t, 𝒮_t, traj)
    @test norm(δ_m - δ_t, Inf) < 1e-6
    @test norm(δ_m, Inf) < 1e-5  # rollout-consistent

    # Magnus + explicit ket_sensitivity: the forward caches the complex Φ built
    # from the Magnus real-iso propagator's Re/Im blocks (the conversion lane),
    # and the per-knot ket-sensitivity solve still agrees with the Tsit5 cell's.
    𝒮_mk = SplineIntegrator(eq, traj; alg = MagnusGL4Alg(n_steps = 40), ket_sensitivity = true)
    @test 𝒮_mk.use_ket_sensitivity
    δ_mk = zeros(𝒮_mk.dim)
    evaluate!(δ_mk, 𝒮_mk, traj)
    @test δ_mk ≈ δ_m atol = 1e-12

    # After a forward pass the Φ_vec caches carry the complex propagator
    for k = 1:(traj.N-1)
        @test !iszero(𝒮_mk.prop_results[k].Φ_vec)
    end

    # The ket-sensitivity per-knot solve works under the Magnus alg too
    # (the sensitivity ODE is always Tsit5-solved regardless of the forward alg)
    compute_ode_jacobian!(𝒮_mk, traj[1], traj[2], 1, nothing)
    @test !isnothing(𝒮_mk.ket_sens_results[1])

    # MagnusGL4Alg buffers are per-knot (thread-safe copies were built)
    @test length(𝒮_m.alg_data.G_drift_copies) == traj.N - 1
end

@testitem "E1: MultiKet constructor fallback lanes: pulse/order inference + error gates" begin
    using DirectTrajOpt
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators:
        _SamplingMemberQTraj, multiket_sens_refresh_count, _mk_multiket_qtraj, spline_order

    H_drift = diagm(ComplexF64[0.0, 1.0, 2.0])
    H_x = ComplexF64[0 1 0; 1 0 1; 0 1 0] / sqrt(2)
    sys = QuantumSystem(H_drift, [H_x], [1.0])

    T = 1.0
    N = 5
    ψ0 = ComplexF64[1.0, 0.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0, 0.0]
    initials = [ψ0, ψ1]
    goals = [ψ1, ψ0]
    times = collect(range(0.0, T, N))

    # Non-spline pulse + no explicit order → defaults to linear (order 1)
    zo = ZeroOrderPulse(0.3 * fill(1.0, 1, N), times)
    eq_zo = MultiKetTrajectory(sys, zo, initials, goals)
    traj_zo = NamedTrajectory(eq_zo, N)
    𝒮 = SplineIntegrator(eq_zo, traj_zo)
    @test spline_order(𝒮) == 1
    δ = zeros(𝒮.dim)
    evaluate!(δ, 𝒮, traj_zo)
    @test norm(δ, Inf) < 1e-5

    # The sampling-member shim: length reports the ket count
    ket_names = state_names(eq_zo)
    shim = _mk_multiket_qtraj(sys, zo, ket_names, :u, traj_zo)
    @test shim isa _SamplingMemberQTraj
    @test length(shim) == 2
    @test state_names(shim) == ket_names
    @test drive_name(shim) == :u
    @test get_system(shim) === sys
    @test get_pulse(shim) === zo

    # Drive-free system: the explicit-drives ArgumentError fires from the
    # shim-form constructor (the drive check precedes everything else).
    sys_free = QuantumSystem(GATES.Z)
    shim_free = _mk_multiket_qtraj(sys_free, zo, ket_names, :u, traj_zo)
    @test_throws ArgumentError SplineIntegrator(shim_free, traj_zo)

    # The refresh-count accessor is loadable and non-decreasing
    c0 = multiket_sens_refresh_count()
    @test c0 >= 0
    @test multiket_sens_refresh_count() >= c0
end

@testitem "E1: MultiKet cubic without du bounds zero-fills the derivative seed" begin
    using DirectTrajOpt
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators: spline_order

    H_drift = diagm(ComplexF64[0.0, 1.0, 2.0])
    H_x = ComplexF64[0 1 0; 1 0 1; 0 1 0] / sqrt(2)
    sys = QuantumSystem(H_drift, [H_x], [1.0])

    T = 1.0
    N = 5
    times = collect(range(0.0, T, N))
    ψ0 = ComplexF64[1.0, 0.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0, 0.0]
    eq = MultiKetTrajectory(
        sys,
        CubicSplinePulse(fill(0.3, 1, N), fill(0.0, 1, N), times),
        [ψ0, ψ1],
        [ψ1, ψ0],
    )

    # Rebuild the combined traj WITHOUT du bounds: the ctor seeds du from zeros
    # (the no-bounds lane) instead of the upper bound.
    traj_full = NamedTrajectory(eq, N)
    @test haskey(traj_full.bounds, :du)
    boundfree = NamedTrajectory(traj_full; bounds = (u = 1.0,))
    @test !haskey(boundfree.bounds, :du)

    𝒮 = SplineIntegrator(eq, boundfree)
    @test spline_order(𝒮) == 3
    δ = zeros(𝒮.dim)
    evaluate!(δ, 𝒮, boundfree)
    @test norm(δ, Inf) < 1e-5

    # With du bounds present the seed is the upper bound instead (covered lane,
    # re-asserted for the pair)
    𝒮_b = SplineIntegrator(eq, traj_full)
    δ_b = zeros(𝒮_b.dim)
    evaluate!(δ_b, 𝒮_b, traj_full)
    @test norm(δ_b, Inf) < 1e-5
end

@testitem "E1: MultiKet standalone jacobian_structure: blocks, du/dθ, globals" begin
    using DirectTrajOpt
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using SparseArrays
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators: jacobian_structure

    H_drift = diagm(ComplexF64[0.0, 1.0, 2.0])
    H_x = ComplexF64[0 1 0; 1 0 1; 0 1 0] / sqrt(2)
    H_y = ComplexF64[0 -im 0; im 0 -im; 0 im 0] / sqrt(2)
    sys = QuantumSystem(H_drift, [H_x, H_y], [1.0, 1.0])

    T = 1.0
    N = 5
    times = collect(range(0.0, T, N))
    ψ0 = ComplexF64[1.0, 0.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0, 0.0]
    eq = MultiKetTrajectory(
        sys,
        CubicSplinePulse(fill(0.3, 2, N), fill(0.0, 2, N), times),
        [ψ0, ψ1],
        [ψ1, ψ0],
    )
    traj = NamedTrajectory(eq, N)

    x_names = state_names(eq)
    ketdim = sys.levels
    # Φ_structure is per-ket in the ISO (real-isomorphic) basis: 2·ketdim square
    Φ_structure = sparse(ones(2 * ketdim, 2 * ketdim))

    # Cubic, du in traj: total_x_dim × (2 z_dim + global_dim) with du columns
    J = jacobian_structure(MultiKetTrajectory, x_names, :u, ketdim, Φ_structure, 3, traj)
    @test size(J) == (2 * 2 * ketdim, 2 * traj.dim)
    # Each ket's state block couples with: its own x cols (Φ structure),
    # k+1 x cols (identity), and all parameter columns
    x_comps_1 = traj.components[x_names[1]]
    @test nnz(J[1:(2ketdim), x_comps_1]) > 0
    @test nnz(J[1:(2ketdim), traj.dim .+ x_comps_1]) == 2 * ketdim

    # Cubic, traj carries :dθ instead of :du: the fallback name picks it up.
    # The state comps are built under the trajectory's OWN state names (the
    # combining-tilde symbols from state_names) — the structure lookup keys on
    # x_names verbatim.
    comps = Pair{Symbol,Matrix{Float64}}[]
    for name in x_names
        push!(comps, name => 0.1 * randn(2 * ketdim, N))
    end
    push!(comps, :u => 0.1 * randn(2, N))
    push!(comps, :dθ => 0.1 * randn(2, N))
    push!(comps, :Δt => fill(T / (N - 1), 1, N))
    push!(comps, :t => reshape(collect(range(0.0, T, N)), 1, :))
    θ_traj = NamedTrajectory(
        (; comps...);
        controls = :u,
        timestep = :Δt,
        bounds = (u = 1.0,),
    )
    J_dθ = jacobian_structure(
        MultiKetTrajectory,
        x_names,
        :u,
        ketdim,
        Φ_structure,
        3,
        θ_traj;
        global_names = Symbol[],
    )
    @test size(J_dθ) == (2 * 2 * ketdim, 2 * θ_traj.dim)

    # Linear: no du columns; parameter cols are [u at k; u at k+1; Δt; t]
    lin_eq = MultiKetTrajectory(
        sys,
        LinearSplinePulse(fill(0.3, 2, N), times),
        [ψ0, ψ1],
        [ψ1, ψ0],
    )
    lin_traj = NamedTrajectory(lin_eq, N)
    J_lin = jacobian_structure(
        MultiKetTrajectory,
        x_names,
        :u,
        ketdim,
        Φ_structure,
        1,
        lin_traj;
        global_names = Symbol[],
    )
    @test size(J_lin) == (2 * 2 * ketdim, 2 * lin_traj.dim)

    # Globals widen the structure and add per-ket dense global columns
    gsys = QuantumSystem(H_drift, [H_x, H_y], [1.0, 1.0]; global_params = (δ = 0.1,))
    g_eq = MultiKetTrajectory(
        gsys,
        LinearSplinePulse(fill(0.3, 2, N), times),
        [ψ0, ψ1],
        [ψ1, ψ0],
    )
    g_traj = NamedTrajectory(g_eq, N)
    J_g = jacobian_structure(
        MultiKetTrajectory,
        x_names,
        :u,
        ketdim,
        Φ_structure,
        1,
        g_traj;
        global_names = [:δ],
    )
    @test size(J_g) == (2 * 2 * ketdim, 2 * g_traj.dim + 1)
    # The global column is dense over both kets' state rows
    @test all(J_g[:, 2 * g_traj.dim + 1] .== 1.0)
end

@testitem "E1: MultiKet seam gates: dense Jacobian/Hessian need the proprietary layout" begin
    using DirectTrajOpt
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra

    H_drift = diagm(ComplexF64[0.0, 1.0, 2.0])
    H_x = ComplexF64[0 1 0; 1 0 1; 0 1 0] / sqrt(2)
    sys = QuantumSystem(H_drift, [H_x], [1.0])
    N = 5
    times = collect(range(0.0, 1.0, N))
    ψ0 = ComplexF64[1.0, 0.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0, 0.0]
    eq = MultiKetTrajectory(sys, LinearSplinePulse(fill(0.3, 1, N), times), [ψ0, ψ1], [ψ1, ψ0])
    traj = NamedTrajectory(eq, N)

    𝒮 = SplineIntegrator(eq, traj)

    # eval_jacobian reads matrix_free_layout(𝒮, traj) — no open-core method
    @test_throws MethodError eval_jacobian(𝒮, traj)

    # eval_hessian_of_lagrangian probes _multiket_directional_hvp_probs first
    # — no open-core method
    μ = zeros(𝒮.dim)
    @test_throws MethodError eval_hessian_of_lagrangian(𝒮, traj, μ)

    # ChebyshevAlg on MultiKet: the matrix-free alg data builder is proprietary,
    # so construction itself is the gate (cubic pulse passes the cubic-only
    # refusal, then MethodErrors at _build_alg_data).
    ceq = MultiKetTrajectory(
        sys,
        CubicSplinePulse(fill(0.3, 1, N), fill(0.0, 1, N), times),
        [ψ0, ψ1],
        [ψ1, ψ0],
    )
    @test_throws MethodError SplineIntegrator(
        ceq,
        traj;
        alg = ChebyshevAlg(bracket = (-8.0, 8.0), n_sub = 8),
    )
end
