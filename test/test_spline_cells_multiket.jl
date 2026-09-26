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
