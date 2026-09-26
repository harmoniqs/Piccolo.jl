# ============================================================================ #
# Cluster E1 (#347) — SplineIntegrator SamplingTrajectory dispatch tests.
#
# Ported from the pre-split embedded suite (Piccolissimo
# src/integrators/spline/spline_integrator_sampling.jl). The sampling
# constructor (one integrator per ensemble member, shared-propagator per
# member) sat at 0%: the SplinePulseProblem sampling template testitem only
# wraps a SamplingProblem around the BASE problem and never constructs the
# spline sampling cells.
#
# ADAPTATION (open-core seam): the MultiKet sampling members are
# SplineIntegrator{MultiKetTrajectory} cells whose assembled-Jacobian path
# runs through the proprietary `matrix_free_layout` hook — so for THOSE the
# assertions are construction + forward (evaluate!); the Ket/Unitary/Density/
# MultiDensity members keep the full test_integrator conformance.
# ============================================================================ #

@testitem "E1: SplineIntegrator dispatch on SamplingTrajectory (Ket)" begin
    using DirectTrajOpt
    using NamedTrajectories
    using Piccolo

    sys1 = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    sys2 = QuantumSystem(1.1 * GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]

    N = 11
    times = collect(range(0, 1.0, length = N))
    pulse = LinearSplinePulse(zeros(2, N), times)

    base_qtraj = KetTrajectory(sys1, pulse, ψ_init, ψ_goal)
    sampling_qtraj = SamplingTrajectory(base_qtraj, [sys1, sys2])
    expanded_traj = NamedTrajectory(sampling_qtraj, N)

    integrators = SplineIntegrator(sampling_qtraj, N)

    # One integrator per ensemble member, each owning its member state
    @test integrators isa Vector{<:SplineIntegrator}
    @test length(integrators) == 2
    @test [𝒮.x_names for 𝒮 in integrators] == [[:ψ̃1], [:ψ̃2]]

    for 𝒮 in integrators
        test_integrator(𝒮, expanded_traj; atol = 1e-3, gauss_newton = true)
    end
end

@testitem "E1: SplineIntegrator dispatch on SamplingTrajectory (Unitary)" begin
    using DirectTrajOpt
    using NamedTrajectories
    using Piccolo

    sys1 = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    sys2 = QuantumSystem(1.1 * GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

    N = 11
    times = collect(range(0, 1.0, length = N))
    pulse = LinearSplinePulse(zeros(2, N), times)

    base_qtraj = UnitaryTrajectory(sys1, pulse, GATES[:H])
    sampling_qtraj = SamplingTrajectory(base_qtraj, [sys1, sys2])
    expanded_traj = NamedTrajectory(sampling_qtraj, N)

    integrators = SplineIntegrator(sampling_qtraj, N)

    @test integrators isa Vector{<:SplineIntegrator}
    @test length(integrators) == 2
    @test [𝒮.x_names for 𝒮 in integrators] == [[:Ũ⃗1], [:Ũ⃗2]]

    for 𝒮 in integrators
        test_integrator(𝒮, expanded_traj; atol = 1e-3, gauss_newton = true)
    end
end

@testitem "E1: SplineIntegrator dispatch on SamplingTrajectory (MultiKet)" begin
    using DirectTrajOpt
    using NamedTrajectories
    using Piccolo

    sys1 = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    sys2 = QuantumSystem(1.1 * GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]

    N = 11
    times = collect(range(0, 1.0, length = N))
    pulse = LinearSplinePulse(zeros(2, N), times)

    base_qtraj = MultiKetTrajectory(sys1, pulse, [ψ0, ψ1], [ψ1, ψ0])
    sampling_qtraj = SamplingTrajectory(base_qtraj, [sys1, sys2])
    expanded_traj = NamedTrajectory(sampling_qtraj, N)

    integrators = SplineIntegrator(sampling_qtraj, N)

    # ONE integrator per member (shared-propagator multiket), not one per ket —
    # the deliberate departure from the bilinear sampling flattening.
    @test integrators isa Vector{<:SplineIntegrator}
    @test length(integrators) == 2
    @test [𝒮.x_names for 𝒮 in integrators] == [[:ψ̃1, :ψ̃2], [:ψ̃3, :ψ̃4]]

    # Forward-only assertions for the multiket members (their assembled-Jacobian
    # path goes through the proprietary matrix_free_layout seam; see the file
    # header). Each member's forward residual must be well-defined, and the
    # MEMBER-1 residual must match the non-sampling multiket cell on sys1.
    for 𝒮 in integrators
        F = zeros(𝒮.dim)
        evaluate!(F, 𝒮, expanded_traj)
        @test !all(iszero, F)
        F2 = zeros(𝒮.dim)
        evaluate!(F2, 𝒮, expanded_traj)
        @test F ≈ F2
    end

    base_integrator = SplineIntegrator(base_qtraj, N)
    base_traj = NamedTrajectory(base_qtraj, N)
    F_base = zeros(base_integrator.dim)
    evaluate!(F_base, base_integrator, base_traj)
    F_m1 = zeros(integrators[1].dim)
    evaluate!(F_m1, integrators[1], expanded_traj)
    @test F_m1 ≈ F_base atol = 1e-10
end

@testitem "E1: SplineIntegrator sampling Jacobian parity (Ket)" begin
    using DirectTrajOpt
    using NamedTrajectories
    using Piccolo

    sys = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    N = 11
    times = range(0, 1.0, length = N)
    pulse = LinearSplinePulse(zeros(2, N), times)

    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)
    direct_𝒮 = SplineIntegrator(qtraj, N)
    direct_traj = NamedTrajectory(qtraj, N)

    sampling_qtraj = SamplingTrajectory(qtraj, [sys])
    sampling_𝒮 = SplineIntegrator(sampling_qtraj, N)[1]
    expanded_traj = NamedTrajectory(sampling_qtraj, N)

    test_integrator(direct_𝒮, direct_traj; atol = 1e-3, gauss_newton = true)
    test_integrator(sampling_𝒮, expanded_traj; atol = 1e-3, gauss_newton = true)
end

@testitem "E1: SplineIntegrator sampling Jacobian parity (Unitary)" begin
    using DirectTrajOpt
    using NamedTrajectories
    using Piccolo

    sys = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    N = 11
    times = range(0, 1.0, length = N)
    pulse = LinearSplinePulse(zeros(2, N), times)

    qtraj = UnitaryTrajectory(sys, pulse, GATES[:H])
    direct_𝒮 = SplineIntegrator(qtraj, N)
    direct_traj = NamedTrajectory(qtraj, N)

    sampling_qtraj = SamplingTrajectory(qtraj, [sys])
    sampling_𝒮 = SplineIntegrator(sampling_qtraj, N)[1]
    expanded_traj = NamedTrajectory(sampling_qtraj, N)

    test_integrator(direct_𝒮, direct_traj; atol = 1e-3, gauss_newton = true)
    test_integrator(sampling_𝒮, expanded_traj; atol = 1e-3, gauss_newton = true)
end

@testitem "E1: SplineIntegrator sampling Jacobian parity (MultiKet: forward)" begin
    using DirectTrajOpt
    using NamedTrajectories
    using Piccolo

    sys = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    N = 11
    times = range(0, 1.0, length = N)
    pulse = LinearSplinePulse(zeros(2, N), times)

    qtraj = MultiKetTrajectory(sys, pulse, [ψ0, ψ1], [ψ1, ψ0])
    direct_𝒮 = SplineIntegrator(qtraj, N)
    direct_traj = NamedTrajectory(qtraj, N)

    sampling_qtraj = SamplingTrajectory(qtraj, [sys])
    sampling_𝒮 = SplineIntegrator(sampling_qtraj, N)[1]
    expanded_traj = NamedTrajectory(sampling_qtraj, N)

    # The single-member sampling cell IS the base cell (forward parity; the
    # assembled Jacobian goes through the proprietary multiket seam).
    F_direct = zeros(direct_𝒮.dim)
    evaluate!(F_direct, direct_𝒮, direct_traj)
    F_sampling = zeros(sampling_𝒮.dim)
    evaluate!(F_sampling, sampling_𝒮, expanded_traj)
    @test F_direct ≈ F_sampling atol = 1e-10
end

@testitem "E1: SplineIntegrator registered in specs layer (#401)" begin
    using Piccolo

    entry = Piccolo.lookup_integrator(:spline)
    @test entry ≢ nothing

    sys = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    N = 11
    times = range(0, 1.0, length = N)
    pulse = LinearSplinePulse(zeros(2, N), times)
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)
    sampling_qtraj = SamplingTrajectory(qtraj, [sys])

    𝒮s = entry.factory(sampling_qtraj, N)
    @test 𝒮s isa Vector{<:SplineIntegrator}
    @test length(𝒮s) == 1
end

@testitem "E1: SplineIntegrator sampling per-member buffer isolation (#401)" begin
    using DirectTrajOpt
    using NamedTrajectories
    using Piccolo

    sys1 = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    sys2 = QuantumSystem(1.05 * GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

    N = 11
    times = range(0, 1.0, length = N)
    pulse = LinearSplinePulse(zeros(2, N), times)
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]

    qtraj = KetTrajectory(sys1, pulse, ψ_init, ψ_goal)
    sampling_qtraj = SamplingTrajectory(qtraj, [sys1, sys2])
    expanded_traj = NamedTrajectory(sampling_qtraj, N)

    𝒮s = SplineIntegrator(sampling_qtraj, N)
    @test length(𝒮s) == 2

    F1 = zeros(𝒮s[1].x_dim * (N - 1))
    F2 = zeros(𝒮s[2].x_dim * (N - 1))

    evaluate!(F1, 𝒮s[1], expanded_traj)
    evaluate!(F2, 𝒮s[2], expanded_traj)
    F1_saved = copy(F1)
    F2_saved = copy(F2)

    evaluate!(F1, 𝒮s[1], expanded_traj)
    evaluate!(F2, 𝒮s[2], expanded_traj)
    @test F1 ≈ F1_saved

    evaluate!(F2, 𝒮s[2], expanded_traj)
    @test F1 ≈ F1_saved
    @test F2 ≈ F2_saved
end

@testitem "E1: SplineIntegrator dispatch on SamplingTrajectory (Density) (AC1)" begin
    using DirectTrajOpt, NamedTrajectories, Piccolo

    L = ComplexF64[0 0.1; 0 0]
    sys1 = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    sys2 =
        OpenQuantumSystem(0.95 * PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    ρ0 = ComplexF64[1 0; 0 0]
    ρg = ComplexF64[0 0; 0 1]
    N = 11
    times = collect(range(0, 1.0, length = N))
    pulse = LinearSplinePulse(randn(1, N) .* 0.1, times)
    base_qtraj = DensityTrajectory(sys1, pulse, ρ0, ρg)
    sampling_qtraj = SamplingTrajectory(base_qtraj, [sys1, sys2])
    expanded_traj = NamedTrajectory(sampling_qtraj, N)

    integrators = SplineIntegrator(sampling_qtraj, N)
    @test integrators isa Vector{<:SplineIntegrator}
    @test length(integrators) == 2
    @test [𝒮.x_names for 𝒮 in integrators] == [[:ρ⃗̃1], [:ρ⃗̃2]]

    for 𝒮 in integrators
        test_integrator(𝒮, expanded_traj; atol = 1e-3, gauss_newton = true)
    end
end

@testitem "E1: SplineIntegrator dispatch on SamplingTrajectory (MultiDensity) (AC2)" begin
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
    pulse = LinearSplinePulse(randn(1, N) .* 0.1, times)
    base_qtraj = MultiDensityTrajectory(sys1, pulse, [ρ0₁, ρ0₂], [ρg₁, ρg₂])
    sampling_qtraj = SamplingTrajectory(base_qtraj, [sys1, sys2])
    expanded_traj = NamedTrajectory(sampling_qtraj, N)

    integrators = SplineIntegrator(sampling_qtraj, N)
    # 2 systems × 2 densities → 2 integrators (shared-propagator per member)
    @test integrators isa Vector{<:SplineIntegrator}
    @test length(integrators) == 2
    @test [𝒮.x_names for 𝒮 in integrators] == [[:ρ⃗̃1, :ρ⃗̃2], [:ρ⃗̃3, :ρ⃗̃4]]
    for 𝒮 in integrators
        test_integrator(𝒮, expanded_traj; atol = 1e-3, gauss_newton = true)
    end
end

@testitem "E1: SplineIntegrator sampling Jacobian parity (Density) (AC3)" begin
    using DirectTrajOpt, NamedTrajectories, Piccolo

    L = ComplexF64[0 0.1; 0 0]
    sys = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    ρ0 = ComplexF64[1 0; 0 0]
    ρg = ComplexF64[0 0; 0 1]
    N = 11
    times = range(0, 1.0, length = N)
    pulse = LinearSplinePulse(randn(1, N) .* 0.1, times)

    qtraj = DensityTrajectory(sys, pulse, ρ0, ρg)
    direct_𝒮 = SplineIntegrator(qtraj, N)
    direct_traj = NamedTrajectory(qtraj, N)

    sampling_qtraj = SamplingTrajectory(qtraj, [sys])
    sampling_𝒮 = SplineIntegrator(sampling_qtraj, N)[1]
    expanded_traj = NamedTrajectory(sampling_qtraj, N)

    test_integrator(direct_𝒮, direct_traj; atol = 1e-3, gauss_newton = true)
    test_integrator(sampling_𝒮, expanded_traj; atol = 1e-3, gauss_newton = true)
end

@testitem "E1: SplineIntegrator sampling Jacobian parity (MultiDensity) (AC3)" begin
    using DirectTrajOpt, NamedTrajectories, Piccolo

    L = ComplexF64[0 0.1; 0 0]
    sys = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    ρ0₁ = ComplexF64[1 0; 0 0]
    ρg₁ = ComplexF64[0 0; 0 1]
    ρ0₂ = ComplexF64[0 0; 0 1]
    ρg₂ = ComplexF64[1 0; 0 0]
    N = 11
    times = range(0, 1.0, length = N)
    pulse = LinearSplinePulse(randn(1, N) .* 0.1, times)

    qtraj = MultiDensityTrajectory(sys, pulse, [ρ0₁, ρ0₂], [ρg₁, ρg₂])
    # Non-sampling constructor returns a single shared-propagator integrator
    direct_𝒮 = SplineIntegrator(qtraj, N)
    direct_traj = NamedTrajectory(qtraj, N)

    # Sampling constructor returns a Vector (one per member)
    sampling_qtraj = SamplingTrajectory(qtraj, [sys])
    sampling_𝒮s = SplineIntegrator(sampling_qtraj, N)
    expanded_traj = NamedTrajectory(sampling_qtraj, N)

    @test direct_𝒮 isa SplineIntegrator
    @test sampling_𝒮s isa AbstractVector
    @test length(sampling_𝒮s) == 1
    sampling_𝒮 = sampling_𝒮s[1]
    @test sampling_𝒮 isa SplineIntegrator
    test_integrator(direct_𝒮, direct_traj; atol = 1e-3, gauss_newton = true)
    test_integrator(sampling_𝒮, expanded_traj; atol = 1e-3, gauss_newton = true)
end

@testitem "E1: SamplingProblem short solve with SplineIntegrator (Ket)" begin
    using LinearAlgebra, Random, DirectTrajOpt, Piccolo

    Random.seed!(42)

    sys_nom = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    sys_pert = QuantumSystem(1.05 * GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

    T = 1.0
    N = 11
    times = range(0, T, length = N)
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    pulse = ZeroOrderPulse(0.1 * randn(2, N), times)
    qtraj = KetTrajectory(sys_nom, pulse, ψ_init, ψ_goal)

    prob = SmoothPulseProblem(qtraj, N; Q = 100.0)
    sampling_bilin = SamplingProblem(prob, [sys_nom, sys_pert])
    solve!(sampling_bilin; max_iter = 50, print_level = 0)

    prob2 = SmoothPulseProblem(qtraj, N; Q = 100.0)
    sampling_spl = SamplingProblem(
        prob2,
        [sys_nom, sys_pert];
        integrator = (sq, n) -> SplineIntegrator(sq, n),
    )
    solve!(sampling_spl; max_iter = 50, print_level = 0)
end

@testitem "E1: SamplingProblem with a DensityTrajectory base errors loudly at the extension point" begin
    using LinearAlgebra, Random, DirectTrajOpt, Piccolo

    Random.seed!(42)

    L = ComplexF64[0 0.1; 0 0]
    sys_nom = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    sys_pert =
        OpenQuantumSystem(0.95 * PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    ρ0 = ComplexF64[1 0; 0 0]
    ρg = ComplexF64[0 0; 0 1]
    N = 11
    T = 1.0
    times = range(0, T, length = N)
    pulse = ZeroOrderPulse(0.1 * randn(1, N), times)
    qtraj = DensityTrajectory(sys_nom, pulse, ρ0, ρg)

    prob = SmoothPulseProblem(qtraj, N; Q = 100.0)

    # Piccolo deliberately provides NO public density fidelity objective: the
    # density sampling objective is a downstream extension (Piccolissimo
    # registers it through the `sampling_state_objective` hook). The gate must
    # error LOUDLY — a null objective would silently solve the wrong problem.
    err = try
        SamplingProblem(
            prob,
            [sys_nom, sys_pert];
            integrator = (sq, n) -> SplineIntegrator(sq, n),
        )
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("sampling_state_objective", err.msg)
    @test occursin("downstream package", err.msg)
end

@testitem "E1: SplineIntegrator sampling auto-detects nominal-system globals" begin
    using DirectTrajOpt
    using NamedTrajectories
    using Piccolo
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators: has_global_dependence

    # The sampling conversion does not attach globals; the sampling constructor
    # re-attaches them from the NOMINAL system's global_params (members share
    # the names) — mirroring the non-sampling conversions.
    sys1 = QuantumSystem(GATES[:Z], [GATES[:X]], [1.0]; global_params = (δ = 0.01,))
    sys2 = QuantumSystem(1.05 * GATES[:Z], [GATES[:X]], [1.0]; global_params = (δ = 0.01,))

    N = 11
    times = collect(range(0, 1.0, length = N))
    pulse = LinearSplinePulse(zeros(1, N), times)
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]

    base_qtraj = KetTrajectory(sys1, pulse, ψ_init, ψ_goal)
    sampling_qtraj = SamplingTrajectory(base_qtraj, [sys1, sys2])
    expanded_traj = NamedTrajectory(sampling_qtraj, N)

    𝒮s = SplineIntegrator(sampling_qtraj, N)
    @test length(𝒮s) == 2
    for 𝒮 in 𝒮s
        @test has_global_dependence(𝒮)
        @test 𝒮.global_names == [:δ]
        @test 𝒮.global_dim == 1
        # u_dim counts controls PLUS the threaded global
        @test 𝒮.u_dim == 1 + 1
    end

    # The member cells run with the globals threaded through the packed ODE
    # parameters: the per-knot call operator takes them explicitly (the
    # expanded conversion trajectory carries no global_data by design — the
    # globals live on the integrator's parameter vector).
    g = [0.01]
    δ = zeros(𝒮s[1].x_dim)
    𝒮s[1](δ, expanded_traj[1], expanded_traj[2], 1, g)
    @test !all(iszero, δ)
    # The member-1 knot-1 residual matches the NON-sampling ket cell on the
    # nominal system with the same globals.
    base_𝒮 = SplineIntegrator(base_qtraj, N)
    base_traj = NamedTrajectory(base_qtraj, N)
    δ_base = zeros(base_𝒮.x_dim)
    base_𝒮(δ_base, base_traj[1], base_traj[2], 1, g)
    @test δ ≈ δ_base atol = 1e-10
end
