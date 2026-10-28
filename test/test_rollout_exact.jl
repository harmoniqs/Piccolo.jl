# ============================================================================ #
# Rollout accuracy (#379)
#
# A zero-order-hold pulse rolls out with the exact piecewise exponential, so the
# rollout IS the exponential integrators' collocation propagation (to roundoff);
# spline pulses integrate the ODE with a default whose error tracks the requested
# tolerance.
# ============================================================================ #

@testitem "ZOH rollout matches the exponential integrators' collocation states (#379)" begin
    using LinearAlgebra
    using NamedTrajectories
    using DirectTrajOpt

    N, T = 11, 2.0
    times = collect(range(0.0, T, length = N))
    s = range(0.0, π, length = N)
    u = 0.4 .* vcat(sin.(s)', cos.(s)')
    zoh = ZeroOrderPulse(u, times)
    sys2 = QuantumSystem(0.5 * PAULIS.Z, [PAULIS.X, PAULIS.Y], [1.0, 1.0])
    a = annihilate(3)
    sys3 = QuantumSystem(-0.2 * a' * a' * a * a, [a + a', im * (a - a')], [1.0, 1.0])
    open_sys = OpenQuantumSystem(
        0.5 * PAULIS.Z,
        [PAULIS.X, PAULIS.Y],
        [1.0, 1.0];
        dissipation_operators = [ComplexF64[0 0.1; 0 0]],
    )
    e(n, i) = (v = zeros(ComplexF64, n); v[i] = 1; v)

    # Largest per-interval collocation defect of the integrator(s), evaluated on the
    # rollout's own states sampled at the knots: x_{k+1} - exp(G(u_k) Δt_k) x_k.
    function max_defect(qtraj, integ)
        traj = NamedTrajectory(qtraj, N)
        return maximum(integ isa AbstractVector ? integ : [integ]) do i
            δ = zeros(i.dim)
            DirectTrajOpt.evaluate!(δ, i, traj)
            maximum(abs, δ)
        end
    end

    fixtures = [
        "unitary 2-level" =>
            (UnitaryTrajectory(sys2, zoh, GATES[:X]), HermitianExponentialIntegrator),
        "unitary 3-level" => (
            UnitaryTrajectory(sys3, zoh, Matrix{ComplexF64}(I, 3, 3)),
            HermitianExponentialIntegrator,
        ),
        "ket 3-level" => (
            KetTrajectory(sys3, zoh, e(3, 1), e(3, 2)),
            HermitianExponentialIntegrator,
        ),
        "multi-ket" => (
            MultiKetTrajectory(sys3, zoh, [e(3, 1), e(3, 2)], [e(3, 2), e(3, 1)]),
            HermitianExponentialIntegrator,
        ),
        "density" => (
            DensityTrajectory(open_sys, zoh, ComplexF64[1 0; 0 0], ComplexF64[0 0; 0 1]),
            NonHermitianExponentialIntegrator,
        ),
        "multi-density" => (
            MultiDensityTrajectory(
                open_sys,
                zoh,
                [ComplexF64[1 0; 0 0], ComplexF64[0 0; 0 1]],
                [ComplexF64[0 0; 0 1], ComplexF64[1 0; 0 0]],
            ),
            NonHermitianExponentialIntegrator,
        ),
    ]
    for (label, (qtraj, Integrator)) in fixtures
        d = max_defect(qtraj, Integrator(qtraj, N))
        @test d ≤ 1e-10
        d ≤ 1e-10 || @info "ZOH collocation defect" label d
    end

    # Through `rollout!` too (the `sync_trajectory!` path) and `rollout`.
    qtraj = UnitaryTrajectory(sys3, zoh, Matrix{ComplexF64}(I, 3, 3))
    rollout!(qtraj, ZeroOrderPulse(0.5 .* u, times))
    @test max_defect(qtraj, HermitianExponentialIntegrator(qtraj, N)) ≤ 1e-10
    q2 = rollout(qtraj, zoh)
    @test max_defect(q2, HermitianExponentialIntegrator(q2, N)) ≤ 1e-10
end

@testitem "ZOH rollout: saveat / n_save preserved, partial-interval saves exact (#379)" begin
    using LinearAlgebra
    using OrdinaryDiffEqTsit5: Tsit5

    # 4 knots, 7 save points: most save times fall strictly inside a knot interval.
    times = [0.0, 0.7, 1.3, 2.0]
    u = [0.3 -0.5 0.8 0.1; 0.2 0.4 -0.6 0.0]
    zoh = ZeroOrderPulse(u, times)
    sys = QuantumSystem(0.5 * PAULIS.Z, [PAULIS.X, PAULIS.Y], [1.0, 1.0])
    open_sys = OpenQuantumSystem(
        0.5 * PAULIS.Z,
        [PAULIS.X, PAULIS.Y],
        [1.0, 1.0];
        dissipation_operators = [ComplexF64[0 0.2; 0 0]],
    )
    ref = (algorithm = Tsit5(), abstol = 1e-13, reltol = 1e-13, n_save = 7)
    save_times = collect(range(0.0, 2.0, length = 7))

    qU = UnitaryTrajectory(sys, zoh, GATES[:X]; n_save = 7)
    qU_ref = UnitaryTrajectory(sys, zoh, GATES[:X]; ref...)
    @test qU.solution.t == save_times
    @test length(qU.solution.u) == 7
    @test maximum(norm.(qU.solution.u .- qU_ref.solution.u)) ≤ 1e-10
    @test qU.solution[:U] == qU.solution.u
    @test all(U -> norm(U' * U - I) ≤ 1e-13, qU.solution.u)

    ψ0, ψg = ComplexF64[1, 0], ComplexF64[0, 1]
    qψ = KetTrajectory(sys, zoh, ψ0, ψg; n_save = 7)
    qψ_ref = KetTrajectory(sys, zoh, ψ0, ψg; ref...)
    @test maximum(norm.(qψ.solution.u .- qψ_ref.solution.u)) ≤ 1e-10

    qm = MultiKetTrajectory(sys, zoh, [ψ0, ψg], [ψg, ψ0]; n_save = 7)
    qm_ref = MultiKetTrajectory(sys, zoh, [ψ0, ψg], [ψg, ψ0]; ref...)
    @test length(qm.solution.u) == 2
    for k = 1:2
        @test qm.solution.u[k].t == save_times
        @test maximum(norm.(qm.solution.u[k].u .- qm_ref.solution.u[k].u)) ≤ 1e-10
    end
    @test fidelity(qm) ≈ fidelity(qm_ref) atol = 1e-12

    ρ0, ρg = ComplexF64[1 0; 0 0], ComplexF64[0 0; 0 1]
    qρ = DensityTrajectory(open_sys, zoh, ρ0, ρg; n_save = 7)
    qρ_ref = DensityTrajectory(open_sys, zoh, ρ0, ρg; ref...)
    @test qρ.solution.t == save_times
    @test maximum(norm.(qρ.solution.u .- qρ_ref.solution.u)) ≤ 1e-10
    @test fidelity(qρ) ≈ fidelity(qρ_ref) atol = 1e-12

    # A non-Hermitian closed system takes the general `exp` branch.
    nh = QuantumSystem(
        0.5 * PAULIS.Z - 0.05im * ComplexF64[0 0; 0 1],
        [PAULIS.X, PAULIS.Y],
        [1.0, 1.0];
        hermitian = false,
    )
    q_nh = KetTrajectory(nh, zoh, ψ0, ψg; n_save = 7)
    q_nh_ref = KetTrajectory(nh, zoh, ψ0, ψg; ref...)
    @test maximum(norm.(q_nh.solution.u .- q_nh_ref.solution.u)) ≤ 1e-10

    # The default `n_save` and the callable are unchanged in shape.
    q101 = UnitaryTrajectory(sys, zoh, GATES[:X])
    @test length(q101.solution.u) == 101
    @test q101(1.0) isa Matrix{ComplexF64}
end

@testitem "ZOH exact rollout: algorithm selection and in-place type stability (#379)" begin
    using LinearAlgebra
    using OrdinaryDiffEqTsit5: Tsit5
    using OrdinaryDiffEqLinear: MagnusAdapt4

    times = collect(range(0.0, 1.0, length = 6))
    u = 0.3 .* ones(1, 6)
    sys = QuantumSystem(0.5 * PAULIS.Z, [PAULIS.X], [1.0])
    zoh = ZeroOrderPulse(u, times)
    lin = LinearSplinePulse(u, times)

    @test UnitaryTrajectory(sys, zoh, GATES[:X]).solution.alg isa PiecewiseExponential
    @test UnitaryTrajectory(sys, zoh, GATES[:X]; algorithm = Tsit5()).solution.alg isa Tsit5
    @test UnitaryTrajectory(sys, lin, GATES[:X]).solution.alg isa Tsit5
    @test UnitaryTrajectory(sys, zoh, GATES[:X]; algorithm = PiecewiseExponential()).solution.alg isa
          PiecewiseExponential
    @test_throws ArgumentError UnitaryTrajectory(
        sys,
        lin,
        GATES[:X];
        algorithm = PiecewiseExponential(),
    )

    cub = CubicSplinePulse(u, zeros(1, 6), times)
    @test UnitaryTrajectory(sys, cub, GATES[:X]).solution.alg isa Tsit5

    # A time-dependent system is not piecewise constant under ZOH: it integrates the ODE
    # with the system's default algorithm, as before #379.
    td = QuantumSystem(
        (u, t) -> 0.5 * PAULIS.Z + u[1] * cos(t) * PAULIS.X,
        [(-1.0, 1.0)];
        time_dependent = true,
    )
    @test UnitaryTrajectory(td, zoh, GATES[:X]).solution.alg isa MagnusAdapt4
    @test UnitaryTrajectory(td, lin, GATES[:X]).solution.alg isa MagnusAdapt4
    # A bare H(u, t) function's `time_dependent` flag is unchecked (default false even
    # when H reads t), so function-built systems keep the ODE path too.
    fn = QuantumSystem((u, t) -> 0.5 * PAULIS.Z + u[1] * cos(t) * PAULIS.X, [(-1.0, 1.0)])
    @test !fn.time_dependent
    @test UnitaryTrajectory(fn, zoh, GATES[:X]).solution.alg isa MagnusAdapt4

    # In-place `rollout!` with `algorithm = nothing` keeps the algorithm the trajectory
    # was built with, so the typed `solution` field accepts the new solve.
    for alg in (Tsit5(), MagnusAdapt4(), nothing)
        qtraj = UnitaryTrajectory(sys, zoh, GATES[:X]; algorithm = alg)
        stored = typeof(qtraj.solution.alg)
        rollout!(qtraj, ZeroOrderPulse(0.5 .* u, times))
        @test qtraj.solution.alg isa stored
        rollout!(qtraj; n_save = 11)
        @test length(qtraj.solution.u) == 11
    end
    open_sys = OpenQuantumSystem(
        0.5 * PAULIS.Z,
        [PAULIS.X],
        [1.0];
        dissipation_operators = [ComplexF64[0 0.1; 0 0]],
    )
    qρ = DensityTrajectory(
        open_sys,
        zoh,
        ComplexF64[1 0; 0 0],
        ComplexF64[0 0; 0 1];
        algorithm = Tsit5(),
    )
    rollout!(qρ, ZeroOrderPulse(0.5 .* u, times))
    @test qρ.solution.alg isa Tsit5
end

@testitem "NamedTrajectory :constant rollout helpers are exact (#379)" begin
    using LinearAlgebra
    using NamedTrajectories

    N, T = 8, 1.5
    times = collect(range(0.0, T, length = N))
    u = 0.4 .* vcat(sin.(range(0, π, length = N))', cos.(range(0, π, length = N))')
    sys = QuantumSystem(0.5 * PAULIS.Z, [PAULIS.X, PAULIS.Y], [1.0, 1.0])
    zoh = ZeroOrderPulse(u, times)
    qU = UnitaryTrajectory(sys, zoh, GATES[:X]; n_save = N)
    qψ = KetTrajectory(sys, zoh, ComplexF64[1, 0], ComplexF64[0, 1]; n_save = N)

    trajU = NamedTrajectory(qU, N)
    Ũ⃗ = unitary_rollout(trajU, sys; interpolation = :constant)
    @test maximum(norm.(eachcol(Ũ⃗ .- trajU.Ũ⃗))) ≤ 1e-12
    @test unitary_rollout_fidelity(trajU, sys; interpolation = :constant) ≈
          unitary_fidelity(qU.solution.u[end], GATES[:X]) atol = 1e-12

    trajψ = NamedTrajectory(qψ, N)
    ψ̃ = ket_rollout(trajψ, sys; interpolation = :constant)
    @test maximum(norm.(eachcol(ψ̃ .- trajψ.ψ̃))) ≤ 1e-12
    @test rollout_fidelity(trajψ, sys; interpolation = :constant) ≈ fidelity(qψ) atol =
        1e-12
end

@testitem "spline rollouts meet the requested tolerance vs a Tsit5 1e-14 reference (#379)" begin
    using LinearAlgebra
    using OrdinaryDiffEqTsit5: Tsit5

    N, T = 11, 2.0
    times = collect(range(0.0, T, length = N))
    s = range(0.0, π, length = N)
    u = 0.4 .* vcat(sin.(s)', cos.(s)')
    du = 0.3 .* vcat(cos.(s)', -sin.(s)')
    a = annihilate(3)
    systems = [
        QuantumSystem(0.5 * PAULIS.Z, [PAULIS.X, PAULIS.Y], [1.0, 1.0]),
        QuantumSystem(-0.1 * a' * a' * a * a, [a + a', im * (a - a')], [1.0, 1.0]),
    ]
    pulses = [LinearSplinePulse(u, times), CubicSplinePulse(u, du, times)]
    for sys in systems, p in pulses
        Id = Matrix{ComplexF64}(I, sys.levels, sys.levels)
        Uref = UnitaryTrajectory(
            sys,
            p,
            Id;
            algorithm = Tsit5(),
            abstol = 1e-14,
            reltol = 1e-14,
        ).solution.u[end]
        for tol in (1e-8, 1e-10)
            U = UnitaryTrajectory(sys, p, Id; abstol = tol, reltol = tol).solution.u[end]
            @test norm(U - Uref) ≤ tol
        end
    end
end
