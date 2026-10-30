# ============================================================================ #
# Pulse basis ↔ integrator pairing (#375)
#
# The table-driven pairing test: every (pulse basis, integrator family) pair a
# problem constructor ACCEPTS must integrate the pulse exactly (the collocation
# defect against the pulse's own reference rollout is at ODE-tolerance level), and
# every pair it REFUSES must be refused with an `ArgumentError` at construction.
# ============================================================================ #

@testitem "basis pairing table: accepted pairs are exact, the rest are refused (#375)" begin
    using LinearAlgebra
    using NamedTrajectories
    using DirectTrajOpt
    using OrdinaryDiffEqTsit5: Tsit5

    N, T = 11, 2.0
    times = collect(range(0.0, T, length = N))
    s = range(0.0, π, length = N)
    u = 0.4 .* vcat(sin.(s)', cos.(s)')
    du = 0.3 .* vcat(cos.(s)', -sin.(s)')
    pulses = Dict(
        :zoh => ZeroOrderPulse(u, times),
        :linear => LinearSplinePulse(u, times),
        :cubic => CubicSplinePulse(u, du, times),
    )
    closed = QuantumSystem(0.5 * PAULIS.Z, [PAULIS.X, PAULIS.Y], [1.0, 1.0])
    open_sys = OpenQuantumSystem(
        0.5 * PAULIS.Z,
        [PAULIS.X, PAULIS.Y],
        [1.0, 1.0];
        dissipation_operators = [ComplexF64[0 0.1; 0 0]],
    )
    # The reference rollout: the trajectory's own ODE solve of the pulse, saved
    # exactly at the knots, at tolerances far below the 1e-8 bar. Tsit5, because
    # the closed-system default (MagnusAdapt4) does not reach 1e-8 here even at
    # 1e-12 tolerances (measured 1.3e-6 on the ZOH fixture, 2.1e-5 on the linear one).
    ref = (algorithm = Tsit5(), abstol = 1e-12, reltol = 1e-12, n_save = N)
    closed_qtraj(p) = UnitaryTrajectory(closed, p, GATES[:X]; ref...)
    open_qtraj(p) =
        DensityTrajectory(open_sys, p, ComplexF64[1 0; 0 0], ComplexF64[0 0; 0 1]; ref...)
    families = Dict(
        :hermitian_exponential =>
            (closed_qtraj, q -> HermitianExponentialIntegrator(q, N)),
        :nonhermitian_exponential =>
            (open_qtraj, q -> NonHermitianExponentialIntegrator(q, N)),
        :bilinear => (closed_qtraj, q -> BilinearIntegrator(q, N)),
        :spline =>
            (closed_qtraj, q -> SplineIntegrator(q, N; alg = Tsit5Alg(tol = 1e-12))),
    )
    template(p) = p isa ZeroOrderPulse ? SmoothPulseProblem : SplinePulseProblem

    # Collocation defect against the reference rollout at the fixed controls. Each
    # propagator is (near-)norm-preserving, so the summed per-interval defect
    # bounds the terminal-state error |x_N(collocation) - x_N(rollout)|.
    function defect(qcp)
        integ = only(filter(i -> !(i isa DerivativeIntegrator), qcp.prob.integrators))
        δ = zeros(integ.dim)
        DirectTrajOpt.evaluate!(δ, integ, get_trajectory(qcp))
        return sum(norm, Iterators.partition(δ, integ.dim ÷ (N - 1)))
    end

    accepted = Set{Tuple{Symbol,Symbol}}()
    for (pk, pulse) in pulses, (fk, (mk_qtraj, mk_integ)) in families
        qtraj = mk_qtraj(pulse)
        qcp = try
            template(pulse)(qtraj, N; integrator = mk_integ(qtraj))
        catch err
            # A refusal is a construction-time ArgumentError (from the template
            # guard or from the spline integrator itself) — or, for cubic + an
            # explicit BilinearIntegrator, the older H1 error that predates #375.
            @test err isa ArgumentError ||
                  (err isa ErrorException && occursin("H1", err.msg))
            # Where the refused pair can still be forced (declared approximation),
            # it is visibly NOT the pulse's physics: the guard is not cosmetic.
            if pk !== :zoh &&
               fk in (:hermitian_exponential, :bilinear) &&
               !(pk === :cubic && fk === :bilinear)
                forced = template(pulse)(
                    qtraj,
                    N;
                    integrator = mk_integ(qtraj),
                    allow_basis_approximation = true,
                )
                @test defect(forced) > 1e-3
            end
            continue
        end
        push!(accepted, (pk, fk))
        @test defect(qcp) ≤ 1e-8
    end

    @test accepted == Set([
        (:zoh, :hermitian_exponential),
        (:zoh, :nonhermitian_exponential),
        (:zoh, :bilinear),
        (:linear, :spline),
        (:cubic, :spline),
    ])
end

@testitem "interpolation_basis: pulses, integrators, trajectories, problems (#375)" begin
    N = 6
    times = collect(range(0.0, 1.0, length = N))
    sys = QuantumSystem(0.5 * PAULIS.Z, [PAULIS.X], [1.0])
    zoh = ZeroOrderPulse(zeros(1, N), times)
    lin = LinearSplinePulse(zeros(1, N), times)
    cub = CubicSplinePulse(zeros(1, N), zeros(1, N), times)
    gauss = GaussianPulse([0.1], [0.2], [0.5], 1.0)

    @test interpolation_basis(zoh) === :zoh
    @test interpolation_basis(lin) === :linear
    @test interpolation_basis(cub) === :cubic
    @test interpolation_basis(gauss) === :analytic

    qz = UnitaryTrajectory(sys, zoh, GATES[:X])
    ql = UnitaryTrajectory(sys, lin, GATES[:X])
    qc = UnitaryTrajectory(sys, cub, GATES[:X])
    @test interpolation_basis(qz) === :zoh
    @test interpolation_basis(SamplingTrajectory(ql, [sys, sys])) === :linear

    @test interpolation_basis(HermitianExponentialIntegrator(qz, N)) === :zoh
    @test interpolation_basis(BilinearIntegrator(qz, N)) === :zoh
    @test interpolation_basis(SplineIntegrator(ql, N)) === :linear
    @test interpolation_basis(SplineIntegrator(qc, N)) === :cubic
    traj = add_control_derivatives(NamedTrajectory(ql, N), 1)
    @test isnothing(interpolation_basis(DerivativeIntegrator(:u, :du, traj)))

    @test interpolation_basis(SmoothPulseProblem(qz, N)) === :zoh
    @test interpolation_basis(SplinePulseProblem(qc, N)) === :cubic
end

@testitem "template guards refuse a mismatched basis unless declared (#375)" begin
    N = 8
    times = collect(range(0.0, 2.0, length = N))
    sys = QuantumSystem(0.5 * PAULIS.Z, [PAULIS.X, PAULIS.Y], [1.0, 1.0])
    u = 0.2 * ones(2, N)
    zoh = UnitaryTrajectory(sys, ZeroOrderPulse(u, times), GATES[:X])
    lin = UnitaryTrajectory(sys, LinearSplinePulse(u, times), GATES[:X])
    cub = UnitaryTrajectory(sys, CubicSplinePulse(u, zeros(2, N), times), GATES[:X])
    refusal(f) =
        try
            f()
            nothing
        catch e
            e
        end

    # ZOH templates: a spline integrator (built as a declared approximation) is
    # still refused by the template — the template must be told too.
    spl_on_zoh =
        SplineIntegrator(zoh, N; spline_order = 1, allow_basis_approximation = true)
    for template in (SmoothPulseProblem, BangBangPulseProblem)
        e = refusal(() -> template(zoh, N; integrator = spl_on_zoh))
        @test e isa ArgumentError
        @test occursin(":zoh", e.msg) && occursin(":linear", e.msg)
        @test occursin("allow_basis_approximation", e.msg)
        @test template(
            zoh,
            N;
            integrator = spl_on_zoh,
            allow_basis_approximation = true,
        ) isa QuantumControlProblem
        # the correct (default) pairing is untouched
        @test template(zoh, N) isa QuantumControlProblem
    end

    # Spline template: exponential / bilinear integrators and `:pwc` are refused.
    for (qtraj, integ) in (
        (lin, HermitianExponentialIntegrator(lin, N)),
        (lin, BilinearIntegrator(lin, N)),
        (cub, HermitianExponentialIntegrator(cub, N)),
    )
        e = refusal(() -> SplinePulseProblem(qtraj, N; integrator = integ))
        @test e isa ArgumentError
        @test occursin("SplineIntegrator", e.msg)
        @test SplinePulseProblem(
            qtraj,
            N;
            integrator = integ,
            allow_basis_approximation = true,
        ) isa QuantumControlProblem
    end
    for qtraj in (lin, cub)
        @test refusal(() -> SplinePulseProblem(qtraj, N; integrator_type = :pwc)) isa
              ArgumentError
        @test SplinePulseProblem(
            qtraj,
            N;
            integrator_type = :pwc,
            allow_basis_approximation = true,
        ) isa QuantumControlProblem
        @test SplinePulseProblem(qtraj, N) isa QuantumControlProblem
    end

    # SamplingProblem guards a supplied integrator the same way.
    base = SplinePulseProblem(lin, N)
    sq = SamplingTrajectory(lin, [sys, sys])
    @test refusal(
        () -> SamplingProblem(
            base,
            [sys, sys];
            integrator = HermitianExponentialIntegrator(sq, N),
        ),
    ) isa ArgumentError

    # check_basis_pairing is the shared, exported entry point.
    @test isnothing(check_basis_pairing(:X, lin, [SplineIntegrator(lin, N)]))
    @test_throws ArgumentError check_basis_pairing(:X, lin, [BilinearIntegrator(lin, N)])
end

@testitem "SamplingProblem members inherit the base integrator type and order (#375)" begin
    using DirectTrajOpt
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators: spline_order

    N = 8
    times = collect(range(0.0, 2.0, length = N))
    sys = QuantumSystem(0.5 * PAULIS.Z, [PAULIS.X, PAULIS.Y], [1.0, 1.0])
    sys2 = QuantumSystem(0.55 * PAULIS.Z, [PAULIS.X, PAULIS.Y], [1.0, 1.0])
    u = 0.2 * ones(2, N)
    dynamics(p) = filter(i -> !(i isa DerivativeIntegrator), p.prob.integrators)

    for (pulse, order) in
        ((LinearSplinePulse(u, times), 1), (CubicSplinePulse(u, zeros(2, N), times), 3))
        base = SplinePulseProblem(UnitaryTrajectory(sys, pulse, GATES[:X]), N)
        members = dynamics(SamplingProblem(base, [sys, sys2]))
        @test length(members) == 2
        @test all(i -> i isa SplineIntegrator, members)
        @test all(i -> spline_order(i) == order, members)
        @test all(i -> interpolation_basis(i) === interpolation_basis(pulse), members)
    end

    # A ZOH base stays on the exact exponential tier (unchanged #334 default).
    zoh = UnitaryTrajectory(sys, ZeroOrderPulse(u, times), GATES[:X])
    members = dynamics(SamplingProblem(SmoothPulseProblem(zoh, N), [sys, sys2]))
    @test all(i -> i isa HermitianExponentialIntegrator, members)

    # A declared-PWC base stays PWC: the members inherit the declared approximation.
    lin = UnitaryTrajectory(sys, LinearSplinePulse(u, times), GATES[:X])
    pwc =
        SplinePulseProblem(lin, N; integrator_type = :pwc, allow_basis_approximation = true)
    members = dynamics(SamplingProblem(pwc, [sys, sys2]))
    @test all(i -> i isa BilinearIntegrator, members)
end

@testitem "rollout helpers default to the pulse's own basis (#375)" begin
    using DirectTrajOpt
    using NamedTrajectories
    using Logging
    using OrdinaryDiffEqTsit5: Tsit5

    N = 21
    T = 10.0
    times = collect(range(0.0, T, length = N))
    sys = QuantumSystem(0.1 * PAULIS.Z, [PAULIS.X, PAULIS.Y], [1.0, 1.0])
    pulse = ZeroOrderPulse(0.1 * randn(2, N), times)
    qcp = SmoothPulseProblem(UnitaryTrajectory(sys, pulse, GATES[:X]), N; Q = 100.0)
    solve!(qcp; max_iter = 30, print_level = 0)

    # The problem's fidelity — `fidelity(qcp)`, the rollout of the extracted
    # zero-order-hold pulse — evaluated at tolerances far below the 1e-8 bar.
    tight = (algorithm = Tsit5(), abstol = 1e-12, reltol = 1e-12)
    F = fidelity(UnitaryTrajectory(sys, get_pulse(qcp.qtraj), GATES[:X]; tight...))
    @test rollout_interpolation(qcp) === :constant
    # Problem-level default: the zero-order hold the problem was optimized in.
    @test abs(unitary_rollout_fidelity(qcp; tight...) - F) ≤ 1e-8
    # NamedTrajectory default: a :ddu chain is the SmoothPulseProblem (ZOH) signature.
    traj = get_trajectory(qcp)
    @test abs(unitary_rollout_fidelity(traj, sys; tight...) - F) ≤ 1e-8
    # The old default (:linear) rolls out a different waveform.
    @test abs(unitary_rollout_fidelity(traj, sys; interpolation = :linear, tight...) - F) >
          1e-6

    # :cubic would read the derivative-integrator :du as Hermite tangents.
    @test_throws ArgumentError unitary_rollout_fidelity(qcp; interpolation = :cubic)
    @test_throws ArgumentError unitary_rollout_fidelity(traj, sys; interpolation = :cubic)

    # Linear spline problem: the default is :linear; :cubic is refused (its :du is a
    # DerivativeIntegrator variable too), while a cubic spline problem rolls out cubic.
    lin = SplinePulseProblem(
        UnitaryTrajectory(sys, LinearSplinePulse(0.1 * randn(2, N), times), GATES[:X]),
        N,
    )
    @test rollout_interpolation(lin) === :linear
    @test_throws ArgumentError rollout_interpolation(lin, :cubic)
    cub = SplinePulseProblem(
        UnitaryTrajectory(
            sys,
            CubicSplinePulse(0.1 * randn(2, N), zeros(2, N), times),
            GATES[:X],
        ),
        N,
    )
    @test rollout_interpolation(cub) === :cubic
    @test 0 ≤ unitary_rollout_fidelity(cub) ≤ 1

    # A bare NamedTrajectory with no basis evidence keeps :linear, and says so.
    lin_traj = get_trajectory(lin)
    @test_logs (:warn, r"no basis tag") match_mode = :any unitary_rollout_fidelity(
        lin_traj,
        sys,
    )
    @test resolve_rollout_interpolation(lin_traj, :u, nothing) === :linear
end

@testitem "post-solve divergence check: fidelity gap, sampling and ensembles (#375)" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra
    using Logging

    N = 11
    times = collect(range(0.0, 5.0, length = N))
    sys = QuantumSystem(0.1 * PAULIS.Z, [PAULIS.X], [2.0])
    sys2 = QuantumSystem(0.12 * PAULIS.Z, [PAULIS.X], [2.0])
    ψ0, ψ1 = ComplexF64[1, 0], ComplexF64[0, 1]

    # Single ket, correct pairing: the gap is at ODE-tolerance level.
    qtraj = KetTrajectory(sys, ZeroOrderPulse(0.15 * ones(1, N), times), ψ0, ψ1)
    qcp = SmoothPulseProblem(qtraj, N)
    sync_trajectory!(qcp; check_divergence = false)
    @test rollout_fidelity_gap(qcp) < 1e-8

    # A collocation error the RELATIVE state check misses (ε < 5e-3) but the
    # ABSOLUTE fidelity check catches (ΔF > 1e-5). Rotate the terminal state a
    # little toward the goal, so only the collocation side moves.
    x = get_trajectory(qcp).ψ̃[:, end]
    ψ = ComplexF64[x[1]+im*x[3], x[2]+im*x[4]]
    θ = 1e-3
    c = ψ1' * ψ
    ψ′ = normalize(ψ + θ * (c / abs(c)) * ψ1)   # phase-aligned: ΔF is first order in θ
    get_trajectory(qcp).ψ̃[:, end] .= [real(ψ′); imag(ψ′)]
    @test rollout_divergence(qcp) < ROLLOUT_DIVERGENCE_RTOL[]
    @test rollout_fidelity_gap(qcp) > ROLLOUT_FIDELITY_ATOL[]
    @test_logs (:warn, r"F_collocation - F_rollout") match_mode = :any _ =
        Piccolo.Control.QuantumControlProblems._warn_on_rollout_divergence(qcp)
    @test_logs min_level = Logging.Warn Piccolo.Control.QuantumControlProblems._warn_on_rollout_divergence(
        qcp;
        atol = 1.0,
    )

    # SamplingProblem: both checks are now defined (they used to be `nothing`).
    # Identical members: the (nominal-tiled) initial collocation states are exact.
    sp = SamplingProblem(SmoothPulseProblem(qtraj, N), [sys, sys])
    sync_trajectory!(sp; check_divergence = false)
    ε = rollout_divergence(sp)
    @test ε isa Float64 && ε < 1e-6
    @test rollout_fidelity_gap(sp) < 1e-8
    get_trajectory(sp).ψ̃2[:, end] .= [0.0, 1.0, 0.0, 0.0]   # break ONE member
    @test rollout_divergence(sp) > ROLLOUT_DIVERGENCE_RTOL[]
    @test rollout_fidelity_gap(sp) > 0.1
    # Distinct members: each member is rolled out with ITS OWN system, so the
    # unsolved problem's nominal-tiled guess for the perturbed member shows up.
    sp2 = SamplingProblem(SmoothPulseProblem(qtraj, N), [sys, sys2])
    sync_trajectory!(sp2; check_divergence = false)
    @test rollout_divergence(sp2) > 1e-6
    @test rollout_fidelity_gap(sp2) > 1e-6

    # Ensemble (multi-ket) trajectories: the fidelity gap is defined too.
    mk = MultiKetTrajectory(
        sys,
        ZeroOrderPulse(0.15 * ones(1, N), times),
        [ψ0, ψ1],
        [ψ1, ψ0],
    )
    mqcp = SmoothPulseProblem(mk, N)
    sync_trajectory!(mqcp; check_divergence = false)
    @test rollout_divergence(mqcp) < 1e-6
    @test rollout_fidelity_gap(mqcp) < 1e-8
end
