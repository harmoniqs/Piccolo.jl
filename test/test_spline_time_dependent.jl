# #386: the SplineIntegrator family evaluates drive coefficients and the drift at a
# frozen t = 0, so explicitly time-dependent systems are refused at construction
# rather than silently integrated as if they were time-independent.

@testitem "SplineIntegrator refuses time-dependent systems (#386)" begin
    using Piccolo

    N = 6
    times = collect(range(0, 1.0, length = N))
    u = 0.3 .* ones(1, N)
    lin = LinearSplinePulse(u, times)
    cub = CubicSplinePulse(u, zeros(1, N), times)
    ψ0 = ComplexF64[1, 0]
    ψ1 = ComplexF64[0, 1]
    ρ0 = ComplexF64[1 0; 0 0]
    ρ1 = ComplexF64[0 0; 0 1]
    L = ComplexF64[0 0.1; 0 0]
    ω = 2π * 1.3

    modulated_drive = QuantumSystem(PAULIS.Z, [PAULIS.X => t -> cos(ω * t)], [1.0])
    modulated_drift = QuantumSystem(PAULIS.Z => t -> cos(ω * t), [PAULIS.X], [1.0])
    declared = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; time_dependent = true)
    ti = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])
    @test modulated_drive.time_dependent && modulated_drift.time_dependent

    closed_trajs(sys, pulse) = [
        KetTrajectory(sys, pulse, ψ0, ψ1),
        UnitaryTrajectory(sys, pulse, GATES[:X]),
        MultiKetTrajectory(sys, pulse, [ψ0, ψ1], [ψ1, ψ0]),
    ]

    refused(qtraj) =
        try
            SplineIntegrator(qtraj, N)
            ""
        catch e
            e isa ArgumentError ? e.msg : "not an ArgumentError: $(typeof(e))"
        end

    for sys in (modulated_drive, modulated_drift, declared), pulse in (lin, cub)
        for qtraj in [
            closed_trajs(sys, pulse)
            [SamplingTrajectory(q, [ti, sys]) for q in closed_trajs(ti, pulse)]
        ]
            msg = refused(qtraj)
            @test occursin("time-dependent", msg)
            @test occursin("t = 0", msg)
            @test occursin("BilinearIntegrator", msg)
        end
    end

    # The SplinePulseProblem default backend inherits the refusal.
    @test_throws "time-dependent" SplinePulseProblem(
        KetTrajectory(modulated_drive, lin, ψ0, ψ1),
        N,
    )

    # Open systems: an explicit flag, and a ModulatedDrive even without one.
    open_td = OpenQuantumSystem(
        PAULIS.Z,
        [PAULIS.X],
        [1.0];
        dissipation_operators = [L],
        time_dependent = true,
    )
    open_mod = OpenQuantumSystem(
        PAULIS.Z,
        AbstractDrive[ModulatedDrive(LinearDrive(PAULIS.X, 1), t -> cos(ω * t))],
        [1.0];
        dissipation_operators = [L],
    )
    for sys in (open_td, open_mod), pulse in (lin, cub)
        for qtraj in [
            DensityTrajectory(sys, pulse, ρ0, ρ1),
            MultiDensityTrajectory(sys, pulse, [ρ0, ρ1], [ρ1, ρ0]),
        ]
            msg = refused(qtraj)
            @test occursin("time-dependent", msg)
            @test occursin("BilinearIntegrator", msg)
        end
    end

    # Time-independent systems are unaffected.
    open_ti = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    for pulse in (lin, cub)
        for qtraj in closed_trajs(ti, pulse)
            @test SplineIntegrator(qtraj, N) !== nothing
        end
        @test SplineIntegrator(DensityTrajectory(open_ti, pulse, ρ0, ρ1), N) !== nothing
    end
end

@testitem "the recommended route for a modulated drive matches rollout (#386)" begin
    using Piccolo
    using DirectTrajOpt
    using NamedTrajectories
    using OrdinaryDiffEqTsit5: Tsit5

    # The refusal points at BilinearIntegrator(qtraj, N), which builds a
    # TimeDependentBilinearIntegrator that evaluates G(u, t) at absolute time
    # tₖ + τ·Δtₖ: on a modulated drive its collocation residual vanishes on a
    # Tsit5 1e-12 rollout, to its own (default) Tsit5 tolerance.
    N = 11
    T = 2.0
    times = collect(range(0, T, length = N))
    ω = 2π * 0.8
    sys = QuantumSystem(0.5 * PAULIS.Z, [PAULIS.X => t -> cos(ω * t)], [1.0])
    pulse = LinearSplinePulse(0.4 .* sin.(π .* times' ./ T), times)
    qtraj = KetTrajectory(
        sys,
        pulse,
        ComplexF64[1, 0],
        ComplexF64[0, 1];
        algorithm = Tsit5(),
        abstol = 1e-12,
        reltol = 1e-12,
    )

    @test_throws "time-dependent" SplineIntegrator(qtraj, N)
    integrator = BilinearIntegrator(qtraj, N)
    @test integrator isa TimeDependentBilinearIntegrator

    traj = NamedTrajectory(qtraj, N)
    δ = zeros(integrator.dim)
    DirectTrajOpt.evaluate!(δ, integrator, traj)
    @test maximum(abs, δ) < 1e-5
end

@testitem "function-built systems that read t warn unless declared (#386)" begin
    using Piccolo

    ω = 2π * 1.3
    L = ComplexF64[0 0.1; 0 0]
    warning = (:warn, r"reads `t`.*time_dependent = true"s)

    # `time_dependent` is the caller's declaration: an undeclared dependence is
    # detected and warned about, and the flag is left as given.
    reads_t = (u, t) -> 0.5 * PAULIS.Z + u[1] * cos(ω * t) * PAULIS.X
    # Time enters only through a control-weighted term: a probe at u = 0 (the
    # midpoint of symmetric bounds) alone would miss it.
    sys = @test_logs warning QuantumSystem(reads_t, [1.0])
    @test !sys.time_dependent
    @test_logs warning QuantumSystem(reads_t, [(-1.0, 1.0)])
    @test_logs warning QuantumSystem((u, t) -> cos(t) * PAULIS.Z + u[1] * PAULIS.X, [1.0])
    open_sys =
        @test_logs warning OpenQuantumSystem(reads_t, [1.0]; dissipation_operators = [L])
    @test !open_sys.time_dependent

    # Globals ride the probe vector: a global-scaled time dependence is seen too.
    with_global = (u, t) -> 0.5 * PAULIS.Z + u[1] * PAULIS.X + u[2] * sin(ω * t) * PAULIS.Y
    @test_logs warning QuantumSystem(with_global, [1.0]; global_params = (δ = 0.2,))

    # Declaring it, or a function that ignores t, is silent.
    @test (@test_logs QuantumSystem(reads_t, [1.0]; time_dependent = true)).time_dependent
    @test_logs OpenQuantumSystem(
        reads_t,
        [1.0];
        dissipation_operators = [L],
        time_dependent = true,
    )
    ignores_t = (u, t) -> 0.5 * PAULIS.Z + u[1] * PAULIS.X
    @test !(@test_logs QuantumSystem(ignores_t, [1.0])).time_dependent
    @test_logs OpenQuantumSystem(ignores_t, [1.0]; dissipation_operators = [L])

    # A drive-built OpenQuantumSystem with a ModulatedDrive is flagged like the
    # closed-system Pair constructor.
    @test OpenQuantumSystem(
        PAULIS.Z,
        AbstractDrive[ModulatedDrive(LinearDrive(PAULIS.X, 1), t -> cos(ω * t))],
        [1.0];
        dissipation_operators = [L],
    ).time_dependent
end
