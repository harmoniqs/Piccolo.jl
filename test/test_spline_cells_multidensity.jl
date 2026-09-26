# ============================================================================ #
# Cluster E1 (#347) — SplineIntegrator MultiDensityTrajectory cell behavior
# tests.
#
# Ported from the pre-split embedded suite (Piccolissimo
# src/integrators/spline/spline_integrator_multidensity.jl). The whole cell sat
# at 0% — the multi-density constructor, call operator, Jacobian assembly and
# Hessian assembly had no tests in the open-core repo at all.
# ============================================================================ #

@testitem "E1: testing SplineIntegrator{MultiDensityTrajectory} with LinearSpline" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra

    # Simple 2-level open system with decay
    L = ComplexF64[0.1 0.0; 0.0 0.0]
    sys = OpenQuantumSystem(
        PAULIS.Z,
        [PAULIS.X, PAULIS.Y],
        [1.0, 1.0];
        dissipation_operators = [L],
    )

    T = 1.0
    N = 5

    ρ0s = [ComplexF64[1.0 0.0; 0.0 0.0], ComplexF64[0.0 0.0; 0.0 1.0]]
    ρgs = [ComplexF64[0.0 0.0; 0.0 1.0], ComplexF64[1.0 0.0; 0.0 0.0]]

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.5, 2, N), times)
    qtraj = MultiDensityTrajectory(sys, pulse, ρ0s, ρgs)

    𝒮 = SplineIntegrator(qtraj, N)
    @test 𝒮 isa SplineIntegrator{MultiDensityTrajectory,LinearSpline}
    @test 𝒮.ketdim == 2
    @test 𝒮.x_dim == 8  # n² * K = 4 * 2

    traj = NamedTrajectory(qtraj, N)
    test_integrator(𝒮, traj; atol = 1e-4, gauss_newton = true)
end

@testitem "E1: testing SplineIntegrator{MultiDensityTrajectory} with CubicSpline" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra

    L = ComplexF64[0.1 0.0; 0.0 0.0]
    sys = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])

    T = 1.0
    N = 5

    ρ0s = [ComplexF64[1.0 0.0; 0.0 0.0], ComplexF64[0.0 0.0; 0.0 1.0]]
    ρgs = [ComplexF64[0.0 0.0; 0.0 1.0], ComplexF64[1.0 0.0; 0.0 0.0]]

    times = collect(range(0.0, T, N))
    pulse = CubicSplinePulse(fill(0.5, 1, N), fill(0.0, 1, N), times)
    qtraj = MultiDensityTrajectory(sys, pulse, ρ0s, ρgs)

    𝒮 = SplineIntegrator(qtraj, N)
    @test 𝒮 isa SplineIntegrator{MultiDensityTrajectory,CubicSpline}

    traj = NamedTrajectory(qtraj, N)
    test_integrator(𝒮, traj; atol = 1e-4, gauss_newton = true)
end

@testitem "E1: SplineIntegrator{MultiDensityTrajectory} constraint satisfaction" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra

    # 3-level system to test non-trivial dimensions (n²=9, K=2 → x_dim=18)
    H_drift = diagm(ComplexF64[1.0, 0.0, -1.0])
    H_drive = zeros(ComplexF64, 3, 3)
    H_drive[1, 2] = H_drive[2, 1] = 1.0
    L = zeros(ComplexF64, 3, 3)
    L[1, 2] = 0.05  # |2⟩ → |1⟩ decay
    sys = OpenQuantumSystem(H_drift, [H_drive], [1.0]; dissipation_operators = [L])

    T = 1.0
    N = 5

    ρ0_1 = zeros(ComplexF64, 3, 3)
    ρ0_1[1, 1] = 1.0
    ρ0_2 = zeros(ComplexF64, 3, 3)
    ρ0_2[2, 2] = 1.0
    ρg_1 = zeros(ComplexF64, 3, 3)
    ρg_1[2, 2] = 1.0
    ρg_2 = zeros(ComplexF64, 3, 3)
    ρg_2[1, 1] = 1.0

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.3, 1, N), times)
    qtraj = MultiDensityTrajectory(sys, pulse, [ρ0_1, ρ0_2], [ρg_1, ρg_2])

    𝒮 = SplineIntegrator(qtraj, N)
    @test 𝒮 isa SplineIntegrator{MultiDensityTrajectory,LinearSpline}
    @test 𝒮.ketdim == 3
    @test 𝒮.x_dim == 18  # n² * K = 9 * 2

    # Test that constraint is satisfied for the rollout trajectory
    traj = NamedTrajectory(qtraj, N)
    δ = zeros(𝒮.dim)
    DirectTrajOpt.evaluate!(δ, 𝒮, traj)
    @test norm(δ, Inf) < 1e-2
end

@testitem "E1: SplineIntegrator{MultiDensityTrajectory} Jacobian with 3-level system" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra

    H_drift = diagm(ComplexF64[1.0, 0.0, -1.0])
    H_drive = zeros(ComplexF64, 3, 3)
    H_drive[1, 2] = H_drive[2, 1] = 1.0
    L = zeros(ComplexF64, 3, 3)
    L[1, 2] = 0.05
    sys = OpenQuantumSystem(H_drift, [H_drive], [1.0]; dissipation_operators = [L])

    T = 1.0
    N = 5

    ρ0_1 = zeros(ComplexF64, 3, 3)
    ρ0_1[1, 1] = 1.0
    ρ0_2 = zeros(ComplexF64, 3, 3)
    ρ0_2[2, 2] = 1.0
    ρg_1 = zeros(ComplexF64, 3, 3)
    ρg_1[2, 2] = 1.0
    ρg_2 = zeros(ComplexF64, 3, 3)
    ρg_2[1, 1] = 1.0

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.3, 1, N), times)
    qtraj = MultiDensityTrajectory(sys, pulse, [ρ0_1, ρ0_2], [ρg_1, ρg_2])

    𝒮 = SplineIntegrator(qtraj, N)
    traj = NamedTrajectory(qtraj, N)
    test_integrator(𝒮, traj; atol = 1e-4, gauss_newton = true)
end

@testitem "E1: SplineIntegrator{MultiDensityTrajectory} NonlinearDrive linear spline" begin
    using DirectTrajOpt
    using Piccolo
    using SparseArrays
    using LinearAlgebra

    # 2-level open system: H = σz + u₁·σx + u₁²·σz  with decay
    L = ComplexF64[0.1 0.0; 0.0 0.0]
    drives = AbstractDrive[
        LinearDrive(sparse(ComplexF64.(PAULIS.X)), 1),
        NonlinearDrive(PAULIS.Z, u -> u[1]^2),
    ]
    sys = OpenQuantumSystem(PAULIS.Z, drives, [1.0]; dissipation_operators = [L])

    T = 1.0
    N = 5

    ρ0_1 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρ0_2 = ComplexF64[0.0 0.0; 0.0 1.0]
    ρg_1 = ComplexF64[0.0 0.0; 0.0 1.0]
    ρg_2 = ComplexF64[1.0 0.0; 0.0 0.0]

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.3, 1, N), times)
    qtraj = MultiDensityTrajectory(sys, pulse, [ρ0_1, ρ0_2], [ρg_1, ρg_2])

    𝒮 = SplineIntegrator(qtraj, N)
    @test 𝒮 isa SplineIntegrator{MultiDensityTrajectory,LinearSpline}

    traj = NamedTrajectory(qtraj, N)
    test_integrator(𝒮, traj; atol = 1e-4, gauss_newton = true)
end

@testitem "E1: SplineIntegrator{MultiDensityTrajectory} NonlinearDrive cubic spline" begin
    using DirectTrajOpt
    using Piccolo
    using SparseArrays
    using LinearAlgebra

    # 2-level open system with bilinear NonlinearDrive:
    # H = σz + u₁·σx + u₂·σy + (u₁·u₂)·σz  with decay
    L = ComplexF64[0.1 0.0; 0.0 0.0]
    drives = AbstractDrive[
        LinearDrive(sparse(ComplexF64.(PAULIS.X)), 1),
        LinearDrive(sparse(ComplexF64.(PAULIS.Y)), 2),
        NonlinearDrive(PAULIS.Z, u -> u[1] * u[2]),
    ]
    sys = OpenQuantumSystem(PAULIS.Z, drives, [1.0, 1.0]; dissipation_operators = [L])

    T = 1.0
    N = 5

    ρ0_1 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρ0_2 = ComplexF64[0.0 0.0; 0.0 1.0]
    ρg_1 = ComplexF64[0.0 0.0; 0.0 1.0]
    ρg_2 = ComplexF64[1.0 0.0; 0.0 0.0]

    times = collect(range(0.0, T, N))
    pulse = CubicSplinePulse(fill(0.3, 2, N), fill(0.0, 2, N), times)
    qtraj = MultiDensityTrajectory(sys, pulse, [ρ0_1, ρ0_2], [ρg_1, ρg_2])

    𝒮 = SplineIntegrator(qtraj, N)
    @test 𝒮 isa SplineIntegrator{MultiDensityTrajectory,CubicSpline}

    traj = NamedTrajectory(qtraj, N)
    test_integrator(𝒮, traj; atol = 1e-4, gauss_newton = true)
end

@testitem "E1: SplineIntegrator{MultiDensityTrajectory} refuses Magnus/Chebyshev (Lindblad lane)" begin
    using Piccolo
    using LinearAlgebra

    L = ComplexF64[0.1 0.0; 0.0 0.0]
    sys = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    ρ0 = ComplexF64[1 0; 0 0]
    ρg = ComplexF64[0 0; 0 1]
    N = 5
    times = collect(range(0.0, 1.0, N))
    pulse = LinearSplinePulse(fill(0.3, 1, N), times)
    qtraj = MultiDensityTrajectory(sys, pulse, [ρ0, ρg], [ρg, ρ0])

    # Both (qtraj, N) and the inner (qtraj, traj) constructor refuse loudly
    @test_throws ErrorException SplineIntegrator(qtraj, N; alg = MagnusGL4Alg())
    traj = NamedTrajectory(qtraj, N)
    @test_throws ErrorException SplineIntegrator(qtraj, traj; alg = MagnusGL4Alg())
    @test_throws ErrorException SplineIntegrator(
        qtraj,
        N;
        alg = ChebyshevAlg(bracket = (-8.0, 8.0), n_sub = 8),
    )
end
