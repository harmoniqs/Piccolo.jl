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

# ── E1 resume fill (#347): inner-ctor lanes, fixed-step probe, globals, and ─── #
# ── the standalone MultiDensity Jacobian structure. ────────────────────────── #

@testitem "E1: multidensity inner-ctor lanes: Chebyshev gate, pulse, du, fixed-step" begin
    using DirectTrajOpt
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators:
        spline_order, get_state_vectors

    L = ComplexF64[0 0.1; 0 0]
    sys = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    ρ0₁ = ComplexF64[1 0; 0 0]
    ρg₁ = ComplexF64[0 0; 0 1]
    ρ0₂ = ComplexF64[0 0; 0 1]
    ρg₂ = ComplexF64[1 0; 0 0]
    N = 5
    times = collect(range(0.0, 1.0, N))

    # The (qtraj, traj) inner form carries its own Magnus/Chebyshev refusals
    qtraj = MultiDensityTrajectory(
        sys,
        LinearSplinePulse(fill(0.3, 1, N), times),
        [ρ0₁, ρ0₂],
        [ρg₁, ρg₂],
    )
    traj = NamedTrajectory(qtraj, N)
    @test_throws ErrorException SplineIntegrator(qtraj, traj; alg = MagnusGL4Alg())
    @test_throws ErrorException SplineIntegrator(
        qtraj,
        traj;
        alg = ChebyshevAlg(bracket = (-8.0, 8.0)),
    )

    # Non-spline pulse defaults to linear
    zo_qtraj = MultiDensityTrajectory(
        sys,
        ZeroOrderPulse(0.3 * fill(1.0, 1, N), times),
        [ρ0₁, ρ0₂],
        [ρg₁, ρg₂],
    )
    zo_traj = NamedTrajectory(zo_qtraj, N)
    𝒮 = SplineIntegrator(zo_qtraj, zo_traj)
    @test spline_order(𝒮) == 1
    δ = zeros(𝒮.dim)
    evaluate!(δ, 𝒮, zo_traj)
    @test norm(δ, Inf) < 1e-4

    # Cubic without du bounds: the derivative seed zero-fills
    cq = MultiDensityTrajectory(
        sys,
        CubicSplinePulse(fill(0.3, 1, N), fill(0.0, 1, N), times),
        [ρ0₁, ρ0₂],
        [ρg₁, ρg₂],
    )
    ctraj = NamedTrajectory(cq, N)
    boundfree = NamedTrajectory(ctraj; bounds = (u = 1.0,))
    𝒮_c = SplineIntegrator(cq, boundfree)
    @test spline_order(𝒮_c) == 3
    δ_c = zeros(𝒮_c.dim)
    evaluate!(δ_c, 𝒮_c, boundfree)
    @test norm(δ_c, Inf) < 1e-4

    # Fixed-step Tsit5: the multidensity Φ-probe sparsity lane
    𝒮_f = SplineIntegrator(qtraj, traj; alg = Tsit5Alg(adaptive = false))
    δ_f = zeros(𝒮_f.dim)
    evaluate!(δ_f, 𝒮_f, traj)
    @test norm(δ_f, Inf) < 1e-4

    # Globals ride through the global-seed lane
    gsys = OpenQuantumSystem(
        PAULIS.Z,
        [PAULIS.X],
        [1.0];
        dissipation_operators = [L],
        global_params = (δ = 0.01,),
    )
    gq = MultiDensityTrajectory(
        gsys,
        LinearSplinePulse(fill(0.3, 1, N), times),
        [ρ0₁, ρ0₂],
        [ρg₁, ρg₂],
    )
    gtraj = NamedTrajectory(gq, N)
    # The inner form does not auto-detect globals — names passed explicitly
    𝒮_g = SplineIntegrator(gq, gtraj; global_names = [:δ])
    @test 𝒮_g.global_names == [:δ]
    @test 𝒮_g.global_dim == 1
    δ_g = zeros(𝒮_g.dim)
    evaluate!(δ_g, 𝒮_g, gtraj)
    @test norm(δ_g, Inf) < 1e-4

    # get_state_vectors returns per-density state vectors
    zₖ = traj[1]
    xs = get_state_vectors(𝒮, zₖ)
    @test length(xs) == 2
end

@testitem "E1: multidensity standalone jacobian_structure: blocks and globals" begin
    using DirectTrajOpt
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using SparseArrays
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators: jacobian_structure

    L = ComplexF64[0 0.1; 0 0]
    sys = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    ρ0₁ = ComplexF64[1 0; 0 0]
    ρg₁ = ComplexF64[0 0; 0 1]
    ρ0₂ = ComplexF64[0 0; 0 1]
    ρg₂ = ComplexF64[1 0; 0 0]
    N = 5
    times = collect(range(0.0, 1.0, N))
    qtraj = MultiDensityTrajectory(
        sys,
        CubicSplinePulse(fill(0.3, 1, N), fill(0.0, 1, N), times),
        [ρ0₁, ρ0₂],
        [ρg₁, ρg₂],
    )
    traj = NamedTrajectory(qtraj, N)

    n = sys.levels
    Φ_structure = sparse(ones(n^2, n^2))

    # Cubic with du: total_x_dim = n_densities · n² rows, both-knot cols + params
    J = jacobian_structure(
        MultiDensityTrajectory,
        state_names(qtraj),
        :u,
        n,
        Φ_structure,
        3,
        traj;
        global_names = Symbol[],
    )
    @test size(J) == (2 * n^2, 2 * traj.dim)
    x_comps_1 = traj.components[state_names(qtraj)[1]]
    @test nnz(J[1:(n^2), x_comps_1]) > 0
    @test nnz(J[1:(n^2), traj.dim .+ x_comps_1]) == n^2  # identity block

    # Linear member: no du columns
    lin_qtraj = MultiDensityTrajectory(
        sys,
        LinearSplinePulse(fill(0.3, 1, N), times),
        [ρ0₁, ρ0₂],
        [ρg₁, ρg₂],
    )
    lin_traj = NamedTrajectory(lin_qtraj, N)
    J_lin = jacobian_structure(
        MultiDensityTrajectory,
        state_names(lin_qtraj),
        :u,
        n,
        Φ_structure,
        1,
        lin_traj;
        global_names = Symbol[],
    )
    @test size(J_lin) == (2 * n^2, 2 * lin_traj.dim)

    # Globals add dense per-density global columns
    gsys = OpenQuantumSystem(
        PAULIS.Z,
        [PAULIS.X],
        [1.0];
        dissipation_operators = [L],
        global_params = (δ = 0.01,),
    )
    gq = MultiDensityTrajectory(
        gsys,
        LinearSplinePulse(fill(0.3, 1, N), times),
        [ρ0₁, ρ0₂],
        [ρg₁, ρg₂],
    )
    g_traj = NamedTrajectory(gq, N)
    J_g = jacobian_structure(
        MultiDensityTrajectory,
        state_names(gq),
        :u,
        n,
        Φ_structure,
        1,
        g_traj;
        global_names = [:δ],
    )
    @test size(J_g) == (2 * n^2, 2 * g_traj.dim + 1)
    @test all(J_g[:, 2*g_traj.dim+1] .== 1.0)
end
