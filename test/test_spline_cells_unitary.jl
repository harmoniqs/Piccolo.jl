# ============================================================================ #
# Cluster E1 (#347) — SplineIntegrator UnitaryTrajectory cell behavior tests.
#
# Ported from the pre-split embedded suite (Piccolissimo
# src/integrators/spline/spline_integrator_unitary.jl, the tests-only companion
# that outlived the 3b open-core move). The open-core split left the DENSE
# unitary cell in Piccolo but its tests behind, which is why this surface sat at
# 77.6% (forward-only coverage via the dense-parity harness and the
# SplinePulseProblem template testitems). These restore the analytic
# Jacobian/Hessian conformance (`test_integrator`), the MagnusGL4 lanes, the
# global-variable lanes and the fixed-step lane.
#
# Adapted for the Piccolo-only environment: `using Piccolissimo` references
# resolve through Piccolo's own re-export seam
# (`Piccolo.Control.QuantumIntegrators.SplineIntegrators`); no ported assertion
# touches the matrix-free/ChebyshevAlg surface that stays proprietary.
# ============================================================================ #

@testitem "E1: testing SplineIntegrator with UnitaryTrajectory" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra

    T = 1.0
    N = 5
    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])

    U_goal = GATES.X

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.5, 2, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)
    traj = NamedTrajectory(qtraj, N)

    𝒮 = SplineIntegrator(qtraj, N)
    @test 𝒮 isa SplineIntegrator{UnitaryTrajectory,LinearSpline}

    test_integrator(𝒮, traj; atol = 1e-4, gauss_newton = true, show_jacobian_diff = true)
end

@testitem "E1: testing SplineIntegrator cubic spline with UnitaryTrajectory" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra

    T = 1.0
    N = 5
    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])

    U_goal = GATES.X

    times = collect(range(0.0, T, N))
    pulse = CubicSplinePulse(fill(0.5, 2, N), fill(0.0, 2, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)
    traj = NamedTrajectory(qtraj, N)

    𝒮 = SplineIntegrator(qtraj, N)
    @test 𝒮 isa SplineIntegrator{UnitaryTrajectory,CubicSpline}

    test_integrator(𝒮, traj; atol = 1e-4, gauss_newton = true, show_jacobian_diff = true)
end

@testitem "E1: SplineIntegrator with UnitaryTrajectory and SplinePulseProblem" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra

    T = 1.0
    N = 10
    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    U_goal = GATES.H

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.5, 2, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)
    integrator = SplineIntegrator(qtraj, N; spline_order = 1)

    @test integrator isa SplineIntegrator{UnitaryTrajectory}

    qcp = SplinePulseProblem(qtraj, N; Q = 100.0, R = 1e-2, integrator = integrator)

    @test qcp isa QuantumControlProblem
    @test length(qcp.prob.integrators) == 2  # dynamics + derivative integrator on du
    @test qcp.prob.integrators[1] isa SplineIntegrator{UnitaryTrajectory}

    solve!(qcp; max_iter = 10, print_level = 0)
end

@testitem "E1: SplineIntegrator cubic splines with SplinePulseProblem" begin
    using DirectTrajOpt
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra

    T = 1.0
    N = 10
    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    U_goal = GATES.H

    times = collect(range(0.0, T, N))
    pulse = CubicSplinePulse(fill(0.5, 2, N), fill(0.0, 2, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)
    integrator = SplineIntegrator(qtraj, N; spline_order = 3)

    @test integrator isa SplineIntegrator{UnitaryTrajectory}
    @test integrator isa SplineIntegrator{UnitaryTrajectory,CubicSpline}

    qcp = SplinePulseProblem(qtraj, N; Q = 100.0, R = 1e-2, integrator = integrator)

    @test qcp isa QuantumControlProblem
    @test length(qcp.prob.integrators) == 1  # dynamics
    @test qcp.prob.integrators[1] isa SplineIntegrator{UnitaryTrajectory,CubicSpline}

    solve!(qcp; max_iter = 10, print_level = 0)

    traj = get_trajectory(qcp)
    @test traj isa NamedTrajectory
end

@testitem "E1: Cubic vs linear splines comparison" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators: spline_order

    T = 1.0
    N = 8
    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    U_goal = GATES.H

    times = collect(range(0.0, T, N))

    pulse_linear = LinearSplinePulse(fill(0.5, 2, N), times)
    qtraj_linear = UnitaryTrajectory(sys, pulse_linear, U_goal)
    integrator_linear = SplineIntegrator(qtraj_linear, N; spline_order = 1)

    @test integrator_linear isa SplineIntegrator{UnitaryTrajectory,LinearSpline}

    pulse_cubic = CubicSplinePulse(fill(0.5, 2, N), fill(0.0, 2, N), times)
    qtraj_cubic = UnitaryTrajectory(sys, pulse_cubic, U_goal)
    integrator_cubic = SplineIntegrator(qtraj_cubic, N; spline_order = 3)

    @test integrator_cubic isa SplineIntegrator{UnitaryTrajectory,CubicSpline}

    # Verify cubic uses more ODE parameters (4*u_dim vs 2*u_dim for controls)
    @test spline_order(integrator_cubic) == 3
    @test spline_order(integrator_linear) == 1

    # Both should work with SplinePulseProblem
    qcp_linear = SplinePulseProblem(
        qtraj_linear,
        N;
        Q = 100.0,
        R = 1e-2,
        integrator = integrator_linear,
    )
    qcp_cubic = SplinePulseProblem(
        qtraj_cubic,
        N;
        Q = 100.0,
        R = 1e-2,
        integrator = integrator_cubic,
    )

    @test qcp_linear isa QuantumControlProblem
    @test qcp_cubic isa QuantumControlProblem
end

@testitem "E1: SplineIntegrator auto-detects global variables" begin
    using Piccolo
    using LinearAlgebra
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators:
        has_global_dependence, global_variables

    T = 1.0
    N = 5

    sys_global =
        QuantumSystem(GATES.Z, [GATES.X], [1.0]; global_params = (δ = 0.5, Ω = 1.0))
    U_goal = GATES.X
    times = collect(range(0, T, N))
    pulse = CubicSplinePulse(fill(0.5, 1, N), fill(0.0, 1, N), times)
    qtraj_global = UnitaryTrajectory(sys_global, pulse, U_goal)

    # Auto-detect: should pick up δ and Ω from sys.global_params
    integrator_auto = SplineIntegrator(qtraj_global, N; spline_order = 3)

    @test has_global_dependence(integrator_auto)
    @test integrator_auto.global_dim == 2
    @test Set(integrator_auto.global_names) == Set([:δ, :Ω])

    # Manual override: specify only one global
    integrator_manual =
        SplineIntegrator(qtraj_global, N; spline_order = 3, global_names = [:δ])

    @test has_global_dependence(integrator_manual)
    @test integrator_manual.global_dim == 2
    @test integrator_manual.global_names == [:δ]

    # System without globals: no global dependence
    sys_no_global = QuantumSystem(GATES.Z, [GATES.X], [1.0])
    qtraj_no_global = UnitaryTrajectory(sys_no_global, pulse, U_goal)
    integrator_no_global = SplineIntegrator(qtraj_no_global, N; spline_order = 3)

    @test !has_global_dependence(integrator_no_global)
    @test integrator_no_global.global_dim == 0
    @test isempty(integrator_no_global.global_names)
    @test isempty(global_variables(integrator_no_global))
end

@testitem "E1: SplineIntegrator with global variables - UnitaryTrajectory" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators:
        has_global_dependence, global_variables

    T = 1.0
    N = 5

    sys = QuantumSystem(GATES.Z, [GATES.X], [1.0]; global_params = (g = 0.5,))
    U_goal = GATES.X
    times = collect(range(0, T, N))
    pulse = CubicSplinePulse(fill(0.5, 1, N), fill(0.0, 1, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)

    traj = NamedTrajectory(qtraj, N)

    @test haskey(traj.global_components, :g)
    @test traj.global_dim == 1

    integrator = SplineIntegrator(qtraj, N, spline_order = 3; global_names = [:g])

    traj = NamedTrajectory(qtraj, N)

    @test integrator.global_dim == 1
    @test integrator.global_names == [:g]
    @test has_global_dependence(integrator)
    @test global_variables(integrator) == [:g]

    test_integrator(
        integrator,
        traj;
        atol = 1e-4,
        gauss_newton = true,
        show_hessian_diff = true,
        show_jacobian_diff = true,
    )
end

@testitem "E1: SplineIntegrator with multiple global variables" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra

    T = 1.0
    N = 4

    sys = QuantumSystem(
        GATES.Z,
        [GATES.X, GATES.Y],
        [1.0, 1.0];
        global_params = (α = 1.0, β = 0.5),
    )
    U_goal = GATES.X
    times = collect(range(0, T, N))
    pulse = CubicSplinePulse(fill(0.5, 2, N), fill(0.0, 2, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)

    traj = NamedTrajectory(qtraj, N)

    @test traj.global_dim == 2

    integrator = SplineIntegrator(qtraj, N, spline_order = 3; global_names = [:α, :β])

    traj = NamedTrajectory(qtraj, N)

    @test integrator.global_dim == 2
    @test length(integrator.global_names) == 2

    test_integrator(integrator, traj; atol = 1e-4, gauss_newton = true)
end

@testitem "E1: SplineIntegrator backward compatibility - no globals" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators: has_global_dependence

    T = 1.0
    N = 5

    sys = QuantumSystem(GATES.Z, [GATES.X], [1.0])
    U_goal = GATES.X
    times = collect(range(0, T, N))
    pulse = CubicSplinePulse(fill(0.5, 1, N), fill(0.0, 1, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)

    traj = NamedTrajectory(qtraj, N)
    @test traj.global_dim == 0

    integrator = SplineIntegrator(qtraj, N, spline_order = 3)

    traj = NamedTrajectory(qtraj, N)

    @test integrator.global_dim == 0
    @test isempty(integrator.global_names)
    @test !has_global_dependence(integrator)

    test_integrator(integrator, traj; atol = 1e-4, gauss_newton = true)
end

@testitem "E1: SplinePulseProblem with SplineIntegrator global support" begin
    using NamedTrajectories
    using DirectTrajOpt
    using DirectTrajOpt: BoundsConstraint
    using Piccolo
    using LinearAlgebra
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators: has_global_dependence

    # System with global parameter: detuning δ
    T = 5.0
    N = 5

    δ_init = 0.01

    sys = QuantumSystem(GATES.Z, [GATES.X], [1.0]; global_params = (δ = δ_init,))
    U_goal = GATES.X

    times = collect(range(0, T, N))
    pulse = CubicSplinePulse(fill(0.5, 1, N), fill(0.0, 1, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)

    integrator = SplineIntegrator(qtraj, N; spline_order = 3, global_names = [:δ])

    δ_bound = 0.5  # Symmetric: [-0.5, 0.5]

    qcp = SplinePulseProblem(
        qtraj,
        N;
        Q = 100.0,
        R = 1e-2,
        integrator = integrator,
        global_bounds = Dict{Symbol,Union{Float64,Tuple{Float64,Float64}}}(:δ => δ_bound),
    )

    bounds_constraints =
        filter(c -> c isa BoundsConstraint && c.is_global, qcp.prob.constraints)
    @test length(bounds_constraints) == 1

    @test qcp isa QuantumControlProblem
    @test length(qcp.prob.integrators) >= 1

    spline_int = first(filter(i -> i isa SplineIntegrator, qcp.prob.integrators))
    @test has_global_dependence(spline_int)
    @test spline_int.global_dim == 1

    traj = get_trajectory(qcp)
    @test haskey(traj.global_components, :δ)
    @test traj.global_dim == 1

    test_integrator(
        spline_int,
        traj;
        atol = 1e-4,
        gauss_newton = true,
        show_jacobian_diff = true,
        show_hessian_diff = true,
    )
end

@testitem "E1: SplineIntegrator Jacobian sparsity for UnitaryTrajectory" begin
    using DirectTrajOpt
    using Piccolo
    using NamedTrajectories
    using SparseArrays
    using LinearAlgebra

    T = 1.0
    N = 10
    sys = QuantumSystem(GATES.Z, [GATES.X], [1.0])
    U_goal = GATES.H
    ketdim = 2
    x_dim = 2 * ketdim * ketdim  # 8 (isomorphic unitary)
    u_dim = 1

    times = collect(range(0.0, T, N))
    pulse = CubicSplinePulse(fill(0.5, u_dim, N), fill(0.0, u_dim, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)
    integrator = SplineIntegrator(qtraj, N; spline_order = 3)

    traj = NamedTrajectory(qtraj, N)
    ∂F = get_jacobian_structure(integrator, traj)

    # Expected sparsity for UnitaryTrajectory with block-diagonal structure:
    # - Each of ketdim columns evolves independently → block_size = 2*ketdim
    # - ∂xₖ: ketdim blocks of (block_size × block_size) = 2 * (2ketdim)² per knot pair
    # - Parameters: x_dim * n_params where n_params = 4*u_dim + 2 (cubic)
    n_params_cubic = 4 * u_dim + 2  # 6
    state_nnz_per_k = 2 * (ketdim * (2 * ketdim)^2)  # 2 * 2 * 16 = 64
    param_nnz_per_k = x_dim * n_params_cubic  # 8 * 6 = 48
    expected_max_nnz = (N - 1) * (state_nnz_per_k + param_nnz_per_k)

    actual_nnz = nnz(∂F)
    println("  UnitaryTrajectory Jacobian nnz: $actual_nnz (expected ≤ $expected_max_nnz)")
    @test actual_nnz ≤ expected_max_nnz

    # Verify the sparsity is actually exploited (should be much less than dense)
    dense_nnz = (N - 1) * x_dim * (2 * x_dim + n_params_cubic)
    @test actual_nnz < dense_nnz
end

@testitem "E1: SplineIntegrator Jacobian sparsity for KetTrajectory" begin
    using DirectTrajOpt
    using Piccolo
    using NamedTrajectories
    using SparseArrays
    using LinearAlgebra

    T = 1.0
    N = 10
    sys = QuantumSystem(GATES.Z, [GATES.X], [1.0])
    ketdim = 2
    x_dim = 2 * ketdim  # 4 (isomorphic ket)
    u_dim = 1
    ψ_init = [1.0 + 0.0im, 0.0 + 0.0im]
    ψ_goal = [0.0 + 1.0im, 1.0 + 0.0im]

    times = collect(range(0.0, T, N))
    pulse = CubicSplinePulse(fill(0.5, u_dim, N), fill(0.0, u_dim, N), times)
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)
    integrator = SplineIntegrator(qtraj, N; spline_order = 3)

    traj = NamedTrajectory(qtraj, N)
    ∂F = get_jacobian_structure(integrator, traj)

    # Expected sparsity for KetTrajectory (dense state blocks):
    # - ∂xₖ and ∂xₖ₊₁: x_dim × x_dim each
    # - Parameters: x_dim * n_params where n_params = 4*u_dim + 2 = 6
    n_params_cubic = 4 * u_dim + 2
    state_nnz_per_k = 2 * x_dim * x_dim  # 32
    param_nnz_per_k = x_dim * n_params_cubic  # 24
    expected_max_nnz = (N - 1) * (state_nnz_per_k + param_nnz_per_k)

    actual_nnz = nnz(∂F)
    println("  KetTrajectory Jacobian nnz: $actual_nnz (expected ≤ $expected_max_nnz)")
    @test actual_nnz ≤ expected_max_nnz
end

@testitem "E1: SplineIntegrator MagnusGL4 basic test" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra

    T = 1.0
    N = 5
    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.5, 2, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, GATES.X)

    𝒮 = SplineIntegrator(qtraj, N; alg = MagnusGL4Alg(n_steps = 20))
    @test 𝒮 isa SplineIntegrator{UnitaryTrajectory,LinearSpline,ComplexF64,MagnusGL4Alg}
    @test 𝒮.ketdim == 2
    @test 𝒮.x_dim == 8  # 2*(2^2)

    traj = NamedTrajectory(qtraj, N)

    test_integrator(𝒮, traj; atol = 1e-4, gauss_newton = true, show_jacobian_diff = true)
end

@testitem "E1: SplineIntegrator MagnusGL4 cubic spline" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra

    T = 1.0
    N = 5
    sys = QuantumSystem(GATES.Z, [GATES.X], [1.0])

    times = collect(range(0.0, T, N))
    pulse = CubicSplinePulse(fill(0.5, 1, N), fill(0.0, 1, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, GATES.X)

    𝒮 = SplineIntegrator(qtraj, N; alg = MagnusGL4Alg(n_steps = 20))
    @test 𝒮 isa SplineIntegrator{UnitaryTrajectory,CubicSpline,ComplexF64,MagnusGL4Alg}

    traj = NamedTrajectory(qtraj, N)
    test_integrator(𝒮, traj; atol = 1e-4, gauss_newton = true, show_jacobian_diff = true)
end

@testitem "E1: SplineIntegrator MagnusGL4 with detuning" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra

    T = 1.0
    N = 5

    δ_init = 0.1
    sys = QuantumSystem(δ_init * GATES.Z, [GATES.X], [1.0])

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.5, 1, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, GATES.X)

    𝒮 = SplineIntegrator(qtraj, N; alg = MagnusGL4Alg(n_steps = 20))
    @test 𝒮 isa SplineIntegrator{UnitaryTrajectory,LinearSpline,ComplexF64,MagnusGL4Alg}

    traj = NamedTrajectory(qtraj, N)
    test_integrator(𝒮, traj; atol = 1e-4, gauss_newton = true)
end

@testitem "E1: Tsit5 vs MagnusGL4 residual agreement" begin
    using DirectTrajOpt
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra

    T = 1.0
    N = 5
    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.5, 2, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, GATES.X)

    𝒮_tsit5 = SplineIntegrator(qtraj, N; alg = Tsit5Alg())
    𝒮_magnus = SplineIntegrator(qtraj, N; alg = MagnusGL4Alg(n_steps = 20))

    traj = NamedTrajectory(qtraj, N)

    δ_tsit5 = zeros(𝒮_tsit5.dim)
    δ_magnus = zeros(𝒮_magnus.dim)

    evaluate!(δ_tsit5, 𝒮_tsit5, traj)
    evaluate!(δ_magnus, 𝒮_magnus, traj)

    println("  ‖δ_tsit5‖  = $(norm(δ_tsit5))")
    println("  ‖δ_magnus‖ = $(norm(δ_magnus))")
    println("  ‖δ_tsit5 - δ_magnus‖ = $(norm(δ_tsit5 - δ_magnus))")

    # Both methods should produce small residuals (trajectory is consistent with dynamics)
    @test norm(δ_tsit5) < 1e-6
    @test norm(δ_magnus) < 1e-6
    # And agree with each other to within solver tolerance
    @test norm(δ_tsit5 - δ_magnus) < 1e-6
end

@testitem "E1: SplineIntegrator with NonlinearDrive - Jacobian test" begin
    using DirectTrajOpt
    using Piccolo
    using SparseArrays
    using LinearAlgebra

    # 2-qubit system with 2 controls (u₁, u₂) and 3 drives:
    #   H = H_drift + u₁·σx⊗I + u₂·I⊗σx + (u₁·u₂)·σz⊗σz
    H_drift = kron(PAULIS.Z, PAULIS.I)
    drives = AbstractDrive[
        LinearDrive(sparse(ComplexF64.(kron(PAULIS.X, PAULIS.I))), 1),
        LinearDrive(sparse(ComplexF64.(kron(PAULIS.I, PAULIS.X))), 2),
        NonlinearDrive(kron(PAULIS.Z, PAULIS.Z), u -> u[1] * u[2]),
    ]
    drive_bounds = [1.0, 1.0]

    sys = QuantumSystem(H_drift, drives, drive_bounds)

    @test length(sys.H_drives) == 3
    @test has_nonlinear_drives(sys.H_drives)

    T = 1.0
    N = 5
    U_goal = kron(GATES.X, GATES.I)

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.3, 2, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)

    𝒮 = SplineIntegrator(qtraj, N)
    @test 𝒮 isa SplineIntegrator{UnitaryTrajectory}

    traj = NamedTrajectory(qtraj, N)
    test_integrator(𝒮, traj; atol = 1e-4, gauss_newton = true, show_jacobian_diff = true)
end

@testitem "E1: SplineIntegrator NonlinearDrive cubic splines" begin
    using DirectTrajOpt
    using Piccolo
    using SparseArrays
    using LinearAlgebra

    # 1-qubit system with a quadratic drive: H = u₁·σx + u₁²·σz
    drives = AbstractDrive[
        LinearDrive(sparse(ComplexF64.(PAULIS.X)), 1),
        NonlinearDrive(PAULIS.Z, u -> u[1]^2),
    ]
    drive_bounds = [1.0]

    sys = QuantumSystem(zeros(ComplexF64, 2, 2), drives, drive_bounds)

    T = 1.0
    N = 5
    U_goal = GATES.X

    times = collect(range(0.0, T, N))
    pulse = CubicSplinePulse(fill(0.3, 1, N), fill(0.0, 1, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)

    𝒮 = SplineIntegrator(qtraj, N)
    @test 𝒮 isa SplineIntegrator{UnitaryTrajectory,CubicSpline}

    traj = NamedTrajectory(qtraj, N)
    test_integrator(𝒮, traj; atol = 1e-4, gauss_newton = true, show_jacobian_diff = true)
end

@testitem "E1: SplineIntegrator NonlinearDrive MagnusGL4" begin
    using DirectTrajOpt
    using Piccolo
    using SparseArrays
    using LinearAlgebra

    # Same system as above: u₁·σx + u₁²·σz
    drives = AbstractDrive[
        LinearDrive(sparse(ComplexF64.(PAULIS.X)), 1),
        NonlinearDrive(PAULIS.Z, u -> u[1]^2),
    ]
    drive_bounds = [1.0]

    sys = QuantumSystem(zeros(ComplexF64, 2, 2), drives, drive_bounds)

    T = 1.0
    N = 5

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.3, 1, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, GATES.X)

    𝒮 = SplineIntegrator(qtraj, N; alg = MagnusGL4Alg(n_steps = 20))
    @test 𝒮 isa SplineIntegrator{UnitaryTrajectory,LinearSpline,ComplexF64,MagnusGL4Alg}

    traj = NamedTrajectory(qtraj, N)
    test_integrator(𝒮, traj; atol = 1e-4, gauss_newton = true, show_jacobian_diff = true)
end

@testitem "E1: SplineIntegrator UnitaryTrajectory with NonlinearDrive and global variables" begin
    using DirectTrajOpt
    using DirectTrajOpt: BoundsConstraint
    using Piccolo
    using NamedTrajectories
    using SparseArrays
    using LinearAlgebra
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators: has_global_dependence

    # System with nonlinear drive AND global parameter:
    #   H = δ·σz + u₁·σx + u₁²·σz
    δ_init = 0.01

    drives = AbstractDrive[
        LinearDrive(sparse(ComplexF64.(PAULIS.X)), 1),
        NonlinearDrive(PAULIS.Z, u -> u[1]^2),
    ]
    drive_bounds = [1.0]

    sys = QuantumSystem(PAULIS.Z, drives, drive_bounds; global_params = (δ = δ_init,))

    @test has_nonlinear_drives(sys.H_drives)

    T = 5.0
    N = 5
    U_goal = GATES.X

    times = collect(range(0.0, T, N))
    pulse = CubicSplinePulse(fill(0.5, 1, N), fill(0.0, 1, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)

    integrator = SplineIntegrator(qtraj, N; global_names = [:δ])

    @test has_global_dependence(integrator)
    @test integrator.global_dim == 1

    traj = NamedTrajectory(qtraj, N)
    @test haskey(traj.global_components, :δ)

    test_integrator(
        integrator,
        traj;
        atol = 1e-4,
        gauss_newton = true,
        show_jacobian_diff = true,
        show_hessian_diff = true,
    )

    # Full optimization with global bounds
    δ_bound = 0.5
    qcp = SplinePulseProblem(
        qtraj,
        N;
        Q = 100.0,
        R = 1e-2,
        integrator = integrator,
        global_bounds = Dict{Symbol,Union{Float64,Tuple{Float64,Float64}}}(:δ => δ_bound),
    )

    solve!(qcp; max_iter = 50, print_level = 0)

    result_traj = get_trajectory(qcp)
    δ_opt = result_traj.global_data[result_traj.global_components[:δ]][1]
    @test isfinite(δ_opt)
    @test δ_opt >= -δ_bound - 1e-5
    @test δ_opt <= δ_bound + 1e-5
    println("  Unitary NonlinearDrive + global: δ_init=$δ_init, δ_opt=$δ_opt")
end

@testitem "E1: SplineIntegrator{UnitaryTrajectory} fixed-step Tsit5 (adaptive=false)" begin
    using DirectTrajOpt
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra

    T = 1.0
    N = 5
    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    U_goal = GATES.X

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.5, 2, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)

    # Issue #180: adaptive=false crashed at construction (complex Φ_structure
    # vs the Float64-typed Tsit5Data field).
    𝒮 = SplineIntegrator(qtraj, N; alg = Tsit5Alg(adaptive = false, ode_h = 0.01))
    @test 𝒮 isa SplineIntegrator{UnitaryTrajectory,LinearSpline}

    traj = NamedTrajectory(qtraj, N)
    test_integrator(𝒮, traj; atol = 1e-4, gauss_newton = true)

    # ode_h is honored: a fine fixed step agrees with the adaptive solve,
    # while a single coarse step over the interval must differ.
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

@testitem "E1: ChebyshevAlg has no UnitaryTrajectory cell (ADR-0003 construction gate)" begin
    using Piccolo
    using LinearAlgebra

    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    U_goal = GATES.H
    N = 5
    times = collect(range(0.0, 1.0, N))
    pulse = LinearSplinePulse(fill(0.3, 2, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)

    # The matrix-free Chebyshev cell is Ket-only; the Unitary path errors at
    # construction instead of silently falling back to a Tsit5 forward.
    @test_throws ErrorException SplineIntegrator(qtraj, N; alg = ChebyshevAlg())
end

# ── E1 resume fill (#347): Magnus lanes (Adapt4 + Dual buffers), inner-ctor ── #
# ── fallback lanes, and the control-before-state Hessian symmetric branch. ─── #

@testitem "E1: Unitary MagnusAdapt4: construction + adaptive Magnus forward parity" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators: compute_ode_jacobian!

    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    U_init = Matrix{ComplexF64}(1.0I, 2, 2)
    N = 5
    times = collect(range(0.0, 1.0, N))
    pulse = LinearSplinePulse(fill(0.3, 2, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, U_init)
    traj = NamedTrajectory(qtraj, N)

    # MagnusAdapt4 construction builds the shared buffer alg-data
    𝒮_a = SplineIntegrator(qtraj, N; alg = MagnusAdapt4Alg(tol = 1e-10))
    @test 𝒮_a isa SplineIntegrator{UnitaryTrajectory,LinearSpline}
    @test 𝒮_a.alg isa MagnusAdapt4Alg

    # Adaptive-Magnus forward: Φ stays unitary and agrees with both the Tsit5
    # cell and the fixed-step GL4 cell.
    𝒮_t = SplineIntegrator(qtraj, N)
    𝒮_g = SplineIntegrator(qtraj, N; alg = MagnusGL4Alg(n_steps = 40))
    δ_a = zeros(𝒮_a.dim)
    δ_t = zeros(𝒮_t.dim)
    δ_g = zeros(𝒮_g.dim)
    evaluate!(δ_a, 𝒮_a, traj)
    evaluate!(δ_t, 𝒮_t, traj)
    evaluate!(δ_g, 𝒮_g, traj)
    @test norm(δ_a, Inf) < 1e-5
    @test norm(δ_a - δ_t, Inf) < 1e-5
    @test norm(δ_a - δ_g, Inf) < 1e-5

    # Lie-group structure: the per-knot propagator is unitary (read through the
    # Jacobian's cached complex Φ — solved by Tsit5 at the cell's tol, so the
    # deviation is solver-tolerance-bound, not machine precision)
    compute_ode_jacobian!(𝒮_a, traj[1], traj[2], 1, nothing)
    Φ = get_propagator(𝒮_a.prop_results[1], 2)
    @test norm(Φ'Φ - I, Inf) < 1e-8
end

@testitem "E1: Unitary MagnusGL4 forward propagates ForwardDiff Duals" begin
    using DirectTrajOpt
    using Piccolo
    using ForwardDiff
    using LinearAlgebra
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators: _forward_propagate

    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    U_init = Matrix{ComplexF64}(1.0I, 2, 2)
    N = 5
    times = collect(range(0.0, 1.0, N))
    pulse = LinearSplinePulse(fill(0.3, 2, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, U_init)
    traj = NamedTrajectory(qtraj, N)

    𝒮 = SplineIntegrator(qtraj, N; alg = MagnusGL4Alg(n_steps = 10))
    pₖ = [fill(0.3, 4); 1.0 / (N - 1); 0.0]  # [uₖ; uₖ₊₁; Δt; t]

    # Float64 lane: pre-allocated per-knot buffers
    Φ_real = _forward_propagate(𝒮, 1, pₖ)
    @test size(Φ_real) == (4, 4)

    # Dual lane: fresh type-matched buffers, values identical to the Float64
    # path (the tangent machinery rides the same fixed-step Magnus product)
    p_dual = ForwardDiff.Dual.(pₖ, one.(pₖ))
    Φ_dual = _forward_propagate(𝒮, 1, p_dual)
    @test eltype(Φ_dual) <: ForwardDiff.Dual
    @test isapprox(ForwardDiff.value.(Φ_dual), Φ_real; atol = 1e-12)
end

@testitem "E1: _spline_unitary inner lanes: pulse inference, globals, drive gates" begin
    using DirectTrajOpt
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators: _spline_unitary, spline_order

    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    U_init = Matrix{ComplexF64}(1.0I, 2, 2)
    N = 5
    times = collect(range(0.0, 1.0, N))

    # Non-spline pulse through the inner constructor: defaults to linear
    zo = ZeroOrderPulse(0.3 * fill(1.0, 2, N), times)
    qtraj_zo = UnitaryTrajectory(sys, zo, U_init)
    traj_zo = NamedTrajectory(qtraj_zo, N)
    𝒮 = _spline_unitary(sys, zo, :Ũ⃗, :u, traj_zo)
    @test spline_order(𝒮) == 1
    δ = zeros(𝒮.dim)
    evaluate!(δ, 𝒮, traj_zo)
    @test norm(δ, Inf) < 1e-5

    # Explicit global_names on a traj WITHOUT globals: the fallback dims the
    # u-vector for globals, then the global seed read fails loudly (a traj must
    # carry the global data it declares — the empty global-component NamedTuple
    # throws FieldError on :δ access).
    @test_throws FieldError _spline_unitary(
        sys,
        zo,
        :Ũ⃗,
        :u,
        traj_zo;
        global_names = [:δ],
    )

    # Drive-free system: the explicit-drives ArgumentError
    sys_free = QuantumSystem(GATES.Z)
    @test_throws ArgumentError _spline_unitary(sys_free, zo, :Ũ⃗, :u, traj_zo)

    # ChebyshevAlg refusal: unitary states buy nothing from matrix-free
    @test_throws ErrorException _spline_unitary(
        sys,
        zo,
        :Ũ⃗,
        :u,
        traj_zo;
        alg = ChebyshevAlg(bracket = (-4.0, 4.0)),
    )
end

@testitem "E1: Unitary cubic without du bounds zero-fills the derivative seed" begin
    using DirectTrajOpt
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators: _spline_unitary, spline_order

    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    U_init = Matrix{ComplexF64}(1.0I, 2, 2)
    N = 5
    times = collect(range(0.0, 1.0, N))
    cpulse = CubicSplinePulse(fill(0.3, 2, N), fill(0.0, 2, N), times)
    qtraj = UnitaryTrajectory(sys, cpulse, U_init)
    traj_full = NamedTrajectory(qtraj, N)
    @test haskey(traj_full.bounds, :du)

    # Strip the du bounds (keep the du component): the ctor seeds du from zeros
    boundfree = NamedTrajectory(traj_full; bounds = (u = 1.0,))
    @test !haskey(boundfree.bounds, :du)

    𝒮 = _spline_unitary(sys, cpulse, :Ũ⃗, :u, boundfree)
    @test spline_order(𝒮) == 3
    δ = zeros(𝒮.dim)
    evaluate!(δ, 𝒮, boundfree)
    @test norm(δ, Inf) < 1e-5
end

@testitem "E1: control-before-state layout: unitary Hessian fills both triangles" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra
    using SparseArrays
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators: _spline_unitary, get_pulse

    # A trajectory whose control component PRECEDES the state: parameter
    # columns land BEFORE x columns, so the (x,p) Hessian fill takes the
    # mirrored triangle branch for every parameter block. The state columns
    # are a genuine constant-control rollout (exact matrix exponentials).
    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    N = 5
    T = 1.0
    Δt = T / (N - 1)
    H = GATES.Z + 0.3 * (GATES.X + GATES.Y)
    Us = [exp(-im * H * (k * Δt)) for k = 0:(N-1)]
    traj = NamedTrajectory(
        (
            u = 0.3 * fill(1.0, 2, N),
            Ũ⃗ = hcat(operator_to_iso_vec.(Us)...),
            Δt = fill(Δt, N),
            t = collect(range(0.0, T, N)),
        );
        controls = :u,
        timestep = :Δt,
        bounds = (u = 1.0,),
    )

    qtraj = UnitaryTrajectory(
        sys,
        LinearSplinePulse(0.3 * fill(1.0, 2, N), collect(range(0.0, T, N))),
        Matrix{ComplexF64}(1.0I, 2, 2),
    )
    𝒮 = _spline_unitary(sys, get_pulse(qtraj), :Ũ⃗, :u, traj)

    # Forward is dynamically consistent on the mirrored layout
    δ = zeros(𝒮.dim)
    evaluate!(δ, 𝒮, traj)
    @test norm(δ, Inf) < 1e-5

    # Jacobian conformance on the mirrored layout
    test_integrator(𝒮, traj; atol = 1e-4, gauss_newton = true)

    # The Hessian fill populates BOTH the (x<p) and mirrored (x>p) triangles
    μ = ones(𝒮.dim)
    H = eval_hessian_of_lagrangian(𝒮, traj, μ)
    @test nnz(H) > 0
end
