# ============================================================================ #
# Cluster E1 (#347) — SplineIntegrator KetTrajectory cell behavior tests.
#
# Ported from the pre-split embedded suite (Piccolissimo
# src/integrators/spline/spline_integrator_ket.jl). The Ket DENSE cell's
# Jacobian conformance, exact-Hessian lanes, ket-level-sensitivity lane, and
# constructor error lanes sat at 40.6% because only the forward path was
# exercised (dense parity + templates).
#
# NOT ported (proprietary surface, exercised by Piccolissimo's own suite):
# ChebyshevAlg/Magnus matrix-free VVP/JVP/HVP items, step_jvp!/step_vjp!/
# step_hvp!, the Rodas5PAlg stiff forward (the `_stiff_rodas5p_solve` hook has
# no Piccolo method), and the ChebyshevAlg forward (`_chebyshev_forward!` hook).
# ============================================================================ #

@testitem "E1: testing SplineIntegrator with KetTrajectory" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra

    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])

    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]

    T = 1.0
    N = 5

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.5, 2, N), times)
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)
    traj = NamedTrajectory(qtraj, N)

    𝒮 = SplineIntegrator(qtraj, N)
    @test 𝒮 isa SplineIntegrator{KetTrajectory,LinearSpline}

    test_integrator(𝒮, traj; atol = 1e-4, gauss_newton = true, show_hessian_diff = true)
end

@testitem "E1: testing SplineIntegrator cubic spline with KetTrajectory" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra

    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])

    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]

    T = 1.0
    N = 5

    times = collect(range(0.0, T, N))
    pulse = CubicSplinePulse(fill(0.5, 2, N), fill(0.0, 2, N), times)
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)
    traj = NamedTrajectory(qtraj, N)

    𝒮 = SplineIntegrator(qtraj, N)
    @test 𝒮 isa SplineIntegrator{KetTrajectory,CubicSpline}

    test_integrator(
        𝒮,
        traj;
        atol = 1e-4,
        gauss_newton = true,
        show_jacobian_diff = true,
        show_hessian_diff = true,
    )
end

@testitem "E1: SplineIntegrator with KetTrajectory and SplinePulseProblem" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra

    T = 1.0
    u_bounds = [1.0, 1.0]

    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], u_bounds)

    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]

    N = 10

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.5, 2, N), times)
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

    integrator = SplineIntegrator(qtraj, N; spline_order = 1)

    @test integrator isa SplineIntegrator{KetTrajectory}

    qcp = SplinePulseProblem(qtraj, N; Q = 50.0, R = 1e-3, integrator = integrator)

    @test qcp isa QuantumControlProblem
    @test length(qcp.prob.integrators) == 2  # dynamics + derivative integrator on du
    @test qcp.prob.integrators[1] isa SplineIntegrator{KetTrajectory}

    solve!(qcp; max_iter = 10, print_level = 0)
end

@testitem "E1: SplineIntegrator cubic splines KetTrajectory with SplinePulseProblem" begin
    using DirectTrajOpt
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra

    T = 1.0
    u_bounds = [1.0, 1.0]

    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], u_bounds)

    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]

    N = 10

    times = collect(range(0.0, T, N))
    pulse = CubicSplinePulse(fill(0.5, 2, N), fill(0.0, 2, N), times)
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

    integrator = SplineIntegrator(qtraj, N; spline_order = 3)

    @test integrator isa SplineIntegrator{KetTrajectory}
    @test integrator isa SplineIntegrator{KetTrajectory,CubicSpline}

    qcp = SplinePulseProblem(qtraj, N; Q = 50.0, R = 1e-3, integrator = integrator)

    @test qcp isa QuantumControlProblem
    @test length(qcp.prob.integrators) == 1  # dynamics
    @test qcp.prob.integrators[1] isa SplineIntegrator{KetTrajectory,CubicSpline}

    solve!(qcp; max_iter = 10, print_level = 0)

    traj = get_trajectory(qcp)
    @test traj isa NamedTrajectory
end

@testitem "E1: SplineIntegrator KetTrajectory with global parameter optimization" begin
    using DirectTrajOpt
    using DirectTrajOpt: BoundsConstraint
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators: has_global_dependence

    T = 5.0
    N = 5

    δ_init = 0.01

    sys = QuantumSystem(GATES.Z, [GATES.X], [1.0]; global_params = (δ = δ_init,))

    ψ_init = ComplexF64[1.0, 0.0]  # |0⟩
    ψ_goal = ComplexF64[0.0, 1.0]  # |1⟩

    times = collect(range(0, T, N))
    pulse = CubicSplinePulse(fill(0.5, 1, N), fill(0.0, 1, N), times)
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

    integrator = SplineIntegrator(qtraj, N; spline_order = 3, global_names = [:δ])

    @test has_global_dependence(integrator)
    @test integrator.global_dim == 1

    traj = NamedTrajectory(qtraj, N)
    @test haskey(traj.global_components, :δ)
    @test traj.global_dim == 1

    test_integrator(
        integrator,
        traj;
        atol = 1e-4,
        gauss_newton = true,
        show_jacobian_diff = true,
        show_hessian_diff = true,
    )

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

    result_traj = get_trajectory(qcp)
    @test haskey(result_traj.global_components, :δ)
    @test result_traj.global_dim == 1
end

@testitem "E1: SplineIntegrator rejects function-based QuantumSystem (KetTrajectory)" begin
    using Piccolo
    using LinearAlgebra

    # Function-based system: H_drives will be empty
    H = (u, t) -> GATES.Z + u[1] * GATES.X + u[2] * GATES.Y
    sys = QuantumSystem(H, [1.0, 1.0])
    @test isempty(sys.H_drives)

    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]

    T = 1.0
    N = 5
    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.5, 2, N), times)
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

    # Should throw ArgumentError with helpful message
    err = try
        SplineIntegrator(qtraj, N)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("explicit drive operators", err.msg)
    @test occursin("Function-based systems", err.msg)
end

@testitem "E1: SplineIntegrator KetTrajectory with NonlinearDrive" begin
    using DirectTrajOpt
    using Piccolo
    using SparseArrays
    using LinearAlgebra

    # 1-qubit system with 1 control and 2 drives:
    #   H = σz + u₁·σx + u₁²·σz
    drives = AbstractDrive[
        LinearDrive(sparse(ComplexF64.(PAULIS.X)), 1),
        NonlinearDrive(PAULIS.Z, u -> u[1]^2),
    ]
    drive_bounds = [1.0]

    sys = QuantumSystem(PAULIS.Z, drives, drive_bounds)

    @test length(sys.H_drives) == 2
    @test has_nonlinear_drives(sys.H_drives)

    T = 1.0
    N = 5

    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.3, 1, N), times)
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

    𝒮 = SplineIntegrator(qtraj, N)
    @test 𝒮 isa SplineIntegrator{KetTrajectory}

    traj = NamedTrajectory(qtraj, N)
    test_integrator(𝒮, traj; atol = 1e-4, gauss_newton = true, show_jacobian_diff = true)
end

@testitem "E1: SplineIntegrator KetTrajectory NonlinearDrive cubic splines" begin
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

    ψ_init = ComplexF64[1.0, 0.0, 0.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0, 0.0, 0.0]

    times = collect(range(0.0, T, N))
    pulse = CubicSplinePulse(fill(0.3, 2, N), fill(0.0, 2, N), times)
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

    𝒮 = SplineIntegrator(qtraj, N)
    @test 𝒮 isa SplineIntegrator{KetTrajectory,CubicSpline}

    traj = NamedTrajectory(qtraj, N)
    test_integrator(𝒮, traj; atol = 1e-4, gauss_newton = true, show_jacobian_diff = true)
end

@testitem "E1: SplineIntegrator KetTrajectory with NonlinearDrive and global variables" begin
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

    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]

    times = collect(range(0.0, T, N))
    pulse = CubicSplinePulse(fill(0.5, 1, N), fill(0.0, 1, N), times)
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

    integrator = SplineIntegrator(qtraj, N; spline_order = 3, global_names = [:δ])

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
    println("  NonlinearDrive + global: δ_init=$δ_init, δ_opt=$δ_opt")
end

@testitem "E1: SplineIntegrator exact Hessian linear spline NonlinearDrive" begin
    using DirectTrajOpt
    using Piccolo
    using SparseArrays
    using LinearAlgebra

    # 1-qubit system with quadratic nonlinear drive: H = σz + u₁·σx + u₁²·σz
    drives = AbstractDrive[
        LinearDrive(sparse(ComplexF64.(PAULIS.X)), 1),
        NonlinearDrive(PAULIS.Z, u -> u[1]^2),
    ]
    drive_bounds = [1.0]

    sys = QuantumSystem(PAULIS.Z, drives, drive_bounds)

    T = 1.0
    N = 5

    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.3, 1, N), times)
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

    𝒮 = SplineIntegrator(qtraj, N; exact_hessian = true)
    @test 𝒮 isa SplineIntegrator{KetTrajectory}
    @test 𝒮.exact_hessian
    @test !isnothing(𝒮.hess_probs)
    @test !isnothing(𝒮.hess_active_pairs)
    @test length(𝒮.hess_active_pairs) > 0

    traj = NamedTrajectory(qtraj, N)

    # gauss_newton=false → compare FULL Hessian against finite differences
    test_integrator(𝒮, traj; atol = 1e-4, gauss_newton = false, show_hessian_diff = true)
end

@testitem "E1: SplineIntegrator exact Hessian cubic spline NonlinearDrive" begin
    using DirectTrajOpt
    using Piccolo
    using SparseArrays
    using LinearAlgebra

    # 2-control system with product nonlinearity:
    #   H = σz⊗I + u₁·σx⊗I + u₂·I⊗σx + (u₁·u₂)·σz⊗σz
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

    ψ_init = ComplexF64[1.0, 0.0, 0.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0, 0.0, 0.0]

    times = collect(range(0.0, T, N))
    pulse = CubicSplinePulse(fill(0.3, 2, N), fill(0.0, 2, N), times)
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

    𝒮 = SplineIntegrator(qtraj, N; exact_hessian = true)
    @test 𝒮 isa SplineIntegrator{KetTrajectory,CubicSpline}
    @test 𝒮.exact_hessian

    traj = NamedTrajectory(qtraj, N)

    # Full Hessian comparison
    test_integrator(𝒮, traj; atol = 1e-4, gauss_newton = false, show_hessian_diff = true)
end

@testitem "E1: SplineIntegrator exact Hessian with LinearDrive-only matches FiniteDiff" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra

    # All-linear system: exact Hessian has (p,p) blocks from ODE nonlinearity
    # (even though drive coefficients are linear, the exponential is not)
    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])

    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]

    T = 1.0
    N = 5

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.5, 2, N), times)
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)
    traj = NamedTrajectory(qtraj, N)

    𝒮 = SplineIntegrator(qtraj, N; exact_hessian = true)
    @test 𝒮.exact_hessian

    # Full Hessian should match finite differences
    test_integrator(𝒮, traj; atol = 1e-4, gauss_newton = false)
end

@testitem "E1: use_ket_sensitivity=true on KetTrajectory: Jacobian agrees with default" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra
    using Random
    Random.seed!(123)

    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]

    T = 1.0
    N = 5
    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.5, 2, N), times)
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)
    traj = NamedTrajectory(qtraj, N)

    # Default: use_ket_sensitivity = false (propagator-level path)
    𝒮_default = SplineIntegrator(qtraj, N; use_ket_sensitivity = false)

    # Ket-level sensitivity path
    𝒮_ket = SplineIntegrator(qtraj, N; use_ket_sensitivity = true)
    @test 𝒮_ket.use_ket_sensitivity
    @test !isnothing(𝒮_ket.ket_sens_results)

    # Compute Jacobian via both paths; should agree to ~1e-7
    J_default = eval_jacobian(𝒮_default, traj)
    J_ket = eval_jacobian(𝒮_ket, traj)

    @test size(J_default) == size(J_ket)
    @test norm(J_default - J_ket) / max(norm(J_default), 1e-12) < 1e-7
end

@testitem "E1: use_ket_sensitivity incompatibility gates on KetTrajectory" begin
    using Piccolo
    using LinearAlgebra

    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]

    N = 5
    times = collect(range(0.0, 1.0, N))
    pulse = LinearSplinePulse(fill(0.3, 2, N), times)
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

    # Ket-level sensitivities cannot express the second-order T_{ij} matrices
    @test_throws ErrorException SplineIntegrator(
        qtraj,
        N;
        use_ket_sensitivity = true,
        exact_hessian = true,
    )

    # The matrix-state forward ODE of the ket-level path has no closed-form
    # Jacobian shape for Rodas5P yet
    @test_throws ErrorException SplineIntegrator(
        qtraj,
        N;
        use_ket_sensitivity = true,
        alg = Rodas5PAlg(),
    )
end

@testitem "E1: ChebyshevAlg sensitivity requests error on KetTrajectory (ADR-0003 Decision 2 gate)" begin
    using Piccolo
    using LinearAlgebra

    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]

    N = 5
    times = collect(range(0.0, 1.0, N))
    pulse = LinearSplinePulse(fill(0.3, 2, N), times)
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

    # Construction-time gates: explicit sensitivity requests error before any
    # matrix-free machinery is reached (the concrete ChebyshevAlg cell itself
    # is proprietary; the GATES are the open-core contract).
    @test_throws ErrorException SplineIntegrator(
        qtraj,
        N;
        alg = ChebyshevAlg(),
        exact_hessian = true,
    )
    @test_throws ErrorException SplineIntegrator(
        qtraj,
        N;
        alg = ChebyshevAlg(),
        use_ket_sensitivity = true,
    )

    # Trajectory types without the Chebyshev cell error at construction
    U_goal = GATES.H
    uqtraj = UnitaryTrajectory(sys, pulse, U_goal)
    @test_throws ErrorException SplineIntegrator(uqtraj, N; alg = ChebyshevAlg())

    initials = [ComplexF64[1.0, 0.0], ComplexF64[0.0, 1.0]]
    goals = [ComplexF64[0.0, 1.0], ComplexF64[1.0, 0.0]]
    mqtraj = MultiKetTrajectory(sys, pulse, initials, goals)
    # MultiKet ChebyshevAlg is cubic-only: the linear pulse errors at
    # construction rather than deep in the first VJP/HVP sweep.
    @test_throws ErrorException SplineIntegrator(mqtraj, N; alg = ChebyshevAlg())
end

@testitem "E1: SplineIntegrator{KetTrajectory} fixed-step Tsit5 (adaptive=false)" begin
    using DirectTrajOpt
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra

    T = 1.0
    N = 5
    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.5, 2, N), times)
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

    # Issue #180: adaptive=false crashed at construction (complex Φ_structure
    # vs the Float64-typed Tsit5Data field).
    𝒮 = SplineIntegrator(qtraj, N; alg = Tsit5Alg(adaptive = false, ode_h = 0.01))
    @test 𝒮 isa SplineIntegrator{KetTrajectory,LinearSpline}

    traj = NamedTrajectory(qtraj, N)
    test_integrator(𝒮, traj; atol = 1e-4, gauss_newton = true)

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

# ── E1 resume fill (#347): inner-constructor lanes, Rodas5P/Magnus/Chebyshev ── #
# ── gates, the ket-sensitivity forward, and the Hessian solve paths. ──────── #

@testitem "E1: _spline_ket inner lanes: pulse inference + globals auto-detect" begin
    using DirectTrajOpt
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators: _spline_ket, spline_order

    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    N = 5
    times = collect(range(0.0, 1.0, N))

    # Non-spline pulse through the inner constructor: defaults to order 1
    zo = ZeroOrderPulse(0.3 * fill(1.0, 2, N), times)
    qtraj_zo = KetTrajectory(sys, zo, ψ_init, ψ_goal)
    traj_zo = NamedTrajectory(qtraj_zo, N)
    𝒮 = _spline_ket(sys, zo, state_name(qtraj_zo), drive_name(qtraj_zo), traj_zo)
    @test spline_order(𝒮) == 1
    δ = zeros(𝒮.dim)
    evaluate!(δ, 𝒮, traj_zo)
    @test norm(δ, Inf) < 1e-5

    # Globals auto-detect via the inner constructor: a system WITH global_params
    # yields the collected names without being asked. The pulse's drive count
    # must match the system (one drive here).
    gsys = QuantumSystem(GATES.Z, [GATES.X], [1.0]; global_params = (δ = 0.01,))
    zo_g = ZeroOrderPulse(0.3 * fill(1.0, 1, N), times)
    gq = KetTrajectory(gsys, zo_g, ψ_init, ψ_goal)
    gtraj = NamedTrajectory(gq, N)
    𝒮_g = _spline_ket(gsys, zo_g, state_name(gq), drive_name(gq), gtraj)
    @test 𝒮_g.global_names == [:δ]
    @test 𝒮_g.global_dim == 1

    # …and a drive-global-free system resolves to an empty name list
    𝒮_0 = _spline_ket(sys, zo, state_name(qtraj_zo), drive_name(qtraj_zo), traj_zo)
    @test isempty(𝒮_0.global_names)
    @test 𝒮_0.global_dim == 0
end

@testitem "E1: Rodas5P ket cell: closed-form Jacobians wired, solves gated" begin
    using DirectTrajOpt
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators: compute_ode_jacobian!

    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    N = 5
    times = collect(range(0.0, 1.0, N))
    pulse = LinearSplinePulse(fill(0.3, 2, N), times)
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)
    traj = NamedTrajectory(qtraj, N)

    # Construction succeeds in the open core: the forward ODE problems carry the
    # closed-form complex Jacobian (built by build_forward_ode_jacobian), and
    # the sensitivity problems carry the block-diagonal one.
    𝒮 = SplineIntegrator(qtraj, N; alg = Rodas5PAlg())
    @test 𝒮 isa SplineIntegrator{KetTrajectory,LinearSpline}
    @test 𝒮.alg isa Rodas5PAlg
    @test !isnothing(𝒮.sens_probs)
    # The sensitivity problems were rebuilt with the closed-form jac attached
    @test 𝒮.sens_probs[1].f.jac !== nothing

    # Forward solve: the concrete stiff solver is proprietary (slice 3b
    # de-scope) — _stiff_rodas5p_solve has no open-core method, so the forward
    # errors loudly instead of silently falling back. The per-knot cells run
    # under Threads.@threads, so the nested MethodError surfaces wrapped in a
    # CompositeException (one TaskFailedException per failing knot task).
    δ = zeros(𝒮.dim)
    @test_throws CompositeException evaluate!(δ, 𝒮, traj)

    # The analytic-Jacobian machinery is still exercised at construction (the
    # Hamiltonian jacobian + ODE jacobian builders); the sensitivity solve hits
    # the same proprietary hook — a DIRECT call here is unwrapped, so the bare
    # MethodError is observable.
    @test_throws MethodError compute_ode_jacobian!(𝒮, traj[1], traj[2], 1, nothing)
end

@testitem "E1: ket algorithm gates: Magnus and Chebyshev construction is proprietary" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra

    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    N = 5
    times = collect(range(0.0, 1.0, N))
    pulse = LinearSplinePulse(fill(0.3, 2, N), times)
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

    # The Magnus ket cell's ChebyshevData (matrix-free midpoint-Magnus cores)
    # is built by the proprietary _build_magnus_ket_data hook
    @test_throws MethodError SplineIntegrator(qtraj, N; alg = MagnusGL4Alg())
    @test_throws MethodError SplineIntegrator(qtraj, N; alg = MagnusAdapt4Alg())

    # ChebyshevAlg alg-data construction routes through the shared (undefined
    # in open core) _build_alg_data dispatch
    @test_throws MethodError SplineIntegrator(
        qtraj,
        N;
        alg = ChebyshevAlg(bracket = (-4.0, 4.0)),
    )
end

@testitem "E1: use_ket_sensitivity forward materializes and caches the propagator" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra

    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    N = 5
    times = collect(range(0.0, 1.0, N))
    pulse = LinearSplinePulse(fill(0.3, 2, N), times)
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)
    traj = NamedTrajectory(qtraj, N)

    𝒮_k = SplineIntegrator(qtraj, N; use_ket_sensitivity = true)
    𝒮_d = SplineIntegrator(qtraj, N)

    # The ket-sensitivity forward solves the identity-initialized matrix ODE
    # (returning Φₖ itself), caches it, and applies it to ψₖ — agreeing with
    # the default vector-state forward.
    δ_k = zeros(𝒮_k.dim)
    δ_d = zeros(𝒮_d.dim)
    evaluate!(δ_k, 𝒮_k, traj)
    evaluate!(δ_d, 𝒮_d, traj)
    @test norm(δ_k - δ_d, Inf) < 1e-10
    @test norm(δ_k, Inf) < 1e-5
    for k = 1:(N-1)
        @test !iszero(𝒮_k.prop_results[k].Φ_vec)
    end
end

@testitem "E1: ket compute_ode_hessian! lanes: exact Hessian solve + no-hessian gate" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators: compute_ode_hessian!

    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    N = 5
    times = collect(range(0.0, 1.0, N))
    pulse = LinearSplinePulse(fill(0.3, 2, N), times)
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)
    traj = NamedTrajectory(qtraj, N)

    # A non-exact cell has no second-order problems: compute_ode_hessian!
    # errors loudly instead of silently degrading to Gauss-Newton.
    𝒮_gn = SplineIntegrator(qtraj, N)
    @test_throws ErrorException compute_ode_hessian!(𝒮_gn, traj[1], traj[2], 1, nothing)

    # The exact cell solves the second-order sensitivity ODE per knot and the
    # result agrees with the Jacobian's sensitivities at the same knot.
    𝒮_ex = SplineIntegrator(qtraj, N; exact_hessian = true)
    compute_ode_hessian!(𝒮_ex, traj[1], traj[2], 1, nothing)
    @test !isnothing(𝒮_ex.hess_probs)

    # Jacobian/Hessian conformance for the SAME cell (dense exact Hessian)
    test_integrator(𝒮_ex, traj; atol = 1e-4, gauss_newton = false)
end

@testitem "E1: exact-Hessian cubic spline matches FiniteDiff (SOSE du/Δt chains)" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra

    # Cubic + exact Hessian drives the second-order sensitivity ODE's full
    # order-3 closure: the Hermite basis du blocks, the Δt↔du chain rules and
    # the (du, Δt) mixed pairs that the linear-spline exact tests never touch.
    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    N = 5
    T = 1.0
    times = collect(range(0.0, T, N))
    pulse = CubicSplinePulse(fill(0.4, 2, N), 0.2 * fill(1.0, 2, N), times)
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)
    traj = NamedTrajectory(qtraj, N)

    𝒮 = SplineIntegrator(qtraj, N; exact_hessian = true)
    @test 𝒮.exact_hessian
    @test !isnothing(𝒮.hess_probs)
    test_integrator(𝒮, traj; atol = 1e-4, gauss_newton = false)
end
