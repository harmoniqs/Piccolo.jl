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
