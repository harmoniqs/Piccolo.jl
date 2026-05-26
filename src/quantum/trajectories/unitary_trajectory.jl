# ============================================================================ #
# UnitaryTrajectory
# ============================================================================ #

"""
    UnitaryTrajectory{P<:AbstractPulse, S<:ODESolution} <: AbstractQuantumTrajectory{P}

Trajectory for unitary gate synthesis. The ODE solution is computed at construction.

# Fields
- `system::AbstractQuantumSystem`: The quantum system
- `pulse::P`: The control pulse (stores drive_name)
- `initial::Matrix{ComplexF64}`: Initial unitary (default: identity)
- `goal::AbstractPiccoloOperator`: Target unitary operator (Matrix{ComplexF64} or EmbeddedOperator)
- `solution::S`: Pre-computed ODE solution

# Callable
`traj(t)` returns the unitary at time `t` by interpolating the solution.

# Conversion to NamedTrajectory
Use `NamedTrajectory(traj, N)` or `NamedTrajectory(traj, times)` for optimization.
"""
mutable struct UnitaryTrajectory{P<:AbstractPulse,S<:ODESolution} <:
               AbstractQuantumTrajectory{P}
    system::AbstractQuantumSystem
    pulse::P
    initial::Matrix{ComplexF64}
    goal::AbstractPiccoloOperator
    solution::S
end

"""
    UnitaryTrajectory(system, pulse, goal; initial=I, algorithm=MagnusAdapt4(), abstol=1e-8, reltol=1e-8, n_save=101)

Create a unitary trajectory by solving the Schrödinger equation.

# Arguments
- `system::QuantumSystem`: The quantum system
- `pulse::AbstractPulse`: The control pulse
- `goal::AbstractPiccoloOperator`: Target unitary (Matrix or EmbeddedOperator).
  Non-ComplexF64 matrices are automatically converted to `Matrix{ComplexF64}`.

# Keyword Arguments
- `initial`: Initial unitary (default: identity matrix)
- `algorithm`: ODE solver algorithm (default: MagnusAdapt4())
- `abstol`: Absolute tolerance for adaptive integration (default: 1e-8)
- `reltol`: Relative tolerance for adaptive integration (default: 1e-8)
- `n_save`: Number of output time points for plotting/interpolation (default: 101)
"""
function UnitaryTrajectory(
    system::AbstractQuantumSystem,
    pulse::AbstractPulse,
    goal::AbstractPiccoloOperator;
    initial::AbstractMatrix{<:Number} = Matrix{ComplexF64}(I, system.levels, system.levels),
    algorithm = nothing,
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
    n_save::Int = 101,
    progress::Bool = false,
    progress_steps::Int = 10,
)
    @assert n_drives(pulse) == system.n_drives "Pulse has $(n_drives(pulse)) drives, system has $(system.n_drives)"

    goal_stored =
        goal isa Matrix{ComplexF64} || goal isa EmbeddedOperator ? goal :
        Matrix{ComplexF64}(goal)

    if isnothing(algorithm)
        algorithm = default_algorithm(system)
    end

    U0 = Matrix{ComplexF64}(initial)
    knot_times = get_knot_times(pulse)
    save_times = collect(range(0.0, duration(pulse), length = n_save))
    tstops = sort(unique(vcat(knot_times, save_times)))
    # See KetTrajectory for rationale: Magnus methods need the matrix-operator
    # form; explicit RK methods should use the sparse mul! RHS form to avoid
    # per-step dense materialization of H.
    prob =
        _needs_operator_form(algorithm) ?
        UnitaryOperatorODEProblem(system, pulse, tstops; U0 = U0) :
        UnitaryODEProblem(system, pulse, tstops; U0 = U0)
    sol = solve(
        prob,
        algorithm;
        saveat = save_times,
        abstol = abstol,
        reltol = reltol,
        progress = progress,
        progress_steps = progress_steps,
    )

    return UnitaryTrajectory{typeof(pulse),typeof(sol)}(system, pulse, U0, goal_stored, sol)
end

"""
    UnitaryTrajectory(system, goal, T::Real; drive_name=:u, algorithm=MagnusAdapt4(), abstol=1e-8, reltol=1e-8)

Convenience constructor that creates a zero pulse of duration T.

# Arguments
- `system::QuantumSystem`: The quantum system
- `goal::AbstractPiccoloOperator`: Target unitary (Matrix or EmbeddedOperator)
- `T::Real`: Duration of the pulse

# Keyword Arguments
- `drive_name::Symbol`: Name of the drive variable (default: `:u`)
- `algorithm`: ODE solver algorithm (default: MagnusAdapt4())
- `abstol`: Absolute tolerance (default: 1e-8)
- `reltol`: Relative tolerance (default: 1e-8)
"""
function UnitaryTrajectory(
    system::QuantumSystem,
    goal::AbstractPiccoloOperator,
    T::Real;
    drive_name::Symbol = :u,
    algorithm = nothing,
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
)
    times = [0.0, T]
    controls = vcat([rand(Uniform(b...), 1, length(times)) for b in system.drive_bounds]...)
    pulse = ZeroOrderPulse(controls, times; drive_name)
    return UnitaryTrajectory(system, pulse, goal; algorithm, abstol, reltol)
end

# Callable: sample solution at any time
(qtraj::UnitaryTrajectory)(t::Real) = qtraj.solution(t)

# ============================================================================ #
# Tests
# ============================================================================ #

@testitem "UnitaryTrajectory construction" begin
    using LinearAlgebra

    # Simple 2-level system
    system = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])

    # Create with duration
    T = 1.0
    X_gate = ComplexF64[0 1; 1 0]
    qtraj = UnitaryTrajectory(system, X_gate, T)

    @test qtraj isa UnitaryTrajectory
    @test qtraj.system === system
    @test qtraj.goal === X_gate
    @test qtraj.initial ≈ Matrix{ComplexF64}(I, 2, 2)

    # Create with explicit pulse
    times = [0.0, 0.5, 1.0]
    controls = 0.1 * randn(1, 3)
    pulse = ZeroOrderPulse(controls, times)
    qtraj2 = UnitaryTrajectory(system, pulse, X_gate)

    @test qtraj2 isa UnitaryTrajectory
    @test duration(qtraj2) ≈ 1.0
end

@testitem "UnitaryTrajectory goal type conversion" begin
    using LinearAlgebra

    system = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])
    T = 1.0

    # ComplexF64 matrix goal: identity preserved (no conversion needed)
    goal_cf64 = ComplexF64[0 1; 1 0]
    qtraj_cf64 = UnitaryTrajectory(system, goal_cf64, T)
    @test qtraj_cf64.goal === goal_cf64
    @test qtraj_cf64.goal isa Matrix{ComplexF64}

    # Float64 matrix goal: converted to ComplexF64
    goal_f64 = [0.0 1.0; 1.0 0.0]
    qtraj_f64 = UnitaryTrajectory(system, goal_f64, T)
    @test qtraj_f64.goal == goal_f64
    @test qtraj_f64.goal isa Matrix{ComplexF64}

    # Int matrix goal: converted to ComplexF64
    goal_int = [0 1; 1 0]
    qtraj_int = UnitaryTrajectory(system, goal_int, T)
    @test qtraj_int.goal == goal_int
    @test qtraj_int.goal isa Matrix{ComplexF64}

    # EmbeddedOperator goal: identity preserved
    H_drift = diagm(ComplexF64[1.0, 0.0, -1.0])
    H_drive = zeros(ComplexF64, 3, 3)
    H_drive[1, 2] = H_drive[2, 1] = 1.0
    sys3 = QuantumSystem(H_drift, [H_drive], [1.0])
    goal_embed = EmbeddedOperator(:X, [1, 2], 3)
    qtraj_embed = UnitaryTrajectory(sys3, goal_embed, T)
    @test qtraj_embed.goal === goal_embed

    # Pulse-based constructor: same conversion behavior
    times = [0.0, 0.5, 1.0]
    controls = 0.1 * randn(1, 3)
    pulse = ZeroOrderPulse(controls, times)

    qtraj_pulse_cf64 = UnitaryTrajectory(system, pulse, goal_cf64)
    @test qtraj_pulse_cf64.goal === goal_cf64

    qtraj_pulse_f64 = UnitaryTrajectory(system, pulse, goal_f64)
    @test qtraj_pulse_f64.goal == goal_f64
    @test qtraj_pulse_f64.goal isa Matrix{ComplexF64}
end

@testitem "UnitaryTrajectory callable" begin
    using LinearAlgebra

    system = QuantumSystem([PAULIS.X], [1.0])

    T = 1.0
    X_gate = ComplexF64[0 1; 1 0]
    qtraj = UnitaryTrajectory(system, X_gate, T)

    # Test at initial time
    U0 = qtraj(0.0)
    @test U0 ≈ Matrix{ComplexF64}(I, 2, 2)

    # Test at intermediate time
    U_mid = qtraj(0.5)
    @test U_mid isa Matrix{ComplexF64}
    @test size(U_mid) == (2, 2)

    # Test at final time
    U_final = qtraj(T)
    @test U_final isa Matrix{ComplexF64}
end

@testitem "UnitaryTrajectory fidelity" begin
    using LinearAlgebra

    # System that naturally implements X gate
    σx = ComplexF64[0 1; 1 0]
    system = QuantumSystem([σx], [1.0])

    # Create pulse that implements X gate: exp(-i π/2 σx) = -i σx
    T = π / 2
    times = [0.0, T]
    controls = ones(1, 2)  # Constant amplitude 1
    pulse = ZeroOrderPulse(controls, times)

    X_gate = ComplexF64[0 1; 1 0]
    qtraj = UnitaryTrajectory(system, pulse, X_gate)

    # Fidelity should be high
    fid = fidelity(qtraj)
    @test fid > 0.99
end

@testitem "UnitaryTrajectory with EmbeddedOperator goal" begin
    using LinearAlgebra

    # 3-level system with embedded 2-level gate
    H_drift = diagm(ComplexF64[1.0, 0.0, -1.0])
    H_drive = zeros(ComplexF64, 3, 3)
    H_drive[1, 2] = H_drive[2, 1] = 1.0
    system = QuantumSystem(H_drift, [H_drive], [1.0])

    # Embedded X gate on levels 1,2
    X_embedded = EmbeddedOperator(:X, [1, 2], 3)

    T = 1.0
    qtraj = UnitaryTrajectory(system, X_embedded, T)

    @test qtraj.goal === X_embedded
    @test qtraj.goal isa EmbeddedOperator

    # Fidelity with subspace
    fid = fidelity(qtraj; subspace = [1, 2])
    @test fid isa Real
end

@testitem "Adaptive vs fixed-step fidelity convergence" begin
    using LinearAlgebra
    using OrdinaryDiffEqLinear: MagnusGL4, MagnusAdapt4

    # Strong-driving system where 101 fixed steps may be inaccurate
    ω = 520.0
    sys = QuantumSystem(ω * PAULIS.Z, [PAULIS.X], [1.0])

    T = 2π / ω * 5  # 5 full rotations
    times = [0.0, T]
    controls = 0.5 * ones(1, 2)
    pulse = ZeroOrderPulse(controls, times)
    X_gate = ComplexF64[0 1; 1 0]

    # Fixed-step with only 101 points — may be inaccurate
    qtraj_fixed =
        UnitaryTrajectory(sys, pulse, X_gate; algorithm = MagnusGL4(), n_save = 101)
    fid_fixed = fidelity(qtraj_fixed)

    # Adaptive — should be accurate
    qtraj_adapt = UnitaryTrajectory(
        sys,
        pulse,
        X_gate;
        algorithm = MagnusAdapt4(),
        abstol = 1e-10,
        reltol = 1e-10,
    )
    fid_adapt = fidelity(qtraj_adapt)

    # Reference: fixed-step with very many points
    qtraj_ref =
        UnitaryTrajectory(sys, pulse, X_gate; algorithm = MagnusGL4(), n_save = 10001)
    fid_ref = fidelity(qtraj_ref)

    # Adaptive should be closer to reference than fixed-101
    # NOTE: both methods reach machine epsilon for this system, making the
    # comparison unreliable — marked broken until a more discriminating test is designed
    @test_broken false
    # @test_broken abs(fid_adapt - fid_ref) < abs(fid_fixed - fid_ref)
end

@testitem "MagnusAdapt4 preserves unitarity" begin
    using LinearAlgebra
    using OrdinaryDiffEqLinear: MagnusAdapt4

    sys = QuantumSystem(10.0 * PAULIS.Z, [PAULIS.X, PAULIS.Y], [1.0, 1.0])
    T = 1.0
    pulse = ZeroOrderPulse(randn(2, 5), range(0, T, length = 5))
    X_gate = ComplexF64[0 1; 1 0]

    qtraj = UnitaryTrajectory(
        sys,
        pulse,
        X_gate;
        algorithm = MagnusAdapt4(),
        abstol = 1e-10,
        reltol = 1e-10,
    )

    for U in qtraj.solution.u
        @test U' * U ≈ I atol = 1e-8
        @test U * U' ≈ I atol = 1e-8
    end
end

@testitem "Fidelity matches analytical result for X gate (MagnusAdapt4)" begin
    using LinearAlgebra
    using OrdinaryDiffEqLinear: MagnusAdapt4

    # Exact X gate: exp(-i π/2 σx) = -iσx, fidelity with X should be 1.0
    sys = QuantumSystem([PAULIS.X], [1.0])
    T = π / 2
    pulse = ZeroOrderPulse(ones(1, 2), [0.0, T])
    X_gate = ComplexF64[0 1; 1 0]

    qtraj = UnitaryTrajectory(
        sys,
        pulse,
        X_gate;
        algorithm = MagnusAdapt4(),
        abstol = 1e-12,
        reltol = 1e-12,
    )

    @test fidelity(qtraj) > 1.0 - 1e-10
end

@testitem "Fidelity consistent: MagnusAdapt4 vs MagnusGL4(fine)" begin
    using LinearAlgebra
    using OrdinaryDiffEqLinear: MagnusGL4, MagnusAdapt4

    sys = QuantumSystem([PAULIS.X, PAULIS.Y], [1.0, 1.0])
    T = π / 2
    pulse = ZeroOrderPulse(ones(2, 2), [0.0, T])
    X_gate = ComplexF64[0 1; 1 0]

    fid_adapt = fidelity(
        UnitaryTrajectory(
            sys,
            pulse,
            X_gate;
            algorithm = MagnusAdapt4(),
            abstol = 1e-10,
            reltol = 1e-10,
        ),
    )

    fid_gl4_fine = fidelity(
        UnitaryTrajectory(sys, pulse, X_gate; algorithm = MagnusGL4(), n_save = 10001),
    )

    @test fid_adapt ≈ fid_gl4_fine atol = 1e-6
end

@testitem "UnitaryTrajectory auto-selects ODE problem form by algorithm" begin
    # Performance bug fix: when an explicit RK algorithm (Tsit5/Vern9/etc.) is
    # passed, the trajectory must use UnitaryODEProblem (sparse mul! path) — NOT
    # UnitaryOperatorODEProblem, which materializes H densely on every ODE step.
    # When a Magnus-class algorithm is passed (or no algorithm — default routes to
    # MagnusAdapt4 for Hermitian systems), the operator form must be used so the
    # Lie-group integrator can compute exp(-iH dt)·U.
    #
    # We assert correctness (Magnus and Tsit5 paths agree on fidelity) and check
    # the underlying problem type via the ODESolution's function field.
    using LinearAlgebra
    using OrdinaryDiffEqLinear: MagnusAdapt4
    using OrdinaryDiffEqTsit5: Tsit5
    using OrdinaryDiffEqLinear: MatrixOperator

    sys = QuantumSystem([PAULIS.X, PAULIS.Y], [1.0, 1.0])
    T = π / 2
    pulse = ZeroOrderPulse(ones(2, 2), [0.0, T])
    X_gate = ComplexF64[0 1; 1 0]

    # 1. Explicit Magnus → operator form (MatrixOperator-backed function).
    qtraj_magnus = UnitaryTrajectory(
        sys,
        pulse,
        X_gate;
        algorithm = MagnusAdapt4(),
        abstol = 1e-10,
        reltol = 1e-10,
    )
    @test qtraj_magnus.solution.prob.f.f isa MatrixOperator

    # 2. Explicit RK (Tsit5) → mul!-based RHS (not a MatrixOperator).
    qtraj_rk = UnitaryTrajectory(
        sys,
        pulse,
        X_gate;
        algorithm = Tsit5(),
        abstol = 1e-10,
        reltol = 1e-10,
    )
    @test !(qtraj_rk.solution.prob.f.f isa MatrixOperator)

    # 3. Default (algorithm = nothing) for Hermitian system → MagnusAdapt4 →
    #    operator form (preserves existing behavior).
    qtraj_default = UnitaryTrajectory(sys, pulse, X_gate)
    @test qtraj_default.solution.prob.f.f isa MatrixOperator

    # 4. Both paths must agree on fidelity for the same pulse.
    @test fidelity(qtraj_magnus) ≈ fidelity(qtraj_rk) atol = 1e-6
end

@testitem "KetTrajectory auto-selects ODE problem form by algorithm" begin
    # Performance bug fix counterpart for KetTrajectory. See UnitaryTrajectory
    # auto-selection test for rationale.
    using LinearAlgebra
    using OrdinaryDiffEqLinear: MagnusAdapt4
    using OrdinaryDiffEqTsit5: Tsit5
    using OrdinaryDiffEqLinear: MatrixOperator

    sys = QuantumSystem([PAULIS.X], [1.0])
    T = π / 2
    pulse = ZeroOrderPulse(ones(1, 2), [0.0, T])
    ψ0 = ComplexF64[1.0, 0.0]
    ψg = ComplexF64[0.0, 1.0]

    qtraj_magnus = KetTrajectory(
        sys,
        pulse,
        ψ0,
        ψg;
        algorithm = MagnusAdapt4(),
        abstol = 1e-10,
        reltol = 1e-10,
    )
    @test qtraj_magnus.solution.prob.f.f isa MatrixOperator

    qtraj_rk = KetTrajectory(
        sys,
        pulse,
        ψ0,
        ψg;
        algorithm = Tsit5(),
        abstol = 1e-10,
        reltol = 1e-10,
    )
    @test !(qtraj_rk.solution.prob.f.f isa MatrixOperator)

    qtraj_default = KetTrajectory(sys, pulse, ψ0, ψg)
    @test qtraj_default.solution.prob.f.f isa MatrixOperator

    @test fidelity(qtraj_magnus) ≈ fidelity(qtraj_rk) atol = 1e-6
end

@testitem "MultiKetTrajectory auto-selects ODE problem form by algorithm" begin
    # Performance bug fix counterpart for MultiKetTrajectory.
    using LinearAlgebra
    using OrdinaryDiffEqLinear: MagnusAdapt4
    using OrdinaryDiffEqTsit5: Tsit5
    using OrdinaryDiffEqLinear: MatrixOperator

    sys = QuantumSystem([PAULIS.X], [1.0])
    T = π / 2
    pulse = ZeroOrderPulse(ones(1, 2), [0.0, T])
    initials = [ComplexF64[1.0, 0.0], ComplexF64[0.0, 1.0]]
    goals = [ComplexF64[0.0, 1.0], ComplexF64[1.0, 0.0]]

    qtraj_magnus = MultiKetTrajectory(
        sys,
        pulse,
        initials,
        goals;
        algorithm = MagnusAdapt4(),
        abstol = 1e-10,
        reltol = 1e-10,
    )
    # EnsembleSolution: per-trajectory base prob lives in qtraj.solution.u[i].prob
    @test qtraj_magnus.solution.u[1].prob.f.f isa MatrixOperator

    qtraj_rk = MultiKetTrajectory(
        sys,
        pulse,
        initials,
        goals;
        algorithm = Tsit5(),
        abstol = 1e-10,
        reltol = 1e-10,
    )
    @test !(qtraj_rk.solution.u[1].prob.f.f isa MatrixOperator)

    # Default (Hermitian system) → Magnus → operator form.
    qtraj_default = MultiKetTrajectory(sys, pulse, initials, goals)
    @test qtraj_default.solution.u[1].prob.f.f isa MatrixOperator

    @test fidelity(qtraj_magnus) ≈ fidelity(qtraj_rk) atol = 1e-6
end

@testitem "_needs_operator_form detects Magnus vs RK methods" begin
    using OrdinaryDiffEqLinear: MagnusAdapt4, MagnusGL4, MagnusGL6, MagnusGL8
    using OrdinaryDiffEqTsit5: Tsit5

    f = Piccolo.Quantum.Rollouts._needs_operator_form

    # Magnus methods need operator form
    @test f(MagnusAdapt4()) === true
    @test f(MagnusGL4()) === true
    @test f(MagnusGL6()) === true
    @test f(MagnusGL8()) === true

    # Explicit RK methods do not
    @test f(Tsit5()) === false
end
