# ============================================================================ #
# UnitaryTrajectory
# ============================================================================ #

"""
    UnitaryTrajectory{P<:AbstractPulse, S<:ODESolution, G} <: AbstractQuantumTrajectory{P}

Trajectory for unitary gate synthesis. The ODE solution is computed at construction.

# Fields
- `system::QuantumSystem`: The quantum system
- `pulse::P`: The control pulse (stores drive_name)
- `initial::Matrix{ComplexF64}`: Initial unitary (default: identity)
- `goal::G`: Target unitary operator (AbstractPiccoloOperator or Matrix)
- `solution::S`: Pre-computed ODE solution

# Callable
`traj(t)` returns the unitary at time `t` by interpolating the solution.

# Conversion to NamedTrajectory
Use `NamedTrajectory(traj, N)` or `NamedTrajectory(traj, times)` for optimization.
"""
mutable struct UnitaryTrajectory{P<:AbstractPulse,S<:ODESolution,G} <:
               AbstractQuantumTrajectory{P}
    system::QuantumSystem
    pulse::P
    initial::Matrix{ComplexF64}
    goal::G
    solution::S
end

"""
    UnitaryTrajectory(system, pulse, goal; initial=I, algorithm=MagnusAdapt4(), abstol=1e-8, reltol=1e-8, n_save=101)

Create a unitary trajectory by solving the Schrödinger equation.

# Arguments
- `system::QuantumSystem`: The quantum system
- `pulse::AbstractPulse`: The control pulse
- `goal`: Target unitary (Matrix or AbstractPiccoloOperator)

# Keyword Arguments
- `initial`: Initial unitary (default: identity matrix)
- `algorithm`: ODE solver algorithm (default: MagnusAdapt4())
- `abstol`: Absolute tolerance for adaptive integration (default: 1e-8)
- `reltol`: Relative tolerance for adaptive integration (default: 1e-8)
- `n_save`: Number of output time points for plotting/interpolation (default: 101)
"""
function UnitaryTrajectory(
    system::QuantumSystem,
    pulse::AbstractPulse,
    goal::G;
    initial::AbstractMatrix{<:Number} = Matrix{ComplexF64}(I, system.levels, system.levels),
    algorithm = MagnusAdapt4(),
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
    n_save::Int = 101,
) where {G}
    @assert n_drives(pulse) == system.n_drives "Pulse has $(n_drives(pulse)) drives, system has $(system.n_drives)"

    U0 = Matrix{ComplexF64}(initial)
    save_times = collect(range(0.0, duration(pulse), length = n_save))
    prob = UnitaryOperatorODEProblem(system, pulse, save_times; U0 = U0)
    sol = solve(prob, algorithm; saveat = save_times, abstol = abstol, reltol = reltol)

    return UnitaryTrajectory{typeof(pulse),typeof(sol),G}(system, pulse, U0, goal, sol)
end

"""
    UnitaryTrajectory(system, goal, T::Real; drive_name=:u, algorithm=MagnusAdapt4(), abstol=1e-8, reltol=1e-8)

Convenience constructor that creates a zero pulse of duration T.

# Arguments
- `system::QuantumSystem`: The quantum system
- `goal`: Target unitary (Matrix or AbstractPiccoloOperator)
- `T::Real`: Duration of the pulse

# Keyword Arguments
- `drive_name::Symbol`: Name of the drive variable (default: `:u`)
- `algorithm`: ODE solver algorithm (default: MagnusAdapt4())
- `abstol`: Absolute tolerance (default: 1e-8)
- `reltol`: Relative tolerance (default: 1e-8)
"""
function UnitaryTrajectory(
    system::QuantumSystem,
    goal::G,
    T::Real;
    drive_name::Symbol = :u,
    algorithm = MagnusAdapt4(),
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
) where {G}
    times = [0.0, T]
    controls = vcat([rand(Uniform(b...), 1, length(times)) for b in system.drive_bounds]...)
    pulse = ZeroOrderPulse(controls, times; drive_name)
    return UnitaryTrajectory(system, pulse, goal; algorithm, abstol, reltol)
end

# Callable: sample solution at any time
(traj::UnitaryTrajectory)(t::Real) = traj.solution(t)

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
    ω = 50.0
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

    # Adaptive should match reference; fixed-101 may not
    @test abs(fid_adapt - fid_ref) < 1e-6
    @test abs(fid_fixed - fid_ref) > abs(fid_adapt - fid_ref)
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
        @test U' * U ≈ I atol=1e-8
        @test U * U' ≈ I atol=1e-8
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

    @test fid_adapt ≈ fid_gl4_fine atol=1e-6
end
