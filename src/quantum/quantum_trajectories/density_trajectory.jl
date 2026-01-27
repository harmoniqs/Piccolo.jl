# ============================================================================ #
# DensityTrajectory
# ============================================================================ #

"""
    DensityTrajectory{P<:AbstractPulse, S<:ODESolution} <: AbstractQuantumTrajectory{P}

Trajectory for open quantum systems (Lindblad dynamics).

# Fields
- `system::OpenQuantumSystem`: The open quantum system
- `pulse::P`: The control pulse
- `initial::Matrix{ComplexF64}`: Initial density matrix ρ₀
- `goal::Matrix{ComplexF64}`: Target density matrix ρ_goal
- `solution::S`: Pre-computed ODE solution

# Callable
`traj(t)` returns the density matrix at time `t`.
"""
mutable struct DensityTrajectory{P<:AbstractPulse,S<:ODESolution} <:
               AbstractQuantumTrajectory{P}
    system::OpenQuantumSystem
    pulse::P
    initial::Matrix{ComplexF64}
    goal::Matrix{ComplexF64}
    solution::S
end

"""
    DensityTrajectory(system, pulse, initial, goal; algorithm=Tsit5())

Create a density matrix trajectory by solving the Lindblad master equation.

# Arguments
- `system::OpenQuantumSystem`: The open quantum system
- `pulse::AbstractPulse`: The control pulse
- `initial::Matrix`: Initial density matrix ρ₀
- `goal::Matrix`: Target density matrix ρ_goal

# Keyword Arguments
- `algorithm`: ODE solver algorithm (default: Tsit5())
"""
function DensityTrajectory(
    system::OpenQuantumSystem,
    pulse::AbstractPulse,
    initial::AbstractMatrix{<:Number},
    goal::AbstractMatrix{<:Number};
    algorithm = Tsit5(),
)
    @assert n_drives(pulse) == system.n_drives "Pulse has $(n_drives(pulse)) drives, system has $(system.n_drives)"

    ρ0 = Matrix{ComplexF64}(initial)
    ρg = Matrix{ComplexF64}(goal)
    times = collect(range(0.0, duration(pulse), length = 101))
    prob = DensityODEProblem(system, pulse, ρ0, times)
    sol = solve(prob, algorithm; saveat = times)

    return DensityTrajectory{typeof(pulse),typeof(sol)}(system, pulse, ρ0, ρg, sol)
end

"""
    DensityTrajectory(system, initial, goal, T::Real; drive_name=:u, algorithm=Tsit5())

Convenience constructor that creates a zero pulse of duration T.
"""
function DensityTrajectory(
    system::OpenQuantumSystem,
    initial::AbstractMatrix{<:Number},
    goal::AbstractMatrix{<:Number},
    T::Real;
    drive_name::Symbol = :u,
    algorithm = Tsit5(),
)
    times = [0.0, T]
    controls = zeros(system.n_drives, 2)
    pulse = ZeroOrderPulse(controls, times; drive_name)
    return DensityTrajectory(system, pulse, initial, goal; algorithm)
end

# Callable: sample solution at any time
(traj::DensityTrajectory)(t::Real) = traj.solution(t)

# ============================================================================ #
# Tests
# ============================================================================ #

@testitem "DensityTrajectory construction" begin
    using LinearAlgebra

    # Simple 2-level open system
    L = ComplexF64[0.1 0.0; 0.0 0.0]  # Decay operator
    system = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])

    # Create with duration
    T = 1.0
    ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρg = ComplexF64[0.0 0.0; 0.0 1.0]

    qtraj = DensityTrajectory(system, ρ0, ρg, T)

    @test qtraj isa DensityTrajectory
    @test qtraj.system === system
    @test qtraj.initial ≈ ρ0
    @test qtraj.goal ≈ ρg

    # Create with explicit pulse
    times = [0.0, 0.5, 1.0]
    controls = 0.1 * randn(1, 3)
    pulse = ZeroOrderPulse(controls, times)

    qtraj2 = DensityTrajectory(system, pulse, ρ0, ρg)

    @test qtraj2 isa DensityTrajectory
    @test duration(qtraj2) ≈ 1.0
end

@testitem "DensityTrajectory callable" begin
    using LinearAlgebra

    L = ComplexF64[0.1 0.0; 0.0 0.0]
    system = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])

    T = 1.0
    ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρg = ComplexF64[0.0 0.0; 0.0 1.0]

    qtraj = DensityTrajectory(system, ρ0, ρg, T)

    # Test at initial time
    ρ_init = qtraj(0.0)
    @test ρ_init ≈ ρ0

    # Test at intermediate time
    ρ_mid = qtraj(0.5)
    @test ρ_mid isa Matrix{ComplexF64}
    @test size(ρ_mid) == (2, 2)
    @test real(tr(ρ_mid)) ≈ 1.0 atol=1e-6  # Trace should be preserved

    # Test at final time
    ρ_final = qtraj(T)
    @test real(tr(ρ_final)) ≈ 1.0 atol=1e-6
end

@testitem "DensityTrajectory fidelity" begin
    using LinearAlgebra

    # Open system without dissipation (to test fidelity calculation)
    σx = ComplexF64[0 1; 1 0]
    system = OpenQuantumSystem([σx], [1.0])

    # Start in pure state |0⟩, target pure state |1⟩
    ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρg = ComplexF64[0.0 0.0; 0.0 1.0]

    # Pulse to rotate |0⟩ → |1⟩
    T = π / 2
    times = [0.0, T]
    controls = ones(1, 2)
    pulse = ZeroOrderPulse(controls, times)

    qtraj = DensityTrajectory(system, pulse, ρ0, ρg)

    # Fidelity should be high
    fid = fidelity(qtraj)
    @test fid > 0.99
end
