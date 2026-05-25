# ============================================================================ #
# KetTrajectory
# ============================================================================ #

"""
    KetTrajectory{P<:AbstractPulse, S<:ODESolution} <: AbstractQuantumTrajectory{P}

Trajectory for quantum state transfer. The ODE solution is computed at construction.

# Fields
- `system::QuantumSystem`: The quantum system
- `pulse::P`: The control pulse
- `initial::Vector{ComplexF64}`: Initial state |ψ₀⟩
- `goal::Vector{ComplexF64}`: Target state |ψ_goal⟩
- `solution::S`: Pre-computed ODE solution

# Callable
`traj(t)` returns the state at time `t` by interpolating the solution.
"""
mutable struct KetTrajectory{P<:AbstractPulse,S<:ODESolution} <:
               AbstractQuantumTrajectory{P}
    system::QuantumSystem
    pulse::P
    initial::Vector{ComplexF64}
    goal::Vector{ComplexF64}
    solution::S
end

"""
    KetTrajectory(system, pulse, initial, goal; algorithm=MagnusAdapt4(), abstol=1e-8, reltol=1e-8, n_save=101, progress=false, progress_steps=10)

Create a ket trajectory by solving the Schrödinger equation.

# Arguments
- `system::QuantumSystem`: The quantum system
- `pulse::AbstractPulse`: The control pulse
- `initial::Vector`: Initial state |ψ₀⟩
- `goal::Vector`: Target state |ψ_goal⟩

# Keyword Arguments
- `algorithm`: ODE solver algorithm (default: MagnusAdapt4())
- `abstol`: Absolute tolerance for adaptive integration (default: 1e-8)
- `reltol`: Relative tolerance for adaptive integration (default: 1e-8)
- `n_save`: Number of output time points (default: 101)
- `progress`: If true, emit ProgressLogging events during the ODE rollout (default: false).
  To see a progress bar in the terminal, install `TerminalLoggers` and set
  `global_logger(TerminalLogger())` before calling.
- `progress_steps`: Granularity of progress reports (default: 10, i.e. ~10% increments)
"""
function KetTrajectory(
    system::QuantumSystem,
    pulse::AbstractPulse,
    initial::AbstractVector{<:Number},
    goal::AbstractVector{<:Number};
    algorithm = nothing,
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
    n_save::Int = 101,
    progress::Bool = false,
    progress_steps::Int = 10,
)
    @assert n_drives(pulse) == system.n_drives "Pulse has $(n_drives(pulse)) drives, system has $(system.n_drives)"

    if isnothing(algorithm)
        algorithm = default_algorithm(system)
    end

    ψ0 = Vector{ComplexF64}(initial)
    ψg = Vector{ComplexF64}(goal)
    knot_times = get_knot_times(pulse)
    save_times = collect(range(0.0, duration(pulse), length = n_save))
    tstops = sort(unique(vcat(knot_times, save_times)))
    # Choose the ODE problem form to match the algorithm. Magnus-class methods
    # need H(t) materialized as a MatrixOperator (KetOperatorODEProblem). Explicit
    # RK methods (Tsit5, Vern9, DP8, ...) only need the dψ/dt = -iH(t)·ψ matvec
    # and should use KetODEProblem to avoid a per-step dense materialization
    # of H — a major hot-path penalty on high-dimensional sparse systems.
    prob = _needs_operator_form(algorithm) ?
        KetOperatorODEProblem(system, pulse, ψ0, tstops) :
        KetODEProblem(system, pulse, ψ0, tstops)
    sol = solve(prob, algorithm;
                saveat = save_times, abstol = abstol, reltol = reltol,
                progress = progress, progress_steps = progress_steps)

    return KetTrajectory{typeof(pulse),typeof(sol)}(system, pulse, ψ0, ψg, sol)
end

"""
    KetTrajectory(system, initial, goal, T::Real; drive_name=:u, algorithm=MagnusAdapt4(), abstol=1e-8, reltol=1e-8)

Convenience constructor that creates a zero pulse of duration T.

# Arguments
- `system::QuantumSystem`: The quantum system
- `initial::Vector`: Initial state |ψ₀⟩
- `goal::Vector`: Target state |ψ_goal⟩
- `T::Real`: Duration of the pulse

# Keyword Arguments
- `drive_name::Symbol`: Name of the drive variable (default: `:u`)
- `algorithm`: ODE solver algorithm (default: MagnusAdapt4())
- `abstol`: Absolute tolerance (default: 1e-8)
- `reltol`: Relative tolerance (default: 1e-8)
"""
function KetTrajectory(
    system::QuantumSystem,
    initial::AbstractVector{<:Number},
    goal::AbstractVector{<:Number},
    T::Real;
    drive_name::Symbol = :u,
    algorithm = nothing,
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
)
    times = [0.0, T]
    controls = vcat([rand(Uniform(b...), 1, length(times)) for b in system.drive_bounds]...)
    pulse = ZeroOrderPulse(controls, times; drive_name)
    return KetTrajectory(system, pulse, initial, goal; algorithm, abstol, reltol)
end

# Callable: sample solution at any time
(qtraj::KetTrajectory)(t::Real) = qtraj.solution(t)

# ============================================================================ #
# Tests
# ============================================================================ #

@testitem "KetTrajectory construction" begin
    using LinearAlgebra

    # Simple 2-level system
    system = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])

    # Create with duration
    T = 1.0
    ψ0 = ComplexF64[1.0, 0.0]
    ψg = ComplexF64[0.0, 1.0]
    qtraj = KetTrajectory(system, ψ0, ψg, T)

    @test qtraj isa KetTrajectory
    @test qtraj.system === system
    @test qtraj.initial ≈ ψ0
    @test qtraj.goal ≈ ψg

    # Create with explicit pulse
    times = [0.0, 0.5, 1.0]
    controls = 0.1 * randn(1, 3)
    pulse = ZeroOrderPulse(controls, times)
    qtraj2 = KetTrajectory(system, pulse, ψ0, ψg)

    @test qtraj2 isa KetTrajectory
    @test duration(qtraj2) ≈ 1.0
end

@testitem "KetTrajectory callable" begin
    using LinearAlgebra

    system = QuantumSystem([PAULIS.X], [1.0])

    T = 1.0
    ψ0 = ComplexF64[1.0, 0.0]
    ψg = ComplexF64[0.0, 1.0]
    qtraj = KetTrajectory(system, ψ0, ψg, T)

    # Test at initial time
    ψ_init = qtraj(0.0)
    @test ψ_init ≈ ψ0

    # Test at intermediate time
    ψ_mid = qtraj(0.5)
    @test ψ_mid isa Vector{ComplexF64}
    @test length(ψ_mid) == 2
    @test norm(ψ_mid) ≈ 1.0  # Should preserve normalization

    # Test at final time
    ψ_final = qtraj(T)
    @test norm(ψ_final) ≈ 1.0
end

@testitem "KetTrajectory fidelity" begin
    using LinearAlgebra

    # System with X drive
    σx = ComplexF64[0 1; 1 0]
    system = QuantumSystem([σx], [1.0])

    # Pulse to transfer |0⟩ → |1⟩: exp(-i π/2 σx) |0⟩ = -i|1⟩
    T = π / 2
    times = [0.0, T]
    controls = ones(1, 2)
    pulse = ZeroOrderPulse(controls, times)

    ψ0 = ComplexF64[1.0, 0.0]
    ψg = ComplexF64[0.0, 1.0]
    qtraj = KetTrajectory(system, pulse, ψ0, ψg)

    # Fidelity should be high
    fid = fidelity(qtraj)
    @test fid > 0.99
end
