# ============================================================================ #
# MultiKetTrajectory
# ============================================================================ #

"""
    MultiKetTrajectory{P<:AbstractPulse, S} <: AbstractQuantumTrajectory{P}

Trajectory for multi-state transfer with a shared pulse. Useful for state-to-state
problems with multiple initial/goal pairs.

# Fields
- `system::QuantumSystem`: The quantum system
- `pulse::P`: The shared control pulse
- `initials::Vector{Vector{ComplexF64}}`: Initial states
- `goals::Vector{Vector{ComplexF64}}`: Target states
- `weights::Vector{Float64}`: Weights for fidelity calculation
- `solution::S`: Pre-computed ensemble solution

# Callable
`traj(t)` returns a vector of states at time `t`.
`traj[i]` returns the i-th trajectory's solution.
"""
mutable struct MultiKetTrajectory{P<:AbstractPulse,S} <: AbstractQuantumTrajectory{P}
    system::QuantumSystem
    pulse::P
    initials::Vector{Vector{ComplexF64}}
    goals::Vector{Vector{ComplexF64}}
    weights::Vector{Float64}
    solution::S
end

"""
    MultiKetTrajectory(system, pulse, initials, goals; weights=..., algorithm=MagnusAdapt4(), abstol=1e-8, reltol=1e-8, n_save=101)

Create a multi-ket trajectory by solving multiple Schrödinger equations.

# Arguments
- `system::QuantumSystem`: The quantum system
- `pulse::AbstractPulse`: The shared control pulse
- `initials::Vector{Vector}`: Initial states
- `goals::Vector{Vector}`: Target states

# Keyword Arguments
- `weights`: Weights for fidelity (default: uniform)
- `algorithm`: ODE solver algorithm (default: MagnusAdapt4())
- `abstol`: Absolute tolerance for adaptive integration (default: 1e-8)
- `reltol`: Relative tolerance for adaptive integration (default: 1e-8)
- `n_save`: Number of output time points (default: 101)
"""
function MultiKetTrajectory(
    system::QuantumSystem,
    pulse::AbstractPulse,
    initials::Vector{<:AbstractVector{<:Number}},
    goals::Vector{<:AbstractVector{<:Number}};
    weights::AbstractVector{<:Real} = fill(1.0 / length(initials), length(initials)),
    algorithm = MagnusAdapt4(),
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
    n_save::Int = 101,
)
    @assert n_drives(pulse) == system.n_drives "Pulse has $(n_drives(pulse)) drives, system has $(system.n_drives)"
    @assert length(initials) == length(goals) == length(weights) "initials, goals, and weights must have same length"

    ψ0s = [Vector{ComplexF64}(ψ) for ψ in initials]
    ψgs = [Vector{ComplexF64}(ψ) for ψ in goals]
    ws = Vector{Float64}(weights)

    knot_times = get_knot_times(pulse)
    save_times = collect(range(0.0, duration(pulse), length = n_save))
    tstops = sort(unique(vcat(knot_times, save_times)))

    # Build ensemble problem
    dummy = zeros(ComplexF64, system.levels)
    base_prob = KetOperatorODEProblem(system, pulse, dummy, tstops)
    prob_func(prob, i, repeat) = remake(prob, u0 = ψ0s[i])
    ensemble_prob = EnsembleProblem(base_prob; prob_func = prob_func)
    sol = solve(
        ensemble_prob,
        algorithm;
        trajectories = length(initials),
        saveat = save_times,
        abstol = abstol,
        reltol = reltol,
    )

    return MultiKetTrajectory{typeof(pulse),typeof(sol)}(system, pulse, ψ0s, ψgs, ws, sol)
end

"""
    MultiKetTrajectory(system, initials, goals, T::Real; weights=..., drive_name=:u, algorithm=MagnusAdapt4(), abstol=1e-8, reltol=1e-8)

Convenience constructor that creates a zero pulse of duration T.

# Arguments
- `system::QuantumSystem`: The quantum system
- `initials::Vector{Vector}`: Initial states
- `goals::Vector{Vector}`: Target states
- `T::Real`: Duration of the pulse

# Keyword Arguments
- `weights`: Weights for fidelity (default: uniform)
- `drive_name::Symbol`: Name of the drive variable (default: `:u`)
- `algorithm`: ODE solver algorithm (default: MagnusAdapt4())
- `abstol`: Absolute tolerance (default: 1e-8)
- `reltol`: Relative tolerance (default: 1e-8)
"""
function MultiKetTrajectory(
    system::QuantumSystem,
    initials::Vector{<:AbstractVector{<:Number}},
    goals::Vector{<:AbstractVector{<:Number}},
    T::Real;
    weights::AbstractVector{<:Real} = fill(1.0 / length(initials), length(initials)),
    drive_name::Symbol = :u,
    algorithm = MagnusAdapt4(),
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
)
    times = [0.0, T]
    controls = vcat([rand(Uniform(b...), 1, length(times)) for b in system.drive_bounds]...)
    pulse = ZeroOrderPulse(controls, times; drive_name)
    return MultiKetTrajectory(
        system,
        pulse,
        initials,
        goals;
        weights,
        algorithm,
        abstol,
        reltol,
    )
end

# Callable: sample all solutions at time t
(traj::MultiKetTrajectory)(t::Real) = [sol(t) for sol in traj.solution]

# Indexing: get individual trajectory solution
Base.getindex(traj::MultiKetTrajectory, i::Int) = traj.solution[i]
Base.length(traj::MultiKetTrajectory) = length(traj.initials)

# ============================================================================ #
# Tests
# ============================================================================ #

@testitem "MultiKetTrajectory construction" begin
    using LinearAlgebra

    # Simple 2-level system
    system = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])

    # Create with duration
    T = 1.0
    initials = [ComplexF64[1.0, 0.0], ComplexF64[0.0, 1.0]]
    goals = [ComplexF64[0.0, 1.0], ComplexF64[1.0, 0.0]]

    qtraj = MultiKetTrajectory(system, initials, goals, T)

    @test qtraj isa MultiKetTrajectory
    @test qtraj.system === system
    @test length(qtraj.initials) == 2
    @test length(qtraj.goals) == 2
    @test length(qtraj.weights) == 2
    @test sum(qtraj.weights) ≈ 1.0  # Default uniform weights

    # Create with explicit pulse and weights
    times = [0.0, 0.5, 1.0]
    controls = 0.1 * randn(1, 3)
    pulse = ZeroOrderPulse(controls, times)
    weights = [0.7, 0.3]

    qtraj2 = MultiKetTrajectory(system, pulse, initials, goals; weights = weights)

    @test qtraj2 isa MultiKetTrajectory
    @test qtraj2.weights ≈ weights
end

@testitem "MultiKetTrajectory callable and indexing" begin
    using LinearAlgebra
    using OrdinaryDiffEqLinear

    system = QuantumSystem([PAULIS.X], [1.0])

    T = 1.0
    initials = [ComplexF64[1.0, 0.0], ComplexF64[0.0, 1.0]]
    goals = [ComplexF64[0.0, 1.0], ComplexF64[1.0, 0.0]]

    qtraj = MultiKetTrajectory(system, initials, goals, T)

    # Test length
    @test length(qtraj) == 2

    # Test callable - returns all states at time t
    states_0 = qtraj(0.0)
    @test length(states_0) == 2
    @test states_0[1] ≈ initials[1]
    @test states_0[2] ≈ initials[2]

    # Test indexing - returns individual solution
    sol1 = qtraj[1]
    @test sol1 isa ODESolution
    @test sol1(0.0) ≈ initials[1]

    sol2 = qtraj[2]
    @test sol2(0.0) ≈ initials[2]
end

@testitem "MultiKetTrajectory fidelity" begin
    using LinearAlgebra

    # System with X drive
    σx = ComplexF64[0 1; 1 0]
    system = QuantumSystem([σx], [1.0])

    # Pulse that swaps |0⟩ ↔ |1⟩
    T = π / 2
    times = [0.0, T]
    controls = ones(1, 2)
    pulse = ZeroOrderPulse(controls, times)

    initials = [ComplexF64[1.0, 0.0], ComplexF64[0.0, 1.0]]
    goals = [ComplexF64[0.0, 1.0], ComplexF64[1.0, 0.0]]

    qtraj = MultiKetTrajectory(system, pulse, initials, goals)

    # Both transfers should have high fidelity
    fid = fidelity(qtraj)
    @test fid > 0.99
end

@testitem "MultiKetTrajectory state_names" begin
    using LinearAlgebra

    system = QuantumSystem([PAULIS.X], [1.0])

    initials = [ComplexF64[1.0, 0.0], ComplexF64[0.0, 1.0], ComplexF64[1.0, 1.0] / √2]
    goals = [ComplexF64[0.0, 1.0], ComplexF64[1.0, 0.0], ComplexF64[1.0, -1.0] / √2]

    qtraj = MultiKetTrajectory(system, initials, goals, 1.0)

    @test state_name(qtraj) == :ψ̃
    @test state_names(qtraj) == [:ψ̃1, :ψ̃2, :ψ̃3]
end
