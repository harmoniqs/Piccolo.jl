# ============================================================================ #
# MultiDensityTrajectory
# ============================================================================ #

"""
    MultiDensityTrajectory{P<:AbstractPulse, S} <: AbstractQuantumTrajectory{P}

Trajectory for multi-state open quantum system optimization with a shared pulse.
Evolves multiple density matrices under the same Lindblad master equation.

# Fields
- `system::OpenQuantumSystem`: The open quantum system
- `pulse::P`: The shared control pulse
- `initials::Vector{Matrix{ComplexF64}}`: Initial density matrices
- `goals::Vector{Matrix{ComplexF64}}`: Target density matrices
- `weights::Vector{Float64}`: Weights for fidelity calculation
- `solution::S`: Pre-computed ensemble solution

# Callable
`traj(t)` returns a vector of density matrices at time `t`.
`traj[i]` returns the i-th trajectory's solution.
"""
mutable struct MultiDensityTrajectory{P<:AbstractPulse,S} <: AbstractQuantumTrajectory{P}
    system::OpenQuantumSystem
    pulse::P
    initials::Vector{Matrix{ComplexF64}}
    goals::Vector{Matrix{ComplexF64}}
    weights::Vector{Float64}
    solution::S
end

"""
    MultiDensityTrajectory(system, pulse, initials, goals; weights=..., algorithm=Tsit5(), abstol=1e-8, reltol=1e-8, n_save=101)

Create a multi-density trajectory by solving multiple Lindblad master equations.

# Arguments
- `system::OpenQuantumSystem`: The open quantum system
- `pulse::AbstractPulse`: The shared control pulse
- `initials::Vector{<:AbstractMatrix}`: Initial density matrices
- `goals::Vector{<:AbstractMatrix}`: Target density matrices

# Keyword Arguments
- `weights`: Weights for fidelity (default: uniform)
- `algorithm`: ODE solver algorithm (default: Tsit5())
- `abstol`: Absolute tolerance for adaptive integration (default: 1e-8)
- `reltol`: Relative tolerance for adaptive integration (default: 1e-8)
- `n_save`: Number of output time points (default: 101)
"""
function MultiDensityTrajectory(
    system::OpenQuantumSystem,
    pulse::AbstractPulse,
    initials::Vector{<:AbstractMatrix{<:Number}},
    goals::Vector{<:AbstractMatrix{<:Number}};
    weights::AbstractVector{<:Real} = fill(1.0 / length(initials), length(initials)),
    algorithm = Tsit5(),
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
    n_save::Int = 101,
)
    @assert n_drives(pulse) == system.n_drives "Pulse has $(n_drives(pulse)) drives, system has $(system.n_drives)"
    @assert length(initials) == length(goals) == length(weights) "initials, goals, and weights must have same length"

    ρ0s = [Matrix{ComplexF64}(ρ) for ρ in initials]
    ρgs = [Matrix{ComplexF64}(ρ) for ρ in goals]
    ws = Vector{Float64}(weights)

    knot_times = get_knot_times(pulse)
    save_times = collect(range(0.0, duration(pulse), length = n_save))
    tstops = sort(unique(vcat(knot_times, save_times)))

    # Build ensemble problem
    dummy = zeros(ComplexF64, system.levels, system.levels)
    base_prob = DensityODEProblem(system, pulse, dummy, tstops)
    prob_func(prob, i, repeat) = remake(prob, u0 = ρ0s[i])
    ensemble_prob = EnsembleProblem(base_prob; prob_func = prob_func)
    sol = solve(
        ensemble_prob,
        algorithm;
        trajectories = length(initials),
        saveat = save_times,
        abstol = abstol,
        reltol = reltol,
    )

    return MultiDensityTrajectory{typeof(pulse),typeof(sol)}(system, pulse, ρ0s, ρgs, ws, sol)
end

"""
    MultiDensityTrajectory(system, initials, goals, T::Real; weights=..., drive_name=:u, algorithm=Tsit5(), abstol=1e-8, reltol=1e-8)

Convenience constructor that creates a zero pulse of duration T.

# Arguments
- `system::OpenQuantumSystem`: The open quantum system
- `initials::Vector{<:AbstractMatrix}`: Initial density matrices
- `goals::Vector{<:AbstractMatrix}`: Target density matrices
- `T::Real`: Duration of the pulse

# Keyword Arguments
- `weights`: Weights for fidelity (default: uniform)
- `drive_name::Symbol`: Name of the drive variable (default: `:u`)
- `algorithm`: ODE solver algorithm (default: Tsit5())
- `abstol`: Absolute tolerance (default: 1e-8)
- `reltol`: Relative tolerance (default: 1e-8)
"""
function MultiDensityTrajectory(
    system::OpenQuantumSystem,
    initials::Vector{<:AbstractMatrix{<:Number}},
    goals::Vector{<:AbstractMatrix{<:Number}},
    T::Real;
    weights::AbstractVector{<:Real} = fill(1.0 / length(initials), length(initials)),
    drive_name::Symbol = :u,
    algorithm = Tsit5(),
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
)
    times = [0.0, T]
    controls = vcat([rand(Uniform(b...), 1, length(times)) for b in system.drive_bounds]...)
    pulse = ZeroOrderPulse(controls, times; drive_name)
    return MultiDensityTrajectory(
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
(qtraj::MultiDensityTrajectory)(t::Real) = [sol(t) for sol in qtraj.solution]

# Indexing: get individual trajectory solution
Base.getindex(qtraj::MultiDensityTrajectory, i::Int) = qtraj.solution[i]
Base.length(qtraj::MultiDensityTrajectory) = length(qtraj.initials)

# ============================================================================ #
# Tests
# ============================================================================ #

@testitem "MultiDensityTrajectory construction" begin
    using LinearAlgebra

    # Simple 2-level open system
    L = ComplexF64[0.1 0.0; 0.0 0.0]
    system = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])

    # Create with duration
    T = 1.0
    ρ0s = [ComplexF64[1.0 0.0; 0.0 0.0], ComplexF64[0.0 0.0; 0.0 1.0]]
    ρgs = [ComplexF64[0.0 0.0; 0.0 1.0], ComplexF64[1.0 0.0; 0.0 0.0]]

    qtraj = MultiDensityTrajectory(system, ρ0s, ρgs, T)

    @test qtraj isa MultiDensityTrajectory
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

    qtraj2 = MultiDensityTrajectory(system, pulse, ρ0s, ρgs; weights = weights)

    @test qtraj2 isa MultiDensityTrajectory
    @test qtraj2.weights ≈ weights
end

@testitem "MultiDensityTrajectory callable and indexing" begin
    using LinearAlgebra

    L = ComplexF64[0.1 0.0; 0.0 0.0]
    system = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])

    T = 1.0
    ρ0s = [ComplexF64[1.0 0.0; 0.0 0.0], ComplexF64[0.0 0.0; 0.0 1.0]]
    ρgs = [ComplexF64[0.0 0.0; 0.0 1.0], ComplexF64[1.0 0.0; 0.0 0.0]]

    qtraj = MultiDensityTrajectory(system, ρ0s, ρgs, T)

    # Test length
    @test length(qtraj) == 2

    # Test callable - returns all states at time t
    states_0 = qtraj(0.0)
    @test length(states_0) == 2
    @test states_0[1] ≈ ρ0s[1]
    @test states_0[2] ≈ ρ0s[2]

    # Test indexing - returns individual solution
    sol1 = qtraj[1]
    @test sol1 isa ODESolution
    @test sol1(0.0) ≈ ρ0s[1]

    sol2 = qtraj[2]
    @test sol2(0.0) ≈ ρ0s[2]
end

@testitem "MultiDensityTrajectory fidelity" begin
    using LinearAlgebra

    # Open system without dissipation (to test fidelity calculation)
    σx = ComplexF64[0 1; 1 0]
    system = OpenQuantumSystem([σx], [1.0])

    # Pulse to rotate |0⟩ → |1⟩
    T = π / 2
    times = [0.0, T]
    controls = ones(1, 2)
    pulse = ZeroOrderPulse(controls, times)

    ρ0s = [ComplexF64[1.0 0.0; 0.0 0.0], ComplexF64[0.0 0.0; 0.0 1.0]]
    ρgs = [ComplexF64[0.0 0.0; 0.0 1.0], ComplexF64[1.0 0.0; 0.0 0.0]]

    qtraj = MultiDensityTrajectory(system, pulse, ρ0s, ρgs)

    # Both transfers should have high fidelity
    fid = fidelity(qtraj)
    @test fid > 0.99
end

@testitem "MultiDensityTrajectory state_names" begin
    using LinearAlgebra

    system = OpenQuantumSystem([PAULIS.X], [1.0])

    ρ0s = [
        ComplexF64[1.0 0.0; 0.0 0.0],
        ComplexF64[0.0 0.0; 0.0 1.0],
        ComplexF64[0.5 0.5; 0.5 0.5],
    ]
    ρgs = [
        ComplexF64[0.0 0.0; 0.0 1.0],
        ComplexF64[1.0 0.0; 0.0 0.0],
        ComplexF64[0.5 -0.5; -0.5 0.5],
    ]

    qtraj = MultiDensityTrajectory(system, ρ0s, ρgs, 1.0)

    @test state_name(qtraj) == :ρ⃗̃
    @test state_names(qtraj) == [:ρ⃗̃1, :ρ⃗̃2, :ρ⃗̃3]
end
