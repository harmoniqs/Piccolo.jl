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
    PerKetStates

Per-ket protocol shim for a `RolloutStates` solution. Implements the two protocols the
`MultiKetTrajectory` solution consumers rely on:
- `(pk::PerKetStates)(t)` — the state at the nearest stored knot time ≤ `t`
  (used by `qtraj(t)` and `NamedTrajectory` conversion sampling)
- `pk.u` — vector of per-knot states, so `pk.u[end]` is the terminal state
  (used by `fidelity`)
"""
struct PerKetStates{M<:AbstractMatrix{ComplexF64}}
    mat::M                    # (ketdim, N_knots) — a view into the shared 3D array
    times::Vector{Float64}
end

(pk::PerKetStates)(t::Real) =
    getfield(pk, :mat)[:, max(searchsortedlast(getfield(pk, :times), t), 1)]

# NOTE: `.u` materializes copies of all N columns per access — fine at current call
# rates (fidelity reads one terminal state per call); do not use in hot loops.
Base.getproperty(pk::PerKetStates, s::Symbol) =
    s === :u ? [getfield(pk, :mat)[:, j] for j in axes(getfield(pk, :mat), 2)] :
    getfield(pk, s)

"""
    RolloutStates

No-ODE solution for `MultiKetTrajectory(...; rollout = :none)`. Holds knot states directly
(initialized by tiling the initial kets — a legitimate warm-start guess that downstream
`NamedTrajectory` sampling may consume) with a MUTABLE interior: a GPU rollout (e.g.
Piccolissimo's `gpu_rollout!`) overwrites `states` in place and sets `real_states = true`,
so the trajectory's solution type parameter never changes and `rollout!`'s
solution-reassignment path is never exercised.

`fidelity` on a trajectory whose `real_states` is still `false` warns (once) — a tiled
placeholder must never be silently reported as a physical fidelity. Direct `rollout!` on a
`RolloutStates` trajectory is an error (see `rollouts_extensions.jl`).
"""
# Concrete view type of one ket's (ketdim × N_knots) slice of the shared 3D array.
# Computed (not hand-written) so the shims are guaranteed to be stored as VIEWS —
# a Matrix-typed field would silently convert-and-COPY, leaving the shims stale
# after an in-place `states .= ...` refresh (review-caught failure mode).
const _KetStatesView = typeof(view(zeros(ComplexF64, 1, 1, 1), :, 1, :))

mutable struct RolloutStates
    states::Array{ComplexF64,3}   # (ketdim, K, N_knots) — shared storage
    times::Vector{Float64}
    real_states::Bool             # false until a real (GPU) rollout refreshes `states`
    u::Vector{PerKetStates{_KetStatesView}}

    function RolloutStates(
        states::Array{ComplexF64,3},
        times::Vector{Float64},
        real_states::Bool,
    )
        shims = [PerKetStates(view(states, :, k, :), times) for k in axes(states, 2)]
        return new(states, times, real_states, shims)
    end
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
    algorithm = nothing,
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
    n_save::Int = 101,
    rollout::Symbol = :cpu,
)
    @assert n_drives(pulse) == system.n_drives "Pulse has $(n_drives(pulse)) drives, system has $(system.n_drives)"
    @assert length(initials) == length(goals) == length(weights) "initials, goals, and weights must have same length"

    if rollout === :none
        # No-ODE fast path: knot states tiled from the initial kets (a readable
        # warm-start guess); refresh later with a GPU rollout. See `RolloutStates`.
        ψ0s = [Vector{ComplexF64}(ψ) for ψ in initials]
        ψgs = [Vector{ComplexF64}(ψ) for ψ in goals]
        ws = Vector{Float64}(weights)
        kt = collect(Float64, get_knot_times(pulse))
        states = Array{ComplexF64,3}(undef, system.levels, length(ψ0s), length(kt))
        for (k, ψ) in enumerate(ψ0s), j in eachindex(kt)
            states[:, k, j] = ψ
        end
        rs = RolloutStates(states, kt, false)
        return MultiKetTrajectory{typeof(pulse),RolloutStates}(
            system,
            pulse,
            ψ0s,
            ψgs,
            ws,
            rs,
        )
    elseif rollout !== :cpu
        throw(ArgumentError("rollout must be :cpu or :none, got $(repr(rollout))"))
    end

    if isnothing(algorithm)
        algorithm = default_algorithm(system)
    end

    ψ0s = [Vector{ComplexF64}(ψ) for ψ in initials]
    ψgs = [Vector{ComplexF64}(ψ) for ψ in goals]
    ws = Vector{Float64}(weights)

    knot_times = get_knot_times(pulse)
    save_times = collect(range(0.0, duration(pulse), length = n_save))
    tstops = sort(unique(vcat(knot_times, save_times)))

    # Build ensemble problem
    dummy = zeros(ComplexF64, system.levels)
    base_prob = KetOperatorODEProblem(system, pulse, dummy, tstops)
    prob_func(prob, i_or_ctx, _repeat = nothing) =
        remake(prob, u0 = ψ0s[_sim_index(i_or_ctx)])
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

Convenience constructor that creates a random pulse of duration T
(each drive sampled uniformly within its `drive_bounds`; pass `rng` for reproducibility).

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
    algorithm = nothing,
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
    rollout::Symbol = :cpu,
    rng::AbstractRNG = default_rng(),
)
    times = [0.0, T]
    controls =
        vcat([rand(rng, Uniform(b...), 1, length(times)) for b in system.drive_bounds]...)
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
        rollout,
    )
end

# Callable: sample all solutions at time t
(qtraj::MultiKetTrajectory)(t::Real) = [sol(t) for sol in qtraj.solution.u]

# Indexing: get individual trajectory solution
Base.getindex(qtraj::MultiKetTrajectory, i::Int) = qtraj.solution.u[i]
Base.length(qtraj::MultiKetTrajectory) = length(qtraj.initials)

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

@testitem "EnsembleSolution .u access returns ODESolution not scalar" begin
    using LinearAlgebra
    using SciMLBase: AbstractTimeseriesSolution

    sys = QuantumSystem(GATES[:Z], [GATES[:X]], [1.0])
    psi0 = ComplexF64[1.0, 0.0]
    psi1 = ComplexF64[0.0, 1.0]

    T = 1.0
    pulse = ZeroOrderPulse(0.5 * ones(1, 10), collect(range(0.0, T, length = 10)))
    qtraj = MultiKetTrajectory(sys, pulse, [psi0, psi1], [psi1, psi0])

    # .solution.u[i] must return an ODESolution-like object, not a scalar
    sol_1 = qtraj.solution.u[1]
    @test sol_1 isa AbstractTimeseriesSolution
    @test length(sol_1.u) > 1

    sol_2 = qtraj.solution.u[2]
    @test sol_2 isa AbstractTimeseriesSolution
    @test length(sol_2.u) > 1

    # length(.solution.u) should be trajectory count, not scalar element count
    @test length(qtraj.solution.u) == 2
end

@testitem "MultiKetTrajectory: rollout=:none tiles initials, readable, guarded" begin
    sys = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])
    ψ0 = [ComplexF64[1, 0], ComplexF64[0, 1]]
    ψg = [ComplexF64[0, 1], ComplexF64[1, 0]]
    qt = MultiKetTrajectory(sys, ψ0, ψg, 1.0; rollout = :none)
    @test qt isa MultiKetTrajectory
    @test qt.solution isa RolloutStates
    @test !qt.solution.real_states

    # placeholder reads SUCCEED (bootstrap requirement): tiled initial at every t
    for (k, ψ) in enumerate(ψ0)
        @test all(qt[k](t) ≈ ψ for t in (0.0, 0.5, 1.0))
    end
    @test qt(0.5) ≈ ψ0                      # trajectory callable samples all kets
    @test qt[1].u[end] ≈ ψ0[1]              # fidelity's terminal-state protocol

    # shims are VIEWS into the shared array — in-place refresh must be visible
    # (guards the silent view→Matrix copy failure mode)
    qt.solution.states[:, 1, :] .= ComplexF64[0, 1]
    @test qt[1](0.5) ≈ ComplexF64[0, 1]
    @test qt[1].u[end] ≈ ComplexF64[0, 1]

    # rollout! is a hard error on BOTH methods
    @test_throws ErrorException rollout!(qt, get_pulse(qt))
    @test_throws ErrorException rollout!(qt)

    # fidelity on a placeholder warns (once) but still returns a number
    @test_logs (:warn, r"placeholder") match_mode = :any begin
        f = fidelity(qt)
        @test f isa Float64
    end

    # bad kwarg value rejected; legacy default unchanged
    @test_throws ArgumentError MultiKetTrajectory(sys, ψ0, ψg, 1.0; rollout = :fast)
    qt2 = MultiKetTrajectory(sys, ψ0, ψg, 1.0)
    @test qt2 isa MultiKetTrajectory
    @test !(qt2.solution isa RolloutStates)
end
