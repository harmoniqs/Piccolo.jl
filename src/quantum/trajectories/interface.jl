# ============================================================================ #
# Common Interface
# ============================================================================ #

"""
    get_system(traj)

Get the quantum system from a trajectory.
"""
get_system(traj::AbstractQuantumTrajectory) = traj.system

"""
    get_pulse(traj)

Get the control pulse from a trajectory.
"""
get_pulse(traj::AbstractQuantumTrajectory) = traj.pulse

"""
    get_initial(traj)

Get the initial state/operator from a trajectory.
"""
get_initial(traj::UnitaryTrajectory) = traj.initial
get_initial(traj::KetTrajectory) = traj.initial
get_initial(traj::MultiKetTrajectory) = traj.initials
get_initial(traj::DensityTrajectory) = traj.initial

"""
    get_goal(traj)

Get the goal state/operator from a trajectory.
"""
get_goal(traj::UnitaryTrajectory) = traj.goal
get_goal(traj::KetTrajectory) = traj.goal
get_goal(traj::MultiKetTrajectory) = traj.goals
get_goal(traj::DensityTrajectory) = traj.goal

"""
    get_solution(traj)

Get the ODE solution from a trajectory.
"""
get_solution(traj::AbstractQuantumTrajectory) = traj.solution

# ============================================================================ #
# Fixed Name Accessors (for NamedTrajectory conversion)
# ============================================================================ #

"""
    state_name(::AbstractQuantumTrajectory)

Get the fixed state variable name for a trajectory type.
- `UnitaryTrajectory` → `:Ũ⃗`
- `KetTrajectory` → `:ψ̃`
- `MultiKetTrajectory` → `:ψ̃` (with index appended: `:ψ̃1`, `:ψ̃2`, etc.)
- `DensityTrajectory` → `:ρ⃗̃`
"""
state_name(::UnitaryTrajectory) = :Ũ⃗
state_name(::KetTrajectory) = :ψ̃
state_name(::MultiKetTrajectory) = :ψ̃  # prefix for :ψ̃1, :ψ̃2, etc.
state_name(::DensityTrajectory) = :ρ⃗̃

"""
    state_names(traj::MultiKetTrajectory)

Get all state names for an ensemble trajectory (`:ψ̃1`, `:ψ̃2`, etc.)
"""
function state_names(traj::MultiKetTrajectory)
    prefix = state_name(traj)
    return [Symbol(prefix, i) for i = 1:length(traj.initials)]
end

"""
    drive_name(traj::AbstractQuantumTrajectory)

Get the drive/control variable name from the trajectory's pulse.
"""
drive_name(traj::AbstractQuantumTrajectory) = drive_name(traj.pulse)

"""
    time_name(::AbstractQuantumTrajectory)

Get the time variable name (always `:t`).
"""
time_name(::AbstractQuantumTrajectory) = :t

"""
    timestep_name(::AbstractQuantumTrajectory)

Get the timestep variable name (always `:Δt`).
"""
timestep_name(::AbstractQuantumTrajectory) = :Δt

"""
    duration(traj)

Get the duration of a trajectory (from its pulse).
"""
duration(traj::AbstractQuantumTrajectory) = duration(traj.pulse)

# ============================================================================ #
# Rollout - Re-solve ODE with new pulse or different ODE parameters
# ============================================================================ #

"""
    rollout!(qtraj::UnitaryTrajectory, pulse::AbstractPulse; algorithm=MagnusGL4(), n_points=101)

Update quantum trajectory in-place with a new pulse by re-solving the ODE.
Mutates `qtraj.pulse` and `qtraj.solution`.

# Arguments
- `qtraj::UnitaryTrajectory`: The trajectory to update
- `pulse::AbstractPulse`: The new control pulse

# Keyword Arguments
- `algorithm`: ODE solver algorithm (default: `MagnusGL4()`)
- `n_points::Int`: Number of time points to sample (default: 101)

# Example
```julia
qtraj = UnitaryTrajectory(sys, old_pulse, goal)
rollout!(qtraj, new_pulse)  # Updates qtraj in-place
fid = fidelity(qtraj)  # Uses new solution
```

See also: [`rollout`](@ref)
"""
function Rollouts.rollout!(
    qtraj::UnitaryTrajectory,
    pulse::AbstractPulse;
    algorithm = MagnusGL4(),
    n_points::Int = 101,
)
    times = collect(range(0.0, duration(pulse), length = n_points))
    prob = UnitaryOperatorODEProblem(qtraj.system, pulse, times; U0 = qtraj.initial)
    sol = solve(prob, algorithm; saveat = times)

    qtraj.pulse = pulse
    qtraj.solution = sol
    return nothing
end

"""
    rollout!(qtraj::UnitaryTrajectory; algorithm=MagnusGL4(), n_points=101, kwargs...)

Update quantum trajectory in-place by re-solving with same pulse but different ODE parameters.
Mutates `qtraj.solution`.

Useful for comparing different solvers or tolerances.

# Keyword Arguments
- `algorithm`: ODE solver algorithm (default: `MagnusGL4()`)
- `n_points::Int`: Number of time points to sample (default: 101)
- Additional kwargs passed to `solve` (e.g., `abstol`, `reltol`)

# Example
```julia
qtraj = UnitaryTrajectory(sys, pulse, goal)

# Compare Magnus vs Runge-Kutta
rollout!(qtraj; algorithm=Tsit5(), abstol=1e-10)
fid_rk = fidelity(qtraj)
```

See also: [`rollout`](@ref)
"""
function Rollouts.rollout!(
    qtraj::UnitaryTrajectory;
    algorithm = MagnusGL4(),
    n_points::Int = 101,
    kwargs...,
)
    times = collect(range(0.0, duration(qtraj.pulse), length = n_points))
    prob = UnitaryOperatorODEProblem(qtraj.system, qtraj.pulse, times; U0 = qtraj.initial)
    sol = solve(prob, algorithm; saveat = times, kwargs...)

    qtraj.solution = sol
    return nothing
end

"""
    rollout!(qtraj::KetTrajectory, pulse::AbstractPulse; algorithm=MagnusGL4(), n_points=101)

Update ket trajectory in-place with a new pulse.
See [`rollout!(::UnitaryTrajectory, ::AbstractPulse)`](@ref) for details.
"""
function Rollouts.rollout!(
    qtraj::KetTrajectory,
    pulse::AbstractPulse;
    algorithm = MagnusGL4(),
    n_points::Int = 101,
)
    times = collect(range(0.0, duration(pulse), length = n_points))
    prob = KetOperatorODEProblem(qtraj.system, pulse, qtraj.initial, times)
    sol = solve(prob, algorithm; saveat = times)

    qtraj.pulse = pulse
    qtraj.solution = sol
    return nothing
end

"""
    rollout!(qtraj::KetTrajectory; algorithm=MagnusGL4(), n_points=101, kwargs...)

Update ket trajectory in-place with same pulse but different ODE parameters.
See [`rollout!(::UnitaryTrajectory; kwargs...)`](@ref) for details.
"""
function Rollouts.rollout!(
    qtraj::KetTrajectory;
    algorithm = MagnusGL4(),
    n_points::Int = 101,
    kwargs...,
)
    times = collect(range(0.0, duration(qtraj.pulse), length = n_points))
    prob = KetOperatorODEProblem(qtraj.system, qtraj.pulse, qtraj.initial, times)
    sol = solve(prob, algorithm; saveat = times, kwargs...)

    qtraj.solution = sol
    return nothing
end

"""
    rollout!(qtraj::MultiKetTrajectory, pulse::AbstractPulse; algorithm=MagnusGL4(), n_points=101)

Update multi-ket trajectory in-place with a new pulse.
See [`rollout!(::UnitaryTrajectory, ::AbstractPulse)`](@ref) for details.
"""
function Rollouts.rollout!(
    qtraj::MultiKetTrajectory,
    pulse::AbstractPulse;
    algorithm = MagnusGL4(),
    n_points::Int = 101,
)
    times = collect(range(0.0, duration(pulse), length = n_points))

    # Build ensemble problem
    dummy = zeros(ComplexF64, qtraj.system.levels)
    base_prob = KetOperatorODEProblem(qtraj.system, pulse, dummy, times)
    prob_func(prob, i, repeat) = remake(prob, u0 = qtraj.initials[i])
    ensemble_prob = EnsembleProblem(base_prob; prob_func = prob_func)
    sol = solve(
        ensemble_prob,
        algorithm;
        trajectories = length(qtraj.initials),
        saveat = times,
    )

    qtraj.pulse = pulse
    qtraj.solution = sol
    return nothing
end

"""
    rollout!(qtraj::MultiKetTrajectory; algorithm=MagnusGL4(), n_points=101, kwargs...)

Update multi-ket trajectory in-place with same pulse but different ODE parameters.
See [`rollout!(::UnitaryTrajectory; kwargs...)`](@ref) for details.
"""
function Rollouts.rollout!(
    qtraj::MultiKetTrajectory;
    algorithm = MagnusGL4(),
    n_points::Int = 101,
    kwargs...,
)
    times = collect(range(0.0, duration(qtraj.pulse), length = n_points))

    # Build ensemble problem
    dummy = zeros(ComplexF64, qtraj.system.levels)
    base_prob = KetOperatorODEProblem(qtraj.system, qtraj.pulse, dummy, times)
    prob_func(prob, i, repeat) = remake(prob, u0 = qtraj.initials[i])
    ensemble_prob = EnsembleProblem(base_prob; prob_func = prob_func)
    sol = solve(
        ensemble_prob,
        algorithm;
        trajectories = length(qtraj.initials),
        saveat = times,
        kwargs...,
    )

    qtraj.solution = sol
    return nothing
end

"""
    rollout!(qtraj::DensityTrajectory, pulse::AbstractPulse; algorithm=Tsit5(), n_points=101)

Update density trajectory in-place with a new pulse.
Note: Default algorithm is `Tsit5()` since density evolution uses standard ODE solvers.
See [`rollout!(::UnitaryTrajectory, ::AbstractPulse)`](@ref) for details.
"""
function Rollouts.rollout!(
    qtraj::DensityTrajectory,
    pulse::AbstractPulse;
    algorithm = Tsit5(),
    n_points::Int = 101,
)
    times = collect(range(0.0, duration(pulse), length = n_points))
    prob = DensityODEProblem(qtraj.system, pulse, qtraj.initial, times)
    sol = solve(prob, algorithm; saveat = times)

    qtraj.pulse = pulse
    qtraj.solution = sol
    return nothing
end

"""
    rollout!(qtraj::DensityTrajectory; algorithm=Tsit5(), n_points=101, kwargs...)

Update density trajectory in-place with same pulse but different ODE parameters.
Note: Default algorithm is `Tsit5()` since density evolution uses standard ODE solvers.
See [`rollout!(::UnitaryTrajectory; kwargs...)`](@ref) for details.
"""
function Rollouts.rollout!(
    qtraj::DensityTrajectory;
    algorithm = Tsit5(),
    n_points::Int = 101,
    kwargs...,
)
    times = collect(range(0.0, duration(qtraj.pulse), length = n_points))
    prob = DensityODEProblem(qtraj.system, qtraj.pulse, qtraj.initial, times)
    sol = solve(prob, algorithm; saveat = times, kwargs...)

    qtraj.solution = sol
    return nothing
end

"""
    rollout!(qtraj::SamplingTrajectory, pulse::AbstractPulse; algorithm=MagnusGL4(), n_points=101)

Update sampling trajectory's base trajectory in-place with a new pulse.
Delegates to the base trajectory's rollout! method.
"""
function Rollouts.rollout!(
    qtraj::SamplingTrajectory,
    pulse::AbstractPulse;
    algorithm = MagnusGL4(),
    n_points::Int = 101,
)
    rollout!(qtraj.base_trajectory, pulse; algorithm = algorithm, n_points = n_points)
    return nothing
end

"""
    rollout!(qtraj::SamplingTrajectory; algorithm=MagnusGL4(), n_points=101, kwargs...)

Update sampling trajectory's base trajectory in-place with new ODE parameters.
Delegates to the base trajectory's rollout! method.
"""
function Rollouts.rollout!(
    qtraj::SamplingTrajectory;
    algorithm = MagnusGL4(),
    n_points::Int = 101,
    kwargs...,
)
    rollout!(qtraj.base_trajectory; algorithm = algorithm, n_points = n_points, kwargs...)
    return nothing
end

"""
    rollout(qtraj::UnitaryTrajectory, pulse::AbstractPulse; algorithm=MagnusGL4(), n_points=101)

Create a new quantum trajectory by rolling out a new pulse through the system.
Returns a new UnitaryTrajectory with the updated pulse and solution.

# Arguments
- `qtraj::UnitaryTrajectory`: The base trajectory (provides system, initial, goal)
- `pulse::AbstractPulse`: The new control pulse to roll out

# Keyword Arguments
- `algorithm`: ODE solver algorithm (default: `MagnusGL4()`)
- `n_points::Int`: Number of time points to sample (default: 101)

# Example
```julia
qtraj = UnitaryTrajectory(sys, old_pulse, goal)

# Roll out a new pulse
qtraj_new = rollout(qtraj, new_pulse)

# Check fidelity
fid = fidelity(qtraj_new)
```

See also: [`extract_pulse`](@ref), [`rollout!`](@ref), [`fidelity`](@ref)
"""
function Rollouts.rollout(
    qtraj::UnitaryTrajectory,
    pulse::AbstractPulse;
    algorithm = MagnusGL4(),
    n_points::Int = 101,
)
    times = collect(range(0.0, duration(pulse), length = n_points))
    prob = UnitaryOperatorODEProblem(qtraj.system, pulse, times; U0 = qtraj.initial)
    sol = solve(prob, algorithm; saveat = times)
    return UnitaryTrajectory(qtraj.system, pulse, qtraj.initial, qtraj.goal, sol)
end

"""
    rollout(qtraj::KetTrajectory, pulse::AbstractPulse; algorithm=MagnusGL4(), n_points=101)

Create a new ket trajectory by rolling out a new pulse.
See [`rollout(::UnitaryTrajectory, ::AbstractPulse)`](@ref) for details.
"""
function Rollouts.rollout(
    qtraj::KetTrajectory,
    pulse::AbstractPulse;
    algorithm = MagnusGL4(),
    n_points::Int = 101,
)
    times = collect(range(0.0, duration(pulse), length = n_points))
    prob = KetOperatorODEProblem(qtraj.system, pulse, qtraj.initial, times)
    sol = solve(prob, algorithm; saveat = times)
    return KetTrajectory(qtraj.system, pulse, qtraj.initial, qtraj.goal, sol)
end

"""
    rollout(qtraj::MultiKetTrajectory, pulse::AbstractPulse; algorithm=MagnusGL4(), n_points=101)

Create a new multi-ket trajectory by rolling out a new pulse.
See [`rollout(::UnitaryTrajectory, ::AbstractPulse)`](@ref) for details.
"""
function Rollouts.rollout(
    qtraj::MultiKetTrajectory,
    pulse::AbstractPulse;
    algorithm = MagnusGL4(),
    n_points::Int = 101,
)
    times = collect(range(0.0, duration(pulse), length = n_points))

    # Build ensemble problem
    dummy = zeros(ComplexF64, qtraj.system.levels)
    base_prob = KetOperatorODEProblem(qtraj.system, pulse, dummy, times)
    prob_func(prob, i, repeat) = remake(prob, u0 = qtraj.initials[i])
    ensemble_prob = EnsembleProblem(base_prob; prob_func = prob_func)
    sol = solve(
        ensemble_prob,
        algorithm;
        trajectories = length(qtraj.initials),
        saveat = times,
    )

    return MultiKetTrajectory(
        qtraj.system,
        pulse,
        qtraj.initials,
        qtraj.goals,
        qtraj.weights,
        sol,
    )
end

"""
    rollout(qtraj::DensityTrajectory, pulse::AbstractPulse; algorithm=Tsit5(), n_points=101)

Create a new density trajectory by rolling out a new pulse.
Note: Default algorithm is `Tsit5()` since density evolution uses standard ODE solvers.
See [`rollout(::UnitaryTrajectory, ::AbstractPulse)`](@ref) for details.
"""
function Rollouts.rollout(
    qtraj::DensityTrajectory,
    pulse::AbstractPulse;
    algorithm = Tsit5(),
    n_points::Int = 101,
)
    times = collect(range(0.0, duration(pulse), length = n_points))
    prob = DensityODEProblem(qtraj.system, pulse, qtraj.initial, times)
    sol = solve(prob, algorithm; saveat = times)
    return DensityTrajectory(qtraj.system, pulse, qtraj.initial, qtraj.goal, sol)
end

# Rollout with same pulse, different ODE parameters (non-mutating)

"""
    rollout(qtraj::UnitaryTrajectory; algorithm=MagnusGL4(), n_points=101, kwargs...)

Re-solve the trajectory with the same pulse but different ODE parameters.
Returns a new UnitaryTrajectory with the updated solution.

Useful for comparing different solvers or tolerances.

# Keyword Arguments
- `algorithm`: ODE solver algorithm (default: `MagnusGL4()`)
- `n_points::Int`: Number of time points to sample (default: 101)
- Additional kwargs passed to `solve` (e.g., `abstol`, `reltol`)

# Example
```julia
qtraj = UnitaryTrajectory(sys, pulse, goal)

# Compare Magnus vs Runge-Kutta
qtraj_rk = rollout(qtraj; algorithm=Tsit5(), abstol=1e-10)
fid_magnus = fidelity(qtraj)
fid_rk = fidelity(qtraj_rk)
```

See also: [`rollout!`](@ref)
"""
function Rollouts.rollout(
    qtraj::UnitaryTrajectory;
    algorithm = MagnusGL4(),
    n_points::Int = 101,
    kwargs...,
)
    times = collect(range(0.0, duration(qtraj.pulse), length = n_points))
    prob = UnitaryOperatorODEProblem(qtraj.system, qtraj.pulse, times; U0 = qtraj.initial)
    sol = solve(prob, algorithm; saveat = times, kwargs...)
    return UnitaryTrajectory(qtraj.system, qtraj.pulse, qtraj.initial, qtraj.goal, sol)
end

"""
    rollout(qtraj::KetTrajectory; algorithm=MagnusGL4(), n_points=101, kwargs...)

Re-solve ket trajectory with same pulse but different ODE parameters.
See [`rollout(::UnitaryTrajectory; kwargs...)`](@ref) for details.
"""
function Rollouts.rollout(
    qtraj::KetTrajectory;
    algorithm = MagnusGL4(),
    n_points::Int = 101,
    kwargs...,
)
    times = collect(range(0.0, duration(qtraj.pulse), length = n_points))
    prob = KetOperatorODEProblem(qtraj.system, qtraj.pulse, qtraj.initial, times)
    sol = solve(prob, algorithm; saveat = times, kwargs...)
    return KetTrajectory(qtraj.system, qtraj.pulse, qtraj.initial, qtraj.goal, sol)
end

"""
    rollout(qtraj::MultiKetTrajectory; algorithm=MagnusGL4(), n_points=101, kwargs...)

Re-solve multi-ket trajectory with same pulse but different ODE parameters.
See [`rollout(::UnitaryTrajectory; kwargs...)`](@ref) for details.
"""
function Rollouts.rollout(
    qtraj::MultiKetTrajectory;
    algorithm = MagnusGL4(),
    n_points::Int = 101,
    kwargs...,
)
    times = collect(range(0.0, duration(qtraj.pulse), length = n_points))

    # Build ensemble problem
    dummy = zeros(ComplexF64, qtraj.system.levels)
    base_prob = KetOperatorODEProblem(qtraj.system, qtraj.pulse, dummy, times)
    prob_func(prob, i, repeat) = remake(prob, u0 = qtraj.initials[i])
    ensemble_prob = EnsembleProblem(base_prob; prob_func = prob_func)
    sol = solve(
        ensemble_prob,
        algorithm;
        trajectories = length(qtraj.initials),
        saveat = times,
        kwargs...,
    )

    return MultiKetTrajectory(
        qtraj.system,
        qtraj.pulse,
        qtraj.initials,
        qtraj.goals,
        qtraj.weights,
        sol,
    )
end

"""
    rollout(qtraj::DensityTrajectory; algorithm=Tsit5(), n_points=101, kwargs...)

Re-solve density trajectory with same pulse but different ODE parameters.
Note: Default algorithm is `Tsit5()` since density evolution uses standard ODE solvers.
See [`rollout(::UnitaryTrajectory; kwargs...)`](@ref) for details.
"""
function Rollouts.rollout(
    qtraj::DensityTrajectory;
    algorithm = Tsit5(),
    n_points::Int = 101,
    kwargs...,
)
    times = collect(range(0.0, duration(qtraj.pulse), length = n_points))
    prob = DensityODEProblem(qtraj.system, qtraj.pulse, qtraj.initial, times)
    sol = solve(prob, algorithm; saveat = times, kwargs...)
    return DensityTrajectory(qtraj.system, qtraj.pulse, qtraj.initial, qtraj.goal, sol)
end

# ============================================================================ #
# Fidelity (extending Rollouts.fidelity)
# ============================================================================ #

"""
    fidelity(traj::UnitaryTrajectory; subspace=nothing)

Compute the fidelity between the final unitary and the goal.
"""
function Rollouts.fidelity(
    traj::UnitaryTrajectory;
    subspace::Union{Nothing,AbstractVector{Int}} = nothing,
)
    U_final = traj.solution.u[end]
    U_goal = traj.goal isa EmbeddedOperator ? traj.goal.operator : traj.goal
    if isnothing(subspace)
        return unitary_fidelity(U_final, U_goal)
    else
        return unitary_fidelity(U_final, U_goal; subspace = subspace)
    end
end

"""
    fidelity(traj::KetTrajectory)

Compute the fidelity between the final state and the goal.
"""
function Rollouts.fidelity(traj::KetTrajectory)
    ψ_final = traj.solution.u[end]
    return abs2(ψ_final' * traj.goal)
end

"""
    fidelity(traj::MultiKetTrajectory)

Compute the weighted average fidelity across all state transfers.
"""
function Rollouts.fidelity(traj::MultiKetTrajectory)
    fids = map(zip(traj.solution, traj.goals)) do (sol, goal)
        abs2(sol.u[end]' * goal)
    end
    return sum(traj.weights .* fids)
end

"""
    fidelity(traj::DensityTrajectory)

Compute the fidelity between the final density matrix and the goal.
Uses trace fidelity: F = tr(ρ_final * ρ_goal)
"""
function Rollouts.fidelity(traj::DensityTrajectory)
    ρ_final = traj.solution.u[end]
    return real(tr(ρ_final * traj.goal))
end

# ============================================================================ #
# Tests
# ============================================================================ #

@testitem "Common interface - getters" begin
    using LinearAlgebra

    # UnitaryTrajectory
    system = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])

    X_gate = ComplexF64[0 1; 1 0]
    qtraj = UnitaryTrajectory(system, X_gate, 1.0)

    @test get_system(qtraj) === system
    @test get_pulse(qtraj) isa AbstractPulse
    @test get_initial(qtraj) ≈ Matrix{ComplexF64}(I, 2, 2)
    @test get_goal(qtraj) === X_gate
    @test duration(qtraj) ≈ 1.0

    # KetTrajectory
    ψ0 = ComplexF64[1.0, 0.0]
    ψg = ComplexF64[0.0, 1.0]
    qtraj_ket = KetTrajectory(system, ψ0, ψg, 1.0)

    @test get_system(qtraj_ket) === system
    @test get_initial(qtraj_ket) ≈ ψ0
    @test get_goal(qtraj_ket) ≈ ψg

    # MultiKetTrajectory
    initials = [ψ0, ψg]
    goals = [ψg, ψ0]
    qtraj_ens = MultiKetTrajectory(system, initials, goals, 1.0)

    @test get_system(qtraj_ens) === system
    @test get_initial(qtraj_ens) == qtraj_ens.initials
    @test get_goal(qtraj_ens) == qtraj_ens.goals
end

@testitem "Common interface - name accessors" begin
    using LinearAlgebra

    system = QuantumSystem([PAULIS.X], [1.0])

    # Test state_name for each trajectory type
    qtraj_u = UnitaryTrajectory(system, ComplexF64[0 1; 1 0], 1.0)
    @test state_name(qtraj_u) == :Ũ⃗

    qtraj_k = KetTrajectory(system, ComplexF64[1, 0], ComplexF64[0, 1], 1.0)
    @test state_name(qtraj_k) == :ψ̃

    qtraj_e = MultiKetTrajectory(system, [ComplexF64[1, 0]], [ComplexF64[0, 1]], 1.0)
    @test state_name(qtraj_e) == :ψ̃
    @test state_names(qtraj_e) == [:ψ̃1]

    # Test drive_name propagation
    times = [0.0, 1.0]
    pulse = ZeroOrderPulse(zeros(1, 2), times; drive_name = :control)
    qtraj_custom = UnitaryTrajectory(system, pulse, ComplexF64[0 1; 1 0])
    @test drive_name(qtraj_custom) == :control

    # Test time_name and timestep_name (always fixed)
    @test time_name(qtraj_u) == :t
    @test timestep_name(qtraj_u) == :Δt
end

@testitem "Interface - DensityTrajectory getters" begin
    using LinearAlgebra

    L = ComplexF64[0.1 0.0; 0.0 0.0]
    system = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])

    ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρg = ComplexF64[0.0 0.0; 0.0 1.0]
    qtraj = DensityTrajectory(system, ρ0, ρg, 1.0)

    @test get_system(qtraj) === system
    @test get_initial(qtraj) ≈ ρ0
    @test get_goal(qtraj) ≈ ρg
    @test state_name(qtraj) == :ρ⃗̃
    @test duration(qtraj) ≈ 1.0
end
@testitem "rollout - UnitaryTrajectory" begin
    using LinearAlgebra
    using OrdinaryDiffEqLinear: MagnusGL4
    using OrdinaryDiffEqTsit5: Tsit5

    system = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])
    X_gate = ComplexF64[0 1; 1 0]

    # Create trajectory with initial pulse (1 drive × 2 timesteps)
    pulse1 = ZeroOrderPulse([0.5 0.5], [0.0, 1.0])
    qtraj = UnitaryTrajectory(system, pulse1, X_gate)
    @test length(qtraj.solution.u) == 101
    fid1 = fidelity(qtraj)

    # Roll out a new pulse
    pulse2 = ZeroOrderPulse([0.8 0.8], [0.0, 1.0])
    qtraj_new = rollout(qtraj, pulse2)

    @test length(qtraj_new.solution.u) == 101
    @test qtraj_new.system === qtraj.system
    @test qtraj_new.pulse === pulse2
    @test qtraj_new.goal === qtraj.goal

    # Should have different fidelity (different pulse)
    fid2 = fidelity(qtraj_new)
    @test fid2 != fid1

    # Roll out with higher resolution
    qtraj_fine = rollout(qtraj, pulse2; n_points = 501)
    @test length(qtraj_fine.solution.u) == 501

    # Roll out with different algorithm
    qtraj_rk = rollout(qtraj, pulse2; algorithm = Tsit5())
    @test length(qtraj_rk.solution.u) == 101
end

@testitem "rollout - KetTrajectory" begin
    using LinearAlgebra
    using OrdinaryDiffEqLinear: MagnusGL4

    system = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])
    ψ0 = ComplexF64[1.0, 0.0]
    ψg = ComplexF64[0.0, 1.0]

    # Create trajectory
    pulse1 = ZeroOrderPulse([0.5 0.5], [0.0, 1.0])
    qtraj = KetTrajectory(system, pulse1, ψ0, ψg)

    # Roll out new pulse
    pulse2 = ZeroOrderPulse([0.8 0.8], [0.0, 1.0])
    qtraj_new = rollout(qtraj, pulse2; n_points = 201)

    @test length(qtraj_new.solution.u) == 201
    @test qtraj_new.pulse === pulse2
end

@testitem "rollout - MultiKetTrajectory" begin
    using LinearAlgebra

    system = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])
    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]

    initials = [ψ0, ψ1]
    goals = [ψ1, ψ0]

    pulse1 = ZeroOrderPulse([0.5 0.5], [0.0, 1.0])
    qtraj = MultiKetTrajectory(system, pulse1, initials, goals)

    # Roll out new pulse
    pulse2 = ZeroOrderPulse([0.8 0.8], [0.0, 1.0])
    qtraj_new = rollout(qtraj, pulse2)

    @test length(qtraj_new.solution) == 2
    @test qtraj_new.pulse === pulse2
end

@testitem "rollout - DensityTrajectory" begin
    using LinearAlgebra
    using OrdinaryDiffEqTsit5: Tsit5

    L = ComplexF64[0.1 0.0; 0.0 0.0]
    system = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])

    ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρg = ComplexF64[0.0 0.0; 0.0 1.0]

    pulse1 = ZeroOrderPulse([0.5 0.5], [0.0, 1.0])
    qtraj = DensityTrajectory(system, pulse1, ρ0, ρg)

    # Roll out new pulse
    pulse2 = ZeroOrderPulse([0.8 0.8], [0.0, 1.0])
    qtraj_new = rollout(qtraj, pulse2; n_points = 301)

    @test length(qtraj_new.solution.u) == 301
    @test qtraj_new.pulse === pulse2
end

# ============================================================================ #
# Update system with optimized global parameters
# ============================================================================ #

"""
    Rollouts._update_system!(qtraj::UnitaryTrajectory, sys::QuantumSystem)

Update the system field in a UnitaryTrajectory with a new QuantumSystem
(typically with updated global parameters after optimization).
"""
function Rollouts._update_system!(qtraj::UnitaryTrajectory, sys::QuantumSystem)
    qtraj.system = sys
    return nothing
end

"""
    Rollouts._update_system!(qtraj::KetTrajectory, sys::QuantumSystem)

Update the system field in a KetTrajectory with a new QuantumSystem
(typically with updated global parameters after optimization).
"""
function Rollouts._update_system!(qtraj::KetTrajectory, sys::QuantumSystem)
    qtraj.system = sys
    return nothing
end

"""
    Rollouts._update_system!(qtraj::MultiKetTrajectory, sys::QuantumSystem)

Update the system field in a MultiKetTrajectory with a new QuantumSystem
(typically with updated global parameters after optimization).
"""
function Rollouts._update_system!(qtraj::MultiKetTrajectory, sys::QuantumSystem)
    qtraj.system = sys
    return nothing
end

"""
    Rollouts._update_system!(qtraj::DensityTrajectory, sys::OpenQuantumSystem)

Update the system field in a DensityTrajectory with a new OpenQuantumSystem
(typically with updated global parameters after optimization).
"""
function Rollouts._update_system!(qtraj::DensityTrajectory, sys::OpenQuantumSystem)
    qtraj.system = sys
    return nothing
end

"""
    Rollouts._update_system!(qtraj::SamplingTrajectory, sys::QuantumSystem)

Update the system in the base_trajectory of a SamplingTrajectory.
Note: This only updates the base trajectory's system, not the systems array.
For updating parameter variations in the systems array, that should be done
through the SamplingTrajectory constructor or direct modification.
"""
function Rollouts._update_system!(qtraj::SamplingTrajectory, sys::QuantumSystem)
    _update_system!(qtraj.base_trajectory, sys)
    return nothing
end
