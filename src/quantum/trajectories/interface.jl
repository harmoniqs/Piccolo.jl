# ============================================================================ #
# Common Interface
# ============================================================================ #

"""
    get_system(qtraj)

Get the quantum system from a trajectory.
"""
get_system(qtraj::AbstractQuantumTrajectory) = qtraj.system

"""
    get_pulse(qtraj)

Get the control pulse from a trajectory.
"""
get_pulse(qtraj::AbstractQuantumTrajectory) = qtraj.pulse

"""
    get_initial(qtraj)

Get the initial state/operator from a trajectory.
"""
get_initial(qtraj::UnitaryTrajectory) = qtraj.initial
get_initial(qtraj::KetTrajectory) = qtraj.initial
get_initial(qtraj::MultiKetTrajectory) = qtraj.initials
get_initial(qtraj::DensityTrajectory) = qtraj.initial
get_initial(qtraj::MultiDensityTrajectory) = qtraj.initials

"""
    get_goal(qtraj)

Get the goal state/operator from a trajectory.
"""
get_goal(qtraj::UnitaryTrajectory) = qtraj.goal
get_goal(qtraj::KetTrajectory) = qtraj.goal
get_goal(qtraj::MultiKetTrajectory) = qtraj.goals
get_goal(qtraj::DensityTrajectory) = qtraj.goal
get_goal(qtraj::MultiDensityTrajectory) = qtraj.goals

"""
    get_solution(qtraj)

Get the ODE solution from a trajectory.
"""
get_solution(qtraj::AbstractQuantumTrajectory) = qtraj.solution

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
state_name(::MultiDensityTrajectory) = :ρ⃗̃  # prefix for :ρ⃗̃1, :ρ⃗̃2, etc.

"""
    state_names(qtraj::MultiKetTrajectory)

Get all state names for an ensemble trajectory (`:ψ̃1`, `:ψ̃2`, etc.)
"""
function state_names(qtraj::MultiKetTrajectory)
    prefix = state_name(qtraj)
    return [Symbol(prefix, i) for i = 1:length(qtraj.initials)]
end

"""
    state_names(qtraj::MultiDensityTrajectory)

Get all state names for a multi-density trajectory (`:ρ⃗̃1`, `:ρ⃗̃2`, etc.)
"""
function state_names(qtraj::MultiDensityTrajectory)
    prefix = state_name(qtraj)
    return [Symbol(prefix, i) for i = 1:length(qtraj.initials)]
end

"""
    drive_name(qtraj::AbstractQuantumTrajectory)

Get the drive/control variable name from the trajectory's pulse.
"""
drive_name(qtraj::AbstractQuantumTrajectory) = drive_name(qtraj.pulse)

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
    duration(qtraj)

Get the duration of a trajectory (from its pulse).
"""
duration(qtraj::AbstractQuantumTrajectory) = duration(qtraj.pulse)

# ============================================================================ #
# Rollout - Re-solve ODE with new pulse or different ODE parameters
# ============================================================================ #

"""
    rollout!(qtraj::UnitaryTrajectory, pulse::AbstractPulse; algorithm=MagnusAdapt4(), n_save=101, abstol=1e-8, reltol=1e-8)

Update quantum trajectory in-place with a new pulse by re-solving the ODE.
Mutates `qtraj.pulse` and `qtraj.solution`.

# Arguments
- `qtraj::UnitaryTrajectory`: The trajectory to update
- `pulse::AbstractPulse`: The new control pulse

# Keyword Arguments
- `algorithm`: ODE solver algorithm (default: `MagnusAdapt4()`)
- `n_save::Int`: Number of output time points (default: 101)
- `abstol`: Absolute tolerance (default: 1e-8)
- `reltol`: Relative tolerance (default: 1e-8)

# Example
```julia
qtraj = UnitaryTrajectory(sys, old_pulse, goal)
rollout!(qtraj, new_pulse)  # Updates qtraj in-place
fid = fidelity(qtraj)  # Uses new solution
```

See also: `rollout`
"""
function Rollouts.rollout!(
    qtraj::UnitaryTrajectory,
    pulse::AbstractPulse;
    algorithm = nothing,
    n_save::Int = 101,
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
)
    knot_times = get_knot_times(pulse)
    save_times = collect(range(0.0, duration(pulse), length = n_save))
    tstops = sort(unique(vcat(knot_times, save_times)))
    prob = UnitaryOperatorODEProblem(qtraj.system, pulse, tstops; U0 = qtraj.initial)
    if isnothing(algorithm)
        algorithm = default_algorithm(qtraj.system)
    end
    sol = solve(prob, algorithm; saveat = save_times, abstol = abstol, reltol = reltol)

    qtraj.pulse = pulse
    qtraj.solution = sol
    return nothing
end

"""
    rollout!(qtraj::UnitaryTrajectory; algorithm=MagnusAdapt4(), n_save=101, abstol=1e-8, reltol=1e-8, kwargs...)

Update quantum trajectory in-place by re-solving with same pulse but different ODE parameters.
Mutates `qtraj.solution`.

Useful for comparing different solvers or tolerances.

# Keyword Arguments
- `algorithm`: ODE solver algorithm (default: `MagnusAdapt4()`)
- `n_save::Int`: Number of output time points (default: 101)
- `abstol`: Absolute tolerance (default: 1e-8)
- `reltol`: Relative tolerance (default: 1e-8)
- Additional kwargs passed to `solve`

# Example
```julia
qtraj = UnitaryTrajectory(sys, pulse, goal)

# Compare Magnus vs Runge-Kutta
rollout!(qtraj; algorithm=Tsit5(), abstol=1e-10)
fid_rk = fidelity(qtraj)
```

See also: `rollout`
"""
function Rollouts.rollout!(
    qtraj::UnitaryTrajectory;
    algorithm = nothing,
    n_save::Int = 101,
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
    kwargs...,
)
    knot_times = get_knot_times(qtraj.pulse)
    save_times = collect(range(0.0, duration(qtraj.pulse), length = n_save))
    tstops = sort(unique(vcat(knot_times, save_times)))
    prob = UnitaryOperatorODEProblem(qtraj.system, qtraj.pulse, tstops; U0 = qtraj.initial)
    if isnothing(algorithm)
        algorithm = default_algorithm(qtraj.system)
    end
    sol = solve(
        prob,
        algorithm;
        saveat = save_times,
        abstol = abstol,
        reltol = reltol,
        kwargs...,
    )

    qtraj.solution = sol
    return nothing
end

"""
    rollout!(qtraj::KetTrajectory, pulse::AbstractPulse; algorithm=MagnusAdapt4(), n_save=101, abstol=1e-8, reltol=1e-8)

Update ket trajectory in-place with a new pulse.
See `rollout!(::UnitaryTrajectory, ::AbstractPulse)` for details.
"""
function Rollouts.rollout!(
    qtraj::KetTrajectory,
    pulse::AbstractPulse;
    algorithm = nothing,
    n_save::Int = 101,
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
)
    knot_times = get_knot_times(pulse)
    save_times = collect(range(0.0, duration(pulse), length = n_save))
    tstops = sort(unique(vcat(knot_times, save_times)))
    prob = KetOperatorODEProblem(qtraj.system, pulse, qtraj.initial, tstops)
    if isnothing(algorithm)
        algorithm = default_algorithm(qtraj.system)
    end
    sol = solve(prob, algorithm; saveat = save_times, abstol = abstol, reltol = reltol)

    qtraj.pulse = pulse
    qtraj.solution = sol
    return nothing
end

"""
    rollout!(qtraj::KetTrajectory; algorithm=MagnusAdapt4(), n_save=101, abstol=1e-8, reltol=1e-8, kwargs...)

Update ket trajectory in-place with same pulse but different ODE parameters.
See `rollout!(::UnitaryTrajectory; kwargs...)` for details.
"""
function Rollouts.rollout!(
    qtraj::KetTrajectory;
    algorithm = nothing,
    n_save::Int = 101,
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
    kwargs...,
)
    knot_times = get_knot_times(qtraj.pulse)
    save_times = collect(range(0.0, duration(qtraj.pulse), length = n_save))
    tstops = sort(unique(vcat(knot_times, save_times)))
    prob = KetOperatorODEProblem(qtraj.system, qtraj.pulse, qtraj.initial, tstops)
    if isnothing(algorithm)
        algorithm = default_algorithm(qtraj.system)
    end
    sol = solve(
        prob,
        algorithm;
        saveat = save_times,
        abstol = abstol,
        reltol = reltol,
        kwargs...,
    )

    qtraj.solution = sol
    return nothing
end

"""
    rollout!(qtraj::MultiKetTrajectory, pulse::AbstractPulse; algorithm=MagnusAdapt4(), n_save=101, abstol=1e-8, reltol=1e-8)

Update multi-ket trajectory in-place with a new pulse.
See `rollout!(::UnitaryTrajectory, ::AbstractPulse)` for details.
"""
function Rollouts.rollout!(
    qtraj::MultiKetTrajectory,
    pulse::AbstractPulse;
    algorithm = nothing,
    n_save::Int = 101,
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
)
    knot_times = get_knot_times(pulse)
    save_times = collect(range(0.0, duration(pulse), length = n_save))
    tstops = sort(unique(vcat(knot_times, save_times)))

    # Build ensemble problem
    dummy = zeros(ComplexF64, qtraj.system.levels)
    base_prob = KetOperatorODEProblem(qtraj.system, pulse, dummy, tstops)
    prob_func(prob, i, repeat) = remake(prob, u0 = qtraj.initials[i])
    ensemble_prob = EnsembleProblem(base_prob; prob_func = prob_func)
    if isnothing(algorithm)
        algorithm = default_algorithm(qtraj.system)
    end
    sol = solve(
        ensemble_prob,
        algorithm;
        trajectories = length(qtraj.initials),
        saveat = save_times,
        abstol = abstol,
        reltol = reltol,
    )

    qtraj.pulse = pulse
    qtraj.solution = sol
    return nothing
end

"""
    rollout!(qtraj::MultiKetTrajectory; algorithm=MagnusAdapt4(), n_save=101, abstol=1e-8, reltol=1e-8, kwargs...)

Update multi-ket trajectory in-place with same pulse but different ODE parameters.
See `rollout!(::UnitaryTrajectory; kwargs...)` for details.
"""
function Rollouts.rollout!(
    qtraj::MultiKetTrajectory;
    algorithm = nothing,
    n_save::Int = 101,
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
    kwargs...,
)
    knot_times = get_knot_times(qtraj.pulse)
    save_times = collect(range(0.0, duration(qtraj.pulse), length = n_save))
    tstops = sort(unique(vcat(knot_times, save_times)))

    # Build ensemble problem
    dummy = zeros(ComplexF64, qtraj.system.levels)
    base_prob = KetOperatorODEProblem(qtraj.system, qtraj.pulse, dummy, tstops)
    prob_func(prob, i, repeat) = remake(prob, u0 = qtraj.initials[i])
    ensemble_prob = EnsembleProblem(base_prob; prob_func = prob_func)
    if isnothing(algorithm)
        algorithm = default_algorithm(qtraj.system)
    end
    sol = solve(
        ensemble_prob,
        algorithm;
        trajectories = length(qtraj.initials),
        saveat = save_times,
        abstol = abstol,
        reltol = reltol,
        kwargs...,
    )

    qtraj.solution = sol
    return nothing
end

"""
    rollout!(qtraj::DensityTrajectory, pulse::AbstractPulse; algorithm=Tsit5(), n_save=101, abstol=1e-8, reltol=1e-8)

Update density trajectory in-place with a new pulse.
Note: Default algorithm is `Tsit5()` since density evolution uses standard ODE solvers.
See `rollout!(::UnitaryTrajectory, ::AbstractPulse)` for details.
"""
function Rollouts.rollout!(
    qtraj::DensityTrajectory,
    pulse::AbstractPulse;
    algorithm = Tsit5(),
    n_save::Int = 101,
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
)
    knot_times = get_knot_times(pulse)
    save_times = collect(range(0.0, duration(pulse), length = n_save))
    tstops = sort(unique(vcat(knot_times, save_times)))
    prob = DensityODEProblem(qtraj.system, pulse, qtraj.initial, tstops)
    if isnothing(algorithm)
        algorithm = default_algorithm(qtraj.system)
    end
    sol = solve(prob, algorithm; saveat = save_times, abstol = abstol, reltol = reltol)

    qtraj.pulse = pulse
    qtraj.solution = sol
    return nothing
end

"""
    rollout!(qtraj::DensityTrajectory; algorithm=Tsit5(), n_save=101, abstol=1e-8, reltol=1e-8, kwargs...)

Update density trajectory in-place with same pulse but different ODE parameters.
Note: Default algorithm is `Tsit5()` since density evolution uses standard ODE solvers.
See `rollout!(::UnitaryTrajectory; kwargs...)` for details.
"""
function Rollouts.rollout!(
    qtraj::DensityTrajectory;
    algorithm = Tsit5(),
    n_save::Int = 101,
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
    kwargs...,
)
    knot_times = get_knot_times(qtraj.pulse)
    save_times = collect(range(0.0, duration(qtraj.pulse), length = n_save))
    tstops = sort(unique(vcat(knot_times, save_times)))
    prob = DensityODEProblem(qtraj.system, qtraj.pulse, qtraj.initial, tstops)
    if isnothing(algorithm)
        algorithm = default_algorithm(qtraj.system)
    end
    sol = solve(
        prob,
        algorithm;
        saveat = save_times,
        abstol = abstol,
        reltol = reltol,
        kwargs...,
    )

    qtraj.solution = sol
    return nothing
end

"""
    rollout!(qtraj::MultiDensityTrajectory, pulse::AbstractPulse; algorithm=Tsit5(), n_save=101, abstol=1e-8, reltol=1e-8)

Update multi-density trajectory in-place with a new pulse.
See `rollout!(::UnitaryTrajectory, ::AbstractPulse)` for details.
"""
function Rollouts.rollout!(
    qtraj::MultiDensityTrajectory,
    pulse::AbstractPulse;
    algorithm = Tsit5(),
    n_save::Int = 101,
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
)
    knot_times = get_knot_times(pulse)
    save_times = collect(range(0.0, duration(pulse), length = n_save))
    tstops = sort(unique(vcat(knot_times, save_times)))

    dummy = zeros(ComplexF64, qtraj.system.levels, qtraj.system.levels)
    base_prob = DensityODEProblem(qtraj.system, pulse, dummy, tstops)
    prob_func(prob, i, repeat) = remake(prob, u0 = qtraj.initials[i])
    ensemble_prob = EnsembleProblem(base_prob; prob_func = prob_func)
    if isnothing(algorithm)
        algorithm = default_algorithm(qtraj.system)
    end
    sol = solve(
        ensemble_prob,
        algorithm;
        trajectories = length(qtraj.initials),
        saveat = save_times,
        abstol = abstol,
        reltol = reltol,
    )

    qtraj.pulse = pulse
    qtraj.solution = sol
    return nothing
end

"""
    rollout!(qtraj::MultiDensityTrajectory; algorithm=Tsit5(), n_save=101, abstol=1e-8, reltol=1e-8, kwargs...)

Update multi-density trajectory in-place with same pulse but different ODE parameters.
See `rollout!(::UnitaryTrajectory; kwargs...)` for details.
"""
function Rollouts.rollout!(
    qtraj::MultiDensityTrajectory;
    algorithm = Tsit5(),
    n_save::Int = 101,
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
    kwargs...,
)
    knot_times = get_knot_times(qtraj.pulse)
    save_times = collect(range(0.0, duration(qtraj.pulse), length = n_save))
    tstops = sort(unique(vcat(knot_times, save_times)))

    dummy = zeros(ComplexF64, qtraj.system.levels, qtraj.system.levels)
    base_prob = DensityODEProblem(qtraj.system, qtraj.pulse, dummy, tstops)
    prob_func(prob, i, repeat) = remake(prob, u0 = qtraj.initials[i])
    ensemble_prob = EnsembleProblem(base_prob; prob_func = prob_func)
    if isnothing(algorithm)
        algorithm = default_algorithm(qtraj.system)
    end
    sol = solve(
        ensemble_prob,
        algorithm;
        trajectories = length(qtraj.initials),
        saveat = save_times,
        abstol = abstol,
        reltol = reltol,
        kwargs...,
    )

    qtraj.solution = sol
    return nothing
end

"""
    rollout!(qtraj::SamplingTrajectory, pulse::AbstractPulse; algorithm=MagnusAdapt4(), n_save=101, abstol=1e-8, reltol=1e-8)

Update sampling trajectory's base trajectory in-place with a new pulse.
Delegates to the base trajectory's rollout! method.
"""
function Rollouts.rollout!(
    qtraj::SamplingTrajectory,
    pulse::AbstractPulse;
    algorithm = nothing,
    n_save::Int = 101,
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
)
    rollout!(
        qtraj.base_trajectory,
        pulse;
        algorithm = algorithm,
        n_save = n_save,
        abstol = abstol,
        reltol = reltol,
    )
    return nothing
end

"""
    rollout!(qtraj::SamplingTrajectory; algorithm=MagnusAdapt4(), n_save=101, abstol=1e-8, reltol=1e-8, kwargs...)

Update sampling trajectory's base trajectory in-place with new ODE parameters.
Delegates to the base trajectory's rollout! method.
"""
function Rollouts.rollout!(
    qtraj::SamplingTrajectory;
    algorithm = nothing,
    n_save::Int = 101,
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
    kwargs...,
)
    rollout!(
        qtraj.base_trajectory;
        algorithm = algorithm,
        n_save = n_save,
        abstol = abstol,
        reltol = reltol,
        kwargs...,
    )
    return nothing
end

"""
    rollout(qtraj::UnitaryTrajectory, pulse::AbstractPulse; algorithm=MagnusAdapt4(), n_save=101, abstol=1e-8, reltol=1e-8)

Create a new quantum trajectory by rolling out a new pulse through the system.
Returns a new UnitaryTrajectory with the updated pulse and solution.

# Arguments
- `qtraj::UnitaryTrajectory`: The base trajectory (provides system, initial, goal)
- `pulse::AbstractPulse`: The new control pulse to roll out

# Keyword Arguments
- `algorithm`: ODE solver algorithm (default: `MagnusAdapt4()`)
- `n_save::Int`: Number of output time points (default: 101)
- `abstol`: Absolute tolerance (default: 1e-8)
- `reltol`: Relative tolerance (default: 1e-8)

# Example
```julia
qtraj = UnitaryTrajectory(sys, old_pulse, goal)

# Roll out a new pulse
qtraj_new = rollout(qtraj, new_pulse)

# Check fidelity
fid = fidelity(qtraj_new)
```

See also: `extract_pulse`, `rollout!`, `fidelity`
"""
function Rollouts.rollout(
    qtraj::UnitaryTrajectory,
    pulse::AbstractPulse;
    algorithm = nothing,
    n_save::Int = 101,
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
)
    knot_times = get_knot_times(pulse)
    save_times = collect(range(0.0, duration(pulse), length = n_save))
    tstops = sort(unique(vcat(knot_times, save_times)))
    prob = UnitaryOperatorODEProblem(qtraj.system, pulse, tstops; U0 = qtraj.initial)
    if isnothing(algorithm)
        algorithm = default_algorithm(qtraj.system)
    end
    sol = solve(prob, algorithm; saveat = save_times, abstol = abstol, reltol = reltol)
    return UnitaryTrajectory(qtraj.system, pulse, qtraj.initial, qtraj.goal, sol)
end

"""
    rollout(qtraj::KetTrajectory, pulse::AbstractPulse; algorithm=MagnusAdapt4(), n_save=101, abstol=1e-8, reltol=1e-8)

Create a new ket trajectory by rolling out a new pulse.
See `rollout(::UnitaryTrajectory, ::AbstractPulse)` for details.
"""
function Rollouts.rollout(
    qtraj::KetTrajectory,
    pulse::AbstractPulse;
    algorithm = nothing,
    n_save::Int = 101,
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
)
    knot_times = get_knot_times(pulse)
    save_times = collect(range(0.0, duration(pulse), length = n_save))
    tstops = sort(unique(vcat(knot_times, save_times)))
    prob = KetOperatorODEProblem(qtraj.system, pulse, qtraj.initial, tstops)
    if isnothing(algorithm)
        algorithm = default_algorithm(qtraj.system)
    end
    sol = solve(prob, algorithm; saveat = save_times, abstol = abstol, reltol = reltol)
    return KetTrajectory(qtraj.system, pulse, qtraj.initial, qtraj.goal, sol)
end

"""
    rollout(qtraj::MultiKetTrajectory, pulse::AbstractPulse; algorithm=MagnusAdapt4(), n_save=101, abstol=1e-8, reltol=1e-8)

Create a new multi-ket trajectory by rolling out a new pulse.
See `rollout(::UnitaryTrajectory, ::AbstractPulse)` for details.
"""
function Rollouts.rollout(
    qtraj::MultiKetTrajectory,
    pulse::AbstractPulse;
    algorithm = nothing,
    n_save::Int = 101,
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
)
    knot_times = get_knot_times(pulse)
    save_times = collect(range(0.0, duration(pulse), length = n_save))
    tstops = sort(unique(vcat(knot_times, save_times)))

    # Build ensemble problem
    dummy = zeros(ComplexF64, qtraj.system.levels)
    base_prob = KetOperatorODEProblem(qtraj.system, pulse, dummy, tstops)
    prob_func(prob, i, repeat) = remake(prob, u0 = qtraj.initials[i])
    ensemble_prob = EnsembleProblem(base_prob; prob_func = prob_func)
    if isnothing(algorithm)
        algorithm = default_algorithm(qtraj.system)
    end
    sol = solve(
        ensemble_prob,
        algorithm;
        trajectories = length(qtraj.initials),
        saveat = save_times,
        abstol = abstol,
        reltol = reltol,
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
    rollout(qtraj::DensityTrajectory, pulse::AbstractPulse; algorithm=Tsit5(), n_save=101, abstol=1e-8, reltol=1e-8)

Create a new density trajectory by rolling out a new pulse.
Note: Default algorithm is `Tsit5()` since density evolution uses standard ODE solvers.
See `rollout(::UnitaryTrajectory, ::AbstractPulse)` for details.
"""
function Rollouts.rollout(
    qtraj::DensityTrajectory,
    pulse::AbstractPulse;
    algorithm = Tsit5(),
    n_save::Int = 101,
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
)
    knot_times = get_knot_times(pulse)
    save_times = collect(range(0.0, duration(pulse), length = n_save))
    tstops = sort(unique(vcat(knot_times, save_times)))
    prob = DensityODEProblem(qtraj.system, pulse, qtraj.initial, tstops)
    if isnothing(algorithm)
        algorithm = default_algorithm(qtraj.system)
    end
    sol = solve(prob, algorithm; saveat = save_times, abstol = abstol, reltol = reltol)
    return DensityTrajectory(qtraj.system, pulse, qtraj.initial, qtraj.goal, sol)
end

"""
    rollout(qtraj::MultiDensityTrajectory, pulse::AbstractPulse; algorithm=Tsit5(), n_save=101, abstol=1e-8, reltol=1e-8)

Create a new multi-density trajectory by rolling out a new pulse.
See `rollout(::UnitaryTrajectory, ::AbstractPulse)` for details.
"""
function Rollouts.rollout(
    qtraj::MultiDensityTrajectory,
    pulse::AbstractPulse;
    algorithm = Tsit5(),
    n_save::Int = 101,
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
)
    knot_times = get_knot_times(pulse)
    save_times = collect(range(0.0, duration(pulse), length = n_save))
    tstops = sort(unique(vcat(knot_times, save_times)))

    dummy = zeros(ComplexF64, qtraj.system.levels, qtraj.system.levels)
    base_prob = DensityODEProblem(qtraj.system, pulse, dummy, tstops)
    prob_func(prob, i, repeat) = remake(prob, u0 = qtraj.initials[i])
    ensemble_prob = EnsembleProblem(base_prob; prob_func = prob_func)
    if isnothing(algorithm)
        algorithm = default_algorithm(qtraj.system)
    end
    sol = solve(
        ensemble_prob,
        algorithm;
        trajectories = length(qtraj.initials),
        saveat = save_times,
        abstol = abstol,
        reltol = reltol,
    )

    return MultiDensityTrajectory(
        qtraj.system,
        pulse,
        qtraj.initials,
        qtraj.goals,
        qtraj.weights,
        sol,
    )
end

# Rollout with same pulse, different ODE parameters (non-mutating)

"""
    rollout(qtraj::UnitaryTrajectory; algorithm=MagnusAdapt4(), n_save=101, abstol=1e-8, reltol=1e-8, kwargs...)

Re-solve the trajectory with the same pulse but different ODE parameters.
Returns a new UnitaryTrajectory with the updated solution.

Useful for comparing different solvers or tolerances.

# Keyword Arguments
- `algorithm`: ODE solver algorithm (default: `MagnusAdapt4()`)
- `n_save::Int`: Number of output time points (default: 101)
- `abstol`: Absolute tolerance (default: 1e-8)
- `reltol`: Relative tolerance (default: 1e-8)
- Additional kwargs passed to `solve`

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
    algorithm = nothing,
    n_save::Int = 101,
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
    kwargs...,
)
    knot_times = get_knot_times(qtraj.pulse)
    save_times = collect(range(0.0, duration(qtraj.pulse), length = n_save))
    tstops = sort(unique(vcat(knot_times, save_times)))
    prob = UnitaryOperatorODEProblem(qtraj.system, qtraj.pulse, tstops; U0 = qtraj.initial)
    if isnothing(algorithm)
        algorithm = default_algorithm(qtraj.system)
    end
    sol = solve(
        prob,
        algorithm;
        saveat = save_times,
        abstol = abstol,
        reltol = reltol,
        kwargs...,
    )
    return UnitaryTrajectory(qtraj.system, qtraj.pulse, qtraj.initial, qtraj.goal, sol)
end

"""
    rollout(qtraj::KetTrajectory; algorithm=MagnusAdapt4(), n_save=101, abstol=1e-8, reltol=1e-8, kwargs...)

Re-solve ket trajectory with same pulse but different ODE parameters.
See `rollout(::UnitaryTrajectory; kwargs...)` for details.
"""
function Rollouts.rollout(
    qtraj::KetTrajectory;
    algorithm = nothing,
    n_save::Int = 101,
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
    kwargs...,
)
    knot_times = get_knot_times(qtraj.pulse)
    save_times = collect(range(0.0, duration(qtraj.pulse), length = n_save))
    tstops = sort(unique(vcat(knot_times, save_times)))
    prob = KetOperatorODEProblem(qtraj.system, qtraj.pulse, qtraj.initial, tstops)
    if isnothing(algorithm)
        algorithm = default_algorithm(qtraj.system)
    end
    sol = solve(
        prob,
        algorithm;
        saveat = save_times,
        abstol = abstol,
        reltol = reltol,
        kwargs...,
    )
    return KetTrajectory(qtraj.system, qtraj.pulse, qtraj.initial, qtraj.goal, sol)
end

"""
    rollout(qtraj::MultiKetTrajectory; algorithm=MagnusAdapt4(), n_save=101, abstol=1e-8, reltol=1e-8, kwargs...)

Re-solve multi-ket trajectory with same pulse but different ODE parameters.
See `rollout(::UnitaryTrajectory; kwargs...)` for details.
"""
function Rollouts.rollout(
    qtraj::MultiKetTrajectory;
    algorithm = nothing,
    n_save::Int = 101,
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
    kwargs...,
)
    knot_times = get_knot_times(qtraj.pulse)
    save_times = collect(range(0.0, duration(qtraj.pulse), length = n_save))
    tstops = sort(unique(vcat(knot_times, save_times)))

    # Build ensemble problem
    dummy = zeros(ComplexF64, qtraj.system.levels)
    base_prob = KetOperatorODEProblem(qtraj.system, qtraj.pulse, dummy, tstops)
    prob_func(prob, i, repeat) = remake(prob, u0 = qtraj.initials[i])
    ensemble_prob = EnsembleProblem(base_prob; prob_func = prob_func)
    if isnothing(algorithm)
        algorithm = default_algorithm(qtraj.system)
    end
    sol = solve(
        ensemble_prob,
        algorithm;
        trajectories = length(qtraj.initials),
        saveat = save_times,
        abstol = abstol,
        reltol = reltol,
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
    rollout(qtraj::DensityTrajectory; algorithm=Tsit5(), n_save=101, abstol=1e-8, reltol=1e-8, kwargs...)

Re-solve density trajectory with same pulse but different ODE parameters.
Note: Default algorithm is `Tsit5()` since density evolution uses standard ODE solvers.
See `rollout(::UnitaryTrajectory; kwargs...)` for details.
"""
function Rollouts.rollout(
    qtraj::DensityTrajectory;
    algorithm = Tsit5(),
    n_save::Int = 101,
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
    kwargs...,
)
    knot_times = get_knot_times(qtraj.pulse)
    save_times = collect(range(0.0, duration(qtraj.pulse), length = n_save))
    tstops = sort(unique(vcat(knot_times, save_times)))
    prob = DensityODEProblem(qtraj.system, qtraj.pulse, qtraj.initial, tstops)
    if isnothing(algorithm)
        algorithm = default_algorithm(qtraj.system)
    end
    sol = solve(
        prob,
        algorithm;
        saveat = save_times,
        abstol = abstol,
        reltol = reltol,
        kwargs...,
    )
    return DensityTrajectory(qtraj.system, qtraj.pulse, qtraj.initial, qtraj.goal, sol)
end

"""
    rollout(qtraj::MultiDensityTrajectory; algorithm=Tsit5(), n_save=101, abstol=1e-8, reltol=1e-8, kwargs...)

Re-solve multi-density trajectory with same pulse but different ODE parameters.
See `rollout(::UnitaryTrajectory; kwargs...)` for details.
"""
function Rollouts.rollout(
    qtraj::MultiDensityTrajectory;
    algorithm = Tsit5(),
    n_save::Int = 101,
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
    kwargs...,
)
    knot_times = get_knot_times(qtraj.pulse)
    save_times = collect(range(0.0, duration(qtraj.pulse), length = n_save))
    tstops = sort(unique(vcat(knot_times, save_times)))

    dummy = zeros(ComplexF64, qtraj.system.levels, qtraj.system.levels)
    base_prob = DensityODEProblem(qtraj.system, qtraj.pulse, dummy, tstops)
    prob_func(prob, i, repeat) = remake(prob, u0 = qtraj.initials[i])
    ensemble_prob = EnsembleProblem(base_prob; prob_func = prob_func)
    if isnothing(algorithm)
        algorithm = default_algorithm(qtraj.system)
    end
    sol = solve(
        ensemble_prob,
        algorithm;
        trajectories = length(qtraj.initials),
        saveat = save_times,
        abstol = abstol,
        reltol = reltol,
        kwargs...,
    )

    return MultiDensityTrajectory(
        qtraj.system,
        qtraj.pulse,
        qtraj.initials,
        qtraj.goals,
        qtraj.weights,
        sol,
    )
end

# ============================================================================ #
# Fidelity (extending Rollouts.fidelity)
# ============================================================================ #

"""
    fidelity(qtraj::UnitaryTrajectory; subspace=nothing, phases=nothing)

Compute the fidelity between the final unitary and the goal.

For `EmbeddedOperator` goals, uses the Pedersen average gate fidelity on the
computational subspace, matching the metric the optimizer minimizes:

    F = 1/(n(n+1)) * (|tr(M'M)| + |tr(M)|²),   M = U_goal' * U_sub

For standard goals (or when an explicit `subspace` is given), uses the standard
unitary fidelity |tr(U'U_goal)|²/N².

When `phases` is provided (a vector of per-qubit Z-phases of length `n_qubits`),
the goal on the computational subspace is adjusted by `Diagonal(phase_vec) *
U_goal_sub` before computing fidelity, where `phase_vec` is a length-`n_sub`
vector whose entries are products of `exp(im * θⱼ)` for each qubit `j` that is
in the `|1⟩` state of the corresponding computational basis vector (via binary
decomposition of the basis index). This matches the free-phase objective used
during optimization.
"""
function Rollouts.fidelity(
    qtraj::UnitaryTrajectory;
    subspace::Union{Nothing,AbstractVector{Int}} = nothing,
    phases::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    U_final = qtraj.solution.u[end]
    if qtraj.goal isa EmbeddedOperator && isnothing(subspace)
        U_goal_sub = if isnothing(phases)
            unembed(qtraj.goal)
        else
            # Apply free phases: same convention as _make_free_phase_goal
            U_base = unembed(qtraj.goal)
            n_sub = size(U_base, 1)
            n_qubits = length(phases)
            phase_diag = map(1:n_sub) do i
                bits = i - 1
                phase = sum(
                    phases[j] for j = 1:n_qubits if (bits >> (n_qubits - j)) & 1 == 1;
                    init = 0.0,
                )
                return exp(im * phase)
            end
            Diagonal(phase_diag) * U_base
        end
        # Use Pedersen formula, consistent with unitary_fidelity_loss(Ũ⃗, ::EmbeddedOperator)
        U_sub = U_final[qtraj.goal.subspace, qtraj.goal.subspace]
        n = length(qtraj.goal.subspace)
        M = U_goal_sub' * U_sub
        return 1 / (n * (n + 1)) * (abs(tr(M' * M)) + abs2(tr(M)))
    else
        if !isnothing(phases)
            @warn "`phases` kwarg is ignored when goal is not an EmbeddedOperator or when `subspace` is provided"
        end
        U_goal = qtraj.goal isa EmbeddedOperator ? qtraj.goal.operator : qtraj.goal
        if isnothing(subspace)
            return unitary_fidelity(U_final, U_goal)
        else
            return unitary_fidelity(U_final, U_goal; subspace = subspace)
        end
    end
end

"""
    fidelity(qtraj::KetTrajectory)

Compute the fidelity between the final state and the goal.
"""
function Rollouts.fidelity(qtraj::KetTrajectory)
    ψ_final = qtraj.solution.u[end]
    return abs2(ψ_final' * qtraj.goal)
end

"""
    fidelity(qtraj::MultiKetTrajectory; phases=nothing, subsystem_levels=nothing)

Compute the coherent fidelity across all state transfers:

    F = |1/n ∑ᵢ ⟨ψᵢ_goal(θ)|ψᵢ⟩|²

When `phases` is provided (a vector of per-qubit Z-phases), each goal is
phase-rotated before computing the overlap. This matches the
`CoherentKetFreePhaseInfidelityObjective` used during optimization.

`subsystem_levels` specifies the Hilbert space dimensions of each subsystem
(e.g., [2, 2, 2, 2] for 4 qubits). Required when `phases` is provided.
"""
function Rollouts.fidelity(
    qtraj::MultiKetTrajectory;
    phases::Union{Nothing,AbstractVector{<:Real}} = nothing,
    subsystem_levels::Union{Nothing,Vector{Int}} = nothing,
)
    n = length(qtraj.goals)
    goals = if !isnothing(phases)
        @assert !isnothing(subsystem_levels) "subsystem_levels required when phases is provided"
        n_qubits = length(phases)
        dim = prod(subsystem_levels)
        # Build phase diagonal: for each basis state, accumulate phases for qubits in |1⟩
        phase_diag = map(1:dim) do i
            bits = i - 1
            phase = sum(
                phases[j] for j in 1:n_qubits
                if (bits >> (n_qubits - j)) & 1 == 1;
                init = 0.0
            )
            return exp(im * phase)
        end
        [phase_diag .* g for g in qtraj.goals]
    else
        qtraj.goals
    end
    overlap_sum = sum(goals[i]' * qtraj.solution[i].u[end] for i = 1:n)
    return abs2(overlap_sum / n)
end

"""
    fidelity(qtraj::DensityTrajectory)

Compute the fidelity between the final density matrix and the goal.
Uses trace fidelity: F = tr(ρ_final * ρ_goal)
"""
function Rollouts.fidelity(qtraj::DensityTrajectory)
    ρ_final = qtraj.solution.u[end]
    return real(tr(ρ_final * qtraj.goal))
end

"""
    fidelity(qtraj::MultiDensityTrajectory)

Compute the weighted average trace fidelity across all density matrix transfers:

    F = ∑ᵢ wᵢ tr(ρᵢ_final * ρᵢ_goal)

Unlike MultiKetTrajectory (which uses a coherent sum), density matrices have no
global phase ambiguity, so we use an incoherent weighted average.
"""
function Rollouts.fidelity(qtraj::MultiDensityTrajectory)
    n = length(qtraj.goals)
    return sum(
        qtraj.weights[i] * real(tr(qtraj.solution[i].u[end] * qtraj.goals[i])) for i = 1:n
    )
end

"""
    fidelity(qtraj::SamplingTrajectory; kwargs...)

Compute the fidelity for each system in the sampling trajectory.

Returns a vector of fidelities, one per system, by rolling out the current pulse
with each system and computing the fidelity against the goal.
"""
function Rollouts.fidelity(qtraj::SamplingTrajectory; kwargs...)
    base = qtraj.base_trajectory
    return [Rollouts.fidelity(_swap_system(base, sys); kwargs...) for sys in qtraj.systems]
end

# Helpers to create a trajectory with a different system (for per-system fidelity evaluation)
_swap_system(qtraj::UnitaryTrajectory, sys::AbstractQuantumSystem) =
    UnitaryTrajectory(sys, qtraj.pulse, get_goal(qtraj))

_swap_system(qtraj::KetTrajectory, sys::AbstractQuantumSystem) =
    KetTrajectory(sys, qtraj.pulse, qtraj.initial, qtraj.goal)

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
    qtraj_fine = rollout(qtraj, pulse2; n_save = 501)
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
    qtraj_new = rollout(qtraj, pulse2; n_save = 201)

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
    qtraj_new = rollout(qtraj, pulse2; n_save = 301)

    @test length(qtraj_new.solution.u) == 301
    @test qtraj_new.pulse === pulse2
end

@testitem "fidelity with EmbeddedOperator uses Pedersen formula" begin
    using LinearAlgebra

    # 2-level system embedded in 3 levels
    H_drift_3 = ComplexF64[0 0 0; 0 1 0; 0 0 2]
    H_drive_3 = ComplexF64[0 1 0; 1 0 1; 0 1 0] / √2
    sys = QuantumSystem(H_drift_3, [H_drive_3], [1.0])

    σx = ComplexF64[0 1; 1 0]
    subspace = [1, 2]
    levels = [3]
    U_goal = EmbeddedOperator(σx, subspace, levels)

    # Create trajectory — fidelity depends on the rollout result
    pulse = ZeroOrderPulse(0.1 * randn(1, 20), collect(range(0.0, 5.0, length = 20)))
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)

    # Should not error — uses Pedersen formula for EmbeddedOperator
    fid = fidelity(qtraj)
    @test 0.0 <= fid <= 1.0
end

@testitem "fidelity with phases kwarg" begin
    using LinearAlgebra

    H_drift_3 = ComplexF64[0 0 0; 0 1 0; 0 0 2]
    H_drive_3 = ComplexF64[0 1 0; 1 0 1; 0 1 0] / √2
    sys = QuantumSystem(H_drift_3, [H_drive_3], [1.0])

    σx = ComplexF64[0 1; 1 0]
    subspace = [1, 2]
    levels = [3]
    U_goal = EmbeddedOperator(σx, subspace, levels)

    pulse = ZeroOrderPulse(0.1 * randn(1, 20), collect(range(0.0, 5.0, length = 20)))
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)

    # phases=zeros should equal no-phase fidelity
    fid_nophase = fidelity(qtraj)
    fid_zerophase = fidelity(qtraj; phases = [0.0])
    @test fid_nophase ≈ fid_zerophase atol = 1e-12

    # phases=nonzero should give a different fidelity
    fid_phase = fidelity(qtraj; phases = [π / 4])
    # Not necessarily different (depends on state) but should not error
    @test 0.0 <= fid_phase <= 1.0
end

@testitem "fidelity with plain matrix goal and subspace" begin
    using LinearAlgebra

    sys = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])
    X_gate = ComplexF64[0 1; 1 0]

    pulse = ZeroOrderPulse([0.5 0.5], [0.0, 1.0])
    qtraj = UnitaryTrajectory(sys, pulse, X_gate)

    # Standard fidelity (no EmbeddedOperator)
    fid = fidelity(qtraj)
    @test 0.0 <= fid <= 1.0

    # With explicit subspace
    fid_sub = fidelity(qtraj; subspace = [1, 2])
    @test 0.0 <= fid_sub <= 1.0
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
    Rollouts._update_system!(qtraj::MultiDensityTrajectory, sys::OpenQuantumSystem)

Update the system field in a MultiDensityTrajectory with a new OpenQuantumSystem
(typically with updated global parameters after optimization).
"""
function Rollouts._update_system!(qtraj::MultiDensityTrajectory, sys::OpenQuantumSystem)
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
