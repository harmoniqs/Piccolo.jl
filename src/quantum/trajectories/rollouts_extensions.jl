# ============================================================================ #
# Rollouts Extensions
#
# This file contains the user-facing entry points for testing a quantum
# trajectory against a system (noisy, shifted, or alternate). It extends the
# `Rollouts` module with methods for each `AbstractQuantumTrajectory` subtype:
#
#   - `rollout!`  — update a trajectory in place by re-solving the ODE
#   - `rollout`   — construct a new trajectory with an updated pulse or tol
#   - `fidelity`  — evaluate fidelity from the stored ODE solution
#   - `_update_system!` — swap in a system with optimized global parameters
#
# This is kept separate from the common interface so users searching for
# "rollout" or "fidelity" find them in an aptly named file.
# ============================================================================ #

# ============================================================================ #
# rollout!  (in-place, pulse + tol variants)
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
    prob_func(prob, i_or_ctx, _repeat = nothing) =
        remake(prob, u0 = qtraj.initials[_sim_index(i_or_ctx)])
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
    prob_func(prob, i_or_ctx, _repeat = nothing) =
        remake(prob, u0 = qtraj.initials[_sim_index(i_or_ctx)])
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
    prob_func(prob, i_or_ctx, _repeat = nothing) =
        remake(prob, u0 = qtraj.initials[_sim_index(i_or_ctx)])
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
    prob_func(prob, i_or_ctx, _repeat = nothing) =
        remake(prob, u0 = qtraj.initials[_sim_index(i_or_ctx)])
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

# ============================================================================ #
# rollout  (non-mutating, pulse + tol variants)
# ============================================================================ #

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
    prob_func(prob, i_or_ctx, _repeat = nothing) =
        remake(prob, u0 = qtraj.initials[_sim_index(i_or_ctx)])
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
    prob_func(prob, i_or_ctx, _repeat = nothing) =
        remake(prob, u0 = qtraj.initials[_sim_index(i_or_ctx)])
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
    prob_func(prob, i_or_ctx, _repeat = nothing) =
        remake(prob, u0 = qtraj.initials[_sim_index(i_or_ctx)])
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
    prob_func(prob, i_or_ctx, _repeat = nothing) =
        remake(prob, u0 = qtraj.initials[_sim_index(i_or_ctx)])
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

export number_operator_phase_diag

@doc raw"""
    number_operator_phase_diag(θ, subsystem_levels) -> Vector

Diagonal of the free-phase rotation ``\bigotimes_j e^{i θ_j \hat n_j}`` in the tensor-product
basis: for basis index `idx` (0-based), level ``s_j`` of subsystem `j` contributes ``s_j θ_j``.

This is the **number-operator** convention, the one `free_phase = true` applies to *ket* goals.
It is deliberately not the same expression as the qubit/binary decomposition that `fidelity`
applies to an `EmbeddedOperator` goal — the two agree when every subsystem is a qubit
(``s_j ∈ \{0,1\}``) and differ above two levels. Kept type-generic so ForwardDiff `Dual`s pass
through.
"""
function number_operator_phase_diag(θ, subsystem_levels::AbstractVector{Int})
    dim = prod(subsystem_levels)
    n_sub = length(subsystem_levels)
    return map(0:(dim-1)) do idx
        phase = zero(eltype(θ))
        remaining = idx
        for j = 1:n_sub
            stride = prod(subsystem_levels[k] for k = (j+1):n_sub; init = 1)
            sj = remaining ÷ stride
            remaining = remaining % stride
            phase += sj * θ[j]
        end
        return exp(im * phase)
    end
end

"""
    fidelity(qtraj::KetTrajectory; phases=nothing, subsystem_levels=nothing)

Compute the fidelity between the final state and the goal.

`phases` applies the [`number_operator_phase_diag`](@ref) rotation to the goal before comparing —
the convention `free_phase = true` optimizes against. `subsystem_levels` gives the tensor-product
factorization; it may be omitted only for a single phase, since with more than one the
factorization cannot be inferred from the state dimension alone.
"""
function Rollouts.fidelity(
    qtraj::KetTrajectory;
    phases::Union{Nothing,AbstractVector{<:Real}} = nothing,
    subsystem_levels::Union{Nothing,AbstractVector{Int}} = nothing,
)
    ψ_final = qtraj.solution.u[end]
    isnothing(phases) && return abs2(ψ_final' * qtraj.goal)

    levels = if !isnothing(subsystem_levels)
        subsystem_levels
    elseif length(phases) == 1
        [length(qtraj.goal)]
    else
        error(
            "fidelity(::KetTrajectory; phases) with $(length(phases)) phases needs " *
            "`subsystem_levels`: the tensor-product factorization cannot be inferred from " *
            "the state dimension $(length(qtraj.goal)) alone.",
        )
    end

    if prod(levels) != length(qtraj.goal)
        error(
            "subsystem_levels = $levels implies dimension $(prod(levels)), but the goal has " *
            "$(length(qtraj.goal)) components.",
        )
    end

    goal_phased = number_operator_phase_diag(phases, collect(Int, levels)) .* qtraj.goal
    return abs2(ψ_final' * goal_phased)
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
# ── rollout = :none (RolloutStates) guards ─────────────────────────────────────
# A `RolloutStates` trajectory must never be CPU-rolled: `rollout!` REASSIGNS
# `qtraj.solution`, and the solution type parameter S is fixed at construction —
# the assignment would convert-error (worse: only AFTER the ODE ran). Placeholder
# READS are fine (they are the bootstrap warm-start guess); only re-rollout and
# placeholder fidelity are guarded.
const _ROLLOUT_NONE_ERR =
    "this MultiKetTrajectory was built with `rollout = :none` (placeholder states); " *
    "CPU `rollout!` is not supported on it — refresh it with a GPU rollout " *
    "(Piccolissimo `gpu_rollout!`) or reconstruct with `rollout = :cpu`"

Rollouts.rollout!(
    ::MultiKetTrajectory{<:AbstractPulse,RolloutStates},
    ::AbstractPulse;
    kwargs...,
) = error(_ROLLOUT_NONE_ERR)

Rollouts.rollout!(::MultiKetTrajectory{<:AbstractPulse,RolloutStates}; kwargs...) =
    error(_ROLLOUT_NONE_ERR)

function Rollouts.fidelity(
    qtraj::MultiKetTrajectory{<:AbstractPulse,RolloutStates};
    kwargs...,
)
    if !qtraj.solution.real_states
        @warn "fidelity() on a rollout=:none placeholder trajectory — the states are " *
              "the tiled initial-ket guess, NOT solved dynamics. Refresh with " *
              "gpu_rollout! first." maxlog = 1
    end
    return invoke(Rollouts.fidelity, Tuple{MultiKetTrajectory}, qtraj; kwargs...)
end

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
                phases[j] for j = 1:n_qubits if (bits >> (n_qubits - j)) & 1 == 1;
                init = 0.0,
            )
            return exp(im * phase)
        end
        [phase_diag .* g for g in qtraj.goals]
    else
        qtraj.goals
    end
    overlap_sum = sum(goals[i]' * qtraj.solution.u[i].u[end] for i = 1:n)
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
        qtraj.weights[i] * real(tr(qtraj.solution.u[i].u[end] * qtraj.goals[i])) for i = 1:n
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

_swap_system(qtraj::DensityTrajectory, sys::OpenQuantumSystem) =
    DensityTrajectory(sys, qtraj.pulse, qtraj.initial, qtraj.goal)

_swap_system(qtraj::MultiKetTrajectory, sys::AbstractQuantumSystem) = MultiKetTrajectory(
    sys,
    qtraj.pulse,
    qtraj.initials,
    qtraj.goals;
    weights = qtraj.weights,
)

_swap_system(qtraj::MultiDensityTrajectory, sys::OpenQuantumSystem) =
    MultiDensityTrajectory(
        sys,
        qtraj.pulse,
        qtraj.initials,
        qtraj.goals;
        weights = qtraj.weights,
    )

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
    Rollouts._update_system!(qtraj.base_trajectory, sys)
    return nothing
end

# ============================================================================ #
# Tests
# ============================================================================ #

@testitem "number_operator_phase_diag pins the ket free-phase convention" begin
    using LinearAlgebra

    # This is the NUMBER-OPERATOR convention: for basis index idx, level s_j of subsystem j
    # contributes s_j * θ_j. It is deliberately NOT the qubit/binary decomposition that
    # `fidelity(::UnitaryTrajectory; phases)` applies to an EmbeddedOperator goal. Conflating the
    # two is one of the incompatible free-phase conventions this cleanup exists to separate, so
    # both the agreement and the divergence are pinned here.
    θ = 0.3

    # One 3-level subsystem: levels 0,1,2 acquire 0, θ, 2θ.
    d3 = number_operator_phase_diag([θ], [3])
    @test length(d3) == 3
    @test d3 ≈ [1, exp(im * θ), exp(2im * θ)]

    # The qubit/binary form would give 0, θ, θ — it cannot express the 2θ. This is exactly why
    # the ket path needs its own function.
    @test !isapprox(d3[3], exp(im * θ); atol = 1e-8)

    # Two qubits: (s1,s2) = (0,0),(0,1),(1,0),(1,1) over idx 0..3.
    θ1, θ2 = 0.2, -0.5
    d22 = number_operator_phase_diag([θ1, θ2], [2, 2])
    @test d22 ≈ [1, exp(im * θ2), exp(im * θ1), exp(im * (θ1 + θ2))]

    # On qubits ONLY, the number-operator and binary conventions coincide (s ∈ {0,1}).
    binary = map(0:3) do bits
        φ = 0.0
        ((bits >> 1) & 1 == 1) && (φ += θ1)
        ((bits >> 0) & 1 == 1) && (φ += θ2)
        exp(im * φ)
    end
    @test d22 ≈ binary

    # Zero phases are the identity, and the result is always unit-modulus.
    @test number_operator_phase_diag([0.0, 0.0], [2, 2]) ≈ ones(4)
    @test all(x -> abs(x) ≈ 1, d3)
end

@testitem "fidelity(::KetTrajectory; phases) applies the rotation and guards its inputs" begin
    using LinearAlgebra

    # A superposition goal is required for the rotation to be observable: a single basis-state
    # goal only picks up a global phase, which abs2 discards.
    sys = OpenQuantumSystem(zeros(ComplexF64, 2, 2), [ComplexF64[0 1; 1 0]], [(-2.0, 2.0)])
    N, T = 11, 1.0
    times = collect(range(0, T, length = N))
    pulse = LinearSplinePulse(fill(0.4, 1, N), times)
    goal = ComplexF64[1, 1] / √2
    qtraj = KetTrajectory(sys, pulse, ComplexF64[1, 0], goal)

    F_plain = fidelity(qtraj)

    # phases = 0 must be EXACTLY the unrotated value — the rotation is the identity there.
    @test fidelity(qtraj; phases = [0.0]) ≈ F_plain

    # A nonzero phase must actually change the number, or the kwarg is inert.
    @test !isapprox(fidelity(qtraj; phases = [1.1]), F_plain; atol = 1e-8)

    # subsystem_levels may be omitted for a single phase (whole space = one subsystem) and must
    # agree with passing it explicitly.
    @test fidelity(qtraj; phases = [1.1]) ≈
          fidelity(qtraj; phases = [1.1], subsystem_levels = [2])

    # With MORE than one phase the factorization is not inferable from the state dimension.
    err = try
        fidelity(qtraj; phases = [0.1, 0.2])
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("subsystem_levels", err.msg)

    # A factorization inconsistent with the goal dimension must be refused, not silently reshaped.
    err2 = try
        fidelity(qtraj; phases = [0.1, 0.2], subsystem_levels = [3, 3])
        nothing
    catch e
        e
    end
    @test err2 isa ErrorException
    @test occursin("implies dimension", err2.msg)
end

@testitem "rollout - UnitaryTrajectory" begin
    using LinearAlgebra
    using OrdinaryDiffEqLinear: MagnusGL4
    using OrdinaryDiffEqTsit5: Tsit5

    system = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])
    X_gate = ComplexF64[0 1; 1 0]

    # Create trajectory with initial pulse (1 drive × 2 knot points)
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

    # In-place variants
    rollout!(qtraj; n_save = 41)
    @test length(qtraj.solution.u) == 41
    rollout!(qtraj, pulse2)
    @test qtraj.pulse === pulse2
    @test length(qtraj.solution.u) == 101
end

@testitem "rollout - KetTrajectory" begin
    using LinearAlgebra
    using OrdinaryDiffEqLinear: MagnusGL4

    system = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])
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

    # In-place variants
    rollout!(qtraj, pulse2)
    @test qtraj.pulse === pulse2
    @test length(qtraj.solution.u) == 101
    rollout!(qtraj; n_save = 41)
    @test length(qtraj.solution.u) == 41

    # Same-pulse re-solve returns a NEW trajectory
    qtraj_resolved = rollout(qtraj; n_save = 41)
    @test qtraj_resolved isa KetTrajectory
    @test qtraj_resolved !== qtraj
    @test length(qtraj_resolved.solution.u) == 41
    @test qtraj_resolved.pulse === qtraj.pulse
end

@testitem "rollout preserves MultiKetTrajectory solution structure" begin
    using LinearAlgebra

    sys = OpenQuantumSystem(GATES[:Z], [GATES[:X]], [1.0])
    psi0 = ComplexF64[1.0, 0.0]
    psi1 = ComplexF64[0.0, 1.0]

    pulse1 = ZeroOrderPulse([0.5 0.5], [0.0, 1.0])
    qtraj = MultiKetTrajectory(sys, pulse1, [psi0, psi1], [psi1, psi0])

    # Roll out with a new pulse
    pulse2 = ZeroOrderPulse([0.8 0.8], [0.0, 1.0])
    qtraj_new = rollout(qtraj, pulse2)

    # Solution structure must be preserved
    @test length(qtraj_new.solution.u) == 2
    @test qtraj_new.pulse === pulse2

    # New solution must not be all zeros
    @test !all(x -> x == zero(x), qtraj_new.solution.u[1].u[end])
    @test !all(x -> x == zero(x), qtraj_new.solution.u[2].u[end])

    # Fidelity must not be exactly 0.0
    @test fidelity(qtraj_new) != 0.0
end

@testitem "rollout - MultiKetTrajectory" begin
    using LinearAlgebra

    system = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])
    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]

    initials = [ψ0, ψ1]
    goals = [ψ1, ψ0]

    pulse1 = ZeroOrderPulse([0.5 0.5], [0.0, 1.0])
    qtraj = MultiKetTrajectory(system, pulse1, initials, goals)

    # Roll out new pulse
    pulse2 = ZeroOrderPulse([0.8 0.8], [0.0, 1.0])
    qtraj_new = rollout(qtraj, pulse2)

    @test length(qtraj_new.solution.u) == 2
    @test qtraj_new.pulse === pulse2

    # In-place variants
    rollout!(qtraj, pulse2)
    @test qtraj.pulse === pulse2
    @test length(qtraj.solution.u) == 2
    rollout!(qtraj; n_save = 41)
    @test length(qtraj.solution.u) == 2
    @test length(qtraj.solution.u[1].t) == 41

    # Same-pulse re-solve returns a NEW trajectory
    qtraj_resolved = rollout(qtraj; n_save = 41)
    @test qtraj_resolved isa MultiKetTrajectory
    @test qtraj_resolved !== qtraj
    @test length(qtraj_resolved.solution.u) == 2
    @test qtraj_resolved.weights == qtraj.weights
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

    # In-place variants
    rollout!(qtraj, pulse2)
    @test qtraj.pulse === pulse2
    @test length(qtraj.solution.u) == 101
    rollout!(qtraj; n_save = 41)
    @test length(qtraj.solution.u) == 41
    rollout!(qtraj; algorithm = nothing)  # falls back to default_algorithm(::OpenQuantumSystem)
    @test length(qtraj.solution.u) == 101

    # Same-pulse re-solve returns a NEW trajectory
    qtraj_resolved = rollout(qtraj; n_save = 41)
    @test qtraj_resolved isa DensityTrajectory
    @test qtraj_resolved !== qtraj
    @test length(qtraj_resolved.solution.u) == 41
end

@testitem "fidelity with EmbeddedOperator uses Pedersen formula" begin
    using LinearAlgebra

    # 2-level system embedded in 3 levels
    H_drift_3 = ComplexF64[0 0 0; 0 1 0; 0 0 2]
    H_drive_3 = ComplexF64[0 1 0; 1 0 1; 0 1 0] / √2
    sys = OpenQuantumSystem(H_drift_3, [H_drive_3], [1.0])

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
    sys = OpenQuantumSystem(H_drift_3, [H_drive_3], [1.0])

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

    sys = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])
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

@testitem "rollout - MultiDensityTrajectory" begin
    using LinearAlgebra

    L = ComplexF64[0.1 0.0; 0.0 0.0]
    system = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])

    ρ0s = [ComplexF64[1.0 0.0; 0.0 0.0], ComplexF64[0.0 0.0; 0.0 1.0]]
    ρgs = [ComplexF64[0.0 0.0; 0.0 1.0], ComplexF64[1.0 0.0; 0.0 0.0]]

    pulse1 = ZeroOrderPulse([0.5 0.5], [0.0, 1.0])
    qtraj = MultiDensityTrajectory(system, pulse1, ρ0s, ρgs)

    # Roll out a new pulse (non-mutating)
    pulse2 = ZeroOrderPulse([0.8 0.8], [0.0, 1.0])
    qtraj_new = rollout(qtraj, pulse2; n_save = 201)

    @test qtraj_new isa MultiDensityTrajectory
    @test qtraj_new !== qtraj
    @test length(qtraj_new.solution.u) == 2
    @test length(qtraj_new.solution.u[1].t) == 201
    @test qtraj_new.pulse === pulse2
    @test qtraj_new.system === qtraj.system
    @test qtraj_new.initials == ρ0s
    @test qtraj_new.goals == ρgs
    @test qtraj_new.weights == qtraj.weights
    @test 0.0 <= fidelity(qtraj_new) <= 1.0

    # Same-pulse re-solve with new ODE parameters (non-mutating)
    qtraj_fine = rollout(qtraj; n_save = 201)
    @test qtraj_fine isa MultiDensityTrajectory
    @test qtraj_fine.pulse === pulse1
    @test length(qtraj_fine.solution.u) == 2

    # In-place variants
    rollout!(qtraj, pulse2)
    @test qtraj.pulse === pulse2
    @test length(qtraj.solution.u) == 2
    rollout!(qtraj; n_save = 41)
    @test length(qtraj.solution.u) == 2
    @test length(qtraj.solution.u[1].t) == 41
    rollout!(qtraj; algorithm = nothing)  # falls back to default_algorithm(::OpenQuantumSystem)
    @test length(qtraj.solution.u) == 2
end

@testitem "fidelity warns and ignores phases for non-EmbeddedOperator goals" begin
    using LinearAlgebra

    # Plain matrix goal: `phases` cannot be applied — the kwarg must warn, not silently apply.
    sys = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])
    X_gate = ComplexF64[0 1; 1 0]
    pulse = ZeroOrderPulse([0.5 0.5], [0.0, 1.0])
    qtraj = UnitaryTrajectory(sys, pulse, X_gate)

    f_unphased = fidelity(qtraj)
    f_phased = @test_logs (:warn, r"phases` kwarg is ignored") match_mode = :any fidelity(
        qtraj;
        phases = [0.3],
    )
    @test f_phased ≈ f_unphased

    # EmbeddedOperator goal with an explicit `subspace`: same warning, standard fidelity path.
    H_drift_3 = ComplexF64[0 0 0; 0 1 0; 0 0 2]
    H_drive_3 = ComplexF64[0 1 0; 1 0 1; 0 1 0] / √2
    sys3 = OpenQuantumSystem(H_drift_3, [H_drive_3], [1.0])
    U_goal = EmbeddedOperator(X_gate, [1, 2], [3])
    qt3 = UnitaryTrajectory(sys3, pulse, U_goal)

    f_sub = @test_logs (:warn, r"phases` kwarg is ignored") match_mode = :any fidelity(
        qt3;
        subspace = [1, 2],
        phases = [0.0],
    )
    @test 0.0 <= f_sub <= 1.0
end

@testitem "fidelity(::MultiKetTrajectory; phases) rotates goals coherently" begin
    using LinearAlgebra

    # Zero drift and zero controls: the states never move, so the phased coherent
    # fidelity is hand-computable. With goals == initials only the |01⟩ goal picks up
    # a phase, so F = |(1 + e^{-iθ₂})/2|² = cos²(θ₂/2) — independent of θ₁.
    sys = OpenQuantumSystem(zeros(ComplexF64, 4, 4), [kron(PAULIS[:X], PAULIS[:I])], [1.0])
    ψa = ComplexF64[1, 0, 0, 0]
    ψb = ComplexF64[0, 1, 0, 0]
    pulse = ZeroOrderPulse(zeros(1, 2), [0.0, 1.0])
    qtraj = MultiKetTrajectory(sys, pulse, [ψa, ψb], [ψa, ψb])

    @test fidelity(qtraj) ≈ 1.0 atol = 1e-10
    @test fidelity(qtraj; phases = [0.0, 0.0], subsystem_levels = [2, 2]) ≈ 1.0 atol = 1e-10
    @test fidelity(qtraj; phases = [0.3, π / 2], subsystem_levels = [2, 2]) ≈ 0.5 atol =
        1e-8
    @test fidelity(qtraj; phases = [0.9, π / 2], subsystem_levels = [2, 2]) ≈ 0.5 atol =
        1e-8

    # phases without subsystem_levels is an assertion error
    @test_throws AssertionError fidelity(qtraj; phases = [0.1, 0.2])
end

@testitem "Rollouts._update_system! swaps the system across trajectory types" begin
    using LinearAlgebra

    sys = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])
    sys2 = OpenQuantumSystem(0.9 * PAULIS.Z, [PAULIS.X], [1.0])
    pulse = ZeroOrderPulse([0.5 0.5], [0.0, 1.0])

    qtraj_u = UnitaryTrajectory(sys, pulse, GATES[:X])
    Rollouts._update_system!(qtraj_u, sys2)
    @test qtraj_u.system === sys2

    ψ0 = ComplexF64[1.0, 0.0]
    ψg = ComplexF64[0.0, 1.0]
    qtraj_k = KetTrajectory(sys, pulse, ψ0, ψg)
    Rollouts._update_system!(qtraj_k, sys2)
    @test qtraj_k.system === sys2

    qtraj_mk = MultiKetTrajectory(sys, pulse, [ψ0, ψg], [ψg, ψ0])
    Rollouts._update_system!(qtraj_mk, sys2)
    @test qtraj_mk.system === sys2

    L = ComplexF64[0.1 0.0; 0.0 0.0]
    osys = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    osys2 =
        OpenQuantumSystem(0.9 * PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρg = ComplexF64[0.0 0.0; 0.0 1.0]

    qtraj_d = DensityTrajectory(osys, pulse, ρ0, ρg)
    Rollouts._update_system!(qtraj_d, osys2)
    @test qtraj_d.system === osys2

    qtraj_md = MultiDensityTrajectory(osys, pulse, [ρ0, ρg], [ρg, ρ0])
    Rollouts._update_system!(qtraj_md, osys2)
    @test qtraj_md.system === osys2

    # SamplingTrajectory delegates to its base trajectory only — the systems
    # array (the sampled variations) is deliberately untouched.
    base = UnitaryTrajectory(sys, pulse, GATES[:X])
    sampling = SamplingTrajectory(base, [sys, sys2])
    Rollouts._update_system!(sampling, sys2)
    @test sampling.base_trajectory.system === sys2
    @test sampling.systems == [sys, sys2]
end
