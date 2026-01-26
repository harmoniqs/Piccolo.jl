# ============================================================================ #
# SamplingTrajectory for Robust Optimization
# ============================================================================ #

export get_systems, get_weights

"""
    SamplingTrajectory{QT<:AbstractQuantumTrajectory} <: AbstractQuantumTrajectory

Wrapper for robust optimization over multiple systems with shared controls.

Used for sampling-based robust optimization where:
- All systems share the same control pulse
- Each system has different dynamics (e.g., parameter variations)
- Optimization minimizes weighted fidelity across all systems

This type does NOT store a NamedTrajectory - use `NamedTrajectory(sampling, N)` for conversion.

# Fields
- `base_trajectory::QT`: Base quantum trajectory (defines pulse, initial, goal)
- `systems::Vector{<:AbstractQuantumSystem}`: Multiple systems to optimize over
- `weights::Vector{Float64}`: Weights for each system in objective

# Example
```julia
sys_nom = QuantumSystem(...)
sys_variations = [QuantumSystem(...) for _ in 1:3]  # Parameter variations
qtraj = UnitaryTrajectory(sys_nom, pulse, U_goal)
sampling = SamplingTrajectory(qtraj, sys_variations, [0.5, 0.3, 0.2])

# Convert to NamedTrajectory for optimization
traj = NamedTrajectory(sampling, 51)  # Creates :Ũ⃗1, :Ũ⃗2, :Ũ⃗3
```
"""
struct SamplingTrajectory{P<:AbstractPulse, QT<:AbstractQuantumTrajectory{P}} <: AbstractQuantumTrajectory{P}
    base_trajectory::QT
    systems::Vector{<:AbstractQuantumSystem}
    weights::Vector{Float64}
end

"""
    SamplingTrajectory(base_trajectory, systems; weights=nothing)

Create a SamplingTrajectory for robust optimization.

# Arguments
- `base_trajectory`: Base quantum trajectory (defines pulse, initial, goal)
- `systems`: Vector of systems with parameter variations

# Keyword Arguments
- `weights`: Optional weights for each system (default: equal weights)
"""
function SamplingTrajectory(
    base_trajectory::QT,
    systems::Vector{<:AbstractQuantumSystem};
    weights::Union{Nothing, Vector{Float64}}=nothing
) where {P<:AbstractPulse, QT<:AbstractQuantumTrajectory{P}}
    n = length(systems)
    if isnothing(weights)
        weights = fill(1.0 / n, n)
    end
    @assert length(weights) == n "Number of weights must match number of systems"
    return SamplingTrajectory{P, QT}(base_trajectory, systems, weights)
end

# Interface implementations for SamplingTrajectory
get_system(traj::SamplingTrajectory) = get_system(traj.base_trajectory)  # Nominal system
get_pulse(traj::SamplingTrajectory) = get_pulse(traj.base_trajectory)
get_initial(traj::SamplingTrajectory) = get_initial(traj.base_trajectory)
get_goal(traj::SamplingTrajectory) = get_goal(traj.base_trajectory)
get_solution(traj::SamplingTrajectory) = get_solution(traj.base_trajectory)
duration(traj::SamplingTrajectory) = duration(traj.base_trajectory)

# Name accessors
state_name(traj::SamplingTrajectory) = state_name(traj.base_trajectory)
drive_name(traj::SamplingTrajectory) = drive_name(traj.base_trajectory)
time_name(traj::SamplingTrajectory) = time_name(traj.base_trajectory)
timestep_name(traj::SamplingTrajectory) = timestep_name(traj.base_trajectory)

"""
    state_names(sampling::SamplingTrajectory)

Get the state variable names for all systems (e.g., [:Ũ⃗1, :Ũ⃗2, :Ũ⃗3]).
"""
function state_names(traj::SamplingTrajectory)
    base = state_name(traj)
    return [Symbol(base, i) for i in 1:length(traj.systems)]
end

"""
    get_systems(sampling::SamplingTrajectory)

Get all systems in the sampling trajectory.
"""
get_systems(traj::SamplingTrajectory) = traj.systems

"""
    get_weights(sampling::SamplingTrajectory)

Get the weights for each system.
"""
get_weights(traj::SamplingTrajectory) = traj.weights

# Length for iteration
Base.length(traj::SamplingTrajectory) = length(traj.systems)

# Callable - sample base trajectory at time t
(traj::SamplingTrajectory)(t::Real) = traj.base_trajectory(t)

# ============================================================================ #
# SamplingTrajectory NamedTrajectory Conversion
# ============================================================================ #

"""
    NamedTrajectory(sampling::SamplingTrajectory, N::Int)
    NamedTrajectory(sampling::SamplingTrajectory, times::AbstractVector)

Convert a SamplingTrajectory to a NamedTrajectory for optimization.

Creates a trajectory with multiple state variables (one per system), 
all sharing the same control pulse. Each state gets a numeric suffix:
- UnitaryTrajectory base → `:Ũ⃗1`, `:Ũ⃗2`, ...
- KetTrajectory base → `:ψ̃1`, `:ψ̃2`, ...

For robust optimization, each state variable represents the evolution under
a different system (e.g., parameter variations), but all share the same controls.

# Example
```julia
# Create sampling trajectory with 3 system variations
sampling = SamplingTrajectory(base_qtraj, [sys1, sys2, sys3])

# Convert to NamedTrajectory with 51 timesteps
traj = NamedTrajectory(sampling, 51)
# Result has: :Ũ⃗1, :Ũ⃗2, :Ũ⃗3, :u, :Δt, :t
```

# Keyword Arguments
- `Δt_bounds`: Optional tuple `(lower, upper)` for timestep bounds. If provided,
  enables free-time optimization (minimum-time problems). Default: `nothing` (no bounds).
"""
function NamedTrajectory(
    sampling::SamplingTrajectory{P, <:UnitaryTrajectory{P}},
    N_or_times::Union{Int, AbstractVector{<:Real}};
    Δt_bounds::Union{Nothing, Tuple{Float64, Float64}}=nothing
) where {P<:AbstractPulse}
    base = sampling.base_trajectory
    times = _sample_times(base, N_or_times)
    T = length(times)
    n_systems = length(sampling.systems)
    snames = state_names(sampling)
    
    # Sample base trajectory for initial state data
    base_states = [base(t) for t in times]
    Ũ⃗_base = hcat([operator_to_iso_vec(U) for U in base_states]...)
    state_dim = size(Ũ⃗_base, 1)
    
    # Build state data for each system (initially all same, dynamics will differ)
    state_data = NamedTuple()
    initial_nt = NamedTuple()
    goal_nt = NamedTuple()
    bounds = NamedTuple()
    
    # All systems share initial and goal (from base trajectory)
    U_init_iso = operator_to_iso_vec(get_initial(base))
    U_goal_iso = operator_to_iso_vec(get_goal(base))
    
    for (i, name) in enumerate(snames)
        state_data = merge(state_data, _named_tuple(name => copy(Ũ⃗_base)))
        initial_nt = merge(initial_nt, _named_tuple(name => U_init_iso))
        goal_nt = merge(goal_nt, _named_tuple(name => U_goal_iso))
        bounds = merge(bounds, _named_tuple(name => (-ones(state_dim), ones(state_dim))))
    end
    
    # Get control data from base pulse
    control_data, control_names, control_bounds = _get_control_data(get_pulse(base), times, get_system(base))
    
    # Compute Δt
    Δt_diff = diff(times)
    Δt = [Δt_diff; Δt_diff[end]]
    
    # Build data
    data = merge(state_data, (; Δt = Δt, t = collect(times)), control_data)
    bounds = merge(bounds, control_bounds)
    # Add Δt bounds if provided
    if !isnothing(Δt_bounds)
        bounds = merge(bounds, (Δt = ([Δt_bounds[1]], [Δt_bounds[2]]),))
    end
    
    return NamedTrajectory(
        data;
        timestep=:Δt,
        controls=(:Δt, control_names...),
        bounds=bounds,
        initial=initial_nt,
        goal=goal_nt
    )
end

function NamedTrajectory(
    sampling::SamplingTrajectory{P, <:KetTrajectory{P}},
    N_or_times::Union{Int, AbstractVector{<:Real}};
    Δt_bounds::Union{Nothing, Tuple{Float64, Float64}}=nothing
) where {P<:AbstractPulse}
    base = sampling.base_trajectory
    times = _sample_times(base, N_or_times)
    T = length(times)
    n_systems = length(sampling.systems)
    snames = state_names(sampling)
    
    # Sample base trajectory for initial state data
    base_states = [base(t) for t in times]
    ψ̃_base = hcat([ket_to_iso(ψ) for ψ in base_states]...)
    state_dim = size(ψ̃_base, 1)
    
    # Build state data for each system
    state_data = NamedTuple()
    initial_nt = NamedTuple()
    goal_nt = NamedTuple()
    bounds = NamedTuple()
    
    # All systems share initial and goal (from base trajectory)
    ψ_init_iso = ket_to_iso(get_initial(base))
    ψ_goal_iso = ket_to_iso(get_goal(base))
    
    for (i, name) in enumerate(snames)
        state_data = merge(state_data, _named_tuple(name => copy(ψ̃_base)))
        initial_nt = merge(initial_nt, _named_tuple(name => ψ_init_iso))
        goal_nt = merge(goal_nt, _named_tuple(name => ψ_goal_iso))
        bounds = merge(bounds, _named_tuple(name => (-ones(state_dim), ones(state_dim))))
    end
    
    # Get control data from base pulse
    control_data, control_names, control_bounds = _get_control_data(get_pulse(base), times, get_system(base))
    
    # Compute Δt
    Δt_diff = diff(times)
    Δt = [Δt_diff; Δt_diff[end]]
    
    # Build data
    data = merge(state_data, (; Δt = Δt, t = collect(times)), control_data)
    bounds = merge(bounds, control_bounds)
    # Add Δt bounds if provided
    if !isnothing(Δt_bounds)
        bounds = merge(bounds, (Δt = ([Δt_bounds[1]], [Δt_bounds[2]]),))
    end
    
    return NamedTrajectory(
        data;
        timestep=:Δt,
        controls=(:Δt, control_names...),
        bounds=bounds,
        initial=initial_nt,
        goal=goal_nt
    )
end

# ============================================================================ #
# Tests for SamplingTrajectory
# ============================================================================ #

@testitem "SamplingTrajectory with UnitaryTrajectory" begin
    using LinearAlgebra
    using NamedTrajectories: NamedTrajectory
    
    # Create base system and variations
    sys_nom = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])
    sys_var1 = QuantumSystem(0.95 * PAULIS.Z, [PAULIS.X], [1.0])
    sys_var2 = QuantumSystem(1.05 * PAULIS.Z, [PAULIS.X], [1.0])
    
    # Create pulse
    T = 1.0
    times = range(0, T, length=11)
    controls = zeros(1, 11)
    pulse = LinearSplinePulse(controls, collect(times))
    
    # Create base trajectory
    U_goal = GATES[:X]
    base_qtraj = UnitaryTrajectory(sys_nom, pulse, U_goal)
    
    # Create sampling trajectory
    systems = [sys_nom, sys_var1, sys_var2]
    weights = [0.5, 0.25, 0.25]
    
    sampling = SamplingTrajectory(base_qtraj, systems; weights=weights)
    
    # Test type and accessors
    @test sampling isa AbstractQuantumTrajectory
    @test sampling isa SamplingTrajectory{<:AbstractPulse, <:UnitaryTrajectory}
    @test get_system(sampling) === sys_nom
    @test length(sampling) == 3
    @test get_systems(sampling) === systems
    @test get_weights(sampling) == weights
    @test state_names(sampling) == [:Ũ⃗1, :Ũ⃗2, :Ũ⃗3]
    @test state_name(sampling) == :Ũ⃗
    @test drive_name(sampling) == :u
    
    # Test NamedTrajectory conversion
    traj = NamedTrajectory(sampling, 11)
    @test :Ũ⃗1 ∈ traj.names
    @test :Ũ⃗2 ∈ traj.names
    @test :Ũ⃗3 ∈ traj.names
    @test :u ∈ traj.names
    @test :Δt ∈ traj.names
    @test :t ∈ traj.names
    
    # Check initial/goal propagated for each state
    for sn in state_names(sampling)
        @test haskey(traj.initial, sn)
        @test haskey(traj.goal, sn)
        @test haskey(traj.bounds, sn)
    end
end

@testitem "SamplingTrajectory with KetTrajectory" begin
    using LinearAlgebra
    using NamedTrajectories: NamedTrajectory
    
    # Create base system and variations
    sys_nom = QuantumSystem([PAULIS.X, PAULIS.Y], [1.0, 1.0])
    sys_var1 = QuantumSystem([0.95 * PAULIS.X, PAULIS.Y], [1.0, 1.0])
    sys_var2 = QuantumSystem([1.05 * PAULIS.X, PAULIS.Y], [1.0, 1.0])
    
    # Create pulse
    T = 1.0
    times = range(0, T, length=11)
    controls = zeros(2, 11)
    pulse = LinearSplinePulse(controls, collect(times))
    
    # Create base ket trajectory
    ψ_init = ComplexF64[1, 0]
    ψ_goal = ComplexF64[0, 1]
    base_qtraj = KetTrajectory(sys_nom, pulse, ψ_init, ψ_goal)
    
    # Create sampling trajectory with default weights
    systems = [sys_nom, sys_var1, sys_var2]
    sampling = SamplingTrajectory(base_qtraj, systems)
    
    # Test type and accessors
    @test sampling isa SamplingTrajectory{<:AbstractPulse, <:KetTrajectory}
    @test length(sampling) == 3
    @test get_weights(sampling) ≈ [1/3, 1/3, 1/3]  # Default equal weights
    @test state_names(sampling) == [:ψ̃1, :ψ̃2, :ψ̃3]
    @test state_name(sampling) == :ψ̃
    
    # Test NamedTrajectory conversion
    traj = NamedTrajectory(sampling, 11)
    @test :ψ̃1 ∈ traj.names
    @test :ψ̃2 ∈ traj.names
    @test :ψ̃3 ∈ traj.names
    @test :u ∈ traj.names
end

@testitem "SamplingTrajectory extract_pulse and rollout" begin
    using LinearAlgebra
    using NamedTrajectories: NamedTrajectory
    
    # Create base system and variations
    sys_nom = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])
    sys_var = QuantumSystem(0.95 * PAULIS.Z, [PAULIS.X], [1.0])
    
    # Create pulse
    T = 1.0
    times = range(0, T, length=11)
    controls = zeros(1, 11)
    pulse = LinearSplinePulse(controls, collect(times))
    
    # Create sampling trajectory
    base_qtraj = UnitaryTrajectory(sys_nom, pulse, GATES[:X])
    sampling = SamplingTrajectory(base_qtraj, [sys_nom, sys_var])
    
    # Convert to NamedTrajectory, modify, extract pulse and rollout
    traj = NamedTrajectory(sampling, 11)
    
    # Modify control values
    new_u = fill(0.5, size(traj.u))
    new_traj = NamedTrajectory(
        (; Ũ⃗1=traj.Ũ⃗1, Ũ⃗2=traj.Ũ⃗2, t=traj.t, Δt=traj.Δt, u=new_u);
        timestep=:Δt,
        controls=(:Δt, :u),
        bounds=traj.bounds,
        initial=traj.initial,
        goal=traj.goal
    )
    
    # Extract pulse and rollout
    new_pulse = extract_pulse(sampling.base_trajectory, new_traj)
    new_base_qtraj = rollout(sampling.base_trajectory, new_pulse)
    new_sampling = SamplingTrajectory(new_base_qtraj, sampling.systems; weights=sampling.weights)
    
    @test new_sampling isa SamplingTrajectory{<:AbstractPulse, <:UnitaryTrajectory}
    @test length(new_sampling) == 2
    @test get_weights(new_sampling) == sampling.weights
    
    # Check pulse was updated
    test_pulse = get_pulse(new_sampling)
    @test test_pulse(0.5)[1] ≈ 0.5
end
@testitem "SamplingTrajectory rollout!" begin
    using LinearAlgebra
    using NamedTrajectories: NamedTrajectory
    
    # Create base system and variations
    sys_nom = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])
    sys_var = QuantumSystem(0.95 * PAULIS.Z, [PAULIS.X], [1.0])
    
    # Create pulse
    T = 1.0
    times = range(0, T, length=11)
    controls = zeros(1, 11)
    pulse = LinearSplinePulse(controls, collect(times))
    
    # Create sampling trajectory
    base_qtraj = UnitaryTrajectory(sys_nom, pulse, GATES[:X])
    sampling = SamplingTrajectory(base_qtraj, [sys_nom, sys_var])
    
    # Convert to NamedTrajectory, modify controls
    traj = NamedTrajectory(sampling, 11)
    new_u = fill(0.5, size(traj.u))
    new_traj = NamedTrajectory(
        (; Ũ⃗1=traj.Ũ⃗1, Ũ⃗2=traj.Ũ⃗2, t=traj.t, Δt=traj.Δt, u=new_u);
        timestep=:Δt,
        controls=(:Δt, :u),
        bounds=traj.bounds,
        initial=traj.initial,
        goal=traj.goal
    )
    
    # Test rollout! with new pulse
    new_pulse = extract_pulse(sampling, new_traj)
    old_solution = get_solution(sampling.base_trajectory)
    rollout!(sampling, new_pulse)
    new_solution = get_solution(sampling.base_trajectory)
    
    @test new_solution !== old_solution  # Solution was updated
    @test get_pulse(sampling)(0.5)[1] ≈ 0.5  # Pulse was updated
    
    # Test rollout! with ODE parameters
    rollout!(sampling; n_points=21)
    @test length(get_solution(sampling.base_trajectory).t) == 21
end