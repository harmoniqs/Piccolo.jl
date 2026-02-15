# ============================================================================ #
# NamedTrajectory Conversion
# ============================================================================ #

"""
    _named_tuple(pairs...)

Create a NamedTuple from pairs of (Symbol, value). This is needed when keys are 
dynamic (stored in variables).

Example:
    name = :x
    _named_tuple(name => 1, :y => 2)  # Returns (x = 1, y = 2)
"""
function _named_tuple(pairs::Pair{Symbol}...)
    keys = Tuple(p.first for p in pairs)
    vals = Tuple(p.second for p in pairs)
    return NamedTuple{keys}(vals)
end

"""
    _sample_times(traj, N::Int)

Generate N uniformly spaced times for sampling.
"""
_sample_times(traj, N::Int) = collect(range(0.0, duration(traj), length = N))

"""
    _sample_times(traj, times::AbstractVector)

Return times as a Float64 vector.
"""
_sample_times(traj, times::AbstractVector{<:Real}) = collect(Float64, times)

"""
    _sample_times(traj, ::Nothing)

For spline pulses, extract native knot times. For other pulses, error.
"""
function _sample_times(traj, ::Nothing)
    pulse = traj.pulse
    if pulse isa AbstractSplinePulse
        return collect(Float64, get_knot_times(pulse))
    else
        error(
            "Cannot infer times for $(typeof(pulse)). " *
            "Provide N::Int for uniform sampling or times::AbstractVector.",
        )
    end
end

"""
    _add_global_data_to_kwargs(nt_kwargs, global_data)

Helper function to process global variables and add them to NamedTrajectory kwargs.
Converts Dict{Symbol, Vector} to flat vector and components NamedTuple.


# Arguments
- `nt_kwargs`: Existing NamedTuple of kwargs to merge with
- `global_data`: Dict mapping global variable names to vectors of values

# Returns
Merged NamedTuple with global_data and global_components added
"""
function _add_global_data_to_kwargs(nt_kwargs, global_data)
    if isnothing(global_data)
        return nt_kwargs
    end

    # Convert Dict{Symbol, Vector} to flat vector and components NamedTuple
    global_names = sort(collect(keys(global_data)))  # Consistent ordering
    global_vec = vcat([global_data[name] for name in global_names]...)
    offset = 0
    global_comps_list = []
    for name in global_names
        len = length(global_data[name])
        push!(global_comps_list, name => ((offset+1):(offset+len)))
        offset += len
    end

    return merge(
        nt_kwargs,
        (global_data = global_vec, global_components = NamedTuple(global_comps_list)),
    )
end

"""
    _get_drive_bounds(sys::QuantumSystem)

Extract drive bounds from system as tuple of (lower, upper) vectors.
"""
function _get_drive_bounds(sys::AbstractQuantumSystem)
    n = sys.n_drives
    lower = [b[1] for b in sys.drive_bounds]
    upper = [b[2] for b in sys.drive_bounds]
    return (lower, upper)
end

"""
    _get_control_data(pulse::Union{ZeroOrderPulse, LinearSplinePulse}, times, sys)

For ZeroOrderPulse and LinearSplinePulse: return `u` data with system bounds and boundary conditions.
Uses the pulse's drive_name to determine variable naming.

# Returns
- `data`: NamedTuple with control data
- `control_names`: Tuple of control variable names
- `bounds`: NamedTuple with control bounds
- `initial_constraints`: NamedTuple with initial value constraints
- `final_constraints`: NamedTuple with final value constraints
"""
function _get_control_data(
    pulse::Union{ZeroOrderPulse,LinearSplinePulse},
    times::AbstractVector,
    sys::AbstractQuantumSystem,
)
    u_name = drive_name(pulse)
    u = hcat([pulse(t) for t in times]...)
    u_bounds = _get_drive_bounds(sys)

    # Extract boundary conditions from pulse (default zeros if not specified)
    initial_u = pulse.initial_value
    final_u = pulse.final_value

    initial_constraints = _named_tuple(u_name => initial_u)
    final_constraints = _named_tuple(u_name => final_u)

    return _named_tuple(u_name => u),
    (u_name,),
    _named_tuple(u_name => u_bounds),
    initial_constraints,
    final_constraints
end

"""
    _get_control_data(pulse::CubicSplinePulse, times, sys)

For CubicSplinePulse: return `u` and `du` data with system bounds and boundary conditions.
Uses the pulse's drive_name to determine variable naming.

When `times` matches the pulse's native knot times, extracts stored `u` and `du` directly.
When resampling to different times, samples `u` via interpolation and computes `du` via
ForwardDiff to get the true spline derivative.

# Returns
- `data`: NamedTuple with control data
- `control_names`: Tuple of control variable names
- `bounds`: NamedTuple with control bounds
- `initial_constraints`: NamedTuple with initial value constraints
- `final_constraints`: NamedTuple with final value constraints
"""
function _get_control_data(
    pulse::CubicSplinePulse,
    times::AbstractVector,
    sys::AbstractQuantumSystem,
)
    u_name = drive_name(pulse)
    du_name = Symbol(:d, u_name)
    n = n_drives(pulse)

    # Get native knot data
    knot_times = get_knot_times(pulse)
    knot_u = get_knot_values(pulse)
    knot_du = get_knot_derivatives(pulse)

    # Check if times match native knots (same length and values)
    is_native =
        length(times) == length(knot_times) &&
        all(isapprox.(times, knot_times; atol = 1e-12))

    if is_native
        # Use stored data directly - preserves Hermite tangents exactly
        u = knot_u
        du = knot_du
    else
        # Resampling: interpolate u and compute du via ForwardDiff
        # Need to clamp times slightly inward at boundaries to avoid extrapolation errors
        t_min, t_max = first(knot_times), last(knot_times)
        eps = 1e-10 * (t_max - t_min)

        u = hcat([pulse(t) for t in times]...)
        du = hcat(
            [
                ForwardDiff.derivative(s -> pulse(s), clamp(t, t_min + eps, t_max - eps)) for t in times
            ]...,
        )
    end

    u_bounds = _get_drive_bounds(sys)
    # du bounds are typically unbounded (controlled by regularization)
    du_bounds = (-Inf * ones(n), Inf * ones(n))

    # Extract boundary conditions from pulse (default zeros if not specified)
    initial_u = pulse.initial_value
    final_u = pulse.final_value
    # For spline pulses, also constrain derivatives to zero at boundaries
    initial_du = zeros(n)
    final_du = zeros(n)

    initial_constraints = _named_tuple(u_name => initial_u, du_name => initial_du)
    final_constraints = _named_tuple(u_name => final_u, du_name => final_du)

    return _named_tuple(u_name => u, du_name => du),
    (u_name, du_name),
    _named_tuple(u_name => u_bounds, du_name => du_bounds),
    initial_constraints,
    final_constraints
end

"""
    _get_control_data(pulse::GaussianPulse, times, sys)

For GaussianPulse: sample as u values with system bounds and boundary conditions.
Uses the pulse's drive_name to determine variable naming.

# Returns
- `data`: NamedTuple with control data
- `control_names`: Tuple of control variable names
- `bounds`: NamedTuple with control bounds
- `initial_constraints`: NamedTuple with initial value constraints (empty for parametric pulses)
- `final_constraints`: NamedTuple with final value constraints (empty for parametric pulses)
"""
function _get_control_data(
    pulse::GaussianPulse,
    times::AbstractVector,
    sys::AbstractQuantumSystem,
)
    u_name = drive_name(pulse)
    u = hcat([pulse(t) for t in times]...)
    u_bounds = _get_drive_bounds(sys)
    # GaussianPulse is parametric, no boundary conditions enforced
    return _named_tuple(u_name => u),
    (u_name,),
    _named_tuple(u_name => u_bounds),
    NamedTuple(),
    NamedTuple()
end

# ============================================================================ #
# Public NamedTrajectory Conversion
# ============================================================================ #

"""
    NamedTrajectory(qtraj::UnitaryTrajectory; kwargs...)
    NamedTrajectory(qtraj::UnitaryTrajectory, N::Int; kwargs...)
    NamedTrajectory(qtraj::UnitaryTrajectory, times::AbstractVector; kwargs...)

Convert a UnitaryTrajectory to a NamedTrajectory for optimization.

The trajectory stores actual times `:t` (not timesteps `:Δt`), which is required
for time-dependent integrators used with `SplinePulseProblem`.

# Stored Variables
- `Ũ⃗`: Isomorphism of unitary (vectorized real representation)
- `u` (or custom drive_name): Control values sampled at times
- `du`: Control derivatives (only for CubicSplinePulse)
- `t`: Times

# Arguments
- `qtraj`: The quantum trajectory to convert
- `N_or_times`: One of:
  - `nothing` (default): Use native knot times from spline pulse (error for non-spline pulses)
  - `N::Int`: Number of uniformly spaced time points
  - `times::AbstractVector`: Specific times to sample at

# Keyword Arguments
- `Δt_bounds`: Optional tuple `(lower, upper)` for timestep bounds. If provided,
  enables free-time optimization (minimum-time problems). Default: `nothing` (no bounds).
- `global_data`: Optional Dict mapping global variable names to initial values (as vectors).
  Note: global variables are optimization variables without explicit box constraints.

# Returns
A NamedTrajectory suitable for direct collocation optimization.
"""
function NamedTrajectory(
    qtraj::UnitaryTrajectory,
    N_or_times::Union{Nothing,Int,AbstractVector{<:Real}} = nothing;
    Δt_bounds::Union{Nothing,Tuple{Float64,Float64}} = nothing,
    global_data::Union{Nothing,Dict{Symbol,<:AbstractVector}} = nothing,
)
    times = _sample_times(qtraj, N_or_times)
    T = length(times)
    s_name = state_name(qtraj)

    # Auto-populate global_data from system if not provided
    if isnothing(global_data) && !isempty(qtraj.system.global_params)
        global_data =
            Dict(name => [val] for (name, val) in pairs(qtraj.system.global_params))
    end

    # Sample unitary states
    states = [qtraj(t) for t in times]
    Ũ⃗ = hcat([operator_to_iso_vec(U) for U in states]...)

    # Get control data based on pulse type
    control_data, control_names, control_bounds, control_initial, control_final =
        _get_control_data(qtraj.pulse, times, qtraj.system)

    # State dimension
    state_dim = size(Ũ⃗, 1)

    # Compute Δt from times (pad to length T by repeating last value)
    Δt_diff = diff(times)
    Δt = [Δt_diff; Δt_diff[end]]

    # Build data NamedTuple with Δt as timestep and t for reference
    data = merge(_named_tuple(s_name => Ũ⃗, :Δt => Δt, :t => collect(times)), control_data)

    # Initial and final conditions (merge state and control boundaries)
    initial =
        merge(_named_tuple(s_name => operator_to_iso_vec(qtraj.initial)), control_initial)
    final_nt = merge(control_final)  # Control boundaries go to final (hard constraints)
    U_goal = qtraj.goal isa EmbeddedOperator ? qtraj.goal.operator : qtraj.goal
    goal_nt = _named_tuple(s_name => operator_to_iso_vec(U_goal))  # State goal (soft constraint)

    # Bounds (state bounded, controls bounded by system, optionally timestep bounded)
    bounds =
        merge(_named_tuple(s_name => (-ones(state_dim), ones(state_dim))), control_bounds)
    # Add Δt bounds if provided
    if !isnothing(Δt_bounds)
        bounds = merge(bounds, (Δt = ([Δt_bounds[1]], [Δt_bounds[2]]),))
    end

    # Build kwargs for NamedTrajectory constructor
    nt_kwargs = (
        timestep = :Δt,
        controls = (:Δt, control_names...),
        bounds = bounds,
        initial = initial,
        final = final_nt,
        goal = goal_nt,
    )

    # Add global variables if provided
    nt_kwargs = _add_global_data_to_kwargs(nt_kwargs, global_data)

    return NamedTrajectory(data; nt_kwargs...)
end

"""
    NamedTrajectory(qtraj::KetTrajectory; kwargs...)
    NamedTrajectory(qtraj::KetTrajectory, N::Int; kwargs...)
    NamedTrajectory(qtraj::KetTrajectory, times::AbstractVector; kwargs...)

Convert a KetTrajectory to a NamedTrajectory for optimization.

# Stored Variables
- `ψ̃`: Isomorphism of ket state (real representation)
- `u` (or custom drive_name): Control values sampled at times
- `du`: Control derivatives (only for CubicSplinePulse)
- `t`: Times

# Arguments
- `N_or_times`: One of:
  - `nothing` (default): Use native knot times from spline pulse
  - `N::Int`: Number of uniformly spaced time points
  - `times::AbstractVector`: Specific times to sample at

# Keyword Arguments
- `Δt_bounds`: Optional tuple `(lower, upper)` for timestep bounds. If provided,
  enables free-time optimization (minimum-time problems). Default: `nothing` (no bounds).
- `global_data`: Optional Dict mapping global variable names to initial values (as vectors).
  Note: global variables are optimization variables without explicit box constraints.
"""
function NamedTrajectory(
    qtraj::KetTrajectory,
    N_or_times::Union{Nothing,Int,AbstractVector{<:Real}} = nothing;
    Δt_bounds::Union{Nothing,Tuple{Float64,Float64}} = nothing,
    global_data::Union{Nothing,Dict{Symbol,<:AbstractVector}} = nothing,
)
    times = _sample_times(qtraj, N_or_times)
    T = length(times)
    s_name = state_name(qtraj)

    # Auto-populate global_data from system if not provided
    if isnothing(global_data) && !isempty(qtraj.system.global_params)
        global_data =
            Dict(name => [val] for (name, val) in pairs(qtraj.system.global_params))
    end

    # Sample ket states
    states = [qtraj(t) for t in times]
    ψ̃ = hcat([ket_to_iso(ψ) for ψ in states]...)

    # Get control data based on pulse type
    control_data, control_names, control_bounds, control_initial, control_final =
        _get_control_data(qtraj.pulse, times, qtraj.system)

    # State dimension
    state_dim = size(ψ̃, 1)

    # Compute Δt from times (pad to length T by repeating last value)
    Δt_diff = diff(times)
    Δt = [Δt_diff; Δt_diff[end]]

    # Build data with Δt as timestep and t for reference
    data = merge(_named_tuple(s_name => ψ̃, :Δt => Δt, :t => collect(times)), control_data)

    # Initial, goal, bounds (merge state and control boundaries)
    initial = merge(_named_tuple(s_name => ket_to_iso(qtraj.initial)), control_initial)
    final_nt = merge(control_final)  # Control boundaries go to final (hard constraints)
    goal_nt = _named_tuple(s_name => ket_to_iso(qtraj.goal))  # State goal (soft constraint)
    bounds =
        merge(_named_tuple(s_name => (-ones(state_dim), ones(state_dim))), control_bounds)
    # Add Δt bounds if provided
    if !isnothing(Δt_bounds)
        bounds = merge(bounds, (Δt = ([Δt_bounds[1]], [Δt_bounds[2]]),))
    end

    # Build kwargs for NamedTrajectory constructor
    nt_kwargs = (
        timestep = :Δt,
        controls = (:Δt, control_names...),
        bounds = bounds,
        initial = initial,
        final = final_nt,
        goal = goal_nt,
    )

    # Add global variables if provided
    nt_kwargs = _add_global_data_to_kwargs(nt_kwargs, global_data)

    return NamedTrajectory(data; nt_kwargs...)
end

"""
    NamedTrajectory(qtraj::MultiKetTrajectory; kwargs...)
    NamedTrajectory(qtraj::MultiKetTrajectory, N::Int; kwargs...)
    NamedTrajectory(qtraj::MultiKetTrajectory, times::AbstractVector; kwargs...)

Convert an MultiKetTrajectory to a NamedTrajectory for optimization.

# Stored Variables
- `ψ̃1`, `ψ̃2`, ...: Isomorphisms of each ket state
- `u` (or custom drive_name): Control values sampled at times
- `du`: Control derivatives (only for CubicSplinePulse)
- `t`: Times

# Arguments
- `N_or_times`: One of:
  - `nothing` (default): Use native knot times from spline pulse
  - `N::Int`: Number of uniformly spaced time points
  - `times::AbstractVector`: Specific times to sample at

# Keyword Arguments
- `Δt_bounds`: Optional tuple `(lower, upper)` for timestep bounds. If provided,
  enables free-time optimization (minimum-time problems). Default: `nothing` (no bounds).
- `global_data`: Optional Dict mapping global variable names to initial values (as vectors).
  Note: global variables are optimization variables without explicit box constraints.
"""
function NamedTrajectory(
    qtraj::MultiKetTrajectory,
    N_or_times::Union{Nothing,Int,AbstractVector{<:Real}} = nothing;
    Δt_bounds::Union{Nothing,Tuple{Float64,Float64}} = nothing,
    global_data::Union{Nothing,Dict{Symbol,<:AbstractVector}} = nothing,
)
    times = _sample_times(qtraj, N_or_times)
    T = length(times)
    n_states = length(qtraj)
    state_prefix = state_name(qtraj)

    # Auto-populate global_data from system if not provided
    if isnothing(global_data) && !isempty(qtraj.system.global_params)
        global_data =
            Dict(name => [val] for (name, val) in pairs(qtraj.system.global_params))
    end

    # Sample all ket states
    state_data = NamedTuple()
    initial_nt = NamedTuple()
    goal_nt = NamedTuple()
    bounds = NamedTuple()

    for i = 1:n_states
        name = Symbol(state_prefix, i)
        sol = qtraj[i]
        states = [sol(t) for t in times]
        ψ̃ = hcat([ket_to_iso(ψ) for ψ in states]...)
        state_dim = size(ψ̃, 1)

        state_data = merge(state_data, _named_tuple(name => ψ̃))
        initial_nt = merge(initial_nt, _named_tuple(name => ket_to_iso(qtraj.initials[i])))
        goal_nt = merge(goal_nt, _named_tuple(name => ket_to_iso(qtraj.goals[i])))
        bounds = merge(bounds, _named_tuple(name => (-ones(state_dim), ones(state_dim))))
    end

    # Get control data
    control_data, control_names, control_bounds, control_initial, control_final =
        _get_control_data(qtraj.pulse, times, qtraj.system)

    # Merge control boundaries with state boundaries
    initial_nt = merge(initial_nt, control_initial)
    final_nt = merge(control_final)  # Control boundaries go to final (hard constraints)
    # goal_nt contains only state goals (soft constraints)

    # Compute Δt from times (pad to length T by repeating last value)
    Δt_diff = diff(times)
    Δt = [Δt_diff; Δt_diff[end]]

    # Build data with Δt as timestep and t for reference
    data = merge(state_data, (; Δt = Δt, t = collect(times)), control_data)
    bounds = merge(bounds, control_bounds)
    # Add Δt bounds if provided
    if !isnothing(Δt_bounds)
        bounds = merge(bounds, (Δt = ([Δt_bounds[1]], [Δt_bounds[2]]),))
    end

    # Build kwargs for NamedTrajectory constructor
    nt_kwargs = (
        timestep = :Δt,
        controls = (:Δt, control_names...),
        bounds = bounds,
        initial = initial_nt,
        final = final_nt,
        goal = goal_nt,
    )

    # Add global variables if provided
    nt_kwargs = _add_global_data_to_kwargs(nt_kwargs, global_data)

    return NamedTrajectory(data; nt_kwargs...)
end

"""
    NamedTrajectory(qtraj::DensityTrajectory; kwargs...)
    NamedTrajectory(qtraj::DensityTrajectory, N::Int; kwargs...)
    NamedTrajectory(qtraj::DensityTrajectory, times::AbstractVector; kwargs...)

Convert a DensityTrajectory to a NamedTrajectory for optimization.

Uses the compact density isomorphism (n² real parameters) which exploits the
Hermiticity of density matrices. This halves the state dimension compared to
the full [vec(Re(ρ)); vec(Im(ρ))] representation.

# Stored Variables
- `ρ⃗̃`: Compact isomorphism of the density matrix (n² real parameters)
- `u` (or custom drive_name): Control values sampled at times
- `du`: Control derivatives (only for CubicSplinePulse)
- `t`: Times

# Arguments
- `N_or_times`: One of:
  - `nothing` (default): Use native knot times from spline pulse
  - `N::Int`: Number of uniformly spaced time points
  - `times::AbstractVector`: Specific times to sample at

# Keyword Arguments
- `Δt_bounds`: Optional tuple `(lower, upper)` for timestep bounds. If provided,
  enables free-time optimization (minimum-time problems). Default: `nothing` (no bounds).
"""
function NamedTrajectory(
    qtraj::DensityTrajectory,
    N_or_times::Union{Nothing,Int,AbstractVector{<:Real}} = nothing;
    Δt_bounds::Union{Nothing,Tuple{Float64,Float64}} = nothing,
)
    times = _sample_times(qtraj, N_or_times)
    T = length(times)
    sname = state_name(qtraj)

    # Sample density matrices and convert to compact isomorphism (n² real params)
    states = [qtraj(t) for t in times]
    ρ̃ = hcat([density_to_compact_iso(ρ) for ρ in states]...)
    state_dim = size(ρ̃, 1)

    # Get control data
    control_data, control_names, control_bounds =
        _get_control_data(qtraj.pulse, times, qtraj.system)

    # Compute Δt from times (pad to length T by repeating last value)
    Δt_diff = diff(times)
    Δt = [Δt_diff; Δt_diff[end]]

    # Build data with Δt as timestep and t for reference
    data = merge(_named_tuple(sname => ρ̃), (; Δt = Δt, t = collect(times)), control_data)

    # Bounds on compact iso components: diagonal Re entries in [0,1],
    # off-diagonal Re/Im entries in [-1,1]. Use generous [-1,1] for all.
    bounds =
        merge(_named_tuple(sname => (-ones(state_dim), ones(state_dim))), control_bounds)

    # Add Δt bounds if provided
    if !isnothing(Δt_bounds)
        bounds = merge(bounds, (Δt = ([Δt_bounds[1]], [Δt_bounds[2]]),))
    end

    # Initial and goal in compact isomorphism
    initial = _named_tuple(sname => density_to_compact_iso(qtraj.initial))
    goal_nt = _named_tuple(sname => density_to_compact_iso(qtraj.goal))

    return NamedTrajectory(
        data;
        timestep = :Δt,
        controls = (:Δt, control_names...),
        bounds = bounds,
        initial = initial,
        goal = goal_nt,
    )
end

# ============================================================================ #
# Tests
# ============================================================================ #

@testitem "UnitaryTrajectory" begin
    using LinearAlgebra
    using NamedTrajectories

    # Simple 2-level system
    system = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])

    # Create trajectory with zero controls
    T = 1.0
    X_gate = [0.0 1.0; 1.0 0.0]
    qtraj = UnitaryTrajectory(system, X_gate, T)

    # Test interface
    @test get_system(qtraj) === system
    @test get_goal(qtraj) === X_gate
    @test state_name(qtraj) == :Ũ⃗
    @test drive_name(qtraj) == :u
    @test time_name(qtraj) == :t
    @test timestep_name(qtraj) == :Δt
    @test duration(qtraj) ≈ T

    # Test callable
    U0 = qtraj(0.0)
    @test U0 ≈ Matrix{ComplexF64}(I, 2, 2)

    # Convert to NamedTrajectory
    N = 11
    traj = NamedTrajectory(qtraj, N)
    @test length(traj.t) == N
    @test haskey(traj.components, :Ũ⃗)
    @test haskey(traj.components, :u)
end

@testitem "UnitaryTrajectory with EmbeddedOperator" begin
    using LinearAlgebra
    using NamedTrajectories

    # Create embedded gate in 3-level system
    X_embedded = EmbeddedOperator(:X, [1, 2], 3)

    # 3-level system
    H_drift = diagm(ComplexF64[1.0, 0.0, -1.0])
    H_drive1 = zeros(ComplexF64, 3, 3)
    H_drive2 = zeros(ComplexF64, 3, 3)
    H_drive1[1, 2] = H_drive1[2, 1] = 1.0
    H_drive2[2, 3] = H_drive2[3, 2] = 1.0
    system = QuantumSystem(H_drift, [H_drive1, H_drive2], [1.0, 1.0])

    # Create trajectory
    T = 1.0
    qtraj = UnitaryTrajectory(system, X_embedded, T)

    # Test that goal is preserved
    @test qtraj.goal === X_embedded

    # Convert to NamedTrajectory
    N = 11
    traj = NamedTrajectory(qtraj, N)
    @test haskey(traj.components, :Ũ⃗)
end

@testitem "KetTrajectory" begin
    using LinearAlgebra
    using NamedTrajectories

    # Simple 2-level system
    system = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])

    # Create trajectory
    T = 1.0
    ψ0 = [1.0 + 0im, 0.0]
    ψg = [0.0, 1.0 + 0im]
    qtraj = KetTrajectory(system, ψ0, ψg, T)

    # Test interface
    @test get_system(qtraj) === system
    @test get_initial(qtraj) == ψ0
    @test get_goal(qtraj) == ψg
    @test state_name(qtraj) == :ψ̃
    @test duration(qtraj) ≈ T

    # Convert to NamedTrajectory
    N = 11
    traj = NamedTrajectory(qtraj, N)
    @test length(traj.t) == N
    @test haskey(traj.components, :ψ̃)
end

@testitem "MultiKetTrajectory" begin
    using LinearAlgebra
    using NamedTrajectories

    # Simple 2-level system
    system = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])

    # Create trajectory with multiple states
    T = 1.0
    initials = [[1.0 + 0im, 0.0], [0.0, 1.0 + 0im]]
    goals = [[0.0, 1.0 + 0im], [1.0 + 0im, 0.0]]
    qtraj = MultiKetTrajectory(system, initials, goals, T)

    # Test interface
    @test length(qtraj) == 2
    @test state_name(qtraj) == :ψ̃
    @test state_names(qtraj) == [:ψ̃1, :ψ̃2]
    @test get_initial(qtraj) == qtraj.initials
    @test get_goal(qtraj) == qtraj.goals

    # Test callable and indexing
    states_0 = qtraj(0.0)
    @test length(states_0) == 2

    # Convert to NamedTrajectory
    N = 11
    traj = NamedTrajectory(qtraj, N)
    @test haskey(traj.components, :ψ̃1)
    @test haskey(traj.components, :ψ̃2)
end

@testitem "DensityTrajectory" begin
    using LinearAlgebra
    using NamedTrajectories

    # Simple 2-level open system
    L = ComplexF64[0.1 0.0; 0.0 0.0]  # Decay operator
    system = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])

    # Create trajectory
    T = 1.0
    ρ0 = [1.0 0.0; 0.0 0.0] .+ 0im
    ρg = [0.0 0.0; 0.0 1.0] .+ 0im
    qtraj = DensityTrajectory(system, ρ0, ρg, T)

    # Test interface
    @test get_system(qtraj) === system
    @test get_initial(qtraj) == ρ0
    @test get_goal(qtraj) == ρg
    @test state_name(qtraj) == :ρ⃗̃
    @test duration(qtraj) ≈ T

    # Convert to NamedTrajectory
    N = 11
    traj = NamedTrajectory(qtraj, N)
    @test haskey(traj.components, :ρ⃗̃)

    # Compact isomorphism: state dim should be n² (not 2n²)
    n = system.levels
    @test traj.dims[:ρ⃗̃] == n^2

    # Round-trip: compact iso → density → compact iso should be consistent
    ρ̃_init = traj[:ρ⃗̃][:, 1]
    ρ_reconstructed = compact_iso_to_density(ρ̃_init)
    @test ρ_reconstructed ≈ ρ0

    # Goal should also round-trip
    ρ̃_goal = traj.goal[:ρ⃗̃]
    ρg_reconstructed = compact_iso_to_density(ρ̃_goal)
    @test ρg_reconstructed ≈ ρg
end

@testitem "NamedTrajectory conversion with specific times" begin
    using LinearAlgebra
    using NamedTrajectories

    # System setup
    system = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])

    # Create trajectory
    T = 1.0
    X_gate = [0.0 1.0; 1.0 0.0]
    qtraj = UnitaryTrajectory(system, X_gate, T)

    # Convert with specific times
    times = [0.0, 0.2, 0.5, 0.8, 1.0]
    traj = NamedTrajectory(qtraj, times)

    @test length(traj.t) == 5
    @test all(vec(traj.t) .≈ times)
end

@testitem "NamedTrajectory with custom drive_name" begin
    using LinearAlgebra
    using NamedTrajectories

    # System setup
    system = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])

    # Create pulse with custom name
    T = 1.0
    times = [0.0, T]
    controls = zeros(1, 2)
    pulse = ZeroOrderPulse(controls, times; drive_name = :a)

    X_gate = [0.0 1.0; 1.0 0.0]
    qtraj = UnitaryTrajectory(system, pulse, X_gate)

    # Check that custom name propagates
    @test drive_name(qtraj) == :a

    traj = NamedTrajectory(qtraj, 11)
    @test haskey(traj.components, :a)
    @test !haskey(traj.components, :u)
end

@testitem "extract_pulse and rollout trajectory from NamedTrajectory" begin
    using LinearAlgebra
    using NamedTrajectories

    # System setup
    system = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])

    # Create trajectory
    T = 1.0
    X_gate = [0.0 1.0; 1.0 0.0]
    qtraj = UnitaryTrajectory(system, X_gate, T)

    # Convert to NamedTrajectory
    N = 11
    traj = NamedTrajectory(qtraj, N)

    # Modify controls
    traj.u .= 0.5 * randn(1, N)

    # Extract pulse and rollout
    pulse = extract_pulse(qtraj, traj)
    qtraj2 = rollout(qtraj, pulse)

    @test qtraj2 isa UnitaryTrajectory
    @test qtraj2.system === system
    @test qtraj2.goal === X_gate
end

@testitem "NamedTrajectory with Δt_bounds" begin
    using LinearAlgebra
    using NamedTrajectories

    # System setup
    system = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])

    # Create trajectory
    T = 1.0
    X_gate = [0.0 1.0; 1.0 0.0]
    qtraj = UnitaryTrajectory(system, X_gate, T)

    # Convert with Δt bounds
    N = 11
    traj = NamedTrajectory(qtraj, N; Δt_bounds = (0.05, 0.2))

    # Check that Δt bounds are set
    @test haskey(traj.bounds, :Δt)
    @test traj.bounds[:Δt][1] == [0.05]
    @test traj.bounds[:Δt][2] == [0.2]
end

@testitem "_add_global_data_to_kwargs helper" begin
    using .QuantumTrajectories: _add_global_data_to_kwargs

    # Test with nothing (should return unchanged)
    nt_kwargs = (timestep = :Δt, controls = (:u,))
    result = _add_global_data_to_kwargs(nt_kwargs, nothing)
    @test result === nt_kwargs
    @test !haskey(result, :global_data)
    @test !haskey(result, :global_components)

    # Test with single-dimensional globals
    global_data = Dict(:δ => [0.5], :Ω => [1.0])
    result = _add_global_data_to_kwargs(nt_kwargs, global_data)
    @test haskey(result, :global_data)
    @test haskey(result, :global_components)
    @test result.global_data == [1.0, 0.5]  # Sorted by Symbol (Unicode codepoint): Ω, δ
    @test result.global_components.Ω == 1:1
    @test result.global_components.δ == 2:2

    # Test with multi-dimensional globals
    global_data_multi = Dict(:δ => [0.5], :α => [1.0, 2.0, 3.0], :Ω => [1.5])
    result = _add_global_data_to_kwargs(nt_kwargs, global_data_multi)
    @test result.global_data == [1.5, 1.0, 2.0, 3.0, 0.5]  # Sorted by Symbol (Unicode codepoint): Ω, α, δ
    @test result.global_components.Ω == 1:1
    @test result.global_components.α == 2:4
    @test result.global_components.δ == 5:5

    # Verify original kwargs are preserved
    @test result.timestep == :Δt
    @test result.controls == (:u,)
end

@testitem "NamedTrajectory with global_data parameter" begin
    using NamedTrajectories

    # Create a system with global parameters
    H_drives = [PAULIS[:X], PAULIS[:Y]]
    sys = QuantumSystem(H_drives, [1.0, 1.0]; global_params = (δ = 0.5, Ω = 1.0))

    # Create trajectory
    pulse = ZeroOrderPulse([0.5 0.3; 0.5 0.3], [0.0, 1.0])
    U_goal = PAULIS[:X]
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)

    # Convert to NamedTrajectory - should auto-populate from system.global_params
    traj_auto = NamedTrajectory(qtraj, 5)
    @test hasfield(typeof(traj_auto), :global_components)
    @test haskey(traj_auto.global_components, :δ)
    @test haskey(traj_auto.global_components, :Ω)
    @test traj_auto.global_data[traj_auto.global_components.δ][1] == 0.5
    @test traj_auto.global_data[traj_auto.global_components.Ω][1] == 1.0

    # Convert with explicit global_data (should override system.global_params)
    custom_global = Dict(:δ => [0.8], :Ω => [1.5])
    traj_custom = NamedTrajectory(qtraj, 5; global_data = custom_global)
    @test traj_custom.global_data[traj_custom.global_components.δ][1] == 0.8
    @test traj_custom.global_data[traj_custom.global_components.Ω][1] == 1.5

    # Test with KetTrajectory
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    qtraj_ket = KetTrajectory(sys, pulse, ψ_init, ψ_goal)
    traj_ket = NamedTrajectory(qtraj_ket, 5)
    @test hasfield(typeof(traj_ket), :global_components)
    @test traj_ket.global_data[traj_ket.global_components.δ][1] == 0.5

    # Test with system without global_params
    sys_no_globals = QuantumSystem(H_drives, [1.0, 1.0])
    qtraj_no_g = UnitaryTrajectory(sys_no_globals, pulse, U_goal)
    traj_no_g = NamedTrajectory(qtraj_no_g, 5)
    # Should work without error, just no global components
    @test !hasfield(typeof(traj_no_g), :global_components) ||
          isempty(traj_no_g.global_components)
end
