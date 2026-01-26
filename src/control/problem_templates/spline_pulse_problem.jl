export SplinePulseProblem

# Helper function to determine spline order from pulse type
_get_spline_order(::LinearSplinePulse) = 1
_get_spline_order(::CubicSplinePulse) = 3

@doc raw"""
    SplinePulseProblem(qtraj::AbstractQuantumTrajectory{<:AbstractSplinePulse}; kwargs...)
    SplinePulseProblem(qtraj::AbstractQuantumTrajectory{<:AbstractSplinePulse}, N::Int; kwargs...)
    SplinePulseProblem(qtraj::AbstractQuantumTrajectory{<:AbstractSplinePulse}, times::AbstractVector; kwargs...)

Construct a `QuantumControlProblem` for spline-based pulse optimization.

Unlike `SmoothPulseProblem` (which uses piecewise constant controls with discrete smoothing 
variables), this problem template is designed for spline-based pulses where the derivative 
variables (`du`) are the actual spline coefficients or slopes.

## Pulse Type Semantics

**LinearSplinePulse**: The `du` variable represents the slope between knots. A `DerivativeIntegrator`
constraint enforces `du[k] = (u[k+1] - u[k]) / Δt`, making the slopes consistent with the linear
interpolation. This constraint ensures mathematical rigor while allowing slope regularization/bounds.

**CubicSplinePulse** (Hermite spline): The `du` variable is the tangent/derivative at each 
knot point, which is a true independent degree of freedom in Hermite interpolation. No 
`DerivativeIntegrator` is added - the optimizer can adjust both `:u` and `:du` independently.

## Mathematical Notes

- **LinearSplinePulse**: Always adds `:du` and `DerivativeIntegrator` to enforce slope consistency
- **CubicSplinePulse**: `:du` values are Hermite tangents (unconstrained, only regularized)

Both pulse types always have `:du` components in the trajectory, simplifying integrator implementations.

# Arguments
- `qtraj::AbstractQuantumTrajectory{<:AbstractSplinePulse}`: Quantum trajectory with spline pulse
- `N_or_times`: One of:
  - `nothing` (default): Use native knot times from spline pulse (ideal for warm-starting)
  - `N::Int`: Number of uniformly spaced timesteps
  - `times::AbstractVector`: Specific sample times

# Keyword Arguments
- `integrator::Union{Nothing, AbstractIntegrator, Vector{<:AbstractIntegrator}}=nothing`: Optional custom integrator(s). If not provided, uses `BilinearIntegrator` (which does not support global variables). A custom integrator is required when `global_names` is specified.
- `global_names::Union{Nothing, Vector{Symbol}}=nothing`: Names of global variables to optimize. Requires a custom integrator (e.g., `SplineIntegrator` from Piccolissimo) that supports global variables.
- `global_bounds::Union{Nothing, Dict{Symbol, Union{Float64, Tuple{Float64, Float64}}}}=nothing`: Bounds for global variables. Keys are variable names, values are either a scalar (symmetric bounds ±value) or a tuple (lower, upper).
- `du_bound::Float64=Inf`: Bound on derivative (slope) magnitude
- `Q::Float64=100.0`: Weight on infidelity/objective
- `R::Float64=1e-2`: Weight on regularization terms
- `R_u::Union{Float64, Vector{Float64}}=R`: Weight on control regularization
- `R_du::Union{Float64, Vector{Float64}}=R`: Weight on derivative regularization  
- `constraints::Vector{<:AbstractConstraint}=AbstractConstraint[]`: Additional constraints
- `piccolo_options::PiccoloOptions=PiccoloOptions()`: Piccolo solver options

# Returns
- `QuantumControlProblem{<:AbstractQuantumTrajectory}`: Wrapper containing trajectory and optimization problem

# Examples
```julia
# Create system and initial pulse
sys = QuantumSystem(H_drift, H_drives, drive_bounds)
pulse = CubicSplinePulse(u_init, du_init, times)
qtraj = UnitaryTrajectory(sys, pulse, U_goal)

# Use native knot structure (best for warm-starting from saved pulse)
qcp = SplinePulseProblem(qtraj; Q=100.0, du_bound=10.0)

# Or resample to different number of knots
qcp = SplinePulseProblem(qtraj, 50; Q=100.0, du_bound=10.0)

solve!(qcp; max_iter=100)
```

See also: [`SmoothPulseProblem`](@ref) for piecewise constant pulses with discrete smoothing.
"""
function SplinePulseProblem(
    qtraj::AbstractQuantumTrajectory{<:AbstractSplinePulse},
    N_or_times::Union{Nothing, Int, AbstractVector{<:Real}}=nothing;
    integrator::Union{Nothing,AbstractIntegrator,Vector{<:AbstractIntegrator}}=nothing,
    global_names::Union{Nothing,Vector{Symbol}}=nothing,
    global_bounds::Union{Nothing,Dict{Symbol,<:Union{Float64,Tuple{Float64,Float64}}}}=nothing,
    du_bound::Float64=Inf,
    Δt_bounds::Union{Nothing,Tuple{Float64,Float64}}=nothing,
    Q::Float64=100.0,
    R::Float64=1e-2,
    R_u::Union{Float64,Vector{Float64}}=R,
    R_du::Union{Float64,Vector{Float64}}=R,
    constraints::Vector{<:AbstractConstraint}=AbstractConstraint[],
    piccolo_options::PiccoloOptions=PiccoloOptions(),
)
    sys = get_system(qtraj)
    control_sym = drive_name(qtraj)
    state_sym = state_name(qtraj)

    if piccolo_options.verbose
        pulse_type = typeof(qtraj.pulse)
        println("    constructing SplinePulseProblem with $(pulse_type)...")
    end

    # Build global_data from system's global_params if present
    global_data = if !isempty(sys.global_params)
        Dict(name => [val] for (name, val) in pairs(sys.global_params))
    else
        nothing
    end

    # Convert quantum trajectory to NamedTrajectory
    # N_or_times=nothing uses native pulse knot times (preserves warm-start exactly)
    base_traj = NamedTrajectory(qtraj, N_or_times; Δt_bounds=Δt_bounds, global_data=global_data)
    N = base_traj.N  # Get actual number of timesteps

    # Always add control derivatives to trajectory
    # For CubicSplinePulse, :du is already included in the base trajectory (Hermite tangents)
    # For LinearSplinePulse, we add :du explicitly (will be constrained by DerivativeIntegrator)
    du_sym = Symbol(:d, control_sym)
    is_linear_spline = !haskey(base_traj.components, du_sym)
    
    traj = if haskey(base_traj.components, du_sym)
        # CubicSplinePulse already has derivative DOFs, but bounds default to (-Inf, Inf)
        # We need to update them if du_bound is specified
        if isfinite(du_bound)
            # Use update_bound! to set the bounds properly
            update_bound!(base_traj, du_sym, (-du_bound, du_bound))
            if piccolo_options.verbose
                println("    set du bounds to ±$du_bound for CubicSplinePulse")
            end
        end
        base_traj
    else
        # LinearSplinePulse: always add derivatives
        du_bounds_vec = isfinite(du_bound) ? fill(du_bound, sys.n_drives) : Float64[]
        if !isempty(du_bounds_vec)
            add_control_derivatives(
                base_traj,
                1;  # Only 1 derivative for spline pulses
                control_name=control_sym,
                derivative_bounds=(du_bounds_vec,)
            )
        else
            add_control_derivatives(
                base_traj,
                1;
                control_name=control_sym
            )
        end
    end

    # Initialize dynamics integrators
    if isnothing(integrator)
        if !isnothing(global_names) && !isempty(global_names)
            error(
                "global_names requires a custom integrator that supports global variables. " *
                "Use SplineIntegrator from Piccolissimo:\n" *
                "  using Piccolissimo\n" *
                "  integrator = SplineIntegrator(qtraj, N; spline_order=$(_get_spline_order(qtraj.pulse)), global_names=$global_names)\n" *
                "  qcp = SplinePulseProblem(qtraj, N; integrator=integrator, ...)"
            )
        end
        # Default to BilinearIntegrator
        default_int = BilinearIntegrator(qtraj, N)
        
        if default_int isa AbstractVector
            dynamics_integrators = AbstractIntegrator[default_int...]
        else
            dynamics_integrators = AbstractIntegrator[default_int]
        end
    elseif integrator isa AbstractIntegrator
        dynamics_integrators = AbstractIntegrator[integrator]
    else
        dynamics_integrators = AbstractIntegrator[integrator...]
    end

    # Get control names
    du_sym = Symbol(:d, control_sym)

    # Build objective: type-specific infidelity + regularization
    J = _state_objective(qtraj, traj, state_sym, Q)

    # Add regularization for control and derivative
    J += QuadraticRegularizer(control_sym, traj, R_u)
    J += QuadraticRegularizer(du_sym, traj, R_du)

    # Apply piccolo options
    J += _apply_piccolo_options(qtraj, piccolo_options, constraints, traj, state_sym)

    # Start with dynamics integrators
    integrators = copy(dynamics_integrators)

    # Add DerivativeIntegrator for LinearSplinePulse to enforce du[k] = (u[k+1] - u[k]) / Δt
    # For CubicSplinePulse, :du values are Hermite tangents (independent DOFs), not constrained
    if is_linear_spline
        push!(integrators, DerivativeIntegrator(control_sym, du_sym, traj))
        if piccolo_options.verbose
            println("    added DerivativeIntegrator for LinearSplinePulse")
        end
    end

    # Add global bounds constraints if specified
    all_constraints = copy(constraints)
    add_global_bounds_constraints!(all_constraints, global_bounds, traj; verbose=piccolo_options.verbose)

    prob = DirectTrajOptProblem(
        traj,
        J,
        integrators;
        constraints=all_constraints
    )

    return QuantumControlProblem(qtraj, prob)
end

# ============================================================================= #
# MultiKetTrajectory Method
# ============================================================================= #

"""
    SplinePulseProblem(qtraj::MultiKetTrajectory{<:AbstractSplinePulse}; kwargs...)
    SplinePulseProblem(qtraj::MultiKetTrajectory{<:AbstractSplinePulse}, N::Int; kwargs...)
    SplinePulseProblem(qtraj::MultiKetTrajectory{<:AbstractSplinePulse}, times::AbstractVector; kwargs...)

Create a spline-based trajectory optimization problem for ensemble ket state transfers.

Uses coherent fidelity objective (phases must align) for gate implementation.

# Arguments  
- `qtraj::MultiKetTrajectory{<:AbstractSplinePulse}`: Ensemble trajectory with spline pulse
- `N_or_times`: One of:
  - `nothing` (default): Use native knot times from spline pulse
  - `N::Int`: Number of uniformly spaced timesteps
  - `times::AbstractVector`: Specific sample times

# Keyword Arguments
Same as the base `SplinePulseProblem` method.
"""
function SplinePulseProblem(
    qtraj::MultiKetTrajectory{<:AbstractSplinePulse},
    N_or_times::Union{Nothing, Int, AbstractVector{<:Real}}=nothing;
    integrator::Union{Nothing,AbstractIntegrator,Vector{<:AbstractIntegrator}}=nothing,
    integrator_type::Symbol=:spline,  # :spline or :ensemble
    parallel_backend::Symbol=:manual,  # :manual (default), :threads, :gpu
    global_names::Union{Nothing,Vector{Symbol}}=nothing,
    global_bounds::Union{Nothing,Dict{Symbol,<:Union{Float64,Tuple{Float64,Float64}}}}=nothing,
    du_bound::Float64=Inf,
    Δt_bounds::Union{Nothing,Tuple{Float64,Float64}}=nothing,
    Q::Float64=100.0,
    R::Float64=1e-2,
    R_u::Union{Float64,Vector{Float64}}=R,
    R_du::Union{Float64,Vector{Float64}}=R,
    constraints::Vector{<:AbstractConstraint}=AbstractConstraint[],
    piccolo_options::PiccoloOptions=PiccoloOptions(),
)
    sys = get_system(qtraj)
    control_sym = drive_name(qtraj)
    snames = state_names(qtraj)
    weights = qtraj.weights
    goals = qtraj.goals

    if piccolo_options.verbose
        pulse_type = typeof(qtraj.pulse)
        println("    constructing SplinePulseProblem for MultiKetTrajectory with $(pulse_type)...")
        println("\twith $(length(qtraj.initials)) state transfers")
    end

    # Build global_data explicitly from system global_params
    global_data = if !isempty(sys.global_params)
        Dict(name => [val] for (name, val) in pairs(sys.global_params))
    else
        nothing
    end

    # Convert quantum trajectory to NamedTrajectory
    # N_or_times=nothing uses native pulse knot times (preserves warm-start exactly)
    base_traj = NamedTrajectory(qtraj, N_or_times; Δt_bounds=Δt_bounds, global_data=global_data)
    N = base_traj.N  # Get actual number of timesteps

    # Always add control derivatives to trajectory
    # For CubicSplinePulse, :du is already included in the base trajectory (Hermite tangents)
    # For LinearSplinePulse, we add :du explicitly (will be constrained by DerivativeIntegrator)
    du_sym = Symbol(:d, control_sym)
    is_linear_spline = !haskey(base_traj.components, du_sym)
    
    traj = if haskey(base_traj.components, du_sym)
        # CubicSplinePulse already has derivative DOFs, but bounds default to (-Inf, Inf)
        # We need to update them if du_bound is specified
        if isfinite(du_bound)
            # Use update_bound! to set the bounds properly
            update_bound!(base_traj, du_sym, (-du_bound, du_bound))
            if piccolo_options.verbose
                println("    set du bounds to ±$du_bound for CubicSplinePulse")
            end
        end
        base_traj
    else
        # LinearSplinePulse: always add derivatives
        du_bounds_vec = isfinite(du_bound) ? fill(du_bound, sys.n_drives) : Float64[]
        if !isempty(du_bounds_vec)
            add_control_derivatives(
                base_traj,
                1;  # Only 1 derivative for spline pulses
                control_name=control_sym,
                derivative_bounds=(du_bounds_vec,)
            )
        else
            add_control_derivatives(
                base_traj,
                1;
                control_name=control_sym
            )
        end
    end

    # Initialize dynamics integrators
    if isnothing(integrator)
        # Check for global_names without integrator
        if !isnothing(global_names) && !isempty(global_names)
            error(
                "global_names requires a custom integrator that supports global variables. " *
                "Use SplineIntegrator from Piccolissimo:\n" *
                "  using Piccolissimo\n" *
                "  integrator = SplineIntegrator(qtraj, N; spline_order=$(_get_spline_order(qtraj.pulse)), global_names=$global_names)\n" *
                "  qcp = SplinePulseProblem(qtraj, N; integrator=integrator, ...)"
            )
        end
        # Choose integrator type based on integrator_type parameter
        if integrator_type == :ensemble
            dynamics_integrators = EnsembleSplineIntegrator(
                qtraj, N;
                spline_order=_get_spline_order(qtraj.pulse),
                parallel_backend=parallel_backend
            )
        else
            dynamics_integrators = BilinearIntegrator(qtraj, N)
        end
        
        if !(dynamics_integrators isa AbstractVector)
            dynamics_integrators = AbstractIntegrator[dynamics_integrators]
        else
            dynamics_integrators = AbstractIntegrator[dynamics_integrators...]
        end
    elseif integrator isa AbstractIntegrator
        dynamics_integrators = AbstractIntegrator[integrator]
    else
        dynamics_integrators = AbstractIntegrator[integrator...]
    end

    # Get control names
    du_sym = Symbol(:d, control_sym)

    # Build objective: coherent fidelity for ensemble
    J = _ensemble_ket_objective(qtraj, traj, snames, weights, goals, Q)

    # Add regularization for control and derivative
    J += QuadraticRegularizer(control_sym, traj, R_u)
    J += QuadraticRegularizer(du_sym, traj, R_du)

    # Apply piccolo options for each state
    J += _apply_piccolo_options(qtraj, piccolo_options, constraints, traj, snames)

    # Start with dynamics integrators
    integrators = copy(dynamics_integrators)

    # Add DerivativeIntegrator for LinearSplinePulse to enforce du[k] = (u[k+1] - u[k]) / Δt
    # For CubicSplinePulse, :du values are Hermite tangents (independent DOFs), not constrained
    if is_linear_spline
        push!(integrators, DerivativeIntegrator(control_sym, du_sym, traj))
        if piccolo_options.verbose
            println("    added DerivativeIntegrator for LinearSplinePulse")
        end
    end

    # Add global bounds constraints if specified
    all_constraints = copy(constraints)
    add_global_bounds_constraints!(all_constraints, global_bounds, traj; verbose=piccolo_options.verbose)

    prob = DirectTrajOptProblem(
        traj,
        J,
        integrators;
        constraints=all_constraints
    )

    return QuantumControlProblem(qtraj, prob)
end

# ============================================================================= #
# Fallback Error Method
# ============================================================================= #

"""
    SplinePulseProblem(qtraj::AbstractQuantumTrajectory, N_or_times; kwargs...)

Fallback method that provides helpful error for non-spline pulse types.
"""
function SplinePulseProblem(
    qtraj::AbstractQuantumTrajectory{P},
    N_or_times::Union{Nothing, Int, AbstractVector{<:Real}}=nothing;
    kwargs...
) where P <: AbstractPulse
    pulse_type = P
    error("""
    SplinePulseProblem is only for spline-based pulses (LinearSplinePulse, CubicSplinePulse).
    
    You provided a trajectory with pulse type: $(pulse_type)
    
    For piecewise constant pulses (ZeroOrderPulse), use SmoothPulseProblem instead:
        qcp = SmoothPulseProblem(qtraj, N; ...)
    """)
end

# ============================================================================= #
# TestItems
# ============================================================================= #

@testitem "SplinePulseProblem with LinearSplinePulse" begin
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    # Simple 2-level system
    σx = ComplexF64[0 1; 1 0]
    σz = ComplexF64[1 0; 0 -1]
    
    H_drift = 0.01 * σz
    H_drives = [σx]
    T = 10.0
    N = 51
    n_drives = 1
    
    # Create system and pulse
    sys = QuantumSystem(H_drift, H_drives, [1.0])
    
    times = collect(range(0.0, T, length=N))
    amps = 0.1 * randn(n_drives, N)
    pulse = LinearSplinePulse(amps, times)
    
    # Goal: X gate
    U_goal = ComplexF64[0 1; 1 0]
    
    # Create trajectory and problem
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)
    qcp = SplinePulseProblem(qtraj, N; Q=100.0, R=1e-2)
    
    @test qcp isa QuantumControlProblem
    @test get_trajectory(qcp) isa NamedTrajectory
    
    # Check that we only have 1 derivative level (du, not ddu)
    traj = get_trajectory(qcp)
    @test haskey(traj.components, :u) || haskey(traj.components, :θ)
    @test !haskey(traj.components, :ddu)  # No second derivative for splines
end

@testitem "SplinePulseProblem with CubicSplinePulse" begin
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    # Simple 2-level system  
    σx = ComplexF64[0 1; 1 0]
    σz = ComplexF64[1 0; 0 -1]
    
    H_drift = 0.01 * σz
    H_drives = [σx]
    T = 10.0
    N = 51
    n_drives = 1
    
    # Create system and pulse
    sys = QuantumSystem(H_drift, H_drives, [1.0])
    
    times = collect(range(0.0, T, length=N))
    amps = 0.1 * randn(n_drives, N)
    derivs = zeros(n_drives, N)  # Hermite spline with derivative DOFs
    pulse = CubicSplinePulse(amps, derivs, times)
    
    # Goal: X gate
    U_goal = ComplexF64[0 1; 1 0]
    
    # Create trajectory and problem
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)
    qcp = SplinePulseProblem(qtraj, N; Q=100.0, R=1e-2)
    
    @test qcp isa QuantumControlProblem
    @test get_trajectory(qcp) isa NamedTrajectory
end

@testitem "SplinePulseProblem du_bound enforcement for CubicSplinePulse" begin
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    # Simple 2-level system  
    σx = ComplexF64[0 1; 1 0]
    σz = ComplexF64[1 0; 0 -1]
    
    H_drift = 0.01 * σz
    H_drives = [σx]
    T = 10.0
    N = 51
    n_drives = 1
    
    # Create system and pulse
    sys = QuantumSystem(H_drift, H_drives, [1.0])
    
    times = collect(range(0.0, T, length=N))
    amps = 0.1 * randn(n_drives, N)
    derivs = zeros(n_drives, N)
    pulse = CubicSplinePulse(amps, derivs, times)
    
    U_goal = ComplexF64[0 1; 1 0]
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)
    
    # Test with du_bound specified
    du_bound = 5.0
    qcp = SplinePulseProblem(qtraj, N; Q=100.0, R=1e-2, du_bound=du_bound)
    
    traj = get_trajectory(qcp)
    
    # Verify du bounds are set correctly
    @test haskey(traj.bounds, :du)
    du_bounds = traj.bounds[:du]
    
    # Bounds are stored as (lower_vector, upper_vector) tuple
    @test length(du_bounds) == 2  # (lower, upper) tuple
    lower_bounds, upper_bounds = du_bounds
    @test length(lower_bounds) == n_drives
    @test length(upper_bounds) == n_drives
    @test all(lower_bounds .≈ -du_bound)
    @test all(upper_bounds .≈ du_bound)
    
    # Test without du_bound (should default to Inf)
    qcp_unbounded = SplinePulseProblem(qtraj, N; Q=100.0, R=1e-2)
    traj_unbounded = get_trajectory(qcp_unbounded)
    
    # Without explicit du_bound, bounds should still be set to Inf (not throw error)
    @test haskey(traj_unbounded.bounds, :du)
end

@testitem "SplinePulseProblem rejects ZeroOrderPulse" begin
    using LinearAlgebra

    σx = ComplexF64[0 1; 1 0]
    σz = ComplexF64[1 0; 0 -1]
    
    H_drift = 0.01 * σz
    H_drives = [σx]
    T = 10.0
    N = 51
    
    sys = QuantumSystem(H_drift, H_drives, [1.0])
    
    times = collect(range(0.0, T, length=N))
    pulse = ZeroOrderPulse(0.1 * randn(1, N), times)
    
    U_goal = ComplexF64[0 1; 1 0]
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)
    
    # Should error with helpful message
    @test_throws ErrorException SplinePulseProblem(qtraj, N)
end

@testitem "SplinePulseProblem with KetTrajectory" begin
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    # Simple 2-level system
    σx = ComplexF64[0 1; 1 0]
    σz = ComplexF64[1 0; 0 -1]
    
    H_drift = 0.01 * σz
    H_drives = [σx]
    T = 10.0
    N = 51
    n_drives = 1
    
    # Create system and pulse
    sys = QuantumSystem(H_drift, H_drives, [1.0])
    
    times = collect(range(0.0, T, length=N))
    amps = 0.1 * randn(n_drives, N)
    pulse = LinearSplinePulse(amps, times)
    
    # State transfer: |0⟩ → |1⟩
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    
    # Create trajectory and problem
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)
    qcp = SplinePulseProblem(qtraj, N; Q=100.0, R=1e-2)
    
    @test qcp isa QuantumControlProblem
    @test qcp.qtraj isa KetTrajectory
    @test get_trajectory(qcp) isa NamedTrajectory
    
    # Check trajectory has proper components
    traj = get_trajectory(qcp)
    @test haskey(traj.components, :ψ̃)
    @test !haskey(traj.components, :ddu)  # No second derivative for splines
end

@testitem "SplinePulseProblem with MultiKetTrajectory" begin
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    # Simple 2-level system
    σx = ComplexF64[0 1; 1 0]
    σz = ComplexF64[1 0; 0 -1]
    
    H_drift = 0.01 * σz
    H_drives = [σx]
    T = 10.0
    N = 51
    n_drives = 1
    
    # Create system and pulse
    sys = QuantumSystem(H_drift, H_drives, [1.0])
    
    times = collect(range(0.0, T, length=N))
    amps = 0.1 * randn(n_drives, N)
    pulse = LinearSplinePulse(amps, times)
    
    # Create ensemble: |0⟩ → |1⟩ and |1⟩ → |0⟩
    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    
    # Create trajectory and problem
    qtraj = MultiKetTrajectory(sys, pulse, [ψ0, ψ1], [ψ1, ψ0])
    qcp = SplinePulseProblem(qtraj, N; Q=100.0, R=1e-2)
    
    @test qcp isa QuantumControlProblem
    @test qcp.qtraj isa MultiKetTrajectory
    @test get_trajectory(qcp) isa NamedTrajectory
    
    # Check trajectory has proper components for both ensemble states
    traj = get_trajectory(qcp)
    @test haskey(traj.components, :ψ̃1)
    @test haskey(traj.components, :ψ̃2)
    @test !haskey(traj.components, :ddu)  # No second derivative for splines
    
    # Should have 2 dynamics integrators (one per state)
    dynamics_integrators = filter(i -> i isa BilinearIntegrator, qcp.prob.integrators)
    @test length(dynamics_integrators) == 2
end

@testitem "SplinePulseProblem with SamplingTrajectory" begin
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    # Simple 2-level system with parameter variation
    σx = ComplexF64[0 1; 1 0]
    σz = ComplexF64[1 0; 0 -1]
    
    H_drift = 0.01 * σz
    H_drives = [σx]
    T = 10.0
    N = 51
    n_drives = 1
    
    # Create nominal and perturbed systems
    sys_nominal = QuantumSystem(H_drift, H_drives, [1.0])
    sys_perturbed = QuantumSystem(1.1 * H_drift, H_drives, [1.0])
    
    times = collect(range(0.0, T, length=N))
    amps = 0.1 * randn(n_drives, N)
    pulse = LinearSplinePulse(amps, times)
    
    # Goal: X gate
    U_goal = ComplexF64[0 1; 1 0]
    
    # Create base trajectory and sampling trajectory
    base_qtraj = UnitaryTrajectory(sys_nominal, pulse, U_goal)
    
    # First create a SplinePulseProblem with base trajectory
    base_qcp = SplinePulseProblem(base_qtraj, N; Q=100.0, R=1e-2)
    
    # Then create SamplingProblem
    sampling_qcp = SamplingProblem(base_qcp, [sys_nominal, sys_perturbed]; Q=100.0)
    
    @test sampling_qcp isa QuantumControlProblem
    @test sampling_qcp.qtraj isa SamplingTrajectory{<:AbstractPulse, <:UnitaryTrajectory}
    
    # Check trajectory has sample states
    traj = get_trajectory(sampling_qcp)
    @test haskey(traj.components, :Ũ⃗1)
    @test haskey(traj.components, :Ũ⃗2)
    @test !haskey(traj.components, :ddu)  # No second derivative for splines
end

@testitem "SplinePulseProblem with global_bounds error handling" begin
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra
    
    # Test that global_bounds throws an informative error when global doesn't exist
    # (global_data must come from integrator - e.g., SplineIntegrator from Piccolissimo)
    
    T = 2.0
    N = 10
    
    sys = QuantumSystem(0.1 * GATES.Z, [GATES.X], [1.0])
    U_goal = GATES.X
    
    # Create pulse
    times = collect(range(0, T, N))
    pulse = CubicSplinePulse(fill(0.5, 1, N), fill(0.0, 1, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)
    
    # Attempting to use global_bounds without globals in trajectory should error
    @test_throws "Global variable :δ not found" SplinePulseProblem(
        qtraj, N;
        Q=100.0, R=1e-2,
        global_bounds=Dict(:δ => 0.5)  # δ doesn't exist in trajectory
    )
end
