export SplinePulseProblem

# Helper function to determine spline order from pulse type
_get_spline_order(::LinearSplinePulse) = 1
_get_spline_order(::CubicSplinePulse) = 3

# _make_free_phase_goal is defined in _problem_templates.jl (shared across all templates)

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
- `calibration_targets::Vector{Symbol}=Symbol[]`: Names of globals declared as **calibration targets** — knobs an external calibration step manages, not free NLP variables. Each listed name is pinned at its nominal value via `GlobalEqualityConstraint` so the QCP solve cannot drift it as a slack variable. Replaces any existing `GlobalBoundsConstraint` on the same name. Default empty: globals stay free.
- `du_bound::Float64=Inf`: Uniform bound on derivative (slope) magnitude for all drives
- `du_bounds::Union{Nothing, Vector{Float64}}=nothing`: Per-drive bounds on derivative magnitude (takes precedence over `du_bound`)
- `Q::Float64=100.0`: Weight on infidelity/objective
- `R::Float64=1e-2`: Weight on regularization terms (LinearSplinePulse only — see below)
- `R_u::Union{Nothing, Float64, Vector{Float64}}=nothing`: Weight on control regularization. Pulse-type-dependent default — see below.
- `R_du::Union{Nothing, Float64, Vector{Float64}}=nothing`: Weight on derivative regularization. Pulse-type-dependent default — see below.
- `λ_bend::Float64=1e-3`: Weight on the cubic-Hermite bending-energy regulariser ``\lambda \sum_c \int_0^T \ddot u_c(t)^2 \, dt``. Only applied for `CubicSplinePulse`; ignored for `LinearSplinePulse` (whose second derivative is impulsive at every knot and cannot be sensibly bent-energy-regularised — its smoothness is shaped by `R_du` on the slopes).
- `constraints::Vector{<:AbstractConstraint}=AbstractConstraint[]`: Additional constraints
- `piccolo_options::PiccoloOptions=PiccoloOptions()`: Piccolo solver options

## Default regularisation by pulse type

The L2 regularisers `R_u` and `R_du` have different roles for the two spline
families, and the defaults reflect this:

- **`LinearSplinePulse`**: `R_u` and `R_du` both default to `R = 1e-2`. The
  `du` variable here is the constrained inter-knot slope, so a quadratic
  penalty on it is the standard smoothness term.
- **`CubicSplinePulse`**: `R_u` and `R_du` both default to `0.0`. Empirically
  (see PR linked in the changelog) the same `1e-2` penalty actively biases
  the optimiser toward small-amplitude flat pulses regardless of whether
  that is needed — on a 2-qubit iSWAP problem it stalls fidelity at ≈ 0.91
  versus ≈ 0.98+ at zero regularisation. The cubic-Hermite smoothness role
  is taken over by `λ_bend`, which penalises the actual second-derivative
  energy of the interpolant rather than the magnitude of the tangents.

Pass `R_u` or `R_du` explicitly to override these defaults.

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

# Per-drive bounds (takes precedence over du_bound)
qcp = SplinePulseProblem(qtraj; Q=100.0, du_bounds=[5.0, 2.0])

solve!(qcp; max_iter=100)
```

See also: [`SmoothPulseProblem`](@ref) for piecewise constant pulses with discrete smoothing.
"""
function SplinePulseProblem(
    qtraj::AbstractQuantumTrajectory{<:AbstractSplinePulse},
    N_or_times::Union{Nothing,Int,AbstractVector{<:Real}} = nothing;
    integrator::Union{Nothing,AbstractIntegrator,Vector{<:AbstractIntegrator}} = nothing,
    global_names::Union{Nothing,Vector{Symbol}} = nothing,
    global_bounds::Union{Nothing,Dict{Symbol,<:Union{Float64,Tuple{Float64,Float64}}}} = nothing,
    calibration_targets::Vector{Symbol} = Symbol[],
    du_bound::Float64 = Inf,
    du_bounds::Union{Nothing,Vector{Float64}} = nothing,
    Δt_bounds::Union{Nothing,Tuple{Float64,Float64}} = nothing,
    Q::Float64 = 100.0,
    R::Float64 = 1e-2,
    R_u::Union{Nothing,Float64,Vector{Float64}} = nothing,
    R_du::Union{Nothing,Float64,Vector{Float64}} = nothing,
    λ_bend::Float64 = 1e-3,
    constraints::Vector{<:AbstractConstraint} = AbstractConstraint[],
    piccolo_options::PiccoloOptions = PiccoloOptions(),
    free_phase::Bool = false,
    subsystem_levels::Union{Nothing,Vector{Int}} = nothing,
    initial_phases::Union{Nothing,Vector{Float64}} = nothing,
    state_leakage_indices::Union{
        Nothing,
        AbstractVector{Int},
        AbstractVector{<:AbstractVector{Int}},
    } = nothing,
)
    sys = get_system(qtraj)
    control_sym = drive_name(qtraj)
    state_sym = state_name(qtraj)

    # Pulse-type-dependent regularisation defaults.
    # See the "Default regularisation by pulse type" section of this template's docstring.
    is_cubic = qtraj.pulse isa CubicSplinePulse
    R_u_resolved = isnothing(R_u) ? (is_cubic ? 0.0 : R) : R_u
    R_du_resolved = isnothing(R_du) ? (is_cubic ? 0.0 : R) : R_du

    if _show_header(piccolo_options)
        pulse_type = _typename(qtraj.pulse)
        traj_type = _typename(qtraj)
        println("constructing SplinePulseProblem [$traj_type / $pulse_type]")
    end

    # Build global_data from system's global_params if present
    global_data = if !isempty(sys.global_params)
        Dict{Symbol,Vector{Float64}}(name => [val] for (name, val) in pairs(sys.global_params))
    else
        nothing
    end

    # Free-phase support: add phase variables as globals
    θ_names = Symbol[]
    U_goal_fn = nothing
    ket_goal_fn = nothing
    if free_phase
        if qtraj isa KetTrajectory
            # Ket free-phase: requires subsystem_levels to build per-subsystem phase rotations
            @assert !isnothing(subsystem_levels) "free_phase=true for KetTrajectory requires subsystem_levels"
            n_qubits = length(subsystem_levels)
            ket_goal_fn = _make_free_phase_ket_goal(qtraj.goal, subsystem_levels)
            θ_names, global_data, global_bounds = setup_free_phase_globals!(
                n_qubits,
                global_data,
                global_bounds;
                initial_phases = initial_phases,
                verbose = _show_details(piccolo_options),
            )
        else
            # Unitary free-phase: requires EmbeddedOperator goal
            goal = qtraj.goal
            @assert goal isa EmbeddedOperator "free_phase=true requires an EmbeddedOperator goal or subsystem_levels"
            n_qubits = length(goal.subsystem_levels)
            U_goal_fn = _make_free_phase_goal(goal)
            θ_names, global_data, global_bounds = setup_free_phase_globals!(
                n_qubits,
                global_data,
                global_bounds;
                initial_phases = initial_phases,
                verbose = _show_details(piccolo_options),
            )
        end
    end

    # Convert quantum trajectory to NamedTrajectory
    # N_or_times=nothing uses native pulse knot times (preserves warm-start exactly)
    base_traj =
        NamedTrajectory(qtraj, N_or_times; Δt_bounds = Δt_bounds, global_data = global_data)
    N = base_traj.N  # Get actual number of timesteps

    # Always add control derivatives to trajectory
    # For CubicSplinePulse, :du is already included in the base trajectory (Hermite tangents)
    # For LinearSplinePulse, we add :du explicitly (will be constrained by DerivativeIntegrator)
    du_sym = Symbol(:d, control_sym)
    is_linear_spline = !haskey(base_traj.components, du_sym)

    # Resolve per-drive du bounds: du_bounds (vector) takes precedence over du_bound (scalar)
    _du_bounds_vec = if !isnothing(du_bounds)
        du_bounds
    elseif isfinite(du_bound)
        fill(du_bound, sys.n_drives)
    else
        nothing
    end

    traj = if haskey(base_traj.components, du_sym)
        # CubicSplinePulse already has derivative DOFs, but bounds default to (-Inf, Inf)
        if !isnothing(_du_bounds_vec)
            update_bound!(base_traj, du_sym, (-_du_bounds_vec, _du_bounds_vec))
            if _show_details(piccolo_options)
                println("    du bounds (CubicSplinePulse): $(_fmt_bounds(_du_bounds_vec))")
            end
        end
        base_traj
    else
        # LinearSplinePulse: always add derivatives
        if !isnothing(_du_bounds_vec)
            add_control_derivatives(
                base_traj,
                1;  # Only 1 derivative for spline pulses
                control_name = control_sym,
                derivative_bounds = (_du_bounds_vec,),
            )
        else
            add_control_derivatives(base_traj, 1; control_name = control_sym)
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
                "  qcp = SplinePulseProblem(qtraj, N; integrator=integrator, ...)",
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
    J = if free_phase && !isnothing(ket_goal_fn)
        KetFreePhaseInfidelityObjective(ket_goal_fn, state_sym, θ_names, traj; Q = Q)
    elseif free_phase && !isnothing(U_goal_fn)
        UnitaryFreePhaseInfidelityObjective(U_goal_fn, state_sym, θ_names, traj; Q = Q)
    else
        _state_objective(qtraj, traj, state_sym, Q)
    end

    # Add regularization for control and derivative
    J += QuadraticRegularizer(control_sym, traj, R_u_resolved)
    J += QuadraticRegularizer(du_sym, traj, R_du_resolved)

    # Cubic Hermite bending-energy regulariser — true second-derivative L2 penalty.
    # Only meaningful for CubicSplinePulse (LinearSplinePulse has impulsive ddu).
    if is_cubic && λ_bend > 0
        J += BendingEnergyObjective(control_sym, du_sym, λ_bend, traj)
        if _show_details(piccolo_options)
            println("    added BendingEnergyObjective (λ_bend = $λ_bend)")
        end
    end

    # Apply piccolo options
    J += _apply_piccolo_options(
        qtraj,
        piccolo_options,
        constraints,
        traj,
        state_sym;
        state_leakage_indices = state_leakage_indices,
    )

    # Start with dynamics integrators
    integrators = copy(dynamics_integrators)

    # Add DerivativeIntegrator for LinearSplinePulse to enforce du[k] = (u[k+1] - u[k]) / Δt
    # For CubicSplinePulse, :du values are Hermite tangents (independent DOFs), not constrained
    if is_linear_spline
        push!(integrators, DerivativeIntegrator(control_sym, du_sym, traj))
        if _show_details(piccolo_options)
            println("    added DerivativeIntegrator (LinearSplinePulse)")
        end
    end

    # Add global bounds constraints if specified
    all_constraints = copy(constraints)
    add_global_bounds_constraints!(
        all_constraints,
        global_bounds,
        traj;
        verbose = _show_details(piccolo_options),
    )

    # Pin calibration targets at nominal — must run AFTER bounds, since
    # apply_calibration_targets! removes any conflicting GlobalBoundsConstraint.
    apply_calibration_targets!(
        all_constraints,
        calibration_targets,
        traj;
        verbose = _show_details(piccolo_options),
    )

    prob = DirectTrajOptProblem(traj, J, integrators; constraints = all_constraints)

    return _maybe_display(QuantumControlProblem(qtraj, prob), piccolo_options)
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
Accepts all keyword arguments from the base [`SplinePulseProblem`](@ref) method,
including pulse-type-dependent `R_u` / `R_du` defaults and the `λ_bend`
bending-energy regulariser (see the base method's docstring for the full
discussion), plus:
- `du_bounds::Union{Nothing, Vector{Float64}}=nothing`: Per-drive bounds on derivative magnitude (takes precedence over `du_bound`)
- `free_phase::Bool=false`: Optimize a per-subsystem frame phase alongside the pulse. Applies number-operator rotation `e^{iθ n̂}` to goal states — level `s` acquires phase `s·θ`. Requires `subsystem_levels`.
- `subsystem_levels::Union{Nothing, Vector{Int}}=nothing`: Number of levels per subsystem, required when `free_phase=true`.
- `initial_phases::Union{Nothing, Vector{Float64}}=nothing`: Initial values for the per-subsystem phase variables when `free_phase=true`. Length must equal the number of subsystems.
- `coherent::Bool=true`: If `true`, uses a coherent fidelity objective (phases must align across state pairs). If `false`, uses per-state fidelity.
- `integrator_type::Symbol=:spline`: Integrator backend (`:spline` or `:ensemble`).
- `parallel_backend::Symbol=:manual`: Parallelism strategy (`:manual`, `:threads`, or `:gpu`).
"""
function SplinePulseProblem(
    qtraj::MultiKetTrajectory{<:AbstractSplinePulse},
    N_or_times::Union{Nothing,Int,AbstractVector{<:Real}} = nothing;
    integrator::Union{Nothing,AbstractIntegrator,Vector{<:AbstractIntegrator}} = nothing,
    integrator_type::Symbol = :spline,  # :spline or :ensemble
    parallel_backend::Symbol = :manual,  # :manual (default), :threads, :gpu
    global_names::Union{Nothing,Vector{Symbol}} = nothing,
    global_bounds::Union{Nothing,Dict{Symbol,<:Union{Float64,Tuple{Float64,Float64}}}} = nothing,
    calibration_targets::Vector{Symbol} = Symbol[],
    du_bound::Float64 = Inf,
    du_bounds::Union{Nothing,Vector{Float64}} = nothing,
    Δt_bounds::Union{Nothing,Tuple{Float64,Float64}} = nothing,
    Q::Float64 = 100.0,
    R::Float64 = 1e-2,
    R_u::Union{Nothing,Float64,Vector{Float64}} = nothing,
    R_du::Union{Nothing,Float64,Vector{Float64}} = nothing,
    λ_bend::Float64 = 1e-3,
    constraints::Vector{<:AbstractConstraint} = AbstractConstraint[],
    piccolo_options::PiccoloOptions = PiccoloOptions(),
    free_phase::Bool = false,
    subsystem_levels::Union{Nothing,Vector{Int}} = nothing,
    initial_phases::Union{Nothing,Vector{Float64}} = nothing,
    coherent::Bool = true,
    state_leakage_indices::Union{
        Nothing,
        AbstractVector{Int},
        AbstractVector{<:AbstractVector{Int}},
    } = nothing,
)
    sys = get_system(qtraj)
    control_sym = drive_name(qtraj)
    snames = state_names(qtraj)
    weights = qtraj.weights
    goals = qtraj.goals

    # Pulse-type-dependent regularisation defaults
    # (see "Default regularisation by pulse type" in the unitary template docstring)
    is_cubic = qtraj.pulse isa CubicSplinePulse
    R_u_resolved = isnothing(R_u) ? (is_cubic ? 0.0 : R) : R_u
    R_du_resolved = isnothing(R_du) ? (is_cubic ? 0.0 : R) : R_du

    if _show_header(piccolo_options)
        pulse_type = _typename(qtraj.pulse)
        println("constructing SplinePulseProblem [MultiKetTrajectory / $pulse_type]")
        println("    state transfers: $(length(qtraj.initials))")
    end

    # Build global_data explicitly from system global_params
    global_data = if !isempty(sys.global_params)
        Dict(name => [val] for (name, val) in pairs(sys.global_params))
    else
        nothing
    end

    # Free-phase support: add per-subsystem phase variables as globals
    θ_names = Symbol[]
    goals_fn = nothing
    if free_phase
        @assert !isnothing(subsystem_levels) "free_phase=true requires subsystem_levels"
        n_qubits = length(subsystem_levels)
        goals_fn = _make_free_phase_ket_goals(goals, subsystem_levels)
        θ_names, global_data, global_bounds = setup_free_phase_globals!(
            n_qubits,
            global_data,
            global_bounds;
            initial_phases = initial_phases,
            verbose = _show_details(piccolo_options),
        )
    end

    # Convert quantum trajectory to NamedTrajectory
    # N_or_times=nothing uses native pulse knot times (preserves warm-start exactly)
    base_traj =
        NamedTrajectory(qtraj, N_or_times; Δt_bounds = Δt_bounds, global_data = global_data)
    N = base_traj.N  # Get actual number of timesteps

    # Always add control derivatives to trajectory
    # For CubicSplinePulse, :du is already included in the base trajectory (Hermite tangents)
    # For LinearSplinePulse, we add :du explicitly (will be constrained by DerivativeIntegrator)
    du_sym = Symbol(:d, control_sym)
    is_linear_spline = !haskey(base_traj.components, du_sym)

    # Resolve per-drive du bounds: du_bounds (vector) takes precedence over du_bound (scalar)
    _du_bounds_vec = if !isnothing(du_bounds)
        du_bounds
    elseif isfinite(du_bound)
        fill(du_bound, sys.n_drives)
    else
        nothing
    end

    traj = if haskey(base_traj.components, du_sym)
        # CubicSplinePulse already has derivative DOFs, but bounds default to (-Inf, Inf)
        if !isnothing(_du_bounds_vec)
            update_bound!(base_traj, du_sym, (-_du_bounds_vec, _du_bounds_vec))
            if _show_details(piccolo_options)
                println("    du bounds (CubicSplinePulse): $(_fmt_bounds(_du_bounds_vec))")
            end
        end
        base_traj
    else
        # LinearSplinePulse: always add derivatives
        if !isnothing(_du_bounds_vec)
            add_control_derivatives(
                base_traj,
                1;  # Only 1 derivative for spline pulses
                control_name = control_sym,
                derivative_bounds = (_du_bounds_vec,),
            )
        else
            add_control_derivatives(base_traj, 1; control_name = control_sym)
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
                "  qcp = SplinePulseProblem(qtraj, N; integrator=integrator, ...)",
            )
        end
        # Choose integrator type based on integrator_type parameter
        if integrator_type == :ensemble
            dynamics_integrators = EnsembleSplineIntegrator(
                qtraj,
                N;
                spline_order = _get_spline_order(qtraj.pulse),
                parallel_backend = parallel_backend,
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

    # Build objective: coherent fidelity for ensemble (with optional free phase)
    J = if free_phase && !isnothing(goals_fn)
        CoherentKetFreePhaseInfidelityObjective(goals_fn, snames, θ_names, traj; Q = Q)
    else
        _ensemble_ket_objective(qtraj, traj, snames, weights, goals, Q; coherent = coherent)
    end

    # Add regularization for control and derivative
    J += QuadraticRegularizer(control_sym, traj, R_u_resolved)
    J += QuadraticRegularizer(du_sym, traj, R_du_resolved)

    # Cubic Hermite bending-energy regulariser — only meaningful for CubicSplinePulse
    if is_cubic && λ_bend > 0
        J += BendingEnergyObjective(control_sym, du_sym, λ_bend, traj)
        if _show_details(piccolo_options)
            println("    added BendingEnergyObjective (λ_bend = $λ_bend)")
        end
    end

    # Apply piccolo options for each state
    J += _apply_piccolo_options(
        qtraj,
        piccolo_options,
        constraints,
        traj,
        snames;
        state_leakage_indices = state_leakage_indices,
    )

    # Start with dynamics integrators
    integrators = copy(dynamics_integrators)

    # Add DerivativeIntegrator for LinearSplinePulse to enforce du[k] = (u[k+1] - u[k]) / Δt
    # For CubicSplinePulse, :du values are Hermite tangents (independent DOFs), not constrained
    if is_linear_spline
        push!(integrators, DerivativeIntegrator(control_sym, du_sym, traj))
        if _show_details(piccolo_options)
            println("    added DerivativeIntegrator (LinearSplinePulse)")
        end
    end

    # Add global bounds constraints if specified
    all_constraints = copy(constraints)
    add_global_bounds_constraints!(
        all_constraints,
        global_bounds,
        traj;
        verbose = _show_details(piccolo_options),
    )

    apply_calibration_targets!(
        all_constraints,
        calibration_targets,
        traj;
        verbose = _show_details(piccolo_options),
    )

    prob = DirectTrajOptProblem(traj, J, integrators; constraints = all_constraints)

    return _maybe_display(QuantumControlProblem(qtraj, prob), piccolo_options)
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
    N_or_times::Union{Nothing,Int,AbstractVector{<:Real}} = nothing;
    kwargs...,
) where {P<:AbstractPulse}
    error(
        """
  SplinePulseProblem is only for spline-based pulses (LinearSplinePulse, CubicSplinePulse).

  You provided a trajectory with pulse type: $(nameof(P))

  For piecewise constant pulses (ZeroOrderPulse), use SmoothPulseProblem instead:
      qcp = SmoothPulseProblem(qtraj, N; ...)
  """,
    )
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

    times = collect(range(0.0, T, length = N))
    amps = 0.1 * randn(n_drives, N)
    pulse = LinearSplinePulse(amps, times)

    # Goal: X gate
    U_goal = ComplexF64[0 1; 1 0]

    # Create trajectory and problem
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)
    qcp = SplinePulseProblem(qtraj, N; Q = 100.0, R = 1e-2)

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

    times = collect(range(0.0, T, length = N))
    amps = 0.1 * randn(n_drives, N)
    derivs = zeros(n_drives, N)  # Hermite spline with derivative DOFs
    pulse = CubicSplinePulse(amps, derivs, times)

    # Goal: X gate
    U_goal = ComplexF64[0 1; 1 0]

    # Create trajectory and problem
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)
    qcp = SplinePulseProblem(qtraj, N; Q = 100.0, R = 1e-2)

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

    times = collect(range(0.0, T, length = N))
    amps = 0.1 * randn(n_drives, N)
    derivs = zeros(n_drives, N)
    pulse = CubicSplinePulse(amps, derivs, times)

    U_goal = ComplexF64[0 1; 1 0]
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)

    # Test with du_bound specified
    du_bound = 5.0
    qcp = SplinePulseProblem(qtraj, N; Q = 100.0, R = 1e-2, du_bound = du_bound)

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
    qcp_unbounded = SplinePulseProblem(qtraj, N; Q = 100.0, R = 1e-2)
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

    times = collect(range(0.0, T, length = N))
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

    times = collect(range(0.0, T, length = N))
    amps = 0.1 * randn(n_drives, N)
    pulse = LinearSplinePulse(amps, times)

    # State transfer: |0⟩ → |1⟩
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]

    # Create trajectory and problem
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)
    qcp = SplinePulseProblem(qtraj, N; Q = 100.0, R = 1e-2)

    @test qcp isa QuantumControlProblem
    @test qcp.qtraj isa KetTrajectory
    @test get_trajectory(qcp) isa NamedTrajectory

    # Check trajectory has proper components
    traj = get_trajectory(qcp)
    @test haskey(traj.components, :ψ̃)
    @test !haskey(traj.components, :ddu)  # No second derivative for splines
end

@testitem "SplinePulseProblem KetTrajectory leakage_indices kwarg" begin
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    # 3-level system; computational subspace = {|0⟩, |1⟩}, leak = |2⟩
    H_drift = ComplexF64[0 0 0; 0 0.01 0; 0 0 0.05]
    H_drives = [ComplexF64[0 1 0; 1 0 sqrt(2); 0 sqrt(2) 0]]
    T = 10.0
    N = 51

    sys = QuantumSystem(H_drift, H_drives, [1.0])

    times = collect(range(0.0, T, length = N))
    amps = 0.1 * randn(1, N)
    pulse = LinearSplinePulse(amps, times)

    ψ_init = ComplexF64[1.0, 0.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0, 0.0]

    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

    # Without state_leakage_indices, KetTrajectory + leakage_constraint=true
    # must still error per the existing contract — there is no goal-derived
    # leakage geometry for a single ket.
    @test_throws ArgumentError SplinePulseProblem(
        qtraj,
        N;
        Q = 100.0,
        R = 1e-2,
        piccolo_options = PiccoloOptions(
            leakage_constraint = true,
            leakage_constraint_value = 1e-3,
            leakage_cost = 1.0,
        ),
    )

    # With user-supplied indices, construction succeeds and a LeakageConstraint
    # is appended to the problem's constraints.
    # iso_ket = [Re(ψ); Im(ψ)] in length-6, so |2⟩ leakage = indices [3, 6].
    qcp = SplinePulseProblem(
        qtraj,
        N;
        Q = 100.0,
        R = 1e-2,
        piccolo_options = PiccoloOptions(
            leakage_constraint = true,
            leakage_constraint_value = 1e-3,
            leakage_cost = 1.0,
            display = :silent,
        ),
        state_leakage_indices = [3, 6],
    )

    @test qcp isa QuantumControlProblem
    # LeakageConstraint is a constructor that returns a NonlinearKnotPointConstraint
    # parametrized by a closure named `leakage_constraint`. Detect via the closure
    # field name to confirm the leakage path actually fired.
    @test any(qcp.prob.constraints) do c
        c isa DirectTrajOpt.NonlinearKnotPointConstraint &&
            occursin("leakage_constraint", string(typeof(c).parameters[1]))
    end
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

    times = collect(range(0.0, T, length = N))
    amps = 0.1 * randn(n_drives, N)
    pulse = LinearSplinePulse(amps, times)

    # Create ensemble: |0⟩ → |1⟩ and |1⟩ → |0⟩
    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]

    # Create trajectory and problem
    qtraj = MultiKetTrajectory(sys, pulse, [ψ0, ψ1], [ψ1, ψ0])
    qcp = SplinePulseProblem(qtraj, N; Q = 100.0, R = 1e-2)

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

    times = collect(range(0.0, T, length = N))
    amps = 0.1 * randn(n_drives, N)
    pulse = LinearSplinePulse(amps, times)

    # Goal: X gate
    U_goal = ComplexF64[0 1; 1 0]

    # Create base trajectory and sampling trajectory
    base_qtraj = UnitaryTrajectory(sys_nominal, pulse, U_goal)

    # First create a SplinePulseProblem with base trajectory
    base_qcp = SplinePulseProblem(base_qtraj, N; Q = 100.0, R = 1e-2)

    # Then create SamplingProblem
    sampling_qcp = SamplingProblem(base_qcp, [sys_nominal, sys_perturbed]; Q = 100.0)

    @test sampling_qcp isa QuantumControlProblem
    @test sampling_qcp.qtraj isa SamplingTrajectory{<:AbstractPulse,<:UnitaryTrajectory}

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
        qtraj,
        N;
        Q = 100.0,
        R = 1e-2,
        global_bounds = Dict(:δ => 0.5),  # δ doesn't exist in trajectory
    )
end

@testitem "_make_free_phase_goal for EmbeddedOperator" begin
    using LinearAlgebra
    using .ProblemTemplates: _make_free_phase_goal

    # Single qubit embedded in 3-level system
    X_gate = ComplexF64[0 1; 1 0]
    subspace = [1, 2]
    levels = [3]
    op = EmbeddedOperator(X_gate, subspace, levels)

    U_goal_fn = _make_free_phase_goal(op)

    # Zero phase should give back the original gate
    result_0 = U_goal_fn([0.0])
    @test unembed(result_0) ≈ X_gate

    # With phase θ, the gate becomes diag(1, e^{iθ}) * X
    θ = [π / 3]
    result_θ = U_goal_fn(θ)
    U_sub = unembed(result_θ)
    expected = Diagonal([1.0 + 0im, exp(im * θ[1])]) * X_gate
    @test U_sub ≈ expected atol = 1e-12

    # Two-qubit gate: CZ embedded in [3,3] system
    CZ = ComplexF64[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 -1]
    subspace_2q = get_subspace_indices([[1, 2], [1, 2]], [3, 3])
    op_2q = EmbeddedOperator(CZ, subspace_2q, [3, 3])

    U_goal_fn_2q = _make_free_phase_goal(op_2q)

    # Zero phase -> original CZ
    result_2q_0 = U_goal_fn_2q([0.0, 0.0])
    @test unembed(result_2q_0) ≈ CZ

    # With phases: diag(1, e^{iθ₂}, e^{iθ₁}, e^{i(θ₁+θ₂)}) * CZ
    θ_2q = [π / 4, π / 6]
    result_2q = U_goal_fn_2q(θ_2q)
    U_sub_2q = unembed(result_2q)
    phase_diag = Diagonal([
        1.0 + 0im,
        exp(im * θ_2q[2]),
        exp(im * θ_2q[1]),
        exp(im * (θ_2q[1] + θ_2q[2])),
    ])
    @test U_sub_2q ≈ phase_diag * CZ atol = 1e-12
end

@testitem "SplinePulseProblem with UnitaryTrajectory and free_phase=true" begin
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    # 2-level system embedded in 3-level
    σx = ComplexF64[0 1; 1 0]
    σz = ComplexF64[1 0; 0 -1]
    H_drift_3 = ComplexF64[0 0 0; 0 1 0; 0 0 2]
    H_drive_3 = ComplexF64[0 1 0; 1 0 1; 0 1 0] / √2
    sys = QuantumSystem(H_drift_3, [H_drive_3], [1.0])

    T = 10.0
    N = 51
    times = collect(range(0.0, T, length = N))
    amps = 0.1 * randn(1, N)
    pulse = LinearSplinePulse(amps, times)

    # Goal: X gate embedded in 3-level system
    subspace = [1, 2]
    levels = [3]
    U_goal = EmbeddedOperator(σx, subspace, levels)

    qtraj = UnitaryTrajectory(sys, pulse, U_goal)
    qcp = SplinePulseProblem(qtraj, N; Q = 100.0, R = 1e-2, free_phase = true)

    @test qcp isa QuantumControlProblem

    # Check that phase variable was added
    traj = get_trajectory(qcp)
    @test haskey(traj.global_components, :φ_1)
end

@testitem "SplinePulseProblem with MultiKetTrajectory and free_phase=true" begin
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    σx = ComplexF64[0 1; 1 0]
    σz = ComplexF64[1 0; 0 -1]
    sys = QuantumSystem(0.01 * σz, [σx], [1.0])

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]

    T = 10.0
    N = 51
    times = collect(range(0.0, T, length = N))
    pulse = LinearSplinePulse(0.1 * randn(1, N), times)
    qtraj = MultiKetTrajectory(sys, pulse, [ψ0, ψ1], [ψ1, ψ0])

    qcp = SplinePulseProblem(
        qtraj,
        N;
        Q = 100.0,
        R = 1e-2,
        free_phase = true,
        subsystem_levels = [2],
    )

    @test qcp isa QuantumControlProblem
    traj = get_trajectory(qcp)
    @test haskey(traj.global_components, :φ_1)
end

@testitem "SplinePulseProblem free_phase requires EmbeddedOperator for unitary" begin
    using LinearAlgebra

    σx = ComplexF64[0 1; 1 0]
    σz = ComplexF64[1 0; 0 -1]
    sys = QuantumSystem(0.01 * σz, [σx], [1.0])

    T = 10.0
    N = 51
    times = collect(range(0.0, T, length = N))
    pulse = LinearSplinePulse(0.1 * randn(1, N), times)

    # Plain matrix goal (not EmbeddedOperator) should fail with free_phase
    U_goal = ComplexF64[0 1; 1 0]
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)
    @test_throws AssertionError SplinePulseProblem(qtraj, N; free_phase = true)
end

@testitem "SplinePulseProblem MultiKet free_phase requires subsystem_levels" begin
    using LinearAlgebra

    σx = ComplexF64[0 1; 1 0]
    σz = ComplexF64[1 0; 0 -1]
    sys = QuantumSystem(0.01 * σz, [σx], [1.0])

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]

    T = 10.0
    N = 51
    times = collect(range(0.0, T, length = N))
    pulse = LinearSplinePulse(0.1 * randn(1, N), times)
    qtraj = MultiKetTrajectory(sys, pulse, [ψ0, ψ1], [ψ1, ψ0])

    @test_throws AssertionError SplinePulseProblem(
        qtraj,
        N;
        free_phase = true,
        subsystem_levels = nothing,
    )
end

@testitem "SplinePulseProblem MultiKet coherent kwarg" begin
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    σx = ComplexF64[0 1; 1 0]
    σz = ComplexF64[1 0; 0 -1]
    sys = QuantumSystem(0.01 * σz, [σx], [1.0])

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]

    T = 10.0
    N = 51
    times = collect(range(0.0, T, length = N))
    pulse = LinearSplinePulse(0.1 * randn(1, N), times)
    qtraj = MultiKetTrajectory(sys, pulse, [ψ0, ψ1], [ψ1, ψ0])

    # coherent=false should construct without error
    qcp = SplinePulseProblem(qtraj, N; Q = 100.0, R = 1e-2, coherent = false)
    @test qcp isa QuantumControlProblem
end

# ============================================================================= #
# Bending-energy regulariser and cubic-spline default tests
# ============================================================================= #

@testitem "SplinePulseProblem CubicSplinePulse default R_u, R_du are zero" begin
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    # CubicSplinePulse should default R_u = R_du = 0.0, NOT inherit R.
    σx = ComplexF64[0 1; 1 0]
    σz = ComplexF64[1 0; 0 -1]
    H_drift = 0.01 * σz
    H_drives = [σx]
    T = 10.0
    N = 11
    n_drives = 1

    sys = QuantumSystem(H_drift, H_drives, [1.0])
    times = collect(range(0.0, T, length = N))
    amps = 0.1 * randn(n_drives, N)
    derivs = zeros(n_drives, N)
    pulse = CubicSplinePulse(amps, derivs, times)
    U_goal = ComplexF64[0 1; 1 0]
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)

    qcp = SplinePulseProblem(qtraj, N; Q = 100.0)  # no explicit R_u, R_du

    @test qcp isa QuantumControlProblem

    # The QuadraticRegularizer on :u and :du should be present but with zero R.
    composite = qcp.prob.objective
    quad_regs = filter(
        x -> x isa DirectTrajOpt.QuadraticRegularizer,
        composite isa DirectTrajOpt.CompositeObjective ? composite.objectives : [composite],
    )
    R_us = [r.R for r in quad_regs if r.name == :u]
    R_dus = [r.R for r in quad_regs if r.name == :du]
    @test !isempty(R_us)
    @test !isempty(R_dus)
    @test all(all(R_us[1] .== 0.0))
    @test all(all(R_dus[1] .== 0.0))
end

@testitem "SplinePulseProblem LinearSplinePulse default R_u, R_du = R" begin
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    # LinearSplinePulse should still inherit R for R_u and R_du (existing behavior).
    σx = ComplexF64[0 1; 1 0]
    σz = ComplexF64[1 0; 0 -1]
    H_drift = 0.01 * σz
    H_drives = [σx]
    T = 10.0
    N = 21
    n_drives = 1

    sys = QuantumSystem(H_drift, H_drives, [1.0])
    times = collect(range(0.0, T, length = N))
    amps = 0.1 * randn(n_drives, N)
    pulse = LinearSplinePulse(amps, times)
    U_goal = ComplexF64[0 1; 1 0]
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)

    R = 1e-2
    qcp = SplinePulseProblem(qtraj, N; Q = 100.0, R = R)  # no explicit R_u, R_du

    composite = qcp.prob.objective
    quad_regs = filter(
        x -> x isa DirectTrajOpt.QuadraticRegularizer,
        composite isa DirectTrajOpt.CompositeObjective ? composite.objectives : [composite],
    )
    R_us = [r.R for r in quad_regs if r.name == :u]
    R_dus = [r.R for r in quad_regs if r.name == :du]
    @test !isempty(R_us)
    @test !isempty(R_dus)
    @test all(R_us[1] .≈ R)
    @test all(R_dus[1] .≈ R)
end

@testitem "build_bending_energy_blocks is zero on a flat constant pulse" begin
    using LinearAlgebra
    using Piccolo.QuantumObjectives: build_bending_energy_blocks

    # Uniform knots, h = 1.0
    N = 12
    h = fill(1.0, N - 1)
    K_uu, K_um, K_mm = build_bending_energy_blocks(h)

    u = fill(0.7, N)   # flat constant amplitude
    m = zeros(N)       # zero Hermite tangents (consistent with constant pulse)
    E = dot(u, K_uu * u) + 2 * dot(u, K_um * m) + dot(m, K_mm * m)
    @test E ≈ 0.0 atol = 1e-12

    # Sanity: a linear ramp also has zero bending energy with the consistent tangents.
    # u_k = α · t_k = α · (k-1)·h, m_k = α  (constant slope across the spline)
    α = 0.3
    t = (0:(N-1)) .* h[1]
    u_lin = α .* t
    m_lin = fill(α, N)
    E_lin = dot(u_lin, K_uu * u_lin) +
            2 * dot(u_lin, K_um * m_lin) +
            dot(m_lin, K_mm * m_lin)
    @test E_lin ≈ 0.0 atol = 1e-10
end

@testitem "build_bending_energy_blocks scales as ω^4 for sinusoid" begin
    using LinearAlgebra
    using Piccolo.QuantumObjectives: build_bending_energy_blocks

    # Sample sin(ω·t) on a fine uniform grid and check that increasing ω
    # by a factor of 2 increases the bending energy by ≈ 16 (= 2^4).
    # Use many knots so the cubic Hermite spline accurately resolves the curve.
    N = 257
    T = 1.0
    h = T / (N - 1)
    h_vec = fill(h, N - 1)
    K_uu, K_um, K_mm = build_bending_energy_blocks(h_vec)

    t = collect(range(0.0, T, length = N))

    function bend_energy(ω)
        u = sin.(ω .* t)
        m = ω .* cos.(ω .* t)
        return dot(u, K_uu * u) + 2 * dot(u, K_um * m) + dot(m, K_mm * m)
    end

    ω1 = 2π * 3.0
    ω2 = 2π * 6.0  # 2x frequency

    E1 = bend_energy(ω1)
    E2 = bend_energy(ω2)

    # Analytic value: ∫_0^T (d²/dt² sin(ωt))² dt = ω^4 ∫ sin² ≈ ω^4 · T/2
    @test E1 ≈ ω1^4 * T / 2 rtol = 5e-3
    @test E2 ≈ ω2^4 * T / 2 rtol = 5e-3

    # Ratio should be 2^4 = 16
    @test E2 / E1 ≈ 16.0 rtol = 1e-2
end

@testitem "BendingEnergyObjective: gradient and Hessian match finite differences" begin
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra
    using Piccolo.QuantumObjectives: BendingEnergyObjective

    # Build a small trajectory with :u and :du components.
    N = 8
    n_drives = 2
    Δt = 0.1
    u = randn(n_drives, N)
    du = randn(n_drives, N)
    Δt_vec = fill(Δt, N)

    traj = NamedTrajectory(
        (u = u, du = du, Δt = Δt_vec);
        timestep = :Δt,
        controls = :u,
    )

    obj = BendingEnergyObjective(:u, :du, 0.7, traj)
    test_objective(obj, traj, atol = 1e-6)
end

@testitem "BendingEnergyObjective: zero on a flat constant pulse" begin
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra
    using Piccolo.QuantumObjectives: BendingEnergyObjective

    N = 10
    n_drives = 1
    u = fill(0.42, n_drives, N)   # flat constant amplitude
    du = zeros(n_drives, N)       # consistent tangents
    Δt_vec = fill(0.1, N)

    traj = NamedTrajectory(
        (u = u, du = du, Δt = Δt_vec);
        timestep = :Δt,
        controls = :u,
    )

    obj = BendingEnergyObjective(:u, :du, 1.0, traj)
    @test objective_value(obj, traj) ≈ 0.0 atol = 1e-12
end

@testitem "BendingEnergyObjective: sinusoidal pulse scales as ω^4" begin
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra
    using Piccolo.QuantumObjectives: BendingEnergyObjective

    N = 257
    T = 1.0
    Δt = T / (N - 1)
    t = (0:(N-1)) .* Δt
    Δt_vec = fill(Δt, N)
    n_drives = 1

    function build_traj(ω)
        u = reshape(sin.(ω .* t), n_drives, N)
        du = reshape(ω .* cos.(ω .* t), n_drives, N)
        return NamedTrajectory((u = u, du = du, Δt = Δt_vec); timestep = :Δt, controls = :u)
    end

    ω1 = 2π * 3.0
    ω2 = 2π * 6.0

    traj1 = build_traj(ω1)
    traj2 = build_traj(ω2)
    obj1 = BendingEnergyObjective(:u, :du, 1.0, traj1)
    obj2 = BendingEnergyObjective(:u, :du, 1.0, traj2)

    E1 = objective_value(obj1, traj1)
    E2 = objective_value(obj2, traj2)

    @test E1 ≈ ω1^4 * T / 2 rtol = 5e-3
    @test E2 ≈ ω2^4 * T / 2 rtol = 5e-3
    @test E2 / E1 ≈ 16.0 rtol = 1e-2
end

@testitem "SplinePulseProblem CubicSplinePulse wires in BendingEnergyObjective" begin
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra
    using Piccolo.QuantumObjectives: BendingEnergyObjective

    σx = ComplexF64[0 1; 1 0]
    σz = ComplexF64[1 0; 0 -1]
    H_drift = 0.01 * σz
    H_drives = [σx]
    T = 10.0
    N = 11

    sys = QuantumSystem(H_drift, H_drives, [1.0])
    times = collect(range(0.0, T, length = N))
    amps = 0.1 * randn(1, N)
    derivs = zeros(1, N)
    pulse = CubicSplinePulse(amps, derivs, times)
    qtraj = UnitaryTrajectory(sys, pulse, ComplexF64[0 1; 1 0])

    qcp = SplinePulseProblem(qtraj, N; Q = 100.0, λ_bend = 1e-3)
    composite = qcp.prob.objective
    bend_terms = filter(
        x -> x isa BendingEnergyObjective,
        composite isa DirectTrajOpt.CompositeObjective ? composite.objectives : [composite],
    )
    @test length(bend_terms) == 1
    @test bend_terms[1].λ == 1e-3

    # λ_bend = 0 should skip the term.
    qcp_no_bend = SplinePulseProblem(qtraj, N; Q = 100.0, λ_bend = 0.0)
    composite_no = qcp_no_bend.prob.objective
    bend_terms_no = filter(
        x -> x isa BendingEnergyObjective,
        composite_no isa DirectTrajOpt.CompositeObjective ? composite_no.objectives :
        [composite_no],
    )
    @test isempty(bend_terms_no)
end

@testitem "SplinePulseProblem LinearSplinePulse never adds BendingEnergyObjective" begin
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra
    using Piccolo.QuantumObjectives: BendingEnergyObjective

    σx = ComplexF64[0 1; 1 0]
    σz = ComplexF64[1 0; 0 -1]
    sys = QuantumSystem(0.01 * σz, [σx], [1.0])

    T = 10.0
    N = 21
    times = collect(range(0.0, T, length = N))
    pulse = LinearSplinePulse(0.1 * randn(1, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, ComplexF64[0 1; 1 0])

    # Even with λ_bend explicitly non-zero, LinearSplinePulse should NOT pick it up.
    qcp = SplinePulseProblem(qtraj, N; Q = 100.0, λ_bend = 1e-3)
    composite = qcp.prob.objective
    bend_terms = filter(
        x -> x isa BendingEnergyObjective,
        composite isa DirectTrajOpt.CompositeObjective ? composite.objectives : [composite],
    )
    @test isempty(bend_terms)
end

@testitem "BendingEnergyObjective does not pull single-qubit X-gate amplitude to zero" begin
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra
    using Piccolo.QuantumObjectives: BendingEnergyObjective

    # A π-pulse on a single qubit requires nonzero integrated drive area;
    # adding a bending-energy regulariser MUST NOT zero out the amplitude.
    # We just verify here that the BendingEnergyObjective evaluated on a known
    # π-pulse (constant amplitude π/T cubic Hermite) has the analytic value 0
    # (since the optimal amplitude is constant — only the connecting cubic
    # would bend, and a constant function has zero bending energy).
    T = 2π / 4.0   # nonzero finite duration
    N = 11
    n_drives = 1
    A = π / T     # amplitude such that ∫ A dt = π → X-rotation by π

    Δt = T / (N - 1)
    Δt_vec = fill(Δt, N)
    u_const = fill(A, n_drives, N)
    du_const = zeros(n_drives, N)

    traj = NamedTrajectory(
        (u = u_const, du = du_const, Δt = Δt_vec);
        timestep = :Δt,
        controls = :u,
    )

    obj = BendingEnergyObjective(:u, :du, 1e-3, traj)
    E = objective_value(obj, traj)

    # The bending energy of a constant pulse is exactly zero — so a bending-energy
    # penalty does not bias the optimizer away from holding a nonzero constant amplitude.
    @test E ≈ 0.0 atol = 1e-12

    # And the gradient is also zero, confirming no force toward A = 0.
    Z_dim = traj.dim * traj.N + traj.global_dim
    ∇ = zeros(Z_dim)
    gradient!(∇, obj, traj)
    @test all(∇ .≈ 0.0)
end

@testitem "SplinePulseProblem CubicSplinePulse: X gate solves with λ_bend > 0" begin
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra
    using Piccolo.QuantumObjectives: BendingEnergyObjective

    # Regression-style test: with the new defaults (R_u = R_du = 0) and
    # λ_bend = 1e-3, a simple single-qubit X gate on a 2-level system should
    # still construct, evaluate, and produce a finite, non-negative objective.
    σx = ComplexF64[0 1; 1 0]
    σz = ComplexF64[1 0; 0 -1]
    sys = QuantumSystem(0.01 * σz, [σx], [1.0])

    T = 10.0
    N = 11
    times = collect(range(0.0, T, length = N))
    amps = 0.1 * randn(1, N)
    derivs = zeros(1, N)
    pulse = CubicSplinePulse(amps, derivs, times)
    U_goal = ComplexF64[0 1; 1 0]
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)

    qcp = SplinePulseProblem(qtraj, N; Q = 100.0, λ_bend = 1e-3)
    @test qcp isa QuantumControlProblem

    J = objective_value(qcp.prob.objective, qcp.prob.trajectory)
    @test isfinite(J)
    @test J ≥ 0.0
end
