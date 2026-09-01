# ----------------------------------------------------------------------------- #
# Template declaration
# ----------------------------------------------------------------------------- #

@doc raw"""
    SplinePulseParams <: AbstractTemplateParams

The typed keyword surface of [`SplinePulseProblem`](@ref). `du_bounds` (per-drive
Hermite-tangent bounds), `integrator_type` and `parallel_backend` live **only**
here; there is no `R_ddu`/`ddu_bound` because spline knot tangents are independent
degrees of freedom rather than discrete derivatives.

`R_u`/`R_du` default to `nothing`, not to `R`: the resolved default is
pulse-type-dependent (0 for cubic Hermite, `R` for linear). The params struct
records the *declared* keyword, so `nothing` here means "let the template resolve
it from the pulse type".
"""
Base.@kwdef struct SplinePulseParams <: AbstractTemplateParams
    Q::Float64 = 100.0
    R::Float64 = 1e-2
    R_u::Union{Nothing,Float64,Vector{Float64}} = nothing
    R_du::Union{Nothing,Float64,Vector{Float64}} = nothing
    R_bend::Union{Nothing,Real,AbstractVector{<:Real}} = nothing
    du_bound::Float64 = Inf
    du_bounds::Union{Nothing,Vector{Float64}} = nothing
    Δt_bounds::Union{Nothing,Tuple{Float64,Float64}} = nothing
    free_phase::Bool = false
    subsystem_levels::Union{Nothing,Vector{Int}} = nothing
    initial_phases::Union{Nothing,Vector{Float64}} = nothing
    coherent::Bool = true
    integrator_type::Symbol = :spline
    spline_interior_bound_constraints::Bool = false
    parallel_backend::Symbol = :manual
    global_names::Union{Nothing,Vector{Symbol}} = nothing
    global_bounds::Union{Nothing,AbstractDict} = nothing
    calibration_targets::Vector{Symbol} = Symbol[]
    state_leakage_indices::Union{
        Nothing,
        AbstractVector{Int},
        AbstractVector{<:AbstractVector{Int}},
    } = nothing
end

@problem_template SplinePulseTemplate begin
    julia_name = SplinePulseProblem
    pulse = AbstractSplinePulse
    trajectories = (UnitaryTrajectory, KetTrajectory, MultiKetTrajectory, DensityTrajectory)
    pulse_kinds = (:linear_spline, :cubic_spline)
    trajectory_kinds = (:unitary, :ket)
    ket_free_phase = true
    params = SplinePulseParams
    passthrough = (:integrator, :constraints, :extra_objectives, :piccolo_options)
    builder = _spline_pulse_problem
    requires_N = false
    hint = """For piecewise constant pulses (ZeroOrderPulse), use SmoothPulseProblem instead:
      qcp = SmoothPulseProblem(qtraj, N; ...)"""
end


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
- `R_bend::Union{Nothing, Real, AbstractVector{<:Real}}=nothing`: Weight on bending-energy regularization ∫u″²dt (`HermiteBendingEnergyRegularizer`). **CubicSplinePulse defaults to `1e-3` (ON)** — a deliberately gentle weight; the fidelity objective stays dominant. Pass `R_bend = 0` to opt out. LinearSplinePulse: bend is undefined on C⁰ families — a nonzero value errors.
- `constraints::Vector{<:AbstractConstraint}=AbstractConstraint[]`: Additional constraints
- `extra_objectives::Vector{<:AbstractObjective}=AbstractObjective[]`: Additional objective terms to compose into the problem's total objective (e.g., `Piccolissimo.HermiteBendingEnergyRegularizer`). Each entry is summed into the objective before constructing the underlying `DirectTrajOptProblem`. Default empty: no extra terms.
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
  versus ≈ 0.98+ at zero regularisation.

Pass `R_u` or `R_du` explicitly to override these defaults.

## Cubic-spline smoothness regularization

Cubic-Hermite-spline smoothness regularization (bending energy, ∫u″²dt) is
**built in and ON by default** for `CubicSplinePulse` via the `R_bend` kwarg
(default `1e-3`; `R_bend = 0` opts out). The closed-form
`HermiteBendingEnergyRegularizer` now lives in Piccolo itself (#309).

The previous long-form — constructing the regularizer against the named
trajectory and passing it via `extra_objectives` — still works and remains
the right tool when you need per-drive weights or a custom weight the kwarg
doesn't express:

```julia
qtraj = UnitaryTrajectory(sys, pulse, U_target)
traj = NamedTrajectory(qtraj, K)
bending_reg = HermiteBendingEnergyRegularizer(traj; R = [0.01, 0.2])

qcp = SplinePulseProblem(qtraj, K;
    Q = 100.0,
    R_bend = 0,                          # don't double-count with the kwarg
    extra_objectives = [bending_reg],
    free_phase = true,
    subsystem_levels = [2, 2],
)
solve!(qcp; max_iter = 500, eval_hessian = true)
```

Any `AbstractObjective` defined against the same trajectory variables can be
injected via `extra_objectives`.

# Returns
- `SplinePulseProblem` (= `QuantumControlProblem{SplinePulseTemplate, QT}`): Wrapper containing trajectory and optimization problem

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
""" SplinePulseProblem

# The construction logic the generated constructor delegates to. Returns the
# untagged problem; `@problem_template`'s constructor stamps on the tag + params.
function _spline_pulse_problem(
    qtraj::AbstractQuantumTrajectory{<:AbstractSplinePulse},
    N_or_times::Union{Nothing,Int,AbstractVector{<:Real}} = nothing;
    integrator::Union{Nothing,AbstractIntegrator,Vector{<:AbstractIntegrator}} = nothing,
    integrator_type::Union{Nothing,Symbol} = nothing,
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
    R_bend::Union{Nothing,Real,AbstractVector{<:Real}} = nothing,  # resolved per pulse type below
    constraints::Vector{<:AbstractConstraint} = AbstractConstraint[],
    extra_objectives::Vector{<:AbstractObjective} = AbstractObjective[],
    piccolo_options::PiccoloOptions = PiccoloOptions(),
    free_phase::Bool = false,
    subsystem_levels::Union{Nothing,Vector{Int}} = nothing,
    initial_phases::Union{Nothing,Vector{Float64}} = nothing,
    state_leakage_indices::Union{
        Nothing,
        AbstractVector{Int},
        AbstractVector{<:AbstractVector{Int}},
    } = nothing,
    spline_interior_bound_constraints::Bool = false,
    n_interior_bound_points::Int = 3,
)
    sys = get_system(qtraj)
    control_sym = drive_name(qtraj)
    state_sym = state_name(qtraj)

    # Pulse-type-dependent regularisation defaults.
    # See the "Default regularisation by pulse type" section of this template's docstring.
    is_cubic = qtraj.pulse isa CubicSplinePulse
    R_u_resolved = isnothing(R_u) ? (is_cubic ? 0.0 : R) : R_u
    R_du_resolved = isnothing(R_du) ? (is_cubic ? 0.0 : R) : R_du

    # Bending energy (#309): default ON (1e-3) for cubic, undefined on C⁰ families.
    if isnothing(R_bend)
        R_bend_resolved = is_cubic ? 1e-3 : nothing
    else
        _bend_nonzero = R_bend isa Real ? !iszero(R_bend) : any(!iszero, R_bend)
        if !is_cubic && _bend_nonzero
            error(
                "R_bend is only defined for CubicSplinePulse (bending energy ∫u″²dt " *
                "requires a C¹ Hermite parameterization). For LinearSplinePulse the " *
                "second derivative is distributional — use R_du/du_bound instead.",
            )
        end
        R_bend_resolved = _bend_nonzero ? R_bend : nothing
    end

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
    N = base_traj.N  # Get actual number of knot points

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
        # Spline pulses must not silently land on PWC dynamics: BilinearIntegrator
        # never reads :du, so a spline pulse would optimize a different waveform
        # than its name promises (issue #275). `integrator_type = :pwc` is the
        # explicit, acknowledged escape hatch.
        isnothing(integrator_type) ||
            integrator_type === :pwc ||
            error(
                "unknown `integrator_type = :$integrator_type`. " *
                "Piccolo ships one backend: `:pwc` (`BilinearIntegrator`).",
            )
        if qtraj.pulse isa CubicSplinePulse
            if integrator_type === :pwc
                @warn "CubicSplinePulse with the PWC backend (`integrator_type = :pwc`): the " *
                      "dynamics treat the drive as piecewise constant and ignore :du — the " *
                      "optimized waveform differs from the cubic spline the pulse object " *
                      "describes. Acknowledged because you asked for it explicitly." maxlog =
                    1
            else
                error(
                    "CubicSplinePulse defaults are not allowed: the default PWC backend " *
                    "silently drops :du (a cubic problem would optimize a piecewise-constant " *
                    "waveform, not a spline — issue #275). Either pass a spline-faithful " *
                    "integrator (Piccolissimo.SplineIntegrator) or explicitly request the PWC " *
                    "backend with `integrator_type = :pwc`.",
                )
            end
        elseif qtraj.pulse isa AbstractSplinePulse
            @warn "SplinePulseProblem default BilinearIntegrator with $(typeof(qtraj.pulse).name.name): " *
                  "Bilinear is PWC and ignores :du. Pass a SplineIntegrator (Piccolissimo) for " *
                  "spline-faithful dynamics, or `integrator_type = :pwc` to acknowledge." maxlog =
                1
        end
        default_int = BilinearIntegrator(qtraj, N)

        if default_int isa AbstractVector
            dynamics_integrators = AbstractIntegrator[default_int...]
        else
            dynamics_integrators = AbstractIntegrator[default_int]
        end
    elseif integrator isa AbstractIntegrator
        # Guard against PWC-vs-spline mismatch (H1)
        if qtraj.pulse isa CubicSplinePulse && integrator isa BilinearIntegrator
            error(
                "CubicSplinePulse with BilinearIntegrator: Bilinear never reads :du (H1). Use SplineIntegrator.",
            )
        end
        dynamics_integrators = AbstractIntegrator[integrator]
    else
        for integ in integrator
            if qtraj.pulse isa CubicSplinePulse && integ isa BilinearIntegrator
                error(
                    "CubicSplinePulse with BilinearIntegrator: Bilinear never reads :du (H1). Use SplineIntegrator.",
                )
            end
        end
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

    # Bending energy (#309): default smoothness term for cubic-Hermite pulses.
    if !isnothing(R_bend_resolved)
        J += HermiteBendingEnergyRegularizer(
            traj;
            R = R_bend_resolved,
            control_name = control_sym,
            derivative_name = du_sym,
            timestep_name = traj.timestep,
        )
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

    # Compose user-supplied extra objective terms (e.g.,
    # Piccolissimo.HermiteBendingEnergyRegularizer). Done after the built-in
    # regularizers so the composite ordering reads: fidelity → built-in regs
    # → piccolo-options terms → user extras.
    for extra_obj in extra_objectives
        J += extra_obj
    end

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

    # Spline interior bounds (H10) — CubicSplinePulse can exceed knot bounds in interior.
    # Slice 3c (#431): CubicSplineBoundConstraint lives in Piccolo
    # (Control.QuantumConstraints.SplineConstraints) — a hard import; the old
    # isdefined(Piccolissimo, ...) soft-dependency sniff is gone.
    if spline_interior_bound_constraints
        if qtraj.pulse isa CubicSplinePulse
            for (drive_idx, (lb, ub)) in enumerate(sys.drive_bounds)
                if isfinite(lb) && isfinite(ub)
                    push!(
                        constraints,
                        CubicSplineBoundConstraint(
                            traj,
                            control_sym,
                            lb,
                            ub;
                            n_interior_points = n_interior_bound_points,
                        ),
                    )
                end
            end
            if _show_details(piccolo_options)
                println(
                    "    added CubicSplineBoundConstraint (n=$(n_interior_bound_points) per segment, H10)",
                )
            end
        else
            @warn "spline_interior_bound_constraints=true only for CubicSplinePulse (LinearSplinePulse interior is linear, knot bounds suffice)."
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
including pulse-type-dependent `R_u` / `R_du` defaults (see the base method's
docstring for the full discussion, including how to attach
`Piccolissimo.HermiteBendingEnergyRegularizer` for cubic-spline smoothness),
plus:
- `du_bounds::Union{Nothing, Vector{Float64}}=nothing`: Per-drive bounds on derivative magnitude (takes precedence over `du_bound`)
- `free_phase::Bool=false`: Optimize a per-subsystem frame phase alongside the pulse. Applies number-operator rotation `e^{iθ n̂}` to goal states — level `s` acquires phase `s·θ`. Requires `subsystem_levels`.
- `subsystem_levels::Union{Nothing, Vector{Int}}=nothing`: Number of levels per subsystem, required when `free_phase=true`.
- `initial_phases::Union{Nothing, Vector{Float64}}=nothing`: Initial values for the per-subsystem phase variables when `free_phase=true`. Length must equal the number of subsystems.
- `coherent::Bool=true`: If `true`, uses a coherent fidelity objective (phases must align across state pairs). If `false`, uses per-state fidelity.
- `integrator_type::Union{Nothing,Symbol}=nothing`: Integrator backend choice. `nothing` (default)
  infers by pulse kind — `ZeroOrderPulse` is silent PWC; spline pulses guard against the silent-PWC
  trap (cubic errors, linear warns; issue #275). `:pwc` explicitly requests the **piecewise-constant**
  `BilinearIntegrator` — allowed with any pulse, acknowledged by warning for splines. The former
  `:spline` and `:ensemble` values raise informative errors: `:spline` silently returned the PWC
  integrator, and `:ensemble` referenced a type that was never defined. For spline-faithful dynamics
  pass `integrator = Piccolissimo.SplineIntegrator(...)` explicitly.
- `parallel_backend::Symbol=:manual`: **Inert.** Its only consumer was the removed `:ensemble`
  branch; setting it to anything other than `:manual` warns and has no effect. Pass a
  parallel integrator via `integrator` instead.
"""
function _spline_pulse_problem(
    qtraj::MultiKetTrajectory{<:AbstractSplinePulse},
    N_or_times::Union{Nothing,Int,AbstractVector{<:Real}} = nothing;
    integrator::Union{Nothing,AbstractIntegrator,Vector{<:AbstractIntegrator}} = nothing,
    integrator_type::Union{Nothing,Symbol} = nothing,  # nothing = infer (see guards); :pwc = acknowledged PWC
    parallel_backend::Symbol = :manual,  # inert — see the docstring
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
    R_bend::Union{Nothing,Real,AbstractVector{<:Real}} = nothing,  # resolved per pulse type below
    constraints::Vector{<:AbstractConstraint} = AbstractConstraint[],
    extra_objectives::Vector{<:AbstractObjective} = AbstractObjective[],
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
    spline_interior_bound_constraints::Bool = false,
    n_interior_bound_points::Int = 3,
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

    # Bending energy (#309): default ON (1e-3) for cubic, undefined on C⁰ families.
    if isnothing(R_bend)
        R_bend_resolved = is_cubic ? 1e-3 : nothing
    else
        _bend_nonzero = R_bend isa Real ? !iszero(R_bend) : any(!iszero, R_bend)
        if !is_cubic && _bend_nonzero
            error(
                "R_bend is only defined for CubicSplinePulse (bending energy ∫u″²dt " *
                "requires a C¹ Hermite parameterization). For LinearSplinePulse the " *
                "second derivative is distributional — use R_du/du_bound instead.",
            )
        end
        R_bend_resolved = _bend_nonzero ? R_bend : nothing
    end

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
    N = base_traj.N  # Get actual number of knot points

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
        if parallel_backend !== :manual
            @warn "`parallel_backend = :$parallel_backend` has no effect: its only consumer \
                   was the `:ensemble` integrator branch, which was never implemented. Pass a \
                   parallel integrator via `integrator = ...` instead." maxlog = 1
        end

        # Spline pulses must not silently land on the PWC default: BilinearIntegrator
        # never reads :du, so a CubicSplinePulse would optimize a piecewise-constant
        # waveform while the name promises a spline (issue #275). Explicit
        # `integrator_type = :pwc` is the acknowledged escape hatch.
        if qtraj.pulse isa CubicSplinePulse && integrator_type !== :pwc
            error(
                "CubicSplinePulse defaults are not allowed: the default PWC backend " *
                "silently drops :du (a cubic problem would optimize a piecewise-constant " *
                "waveform, not a spline — issue #275). Either pass a spline-faithful " *
                "integrator (Piccolissimo.SplineIntegrator) or explicitly request the PWC " *
                "backend with `integrator_type = :pwc`.",
            )
        end
        if qtraj.pulse isa CubicSplinePulse && integrator_type === :pwc
            @warn "CubicSplinePulse with the PWC backend (`integrator_type = :pwc`): the " *
                  "dynamics treat the drive as piecewise constant and ignore :du — the " *
                  "optimized waveform differs from the cubic spline the pulse object " *
                  "describes. Acknowledged because you asked for it explicitly." maxlog = 1
        elseif qtraj.pulse isa AbstractSplinePulse && integrator_type !== :pwc
            @warn "SplinePulseProblem default BilinearIntegrator with $(typeof(qtraj.pulse).name.name): use SplineIntegrator from Piccolissimo for correct spline physics (Bilinear is PWC, ignores :du)." maxlog =
                1
        end

        # `integrator_type` names what you actually get. `:pwc` is the only backend Piccolo
        # ships: `BilinearIntegrator`, which models the drive as piecewise constant on each
        # interval. There is deliberately no `:spline` value — see the errors below.
        if isnothing(integrator_type) || integrator_type === :pwc
            dynamics_integrators = BilinearIntegrator(qtraj, N)
        elseif integrator_type === :spline
            error(
                """
                `integrator_type = :spline` is not available in Piccolo.

                It previously accepted this value and silently returned a
                **piecewise-constant** `BilinearIntegrator` instead — so a spline pulse was
                optimized against PWC dynamics while the name claimed otherwise. Measured
                cost of that mismatch: the optimizer reports ~1e-8 infidelity for a pulse
                that actually achieves ~1e-3 (see `rollout_divergence`).

                Piccolo ships no spline integrator. Either:
                  • pass one explicitly (requires Piccolissimo):
                      using Piccolissimo
                      integrator = SplineIntegrator(qtraj, N; spline_order = $(_get_spline_order(qtraj.pulse)))
                      SplinePulseProblem(qtraj, N; integrator = integrator, ...)
                  • or request the PWC backend by its real name, `integrator_type = :pwc`,
                    accepting that the pulse is integrated as piecewise constant.
                """,
            )
        elseif integrator_type === :ensemble
            error("""
                  `integrator_type = :ensemble` was never implemented.

                  It referenced an `EnsembleSplineIntegrator` that is defined nowhere in
                  Piccolo or DirectTrajOpt, so this value could only ever throw an
                  `UndefVarError` — which is why `test/jet.jl` ran with `broken = true`.

                  For parallel multi-ket dynamics, pass an integrator explicitly from
                  Piccolissimo. For the shipped backend, use `integrator_type = :pwc`.
                  """)
        else
            error(
                "unknown `integrator_type = :$integrator_type`. " *
                "Piccolo ships one backend: `:pwc` (`BilinearIntegrator`). " *
                "Pass `integrator = ...` for anything else.",
            )
        end

        if !(dynamics_integrators isa AbstractVector)
            dynamics_integrators = AbstractIntegrator[dynamics_integrators]
        else
            dynamics_integrators = AbstractIntegrator[dynamics_integrators...]
        end
    elseif integrator isa AbstractIntegrator
        if qtraj.pulse isa CubicSplinePulse && integrator isa BilinearIntegrator
            error(
                "CubicSplinePulse with BilinearIntegrator: Bilinear never reads :du (H1). Use SplineIntegrator.",
            )
        end
        dynamics_integrators = AbstractIntegrator[integrator]
    else
        for integ in integrator
            if qtraj.pulse isa CubicSplinePulse && integ isa BilinearIntegrator
                error(
                    "CubicSplinePulse with BilinearIntegrator: Bilinear never reads :du (H1). Use SplineIntegrator.",
                )
            end
        end
        dynamics_integrators = AbstractIntegrator[integrator...]
    end

    # Get control names
    du_sym = Symbol(:d, control_sym)

    # Build objective: coherent fidelity for ensemble (with optional free phase)
    J = if free_phase && !isnothing(goals_fn)
        CoherentKetFreePhaseInfidelityObjective(
            goals_fn,
            snames,
            θ_names,
            traj;
            Q = Q,
            weights = weights,
        )
    else
        _ensemble_ket_objective(qtraj, traj, snames, weights, goals, Q; coherent = coherent)
    end

    # Add regularization for control and derivative
    J += QuadraticRegularizer(control_sym, traj, R_u_resolved)
    J += QuadraticRegularizer(du_sym, traj, R_du_resolved)

    # Bending energy (#309): default smoothness term for cubic-Hermite pulses.
    if !isnothing(R_bend_resolved)
        J += HermiteBendingEnergyRegularizer(
            traj;
            R = R_bend_resolved,
            control_name = control_sym,
            derivative_name = du_sym,
            timestep_name = traj.timestep,
        )
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

    # Compose user-supplied extra objective terms (e.g.,
    # Piccolissimo.HermiteBendingEnergyRegularizer). Same ordering convention
    # as the unitary method: fidelity → built-in regs → piccolo-options → extras.
    for extra_obj in extra_objectives
        J += extra_obj
    end

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

    # Spline interior bounds (H10) — CubicSplinePulse can exceed knot bounds in interior.
    # Slice 3c (#431): CubicSplineBoundConstraint lives in Piccolo
    # (Control.QuantumConstraints.SplineConstraints) — a hard import; the old
    # isdefined(Piccolissimo, ...) soft-dependency sniff is gone.
    if spline_interior_bound_constraints
        if qtraj.pulse isa CubicSplinePulse
            for (drive_idx, (lb, ub)) in enumerate(sys.drive_bounds)
                if isfinite(lb) && isfinite(ub)
                    push!(
                        constraints,
                        CubicSplineBoundConstraint(
                            traj,
                            control_sym,
                            lb,
                            ub;
                            n_interior_points = n_interior_bound_points,
                        ),
                    )
                end
            end
            if _show_details(piccolo_options)
                println(
                    "    added CubicSplineBoundConstraint (n=$(n_interior_bound_points) per segment, H10)",
                )
            end
        else
            @warn "spline_interior_bound_constraints=true only for CubicSplinePulse (LinearSplinePulse interior is linear, knot bounds suffice)."
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

# The wrong-pulse fallback ("use SmoothPulseProblem instead") is now *generated*
# by `@problem_template` from the declaration's `hint`.

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
    qcp = SplinePulseProblem(qtraj, N; Q = 100.0, R = 1e-2, integrator_type = :pwc)

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
    qcp = SplinePulseProblem(
        qtraj,
        N;
        Q = 100.0,
        R = 1e-2,
        du_bound = du_bound,
        integrator_type = :pwc,
    )

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
    qcp_unbounded =
        SplinePulseProblem(qtraj, N; Q = 100.0, R = 1e-2, integrator_type = :pwc)
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

@testitem "integrator_type names the backend it actually returns" begin
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    σx = ComplexF64[0 1; 1 0]
    sys = QuantumSystem(0.01 * ComplexF64[1 0; 0 -1], [σx], [1.0])
    N, T = 21, 10.0
    times = collect(range(0.0, T, length = N))
    pulse = LinearSplinePulse(0.1 * randn(1, N), times)
    ψ0, ψ1 = ComplexF64[1.0, 0.0], ComplexF64[0.0, 1.0]
    mk() = MultiKetTrajectory(sys, pulse, [ψ0, ψ1], [ψ1, ψ0])

    # The default is `:pwc`, and it means what it says: BilinearIntegrator.
    qcp = SplinePulseProblem(mk(), N)
    @test length(filter(i -> i isa BilinearIntegrator, qcp.prob.integrators)) == 2
    @test SplinePulseProblem(mk(), N; integrator_type = :pwc) isa QuantumControlProblem

    # `:spline` used to return the PWC integrator silently. It must now say so instead:
    # optimizing a spline against PWC dynamics is the documented 5-orders-of-magnitude
    # misreporting hazard, so a wrong answer is worse than a refusal.
    err = try
        SplinePulseProblem(mk(), N; integrator_type = :spline)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("not available in Piccolo", err.msg)
    @test occursin("SplineIntegrator", err.msg)   # names the way to actually get one

    # `:ensemble` referenced `EnsembleSplineIntegrator`, which is defined nowhere — the
    # reason `test/jet.jl` carried `broken = true`. It must fail with an explanation, not
    # an UndefVarError.
    err2 = try
        SplinePulseProblem(mk(), N; integrator_type = :ensemble)
        nothing
    catch e
        e
    end
    @test err2 isa ErrorException
    @test occursin("never implemented", err2.msg)

    # An unknown value lists what is valid.
    @test_throws ErrorException SplinePulseProblem(mk(), N; integrator_type = :nonsense)

    # `parallel_backend` is inert now that the ensemble branch is gone; say so.
    # `match_mode = :any` because DirectTrajOpt also warns about missing Δt bounds here,
    # and the default exact-match mode would fail on that unrelated record.
    @test_logs (:warn, r"has no effect") match_mode = :any SplinePulseProblem(
        mk(),
        N;
        parallel_backend = :threads,
    )
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

    @test sampling_qcp isa SamplingProblem
    @test sampling_qcp isa AbstractQuantumControlProblem
    @test inner(sampling_qcp) isa QuantumControlProblem
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
    # (declare the PWC backend so the #275 guard doesn't fire first)
    @test_throws "Global variable :δ not found" SplinePulseProblem(
        qtraj,
        N;
        Q = 100.0,
        R = 1e-2,
        global_bounds = Dict(:δ => 0.5),  # δ doesn't exist in trajectory
        integrator_type = :pwc,
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

    # The free-phase branch must honor trajectory weights (issue #263)
    ψp = ComplexF64[1.0, 1.0] / √2
    ψm = ComplexF64[1.0, -1.0] / √2
    pulse_det =
        LinearSplinePulse(0.1 * reshape(cos.(2π .* (0:(N-1)) ./ (N - 1)), 1, N), times)

    function free_phase_objective_value(ws)
        qt = MultiKetTrajectory(sys, pulse_det, [ψ0, ψ1, ψp], [ψ1, ψ0, ψm]; weights = ws)
        p = SplinePulseProblem(
            qt,
            N;
            Q = 100.0,
            R = 1e-2,
            free_phase = true,
            subsystem_levels = [2],
        )
        objective_value(p.prob.objective, p.prob.trajectory)
    end

    @test free_phase_objective_value([0.8, 0.1, 0.1]) !=
          free_phase_objective_value([0.1, 0.1, 0.8])
    @test free_phase_objective_value(fill(1 / 3, 3)) ===
          free_phase_objective_value(fill(1.0, 3))
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
# Cubic-spline default-regularization tests
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

    qcp = SplinePulseProblem(qtraj, N; Q = 100.0, integrator_type = :pwc)  # no explicit R_u, R_du

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

@testitem "SplinePulseProblem CubicSplinePulse: defaults do not over-regularize X gate" begin
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    # With the new defaults (R_u = R_du = 0 for CubicSplinePulse), a simple
    # single-qubit X gate problem should construct, evaluate, and produce a
    # finite, non-negative objective whose regularisation term is exactly zero
    # — confirming the empirically-validated headline of PR #214 (the
    # cubic-spline `R = 1e-2` default was stalling fidelity at ≈ 0.91 vs 0.98+).
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

    qcp = SplinePulseProblem(qtraj, N; Q = 100.0, integrator_type = :pwc)
    @test qcp isa QuantumControlProblem

    # The QuadraticRegularizer contribution must be exactly zero on any
    # trajectory state — the defaults zero out R_u and R_du.
    composite = qcp.prob.objective
    quad_regs = filter(
        x -> x isa DirectTrajOpt.QuadraticRegularizer,
        composite isa DirectTrajOpt.CompositeObjective ? composite.objectives : [composite],
    )
    reg_J = sum(objective_value(r, qcp.prob.trajectory) for r in quad_regs; init = 0.0)
    @test reg_J == 0.0

    J = objective_value(qcp.prob.objective, qcp.prob.trajectory)
    @test isfinite(J)
    @test J ≥ 0.0
end

# ============================================================================= #
# extra_objectives wiring tests
# ============================================================================= #

@testitem "SplinePulseProblem extra_objectives default is identity (unitary)" begin
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    # Passing `extra_objectives = AbstractObjective[]` (the default) must
    # produce the same composite-objective term count as omitting the kwarg.
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

    qcp_default = SplinePulseProblem(qtraj, N; Q = 100.0, integrator_type = :pwc)
    qcp_empty = SplinePulseProblem(
        qtraj,
        N;
        Q = 100.0,
        extra_objectives = AbstractObjective[],
        integrator_type = :pwc,
    )

    n_default =
        qcp_default.prob.objective isa DirectTrajOpt.CompositeObjective ?
        length(qcp_default.prob.objective.objectives) : 1
    n_empty =
        qcp_empty.prob.objective isa DirectTrajOpt.CompositeObjective ?
        length(qcp_empty.prob.objective.objectives) : 1
    @test n_default == n_empty

    J_default = objective_value(qcp_default.prob.objective, qcp_default.prob.trajectory)
    J_empty = objective_value(qcp_empty.prob.objective, qcp_empty.prob.trajectory)
    @test J_default ≈ J_empty
end

@testitem "SplinePulseProblem extra_objectives appends and evaluates (unitary)" begin
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    # Passing one extra QuadraticRegularizer (on an already-present trajectory
    # variable) must (a) increase the composite term count by one and
    # (b) increase the total objective value by exactly that regularizer's
    # contribution.
    σx = ComplexF64[0 1; 1 0]
    σz = ComplexF64[1 0; 0 -1]
    sys = QuantumSystem(0.01 * σz, [σx], [1.0])

    T = 10.0
    N = 11
    times = collect(range(0.0, T, length = N))
    amps = 0.5 * randn(1, N)  # nonzero so the regularizer actually contributes
    derivs = zeros(1, N)
    pulse = CubicSplinePulse(amps, derivs, times)
    U_goal = ComplexF64[0 1; 1 0]
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)

    qcp_baseline = SplinePulseProblem(qtraj, N; Q = 100.0, integrator_type = :pwc)

    # Build the extra regularizer against the actual trajectory we will use.
    # The unitary template emits this same trajectory via NamedTrajectory(qtraj, N).
    extra_R = 1e-1
    traj_for_extra = qcp_baseline.prob.trajectory
    extra_reg = QuadraticRegularizer(:u, traj_for_extra, extra_R)

    qcp_extra = SplinePulseProblem(
        qtraj,
        N;
        Q = 100.0,
        extra_objectives = AbstractObjective[extra_reg],
        integrator_type = :pwc,
    )

    n_baseline =
        qcp_baseline.prob.objective isa DirectTrajOpt.CompositeObjective ?
        length(qcp_baseline.prob.objective.objectives) : 1
    n_extra =
        qcp_extra.prob.objective isa DirectTrajOpt.CompositeObjective ?
        length(qcp_extra.prob.objective.objectives) : 1
    @test n_extra == n_baseline + 1

    # Sanity-check the value: extra objective is summed in linearly.
    J_baseline = objective_value(qcp_baseline.prob.objective, qcp_baseline.prob.trajectory)
    J_extra = objective_value(qcp_extra.prob.objective, qcp_extra.prob.trajectory)
    extra_contribution = objective_value(extra_reg, qcp_extra.prob.trajectory)
    @test J_extra ≈ J_baseline + extra_contribution
    @test extra_contribution > 0.0  # nonzero pulse + nonzero R
end

@testitem "SplinePulseProblem extra_objectives appends and evaluates (MultiKet)" begin
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    # Mirror of the unitary test for the MultiKetTrajectory method signature.
    σx = ComplexF64[0 1; 1 0]
    σz = ComplexF64[1 0; 0 -1]
    sys = QuantumSystem(0.01 * σz, [σx], [1.0])

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]

    T = 10.0
    N = 11
    times = collect(range(0.0, T, length = N))
    pulse = LinearSplinePulse(0.5 * randn(1, N), times)
    qtraj = MultiKetTrajectory(sys, pulse, [ψ0, ψ1], [ψ1, ψ0])

    qcp_baseline = SplinePulseProblem(qtraj, N; Q = 100.0, R = 1e-2)

    extra_R = 1e-1
    traj_for_extra = qcp_baseline.prob.trajectory
    extra_reg = QuadraticRegularizer(:u, traj_for_extra, extra_R)

    qcp_extra = SplinePulseProblem(
        qtraj,
        N;
        Q = 100.0,
        R = 1e-2,
        extra_objectives = AbstractObjective[extra_reg],
    )

    n_baseline =
        qcp_baseline.prob.objective isa DirectTrajOpt.CompositeObjective ?
        length(qcp_baseline.prob.objective.objectives) : 1
    n_extra =
        qcp_extra.prob.objective isa DirectTrajOpt.CompositeObjective ?
        length(qcp_extra.prob.objective.objectives) : 1
    @test n_extra == n_baseline + 1

    J_baseline = objective_value(qcp_baseline.prob.objective, qcp_baseline.prob.trajectory)
    J_extra = objective_value(qcp_extra.prob.objective, qcp_extra.prob.trajectory)
    extra_contribution = objective_value(extra_reg, qcp_extra.prob.trajectory)
    @test J_extra ≈ J_baseline + extra_contribution
end

@testitem "SplinePulseProblem #275 guards: cubic defaults, explicit integrator paths" begin
    using Piccolo

    σx = ComplexF64[0 1; 1 0]
    σz = ComplexF64[1 0; 0 -1]
    T, N = 10.0, 11
    times = collect(range(0.0, T, length = N))
    sys = QuantumSystem(0.01 * σz, [σx], [1.0])

    lin = LinearSplinePulse(0.1 * randn(1, N), times)
    cub = CubicSplinePulse(0.1 * randn(1, N), zeros(1, N), times)
    U_goal = σx
    uq_lin = UnitaryTrajectory(sys, lin, U_goal)
    uq_cub = UnitaryTrajectory(sys, cub, U_goal)

    # CubicSplinePulse + no integrator declared: ERROR (issue #275 — the default
    # PWC backend would silently drop :du and optimize a different waveform).
    err = try
        SplinePulseProblem(uq_cub, N)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("CubicSplinePulse defaults are not allowed", err.msg)
    @test occursin("integrator_type = :pwc", err.msg)

    # global_names without an integrator is an error on both trajectory methods
    err = try
        SplinePulseProblem(uq_lin, N; global_names = [:δ])
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("global_names requires a custom integrator", err.msg)

    ψ0, ψ1 = ComplexF64[1.0, 0.0], ComplexF64[0.0, 1.0]
    mk = MultiKetTrajectory(sys, lin, [ψ0, ψ1], [ψ1, ψ0])
    err = try
        SplinePulseProblem(mk, N; global_names = [:δ])
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("global_names requires a custom integrator", err.msg)

    # explicit single-integrator kwarg: fine for a linear spline (warns on the
    # PWC default elsewhere), H1 error for cubic + BilinearIntegrator
    integ = BilinearIntegrator(uq_lin, N)
    qcp = SplinePulseProblem(uq_lin, N; integrator = integ)
    @test qcp isa QuantumControlProblem

    err = try
        SplinePulseProblem(uq_cub, N; integrator = BilinearIntegrator(uq_cub, N))
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("H1", err.msg)

    # explicit vector-of-integrators kwarg: same contract per element
    qcp_vec = SplinePulseProblem(uq_lin, N; integrator = [integ])
    @test qcp_vec isa QuantumControlProblem

    err = try
        SplinePulseProblem(uq_cub, N; integrator = [BilinearIntegrator(uq_cub, N)])
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("H1", err.msg)
end

@testitem "SplinePulseProblem du_bounds vectors, verbose details, KetTrajectory free_phase" begin
    using Piccolo

    σx = ComplexF64[0 1; 1 0]
    σz = ComplexF64[1 0; 0 -1]
    T, N = 10.0, 11
    times = collect(range(0.0, T, length = N))
    # two drives so per-drive du_bounds differ from the scalar fill
    sys = QuantumSystem(0.01 * σz, [σx, σx], [1.0, 1.0])
    U_goal = σx

    lin2 = LinearSplinePulse(0.1 * randn(2, N), times)
    cub2 = CubicSplinePulse(0.1 * randn(2, N), zeros(2, N), times)
    uq_lin2 = UnitaryTrajectory(sys, lin2, U_goal)
    uq_cub2 = UnitaryTrajectory(sys, cub2, U_goal)

    # per-drive du_bounds: LinearSplinePulse derivative-add path
    qcp = SplinePulseProblem(uq_lin2, N; du_bounds = [5.0, 2.0])
    lo, hi = get_trajectory(qcp).bounds[:du]
    @test lo ≈ [-5.0, -2.0]
    @test hi ≈ [5.0, 2.0]

    # per-drive du_bounds: CubicSplinePulse update_bound! path
    qcp2 = SplinePulseProblem(uq_cub2, N; du_bounds = [5.0, 2.0], integrator_type = :pwc)
    lo2, hi2 = get_trajectory(qcp2).bounds[:du]
    @test lo2 ≈ [-5.0, -2.0]
    @test hi2 ≈ [5.0, 2.0]

    # detailed display exercises the du-bounds and DerivativeIntegrator printlns
    opts = PiccoloOptions(; display = :detailed)
    @test SplinePulseProblem(
        uq_cub2,
        N;
        du_bound = 4.0,
        integrator_type = :pwc,
        piccolo_options = opts,
    ) isa QuantumControlProblem
    @test SplinePulseProblem(uq_lin2, N; piccolo_options = opts) isa QuantumControlProblem

    # same printlns on the MultiKetTrajectory method
    ψ0, ψ1 = ComplexF64[1.0, 0.0], ComplexF64[0.0, 1.0]
    mk_cub2 = MultiKetTrajectory(sys, cub2, [ψ0, ψ1], [ψ1, ψ0])
    mk_lin2 = MultiKetTrajectory(sys, lin2, [ψ0, ψ1], [ψ1, ψ0])
    @test SplinePulseProblem(
        mk_cub2,
        N;
        du_bound = 4.0,
        integrator_type = :pwc,
        piccolo_options = opts,
    ) isa QuantumControlProblem
    @test SplinePulseProblem(mk_lin2, N; piccolo_options = opts) isa QuantumControlProblem

    # single KetTrajectory free-phase: subsystem_levels gate + φ globals +
    # KetFreePhaseInfidelityObjective
    sys1 = QuantumSystem(0.01 * σz, [σx], [1.0])
    lin1 = LinearSplinePulse(0.1 * randn(1, N), times)
    ψ0, ψ1 = ComplexF64[1.0, 0.0], ComplexF64[0.0, 1.0]
    kq = KetTrajectory(sys1, lin1, ψ0, ψ1)
    qcp3 = SplinePulseProblem(kq, N; free_phase = true, subsystem_levels = [2])
    @test qcp3 isa QuantumControlProblem
    @test haskey(get_trajectory(qcp3).global_components, :φ_1)

    err = try
        SplinePulseProblem(kq, N; free_phase = true)
        nothing
    catch e
        e
    end
    @test err isa AssertionError
end

@testitem "SplinePulseProblem global_params system data and spline_interior_bound_constraints" begin
    using Piccolo

    σx = ComplexF64[0 1; 1 0]
    σz = ComplexF64[1 0; 0 -1]
    T, N = 10.0, 11
    times = collect(range(0.0, T, length = N))
    U_goal = σx
    ψ0, ψ1 = ComplexF64[1.0, 0.0], ComplexF64[0.0, 1.0]

    # a system carrying global_params seeds global_data on BOTH trajectory
    # methods (the φ free-phase globals merge on top of it when asked for).
    sysg = QuantumSystem(0.01 * σz, [σx], [1.0]; global_params = (δ = 0.2,))
    lin = LinearSplinePulse(0.1 * randn(1, N), times)

    uq = UnitaryTrajectory(sysg, lin, U_goal)
    qcp = SplinePulseProblem(uq, N)
    @test haskey(get_trajectory(qcp).global_components, :δ)

    mk = MultiKetTrajectory(sysg, lin, [ψ0, ψ1], [ψ1, ψ0])
    qcp2 = SplinePulseProblem(mk, N)
    @test haskey(get_trajectory(qcp2).global_components, :δ)

    sys1 = QuantumSystem(0.01 * σz, [σx], [1.0])
    cub = CubicSplinePulse(0.1 * randn(1, N), zeros(1, N), times)
    lin1 = LinearSplinePulse(0.1 * randn(1, N), times)

    # spline_interior_bound_constraints binds via Piccolo's own
    # SplineConstraints submodule (open-core slice 3c, #431): a HARD import,
    # no private-package sniff. Piccolissimo is not — and cannot be — a
    # dependency of this package, so this test env IS the
    # no-private-package proof (AC 2): the constraint must land and no
    # fallback warning may fire.
    bounded_sys = QuantumSystem(0.01 * σz, [σx], [(-1.0, 1.0)])
    qcp_bounded = SplinePulseProblem(
        UnitaryTrajectory(bounded_sys, cub, U_goal),
        N;
        integrator_type = :pwc,
        spline_interior_bound_constraints = true,
    )
    @test any(c -> c isa CubicSplineBoundConstraint, qcp_bounded.prob.constraints)
    # ...and silently (no fallback warning) for the multi-ket method too.
    mk_cub_bounded = MultiKetTrajectory(bounded_sys, cub, [ψ0, ψ1], [ψ1, ψ0])
    qcp_bounded_mk = SplinePulseProblem(
        mk_cub_bounded,
        N;
        integrator_type = :pwc,
        spline_interior_bound_constraints = true,
    )
    @test any(c -> c isa CubicSplineBoundConstraint, qcp_bounded_mk.prob.constraints)

    # LinearSplinePulse: knot bounds suffice; the kwarg warns and is inert
    @test_logs (:warn, r"only for CubicSplinePulse") match_mode = :any SplinePulseProblem(
        UnitaryTrajectory(sys1, lin1, U_goal),
        N;
        spline_interior_bound_constraints = true,
    )

    mk_lin = MultiKetTrajectory(sys1, lin1, [ψ0, ψ1], [ψ1, ψ0])
    @test_logs (:warn, r"only for CubicSplinePulse") match_mode = :any SplinePulseProblem(
        mk_lin,
        N;
        spline_interior_bound_constraints = true,
    )
end

@testitem "SplinePulseProblem MultiKet cubic guards, du bounds, explicit integrators" begin
    using Piccolo

    σx = ComplexF64[0 1; 1 0]
    σz = ComplexF64[1 0; 0 -1]
    T, N = 10.0, 11
    times = collect(range(0.0, T, length = N))
    sys = QuantumSystem(0.01 * σz, [σx], [1.0])
    ψ0, ψ1 = ComplexF64[1.0, 0.0], ComplexF64[0.0, 1.0]

    lin = LinearSplinePulse(0.1 * randn(1, N), times)
    cub = CubicSplinePulse(0.1 * randn(1, N), zeros(1, N), times)
    mk_lin = MultiKetTrajectory(sys, lin, [ψ0, ψ1], [ψ1, ψ0])
    mk_cub = MultiKetTrajectory(sys, cub, [ψ0, ψ1], [ψ1, ψ0])

    # cubic + no integrator_type on the multi-ket method: same #275 error
    err = try
        SplinePulseProblem(mk_cub, N)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("CubicSplinePulse defaults are not allowed", err.msg)

    # cubic + acknowledged :pwc backend: warns, constructs
    @test_logs (:warn, r"PWC backend") match_mode = :any SplinePulseProblem(
        mk_cub,
        N;
        integrator_type = :pwc,
    )

    # scalar du_bound: linear add-derivatives path and cubic update_bound! path
    qcp = SplinePulseProblem(mk_lin, N; du_bound = 5.0)
    lo, hi = get_trajectory(qcp).bounds[:du]
    @test lo ≈ [-5.0] && hi ≈ [5.0]

    qcp2 = SplinePulseProblem(mk_cub, N; du_bound = 5.0, integrator_type = :pwc)
    lo2, hi2 = get_trajectory(qcp2).bounds[:du]
    @test lo2 ≈ [-5.0] && hi2 ≈ [5.0]

    # per-drive du_bounds vector on the multi-ket linear path
    sys2 = QuantumSystem(0.01 * σz, [σx, σx], [1.0, 1.0])
    lin2 = LinearSplinePulse(0.1 * randn(2, N), times)
    mk_lin2 = MultiKetTrajectory(sys2, lin2, [ψ0, ψ1], [ψ1, ψ0])
    qcp3 = SplinePulseProblem(mk_lin2, N; du_bounds = [5.0, 2.0])
    lo3, hi3 = get_trajectory(qcp3).bounds[:du]
    @test lo3 ≈ [-5.0, -2.0] && hi3 ≈ [5.0, 2.0]

    # single-integrator kwarg branch: linear + a single integrator constructs
    # (the vector-returning BilinearIntegrator(multiket, N) takes the vector
    # branch; the single branch is the parallel-integrator contract).
    uq_lin = UnitaryTrajectory(sys, lin, σx)
    qcp4 = SplinePulseProblem(mk_lin, N; integrator = BilinearIntegrator(uq_lin, N))
    @test qcp4 isa QuantumControlProblem

    # cubic + a single BilinearIntegrator: the H1 guard fires before anything
    # is built
    uq_cub = UnitaryTrajectory(sys, cub, σx)
    err = try
        SplinePulseProblem(mk_cub, N; integrator = BilinearIntegrator(uq_cub, N))
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("H1", err.msg)

    # 2-state ensemble: BilinearIntegrator returns a VECTOR — the vector branch
    integs = BilinearIntegrator(mk_lin, N)
    @test integs isa AbstractVector
    qcp5 = SplinePulseProblem(mk_lin, N; integrator = integs)
    @test qcp5 isa QuantumControlProblem

    err = try
        SplinePulseProblem(mk_cub, N; integrator = BilinearIntegrator(mk_cub, N))
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("H1", err.msg)
end

@testitem "coverage: _get_spline_order + verbose construction path" begin
    using Piccolo

    @test Piccolo.Control.ProblemTemplates._get_spline_order(
        LinearSplinePulse(randn(1, 5), collect(0.0:0.1:0.4)),
    ) == 1
    @test Piccolo.Control.ProblemTemplates._get_spline_order(
        CubicSplinePulse(randn(1, 5), randn(1, 5), collect(0.0:0.1:0.4)),
    ) == 3

    # the :detailed display level exercises the verbose println branches
    sys = QuantumSystem(0.1 * PAULIS[:Z], [PAULIS[:X]], [(-1.0, 1.0)])
    times = collect(range(0.0, 1.0; length = 11))
    pulse = LinearSplinePulse(0.1 .* randn(1, 11), times)
    qtraj = UnitaryTrajectory(sys, pulse, GATES[:X])
    qcp = SplinePulseProblem(
        qtraj,
        11;
        Q = 100.0,
        R = 1e-2,
        piccolo_options = PiccoloOptions(display = :detailed),
    )
    @test qcp isa QuantumControlProblem
end

# ============================================================================= #
# Bending-energy default tests (#309)
# ============================================================================= #

@testitem "SplinePulseProblem R_bend default ON for cubic, absent when opted out" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

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

    terms(qcp) =
        let obj = qcp.prob.objective
            obj isa DirectTrajOpt.CompositeObjective ? obj.objectives : [obj]
        end

    # Default: bending term present with R = [1e-3]
    qcp_default = SplinePulseProblem(qtraj, N; Q = 100.0, integrator_type = :pwc)
    bend_default = filter(x -> x isa HermiteBendingEnergyRegularizer, terms(qcp_default))
    @test length(bend_default) == 1
    @test bend_default[1].R == [1e-3]

    # Opt-out: R_bend = 0 → no bending term
    qcp_off = SplinePulseProblem(qtraj, N; Q = 100.0, integrator_type = :pwc, R_bend = 0)
    bend_off = filter(x -> x isa HermiteBendingEnergyRegularizer, terms(qcp_off))
    @test isempty(bend_off)

    # Linear spline + explicit nonzero R_bend → error
    lin_pulse = LinearSplinePulse(0.1 * randn(1, N), times)
    lin_qtraj = UnitaryTrajectory(sys, lin_pulse, U_goal)
    err = try
        SplinePulseProblem(lin_qtraj, N; Q = 100.0, R_bend = 0.5)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("CubicSplinePulse", err.msg)
end

@testitem "SplinePulseProblem MultiKet R_bend default ON for cubic" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    σx = ComplexF64[0 1; 1 0]
    σz = ComplexF64[1 0; 0 -1]
    sys = QuantumSystem(0.01 * σz, [σx], [1.0])

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]

    T = 10.0
    N = 11
    times = collect(range(0.0, T, length = N))
    pulse = CubicSplinePulse(0.1 * randn(1, N), zeros(1, N), times)
    qtraj = MultiKetTrajectory(sys, pulse, [ψ0, ψ1], [ψ1, ψ0])

    qcp = SplinePulseProblem(qtraj, N; Q = 100.0, integrator_type = :pwc)
    terms_ = let obj = qcp.prob.objective
        obj isa DirectTrajOpt.CompositeObjective ? obj.objectives : [obj]
    end
    bend = filter(x -> x isa HermiteBendingEnergyRegularizer, terms_)
    @test length(bend) == 1
    @test bend[1].R == [1e-3]
end

@testitem "bending closed form ≈ fine-mesh FD; grid-refinement invariance (#309)" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using Random

    # One fixed continuous cubic pulse discretized at several knot counts:
    # the closed-form bending energy must be a property of the PULSE, not the
    # grid (Riemann property), and must agree with an independent fine-mesh
    # finite-difference evaluation of the same Hermite spline.
    Random.seed!(309)
    f(t) = sin(2π * t) + 0.3 * cos(6π * t)   # smooth, fixed
    df(t) = 2π * cos(2π * t) - 1.8π * sin(6π * t)
    T = 1.0

    closed_form_bend(N) = begin
        times = collect(range(0.0, T, length = N))
        u = reshape(f.(times), 1, N)
        du = reshape(df.(times), 1, N)
        Δt = fill(T / (N - 1), N)
        traj = NamedTrajectory((u = u, du = du, Δt = Δt); timestep = :Δt, controls = :u)
        reg = HermiteBendingEnergyRegularizer(traj; R = 1.0)
        DirectTrajOpt.objective_value(reg, traj)
    end

    # J = (1/2)∫u″² — the regularizer's R=1 convention
    # (J = (R/2)·∫f''² with R=1 → J = 0.5·∫f''²)
    d2f(t) = -(2π)^2 * sin(2π * t) - 10.8π^2 * cos(6π * t)  # 0.3·(6π)² = 10.8π²
    N_ref = 200_000
    h = T / (N_ref - 1)
    J_exact = 0.5 * h * sum(abs2, d2f.(collect(0:h:(T-h))))

    Js = [closed_form_bend(N) for N in (51, 201, 801)]
    rel_drift = maximum(Js) / minimum(Js) - 1.0
    rel_err = abs(Js[2] - J_exact) / J_exact

    @test rel_drift < 0.02   # grid-refinement invariance (Riemann property)
    @test rel_err < 0.01     # closed form ≈ independent quadrature
end

@testitem "shape_metrics quartet sanity (#309)" begin
    using Piccolo
    using Random
    Random.seed!(309)

    times = collect(range(0.0, 1.0, length = 21))
    amps = 0.1 * randn(2, 21)
    derivs = zeros(2, 21)
    pulse = CubicSplinePulse(amps, derivs, times)

    m = shape_metrics(pulse)
    @test length(m.bend) == 2
    @test all(m.bend .> 0)
    @test all(m.int_u2 .> 0)
    @test all(m.max_du .> 0)
    @test m.T ≈ 1.0
    @test m.parameterization == 3

    lin = LinearSplinePulse(0.1 * randn(2, 21), times)
    ml = shape_metrics(lin)
    @test ml.parameterization == 1

    # bend on a C² trajectory: linear u with constant slope → ~0.
    # (Construct via cubic knots that reproduce y = t exactly.)
    lin_times = collect(range(0.0, 2.0, length = 6))
    lin_u = reshape(collect(lin_times), 1, 6)
    lin_du = ones(1, 6)
    lin_pulse = CubicSplinePulse(lin_u, lin_du, lin_times)
    m_c2 = shape_metrics(lin_pulse; mesh = 2^12)
    @test m_c2.bend[1] < 1e-6
end
