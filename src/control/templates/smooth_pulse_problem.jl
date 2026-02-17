export SmoothPulseProblem

@doc raw"""
    SmoothPulseProblem(qtraj::AbstractQuantumTrajectory{<:ZeroOrderPulse}, N::Int; kwargs...)

Construct a `QuantumControlProblem` for smooth pulse optimization with piecewise constant controls.

**Note**: This problem template is for `ZeroOrderPulse` only. For spline-based pulses
(`LinearSplinePulse`, `CubicSplinePulse`), use `SplinePulseProblem` instead.

The problem adds discrete derivative variables (du, ddu) that:
- Regularize control changes between timesteps
- Enforce smoothness via `DerivativeIntegrator` constraints

# Arguments
- `qtraj::AbstractQuantumTrajectory{<:ZeroOrderPulse}`: Quantum trajectory with piecewise constant pulse
- `N::Int`: Number of timesteps for discretization

# Keyword Arguments
- `integrator::Union{Nothing, AbstractIntegrator, Vector{<:AbstractIntegrator}}=nothing`: Optional custom integrator(s). If not provided, uses BilinearIntegrator (which does not support global variables). A custom integrator is required when `global_names` is specified.
- `global_names::Union{Nothing, Vector{Symbol}}=nothing`: Names of global variables to optimize. Requires a custom integrator (e.g., HermitianExponentialIntegrator from Piccolissimo) that supports global variables.
- `global_bounds::Union{Nothing, Dict{Symbol, Union{Float64, Tuple{Float64, Float64}}}}=nothing`: Bounds for global variables. Keys are variable names, values are either a scalar (symmetric bounds ±value) or a tuple (lower, upper).
- `du_bound::Float64=Inf`: Bound on discrete first derivative (controls jump rate)
- `ddu_bound::Float64=1.0`: Bound on discrete second derivative (controls acceleration)
- `Q::Float64=100.0`: Weight on infidelity/objective
- `R::Float64=1e-2`: Weight on regularization terms (u, u̇, ü)
- `R_u::Union{Float64, Vector{Float64}}=R`: Weight on control regularization
- `R_du::Union{Float64, Vector{Float64}}=R`: Weight on first derivative regularization
- `R_ddu::Union{Float64, Vector{Float64}}=R`: Weight on second derivative regularization
- `constraints::Vector{<:AbstractConstraint}=AbstractConstraint[]`: Additional constraints
- `piccolo_options::PiccoloOptions=PiccoloOptions()`: Piccolo solver options

# Returns
- `QuantumControlProblem`: Wrapper containing quantum trajectory and optimization problem

# Examples
```julia
# Unitary gate synthesis with piecewise constant pulse
sys = QuantumSystem(H_drift, H_drives, drive_bounds)
pulse = ZeroOrderPulse(0.1 * randn(n_drives, N), collect(range(0.0, T, length=N)))
qtraj = UnitaryTrajectory(sys, pulse, U_goal)
qcp = SmoothPulseProblem(qtraj, N; Q=100.0, R=1e-2)
solve!(qcp; max_iter=100)

# Quantum state transfer
pulse = ZeroOrderPulse(0.1 * randn(n_drives, N), collect(range(0.0, T, length=N)))
qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)
qcp = SmoothPulseProblem(qtraj, N; Q=50.0, R=1e-3)
solve!(qcp)
```

See also: [`SplinePulseProblem`](@ref) for spline-based pulses.
"""
function SmoothPulseProblem(
    qtraj::AbstractQuantumTrajectory{<:ZeroOrderPulse},
    N::Int;
    integrator::Union{Nothing,AbstractIntegrator,Vector{<:AbstractIntegrator}} = nothing,
    global_names::Union{Nothing,Vector{Symbol}} = nothing,
    global_bounds::Union{Nothing,Dict{Symbol,<:Union{Float64,Tuple{Float64,Float64}}}} = nothing,
    du_bound::Float64 = Inf,
    ddu_bound::Float64 = 1.0,
    Δt_bounds::Union{Nothing,Tuple{Float64,Float64}} = nothing,
    Q::Float64 = 100.0,
    R::Float64 = 1e-2,
    R_u::Union{Float64,Vector{Float64}} = R,
    R_du::Union{Float64,Vector{Float64}} = R,
    R_ddu::Union{Float64,Vector{Float64}} = R,
    constraints::Vector{<:AbstractConstraint} = AbstractConstraint[],
    piccolo_options::PiccoloOptions = PiccoloOptions(),
)
    if piccolo_options.verbose
        traj_type = split(string(typeof(qtraj).name.name), ".")[end]
        println("    constructing SmoothPulseProblem for $traj_type...")
    end

    # Extract info from quantum trajectory
    sys = get_system(qtraj)
    state_sym = state_name(qtraj)
    control_sym = drive_name(qtraj)

    # Build global_data from system's global_params if present
    global_data = if hasproperty(sys, :global_params) && !isempty(sys.global_params)
        Dict(name => [val] for (name, val) in pairs(sys.global_params))
    else
        nothing
    end

    # Convert quantum trajectory to NamedTrajectory
    base_traj = if isnothing(global_data)
        NamedTrajectory(qtraj, N; Δt_bounds = Δt_bounds)
    else
        NamedTrajectory(qtraj, N; Δt_bounds = Δt_bounds, global_data = global_data)
    end

    # Add control derivatives to trajectory (always 2 derivatives for smooth pulses)
    du_bounds = fill(du_bound, sys.n_drives)
    ddu_bounds = fill(ddu_bound, sys.n_drives)

    traj_smooth = add_control_derivatives(
        base_traj,
        2;  # Always use 2 derivatives
        control_name = control_sym,
        derivative_bounds = (du_bounds, ddu_bounds),
    )

    # Initialize dynamics integrators - handle both single integrator and vector of integrators
    if isnothing(integrator)
        # Check for global_names without integrator
        if !isnothing(global_names) && !isempty(global_names)
            error(
                "global_names requires a custom integrator that supports global variables. " *
                "Use HermitianExponentialIntegrator from Piccolissimo:\n" *
                "  using Piccolissimo\n" *
                "  integrator = HermitianExponentialIntegrator(qtraj, N; global_names=$global_names)\n" *
                "  qcp = SmoothPulseProblem(qtraj, N; integrator=integrator, ...)",
            )
        end
        # Use default BilinearIntegrator for the trajectory type
        default_int = BilinearIntegrator(qtraj, N)
        if default_int isa AbstractVector
            dynamics_integrators = AbstractIntegrator[default_int...]
        else
            dynamics_integrators = AbstractIntegrator[default_int]
        end
    elseif integrator isa AbstractIntegrator
        # Single custom integrator provided
        dynamics_integrators = AbstractIntegrator[integrator]
    else
        # Vector of custom integrators provided
        dynamics_integrators = AbstractIntegrator[integrator...]
    end

    # Get control derivative names
    control_names =
        [name for name ∈ traj_smooth.names if endswith(string(name), string(control_sym))]

    # Build objective: type-specific infidelity + regularization
    J = _state_objective(qtraj, traj_smooth, state_sym, Q)

    # Add regularization for control and derivatives
    J += QuadraticRegularizer(control_names[1], traj_smooth, R_u)
    J += QuadraticRegularizer(control_names[2], traj_smooth, R_du)
    J += QuadraticRegularizer(control_names[3], traj_smooth, R_ddu)

    # Add optional Piccolo constraints and objectives
    J += _apply_piccolo_options(qtraj, piccolo_options, constraints, traj_smooth, state_sym)

    # Start with dynamics integrators
    integrators = copy(dynamics_integrators)

    # Add derivative integrators (always 2 for smooth pulses)
    push!(
        integrators,
        DerivativeIntegrator(control_names[1], control_names[2], traj_smooth),
    )

    push!(
        integrators,
        DerivativeIntegrator(control_names[2], control_names[3], traj_smooth),
    )

    # Note: TimeConsistencyConstraint is auto-applied by DirectTrajOpt when :t and :Δt present

    # Add global bounds constraints if specified
    all_constraints = copy(constraints)
    add_global_bounds_constraints!(
        all_constraints,
        global_bounds,
        traj_smooth;
        verbose = piccolo_options.verbose,
    )

    prob = DirectTrajOptProblem(traj_smooth, J, integrators; constraints = all_constraints)

    return QuantumControlProblem(qtraj, prob)
end

# ============================================================================= #
# EnsembleTrajectory Constructor
# ============================================================================= #

@doc raw"""
    SmoothPulseProblem(qtraj::MultiKetTrajectory{<:ZeroOrderPulse}, N::Int; kwargs...)

Construct a `QuantumControlProblem` for smooth pulse optimization over an ensemble of ket state transfers
with piecewise constant controls.

This handles the case where you want to optimize a single pulse that achieves multiple 
state transfers simultaneously (e.g., |0⟩→|1⟩ and |1⟩→|0⟩ for an X gate via state transfer).

**Note**: This problem template is for `ZeroOrderPulse` only. For spline-based pulses,
use `SplinePulseProblem` instead.

# Arguments
- `qtraj::MultiKetTrajectory{<:ZeroOrderPulse}`: Ensemble of ket state transfers with piecewise constant pulse
- `N::Int`: Number of timesteps for the discretization

# Keyword Arguments
- `integrator::Union{Nothing, AbstractIntegrator, Vector{<:AbstractIntegrator}}=nothing`: Optional custom integrator(s). If not provided, the default `BilinearIntegrator` is used. When `global_names` is specified, you must supply a custom integrator here (i.e., do not rely on the default `BilinearIntegrator`) that supports global variables.
- `global_names::Union{Nothing, Vector{Symbol}}=nothing`: Names of global variables to optimize. Requires a custom integrator provided via `integrator` (e.g., `HermitianExponentialIntegrator` from Piccolissimo) that supports global variables.
- `global_bounds::Union{Nothing, Dict{Symbol, Union{Float64, Tuple{Float64, Float64}}}}=nothing`: Bounds for global variables. Keys are variable names, values are either a scalar (symmetric bounds ±value) or a tuple (lower, upper).
- `du_bound::Float64=Inf`: Bound on discrete first derivative
- `ddu_bound::Float64=1.0`: Bound on discrete second derivative
- `Q::Float64=100.0`: Weight on infidelity/objective
- `R::Float64=1e-2`: Weight on regularization terms (u, u̇, ü)
- `R_u::Union{Float64, Vector{Float64}}=R`: Weight on control regularization
- `R_du::Union{Float64, Vector{Float64}}=R`: Weight on first derivative regularization
- `R_ddu::Union{Float64, Vector{Float64}}=R`: Weight on second derivative regularization
- `constraints::Vector{<:AbstractConstraint}=AbstractConstraint[]`: Additional constraints
- `piccolo_options::PiccoloOptions=PiccoloOptions()`: Piccolo solver options

# Returns
- `QuantumControlProblem{MultiKetTrajectory}`: Wrapper containing ensemble trajectory and optimization problem

# Examples
```julia
# Create ensemble for X gate via state transfer
sys = QuantumSystem(H_drift, H_drives, drive_bounds)
pulse = ZeroOrderPulse(0.1 * randn(n_drives, N), collect(range(0.0, T, length=N)))

ψ0 = ComplexF64[1.0, 0.0]
ψ1 = ComplexF64[0.0, 1.0]

ensemble_qtraj = MultiKetTrajectory(sys, pulse, [ψ0, ψ1], [ψ1, ψ0])
qcp = SmoothPulseProblem(ensemble_qtraj, N; Q=100.0, R=1e-2)
solve!(qcp; max_iter=100)
```

See also: [`SplinePulseProblem`](@ref) for spline-based pulses.
"""
function SmoothPulseProblem(
    qtraj::MultiKetTrajectory{<:ZeroOrderPulse},
    N::Int;
    integrator::Union{Nothing,AbstractIntegrator,Vector{<:AbstractIntegrator}} = nothing,
    global_names::Union{Nothing,Vector{Symbol}} = nothing,
    global_bounds::Union{Nothing,Dict{Symbol,<:Union{Float64,Tuple{Float64,Float64}}}} = nothing,
    du_bound::Float64 = Inf,
    ddu_bound::Float64 = 1.0,
    Δt_bounds::Union{Nothing,Tuple{Float64,Float64}} = nothing,
    Q::Float64 = 100.0,
    R::Float64 = 1e-2,
    R_u::Union{Float64,Vector{Float64}} = R,
    R_du::Union{Float64,Vector{Float64}} = R,
    R_ddu::Union{Float64,Vector{Float64}} = R,
    constraints::Vector{<:AbstractConstraint} = AbstractConstraint[],
    piccolo_options::PiccoloOptions = PiccoloOptions(),
)
    if piccolo_options.verbose
        println(
            "    constructing SmoothPulseProblem for MultiKetTrajectory ($(length(qtraj.initials)) states)...",
        )
    end

    # Extract info from ensemble trajectory
    sys = get_system(qtraj)
    control_sym = drive_name(qtraj)
    snames = state_names(qtraj)
    weights = qtraj.weights
    goals = qtraj.goals

    # Build global_data from system's global_params if present
    global_data = if !isempty(sys.global_params)
        Dict(name => [val] for (name, val) in pairs(sys.global_params))
    else
        nothing
    end

    # Convert quantum trajectory to NamedTrajectory
    base_traj = NamedTrajectory(qtraj, N; Δt_bounds = Δt_bounds, global_data = global_data)

    # Add control derivatives to trajectory
    du_bounds = fill(du_bound, sys.n_drives)
    ddu_bounds = fill(ddu_bound, sys.n_drives)

    traj_smooth = add_control_derivatives(
        base_traj,
        2;  # Always use 2 derivatives
        control_name = control_sym,
        derivative_bounds = (du_bounds, ddu_bounds),
    )

    # Get control derivative names
    control_names =
        [name for name ∈ traj_smooth.names if endswith(string(name), string(control_sym))]

    # Build objective: weighted sum of infidelities for each state
    J = _ensemble_ket_objective(qtraj, traj_smooth, snames, weights, goals, Q)

    # Add regularization for control and derivatives
    J += QuadraticRegularizer(control_names[1], traj_smooth, R_u)
    J += QuadraticRegularizer(control_names[2], traj_smooth, R_du)
    J += QuadraticRegularizer(control_names[3], traj_smooth, R_ddu)

    # Apply piccolo options for each state
    J += _apply_piccolo_options(qtraj, piccolo_options, constraints, traj_smooth, snames)

    # Build integrators: one dynamics integrator per state
    if isnothing(integrator)
        # Check for global_names without integrator
        if !isnothing(global_names) && !isempty(global_names)
            error(
                "global_names requires a custom integrator that supports global variables. " *
                "Use HermitianExponentialIntegrator from Piccolissimo:\n" *
                "  using Piccolissimo\n" *
                "  integrator = HermitianExponentialIntegrator(qtraj, N; global_names=$global_names)\n" *
                "  qcp = SmoothPulseProblem(qtraj, N; integrator=integrator, ...)",
            )
        end
        dynamics_integrators = BilinearIntegrator(qtraj, N)
    elseif integrator isa AbstractIntegrator
        dynamics_integrators = AbstractIntegrator[integrator]
    else
        dynamics_integrators = AbstractIntegrator[integrator...]
    end

    integrators = AbstractIntegrator[dynamics_integrators...]

    # Add derivative integrators (always 2 for smooth pulses)
    push!(
        integrators,
        DerivativeIntegrator(control_names[1], control_names[2], traj_smooth),
    )
    push!(
        integrators,
        DerivativeIntegrator(control_names[2], control_names[3], traj_smooth),
    )

    # Note: TimeConsistencyConstraint is auto-applied by DirectTrajOpt when :t and :Δt present

    # Add global bounds constraints if specified
    all_constraints = copy(constraints)
    add_global_bounds_constraints!(
        all_constraints,
        global_bounds,
        traj_smooth;
        verbose = piccolo_options.verbose,
    )

    prob = DirectTrajOptProblem(traj_smooth, J, integrators; constraints = all_constraints)

    return QuantumControlProblem(qtraj, prob)
end

# ============================================================================= #
# Type-specific helper functions
# ============================================================================= #

# ----------------------------------------------------------------------------- #
# Single trajectory objectives
# ----------------------------------------------------------------------------- #

# Unitary trajectory: single infidelity objective
function _state_objective(
    qtraj::UnitaryTrajectory,
    traj::NamedTrajectory,
    state_sym::Symbol,
    Q::Float64,
)
    U_goal = qtraj.goal
    return UnitaryInfidelityObjective(U_goal, state_sym, traj; Q = Q)
end

# Ket trajectory: single infidelity objective
function _state_objective(
    qtraj::KetTrajectory,
    traj::NamedTrajectory,
    state_sym::Symbol,
    Q::Float64,
)
    return KetInfidelityObjective(state_sym, traj; Q = Q)
end

# Density trajectory: general mixed state fidelity via compact isomorphism
function _state_objective(
    qtraj::DensityTrajectory,
    traj::NamedTrajectory,
    state_sym::Symbol,
    Q::Float64,
)
    ρ_goal = Matrix{ComplexF64}(qtraj.goal)
    return DensityMatrixInfidelityObjective(state_sym, ρ_goal, traj; Q = Q)
end

# Apply Piccolo options with trajectory-type-specific logic
function _apply_piccolo_options(
    qtraj::UnitaryTrajectory,
    piccolo_options::PiccoloOptions,
    constraints::Vector{<:AbstractConstraint},
    traj::NamedTrajectory,
    state_sym::Symbol,
)
    U_goal = qtraj.goal
    return apply_piccolo_options!(
        piccolo_options,
        constraints,
        traj;
        state_names = state_sym,
        state_leakage_indices = U_goal isa EmbeddedOperator ?
                                get_iso_vec_leakage_indices(U_goal) : nothing,
    )
end

function _apply_piccolo_options(
    qtraj::KetTrajectory,
    piccolo_options::PiccoloOptions,
    constraints::Vector{<:AbstractConstraint},
    traj::NamedTrajectory,
    state_sym::Symbol,
)
    return apply_piccolo_options!(
        piccolo_options,
        constraints,
        traj;
        state_names = state_sym,
    )
end

function _apply_piccolo_options(
    qtraj::DensityTrajectory,
    piccolo_options::PiccoloOptions,
    constraints::Vector{<:AbstractConstraint},
    traj::NamedTrajectory,
    state_sym::Symbol,
)
    return apply_piccolo_options!(
        piccolo_options,
        constraints,
        traj;
        state_names = state_sym,
    )
end

# ----------------------------------------------------------------------------- #
# Ensemble ket trajectory objectives
# ----------------------------------------------------------------------------- #

"""
    _ensemble_ket_objective(qtraj::MultiKetTrajectory, traj, state_names, weights, goals, Q)

Create a coherent fidelity objective for ensemble state transfers.

For ensemble trajectories (implementing a gate via multiple state transfers),
we use coherent fidelity:
    F_coherent = |1/n ∑ᵢ ⟨ψᵢ_goal|ψᵢ⟩|²

This requires all state overlaps to have aligned phases, which is essential
for gate implementation (the gate should have a single global phase).
"""
function _ensemble_ket_objective(
    qtraj::MultiKetTrajectory,
    traj::NamedTrajectory,
    snames::Vector{Symbol},
    weights::Vector{Float64},
    goals::Vector,
    Q::Float64,
)
    # Use coherent fidelity - phases must align for gate implementation
    return CoherentKetInfidelityObjective(goals, snames, traj; Q = Q)
end

# ----------------------------------------------------------------------------- #
# Ensemble piccolo options
# ----------------------------------------------------------------------------- #

function _apply_piccolo_options(
    qtraj::MultiKetTrajectory,
    piccolo_options::PiccoloOptions,
    constraints::Vector{<:AbstractConstraint},
    traj::NamedTrajectory,
    snames::Vector{Symbol},
)
    # Apply piccolo options for all state variables in the ensemble
    return apply_piccolo_options!(piccolo_options, constraints, traj; state_names = snames)
end

# ============================================================================= #
# Fallback Error Method
# ============================================================================= #

"""
    SmoothPulseProblem(qtraj::AbstractQuantumTrajectory, N::Int; kwargs...)

Fallback method that provides helpful error for non-ZeroOrderPulse types.
"""
function SmoothPulseProblem(
    qtraj::AbstractQuantumTrajectory{P},
    N::Int;
    kwargs...,
) where {P<:AbstractPulse}
    pulse_type = P
    error(
        """
  SmoothPulseProblem is only for piecewise constant pulses (ZeroOrderPulse).

  You provided a trajectory with pulse type: $(pulse_type)

  For spline-based pulses (LinearSplinePulse, CubicSplinePulse), use SplinePulseProblem instead:
      qcp = SplinePulseProblem(qtraj, N; ...)
  """,
    )
end

# ============================================================================= #
# Tests
# ============================================================================= #

@testitem "SmoothPulseProblem with UnitaryTrajectory" tags = [:experimental] begin
    using DirectTrajOpt
    using LinearAlgebra

    T = 10.0
    N = 50
    sys = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    U_goal = GATES[:H]

    # Create pulse and quantum trajectory
    pulse = ZeroOrderPulse(0.1 * randn(2, N), collect(range(0.0, T, length = N)))
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)
    qcp = SmoothPulseProblem(qtraj, N; Q = 100.0, R = 1e-2)

    @test qcp isa QuantumControlProblem
    @test length(qcp.prob.integrators) == 3  # dynamics + du + ddu
    @test haskey(qcp.prob.trajectory.components, :u)
    @test haskey(qcp.prob.trajectory.components, :du)
    @test haskey(qcp.prob.trajectory.components, :ddu)

    # Test accessors
    @test get_system(qcp) === sys
    @test qtraj.goal === U_goal
    @test get_trajectory(qcp) === qcp.prob.trajectory

    # Solve and verify
    solve!(qcp; max_iter = 100, print_level = 1, verbose = false)

    # Test fidelity after solve
    traj = get_trajectory(qcp)
    Ũ⃗_final = traj[end][state_name(qtraj)]
    U_final = iso_vec_to_operator(Ũ⃗_final)
    fid = unitary_fidelity(U_final, U_goal)
    @test fid > 0.9

    # Test dynamics constraints are satisfied
    dynamics_integrator = qcp.prob.integrators[1]
    δ = zeros(dynamics_integrator.dim)
    DirectTrajOpt.evaluate!(δ, dynamics_integrator, traj)
    @test norm(δ, Inf) < 1e-3
end

@testitem "SmoothPulseProblem rejects spline pulses" begin
    using LinearAlgebra

    T = 10.0
    N = 50
    sys = QuantumSystem(GATES[:Z], [GATES[:X]], [1.0])
    U_goal = GATES[:X]

    times = collect(range(0.0, T, length = N))

    # LinearSplinePulse should be rejected
    pulse_linear = LinearSplinePulse(0.1 * randn(1, N), times)
    qtraj_linear = UnitaryTrajectory(sys, pulse_linear, U_goal)
    @test_throws ErrorException SmoothPulseProblem(qtraj_linear, N)

    # CubicSplinePulse should be rejected
    pulse_cubic = CubicSplinePulse(0.1 * randn(1, N), zeros(1, N), times)
    qtraj_cubic = UnitaryTrajectory(sys, pulse_cubic, U_goal)
    @test_throws ErrorException SmoothPulseProblem(qtraj_cubic, N)
end

@testitem "SmoothPulseProblem with KetTrajectory" begin
    using DirectTrajOpt
    using LinearAlgebra

    T = 10.0
    N = 50
    sys = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]

    # Create pulse and quantum trajectory
    pulse = ZeroOrderPulse(randn(2, N), collect(range(0.0, T, length = N)))
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)
    qcp = SmoothPulseProblem(qtraj, N; Q = 50.0, R = 1e-3)

    @test qcp isa QuantumControlProblem
    @test length(qcp.prob.integrators) == 3
    @test haskey(qcp.prob.trajectory.components, :du)
    @test haskey(qcp.prob.trajectory.components, :ddu)

    # Solve and verify
    solve!(qcp; max_iter = 100, print_level = 5, verbose = false)

    # Test fidelity after solve
    traj = get_trajectory(qcp)
    ψ̃_final = traj[end][state_name(qtraj)]
    ψ_final = iso_to_ket(ψ̃_final)
    fid = fidelity(ψ_final, ψ_goal)
    @test fid > 0.9

    # Test dynamics constraints are satisfied
    dynamics_integrator = qcp.prob.integrators[1]
    δ = zeros(dynamics_integrator.dim)
    DirectTrajOpt.evaluate!(δ, dynamics_integrator, traj)
    @test norm(δ, Inf) < 1e-3
end

@testitem "SmoothPulseProblem with DensityTrajectory" begin
    using DirectTrajOpt
    using LinearAlgebra

    # Simple 2-level open system with weak dissipation
    L = ComplexF64[0.1 0.0; 0.0 0.0]
    sys = OpenQuantumSystem(
        PAULIS.Z,
        [PAULIS.X, PAULIS.Y],
        [1.0, 1.0];
        dissipation_operators = [L],
    )

    ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρg = ComplexF64[0.0 0.0; 0.0 1.0]

    T = 10.0
    N = 50

    pulse = ZeroOrderPulse(0.1 * randn(2, N), collect(range(0.0, T, length = N)))
    qtraj = DensityTrajectory(sys, pulse, ρ0, ρg)

    qcp = SmoothPulseProblem(qtraj, N; Q = 100.0, R = 1e-2)

    @test qcp isa QuantumControlProblem
    @test qcp.qtraj isa DensityTrajectory

    # Check trajectory components
    @test haskey(qcp.prob.trajectory.components, :ρ⃗̃)
    @test haskey(qcp.prob.trajectory.components, :du)
    @test haskey(qcp.prob.trajectory.components, :ddu)

    # Compact iso: state dim should be n² (not 2n²)
    n = sys.levels
    @test qcp.prob.trajectory.dims[:ρ⃗̃] == n^2

    # Integrators: 1 dynamics + 2 derivatives = 3
    @test length(qcp.prob.integrators) == 3

    # Solve and verify
    solve!(qcp; max_iter = 150, print_level = 1, verbose = false)

    # Check fidelity via compact iso
    traj = get_trajectory(qcp)
    ρ̃_final = traj[end][:ρ⃗̃]
    ρ_final = compact_iso_to_density(ρ̃_final)
    fid = real(tr(ρ_final * ρg))
    @test fid > 0.9

    # Dynamics constraints should be satisfied
    dynamics_integrator = qcp.prob.integrators[1]
    δ = zeros(dynamics_integrator.dim)
    DirectTrajOpt.evaluate!(δ, dynamics_integrator, traj)
    @test norm(δ, Inf) < 1e-3
end

@testitem "SmoothPulseProblem with MultiKetTrajectory" tags = [:experimental] begin
    using DirectTrajOpt
    using LinearAlgebra

    # Use case: Optimize a single pulse that transfers multiple initial states 
    # to their respective goal states (gate-like behavior for kets)

    T = 10.0
    N = 50
    sys = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

    # Create initial and goal states for different state transfers
    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]

    # Create ensemble ket trajectory for X gate via state transfer
    # |0⟩ → |1⟩ and |1⟩ → |0⟩
    pulse = ZeroOrderPulse(randn(2, N), collect(range(0.0, T, length = N)))
    ensemble_qtraj = MultiKetTrajectory(sys, pulse, [ψ0, ψ1], [ψ1, ψ0])
    goals = ensemble_qtraj.goals
    snames = state_names(ensemble_qtraj)

    # Create problem using the new constructor
    qcp = SmoothPulseProblem(ensemble_qtraj, N; Q = 100.0, R = 1e-2)

    @test qcp isa QuantumControlProblem
    @test qcp.qtraj isa MultiKetTrajectory

    # Check trajectory components: 2 states + controls + derivatives
    @test haskey(qcp.prob.trajectory.components, :ψ̃1)
    @test haskey(qcp.prob.trajectory.components, :ψ̃2)
    @test haskey(qcp.prob.trajectory.components, :u)
    @test haskey(qcp.prob.trajectory.components, :du)
    @test haskey(qcp.prob.trajectory.components, :ddu)

    # Check integrators: 2 dynamics + 2 derivatives = 4
    @test length(qcp.prob.integrators) == 4

    # Solve and verify
    solve!(qcp; max_iter = 150, print_level = 1, verbose = true)

    # Test fidelity after solve for both states
    traj = get_trajectory(qcp)
    for (i, (name, goal)) in enumerate(zip(snames, goals))
        ψ̃_final = traj[end][name]
        ψ_final = iso_to_ket(ψ̃_final)
        fid = fidelity(ψ_final, goal)
        @test fid > 0.9
    end

    # Test dynamics constraints are satisfied for all integrators
    for integrator in qcp.prob.integrators[1:2]  # First 2 are dynamics
        δ = zeros(integrator.dim)
        DirectTrajOpt.evaluate!(δ, integrator, traj)
        @test norm(δ, Inf) < 1e-3
    end
end

# ============================================================================= #
# MultiKetTrajectory Tests (manual setup)
# ============================================================================= #

@testitem "MultiKetTrajectory manual setup" begin
    using DirectTrajOpt
    using NamedTrajectories

    # Use case: Optimize a single pulse that transfers multiple initial states 
    # to their respective goal states (gate-like behavior for kets)

    T = 1.0
    N = 50
    sys = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

    # Create initial and goal states
    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]

    # Create ensemble ket trajectory
    pulse = ZeroOrderPulse(0.1 * randn(2, N), collect(range(0.0, T, length = N)))
    ensemble_qtraj = MultiKetTrajectory(sys, pulse, [ψ0, ψ1], [ψ1, ψ0])

    @test ensemble_qtraj isa MultiKetTrajectory
    @test state_names(ensemble_qtraj) == [:ψ̃1, :ψ̃2]

    # Convert to NamedTrajectory
    new_traj = NamedTrajectory(ensemble_qtraj, N)
    snames = state_names(ensemble_qtraj)

    @test haskey(new_traj.components, :ψ̃1)
    @test haskey(new_traj.components, :ψ̃2)
    @test haskey(new_traj.components, :u)

    # Verify both initial conditions are set correctly
    @test new_traj.initial[:ψ̃1] ≈ ket_to_iso(ψ0)
    @test new_traj.initial[:ψ̃2] ≈ ket_to_iso(ψ1)

    # Build objective: coherent fidelity for gate implementation
    goals = ensemble_qtraj.goals
    J = CoherentKetInfidelityObjective(goals, snames, new_traj; Q = 50.0)

    # Add regularization
    J += QuadraticRegularizer(:u, new_traj, 1e-3)

    # Build integrators: one BilinearIntegrator per state (using dispatch)
    integrators = BilinearIntegrator(ensemble_qtraj, N)

    # Create and solve problem
    prob = DirectTrajOptProblem(new_traj, J, integrators)
    qcp = QuantumControlProblem(ensemble_qtraj, prob)

    solve!(qcp; max_iter = 50, verbose = false, print_level = 1)

    # Both state transfers should have reasonable fidelity
end

@testitem "EnsembleTrajectory vs SamplingTrajectory distinction" begin
    using DirectTrajOpt
    using NamedTrajectories

    # This test documents the key difference:
    # - MultiKetTrajectory: SAME system, DIFFERENT initial/goal states
    # - SamplingTrajectory: DIFFERENT systems, SAME goal

    T = 1.0
    N = 50
    sys = QuantumSystem(GATES[:Z], [GATES[:X]], [1.0])

    # ===== MultiKetTrajectory setup =====
    # Multiple state transfers on the SAME system
    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]

    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    ensemble_qtraj = MultiKetTrajectory(sys, pulse, [ψ0, ψ1], [ψ1, ψ0])

    @test get_system(ensemble_qtraj) === sys  # Single system
    @test length(ensemble_qtraj.initials) == 2  # Multiple state transfers
    @test state_names(ensemble_qtraj) == [:ψ̃1, :ψ̃2]

    # ===== SamplingTrajectory setup =====
    # Same goal, different systems (robust optimization)
    sys_perturbed = QuantumSystem(1.1 * GATES[:Z], [GATES[:X]], [1.0])

    pulse_unitary = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    qtraj_unitary = UnitaryTrajectory(sys, pulse_unitary, GATES[:X])

    sampling_qtraj = SamplingTrajectory(qtraj_unitary, [sys, sys_perturbed])

    @test get_system(sampling_qtraj) === sys  # Nominal system
    @test length(sampling_qtraj.systems) == 2  # Multiple systems
    @test state_names(sampling_qtraj) == [:Ũ⃗1, :Ũ⃗2]

    # Key differences:
    # 1. MultiKetTrajectory has `initials`/`goals` fields (multiple init/goals)
    # 2. SamplingTrajectory has `systems` field (multiple system params)
    # 3. State naming: numbered suffix (ψ̃1, ψ̃2 vs Ũ⃗1, Ũ⃗2)
end

# ============================================================================= #
# Time-Dependent System Tests
# ============================================================================= #

@testitem "SmoothPulseProblem with time-dependent UnitaryTrajectory" begin
    using DirectTrajOpt
    using LinearAlgebra

    # Time-dependent Hamiltonian with oscillating drive
    ω = 2π * 1.0
    H(u, t) = GATES[:Z] + u[1] * cos(ω * t) * GATES[:X] + u[2] * sin(ω * t) * GATES[:Y]

    T = 5.0
    N = 50
    sys = QuantumSystem(H, [1.0, 1.0])

    U_goal = GATES[:H]
    times = collect(range(0.0, T, length = N))
    pulse = ZeroOrderPulse(0.1 * randn(2, N), times)
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)

    qcp = SmoothPulseProblem(qtraj, N; Q = 100.0, R = 1e-2)

    @test qcp isa QuantumControlProblem

    # TimeConsistencyConstraint is auto-applied via get_trajectory_constraints
    # Should have: 1 dynamics + 2 derivatives = 3 integrators
    @test length(qcp.prob.integrators) == 3

    # Solve and verify
    solve!(qcp; max_iter = 50, print_level = 5, verbose = false)

    # Test fidelity after solve
    traj = get_trajectory(qcp)
    Ũ⃗_final = traj[end][state_name(qtraj)]
    U_final = iso_vec_to_operator(Ũ⃗_final)
    fid = unitary_fidelity(U_final, U_goal)
    @test fid > 0.85

    # Test dynamics constraints are satisfied
    dynamics_integrator = qcp.prob.integrators[1]
    δ = zeros(dynamics_integrator.dim)
    DirectTrajOpt.evaluate!(δ, dynamics_integrator, traj)
    @test norm(δ, Inf) < 1e-3
end

@testitem "SmoothPulseProblem with time-dependent KetTrajectory" begin
    using DirectTrajOpt
    using LinearAlgebra

    # Time-dependent Hamiltonian with oscillating drive
    ω = 2π * 1.0e-1
    H(u, t) = GATES[:Z] + u[1] * cos(ω * t) * GATES[:X]

    T = 10.0
    N = 50
    sys = QuantumSystem(H, [1.0])

    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

    qcp = SmoothPulseProblem(qtraj, N; Q = 50.0, R = 1e-3)

    @test qcp isa QuantumControlProblem

    # TimeConsistencyConstraint is auto-applied via get_trajectory_constraints
    # Should have: 1 dynamics + 2 derivatives = 3 integrators
    @test length(qcp.prob.integrators) == 3

    # Solve and verify
    solve!(qcp; max_iter = 100, print_level = 1, verbose = false)

    # Test fidelity after solve
    traj = get_trajectory(qcp)
    ψ̃_final = traj[end][state_name(qtraj)]
    ψ_final = iso_to_ket(ψ̃_final)
    fid = fidelity(ψ_final, ψ_goal)
    @test fid > 0.85

    # Test dynamics constraints are satisfied
    dynamics_integrator = qcp.prob.integrators[1]
    δ = zeros(dynamics_integrator.dim)
    DirectTrajOpt.evaluate!(δ, dynamics_integrator, traj)
    @test norm(δ, Inf) < 1e-3
end

@testitem "SmoothPulseProblem with SamplingTrajectory (Unitary)" begin
    using DirectTrajOpt
    using LinearAlgebra

    T = 10.0
    N = 50

    # System with uncertainty in drift
    sys_nominal = QuantumSystem(GATES[:Z], [GATES[:X]], [1.0])
    sys_perturbed = QuantumSystem(1.1 * GATES[:Z], [GATES[:X]], [1.0])

    pulse = ZeroOrderPulse(0.5 * randn(1, N), collect(range(0.0, T, length = N)))
    qtraj = UnitaryTrajectory(sys_nominal, pulse, GATES[:X])
    qcp = SmoothPulseProblem(qtraj, N; Q = 100.0, R = 1e-2)

    sampling_prob = SamplingProblem(qcp, [sys_nominal, sys_perturbed]; Q = 100.0)

    @test sampling_prob isa QuantumControlProblem
    @test sampling_prob.qtraj isa SamplingTrajectory{<:AbstractPulse,<:UnitaryTrajectory}

    # Check trajectory has sample states
    traj = get_trajectory(sampling_prob)
    @test haskey(traj.components, :Ũ⃗1)
    @test haskey(traj.components, :Ũ⃗2)

    # Solve
    solve!(sampling_prob; max_iter = 150, verbose = false, print_level = 5)

    # Test dynamics constraints are satisfied
    for integrator in sampling_prob.prob.integrators
        if integrator isa BilinearIntegrator
            δ = zeros(integrator.dim)
            DirectTrajOpt.evaluate!(δ, integrator, traj)
            @test norm(δ, Inf) < 1e-2
        end
    end
end

@testitem "SmoothPulseProblem with SamplingTrajectory (Ket)" begin
    using DirectTrajOpt
    using LinearAlgebra

    T = 10.0
    N = 50

    # System with uncertainty in drift
    sys_nominal = QuantumSystem(GATES[:Z], [GATES[:X]], [1.0])
    sys_perturbed = QuantumSystem(1.1 * GATES[:Z], [GATES[:X]], [1.0])

    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]

    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    qtraj = KetTrajectory(sys_nominal, pulse, ψ_init, ψ_goal)
    qcp = SmoothPulseProblem(qtraj, N; Q = 100.0, R = 1e-2)

    sampling_prob = SamplingProblem(qcp, [sys_nominal, sys_perturbed]; Q = 100.0)

    @test sampling_prob isa QuantumControlProblem
    @test sampling_prob.qtraj isa SamplingTrajectory{<:AbstractPulse,<:KetTrajectory}

    # Check trajectory has sample states
    traj = get_trajectory(sampling_prob)
    @test haskey(traj.components, :ψ̃1)
    @test haskey(traj.components, :ψ̃2)

    # Solve
    solve!(sampling_prob; max_iter = 50, verbose = false, print_level = 1)

    # Test dynamics constraints are satisfied
    for integrator in sampling_prob.prob.integrators
        if integrator isa BilinearIntegrator
            δ = zeros(integrator.dim)
            DirectTrajOpt.evaluate!(δ, integrator, traj)
            @test norm(δ, Inf) < 1e-3
        end
    end
end

@testitem "SmoothPulseProblem with time-dependent MultiKetTrajectory" begin
    using DirectTrajOpt
    using LinearAlgebra

    # Time-dependent Hamiltonian with oscillating drive
    ω = 2π * 5.0
    H(u, t) = GATES[:Z] + u[1] * cos(ω * t) * GATES[:X] + u[2] * sin(ω * t) * GATES[:Y]

    T = 10.0
    N = 50
    sys = QuantumSystem(H, [1.0, 1.0])

    # Create ensemble ket trajectory for X gate
    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]

    pulse = ZeroOrderPulse(0.1 * randn(2, N), collect(range(0.0, T, length = N)))
    ensemble_qtraj = MultiKetTrajectory(sys, pulse, [ψ0, ψ1], [ψ1, ψ0])

    qcp = SmoothPulseProblem(ensemble_qtraj, N; Q = 50.0, R = 1e-3)

    @test qcp isa QuantumControlProblem
    @test qcp.qtraj isa MultiKetTrajectory

    # TimeConsistencyConstraint is auto-applied via get_trajectory_constraints
    # Should have: 2 dynamics + 2 derivatives = 4 integrators
    @test length(qcp.prob.integrators) == 4

    # Solve and verify
    solve!(qcp; max_iter = 50, print_level = 1, verbose = false)

    # Test fidelity for both states after solve
    traj = get_trajectory(qcp)
    snames = state_names(ensemble_qtraj)
    goals = ensemble_qtraj.goals

    for (name, goal) in zip(snames, goals)
        ψ̃_final = traj[end][name]
        ψ_final = iso_to_ket(ψ̃_final)
        fid = fidelity(ψ_final, goal)
        @test fid > 0.80
    end
end

@testitem "SmoothPulseProblem with time-dependent SamplingTrajectory (Unitary)" begin
    using DirectTrajOpt
    using LinearAlgebra

    # Time-dependent Hamiltonian with oscillating drive
    ω = 2π * 5.0
    H1(u, t) = GATES[:Z] + u[1] * cos(ω * t) * GATES[:X]
    H2(u, t) = 1.1 * GATES[:Z] + u[1] * cos(ω * t) * GATES[:X]  # Perturbed

    T = 1.0
    N = 50
    sys_nominal = QuantumSystem(H1, [1.0])
    sys_perturbed = QuantumSystem(H2, [1.0])

    U_goal = GATES[:X]
    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    qtraj = UnitaryTrajectory(sys_nominal, pulse, U_goal)

    qcp = SmoothPulseProblem(qtraj, N; Q = 100.0, R = 1e-2)

    # Create sampling problem
    sampling_prob = SamplingProblem(qcp, [sys_nominal, sys_perturbed]; Q = 100.0)

    @test sampling_prob isa QuantumControlProblem
    @test sampling_prob.qtraj isa SamplingTrajectory{<:AbstractPulse,<:UnitaryTrajectory}

    # Check trajectory has sample states
    traj = get_trajectory(sampling_prob)
    @test haskey(traj.components, :Ũ⃗1)
    @test haskey(traj.components, :Ũ⃗2)

    # TimeConsistencyConstraint is auto-applied
    # Integrators: 2 dynamics (samples) + 2 derivatives = 4
    # (depending on SamplingProblem implementation)

    # Solve
    solve!(sampling_prob; max_iter = 50, verbose = false, print_level = 1)

    # Test dynamics constraints are satisfied
    for integrator in sampling_prob.prob.integrators
        if integrator isa BilinearIntegrator
            δ = zeros(integrator.dim)
            DirectTrajOpt.evaluate!(δ, integrator, traj)
            @test norm(δ, Inf) < 1e-2
        end
    end
end

@testitem "SmoothPulseProblem with global_names requires custom integrator" begin

    # System with global parameters
    T = 2.0
    N = 10

    H = (u, t) -> begin
        δ = u[2]  # Global detuning
        δ * GATES.Z + u[1] * GATES.X
    end

    δ_init = 0.1
    sys = QuantumSystem(H, [1.0]; time_dependent = true, global_params = (δ = δ_init,))
    U_goal = GATES.X

    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)

    # Should error when global_names provided without custom integrator
    @test_throws ErrorException SmoothPulseProblem(
        qtraj,
        N;
        Q = 100.0,
        R = 1e-2,
        global_names = [:δ],
    )
end

@testitem "SmoothPulseProblem with global_bounds error handling" begin
    using NamedTrajectories
    using DirectTrajOpt

    # Test that global_bounds throws an informative error when global doesn't exist
    # (global_data must come from integrator - e.g., HermitianExponentialIntegrator from Piccolissimo)

    T = 5.0
    N = 10

    sys = QuantumSystem(0.1 * GATES.Z, [GATES.X], [1.0])
    U_goal = GATES.X

    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)

    # Attempting to use global_bounds without globals in trajectory should error
    @test_throws "Global variable :δ not found" SmoothPulseProblem(
        qtraj,
        N;
        Q = 100.0,
        R = 1e-2,
        global_bounds = Dict(:δ => 0.5),  # δ doesn't exist in trajectory
    )
end
