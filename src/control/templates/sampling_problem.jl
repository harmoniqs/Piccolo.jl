export SamplingProblem

# Note: SamplingTrajectory is now exported from PiccoloQuantumObjects

# ============================================================================= #
# Sampling Objective Utilities
# ============================================================================= #

"""
    extract_regularization(objective, state_sym::Symbol, new_traj::NamedTrajectory) -> AbstractObjective

Extract regularization terms (non-state-dependent objectives) from a composite objective,
filtering to only include terms for variables that exist in the new trajectory.

Used by `SamplingProblem` to extract shared regularizers (e.g., control penalty) from
the base problem while excluding regularizers for variables that don't exist in the
sampling trajectory (e.g., `:du`, `:ddu` which are added by `SmoothPulseProblem`).
"""
function extract_regularization(objective, state_sym::Symbol, new_traj::NamedTrajectory)
    objs = hasproperty(objective, :objectives) ? objective.objectives : [objective]

    regularizers = filter(objs) do term
        # Get variable names this term depends on
        term_syms = if hasproperty(term, :syms)
            term.syms
        elseif hasproperty(term, :var_names)
            term.var_names
        elseif hasproperty(term, :name) && term.name isa Symbol
            # QuadraticRegularizer has a single :name field
            (term.name,)
        else
            Symbol[]
        end
        # Only include if:
        # 1. It doesn't depend on the state symbol
        # 2. All its variables exist in the new trajectory
        state_sym ∉ term_syms && all(s -> s ∈ new_traj.names, term_syms)
    end

    return isempty(regularizers) ? NullObjective() : reduce(+, regularizers)
end

# ============================================================================= #
# Sampling State Objective (dispatch-based)
# ============================================================================= #

"""
    sampling_state_objective(qtraj, traj, state_sym, Q)

Create the state-dependent objective for a sampling member.
Dispatches on quantum trajectory type.
"""
function sampling_state_objective(
    qtraj::UnitaryTrajectory,
    traj::NamedTrajectory,
    state_sym::Symbol,
    Q::Float64,
)
    goal = get_goal(qtraj)
    # Use free-phase objective when globals exist and goal is EmbeddedOperator
    if traj.global_dim > 0 && goal isa EmbeddedOperator
        # `Symbol[traj.global_names...]` rather than `collect(traj.global_names)` so
        # the eltype is `Symbol` even when the source tuple is empty (JET infers
        # `Vector{Union{}}` from `collect` on a possibly-empty tuple, which then
        # fails to match `θ_names::AbstractVector{Symbol}` downstream).
        θ_names = Symbol[traj.global_names...]
        U_goal_fn = _make_free_phase_goal(goal)
        return UnitaryFreePhaseInfidelityObjective(
            U_goal_fn,
            state_sym,
            θ_names,
            traj;
            Q = Q,
        )
    else
        return UnitaryInfidelityObjective(goal, state_sym, traj; Q = Q)
    end
end

function sampling_state_objective(
    qtraj::KetTrajectory,
    traj::NamedTrajectory,
    state_sym::Symbol,
    Q::Float64,
)
    ψ_goal = get_goal(qtraj)
    return KetInfidelityObjective(ψ_goal, state_sym, traj; Q = Q)
end

function sampling_state_objective(
    qtraj::DensityTrajectory,
    traj::NamedTrajectory,
    state_sym::Symbol,
    Q::Float64,
)
    # DensityTrajectory doesn't have a fidelity objective yet
    # Return NullObjective for now
    return NullObjective(traj)
end

# ============================================================================= #
# SamplingProblem Constructor
# ============================================================================= #

@doc raw"""
    SamplingProblem(qcp::QuantumControlProblem, systems::Vector{<:AbstractQuantumSystem}; kwargs...)

Construct a `SamplingProblem` from an existing `QuantumControlProblem` and a list of systems.

This creates a robust optimization problem where the controls are shared across all systems,
but each system evolves according to its own dynamics. The objective is the weighted sum of
fidelity objectives for each system.

# Arguments
- `qcp::QuantumControlProblem`: The base problem (defines nominal trajectory, objective, etc.)
- `systems::Vector{<:AbstractQuantumSystem}`: List of systems to optimize over

# Keyword Arguments
- `weights::Vector{Float64}=fill(1.0, length(systems))`: Weights for each system
- `Q::Float64=100.0`: Weight on infidelity objective (explicit, not extracted from base problem)
- `integrator::Union{Nothing, Function}=nothing`: Optional integrator factory function. When
  provided, it is called as `integrator(sampling_qtraj, N)` and must return an integrator or
  vector of integrators. When `nothing` (default), `BilinearIntegrator` is used.
- `calibration_targets::Vector{Symbol}=Symbol[]`: Names of globals declared as **calibration targets** — knobs an external calibration step manages, not free NLP variables. SamplingProblem builds a fresh constraint list (rather than inheriting from the base `qcp`), so calibration_target pins set on the base `qcp` are *not* automatically carried over — pass them here explicitly. Default empty: globals stay free.
- `piccolo_options::PiccoloOptions=PiccoloOptions()`: Options for the solver

# Returns
- `QuantumControlProblem{SamplingTrajectory}`: A new problem with the sampling trajectory
"""
function SamplingProblem(
    qcp::QuantumControlProblem,
    systems::Vector{<:AbstractQuantumSystem};
    weights::Vector{Float64} = fill(1.0, length(systems)),
    Q::Float64 = 100.0,
    integrator::Union{Nothing,Function} = nothing,
    calibration_targets::Vector{Symbol} = Symbol[],
    piccolo_options::PiccoloOptions = PiccoloOptions(),
)
    if _show_header(piccolo_options)
        println("constructing SamplingProblem")
        println("    systems: $(length(systems))")
    end

    base_qtraj = qcp.qtraj
    state_sym = state_name(base_qtraj)
    base_traj = get_trajectory(qcp)

    # 1. Create SamplingTrajectory wrapper (new API: no stored trajectory)
    sampling_qtraj = SamplingTrajectory(base_qtraj, systems; weights)

    # 2. Build trajectory from sampling trajectory (this creates duplicated states)
    #    Propagate Δt bounds and global variables from base problem
    N = base_traj.N
    Δt_bounds = if haskey(base_traj.bounds, :Δt)
        (base_traj.bounds[:Δt][1][1], base_traj.bounds[:Δt][2][1])
    else
        nothing
    end

    new_traj = NamedTrajectory(sampling_qtraj, N; Δt_bounds = Δt_bounds)

    # Propagate global variables (e.g., free-phase φ_1, φ_2) from base trajectory.
    # Directly mutate the struct fields — NamedTrajectory doesn't have a public API
    # for adding globals after construction.
    if base_traj.global_dim > 0
        new_traj = NamedTrajectory(
            new_traj.datavec,
            new_traj.components,
            new_traj.N;
            timestep = new_traj.timestep,
            controls = new_traj.control_names,
            bounds = new_traj.bounds,
            initial = new_traj.initial,
            final = isnothing(new_traj.final_) ? NamedTuple() : new_traj.final_,
            goal = new_traj.goal,
            global_data = copy(base_traj.global_data),
            global_components = base_traj.global_components,
        )
    end
    snames = state_names(sampling_qtraj)

    # 3. Build objective: weighted state objectives + shared regularization
    J_state = sum(
        sampling_state_objective(base_qtraj, new_traj, name, w * Q) for
        (name, w) in zip(snames, weights)
    )
    J_reg = extract_regularization(qcp.prob.objective, state_sym, new_traj)
    J_total = J_state + J_reg

    # 4. Build integrators: dynamics for each system
    #    Note: We don't carry over DerivativeIntegrators from the base problem
    #    because they operate on :du, :ddu which don't exist in the sampling trajectory.
    #    For now, SamplingProblem operates on the raw controls without derivative smoothing.
    #    TODO: Consider adding an option to preserve smoothness constraints.

    # Use BilinearIntegrator by default, or a custom integrator factory via the
    # `integrator` kwarg (must accept (sampling_qtraj, N) and return integrator(s)).
    if isnothing(integrator)
        dynamics_integrators = BilinearIntegrator(sampling_qtraj, N)
    else
        dynamics_integrators = integrator(sampling_qtraj, N)
    end

    all_integrators =
        dynamics_integrators isa AbstractVector ? dynamics_integrators :
        [dynamics_integrators]

    # 5. Construct problem (TimeConsistencyConstraint auto-applied)
    constraints = AbstractConstraint[]

    # Pin calibration targets at nominal — knobs the user has declared an
    # external calibration step will manage, not free NLP variables. SamplingProblem
    # builds a fresh constraint list (cf. base qcp's), so this is not inherited
    # automatically; the user passes the list explicitly.
    apply_calibration_targets!(
        constraints,
        calibration_targets,
        new_traj;
        verbose = _show_details(piccolo_options),
    )

    prob =
        DirectTrajOptProblem(new_traj, J_total, all_integrators; constraints = constraints)

    return _maybe_display(QuantumControlProblem(sampling_qtraj, prob), piccolo_options)
end

# ============================================================================= #
# Composability with MinimumTimeProblem
# ============================================================================= #

function _update_goal(qtraj::SamplingTrajectory, new_goal)
    new_base = _update_goal(qtraj.base_trajectory, new_goal)
    # Reconstruct the SamplingTrajectory around the updated base — a thin
    # SamplingTrajectory only owns (base_trajectory, systems, weights).
    return SamplingTrajectory(new_base, qtraj.systems; weights = qtraj.weights)
end

function _final_fidelity_constraint(
    qtraj::SamplingTrajectory,
    final_fidelity::Float64,
    traj::NamedTrajectory;
    subsystem_levels::Union{Nothing,Vector{Int}} = nothing,
)
    constraints = [
        _sampling_fidelity_constraint(qtraj.base_trajectory, name, final_fidelity, traj) for name in state_names(qtraj)
    ]
    return constraints
end

# Dispatch on trajectory type for fidelity constraint
function _sampling_fidelity_constraint(
    qtraj::UnitaryTrajectory,
    state_sym::Symbol,
    final_fidelity::Float64,
    traj::NamedTrajectory,
)
    return FinalUnitaryFidelityConstraint(qtraj.goal, state_sym, final_fidelity, traj)
end

function _sampling_fidelity_constraint(
    qtraj::KetTrajectory,
    state_sym::Symbol,
    final_fidelity::Float64,
    traj::NamedTrajectory,
)
    return FinalKetFidelityConstraint(qtraj.goal, state_sym, final_fidelity, traj)
end

# Tests
@testitem "SamplingProblem Construction" begin
    using DirectTrajOpt

    T = 10.0
    N = 50

    # Define system
    sys = QuantumSystem(GATES[:Z], [GATES[:X]], [1.0])

    # Create base problem
    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    qtraj = UnitaryTrajectory(sys, pulse, GATES[:H])
    qcp = SmoothPulseProblem(qtraj, N; Q = 100.0)

    # Create sampling problem with 2 systems
    systems = [sys, sys] # Identical systems for testing
    sampling_prob = SamplingProblem(qcp, systems)

    @test sampling_prob isa QuantumControlProblem
    @test sampling_prob.qtraj isa SamplingTrajectory
    @test length(sampling_prob.qtraj.systems) == 2

    # Check trajectory components (now use numbered suffix :Ũ⃗1, :Ũ⃗2, etc.)
    traj = get_trajectory(sampling_prob)
    @test haskey(traj.components, :Ũ⃗1)
    @test haskey(traj.components, :Ũ⃗2)
    @test haskey(traj.components, :u)

    # Check integrators
    # Should have 2 dynamics integrators (one per system)
    # SamplingProblem doesn't carry derivative integrators from base problem
    @test length(sampling_prob.prob.integrators) == 2
end

@testitem "SamplingProblem Solving" tags = [:sampling_problem] begin
    using DirectTrajOpt

    T = 10.0
    N = 50

    # Simple robust optimization
    # System with uncertainty in drift
    sys_nominal = QuantumSystem(GATES[:Z], [GATES[:X]], [1.0])
    sys_perturbed = QuantumSystem(1.1 * GATES[:Z], [GATES[:X]], [1.0])

    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    qtraj = UnitaryTrajectory(sys_nominal, pulse, GATES[:X])
    qcp = SmoothPulseProblem(qtraj, N; Q = 100.0)

    sampling_prob = SamplingProblem(qcp, [sys_nominal, sys_perturbed])

    # Solve
    solve!(sampling_prob; max_iter = 10, verbose = false)

    # Check that we have a solution
    @test sampling_prob.prob.objective(sampling_prob.trajectory) < 1e10 # Just check it didn't blow up
end

@testitem "SamplingProblem with KetTrajectory" begin
    using DirectTrajOpt

    T = 1.0
    N = 50

    # Robust state transfer over parameter uncertainty
    sys_nominal = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    sys_perturbed = QuantumSystem(1.1 * GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    pulse = ZeroOrderPulse(0.1 * randn(2, N), collect(range(0.0, T, length = N)))
    qtraj = KetTrajectory(sys_nominal, pulse, ψ_init, ψ_goal)

    qcp = SmoothPulseProblem(qtraj, N; Q = 50.0, R = 1e-3)

    # Create sampling problem
    sampling_prob = SamplingProblem(qcp, [sys_nominal, sys_perturbed]; Q = 50.0)

    @test sampling_prob isa QuantumControlProblem
    @test sampling_prob.qtraj isa SamplingTrajectory

    # Check trajectory has sample states (now use numbered suffix)
    traj = get_trajectory(sampling_prob)
    @test haskey(traj.components, :ψ̃1)
    @test haskey(traj.components, :ψ̃2)

    # Solve
    solve!(sampling_prob; max_iter = 10, verbose = false, print_level = 1)
end

@testitem "SamplingProblem with custom weights" begin
    using DirectTrajOpt

    T = 10.0
    N = 50

    sys1 = QuantumSystem(GATES[:Z], [GATES[:X]], [1.0])
    sys2 = QuantumSystem(1.1 * GATES[:Z], [GATES[:X]], [1.0])
    sys3 = QuantumSystem(0.9 * GATES[:Z], [GATES[:X]], [1.0])

    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    qtraj = UnitaryTrajectory(sys1, pulse, GATES[:X])
    qcp = SmoothPulseProblem(qtraj, N; Q = 100.0)

    # Non-uniform weights - emphasize nominal system
    weights = [0.6, 0.2, 0.2]
    sampling_prob = SamplingProblem(qcp, [sys1, sys2, sys3]; weights = weights, Q = 100.0)

    @test sampling_prob.qtraj.weights == weights
    @test length(sampling_prob.qtraj.systems) == 3

    # Should have 3 sample states (numbered suffix)
    traj = get_trajectory(sampling_prob)
    @test haskey(traj.components, :Ũ⃗1)
    @test haskey(traj.components, :Ũ⃗2)
    @test haskey(traj.components, :Ũ⃗3)

    solve!(sampling_prob; max_iter = 5, verbose = false, print_level = 1)
end

@testitem "SamplingProblem + MinimumTimeProblem composition" begin
    using DirectTrajOpt

    T = 1.0
    N = 50

    # Robust minimum-time optimization
    sys_nominal = QuantumSystem(0.1 * GATES[:Z], [GATES[:X]], [1.0])
    sys_perturbed = QuantumSystem(0.11 * GATES[:Z], [GATES[:X]], [1.0])

    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    qtraj = UnitaryTrajectory(sys_nominal, pulse, GATES[:X])
    qcp = SmoothPulseProblem(qtraj, N; Q = 100.0, R = 1e-2, Δt_bounds = (0.01, 0.5))

    # First create sampling problem
    sampling_prob = SamplingProblem(qcp, [sys_nominal, sys_perturbed]; Q = 100.0)
    solve!(sampling_prob; max_iter = 20, verbose = false, print_level = 1)

    # Then convert to minimum-time
    mintime_prob = MinimumTimeProblem(sampling_prob; final_fidelity = 0.90, D = 50.0)

    @test mintime_prob isa QuantumControlProblem
    @test mintime_prob.qtraj isa SamplingTrajectory

    # Solve minimum-time
    solve!(mintime_prob; max_iter = 20, verbose = false, print_level = 1)
end

@testitem "SamplingProblem with EmbeddedOperator" begin
    using DirectTrajOpt

    # Minimal setup (reproducing the bug from main.jl)
    T = 1.0
    N = 10
    times = collect(range(0, T, length = N))
    initial_controls = zeros(2, N)
    pulse = ZeroOrderPulse(initial_controls, times)

    # Create systems
    sys1 = TransmonSystem(levels = 2, δ = -0.18)
    sys2 = TransmonSystem(levels = 2, δ = -0.20)

    # Create embedded operator (this is what caused the bug)
    X_embedded = EmbeddedOperator(GATES[:X], sys2)

    # Create trajectory with embedded operator as goal
    qtraj = UnitaryTrajectory(sys2, pulse, X_embedded)
    qcp = SmoothPulseProblem(qtraj, N; Q = 100.0)

    # This should not fail - it's the bug we're fixing
    sampling_prob = SamplingProblem(qcp, [sys1, sys2])

    @test sampling_prob isa QuantumControlProblem
    @test sampling_prob.qtraj isa SamplingTrajectory
    @test length(sampling_prob.qtraj.systems) == 2
end

@testitem "SamplingProblem with DensityTrajectory" tags = [:density, :skip] begin
    # TODO: DensityTrajectory support for SamplingProblem is not yet complete
    # Needs: BilinearIntegrator dispatch, SamplingTrajectory NamedTrajectory conversion
    @test_skip "DensityTrajectory support not yet implemented"
end

@testitem "SamplingProblem with custom integrator factory" begin
    using DirectTrajOpt
    using LinearAlgebra

    T = 10.0
    N = 50

    sys_nominal = QuantumSystem(GATES[:Z], [GATES[:X]], [1.0])
    sys_perturbed = QuantumSystem(1.1 * GATES[:Z], [GATES[:X]], [1.0])

    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    qtraj = UnitaryTrajectory(sys_nominal, pulse, GATES[:X])
    qcp = SmoothPulseProblem(qtraj, N; Q = 100.0)

    # Custom integrator factory — reimplements default BilinearIntegrator logic
    custom_factory(sqtraj, n) = BilinearIntegrator(sqtraj, n)

    sampling_prob =
        SamplingProblem(qcp, [sys_nominal, sys_perturbed]; integrator = custom_factory)

    @test sampling_prob isa QuantumControlProblem
    @test length(sampling_prob.prob.integrators) == 2

    solve!(sampling_prob; max_iter = 5, verbose = false, print_level = 1)
end
