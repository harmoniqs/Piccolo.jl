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
    return UnitaryInfidelityObjective(get_goal(qtraj), state_sym, traj; Q = Q)
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
- `piccolo_options::PiccoloOptions=PiccoloOptions()`: Options for the solver

# Returns
- `QuantumControlProblem{SamplingTrajectory}`: A new problem with the sampling trajectory
"""
function SamplingProblem(
    qcp::QuantumControlProblem,
    systems::Vector{<:AbstractQuantumSystem};
    weights::Vector{Float64} = fill(1.0, length(systems)),
    Q::Float64 = 100.0,
    piccolo_options::PiccoloOptions = PiccoloOptions(),
)
    if piccolo_options.verbose
        println("    constructing SamplingProblem...")
        println("\tusing $(length(systems)) systems")
    end

    base_qtraj = qcp.qtraj
    state_sym = state_name(base_qtraj)
    base_traj = get_trajectory(qcp)

    # 1. Create SamplingTrajectory wrapper (new API: no stored trajectory)
    sampling_qtraj = SamplingTrajectory(base_qtraj, systems; weights)

    # 2. Build trajectory from sampling trajectory (this creates duplicated states)
    #    Propagate Δt bounds from base problem if they exist
    N = base_traj.N
    Δt_bounds = if haskey(base_traj.bounds, :Δt)
        (base_traj.bounds[:Δt][1][1], base_traj.bounds[:Δt][2][1])
    else
        nothing
    end
    new_traj = NamedTrajectory(sampling_qtraj, N; Δt_bounds = Δt_bounds)
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

    # Use BilinearIntegrator dispatch on SamplingTrajectory
    dynamics_integrators = BilinearIntegrator(sampling_qtraj, N)

    all_integrators =
        dynamics_integrators isa AbstractVector ? dynamics_integrators :
        [dynamics_integrators]

    # 5. Construct problem (TimeConsistencyConstraint auto-applied)
    prob = DirectTrajOptProblem(
        new_traj,
        J_total,
        all_integrators;
        constraints = AbstractConstraint[],
    )

    return QuantumControlProblem(sampling_qtraj, prob)
end

# ============================================================================= #
# Composability with MinimumTimeProblem
# ============================================================================= #

function _update_goal(qtraj::SamplingTrajectory, new_goal)
    new_base = _update_goal(qtraj.base_trajectory, new_goal)
    return update_base_trajectory(qtraj, new_base)
end

function _final_fidelity_constraint(
    qtraj::SamplingTrajectory,
    final_fidelity::Float64,
    traj::NamedTrajectory,
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

@testitem "SamplingProblem with DensityTrajectory" tags = [:density, :skip] begin
    # TODO: DensityTrajectory support for SamplingProblem is not yet complete
    # Needs: BilinearIntegrator dispatch, SamplingTrajectory NamedTrajectory conversion
    @test_skip "DensityTrajectory support not yet implemented"
end
