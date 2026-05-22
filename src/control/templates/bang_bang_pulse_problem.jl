export BangBangPulseProblem

@doc raw"""
    BangBangPulseProblem(qtraj::AbstractQuantumTrajectory{<:ZeroOrderPulse}, N::Int; kwargs...)

Construct a `QuantumControlProblem` that promotes **bang-bang** (piecewise-constant,
few-switch) controls by penalizing ``\|du\|_1`` via an **exact slack reformulation**.

Unlike `SmoothPulseProblem` (which uses 2 derivative levels with L2 regularization),
this stores only 1 derivative (`du`) and uses slack variables to impose an exact L1
penalty on it, promoting sparsity in `du` and thus fewer switches.

# L1 penalty via slack variables

Introduces slack variables ``s \geq 0`` (same dimension as `du`) and enforces:
```math
|du_{k,i}| \leq s_{k,i}
```
Then minimizes the linear cost ``R_{du} \sum_k \sum_i s_{k,i} \Delta t_k``.
At optimality, ``s = |du|``, giving the exact L1 norm.

# Arguments
- `qtraj::AbstractQuantumTrajectory{<:ZeroOrderPulse}`: Quantum trajectory with piecewise constant pulse
- `N::Int`: Number of timesteps for discretization

# Keyword Arguments
- `integrator::Union{Nothing, AbstractIntegrator, Vector{<:AbstractIntegrator}}=nothing`: Optional custom integrator(s). If not provided, uses BilinearIntegrator.
- `global_names::Union{Nothing, Vector{Symbol}}=nothing`: Names of global variables to optimize. Requires a custom integrator.
- `global_bounds::Union{Nothing, Dict{Symbol, Union{Float64, Tuple{Float64, Float64}}}}=nothing`: Bounds for global variables.
- `calibration_targets::Vector{Symbol}=Symbol[]`: Names of globals declared as **calibration targets** — knobs an external calibration step manages, not free NLP variables. Each listed name is pinned at its nominal value via `GlobalEqualityConstraint` so the QCP solve cannot drift it as a slack variable. Default empty: globals stay free.
- `du_bound::Float64=Inf`: Bound on discrete first derivative
- `Δt_bounds::Union{Nothing, Tuple{Float64, Float64}}=nothing`: Timestep bounds
- `Q::Float64=100.0`: Weight on infidelity/objective
- `R::Float64=1e-2`: Default regularization weight
- `R_u::Union{Float64, Vector{Float64}}=0.0`: L2 weight on control amplitude (defaults to 0 — bang-bang needs no amplitude regularization)
- `R_du::Union{Float64, Vector{Float64}}=R`: L1 weight on `du` (applied to slacks)
- `constraints::Vector{<:AbstractConstraint}=AbstractConstraint[]`: Additional constraints
- `piccolo_options::PiccoloOptions=PiccoloOptions()`: Piccolo solver options
- `free_phase::Bool=false`: Optimize per-qubit Z-phase rotations alongside the pulse. For `UnitaryTrajectory`, requires an `EmbeddedOperator` goal. For `KetTrajectory`, requires `subsystem_levels`.
- `subsystem_levels::Union{Nothing, Vector{Int}}=nothing`: Number of levels per subsystem (required for `free_phase=true` with `KetTrajectory` and `MultiKetTrajectory`).
- `initial_phases::Union{Nothing, Vector{Float64}}=nothing`: Initial values for per-qubit phase variables.

# Returns
- `QuantumControlProblem`: Wrapper containing quantum trajectory and optimization problem

# Examples
```julia
sys = QuantumSystem(H_drift, H_drives, drive_bounds)
pulse = ZeroOrderPulse(0.1 * randn(n_drives, N), collect(range(0.0, T, length=N)))
qtraj = UnitaryTrajectory(sys, pulse, U_goal)
qcp = BangBangPulseProblem(qtraj, N; Q=100.0, R_du=1e-1)
solve!(qcp; max_iter=200)
```

See also: [`SmoothPulseProblem`](@ref) for smooth (L2-regularized) controls.
"""
function BangBangPulseProblem(
    qtraj::AbstractQuantumTrajectory{<:ZeroOrderPulse},
    N::Int;
    integrator::Union{Nothing,AbstractIntegrator,Vector{<:AbstractIntegrator}} = nothing,
    global_names::Union{Nothing,Vector{Symbol}} = nothing,
    global_bounds::Union{Nothing,Dict{Symbol,<:Union{Float64,Tuple{Float64,Float64}}}} = nothing,
    calibration_targets::Vector{Symbol} = Symbol[],
    du_bound::Float64 = Inf,
    Δt_bounds::Union{Nothing,Tuple{Float64,Float64}} = nothing,
    Q::Float64 = 100.0,
    R::Float64 = 1e-2,
    R_u::Union{Float64,Vector{Float64}} = 0.0,
    R_du::Union{Float64,Vector{Float64}} = R,
    constraints::Vector{<:AbstractConstraint} = AbstractConstraint[],
    piccolo_options::PiccoloOptions = PiccoloOptions(),
    free_phase::Bool = false,
    subsystem_levels::Union{Nothing,Vector{Int}} = nothing,
    initial_phases::Union{Nothing,Vector{Float64}} = nothing,
)
    if _show_header(piccolo_options)
        println("constructing BangBangPulseProblem [$(_typename(qtraj))]")
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

    # Free-phase support: add per-qubit phase variables as globals
    θ_names = Symbol[]
    U_goal_fn = nothing
    ket_goal_fn = nothing
    if free_phase
        if qtraj isa KetTrajectory
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
    base_traj = if isnothing(global_data)
        NamedTrajectory(qtraj, N; Δt_bounds = Δt_bounds)
    else
        NamedTrajectory(qtraj, N; Δt_bounds = Δt_bounds, global_data = global_data)
    end

    # Add 1 control derivative (not 2 like SmoothPulseProblem)
    du_bounds = fill(du_bound, sys.n_drives)

    traj_bb = add_control_derivatives(
        base_traj,
        1;
        control_name = control_sym,
        derivative_bounds = (du_bounds,),
    )

    # Get control derivative names
    control_names =
        [name for name ∈ traj_bb.names if endswith(string(name), string(control_sym))]

    du_sym = control_names[2]  # e.g. :du
    s_du_sym = Symbol(:s_, du_sym)  # e.g. :s_du

    # Add slack variable: s_du ≥ 0, initialized at |du|
    n_drives = sys.n_drives
    s_du_data = abs.(traj_bb[du_sym])
    traj_bb = add_component(
        traj_bb,
        s_du_sym,
        s_du_data;
        type = :control,
        bounds = merge(
            traj_bb.bounds,
            NamedTuple{(s_du_sym,)}(((zeros(n_drives), fill(Inf, n_drives)),)),
        ),
    )

    # Initialize dynamics integrators
    # Combine user-specified global_names with free_phase θ_names for integrator
    all_global_names = if !isempty(θ_names)
        gn = isnothing(global_names) ? Symbol[] : copy(global_names)
        append!(gn, θ_names)
        gn
    else
        global_names
    end

    if isnothing(integrator)
        if !isnothing(all_global_names) && !isempty(all_global_names)
            error(
                "free_phase=true or global_names requires a custom integrator that supports global variables. " *
                "Use HermitianExponentialIntegrator from Piccolissimo:\n" *
                "  using Piccolissimo\n" *
                "  integrator = HermitianExponentialIntegrator(qtraj, N; global_names=$all_global_names)\n" *
                "  qcp = BangBangPulseProblem(qtraj, N; integrator=integrator, ...)",
            )
        end
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

    # Build objective: type-specific infidelity + regularization
    J = if free_phase && !isnothing(ket_goal_fn)
        KetFreePhaseInfidelityObjective(ket_goal_fn, state_sym, θ_names, traj_bb; Q = Q)
    elseif free_phase && !isnothing(U_goal_fn)
        UnitaryFreePhaseInfidelityObjective(U_goal_fn, state_sym, θ_names, traj_bb; Q = Q)
    else
        _state_objective(qtraj, traj_bb, state_sym, Q)
    end

    # L2 on control amplitude
    J += QuadraticRegularizer(control_names[1], traj_bb, R_u)

    # L1 on du via slack variables
    J += LinearRegularizer(s_du_sym, traj_bb, R_du)

    # Add optional Piccolo constraints and objectives
    J += _apply_piccolo_options(qtraj, piccolo_options, constraints, traj_bb, state_sym)

    # Build integrators: dynamics + 1 derivative integrator (not 2)
    integrators = copy(dynamics_integrators)
    push!(integrators, DerivativeIntegrator(control_names[1], control_names[2], traj_bb))

    # Build constraints: L1 slack constraint ties |du| ≤ s_du
    all_constraints = copy(constraints)
    push!(all_constraints, L1SlackConstraint(du_sym, s_du_sym, traj_bb))

    add_global_bounds_constraints!(
        all_constraints,
        global_bounds,
        traj_bb;
        verbose = _show_details(piccolo_options),
    )

    apply_calibration_targets!(
        all_constraints,
        calibration_targets,
        traj_bb;
        verbose = _show_details(piccolo_options),
    )

    prob = DirectTrajOptProblem(traj_bb, J, integrators; constraints = all_constraints)

    return _maybe_display(QuantumControlProblem(qtraj, prob), piccolo_options)
end

# ============================================================================= #
# MultiKetTrajectory Constructor
# ============================================================================= #

@doc raw"""
    BangBangPulseProblem(qtraj::MultiKetTrajectory{<:ZeroOrderPulse}, N::Int; kwargs...)

Construct a `QuantumControlProblem` for bang-bang pulse optimization over an ensemble
of ket state transfers with piecewise constant controls.

See the single-trajectory method for full documentation of keyword arguments.
"""
function BangBangPulseProblem(
    qtraj::MultiKetTrajectory{<:ZeroOrderPulse},
    N::Int;
    integrator::Union{Nothing,AbstractIntegrator,Vector{<:AbstractIntegrator}} = nothing,
    global_names::Union{Nothing,Vector{Symbol}} = nothing,
    global_bounds::Union{Nothing,Dict{Symbol,<:Union{Float64,Tuple{Float64,Float64}}}} = nothing,
    calibration_targets::Vector{Symbol} = Symbol[],
    du_bound::Float64 = Inf,
    Δt_bounds::Union{Nothing,Tuple{Float64,Float64}} = nothing,
    Q::Float64 = 100.0,
    R::Float64 = 1e-2,
    R_u::Union{Float64,Vector{Float64}} = 0.0,
    R_du::Union{Float64,Vector{Float64}} = R,
    constraints::Vector{<:AbstractConstraint} = AbstractConstraint[],
    piccolo_options::PiccoloOptions = PiccoloOptions(),
    free_phase::Bool = false,
    subsystem_levels::Union{Nothing,Vector{Int}} = nothing,
    initial_phases::Union{Nothing,Vector{Float64}} = nothing,
    coherent::Bool = true,
)
    if _show_header(piccolo_options)
        println(
            "    constructing BangBangPulseProblem for $(nameof(typeof(qtraj))) ($(length(qtraj.initials)) states)...",
        )
    end

    # Extract info from ensemble trajectory
    sys = get_system(qtraj)
    control_sym = drive_name(qtraj)
    snames = state_names(qtraj)
    weights = qtraj.weights
    goals = qtraj.goals

    # Build global_data from system's global_params if present
    global_data = if hasproperty(sys, :global_params) && !isempty(sys.global_params)
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
    base_traj = if isnothing(global_data)
        NamedTrajectory(qtraj, N; Δt_bounds = Δt_bounds)
    else
        NamedTrajectory(qtraj, N; Δt_bounds = Δt_bounds, global_data = global_data)
    end

    # Add 1 control derivative
    du_bounds = fill(du_bound, sys.n_drives)

    traj_bb = add_control_derivatives(
        base_traj,
        1;
        control_name = control_sym,
        derivative_bounds = (du_bounds,),
    )

    # Get control derivative names
    control_names =
        [name for name ∈ traj_bb.names if endswith(string(name), string(control_sym))]

    du_sym = control_names[2]
    s_du_sym = Symbol(:s_, du_sym)

    # Add slack variable
    n_drives = sys.n_drives
    s_du_data = abs.(traj_bb[du_sym])
    traj_bb = add_component(
        traj_bb,
        s_du_sym,
        s_du_data;
        type = :control,
        bounds = merge(
            traj_bb.bounds,
            NamedTuple{(s_du_sym,)}(((zeros(n_drives), fill(Inf, n_drives)),)),
        ),
    )

    # Build objective: weighted sum of infidelities for each state
    J = if free_phase && !isnothing(goals_fn)
        CoherentKetFreePhaseInfidelityObjective(goals_fn, snames, θ_names, traj_bb; Q = Q)
    else
        _ensemble_ket_objective(qtraj, traj_bb, snames, weights, goals, Q; coherent = coherent)
    end

    # L2 on control amplitude
    J += QuadraticRegularizer(control_names[1], traj_bb, R_u)

    # L1 on du via slack variables
    J += LinearRegularizer(s_du_sym, traj_bb, R_du)

    # Apply piccolo options for each state
    J += _apply_piccolo_options(qtraj, piccolo_options, constraints, traj_bb, snames)

    # Combine user-specified global_names with free_phase θ_names for integrator
    all_global_names = if !isempty(θ_names)
        gn = isnothing(global_names) ? Symbol[] : copy(global_names)
        append!(gn, θ_names)
        gn
    else
        global_names
    end

    # Build integrators
    if isnothing(integrator)
        if !isnothing(all_global_names) && !isempty(all_global_names)
            error(
                "free_phase=true or global_names requires a custom integrator that supports global variables. " *
                "Use HermitianExponentialIntegrator from Piccolissimo:\n" *
                "  using Piccolissimo\n" *
                "  integrator = HermitianExponentialIntegrator(qtraj, N; global_names=$all_global_names)\n" *
                "  qcp = BangBangPulseProblem(qtraj, N; integrator=integrator, ...)",
            )
        end
        dynamics_integrators = BilinearIntegrator(qtraj, N)
    elseif integrator isa AbstractIntegrator
        dynamics_integrators = AbstractIntegrator[integrator]
    else
        dynamics_integrators = AbstractIntegrator[integrator...]
    end

    integrators = AbstractIntegrator[dynamics_integrators...]

    # 1 derivative integrator (not 2)
    push!(integrators, DerivativeIntegrator(control_names[1], control_names[2], traj_bb))

    # Build constraints
    all_constraints = copy(constraints)
    push!(all_constraints, L1SlackConstraint(du_sym, s_du_sym, traj_bb))

    add_global_bounds_constraints!(
        all_constraints,
        global_bounds,
        traj_bb;
        verbose = _show_details(piccolo_options),
    )

    apply_calibration_targets!(
        all_constraints,
        calibration_targets,
        traj_bb;
        verbose = _show_details(piccolo_options),
    )

    prob = DirectTrajOptProblem(traj_bb, J, integrators; constraints = all_constraints)

    return _maybe_display(QuantumControlProblem(qtraj, prob), piccolo_options)
end

# ============================================================================= #
# Fallback Error Method
# ============================================================================= #

"""
    BangBangPulseProblem(qtraj::AbstractQuantumTrajectory, N::Int; kwargs...)

Fallback method that provides helpful error for non-ZeroOrderPulse types.
"""
function BangBangPulseProblem(
    qtraj::AbstractQuantumTrajectory{P},
    N::Int;
    kwargs...,
) where {P<:AbstractPulse}
    pulse_type = P
    error(
        """
  BangBangPulseProblem is only for piecewise constant pulses (ZeroOrderPulse).

  You provided a trajectory with pulse type: $(pulse_type)

  For spline-based pulses (LinearSplinePulse, CubicSplinePulse), use SplinePulseProblem instead:
      qcp = SplinePulseProblem(qtraj, N; ...)
  """,
    )
end

# ============================================================================= #
# Tests
# ============================================================================= #

@testitem "BangBangPulseProblem with UnitaryTrajectory" tags = [:experimental] begin
    using DirectTrajOpt
    using LinearAlgebra

    T = 10.0
    N = 50
    sys = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    U_goal = GATES[:H]

    pulse = ZeroOrderPulse(0.1 * randn(2, N), collect(range(0.0, T, length = N)))
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)
    qcp = BangBangPulseProblem(qtraj, N; Q = 100.0, R_du = 1e-1)

    @test qcp isa QuantumControlProblem
    @test length(qcp.prob.integrators) == 2  # dynamics + du (not ddu)
    @test haskey(qcp.prob.trajectory.components, :u)
    @test haskey(qcp.prob.trajectory.components, :du)
    @test haskey(qcp.prob.trajectory.components, :s_du)

    # Solve and verify
    solve!(qcp; max_iter = 200, print_level = 1, verbose = false)

    # Test fidelity after solve
    traj = get_trajectory(qcp)
    Ũ⃗_final = traj[end][state_name(qtraj)]
    U_final = iso_vec_to_operator(Ũ⃗_final)
    fid = unitary_fidelity(U_final, U_goal)
    @test fid > 0.9

    # Verify slack constraint: s_du >= |du|
    for t = 1:N
        du = traj[t][:du]
        s = traj[t][:s_du]
        @test all(s .>= abs.(du) .- 1e-6)
    end
end

@testitem "BangBangPulseProblem with KetTrajectory" begin
    using DirectTrajOpt
    using LinearAlgebra

    T = 10.0
    N = 50
    sys = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]

    pulse = ZeroOrderPulse(randn(2, N), collect(range(0.0, T, length = N)))
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)
    qcp = BangBangPulseProblem(qtraj, N; Q = 50.0, R_du = 1e-1)

    @test qcp isa QuantumControlProblem
    @test length(qcp.prob.integrators) == 2
    @test haskey(qcp.prob.trajectory.components, :du)
    @test haskey(qcp.prob.trajectory.components, :s_du)

    # Solve and verify
    solve!(qcp; max_iter = 200, print_level = 1, verbose = false)

    traj = get_trajectory(qcp)
    ψ̃_final = traj[end][state_name(qtraj)]
    ψ_final = iso_to_ket(ψ̃_final)
    fid = fidelity(ψ_final, ψ_goal)
    @test fid > 0.9
end

@testitem "BangBangPulseProblem rejects spline pulses" begin
    using LinearAlgebra

    T = 10.0
    N = 50
    sys = QuantumSystem(GATES[:Z], [GATES[:X]], [1.0])
    U_goal = GATES[:X]

    times = collect(range(0.0, T, length = N))

    pulse_linear = LinearSplinePulse(0.1 * randn(1, N), times)
    qtraj_linear = UnitaryTrajectory(sys, pulse_linear, U_goal)
    @test_throws ErrorException BangBangPulseProblem(qtraj_linear, N)
end

@testitem "BangBangPulseProblem with MultiKetTrajectory" tags = [:experimental] begin
    using DirectTrajOpt
    using LinearAlgebra

    T = 10.0
    N = 50
    sys = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]

    # Deterministic small smooth init. Random init was hitting two issues:
    # (1) ODE integrator querying t = t_final + ε (now structurally fixed by
    # ZeroOrderPulse's constant-extrapolation default), and (2) randn-stream
    # variability across Julia versions producing different convergence
    # outcomes — replaced with a smooth deterministic warmstart so the BangBang
    # pipeline is exercised reproducibly.
    times_arr = (0:(N-1)) ./ (N - 1)
    u_init =
        0.1 *
        vcat(reshape(cos.(2π .* times_arr), 1, N), reshape(sin.(2π .* times_arr), 1, N))
    pulse = ZeroOrderPulse(u_init, collect(range(0.0, T, length = N)))
    ensemble_qtraj = MultiKetTrajectory(sys, pulse, [ψ0, ψ1], [ψ1, ψ0])

    qcp = BangBangPulseProblem(ensemble_qtraj, N; Q = 100.0, R_du = 1e-1)

    @test qcp isa QuantumControlProblem
    @test qcp.qtraj isa MultiKetTrajectory

    @test haskey(qcp.prob.trajectory.components, :ψ̃1)
    @test haskey(qcp.prob.trajectory.components, :ψ̃2)
    @test haskey(qcp.prob.trajectory.components, :du)
    @test haskey(qcp.prob.trajectory.components, :s_du)

    # 2 dynamics + 1 derivative = 3 (not 4 like SmoothPulseProblem)
    @test length(qcp.prob.integrators) == 3

    solve!(qcp; max_iter = 200, print_level = 1, verbose = false)

    traj = get_trajectory(qcp)
    snames = state_names(ensemble_qtraj)
    goals = ensemble_qtraj.goals

    for (name, goal) in zip(snames, goals)
        ψ̃_final = traj[end][name]
        ψ_final = iso_to_ket(ψ̃_final)
        fid = fidelity(ψ_final, goal)
        @test fid > 0.9
    end
end

@testitem "BangBangPulseProblem with UnitaryTrajectory and free_phase=true" begin
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    T = 10.0
    N = 50

    # 3-level system with X gate on computational subspace
    σx = ComplexF64[0 1; 1 0]
    H_drift = Diagonal(ComplexF64[0, 0, 1.0])
    H_drive = ComplexF64[0 1 0; 1 0 0; 0 0 0]
    sys = QuantumSystem(H_drift, [H_drive], [1.0])

    subspace = get_subspace_indices([1, 2], 3)
    U_goal = EmbeddedOperator(σx, subspace, [3])

    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)
    @test_throws ErrorException BangBangPulseProblem(
        qtraj,
        N;
        Q = 100.0,
        R_du = 1e-1,
        free_phase = true,
    )
end

@testitem "BangBangPulseProblem with MultiKetTrajectory and free_phase=true" begin
    using NamedTrajectories
    using DirectTrajOpt
    using LinearAlgebra

    T = 10.0
    N = 50
    sys = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]

    pulse = ZeroOrderPulse(randn(2, N), collect(range(0.0, T, length = N)))
    qtraj = MultiKetTrajectory(sys, pulse, [ψ0, ψ1], [ψ1, ψ0])

    @test_throws ErrorException BangBangPulseProblem(
        qtraj,
        N;
        Q = 100.0,
        R_du = 1e-1,
        free_phase = true,
        subsystem_levels = [2],
    )
end

@testitem "BangBangPulseProblem free_phase requires EmbeddedOperator for unitary" begin
    using LinearAlgebra

    T = 10.0
    N = 50
    sys = QuantumSystem(GATES[:Z], [GATES[:X]], [1.0])

    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    U_goal = GATES[:X]
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)
    @test_throws AssertionError BangBangPulseProblem(qtraj, N; free_phase = true)
end

@testitem "BangBangPulseProblem MultiKet free_phase requires subsystem_levels" begin
    using LinearAlgebra

    T = 10.0
    N = 50
    sys = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    pulse = ZeroOrderPulse(randn(2, N), collect(range(0.0, T, length = N)))
    qtraj = MultiKetTrajectory(sys, pulse, [ψ0, ψ1], [ψ1, ψ0])

    @test_throws AssertionError BangBangPulseProblem(
        qtraj,
        N;
        free_phase = true,
        subsystem_levels = nothing,
    )
end
