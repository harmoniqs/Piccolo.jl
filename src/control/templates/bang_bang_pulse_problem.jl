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
- `du_bound::Float64=Inf`: Bound on discrete first derivative
- `Δt_bounds::Union{Nothing, Tuple{Float64, Float64}}=nothing`: Timestep bounds
- `Q::Float64=100.0`: Weight on infidelity/objective
- `R::Float64=1e-2`: Default regularization weight
- `R_u::Union{Float64, Vector{Float64}}=0.0`: L2 weight on control amplitude (defaults to 0 — bang-bang needs no amplitude regularization)
- `R_du::Union{Float64, Vector{Float64}}=R`: L1 weight on `du` (applied to slacks)
- `constraints::Vector{<:AbstractConstraint}=AbstractConstraint[]`: Additional constraints
- `piccolo_options::PiccoloOptions=PiccoloOptions()`: Piccolo solver options

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
    du_bound::Float64 = Inf,
    Δt_bounds::Union{Nothing,Tuple{Float64,Float64}} = nothing,
    Q::Float64 = 100.0,
    R::Float64 = 1e-2,
    R_u::Union{Float64,Vector{Float64}} = 0.0,
    R_du::Union{Float64,Vector{Float64}} = R,
    constraints::Vector{<:AbstractConstraint} = AbstractConstraint[],
    piccolo_options::PiccoloOptions = PiccoloOptions(),
)
    if piccolo_options.verbose
        traj_type = split(string(typeof(qtraj).name.name), ".")[end]
        println("    constructing BangBangPulseProblem for $traj_type...")
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
    if isnothing(integrator)
        if !isnothing(global_names) && !isempty(global_names)
            error(
                "global_names requires a custom integrator that supports global variables. " *
                "Use HermitianExponentialIntegrator from Piccolissimo:\n" *
                "  using Piccolissimo\n" *
                "  integrator = HermitianExponentialIntegrator(qtraj, N; global_names=$global_names)\n" *
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
    J = _state_objective(qtraj, traj_bb, state_sym, Q)

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
        verbose = piccolo_options.verbose,
    )

    prob = DirectTrajOptProblem(traj_bb, J, integrators; constraints = all_constraints)

    return QuantumControlProblem(qtraj, prob)
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
    du_bound::Float64 = Inf,
    Δt_bounds::Union{Nothing,Tuple{Float64,Float64}} = nothing,
    Q::Float64 = 100.0,
    R::Float64 = 1e-2,
    R_u::Union{Float64,Vector{Float64}} = 0.0,
    R_du::Union{Float64,Vector{Float64}} = R,
    constraints::Vector{<:AbstractConstraint} = AbstractConstraint[],
    piccolo_options::PiccoloOptions = PiccoloOptions(),
)
    if piccolo_options.verbose
        println(
            "    constructing BangBangPulseProblem for MultiKetTrajectory ($(length(qtraj.initials)) states)...",
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
    J = _ensemble_ket_objective(qtraj, traj_bb, snames, weights, goals, Q)

    # L2 on control amplitude
    J += QuadraticRegularizer(control_names[1], traj_bb, R_u)

    # L1 on du via slack variables
    J += LinearRegularizer(s_du_sym, traj_bb, R_du)

    # Apply piccolo options for each state
    J += _apply_piccolo_options(qtraj, piccolo_options, constraints, traj_bb, snames)

    # Build integrators
    if isnothing(integrator)
        if !isnothing(global_names) && !isempty(global_names)
            error(
                "global_names requires a custom integrator that supports global variables. " *
                "Use HermitianExponentialIntegrator from Piccolissimo:\n" *
                "  using Piccolissimo\n" *
                "  integrator = HermitianExponentialIntegrator(qtraj, N; global_names=$global_names)\n" *
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
        verbose = piccolo_options.verbose,
    )

    prob = DirectTrajOptProblem(traj_bb, J, integrators; constraints = all_constraints)

    return QuantumControlProblem(qtraj, prob)
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

@testitem "BangBangPulseProblem with MultiKetTrajectory" tags=[:experimental] begin
    using DirectTrajOpt
    using LinearAlgebra

    T = 10.0
    N = 50
    sys = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]

    pulse = ZeroOrderPulse(randn(2, N), collect(range(0.0, T, length = N)))
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
