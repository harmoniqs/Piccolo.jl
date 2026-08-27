# ============================================================================ #
# HermitianExponentialIntegrator dispatch over SamplingTrajectory (closed-system)
#
# Piccolissimo.jl#401 / parent #395. One integrator per ensemble member, built
# from the expanded sampling trajectory — mirroring Piccolo's
# `BilinearIntegrator(qtraj::SamplingTrajectory, N)` precedent. Per-member
# construction shares NO mutable scratch across members: each member integrator
# allocates its own ∂ℰs/μ∂²ℰs/exp_eigen!/DK buffer sets inside the per-base
# cores (the #354 race-fix discipline).
# ============================================================================ #

"""
    HermitianExponentialIntegrator(qtraj::SamplingTrajectory, N::Int; kwargs...)

Create a vector of HermitianExponentialIntegrators for each system in a
SamplingTrajectory — one integrator per ensemble member, built from the expanded
trajectory. Mirrors `BilinearIntegrator(qtraj::SamplingTrajectory, N)`; unlike the
bilinear path, a multi-ket member gets ONE integrator over its ket sub-states
(shared-propagator), not one per ket.

Supported closed-system bases: `UnitaryTrajectory`, `KetTrajectory`,
`MultiKetTrajectory`. Density bases route to the NonHermitianExponentialIntegrator
family (out of scope here).

# Keyword Arguments
- `global_names::Union{Vector{Symbol}, Nothing}=nothing`: global (time-invariant)
  variables in the dynamics. If `nothing`, auto-detects from the NOMINAL system's
  `global_params` (members share the global names; values are decision variables
  read from the trajectory at evaluation time).
- `gauss_newton`, `matrix_free`, `use_analytical`: as in the per-base constructors.
"""
function HermitianExponentialIntegrator(
    qtraj::SamplingTrajectory,
    N::Int;
    global_names::Union{Vector{Symbol},Nothing} = nothing,
    gauss_newton::Bool = false,
    matrix_free::Bool = false,
    use_analytical::Bool = true,
)
    base = qtraj.base_trajectory
    nominal_sys = get_system(base)

    # Resolve globals from the nominal system (members share the names), matching
    # the per-base constructors' auto-detection.
    resolved_global_names = if !isnothing(global_names)
        global_names
    elseif !isempty(nominal_sys.global_params)
        collect(keys(nominal_sys.global_params))
    else
        Symbol[]
    end

    # Expanded sampling trajectory. The sampling conversion does not attach
    # globals, so attach them here — mirroring the non-sampling conversions (which
    # auto-attach sys.global_params) and SamplingProblem's own propagation.
    traj = NamedTrajectory(qtraj, N)
    if !isempty(resolved_global_names)
        global_data = Dict{Symbol,Vector{Float64}}(
            name => [
                haskey(nominal_sys.global_params, name) ?
                nominal_sys.global_params[name] : 0.0,
            ] for name in resolved_global_names
        )
        traj = NamedTrajectory(
            traj.datavec,
            traj.components,
            traj.N;
            timestep = traj.timestep,
            controls = traj.control_names,
            bounds = traj.bounds,
            initial = traj.initial,
            final = isnothing(traj.final_) ? NamedTuple() : traj.final_,
            goal = traj.goal,
            global_data = global_data,
        )
    end

    control_sym = drive_name(qtraj)
    return map(zip(qtraj.systems, sampling_member_states(qtraj))) do (sys, states)
        _sampling_hermitian_integrator(
            base,
            sys,
            traj,
            states,
            control_sym,
            resolved_global_names;
            gauss_newton = gauss_newton,
            matrix_free = matrix_free,
            use_analytical = use_analytical,
        )
    end
end

# Per-member construction helpers — dispatch on the base trajectory type, one
# integrator per member built from the expanded trajectory (the bilinear
# `_sampling_integrator` precedent).
function _sampling_hermitian_integrator(
    base::KetTrajectory,
    sys,
    traj::NamedTrajectory,
    x::Symbol,
    u::Symbol,
    global_names::Vector{Symbol};
    kwargs...,
)
    return _hermitian_exp_ket(sys, x, u, traj, global_names; kwargs...)
end

function _sampling_hermitian_integrator(
    base::UnitaryTrajectory,
    sys,
    traj::NamedTrajectory,
    x::Symbol,
    u::Symbol,
    global_names::Vector{Symbol};
    kwargs...,
)
    return _hermitian_exp_unitary(sys, x, u, traj, global_names; kwargs...)
end

function _sampling_hermitian_integrator(
    base::MultiKetTrajectory,
    sys,
    traj::NamedTrajectory,
    xs::Vector{Symbol},
    u::Symbol,
    global_names::Vector{Symbol};
    kwargs...,
)
    return _hermitian_exp_multiket(sys, xs, u, traj, global_names; kwargs...)
end

@testitem "HermitianExponentialIntegrator dispatch on SamplingTrajectory (Ket)" begin
    using DirectTrajOpt
    using NamedTrajectories
    using Piccolo

    sys1 = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    sys2 = QuantumSystem(1.1 * GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]

    N = 11
    times = collect(range(0, 1.0, length = N))
    pulse = LinearSplinePulse(zeros(2, N), times)

    base_qtraj = KetTrajectory(sys1, pulse, ψ_init, ψ_goal)
    sampling_qtraj = SamplingTrajectory(base_qtraj, [sys1, sys2])
    expanded_traj = NamedTrajectory(sampling_qtraj, N)

    integrators = HermitianExponentialIntegrator(sampling_qtraj, N)

    # One integrator per ensemble member, each owning its member state
    @test integrators isa Vector{<:HermitianExponentialIntegrator}
    @test length(integrators) == 2
    @test [x_name(ℰ) for ℰ in integrators] == [:ψ̃1, :ψ̃2]

    # Each member integrator is consistent with the expanded trajectory:
    # test_integrator checks evaluate!/jacobian/hessian against finite
    # differences through the integrator's public interface.
    for ℰ in integrators
        test_integrator(ℰ, expanded_traj; atol = 1e-3)
    end
end

@testitem "HermitianExponentialIntegrator dispatch on SamplingTrajectory (Unitary)" begin
    using DirectTrajOpt
    using NamedTrajectories
    using Piccolo

    sys1 = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    sys2 = QuantumSystem(1.1 * GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

    N = 11
    times = collect(range(0, 1.0, length = N))
    pulse = LinearSplinePulse(zeros(2, N), times)

    base_qtraj = UnitaryTrajectory(sys1, pulse, GATES[:H])
    sampling_qtraj = SamplingTrajectory(base_qtraj, [sys1, sys2])
    expanded_traj = NamedTrajectory(sampling_qtraj, N)

    integrators = HermitianExponentialIntegrator(sampling_qtraj, N)

    @test integrators isa Vector{<:HermitianExponentialIntegrator}
    @test length(integrators) == 2
    @test [x_name(ℰ) for ℰ in integrators] == [:Ũ⃗1, :Ũ⃗2]

    for ℰ in integrators
        test_integrator(ℰ, expanded_traj; atol = 1e-3)
    end
end

@testitem "HermitianExponentialIntegrator dispatch on SamplingTrajectory (MultiKet)" begin
    using DirectTrajOpt
    using NamedTrajectories
    using Piccolo

    sys1 = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    sys2 = QuantumSystem(1.1 * GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]

    N = 11
    times = collect(range(0, 1.0, length = N))
    pulse = LinearSplinePulse(zeros(2, N), times)

    base_qtraj = MultiKetTrajectory(sys1, pulse, [ψ0, ψ1], [ψ1, ψ0])
    sampling_qtraj = SamplingTrajectory(base_qtraj, [sys1, sys2])
    expanded_traj = NamedTrajectory(sampling_qtraj, N)

    integrators = HermitianExponentialIntegrator(sampling_qtraj, N)

    # ONE integrator per member (shared-propagator multiket), not one per ket —
    # the deliberate departure from the bilinear sampling flattening.
    @test integrators isa Vector{<:HermitianExponentialIntegrator}
    @test length(integrators) == 2
    @test [ℰ.x_names for ℰ in integrators] == [[:ψ̃1, :ψ̃2], [:ψ̃3, :ψ̃4]]

    for ℰ in integrators
        test_integrator(ℰ, expanded_traj; atol = 1e-3)
    end
end

@testitem "HermitianExponentialIntegrator sampling Jacobian parity (Ket)" begin
    using DirectTrajOpt
    using NamedTrajectories
    using Piccolo

    sys = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    N = 11
    times = range(0, 1.0, length = N)
    pulse = LinearSplinePulse(zeros(2, N), times)

    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)
    direct_ℰ = HermitianExponentialIntegrator(qtraj, N)
    direct_traj = NamedTrajectory(qtraj, N)

    sampling_qtraj = SamplingTrajectory(qtraj, [sys])
    sampling_ℰ = HermitianExponentialIntegrator(sampling_qtraj, N)[1]
    expanded_traj = NamedTrajectory(sampling_qtraj, N)

    test_integrator(direct_ℰ, direct_traj; atol = 1e-3)
    test_integrator(sampling_ℰ, expanded_traj; atol = 1e-3)
end

@testitem "HermitianExponentialIntegrator sampling Jacobian parity (Unitary)" begin
    using DirectTrajOpt
    using NamedTrajectories
    using Piccolo

    sys = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    N = 11
    times = range(0, 1.0, length = N)
    pulse = LinearSplinePulse(zeros(2, N), times)

    qtraj = UnitaryTrajectory(sys, pulse, GATES[:H])
    direct_ℰ = HermitianExponentialIntegrator(qtraj, N)
    direct_traj = NamedTrajectory(qtraj, N)

    sampling_qtraj = SamplingTrajectory(qtraj, [sys])
    sampling_ℰ = HermitianExponentialIntegrator(sampling_qtraj, N)[1]
    expanded_traj = NamedTrajectory(sampling_qtraj, N)

    test_integrator(direct_ℰ, direct_traj; atol = 1e-3)
    test_integrator(sampling_ℰ, expanded_traj; atol = 1e-3)
end

@testitem "HermitianExponentialIntegrator sampling Jacobian parity (MultiKet)" begin
    using DirectTrajOpt
    using NamedTrajectories
    using Piccolo

    sys = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    N = 11
    times = range(0, 1.0, length = N)
    pulse = LinearSplinePulse(zeros(2, N), times)

    qtraj = MultiKetTrajectory(sys, pulse, [ψ0, ψ1], [ψ1, ψ0])
    direct_ℰ = HermitianExponentialIntegrator(qtraj, N)
    direct_traj = NamedTrajectory(qtraj, N)

    sampling_qtraj = SamplingTrajectory(qtraj, [sys])
    sampling_ℰ = HermitianExponentialIntegrator(sampling_qtraj, N)[1]
    expanded_traj = NamedTrajectory(sampling_qtraj, N)

    test_integrator(direct_ℰ, direct_traj; atol = 1e-3)
    test_integrator(sampling_ℰ, expanded_traj; atol = 1e-3)
end

@testitem "SamplingProblem short solve with HermitianExponentialIntegrator (CPU)" begin
    using LinearAlgebra
    using Random
    using DirectTrajOpt
    using Piccolo

    Random.seed!(42)

    sys_nom = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    sys_pert = QuantumSystem(1.05 * GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

    N = 11
    T = 1.0
    times = range(0, T, length = N)
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    pulse = ZeroOrderPulse(0.1 * randn(2, N), times)
    qtraj = KetTrajectory(sys_nom, pulse, ψ_init, ψ_goal)

    prob = SmoothPulseProblem(qtraj, N; Q = 100.0)
    sampling_bilin = SamplingProblem(prob, [sys_nom, sys_pert])
    solve!(sampling_bilin; max_iter = 50, print_level = 0)

    prob2 = SmoothPulseProblem(qtraj, N; Q = 100.0)
    sampling_exp = SamplingProblem(
        prob2,
        [sys_nom, sys_pert];
        integrator = (sq, n) -> HermitianExponentialIntegrator(sq, n),
    )
    solve!(sampling_exp; max_iter = 50, print_level = 0)
end

@testitem "HermitianExponentialIntegrator registered in specs layer (#401)" begin
    using Piccolo

    entry = Piccolo.lookup_integrator(:hermitian_exponential)
    @test entry ≢ nothing

    sys = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    N = 11
    times = range(0, 1.0, length = N)
    pulse = LinearSplinePulse(zeros(2, N), times)
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)
    sampling_qtraj = SamplingTrajectory(qtraj, [sys])

    ℰs = entry.factory(sampling_qtraj, N)
    @test ℰs isa Vector{<:HermitianExponentialIntegrator}
    @test length(ℰs) == 1
end

@testitem "HermitianExponentialIntegrator sampling per-member buffer isolation (#401)" begin
    using DirectTrajOpt
    using NamedTrajectories
    using Piccolo

    sys1 = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    sys2 = QuantumSystem(1.05 * GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

    N = 11
    times = range(0, 1.0, length = N)
    pulse = LinearSplinePulse(zeros(2, N), times)
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]

    qtraj = KetTrajectory(sys1, pulse, ψ_init, ψ_goal)
    sampling_qtraj = SamplingTrajectory(qtraj, [sys1, sys2])
    expanded_traj = NamedTrajectory(sampling_qtraj, N)

    ℰs = HermitianExponentialIntegrator(sampling_qtraj, N)
    @test length(ℰs) == 2

    F1 = zeros(ℰs[1].x_dim * (N - 1))
    F2 = zeros(ℰs[2].x_dim * (N - 1))

    evaluate!(F1, ℰs[1], expanded_traj)
    evaluate!(F2, ℰs[2], expanded_traj)
    F1_saved = copy(F1)
    F2_saved = copy(F2)

    evaluate!(F1, ℰs[1], expanded_traj)
    evaluate!(F2, ℰs[2], expanded_traj)
    @test F1 ≈ F1_saved

    evaluate!(F2, ℰs[2], expanded_traj)
    @test F1 ≈ F1_saved
end
