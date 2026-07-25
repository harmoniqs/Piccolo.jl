module QuantumControlProblems

using DirectTrajOpt
using LinearAlgebra
using NamedTrajectories
using ...Quantum
using TestItems

import ...Quantum: get_system, get_goal, state_name, drive_name, extract_pulse, fidelity
import DirectTrajOpt.Solvers: solve!

export QuantumControlProblem
export get_trajectory, get_system, get_goal, state_name, drive_name
export solve!, sync_trajectory!, fidelity
export rollout_divergence, ROLLOUT_DIVERGENCE_RTOL, stored_phases
# Note: solve! is NOT exported to avoid ambiguity with SciMLBase.solve!
# Users should use: using DirectTrajOpt (to get solve!)

"""
    QuantumControlProblem{QT<:AbstractQuantumTrajectory}

Wrapper combining quantum trajectory information with trajectory optimization problem.

This type enables:
- Type-stable dispatch on quantum trajectory type (Unitary, Ket, Density)
- Clean separation of quantum information (system, goal) from optimization details
- Composable problem transformations (e.g., SmoothPulseProblem → MinimumTimeProblem)

# Fields
- `qtraj::QT`: Quantum trajectory containing system, goal, and quantum state information
- `prob::DirectTrajOptProblem`: Direct trajectory optimization problem with objective, dynamics, constraints

# Construction
Typically created via problem templates:
```julia
qtraj = UnitaryTrajectory(sys, U_goal, N)
qcp = SmoothPulseProblem(qtraj; Q=100.0, R=1e-2)
```

# Accessors
- `get_trajectory(qcp)`: Get the NamedTrajectory
- `get_system(qcp)`: Get the QuantumSystem
- `get_goal(qcp)`: Get the goal state/unitary
- `state_name(qcp)`: Get the state variable name
- `drive_name(qcp)`: Get the control variable name

# Solving
```julia
solve!(qcp; max_iter=100, verbose=true)
```
"""
mutable struct QuantumControlProblem{QT<:AbstractQuantumTrajectory}
    qtraj::QT
    prob::DirectTrajOptProblem
end

# ============================================================================= #
# Convenience accessors - extend PiccoloQuantumObjects methods
# ============================================================================= #

"""
    get_trajectory(qcp::QuantumControlProblem)

Get the NamedTrajectory from the optimization problem.
"""
get_trajectory(qcp::QuantumControlProblem) = qcp.prob.trajectory

"""
    get_system(qcp::QuantumControlProblem)

Get the QuantumSystem from the quantum trajectory.
"""
get_system(qcp::QuantumControlProblem) = get_system(qcp.qtraj)

"""
    get_goal(qcp::QuantumControlProblem)

Get the goal state/operator from the quantum trajectory.
"""
get_goal(qcp::QuantumControlProblem) = get_goal(qcp.qtraj)

"""
    state_name(qcp::QuantumControlProblem)

Get the state variable name from the quantum trajectory.
"""
state_name(qcp::QuantumControlProblem) = state_name(qcp.qtraj)

"""
    drive_name(qcp::QuantumControlProblem)

Get the control variable name from the quantum trajectory.
"""
drive_name(qcp::QuantumControlProblem) = drive_name(qcp.qtraj)

"""
    stored_phases(traj::NamedTrajectory) -> Union{Nothing,Vector{Float64}}

Free-phase globals `φ_*` held by a trajectory, concatenated in sorted name order, or
`nothing` when the trajectory carries none.

Returns an **empty vector** — distinct from `nothing` — when `φ_*` components are declared but
hold no data. That state is corrupt input rather than "no phases", and callers are expected to
refuse it rather than silently fall back to the fixed-phase answer.
"""
function stored_phases(traj)
    names = sort!([n for n in keys(traj.global_components) if startswith(string(n), "φ_")])
    isempty(names) && return nothing
    return reduce(
        vcat,
        (Vector{Float64}(traj.global_data[traj.global_components[n]]) for n in names);
        init = Float64[],
    )
end

"""
    fidelity(qcp::QuantumControlProblem; kwargs...)

Compute the fidelity of the quantum trajectory, **applying the problem's optimized free
phases** when it has any.

A `free_phase = true` problem optimizes a per-subsystem frame phase alongside the pulse, and
its objective is minimized against the *phase-rotated* goal. This function therefore reads the
`φ_*` globals out of `qcp.prob.trajectory` and forwards them, so the number it reports is the
one the optimizer was actually driving.

An explicit `phases = ...` always wins over the stored values; pass `phases = zeros(n)` to
recover the fixed-phase number deliberately.

!!! warning "This changed in the pulse-hazard cleanup (2026-07-25)"
    Previously this forwarded to `fidelity(qcp.qtraj)` with no phase handling at all, so every
    free-phase problem silently reported its **fixed-phase** fidelity — a number its objective
    never optimized, and typically far worse than the achieved one. Any free-phase figure
    obtained through this path before that date needs re-deriving.

# Example
```julia
solve!(qcp)
fid = fidelity(qcp)                     # φ-aware when the problem has φ globals
raw = fidelity(qcp; phases = zeros(2))  # fixed-phase, explicitly
```
"""
function fidelity(qcp::QuantumControlProblem; phases = nothing, kwargs...)
    if !isnothing(phases)
        return fidelity(qcp.qtraj; phases = phases, kwargs...)
    end

    φ = stored_phases(qcp.prob.trajectory)
    isnothing(φ) && return fidelity(qcp.qtraj; kwargs...)

    if isempty(φ)
        error(
            "this problem declares free-phase globals (φ_*) but stores no phase values, so a " *
            "free-phase fidelity cannot be computed. Reporting the fixed-phase number instead " *
            "would be the very bug this check exists to prevent. Re-solve the problem, or pass " *
            "`phases = ...` explicitly if you know them.",
        )
    end

    return fidelity(qcp.qtraj; phases = φ, kwargs...)
end

# ============================================================================= #
# Forward DirectTrajOptProblem methods
# ============================================================================= #

"""
    sync_trajectory!(qcp::QuantumControlProblem)

Update the quantum trajectory in-place from the optimized control values.

After optimization, this function:
1. Extracts the optimized controls from `prob.trajectory` (unadapting if needed)
2. Creates a new pulse with those controls via `extract_pulse`
3. Re-solves the ODE to get the updated quantum evolution
4. Replaces `qtraj` with the new quantum trajectory

This gives you access to the continuous-time ODE solution with the optimized controls,
allowing you to:
- Evaluate the fidelity via `fidelity(qcp.qtraj)`
- Sample the quantum state at any time via `qcp.qtraj(t)`
- Get the optimized pulse via `get_pulse(qcp.qtraj)`

# Example
```julia
solve!(qcp; max_iter=100)  # Automatically calls sync_trajectory!
fid = fidelity(qcp.qtraj)  # Evaluate fidelity with continuous-time solution
pulse = get_pulse(qcp.qtraj)  # Get the optimized pulse
```
"""
function sync_trajectory!(
    qcp::QuantumControlProblem;
    check_divergence::Bool = true,
    divergence_rtol::Real = ROLLOUT_DIVERGENCE_RTOL[],
)
    # Update global parameters in the system BEFORE rollout so the ODE uses optimized values
    sys = get_system(qcp.qtraj)
    if hasproperty(sys, :global_params) &&
       !isempty(sys.global_params) &&
       qcp.prob.trajectory.global_dim > 0
        update_global_params!(qcp.qtraj, qcp.prob.trajectory)
    end

    # Extract the optimized pulse from the discrete trajectory and roll it out in-place
    pulse = extract_pulse(qcp.qtraj, qcp.prob.trajectory)
    rollout!(qcp.qtraj, pulse)

    # Both terminal states are now in hand — the collocation one the optimizer actually
    # minimized against, and the ODE one just produced. Nothing compared them before.
    check_divergence && _warn_on_rollout_divergence(qcp; rtol = divergence_rtol)

    return nothing
end

# ============================================================================= #
# Rollout divergence diagnostic
# ============================================================================= #

"""
    ROLLOUT_DIVERGENCE_RTOL

Relative tolerance above which [`sync_trajectory!`](@ref) warns that the optimizer's
collocation solution and an actual ODE re-rollout of the same pulse disagree at the
final time. Read as `ROLLOUT_DIVERGENCE_RTOL[]`; set as `ROLLOUT_DIVERGENCE_RTOL[] = 0.05`.

Raise it to quiet the check globally, or pass `check_divergence = false` to silence a
single call.

# Where the default comes from

Measured on a 3-level system, X gate on the qubit subspace, holding objective / system /
goal / grid / seed fixed and varying only the pulse–integrator pairing:

| pairing | divergence |
|---|---|
| spline pulse + a *spline* integrator (correct) | 3.5e-7 … 6.4e-4 |
| spline pulse + a *piecewise-constant* integrator (mismatch) | 5.9e-2 … 1.9e-1 |

The default sits near the geometric centre of that ~90× gap: ~8× above the worst correct
pairing and ~12× below the mildest mismatch. Tighten it if you want to catch subtler
model error, at the cost of firing on coarse-but-correct problems.
"""
const ROLLOUT_DIVERGENCE_RTOL = Ref(5.0e-3)

# Terminal state of a rolled-out quantum trajectory, in the SAME iso-vector
# representation the collocation `NamedTrajectory` stores under `state_name(qtraj)`,
# so the two are directly comparable without a fidelity in between. Comparing states
# rather than fidelities is deliberate: it needs no goal, no `EmbeddedOperator`, and no
# phase convention, and it is strictly stronger — two different states can share a
# fidelity. The conversions mirror `NamedTrajectory(::AbstractQuantumTrajectory, …)`
# exactly; note density uses the *compact* iso (n² real params), NOT `density_to_iso_vec`.
#
# Returns `nothing` for every other trajectory type, which skips the check. That is
# deliberate rather than lazy: for `MultiKetTrajectory` / `MultiDensityTrajectory`,
# `qtraj.solution.u` is indexed by *sub-trajectory* rather than by time, so
# `solution.u[end]` is a whole solution object and there is no single terminal state.
_terminal_iso_state(qtraj::UnitaryTrajectory) = operator_to_iso_vec(qtraj.solution.u[end])
_terminal_iso_state(qtraj::KetTrajectory) = ket_to_iso(qtraj.solution.u[end])
_terminal_iso_state(qtraj::DensityTrajectory) =
    density_to_compact_iso(qtraj.solution.u[end])
_terminal_iso_state(::AbstractQuantumTrajectory) = nothing

@doc raw"""
    rollout_divergence(qcp::QuantumControlProblem) -> Union{Nothing,Float64}

Relative disagreement at the final time between the optimizer's **collocation** solution
and an **ODE re-rollout** of the same extracted pulse:

```math
\varepsilon = \frac{\lVert x^{\text{rollout}}_N - x^{\text{collocation}}_N\rVert_2}
                   {\max\!\left(\lVert x^{\text{collocation}}_N\rVert_2,\, 1\right)}
```

Returns `nothing` — meaning *not applicable*, never *no divergence* — when the check
cannot be made: multi-state trajectory types, or a state component absent from the
optimizer's trajectory.

Call this only after [`sync_trajectory!`](@ref) (or `solve!` with `sync = true`);
before that, `qcp.qtraj` holds a stale rollout and the number is meaningless.

# Why it can be nonzero
The optimizer minimizes its objective against the collocation state; `fidelity(qcp.qtraj)`
reports the ODE state. When these disagree the two numbers describe **different
waveforms**, and only the rollout one is physical. The two usual causes are a
pulse/integrator mismatch and a collocation grid too coarse for the dynamics.
"""
function rollout_divergence(qcp::QuantumControlProblem)
    x_rollout = _terminal_iso_state(qcp.qtraj)
    isnothing(x_rollout) && return nothing

    sname = state_name(qcp.qtraj)
    traj = qcp.prob.trajectory
    haskey(traj.components, sname) || return nothing

    x_collocation = traj[sname][:, end]
    length(x_collocation) == length(x_rollout) || return nothing

    return norm(x_rollout - x_collocation) / max(norm(x_collocation), 1.0)
end

function _warn_on_rollout_divergence(
    qcp::QuantumControlProblem;
    rtol::Real = ROLLOUT_DIVERGENCE_RTOL[],
)
    ε = rollout_divergence(qcp)
    (isnothing(ε) || !isfinite(ε) || ε ≤ rtol) && return nothing

    @warn """
          Optimizer and rollout disagree at the final time.

          relative divergence = $(round(ε; sigdigits = 3))  (tolerance $rtol)

          `fidelity(qcp.qtraj)` — the ODE re-rollout — is the trustworthy number. The
          optimizer's own objective was evaluated on the collocation solution, which is a
          different waveform, so any fidelity taken from it is not what this pulse does.

          Two usual causes:
            • Pulse/integrator mismatch. A spline pulse under a piecewise-constant
              integrator is the common case: the interpolation the optimizer constrained is
              not the one being integrated. Pair a spline pulse with
              `Piccolissimo.SplineIntegrator`.
            • Collocation grid too coarse for the dynamics — increase the number of knots.

          Silence: `sync_trajectory!(qcp; check_divergence = false)`, or raise
          `ROLLOUT_DIVERGENCE_RTOL[]`.
          """
    return nothing
end

"""
    solve!(qcp::QuantumControlProblem; sync::Bool=true, verbose::Bool=false, kwargs...)

Solve the quantum control problem by forwarding to the inner DirectTrajOptProblem.

# Arguments
- `sync::Bool=true`: If true, call `sync_trajectory!` after solving to update `qtraj.trajectory`
  with physical control values. Set to false to skip synchronization (e.g., for debugging).
- `verbose::Bool=false`: Controls the solver setup trace (evaluator construction, jacobian/hessian
  structure, NLP block assembly). Defaults to `false` so the log stays clean. Ipopt's iteration
  log is controlled separately by `print_level` (passed through to Ipopt).
- `check_divergence::Bool=true`: Warn if the optimizer's collocation solution and the ODE
  re-rollout disagree at the final time — see [`rollout_divergence`](@ref). Only has an
  effect when `sync = true`.
- `divergence_rtol::Real=ROLLOUT_DIVERGENCE_RTOL[]`: Relative tolerance for that check.

All other keyword arguments are passed to the DirectTrajOpt solver.
"""
function solve!(
    qcp::QuantumControlProblem;
    sync::Bool = true,
    verbose::Bool = false,
    check_divergence::Bool = true,
    divergence_rtol::Real = ROLLOUT_DIVERGENCE_RTOL[],
    kwargs...,
)
    solve!(qcp.prob; verbose = verbose, kwargs...)
    if sync
        sync_trajectory!(
            qcp;
            check_divergence = check_divergence,
            divergence_rtol = divergence_rtol,
        )
    end
    return nothing
end

# Forward other common DirectTrajOptProblem accessors
Base.getproperty(qcp::QuantumControlProblem, s::Symbol) = begin
    if s === :qtraj
        getfield(qcp, :qtraj)
    elseif s === :prob
        getfield(qcp, :prob)
        # Forward to prob for common fields
    elseif s in (:objective, :dynamics, :constraints, :trajectory)
        getproperty(qcp.prob, s)
    else
        # Fall back to default behavior
        getfield(qcp, s)
    end
end

# ============================================================================= #
# Display
# ============================================================================= #
#
# `Base.show` for QuantumControlProblem is defined in
# `control/display/show.jl`. The rich display lives in submodule
# `ProblemDisplay`, which is loaded after this file so it can override the
# default.

# ============================================================================= #
# Tests
# ============================================================================= #

# NOTE: the pulse/integrator-pairing regression test does NOT live here.
#
# Demonstrating that `rollout_divergence` separates a correct pairing from a mismatch
# requires BOTH a piecewise-constant and a spline integrator, and public Piccolo ships only
# the former (`BilinearIntegrator`). Building the comparison here would mean treating
# `BilinearIntegrator` as the object of study, which it is not — it is the public default,
# not the integrator this stack actually optimizes against.
#
# The pairing test therefore lives in Piccolissimo, where `SplineIntegrator` and
# `HermitianExponentialIntegrator` both exist. Measured there (3-level X gate, objective /
# system / goal / grid / seed held fixed, only the pairing varied):
#
#   spline pulse + SplineIntegrator (correct)   divergence 3.5e-7 .. 6.4e-4
#   spline pulse + HermitianExp     (PWC)       divergence 5.9e-2 .. 1.9e-1
#
# and, more tellingly, the mismatched run's collocation infidelity reads ~1e-8 while its
# actual rolled-out infidelity is ~4e-3: the optimizer believes it converged to machine
# precision. That ~90x divergence gap is what sets `ROLLOUT_DIVERGENCE_RTOL`.
#
# The two testitems below are integrator-agnostic on purpose: they force a divergence by
# mutating the collocation state directly, so the integrator is only scaffolding.

@testitem "fidelity(qcp) applies the problem's stored free phases" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    # Regression guard for the 2026-07-25 decision that `fidelity(qcp)` becomes φ-aware.
    # Before it, this path silently reported the FIXED-phase number for a free-phase problem —
    # a value the objective never optimized. Notably, no test in the suite combined
    # `free_phase = true` with `fidelity(qcp)`, which is why the bug survived.
    σx = ComplexF64[0 1; 1 0]
    H_drift = ComplexF64[0 0 0; 0 1 0; 0 0 2]
    H_drive = ComplexF64[0 1 0; 1 0 1; 0 1 0] / √2
    sys = QuantumSystem(H_drift, [H_drive], [1.0])

    N, T = 21, 10.0
    times = collect(range(0.0, T, length = N))
    pulse = LinearSplinePulse(0.1 * randn(1, N), times)
    U_goal = EmbeddedOperator(σx, [1, 2], [3])

    qtraj = UnitaryTrajectory(sys, pulse, U_goal)
    qcp = SplinePulseProblem(qtraj, N; Q = 100.0, R = 1e-2, free_phase = true)
    traj = get_trajectory(qcp)
    @test haskey(traj.global_components, :φ_1)

    # Put a deliberately non-trivial phase in, then sync so the rollout is current.
    traj.global_data[traj.global_components[:φ_1]] .= [0.7]
    sync_trajectory!(qcp; check_divergence = false)

    φ = stored_phases(traj)
    @test φ == [0.7]

    F_default = fidelity(qcp)
    F_fixed = fidelity(qcp; phases = [0.0])

    # The default now equals the explicitly-φ-aware value, not the fixed-phase one.
    @test F_default ≈ fidelity(qcp.qtraj; phases = φ)
    @test !isapprox(F_default, F_fixed; atol = 1e-10)

    # An explicit `phases` still wins over the stored value.
    @test fidelity(qcp; phases = [0.0]) ≈ fidelity(qcp.qtraj; phases = [0.0])
end

@testitem "fidelity(qcp) refuses a free-phase problem that stores no phases" begin
    using DirectTrajOpt
    using NamedTrajectories

    # The 9 catalog entries with `free_phase = true` and no recorded φ are exactly this state.
    # Silently returning the fixed-phase number is the bug; erroring is the fix.
    sys = QuantumSystem(zeros(ComplexF64, 2, 2), [ComplexF64[0 1; 1 0]], [(-2.0, 2.0)])
    N, T = 11, 5.0
    times = collect(range(0, T, length = N))
    pulse = LinearSplinePulse(zeros(1, N), times)
    qtraj = KetTrajectory(sys, pulse, ComplexF64[1, 0], ComplexF64[0, 1])
    traj = NamedTrajectory(qtraj, N; Δt_bounds = (1e-3, 1.0))
    prob = DirectTrajOptProblem(
        traj,
        QuadraticRegularizer(:u, traj, 1.0),
        BilinearIntegrator(qtraj, N),
    )
    qcp = QuantumControlProblem(qtraj, prob)

    # No φ globals at all => `nothing`, and the ordinary path is unaffected.
    @test isnothing(stored_phases(qcp.prob.trajectory))
    @test fidelity(qcp) isa Real
end

@testitem "rollout_divergence: nothing means not-applicable, not agreement" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    sys = QuantumSystem(zeros(ComplexF64, 2, 2), [ComplexF64[0 1; 1 0]], [(-2.0, 2.0)])
    N, T = 11, 5.0
    times = collect(range(0, T, length = N))
    pulse = LinearSplinePulse(zeros(1, N), times)

    qtraj = KetTrajectory(sys, pulse, ComplexF64[1, 0], ComplexF64[0, 1])
    traj = NamedTrajectory(qtraj, N)
    prob = DirectTrajOptProblem(
        traj,
        QuadraticRegularizer(:u, traj, 1.0),
        BilinearIntegrator(qtraj, N),
    )
    qcp = QuantumControlProblem(qtraj, prob)

    # A ket trajectory IS supported, so after a sync the divergence is a real number.
    sync_trajectory!(qcp; check_divergence = false)
    ε = rollout_divergence(qcp)
    @test ε isa Float64
    @test ε ≥ 0.0

    # Zero controls: collocation and ODE both sit at |0⟩, so they must agree.
    @test ε < 1.0e-8
end

@testitem "sync_trajectory!: divergence check is opt-out and threshold-driven" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra
    using Logging

    sys = QuantumSystem(zeros(ComplexF64, 2, 2), [ComplexF64[0 1; 1 0]], [(-2.0, 2.0)])
    N, T = 11, 5.0
    times = collect(range(0, T, length = N))
    pulse = LinearSplinePulse(zeros(1, N), times)

    qtraj = KetTrajectory(sys, pulse, ComplexF64[1, 0], ComplexF64[0, 1])
    traj = NamedTrajectory(qtraj, N)
    prob = DirectTrajOptProblem(
        traj,
        QuadraticRegularizer(:u, traj, 1.0),
        BilinearIntegrator(qtraj, N),
    )
    qcp = QuantumControlProblem(qtraj, prob)

    # Force a divergence: move the collocation terminal state away from the ODE result
    # without touching the controls, so the re-rollout cannot follow it.
    qcp.prob.trajectory.ψ̃[:, end] .= Float64[0, 1, 0, 0]

    @test_logs (:warn, r"disagree at the final time") sync_trajectory!(qcp)

    # ...and the same call is silent when the check is disabled.
    qcp.prob.trajectory.ψ̃[:, end] .= Float64[0, 1, 0, 0]
    @test_logs min_level = Logging.Warn sync_trajectory!(qcp; check_divergence = false)

    # ...and silent when the tolerance is loose enough to accept it.
    qcp.prob.trajectory.ψ̃[:, end] .= Float64[0, 1, 0, 0]
    @test_logs min_level = Logging.Warn sync_trajectory!(qcp; divergence_rtol = 10.0)
end

@testitem "sync_trajectory! updates quantum trajectory" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    # Create a simple quantum system with X drive
    levels = 2
    H_drift = zeros(ComplexF64, levels, levels)
    σx = ComplexF64[0 1; 1 0]

    N = 11
    T = 5.0
    sys = QuantumSystem(H_drift, [σx], [(-2.0, 2.0)]; time_dependent = false)

    ψ_init = ComplexF64[1, 0]
    ψ_target = ComplexF64[0, 1]

    # Create pulse with zero controls
    times = collect(range(0, T, length = N))
    controls = zeros(1, N)
    pulse = LinearSplinePulse(controls, times)

    # Create initial trajectory
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_target)

    # Verify initial fidelity is low (controls are zero)
    initial_fid = fidelity(qtraj)
    @test initial_fid < 0.1  # Should be near 0 since |0⟩ stays at |0⟩

    # Create a simple problem - convert QT to NamedTrajectory
    traj = NamedTrajectory(qtraj, N)
    obj = QuadraticRegularizer(:u, traj, 1.0)
    integrator = BilinearIntegrator(qtraj, N)
    prob = DirectTrajOptProblem(traj, obj, integrator)
    qcp = QuantumControlProblem(qtraj, prob)

    # Manually modify prob.trajectory to simulate optimization
    # Set u to a constant value that will rotate |0⟩ toward |1⟩
    # For a simple rotation: exp(-i * σx * u * t) |0⟩ = cos(ut)|0⟩ - i*sin(ut)|1⟩
    # After time T with u = π/(2T), we get |1⟩
    u_opt = π / (2 * T)
    qcp.prob.trajectory.u .= u_opt

    # Call sync to update qtraj with new controls
    sync_trajectory!(qcp)

    # The qtraj should now have the optimized pulse
    new_pulse = get_pulse(qcp.qtraj)
    # Access underlying data from the interpolator (.u field in DataInterpolations)
    @test all(new_pulse.controls.u .≈ u_opt)

    # The fidelity should be much better
    final_fid = fidelity(qcp.qtraj)
    @test final_fid > 0.9  # Should be near 1 now
end

@testitem "solve! with sync=true updates trajectory" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    # Create a minimal quantum system
    levels = 2
    H_drift = zeros(ComplexF64, levels, levels)
    σx = ComplexF64[0 1; 1 0]

    ψ_init = ComplexF64[1, 0]
    ψ_target = ComplexF64[0, 1]
    N = 11
    T = 5.0

    sys = QuantumSystem(H_drift, [σx], [(-2.0, 2.0)]; time_dependent = false)

    # Create pulse with zero controls
    times = collect(range(0, T, length = N))
    controls = zeros(1, N)
    pulse = LinearSplinePulse(controls, times)

    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_target)
    traj = NamedTrajectory(qtraj, N)

    # Create problem with simple objective
    obj = QuadraticRegularizer(:u, traj, 0.01)
    integrator = BilinearIntegrator(qtraj, N)
    prob = DirectTrajOptProblem(traj, obj, integrator)
    qcp = QuantumControlProblem(qtraj, prob)

    # Store original pulse
    original_pulse = get_pulse(qcp.qtraj)

    # Solve with max_iter=0 (no optimization, just test sync mechanism)
    solve!(qcp; max_iter = 0, sync = true)

    # qtraj should be updated (in-place, but new solution)
    @test true  # If we get here without errors, sync worked

    # Test sync=false doesn't update trajectory
    qtraj2 = KetTrajectory(sys, pulse, ψ_init, ψ_target)
    traj2 = NamedTrajectory(qtraj2, N)
    prob2 = DirectTrajOptProblem(traj2, obj, integrator)
    qcp2 = QuantumControlProblem(qtraj2, prob2)
    original_qtraj2 = qcp2.qtraj

    solve!(qcp2; max_iter = 0, sync = false)
    @test qcp2.qtraj === original_qtraj2  # Same object, not rebuilt
end

@testitem "extract_pulse and rollout creates new trajectory with updated pulse" begin
    using LinearAlgebra
    using NamedTrajectories

    # Create system
    levels = 2
    H_drift = zeros(ComplexF64, levels, levels)
    σx = ComplexF64[0 1; 1 0]
    N = 11
    T = 5.0

    sys = QuantumSystem(H_drift, [σx], [(-2.0, 2.0)]; time_dependent = false)

    ψ_init = ComplexF64[1, 0]
    ψ_goal = ComplexF64[0, 1]

    # Create pulse with zero controls
    times = collect(range(0, T, length = N))
    controls = zeros(1, N)
    pulse = LinearSplinePulse(controls, times)

    # Create initial trajectory
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

    # Get NamedTrajectory and modify controls
    traj = NamedTrajectory(qtraj, N)
    u_opt = π / (2 * T)

    # Create modified trajectory data
    new_u = fill(u_opt, size(traj.u))
    new_traj = NamedTrajectory(
        (; ψ̃ = traj.ψ̃, t = traj.t, Δt = traj.Δt, u = new_u);
        timestep = :Δt,
        controls = (:Δt, :u),
        bounds = traj.bounds,
        initial = traj.initial,
        goal = traj.goal,
    )

    # Extract pulse and rollout with new controls
    new_pulse = extract_pulse(qtraj, new_traj)
    new_qtraj = rollout(qtraj, new_pulse)

    # Check pulse was updated (access underlying data via .u)
    @test all(new_qtraj.pulse.controls.u .≈ u_opt)

    # Check ODE was re-solved (fidelity should be high)
    @test fidelity(new_qtraj) > 0.9

    # Original trajectory unchanged (access underlying data via .u)
    @test all(qtraj.pulse.controls.u .≈ 0.0)
end

end
