# ============================================================================ #
# Solution Validation
#
# `validate_solution` self-certifies a solved `QuantumControlProblem` by
# comparing two *independently computed* fidelities:
#
#   1. `F_reported`  — the fidelity the optimizer's objective drove to, read
#                      from the final state stored in the discrete optimized
#                      `NamedTrajectory`. This is the number a user sees on the
#                      solve log / via `fidelity(qcp)`.
#   2. `F_rerolled`  — the fidelity of an *independent* re-rollout: the optimized
#                      pulse is extracted and propagated fresh through the
#                      standard exponential/Magnus (bare) propagator, then scored
#                      with the same fidelity convention the objective used.
#
# When the optimizer's own integrator disagrees with an honest re-integration of
# the same pulse, these two numbers diverge — sometimes catastrophically (an
# observed `F_reported = 0.999965` vs `F_rerolled = 0.000364` on a Stanford
# Fock-|1⟩ state-transfer, from a SplineIntegrator ↔ bare-propagator mismatch).
# `validate_solution` turns that silent wrongness into a loud, machine-checkable
# `pass == false` so unattended multi-run scheduling can self-certify every solve.
# ============================================================================ #

module Validation

export ValidationResult
export validate_solution

using LinearAlgebra
using NamedTrajectories

using ...Quantum
using ...Quantum: get_goal, get_system, state_name, state_names, extract_pulse
using ...Quantum.Rollouts: rollout, fidelity
using ...Quantum.QuantumTrajectories:
    AbstractQuantumTrajectory,
    UnitaryTrajectory,
    KetTrajectory,
    MultiKetTrajectory
using ...Quantum.EmbeddedOperators: EmbeddedOperator, unembed
using ...Quantum.Isomorphisms: iso_vec_to_operator, iso_to_ket

using ..QuantumControlProblems: QuantumControlProblem, get_trajectory
using ..QuantumObjectives: unitary_fidelity_loss, ket_fidelity_loss

using TestItems

# ---------------------------------------------------------------------------- #
# Result type
# ---------------------------------------------------------------------------- #

"""
    ValidationResult

Outcome of [`validate_solution`](@ref): a comparison between the optimizer's
reported fidelity and an independent re-rollout of the same pulse.

# Fields
- `F_reported::Float64`: Fidelity from the optimizer's discrete solution
  (the value the objective minimized against — what `fidelity(qcp)` reports).
- `F_rerolled::Float64`: Fidelity of an independent re-rollout of the extracted
  pulse through the standard exponential/Magnus (bare) propagator, scored with
  the same fidelity convention.
- `ΔF::Float64`: `abs(F_reported - F_rerolled)` — the disagreement.
- `pass::Bool`: `true` iff `ΔF ≤ atol`.
- `atol::Float64`: The tolerance used for the pass/fail decision.
"""
struct ValidationResult
    F_reported::Float64
    F_rerolled::Float64
    ΔF::Float64
    pass::Bool
    atol::Float64
end

function Base.show(io::IO, ::MIME"text/plain", r::ValidationResult)
    status = r.pass ? "PASS" : "FAIL"
    println(io, "ValidationResult ($status)")
    println(io, "  F_reported: ", r.F_reported)
    println(io, "  F_rerolled: ", r.F_rerolled)
    println(io, "  ΔF:         ", r.ΔF)
    print(io, "  atol:       ", r.atol)
end

Base.show(io::IO, r::ValidationResult) = print(
    io,
    "ValidationResult(",
    r.pass ? "pass" : "FAIL",
    ", F_reported=",
    r.F_reported,
    ", F_rerolled=",
    r.F_rerolled,
    ", ΔF=",
    r.ΔF,
    ")",
)

# ---------------------------------------------------------------------------- #
# Reported fidelity — read directly from the discrete optimized trajectory
#
# These recompute exactly the terminal loss the objective used, from the *final*
# discrete state column of the optimized `NamedTrajectory`. They intentionally
# do NOT re-solve any ODE — that is what makes them "what the optimizer saw".
# ---------------------------------------------------------------------------- #

function _reported_fidelity(
    qtraj::UnitaryTrajectory,
    traj::NamedTrajectory,
    phases::Union{Nothing,AbstractVector{<:Real}},
)
    s_name = state_name(qtraj)
    Ũ⃗_final = traj[s_name][:, end]
    goal = get_goal(qtraj)
    if isnothing(phases)
        return unitary_fidelity_loss(Ũ⃗_final, goal)
    else
        # Free-phase goal: apply per-qubit Z-phases to the computational subspace,
        # matching `_make_free_phase_goal` / `fidelity(::UnitaryTrajectory; phases)`.
        goal isa EmbeddedOperator ||
            error("`phases` requires an EmbeddedOperator goal for UnitaryTrajectory")
        U_base = unembed(goal)
        n_sub = size(U_base, 1)
        n_qubits = length(phases)
        phase_diag = _phase_diagonal(phases, n_sub, n_qubits)
        U_goal_sub = Diagonal(phase_diag) * U_base
        U_sub = iso_vec_to_operator(Ũ⃗_final)[goal.subspace, goal.subspace]
        n = length(goal.subspace)
        M = U_goal_sub' * U_sub
        return 1 / (n * (n + 1)) * (abs(tr(M' * M)) + abs2(tr(M)))
    end
end

function _reported_fidelity(
    qtraj::KetTrajectory,
    traj::NamedTrajectory,
    phases::Union{Nothing,AbstractVector{<:Real}},
)
    isnothing(phases) ||
        error("`phases` is not supported for KetTrajectory validation")
    s_name = state_name(qtraj)
    ψ̃_final = traj[s_name][:, end]
    return ket_fidelity_loss(ψ̃_final, ComplexF64.(get_goal(qtraj)))
end

function _reported_fidelity(
    qtraj::MultiKetTrajectory,
    traj::NamedTrajectory,
    phases::Union{Nothing,AbstractVector{<:Real}},
)
    isnothing(phases) ||
        error("`phases` is not supported for MultiKetTrajectory validation")
    goals = get_goal(qtraj)
    n = length(goals)
    snames = state_names(qtraj)
    # Coherent fidelity — same convention as `fidelity(::MultiKetTrajectory)` and
    # `CoherentKetInfidelityObjective`: F = |1/n ∑ᵢ ⟨ψᵢ_goal|ψᵢ⟩|²
    overlap_sum = sum(
        ComplexF64.(goals[i])' * iso_to_ket(traj[snames[i]][:, end]) for i = 1:n
    )
    return abs2(overlap_sum / n)
end

# per-qubit Z-phase diagonal over a computational subspace (binary decomposition)
function _phase_diagonal(phases, n_sub::Int, n_qubits::Int)
    return map(1:n_sub) do i
        bits = i - 1
        phase = sum(
            phases[j] for j = 1:n_qubits if (bits >> (n_qubits - j)) & 1 == 1;
            init = 0.0,
        )
        return exp(im * phase)
    end
end

# ---------------------------------------------------------------------------- #
# validate_solution
# ---------------------------------------------------------------------------- #

"""
    validate_solution(qcp::QuantumControlProblem; atol=1e-4, algorithm=nothing,
                      n_save=101, phases=nothing, verbose=true) -> ValidationResult

Self-certify a solved `QuantumControlProblem` by comparing the optimizer's
reported fidelity against an **independent re-rollout** of the same pulse.

Two fidelities are computed with the *same* fidelity convention the objective
used (standard `|tr(U'U_goal)|²/N²` for matrix goals, the Pedersen average-gate
formula for `EmbeddedOperator` goals, coherent ket overlap for state transfer):

- `F_reported`: read from the **final state** of the discrete optimized
  trajectory (`get_trajectory(qcp)`) — the number the optimizer converged to and
  the one `fidelity(qcp)` reports. No ODE is re-solved for this value.
- `F_rerolled`: the optimized pulse is extracted (`extract_pulse`) and propagated
  **fresh** through the standard exponential/Magnus (bare) propagator
  (`rollout`), then scored. This is the honest, integrator-independent check.

If `abs(F_reported - F_rerolled) > atol` the result is marked `pass = false` and
(when `verbose`) a `@warn` is emitted reporting **both** numbers. A large `ΔF`
signals that the optimizer's own integrator disagrees with an honest
re-integration of the pulse — e.g. a `SplineIntegrator` ↔ bare-propagator
mismatch — which would otherwise pass silently.

# Arguments
- `qcp::QuantumControlProblem`: A solved problem (call `solve!` first). Supported
  trajectory types: `UnitaryTrajectory`, `KetTrajectory`, `MultiKetTrajectory`.

# Keyword Arguments
- `atol::Real=1e-4`: Pass threshold on `ΔF = |F_reported - F_rerolled|`.
- `algorithm=nothing`: ODE algorithm for the re-rollout. Defaults to the system's
  `default_algorithm` (bare exponential/Magnus).
- `n_save::Int=101`: Save points for the re-rollout.
- `phases=nothing`: Optional per-qubit Z-phases (for free-phase gates with an
  `EmbeddedOperator` goal); applied consistently to both fidelities.
- `verbose::Bool=true`: Emit a `@warn` with both fidelities on failure.

# Returns
A [`ValidationResult`](@ref) with `F_reported`, `F_rerolled`, `ΔF`, `pass`, `atol`.

# Example
```julia
solve!(qcp; max_iter=300)
result = validate_solution(qcp; atol=1e-4)
result.pass || @warn "Solve did not self-certify" result
```
"""
function validate_solution(
    qcp::QuantumControlProblem;
    atol::Real = 1e-4,
    algorithm = nothing,
    n_save::Int = 101,
    phases::Union{Nothing,AbstractVector{<:Real}} = nothing,
    verbose::Bool = true,
)
    qtraj = qcp.qtraj
    traj = get_trajectory(qcp)

    # 1. Reported fidelity — from the discrete optimized trajectory (no re-solve).
    F_reported = _reported_fidelity(qtraj, traj, phases)

    # 2. Re-rolled fidelity — extract the pulse and propagate it fresh through the
    #    bare exponential/Magnus propagator, on the same system + goal.
    pulse = extract_pulse(qtraj, traj)
    qtraj_rerolled = if isnothing(algorithm)
        rollout(qtraj, pulse; n_save = n_save)
    else
        rollout(qtraj, pulse; algorithm = algorithm, n_save = n_save)
    end
    F_rerolled = _rerolled_fidelity(qtraj_rerolled, phases)

    ΔF = abs(F_reported - F_rerolled)
    pass = ΔF ≤ atol

    if !pass && verbose
        @warn(
            "validate_solution: reported and re-rolled fidelities disagree — " *
            "the optimizer's integrator may not match an honest re-integration " *
            "of the pulse (e.g. SplineIntegrator ↔ bare-propagator mismatch).",
            F_reported,
            F_rerolled,
            ΔF,
            atol,
        )
    end

    return ValidationResult(F_reported, F_rerolled, ΔF, pass, Float64(atol))
end

# Re-rolled fidelity dispatches through the public `fidelity` API so it uses the
# exact same convention as `F_reported` (Pedersen for EmbeddedOperator, etc.).
_rerolled_fidelity(qtraj::UnitaryTrajectory, phases) = fidelity(qtraj; phases = phases)
_rerolled_fidelity(qtraj::KetTrajectory, ::Nothing) = fidelity(qtraj)
function _rerolled_fidelity(qtraj::MultiKetTrajectory, phases)
    return isnothing(phases) ? fidelity(qtraj) :
           error("`phases` re-rollout for MultiKetTrajectory requires subsystem_levels")
end

# ---------------------------------------------------------------------------- #
# Tests
# ---------------------------------------------------------------------------- #

@testitem "validate_solution: consistent solve passes" begin
    using DirectTrajOpt
    using LinearAlgebra

    # A converged 2-level state-transfer via the standard SmoothPulseProblem
    # machinery (BilinearIntegrator). Because the optimizer's discretization is
    # fine (N=50) the discrete solution and an honest bare-propagator re-rollout
    # of the same pulse agree closely, so the result must pass.
    T = 10.0
    N = 50
    sys = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]

    pulse = ZeroOrderPulse(randn(2, N), collect(range(0.0, T, length = N)))
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)
    qcp = SmoothPulseProblem(qtraj, N; Q = 50.0, R = 1e-3)

    solve!(qcp; max_iter = 100, verbose = false, print_level = 1)

    result = validate_solution(qcp; atol = 1e-2)

    @test result isa ValidationResult
    @test 0.0 <= result.F_reported <= 1.0
    @test 0.0 <= result.F_rerolled <= 1.0
    @test result.ΔF == abs(result.F_reported - result.F_rerolled)
    @test result.ΔF ≤ 1e-2
    @test result.pass

    # Sanity: the solve reached high fidelity AND the honest re-rollout confirms
    # it — so the pass certifies a genuinely good solution, not a self-consistent
    # lie. Both numbers must be high.
    @test result.F_reported > 0.9
    @test result.F_rerolled > 0.9
end

@testitem "validate_solution: tampered trajectory fails" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    # Solve a small problem, then TAMPER with the stored final state so the
    # reported fidelity is a lie while the pulse (hence the re-rollout) is
    # unchanged. This synthesizes exactly the optimizer-vs-reality mismatch the
    # detector exists to catch.
    levels = 2
    H_drift = zeros(ComplexF64, levels, levels)
    σx = ComplexF64[0 1; 1 0]
    N = 21
    T = 5.0

    sys = QuantumSystem(H_drift, [σx], [(-2.0, 2.0)]; time_dependent = false)
    ψ_init = ComplexF64[1, 0]
    ψ_goal = ComplexF64[0, 1]

    times = collect(range(0, T, length = N))
    pulse = LinearSplinePulse(zeros(1, N), times)
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

    traj = NamedTrajectory(qtraj, N)
    obj = KetInfidelityObjective(:ψ̃, traj; Q = 100.0)
    integrator = BilinearIntegrator(qtraj, N)
    prob = DirectTrajOptProblem(traj, obj, integrator)
    qcp = QuantumControlProblem(qtraj, prob)

    # Do NOT optimize meaningfully — start from a near-zero pulse so the *real*
    # (re-rolled) fidelity is low, then overwrite the stored final state with the
    # perfect goal so the *reported* fidelity is spuriously ~1.
    solve!(qcp; max_iter = 0, verbose = false, print_level = 1)

    # Tamper: force the final discrete state column to be exactly the goal ket.
    ψ̃_goal = qcp.prob.trajectory.goal[:ψ̃]
    qcp.prob.trajectory.data[qcp.prob.trajectory.components.ψ̃, end] .= ψ̃_goal

    result = validate_solution(qcp; atol = 1e-4, verbose = false)

    # Reported (from the tampered state) should be ~1; re-rolled (from the real,
    # near-zero pulse) should be far lower ⇒ big ΔF ⇒ pass == false.
    @test result.F_reported > 0.99
    @test result.F_rerolled < 0.5
    @test result.ΔF > 1e-4
    @test !result.pass
end

@testitem "validate_solution: EmbeddedOperator goal uses Pedersen on both sides" begin
    using DirectTrajOpt
    using LinearAlgebra

    # 2-level X gate embedded in a 3-level system, solved with SmoothPulseProblem.
    # Both F_reported and F_rerolled must use the Pedersen average-gate formula on
    # the computational subspace, so a self-consistent solve validates cleanly and
    # both fidelities are high.
    levels = 3
    a = annihilate(levels)
    H_drift = 2π * (-0.2 / 2) * (a' * a' * a * a)   # anharmonicity
    H_drives = [2π * (a + a'), 2π * 1im * (a - a')]
    sys = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

    σx = ComplexF64[0 1; 1 0]
    U_goal = EmbeddedOperator(σx, sys)  # embed X on the |0⟩,|1⟩ subspace

    N = 50
    T = 10.0
    pulse = ZeroOrderPulse(0.1 * randn(2, N), collect(range(0.0, T, length = N)))
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)
    qcp = SmoothPulseProblem(qtraj, N; Q = 100.0, R = 1e-3)

    solve!(qcp; max_iter = 200, verbose = false, print_level = 1)

    result = validate_solution(qcp; atol = 1e-2)

    @test result isa ValidationResult
    # Pedersen average-gate fidelity is not clamped, so allow tiny fp roundoff
    # above 1.0 (same behaviour as the objective / `fidelity(::UnitaryTrajectory)`).
    @test 0.0 <= result.F_reported <= 1.0 + 1e-6
    @test 0.0 <= result.F_rerolled <= 1.0 + 1e-6
    # Reported and re-rolled use the same Pedersen convention → should agree.
    @test result.pass
    # Both should be high on a converged embedded-gate solve.
    @test result.F_reported > 0.9
    @test result.F_rerolled > 0.9
end

end  # module Validation
