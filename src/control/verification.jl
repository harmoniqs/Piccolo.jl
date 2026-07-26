module Verification

using LinearAlgebra          # Diagonal, in the unitary free-phase goal rotation
using NamedTrajectories
using TestItems
using ...Quantum

using ..QuantumControlProblems: QuantumControlProblem, fidelity, stored_phases
# `rollout_divergence` is imported for its docstring cross-reference, not called here: Documenter
# resolves `@ref` against the enclosing module's scope, so a `[`rollout_divergence`](@ref)` in
# `verify`'s docstring fails the docs build ([:cross_references]) unless the name is in scope.
using ..QuantumControlProblems: rollout_divergence
using ..QuantumObjectives: unitary_fidelity_loss, ket_fidelity_loss

export verify

# Fidelity of the OPTIMIZER's own terminal state — the collocation state it minimized
# against — rather than of the rolled-out ODE solution. This lives here, not in
# `problems.jl`, because the `*_fidelity_loss` functions are defined in `objectives.jl`,
# which is included after it.
#
# Returns `nothing` where no comparable definition is available (density, MultiDensity, sampling,
# and free-phase cases needing a factorization we are not given): `nothing` means NOT COMPARABLE,
# never "they agree". MultiKet IS covered, via its coherent-fidelity definition below.
function _optimizer_side_fidelity(qcp::QuantumControlProblem, phases)
    # MultiKet is not a single-state case: its objective is a coherent sum across sub-states.
    qcp.qtraj isa MultiKetTrajectory &&
        return _multiket_optimizer_side_fidelity(qcp, phases)

    traj = qcp.prob.trajectory
    sname = state_name(qcp.qtraj)
    haskey(traj.components, sname) || return nothing
    x_end = traj[sname][:, end]
    return _optimizer_side_fidelity(qcp.qtraj, x_end, phases)
end

function _optimizer_side_fidelity(qtraj::UnitaryTrajectory, Ũ⃗, phases)
    goal = qtraj.goal
    if isnothing(phases)
        return Float64(unitary_fidelity_loss(Ũ⃗, goal))
    end
    # Free-phase goals rotate on the computational subspace, which is only defined for an
    # EmbeddedOperator. Without one there is nothing to rotate, so decline rather than guess.
    goal isa EmbeddedOperator || return nothing
    U_sub = unembed(goal)
    n_sub = size(U_sub, 1)
    n_q = length(phases)
    # Qubit/binary convention — the same one `fidelity(::UnitaryTrajectory; phases)` uses on
    # the unembedded subspace. NOT `number_operator_phase_diag`, which is the ket convention.
    diag = map(1:n_sub) do i
        bits = i - 1
        φ = sum(phases[j] for j = 1:n_q if (bits >> (n_q - j)) & 1 == 1; init = 0.0)
        return exp(im * φ)
    end
    rotated = EmbeddedOperator(Diagonal(diag) * U_sub, goal.subspace, goal.subsystem_levels)
    return Float64(unitary_fidelity_loss(Ũ⃗, rotated))
end

function _optimizer_side_fidelity(qtraj::KetTrajectory, ψ̃, phases)
    goal = qtraj.goal
    isnothing(phases) && return Float64(ket_fidelity_loss(ψ̃, goal))
    length(phases) == 1 || return nothing   # factorization unknown; see fidelity(::KetTrajectory)
    diag = number_operator_phase_diag(phases, [length(goal)])
    return Float64(ket_fidelity_loss(ψ̃, diag .* goal))
end

_optimizer_side_fidelity(::AbstractQuantumTrajectory, _, _) = nothing

# MultiKet needs its own path: the objective is the COHERENT fidelity across all transfers,
# F = |1/n Σᵢ ⟨ψᵢ_goal|ψᵢ⟩|², not a per-state loss, and the states live under one component per
# sub-trajectory (`:ψ̃1`, `:ψ̃2`, …) rather than a single one. Mirrors
# `fidelity(::MultiKetTrajectory)` so both sides of `verify` use the same definition.
#
# This case is why `verify` can replace the hand-rolled optimizer-vs-rollout comparison in
# `fluxonium-2q/scripts/probe/rollout_fidelity_check.jl`, which is a MultiKet problem.
function _multiket_optimizer_side_fidelity(qcp::QuantumControlProblem, phases)
    qtraj = qcp.qtraj
    traj = qcp.prob.trajectory
    names = state_names(qtraj)
    n = length(qtraj.goals)
    length(names) == n || return nothing
    all(nm -> haskey(traj.components, nm), names) || return nothing

    goals = if isnothing(phases)
        qtraj.goals
    else
        # Free-phase MultiKet goals rotate under the same qubit/binary convention that
        # `fidelity(::MultiKetTrajectory)` applies, which needs subsystem_levels we do not have
        # here. Decline rather than guess at the factorization.
        return nothing
    end

    overlap = sum(goals[i]' * iso_to_ket(traj[names[i]][:, end]) for i = 1:n)
    return Float64(abs2(overlap / n))
end

"""
    verify(qcp::QuantumControlProblem; kwargs...) -> (; F_optimizer, F_rollout, Δ, phases)

Recompute fidelity two ways and return both plus their gap:

- `F_rollout` — an independent ODE re-rollout of the extracted pulse. **This is the physical
  number; put it in artifacts, papers and catalogs.**
- `F_optimizer` — the optimizer's own collocation terminal state. Covers unitary, ket and MultiKet.
  `nothing` means *not comparable*, never *they agree*: density, MultiDensity, sampling, free-phase
  MultiKet, or a unitary free-phase problem without an `EmbeddedOperator` goal.
- `Δ` — `abs(F_optimizer - F_rollout)`, or `nothing` when `F_optimizer` is.
- `phases` — the phases actually applied, so callers can record the convention used.

Both sides get the **same** phase treatment, so the gap is never an artifact of one side being
φ-aware. A large `Δ` means the optimizer converged against a model the pulse does not realize.
Compare [`rollout_divergence`](@ref), which compares terminal *states* and so also catches cases
where two different states share a fidelity.

# Example
```julia
solve!(qcp)
v = verify(qcp)
v.Δ > 1e-6 && @warn "optimizer and rollout disagree" v.F_optimizer v.F_rollout
```
"""
function verify(qcp::QuantumControlProblem; phases = nothing, kwargs...)
    φ = isnothing(phases) ? stored_phases(qcp.prob.trajectory) : phases
    if !isnothing(φ) && isempty(φ)
        error(
            "this problem declares free-phase globals (φ_*) but stores no phase values; " *
            "`verify` will not silently report fixed-phase numbers. Pass `phases = ...`.",
        )
    end

    F_rollout = if isnothing(φ)
        Float64(fidelity(qcp.qtraj; kwargs...))
    else
        Float64(fidelity(qcp.qtraj; phases = φ, kwargs...))
    end

    F_optimizer = _optimizer_side_fidelity(qcp, φ)
    Δ = isnothing(F_optimizer) ? nothing : abs(F_optimizer - F_rollout)

    return (; F_optimizer, F_rollout, Δ, phases = φ)
end

# ============================================================================= #
# Tests
# ============================================================================= #

@testitem "verify: reports both fidelities and their gap" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

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
    sync_trajectory!(qcp; check_divergence = false)

    v = verify(qcp)
    @test haskey(v, :F_optimizer) && haskey(v, :F_rollout) && haskey(v, :Δ)
    @test 0 ≤ v.F_rollout ≤ 1
    @test v.Δ ≈ abs(v.F_optimizer - v.F_rollout)
    @test isnothing(v.phases)          # no φ globals on this problem

    # Zero controls: collocation and rollout agree, so the gap is ~0.
    @test v.Δ < 1e-8
end

@testitem "verify: gap is large exactly when the collocation state is wrong" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

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
    sync_trajectory!(qcp; check_divergence = false)

    # Claim the collocation solution reached the goal while the controls stay at zero, i.e.
    # exactly the "optimizer converged against a model the pulse does not realize" failure.
    qcp.prob.trajectory.ψ̃[:, end] .= Float64[0, 1, 0, 0]

    v = verify(qcp)
    @test v.F_optimizer > 0.99      # the optimizer believes it succeeded
    @test v.F_rollout < 0.01        # the pulse does nothing of the sort
    @test v.Δ > 0.9
end

@testitem "verify covers MultiKet on both sides" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    # MultiKet's objective is the COHERENT fidelity across transfers, so `verify` needs its own
    # optimizer-side definition. Without it F_optimizer was `nothing` for exactly the case the
    # hand-rolled check in fluxonium-2q/scripts/probe/rollout_fidelity_check.jl exists to cover,
    # so `verify` could not replace it.
    σx = ComplexF64[0 1; 1 0]
    sys = QuantumSystem(0.01 * ComplexF64[1 0; 0 -1], [σx], [1.0])
    N, T = 21, 10.0
    times = collect(range(0.0, T, length = N))
    pulse = LinearSplinePulse(0.1 * randn(1, N), times)
    ψ0, ψ1 = ComplexF64[1.0, 0.0], ComplexF64[0.0, 1.0]

    qtraj = MultiKetTrajectory(sys, pulse, [ψ0, ψ1], [ψ1, ψ0])
    qcp = SplinePulseProblem(qtraj, N; Q = 100.0, R = 1e-2)
    sync_trajectory!(qcp; check_divergence = false)

    v = verify(qcp)
    @test v.F_optimizer isa Float64        # NOT `nothing` any more
    @test v.F_rollout isa Float64
    @test 0 ≤ v.F_optimizer ≤ 1
    @test v.Δ ≈ abs(v.F_optimizer - v.F_rollout)

    # Corrupting one sub-state's collocation value must move the optimizer-side number, proving it
    # actually reads the collocation states rather than the rollout.
    traj = get_trajectory(qcp)
    before = verify(qcp).F_optimizer
    traj.ψ̃1[:, end] .*= -1
    @test verify(qcp).F_optimizer != before
end

@testitem "verify handles a UNITARY free-phase problem (the catalog's actual case)" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    # The unitary free-phase branch builds a phase-rotated EmbeddedOperator using the
    # QUBIT/BINARY convention — a different expression from the ket path's number-operator form.
    # It was untested, and it is the production case: all 9 free_phase catalog entries are
    # unitary. An untested convention is how incompatible free-phase definitions accumulate.
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

    traj.global_data[traj.global_components[:φ_1]] .= [0.6]
    sync_trajectory!(qcp; check_divergence = false)

    v = verify(qcp)
    @test v.F_optimizer isa Float64        # the EmbeddedOperator rotation path is reached
    @test v.F_rollout isa Float64
    @test v.phases == [0.6]
    @test v.Δ ≈ abs(v.F_optimizer - v.F_rollout)

    # The rollout side must agree with asking for that φ directly — i.e. verify is not applying
    # some other convention on one side.
    @test v.F_rollout ≈ fidelity(qcp.qtraj; phases = [0.6])

    # And φ must actually matter on BOTH sides, or the rotation is inert.
    v0 = verify(qcp; phases = [0.0])
    @test !isapprox(v0.F_optimizer, v.F_optimizer; atol = 1e-10)
    @test !isapprox(v0.F_rollout, v.F_rollout; atol = 1e-10)
end

end # module Verification
