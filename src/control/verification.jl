module Verification

using NamedTrajectories
using TestItems
using ...Quantum

using ..QuantumControlProblems: QuantumControlProblem, fidelity, stored_phases
using ..QuantumObjectives: unitary_fidelity_loss, ket_fidelity_loss

export verify

# Fidelity of the OPTIMIZER's own terminal state — the collocation state it minimized
# against — rather than of the rolled-out ODE solution. This lives here, not in
# `problems.jl`, because the `*_fidelity_loss` functions are defined in `objectives.jl`,
# which is included after it.
#
# Returns `nothing` for trajectory types with no single-state fidelity loss (density,
# multi-state, sampling): `nothing` means NOT COMPARABLE, never "they agree".
function _optimizer_side_fidelity(qcp::QuantumControlProblem, phases)
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

"""
    verify(qcp::QuantumControlProblem; kwargs...) -> (; F_optimizer, F_rollout, Δ, phases)

Recompute fidelity two ways and return both plus their gap:

- `F_rollout` — an independent ODE re-rollout of the extracted pulse. **This is the physical
  number. Put it in artifacts, papers and catalogs.**
- `F_optimizer` — the optimizer's own collocation terminal state. `nothing` when the trajectory
  type has no single-state fidelity loss (density / multi-state / sampling); `nothing` means
  *not comparable*, never *they agree*.
- `Δ` — `abs(F_optimizer - F_rollout)`, or `nothing` when `F_optimizer` is.
- `phases` — the free phases actually applied, so the caller can record which convention
  produced the numbers.

Both sides use the **same** phase treatment, so the comparison is never confounded by one side
being φ-aware and the other not.

A large `Δ` means the optimizer converged against a model the pulse does not realize — usually
a pulse/integrator mismatch. Compare [`rollout_divergence`](@ref), which compares terminal
*states* instead and so catches cases where two different states happen to share a fidelity.

This replaces hand-rolled copies of the same check that had accumulated in demo scripts.

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

end # module Verification
