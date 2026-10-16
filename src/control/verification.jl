module Verification

using LinearAlgebra          # Diagonal, in the unitary free-phase goal rotation
using NamedTrajectories
using TestItems
using ...Quantum

using ...Quantum: extract_pulse
using ...Quantum.Rollouts: rollout

using ..QuantumControlProblems:
    QuantumControlProblem,
    AbstractQuantumControlProblem,
    fidelity,
    stored_phases,
    get_trajectory
using ..QuantumControlProblems: rollout_divergence, ROLLOUT_DIVERGENCE_RTOL
# `sync_trajectory!` is imported for its docstring cross-reference, not called here:
# Documenter resolves `@ref` against the enclosing module's scope, so a
# `[`sync_trajectory!`](@ref)` in `verify`'s docstring fails the docs build
# ([:cross_references]) unless the name is in scope. Same reason `rollout_divergence` was
# imported before it had a caller here.
using ..QuantumControlProblems: sync_trajectory!
using ..QuantumObjectives: unitary_fidelity_loss, ket_fidelity_loss

export verify, VerificationReport

# Fidelity of the OPTIMIZER's own terminal state — the collocation state it minimized
# against — rather than of the rolled-out ODE solution. This lives here, not in
# `problems.jl`, because the `*_fidelity_loss` functions are defined in `objectives.jl`,
# which is included after it.
#
# Returns `nothing` where no comparable definition is available (density, MultiDensity, sampling,
# and free-phase cases needing a factorization we are not given): `nothing` means NOT COMPARABLE,
# never "they agree". MultiKet IS covered, via its coherent-fidelity definition below.
function _optimizer_side_fidelity(qcp::AbstractQuantumControlProblem, phases)
    # MultiKet is not a single-state case: its objective is a coherent sum across sub-states.
    qcp.qtraj isa MultiKetTrajectory &&
        return _multiket_optimizer_side_fidelity(qcp, phases)

    # SamplingTrajectory: per-member collocation readback, weighted aggregate
    qcp.qtraj isa SamplingTrajectory &&
        return _sampling_optimizer_side_fidelity(qcp, phases)

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

# SamplingTrajectory: per-member collocation readback, weighted with the sampling weights.
# Density / MultiDensity bases return `nothing` — no comparable definition.
function _sampling_optimizer_side_fidelity(qcp::AbstractQuantumControlProblem, phases)
    sampling_qtraj = qcp.qtraj
    base = sampling_qtraj.base_trajectory
    traj = qcp.prob.trajectory
    snames = state_names(sampling_qtraj)
    member_states = sampling_member_states(sampling_qtraj)
    weights = sampling_qtraj.weights

    # Density bases: no comparable definition
    base isa DensityTrajectory && return nothing
    base isa MultiDensityTrajectory && return nothing

    # Single-state bases (unitary, ket): one state component per member
    if base isa UnitaryTrajectory || base isa KetTrajectory
        member_fids = Float64[]
        for (i, sname) in enumerate(snames)
            haskey(traj.components, sname) || return nothing
            x_end = traj[sname][:, end]
            fid = _optimizer_side_fidelity(base, x_end, phases)
            isnothing(fid) && return nothing
            push!(member_fids, fid)
        end
        return sum(w * f for (w, f) in zip(weights, member_fids))
    end

    # MultiKet base: per-member coherent fidelity over the member's ket sub-states
    if base isa MultiKetTrajectory
        isnothing(phases) || return nothing  # free-phase MultiKet not supported
        member_fids = Float64[]
        for (i, mstates) in enumerate(member_states)
            # mstates is a Vector{Symbol} for multi-state bases
            n_sub = length(mstates)
            all(nm -> haskey(traj.components, nm), mstates) || return nothing
            overlap =
                sum(base.goals[j]' * iso_to_ket(traj[mstates[j]][:, end]) for j = 1:n_sub)
            push!(member_fids, Float64(abs2(overlap / n_sub)))
        end
        return sum(w * f for (w, f) in zip(weights, member_fids))
    end

    return nothing
end

# MultiKet needs its own path: the objective is the COHERENT fidelity across all transfers,
# F = |1/n Σᵢ ⟨ψᵢ_goal|ψᵢ⟩|², not a per-state loss, and the states live under one component per
# sub-trajectory (`:ψ̃1`, `:ψ̃2`, …) rather than a single one. Mirrors
# `fidelity(::MultiKetTrajectory)` so both sides of `verify` use the same definition.
#
# This case is why `verify` can replace the hand-rolled optimizer-vs-rollout comparison in
# `fluxonium-2q/scripts/probe/rollout_fidelity_check.jl`, which is a MultiKet problem.
function _multiket_optimizer_side_fidelity(qcp::AbstractQuantumControlProblem, phases)
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
    VerificationReport

Outcome of [`verify`](@ref): both fidelities, both disagreement measures, the phase
convention used, and a verdict.

Every numeric field except `F_rollout` is `Union{Nothing,Float64}`, and the absence is
load-bearing — `nothing` means *the comparison could not be made*, never *it agreed*. That
is why the verdict is a three-valued `status` rather than a `Bool`: for a problem where
neither comparison applies, `pass = true` would certify a check that never ran, and
`pass = false` would report a disagreement that was never measured. The pulse catalog's
free-phase entries are all unitary free-phase, one of those cases, so this is the common
path rather than a corner.

# Fields
- `F_rollout::Float64`: fidelity of an independent ODE re-rollout. **The physical number
  — this is what belongs in artifacts, papers and catalogs.**
- `F_optimizer::Union{Nothing,Float64}`: fidelity of the optimizer's own collocation
  terminal state. `nothing` for density, MultiDensity, sampling, free-phase MultiKet, or a
  unitary free-phase problem without an `EmbeddedOperator` goal.
- `ΔF::Union{Nothing,Float64}`: `abs(F_optimizer - F_rollout)`, compared against `atol`.
- `divergence::Union{Nothing,Float64}`: [`rollout_divergence`](@ref) — relative
  disagreement of the terminal *states*, compared against `rtol`. `nothing` when not
  applicable. Strictly stronger than `ΔF` at detecting a mismatch, since two different
  states can share a fidelity; it just yields no fidelity to publish.
- `phases::Union{Nothing,Vector{Float64}}`: the phases actually applied to both sides.
- `status::Symbol`: `:pass` (every applicable comparison was within tolerance), `:fail`
  (at least one exceeded it), or `:unverified` (neither applied — nothing was checked).
- `atol::Float64`, `rtol::Float64`: the thresholds used.
"""
struct VerificationReport
    F_rollout::Float64
    F_optimizer::Union{Nothing,Float64}
    ΔF::Union{Nothing,Float64}
    divergence::Union{Nothing,Float64}
    phases::Union{Nothing,Vector{Float64}}
    status::Symbol
    atol::Float64
    rtol::Float64
end

# The verdict rule, factored out of `verify` so it is testable on its own: reaching every
# branch through a real solve would mean manufacturing a problem for which each comparison
# individually declines, and the `:unverified` branch is the one most worth pinning.
#
# `:pass` means every comparison that COULD be made was, and agreed. When neither could,
# the answer is `:unverified` — `:pass` there would certify a check that never ran, which
# is worse than running no check at all, and `:fail` would report a disagreement nobody
# measured.
function _verdict(ΔF, divergence, atol, rtol)
    isnothing(ΔF) && isnothing(divergence) && return :unverified
    exceeded =
        (!isnothing(ΔF) && ΔF > atol) || (!isnothing(divergence) && divergence > rtol)
    return exceeded ? :fail : :pass
end

# `nothing` prints as an explicit reason, never as a blank or a zero: a report that is
# silent about an unmade comparison reads as a clean bill of health.
_fmt(x::Nothing, why::AbstractString) = "— ($why)"
_fmt(x::Real, ::AbstractString) = string(x)

function Base.show(io::IO, ::MIME"text/plain", r::VerificationReport)
    println(io, "VerificationReport ($(uppercase(string(r.status))))")
    println(io, "  F_rollout   (physical): ", r.F_rollout)
    println(io, "  F_optimizer           : ", _fmt(r.F_optimizer, "not comparable"))
    println(io, "  ΔF          (atol ", r.atol, "): ", _fmt(r.ΔF, "not comparable"))
    println(io, "  divergence  (rtol ", r.rtol, "): ", _fmt(r.divergence, "not applicable"))
    print(io, "  phases                : ", r.phases === nothing ? "none" : r.phases)
end

Base.show(io::IO, r::VerificationReport) = print(
    io,
    "VerificationReport(",
    r.status,
    ", F_rollout=",
    r.F_rollout,
    ", ΔF=",
    r.ΔF === nothing ? "nothing" : string(r.ΔF),
    ", divergence=",
    r.divergence === nothing ? "nothing" : string(r.divergence),
    ")",
)

"""
    verify(qcp::QuantumControlProblem; kwargs...) -> VerificationReport

Self-certify a solved problem: recompute fidelity two ways, compare the terminal states as
well, and return a [`VerificationReport`](@ref) carrying both numbers, both disagreements,
the phase convention used, and a verdict.

The two comparisons answer different questions and neither subsumes the other, so `verify`
runs both and fails if *either* exceeds tolerance:

- `ΔF` (against `atol`) compares **fidelities**, so it exposes "the optimizer thinks it
  converged" directly and yields a number you can publish.
- `divergence` (against `rtol`) compares terminal **states** via
  [`rollout_divergence`](@ref) — strictly stronger at detecting a mismatch, since two
  different states can share a fidelity, but it yields no fidelity.

Both fidelity sides get the **same** phase treatment, so `ΔF` is never an artifact of one
side being φ-aware. A large `ΔF` or `divergence` means the optimizer converged against a
model the pulse does not realize; only `F_rollout` is physical.

Call after [`sync_trajectory!`](@ref) (or `solve!` with `sync = true`) — before that,
`qcp.qtraj` holds a stale rollout and every number here is meaningless.

# Keyword arguments
- `atol::Real = 1e-4`: tolerance on `ΔF`.
- `rtol::Real = ROLLOUT_DIVERGENCE_RTOL[]`: tolerance on `divergence`.
- `algorithm = nothing`: when given, re-extract the pulse and roll it out **fresh** with
  this ODE algorithm instead of reading the stored `qcp.qtraj` solution. Passing a
  different integrator *family* than the solve used is what makes `F_rollout` genuinely
  independent rather than a re-read of the solve's own discretization.
- `n_save::Int = 101`: save points for that fresh rollout (ignored without `algorithm`).
- `phases = nothing`: phase convention; defaults to the problem's stored `φ_*` globals.
- `verbose::Bool = true`: `@warn` with both numbers when the verdict is `:fail`.

# Example
```julia
solve!(qcp)
r = verify(qcp)
r.status === :pass || @warn "solve did not self-certify" r
# genuinely independent check — a different integrator family than the solve used:
r2 = verify(qcp; algorithm = MagnusAdapt4())
```
"""
function verify(
    qcp::AbstractQuantumControlProblem;
    atol::Real = 1e-4,
    rtol::Real = ROLLOUT_DIVERGENCE_RTOL[],
    algorithm = nothing,
    n_save::Int = 101,
    phases = nothing,
    verbose::Bool = true,
    kwargs...,
)
    φ = isnothing(phases) ? stored_phases(qcp.prob.trajectory) : phases
    if !isnothing(φ) && isempty(φ)
        error(
            "this problem declares free-phase globals (φ_*) but stores no phase values; " *
            "`verify` will not silently report fixed-phase numbers. Pass `phases = ...`.",
        )
    end

    # Default reads the stored rollout; `algorithm` forces a fresh propagation of the
    # extracted pulse, which is the only way to get a family different from the solve's.
    qtraj_scored = if isnothing(algorithm)
        qcp.qtraj
    else
        rollout(
            qcp.qtraj,
            extract_pulse(qcp.qtraj, get_trajectory(qcp));
            algorithm = algorithm,
            n_save = n_save,
        )
    end

    F_rollout = if qtraj_scored isa SamplingTrajectory
        members = if isnothing(φ)
            fidelity(qtraj_scored; kwargs...)
        else
            fidelity(qtraj_scored; phases = φ, kwargs...)
        end
        Float64(sum(w * f for (w, f) in zip(qtraj_scored.weights, members)))
    else
        if isnothing(φ)
            Float64(fidelity(qtraj_scored; kwargs...))
        else
            Float64(fidelity(qtraj_scored; phases = φ, kwargs...))
        end
    end

    F_optimizer = _optimizer_side_fidelity(qcp, φ)
    ΔF = isnothing(F_optimizer) ? nothing : abs(F_optimizer - F_rollout)
    divergence = rollout_divergence(qcp)

    report = VerificationReport(
        F_rollout,
        F_optimizer,
        ΔF,
        divergence,
        isnothing(φ) ? nothing : Vector{Float64}(φ),
        _verdict(ΔF, divergence, atol, rtol),
        Float64(atol),
        Float64(rtol),
    )

    if verbose && report.status === :fail
        @warn "verify: the solve did not self-certify — only `F_rollout` is physical" report
    end

    return report
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
    # `verify` returns a `VerificationReport` now, not a NamedTuple, so this is
    # `hasproperty` rather than `haskey` — and it covers the two fields the struct added
    # (`divergence`, `status`) alongside the renamed `ΔF`.
    @test v isa VerificationReport
    @test all(
        f -> hasproperty(v, f),
        (:F_optimizer, :F_rollout, :ΔF, :divergence, :phases, :status, :atol, :rtol),
    )
    @test 0 ≤ v.F_rollout ≤ 1
    @test v.ΔF ≈ abs(v.F_optimizer - v.F_rollout)
    @test isnothing(v.phases)          # no φ globals on this problem

    # Zero controls: collocation and rollout agree, so the gap is ~0.
    @test v.ΔF < 1e-8
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

    v = verify(qcp; verbose = false)   # this case IS the failure; don't log it
    @test v.F_optimizer > 0.99      # the optimizer believes it succeeded
    @test v.F_rollout < 0.01        # the pulse does nothing of the sort
    @test v.ΔF > 0.9
    @test v.status === :fail        # and the verdict says so
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
    @test v.ΔF ≈ abs(v.F_optimizer - v.F_rollout)

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
    @test v.ΔF ≈ abs(v.F_optimizer - v.F_rollout)

    # The rollout side must agree with asking for that φ directly — i.e. verify is not applying
    # some other convention on one side.
    @test v.F_rollout ≈ fidelity(qcp.qtraj; phases = [0.6])

    # And φ must actually matter on BOTH sides, or the rotation is inert.
    v0 = verify(qcp; phases = [0.0])
    @test !isapprox(v0.F_optimizer, v.F_optimizer; atol = 1e-10)
    @test !isapprox(v0.F_rollout, v.F_rollout; atol = 1e-10)
end

# ============================================================================= #
# Tests: sampling trajectory verification
# ============================================================================= #

@testitem "verify: unitary sampling reports weighted agreement" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    # Small 2-qubit unitary problem with two perturbed systems
    T, N = 2.0, 11
    sys_nom = QuantumSystem(GATES[:Z], [GATES[:X]], [1.0])
    sys_var = QuantumSystem(1.05 * GATES[:Z], [GATES[:X]], [1.0])

    pulse = ZeroOrderPulse(0.2 * randn(1, N), collect(range(0.0, T, length = N)))
    qtraj = UnitaryTrajectory(sys_nom, pulse, GATES[:X])
    qcp = SmoothPulseProblem(qtraj, N; Q = 50.0, R = 1e-3)

    sampling_prob = SamplingProblem(qcp, [sys_nom, sys_var]; Q = 50.0, weights = [0.5, 0.5])
    solve!(sampling_prob; max_iter = 20, verbose = false, print_level = 0)

    r = verify(sampling_prob)
    @test r isa VerificationReport
    @test r.F_optimizer isa Float64    # NOT nothing — per-member collocation readback enabled
    @test 0.0 <= r.F_rollout <= 1.0 + 1e-8
    @test r.ΔF !== nothing
    @test r.ΔF ≈ abs(r.F_optimizer - r.F_rollout)

    # Corrupt one member's collocation state — the optimizer-side number must move
    traj = get_trajectory(sampling_prob)
    before = verify(sampling_prob).F_optimizer
    traj.Ũ⃗1[:, end] .= 0.0
    @test verify(sampling_prob).F_optimizer != before
end

@testitem "verify: ket sampling reports weighted agreement" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    T, N = 2.0, 11
    sys_nom = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    sys_var = QuantumSystem(1.05 * GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

    ψ0 = ComplexF64[1.0, 0.0]
    ψg = ComplexF64[0.0, 1.0]
    pulse = ZeroOrderPulse(0.2 * randn(2, N), collect(range(0.0, T, length = N)))
    qtraj = KetTrajectory(sys_nom, pulse, ψ0, ψg)
    qcp = SmoothPulseProblem(qtraj, N; Q = 50.0, R = 1e-3)

    sampling_prob = SamplingProblem(qcp, [sys_nom, sys_var]; Q = 50.0, weights = [0.5, 0.5])
    solve!(sampling_prob; max_iter = 20, verbose = false, print_level = 0)

    r = verify(sampling_prob)
    @test r isa VerificationReport
    @test r.F_optimizer isa Float64
    @test 0.0 <= r.F_rollout <= 1.0 + 1e-8
    @test r.ΔF !== nothing
    @test r.ΔF ≈ abs(r.F_optimizer - r.F_rollout)
end

@testitem "verify: multi-ket sampling reports weighted coherent agreement" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    T, N = 1.0, 11
    sys_nom = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    sys_var = QuantumSystem(1.05 * GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]
    pulse = ZeroOrderPulse(0.1 * randn(2, N), collect(range(0.0, T, length = N)))
    qtraj = MultiKetTrajectory(sys_nom, pulse, [ψ0, ψ1], [ψ1, ψ0])
    qcp = SmoothPulseProblem(qtraj, N; Q = 30.0, R = 1e-3)

    sampling_prob = SamplingProblem(qcp, [sys_nom, sys_var]; Q = 30.0, weights = [0.5, 0.5])
    solve!(sampling_prob; max_iter = 10, verbose = false, print_level = 0)

    r = verify(sampling_prob)
    @test r isa VerificationReport
    @test r.F_optimizer isa Float64
    @test 0.0 <= r.F_rollout <= 1.0 + 1e-8
    @test r.ΔF !== nothing
    @test r.ΔF ≈ abs(r.F_optimizer - r.F_rollout)
end

@testitem "verify: density sampling returns no comparable definition" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    T, N = 1.0, 11
    L = ComplexF64[0.0 0.1; 0.0 0.0]
    sys_nom = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    sys_var =
        OpenQuantumSystem(0.95 * PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])

    ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρg = ComplexF64[0.0 0.0; 0.0 1.0]
    pulse = ZeroOrderPulse(0.2 * randn(1, N), collect(range(0.0, T, length = N)))
    base_qtraj = DensityTrajectory(sys_nom, pulse, ρ0, ρg)

    sampling_qtraj = SamplingTrajectory(base_qtraj, [sys_nom, sys_var])
    traj = NamedTrajectory(sampling_qtraj, N)

    prob = DirectTrajOptProblem(traj, NullObjective(), AbstractIntegrator[])
    qcp = QuantumControlProblem(sampling_qtraj, prob)
    sync_trajectory!(qcp; check_divergence = false)

    r = verify(qcp)
    @test r.F_optimizer === nothing
    @test r.ΔF === nothing
    @test r.divergence === nothing
    @test r.status === :unverified
end

@testitem "verify: multi-density sampling returns no comparable definition" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    T, N = 1.0, 11
    L = ComplexF64[0.0 0.1; 0.0 0.0]
    sys_nom = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    sys_var =
        OpenQuantumSystem(0.95 * PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])

    ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρ1 = ComplexF64[0.0 0.0; 0.0 1.0]
    pulse = ZeroOrderPulse(0.2 * randn(1, N), collect(range(0.0, T, length = N)))
    base_qtraj = MultiDensityTrajectory(sys_nom, pulse, [ρ0, ρ1], [ρ1, ρ0])

    sampling_qtraj = SamplingTrajectory(base_qtraj, [sys_nom, sys_var])
    traj = NamedTrajectory(sampling_qtraj, N)

    prob = DirectTrajOptProblem(traj, NullObjective(), AbstractIntegrator[])
    qcp = QuantumControlProblem(sampling_qtraj, prob)
    sync_trajectory!(qcp; check_divergence = false)

    r = verify(qcp)
    @test r.F_optimizer === nothing
    @test r.ΔF === nothing
    @test r.divergence === nothing
    @test r.status === :unverified
end

end # module Verification

@testitem "verify: the verdict never reports success for a comparison that did not run" begin
    using Piccolo.Control.Verification: _verdict

    # `:unverified` is the branch that matters. A `pass::Bool` cannot express it, so it
    # would have to collapse into one of the other two — and `pass = true` on a problem
    # where neither comparison applied is a worse outcome than running no check at all.
    # The catalog's free-phase entries are all unitary free-phase, which is one of the
    # `F_optimizer === nothing` cases, so this is a live path.
    @test _verdict(nothing, nothing, 1e-4, 5e-3) === :unverified

    # Either comparison alone is enough to reach a real verdict.
    @test _verdict(1e-9, nothing, 1e-4, 5e-3) === :pass
    @test _verdict(nothing, 1e-9, 1e-4, 5e-3) === :pass
    @test _verdict(1e-2, nothing, 1e-4, 5e-3) === :fail
    @test _verdict(nothing, 1e-1, 1e-4, 5e-3) === :fail

    # Both present: agreement requires BOTH within tolerance, since the two measure
    # different things and neither subsumes the other.
    @test _verdict(1e-9, 1e-9, 1e-4, 5e-3) === :pass
    @test _verdict(1e-2, 1e-9, 1e-4, 5e-3) === :fail
    @test _verdict(1e-9, 1e-1, 1e-4, 5e-3) === :fail

    # The state comparison catching what the fidelity comparison misses is the whole
    # reason both run: two different states can share a fidelity.
    @test _verdict(0.0, 1e-1, 1e-4, 5e-3) === :fail

    # Boundaries are inclusive — exactly at tolerance passes.
    @test _verdict(1e-4, nothing, 1e-4, 5e-3) === :pass
    @test _verdict(nothing, 5e-3, 1e-4, 5e-3) === :pass
end

@testitem "verify: a real converged solve self-certifies" begin
    using NamedTrajectories
    using LinearAlgebra
    using OrdinaryDiffEqLinear: MagnusAdapt4

    # Ported from #252. main's other verify tests hand-build a problem and never run the
    # optimizer, so nothing exercised the verdict end-to-end on a genuine solve.
    #
    # The setup is #252's verbatim, and it has to be: a 2-drive (X,Y) system over a Z drift
    # at N = 50 converges tightly enough for BOTH comparisons to agree (measured ΔF ≈ 5e-8,
    # divergence ≈ 2e-4). A single-drive variant does NOT — it lands at ΔF ≈ 0.13 with
    # divergence ≈ 0.36, i.e. a genuine pairing/resolution failure, and would make this a
    # test of the problem rather than of `verify`.
    T, N = 10.0, 50
    sys = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    pulse = ZeroOrderPulse(randn(2, N), collect(range(0.0, T, length = N)))
    qtraj = KetTrajectory(sys, pulse, ComplexF64[1.0, 0.0], ComplexF64[0.0, 1.0])
    qcp = SmoothPulseProblem(qtraj, N; Q = 50.0, R = 1e-3)
    solve!(qcp; max_iter = 100, verbose = false, print_level = 1)

    r = verify(qcp; atol = 1e-2)
    @test r isa VerificationReport
    @test r.status === :pass
    @test 0.0 <= r.F_rollout <= 1.0 + 1e-8
    @test r.F_optimizer !== nothing
    @test r.ΔF ≈ abs(r.F_optimizer - r.F_rollout)
    # A pass must certify a genuinely good solve, not a self-consistent one: both numbers
    # high, not merely in agreement.
    @test r.F_rollout > 0.9
    @test r.F_optimizer > 0.9

    # `algorithm` forces a FRESH rollout of the extracted pulse rather than re-reading the
    # stored solution — the only way to score against a different integrator family. The
    # numbers should agree closely on a well-resolved solve, but must be computed, not copied.
    r2 = verify(qcp; atol = 1e-2, algorithm = MagnusAdapt4(), n_save = 51)
    @test r2 isa VerificationReport
    @test isapprox(r2.F_rollout, r.F_rollout; atol = 1e-3)

    # `show` must name an unmade comparison rather than leaving it blank.
    txt = sprint(show, MIME"text/plain"(), verify(qcp))
    @test occursin("F_rollout", txt)
    @test occursin("PASS", txt)
end
