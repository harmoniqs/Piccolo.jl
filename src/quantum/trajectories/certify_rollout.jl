# ============================================================================ #
# Post-solve rollout certification
# ============================================================================ #
#
# Direct collocation certifies fidelity on the *stored terminal decision
# variable*. A coarse collocation can therefore be gamed: the optimizer finds
# controls whose discrete defects are feasible while the *continuous* evolution
# of the same physical pulse misses the target (observed on B-spline min-time:
# certified F ≥ 0.999 against a true rollout F ≈ 0.5 — see
# stanford-bosonics MINTIME_COLLOCATION_FINDING.md). `certify_rollout`
# re-propagates the returned pulse on a fine grid *independent of the
# collocation* (tight-tolerance adaptive ODE of the SAME physical pulse, as
# rebuilt by `extract_pulse` over the actual ∑Δt) and compares. Minimum-time
# results are certified iff |F_certified − F_collocation| ≤ tol (default 1e-3);
# acceptance must gate on this, never on the solver's primal infeasibility.

export certify_rollout

# ---------------------------------------------------------------------------- #
# Free-phase globals stored in the NamedTrajectory
# ---------------------------------------------------------------------------- #

# Names of the per-subsystem phase globals (φ_1, φ_2, …), sorted for a
# deterministic ordering — same detection used by the problem templates.
function _free_phase_global_names(traj::NamedTrajectory)
    θ_names = Symbol[
        name for name in keys(traj.global_components) if startswith(string(name), "φ_")
    ]
    return sort!(θ_names)
end

_free_phase_global_values(traj::NamedTrajectory, θ_names::AbstractVector{Symbol}) =
    Float64[traj.global_data[traj.global_components[name]][1] for name in θ_names]

# Number-operator phase diagonal e^{iθ·n̂} per subsystem: basis element
# |s₁,…,sₙ⟩ acquires phase exp(i Σⱼ sⱼ θⱼ). Must match the convention of
# `ProblemTemplates._make_free_phase_ket_goals` (control layer), which builds
# the free-phase objectives/constraints the solve optimized against —
# duplicated here (Quantum layer) because Control depends on Quantum, not
# vice versa. Type-generic in θ for ForwardDiff (phase optimization below).
function _number_op_phase_diag(θ, subsystem_levels::Vector{Int})
    dim = prod(subsystem_levels)
    n_subsystems = length(subsystem_levels)
    return map(0:(dim-1)) do idx
        phase = zero(eltype(θ))
        remaining = idx
        for j = 1:n_subsystems
            stride = prod(subsystem_levels[k] for k = (j+1):n_subsystems; init = 1)
            sj = remaining ÷ stride
            remaining = remaining % stride
            phase += sj * θ[j]
        end
        return exp(im * phase)
    end
end

# Per-qubit binary Z-phase diagonal on a computational subspace — the unitary
# free-phase convention (`ProblemTemplates._make_free_phase_goal`).
function _qubit_phase_diag(θ::AbstractVector{<:Real}, n_sub::Int)
    n_qubits = length(θ)
    return map(1:n_sub) do i
        bits = i - 1
        phase = sum(
            θ[j] for j = 1:n_qubits if (bits >> (n_qubits - j)) & 1 == 1;
            init = 0.0,
        )
        return exp(im * phase)
    end
end

# Maximize a smooth scalar F(θ) over the (low-dimensional) phase vector,
# warm-started at the solver's phases: damped Newton ascent with a gradient
# fallback and backtracking. The stored phases are already near-optimal after
# a converged solve; this is a refinement so the certified fidelity is the
# *phase-optimized* one, not an artifact of slightly-stale phases.
function _maximize_over_phases(
    F::Function,
    θ₀::AbstractVector{<:Real};
    max_iter::Int = 50,
    gtol::Float64 = 1e-12,
)
    θ = collect(Float64, θ₀)
    Fθ = F(θ)
    for _ = 1:max_iter
        g = ForwardDiff.gradient(F, θ)
        norm(g) < gtol && break
        Δ = try
            -(ForwardDiff.hessian(F, θ) \ g)   # Newton step (ascent near a max)
        catch
            copy(g)
        end
        if !all(isfinite, Δ) || dot(Δ, g) <= 0
            Δ = copy(g)                        # fall back to plain ascent
        end
        α = 1.0
        improved = false
        while α >= 1e-8
            F_trial = F(θ .+ α .* Δ)
            if F_trial > Fθ
                θ .+= α .* Δ
                Fθ = F_trial
                improved = true
                break
            end
            α /= 2
        end
        improved || break
    end
    return θ, Fθ
end

# ---------------------------------------------------------------------------- #
# Collocation-side fidelity — evaluated on the stored terminal decision
# variable, exactly what the solver's final-fidelity floor constraint sees.
# ---------------------------------------------------------------------------- #

function _collocation_fidelity(
    qtraj::KetTrajectory,
    traj::NamedTrajectory;
    subsystem_levels::Union{Nothing,Vector{Int}} = nothing,
)
    ψ_N = iso_to_ket(collect(traj[traj.N][state_name(qtraj)]))
    θ_names = _free_phase_global_names(traj)
    isempty(θ_names) && return abs2(dot(qtraj.goal, ψ_N))
    isnothing(subsystem_levels) && throw(
        ArgumentError(
            "free-phase globals ($(join(θ_names, ", "))) present but " *
            "subsystem_levels is missing — the free-phase ket fidelity cannot " *
            "be evaluated without the per-subsystem level structure. Pass " *
            "subsystem_levels (e.g. subsystem_levels = [3, 20]).",
        ),
    )
    θ = _free_phase_global_values(traj, θ_names)
    goal_θ = _number_op_phase_diag(θ, subsystem_levels) .* qtraj.goal
    return abs2(dot(goal_θ, ψ_N))
end

function _collocation_fidelity(
    qtraj::MultiKetTrajectory,
    traj::NamedTrajectory;
    subsystem_levels::Union{Nothing,Vector{Int}} = nothing,
)
    snames = state_names(qtraj)
    n = length(snames)
    ψs = [iso_to_ket(collect(traj[traj.N][name])) for name in snames]
    θ_names = _free_phase_global_names(traj)
    goals = if isempty(θ_names)
        qtraj.goals
    else
        isnothing(subsystem_levels) && throw(
            ArgumentError(
                "free-phase globals ($(join(θ_names, ", "))) present but " *
                "subsystem_levels is missing — pass subsystem_levels.",
            ),
        )
        θ = _free_phase_global_values(traj, θ_names)
        d = _number_op_phase_diag(θ, subsystem_levels)
        [d .* g for g in qtraj.goals]
    end
    return abs2(sum(dot(goals[i], ψs[i]) for i = 1:n) / n)
end

function _collocation_fidelity(
    qtraj::UnitaryTrajectory,
    traj::NamedTrajectory;
    subsystem_levels::Union{Nothing,Vector{Int}} = nothing,
)
    U_final = iso_vec_to_operator(collect(traj[traj.N][state_name(qtraj)]))
    goal = qtraj.goal
    θ_names = _free_phase_global_names(traj)
    if goal isa EmbeddedOperator
        # Pedersen average gate fidelity on the computational subspace —
        # consistent with `unitary_fidelity_loss(Ũ⃗, ::EmbeddedOperator)` and
        # `fidelity(::UnitaryTrajectory; phases)`.
        U_base = unembed(goal)
        U_goal_sub = if isempty(θ_names)
            U_base
        else
            θ = _free_phase_global_values(traj, θ_names)
            Diagonal(_qubit_phase_diag(θ, size(U_base, 1))) * U_base
        end
        U_sub = U_final[goal.subspace, goal.subspace]
        n = length(goal.subspace)
        M = U_goal_sub' * U_sub
        return 1 / (n * (n + 1)) * (abs(tr(M' * M)) + abs2(tr(M)))
    else
        return unitary_fidelity(U_final, goal)
    end
end

_collocation_fidelity(qtraj::AbstractQuantumTrajectory, ::NamedTrajectory; kwargs...) =
    throw(
        ArgumentError(
            "certify_rollout is not implemented for $(nameof(typeof(qtraj))); " *
            "supported: KetTrajectory, MultiKetTrajectory, UnitaryTrajectory.",
        ),
    )

# ---------------------------------------------------------------------------- #
# Certified-side fidelity — evaluated on the fine-rollout terminal state of
# the SAME physical pulse (phase-optimized for free-phase kets).
# ---------------------------------------------------------------------------- #

function _certified_fidelity(
    qtraj_fine::KetTrajectory,
    traj::NamedTrajectory;
    subsystem_levels::Union{Nothing,Vector{Int}} = nothing,
    optimize_phases::Bool = true,
)
    ψ_fin = qtraj_fine.solution.u[end]
    θ_names = _free_phase_global_names(traj)
    isempty(θ_names) && return abs2(dot(qtraj_fine.goal, ψ_fin))
    isnothing(subsystem_levels) &&
        throw(ArgumentError("free-phase globals present but subsystem_levels is missing"))
    θ₀ = _free_phase_global_values(traj, θ_names)
    F =
        θ -> abs2(
            dot(_number_op_phase_diag(θ, subsystem_levels) .* qtraj_fine.goal, ψ_fin),
        )
    optimize_phases || return F(θ₀)
    _, F★ = _maximize_over_phases(F, θ₀)
    return max(F★, F(θ₀))
end

function _certified_fidelity(
    qtraj_fine::MultiKetTrajectory,
    traj::NamedTrajectory;
    subsystem_levels::Union{Nothing,Vector{Int}} = nothing,
    optimize_phases::Bool = true,
)
    n = length(qtraj_fine.goals)
    ψs = [qtraj_fine.solution.u[i].u[end] for i = 1:n]
    θ_names = _free_phase_global_names(traj)
    coherent(goals) = abs2(sum(dot(goals[i], ψs[i]) for i = 1:n) / n)
    isempty(θ_names) && return coherent(qtraj_fine.goals)
    isnothing(subsystem_levels) &&
        throw(ArgumentError("free-phase globals present but subsystem_levels is missing"))
    θ₀ = _free_phase_global_values(traj, θ_names)
    F = θ -> begin
        d = _number_op_phase_diag(θ, subsystem_levels)
        abs2(sum(dot(d .* qtraj_fine.goals[i], ψs[i]) for i = 1:n) / n)
    end
    optimize_phases || return F(θ₀)
    _, F★ = _maximize_over_phases(F, θ₀)
    return max(F★, F(θ₀))
end

function _certified_fidelity(
    qtraj_fine::UnitaryTrajectory,
    traj::NamedTrajectory;
    subsystem_levels::Union{Nothing,Vector{Int}} = nothing,
    optimize_phases::Bool = true,
)
    θ_names = _free_phase_global_names(traj)
    if qtraj_fine.goal isa EmbeddedOperator && !isempty(θ_names)
        # Evaluate at the solver's phases (the free-phase unitary convention);
        # phase refinement is currently implemented for kets only.
        θ = _free_phase_global_values(traj, θ_names)
        return Rollouts.fidelity(qtraj_fine; phases = θ)
    end
    return Rollouts.fidelity(qtraj_fine)
end

# ---------------------------------------------------------------------------- #
# Public entry point
# ---------------------------------------------------------------------------- #

"""
    certify_rollout(qtraj, traj::NamedTrajectory; kwargs...) -> NamedTuple

Post-solve verification that an optimized trajectory's certified fidelity is
physically real. Re-propagates the returned pulse on a fine grid *independent
of the collocation* — a tight-tolerance adaptive ODE rollout of the SAME
physical pulse, rebuilt by `extract_pulse` over the actual optimized duration
`∑Δt` — and compares against the fidelity the collocation certified on its
stored terminal decision variable.

A coarse collocation can be gamed by the optimizer (discrete defects feasible
while the continuous evolution misses the target), which historically produced
"minimum-time" B-spline pulses with certified F ≥ 0.999 and true rollout
F ≈ 0.5. Minimum-time results are **certified iff**
`|F_certified − F_collocation| ≤ tol` (default `1e-3`); gate acceptance on
this, never on the solver's primal infeasibility.

Fidelity conventions per trajectory type (matching the objective/constraint
the solve optimized against):
- `KetTrajectory`, free phase (φ_ globals present): number-operator
  phase-rotated goal `e^{iθ·n̂}|ψ_goal⟩`; requires `subsystem_levels`
  (`ArgumentError` otherwise). The certified fidelity is **phase-optimized**
  (Newton refinement warm-started at the solver's phases).
- `KetTrajectory`, fixed phase: `|⟨ψ_goal|ψ⟩|²`.
- `MultiKetTrajectory`: coherent fidelity `|1/n Σᵢ ⟨ψᵢ_goal(θ)|ψᵢ⟩|²`
  (phase-optimized when free-phase globals are present).
- `UnitaryTrajectory`: Pedersen subspace fidelity for `EmbeddedOperator` goals
  (with the solver's phases when free-phase globals are present), standard
  unitary fidelity otherwise.

# Arguments
- `qtraj`: quantum trajectory (provides system, goal, and pulse type)
- `traj`: the solved `NamedTrajectory` (e.g. `get_trajectory(qcp)` post-`solve!`)

# Keyword Arguments
- `tol::Float64 = 1e-3`: certification gap tolerance
- `subsystem_levels`: per-subsystem level structure, required for free-phase kets
- `algorithm = nothing`: ODE algorithm for the fine rollout (default: the
  system's `default_algorithm`)
- `abstol = 1e-10`, `reltol = 1e-10`: fine-rollout tolerances
- `n_save::Int = 101`: number of saved rollout points
- `optimize_phases::Bool = true`: refine free-phase kets' phases on the fine state

# Returns
NamedTuple with fields:
- `fidelity`: the verified (fine-rollout, phase-optimized) fidelity
- `collocation_fidelity`: fidelity of the stored terminal decision variable
- `gap`: `|fidelity − collocation_fidelity|`
- `certified::Bool`: `gap ≤ tol`
- `duration`: physical duration of the extracted pulse (`∑Δt`)
- `pulse`: the extracted pulse that was rolled out
"""
function certify_rollout(
    qtraj::AbstractQuantumTrajectory,
    traj::NamedTrajectory;
    tol::Float64 = 1e-3,
    subsystem_levels::Union{Nothing,Vector{Int}} = nothing,
    algorithm = nothing,
    abstol::Real = 1e-10,
    reltol::Real = 1e-10,
    n_save::Int = 101,
    optimize_phases::Bool = true,
)
    F_colloc = _collocation_fidelity(qtraj, traj; subsystem_levels = subsystem_levels)
    pulse = extract_pulse(qtraj, traj)
    qtraj_fine = Rollouts.rollout(
        qtraj,
        pulse;
        algorithm = algorithm,
        n_save = n_save,
        abstol = abstol,
        reltol = reltol,
    )
    F_cert = _certified_fidelity(
        qtraj_fine,
        traj;
        subsystem_levels = subsystem_levels,
        optimize_phases = optimize_phases,
    )
    gap = abs(F_cert - F_colloc)
    return (
        fidelity = F_cert,
        collocation_fidelity = F_colloc,
        gap = gap,
        certified = gap <= tol,
        duration = duration(pulse),
        pulse = pulse,
    )
end

# ============================================================================ #
# Tests
# ============================================================================ #

@testitem "certify_rollout agrees on an honest fixed-phase ket trajectory" begin
    using LinearAlgebra
    using NamedTrajectories

    system = QuantumSystem(0.0 * PAULIS.Z, [PAULIS.X], [1.0])

    # Constant pulse u = π/(2T): |0⟩ → -i|1⟩ exactly, so both the analytic
    # terminal state and the fine rollout are known.
    T = 2.0
    N = 11
    times = collect(range(0.0, T, length = N))
    u = fill(π / (2T), 1, N)
    pulse = LinearSplinePulse(u, times)

    ψ0 = ComplexF64[1.0, 0.0]
    ψg = ComplexF64[0.0, -im]
    qtraj = KetTrajectory(system, pulse, ψ0, ψg)

    traj = NamedTrajectory(qtraj, times)
    # The NamedTrajectory conversion samples states from the ODE solution, so
    # the stored terminal state is honest — certification must pass with a
    # tiny gap.
    result = certify_rollout(qtraj, traj)
    @test result.certified
    @test result.gap <= 1e-3
    @test result.fidelity >= 0.999
    @test result.collocation_fidelity >= 0.999
    @test result.duration ≈ T
end

@testitem "certify_rollout catches a gamed terminal state" begin
    using LinearAlgebra
    using NamedTrajectories

    system = QuantumSystem(0.0 * PAULIS.Z, [PAULIS.X], [1.0])

    # Zero pulse: the true evolution stays at |0⟩ (rollout fidelity ≈ 0), but
    # we overwrite the stored terminal decision variable with the goal — the
    # discretization-gaming signature. certify_rollout must flag it.
    T = 2.0
    N = 11
    times = collect(range(0.0, T, length = N))
    pulse = LinearSplinePulse(zeros(1, N), times)

    ψ0 = ComplexF64[1.0, 0.0]
    ψg = ComplexF64[0.0, 1.0]
    qtraj = KetTrajectory(system, pulse, ψ0, ψg)

    traj = NamedTrajectory(qtraj, times)
    traj[traj.N][:ψ̃] .= ket_to_iso(ψg)

    result = certify_rollout(qtraj, traj)
    @test !result.certified
    @test result.gap > 0.9
    @test result.collocation_fidelity > 0.99  # the collocation *claims* success
    @test result.fidelity < 0.1               # the physics says otherwise
end

@testitem "certify_rollout free-phase ket requires subsystem_levels" begin
    using LinearAlgebra
    using NamedTrajectories

    system = QuantumSystem(0.1 * PAULIS.Z, [PAULIS.X], [1.0])

    T = 2.0
    N = 11
    times = collect(range(0.0, T, length = N))
    pulse = LinearSplinePulse(0.1 * randn(1, N), times)

    ψ0 = ComplexF64[1.0, 0.0]
    ψg = ComplexF64[1 / sqrt(2), 1 / sqrt(2)]
    qtraj = KetTrajectory(system, pulse, ψ0, ψg)

    # Trajectory carrying a free-phase global φ_1
    traj = NamedTrajectory(qtraj, times; global_data = Dict(:φ_1 => [0.3]))

    @test_throws ArgumentError certify_rollout(qtraj, traj)

    # With subsystem_levels and phases frozen at the stored value, both sides
    # evaluate the same fidelity of the same honest (ODE-sampled) state — the
    # gap must vanish. (With phase optimization on, the certified side may
    # legitimately exceed the collocation side here because the stored φ_1 was
    # never optimized; after a real solve the phases are optimal and the two
    # coincide.)
    result_frozen =
        certify_rollout(qtraj, traj; subsystem_levels = [2], optimize_phases = false)
    @test result_frozen.certified
    @test result_frozen.gap <= 1e-6
    # Phase optimization can only improve on the stored phase value.
    result = certify_rollout(qtraj, traj; subsystem_levels = [2])
    @test result.fidelity >= result_frozen.fidelity - 1e-12
end
