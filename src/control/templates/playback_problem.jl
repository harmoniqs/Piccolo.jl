export PlaybackProblem, PlaybackParams, PlaybackTieConstraint, PlaybackTiedIntegrator
export playback_tie_residual

using DirectTrajOpt: AbstractLinearConstraint, AbstractIntegrator
import DirectTrajOpt.CommonInterface: evaluate!, eval_jacobian, eval_hessian_of_lagrangian
import DirectTrajOpt.Integrators:
    get_jacobian_structure, get_hessian_of_lagrangian_structure
import ..QuantumControlProblems: interpolation_basis
import ..Verification: played_fidelity, design_pulse
using ...Quantum: PlaybackSpec, play, playback_jacobian, design_parameters
using ...Quantum: with_design_parameters, extract_pulse, n_drives, get_knot_times
using ...Quantum.Rollouts: rollout

# ============================================================================= #
# Optimize through playback (#387)
# ============================================================================= #
#
# Formulation (the design choice of record). Away from clipping and quantization,
# playback is linear in the design parameters: the played samples are `y = P·θ`,
# `P = playback_jacobian(spec, design)`. The NLP is built on the PLAYED lattice:
#
#   - the trajectory has one knot per played sample (plus the final knot), a fixed
#     timestep `Δt = dt_out`, and the played samples as its control component `u`;
#   - the dynamics are the exact exponential integrator, which is exact for the held
#     (ZOH) played waveform, so collocation == rollout of the played pulse;
#   - the design parameters `θ` are ONE global block (`:<u>_design`), the problem's
#     only extra decision variables (their count is the design's, independent of
#     the lattice);
#   - a linear equality `PlaybackTieConstraint` ties them: `u_k = (P θ)_k` at every
#     knot. Linear equalities are satisfied to rounding by every interior-point
#     step that is taken in full, so the tie is exact in practice, not a penalty.
#
# A tie constraint was chosen over an integrator that maps θ → u internally because it
# keeps the exponential integrator untouched (its analytic Daleckii–Krein derivatives,
# its sparsity) and keeps `u` a real trajectory component, so `extract_pulse`,
# `rollout`, `verify` and the displays all see the played waveform for free. The only
# adapter needed is `PlaybackTiedIntegrator`, which shows the exponential integrator a
# view of the trajectory without the θ block (the integrator treats every trajectory
# global as a dynamics global) and pads its derivatives with zero θ columns.

@doc raw"""
    PlaybackParams <: AbstractTemplateParams

The typed keyword surface of [`PlaybackProblem`](@ref): the infidelity weight `Q`
and the weight `R_u` of the quadratic regularizer on the played samples.
"""
Base.@kwdef struct PlaybackParams <: AbstractTemplateParams
    Q::Float64 = 100.0
    R_u::Float64 = 1e-2
end

@doc raw"""
    PlaybackProblem{P<:AbstractQuantumControlProblem, QT} <: AbstractProblemWrapper

Optimize the waveform the hardware **plays** rather than the designed one. Built from
a design problem (the `inner` problem: a [`SmoothPulseProblem`](@ref) or
[`SplinePulseProblem`](@ref)) and a [`PlaybackSpec`](@ref):

```math
\min_{\theta}\; Q\,\bigl(1 - F\bigl(U_{N}(y)\bigr)\bigr) + R_u \lVert y \rVert^2
\quad\text{s.t.}\quad y = P\,\theta,\qquad
U_{k+1} = e^{-i H(y_k)\,\Delta t}\, U_k ,
```

where `θ` are the design parameters ([`design_parameters`](@ref)), `P` is the
playback Jacobian ([`playback_jacobian`](@ref)) and `y` the played samples on the
lattice, each held for `Δt = dt_out`. The dynamics are exact for the held played
waveform; `θ` is a single global block of the trajectory and the tie `y = Pθ` is a
linear equality ([`PlaybackTieConstraint`](@ref)).

Usually built through the templates' `playback` keyword:

```julia
qcp = SplinePulseProblem(qtraj; playback = PlaybackSpec(dt_sample = 4.0, filter = RCLowpass(20.0)))
solve!(qcp; max_iter = 200)
played_fidelity(qcp)         # F_design, F_played of the optimized design
design_pulse(qcp)            # the optimized design pulse
inner(qcp)                   # the same problem on the design waveform
```

# Scope
- closed-system `UnitaryTrajectory` / `KetTrajectory` problems without system
  global parameters;
- `:hold` or FIR DAC interpolation (the played waveform is held on the output
  lattice), any causal filter and per-channel gain; `t0 = 0` and `delay = 0`
  (the played waveform starts at `t = 0`);
- no quantization (`bits = nothing`): rounding is not differentiable. Clipping is not
  modelled inside the NLP — keep the played samples inside `full_scale` (the
  system's drive bounds bound the played samples).
"""
mutable struct PlaybackProblem{
    P<:AbstractQuantumControlProblem,
    QT<:AbstractQuantumTrajectory,
} <: AbstractProblemWrapper
    inner::P
    qtraj::QT
    prob::DirectTrajOptProblem
    params::PlaybackParams
    spec::Union{Nothing,AbstractProblemSpec}
    playback::PlaybackSpec
    design::AbstractPulse
    tie::SparseMatrixCSC{Float64,Int}
end

wrapper_kind(::PlaybackProblem) = :playback

rewrap(w::PlaybackProblem, qtraj, prob::DirectTrajOptProblem) =
    PlaybackProblem(inner(w), qtraj, prob, w.params, w.spec, w.playback, w.design, w.tie)

# the template constructors stamp the tag + params onto the inner design problem
_retag(w::PlaybackProblem, tag::AbstractProblemTemplate, params::AbstractTemplateParams) =
    PlaybackProblem(
        _retag(inner(w), tag, params),
        w.qtraj,
        w.prob,
        w.params,
        w.spec,
        w.playback,
        w.design,
        w.tie,
    )

"""
    PlaybackProblem(qcp, spec::PlaybackSpec; Q = 100.0, R_u = 1e-2, piccolo_options = PiccoloOptions())

Build the optimize-through-playback problem for the design problem `qcp`. The design
pulse (and so the design parameters' initial values and the design knots) is
`qcp`'s current pulse; its boundary conditions on the design knots are kept. The
NLP is new: see [`PlaybackProblem`](@ref) for the formulation.
"""
function PlaybackProblem(
    qcp::AbstractQuantumControlProblem,
    spec::PlaybackSpec;
    Q::Real = 100.0,
    R_u::Real = 1e-2,
    piccolo_options::PiccoloOptions = PiccoloOptions(),
)
    qtraj = quantum_trajectory(qcp)
    _check_playback_scope(qtraj, spec)
    design_traj = get_trajectory(qcp)
    design = extract_pulse(qtraj, design_traj)
    u_name = drive_name(qtraj)
    θ_name = _design_global(u_name)

    _, info = play(spec, design)
    any(info.clipped) && @warn "PlaybackProblem: the initial design clips at full scale \
        ($(count(info.clipped)) samples); clipping is not modelled inside the NLP"
    J, _ = playback_jacobian(spec, design)
    n = n_drives(design)
    # the final knot holds the last played sample (it carries no dynamics)
    tie = [J; J[(end-n+1):end, :]]
    N_out = size(info.samples, 2)
    times = collect(info.dt_out .* (0:N_out))
    played = ZeroOrderPulse(
        hcat(info.samples, info.samples[:, end]),
        times;
        drive_name = u_name,
        initial_value = :free,
        final_value = :free,
    )
    played_qtraj = rollout(qtraj, played)
    θ0 = design_parameters(design)
    traj = NamedTrajectory(
        played_qtraj,
        times;
        Δt_bounds = (info.dt_out, info.dt_out),
        global_data = Dict(θ_name => θ0),
    )

    state_sym = state_name(qtraj)
    constraints = AbstractConstraint[
        PlaybackTieConstraint(u_name, θ_name, tie),
        _design_boundary_constraints(design, design_traj, u_name, θ_name)...,
    ]
    obj = _state_objective(played_qtraj, traj, state_sym, Float64(Q))
    obj += QuadraticRegularizer(u_name, traj, Float64(R_u))
    obj +=
        _apply_piccolo_options(played_qtraj, piccolo_options, constraints, traj, state_sym)
    integrator =
        PlaybackTiedIntegrator(HermitianExponentialIntegrator(played_qtraj, traj.N), θ_name)
    prob = DirectTrajOptProblem(traj, obj, AbstractIntegrator[integrator]; constraints)
    return PlaybackProblem(
        qcp,
        played_qtraj,
        prob,
        PlaybackParams(Float64(Q), Float64(R_u)),
        nothing,
        spec,
        design,
        tie,
    )
end

_design_global(u_name::Symbol) = Symbol(u_name, :_design)

function _check_playback_scope(qtraj, spec::PlaybackSpec)
    qtraj isa Union{UnitaryTrajectory,KetTrajectory} || throw(
        ArgumentError(
            "playback: optimize-through-playback supports UnitaryTrajectory and " *
            "KetTrajectory problems, got $(_typename(qtraj))",
        ),
    )
    isempty(get_system(qtraj).global_params) ||
        throw(ArgumentError("playback: systems with global parameters are not supported"))
    spec.dac_interp === :linear && throw(
        ArgumentError(
            "playback: a :linear DAC plays a ramp, not a held waveform; optimize-through-" *
            "playback needs :hold or FIR interpolation (verify a :linear DAC with played_fidelity)",
        ),
    )
    isnothing(spec.bits) || throw(
        ArgumentError(
            "playback: quantization (bits = $(spec.bits)) is not differentiable; optimize " *
            "with a spec without `bits` and verify the quantized spec with played_fidelity",
        ),
    )
    spec.t0 + spec.delay == 0 || throw(
        ArgumentError(
            "playback: the played lattice must start at t = 0 (t0 = $(spec.t0), " *
            "delay = $(spec.delay))",
        ),
    )
    return nothing
end

# The design problem's boundary conditions on its control knots (initial / final
# values of `u`, and of `du` for a cubic design), restated on the θ block.
function _design_boundary_constraints(design, design_traj, u_name, θ_name)
    M = length(get_knot_times(design))
    n = n_drives(design)
    stride = length(design_parameters(design)) ÷ M
    rows = Int[]
    vals = Float64[]
    for (fixed, k) in ((design_traj.initial, 1), (design_traj.final, M)),
        (name, offset) in ((u_name, 0), (Symbol(:d, u_name), n))

        (offset + n ≤ stride && haskey(fixed, name)) || continue
        append!(rows, (k - 1) * stride + offset .+ (1:n))
        append!(vals, fixed[name])
    end
    isempty(rows) && return AbstractConstraint[]
    A = sparse(1:length(rows), rows, 1.0, length(rows), stride * M)
    return AbstractConstraint[GlobalLinearConstraint(
        θ_name,
        A,
        vals;
        label = "design boundary conditions on :$θ_name",
    )]
end

# ----------------------------------------------------------------------------- #
# The linear tie u_k = (P θ)_k
# ----------------------------------------------------------------------------- #

"""
    PlaybackTieConstraint(control, design, P; label)

Linear equality tying the played-sample control component `control` to the design
global block `design`: `vec(u) = P · θ` (rows ordered like `vec(u)`, drive fastest,
then knot). Emitted as one affine `MOI.EqualTo(0)` row per played sample.
"""
struct PlaybackTieConstraint <: AbstractLinearConstraint
    control::Symbol
    design::Symbol
    P::SparseMatrixCSC{Float64,Int}
    label::String
end

PlaybackTieConstraint(
    control::Symbol,
    design::Symbol,
    P::AbstractMatrix;
    label::String = "playback tie :$control = P·:$design",
) = PlaybackTieConstraint(control, design, sparse(P), label)

function (con::PlaybackTieConstraint)(opt, vars::AbstractVector, traj::NamedTrajectory)
    MOI = DirectTrajOpt.Solvers.MOI
    traj.warp === nothing ||
        throw(ArgumentError("PlaybackTieConstraint does not support time warps"))
    u = traj.components[con.control]
    g = traj.global_components[con.design]
    n = length(u)
    size(con.P) == (n * traj.N, length(g)) || throw(
        DimensionMismatch(
            "PlaybackTieConstraint: P is $(size(con.P)), expected $((n * traj.N, length(g)))",
        ),
    )
    base = traj.dim * traj.N
    Pt = sparse(transpose(con.P))      # column r of Pt is row r of P
    for k = 1:traj.N, c = 1:n
        r = (k - 1) * n + c
        terms = [MOI.ScalarAffineTerm(1.0, vars[(k-1)*traj.dim+u[c]])]
        for j in nzrange(Pt, r)
            push!(
                terms,
                MOI.ScalarAffineTerm(-nonzeros(Pt)[j], vars[base+g[rowvals(Pt)[j]]]),
            )
        end
        MOI.add_constraint(opt, MOI.ScalarAffineFunction(terms, 0.0), MOI.EqualTo(0.0))
    end
    return nothing
end

Base.show(io::IO, c::PlaybackTieConstraint) = print(
    io,
    "PlaybackTieConstraint: \"$(c.label)\" ($(size(c.P, 1)) rows, $(size(c.P, 2)) design parameters)",
)

"""
    playback_tie_residual(qcp::PlaybackProblem) -> Float64

`maximum |vec(u) - P·θ|` on the problem's current trajectory: how exactly the played
samples follow the design parameters.
"""
function playback_tie_residual(qcp::PlaybackProblem)
    traj = get_trajectory(qcp)
    u = vec(traj[drive_name(quantum_trajectory(qcp))])
    θ = traj.global_data[traj.global_components[_design_global(drive_name(qcp.qtraj))]]
    return maximum(abs, u .- qcp.tie * θ)
end

# ----------------------------------------------------------------------------- #
# Design pulse, played fidelity
# ----------------------------------------------------------------------------- #

"""
    design_pulse(qcp::PlaybackProblem) -> AbstractPulse

The design pulse at the problem's current design parameters.
"""
function design_pulse(qcp::PlaybackProblem)
    traj = get_trajectory(qcp)
    θ = traj.global_data[traj.global_components[_design_global(drive_name(qcp.qtraj))]]
    return with_design_parameters(qcp.design, θ)
end

"""
    played_fidelity(qcp::PlaybackProblem; kwargs...)

[`played_fidelity`](@ref) of the optimized design pulse under the problem's own
playback spec.
"""
played_fidelity(qcp::PlaybackProblem; kwargs...) =
    played_fidelity(qcp, qcp.playback; kwargs...)

# ----------------------------------------------------------------------------- #
# The integrator adapter
# ----------------------------------------------------------------------------- #

"""
    PlaybackTiedIntegrator(inner, design)

Run the dynamics integrator `inner` on a trajectory that also carries the design
global block `design`, which the dynamics do not depend on. The inner integrator sees
the trajectory without the block (sharing its knot data) and its Jacobian and Hessian
are padded with zero columns for it. The block must be the trajectory's only global.
"""
struct PlaybackTiedIntegrator{I<:AbstractIntegrator} <: AbstractIntegrator
    inner::I
    design::Symbol
    x_names::Vector{Symbol}
    x_dim::Int
    dim::Int
end

PlaybackTiedIntegrator(inner::AbstractIntegrator, design::Symbol) =
    PlaybackTiedIntegrator(inner, design, inner.x_names, inner.x_dim, inner.dim)

interpolation_basis(ℰ::PlaybackTiedIntegrator) = interpolation_basis(ℰ.inner)

function _without_design(ℰ::PlaybackTiedIntegrator, traj::NamedTrajectory)
    keys(traj.global_components) == (ℰ.design,) || throw(
        ArgumentError(
            "PlaybackTiedIntegrator: :$(ℰ.design) must be the trajectory's only global " *
            "block, found $(keys(traj.global_components))",
        ),
    )
    return NamedTrajectory(traj; global_data = Float64[], global_components = NamedTuple())
end

_pad_columns(A, g) = [A spzeros(size(A, 1), g)]
_pad_square(A, g) = blockdiag(A, spzeros(g, g))

evaluate!(δ::AbstractVector, ℰ::PlaybackTiedIntegrator, traj::NamedTrajectory) =
    evaluate!(δ, ℰ.inner, _without_design(ℰ, traj))

eval_jacobian(ℰ::PlaybackTiedIntegrator, traj::NamedTrajectory) =
    _pad_columns(eval_jacobian(ℰ.inner, _without_design(ℰ, traj)), traj.global_dim)

get_jacobian_structure(ℰ::PlaybackTiedIntegrator, traj::NamedTrajectory) =
    _pad_columns(get_jacobian_structure(ℰ.inner, _without_design(ℰ, traj)), traj.global_dim)

eval_hessian_of_lagrangian(
    ℰ::PlaybackTiedIntegrator,
    traj::NamedTrajectory,
    μ::AbstractVector,
) = _pad_square(
    eval_hessian_of_lagrangian(ℰ.inner, _without_design(ℰ, traj), μ),
    traj.global_dim,
)

get_hessian_of_lagrangian_structure(ℰ::PlaybackTiedIntegrator, traj::NamedTrajectory) =
    _pad_square(
        get_hessian_of_lagrangian_structure(ℰ.inner, _without_design(ℰ, traj)),
        traj.global_dim,
    )

# ----------------------------------------------------------------------------- #
# The templates' `playback` keyword
# ----------------------------------------------------------------------------- #

# Keywords whose meaning is tied to the design lattice or to globals; refused rather
# than silently dropped when `playback` is set.
const _PLAYBACK_REFUSED = (
    :integrator,
    :global_names,
    :global_bounds,
    :calibration_targets,
    :free_phase,
    :Δt_bounds,
)

_with_playback(qcp, ::Nothing, kwargs) = qcp

function _with_playback(qcp, spec::PlaybackSpec, kwargs)
    for k in _PLAYBACK_REFUSED
        v = get(kwargs, k, nothing)
        (isnothing(v) || v === false || (v isa AbstractVector && isempty(v))) && continue
        throw(ArgumentError("playback: keyword `$k` is not supported with `playback`"))
    end
    R_u = get(kwargs, :R_u, get(kwargs, :R, 1e-2))
    R_u isa Real || throw(
        ArgumentError("playback: `R_u` must be a scalar (it weights the played samples)"),
    )
    return PlaybackProblem(
        qcp,
        spec;
        Q = get(kwargs, :Q, 100.0),
        R_u = R_u,
        piccolo_options = get(kwargs, :piccolo_options, PiccoloOptions()),
    )
end

# ============================================================================= #
# Tests
# ============================================================================= #

@testitem "PlaybackProblem: optimizing through playback beats design-then-play (4 ns hold + RC)" begin
    using LinearAlgebra
    using NamedTrajectories

    # Detuned qubit (Δ/2π = 20 MHz, ns units), X gate in 40 ns, cubic design with 11
    # knots; the hardware holds 4 ns samples and low-passes them (τ = 6 ns RC).
    sys = QuantumSystem(
        2π * 0.02 / 2 * PAULIS.Z,
        [PAULIS.X / 2, PAULIS.Y / 2],
        [(-0.5, 0.5), (-0.5, 0.5)],
    )
    T, M = 40.0, 11
    times = collect(range(0, T, length = M))
    u0 = [0.08 * sin(π * t / T)^2 * (c == 1 ? 1 : 0.3) for c = 1:2, t in times]
    design_qtraj() =
        UnitaryTrajectory(sys, CubicSplinePulse(u0, zeros(2, M), times), GATES[:X])
    spec = PlaybackSpec(dt_sample = 4.0, filter = RCLowpass(6.0))
    opts = PiccoloOptions(display = :silent)

    # design-then-play: optimize the design waveform (duration fixed), then play it
    dt = T / (M - 1)
    base = SplinePulseProblem(design_qtraj(); Δt_bounds = (dt, dt), piccolo_options = opts)
    solve!(base; max_iter = 150, verbose = false, print_level = 0)
    r_base = played_fidelity(base, spec)
    @test r_base.F_design > 0.999

    # optimize through playback
    qcp = SplinePulseProblem(design_qtraj(); playback = spec, piccolo_options = opts)
    @test qcp isa PlaybackProblem
    @test inner(qcp) isa SplinePulseProblem
    θ0 = copy(get_trajectory(qcp).global_data)
    solve!(qcp; max_iter = 300, verbose = false, print_level = 0)
    r = played_fidelity(qcp)
    @test r.F_played > 0.9999
    @test r.F_played > r_base.F_played + 0.5
    @info "playback benchmark (played fidelity)" design_then_play = r_base.F_played through_playback =
        r.F_played

    traj = get_trajectory(qcp)
    θ = traj.global_data
    @test θ != θ0                                   # the design parameters moved

    # AC3: the decision variables are the design knots, the tie is exact
    @test length(θ) == length(design_parameters(design_pulse(qcp))) == 4M
    @test playback_tie_residual(qcp) ≤ 1e-12
    _, info = play(spec, design_pulse(qcp))
    @test maximum(abs, traj[:u][:, 1:(end-1)] .- info.samples) ≤ 1e-12

    # AC2: collocation == exact propagation of the held played waveform
    steps = [
        exp(-im * Matrix(sys.H(info.samples[:, k], 0.0)) * info.dt_out) for k = 1:(traj.N-1)
    ]
    exact = accumulate((U, E) -> E * U, steps; init = Matrix{ComplexF64}(I, 2, 2))
    @test maximum(
        norm(iso_vec_to_operator(traj[:Ũ⃗][:, k+1]) - exact[k]) for k = 1:(traj.N-1)
    ) ≤ 1e-8
    v = verify(qcp; verbose = false)
    @test v.ΔF ≤ 1e-6
end

@testitem "PlaybackProblem: construction, lattice independence and scope" begin
    using LinearAlgebra
    using NamedTrajectories

    sys = QuantumSystem(0.1 * PAULIS.Z, [PAULIS.X, PAULIS.Y], [(-1.0, 1.0), (-1.0, 1.0)])
    T, M = 20.0, 6
    times = collect(range(0, T, length = M))
    u0 = 0.1 .* [sin(π * t / T) for _ = 1:2, t in times]
    opts = PiccoloOptions(display = :silent)
    cubic() = UnitaryTrajectory(sys, CubicSplinePulse(u0, zeros(2, M), times), GATES[:X])

    # the design-parameter count does not depend on the lattice
    probs = map((4.0, 1.0)) do dt_s
        SplinePulseProblem(
            cubic();
            playback = PlaybackSpec(dt_sample = dt_s, filter = RCLowpass(3.0)),
            piccolo_options = opts,
        )
    end
    @test [p.prob.trajectory.N for p in probs] == [6, 21]
    @test [length(p.prob.trajectory.global_data) for p in probs] == [4M, 4M]
    for p in probs
        traj = get_trajectory(p)
        @test traj.global_names == (:u_design,)
        @test all(traj.bounds.Δt[1] .== traj.bounds.Δt[2])       # Δt fixed to the lattice
        @test playback_tie_residual(p) ≤ 1e-15
        @test template_params(p) isa PlaybackParams
        @test nameof(typeof(template_params(inner(p)))) === :SplinePulseParams
        @test wrapper_kind(p) === :playback
        @test p.qtraj.pulse isa ZeroOrderPulse
        # design boundary conditions are carried to the design knots
        @test any(c -> c isa GlobalLinearConstraint, p.prob.constraints)
    end

    # linear design on SplinePulseProblem, held design on SmoothPulseProblem, ket goals
    spec = PlaybackSpec(dt_sample = 2.0, filter = ButterworthLowpass(2, 0.1))
    lin = SplinePulseProblem(
        KetTrajectory(
            sys,
            LinearSplinePulse(u0, times),
            ComplexF64[1, 0],
            ComplexF64[0, 1],
        );
        playback = spec,
        piccolo_options = opts,
    )
    @test lin isa PlaybackProblem
    @test length(get_trajectory(lin).global_data) == 2M
    zoh = SmoothPulseProblem(
        UnitaryTrajectory(sys, ZeroOrderPulse(u0, times), GATES[:X]),
        M;
        playback = spec,
        piccolo_options = opts,
    )
    @test zoh isa PlaybackProblem
    @test inner(zoh) isa SmoothPulseProblem
    @test playback_tie_residual(zoh) ≤ 1e-15
    @test design_pulse(zoh) isa ZeroOrderPulse

    # out of scope: refused, never silently approximated
    for bad in (
        PlaybackSpec(dt_sample = 2.0, dac_interp = :linear),
        PlaybackSpec(dt_sample = 2.0, bits = 12, full_scale = 1.0),
        PlaybackSpec(dt_sample = 2.0, delay = 3.0),
    )
        @test_throws ArgumentError SplinePulseProblem(
            cubic();
            playback = bad,
            piccolo_options = opts,
        )
    end
    @test_throws ArgumentError SplinePulseProblem(
        cubic();
        playback = spec,
        Δt_bounds = (0.1, 5.0),
        piccolo_options = opts,
    )
end

@testitem "templates without `playback` are unchanged" begin
    using NamedTrajectories
    using Piccolo.Control.ProblemTemplates: _spline_pulse_problem, _smooth_pulse_problem

    sys = QuantumSystem(0.1 * PAULIS.Z, [PAULIS.X, PAULIS.Y], [(-1.0, 1.0), (-1.0, 1.0)])
    T, M = 10.0, 6
    times = collect(range(0, T, length = M))
    u0 = 0.1 .* [sin(π * t / T) for _ = 1:2, t in times]
    opts = PiccoloOptions(display = :silent)
    cubic() = UnitaryTrajectory(sys, CubicSplinePulse(u0, zeros(2, M), times), GATES[:X])
    zoh() = UnitaryTrajectory(sys, ZeroOrderPulse(u0, times), GATES[:X])

    pairs = [
        (
            SplinePulseProblem(cubic(); piccolo_options = opts),
            _spline_pulse_problem(cubic(); piccolo_options = opts),
        ),
        (
            SplinePulseProblem(cubic(); playback = nothing, piccolo_options = opts),
            _spline_pulse_problem(cubic(); piccolo_options = opts),
        ),
        (
            SmoothPulseProblem(zoh(), M; piccolo_options = opts),
            _smooth_pulse_problem(zoh(), M; piccolo_options = opts),
        ),
    ]
    for (qcp, ref) in pairs
        @test qcp isa QuantumControlProblem
        @test !(qcp isa AbstractProblemWrapper)
        @test !hasproperty(template_params(qcp), :playback)
        @test get_trajectory(qcp).datavec == get_trajectory(ref).datavec
        @test typeof.(qcp.prob.integrators) == typeof.(ref.prob.integrators)
        @test typeof.(qcp.prob.constraints) == typeof.(ref.prob.constraints)
        solve!(qcp; max_iter = 5, verbose = false, print_level = 0)
        solve!(ref; max_iter = 5, verbose = false, print_level = 0)
        @test get_trajectory(qcp).datavec == get_trajectory(ref).datavec   # bit-for-bit
    end
end
