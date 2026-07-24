export BSplineDilationBoundConstraint

using SparseArrays
using TrajectoryIndexingUtils
using NamedTrajectories: get_timesteps
import DirectTrajOpt.CommonInterface

@doc raw"""
    BSplineDilationBoundConstraint <: AbstractNonlinearConstraint

Free-time (minimum-time-valid) derivative bound on a B-spline control-point
global block: with `s = (∑ₖ Δtₖ) / T_ref` the uniform time dilation of the
trajectory relative to the construction-time duration `T_ref` that the knot
vector spans, each row `r` of the derivative-difference operator `A` (the
scaled first/second control-point differences from
[`bspline_slew_constraint`](@ref) / [`bspline_accel_constraint`](@ref))
is bounded by the dilation-scaled cap

    |(A c)_r| ≤ b_r · s^p,

emitted as the two one-sided inequality rows (interleaved per base row)

    g⁺_r =  (A c)_r − b_r · s^p ≤ 0,
    g⁻_r = −(A c)_r − b_r · s^p ≤ 0.

Physically: under uniform dilation the pulse's physical derivatives scale as
`1/s` per derivative order, so the knot-parameter derivative rows must be
allowed to grow as `s^p` for the physical bound (`v_max` for `p = 1` slew,
`a_max` for `p = 2` curvature) to stay exactly enforced as the pulse
compresses or stretches. The timestep sum runs over `k = 1 … N-1`, matching
`NamedTrajectories.get_duration` and `MinimumTimeObjective` (the last knot's
`Δt` is a padded dummy).

All derivatives are exact and analytic (no ForwardDiff):

- Jacobian: `±A` on the `:c_<drive>` global columns;
  `−b_r · p · s^{p−1} / T_ref` on every `Δtₖ` column (`k ≤ N-1`). For `p = 2`
  this is `−2 (b_r / T_ref²) ∑ₖ Δtₖ`.
- Hessian of the Lagrangian: zero for `p = 1`;
  for `p = 2` the constant `−2 b_r / T_ref²` on every `Δt × Δt` pair,
  contracted with the multipliers.

!!! note
    The Jacobian/Hessian sparsity is captured once by the `Evaluator` at the
    initial trajectory; the `Δt` Jacobian entries for `p = 2` vanish only when
    `∑ₖ Δtₖ = 0`, which no valid trajectory satisfies.

# Fields
- `A::SparseMatrixCSC{Float64,Int}`: `(n_rows × global_dim)` derivative rows
  over the `:c_<drive>` global block.
- `b::Vector{Float64}`: per-row bound at `s = 1` (length `n_rows`).
- `power::Int`: dilation order `p` (`1` = slew, `2` = acceleration).
- `T_ref::Float64`: reference duration spanned by the pulse's knot vector.
- `global_name::Symbol`: global block to constrain (e.g. `:c_u`).
- `equality::Bool`: always `false` (inequality `g ≤ 0`).
- `dim::Int`: total constraint dimension, `2 · n_rows`.
- `label::String`
"""
struct BSplineDilationBoundConstraint <: AbstractNonlinearConstraint
    A::SparseMatrixCSC{Float64,Int}
    b::Vector{Float64}
    power::Int
    T_ref::Float64
    global_name::Symbol
    equality::Bool
    dim::Int
    label::String
end

"""
    BSplineDilationBoundConstraint(global_name, A, b, power, T_ref; label=...)

Construct the free-time dilation bound `|A·g| ≤ b · ((∑ₖ Δtₖ)/T_ref)^power`
on the global block `global_name`. `A` may be any `AbstractMatrix` (stored
sparse); `b` is the per-row bound at unit dilation.
"""
function BSplineDilationBoundConstraint(
    global_name::Symbol,
    A::AbstractMatrix,
    b::AbstractVector,
    power::Int,
    T_ref::Real;
    label::String = "B-spline free-time dilation bound on :$global_name",
)
    power in (1, 2) || throw(ArgumentError(
        "power must be 1 (slew) or 2 (acceleration); got $power"))
    T_ref > 0 || throw(ArgumentError("T_ref must be positive; got $T_ref"))
    size(A, 1) == length(b) || throw(ArgumentError(
        "row count mismatch: A has $(size(A, 1)) rows, b has $(length(b))"))
    all(>=(0), b) || throw(ArgumentError("bounds b must be nonnegative"))
    As = convert(SparseMatrixCSC{Float64,Int}, A)
    return BSplineDilationBoundConstraint(
        As,
        Vector{Float64}(b),
        power,
        Float64(T_ref),
        global_name,
        false,                # inequality g ≤ 0
        2 * size(As, 1),
        label,
    )
end

function Base.show(io::IO, con::BSplineDilationBoundConstraint)
    kind = con.power == 1 ? "slew" : "acceleration"
    print(
        io,
        "BSplineDilationBoundConstraint: \"$(con.label)\" " *
        "($(size(con.A, 1)) $kind rows on :$(con.global_name), T_ref = $(con.T_ref))",
    )
end

# Uniform dilation s = (∑_{k=1}^{N-1} Δt_k) / T_ref — the sum convention
# matches get_duration / MinimumTimeObjective (last Δt is a padded dummy).
function _dilation(con::BSplineDilationBoundConstraint, traj::NamedTrajectory)
    @assert traj.timestep isa Symbol (
        "BSplineDilationBoundConstraint requires a variable timestep " *
        "(traj.timestep must be a Symbol)")
    Δt = get_timesteps(traj)
    return sum(view(Δt, 1:(traj.N - 1))) / con.T_ref
end

# Global-block view of the constrained control points, with shape validation.
function _dilation_cp_view(con::BSplineDilationBoundConstraint, traj::NamedTrajectory)
    haskey(traj.global_components, con.global_name) || throw(ArgumentError(
        "global block :$(con.global_name) not found in trajectory " *
        "(have: $(collect(keys(traj.global_components))))"))
    gcomps = traj.global_components[con.global_name]
    length(gcomps) == size(con.A, 2) || throw(ArgumentError(
        "global block :$(con.global_name) has length $(length(gcomps)); " *
        "constraint rows expect $(size(con.A, 2))"))
    return view(traj.global_data, gcomps), gcomps
end

# --------------------------------------------------------------------------- #
# CommonInterface implementation (exact analytic derivatives)
# --------------------------------------------------------------------------- #

function CommonInterface.evaluate!(
    values::AbstractVector,
    con::BSplineDilationBoundConstraint,
    traj::NamedTrajectory,
)
    c, _ = _dilation_cp_view(con, traj)
    sᵖ = _dilation(con, traj)^con.power
    Ac = con.A * c
    for r in 1:size(con.A, 1)
        cap = con.b[r] * sᵖ
        values[2r - 1] = Ac[r] - cap
        values[2r] = -Ac[r] - cap
    end
    return nothing
end

function CommonInterface.eval_jacobian(
    con::BSplineDilationBoundConstraint,
    traj::NamedTrajectory,
)
    Z_dim = traj.dim * traj.N + traj.global_dim
    _, gcomps = _dilation_cp_view(con, traj)
    gcols = traj.dim * traj.N .+ gcomps

    n_rows = size(con.A, 1)
    Is, Js, Vs = Int[], Int[], Float64[]
    sizehint = 2 * nnz(con.A) + 2 * n_rows * (traj.N - 1)
    sizehint!(Is, sizehint); sizehint!(Js, sizehint); sizehint!(Vs, sizehint)

    # Control-point block: constant ±A on the global columns.
    Arows = rowvals(con.A)
    Avals = nonzeros(con.A)
    for j in 1:size(con.A, 2)
        for idx in nzrange(con.A, j)
            r = Arows[idx]
            v = Avals[idx]
            push!(Is, 2r - 1); push!(Js, gcols[j]); push!(Vs, v)
            push!(Is, 2r);     push!(Js, gcols[j]); push!(Vs, -v)
        end
    end

    # Timestep block: ∂(−b_r s^p)/∂Δt_k = −b_r · p · s^{p−1} / T_ref,
    # identical for both one-sided rows of a base row.
    s = _dilation(con, traj)
    dcap = con.power * s^(con.power - 1) / con.T_ref
    Δt_comps = traj.components[traj.timestep]
    for k in 1:(traj.N - 1)
        col = slice(k, Δt_comps, traj.dim)[1]
        for r in 1:n_rows
            coeff = -con.b[r] * dcap
            push!(Is, 2r - 1); push!(Js, col); push!(Vs, coeff)
            push!(Is, 2r);     push!(Js, col); push!(Vs, coeff)
        end
    end

    return sparse(Is, Js, Vs, con.dim, Z_dim)
end

function CommonInterface.eval_hessian_of_lagrangian(
    con::BSplineDilationBoundConstraint,
    traj::NamedTrajectory,
    μ::AbstractVector,
)
    Z_dim = traj.dim * traj.N + traj.global_dim

    # Slew (p = 1) is jointly affine in (c, Δt): zero Hessian.
    con.power == 1 && return spzeros(Z_dim, Z_dim)

    # Acceleration (p = 2): ∂²g_i/∂Δt_j∂Δt_k = −2 b_{⌈i/2⌉} / T_ref² for every
    # Δt pair — a constant block, contracted with the multipliers here.
    w = 0.0
    for r in 1:size(con.A, 1)
        w -= (2.0 / con.T_ref^2) * con.b[r] * (μ[2r - 1] + μ[2r])
    end

    Δt_comps = traj.components[traj.timestep]
    cols = [slice(k, Δt_comps, traj.dim)[1] for k in 1:(traj.N - 1)]
    n = length(cols)
    Is = repeat(cols, inner = n)
    Js = repeat(cols, outer = n)
    return sparse(Is, Js, fill(w, n * n), Z_dim, Z_dim)
end

# ============================================================================= #
# TestItems
# ============================================================================= #

@testitem "bspline free-time constraint jacobians/hessian match finite differences" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using DirectTrajOpt: CommonInterface
    using LinearAlgebra
    using SparseArrays
    using Random

    Random.seed!(42)

    # Small randomized Layout-A-style trajectory: dummy state :x, control
    # points in the :c_u global block, variable (nonuniform) timesteps so the
    # dilation s ≠ 1 and the FD check hits a generic point.
    M, order, n_d = 9, 4, 2
    T_ref = 3.0
    N = M - order + 2
    times = collect(range(0.0, T_ref, length = N))
    pulse = BSplinePulse(randn(n_d, M), times; order = order)

    Δt = reshape(0.3 .+ 0.4 * rand(N), 1, N)
    traj = NamedTrajectory(
        (x = randn(2, N), Δt = Δt),
        (c_u = collect(vec(pulse.control_points)),);
        timestep = :Δt,
        controls = (:Δt,),
    )

    slew = bspline_slew_constraint(pulse, [1.3, 0.9]; free_time = true, T_ref = T_ref)
    accel = bspline_accel_constraint(pulse, [2.1, 1.7]; free_time = true, T_ref = T_ref)

    n_traj = traj.dim * traj.N
    Z0 = vcat(collect(traj.datavec), collect(traj.global_data))

    eval_at = (con, Z) -> begin
        wrap = NamedTrajectory(
            traj;
            datavec = Z[1:n_traj],
            global_data = Z[(n_traj + 1):end],
        )
        vals = zeros(con.dim)
        CommonInterface.evaluate!(vals, con, wrap)
        return vals
    end

    # Central-difference Jacobian over the full (datavec, global_data) vector.
    # g is affine in c and polynomial (degree ≤ 2) in Δt, so central FD is
    # exact up to roundoff — the 1e-8 tolerance is comfortably met.
    h = 1e-5
    for con in (slew, accel)
        ∂g = Matrix(CommonInterface.eval_jacobian(con, traj))
        @test size(∂g) == (con.dim, length(Z0))
        ∂g_fd = zeros(con.dim, length(Z0))
        for j in eachindex(Z0)
            Zp = copy(Z0); Zp[j] += h
            Zm = copy(Z0); Zm[j] -= h
            ∂g_fd[:, j] = (eval_at(con, Zp) - eval_at(con, Zm)) / (2h)
        end
        @test maximum(abs.(∂g - ∂g_fd)) <= 1e-8
    end

    # Accel Hessian of the Lagrangian vs central second differences of μ'g.
    μ = rand(accel.dim)
    H = Matrix(CommonInterface.eval_hessian_of_lagrangian(accel, traj, μ))
    φ = Z -> μ' * eval_at(accel, Z)
    h2 = 1e-2
    H_fd = zeros(length(Z0), length(Z0))
    for j in eachindex(Z0), l in j:length(Z0)
        Zpp = copy(Z0); Zpp[j] += h2; Zpp[l] += h2
        Zpm = copy(Z0); Zpm[j] += h2; Zpm[l] -= h2
        Zmp = copy(Z0); Zmp[j] -= h2; Zmp[l] += h2
        Zmm = copy(Z0); Zmm[j] -= h2; Zmm[l] -= h2
        H_fd[j, l] = (φ(Zpp) - φ(Zpm) - φ(Zmp) + φ(Zmm)) / (4h2^2)
        H_fd[l, j] = H_fd[j, l]
    end
    @test maximum(abs.(H - H_fd)) <= 1e-8

    # Slew (p = 1) is jointly affine: exactly zero Hessian.
    Hs = CommonInterface.eval_hessian_of_lagrangian(slew, traj, rand(slew.dim))
    @test nnz(Hs) == 0
end

@testitem "bspline free-time constraints equal fixed-time at s = 1" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using DirectTrajOpt: CommonInterface
    using Random

    Random.seed!(0)

    M, order, n_d = 8, 4, 1
    T_ref = 2.0
    N = M - order + 2
    times = collect(range(0.0, T_ref, length = N))
    pulse = BSplinePulse(randn(n_d, M), times; order = order)

    # Nominal timesteps: ∑_{k<N} Δt_k = T_ref exactly, i.e. s = 1.
    Δt = fill(T_ref / (N - 1), 1, N)
    traj = NamedTrajectory(
        (x = randn(2, N), Δt = Δt),
        (c_u = collect(vec(pulse.control_points)),);
        timestep = :Δt,
        controls = (:Δt,),
    )
    c = traj.global_data[traj.global_components[:c_u]]

    for (make, bound) in ((bspline_slew_constraint, 0.7), (bspline_accel_constraint, 1.1))
        fixed = make(pulse, bound)   # GlobalLinearConstraint: −b ≤ A·c ≤ b
        free = make(pulse, bound; free_time = true, T_ref = T_ref)
        @test free isa BSplineDilationBoundConstraint
        vals = zeros(free.dim)
        CommonInterface.evaluate!(vals, free, traj)
        Ac = fixed.A * c
        # Interleaved rows (g⁺, g⁻) must reproduce the fixed-time margins.
        @test vals[1:2:end] ≈ Ac .- fixed.ub atol = 1e-12
        @test vals[2:2:end] ≈ -Ac .- fixed.ub atol = 1e-12
    end
end

@testitem "bspline free-time constraints tighten under uniform timestep shrink" begin
    using Piccolo
    using NamedTrajectories
    using DirectTrajOpt
    using DirectTrajOpt: CommonInterface
    using Random

    Random.seed!(3)

    M, order, n_d = 10, 4, 1
    T_ref = 2.0
    N = M - order + 2
    times = collect(range(0.0, T_ref, length = N))
    pulse = BSplinePulse(randn(n_d, M), times; order = order)
    c = vec(pulse.control_points)

    make_traj = Δt_val -> NamedTrajectory(
        (x = zeros(1, N), Δt = fill(Δt_val, 1, N)),
        (c_u = collect(c),);
        timestep = :Δt,
        controls = (:Δt,),
    )
    traj_nom = make_traj(T_ref / (N - 1))       # s = 1
    traj_half = make_traj(T_ref / (N - 1) / 2)  # s = 1/2

    for (make, p) in ((bspline_slew_constraint, 1), (bspline_accel_constraint, 2))
        rows = make(pulse, 1.0)                 # unit-bound rows: |A·c| per row
        Ac_abs = abs.(rows.A * c)
        # Feasible at s = 1, and the largest rows violate at s = 1/2.
        bound = 1.05 * maximum(Ac_abs)
        con = make(pulse, bound; free_time = true, T_ref = T_ref)

        vals_nom = zeros(con.dim)
        CommonInterface.evaluate!(vals_nom, con, traj_nom)
        @test all(vals_nom .<= 0)               # feasible at nominal duration

        vals_half = zeros(con.dim)
        CommonInterface.evaluate!(vals_half, con, traj_half)

        # The cap scales as s^p: every row tightens by exactly b·(1 − (1/2)^p).
        @test all(vals_half .- vals_nom .≈ bound * (1 - 0.5^p))

        # Violation appears exactly when |A·c| exceeds the scaled bound.
        viol = [max(vals_half[2r - 1], vals_half[2r]) > 0 for r in 1:length(Ac_abs)]
        @test viol == (Ac_abs .> bound * 0.5^p)
        @test any(viol)                         # the max row always violates
        @test !all(viol)                        # small rows stay feasible
    end
end
