# ============================================================================ #
# SplineType trait hierarchy
#
# Moved from spline_integrator_type.jl (pure move, #132) so the shared
# algorithm-data layer (alg_data.jl) can dispatch on SplineType — it is
# included before alg_data.jl in _spline_integrators.jl.
# ============================================================================ #

"""
    SplineType

Abstract type for spline interpolation methods. Concrete types:
- `LinearSpline`: Linear interpolation between knot points
- `CubicSpline`: Cubic Hermite interpolation with derivative control
"""
abstract type SplineType end

"""
    LinearSpline <: SplineType

Linear spline interpolation. Control derivatives are constrained by knot values.
"""
struct LinearSpline <: SplineType end

"""
    CubicSpline <: SplineType

Cubic Hermite spline interpolation. Control derivatives are independent DOFs.
"""
struct CubicSpline <: SplineType end

# Trait inference from pulse types
"""
    SplineType(pulse::AbstractPulse)

Infer the SplineType from a pulse object.
"""
SplineType(::LinearSplinePulse) = LinearSpline()
SplineType(::CubicSplinePulse) = CubicSpline()

# Convert to integer for backward compatibility (temporary)
@inline spline_order(::LinearSpline) = 1
@inline spline_order(::CubicSpline) = 3
spline_order(pulse::AbstractPulse) = spline_order(SplineType(pulse))

# ── The packed ODE-parameter layout, declared per SplineType (#338) ────────── #
#
# The matrix-free driver builds its per-knot parameter gather table from THIS
# declaration, so a third `SplineType` (B-splines: PR #276, #266/#267) extends
# the interface by adding two methods here and NOTHING in the driver. Before
# #338 the layout was open-coded as `if spline_order(𝒮) == 1 … else cubic …` in
# `get_param_indices`, `build_ode_params` and `build_ode_params!`, which is
# precisely the pattern a third type could not extend.

"""
    AbstractParamBlockRole

Which trajectory quantity a packed ODE-parameter block draws from. Resolved by
`param_block_component` and `param_block_carries_globals`, both
of which DISPATCH — so a new role is a new method beside the new
[`SplineType`](@ref), never a new branch in the driver.
"""
abstract type AbstractParamBlockRole end

"""
    ControlValueBlock <: AbstractParamBlockRole

A block drawn from the control VALUE component (`u_name`). Globals ride these
blocks: a global is constant across knots, so both endpoint blocks point at the
same global column and their contributions sum — which is the
dynamics-entering-global chain rule.
"""
struct ControlValueBlock <: AbstractParamBlockRole end

"""
    ControlDerivBlock <: AbstractParamBlockRole

A block drawn from the control DERIVATIVE component (`du_name`). Global
derivatives are identically zero, so those slots are zero-filled and carry no
trajectory column.
"""
struct ControlDerivBlock <: AbstractParamBlockRole end

"""
    ControlDeriv2Block <: AbstractParamBlockRole

A block drawn from the control SECOND-derivative component (`ddu_name`). No
in-tree spline type uses it yet; it is declared so the role vocabulary is
demonstrably open (a quintic-Hermite or B-spline packing needs it) rather than
a two-case enumeration.
"""
struct ControlDeriv2Block <: AbstractParamBlockRole end

"""
    param_blocks(S::SplineType) -> Tuple{Vararg{Tuple{AbstractParamBlockRole,Int}}}

The packed ODE-parameter block layout this spline type declares.

The parameter vector is `[block₁, block₂, …, Δtₖ, tₖ]`, where each block is
`u_dim`-wide (`u_dim = control_dim + global_dim`). Each entry of the returned
tuple is `(role, knot_offset)` with `knot_offset ∈ {0, 1}` selecting knot `k` or
knot `k+1`.

`LinearSpline` declares `[uₖ, uₖ₊₁]`; `CubicSpline` declares
`[uₖ, uₖ₊₁, duₖ, duₖ₊₁]`. A new spline type declares its own `spline_order` and
`param_blocks` and nothing else changes — `get_param_indices`,
`build_ode_params!` and the matrix-free driver's `fill_knot_params!` are all
generic over this declaration (#338 AC6). B-splines arrive via PR #276 and
#266/#267.
"""
param_blocks(::LinearSpline) = ((ControlValueBlock(), 0), (ControlValueBlock(), 1))
param_blocks(::CubicSpline) = (
    (ControlValueBlock(), 0),
    (ControlValueBlock(), 1),
    (ControlDerivBlock(), 0),
    (ControlDerivBlock(), 1),
)

"""
    n_param_blocks(S::SplineType) -> Int

Number of `u_dim`-wide blocks in the packed ODE parameter layout — `p_dim` is
`n_param_blocks(S) * u_dim`, and `Δt` / `t` follow at `p_dim + 1` / `p_dim + 2`.
"""
n_param_blocks(S::SplineType) = length(param_blocks(S))

"""
    param_block_carries_globals(role::AbstractParamBlockRole) -> Bool

Whether this block's global slots carry the global VALUES (`true`) or are
identically zero (`false`).
"""
param_block_carries_globals(::ControlValueBlock) = true
param_block_carries_globals(::AbstractParamBlockRole) = false

# ── Operator-aware helpers ─────────────────────────────────────────────────

"""
    _drift_matrix(H_drift) -> AbstractMatrix

Materialize H_drift to a matrix for use in `Isomorphisms.G` (Magnus path).
If H_drift is already an AbstractMatrix, return it directly.
If it's an AbstractDynamicsOperator, materialize it.
"""
_drift_matrix(H::AbstractMatrix) = H
_drift_matrix(H::AbstractDynamicsOperator) = materialize(H)
