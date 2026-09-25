# ============================================================================ #
# SplineConstraints — spline-shape constraint family (Piccolo home).
#
# Open-core slice 3c (harmoniqs/Piccolissimo.jl#431). This module owns the
# nonlinear knot/segment constraint types, the spline bound/extrema/slope
# family, the two stencil-table constraints (`CubicSplineBoundConstraint`,
# `HermiteSmoothAccelerationConstraint`), the functional-indexed
# `ConstraintStencilTable` with its host application kernels (ADR-0010), and
# the one-definition Hermite primitives — moved from Piccolissimo per
# docs/moved-file-review.md rows 20–25 of the open-core split's manifest.
#
# What STAYED in Piccolissimo (proprietary): the device stencil table (#333),
# knot-axis sharding (#334), the Altissimo backends that consume this surface,
# and `HermiteC2Regularizer` (row 26 struck by amendment 9 — deprecated there,
# removed on the later C2 arc; see that repo's
# `objectives/hermite_c2_regularizer.jl`). `HermiteBendingEnergyRegularizer`
# moved in slice 3a (#309) and lives in QuantumObjectives.
#
# No extension seams are needed this slice: the moved family is
# self-contained, and the dependency arrow runs Piccolissimo → Piccolo (the
# stays import this surface back — the 3a/3b convention).
# ============================================================================ #

module SplineConstraints

using SparseArrays
using ForwardDiff
using LinearAlgebra
using NamedTrajectories
using TrajectoryIndexingUtils
using TestItemRunner

using DirectTrajOpt
using DirectTrajOpt: AbstractNonlinearConstraint

# Import the common interface functions (defined in CommonInterface).
# These are the generic functions we extend with our concrete implementations.
# Imported WITHOUT export: interface functions come from CommonInterface (the
# drift-guard seam convention, see test/test_spline_reexport_seam.jl).
using DirectTrajOpt.CommonInterface
import DirectTrajOpt.CommonInterface:
    evaluate!, jacobian_structure, jacobian!, hessian_structure, hessian_of_lagrangian
import DirectTrajOpt.CommonInterface: eval_jacobian, eval_hessian_of_lagrangian

# Import constraint-specific functions (defined in DirectTrajOpt.Constraints).
# Note: These are constraint-specific, not part of CommonInterface.
import DirectTrajOpt.Constraints: get_full_jacobian, get_full_hessian
import DirectTrajOpt.Constraints: jacobian!, hessian_of_lagrangian!

# ── Family surface ───────────────────────────────────────────────────────── #
export OptimizedNonlinearKnotPointConstraint
export NonlinearSegmentConstraint
export CubicSplineExtremaConstraint,
    CubicSplineSufficientBoundConstraint, CubicSplineSlopeConstraint
export CubicSplineBoundConstraint
export HermiteSmoothAccelerationConstraint

# The functional-indexed stencil table (ADR-0010) and its generic application
# kernels. Consumers: `HermiteSmoothAccelerationConstraint` (#330, two-interval
# continuity rows) and `CubicSplineBoundConstraint` (#331, single-interval
# sample rows). The kernels were reused UNCHANGED by the second port — the
# table needed no widening.
export ConstraintStencilTable
export stencil_structure,
    stencil_fill_values!,
    stencil_assemble!,
    stencil_scatter_functional!,
    stencil_expand_rows!,
    stencil_coeff_range,
    stencil_functional_rows,
    stencil_n_entries,
    stencil_width

# Coefficient-refresh bookkeeping and the constraint-type gradient trait (#332).
# `supports_matrix_free_constraint_gradient` is the inequality-path sibling of
# the backend's `_supports_matrix_free_eq_gradient`; a constraint opts in by
# defining `constraint_stencil_table` + `refresh_constraint_coefficients!` and
# routability is then DERIVED from the declared stencil width. The application
# kernels' DEVICE half (#333) and the sharded path (#334) stayed in
# Piccolissimo and extend this surface from there.
export stencil_refresh_token,
    stencil_touch!,
    stencil_jvp!,
    stencil_vjp!,
    constraint_stencil_table,
    refresh_constraint_coefficients!,
    supports_matrix_free_constraint_gradient,
    UNBOUNDED_STENCIL_WIDTH

# The exact inequality HVP trait (#458): DECLARED per family (unlike the
# DERIVED gradient trait — no structural property of a table implies a
# second-order action), with the per-functional weight contraction shared by
# every opting-in family.
export supports_matrix_free_constraint_hvp,
    constraint_stencil_hvp!, stencil_functional_weight

# Cubic Hermite spline primitives — one definition, shared by every spline constraint.
export hermite_basis_functions,
    hermite_derivative_basis,
    evaluate_hermite_spline,
    evaluate_hermite_derivative,
    hermite_value_gradient,
    hermite_accel_start,
    hermite_accel_end,
    hermite_accel_start_gradient,
    hermite_accel_end_gradient,
    hermite_accel_jump_gradient

# ── Includes (definition order matters) ──────────────────────────────────── #
include("hermite_primitives.jl")
include("stencil_table.jl")
include("nonlinear_knot_point_constraint.jl")
include("nonlinear_segment_constraint.jl")
include("spline_bound_constraints.jl")
include("cubic_spline_bound_constraint.jl")
include("hermite_smooth_acceleration_constraint.jl")
# `cubic_spline_bound_constraint.jl` and `hermite_smooth_acceleration_constraint.jl`
# come after the primitives + table they build on; they are mutually independent
# (neither references the other's symbols), so their relative order is free.

# Create alias for the optimized version (exported; the bare struct name stays
# internal so `using Piccolo` keeps resolving DirectTrajOpt's
# `NonlinearKnotPointConstraint` without ambiguity).
const OptimizedNonlinearKnotPointConstraint = NonlinearKnotPointConstraint

end
