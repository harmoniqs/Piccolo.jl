"""
    HermiteSmoothAccelerationConstraint <: AbstractNonlinearConstraint

Unified constraint that enforces both acceleration bounds AND C² continuity for Hermite cubic splines.

This constraint combines:
1. **Acceleration bounds**: |a| ≤ a_max at all knot points
2. **C² continuity**: acceleration continuity at interior knots (no jumps)

By enforcing continuity, we only need to bound acceleration once per interior knot (instead of
twice for endpoints of adjacent segments), making this more efficient than separate constraints.

**Constraint structure:**

For N knots and n_drives control drives (all segments constrained):
- At knot 1 (start): bound a_start(segment 1) with 2×n_drives inequalities (±a_max)
- At knots 2 to N-1 (interior):
  - Continuity (double-sided): a_end(k-1) = a_start(k) → 2×(N-2)×n_drives inequalities
  - Bounds: |a| ≤ a_max where a is the common acceleration → 2×(N-2)×n_drives inequalities
- At knot N (end): bound a_end(segment N-1) with 2×n_drives inequalities (±a_max)

Total: `4 × n_segments × n_drives` inequalities, where `n_segments = length(segments)`.

**Acceleration formulas** (`hermite_accel_start` / `hermite_accel_end`):

At start of segment k (s=0, knot k):
    a_start(k) = (6/Δt_k²)(u_{k+1} - u_k) - (2/Δt_k)(2du_k + du_{k+1})

At end of segment k (s=1, knot k+1):
    a_end(k) = -(6/Δt_k²)(u_{k+1} - u_k) + (2/Δt_k)(du_k + 2du_{k+1})

# Constructor
```julia
HermiteSmoothAccelerationConstraint(
    traj::NamedTrajectory;
    a_max::Real,
    control_name::Symbol=:u,
    segments::AbstractVector{Int}=1:traj.N-1,
)
```

# Arguments
- `traj`: NamedTrajectory with cubic spline controls (requires :u, :du, and :Δt components)
- `a_max`: Maximum allowed acceleration magnitude
- `control_name`: Name of the control variable (default: :u)
- `segments`: Which spline segments to constrain (default: all). Continuity is enforced
  at the interior knots of each maximal contiguous run of `segments`; each run's leading
  `a_start` and every segment's `a_end` are bounded.

# Example
```julia
# Smooth acceleration with bounds ±10 MHz/μs²
constraint = HermiteSmoothAccelerationConstraint(
    traj;
    a_max=10.0
)
```

# Representation

The residual and Jacobian are served from a `ConstraintStencilTable` (ADR-0010): one
table row per unique scalar acceleration functional, the constant columns declared
once at construction, and only the coefficients refreshed per iterate. The kernels
read the trajectory's flat data vector at those precomputed columns through a function
barrier — never the named-component accessor, which is inferred as `Any` and boxes
~38.6 B per scalar read.

# Allocation

**Not** zero-allocation, and zero allocation is not the target. The interface hands the
constraint a `NamedTrajectory` whose `datavec` is an abstractly-typed field, so one
dynamic access costs a constant floor of ~80 B; the backend's own per-callback
trajectory wrapper costs ~3.3 kB, roughly forty times a whole kernel here. What IS
guaranteed, and gated, is **knot-flatness**: per-call allocation for both `evaluate!`
and the Jacobian is independent of the knot count. See ADR-0009 decision 3 and
ADR-0010 decision 6.

# Notes
- Replaces both HermiteAccelerationConstraint and HermiteAccelerationContinuityConstraint
- More efficient: fewer total constraints than using both separately
- Results in C² continuous splines with bounded acceleration
"""
struct HermiteSmoothAccelerationConstraint <: AbstractNonlinearConstraint
    control_name::Symbol
    derivative_name::Symbol
    timestep_name::Symbol
    n_drives::Int
    N_segments::Int
    segments::Vector{Int}
    a_max::Float64
    equality::Bool  # Always false (pure inequality formulation)
    g_dim::Int  # Total constraints: 4 * length(segments) * n_drives
    dim::Int  # Same as g_dim for this constraint

    # Functional-indexed stencil table: the declared structure (constant columns,
    # signed row map, per-row constant, stencil width) plus the per-iterate
    # coefficients. Application is generic over it.
    table::ConstraintStencilTable

    # Per-functional kind tag — which closed-form scalar this functional is. The
    # constraint's own math; the table neither reads nor interprets it.
    kinds::Vector{Int8}

    # Jacobian sparsity structure (derived from the table, cached for the interface)
    jac_structure::Vector{Tuple{Int,Int}}

    # Hessian sparsity structure
    hess_structure::Vector{Tuple{Int,Int}}
end

# ----------------------------------------------------------------------------- #
# Functional kinds
# ----------------------------------------------------------------------------- #
#
# Column layouts, fixed by the declaration below and read by the barrier kernels:
#
#   HSA_A_START (5): [u_k, u_{k+1}, du_k, du_{k+1}, Δt_k]
#   HSA_A_END   (5): [u_k, u_{k+1}, du_k, du_{k+1}, Δt_k]
#   HSA_JUMP    (8): [u_{k-1}, u_k, u_{k+1}, du_{k-1}, du_k, du_{k+1}, Δt_{k-1}, Δt_k]

const HSA_A_START = Int8(1)
const HSA_A_END = Int8(2)
const HSA_JUMP = Int8(3)

# Stencil width: a bound/continuity row reads at most knots k-1, k, k+1 around its
# anchor, i.e. one knot either side. DECLARED, never acted on here — the sharded
# driver takes `max` over dynamics and every routed constraint and exchanges once.
const HSA_STENCIL_WIDTH = 1

function HermiteSmoothAccelerationConstraint(
    traj::NamedTrajectory;
    a_max::Real,
    control_name::Symbol = :u,
    segments::AbstractVector{Int} = 1:(traj.N-1),
)
    # Validate trajectory components
    @assert haskey(traj.components, control_name) "Trajectory missing control: $control_name"

    derivative_name = Symbol("d" * string(control_name))
    @assert haskey(traj.components, derivative_name) "Trajectory missing derivative: $derivative_name"

    timestep_name = traj.timestep
    @assert !isnothing(timestep_name) "Trajectory missing timestep"

    n_drives = traj.dims[control_name]
    N_timesteps = traj.N

    segs = sort(unique(collect(Int, segments)))
    @assert !isempty(segs) "segments is empty — nothing to constrain"
    @assert first(segs) >= 1 && last(segs) <= N_timesteps - 1 "segments must lie in 1:$(N_timesteps-1)"
    N_segments = length(segs)
    in_seg = falses(N_timesteps - 1)
    for k in segs
        in_seg[k] = true
    end

    # All constraints as inequalities, `±` row pairs sharing one functional:
    #   continuity at every interior knot of a segment run, plus a bound on each
    #   run's leading a_start and on every segment's a_end.
    g_dim = 4 * N_segments * n_drives

    # ── Column tables. THE FULL DECISION-VECTOR WIDTH, INCLUDING GLOBALS, and it is
    # carried as data on the table rather than recomputed at each use site: sizing the
    # Jacobian without the global columns is what made the backend's row-stacking throw
    # outright on any trajectory carrying free phases.
    n_cols = traj.dim * traj.N + traj.global_dim
    u_comps = traj.components[control_name]
    du_comps = traj.components[derivative_name]
    dt_comps = traj.components[timestep_name]
    u_col = (k, d) -> slice(k, u_comps, traj.dim)[d]
    du_col = (k, d) -> slice(k, du_comps, traj.dim)[d]
    dt_col = k -> slice(k, dt_comps, traj.dim)[1]

    functional_cols = Vector{Vector{Int}}()
    kinds = Int8[]
    row_map = Int[]
    row_offset = Float64[]

    # A functional emitted with its `±` row pair — CONTIGUOUS by construction, which is
    # the invariant `ConstraintStencilTable` enforces.
    push_functional! = (kind, cols, offset) -> begin
        push!(functional_cols, cols)
        push!(kinds, kind)
        f = length(functional_cols)
        push!(row_map, f)
        push!(row_offset, offset)
        push!(row_map, -f)
        push!(row_offset, offset)
        return f
    end

    # ── Continuity: knot k is interior to a run iff segments k-1 and k are both in.
    # Row order matches the pre-port implementation exactly (all-segments default:
    # k = 2:N-1).
    for k in segs
        (k >= 2 && in_seg[k-1]) || continue
        for d = 1:n_drives
            push_functional!(
                HSA_JUMP,
                [
                    u_col(k - 1, d),
                    u_col(k, d),
                    u_col(k + 1, d),
                    du_col(k - 1, d),
                    du_col(k, d),
                    du_col(k + 1, d),
                    dt_col(k - 1),
                    dt_col(k),
                ],
                0.0,
            )
        end
    end

    # ── Bounds: each run's leading a_start (default: knot 1) …
    for k in segs
        (k == 1 || !in_seg[k-1]) || continue
        for d = 1:n_drives
            push_functional!(
                HSA_A_START,
                [u_col(k, d), u_col(k + 1, d), du_col(k, d), du_col(k + 1, d), dt_col(k)],
                -Float64(a_max),
            )
        end
    end

    # … then a_end of every non-final segment of a run (default: interior knots 2:N-1,
    # whose common acceleration is a_end(k-1)) …
    for k in segs
        (k + 1 <= N_timesteps - 1 && in_seg[k+1]) || continue
        for d = 1:n_drives
            push_functional!(
                HSA_A_END,
                [u_col(k, d), u_col(k + 1, d), du_col(k, d), du_col(k + 1, d), dt_col(k)],
                -Float64(a_max),
            )
        end
    end

    # … then a_end of each run's final segment (default: knot N).
    for k in segs
        (k + 1 > N_timesteps - 1 || !in_seg[k+1]) || continue
        for d = 1:n_drives
            push_functional!(
                HSA_A_END,
                [u_col(k, d), u_col(k + 1, d), du_col(k, d), du_col(k + 1, d), dt_col(k)],
                -Float64(a_max),
            )
        end
    end

    @assert length(row_map) == g_dim "internal: emitted $(length(row_map)) rows, expected $g_dim"

    table = ConstraintStencilTable(
        functional_cols,
        row_map,
        row_offset,
        n_cols;
        stencil_width = HSA_STENCIL_WIDTH,
    )

    # Build Hessian sparsity structure
    # Include second derivatives involving Δt for:
    # 1. Within-segment terms (for bound constraints)
    # 2. Cross-segment terms (for continuity constraints)
    hess_structure = Tuple{Int,Int}[]

    # Within-segment terms, for each constrained segment
    for k in segs
        for drive_idx = 1:n_drives
            var_indices = [
                u_col(k, drive_idx),
                u_col(k + 1, drive_idx),
                du_col(k, drive_idx),
                du_col(k + 1, drive_idx),
                dt_col(k),
            ]

            for i = 1:length(var_indices)
                for j = 1:i
                    idx_i = var_indices[i]
                    idx_j = var_indices[j]
                    push!(hess_structure, (max(idx_i, idx_j), min(idx_i, idx_j)))
                end
            end
        end
    end

    # Cross-segment Hessian entries for continuity constraints:
    # Δt_{k-1} with variables from segment k, and Δt_k with variables from segment k-1
    for k in segs
        (k >= 2 && in_seg[k-1]) || continue
        dt_km1 = dt_col(k - 1)
        dt_k = dt_col(k)

        for drive_idx = 1:n_drives
            seg_k_vars = [
                u_col(k, drive_idx),
                u_col(k + 1, drive_idx),
                du_col(k, drive_idx),
                du_col(k + 1, drive_idx),
                dt_k,
            ]
            for var in seg_k_vars
                push!(hess_structure, (max(dt_km1, var), min(dt_km1, var)))
            end

            seg_km1_vars = [
                u_col(k - 1, drive_idx),
                u_col(k, drive_idx),
                du_col(k - 1, drive_idx),
                du_col(k, drive_idx),
                dt_km1,
            ]
            for var in seg_km1_vars
                push!(hess_structure, (max(dt_k, var), min(dt_k, var)))
            end
        end
    end

    # Remove duplicates and sort
    hess_structure = unique(hess_structure)
    sort!(hess_structure)

    return HermiteSmoothAccelerationConstraint(
        control_name,
        derivative_name,
        timestep_name,
        n_drives,
        N_segments,
        segs,
        Float64(a_max),
        false,  # equality field (pure inequality formulation)
        g_dim,
        g_dim,  # dim = g_dim
        table,
        kinds,
        stencil_structure(table),
        hess_structure,
    )
end

# CommonInterface methods
import DirectTrajOpt.CommonInterface: jacobian_structure, hessian_structure
import DirectTrajOpt.CommonInterface: evaluate!, jacobian!, hessian_of_lagrangian

jacobian_structure(c::HermiteSmoothAccelerationConstraint) = c.jac_structure
hessian_structure(c::HermiteSmoothAccelerationConstraint) = c.hess_structure

"""
    stencil_width(constraint::HermiteSmoothAccelerationConstraint) -> Int

The declared stencil width: the maximum number of knots either side of its anchor that
one constraint row reads. 1 here — a continuity row couples `k-1, k, k+1`.

Declared, never acted on: halo exchange belongs to the sharded driver, which takes
`max` over dynamics and every routed constraint and exchanges once per callback.
"""
stencil_width(c::HermiteSmoothAccelerationConstraint) = c.table.stencil_width

# ----------------------------------------------------------------------------- #
# Barrier kernels
# ----------------------------------------------------------------------------- #
#
# `Z` arrives as `traj.datavec`, an abstractly-typed FIELD — reading scalars through it
# at the call site (or through `traj[name]`, which is inferred as `Any`) boxes every
# intermediate and was ~98% of the pre-port residual allocation. These functions are
# the FUNCTION BARRIER: one dynamic dispatch at entry, then a fully concrete body.
# Passing component *views* through a barrier is measurably not sufficient — it leaves
# a constant kilobyte-scale floor.

"""
    _hsa_functional_value(kind, Z, cols, o)

The closed-form scalar for one functional: `kind` selects the formula, `cols[o+1…]`
are its precomputed decision-vector columns. Generic in `eltype(Z)`.
"""
@inline function _hsa_functional_value(
    kind::Int8,
    Z::AbstractVector,
    cols::Vector{Int},
    o::Int,
)
    @inbounds if kind == HSA_JUMP
        u_km1 = Z[cols[o+1]]
        u_k = Z[cols[o+2]]
        u_kp1 = Z[cols[o+3]]
        du_km1 = Z[cols[o+4]]
        du_k = Z[cols[o+5]]
        du_kp1 = Z[cols[o+6]]
        Δt_km1 = Z[cols[o+7]]
        Δt_k = Z[cols[o+8]]
        return hermite_accel_end(u_km1, u_k, du_km1, du_k, Δt_km1) -
               hermite_accel_start(u_k, u_kp1, du_k, du_kp1, Δt_k)
    else
        u_a = Z[cols[o+1]]
        u_b = Z[cols[o+2]]
        du_a = Z[cols[o+3]]
        du_b = Z[cols[o+4]]
        Δt = Z[cols[o+5]]
        return kind == HSA_A_START ? hermite_accel_start(u_a, u_b, du_a, du_b, Δt) :
               hermite_accel_end(u_a, u_b, du_a, du_b, Δt)
    end
end

"""
    _hsa_residual!(g, Z, kinds, table)

Residual barrier: form each functional and scatter it straight into its (contiguous)
`±` row pair. No per-functional buffer, so this is generic in `eltype(g)` /
`eltype(Z)` — forward-mode AD of this very function is the Jacobian oracle.
"""
function _hsa_residual!(
    g::AbstractVector,
    Z::AbstractVector,
    kinds::Vector{Int8},
    table::ConstraintStencilTable,
)
    cols = table.cols
    col_ptr = table.col_ptr
    @inbounds for f in eachindex(kinds)
        val = _hsa_functional_value(kinds[f], Z, cols, col_ptr[f] - 1)
        stencil_scatter_functional!(g, table, f, val)
    end
    return nothing
end

"""
    _hsa_refresh_coefficients!(Z, kinds, table)

Per-iterate coefficient refresh: overwrite `table.coeffs` with ∂functional/∂column, in
the declared column order. The `±` row signs are the table's business, not this
kernel's.
"""
function _hsa_refresh_coefficients!(
    Z::AbstractVector,
    kinds::Vector{Int8},
    table::ConstraintStencilTable,
)
    cols = table.cols
    col_ptr = table.col_ptr
    coeffs = table.coeffs
    @inbounds for f in eachindex(kinds)
        o = col_ptr[f] - 1
        kind = kinds[f]
        if kind == HSA_JUMP
            u_km1 = Z[cols[o+1]]
            u_k = Z[cols[o+2]]
            u_kp1 = Z[cols[o+3]]
            du_km1 = Z[cols[o+4]]
            du_k = Z[cols[o+5]]
            du_kp1 = Z[cols[o+6]]
            i_km1 = 1.0 / Z[cols[o+7]]
            i_k = 1.0 / Z[cols[o+8]]
            i_km1_2 = i_km1 * i_km1
            i_km1_3 = i_km1_2 * i_km1
            i_k_2 = i_k * i_k
            i_k_3 = i_k_2 * i_k
            # ∂/∂(u_{k-1}, u_k, u_{k+1}, du_{k-1}, du_k, du_{k+1}, Δt_{k-1}, Δt_k)
            coeffs[o+1] = 6 * i_km1_2
            coeffs[o+2] = -6 * i_km1_2 + 6 * i_k_2
            coeffs[o+3] = -6 * i_k_2
            coeffs[o+4] = 2 * i_km1
            coeffs[o+5] = 4 * i_km1 + 4 * i_k
            coeffs[o+6] = 2 * i_k
            coeffs[o+7] = 12 * i_km1_3 * (u_k - u_km1) - 2 * i_km1_2 * (du_km1 + 2 * du_k)
            coeffs[o+8] = 12 * i_k_3 * (u_kp1 - u_k) - 2 * i_k_2 * (2 * du_k + du_kp1)
        else
            u_a = Z[cols[o+1]]
            u_b = Z[cols[o+2]]
            du_a = Z[cols[o+3]]
            du_b = Z[cols[o+4]]
            i = 1.0 / Z[cols[o+5]]
            i2 = i * i
            i3 = i2 * i
            if kind == HSA_A_START
                # a_start = 6/Δt²(u_b - u_a) - 2/Δt(2 du_a + du_b)
                coeffs[o+1] = -6 * i2
                coeffs[o+2] = 6 * i2
                coeffs[o+3] = -4 * i
                coeffs[o+4] = -2 * i
                coeffs[o+5] = -12 * i3 * (u_b - u_a) + 2 * i2 * (2 * du_a + du_b)
            else
                # a_end = -6/Δt²(u_b - u_a) + 2/Δt(du_a + 2 du_b)
                coeffs[o+1] = 6 * i2
                coeffs[o+2] = -6 * i2
                coeffs[o+3] = 2 * i
                coeffs[o+4] = 4 * i
                coeffs[o+5] = 12 * i3 * (u_b - u_a) - 2 * i2 * (du_a + 2 * du_b)
            end
        end
    end
    stencil_touch!(table)
    return nothing
end

# ── The constraint-type gradient trait (#332) ─────────────────────────────────── #
#
# Two methods, and `supports_matrix_free_constraint_gradient` follows from them: this
# constraint declares a bounded stencil width (a continuity row couples `k-1, k, k+1`,
# i.e. one knot either side), so the backend's inequality path serves its rows from
# `stencil_jvp!` / `stencil_vjp!` and gives it no block in the assembled matrix.
# `eval_jacobian` below is untouched and stays the assembled contract for Ipopt/MadNLP.

constraint_stencil_table(c::HermiteSmoothAccelerationConstraint) = c.table

refresh_constraint_coefficients!(
    c::HermiteSmoothAccelerationConstraint,
    traj::NamedTrajectory,
) = _hsa_refresh_coefficients!(traj.datavec, c.kinds, c.table)

"""
    evaluate!(g, constraint, traj)

Evaluate acceleration constraints as pure inequalities (≤ 0):
- Continuity (double-sided): a_end(k-1) - a_start(k) ≤ 0 AND a_start(k) - a_end(k-1) ≤ 0
- Bounds (double-sided): a - a_max ≤ 0 AND -a - a_max ≤ 0
"""
function evaluate!(
    g::AbstractVector,
    constraint::HermiteSmoothAccelerationConstraint,
    traj::NamedTrajectory,
)
    _hsa_residual!(g, traj.datavec, constraint.kinds, constraint.table)
    return nothing
end

"""
    jacobian!(jac, constraint, traj)

Compute the Jacobian's nonzero values, in `jacobian_structure` order.
"""
function jacobian!(
    jac::AbstractVector,
    constraint::HermiteSmoothAccelerationConstraint,
    traj::NamedTrajectory,
)
    _hsa_refresh_coefficients!(traj.datavec, constraint.kinds, constraint.table)
    stencil_fill_values!(jac, constraint.table)
    return nothing
end

# The `± a_end(k) - a_max` Hessian block, shared by the two bound emissions in
# `hessian_of_lagrangian!` below (previously written out twice under `k-1`/`k` aliases).
function _hsa_add_a_end_hessian!(
    add_hess!,
    μ,
    constraint_idx0,
    u,
    du,
    Δt,
    k,
    n_drives,
    u_col,
    du_col,
    dt_col,
)
    inv_Δt = 1.0 / Δt[k]
    inv_Δt_2 = inv_Δt * inv_Δt
    inv_Δt_3 = inv_Δt_2 * inv_Δt
    inv_Δt_4 = inv_Δt_3 * inv_Δt

    idx = constraint_idx0
    for drive_idx = 1:n_drives
        @inbounds begin
            u_a = u[drive_idx, k]
            u_b = u[drive_idx, k+1]
            du_a = du[drive_idx, k]
            du_b = du[drive_idx, k+1]
        end

        u_a_v = u_col(k, drive_idx)
        u_b_v = u_col(k + 1, drive_idx)
        du_a_v = du_col(k, drive_idx)
        du_b_v = du_col(k + 1, drive_idx)
        dt_v = dt_col(k)

        # g = a_end - a_max
        @inbounds μ_pos = μ[idx]
        add_hess!(
            dt_v,
            dt_v,
            μ_pos * (-36 * inv_Δt_4 * (u_b - u_a) + 4 * inv_Δt_3 * (du_a + 2 * du_b)),
        )
        add_hess!(u_a_v, dt_v, μ_pos * (-12*inv_Δt_3))
        add_hess!(u_b_v, dt_v, μ_pos * (12*inv_Δt_3))
        add_hess!(du_a_v, dt_v, μ_pos * (-2*inv_Δt_2))
        add_hess!(du_b_v, dt_v, μ_pos * (-4*inv_Δt_2))
        idx += 1

        # g = -a_end - a_max
        @inbounds μ_neg = μ[idx]
        add_hess!(
            dt_v,
            dt_v,
            -μ_neg * (-36 * inv_Δt_4 * (u_b - u_a) + 4 * inv_Δt_3 * (du_a + 2 * du_b)),
        )
        add_hess!(u_a_v, dt_v, -μ_neg * (-12*inv_Δt_3))
        add_hess!(u_b_v, dt_v, -μ_neg * (12*inv_Δt_3))
        add_hess!(du_a_v, dt_v, -μ_neg * (-2*inv_Δt_2))
        add_hess!(du_b_v, dt_v, -μ_neg * (-4*inv_Δt_2))
        idx += 1
    end

    return nothing
end

"""
    hessian_of_lagrangian!(hess, μ, constraint, traj)

Compute Hessian of Lagrangian weighted by multipliers μ.

The acceleration formulas are:
- a_start = (6/Δt²)(u_{k+1} - u_k) - (2/Δt)(2*du_k + du_{k+1})
- a_end = -(6/Δt²)(u_{k+1} - u_k) + (2/Δt)(du_k + 2*du_{k+1})

Second derivatives (all involve Δt since a is linear in u, du):
- ∂²a_start/∂Δt² = 36/Δt⁴ * (u_{k+1} - u_k) - 4/Δt³ * (2*du_k + du_{k+1})
- ∂²a_start/∂u_k∂Δt = 12/Δt³
- ∂²a_start/∂u_{k+1}∂Δt = -12/Δt³
- ∂²a_start/∂du_k∂Δt = 4/Δt²
- ∂²a_start/∂du_{k+1}∂Δt = 2/Δt²

- ∂²a_end/∂Δt² = -36/Δt⁴ * (u_{k+1} - u_k) + 4/Δt³ * (du_k + 2*du_{k+1})
- ∂²a_end/∂u_k∂Δt = -12/Δt³
- ∂²a_end/∂u_{k+1}∂Δt = 12/Δt³
- ∂²a_end/∂du_k∂Δt = -2/Δt²
- ∂²a_end/∂du_{k+1}∂Δt = -4/Δt²

Deliberately NOT ported onto the stencil table: the backend does not use a Hessian for
algebraic constraints (its inequality path passes
`compute_inequality_constraint_hvp! = nothing`), so this stays on the upstream
`hess_structure` path for the interior-point solvers that do.
"""
@views function hessian_of_lagrangian!(
    hess::AbstractVector,
    μ::AbstractVector,
    constraint::HermiteSmoothAccelerationConstraint,
    traj::NamedTrajectory,
)
    u = traj[constraint.control_name]
    du = traj[constraint.derivative_name]
    Δt = traj[constraint.timestep_name]

    n_drives = constraint.n_drives
    hess_struct = constraint.hess_structure
    segs = constraint.segments
    N_timesteps = traj.N
    in_seg = falses(N_timesteps - 1)
    for k in segs
        in_seg[k] = true
    end

    fill!(hess, 0.0)

    # Build lookup table for O(1) Hessian index access
    hess_lookup = Dict{Tuple{Int,Int},Int}()
    for (idx, pair) in enumerate(hess_struct)
        hess_lookup[pair] = idx
    end

    # Helper to add to Hessian (lower triangular format: row >= col)
    add_hess! = (i, j, val) -> begin
        key = (max(i, j), min(i, j))
        idx = get(hess_lookup, key, nothing)
        if !isnothing(idx)
            @inbounds hess[idx] += val
        end
        return nothing
    end

    u_comps = traj.components[constraint.control_name]
    du_comps = traj.components[constraint.derivative_name]
    dt_comps = traj.components[constraint.timestep_name]
    u_col = (k, d) -> slice(k, u_comps, traj.dim)[d]
    du_col = (k, d) -> slice(k, du_comps, traj.dim)[d]
    dt_col = k -> slice(k, dt_comps, traj.dim)[1]

    constraint_idx = 1

    # === Continuity inequalities (double-sided), in the same row order as the
    # residual: interior knots of each segment run ===
    for k in segs
        (k >= 2 && in_seg[k-1]) || continue

        inv_Δt_km1 = 1.0 / Δt[k-1]
        inv_Δt_k = 1.0 / Δt[k]
        inv_Δt_km1_2 = inv_Δt_km1 * inv_Δt_km1
        inv_Δt_km1_3 = inv_Δt_km1_2 * inv_Δt_km1
        inv_Δt_km1_4 = inv_Δt_km1_3 * inv_Δt_km1
        inv_Δt_k_2 = inv_Δt_k * inv_Δt_k
        inv_Δt_k_3 = inv_Δt_k_2 * inv_Δt_k
        inv_Δt_k_4 = inv_Δt_k_3 * inv_Δt_k

        for drive_idx = 1:n_drives
            @inbounds begin
                u_km1_val = u[drive_idx, k-1]
                u_k_val = u[drive_idx, k]
                u_kp1_val = u[drive_idx, k+1]
                du_km1_val = du[drive_idx, k-1]
                du_k_val = du[drive_idx, k]
                du_kp1_val = du[drive_idx, k+1]
            end

            u_km1_v = u_col(k - 1, drive_idx)
            u_k_v = u_col(k, drive_idx)
            u_kp1_v = u_col(k + 1, drive_idx)
            du_km1_v = du_col(k - 1, drive_idx)
            du_k_v = du_col(k, drive_idx)
            du_kp1_v = du_col(k + 1, drive_idx)
            dt_km1_v = dt_col(k - 1)
            dt_k_v = dt_col(k)

            # Forward inequality: g = a_end(k-1) - a_start(k)
            @inbounds μ_fwd = μ[constraint_idx]

            # ∂²a_end(k-1)/∂Δt_{k-1}²
            add_hess!(
                dt_km1_v,
                dt_km1_v,
                μ_fwd * (
                    -36 * inv_Δt_km1_4 * (u_k_val - u_km1_val) +
                    4 * inv_Δt_km1_3 * (du_km1_val + 2*du_k_val)
                ),
            )
            # -∂²a_start(k)/∂Δt_k²
            add_hess!(
                dt_k_v,
                dt_k_v,
                μ_fwd * (-(
                    36 * inv_Δt_k_4 * (u_kp1_val - u_k_val) -
                    4 * inv_Δt_k_3 * (2*du_k_val + du_kp1_val)
                )),
            )

            # ∂²a_end(k-1)/∂var∂Δt_{k-1}
            add_hess!(u_km1_v, dt_km1_v, μ_fwd * (-12*inv_Δt_km1_3))
            add_hess!(u_k_v, dt_km1_v, μ_fwd * (12*inv_Δt_km1_3))
            add_hess!(du_km1_v, dt_km1_v, μ_fwd * (-2*inv_Δt_km1_2))
            add_hess!(du_k_v, dt_km1_v, μ_fwd * (-4*inv_Δt_km1_2))

            # -∂²a_start(k)/∂var∂Δt_k
            add_hess!(u_k_v, dt_k_v, μ_fwd * (-12*inv_Δt_k_3))
            add_hess!(u_kp1_v, dt_k_v, μ_fwd * (12*inv_Δt_k_3))
            add_hess!(du_k_v, dt_k_v, μ_fwd * (-4*inv_Δt_k_2))
            add_hess!(du_kp1_v, dt_k_v, μ_fwd * (-2*inv_Δt_k_2))

            constraint_idx += 1

            # Backward inequality: g = -g_fwd, so Hessian is negated
            @inbounds μ_bwd = μ[constraint_idx]

            add_hess!(
                dt_km1_v,
                dt_km1_v,
                -μ_bwd * (
                    -36 * inv_Δt_km1_4 * (u_k_val - u_km1_val) +
                    4 * inv_Δt_km1_3 * (du_km1_val + 2*du_k_val)
                ),
            )
            add_hess!(
                dt_k_v,
                dt_k_v,
                -μ_bwd * (-(
                    36 * inv_Δt_k_4 * (u_kp1_val - u_k_val) -
                    4 * inv_Δt_k_3 * (2*du_k_val + du_kp1_val)
                )),
            )
            add_hess!(u_km1_v, dt_km1_v, -μ_bwd * (-12*inv_Δt_km1_3))
            add_hess!(u_k_v, dt_km1_v, -μ_bwd * (12*inv_Δt_km1_3))
            add_hess!(du_km1_v, dt_km1_v, -μ_bwd * (-2*inv_Δt_km1_2))
            add_hess!(du_k_v, dt_km1_v, -μ_bwd * (-4*inv_Δt_km1_2))
            add_hess!(u_k_v, dt_k_v, -μ_bwd * (-12*inv_Δt_k_3))
            add_hess!(u_kp1_v, dt_k_v, -μ_bwd * (12*inv_Δt_k_3))
            add_hess!(du_k_v, dt_k_v, -μ_bwd * (-4*inv_Δt_k_2))
            add_hess!(du_kp1_v, dt_k_v, -μ_bwd * (-2*inv_Δt_k_2))

            constraint_idx += 1
        end
    end

    # === Bound inequalities: each run's leading a_start ===
    for k in segs
        (k == 1 || !in_seg[k-1]) || continue

        inv_Δt_k = 1.0 / Δt[k]
        inv_Δt_k_2 = inv_Δt_k * inv_Δt_k
        inv_Δt_k_3 = inv_Δt_k_2 * inv_Δt_k
        inv_Δt_k_4 = inv_Δt_k_3 * inv_Δt_k

        for drive_idx = 1:n_drives
            @inbounds begin
                u_k_val = u[drive_idx, k]
                u_kp1_val = u[drive_idx, k+1]
                du_k_val = du[drive_idx, k]
                du_kp1_val = du[drive_idx, k+1]
            end

            u_k_v = u_col(k, drive_idx)
            u_kp1_v = u_col(k + 1, drive_idx)
            du_k_v = du_col(k, drive_idx)
            du_kp1_v = du_col(k + 1, drive_idx)
            dt_k_v = dt_col(k)

            # g = a_start - a_max, so ∂²g = ∂²a_start
            @inbounds μ_pos = μ[constraint_idx]
            add_hess!(
                dt_k_v,
                dt_k_v,
                μ_pos * (
                    36 * inv_Δt_k_4 * (u_kp1_val - u_k_val) -
                    4 * inv_Δt_k_3 * (2*du_k_val + du_kp1_val)
                ),
            )
            add_hess!(u_k_v, dt_k_v, μ_pos * (12*inv_Δt_k_3))
            add_hess!(u_kp1_v, dt_k_v, μ_pos * (-12*inv_Δt_k_3))
            add_hess!(du_k_v, dt_k_v, μ_pos * (4*inv_Δt_k_2))
            add_hess!(du_kp1_v, dt_k_v, μ_pos * (2*inv_Δt_k_2))
            constraint_idx += 1

            # g = -a_start - a_max, so ∂²g = -∂²a_start
            @inbounds μ_neg = μ[constraint_idx]
            add_hess!(
                dt_k_v,
                dt_k_v,
                -μ_neg * (
                    36 * inv_Δt_k_4 * (u_kp1_val - u_k_val) -
                    4 * inv_Δt_k_3 * (2*du_k_val + du_kp1_val)
                ),
            )
            add_hess!(u_k_v, dt_k_v, -μ_neg * (12*inv_Δt_k_3))
            add_hess!(u_kp1_v, dt_k_v, -μ_neg * (-12*inv_Δt_k_3))
            add_hess!(du_k_v, dt_k_v, -μ_neg * (4*inv_Δt_k_2))
            add_hess!(du_kp1_v, dt_k_v, -μ_neg * (2*inv_Δt_k_2))
            constraint_idx += 1
        end
    end

    # === Bound inequalities: a_end of every non-final segment of a run ===
    for k in segs
        (k + 1 <= N_timesteps - 1 && in_seg[k+1]) || continue
        _hsa_add_a_end_hessian!(
            add_hess!,
            μ,
            constraint_idx,
            u,
            du,
            Δt,
            k,
            n_drives,
            u_col,
            du_col,
            dt_col,
        )
        constraint_idx += 2 * n_drives
    end

    # === Bound inequalities: a_end of each run's final segment ===
    for k in segs
        (k + 1 > N_timesteps - 1 || !in_seg[k+1]) || continue
        _hsa_add_a_end_hessian!(
            add_hess!,
            μ,
            constraint_idx,
            u,
            du,
            Δt,
            k,
            n_drives,
            u_col,
            du_col,
            dt_col,
        )
        constraint_idx += 2 * n_drives
    end

    return nothing
end

# Evaluation wrappers for DirectTrajOpt
"""
    eval_jacobian(constraint, traj) -> SparseMatrixCSC

The ASSEMBLED Jacobian entry point — unchanged in contract (this is what Ipopt/MadNLP
and the backend's row-stacking consume). Served from the stencil table's cached sparse
matrix: the pattern is built once at construction and only `nzval` is refreshed, so the
returned matrix is the same object on every call and the call allocates nothing that
grows with the knot count.

`size(J) == (g_dim, traj.dim * traj.N + traj.global_dim)` — the FULL decision-vector
width. Sizing it without the global columns made the row-stacking throw on any
trajectory carrying free phases.
"""
function DirectTrajOpt.CommonInterface.eval_jacobian(
    constraint::HermiteSmoothAccelerationConstraint,
    traj::NamedTrajectory,
)
    _hsa_refresh_coefficients!(traj.datavec, constraint.kinds, constraint.table)
    return stencil_assemble!(constraint.table)
end

function DirectTrajOpt.CommonInterface.eval_hessian_of_lagrangian(
    constraint::HermiteSmoothAccelerationConstraint,
    traj::NamedTrajectory,
    μ::AbstractVector,
)
    hess = zeros(length(constraint.hess_structure))
    hessian_of_lagrangian!(hess, μ, constraint, traj)

    # Convert to sparse matrix (lower triangular storage)
    I_rows = [pair[1] for pair in constraint.hess_structure]
    J_cols = [pair[2] for pair in constraint.hess_structure]

    # Full decision-vector width, matching `eval_jacobian` and DirectTrajOpt's
    # `Z_dim = traj.dim * traj.N + traj.global_dim` convention.
    n = constraint.table.n_cols
    M_lower = sparse(I_rows, J_cols, hess, n, n)

    # Make symmetric: M = M_lower + M_lower' - diag(M_lower)
    # This properly handles off-diagonal elements
    M_sym = M_lower + M_lower' - spdiagm(diag(M_lower))

    # Return sparse matrix for Ipopt compatibility
    return M_sym
end

# ============================================================================= #
#   Reference (pre-port) residual — the residual-semantics witness               #
# ============================================================================= #
#
# Kept reachable so the port's residual can be checked against the implementation it
# replaces, on randomised trajectories (#330 AC1). Reads the trajectory through the
# named-component accessor exactly as the pre-port code did — which is why it is
# test-only and on no hot path. All-segments only: the pre-port constraint had no
# `segments` option.

"""
    _reference_evaluate!(g, constraint, traj)

The pre-port residual, verbatim in formulation and row order. Test-only witness for
"residual semantics unchanged"; requires `constraint.segments == 1:traj.N-1`.
"""
@views function _reference_evaluate!(
    g::AbstractVector,
    constraint::HermiteSmoothAccelerationConstraint,
    traj::NamedTrajectory,
)
    @assert constraint.segments == collect(1:(traj.N-1)) "reference residual covers all segments only"

    u = traj[constraint.control_name]
    du = traj[constraint.derivative_name]
    Δt = traj[constraint.timestep_name]

    N_timesteps = traj.N
    n_drives = constraint.n_drives
    a_max = constraint.a_max

    idx = 1

    # === Continuity inequalities (double-sided) at interior knots ===
    for k = 2:(N_timesteps-1)
        inv_Δt_km1 = 1.0 / Δt[k-1]
        inv_Δt_k = 1.0 / Δt[k]
        inv_Δt_km1_2 = inv_Δt_km1 * inv_Δt_km1
        inv_Δt_k_2 = inv_Δt_k * inv_Δt_k

        for drive_idx = 1:n_drives
            u_km1 = u[drive_idx, k-1]
            u_k = u[drive_idx, k]
            u_kp1 = u[drive_idx, k+1]
            du_km1 = du[drive_idx, k-1]
            du_k = du[drive_idx, k]
            du_kp1 = du[drive_idx, k+1]

            a_end_km1 =
                -6 * inv_Δt_km1_2 * (u_k - u_km1) + 2 * inv_Δt_km1 * (du_km1 + 2*du_k)
            a_start_k = 6 * inv_Δt_k_2 * (u_kp1 - u_k) - 2 * inv_Δt_k * (2*du_k + du_kp1)

            g[idx] = a_end_km1 - a_start_k  # ≤ 0
            idx += 1
            g[idx] = a_start_k - a_end_km1  # ≤ 0
            idx += 1
        end
    end

    # === Bound inequalities at all knots ===

    # Knot 1: a_start(segment 1)
    inv_Δt_1 = 1.0 / Δt[1]
    inv_Δt_1_2 = inv_Δt_1 * inv_Δt_1
    for drive_idx = 1:n_drives
        u_1 = u[drive_idx, 1]
        u_2 = u[drive_idx, 2]
        du_1 = du[drive_idx, 1]
        du_2 = du[drive_idx, 2]

        a = 6 * inv_Δt_1_2 * (u_2 - u_1) - 2 * inv_Δt_1 * (2*du_1 + du_2)

        g[idx] = a - a_max      # ≤ 0
        idx += 1
        g[idx] = -a - a_max     # ≤ 0
        idx += 1
    end

    # Interior knots: common acceleration (use a_end formula)
    for k = 2:(N_timesteps-1)
        inv_Δt_km1 = 1.0 / Δt[k-1]
        inv_Δt_km1_2 = inv_Δt_km1 * inv_Δt_km1

        for drive_idx = 1:n_drives
            u_km1 = u[drive_idx, k-1]
            u_k = u[drive_idx, k]
            du_km1 = du[drive_idx, k-1]
            du_k = du[drive_idx, k]

            a = -6 * inv_Δt_km1_2 * (u_k - u_km1) + 2 * inv_Δt_km1 * (du_km1 + 2*du_k)

            g[idx] = a - a_max      # ≤ 0
            idx += 1
            g[idx] = -a - a_max     # ≤ 0
            idx += 1
        end
    end

    # Knot N: a_end(segment N-1)
    k = N_timesteps - 1
    inv_Δt_k = 1.0 / Δt[k]
    inv_Δt_k_2 = inv_Δt_k * inv_Δt_k
    for drive_idx = 1:n_drives
        u_k = u[drive_idx, k]
        u_kp1 = u[drive_idx, k+1]
        du_k = du[drive_idx, k]
        du_kp1 = du[drive_idx, k+1]

        a = -6 * inv_Δt_k_2 * (u_kp1 - u_k) + 2 * inv_Δt_k * (du_k + 2*du_kp1)

        g[idx] = a - a_max      # ≤ 0
        idx += 1
        g[idx] = -a - a_max     # ≤ 0
        idx += 1
    end

    return nothing
end

@testitem "HermiteSmoothAccelerationConstraint derivative validation" begin
    using NamedTrajectories
    using Piccolo
    using DirectTrajOpt: test_constraint
    using LinearAlgebra
    using Random

    # Seeded: the μ∂²g FD check sits near the tolerance (Hessian entries carry
    # 1/Δt³ factors, so FD noise on O(10–100) entries approaches the 1e-3 atol);
    # an unseeded fixture made this intermittently red in full-suite runs.
    Random.seed!(1)

    # Create trajectory with cubic spline pulse
    N = 6
    T = 2.0
    times = collect(range(0, T, N))
    n_drives = 2
    u = randn(n_drives, N)
    du = randn(n_drives, N) * 0.5

    pulse = CubicSplinePulse(u, du, times)
    sys = TransmonSystem(levels = 3)
    U_goal = Matrix{ComplexF64}(I, 3, 3)
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)
    traj = NamedTrajectory(qtraj, N)

    # Create constraint
    constraint = HermiteSmoothAccelerationConstraint(traj; a_max = 10.0)

    # Validate Jacobian and Hessian against finite differences
    # Use atol=1e-3 to accommodate finite difference precision on diagonal Hessian terms
    test_constraint(constraint, traj; atol = 1e-3, show_jacobian_diff = true)
end

@testitem "#330 AC1/AC2: ported residual matches the reference; Jacobian matches forward-mode AD" begin
    using NamedTrajectories

    using Piccolo: HermiteSmoothAccelerationConstraint, eval_jacobian
    using DirectTrajOpt.CommonInterface: evaluate!, jacobian_structure
    using DirectTrajOpt.Constraints: jacobian!
    using Piccolo.Control.QuantumConstraints.SplineConstraints: _reference_evaluate!
    using ForwardDiff
    using SparseArrays
    using LinearAlgebra
    using Random

    Random.seed!(0x330AC12)

    # Trajectory factory, with and without global variables. The globals case is the
    # one that matters for the DERIVATIVE oracle — the pre-port implementation sized
    # its Jacobian without the global columns, so parity against it proves nothing
    # there. Forward-mode AD OF THE PORTED RESIDUAL is the oracle: everything here is
    # closed-form, so a tolerance loose enough to absorb finite-difference truncation
    # error is loose enough to hide a wrong stencil coefficient.
    function fixture(N, n_drives; globals::Bool, segments = nothing)
        comps = (
            u = randn(n_drives, N),
            du = randn(n_drives, N) * 0.5,
            Δt = reshape(0.05 .+ 0.1 .* rand(N), 1, N),
        )
        traj =
            globals ?
            NamedTrajectory(
                comps;
                timestep = :Δt,
                controls = :u,
                global_data = [0.3, -0.8],
                global_components = (ϕ = 1:2,),
            ) : NamedTrajectory(comps; timestep = :Δt, controls = :u)
        c =
            isnothing(segments) ? HermiteSmoothAccelerationConstraint(traj; a_max = 7.0) :
            HermiteSmoothAccelerationConstraint(traj; a_max = 7.0, segments = segments)
        return traj, c
    end

    # ── AC1: residual semantics unchanged, to machine precision, on randomised
    # trajectories. The reference is the pre-port implementation, kept reachable.
    for trial = 1:4
        N = rand(4:9)
        nd = rand(1:3)
        traj, c = fixture(N, nd; globals = false)
        g_new = zeros(c.g_dim)
        g_ref = zeros(c.g_dim)
        evaluate!(g_new, c, traj)
        _reference_evaluate!(g_ref, c, traj)
        @test g_new == g_ref        # bit-identical: same formulation, same row order
    end

    # Same, with globals present (neither the formulation nor the reference reads them)
    traj_g, c_g = fixture(7, 2; globals = true)
    g_new = zeros(c_g.g_dim)
    g_ref = zeros(c_g.g_dim)
    evaluate!(g_new, c_g, traj_g)
    _reference_evaluate!(g_ref, c_g, traj_g)
    @test g_new == g_ref

    # ── AC2: the assembled Jacobian matches forward-mode AD of the ported residual.
    function ad_jacobian(traj, c)
        n_data = traj.dim * traj.N
        f =
            Z -> begin
                t = NamedTrajectory(
                    traj;
                    datavec = Z[1:n_data],
                    global_data = Z[(n_data+1):end],
                )
                g = zeros(eltype(Z), c.g_dim)
                evaluate!(g, c, t)
                return g
            end
        return ForwardDiff.jacobian(f, collect(vec(traj)))
    end

    for (label, fx) in (
        "no globals" => fixture(8, 2; globals = false),
        "with globals" => fixture(8, 2; globals = true),
        "segment subset" => fixture(10, 2; globals = true, segments = 3:7),
        "split segment runs" => fixture(12, 1; globals = false, segments = [1, 2, 3, 7, 8]),
    )
        traj, c = fx
        J = eval_jacobian(c, traj)
        J_ad = ad_jacobian(traj, c)
        @test size(J) == size(J_ad)      # includes the global columns
        @test size(J, 2) == traj.dim * traj.N + traj.global_dim
        @test maximum(abs, Matrix(J) - J_ad) < 1e-11 * max(1.0, maximum(abs, J_ad))

        # The flat-value entry point agrees with the assembled one, entry for entry.
        vals = zeros(length(jacobian_structure(c)))
        jacobian!(vals, c, traj)
        for (e, (r, col)) in enumerate(jacobian_structure(c))
            @test vals[e] == J[r, col]
        end
    end

    # A subset of segments really does constrain fewer rows, and the table enumerates
    # its functionals from those actual rows.
    traj_s, c_s = fixture(10, 2; globals = false, segments = 3:7)
    @test c_s.g_dim == 4 * 5 * 2
    @test c_s.table.n_functionals == 2 * 5 * 2
    @test c_s.table.n_rows == c_s.g_dim

    # ── AC5: the assembled entry point's sparsity structure — rows, columns AND their
    # ORDER — is exactly the pre-port one. `jacobian_structure` is the upstream
    # interior-point contract (Ipopt/MadNLP consume it positionally alongside
    # `jacobian!`'s value vector), so a reordering would be a silent miscoupling rather
    # than a test failure. Reconstructed here independently of the constructor, from the
    # pre-port emission order: all continuity rows (±, per interior knot, per drive),
    # then knot 1's a_start, then a_end at interior knots 2:N-1, then knot N's a_end.
    let N = 7, nd = 2
        traj, c = fixture(N, nd; globals = false)
        u_c = traj.components[:u]
        du_c = traj.components[:du]
        dt_c = traj.components[:Δt]
        ucol = (k, d) -> ((k - 1) * traj.dim) + u_c[d]
        dcol = (k, d) -> ((k - 1) * traj.dim) + du_c[d]
        tcol = k -> ((k - 1) * traj.dim) + dt_c[1]

        expected = Tuple{Int,Int}[]
        row = 0
        emit! = cols -> begin
            for _ = 1:2                       # the ± row pair
                row += 1
                for col in cols
                    push!(expected, (row, col))
                end
            end
        end
        for k = 2:(N-1), d = 1:nd
            emit!([
                ucol(k - 1, d),
                ucol(k, d),
                ucol(k + 1, d),
                dcol(k - 1, d),
                dcol(k, d),
                dcol(k + 1, d),
                tcol(k - 1),
                tcol(k),
            ])
        end
        for d = 1:nd
            emit!([ucol(1, d), ucol(2, d), dcol(1, d), dcol(2, d), tcol(1)])
        end
        for k = 2:(N-1), d = 1:nd
            emit!([ucol(k - 1, d), ucol(k, d), dcol(k - 1, d), dcol(k, d), tcol(k - 1)])
        end
        for d = 1:nd
            k = N - 1
            emit!([ucol(k, d), ucol(k + 1, d), dcol(k, d), dcol(k + 1, d), tcol(k)])
        end

        @test jacobian_structure(c) == expected
        @test row == c.g_dim
    end
end

@testitem "#330 AC3: residual and Jacobian allocation are invariant to knot count" begin
    using NamedTrajectories

    using Piccolo: HermiteSmoothAccelerationConstraint, eval_jacobian
    using DirectTrajOpt.CommonInterface: evaluate!
    using DirectTrajOpt.Constraints: jacobian!
    using Random

    # KNOT-FLATNESS, not zero allocation. The interface hands the constraint a
    # trajectory whose `datavec` is an abstractly-typed field, so one dynamic access
    # costs a constant floor no kernel work removes (ADR-0010 decision 6). What is
    # asserted is that the floor is CONSTANT: an implementation that still reads the
    # trajectory through the `Any`-inferred accessor, or rebuilds a hash lookup and a
    # whole sparse matrix per call, grows linearly in N and fails the `==`.
    #
    # Shape copied from the dynamics-HVP gates (`#337 alloc:` in
    # spline_integrator_ket.jl): local-scoped closures only (a callback at global scope
    # captures `Any` and inflates the count), and the MINIMUM over repeats because
    # `@allocated` sums across every thread in the process.
    function measure(N)
        Random.seed!(0x330A11C)
        n_drives = 2
        traj = NamedTrajectory(
            (
                u = randn(n_drives, N),
                du = randn(n_drives, N),
                Δt = reshape(fill(0.1, N), 1, N),
            );
            timestep = :Δt,
            controls = :u,
        )
        c = HermiteSmoothAccelerationConstraint(traj; a_max = 10.0)
        g = zeros(c.g_dim)
        vals = zeros(length(c.jac_structure))
        res = () -> evaluate!(g, c, traj)
        jacv = () -> jacobian!(vals, c, traj)
        asm = () -> eval_jacobian(c, traj)
        best = f -> (f(); f(); minimum(@allocated(f()) for _ = 1:5))
        return (residual = best(res), jacobian = best(jacv), assembled = best(asm))
    end

    N_SMALL, N_LARGE = 33, 129

    # Process-level warmup over BOTH sizes before either is measured: the FIRST
    # `@allocated` site in a process reports a different figure from every subsequent
    # one, and with `==` as the assertion whichever size went first would otherwise
    # carry that one-time step as phantom N-growth.
    measure(N_SMALL)
    measure(N_LARGE)

    a = measure(N_SMALL)
    b = measure(N_LARGE)

    for k in (:residual, :jacobian, :assembled)
        # INVARIANCE, not "less than before": a bound of that form is satisfied by an
        # implementation that still allocates linearly in N.
        @test getproperty(b, k) == getproperty(a, k)
        # …plus a loose absolute ceiling, so the constant floor cannot creep either.
        @test getproperty(b, k) <= 4096
    end
end

@testitem "HermiteSmoothAccelerationConstraint basic" begin
    using NamedTrajectories
    using Piccolo
    using LinearAlgebra

    N = 5
    T = 2.0
    times = collect(range(0, T, N))
    u = randn(2, N)
    du = randn(2, N)

    pulse = CubicSplinePulse(u, du, times)
    sys = TransmonSystem(levels = 3)
    U_goal = Matrix{ComplexF64}(I, 3, 3)
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)
    traj = NamedTrajectory(qtraj, N)

    constraint = HermiteSmoothAccelerationConstraint(traj; a_max = 10.0)

    # Check dimensions - pure inequality formulation:
    # - Continuity (double-sided): 2 * (N-2) * n_drives = 2 * 3 * 2 = 12
    # - Bounds (double-sided): 2 * N * n_drives = 2 * 5 * 2 = 20
    # - Total: 4 * (N-1) * n_drives = 4 * 4 * 2 = 32
    n_drives = 2
    N_segments = N - 1  # 4 segments
    expected_total = 4 * N_segments * n_drives  # 32

    @test constraint.g_dim == expected_total
    @test constraint.equality == false  # Pure inequality
    @test constraint.N_segments == N_segments

    # The stencil table: one functional per unique scalar, two rows each (±).
    @test constraint.table.n_rows == expected_total
    @test constraint.table.n_functionals == expected_total ÷ 2
    @test constraint.table.n_cols == traj.dim * traj.N + traj.global_dim
    @test constraint.table.stencil_width == 1
end

@testitem "HermiteSmoothAccelerationConstraint constraint counts" begin
    using NamedTrajectories
    using Piccolo
    using LinearAlgebra

    # Test different trajectory sizes
    for N in [4, 6, 10]
        T = 2.0
        times = collect(range(0, T, N))
        n_drives = 2
        u = randn(n_drives, N)
        du = randn(n_drives, N)

        pulse = CubicSplinePulse(u, du, times)
        sys = TransmonSystem(levels = 3)
        U_goal = Matrix{ComplexF64}(I, 3, 3)
        qtraj = UnitaryTrajectory(sys, pulse, U_goal)
        traj = NamedTrajectory(qtraj, N)

        constraint = HermiteSmoothAccelerationConstraint(traj; a_max = 10.0)

        # Expected counts - pure inequality formulation:
        # - Continuity (double-sided): 2 * (N-2) * n_drives
        # - Bounds (double-sided): 2 * N * n_drives
        # - Total: 4 * (N-1) * n_drives
        N_segments = N - 1
        expected_total = 4 * N_segments * n_drives

        @test constraint.equality == false
        @test constraint.g_dim == expected_total
        @test constraint.N_segments == N_segments

        # Test evaluation produces correct sized output
        g = zeros(constraint.g_dim)
        evaluate!(g, constraint, traj)
        @test length(g) == expected_total
        @test all(isfinite.(g))
    end
end

# ============================================================================= #
#   Exact inequality HVP on the stencil table (#458)                            #
# ============================================================================= #
#
# The assembled `hessian_of_lagrangian!` above stays on the `hess_structure` path for the
# interior-point solvers. The MATRIX-FREE backend path consumes the per-functional action
# below instead: `Σ_f ω_f ∇²F_f · v` over the table's functionals, the second-order
# analogue of the #332 pairing (a functional's `±` row pair shares ONE functional, so one
# weight contraction `ω_f = Σ_r sign_r w_r` serves both rows before any second-order work).
# Every closed form here is the one `hessian_of_lagrangian!` already carries, re-expressed
# as an ACTION — and the ForwardDiff oracle below is what keeps the two from drifting,
# exactly as the residual witness does for the Jacobian.

"""
    _hsa_functional_hessian_action!(Hv, v, ω, kind, Z, cols, o)

Add `ω · ∇²F_kind(x) · v` to `Hv` for one functional of the given kind, reading the
iterate scalars from `Z` at the declared columns `cols[o+1…]`. ACCUMULATES. Allocation-free.

The accelerations are linear in `u`/`du` and rational in `Δt`, so every second derivative
carries a `Δt` index — exactly the entries `hessian_of_lagrangian!` writes. All closed
forms follow that function (and the docstring above it); the ForwardDiff oracle testitem
is the drift guard.
"""
@inline function _hsa_functional_hessian_action!(
    Hv::AbstractVector{Float64},
    v::AbstractVector{Float64},
    ω::Float64,
    kind::Int8,
    Z::AbstractVector{Float64},
    cols::Vector{Int},
    o::Int,
)
    @inbounds if kind == HSA_JUMP
        # F = a_end(k-1) − a_start(k), columns [u_{k-1},u_k,u_{k+1},du_{k-1},du_k,du_{k+1},Δt_{k-1},Δt_k]
        u_km1 = Z[cols[o+1]]
        u_k = Z[cols[o+2]]
        u_kp1 = Z[cols[o+3]]
        du_km1 = Z[cols[o+4]]
        du_k = Z[cols[o+5]]
        du_kp1 = Z[cols[o+6]]
        i_km1 = 1.0 / Z[cols[o+7]]
        i_k = 1.0 / Z[cols[o+8]]
        i_km1_2 = i_km1 * i_km1
        i_km1_3 = i_km1_2 * i_km1
        i_km1_4 = i_km1_3 * i_km1
        i_k_2 = i_k * i_k
        i_k_3 = i_k_2 * i_k
        i_k_4 = i_k_3 * i_k
        vt_m1 = v[cols[o+7]]
        vt_k = v[cols[o+8]]
        # −a_start(k) part: every ∂²a_start entry negated.
        h_k = -(36 * i_k_4 * (u_kp1 - u_k) - 4 * i_k_3 * (2 * du_k + du_kp1))
        # a_end(k-1) part.
        h_m1 = (-36 * i_km1_4 * (u_k - u_km1) + 4 * i_km1_3 * (du_km1 + 2 * du_k))
        # Δt_{k-1} row/column (the a_end block only)
        Hv[cols[o+7]] +=
            ω * (
                h_m1 * vt_m1 - 12 * i_km1_3 * v[cols[o+1]] + 12 * i_km1_3 * v[cols[o+2]] -
                2 * i_km1_2 * v[cols[o+4]] - 4 * i_km1_2 * v[cols[o+5]]
            )
        # Δt_k row/column (the −a_start block only)
        Hv[cols[o+8]] +=
            ω * (
                h_k * vt_k - 12 * i_k_3 * v[cols[o+2]] + 12 * i_k_3 * v[cols[o+3]] -
                4 * i_k_2 * v[cols[o+5]] - 2 * i_k_2 * v[cols[o+6]]
            )
        # the shared u_k / du_k columns carry one entry from each block
        Hv[cols[o+1]] += ω * (-12 * i_km1_3 * vt_m1)
        Hv[cols[o+2]] += ω * (12 * i_km1_3 * vt_m1 - 12 * i_k_3 * vt_k)
        Hv[cols[o+3]] += ω * (12 * i_k_3 * vt_k)
        Hv[cols[o+4]] += ω * (-2 * i_km1_2 * vt_m1)
        Hv[cols[o+5]] += ω * (-4 * i_km1_2 * vt_m1 - 4 * i_k_2 * vt_k)
        Hv[cols[o+6]] += ω * (-2 * i_k_2 * vt_k)
    else
        # a_start / a_end over [u_a, u_b, du_a, du_b, Δt]
        u_a = Z[cols[o+1]]
        u_b = Z[cols[o+2]]
        du_a = Z[cols[o+3]]
        du_b = Z[cols[o+4]]
        i = 1.0 / Z[cols[o+5]]
        i2 = i * i
        i3 = i2 * i
        i4 = i3 * i
        vt = v[cols[o+5]]
        if kind == HSA_A_START
            h = 36 * i4 * (u_b - u_a) - 4 * i3 * (2 * du_a + du_b)
            c_ua = 12 * i3
            c_ub = -12 * i3
            c_dua = 4 * i2
            c_dub = 2 * i2
        else
            h = -36 * i4 * (u_b - u_a) + 4 * i3 * (du_a + 2 * du_b)
            c_ua = -12 * i3
            c_ub = 12 * i3
            c_dua = -2 * i2
            c_dub = -4 * i2
        end
        Hv[cols[o+5]] +=
            ω * (
                h * vt +
                c_ua * v[cols[o+1]] +
                c_ub * v[cols[o+2]] +
                c_dua * v[cols[o+3]] +
                c_dub * v[cols[o+4]]
            )
        Hv[cols[o+1]] += ω * (c_ua * vt)
        Hv[cols[o+2]] += ω * (c_ub * vt)
        Hv[cols[o+3]] += ω * (c_dua * vt)
        Hv[cols[o+4]] += ω * (c_dub * vt)
    end
    return nothing
end

function constraint_stencil_hvp!(
    Hv::AbstractVector{Float64},
    c::HermiteSmoothAccelerationConstraint,
    w::AbstractVector{Float64},
    v::AbstractVector{Float64},
    x::AbstractVector{Float64},
)
    t = c.table
    length(Hv) == t.n_cols || throw(
        DimensionMismatch(
            "Hv has length $(length(Hv)), expected the full decision width $(t.n_cols)",
        ),
    )
    kinds = c.kinds
    cols = t.cols
    col_ptr = t.col_ptr
    @inbounds for f in eachindex(kinds)
        ω = stencil_functional_weight(t, w, f)
        _hsa_functional_hessian_action!(Hv, v, ω, kinds[f], x, cols, col_ptr[f] - 1)
    end
    return nothing
end

supports_matrix_free_constraint_hvp(::HermiteSmoothAccelerationConstraint) = true

@testitem "#458 HSA exact inequality HVP: FD parity of the weighted residual Hessian" begin
    using NamedTrajectories

    using Piccolo:
        HermiteSmoothAccelerationConstraint,
        constraint_stencil_hvp!,
        supports_matrix_free_constraint_hvp
    using DirectTrajOpt.CommonInterface: evaluate!
    using ForwardDiff
    using LinearAlgebra
    using Random

    Random.seed!(0x45800A1)

    function fixture(N, n_drives; globals::Bool, segments = nothing)
        comps = (
            u = randn(n_drives, N),
            du = randn(n_drives, N) * 0.5,
            Δt = reshape(0.05 .+ 0.1 .* rand(N), 1, N),
        )
        traj =
            globals ?
            NamedTrajectory(
                comps;
                timestep = :Δt,
                controls = :u,
                global_data = [0.3, -0.8],
                global_components = (ϕ = 1:2,),
            ) : NamedTrajectory(comps; timestep = :Δt, controls = :u)
        c =
            isnothing(segments) ? HermiteSmoothAccelerationConstraint(traj; a_max = 7.0) :
            HermiteSmoothAccelerationConstraint(traj; a_max = 7.0, segments = segments)
        return traj, c
    end

    # The AD oracle: ∇²(wᵀg)(x) · v with g the residual, x the flat decision vector.
    function oracle(traj, c, w, v)
        n_data = traj.dim * traj.N
        n_vars = n_data + traj.global_dim
        gfunc =
            Z -> begin
                t = NamedTrajectory(
                    traj;
                    datavec = Z[1:n_data],
                    global_data = Z[(n_data+1):end],
                )
                g = zeros(eltype(Z), c.g_dim)
                evaluate!(g, c, t)
                return dot(w, g)
            end
        x = collect(vec(traj))
        return ForwardDiff.hessian(gfunc, x) * v, x, n_vars
    end

    for (label, fx) in (
        "no globals" => fixture(8, 2; globals = false),
        "with globals" => fixture(8, 2; globals = true),
        "segment subset" => fixture(10, 2; globals = true, segments = 3:7),
        "split runs, one drive" =>
            fixture(9, 1; globals = false, segments = [1, 2, 3, 6, 7]),
    )
        traj, c = fx
        @test supports_matrix_free_constraint_hvp(c) == true

        w = randn(MersenneTwister(0x458A), c.g_dim)
        n_vars = traj.dim * traj.N + traj.global_dim
        v = randn(MersenneTwister(0x458B), n_vars)
        Hv_ref, x, _ = oracle(traj, c, w, v)

        # The action matches the exact Hessian of the weighted residual.
        Hv = zeros(n_vars)
        constraint_stencil_hvp!(Hv, c, w, v, x)
        @test norm(Hv) > 0                    # not vacuous
        @test norm(Hv - Hv_ref) <= 1e-8 * max(1.0, norm(Hv_ref))

        # Linearity in the weights: w₁ + w₂ doubles through the contraction.
        Hv_sum = zeros(n_vars)
        constraint_stencil_hvp!(Hv_sum, c, 2.0 .* w, v, x)
        @test norm(Hv_sum - 2.0 .* Hv) <= 1e-8 * max(1.0, norm(Hv))

        # ACCUMULATES: a second call doubles the output.
        constraint_stencil_hvp!(Hv, c, w, v, x)
        @test norm(Hv - 2.0 .* Hv_ref) <= 1e-8 * max(1.0, norm(Hv_ref))

        # A balanced ± pair on EVERY functional (each functional is a `±` pair here)
        # contracts ω = 0 before any second-order work: exactly zero, everywhere.
        w_bal = fill(0.75, c.g_dim)
        Hv_bal = zeros(n_vars)
        constraint_stencil_hvp!(Hv_bal, c, w_bal, v, x)
        @test all(iszero, Hv_bal)

        # Independent of the Jacobian coefficients: a refresh at a DIFFERENT iterate
        # (rewriting `table.coeffs` for another x) must not move the HVP at this x.
        using Piccolo: refresh_constraint_coefficients!
        refresh_constraint_coefficients!(c, traj)          # coeffs at `x`
        Hv_a = zeros(n_vars)
        constraint_stencil_hvp!(Hv_a, c, w, v, x)
        far = x .+ 0.5 .* randn(MersenneTwister(0x458C), n_vars)
        far_traj = NamedTrajectory(
            traj;
            datavec = far[1:(traj.dim*traj.N)],
            global_data = far[(traj.dim*traj.N+1):end],
        )
        refresh_constraint_coefficients!(c, far_traj)      # coeffs at `far`
        Hv_b = zeros(n_vars)
        constraint_stencil_hvp!(Hv_b, c, w, v, x)
        @test Hv_a == Hv_b                 # bitwise: the HVP never reads `coeffs`
    end
end

@testitem "#458 HSA inequality HVP: allocation is invariant to knot count" begin
    using NamedTrajectories

    using Piccolo: HermiteSmoothAccelerationConstraint, constraint_stencil_hvp!
    using Random

    function measure(N)
        Random.seed!(0x45800A2)
        n_drives = 2
        traj = NamedTrajectory(
            (
                u = randn(n_drives, N),
                du = randn(n_drives, N),
                Δt = reshape(fill(0.1, N), 1, N),
            );
            timestep = :Δt,
            controls = :u,
        )
        c = HermiteSmoothAccelerationConstraint(traj; a_max = 10.0)
        n_vars = traj.dim * traj.N + traj.global_dim
        Hv = zeros(n_vars)
        w = randn(c.g_dim)
        v = randn(n_vars)
        x = collect(vec(traj))
        f = () -> constraint_stencil_hvp!(Hv, c, w, v, x)
        f()
        f()
        return minimum(@allocated(f()) for _ = 1:5)
    end

    measure(33)
    measure(129)
    a = measure(33)
    b = measure(129)
    # KNOT-FLATNESS (ADR-0009 decision 3 / ADR-0010 decision 6): invariance, plus a
    # loose absolute ceiling so the constant floor cannot creep either.
    @test b == a
    @test b <= 4096
end
