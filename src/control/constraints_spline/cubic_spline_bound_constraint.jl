"""
    CubicSplineBoundConstraint <: AbstractNonlinearConstraint

Constraint that ensures a cubic Hermite spline stays within bounds at all points
along the interpolation, not just at knot points.

For each segment between knots k and k+1, the cubic polynomial can have extrema
at interior points. This constraint finds those extrema and ensures they satisfy bounds.

# Constructor
```julia
CubicSplineBoundConstraint(
    traj::NamedTrajectory,
    control_name::Symbol,
    lower_bound::Real,
    upper_bound::Real;
    times::Union{Nothing, AbstractVector{Int}}=nothing,
    n_interior_points::Int=2  # Number of points to check between knots
)
```

# Arguments
- `traj`: NamedTrajectory with cubic spline controls (requires :u and :du components)
- `control_name`: Name of the control variable (e.g., :u)
- `lower_bound`: Lower bound for the spline
- `upper_bound`: Upper bound for the spline
- `times`: Which knot-point segments to constrain (default: all segments 1:N-1)
- `n_interior_points`: Number of interior points to check per segment (default: 2)

# Example
```julia
# Ensure first control component stays between -1 and 1 along entire spline
constraint = CubicSplineBoundConstraint(
    traj, :u, -1.0, 1.0;
    times=1:traj.N-1
)
```

# Representation

The residual and the assembled Jacobian are served from a `ConstraintStencilTable`
(ADR-0010): one table *functional* per unique scalar sample `s(τ)` the constraint forms,
its constant columns declared once at construction, and only the coefficients refreshed
per iterate. Each sample feeds a `±` row pair (`lower - s ≤ 0` and `s - upper ≤ 0`)
sharing that one scalar, which is the structure `#332`'s matrix-free kernels exploit.

The kernels read the trajectory's flat `datavec` at those precomputed columns through a
**function barrier** — never `traj[name]`, which is inferred as `Any` and boxes every
scalar read. That was 1887 µs and megabyte-scale allocation per Jacobian at `N`=801.

**Stencil width 1.** Every row reads a single interval —
`u_k, u_{k+1}, du_k, du_{k+1}, Δt_k` — so no column lies more than one knot from its
anchor `k`. The constructor checks that against the declared columns, because the
declared width is what the sharded backend maximises over and a wrong declaration is
otherwise silent until multi-rank. Contrast `HermiteSmoothAccelerationConstraint`, whose
continuity rows couple two adjacent intervals (`k-1, k, k+1`) at the same width.

# Allocation

**Not** zero-allocation, and zero allocation is not the target: the interface hands the
constraint a `NamedTrajectory` whose `datavec` is an abstractly-typed field, so one
dynamic access costs a constant floor of ~80 B that no kernel work removes. What is
guaranteed, and gated, is **knot-flatness** — per-call allocation for `evaluate!` and for
both Jacobian entry points is independent of the knot count. See ADR-0009 decision 3 and
ADR-0010 decision 6.

# Notes
- Requires both :u and :du components in trajectory (cubic Hermite spline DOFs)
- Checks spline values at uniformly spaced interior points
- For more accuracy, increase `n_interior_points` (but adds more constraints)
- `n_interior_points` is a fixed count per segment, so the functional enumeration is
  static at construction and the table never resizes per iterate
"""
struct CubicSplineBoundConstraint <: AbstractNonlinearConstraint
    control_name::Symbol
    derivative_name::Symbol
    timestep_name::Symbol
    lower_bound::Float64
    upper_bound::Float64
    times::Vector{Int}
    n_interior_points::Int
    n_drives::Int
    equality::Bool
    g_dim::Int
    dim::Int
    var_dim::Int

    # Functional-indexed stencil table: the declared structure (constant columns, signed
    # row map, per-row constant, stencil width) plus the per-iterate coefficients.
    # Application — sparse fill, cached assembly, row scatter — is generic over it.
    table::ConstraintStencilTable

    # Per-functional interior-point parameter τ ∈ (0,1). This constraint's analogue of
    # the smooth-acceleration constraint's `kinds` tag: there is exactly ONE closed-form
    # functional here (the spline value), parameterised by τ, so the per-functional datum
    # is a τ rather than a kind enum. The constraint's own math; the table never reads it.
    τs::Vector{Float64}

    # Jacobian sparsity structure (derived from the table, cached for the interface)
    jac_structure::Vector{Tuple{Int,Int}}

    # Hessian, declared zero. Exact only while Δt is not a decision variable: the sample
    # is linear in `u` and `du`, but `∂²s/∂Δt∂du_k = h10(τ) ≠ 0`. That is a PRE-EXISTING
    # formulation gap, out of scope for #331 (a representation change, not a formulation
    # change) and recorded in the `eval_hessian_of_lagrangian` docstring below. Sized at
    # the FULL decision-vector width from the table, so it stacks against a trajectory
    # carrying globals — the pre-port sizing omitted them.
    μ∂²g_full::SparseMatrixCSC{Float64,Int}
end

# Stencil width: one interval per row (`k, k+1`), i.e. no column further than one knot
# from the anchor. DECLARED, never acted on here — halo exchange belongs to the sharded
# driver, which takes `max` over dynamics and every routed constraint and exchanges once.
const CSB_STENCIL_WIDTH = 1

function CubicSplineBoundConstraint(
    traj::NamedTrajectory,
    control_name::Symbol = :u,
    lower_bound::Real = -Inf,
    upper_bound::Real = Inf;
    times::Union{Nothing,AbstractVector{Int}} = nothing,
    n_interior_points::Int = 2,
)
    @assert control_name ∈ traj.names "Control $control_name not found in trajectory"

    derivative_name = Symbol(:d, control_name)
    @assert derivative_name ∈ traj.names "Derivative $derivative_name not found - cubic spline requires both u and du"

    timestep_name = traj.timestep
    @assert haskey(traj.components, timestep_name) "Trajectory missing timestep component $timestep_name"

    # Default: all segments
    if isnothing(times)
        times = collect(1:(traj.N-1))
    else
        times = collect(times)
    end

    @assert !isempty(times) "times is empty — nothing to constrain"
    @assert all(k -> 1 <= k <= traj.N - 1, times) "times must lie in 1:$(traj.N-1)"
    @assert n_interior_points >= 1 "n_interior_points must be at least 1, got $n_interior_points"

    n_drives = traj.dims[control_name]

    # For each segment k→k+1, check n_interior_points in each drive
    # Total constraints: length(times) * n_interior_points * n_drives * 2 (upper and lower)
    g_dim = length(times) * n_interior_points * n_drives * 2

    # Variable dimension: u and du at two knots per segment
    var_dim = 2 * n_drives * 2  # (u_k, u_{k+1}, du_k, du_{k+1})

    dim = g_dim

    # Interior sample points, uniformly spaced in τ ∈ (0,1). A FIXED count per segment, so
    # the functional enumeration below is static: the table needs no per-iterate resizing.
    # Materialised once here; iterating the range and collecting it give bit-identical
    # Float64s, which is what lets the residual-parity witness assert `==`.
    τ_values = collect(range(0.0, 1.0, length = n_interior_points + 2)[2:(end-1)])

    # ── Structure declaration ─────────────────────────────────────────────────────── #
    #
    # THE FULL DECISION-VECTOR WIDTH, INCLUDING GLOBALS, carried as data on the table
    # rather than recomputed at each use site: sizing the Jacobian without the global
    # columns is what made the backend's row-stacking throw outright on any trajectory
    # carrying free phases.
    n_cols = traj.dim * traj.N + traj.global_dim
    n_knot_cols = traj.dim * traj.N
    u_comps = traj.components[control_name]
    du_comps = traj.components[derivative_name]
    dt_comps = traj.components[timestep_name]
    u_col = (k, d) -> slice(k, u_comps, traj.dim)[d]
    du_col = (k, d) -> slice(k, du_comps, traj.dim)[d]
    dt_col = k -> slice(k, dt_comps, traj.dim)[1]

    lo = Float64(lower_bound)
    hi = Float64(upper_bound)

    functional_cols = Vector{Vector{Int}}()
    τs = Float64[]
    row_map = Int[]
    row_offset = Float64[]

    # Row order is the pre-port one exactly: segment, then drive, then sample point, and
    # each sample's LOWER row before its UPPER row. Both rows of a sample are emitted
    # together, so a functional's rows are CONTIGUOUS by construction — the invariant
    # `ConstraintStencilTable` enforces.
    for k in times, d = 1:n_drives, τ in τ_values
        cols = [u_col(k, d), u_col(k + 1, d), du_col(k, d), du_col(k + 1, d), dt_col(k)]

        # AC6 of #331: the DECLARED stencil width, checked against the declared columns.
        # Every column must sit within `CSB_STENCIL_WIDTH` knots of the anchor `k`. The
        # sharded backend takes `max` over dynamics and every routed constraint and
        # exchanges one halo at that width, so an under-declaration corrupts nothing
        # locally and everything at multi-rank — it has to be caught here.
        for c in cols
            knot = ((c - 1) ÷ traj.dim) + 1
            @assert c <= n_knot_cols "column $c is not a knot column"
            @assert abs(knot - k) <= CSB_STENCIL_WIDTH "declared column $c sits $(abs(knot - k)) knots from anchor $k, exceeding the declared stencil width $CSB_STENCIL_WIDTH"
        end

        push!(functional_cols, cols)
        push!(τs, τ)
        f = length(functional_cols)
        # Lower: `lo - s ≤ 0`  →  sign −1, offset  lo
        push!(row_map, -f)
        push!(row_offset, lo)
        # Upper: `s - hi ≤ 0`  →  sign +1, offset −hi
        push!(row_map, f)
        push!(row_offset, -hi)
    end

    @assert length(row_map) == g_dim "internal: emitted $(length(row_map)) rows, expected $g_dim"

    table = ConstraintStencilTable(
        functional_cols,
        row_map,
        row_offset,
        n_cols;
        stencil_width = CSB_STENCIL_WIDTH,
    )

    μ∂²g_full = spzeros(n_cols, n_cols)

    return CubicSplineBoundConstraint(
        control_name,
        derivative_name,
        timestep_name,
        lo,
        hi,
        times,
        n_interior_points,
        n_drives,
        false,  # inequality constraint (bounds)
        g_dim,
        dim,
        var_dim,
        table,
        τs,
        stencil_structure(table),
        μ∂²g_full,
    )
end

jacobian_structure(c::CubicSplineBoundConstraint) = c.jac_structure

"""
    stencil_width(constraint::CubicSplineBoundConstraint) -> Int

The declared stencil width: the maximum number of knots either side of its anchor that
one constraint row reads. 1 here — every row samples a single interval `k → k+1`.

Declared, never acted on: halo exchange belongs to the sharded driver, which takes `max`
over dynamics and every routed constraint and exchanges once per callback.
"""
stencil_width(c::CubicSplineBoundConstraint) = c.table.stencil_width

# ----------------------------------------------------------------------------- #
# Barrier kernels
# ----------------------------------------------------------------------------- #
#
# `Z` arrives as `traj.datavec`, an abstractly-typed FIELD — reading scalars through it at
# the call site (or through `traj[name]`, inferred as `Any`) boxes every intermediate and
# was the bulk of the pre-port allocation. These functions are the FUNCTION BARRIER: one
# dynamic dispatch at entry, then a fully concrete body. Passing component *views* through
# a barrier is measurably not sufficient — it leaves a constant kilobyte-scale floor.

"""
    _csb_residual!(g, Z, τs, table)

Residual barrier: sample the spline at each functional's τ and scatter the value straight
into its (contiguous) `±` row pair. No per-functional buffer, so this is generic in
`eltype(g)` / `eltype(Z)` — forward-mode AD of this very function is the Jacobian oracle.
"""
function _csb_residual!(
    g::AbstractVector,
    Z::AbstractVector,
    τs::Vector{Float64},
    table::ConstraintStencilTable,
)
    cols = table.cols
    col_ptr = table.col_ptr
    @inbounds for f in eachindex(τs)
        o = col_ptr[f] - 1
        u_k = Z[cols[o+1]]
        u_k1 = Z[cols[o+2]]
        du_k = Z[cols[o+3]]
        du_k1 = Z[cols[o+4]]
        Δt = Z[cols[o+5]]
        val = evaluate_hermite_spline(τs[f], u_k, u_k1, du_k, du_k1, Δt)
        stencil_scatter_functional!(g, table, f, val)
    end
    return nothing
end

"""
    _csb_refresh_coefficients!(Z, τs, table)

Per-iterate coefficient refresh: overwrite `table.coeffs` with ∂s/∂column, in the declared
column order `(u_k, u_{k+1}, du_k, du_{k+1}, Δt_k)`. The `±` row signs and the per-row
bound offsets are the table's business, not this kernel's.
"""
function _csb_refresh_coefficients!(
    Z::AbstractVector,
    τs::Vector{Float64},
    table::ConstraintStencilTable,
)
    cols = table.cols
    col_ptr = table.col_ptr
    coeffs = table.coeffs
    @inbounds for f in eachindex(τs)
        o = col_ptr[f] - 1
        du_k = Z[cols[o+3]]
        du_k1 = Z[cols[o+4]]
        Δt = Z[cols[o+5]]
        ∂u_k, ∂u_k1, ∂du_k, ∂du_k1, ∂Δt = hermite_value_gradient(du_k, du_k1, Δt, τs[f])
        coeffs[o+1] = ∂u_k
        coeffs[o+2] = ∂u_k1
        coeffs[o+3] = ∂du_k
        coeffs[o+4] = ∂du_k1
        coeffs[o+5] = ∂Δt
    end
    stencil_touch!(table)
    return nothing
end

# ── The constraint-type gradient trait (#332) ─────────────────────────────────── #
#
# Two methods, and `supports_matrix_free_constraint_gradient` follows from them: this
# constraint declares a bounded stencil width (one interval per row), so the backend's
# inequality path serves its rows from `stencil_jvp!` / `stencil_vjp!` and gives it no
# block in the assembled matrix. `eval_jacobian` above is untouched and stays the
# assembled contract for Ipopt/MadNLP.

constraint_stencil_table(c::CubicSplineBoundConstraint) = c.table

refresh_constraint_coefficients!(c::CubicSplineBoundConstraint, traj::NamedTrajectory) =
    _csb_refresh_coefficients!(traj.datavec, c.τs, c.table)

# ----------------------------------------------------------------------------- #
# DirectTrajOpt interface — three-line wrappers over the barriers
# ----------------------------------------------------------------------------- #

"""
    evaluate!(values, constraint, traj)

Evaluate the interior-sample bounds as pure inequalities (≤ 0), two rows per sample:
`lower - s(τ) ≤ 0` and `s(τ) - upper ≤ 0`.
"""
function CommonInterface.evaluate!(
    values::AbstractVector,
    constraint::CubicSplineBoundConstraint,
    traj::NamedTrajectory,
)
    _csb_residual!(values, traj.datavec, constraint.τs, constraint.table)
    return nothing
end

"""
    jacobian!(jac, constraint, traj)

Compute the Jacobian's nonzero values, in `jacobian_structure` order.
"""
function jacobian!(
    jac::AbstractVector,
    constraint::CubicSplineBoundConstraint,
    traj::NamedTrajectory,
)
    _csb_refresh_coefficients!(traj.datavec, constraint.τs, constraint.table)
    stencil_fill_values!(jac, constraint.table)
    return nothing
end

"""
    jacobian!(constraint, traj)

Refresh the cached assembled Jacobian in place (the pre-port two-argument form; the
matrix is then read from `get_full_jacobian` / `eval_jacobian`).
"""
function jacobian!(constraint::CubicSplineBoundConstraint, traj::NamedTrajectory)
    _csb_refresh_coefficients!(traj.datavec, constraint.τs, constraint.table)
    stencil_assemble!(constraint.table)
    return nothing
end

function DirectTrajOpt.Constraints.get_full_jacobian(
    constraint::CubicSplineBoundConstraint,
    traj::NamedTrajectory,
)
    _csb_refresh_coefficients!(traj.datavec, constraint.τs, constraint.table)
    return stencil_assemble!(constraint.table)
end

"""
    eval_jacobian(constraint, traj) -> SparseMatrixCSC

The ASSEMBLED Jacobian entry point — unchanged in contract (this is what Ipopt/MadNLP and
the backend's row-stacking consume). Served from the stencil table's cached sparse matrix:
the pattern is built once at construction and only `nzval` is refreshed, so the returned
matrix is the same object on every call and the call allocates nothing that grows with the
knot count.

`size(J) == (g_dim, traj.dim * traj.N + traj.global_dim)` — the FULL decision-vector
width. Sizing it without the global columns made the row-stacking throw on any trajectory
carrying free phases.
"""
function eval_jacobian(constraint::CubicSplineBoundConstraint, traj::NamedTrajectory)
    _csb_refresh_coefficients!(traj.datavec, constraint.τs, constraint.table)
    return stencil_assemble!(constraint.table)
end

function hessian_of_lagrangian!(
    constraint::CubicSplineBoundConstraint,
    traj::NamedTrajectory,
    μ::AbstractVector,
)
    return nothing
end

function DirectTrajOpt.Constraints.get_full_hessian(
    constraint::CubicSplineBoundConstraint,
    traj::NamedTrajectory,
    μ::AbstractVector,
)
    return constraint.μ∂²g_full
end

"""
    eval_hessian_of_lagrangian(constraint, traj, μ) -> SparseMatrixCSC

Declared zero, at the full decision-vector width.

The sample `s(τ) = h00·u_k + h10·Δt·du_k + h01·u_{k+1} + h11·Δt·du_{k+1}` is linear in
`u` and `du`, so with a FIXED timestep the Hessian really is zero. With `Δt` a decision
variable it is not: `∂²s/∂Δt∂du_k = h10(τ)`. That gap is inherited from the pre-port
implementation and deliberately left alone — #331 is a representation change, not a
formulation change — and is the same class of issue as the `CubicSplineExtremaConstraint`
discontinuity recorded in ADR-0010's consequences.
"""
function CommonInterface.eval_hessian_of_lagrangian(
    constraint::CubicSplineBoundConstraint,
    traj::NamedTrajectory,
    μ::AbstractVector,
)
    return constraint.μ∂²g_full
end

# ============================================================================= #
#   Reference (pre-port) residual and Jacobian — the semantics witnesses         #
# ============================================================================= #
#
# Kept reachable so the port can be checked against the implementation it replaces, on
# randomised trajectories (#331 AC1/AC5). Both read the trajectory through the
# named-component accessor exactly as the pre-port code did — which is why they are
# test-only and on no hot path.

"""
    _csb_reference_evaluate!(values, constraint, traj)

The pre-port residual, verbatim in formulation and row order. Test-only witness for
"residual semantics unchanged".
"""
function _csb_reference_evaluate!(
    values::AbstractVector,
    constraint::CubicSplineBoundConstraint,
    traj::NamedTrajectory,
)
    u = traj[constraint.control_name]
    du = traj[constraint.derivative_name]
    Δt = traj[constraint.timestep_name]

    constraint_idx = 1

    for k in constraint.times
        k_next = k + 1
        Δt_k = Δt[k]

        τ_values = range(0.0, 1.0, length = constraint.n_interior_points + 2)[2:(end-1)]

        for drive_idx = 1:constraint.n_drives
            u_k = u[drive_idx, k]
            u_k1 = u[drive_idx, k_next]
            du_k = du[drive_idx, k]
            du_k1 = du[drive_idx, k_next]

            for τ in τ_values
                s = evaluate_hermite_spline(τ, u_k, u_k1, du_k, du_k1, Δt_k)

                # Lower bound: lower_bound ≤ s  →  lower_bound - s ≤ 0
                values[constraint_idx] = constraint.lower_bound - s
                constraint_idx += 1

                # Upper bound: s ≤ upper_bound  →  s - upper_bound ≤ 0
                values[constraint_idx] = s - constraint.upper_bound
                constraint_idx += 1
            end
        end
    end

    return nothing
end

"""
    _csb_reference_jacobian(constraint, traj) -> SparseMatrixCSC

The pre-port assembled Jacobian, verbatim: built fresh, at the pre-port width
`traj.dim * traj.N` (which is exactly the sizing bug — the global columns are missing).
Test-only witness for "the assembled entry point is numerically identical".
"""
function _csb_reference_jacobian(
    constraint::CubicSplineBoundConstraint,
    traj::NamedTrajectory,
)
    du = traj[constraint.derivative_name]
    Δt = traj[constraint.timestep_name]

    I_rows = Int[]
    J_cols = Int[]
    V = Float64[]
    constraint_idx = 1

    for k in constraint.times
        k_next = k + 1
        Δt_k = Δt[k]

        τ_values = range(0.0, 1.0, length = constraint.n_interior_points + 2)[2:(end-1)]

        u_k_slice = slice(k, traj.components[constraint.control_name], traj.dim)
        u_k1_slice = slice(k_next, traj.components[constraint.control_name], traj.dim)
        du_k_slice = slice(k, traj.components[constraint.derivative_name], traj.dim)
        du_k1_slice = slice(k_next, traj.components[constraint.derivative_name], traj.dim)
        Δt_k_idx = slice(k, traj.components[constraint.timestep_name], traj.dim)[1]

        for drive_idx = 1:constraint.n_drives
            var_idxs = [
                u_k_slice[drive_idx],
                u_k1_slice[drive_idx],
                du_k_slice[drive_idx],
                du_k1_slice[drive_idx],
                Δt_k_idx,
            ]
            du_k_val = du[drive_idx, k]
            du_k1_val = du[drive_idx, k_next]

            for τ in τ_values
                # The body of the retired `eval_cubic_hermite_jacobian`, written out here
                # ON PURPOSE: this is the pre-port witness, and calling the shared
                # `hermite_value_gradient` (which the port itself calls) would stop it
                # witnessing the coefficient transcription. Not an AC7 duplication —
                # nothing on a hot path reaches it.
                h00, h10, h01, h11 = hermite_basis_functions(τ)
                ∂s = [h00, h01, Δt_k * h10, Δt_k * h11, h10 * du_k_val + h11 * du_k1_val]

                for (j, v) in zip(var_idxs, ∂s)          # lower: ∂g/∂x = -∂s/∂x
                    push!(I_rows, constraint_idx)
                    push!(J_cols, j)
                    push!(V, -v)
                end
                constraint_idx += 1

                for (j, v) in zip(var_idxs, ∂s)          # upper: ∂g/∂x = +∂s/∂x
                    push!(I_rows, constraint_idx)
                    push!(J_cols, j)
                    push!(V, v)
                end
                constraint_idx += 1
            end
        end
    end

    return sparse(I_rows, J_cols, V, constraint.g_dim, traj.dim * traj.N)
end

# ============================================================================= #
# Tests
# ============================================================================= #

@testitem "CubicSplineBoundConstraint - basic construction" begin
    using NamedTrajectories
    using Piccolo

    N = 10
    n_drives = 2

    # Create trajectory with cubic spline components
    traj = NamedTrajectory(
        (
            x = randn(3, N),
            u = randn(n_drives, N),
            du = randn(n_drives, N),
            Δt = fill(0.1, N),
        );
        timestep = :Δt,
        controls = :u,
    )

    # Create constraint
    constraint = CubicSplineBoundConstraint(traj, :u, -1.0, 1.0)

    @test constraint isa CubicSplineBoundConstraint
    @test constraint.control_name == :u
    @test constraint.derivative_name == :du
    @test constraint.lower_bound == -1.0
    @test constraint.upper_bound == 1.0
    @test constraint.n_drives == n_drives
    @test constraint.times == collect(1:(N-1))
    @test constraint.n_interior_points == 2  # default

    # Check constraint dimensions
    # Each segment has n_interior_points * n_drives * 2 (lower+upper) constraints
    expected_dim = (N-1) * 2 * n_drives * 2
    @test constraint.g_dim == expected_dim
    @test constraint.dim == expected_dim
end

@testitem "CubicSplineBoundConstraint - custom times and interior points" begin
    using NamedTrajectories
    using Piccolo

    N = 20
    traj = NamedTrajectory(
        (u = randn(1, N), du = randn(1, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
    )

    # Custom settings
    times = 5:15
    n_interior = 3
    constraint = CubicSplineBoundConstraint(
        traj,
        :u,
        -2.0,
        2.0;
        times = times,
        n_interior_points = n_interior,
    )

    @test constraint.times == collect(times)
    @test constraint.n_interior_points == n_interior
    @test constraint.g_dim == length(times) * n_interior * 1 * 2
end

@testitem "CubicSplineBoundConstraint - evaluation" begin
    using NamedTrajectories
    using Piccolo
    using LinearAlgebra

    N = 5
    n_drives = 1

    # Create simple trajectory with known values
    u = [0.0 0.5 1.0 0.5 0.0]  # 1×5 matrix
    du = [1.0 0.0 -1.0 -1.0 0.0]  # derivatives
    Δt_vec = fill(0.2, N)

    traj = NamedTrajectory(
        (u = u, du = du, Δt = reshape(Δt_vec, 1, N));
        timestep = :Δt,
        controls = :u,
    )

    # Bounds that should be satisfied
    constraint = CubicSplineBoundConstraint(traj, :u, -0.5, 1.5; n_interior_points = 3)

    values = zeros(constraint.g_dim)
    evaluate!(values, constraint, traj)

    # All constraints should be satisfied (≤ 0 for inequality constraints)
    # Lower bound: lower - s ≤ 0  →  s ≥ lower
    # Upper bound: s - upper ≤ 0  →  s ≤ upper
    @test length(values) == constraint.g_dim

    # With bounds [-0.5, 1.5] and max control at 1.0, all should be feasible
    # Check a few constraints are within reasonable range
    @test all(isfinite.(values))
end

@testitem "CubicSplineBoundConstraint - detects violations" begin
    using NamedTrajectories
    using Piccolo

    N = 5

    # Create trajectory that will violate bounds at interior points
    # Use knot points that are within bounds but derivatives create overshoot
    u = [0.0 0.1 0.0 0.1 0.0]  # Oscillating but within [-0.2, 0.2]
    du = [2.0 -2.0 2.0 -2.0 0.0]  # Large opposing derivatives cause overshoot
    Δt_vec = fill(0.2, N)

    traj = NamedTrajectory(
        (u = u, du = du, Δt = reshape(Δt_vec, 1, N));
        timestep = :Δt,
        controls = :u,
    )

    # Very tight bounds - the cubic interpolation will overshoot these
    constraint = CubicSplineBoundConstraint(traj, :u, -0.15, 0.15; n_interior_points = 5)

    values = zeros(constraint.g_dim)
    evaluate!(values, constraint, traj)

    # Should have violations (positive constraint values)
    # The Hermite spline with large opposing derivatives will exceed bounds
    @test any(values .> 1e-6)  # Use tolerance to avoid numerical issues
end

@testitem "CubicSplineBoundConstraint - jacobian correctness" begin
    using NamedTrajectories
    using Piccolo
    using FiniteDiff
    using SparseArrays
    using LinearAlgebra

    N = 6
    n_drives = 2

    traj = NamedTrajectory(
        (u = 0.1 * randn(n_drives, N), du = 0.1 * randn(n_drives, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
    )

    constraint = CubicSplineBoundConstraint(
        traj,
        :u,
        -1.0,
        1.0;
        times = 2:(N-2),  # Test subset
        n_interior_points = 2,
    )

    # Function to evaluate constraints
    g_func = Z -> begin
        test_traj = NamedTrajectory(traj; datavec = Z)
        vals = zeros(constraint.g_dim)
        evaluate!(vals, constraint, test_traj)
        return vals
    end

    # Analytical Jacobian
    ∂g_analytical = Matrix(eval_jacobian(constraint, traj))

    # Finite difference Jacobian
    Z = copy(traj.datavec)
    ∂g_fd = FiniteDiff.finite_difference_jacobian(g_func, Z)

    # Check they match
    rel_error = norm(∂g_analytical - ∂g_fd) / (norm(∂g_fd) + 1e-10)
    @test rel_error < 1e-5
end

@testitem "CubicSplineBoundConstraint - multiple drives" begin
    using NamedTrajectories
    using Piccolo
    using SparseArrays

    N = 8
    n_drives = 3

    traj = NamedTrajectory(
        (u = randn(n_drives, N), du = randn(n_drives, N), Δt = fill(0.15, N));
        timestep = :Δt,
        controls = :u,
    )

    constraint = CubicSplineBoundConstraint(traj, :u, -2.0, 2.0; n_interior_points = 3)

    # Check dimensions scale correctly with number of drives
    expected_constraints = (N-1) * 3 * n_drives * 2
    @test constraint.g_dim == expected_constraints

    # Evaluate and check structure
    values = zeros(constraint.g_dim)
    evaluate!(values, constraint, traj)

    @test length(values) == expected_constraints
    @test all(isfinite.(values))

    # Check Jacobian sparsity
    J = eval_jacobian(constraint, traj)
    @test nnz(J) > 0
    @test size(J) == (constraint.g_dim, traj.dim * N + traj.global_dim)
end

@testitem "CubicSplineBoundConstraint - hermite evaluation accuracy" begin
    # AC7 of #331: this testitem used to define its own copies of the four Hermite basis
    # polynomials and of the value gradient, so it tested those copies and nothing in the
    # library. It now exercises the CONSOLIDATED primitives — which is what the constraint
    # itself calls.
    using Piccolo: evaluate_hermite_spline, hermite_value_gradient
    using ForwardDiff

    u_k = 0.0
    u_k1 = 1.0
    du_k = 0.5
    du_k1 = -0.5
    Δt = 0.1

    # At τ=0, should equal u_k; at τ=1, u_k1
    @test evaluate_hermite_spline(0.0, u_k, u_k1, du_k, du_k1, Δt) ≈ u_k
    @test evaluate_hermite_spline(1.0, u_k, u_k1, du_k, du_k1, Δt) ≈ u_k1

    # Value gradient at the segment boundaries, in the constraint's declared column order
    # (u_k, u_{k+1}, du_k, du_{k+1}, Δt).
    ∂s0 = hermite_value_gradient(du_k, du_k1, Δt, 0.0)
    @test ∂s0[1] ≈ 1.0  # ∂s/∂u_k at τ=0
    @test ∂s0[2] ≈ 0.0  # ∂s/∂u_k1 at τ=0

    ∂s1 = hermite_value_gradient(du_k, du_k1, Δt, 1.0)
    @test ∂s1[1] ≈ 0.0  # ∂s/∂u_k at τ=1
    @test ∂s1[2] ≈ 1.0  # ∂s/∂u_k1 at τ=1

    # …and the gradient IS the gradient of the value, including the Δt partial that the
    # retired `eval_cubic_hermite_jacobian` copy never checked.
    for τ in (0.0, 0.25, 0.5, 0.75, 1.0)
        ad = ForwardDiff.gradient(
            v -> evaluate_hermite_spline(τ, v[1], v[2], v[3], v[4], v[5]),
            [u_k, u_k1, du_k, du_k1, Δt],
        )
        @test all(isapprox.(collect(hermite_value_gradient(du_k, du_k1, Δt, τ)), ad))
    end
end

@testitem "#331 AC1/AC2/AC5: ported residual matches the reference; Jacobian matches forward-mode AD" begin
    using NamedTrajectories

    using Piccolo: CubicSplineBoundConstraint, eval_jacobian
    using DirectTrajOpt.CommonInterface: evaluate!, jacobian_structure
    using DirectTrajOpt.Constraints: get_full_jacobian, jacobian!
    using Piccolo.Control.QuantumConstraints.SplineConstraints:
        _csb_reference_evaluate!, _csb_reference_jacobian
    using ForwardDiff
    using SparseArrays
    using LinearAlgebra
    using Random

    Random.seed!(0x331AC12)

    # Trajectory factory, with and without global variables. The globals case is the one
    # that matters for the DERIVATIVE oracle — the pre-port implementation sized its
    # Jacobian without the global columns, so parity against it proves nothing there.
    # Forward-mode AD OF THE PORTED RESIDUAL is the oracle: everything here is closed
    # form, so a tolerance loose enough to absorb finite-difference truncation error is
    # loose enough to hide a wrong stencil coefficient.
    function fixture(N, n_drives; globals::Bool, n_interior = 2, times = nothing)
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
        c = CubicSplineBoundConstraint(
            traj,
            :u,
            -0.9,
            1.1;
            n_interior_points = n_interior,
            times = times,
        )
        return traj, c
    end

    # ── AC1: residual semantics unchanged, to machine precision, on randomised
    # trajectories. The reference is the pre-port implementation, kept reachable.
    for trial = 1:5
        N = rand(4:9)
        nd = rand(1:3)
        ni = rand(1:4)                        # includes the non-default interior counts
        traj, c = fixture(N, nd; globals = false, n_interior = ni)
        g_new = zeros(c.g_dim)
        g_ref = zeros(c.g_dim)
        evaluate!(g_new, c, traj)
        _csb_reference_evaluate!(g_ref, c, traj)
        @test g_new == g_ref        # bit-identical: same formulation, same row order
    end

    # Same, with globals present (neither the formulation nor the reference reads them),
    # and on a subset of segments.
    for fx in (
        fixture(7, 2; globals = true, n_interior = 3),
        fixture(9, 2; globals = false, n_interior = 5, times = 3:6),
    )
        traj, c = fx
        g_new = zeros(c.g_dim)
        g_ref = zeros(c.g_dim)
        evaluate!(g_new, c, traj)
        _csb_reference_evaluate!(g_ref, c, traj)
        @test g_new == g_ref
    end

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
        "non-default interior count" => fixture(8, 2; globals = true, n_interior = 5),
        "segment subset" => fixture(10, 3; globals = true, times = 4:8),
        "single interior point" => fixture(6, 1; globals = false, n_interior = 1),
    )
        traj, c = fx
        J = eval_jacobian(c, traj)
        J_ad = ad_jacobian(traj, c)
        @test size(J) == size(J_ad)      # includes the global columns
        @test size(J, 2) == traj.dim * traj.N + traj.global_dim
        @test maximum(abs, Matrix(J) - J_ad) < 1e-12 * max(1.0, maximum(abs, J_ad))

        # The flat-value entry point agrees with the assembled one, entry for entry.
        vals = zeros(length(jacobian_structure(c)))
        jacobian!(vals, c, traj)
        for (e, (r, col)) in enumerate(jacobian_structure(c))
            @test vals[e] == J[r, col]
        end
    end

    # ── AC5: the assembled entry point is still there and NUMERICALLY IDENTICAL to the
    # pre-port one, which is rebuilt from scratch here the way the pre-port code did.
    for (label, fx) in (
        "no globals" => fixture(9, 2; globals = false, n_interior = 3),
        "with globals" => fixture(9, 2; globals = true, n_interior = 3),
        "segment subset" => fixture(11, 2; globals = false, n_interior = 2, times = 2:9),
    )
        traj, c = fx
        J = eval_jacobian(c, traj)
        J_ref = _csb_reference_jacobian(c, traj)
        n_data = traj.dim * traj.N
        @test size(J_ref) == (c.g_dim, n_data)                 # the pre-port width
        @test Matrix(J[:, 1:n_data]) == Matrix(J_ref)          # bit-identical values
        # …and the columns the pre-port sizing dropped are structurally absent, so the
        # only difference is that the block now stacks against a globals trajectory.
        @test iszero(J[:, (n_data+1):end])

        # A CACHED reference refilled in place: same object every call, and the two
        # assembled entry points return that same object.
        @test eval_jacobian(c, traj) === J
        @test get_full_jacobian(c, traj) === J
        @test J === c.table.J
    end
end

@testitem "#331 AC3: residual and Jacobian allocation are invariant to knot count" begin
    using NamedTrajectories

    using Piccolo: CubicSplineBoundConstraint, eval_jacobian
    using DirectTrajOpt.CommonInterface: evaluate!
    using DirectTrajOpt.Constraints: jacobian!
    using Random

    # KNOT-FLATNESS, not zero allocation. The interface hands the constraint a trajectory
    # whose `datavec` is an abstractly-typed field, so one dynamic access costs a constant
    # floor no kernel work removes (ADR-0010 decision 6). What is asserted is that the
    # floor is CONSTANT: an implementation that still reads the trajectory through the
    # `Any`-inferred accessor, allocates a five-element gradient vector per sample point,
    # or rebuilds a whole sparse matrix per call, grows linearly in N and fails the `==`.
    #
    # Shape copied from #330's gate: local-scoped closures only (a callback at global
    # scope captures `Any` and inflates the count), and the MINIMUM over repeats because
    # `@allocated` sums across every thread in the process.
    function measure(N)
        Random.seed!(0x331A11C)
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
        c = CubicSplineBoundConstraint(traj, :u, -1.0, 1.0; n_interior_points = 3)
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
    # `@allocated` site in a process reports a different figure from every subsequent one,
    # and with `==` as the assertion whichever size went first would otherwise carry that
    # one-time step as phantom N-growth.
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

@testitem "#331 AC6: declared stencil width of one, and the declared structure it implies" begin
    using NamedTrajectories

    using Piccolo:
        CubicSplineBoundConstraint,
        eval_jacobian,
        hermite_value_gradient,
        stencil_coeff_range,
        stencil_functional_rows,
        stencil_width
    using Random

    Random.seed!(0x331AC06)

    # A non-default interior-point count, because the functional count scales with it.
    N, n_drives, n_interior = 8, 2, 4
    times = 2:6
    traj = NamedTrajectory(
        (
            u = randn(n_drives, N),
            du = 0.5 * randn(n_drives, N),
            Δt = reshape(0.05 .+ 0.1 .* rand(N), 1, N),
        );
        timestep = :Δt,
        controls = :u,
        global_data = [0.7],
        global_components = (ϕ = 1:1,),
    )
    c = CubicSplineBoundConstraint(
        traj,
        :u,
        -0.8,
        1.3;
        times = times,
        n_interior_points = n_interior,
    )

    # ── The declaration itself. ONE, because every row samples a single interval — unlike
    # the smooth-acceleration constraint's continuity rows, which couple two adjacent
    # intervals. The sharded backend takes `max` over dynamics and every routed
    # constraint, so an under-declaration is silent locally and wrong at multi-rank; this
    # is the only place it can be caught cheaply.
    @test c.table.stencil_width == 1
    @test stencil_width(c) == 1

    # ── …and the declaration is CONSISTENT with the declared columns: every column of
    # every functional sits within `stencil_width` knots of that functional's anchor. This
    # is the assertion that a wrong declaration would fail, as opposed to merely reading
    # the field back.
    t = c.table
    for (f, k) in
        enumerate(Iterators.flatten(fill(k, n_drives * n_interior) for k in times))
        knots = [((col - 1) ÷ traj.dim) + 1 for col in t.cols[stencil_coeff_range(t, f)]]
        @test minimum(knots) == k && maximum(knots) == k + 1
        @test maximum(abs.(knots .- k)) <= stencil_width(c)
    end

    # ── Functional enumeration is static and scales with the interior-point count; each
    # functional is read by exactly one `±` row pair sharing its scalar.
    n_functionals = length(times) * n_drives * n_interior
    @test t.n_functionals == n_functionals
    @test t.n_rows == c.g_dim == 2 * n_functionals
    @test t.n_cols == traj.dim * traj.N + traj.global_dim
    for f = 1:n_functionals
        rows = stencil_functional_rows(t, f)
        # `2 * f`, never `2f`: the formatter tightens `(2f - 1)` to `2f-1`, which re-lexes
        # as the Float32 literal 0.2f0 and silently turns this into a nonsense range.
        @test rows == (2*f-1):(2*f)
        @test t.row_sign[rows[1]] == -1.0 && t.row_offset[rows[1]] == -0.8   # lower - s
        @test t.row_sign[rows[2]] == +1.0 && t.row_offset[rows[2]] == -1.3   # s - upper
    end

    # ── The COEFFICIENTS are this constraint's own math, so check them directly against
    # the shared Hermite primitive rather than only through the assembled matrix.
    eval_jacobian(c, traj)                   # refreshes `t.coeffs`
    u_c, du_c, dt_c = traj.components[:u], traj.components[:du], traj.components[:Δt]
    # A `Ref`, not `f = 0` + `f += 1`: at top level inside a `for`, `f += 1` binds a NEW
    # local and throws an UndefVarError.
    fref = Ref(0)
    for k in times, d = 1:n_drives, i = 1:n_interior
        fref[] += 1
        f = fref[]
        τ = collect(range(0.0, 1.0, length = n_interior + 2)[2:(end-1)])[i]
        expected = hermite_value_gradient(
            traj.datavec[(k-1)*traj.dim+du_c[d]],
            traj.datavec[k*traj.dim+du_c[d]],
            traj.datavec[(k-1)*traj.dim+dt_c[1]],
            τ,
        )
        @test t.coeffs[stencil_coeff_range(t, f)] == collect(expected)
        # …and the columns they sit against are the declared ones, in order.
        @test t.cols[stencil_coeff_range(t, f)] == [
            (k - 1) * traj.dim + u_c[d],
            k * traj.dim + u_c[d],
            (k - 1) * traj.dim + du_c[d],
            k * traj.dim + du_c[d],
            (k - 1) * traj.dim + dt_c[1],
        ]
    end
    @test fref[] == n_functionals
end

@testitem "#458 CSB exact inequality HVP: FD parity of the weighted residual Hessian" begin
    using NamedTrajectories

    using Piccolo:
        CubicSplineBoundConstraint,
        constraint_stencil_hvp!,
        supports_matrix_free_constraint_hvp
    using DirectTrajOpt.CommonInterface: evaluate!
    using ForwardDiff
    using LinearAlgebra
    using Random

    Random.seed!(0x45800B1)

    function fixture(N, n_drives; globals::Bool, n_interior = 2, times = nothing)
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
        c = CubicSplineBoundConstraint(
            traj,
            :u,
            -0.9,
            1.1;
            n_interior_points = n_interior,
            times = times,
        )
        return traj, c
    end

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
        "non-default interior count" => fixture(8, 2; globals = true, n_interior = 4),
        "segment subset" => fixture(10, 3; globals = true, times = 4:8),
    )
        traj, c = fx
        @test supports_matrix_free_constraint_hvp(c) == true

        w = randn(MersenneTwister(0x458D), c.g_dim)
        n_vars = traj.dim * traj.N + traj.global_dim
        v = randn(MersenneTwister(0x458E), n_vars)
        Hv_ref, x, _ = oracle(traj, c, w, v)

        # The action matches the exact Hessian of the weighted residual — INCLUDING the
        # Δt/du cross terms the assembled entry point still declares zero (the documented
        # gap; the matrix-free slot is where it closes).
        Hv = zeros(n_vars)
        constraint_stencil_hvp!(Hv, c, w, v, x)
        @test norm(Hv) > 0
        @test norm(Hv - Hv_ref) <= 1e-9 * max(1.0, norm(Hv_ref))

        # SUPPORT: the only columns the action can touch are the Δt and du columns of
        # each functional's segment — assert a v supported purely on the u columns gets
        # exactly zero everywhere (bilinearity: s has no u² term).
        v_u_only = zeros(n_vars)
        u_comps = traj.components[:u]
        for k = 1:traj.N, d = 1:size(traj[:u], 1)
            v_u_only[(k-1)*traj.dim+u_comps[d]] = 1.0
        end
        Hv_u = zeros(n_vars)
        constraint_stencil_hvp!(Hv_u, c, w, v_u_only, x)
        @test all(iszero, Hv_u)

        # ACCUMULATES + balanced ± pairs cancel exactly (every functional is a `±` pair).
        constraint_stencil_hvp!(Hv, c, w, v, x)
        @test norm(Hv - 2.0 .* Hv_ref) <= 1e-9 * max(1.0, norm(Hv_ref))
        w_bal = fill(0.4, c.g_dim)
        Hv_bal = zeros(n_vars)
        constraint_stencil_hvp!(Hv_bal, c, w_bal, v, x)
        @test all(iszero, Hv_bal)
    end
end

@testitem "#458 CSB inequality HVP: allocation is invariant to knot count" begin
    using NamedTrajectories
    using Piccolo: CubicSplineBoundConstraint, constraint_stencil_hvp!
    using Random

    function measure(N)
        Random.seed!(0x45800B2)
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
        c = CubicSplineBoundConstraint(traj, :u, -1.0, 1.0; n_interior_points = 3)
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
    # KNOT-FLATNESS, invariance + loose ceiling (ADR-0009 decision 3 / ADR-0010 decision 6).
    @test b == a
    @test b <= 4096
end

# The matrix-free backend path closes that gap for its own slot: the sample
# `s(τ) = h00·u_k + h10·Δt·du_k + h01·u_{k+1} + h11·Δt·du_{k+1}` is bilinear in
# `(Δt, du)`, so its only nonzero second derivatives are the CONSTANT cross terms
# `∂²s/∂Δt∂du_k = h10(τ)` and `∂²s/∂Δt∂du_{k+1} = h11(τ)`. The action below supplies
# exactly those; on a trajectory whose Δt is not a decision column the terms vanish
# structurally (no Δt column exists to index).

function constraint_stencil_hvp!(
    Hv::AbstractVector{Float64},
    c::CubicSplineBoundConstraint,
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
    τs = c.τs
    cols = t.cols
    col_ptr = t.col_ptr
    @inbounds for f in eachindex(τs)
        ω = stencil_functional_weight(t, w, f)
        o = col_ptr[f] - 1
        h10, h11 = hermite_basis_functions(τs[f])[2], hermite_basis_functions(τs[f])[4]
        vt = v[cols[o+5]]                       # the Δt column
        v_du_k = v[cols[o+3]]
        v_du_k1 = v[cols[o+4]]
        # symmetric cross terms: (Δt, du_k) with h10(τ), (Δt, du_{k+1}) with h11(τ)
        Hv[cols[o+3]] += ω * h10 * vt
        Hv[cols[o+4]] += ω * h11 * vt
        Hv[cols[o+5]] += ω * (h10 * v_du_k + h11 * v_du_k1)
    end
    return nothing
end

supports_matrix_free_constraint_hvp(::CubicSplineBoundConstraint) = true
