# ============================================================================= #
#   ConstraintStencilTable — the functional-indexed intermediate representation   #
# ============================================================================= #
#
# ADR-0010 decision 4. The table's ROWS ARE FUNCTIONALS, not constraint rows: each
# unique scalar quantity the constraint forms is one table row, its CONSTANT columns
# are stored separately from its PER-ITERATE coefficients, and a SIGNED row-to-
# functional map says which constraint rows read which functional with which sign.
#
# Why not store constraint rows: every double-sided inequality emits a `±` row pair
# sharing one underlying scalar, so one directional-derivative evaluation serves two
# rows. A per-row table (or a sparse matrix) stores both rows' stencils
# independently and forfeits exactly that pairing — which is the entire reason the
# matrix-free kernels (#332) beat the SpMV. See CONTEXT.md, "Paired constraint rows".
#
# The split of labour: THE CONSTRAINT OWNS ITS MATH, THE TABLE OWNS ITS APPLICATION.
# A constraint implements exactly three things —
#
#   (a) a one-time STRUCTURE DECLARATION: which columns each functional touches,
#       the signed row-to-functional map, the per-row constant, and the stencil
#       width. That is this constructor's argument list.
#   (b) a per-iterate COEFFICIENT REFRESH: overwrite `table.coeffs` in place.
#   (c) a RESIDUAL.
#
# Application — sparse fill, sparsity construction, and (in #332) JVP/VJP — is
# generic over the table and lives here.
#
# ALLOCATION. Zero allocation is not the target and is not attainable: the
# `CommonInterface` hands a constraint a `NamedTrajectory` whose `datavec` is an
# abstractly-typed field, so one dynamic access costs a constant floor of ~80 B that
# no kernel work removes, and the backend's own per-callback trajectory wrapper costs
# ~3.3 kB — roughly forty times a whole optimised constraint kernel. The criterion is
# KNOT-FLATNESS: per-call allocation must not grow with the knot count. See ADR-0009
# decision 3 and ADR-0010 decision 6.

"""
    ConstraintStencilTable

Functional-indexed stencil table for an algebraic constraint's Jacobian.

Each **functional** is one unique scalar quantity the constraint forms (e.g. the
acceleration at one knot for one drive). Constraint **rows** are then signed,
offset copies of functionals: `g[r] = row_sign[r] * F[row_functional[r]] +
row_offset[r]`.

# Fields (declared)
- `n_functionals`, `n_rows`: table and constraint dimensions.
- `n_cols`: the FULL decision-vector width, `traj.dim * traj.N + traj.global_dim`.
  Carried as data, never recomputed at a use site — two constraints previously sized
  their Jacobians without the global columns, which made the backend's row-stacking
  throw on any trajectory carrying free phases.
- `stencil_width`: the maximum number of knots either side of its anchor that one row
  reads. A *declared* property; the constraint never exchanges halos itself (the
  sharded driver takes `max` over dynamics and every routed constraint and exchanges
  once, ADR-0010 decision 7).
- `col_ptr`, `cols`: constant per-functional column lists, flat/CSR (functional `f`
  owns `cols[col_ptr[f]:col_ptr[f+1]-1]`).
- `coeffs`: per-iterate partials, parallel to `cols`. **Mutable** — this is what a
  coefficient refresh writes.
- `row_map`: signed row-to-functional map; `sign` is the row's sign, `abs` its
  functional index.
- `row_offset`: per-row constant.

# Derived
`row_sign`, `row_functional`, `functional_row_ptr` (contiguity makes this
well-defined), the flat assembled-Jacobian entry tables, and the cached
`SparseMatrixCSC` the assembled entry point returns.

# Invariant
A functional's constraint rows are **contiguous**, enforced at construction. This is
what lets the residual write a functional's rows the moment its value is formed,
with no per-functional buffer and hence no eltype restriction (ForwardDiff works),
and it is what `#332`'s pair contraction will index by.
"""
struct ConstraintStencilTable
    # ── declared ─────────────────────────────────────────────────────────────── #
    n_functionals::Int
    n_rows::Int
    n_cols::Int
    stencil_width::Int
    col_ptr::Vector{Int}
    cols::Vector{Int}
    coeffs::Vector{Float64}
    row_map::Vector{Int}
    row_offset::Vector{Float64}
    # ── derived: row view ────────────────────────────────────────────────────── #
    row_sign::Vector{Float64}
    row_functional::Vector{Int}
    functional_row_ptr::Vector{Int}
    # ── derived: assembled-Jacobian entry tables (structure order) ───────────── #
    entry_row::Vector{Int}
    entry_col::Vector{Int}
    entry_coeff::Vector{Int}     # index into `coeffs`
    entry_sign::Vector{Float64}
    J::SparseMatrixCSC{Float64,Int}
    entry_nzpos::Vector{Int}     # entry k → index into `J.nzval`
    # ── derived: coefficient-generation counter ──────────────────────────────── #
    # Bumped by EVERY coefficient refresh (`stencil_touch!`, called from each
    # constraint's refresh barrier). #332's routed products cache "the coefficients
    # already match this iterate" on the committed-iterate stamp; but `coeffs` lives on
    # the CONSTRAINT, not on the caching closure, so a second consumer of the same
    # constraint — another callback set, or a plain `eval_jacobian` from the
    # factorisation evaluator — can rewrite it at a different iterate between two
    # products. The stamp alone would then serve a coefficient set belonging to some
    # other `x`: exactly #329's staleness defect one level down. A cached read is valid
    # only when the stamp matches AND this counter is the one the refresh recorded.
    refresh_token::Base.RefValue{Int}
end

function ConstraintStencilTable(
    functional_cols::AbstractVector,
    row_map::AbstractVector{Int},
    row_offset::AbstractVector{<:Real},
    n_cols::Int;
    stencil_width::Int,
)
    n_functionals = length(functional_cols)
    n_rows = length(row_map)

    n_functionals > 0 || throw(ArgumentError("stencil table declares no functionals"))
    n_rows > 0 || throw(ArgumentError("stencil table declares no constraint rows"))
    length(row_offset) == n_rows || throw(
        ArgumentError(
            "row_offset has length $(length(row_offset)), expected $n_rows (one per row)",
        ),
    )
    stencil_width >= 0 ||
        throw(ArgumentError("stencil_width must be non-negative, got $stencil_width"))
    n_cols > 0 || throw(ArgumentError("n_cols must be positive, got $n_cols"))

    # ── flatten the declared columns ─────────────────────────────────────────── #
    col_ptr = Vector{Int}(undef, n_functionals + 1)
    col_ptr[1] = 1
    for f = 1:n_functionals
        nc = length(functional_cols[f])
        nc > 0 || throw(ArgumentError("functional $f declares no columns"))
        col_ptr[f+1] = col_ptr[f] + nc
    end
    cols = Vector{Int}(undef, col_ptr[end] - 1)
    for f = 1:n_functionals
        for (j, c) in enumerate(functional_cols[f])
            1 <= c <= n_cols || throw(
                ArgumentError(
                    "functional $f column $c out of range 1:$n_cols — the column count " *
                    "must include the global variables",
                ),
            )
            cols[col_ptr[f]+j-1] = c
        end
    end
    coeffs = zeros(Float64, length(cols))

    # ── row view + the contiguity invariant ──────────────────────────────────── #
    row_sign = Vector{Float64}(undef, n_rows)
    row_functional = Vector{Int}(undef, n_rows)
    row_count = zeros(Int, n_functionals)
    first_row = zeros(Int, n_functionals)
    last_row = zeros(Int, n_functionals)
    for r = 1:n_rows
        m = row_map[r]
        m != 0 || throw(ArgumentError("row $r has a zero (unsigned) functional map entry"))
        f = abs(m)
        f <= n_functionals || throw(
            ArgumentError(
                "row $r maps to functional $f but only $n_functionals are declared",
            ),
        )
        row_functional[r] = f
        row_sign[r] = m > 0 ? 1.0 : -1.0
        row_count[f] += 1
        first_row[f] == 0 && (first_row[f] = r)
        last_row[f] = r
    end
    for f = 1:n_functionals
        row_count[f] > 0 || throw(
            ArgumentError(
                "functional $f is declared but no constraint row reads it — the table " *
                "must enumerate functionals from actual rows",
            ),
        )
        last_row[f] - first_row[f] + 1 == row_count[f] || throw(
            ArgumentError(
                "functional $f's constraint rows are not contiguous: $(row_count[f]) " *
                "rows spread over $(first_row[f]):$(last_row[f])",
            ),
        )
    end
    functional_row_ptr = Vector{Int}(undef, n_functionals + 1)
    for f = 1:n_functionals
        functional_row_ptr[f] = first_row[f]
    end
    functional_row_ptr[n_functionals+1] = last_row[n_functionals] + 1

    # ── assembled-Jacobian entry tables, in fill (structure) order ───────────── #
    n_entries = 0
    for r = 1:n_rows
        f = row_functional[r]
        n_entries += col_ptr[f+1] - col_ptr[f]
    end
    entry_row = Vector{Int}(undef, n_entries)
    entry_col = Vector{Int}(undef, n_entries)
    entry_coeff = Vector{Int}(undef, n_entries)
    entry_sign = Vector{Float64}(undef, n_entries)
    e = 0
    for r = 1:n_rows
        f = row_functional[r]
        s = row_sign[r]
        for j = col_ptr[f]:(col_ptr[f+1]-1)
            e += 1
            entry_row[e] = r
            entry_col[e] = cols[j]
            entry_coeff[e] = j
            entry_sign[e] = s
        end
    end

    # Cached assembled Jacobian. The pattern never moves, so this is built ONCE and
    # refilled through `entry_nzpos`; rebuilding a sparse matrix per call was ~6.9 MB
    # of the measured 6.88 MB per Jacobian at N=801.
    J = sparse(entry_row, entry_col, ones(Float64, n_entries), n_rows, n_cols)
    entry_nzpos = Vector{Int}(undef, n_entries)
    for e = 1:n_entries
        r, c = entry_row[e], entry_col[e]
        pos = 0
        @inbounds for k = J.colptr[c]:(J.colptr[c+1]-1)
            if J.rowval[k] == r
                pos = k
                break
            end
        end
        pos == 0 && error("internal: entry ($r, $c) missing from the assembled pattern")
        entry_nzpos[e] = pos
    end

    return ConstraintStencilTable(
        n_functionals,
        n_rows,
        n_cols,
        stencil_width,
        col_ptr,
        cols,
        coeffs,
        collect(Int, row_map),
        collect(Float64, row_offset),
        row_sign,
        row_functional,
        functional_row_ptr,
        entry_row,
        entry_col,
        entry_coeff,
        entry_sign,
        J,
        entry_nzpos,
        Ref(0),
    )
end

# ----------------------------------------------------------------------------- #
# Accessors
# ----------------------------------------------------------------------------- #

"""
    stencil_n_entries(table) -> Int

Number of assembled-Jacobian nonzero entries (rows × their functional's columns).
"""
stencil_n_entries(t::ConstraintStencilTable) = length(t.entry_row)

"""
    stencil_coeff_range(table, f) -> UnitRange{Int}

Indices into `table.cols` / `table.coeffs` owned by functional `f`. This is what a
coefficient refresh writes into.
"""
@inline stencil_coeff_range(t::ConstraintStencilTable, f::Int) =
    t.col_ptr[f]:(t.col_ptr[f+1]-1)

"""
    stencil_functional_rows(table, f) -> UnitRange{Int}

The contiguous constraint rows that read functional `f`.
"""
@inline stencil_functional_rows(t::ConstraintStencilTable, f::Int) =
    t.functional_row_ptr[f]:(t.functional_row_ptr[f+1]-1)

"""
    stencil_refresh_token(table) -> Int

The table's coefficient generation: a monotone counter bumped by every coefficient
refresh. A consumer that caches "the coefficients already match iterate `x`" must record
this alongside its iterate stamp and re-check both, because `coeffs` belongs to the
constraint and any other consumer may rewrite it in between. See the field's comment.
"""
@inline stencil_refresh_token(t::ConstraintStencilTable) = t.refresh_token[]

"""
    stencil_touch!(table)

Record that `table.coeffs` has just been rewritten. Called at the end of every
constraint's coefficient-refresh barrier; allocation-free.
"""
@inline function stencil_touch!(t::ConstraintStencilTable)
    t.refresh_token[] += 1
    return nothing
end

# ----------------------------------------------------------------------------- #
# Generic application
# ----------------------------------------------------------------------------- #

"""
    stencil_structure(table) -> Vector{Tuple{Int,Int}}

The assembled Jacobian's `(row, column)` sparsity structure, in the same order
[`stencil_fill_values!`](@ref) writes. Built once; this is what a constraint returns
from `jacobian_structure`.
"""
function stencil_structure(t::ConstraintStencilTable)
    return [(t.entry_row[e], t.entry_col[e]) for e = 1:stencil_n_entries(t)]
end

"""
    stencil_fill_values!(vals, table)

Write the assembled Jacobian's nonzero values into `vals`, in
[`stencil_structure`](@ref) order: `vals[e] = row_sign * coeff`. Allocation-free.
"""
function stencil_fill_values!(vals::AbstractVector, t::ConstraintStencilTable)
    n = stencil_n_entries(t)
    length(vals) >= n ||
        throw(DimensionMismatch("value buffer has length $(length(vals)), need $n"))
    coeffs = t.coeffs
    sgn = t.entry_sign
    src = t.entry_coeff
    @inbounds for e = 1:n
        vals[e] = sgn[e] * coeffs[src[e]]
    end
    return vals
end

"""
    stencil_assemble!(table) -> SparseMatrixCSC

Refresh and return the CACHED assembled Jacobian. The pattern is built once at
construction, so this only overwrites `nzval` — the returned matrix is the same
object on every call (the pre-allocated-reference convention already used by
`CubicSplineBoundConstraint` and DirectTrajOpt's `KnotPointConstraint`).
"""
function stencil_assemble!(t::ConstraintStencilTable)
    J = t.J
    nz = J.nzval
    coeffs = t.coeffs
    sgn = t.entry_sign
    src = t.entry_coeff
    pos = t.entry_nzpos
    @inbounds for e in eachindex(pos)
        nz[pos[e]] = sgn[e] * coeffs[src[e]]
    end
    return J
end

"""
    stencil_scatter_functional!(g, table, f, val)

Write functional `f`'s value into every constraint row that reads it:
`g[r] = row_sign[r] * val + row_offset[r]`.

Generic in `eltype(g)` and in `typeof(val)` (ForwardDiff duals included) and
allocation-free, because the invariant that `f`'s rows are contiguous means no
per-functional buffer is needed — the value is consumed the moment it is formed.
"""
@inline function stencil_scatter_functional!(
    g::AbstractVector,
    t::ConstraintStencilTable,
    f::Int,
    val,
)
    sgn = t.row_sign
    off = t.row_offset
    @inbounds for r = t.functional_row_ptr[f]:(t.functional_row_ptr[f+1]-1)
        g[r] = sgn[r] * val + off[r]
    end
    return nothing
end

"""
    stencil_expand_rows!(g, table, fvals)

Expand a full vector of functional values into constraint rows. Convenience wrapper
over [`stencil_scatter_functional!`](@ref); the hot residual path scatters as it goes
instead.
"""
function stencil_expand_rows!(
    g::AbstractVector,
    t::ConstraintStencilTable,
    fvals::AbstractVector,
)
    length(fvals) == t.n_functionals || throw(
        DimensionMismatch("fvals has length $(length(fvals)), expected $(t.n_functionals)"),
    )
    for f = 1:t.n_functionals
        stencil_scatter_functional!(g, t, f, fvals[f])
    end
    return g
end

# ----------------------------------------------------------------------------- #
# Matrix-free application (#332)
# ----------------------------------------------------------------------------- #
#
# THE REASON THESE EXIST. The assembled Jacobian is cheap to refill, so the naive
# reading is "cache it and pay a SpMV". That is wrong for a structural reason: a
# double-sided inequality's `±` row pair shares ONE underlying scalar, so one
# directional-derivative evaluation serves both rows. A sparse matrix stores the two
# rows' stencils independently and re-does the contraction per row; these kernels
# contract ONCE PER FUNCTIONAL instead.
#
#   JVP: contract the functional's stencil against `v` once  →  write both rows.
#   VJP: contract the pair's WEIGHTS first (`ω = Σ_r sign_r · w_r`)  →  scatter once.
#
# Combine coefficients, apply once — the constraint-side analogue of the MultiKet
# parameter-block reassociation. An implementation that scatters per row forfeits the
# entire win (ADR-0010 decision 3, CONTEXT.md "Paired constraint rows").
#
# Both are range reductions rather than gathers because a functional's rows are
# CONTIGUOUS, an invariant this table enforces at construction.

"""
    stencil_jvp!(Jv, table, v) -> Jv

Matrix-free Jacobian-vector product: `Jv[r] = row_sign[r] · (∂F_f/∂z ⋅ v)` for every
constraint row `r` reading functional `f`.

One contraction per FUNCTIONAL, then every row of its `±` pair is written from that one
scalar — this is where the win over the equivalent SpMV comes from.

`Jv` is **overwritten** over `1:table.n_rows` (a routed block owns a disjoint row range,
so no accumulation is needed and none is done). `v` is the FULL decision vector, so
`length(v)` must be `table.n_cols`. `row_offset` plays no part: it is the residual's
constant and has no derivative.

Allocation-free.
"""
function stencil_jvp!(
    Jv::AbstractVector{Float64},
    t::ConstraintStencilTable,
    v::AbstractVector{Float64},
)
    length(Jv) == t.n_rows ||
        throw(DimensionMismatch("Jv has length $(length(Jv)), expected $(t.n_rows) rows"))
    length(v) == t.n_cols || throw(
        DimensionMismatch(
            "v has length $(length(v)), expected the full decision width $(t.n_cols)",
        ),
    )
    cols = t.cols
    coeffs = t.coeffs
    col_ptr = t.col_ptr
    sgn = t.row_sign
    rp = t.functional_row_ptr
    @inbounds for f = 1:t.n_functionals
        # TWO accumulators, not one. A stencil is five to eight columns long, so the
        # single-accumulator form is a serial chain of that many dependent FMAs with
        # nothing to overlap them — measured, splitting it is worth ~1.2-1.3× on the whole
        # kernel (33.3 → 27.8 µs at N=801 for the continuity constraint). Anything wider
        # stops paying: at this stencil length the gather, not the chain, becomes the
        # limit. The pairwise order is why parity with the assembled matvec is asserted at
        # `rtol`, not bitwise — both are correct summations of the same terms.
        lo = col_ptr[f]
        hi = col_ptr[f+1] - 1
        d1 = 0.0
        d2 = 0.0
        j = lo
        while j + 1 <= hi
            d1 += coeffs[j] * v[cols[j]]
            d2 += coeffs[j+1] * v[cols[j+1]]
            j += 2
        end
        j <= hi && (d1 += coeffs[j] * v[cols[j]])
        d = d1 + d2
        for r = rp[f]:(rp[f+1]-1)
            Jv[r] = sgn[r] * d
        end
    end
    return Jv
end

"""
    stencil_vjp!(JTw, table, w) -> JTw

Matrix-free transpose product: `JTw += Jᵀw`.

The pair's weights are contracted FIRST — `ω_f = Σ_{r reads f} row_sign[r] · w[r]`, a
reduction over a contiguous range — and only then is the functional's stencil scattered,
once. Scattering per row instead would touch every partial twice and forfeit the win.

**ACCUMULATES** into `JTw` (unlike [`stencil_jvp!`](@ref), which overwrites): several
routed constraints share one `n_cols`-wide output, and so does the dense remainder. The
caller zeroes. `length(JTw)` must be `table.n_cols`, `length(w)` must be `table.n_rows`
(pass a `view` of the block's rows).

Allocation-free.
"""
function stencil_vjp!(
    JTw::AbstractVector{Float64},
    t::ConstraintStencilTable,
    w::AbstractVector{Float64},
)
    length(JTw) == t.n_cols || throw(
        DimensionMismatch(
            "JTw has length $(length(JTw)), expected the full decision width $(t.n_cols)",
        ),
    )
    length(w) == t.n_rows ||
        throw(DimensionMismatch("w has length $(length(w)), expected $(t.n_rows) rows"))
    cols = t.cols
    coeffs = t.coeffs
    col_ptr = t.col_ptr
    sgn = t.row_sign
    rp = t.functional_row_ptr
    @inbounds for f = 1:t.n_functionals
        ω = 0.0
        for r = rp[f]:(rp[f+1]-1)
            ω += sgn[r] * w[r]
        end
        for j = col_ptr[f]:(col_ptr[f+1]-1)
            JTw[cols[j]] += ω * coeffs[j]
        end
    end
    return JTw
end

# ----------------------------------------------------------------------------- #
# The constraint-type gradient trait (#332)
# ----------------------------------------------------------------------------- #
#
# ADR-0010 decisions 1 and 2: the gradient method for an algebraic constraint is
# dispatched on the CONSTRAINT TYPE (constraints have no integration algorithm, so
# ADR-0003's algorithm-dispatch rule does not reach them), and the vocabulary stays
# INSIDE the package that owns the constraint types (Piccolo since the open-core
# split's 3c slice — re-exported by Piccolissimo —, so an upstream split buys no
# consumer). This is the inequality-path sibling of the backend's
# `_supports_matrix_free_eq_gradient`.
#
# HOW A CONSTRAINT OPTS IN — exactly two methods:
#
#   constraint_stencil_table(c::MyConstraint) = c.table
#   refresh_constraint_coefficients!(c::MyConstraint, traj) =
#       _my_refresh_coefficients!(traj.datavec, c.tag, c.table)
#
# `supports_matrix_free_constraint_gradient` is then DERIVED and must not be overridden:
# routability is a property of the declared table, not an independent opt-in flag, so a
# constraint cannot claim to be routable without supplying a table whose stencil width is
# bounded. That is what "refused at the trait, not silently routed" means.
#
# #333 (device kernels) and #334 (the sharded path) extend this same trait rather than
# adding a parallel one: they consume `constraint_stencil_table` (the declared structure
# they upload / slab-index) and `refresh_constraint_coefficients!` (the per-iterate
# cadence), and gate on the same predicate.

"""
    UNBOUNDED_STENCIL_WIDTH

The stencil width a constraint declares when its rows are **not** knot-local — a
periodicity row coupling knot 1 with knot `N`, or a global `∫u dt` budget coupling all of
them. Such a constraint is refused by
[`supports_matrix_free_constraint_gradient`](@ref) and stays on the assembled path: its
rows are not owned by one rank, so it needs non-adjacent point-to-point or a reduction
rather than the single bounded halo the driver exchanges (ADR-0010 decision 7).
"""
const UNBOUNDED_STENCIL_WIDTH = typemax(Int)

"""
    constraint_stencil_table(constraint) -> Union{Nothing,ConstraintStencilTable}

The constraint's declared stencil table, or `nothing` (the default) for a constraint that
has none — the three ForwardDiff-backed types, which stay on the assembled path.

One half of the gradient trait; see this section's comment for the two-method opt-in.
"""
constraint_stencil_table(::Any) = nothing

"""
    refresh_constraint_coefficients!(constraint, traj)

Overwrite the constraint's stencil coefficients from `traj` for the current iterate, and
[`stencil_touch!`](@ref) the table. The other half of the gradient trait.

Called on the **once-per-accepted-iterate** cadence, never per product: keeping the
constant columns separate from the per-iterate coefficients is exactly what buys that
cadence (ADR-0010 decision 4), and re-refreshing per product would cost more than the
SpMV the kernels are replacing.
"""
function refresh_constraint_coefficients!(c, ::Any)
    throw(
        ArgumentError(
            "$(typeof(c)) declares a stencil table but no `refresh_constraint_coefficients!` " *
            "method — a routable constraint must supply BOTH halves of the gradient trait",
        ),
    )
end

"""
    supports_matrix_free_constraint_gradient(constraint) -> Bool

Whether the backend should serve this constraint's inequality rows from the matrix-free
[`stencil_jvp!`](@ref) / [`stencil_vjp!`](@ref) kernels instead of from a block of the
assembled inequality Jacobian.

DERIVED, never overridden: `true` exactly when the constraint declares a stencil table
whose stencil width is bounded. A constraint whose rows are not knot-local declares
[`UNBOUNDED_STENCIL_WIDTH`](@ref) and is refused here — at the trait, rather than being
silently routed into a kernel whose halo assumption it violates.
"""
supports_matrix_free_constraint_gradient(c) = _routable_stencil(constraint_stencil_table(c))

_routable_stencil(::Nothing) = false
_routable_stencil(t::ConstraintStencilTable) =
    0 <= t.stencil_width < UNBOUNDED_STENCIL_WIDTH

# The HVP sibling of that trait (#458) lives right below, with the kernels it names.

@testitem "ConstraintStencilTable: structure, ± pairing, cached assembly" begin

    using Piccolo:
        ConstraintStencilTable,
        stencil_assemble!,
        stencil_coeff_range,
        stencil_expand_rows!,
        stencil_fill_values!,
        stencil_functional_rows,
        stencil_n_entries,
        stencil_structure
    using SparseArrays

    # Two functionals; each read by a `+`/`−` row pair — the double-sided shape.
    fcols = [[1, 2, 5], [2, 3]]
    row_map = [1, -1, 2, -2]
    row_offset = [-10.0, -10.0, 0.0, 0.0]
    t = ConstraintStencilTable(fcols, row_map, row_offset, 7; stencil_width = 2)

    @test t.n_functionals == 2
    @test t.n_rows == 4
    @test t.n_cols == 7            # includes the global columns, carried as data
    @test t.stencil_width == 2
    @test stencil_n_entries(t) == 3 + 3 + 2 + 2
    @test collect(stencil_coeff_range(t, 1)) == [1, 2, 3]
    @test collect(stencil_coeff_range(t, 2)) == [4, 5]
    @test stencil_functional_rows(t, 1) == 1:2
    @test stencil_functional_rows(t, 2) == 3:4

    # Structure is row-major, each row carrying its functional's columns in order.
    @test stencil_structure(t) ==
          [(1, 1), (1, 2), (1, 5), (2, 1), (2, 2), (2, 5), (3, 2), (3, 3), (4, 2), (4, 3)]

    # A coefficient refresh writes `coeffs`; application is generic over it.
    t.coeffs .= [1.0, 2.0, 3.0, 4.0, 5.0]
    vals = zeros(stencil_n_entries(t))
    stencil_fill_values!(vals, t)
    @test vals == [1.0, 2.0, 3.0, -1.0, -2.0, -3.0, 4.0, 5.0, -4.0, -5.0]

    J = stencil_assemble!(t)
    @test size(J) == (4, 7)
    @test J[1, 1] == 1.0 && J[1, 2] == 2.0 && J[1, 5] == 3.0
    @test J[2, 1] == -1.0 && J[2, 2] == -2.0 && J[2, 5] == -3.0
    @test J[3, 2] == 4.0 && J[3, 3] == 5.0
    @test J[4, 2] == -4.0 && J[4, 3] == -5.0

    # The assembled matrix is a CACHED reference refilled in place — no per-call
    # sparse rebuild. Same object, new values.
    J_first = J
    t.coeffs .*= 2
    J_second = stencil_assemble!(t)
    @test J_second === J_first
    @test J_second[1, 1] == 2.0

    # Row expansion applies the signed map and the per-row constant.
    g = zeros(4)
    stencil_expand_rows!(g, t, [7.0, -1.5])
    @test g == [7.0 - 10.0, -7.0 - 10.0, -1.5, 1.5]

    # Allocation-free application (both are pure index arithmetic).
    fv = [7.0, -1.5]
    fill_alloc = () -> stencil_fill_values!(vals, t)
    asm_alloc = () -> stencil_assemble!(t)
    exp_alloc = () -> stencil_expand_rows!(g, t, fv)
    for f in (fill_alloc, asm_alloc, exp_alloc)
        f()
        f()
        @test minimum(@allocated(f()) for _ = 1:3) == 0
    end
end

@testitem "ConstraintStencilTable: the contiguity invariant is enforced at construction" begin
    using Piccolo: ConstraintStencilTable

    fcols = [[1, 2], [2, 3]]

    # Deliberately MIS-ORDERED: functional 1's two rows are split by functional 2's.
    # Assuming contiguity instead of enforcing it would silently corrupt every
    # functional-indexed application (the residual writes a functional's rows as one
    # contiguous span).
    @test_throws ArgumentError ConstraintStencilTable(
        fcols,
        [1, 2, -1, -2],
        zeros(4),
        6;
        stencil_width = 1,
    )

    # The correctly-grouped version of the same declaration constructs.
    t = ConstraintStencilTable(fcols, [1, -1, 2, -2], zeros(4), 6; stencil_width = 1)
    @test t.n_rows == 4

    # …and the other declaration errors, each on its own message.
    @test_throws ArgumentError ConstraintStencilTable(
        fcols,
        [1, -1, 0, -2],           # zero map entry: no sign
        zeros(4),
        6;
        stencil_width = 1,
    )
    @test_throws ArgumentError ConstraintStencilTable(
        fcols,
        [1, -1, 3, -3],           # functional out of range
        zeros(4),
        6;
        stencil_width = 1,
    )
    @test_throws ArgumentError ConstraintStencilTable(
        fcols,
        [1, -1],                  # functional 2 declared but unread
        zeros(2),
        6;
        stencil_width = 1,
    )
    @test_throws ArgumentError ConstraintStencilTable(
        fcols,
        [1, -1, 2, -2],
        zeros(4),
        2;                        # a column beyond n_cols — the globals-omitted bug
        stencil_width = 1,
    )
    @test_throws ArgumentError ConstraintStencilTable(
        fcols,
        [1, -1, 2, -2],
        zeros(3),                 # row_offset/row_map length disagreement
        6;
        stencil_width = 1,
    )
end

@testitem "#332 stencil kernels: JVP/VJP match the assembled action, adjoint identity, ± reassociation" begin

    using Piccolo:
        ConstraintStencilTable,
        stencil_assemble!,
        stencil_jvp!,
        stencil_refresh_token,
        stencil_touch!,
        stencil_vjp!
    using LinearAlgebra
    using Random

    Random.seed!(0x332001)

    # A `±` pair on functional 1, an ASYMMETRIC row group on functional 2 (three rows,
    # mixed signs — the shape a per-row kernel would get right by accident and a paired
    # one must get right on purpose), and a single-row functional 3.
    fcols = [[1, 2, 7], [2, 3, 4], [5, 6]]
    row_map = [1, -1, 2, -2, 2, 3]
    row_offset = [-1.0, -2.0, 0.5, 0.25, -3.0, 4.0]
    n_cols = 9
    t = ConstraintStencilTable(fcols, row_map, row_offset, n_cols; stencil_width = 1)
    t.coeffs .= randn(length(t.coeffs))

    J = Matrix(stencil_assemble!(t))          # the ORACLE: the assembled Jacobian
    v = randn(n_cols)
    w = randn(t.n_rows)

    # ── Parity with the assembled action, to machine precision. ──────────────────── #
    Jv = zeros(t.n_rows)
    stencil_jvp!(Jv, t, v)
    @test Jv ≈ J * v atol = 1e-14

    JTw = zeros(n_cols)
    stencil_vjp!(JTw, t, w)
    @test JTw ≈ J' * w atol = 1e-14

    # ── Adjoint identity ⟨Jv, w⟩ == ⟨v, Jᵀw⟩, kernel against kernel. ─────────────── #
    @test dot(Jv, w) ≈ dot(v, JTw) atol = 1e-13

    # ── The `±` reassociation, asserted on VALUES rather than on timing: the pair's
    #    rows must be exact negatives of one another in the JVP (one contraction, two
    #    signed writes), and a pair carrying EQUAL weights must contribute exactly
    #    nothing to the VJP because ω = w₊ - w₋ = 0 cancels BEFORE any partial is
    #    touched. An implementation that scatters per row still gets the first; only one
    #    that contracts the weights first gets the second exactly.
    @test Jv[1] == -Jv[2]
    w_bal = zeros(t.n_rows)
    w_bal[1] = 0.75
    w_bal[2] = 0.75
    JTw_bal = zeros(n_cols)
    stencil_vjp!(JTw_bal, t, w_bal)
    @test all(iszero, JTw_bal)

    # ── VJP ACCUMULATES (several routed blocks share one output). ─────────────────── #
    acc = copy(JTw)
    stencil_vjp!(acc, t, w)
    @test acc ≈ 2 .* JTw atol = 1e-13
    # …and the JVP OVERWRITES.
    Jv_dirty = fill(1e6, t.n_rows)
    stencil_jvp!(Jv_dirty, t, v)
    @test Jv_dirty == Jv

    # ── Allocation-free. Local-scoped closures only. ─────────────────────────────── #
    jvp_f = () -> stencil_jvp!(Jv, t, v)
    vjp_f = () -> stencil_vjp!(JTw, t, w)
    for f in (jvp_f, vjp_f)
        f()
        f()
        @test minimum(@allocated(f()) for _ = 1:3) == 0
    end

    # ── Dimension contracts are enforced, not silently truncated. ────────────────── #
    @test_throws DimensionMismatch stencil_jvp!(zeros(t.n_rows), t, zeros(n_cols - 1))
    @test_throws DimensionMismatch stencil_jvp!(zeros(t.n_rows - 1), t, v)
    @test_throws DimensionMismatch stencil_vjp!(zeros(n_cols - 1), t, w)
    @test_throws DimensionMismatch stencil_vjp!(zeros(n_cols), t, zeros(t.n_rows - 1))

    # ── The refresh token is monotone and bumped only by an explicit touch. ──────── #
    tok = stencil_refresh_token(t)
    stencil_jvp!(Jv, t, v)
    @test stencil_refresh_token(t) == tok      # a product never touches the coefficients
    stencil_touch!(t)
    @test stencil_refresh_token(t) == tok + 1
    touch_f = () -> stencil_touch!(t)
    touch_f()
    @test minimum(@allocated(touch_f()) for _ = 1:3) == 0
end

@testitem "#332 gradient trait: derived from the declared table, unbounded width refused" begin
    using Piccolo

    using Piccolo:
        ConstraintStencilTable,
        UNBOUNDED_STENCIL_WIDTH,
        constraint_stencil_table,
        refresh_constraint_coefficients!,
        supports_matrix_free_constraint_gradient

    fcols = [[1, 2], [2, 3]]
    bounded = ConstraintStencilTable(fcols, [1, -1, 2, -2], zeros(4), 5; stencil_width = 1)
    unbounded = ConstraintStencilTable(
        fcols,
        [1, -1, 2, -2],
        zeros(4),
        5;
        stencil_width = UNBOUNDED_STENCIL_WIDTH,
    )

    # A constraint with no table at all — the three ForwardDiff-backed types' situation.
    struct _NoTableCon end
    @test constraint_stencil_table(_NoTableCon()) === nothing
    @test supports_matrix_free_constraint_gradient(_NoTableCon()) == false

    # Opting in is two methods: the table, and the refresh.
    struct _BoundedCon
        table::ConstraintStencilTable
    end
    Piccolo.Control.QuantumConstraints.SplineConstraints.constraint_stencil_table(
        c::_BoundedCon,
    ) = c.table
    @test supports_matrix_free_constraint_gradient(_BoundedCon(bounded)) == true

    # A NON-knot-local constraint is refused AT THE TRAIT, not silently routed into a
    # kernel whose one-knot-halo assumption it violates.
    @test supports_matrix_free_constraint_gradient(_BoundedCon(unbounded)) == false

    # Declaring a table without a refresh is a loud error, not a wrong answer.
    @test_throws ArgumentError refresh_constraint_coefficients!(
        _BoundedCon(bounded),
        nothing,
    )
end

# ----------------------------------------------------------------------------- #
# The inequality HVP trait (#458)                                               #
# ----------------------------------------------------------------------------- #

"""
    supports_matrix_free_constraint_hvp(constraint) -> Bool

Whether this constraint supplies an exact per-functional Hessian action the backend's
inequality path can compose into `compute_inequality_constraint_hvp!` — the matrix-free
action `Σ_f ω_f ∇²F_f · v` over the constraint's stencil functionals, with
`ω_f = Σ_{r reads f} row_sign[r] · w[r]` ([`stencil_functional_weight`](@ref)).

**DECLARED, never derived, never delegated** (ADR-0010 decision 1's vocabulary): there is
no structural property of a [`ConstraintStencilTable`](@ref) that implies a second-order
action — the table carries first-order coefficient data only — so unlike
[`supports_matrix_free_constraint_gradient`](@ref) this trait is an explicit per-family
opt-in. The default is `false`, and a family opts in by declaring `true` AND supplying a
[`constraint_stencil_hvp!`](@ref) method (whose default throws — declaring support without
the action is a loud error, mirroring the `refresh_constraint_coefficients!` pairing).
Constraints without the trait keep the Gauss-Newton fold (`ρJᵀJ`) as their entire
inequality curvature; nothing infers a Hessian through AD at the backend.
"""
supports_matrix_free_constraint_hvp(::Any) = false

"""
    stencil_functional_weight(table, w, f) -> Float64

Contract the stacked multiplier weights over functional `f`'s (contiguous) constraint
rows: `ω_f = Σ_r row_sign[r] · w[r]`. This is the VJP's weight contraction standing alone,
and it is the ONLY weight information an exact inequality HVP needs: constraint rows are
signed copies of functionals, `g[r] = row_sign[r] · F[row_functional[r]] + row_offset[r]`,
so `Σ_r w_r ∇²g_r = Σ_f ω_f ∇²F_f` — the `row_offset` drops out as an affine constant and
the `±` pair collapses into one scalar weight before any second-order work happens.

`length(w)` must be `table.n_rows` (pass a `view` of the block's rows). Allocation-free.
"""
@inline function stencil_functional_weight(
    t::ConstraintStencilTable,
    w::AbstractVector{Float64},
    f::Int,
)
    length(w) == t.n_rows ||
        throw(DimensionMismatch("w has length $(length(w)), expected $(t.n_rows) rows"))
    ω = 0.0
    sgn = t.row_sign
    rp = t.functional_row_ptr
    @inbounds for r = rp[f]:(rp[f+1]-1)
        ω += sgn[r] * w[r]
    end
    return ω
end

"""
    constraint_stencil_hvp!(Hv, constraint, w, v, x)

Add the constraint's exact inequality Hessian-of-Lagrangian action to `Hv`:
`Hv += Σ_f ω_f ∇²F_f(x) · v` over the constraint's stencil functionals, with the
multiplier weights `w` over the constraint's OWN rows ([`stencil_functional_weight`](@ref)
does the contraction) and `v` the full decision vector.

**ACCUMULATES** into `Hv` (unlike the routed JVP, which overwrites): several routed
constraints share one `n_cols`-wide output, and so does the dense remainder. The callback
that calls this zeroes first. `length(Hv)` must be the full decision width
`table.n_cols`.

One method per opting-in family — the constraint owns its second-order math exactly as it
owns its coefficient refresh. The `::Any` default throws: reaching it means the caller
gated on [`supports_matrix_free_constraint_hvp`](@ref) but the family never supplied its
action, which is a wiring bug, not a runtime contingency.
"""
function constraint_stencil_hvp!(
    Hv::AbstractVector{Float64},
    c::Any,
    w::AbstractVector{Float64},
    v::AbstractVector{Float64},
    x::AbstractVector{Float64},
)
    throw(
        ArgumentError(
            "$(typeof(c)) declares `supports_matrix_free_constraint_hvp` but supplies no " *
            "`constraint_stencil_hvp!` method — an opting-in family must supply BOTH " *
            "halves, exactly like the gradient trait",
        ),
    )
end

@testitem "#458 HVP trait: declared per family, default false; weight contraction" begin

    using Piccolo:
        ConstraintStencilTable,
        stencil_functional_weight,
        supports_matrix_free_constraint_hvp

    # The default is false and is stated on Any: a ForwardDiff-backed constraint (or any
    # type at all) never carries an HVP capability into the backend.
    struct _NoHVPCon end
    @test supports_matrix_free_constraint_hvp(_NoHVPCon()) == false

    # The weight contraction is the VJP's, standing alone: signed, over the contiguous
    # row range, on the FUNCTIONALS — the ± pair collapses into one scalar.
    fcols = [[1, 2, 7], [2, 3, 4], [5, 6]]
    row_map = [1, -1, 2, -2, 2, 3]
    t = ConstraintStencilTable(fcols, row_map, zeros(6), 9; stencil_width = 1)
    w = [0.25, 0.75, -1.5, 2.0, 0.5, -0.5]
    # functional 1: rows 1,2 (signs +,−)  →  0.25 − 0.75
    @test stencil_functional_weight(t, w, 1) == 0.25 - 0.75
    # functional 2: rows 3,4,5 (signs +,−,+ — the asymmetric group)
    #   →  −1.5 − 2.0 + 0.5
    @test stencil_functional_weight(t, w, 2) == -1.5 - 2.0 + 0.5
    # functional 3: row 6 (sign +)  →  −0.5
    @test stencil_functional_weight(t, w, 3) == -0.5

    # A balanced ± pair contracts to EXACTLY zero — the whole pairing win, one level up.
    w_bal = zeros(6)
    w_bal[1] = 0.75
    w_bal[2] = 0.75
    @test stencil_functional_weight(t, w_bal, 1) == 0.0

    # Dimension contract enforced, not silently truncated.
    @test_throws DimensionMismatch stencil_functional_weight(t, zeros(5), 1)

    # Allocation-free. Local-scoped `let` bindings on purpose: a closure at the
    # testitem's module scope captures `Any`-typed globals and boxes the scalar
    # return — the same measurement fault the #332 testitem's comment documents.
    let wt = t, ww = w
        wf = () -> stencil_functional_weight(wt, ww, 2)
        wf()
        wf()
        @test minimum(@allocated(wf()) for _ = 1:3) == 0
    end
end
