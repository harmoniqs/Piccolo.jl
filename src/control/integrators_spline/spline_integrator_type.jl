export SplineIntegrator, canonical_hessian_knot_dim
export SplineType, LinearSpline, CubicSpline
export IntegrationAlgorithm, Tsit5Alg, MagnusGL4Alg, MagnusAdapt4Alg
export Tsit5Data, MagnusGL4Data, MagnusAdapt4Data
export eval_constraint_jvp!, eval_constraint_vjp!, eval_constraint_hvp!
export build_multiket_hvp_cache

# NOTE: the SplineType trait hierarchy and `_drift_matrix` moved to
# ../shared/spline_types.jl (pure move, #132) so alg_data.jl can dispatch on
# SplineType.

"""
    SplineIntegrator{T, S, R, A, D} <: AbstractSplineIntegrator

Unified spline integrator with algorithm dispatch for forward propagation.

# Type Parameters
- `T`: Quantum trajectory type (UnitaryTrajectory, KetTrajectory, MultiKetTrajectory, DensityTrajectory)
- `S`: Spline type (LinearSpline or CubicSpline)
- `R`: Numeric type for propagator (ComplexF64 or Float64)
- `A`: Integration algorithm (Tsit5Alg, MagnusGL4Alg, ...)
- `D`: Algorithm-specific data type (Tsit5Data{R} or MagnusGL4Data)

# Fields
- `x_names::Vector{Symbol}`: Names of state variables (single element for non-ensemble)
- `u_name::Symbol`: Name of control variable
- `x_dim::Int`: Total constraint dimension per knot point
- `u_dim::Int`: Control dimension
- `dim::Int`: Total constraint dimension
- `tol::Float64`: Tolerance for integration
- `prop_results::Vector{PropagatorResult{R}}`: Complex propagator + sensitivities per knot point
- `ketdim::Int`: Dimension of ket space
- `global_names::Vector{Symbol}`: Names of global (time-invariant) variables
- `global_dim::Int`: Dimension of global variables
- `sens_probs`: Sensitivity-augmented ODE problems for analytical Jacobian
- `sens_state`: Preallocated buffer for sensitivity ODE solution
- `alg::A`: Integration algorithm instance
- `alg_data::D`: Algorithm-specific pre-allocated data

# ODE Parameter Layout
Parameters passed to ODE: `[uₖ, uₖ₊₁, (duₖ, duₖ₊₁), Δtₖ, tₖ]`
where `uₖ = [controls..., globals...]` - globals are embedded in u vectors
- Linear spline: 2*(control_dim + global_dim) for uₖ and uₖ₊₁
- Cubic spline: 4*(control_dim + global_dim) for uₖ, uₖ₊₁, duₖ, duₖ₊₁
- Note: duₖ includes zeros for global derivatives since globals don't change

# Internal Representation
The ODE runs in complex domain (ℂ^{n×n}). Propagator Φ and sensitivities S_j are
stored as ComplexF64 in `prop_results`. Conversion to real isomorphic form happens
at the Jacobian/Hessian assembly boundary via dispatch to trajectory-specific `fill_*` methods.
"""
struct SplineIntegrator{
    T<:AbstractQuantumTrajectory,
    S<:SplineType,
    R<:Number,
    A<:IntegrationAlgorithm,
    D,
} <: AbstractSplineIntegrator
    x_names::Vector{Symbol}  # Vector for ensemble support; single element for non-ensemble
    u_name::Symbol
    x_dim::Int
    u_dim::Int
    dim::Int
    tol::Float64
    prop_results::Vector{PropagatorResult{R}}
    ketdim::Int
    global_names::Vector{Symbol}
    global_dim::Int
    # Sensitivity ODE support (analytical Jacobian computation)
    # Nothing for function-based systems that lack explicit H_drives matrices
    sens_probs::Union{Nothing,Vector{ODEProblem}}
    sens_state::Union{Nothing,Vector{R}}
    # Second-order sensitivity ODE support (exact Hessian computation)
    exact_hessian::Bool
    hess_probs::Union{Nothing,Vector{ODEProblem}}
    hess_state::Union{Nothing,Vector{R}}
    hess_active_pairs::Union{Nothing,Vector{Tuple{Int,Int}}}
    hess_n_params::Int
    # Ket-level sensitivity support (MultiKetTrajectory only)
    # When true, sensitivity ODE propagates K n-dim ket vectors instead of n×n matrices.
    # Reduces sensitivity cost from O(n²J) to O(KnJ) per RHS evaluation.
    use_ket_sensitivity::Bool
    ket_sens_results::Union{Nothing,Vector{Matrix{R}}}  # per-knot (K*ketdim, n_params); K = n_kets, not the knot count
    # Algorithm dispatch
    alg::A
    alg_data::D
end

# Convenience accessor for single state name (non-ensemble case)
x_name(𝒮::SplineIntegrator) = 𝒮.x_names[1]

# SplineType trait instance of an integrator (for SplineType-dispatched seams)
_spline_type(::SplineIntegrator{T,S}) where {T,S} = S()
single_state_dim(𝒮::SplineIntegrator) = 𝒮.x_dim ÷ length(𝒮.x_names)

# Propagator side dimension: n for ket/unitary/multiket (complex), n² for density (compact real)
# Used for reshaping PropagatorResult in shared eval_jacobian/eval_hessian
propagator_side_dim(𝒮::SplineIntegrator) = 𝒮.ketdim
propagator_side_dim(𝒮::SplineIntegrator{DensityTrajectory}) = 𝒮.ketdim^2
propagator_side_dim(𝒮::SplineIntegrator{MultiDensityTrajectory}) = 𝒮.ketdim^2

"""
    matrix_free_step_scratch(𝒮, x_dim, n_ode_params) -> cell scratch

The per-**task** scratch a cell's inner kernels take. The driver calls this ONCE
per spawned task (never per knot), so this is where a cell declares its buffer set
without the outer layer knowing its shape — which is what lets ONE outer body serve
a `1`-wide core (Ket) and a `d`-wide one (Unitary, #345).

Every dimension must trace to the constructing trajectory (`ketdim`, `x_dim`,
`n_ode_params` from `matrix_free_layout`) — never to a global `N` (ADR-0009).

Methods live with their cells: `KetStepScratch` in `spline_integrator_ket.jl`,
`UnitaryStepScratch` in `spline_integrator_unitary.jl`.
"""
function matrix_free_step_scratch end

"""
    _refresh_prop_results!(𝒮, traj, globals)

Solve the per-knot first-order sensitivity ODE for every interval of `traj`, in
parallel, filling `𝒮.prop_results`. Disjoint writes per `k` ⇒ thread-safe.

This is the ONE place the refresh loop lives, so a per-iterate cache builder can
reuse it without a representation file having to spawn anything of its own.
"""
function _refresh_prop_results!(𝒮::SplineIntegrator, traj::NamedTrajectory, globals)
    N = traj.N
    Threads.@threads for k = 1:(N-1)
        compute_ode_jacobian!(𝒮, traj[k], traj[k+1], k, globals)
    end
    return nothing
end

"""
    refresh_sensitivities!(𝒮, traj, globals) -> nothing

The refresh half of `_refresh_or_restore_sensitivities!`, as a CELL SEAM.

The generic method is `_refresh_prop_results!` and that is the whole of it
numerically. It exists as a seam because a cell may need to OBSERVE its own refresh
passes: the Ket cell counts them (`ket_sens_refresh_count`, #340's witness that "two
applies at the same decision vector solve once, not twice"), and that witness has to
keep counting once the shared outer body of #338/#345 drives the refresh instead of
the cell's own body. A cell that wants no bookkeeping declares nothing.
"""
refresh_sensitivities!(𝒮::SplineIntegrator, traj::NamedTrajectory, globals) =
    _refresh_prop_results!(𝒮, traj, globals)

# Access spline_order from the type parameter — GENERICALLY, by delegating to the
# `SplineType` trait (#338 AC6). The enumerated `{T,LinearSpline}` / `{T,CubicSpline}`
# pair this replaces meant a third spline type silently had no `spline_order` and
# hit a MethodError deep inside the driver; now it inherits one from its own
# `spline_order(::MySpline)` declaration in `spline_types.jl`.
#
# `@inline` is LOAD-BEARING, not decoration. `_mf_layout_key` folds
# `spline_order(𝒮)` into the layout-cache key tuple, and that tuple is only
# stack-allocated if the whole key construction stays inlinable — the #307 gate
# asserts a cache HIT allocates zero bytes. Replacing the two literal-returning
# methods with a one-hop delegation was enough to push it over the threshold and
# cost 192 B per lookup, i.e. per matrix-free product. Measured.
@inline spline_order(𝒮::SplineIntegrator{T,S}) where {T,S<:SplineType} = spline_order(S())

# Likewise the packed-parameter block layout (#338 AC6).
@inline param_blocks(𝒮::SplineIntegrator{T,S}) where {T,S<:SplineType} = param_blocks(S())
@inline n_param_blocks(𝒮::SplineIntegrator{T,S}) where {T,S<:SplineType} =
    n_param_blocks(S())

# Role → trajectory component. DISPATCHED, so a spline type that needs a role
# beyond these declares one method beside itself; nothing in the driver branches
# on a role symbol (#338 AC6).
param_block_component(𝒮::SplineIntegrator, ::ControlValueBlock) = 𝒮.u_name
param_block_component(𝒮::SplineIntegrator, ::ControlDerivBlock) = du_name(𝒮)
param_block_component(𝒮::SplineIntegrator, ::ControlDeriv2Block) = ddu_name(𝒮)

# Helper functions to get derivative names from control name
du_name(𝒮::SplineIntegrator) = Symbol("d", 𝒮.u_name)
ddu_name(𝒮::SplineIntegrator) = Symbol("dd", 𝒮.u_name)

# Helper to compute canonical hessian knot dimension - dispatch on SplineType
canonical_hessian_knot_dim(𝒮::SplineIntegrator{T,LinearSpline}) where {T} =
    𝒮.x_dim + 𝒮.u_dim + 2  # x, u, Δt, t
canonical_hessian_knot_dim(𝒮::SplineIntegrator{T,CubicSpline}) where {T} =
    𝒮.x_dim + 2 * 𝒮.u_dim + 2  # x, u, du, Δt, t

# Interface methods for global variable support
has_global_dependence(𝒮::SplineIntegrator) = !isempty(𝒮.global_names)
global_variables(𝒮::SplineIntegrator) = 𝒮.global_names

"""
    extract_globals(𝒮::SplineIntegrator, traj::NamedTrajectory)

Extract global variable values from trajectory as a concatenated vector.
Returns `nothing` if no globals are present.
"""
function extract_globals(𝒮::SplineIntegrator, traj::NamedTrajectory)
    if 𝒮.global_dim == 0 || isempty(𝒮.global_names)
        return nothing
    end
    return vcat(
        [traj.global_data[traj.global_components[name]] for name in 𝒮.global_names]...,
    )
end

# ============================================================================ #
# Hessian structure - Shared for all types
# Uses CANONICAL ordering: [x_k, u_k, (du_k), Δt, t, x_{k+1}, u_{k+1}, (du_{k+1}), g]
# This is independent of trajectory component ordering
# 
# Creates structure for a SINGLE constraint (one knot point).
# Use get_hessian_of_lagrangian_structure() to assemble the full problem.
# ============================================================================ #

"""
    canonical_block_hessian_structure(x_dim, u_dim, spline_order, global_dim=0)
    canonical_block_hessian_structure(x_dims, u_dim, spline_order, global_dim=0)

Create sparsity structure for Hessian of Lagrangian for a single constraint block.

Returns a symmetric sparse matrix in CANONICAL ordering representing the sparsity
pattern for μ-weighted Hessian ∂²𝒮/∂z². This structure captures how state variables
couple with control parameters and time variables through the dynamics.

# Arguments
- `x_dim::Int` or `x_dims::Vector{Int}`: State dimension(s)
- `u_dim::Int`: Control dimension
- `spline_order::Int`: 1 (linear) or 3 (cubic Hermite)
- `global_dim::Int=0`: Number of global (time-invariant) variables

# Canonical Layout
Per-knot variables in order:
- States: x₁, x₂, ..., xₙ (for ensemble: all states sequential)
- Controls: u₁, u₂, ..., uₘ
- Derivatives (order 3 only): du₁, du₂, ..., duₘ
- Timestep: Δt
- Time: t

Block structure: [knot_k; knot_k+1; global_vars]

# Coupling Pattern
- States (x) couple with parameters (u, du, Δt, t, globals)
- Parameters couple with each other
- This reflects the structure of ∂²𝒮/∂xᵢ∂pⱼ and ∂²𝒮/∂pᵢ∂pⱼ
"""
function canonical_block_hessian_structure(
    x_dim::Int,
    u_dim::Int,
    spline_order::Int,
    global_dim::Int = 0,
)
    return _hessian_structure_impl([x_dim], u_dim, spline_order, global_dim)
end

# Multi-state version for ensembles - delegates to shared implementation
function canonical_block_hessian_structure(
    x_dims::Vector{Int},
    u_dim::Int,
    spline_order::Int,
    global_dim::Int = 0,
)
    return _hessian_structure_impl(x_dims, u_dim, spline_order, global_dim)
end

# Shared implementation
function _hessian_structure_impl(
    x_dims::Vector{Int},
    u_dim::Int,
    spline_order::Int,
    global_dim::Int,
)
    total_x_dim = sum(x_dims)

    if spline_order == 1
        knot_dim = total_x_dim + u_dim + 2
    elseif spline_order == 3
        knot_dim = total_x_dim + 2 * u_dim + 2
    end

    total_dim = 2 * knot_dim + global_dim
    μ∂²𝒮 = spzeros(total_dim, total_dim)

    # Canonical indices for knot k - states come first sequentially
    offset = 0
    x_comps_k_list = Vector{UnitRange{Int}}()
    for xd in x_dims
        push!(x_comps_k_list, (offset+1):(offset+xd))
        offset += xd
    end

    u_comps_k = (offset+1):(offset+u_dim)
    if spline_order == 3
        du_comps_k = (offset+u_dim+1):(offset+2*u_dim)
        Δt_comp_k = offset + 2 * u_dim + 1
        t_comp_k = offset + 2 * u_dim + 2
    else
        Δt_comp_k = offset + u_dim + 1
        t_comp_k = offset + u_dim + 2
    end

    # Canonical indices for knot k+1
    u_comps_k1 = knot_dim .+ ((offset+1):(offset+u_dim))
    if spline_order == 3
        du_comps_k1 = knot_dim .+ ((offset+u_dim+1):(offset+2*u_dim))
    end

    # Build parameter components
    if spline_order == 1
        p_comps = [collect(u_comps_k); collect(u_comps_k1); Δt_comp_k; t_comp_k]
    elseif spline_order == 3
        p_comps = [
            collect(u_comps_k);
            collect(u_comps_k1);
            collect(du_comps_k);
            collect(du_comps_k1);
            Δt_comp_k;
            t_comp_k
        ]
    end

    # Add global components if present
    if global_dim > 0
        global_comps = (2*knot_dim+1):(2*knot_dim+global_dim)
        p_comps = [p_comps; collect(global_comps)]
    end

    # μ∂ₓₖ∂ₚ𝒮 - states at knot point k couple with parameters
    for x_comps_k in x_comps_k_list
        μ∂²𝒮[x_comps_k, p_comps] .= 1.0
    end

    # μ∂ₓₖ₊₁∂ₚ𝒮 - states at knot point k+1 ALSO couple with parameters
    for x_comps_k in x_comps_k_list
        x_comps_k1 = knot_dim .+ x_comps_k
        μ∂²𝒮[x_comps_k1, p_comps] .= 1.0
    end

    # Gauss-Newton: no (p,p) blocks

    # Return full symmetric structure - triu will be taken in eval_hessian_of_lagrangian
    return sparse(Symmetric(μ∂²𝒮))
end

# ============================================================================ #
# API Methods - Shared implementations
# ============================================================================ #

"""
    get_param_indices(𝒮::SplineIntegrator, traj::NamedTrajectory, k::Int)

Build the mapping between ODE parameter indices and trajectory columns.

Returns `(traj_cols, ode_indices)` where:
- `traj_cols`: Column indices in the full Jacobian/Hessian (trajectory ordering)
- `ode_indices`: Corresponding indices in `∂ₚΦ` (ODE parameter ordering)

# ODE Parameter Layout
`[uₖ, uₖ₊₁, (duₖ, duₖ₊₁), Δt, t]` where `uₖ = [controls..., globals...]`
"""
function get_param_indices(𝒮::SplineIntegrator, traj::NamedTrajectory, k::Int)
    z_dim = traj.dim
    u_dim = 𝒮.u_dim
    control_dim = u_dim - 𝒮.global_dim
    Δt_comp = traj.components[traj.timestep][1]
    t_comp = traj.components[:t][1]
    global_offset = z_dim * traj.N

    # Build parallel arrays: trajectory columns and corresponding ODE param indices
    traj_cols = Int[]
    ode_indices = Int[]

    # ── Blocks, driven by the SplineType's declared `param_blocks` (#338 AC6) ──
    #
    # Each declared block is `u_dim` wide and contributes `control_dim` mapped
    # slots. `:value` blocks additionally map the globals (a global is constant
    # across knots, so both endpoint blocks point at the SAME global column and
    # the two contributions sum — which is exactly the dynamics-entering-global
    # chain rule). `:deriv` blocks leave their global slots unmapped, because a
    # global derivative is identically zero.
    #
    # The `if spline_order(𝒮) == 3` this replaces could not admit a third spline
    # type; the loop below needs no edit for one.
    blocks = param_blocks(𝒮)
    p_dim = 0
    for (b, (role, knot_offset)) in enumerate(blocks)
        base = (b - 1) * u_dim
        comp_name = param_block_component(𝒮, role)
        # A block whose component is absent from the trajectory (a cubic integrator
        # over a trajectory carrying no `du`) contributes no columns — matching the
        # pre-#338 `haskey(traj.components, du_sym)` guard. This is a presence check
        # on the trajectory, not a branch on the spline type or the role.
        haskey(traj.components, comp_name) || continue
        comps = traj.components[comp_name]
        append!(traj_cols, (k - 1 + knot_offset) * z_dim .+ collect(comps))
        append!(ode_indices, (base+1):(base+control_dim))
        if param_block_carries_globals(role) && 𝒮.global_dim > 0
            for (i, name) in enumerate(𝒮.global_names)
                append!(traj_cols, global_offset .+ collect(traj.global_components[name]))
                push!(ode_indices, base + control_dim + i)
            end
        end
        # `p_dim` only advances past a block that is actually present, preserving
        # the pre-#338 behaviour where a cubic integrator over a du-less trajectory
        # kept the LINEAR `p_dim = 2·u_dim` Δt/t slots.
        p_dim = base + u_dim
    end

    # Δt: ODE index p_dim + 1 (no separate global slots)
    push!(traj_cols, (k - 1) * z_dim + Δt_comp)
    push!(ode_indices, p_dim + 1)

    # t: ODE index p_dim + 2
    push!(traj_cols, (k - 1) * z_dim + t_comp)
    push!(ode_indices, p_dim + 2)

    return traj_cols, ode_indices
end

"""
    get_param_indices(𝒮::SplineIntegrator, traj::NamedTrajectory)

Return ODE parameter indices without trajectory column mapping.
Useful when trajectory columns are computed separately (e.g., MultiKetTrajectory).

Returns `(ctrl_indices, Δt_idx, t_idx, global_indices)` where all indices
are into the ODE parameter vector `[uₖ, uₖ₊₁, (duₖ, duₖ₊₁), Δt, t]`
where `uₖ = [controls..., globals...]`.

Note: ctrl_indices includes both controls and embedded globals for all knot values.
"""
function get_param_indices(𝒮::SplineIntegrator, traj::NamedTrajectory)
    u_dim = 𝒮.u_dim

    # Control indices: one `u_dim`-wide slot range per declared block (#338 AC6).
    # Each u includes [controls..., globals...]
    p_dim = n_param_blocks(𝒮) * u_dim
    ctrl_indices = collect(1:p_dim)

    # Δt and t come directly after controls (globals are embedded in u)
    Δt_idx = p_dim + 1
    t_idx = p_dim + 2

    # Global indices are empty since globals are embedded in ctrl_indices
    global_indices = Int[]

    return ctrl_indices, Δt_idx, t_idx, global_indices
end

"""
    build_ode_params(𝒮::SplineIntegrator, zₖ::KnotPoint, zₖ₊₁::KnotPoint, globals)

Construct the ODE parameter vector from knot points.

Layout: `[uₖ, uₖ₊₁, (duₖ, duₖ₊₁), Δtₖ, tₖ]`
where `uₖ = [controls..., globals...]` (globals are appended to controls)

# Arguments
- `globals::Union{Nothing, AbstractVector{<:Real}}`: Pre-extracted global values from `extract_globals`
"""
function build_ode_params(
    𝒮::SplineIntegrator,
    zₖ::KnotPoint,
    zₖ₊₁::KnotPoint,
    globals::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    p = zeros(Float64, n_param_blocks(𝒮) * 𝒮.u_dim + 2)
    return build_ode_params!(p, 𝒮, zₖ, zₖ₊₁, globals)
end

"""
    build_ode_params!(pₖ, 𝒮, zₖ, zₖ₊₁, globals=nothing)

In-place `build_ode_params`, writing into a caller-owned buffer of length
`layout.n_ode_params`. Zero allocations — the matrix-free MultiKet kernels build
the parameter vector once per knot per matvec, where the allocating form was
knot-scaled hot-path allocation (#307).

`globals` must be supplied whenever `𝒮.global_dim > 0`: the packed layout reserves
`u_dim = control_dim + global_dim` slots per `u` block, which is also what
`get_param_indices` assumes.
"""
function build_ode_params!(
    pₖ::AbstractVector{Float64},
    𝒮::SplineIntegrator,
    zₖ::KnotPoint,
    zₖ₊₁::KnotPoint,
    globals::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    u_dim = 𝒮.u_dim
    gd = 𝒮.global_dim
    gd > 0 &&
        globals === nothing &&
        throw(ArgumentError("build_ode_params!: globals required when 𝒮.global_dim > 0"))
    cd = u_dim - gd
    # One `u_dim`-wide slot range per block the SplineType declares (#338 AC6).
    # The `if spline_order(𝒮) == 1 … else cubic …` this replaces was the second of
    # three places a third spline type would have had to edit. The role's component
    # and its global handling are DISPATCHED, not branched on.
    blocks = param_blocks(𝒮)
    p_dim = 0
    @inbounds for (b, (role, knot_offset)) in enumerate(blocks)
        base = (b - 1) * u_dim
        z = knot_offset == 0 ? zₖ : zₖ₊₁
        vals = z[param_block_component(𝒮, role)]
        for j = 1:cd
            pₖ[base+j] = vals[j]
        end
        # Value blocks carry the global VALUES (constant across knots); every other
        # role's global slots are identically zero.
        gval = param_block_carries_globals(role)
        for j = 1:gd
            pₖ[base+cd+j] = gval ? globals[j] : 0.0
        end
        p_dim = base + u_dim
    end
    @inbounds begin
        pₖ[p_dim+1] = zₖ.timestep
        pₖ[p_dim+2] = haskey(zₖ.components, :t) ? zₖ[:t][1] : 0.0
    end
    return pₖ
end

"""
    get_hessian_of_lagrangian_structure(𝒮, traj)

Assemble the sparsity structure for Hessian of Lagrangian across all knot points.

Returns a sparse matrix of size (Z_dim × Z_dim) where Z_dim = z_dim*N + global_dim.
"""
function get_hessian_of_lagrangian_structure(𝒮::SplineIntegrator, traj::NamedTrajectory)
    N = traj.N
    z_dim = traj.dim
    x_dim = 𝒮.x_dim
    global_dim = traj.global_dim
    Z_dim = z_dim * N + global_dim
    μ∂²F = spzeros(Z_dim, Z_dim)

    # Get parameter indices for a representative knot point
    traj_cols, _ = get_param_indices(𝒮, traj, 1)

    # State component indices (relative to knot start)
    x_comps_rel = vcat([collect(traj.components[name]) for name in 𝒮.x_names]...)

    # For each knot point, mark the sparsity pattern
    for k = 1:(N-1)
        traj_cols_k, ode_indices_k = get_param_indices(𝒮, traj, k)

        # State indices for knot k
        x_cols_k = (k - 1) * z_dim .+ x_comps_rel

        # (x,p) cross-terms: state-parameter coupling
        for xi in x_cols_k
            for pj in traj_cols_k
                if xi <= pj
                    μ∂²F[xi, pj] = 1.0
                else
                    μ∂²F[pj, xi] = 1.0
                end
            end
        end

        # (p,p) blocks from exact Hessian
        if 𝒮.exact_hessian && !isnothing(𝒮.hess_active_pairs)
            # Build reverse lookup: ODE param index → trajectory column
            ode_to_traj = Dict{Int,Int}()
            for (tc, oi) in zip(traj_cols_k, ode_indices_k)
                ode_to_traj[oi] = tc
            end

            for (ode_i, ode_j) in 𝒮.hess_active_pairs
                haskey(ode_to_traj, ode_i) || continue
                haskey(ode_to_traj, ode_j) || continue
                traj_i = ode_to_traj[ode_i]
                traj_j = ode_to_traj[ode_j]
                ri, rj = minmax(traj_i, traj_j)
                μ∂²F[ri, rj] = 1.0
            end
        end
    end

    return μ∂²F
end

"""
    fill_state_jacobian_structure!(∂F, rows, x_cols_k, x_cols_k1, 𝒮, traj, k)

Fill the state-dependent Jacobian structure for knot point k.
Dispatches on SplineIntegrator type to exploit trajectory-specific sparsity.
"""
function fill_state_jacobian_structure! end

# UnitaryTrajectory: block-diagonal structure (each unitary column evolves independently)
function fill_state_jacobian_structure!(
    ∂F::SparseMatrixCSC,
    rows::AbstractVector{Int},
    z_dim::Int,
    k::Int,
    𝒮::SplineIntegrator{UnitaryTrajectory},
    traj::NamedTrajectory,
)
    ketdim = 𝒮.ketdim
    block_size = 2 * ketdim
    x_comps_rel = collect(traj.components[𝒮.x_names[1]])
    x_cols_k = (k - 1) * z_dim .+ x_comps_rel
    x_cols_k1 = k * z_dim .+ x_comps_rel

    for i = 1:ketdim
        block_rows = rows[slice(i, block_size)]
        block_x_cols_k = x_cols_k[slice(i, block_size)]
        block_x_cols_k1 = x_cols_k1[slice(i, block_size)]

        # ∂xₖ𝒮: block-diagonal (-Φₖ blocks)
        for r in block_rows, c in block_x_cols_k
            ∂F[r, c] = 1.0
        end

        # ∂xₖ₊₁𝒮: block-diagonal (identity blocks)
        for r in block_rows, c in block_x_cols_k1
            ∂F[r, c] = 1.0
        end
    end
end

# MultiKetTrajectory: block-diagonal structure (each ket evolves independently)
function fill_state_jacobian_structure!(
    ∂F::SparseMatrixCSC,
    rows::AbstractVector{Int},
    z_dim::Int,
    k::Int,
    𝒮::SplineIntegrator{MultiKetTrajectory},
    traj::NamedTrajectory,
)
    ketdim = 𝒮.ketdim
    block_size = 2 * ketdim

    for (ket_idx, x_name) in enumerate(𝒮.x_names)
        x_comps_rel = collect(traj.components[x_name])
        x_cols_k = (k - 1) * z_dim .+ x_comps_rel
        x_cols_k1 = k * z_dim .+ x_comps_rel

        # Rows for this ket
        ket_rows = rows[slice(ket_idx, block_size)]

        # ∂xₖ𝒮: dense block for this ket
        for r in ket_rows, c in x_cols_k
            ∂F[r, c] = 1.0
        end

        # ∂xₖ₊₁𝒮: identity block for this ket
        for r in ket_rows, c in x_cols_k1
            ∂F[r, c] = 1.0
        end
    end
end

# Default (KetTrajectory, etc.): dense state blocks
function fill_state_jacobian_structure!(
    ∂F::SparseMatrixCSC,
    rows::AbstractVector{Int},
    z_dim::Int,
    k::Int,
    𝒮::SplineIntegrator,
    traj::NamedTrajectory,
)
    x_comps_rel = vcat([collect(traj.components[name]) for name in 𝒮.x_names]...)
    x_cols_k = (k - 1) * z_dim .+ x_comps_rel
    x_cols_k1 = k * z_dim .+ x_comps_rel
    ∂F[rows, x_cols_k] .= 1.0
    ∂F[rows, x_cols_k1] .= 1.0
end

"""
    get_jacobian_structure(𝒮::SplineIntegrator, traj::NamedTrajectory)

Get the sparsity structure for the constraint Jacobian across all knot points.

For UnitaryTrajectory: exploits block-diagonal structure - each column of the unitary
evolves independently, so ∂F/∂xₖ is block-diagonal (ketdim blocks of size 2ketdim×2ketdim).

For MultiKetTrajectory: each ket evolves independently, so ∂F/∂xₖ is block-diagonal
with n_kets blocks of size 2ketdim×2ketdim.
"""
function get_jacobian_structure(𝒮::SplineIntegrator, traj::NamedTrajectory)
    N = traj.N
    x_dim = 𝒮.x_dim
    z_dim = traj.dim
    global_dim = traj.global_dim
    F_dim = x_dim * (N - 1)
    Z_dim = z_dim * N + global_dim
    ∂F = spzeros(F_dim, Z_dim)

    for k = 1:(N-1)
        rows = slice(k, x_dim)

        # Fill state Jacobian structure (dispatches on trajectory type)
        fill_state_jacobian_structure!(∂F, rows, z_dim, k, 𝒮, traj)

        # ∂p𝒮: Parameter columns (same for all trajectory types)
        traj_cols, _ = get_param_indices(𝒮, traj, k)
        for col in traj_cols
            ∂F[rows, col] .= 1.0
        end
    end

    return ∂F
end

@views function evaluate!(δ::AbstractVector, 𝒮::SplineIntegrator, traj::NamedTrajectory)
    x_dim = 𝒮.x_dim
    globals = extract_globals(𝒮, traj)
    Threads.@threads for k = 1:(traj.N-1)
        δₖ = δ[slice(k, x_dim)]
        𝒮(δₖ, traj[k], traj[k+1], k, globals)
    end
    return nothing
end

"""
    eval_jacobian(𝒮::SplineIntegrator, traj::NamedTrajectory)

Evaluate the constraint Jacobian by:
1. Computing ODE derivatives in parallel (fills DiffResults)
2. Assembling the full Jacobian directly from DiffResults

The Jacobian for constraint k is:
- ∂𝒮/∂ψₖ = -Φₖ  (or column-wise for unitary)
- ∂𝒮/∂ψₖ₊₁ = I
- ∂𝒮/∂p = -(∂ₚΦₖ) ψₖ
"""
@views function eval_jacobian(𝒮::SplineIntegrator, traj::NamedTrajectory)
    N = traj.N
    x_dim = 𝒮.x_dim
    z_dim = traj.dim
    F_dim = x_dim * (N - 1)
    Z_dim = z_dim * N + traj.global_dim

    # Extract globals once
    globals = extract_globals(𝒮, traj)

    # 1. Compute ODE Jacobians in parallel (fills ode_jac_results)
    Threads.@threads for k = 1:(N-1)
        compute_ode_jacobian!(𝒮, traj[k], traj[k+1], k, globals)
    end

    # 2. Assemble full Jacobian from PropagatorResults
    ∂F = spzeros(F_dim, Z_dim)

    # State component indices (relative to knot start)
    x_comps_rel = vcat([collect(traj.components[name]) for name in 𝒮.x_names]...)

    # Preallocate temp vectors for matrix-vector products
    temp_vec = zeros(x_dim)
    tmp_complex = Vector{ComplexF64}(undef, 𝒮.ketdim)
    ψ_complex = Vector{ComplexF64}(undef, 𝒮.ketdim)

    @inbounds for k = 1:(N-1)
        rows = slice(k, x_dim)

        # Extract propagator and sensitivities from PropagatorResult
        pdim = propagator_side_dim(𝒮)
        Φₖ = get_propagator(𝒮.prop_results[k], pdim)
        ∂ₚΦₖ = get_sensitivities(𝒮.prop_results[k], pdim)

        # Get state at knot k
        ψₖ = get_state_vector(𝒮, traj[k])

        # ∂xₖ𝒮 = -Φₖ (or block-wise for unitary)
        x_cols_k = (k - 1) * z_dim .+ x_comps_rel
        fill_state_jacobian!(∂F, rows, x_cols_k, Φₖ, 𝒮)

        # ∂xₖ₊₁𝒮 = I
        x_cols_k1 = k * z_dim .+ x_comps_rel
        for (i, row) in enumerate(rows)
            ∂F[row, x_cols_k1[i]] = 1.0
        end

        # ∂p𝒮 = -(∂ₚΦₖ) ψₖ
        traj_cols, ode_indices = get_param_indices(𝒮, traj, k)
        for (col, ode_idx) in zip(traj_cols, ode_indices)
            fill_param_jacobian!(
                ∂F,
                rows,
                col,
                ∂ₚΦₖ,
                ode_idx,
                ψₖ,
                temp_vec,
                tmp_complex,
                ψ_complex,
                𝒮,
            )
        end
    end

    return ∂F
end

# ─────────────────────────────────────────────────────────────────────────── #
# Per-operation precompute-needed predicates.
#
# The trajectory-level `eval_constraint_*!` wrappers call `compute_ode_jacobian!`
# per knot before the per-step accumulate loop **only when the per-step body
# requires `prop_results`**.
#
# CORRECTED BY ADR-0009, ENFORCED BY #340: self-containedness is a property of the
# CELL that must be DEMONSTRATED COLD, PER PRODUCT, never inferred from the
# alg-dispatch design. The Chebyshev and Magnus Ket cells are genuinely
# self-contained. For the Tsit5Alg Ket cell, #340 measured each product cold — no
# preceding `eval_jacobian` / `compute_ode_jacobian!`, no warming call before the
# assertion — against the same product after a refresh, and against a central
# finite difference of the integrator's OWN constraint residual:
#
#   JVP  cold-vs-FD 6.51e-1,  cold-vs-warm 6.51e-1   ⇒ NOT self-contained ⇒ true
#   VJP                       cold-vs-warm 7.77e-1   ⇒ NOT self-contained ⇒ true
#   HVP                       cold-vs-warm 0.0       ⇒     self-contained ⇒ false
#
# Root cause of the first two, also found by #340: the `Tsit5Alg`-specialised
# `step_jvp!` / `step_vjp!` in `spline_integrator_ket.jl` — the bodies that DO
# integrate their own ket-level ODEs — are UNREACHABLE. Their signatures use
# unbounded `where {S,R,D}`, which Julia does not rank above the propagator-level
# `SplineIntegrator{KetTrajectory}` method (the bug-#247 shadowing pattern these
# predicates carry bounded typevars to avoid), so a Ket Tsit5 JVP/VJP runs the
# propagator-level body and reads `prop_results`. The old `false` returns were read
# off the shadowed bodies rather than off the dispatched ones. `step_hvp!` has no
# propagator-level Ket method to be shadowed by, so its Tsit5 body does run — and
# it reconstructs `ψ(1)` from its own forward solve, which is why its `false` is
# correct and is now pinned by a cold test rather than by an inference.
#
# A cell that is not self-contained is made correct by the per-iterate sensitivity
# cache (`build_ket_sens_cache`), NOT by a predicate: with a cache the refresh is
# once per accepted iterate, without one it is once per product.
#
# Skipping the precompute eliminates one O((1+n_p)·d²)-state sensitivity ODE
# solve per knot per call, which is the dominant cost when matrix-free HVP
# is invoked many times per outer step (NewtonCG inner CG iterations).
# ─────────────────────────────────────────────────────────────────────────── #
"""
    eval_hessian_of_lagrangian(𝒮::SplineIntegrator, traj::NamedTrajectory, μ::AbstractVector)

Evaluate the Lagrangian-weighted Hessian by:
1. Computing ODE Hessians in parallel (fills ode_hess_results)
2. Assembling directly from DiffResults

The Hessian blocks are:
- ∂²𝒮/∂ψₖ∂p = -(∂ₚΦₖ)ᵀ μ
- ∂²𝒮/∂pᵢ∂pⱼ = -μᵀ (∂²ₚᵢₚⱼΦₖ) ψₖ

Uses the Gauss-Newton approximation:
- Reuses already-computed `ode_jac_results` for cross-terms (∂ψₖ∂p)
- Drops (p,p) blocks (second-order ODE terms)
"""
function eval_hessian_of_lagrangian(
    𝒮::SplineIntegrator,
    traj::NamedTrajectory,
    μ::AbstractVector,
)
    N = traj.N
    x_dim = 𝒮.x_dim
    z_dim = traj.dim
    global_dim = traj.global_dim
    Z_dim = z_dim * N + global_dim

    # Assemble Hessian from prop_results
    # Note: In Ipopt, eval_jac_g is always called before eval_h, so prop_results is populated
    μ∂²F = spzeros(Z_dim, Z_dim)

    # State component indices (relative to knot start)
    x_comps_rel = vcat([collect(traj.components[name]) for name in 𝒮.x_names]...)

    # Pre-compute (in parallel) the per-knot second-order sensitivity ODE
    # solves before the sequential sparse-matrix accumulation pass. The
    # ODE solves are the dominant cost; sparse-matrix writes need to be
    # sequential because adjacent knots' globals share columns. Disjoint
    # writes on `𝒮.prop_results[k]` make the precompute pass thread-safe.
    use_exact = 𝒮.exact_hessian && !isnothing(𝒮.hess_active_pairs)
    T_mats_per_knot = if use_exact
        _T_mats = Vector{Vector{Matrix{ComplexF64}}}(undef, N - 1)
        Threads.@threads for k = 1:(N-1)
            _T_mats[k] = _compute_ode_hessian_from_traj(𝒮, traj, k)
        end
        _T_mats
    else
        nothing
    end

    # Preallocate temp vectors (sequential pass — single set is fine)
    temp_vec = zeros(x_dim)
    tmp1 = Vector{ComplexF64}(undef, 𝒮.ketdim)
    tmp2 = Vector{ComplexF64}(undef, 𝒮.ketdim)

    @inbounds for k = 1:(N-1)
        μₖ = μ[slice(k, x_dim)]

        # Get parameter mapping
        traj_cols, ode_indices = get_param_indices(𝒮, traj, k)

        # State indices for knot k
        x_cols_k = (k - 1) * z_dim .+ x_comps_rel

        # Extract ∂ₚΦₖ from prop_results (already computed for Jacobian or by
        # the parallel `_compute_ode_hessian_from_traj` precompute above)
        pdim = propagator_side_dim(𝒮)
        ∂ₚΦₖ = get_sensitivities(𝒮.prop_results[k], pdim)

        # (x,p) cross-terms: -(∂ₚΦₖ)ᵀ μ
        for (col, ode_idx) in zip(traj_cols, ode_indices)
            fill_state_param_hessian!(
                μ∂²F,
                x_cols_k,
                col,
                ∂ₚΦₖ,
                ode_idx,
                μₖ,
                temp_vec,
                tmp1,
                tmp2,
                𝒮,
            )
        end

        # (p,p) blocks from exact Hessian (second-order sensitivities)
        if use_exact
            ψₖ = traj.data[x_comps_rel, k]

            # Build reverse lookup: ODE param index → trajectory column
            # Note: globals map to the same traj_col from both uₖ and uₖ₊₁ slots
            ode_to_traj = Dict{Int,Int}()
            for (tc, oi) in zip(traj_cols, ode_indices)
                ode_to_traj[oi] = tc
            end

            T_mats = T_mats_per_knot[k]

            for (pair_k, (ode_i, ode_j)) in enumerate(𝒮.hess_active_pairs)
                # Map ODE indices to trajectory columns
                haskey(ode_to_traj, ode_i) || continue
                haskey(ode_to_traj, ode_j) || continue
                traj_i = ode_to_traj[ode_i]
                traj_j = ode_to_traj[ode_j]

                # Compute -real(μ_C' * T_{ij} * ψ_C)
                val = hessian_pp_contraction(T_mats[pair_k], μₖ, ψₖ, tmp1, tmp2)

                # Store in upper triangle
                ri, rj = minmax(traj_i, traj_j)
                μ∂²F[ri, rj] += val
            end
        end
    end

    return μ∂²F
end

"""
    _compute_ode_hessian_from_traj(𝒮::SplineIntegrator, traj, k)

Solve the second-order sensitivity ODE for knot point k using trajectory data.
Builds ODE parameters from the trajectory and calls the analytical solver.
Also re-populates `𝒮.prop_results[k]` with Φ and first-order sensitivities.

Returns a vector of T_{ij} matrices for the active pairs.
"""
function _compute_ode_hessian_from_traj(𝒮::SplineIntegrator, traj::NamedTrajectory, k::Int)
    # Build ODE parameters from trajectory data (same as eval_jacobian path)
    zₖ = traj[k]
    zₖ₊₁ = traj[k+1]
    globals = extract_globals(𝒮, traj)
    pₖ = build_ode_params(𝒮, zₖ, zₖ₊₁, globals)
    return _compute_ode_hessian_analytical!(𝒮, pₖ, k)
end

# ============================================================================ #
# Dispatch helpers - implemented per trajectory type
# ============================================================================ #

"""
    compute_ode_jacobian!(𝒮, zₖ, zₖ₊₁, k, globals)

Compute ODE Jacobian and store in `𝒮.ode_jac_results[k]`.
Implemented per trajectory type.

# Arguments
- `globals::Union{Nothing, AbstractVector{<:Real}}`: Pre-extracted global values
"""
function compute_ode_jacobian! end

"""
    compute_ode_hessian!(𝒮, zₖ, zₖ₊₁, k, globals)

Compute ODE Hessian and store in `𝒮.ode_hess_results[k]`.
Implemented per trajectory type.

# Arguments
- `globals::Union{Nothing, AbstractVector{<:Real}}`: Pre-extracted global values
"""
function compute_ode_hessian! end

"""
    get_state_vector(𝒮, zₖ)

Extract the state vector(s) from a knot point for matrix-vector products.
"""
function get_state_vector end

"""
    fill_state_jacobian!(∂F, rows, x_cols, Φₖ, 𝒮)

Fill the ∂𝒮/∂ψₖ block of the Jacobian.
Default dispatches on element type of Φₖ:
- Real (density): direct -Φₖ
- Complex (ket): -propagator_to_iso(Φₖ)
Overridden for Unitary (column-wise blocks).
"""
function fill_state_jacobian!(
    ∂F,
    rows,
    x_cols,
    Φₖ::AbstractMatrix{<:Real},
    𝒮::SplineIntegrator,
)
    # Real propagator (density): direct assignment
    ∂F[rows, x_cols] = -Φₖ
end

function fill_state_jacobian!(
    ∂F,
    rows,
    x_cols,
    Φₖ::AbstractMatrix{<:Complex},
    𝒮::SplineIntegrator,
)
    # Complex propagator (ket): convert to real iso form
    ∂F[rows, x_cols] = -propagator_to_iso(Φₖ)
end

"""
    fill_param_jacobian!(∂F, rows, col, ∂ₚΦₖ, ode_idx, ψₖ, temp_vec, tmp_complex, ψ_complex, 𝒮)

Fill a single parameter column of the Jacobian: -(∂ₚΦₖ[:,:,ode_idx]) ψₖ
Dispatches on element type: real (density) vs complex (ket).

# Pre-allocated buffers (used by complex path, ignored by real path)
- `tmp_complex`: complex n-vector workspace
- `ψ_complex`: complex n-vector for iso→complex state conversion
"""
function fill_param_jacobian!(
    ∂F,
    rows,
    col,
    ∂ₚΦₖ::AbstractArray{<:Real,3},
    ode_idx,
    ψₖ,
    temp_vec,
    tmp_complex,
    ψ_complex,
    𝒮::SplineIntegrator,
)
    # Real propagator (density): direct matrix-vector product
    mul!(temp_vec, @view(∂ₚΦₖ[:, :, ode_idx]), ψₖ)
    for (i, r) in enumerate(rows)
        ∂F[r, col] += -temp_vec[i]
    end
end

function fill_param_jacobian!(
    ∂F,
    rows,
    col,
    ∂ₚΦₖ::AbstractArray{<:Complex,3},
    ode_idx,
    ψₖ,
    temp_vec,
    tmp_complex,
    ψ_complex,
    𝒮::SplineIntegrator,
)
    # Complex propagator (ket): convert iso state → complex, apply sensitivity
    Sⱼ = @view ∂ₚΦₖ[:, :, ode_idx]
    iso_to_complex_ket!(ψ_complex, ψₖ)
    sensitivity_to_jac_col!(temp_vec, Sⱼ, ψ_complex, tmp_complex)
    # sensitivity_to_jac_col! computes iso(-Sⱼ·ψ), already includes negative sign
    for (i, r) in enumerate(rows)
        ∂F[r, col] += temp_vec[i]
    end
end

"""
    fill_state_param_hessian!(μ∂²F, x_cols, col, ∂ₚΦₖ, ode_idx, μₖ, temp_vec, tmp1, tmp2, 𝒮)

Fill the (x_k, p) block of the Hessian: -(∂ₚΦₖ[:,:,ode_idx])ᵀ μₖ
Dispatches on element type: real (density) vs complex (ket).

# Pre-allocated buffers (used by complex path, ignored by real path)
- `tmp1`: complex n-vector workspace for iso→complex conversion
- `tmp2`: complex n-vector workspace for mul! result
"""
function fill_state_param_hessian!(
    μ∂²F,
    x_cols,
    col,
    ∂ₚΦₖ::AbstractArray{<:Real,3},
    ode_idx,
    μₖ,
    temp_vec,
    tmp1,
    tmp2,
    𝒮::SplineIntegrator,
)
    # Real propagator (density): direct -(∂ₚΦ)ᵀ μ
    mul!(temp_vec, @view(∂ₚΦₖ[:, :, ode_idx])', μₖ)
    for (i, xi) in enumerate(x_cols)
        if xi <= col
            μ∂²F[xi, col] += -temp_vec[i]
        else
            μ∂²F[col, xi] += -temp_vec[i]
        end
    end
end

function fill_state_param_hessian!(
    μ∂²F,
    x_cols,
    col,
    ∂ₚΦₖ::AbstractArray{<:Complex,3},
    ode_idx,
    μₖ,
    temp_vec,
    tmp1,
    tmp2,
    𝒮::SplineIntegrator,
)
    # Complex propagator (ket): convert via sensitivity_to_hess_col! with pre-allocated buffers
    Sⱼ = @view ∂ₚΦₖ[:, :, ode_idx]
    sensitivity_to_hess_col!(temp_vec, Sⱼ, μₖ, tmp1, tmp2)
    # sensitivity_to_hess_col! computes iso(-Sⱼ'·μ_C), already includes negative sign
    for (i, xi) in enumerate(x_cols)
        if xi <= col
            μ∂²F[xi, col] += temp_vec[i]
        else
            μ∂²F[col, xi] += temp_vec[i]
        end
    end
end

# ============================================================================ #
# Shared Sensitivity ODE Builder
# ============================================================================ #

"""
    build_sensitivity_ode(drift_op, drives, u_dim, ketdim, order)

Build sensitivity-augmented ODE from `AbstractDrive`s (generalized version).

This method supports both linear and nonlinear control maps. Each drive
provides its own coefficient function and Jacobian, enabling control mixing
(e.g., `u₁*u₂`) and other nonlinear functions of controls.

# Arguments
- `drift_op::AbstractDynamicsOperator`: Drift Hamiltonian operator
- `drives::AbstractVector{<:AbstractDrive}`: Drive terms with coefficients
- `u_dim::Int`: Control dimension (includes globals)
- `ketdim::Int`: Ket space dimension
- `order::Int`: Spline order (1=linear, 3=cubic)
"""
function build_sensitivity_ode(
    drift_op::AbstractDynamicsOperator,
    drives::AbstractVector{<:AbstractDrive},
    u_dim::Int,
    ketdim::Int,
    order::Int,
)
    Φ_dim = ketdim^2
    n_terms = length(drives)

    if order == 1
        n_params = 2 * u_dim + 2  # [uₖ, uₖ₊₁, Δt, t]
    elseif order == 3
        n_params = 4 * u_dim + 2  # [uₖ, uₖ₊₁, duₖ, duₖ₊₁, Δt, t]
    else
        error("Unsupported spline order: $order")
    end

    if isempty(drives)
        return nothing, n_params
    end

    # Pre-build operator wrappers from drive Hamiltonians
    drive_ops = [_to_operator(d.H) for d in drives]

    # Precompute control→drives mapping for structural sparsity.
    # control_to_drives[j] lists drive term indices whose Jacobian ∂c/∂u_j can be nonzero.
    # When active_controls(d) is empty, the drive depends on all controls.
    control_to_drives = [Int[] for _ = 1:u_dim]
    drives_active = [active_controls(d) for d in drives]
    for (t_idx, ac) in enumerate(drives_active)
        if isempty(ac)
            for j = 1:u_dim
                push!(control_to_drives[j], t_idx)
            end
        else
            for j in ac
                push!(control_to_drives[j], t_idx)
            end
        end
    end

    make_f!_sens = if order == 1
        () -> begin
            u_interp = zeros(u_dim)
            (dx, x, p, τ) -> begin
                uₖ = @view p[1:u_dim]
                uₖ₊₁ = @view p[(u_dim+1):2u_dim]
                Δtₖ = p[end-1]

                onemτ = 1 - τ
                @inbounds for i = 1:u_dim
                    u_interp[i] = onemτ * uₖ[i] + τ * uₖ₊₁[i]
                end

                # Apply H(u(τ)) to a matrix: drift + Σ coeff(d, u_τ) * d_op
                @inline function apply_H!(dM, M)
                    apply!(dM, drift_op, M, -im * Δtₖ, false)
                    @inbounds for t_idx = 1:n_terms
                        c = drive_coeff(drives[t_idx], u_interp)
                        apply!(dM, drive_ops[t_idx], M, -im * Δtₖ * c, true)
                    end
                end

                # dΦ/dτ = -i·Δt·H·Φ
                Φ_vec = @view x[1:Φ_dim]
                Φ_mat = reshape(Φ_vec, ketdim, ketdim)
                dΦ_vec = @view dx[1:Φ_dim]
                dΦ_mat = reshape(dΦ_vec, ketdim, ketdim)
                apply_H!(dΦ_mat, Φ_mat)

                # Sensitivities for uₖ[j]: spline basis = (1-τ)
                @inbounds for j = 1:u_dim
                    offset = Φ_dim + (j - 1) * Φ_dim
                    Sⱼ_mat = reshape(@view(x[(offset+1):(offset+Φ_dim)]), ketdim, ketdim)
                    dSⱼ_mat = reshape(@view(dx[(offset+1):(offset+Φ_dim)]), ketdim, ketdim)
                    apply_H!(dSⱼ_mat, Sⱼ_mat)
                    # Forcing: Σ_drives (∂coeff/∂u_j) * (1-τ) * drive_op * Φ
                    for t_idx in control_to_drives[j]
                        dc = drive_coeff_jac(drives[t_idx], u_interp, j)
                        if dc != 0.0
                            apply!(
                                dSⱼ_mat,
                                drive_ops[t_idx],
                                Φ_mat,
                                -im * Δtₖ * onemτ * dc,
                                true,
                            )
                        end
                    end
                end

                # Sensitivities for uₖ₊₁[j]: spline basis = τ
                @inbounds for j = 1:u_dim
                    offset = Φ_dim + (u_dim + j - 1) * Φ_dim
                    Sⱼ_mat = reshape(@view(x[(offset+1):(offset+Φ_dim)]), ketdim, ketdim)
                    dSⱼ_mat = reshape(@view(dx[(offset+1):(offset+Φ_dim)]), ketdim, ketdim)
                    apply_H!(dSⱼ_mat, Sⱼ_mat)
                    for t_idx in control_to_drives[j]
                        dc = drive_coeff_jac(drives[t_idx], u_interp, j)
                        if dc != 0.0
                            apply!(dSⱼ_mat, drive_ops[t_idx], Φ_mat, -im * Δtₖ * τ * dc, true)
                        end
                    end
                end

                # Sensitivity for Δt: ∂(-iΔtH(u))/∂Δt = -iH(u)
                let j = 2 * u_dim + 1
                    offset = Φ_dim + (j - 1) * Φ_dim
                    Sⱼ_mat = reshape(@view(x[(offset+1):(offset+Φ_dim)]), ketdim, ketdim)
                    dSⱼ_mat = reshape(@view(dx[(offset+1):(offset+Φ_dim)]), ketdim, ketdim)
                    apply_H!(dSⱼ_mat, Sⱼ_mat)
                    apply!(dSⱼ_mat, drift_op, Φ_mat, -im, true)
                    @inbounds for t_idx = 1:n_terms
                        c = drive_coeff(drives[t_idx], u_interp)
                        apply!(dSⱼ_mat, drive_ops[t_idx], Φ_mat, -im * c, true)
                    end
                end

                # Sensitivity for t (zero forcing for time-independent)
                let j = 2 * u_dim + 2
                    offset = Φ_dim + (j - 1) * Φ_dim
                    Sⱼ_mat = reshape(@view(x[(offset+1):(offset+Φ_dim)]), ketdim, ketdim)
                    dSⱼ_mat = reshape(@view(dx[(offset+1):(offset+Φ_dim)]), ketdim, ketdim)
                    apply_H!(dSⱼ_mat, Sⱼ_mat)
                end

                return nothing
            end
        end

    else  # order == 3
        () -> begin
            u_interp = zeros(u_dim)
            jac_cache = zeros(u_dim, n_terms)
            (dx, x, p, τ) -> begin
                uₖ = @view p[1:u_dim]
                uₖ₊₁ = @view p[(u_dim+1):2u_dim]
                duₖ = @view p[(2u_dim+1):3u_dim]
                duₖ₊₁ = @view p[(3u_dim+1):4u_dim]
                Δtₖ = p[end-1]

                τ2 = τ * τ
                τ3 = τ2 * τ
                h00 = 2τ3 - 3τ2 + 1
                h10 = (τ3 - 2τ2 + τ) * Δtₖ
                h01 = -2τ3 + 3τ2
                h11 = (τ3 - τ2) * Δtₖ
                dh10_dΔt = τ3 - 2τ2 + τ
                dh11_dΔt = τ3 - τ2

                # Interpolate full control vector
                @inbounds for i = 1:u_dim
                    u_interp[i] =
                        h00 * uₖ[i] + h10 * duₖ[i] + h01 * uₖ₊₁[i] + h11 * duₖ₊₁[i]
                end

                # Op-fusion: per-substep jac coefficient cache.
                # jac_cache[j, t_idx] = ∂c_t/∂u_j; non-active (j, t_idx) entries left at 0.
                @inbounds for t_idx = 1:n_terms
                    ac = drives_active[t_idx]
                    iter = isempty(ac) ? (1:u_dim) : ac
                    for j = 1:u_dim
                        jac_cache[j, t_idx] = 0.0
                    end
                    for j in iter
                        jac_cache[j, t_idx] = drive_coeff_jac(drives[t_idx], u_interp, j)
                    end
                end

                # Op-fusion: wide-matrix homogeneous A·block.
                # Augmented state x = [Φ; S₁; ...; S_{n_params}] is contiguous in memory.
                # Reshape as ketdim × (ketdim·(1+n_params)) and apply A once per drive_op.
                wide_state = reshape(x, ketdim, ketdim * (1 + n_params))
                wide_dstate = reshape(dx, ketdim, ketdim * (1 + n_params))
                apply!(wide_dstate, drift_op, wide_state, -im * Δtₖ, false)
                @inbounds for t_idx = 1:n_terms
                    c = drive_coeff(drives[t_idx], u_interp)
                    apply!(wide_dstate, drive_ops[t_idx], wide_state, -im * Δtₖ * c, true)
                end

                # Per-block views for inhomogeneous forcing
                Φ_vec = @view x[1:Φ_dim]
                Φ_mat = reshape(Φ_vec, ketdim, ketdim)

                # uₖ sensitivities: spline basis = h00
                @inbounds for j = 1:u_dim
                    offset = Φ_dim + (j - 1) * Φ_dim
                    dSⱼ_mat = reshape(@view(dx[(offset+1):(offset+Φ_dim)]), ketdim, ketdim)
                    for t_idx in control_to_drives[j]
                        dc = jac_cache[j, t_idx]
                        if dc != 0.0
                            apply!(dSⱼ_mat, drive_ops[t_idx], Φ_mat, -im * Δtₖ * h00 * dc, true)
                        end
                    end
                end

                # uₖ₊₁ sensitivities: spline basis = h01
                @inbounds for j = 1:u_dim
                    offset = Φ_dim + (u_dim + j - 1) * Φ_dim
                    dSⱼ_mat = reshape(@view(dx[(offset+1):(offset+Φ_dim)]), ketdim, ketdim)
                    for t_idx in control_to_drives[j]
                        dc = jac_cache[j, t_idx]
                        if dc != 0.0
                            apply!(dSⱼ_mat, drive_ops[t_idx], Φ_mat, -im * Δtₖ * h01 * dc, true)
                        end
                    end
                end

                # duₖ sensitivities: spline basis = h10
                @inbounds for j = 1:u_dim
                    offset = Φ_dim + (2u_dim + j - 1) * Φ_dim
                    dSⱼ_mat = reshape(@view(dx[(offset+1):(offset+Φ_dim)]), ketdim, ketdim)
                    for t_idx in control_to_drives[j]
                        dc = jac_cache[j, t_idx]
                        if dc != 0.0
                            apply!(dSⱼ_mat, drive_ops[t_idx], Φ_mat, -im * Δtₖ * h10 * dc, true)
                        end
                    end
                end

                # duₖ₊₁ sensitivities: spline basis = h11
                @inbounds for j = 1:u_dim
                    offset = Φ_dim + (3u_dim + j - 1) * Φ_dim
                    dSⱼ_mat = reshape(@view(dx[(offset+1):(offset+Φ_dim)]), ketdim, ketdim)
                    for t_idx in control_to_drives[j]
                        dc = jac_cache[j, t_idx]
                        if dc != 0.0
                            apply!(dSⱼ_mat, drive_ops[t_idx], Φ_mat, -im * Δtₖ * h11 * dc, true)
                        end
                    end
                end

                # Δt sensitivity: homogeneous done by wide-matrix; add inhomogeneous forcing.
                let j = 4 * u_dim + 1
                    offset = Φ_dim + (j - 1) * Φ_dim
                    dSⱼ_mat = reshape(@view(dx[(offset+1):(offset+Φ_dim)]), ketdim, ketdim)
                    # -iH(u) contribution
                    apply!(dSⱼ_mat, drift_op, Φ_mat, -im, true)
                    @inbounds for t_idx = 1:n_terms
                        c = drive_coeff(drives[t_idx], u_interp)
                        apply!(dSⱼ_mat, drive_ops[t_idx], Φ_mat, -im * c, true)
                    end
                    # -iΔt · Σ_drives coeff_jac(d, u, i) · ∂u_i/∂Δt · drive_op · Φ
                    @inbounds for t_idx = 1:n_terms
                        op = drive_ops[t_idx]
                        ac = drives_active[t_idx]
                        iter = isempty(ac) ? (1:u_dim) : ac
                        total_du_dΔt = 0.0
                        for i in iter
                            dc = jac_cache[i, t_idx]
                            if dc != 0.0
                                du_i_dΔt = dh10_dΔt * duₖ[i] + dh11_dΔt * duₖ₊₁[i]
                                total_du_dΔt += dc * du_i_dΔt
                            end
                        end
                        if total_du_dΔt != 0.0
                            apply!(dSⱼ_mat, op, Φ_mat, -im * Δtₖ * total_du_dΔt, true)
                        end
                    end
                end

                # t sensitivity: homogeneous done by wide-matrix above; zero forcing.

                return nothing
            end
        end
    end

    return make_f!_sens, n_params
end

"""
    build_ket_jvp_ode(drift_op, drives, u_dim, ketdim, order)
        -> (make_f!_jvp, n_params)

Build an in-place ODE RHS for the **ket-level forward tangent** used by Phase 2's
matrix-free JVP. State layout: `x = [ψ; Y] ∈ ℂ^{2d}` where `d = ketdim`.

Dynamics:
    dψ/dτ = -i·Δt·H(τ; u(τ)) · ψ
    dY/dτ = -i·Δt·H(τ; u(τ)) · Y  -  i·Δt·M_v(τ; u(τ), v_p) · ψ

where `M_v(τ; ·, v_p) = Σⱼ (∂_{p_j}H)(τ; ·) · v_p[j]` is the v_p-contracted
parameter-derivative operator. `v_p` is appended at the end of the ODE
parameter vector after the standard `[uₖ, uₖ₊₁, (duₖ, duₖ₊₁), Δt, t]` block.

This contracts the work that `build_sensitivity_ode` does in parallel for `n_p`
sensitivities into a single forcing term, at the cost of having to rebuild and
re-solve the ODE per JVP-direction `v` (vs amortizing one big solve). The per-
solve state size drops from `(1 + n_p)·d²` complex to `2d` complex — the goal
of Phase 2.

Returns `(make_f!_jvp, n_params)` matching `build_sensitivity_ode`'s shape.
The returned `make_f!_jvp()` returns a closure suitable for `ODEProblem`.
"""
function build_ket_jvp_ode(
    drift_op::AbstractDynamicsOperator,
    drives::AbstractVector{<:AbstractDrive},
    u_dim::Int,
    ketdim::Int,
    order::Int,
)
    n_terms = length(drives)

    if order == 1
        n_params = 2 * u_dim + 2  # [uₖ, uₖ₊₁, Δt, t]
    elseif order == 3
        n_params = 4 * u_dim + 2  # [uₖ, uₖ₊₁, duₖ, duₖ₊₁, Δt, t]
    else
        error("Unsupported spline order: $order")
    end

    if isempty(drives)
        return nothing, n_params
    end

    n = ketdim

    # Pre-build operator wrappers from drive Hamiltonians
    drive_ops = [_to_operator(d.H) for d in drives]

    # Precompute control→drives mapping for structural sparsity (matches
    # build_sensitivity_ode). control_to_drives[j] lists drive indices whose
    # ∂c/∂u_j can be nonzero.
    control_to_drives = [Int[] for _ = 1:u_dim]
    drives_active = [active_controls(d) for d in drives]
    for (t_idx, ac) in enumerate(drives_active)
        if isempty(ac)
            for j = 1:u_dim
                push!(control_to_drives[j], t_idx)
            end
        else
            for j in ac
                push!(control_to_drives[j], t_idx)
            end
        end
    end

    make_f!_jvp = if order == 1
        () -> begin
            u_interp = zeros(u_dim)
            (dx, x, p, τ) -> begin
                # Standard params (first 2u_dim+2 slots), then v_p tail
                uₖ = @view p[1:u_dim]
                uₖ₊₁ = @view p[(u_dim+1):2u_dim]
                Δtₖ = p[2u_dim+1]
                # tₖ   = p[2u_dim + 2]  # unused; time-independent forcing
                v_p_offset = 2u_dim + 2  # v_p at p[v_p_offset+1 : v_p_offset+n_params]

                onemτ = 1 - τ
                @inbounds for i = 1:u_dim
                    u_interp[i] = onemτ * uₖ[i] + τ * uₖ₊₁[i]
                end

                ψ = @view x[1:n]
                Y = @view x[(n+1):2n]
                dψ = @view dx[1:n]
                dY = @view dx[(n+1):2n]

                # ψ̇ = -iΔt·H(τ)·ψ
                apply!(dψ, drift_op, ψ, -im * Δtₖ, false)
                @inbounds for t_idx = 1:n_terms
                    c = drive_coeff(drives[t_idx], u_interp)
                    apply!(dψ, drive_ops[t_idx], ψ, -im * Δtₖ * c, true)
                end

                # Ẏ = -iΔt·H(τ)·Y  (homogeneous part)
                apply!(dY, drift_op, Y, -im * Δtₖ, false)
                @inbounds for t_idx = 1:n_terms
                    c = drive_coeff(drives[t_idx], u_interp)
                    apply!(dY, drive_ops[t_idx], Y, -im * Δtₖ * c, true)
                end

                # Forcing: -iΔt·M_v(τ)·ψ where M_v = Σⱼ (∂_{p_j}H) v_p[j]

                # uₖ[j] sensitivities: spline basis = (1-τ); ODE param j = j (1..u_dim)
                @inbounds for j = 1:u_dim
                    v_pj = p[v_p_offset+j]
                    v_pj == 0.0 && continue
                    for t_idx in control_to_drives[j]
                        dc = drive_coeff_jac(drives[t_idx], u_interp, j)
                        if dc != 0.0
                            apply!(dY, drive_ops[t_idx], ψ, -im * Δtₖ * onemτ * dc * v_pj, true)
                        end
                    end
                end

                # uₖ₊₁[j] sensitivities: spline basis = τ; ODE param = u_dim + j
                @inbounds for j = 1:u_dim
                    v_pj = p[v_p_offset+u_dim+j]
                    v_pj == 0.0 && continue
                    for t_idx in control_to_drives[j]
                        dc = drive_coeff_jac(drives[t_idx], u_interp, j)
                        if dc != 0.0
                            apply!(dY, drive_ops[t_idx], ψ, -im * Δtₖ * τ * dc * v_pj, true)
                        end
                    end
                end

                # Δt sensitivity: ODE param = 2u_dim + 1
                # ∂(-iΔt H)/∂Δt = -i·H, applied to ψ
                let v_pj = p[v_p_offset+2u_dim+1]
                    if v_pj != 0.0
                        apply!(dY, drift_op, ψ, -im * v_pj, true)
                        @inbounds for t_idx = 1:n_terms
                            c = drive_coeff(drives[t_idx], u_interp)
                            apply!(dY, drive_ops[t_idx], ψ, -im * c * v_pj, true)
                        end
                    end
                end

                # t sensitivity: zero forcing for time-independent drives
                # (slot exists at v_p_offset + 2u_dim + 2 for indexing parity)

                return nothing
            end
        end

    else  # order == 3
        () -> begin
            u_interp = zeros(u_dim)
            jac_cache = zeros(u_dim, n_terms)
            (dx, x, p, τ) -> begin
                uₖ = @view p[1:u_dim]
                uₖ₊₁ = @view p[(u_dim+1):2u_dim]
                duₖ = @view p[(2u_dim+1):3u_dim]
                duₖ₊₁ = @view p[(3u_dim+1):4u_dim]
                Δtₖ = p[4u_dim+1]
                # tₖ  = p[4u_dim + 2]  # unused
                v_p_offset = 4u_dim + 2

                τ2 = τ * τ
                τ3 = τ2 * τ
                h00 = 2τ3 - 3τ2 + 1
                h10 = (τ3 - 2τ2 + τ) * Δtₖ
                h01 = -2τ3 + 3τ2
                h11 = (τ3 - τ2) * Δtₖ
                dh10_dΔt = τ3 - 2τ2 + τ
                dh11_dΔt = τ3 - τ2

                @inbounds for i = 1:u_dim
                    u_interp[i] =
                        h00 * uₖ[i] + h10 * duₖ[i] + h01 * uₖ₊₁[i] + h11 * duₖ₊₁[i]
                end

                # Op-fusion: per-substep jac coefficient cache.
                # jac_cache[j, t_idx] = ∂c_t/∂u_j; non-active entries left at 0.
                @inbounds for t_idx = 1:n_terms
                    ac = drives_active[t_idx]
                    iter = isempty(ac) ? (1:u_dim) : ac
                    for j = 1:u_dim
                        jac_cache[j, t_idx] = 0.0
                    end
                    for j in iter
                        jac_cache[j, t_idx] = drive_coeff_jac(drives[t_idx], u_interp, j)
                    end
                end

                # Op-fusion: wide-matrix homogeneous A·block on [ψ Y].
                # [ψ; Y] is contiguous in x; reshape as n × 2 wide matrix.
                wide_x = reshape(@view(x[1:2n]), n, 2)
                wide_dx = reshape(@view(dx[1:2n]), n, 2)
                apply!(wide_dx, drift_op, wide_x, -im * Δtₖ, false)
                @inbounds for t_idx = 1:n_terms
                    c = drive_coeff(drives[t_idx], u_interp)
                    apply!(wide_dx, drive_ops[t_idx], wide_x, -im * Δtₖ * c, true)
                end

                # Per-block views for inhomogeneous forcing (only ψ source and dY dest)
                ψ = @view x[1:n]
                dY = @view dx[(n+1):2n]

                # uₖ[j] sensitivities: basis = h00; ODE param = j
                @inbounds for j = 1:u_dim
                    v_pj = p[v_p_offset+j]
                    v_pj == 0.0 && continue
                    for t_idx in control_to_drives[j]
                        dc = jac_cache[j, t_idx]
                        if dc != 0.0
                            apply!(dY, drive_ops[t_idx], ψ, -im * Δtₖ * h00 * dc * v_pj, true)
                        end
                    end
                end

                # uₖ₊₁[j] sensitivities: basis = h01; ODE param = u_dim + j
                @inbounds for j = 1:u_dim
                    v_pj = p[v_p_offset+u_dim+j]
                    v_pj == 0.0 && continue
                    for t_idx in control_to_drives[j]
                        dc = jac_cache[j, t_idx]
                        if dc != 0.0
                            apply!(dY, drive_ops[t_idx], ψ, -im * Δtₖ * h01 * dc * v_pj, true)
                        end
                    end
                end

                # duₖ[j] sensitivities: basis = h10; ODE param = 2u_dim + j
                @inbounds for j = 1:u_dim
                    v_pj = p[v_p_offset+2u_dim+j]
                    v_pj == 0.0 && continue
                    for t_idx in control_to_drives[j]
                        dc = jac_cache[j, t_idx]
                        if dc != 0.0
                            apply!(dY, drive_ops[t_idx], ψ, -im * Δtₖ * h10 * dc * v_pj, true)
                        end
                    end
                end

                # duₖ₊₁[j] sensitivities: basis = h11; ODE param = 3u_dim + j
                @inbounds for j = 1:u_dim
                    v_pj = p[v_p_offset+3u_dim+j]
                    v_pj == 0.0 && continue
                    for t_idx in control_to_drives[j]
                        dc = jac_cache[j, t_idx]
                        if dc != 0.0
                            apply!(dY, drive_ops[t_idx], ψ, -im * Δtₖ * h11 * dc * v_pj, true)
                        end
                    end
                end

                # Δt sensitivity: ODE param = 4u_dim + 1
                # ∂(-iΔtH(u(τ;Δt)))/∂Δt = -iH(u) + (-iΔt) Σ_drives (∂c/∂u_j)·(∂u_j/∂Δt)·H_d
                # where ∂u_j/∂Δt comes from the Hermite basis (h10, h11 carry an explicit
                # Δt factor: h10 = (τ3-2τ2+τ)·Δt, h11 = (τ3-τ2)·Δt). The cross terms via
                # `total_du_dΔt` are non-trivial — must port to keep Δt-direction tangents
                # correct. See spline_integrator_type.jl:1429-1457.
                let v_pj = p[v_p_offset+4u_dim+1]
                    if v_pj != 0.0
                        # -i·H(u)·ψ contribution
                        apply!(dY, drift_op, ψ, -im * v_pj, true)
                        @inbounds for t_idx = 1:n_terms
                            c = drive_coeff(drives[t_idx], u_interp)
                            apply!(dY, drive_ops[t_idx], ψ, -im * c * v_pj, true)
                        end
                        # -iΔt · Σ_drives coeff_jac · ∂u/∂Δt · H_d · ψ contribution
                        @inbounds for t_idx = 1:n_terms
                            op = drive_ops[t_idx]
                            ac = drives_active[t_idx]
                            iter = isempty(ac) ? (1:u_dim) : ac
                            total_du_dΔt = 0.0
                            for i in iter
                                dc = jac_cache[i, t_idx]
                                if dc != 0.0
                                    du_i_dΔt = dh10_dΔt * duₖ[i] + dh11_dΔt * duₖ₊₁[i]
                                    total_du_dΔt += dc * du_i_dΔt
                                end
                            end
                            if total_du_dΔt != 0.0
                                apply!(dY, op, ψ, -im * Δtₖ * total_du_dΔt * v_pj, true)
                            end
                        end
                    end
                end

                # t sensitivity: zero forcing (slot at v_p_offset + 4u_dim + 2)

                return nothing
            end
        end
    end

    return make_f!_jvp, n_params
end

"""
    build_hvp_forward_ode(drift_op, drives, u_dim, ketdim, order)
        -> (make_f!_hvp_fwd, n_params)

Forward augmented ODE for HVP forward-over-reverse (Phase 4). Dynamics are
**structurally identical** to `build_ket_jvp_ode` — same state `[ψ; δψ] ∈ ℂ^{2d}`,
same forcing `δψ̇ = -iΔt H δψ - iΔt M_v ψ` where `M_v = Σⱼ (∂_{p_j}H) v_{p_j}`.

The only purpose of this thin alias is to reserve a separate `Vector{ODEProblem}`
slot in `Tsit5Data.hvp_fwd_probs` so the HVP forward path can be specialized
in future phases (e.g., dense-output retention) without disturbing Phase 2's
JVP path. v0.1 just delegates to `build_ket_jvp_ode`.

Initial conditions (set by `step_hvp!` caller): `ψ(0) = ψ_k`, `δψ(0) = -v_{x_k}`.
After solve: `δψ(1)` carries the JVP at this knot (Φ_p·v_p − Φ·v_{x_k}).
"""
function build_hvp_forward_ode(
    drift_op::AbstractDynamicsOperator,
    drives::AbstractVector{<:AbstractDrive},
    u_dim::Int,
    ketdim::Int,
    order::Int,
)
    return build_ket_jvp_ode(drift_op, drives, u_dim, ketdim, order)
end

"""
    build_adjoint_ode(drift_op, drives, u_dim, ketdim, order)
        -> (make_f!_adj, n_params)

Build an in-place ODE RHS for the **ket-level adjoint VJP backward solve** used
by Phase 3's matrix-free VJP. State layout: `x = [ψ; λ; g] ∈ ℂ^{2d + n_p}` where
`d = ketdim`.

Dynamics (direction-agnostic; integrate via `tspan` to control direction):

    ψ̇ = -iΔt·H(τ; u(τ)) · ψ
    λ̇ = -iΔt·H(τ; u(τ)) · λ
    ġ_i = ⟨ -iΔt ∂_{p_i} H(τ; u(τ)) ψ, λ ⟩_ℂ   for i = 1..n_p

(Julia convention: `dot(a,b) = Σ conj(a) b`.)

For Hermitian H, the forward ψ-equation reversed in time gives backward
unitary evolution; the adjoint λ uses the same dynamics. Integration backward
(`tspan = (1.0, 0.0)`) with terminal conditions `ψ(1) = Φ ψ_k`, `λ(1) = w`,
`g(1) = 0` produces `ψ(0) = ψ_k`, `λ(0) = Φ' w`, and
`g(0) = -⟨∂_{p_i}ψ_{k+1}^{prop}, w⟩_ℂ` (sign from the backward integration of
the variational identity; see `step_vjp!` for the iso-output convention).

This contracts the work that `build_sensitivity_ode` does in parallel for `n_p`
sensitivities into a single backward solve, accumulating the parameter
contribution as a complex inner product into `dg[j]`. The per-solve state size
drops from `(1 + n_p)·d²` complex to `2d + n_p` complex — the goal of Phase 3.

Returns `(make_f!_adj, n_params)` matching `build_sensitivity_ode`'s and
`build_ket_jvp_ode`'s shape.
"""
function build_adjoint_ode(
    drift_op::AbstractDynamicsOperator,
    drives::AbstractVector{<:AbstractDrive},
    u_dim::Int,
    ketdim::Int,
    order::Int,
)
    n_terms = length(drives)

    if order == 1
        n_params = 2 * u_dim + 2  # [uₖ, uₖ₊₁, Δt, t]
    elseif order == 3
        n_params = 4 * u_dim + 2  # [uₖ, uₖ₊₁, duₖ, duₖ₊₁, Δt, t]
    else
        error("Unsupported spline order: $order")
    end

    if isempty(drives)
        return nothing, n_params
    end

    n = ketdim

    # Pre-build operator wrappers from drive Hamiltonians
    drive_ops = [_to_operator(d.H) for d in drives]

    # Precompute control→drives mapping for structural sparsity (matches
    # build_sensitivity_ode and build_ket_jvp_ode). control_to_drives[j] lists
    # drive indices whose ∂c/∂u_j can be nonzero.
    control_to_drives = [Int[] for _ = 1:u_dim]
    drives_active = [active_controls(d) for d in drives]
    for (t_idx, ac) in enumerate(drives_active)
        if isempty(ac)
            for j = 1:u_dim
                push!(control_to_drives[j], t_idx)
            end
        else
            for j in ac
                push!(control_to_drives[j], t_idx)
            end
        end
    end

    make_f!_adj = if order == 1
        () -> begin
            # Closure-captured scratch — fresh per ODEProblem instance.
            # Each per-knot ODEProblem in vjp_probs has its own closure (and
            # thus its own scratch), so Threads.@threads over knots is race-
            # free — each thread holds a different ODEProblem and uses its
            # closure's scratch. Mirrors Phase 2's build_ket_jvp_ode pattern
            # for `u_interp`.
            u_interp = zeros(u_dim)
            tmp_complex = Vector{ComplexF64}(undef, n)
            (dx, x, p, τ) -> begin
                # Standard params (first 2u_dim+2 slots); no v_p tail (VJP
                # does not contract over params at the parameter level — the
                # contraction happens via the complex inner product in dg).
                uₖ = @view p[1:u_dim]
                uₖ₊₁ = @view p[(u_dim+1):2u_dim]
                Δtₖ = p[2u_dim+1]
                # tₖ  = p[2u_dim + 2]  # unused; time-independent forcing

                onemτ = 1 - τ
                @inbounds for i = 1:u_dim
                    u_interp[i] = onemτ * uₖ[i] + τ * uₖ₊₁[i]
                end

                ψ = @view x[1:n]
                λ = @view x[(n+1):2n]
                # g  = @view x[(2n+1):(2n+n_params)]  # only written via dg
                dψ = @view dx[1:n]
                dλ = @view dx[(n+1):2n]
                dg = @view dx[(2n+1):(2n+n_params)]

                # ψ̇ = -iΔt·H(τ)·ψ
                apply!(dψ, drift_op, ψ, -im * Δtₖ, false)
                @inbounds for t_idx = 1:n_terms
                    c = drive_coeff(drives[t_idx], u_interp)
                    apply!(dψ, drive_ops[t_idx], ψ, -im * Δtₖ * c, true)
                end

                # λ̇ = -iΔt·H(τ)·λ  (same dynamics as ψ; for Hermitian H,
                # this serves as the adjoint when integrated backward)
                apply!(dλ, drift_op, λ, -im * Δtₖ, false)
                @inbounds for t_idx = 1:n_terms
                    c = drive_coeff(drives[t_idx], u_interp)
                    apply!(dλ, drive_ops[t_idx], λ, -im * Δtₖ * c, true)
                end

                # ġ_i = ⟨ -iΔt ∂_{p_i} H(τ) ψ, λ ⟩  for each parameter i
                # = -conj(-iΔt) · Σ_k conj(∂_i H ψ)_k · λ_k
                # = (-iΔt) · Σ_k (∂_i H ψ)_k * λ_k... no, we use Julia dot
                # convention: dot(a,b) = Σ conj(a) b. So
                # dg[i] = dot(-iΔt ∂_i H ψ, λ) = conj(-iΔt) dot(∂_i H ψ, λ)
                #       = (iΔt) dot(∂_i H ψ, λ).
                # We compute tmp = ∂_i H ψ first, then accumulate
                # acc = sum(conj(tmp_k) * λ_k), and write
                # dg[i] = (iΔt) * acc... but the plan code says
                # dg[i] = -iΔt * acc; the per-knot equivalence test below
                # is the source of truth for the sign convention. Plan
                # comment justifies the joint sign cancellation with
                # `out[col] += -real(g_zero[ode_idx])`.

                # uₖ[j] sensitivities: spline basis = (1-τ); ODE param j = j (1..u_dim)
                @inbounds for j = 1:u_dim
                    fill!(tmp_complex, 0)
                    any_nz = false
                    for t_idx in control_to_drives[j]
                        dc = drive_coeff_jac(drives[t_idx], u_interp, j)
                        if dc != 0.0
                            apply!(tmp_complex, drive_ops[t_idx], ψ, onemτ * dc, true)
                            any_nz = true
                        end
                    end
                    if any_nz
                        acc = zero(ComplexF64)
                        for k = 1:n
                            acc += conj(tmp_complex[k]) * λ[k]
                        end
                        dg[j] = -im * Δtₖ * acc
                    else
                        dg[j] = zero(ComplexF64)
                    end
                end

                # uₖ₊₁[j] sensitivities: spline basis = τ; ODE param = u_dim + j
                @inbounds for j = 1:u_dim
                    fill!(tmp_complex, 0)
                    any_nz = false
                    for t_idx in control_to_drives[j]
                        dc = drive_coeff_jac(drives[t_idx], u_interp, j)
                        if dc != 0.0
                            apply!(tmp_complex, drive_ops[t_idx], ψ, τ * dc, true)
                            any_nz = true
                        end
                    end
                    if any_nz
                        acc = zero(ComplexF64)
                        for k = 1:n
                            acc += conj(tmp_complex[k]) * λ[k]
                        end
                        dg[u_dim+j] = -im * Δtₖ * acc
                    else
                        dg[u_dim+j] = zero(ComplexF64)
                    end
                end

                # Δt sensitivity: ODE param = 2u_dim + 1
                # ∂(-iΔtH)/∂Δt = -i·H, applied to ψ
                let
                    fill!(tmp_complex, 0)
                    apply!(tmp_complex, drift_op, ψ, one(ComplexF64), true)
                    @inbounds for t_idx = 1:n_terms
                        c = drive_coeff(drives[t_idx], u_interp)
                        apply!(tmp_complex, drive_ops[t_idx], ψ, ComplexF64(c), true)
                    end
                    acc = zero(ComplexF64)
                    @inbounds for k = 1:n
                        acc += conj(tmp_complex[k]) * λ[k]
                    end
                    # ġ_{Δt} = ⟨-i·H ψ, λ⟩ = (i)·dot(H ψ, λ)
                    # Joint sign convention with output negation: write the
                    # raw -i·acc here, the output side flips sign.
                    dg[2u_dim+1] = -im * acc
                end

                # t sensitivity: zero forcing for time-independent drives
                # (slot exists at 2u_dim + 2 for indexing parity)
                dg[2u_dim+2] = zero(ComplexF64)

                return nothing
            end
        end

    else  # order == 3
        () -> begin
            u_interp = zeros(u_dim)
            tmp_complex = Vector{ComplexF64}(undef, n)
            (dx, x, p, τ) -> begin
                uₖ = @view p[1:u_dim]
                uₖ₊₁ = @view p[(u_dim+1):2u_dim]
                duₖ = @view p[(2u_dim+1):3u_dim]
                duₖ₊₁ = @view p[(3u_dim+1):4u_dim]
                Δtₖ = p[4u_dim+1]
                # tₖ  = p[4u_dim + 2]  # unused

                τ2 = τ * τ
                τ3 = τ2 * τ
                h00 = 2τ3 - 3τ2 + 1
                h10 = (τ3 - 2τ2 + τ) * Δtₖ
                h01 = -2τ3 + 3τ2
                h11 = (τ3 - τ2) * Δtₖ
                dh10_dΔt = τ3 - 2τ2 + τ
                dh11_dΔt = τ3 - τ2

                @inbounds for i = 1:u_dim
                    u_interp[i] =
                        h00 * uₖ[i] + h10 * duₖ[i] + h01 * uₖ₊₁[i] + h11 * duₖ₊₁[i]
                end

                ψ = @view x[1:n]
                λ = @view x[(n+1):2n]
                dψ = @view dx[1:n]
                dλ = @view dx[(n+1):2n]
                dg = @view dx[(2n+1):(2n+n_params)]

                # ψ̇ = -iΔt·H(τ)·ψ
                apply!(dψ, drift_op, ψ, -im * Δtₖ, false)
                @inbounds for t_idx = 1:n_terms
                    c = drive_coeff(drives[t_idx], u_interp)
                    apply!(dψ, drive_ops[t_idx], ψ, -im * Δtₖ * c, true)
                end

                # λ̇ = -iΔt·H(τ)·λ
                apply!(dλ, drift_op, λ, -im * Δtₖ, false)
                @inbounds for t_idx = 1:n_terms
                    c = drive_coeff(drives[t_idx], u_interp)
                    apply!(dλ, drive_ops[t_idx], λ, -im * Δtₖ * c, true)
                end

                # uₖ[j]: basis = h00; ODE param = j
                @inbounds for j = 1:u_dim
                    fill!(tmp_complex, 0)
                    any_nz = false
                    for t_idx in control_to_drives[j]
                        dc = drive_coeff_jac(drives[t_idx], u_interp, j)
                        if dc != 0.0
                            apply!(tmp_complex, drive_ops[t_idx], ψ, h00 * dc, true)
                            any_nz = true
                        end
                    end
                    if any_nz
                        acc = zero(ComplexF64)
                        for k = 1:n
                            acc += conj(tmp_complex[k]) * λ[k]
                        end
                        dg[j] = -im * Δtₖ * acc
                    else
                        dg[j] = zero(ComplexF64)
                    end
                end

                # uₖ₊₁[j]: basis = h01; ODE param = u_dim + j
                @inbounds for j = 1:u_dim
                    fill!(tmp_complex, 0)
                    any_nz = false
                    for t_idx in control_to_drives[j]
                        dc = drive_coeff_jac(drives[t_idx], u_interp, j)
                        if dc != 0.0
                            apply!(tmp_complex, drive_ops[t_idx], ψ, h01 * dc, true)
                            any_nz = true
                        end
                    end
                    if any_nz
                        acc = zero(ComplexF64)
                        for k = 1:n
                            acc += conj(tmp_complex[k]) * λ[k]
                        end
                        dg[u_dim+j] = -im * Δtₖ * acc
                    else
                        dg[u_dim+j] = zero(ComplexF64)
                    end
                end

                # duₖ[j]: basis = h10; ODE param = 2u_dim + j
                @inbounds for j = 1:u_dim
                    fill!(tmp_complex, 0)
                    any_nz = false
                    for t_idx in control_to_drives[j]
                        dc = drive_coeff_jac(drives[t_idx], u_interp, j)
                        if dc != 0.0
                            apply!(tmp_complex, drive_ops[t_idx], ψ, h10 * dc, true)
                            any_nz = true
                        end
                    end
                    if any_nz
                        acc = zero(ComplexF64)
                        for k = 1:n
                            acc += conj(tmp_complex[k]) * λ[k]
                        end
                        dg[2u_dim+j] = -im * Δtₖ * acc
                    else
                        dg[2u_dim+j] = zero(ComplexF64)
                    end
                end

                # duₖ₊₁[j]: basis = h11; ODE param = 3u_dim + j
                @inbounds for j = 1:u_dim
                    fill!(tmp_complex, 0)
                    any_nz = false
                    for t_idx in control_to_drives[j]
                        dc = drive_coeff_jac(drives[t_idx], u_interp, j)
                        if dc != 0.0
                            apply!(tmp_complex, drive_ops[t_idx], ψ, h11 * dc, true)
                            any_nz = true
                        end
                    end
                    if any_nz
                        acc = zero(ComplexF64)
                        for k = 1:n
                            acc += conj(tmp_complex[k]) * λ[k]
                        end
                        dg[3u_dim+j] = -im * Δtₖ * acc
                    else
                        dg[3u_dim+j] = zero(ComplexF64)
                    end
                end

                # Δt sensitivity: ODE param = 4u_dim + 1
                # ∂(-iΔtH(u(τ;Δt)))/∂Δt = -iH(u) + (-iΔt) Σ_drives (∂c/∂u_j)·(∂u_j/∂Δt)·H_d
                # Mirror Phase 2's two-part contribution (drift + drives) plus
                # the cross-term via total_du_dΔt for the cubic basis.
                let
                    # -i·H(u)·ψ contribution: tmp_a = H(u) ψ
                    fill!(tmp_complex, 0)
                    apply!(tmp_complex, drift_op, ψ, one(ComplexF64), true)
                    @inbounds for t_idx = 1:n_terms
                        c = drive_coeff(drives[t_idx], u_interp)
                        apply!(tmp_complex, drive_ops[t_idx], ψ, ComplexF64(c), true)
                    end
                    # Add the cross-term contribution: Δt · Σ_drives total_du_dΔt · H_d ψ
                    @inbounds for t_idx = 1:n_terms
                        op = drive_ops[t_idx]
                        ac = drives_active[t_idx]
                        iter = isempty(ac) ? (1:u_dim) : ac
                        total_du_dΔt = 0.0
                        for i in iter
                            dc = drive_coeff_jac(drives[t_idx], u_interp, i)
                            if dc != 0.0
                                du_i_dΔt = dh10_dΔt * duₖ[i] + dh11_dΔt * duₖ₊₁[i]
                                total_du_dΔt += dc * du_i_dΔt
                            end
                        end
                        if total_du_dΔt != 0.0
                            apply!(tmp_complex, op, ψ, ComplexF64(Δtₖ * total_du_dΔt), true)
                        end
                    end
                    acc = zero(ComplexF64)
                    @inbounds for k = 1:n
                        acc += conj(tmp_complex[k]) * λ[k]
                    end
                    # ġ_{Δt} = ⟨-i·(H ψ + Δt·cross), λ⟩, raw form (sign cancels
                    # with output negation in step_vjp!).
                    dg[4u_dim+1] = -im * acc
                end

                # t sensitivity: zero forcing (slot at 4u_dim + 2)
                dg[4u_dim+2] = zero(ComplexF64)

                return nothing
            end
        end
    end

    return make_f!_adj, n_params
end

"""
    build_second_order_adjoint_ode(drift_op, drives, u_dim, ketdim, order)
        -> (make_f!_2nd_adj, n_params)

Backward second-order adjoint ODE for HVP forward-over-reverse (Phase 4).

State layout: x ∈ ℂ^{4d + n_p}
    x[1:d]               = ψ'(τ)        — round-trip ψ via unitary reversal
    x[d+1:2d]            = δψ'(τ)        — round-trip δψ via tangent ODE
    x[2d+1:3d]           = λ(τ)          — gradient adjoint (Phase 3 dynamics)
    x[3d+1:4d]           = ν(τ)          — second-order adjoint
    x[4d+1:4d+n_p]       = g(τ)          — HVP integrand accumulator

Parameter layout (same as `build_adjoint_ode`/`build_ket_jvp_ode`):
    p[1:n_p]                   = standard params [uₖ, uₖ₊₁, (duₖ, duₖ₊₁), Δt, t]
    p[n_p+1 : 2*n_p]           = v_p tail (HVP direction's parameter block)

Dynamics (direction-agnostic; tspan controls direction):

    ψ'̇  = -iΔt H(τ) ψ'                   (forward; reversible)
    δψ'̇ = -iΔt H(τ) δψ' - iΔt M_v(τ) ψ'   (forward tangent eq, reversible)
    λ̇  = -iΔt H(τ) λ                     (gradient adjoint, Phase 3)
    ν̇  = -iΔt H(τ) ν - iΔt M_v(τ)ᵀ λ      (second-order adjoint forcing)
    ġ_i = ⟨-iΔt ∂_{p_i}H δψ', λ⟩ +
          ⟨-iΔt ∂_{p_i}H ψ', ν⟩ +
          Σⱼ v_{p_j} ⟨-iΔt ∂²_{p_i p_j}H ψ', λ⟩

For Hermitian H, the ν forcing M_v^†λ uses the same `apply!` path as the
forward forcing — H_d is Hermitian, so `apply!(dν, H_d, λ, scale)` computes
-iΔt·dc·bw·v_pj·H_d λ = (-iΔt M_v λ), the same expression as for δψ'.

The third forcing term ġ_i involves `drive_coeff_hess`. For LinearDrive
this is identically zero. For NonlinearDrive (e.g. |α|²), it's nonzero
and is precisely the term Gauss-Newton drops, causing the Stanford-bosonic
catastrophic failure mode. **Including this term is the entire point of Phase 4.**

Solved with `tspan = (1.0, 0.0)`. Terminal conditions set by `step_hvp!`:
    ψ'(1)  = ψ(1)        from forward solve
    δψ'(1) = δψ(1)       from forward solve
    λ(1)   = μ_k         (from optimizer)
    ν(1)   = 0           (linear-in-ψ output, ∇²g = 0)
    g(1)   = 0
"""
function build_second_order_adjoint_ode(
    drift_op::AbstractDynamicsOperator,
    drives::AbstractVector{<:AbstractDrive},
    u_dim::Int,
    ketdim::Int,
    order::Int,
)
    n_terms = length(drives)

    if order == 1
        n_params = 2 * u_dim + 2  # [uₖ, uₖ₊₁, Δt, t]
    elseif order == 3
        n_params = 4 * u_dim + 2  # [uₖ, uₖ₊₁, duₖ, duₖ₊₁, Δt, t]
    else
        error("Unsupported spline order: $order")
    end

    if isempty(drives)
        return nothing, n_params
    end

    n = ketdim

    # Pre-build operator wrappers from drive Hamiltonians
    drive_ops = [_to_operator(d.H) for d in drives]

    # control→drives mapping (matches build_adjoint_ode/build_ket_jvp_ode)
    control_to_drives = [Int[] for _ = 1:u_dim]
    drives_active = [active_controls(d) for d in drives]
    for (t_idx, ac) in enumerate(drives_active)
        if isempty(ac)
            for j = 1:u_dim
                push!(control_to_drives[j], t_idx)
            end
        else
            for j in ac
                push!(control_to_drives[j], t_idx)
            end
        end
    end

    # #353: per-drive control-index list, with the "no declared active controls ⇒
    # relevant to every control" expansion done ONCE here. The cubic branch's RHS
    # needs to sum over a drive's controls (for the Δt chain rule), and the
    # `isempty(ac) ? (1:u_dim) : ac` idiom used elsewhere is a small Union that
    # makes the loop dynamic — this RHS is on the knot-flat allocation contract
    # (#342 / ADR-0009 decision 3), so it iterates a concrete `Vector{Int}`.
    drive_controls = [isempty(ac) ? collect(1:u_dim) : collect(ac) for ac in drives_active]

    # Build the parameter → u-component mapping for the drive_coeff_hess
    # chain rule. drive_coeff_hess(d, u, m, n) is indexed by u-component
    # indices, NOT parameter indices. Each parameter p_i corresponds to a
    # specific u-component m_i with a τ-dependent basis weight.
    #
    # Layout (linear, n_params = 2*u_dim + 2):
    #   p_1..p_{u_dim}         → uₖ[1..u_dim],   bw = (1-τ)
    #   p_{u_dim+1..2*u_dim}    → uₖ₊₁[1..u_dim], bw = τ
    #   p_{2*u_dim+1}           → Δtₖ            (no u-component; m_i = 0)
    #   p_{2*u_dim+2}           → tₖ             (no u-component; m_i = 0)
    #
    # Layout (cubic, n_params = 4*u_dim + 2):
    #   p_1..p_{u_dim}         → uₖ,           bw = h00(τ)
    #   p_{u_dim+1..2*u_dim}   → uₖ₊₁,         bw = h01(τ)
    #   p_{2*u_dim+1..3*u_dim} → duₖ,          bw = h10(τ)·Δtₖ  *
    #   p_{3*u_dim+1..4*u_dim} → duₖ₊₁,        bw = h11(τ)·Δtₖ  *
    #   p_{4*u_dim+1}          → Δtₖ           (m_i = 0)
    #   p_{4*u_dim+2}          → tₖ            (m_i = 0)
    #
    # * Cubic du-parameter weights carry an explicit Δtₖ factor, so for the
    #   cubic cell Δt enters BOTH through the basis weights and through the
    #   interpolated control values u(τ; Δt). #353: the cubic branch below
    #   carries every path of that chain rule. Writing h̃10 = τ³-2τ²+τ and
    #   h̃11 = τ³-τ² (so h10 = h̃10·Δt, h11 = h̃11·Δt) and
    #   D_m ≡ ∂u_m/∂Δt = h̃10·duₖ[m] + h̃11·duₖ₊₁[m]:
    #
    #     ∂A/∂Δt          = -i·H(u) - iΔt·Σ_d (Σ_m ∂c_d/∂u_m · D_m)·H_d
    #     ∂²A/∂p_i∂Δt     = -i·b_i·Σ_d ∂c_d/∂u_{m_i}·H_d                  (a)
    #                       -iΔt·b_i·Σ_d (Σ_m ∂²c_d/∂u_{m_i}∂u_m · D_m)·H_d (b)
    #                       -iΔt·(∂b_i/∂Δt)·Σ_d ∂c_d/∂u_{m_i}·H_d          (c)
    #     ∂²A/∂Δt²        = -2i·Σ_d (Σ_m ∂c_d/∂u_m · D_m)·H_d
    #                       -iΔt·Σ_d (Σ_{m,m'} ∂²c_d/∂u_m∂u_{m'} D_m D_{m'})·H_d
    #
    #   with ∂b_i/∂Δt = h̃10 for the duₖ group, h̃11 for duₖ₊₁, 0 otherwise.
    #   Before #353 the cubic branch had only (a) — the Δt-chain-rule part of
    #   ∂A/∂Δt was missing from the δψ' / ν forcings and from the Δt-slot T1/T2,
    #   (b) and (c) were missing from the (u, Δt) cross-Hessian, and the whole
    #   (Δt, Δt) diagonal was absent. Measured on the Ket Tsit5 cubic cell that
    #   was a 4.5e-1 relative error in the `du` block and 6.0e-2 in `Δt`.

    make_f!_2nd_adj = if order == 1
        () -> begin
            # Closure-captured scratch — fresh per ODEProblem
            u_interp = zeros(u_dim)
            tmp1 = Vector{ComplexF64}(undef, n)
            tmp2 = Vector{ComplexF64}(undef, n)

            (dx, x, p, τ) -> begin
                uₖ = @view p[1:u_dim]
                uₖ₊₁ = @view p[(u_dim+1):2u_dim]
                Δtₖ = p[2u_dim+1]
                # tₖ  = p[2u_dim + 2]  # unused
                v_p_offset = 2u_dim + 2  # v_p at p[v_p_offset+1 : v_p_offset+n_params]

                onemτ = 1 - τ
                @inbounds for i = 1:u_dim
                    u_interp[i] = onemτ * uₖ[i] + τ * uₖ₊₁[i]
                end

                # State views
                ψ′ = @view x[1:n]
                δψ′ = @view x[(n+1):2n]
                λ = @view x[(2n+1):3n]
                ν = @view x[(3n+1):4n]
                # g  in-place via dx writes only
                dψ′ = @view dx[1:n]
                dδψ′ = @view dx[(n+1):2n]
                dλ = @view dx[(2n+1):3n]
                dν = @view dx[(3n+1):4n]
                dg = @view dx[(4n+1):(4n+n_params)]

                # ── ψ'̇ = -iΔt H ψ' (round-trip ψ via unitary reversal) ──
                apply!(dψ′, drift_op, ψ′, -im * Δtₖ, false)
                @inbounds for t_idx = 1:n_terms
                    c = drive_coeff(drives[t_idx], u_interp)
                    apply!(dψ′, drive_ops[t_idx], ψ′, -im * Δtₖ * c, true)
                end

                # ── δψ'̇ = -iΔt H δψ' - iΔt M_v ψ' ──
                # Homogeneous part
                apply!(dδψ′, drift_op, δψ′, -im * Δtₖ, false)
                @inbounds for t_idx = 1:n_terms
                    c = drive_coeff(drives[t_idx], u_interp)
                    apply!(dδψ′, drive_ops[t_idx], δψ′, -im * Δtₖ * c, true)
                end
                # Forcing M_v ψ' : uₖ contribution
                @inbounds for j = 1:u_dim
                    v_pj = p[v_p_offset+j]
                    v_pj == 0.0 && continue
                    for t_idx in control_to_drives[j]
                        dc = drive_coeff_jac(drives[t_idx], u_interp, j)
                        if dc != 0.0
                            apply!(
                                dδψ′,
                                drive_ops[t_idx],
                                ψ′,
                                -im * Δtₖ * onemτ * dc * v_pj,
                                true,
                            )
                        end
                    end
                end
                # Forcing M_v ψ' : uₖ₊₁ contribution
                @inbounds for j = 1:u_dim
                    v_pj = p[v_p_offset+u_dim+j]
                    v_pj == 0.0 && continue
                    for t_idx in control_to_drives[j]
                        dc = drive_coeff_jac(drives[t_idx], u_interp, j)
                        if dc != 0.0
                            apply!(dδψ′, drive_ops[t_idx], ψ′, -im * Δtₖ * τ * dc * v_pj, true)
                        end
                    end
                end
                # Forcing M_v ψ' : Δt contribution (∂(-iΔtH)/∂Δt = -iH)
                let v_pj = p[v_p_offset+2u_dim+1]
                    if v_pj != 0.0
                        apply!(dδψ′, drift_op, ψ′, -im * v_pj, true)
                        @inbounds for t_idx = 1:n_terms
                            c = drive_coeff(drives[t_idx], u_interp)
                            apply!(dδψ′, drive_ops[t_idx], ψ′, -im * c * v_pj, true)
                        end
                    end
                end
                # t slot has zero forcing

                # ── λ̇ = -iΔt H λ ──
                apply!(dλ, drift_op, λ, -im * Δtₖ, false)
                @inbounds for t_idx = 1:n_terms
                    c = drive_coeff(drives[t_idx], u_interp)
                    apply!(dλ, drive_ops[t_idx], λ, -im * Δtₖ * c, true)
                end

                # ── ν̇ = -iΔt H ν - iΔt M_v λ ──
                # (M_v^† λ for Hermitian H equals M_v λ)
                # Homogeneous part
                apply!(dν, drift_op, ν, -im * Δtₖ, false)
                @inbounds for t_idx = 1:n_terms
                    c = drive_coeff(drives[t_idx], u_interp)
                    apply!(dν, drive_ops[t_idx], ν, -im * Δtₖ * c, true)
                end
                # Forcing M_v λ : uₖ contribution
                @inbounds for j = 1:u_dim
                    v_pj = p[v_p_offset+j]
                    v_pj == 0.0 && continue
                    for t_idx in control_to_drives[j]
                        dc = drive_coeff_jac(drives[t_idx], u_interp, j)
                        if dc != 0.0
                            apply!(dν, drive_ops[t_idx], λ, -im * Δtₖ * onemτ * dc * v_pj, true)
                        end
                    end
                end
                # Forcing M_v λ : uₖ₊₁ contribution
                @inbounds for j = 1:u_dim
                    v_pj = p[v_p_offset+u_dim+j]
                    v_pj == 0.0 && continue
                    for t_idx in control_to_drives[j]
                        dc = drive_coeff_jac(drives[t_idx], u_interp, j)
                        if dc != 0.0
                            apply!(dν, drive_ops[t_idx], λ, -im * Δtₖ * τ * dc * v_pj, true)
                        end
                    end
                end
                # Forcing M_v λ : Δt contribution
                let v_pj = p[v_p_offset+2u_dim+1]
                    if v_pj != 0.0
                        apply!(dν, drift_op, λ, -im * v_pj, true)
                        @inbounds for t_idx = 1:n_terms
                            c = drive_coeff(drives[t_idx], u_interp)
                            apply!(dν, drive_ops[t_idx], λ, -im * c * v_pj, true)
                        end
                    end
                end

                # ── ġ_i = three-term integrand ──
                # Helper: param i → (u-component m_i, basis-weight bw_i, is_dt)
                # Linear: uₖ has bw=onemτ, uₖ₊₁ has bw=τ, Δt has m_i=0/dt=true, t has m_i=0/dt=false
                #
                # The three terms:
                # T1: ⟨-iΔt ∂_{p_i}H δψ', λ⟩
                # T2: ⟨-iΔt ∂_{p_i}H ψ', ν⟩
                # T3: Σⱼ v_{p_j} ⟨-iΔt ∂²_{p_i p_j}H ψ', λ⟩
                #     where ∂²_{p_i p_j}H is non-zero only when both p_i, p_j map to
                #     u-components (m_i ≥ 1 and m_j ≥ 1). For Δt cross-terms see below.
                #
                # Loop over each parameter, build ∂_{p_i}H · vec into tmp_complex, dot.

                # Param ordering: 1..u_dim = uₖ (bw=onemτ),
                #                 u_dim+1..2u_dim = uₖ₊₁ (bw=τ),
                #                 2u_dim+1 = Δt, 2u_dim+2 = t

                # #342: `::Int` on the bound, same reason as the order-3 branch
                # below — `n_params` is captured non-concretely and an unannotated
                # `1:n_params` makes the loop dynamic and allocating per RHS call.
                # (The order-1 branch's remaining `@view`-based state blocks are a
                # recorded follow-up: only the CUBIC cell is in the committed
                # allocation artifact.)
                @inbounds for i = 1:(n_params::Int)
                    # Determine (m_i, bw_i, is_dt) for parameter i
                    if i <= u_dim
                        m_i = i
                        bw_i = onemτ
                        is_dt_i = false
                    elseif i <= 2u_dim
                        m_i = i - u_dim
                        bw_i = τ
                        is_dt_i = false
                    elseif i == 2u_dim + 1
                        m_i = 0
                        bw_i = 0.0  # not used; Δt path handled specially
                        is_dt_i = true
                    else  # i == 2u_dim + 2 (t slot)
                        dg[i] = zero(ComplexF64)
                        continue
                    end

                    # Build (∂_{p_i}H) δψ' into tmp1
                    fill!(tmp1, 0)
                    if is_dt_i
                        # ∂(-iΔtH)/∂Δt = -iH; the corresponding integrand prefactor
                        # already incorporates Δt^0 (no Δt scaling on tmp). The
                        # outer "-iΔt" prefactor adds Δt back. To match Phase 3's
                        # sign convention we accumulate H · vec into tmp (not -i·H).
                        apply!(tmp1, drift_op, δψ′, one(ComplexF64), true)
                        for t_idx = 1:n_terms
                            c = drive_coeff(drives[t_idx], u_interp)
                            apply!(tmp1, drive_ops[t_idx], δψ′, ComplexF64(c), true)
                        end
                    else
                        # ∂_{p_i}H = bw_i · Σ_drives c_jac(d, u, m_i) · H_d
                        for t_idx in control_to_drives[m_i]
                            dc = drive_coeff_jac(drives[t_idx], u_interp, m_i)
                            dc == 0.0 && continue
                            apply!(tmp1, drive_ops[t_idx], δψ′, bw_i * dc, true)
                        end
                    end
                    acc1 = zero(ComplexF64)
                    for k = 1:n
                        acc1 += conj(tmp1[k]) * λ[k]
                    end

                    # Build (∂_{p_i}H) ψ' into tmp1, dot with ν → acc2
                    fill!(tmp1, 0)
                    if is_dt_i
                        apply!(tmp1, drift_op, ψ′, one(ComplexF64), true)
                        for t_idx = 1:n_terms
                            c = drive_coeff(drives[t_idx], u_interp)
                            apply!(tmp1, drive_ops[t_idx], ψ′, ComplexF64(c), true)
                        end
                    else
                        for t_idx in control_to_drives[m_i]
                            dc = drive_coeff_jac(drives[t_idx], u_interp, m_i)
                            dc == 0.0 && continue
                            apply!(tmp1, drive_ops[t_idx], ψ′, bw_i * dc, true)
                        end
                    end
                    acc2 = zero(ComplexF64)
                    for k = 1:n
                        acc2 += conj(tmp1[k]) * ν[k]
                    end

                    # Build (Σⱼ v_{p_j} ∂²_{p_i p_j}H) ψ' into tmp2, dot with λ → acc3.
                    #
                    # ∂²H/∂p_i∂p_j has three sources:
                    # (a) Both p_i, p_j map to u-components (m_i, m_j ≥ 1):
                    #     ∂²H/∂p_i∂p_j = bw_i·bw_j·Σ_d c_hess(d, u, m_i, m_j)·H_d
                    # (b) p_i is a u-component, p_j = Δt (or vice versa):
                    #     ∂²(-iΔtH)/∂u·∂Δt = -i·∂_u H = -i·bw·Σ c_jac · H_d
                    #     The outer "-iΔt" prefactor on this term should become
                    #     "-i" (one Δt cancelled by the differentiation). To keep
                    #     all terms multiplied by the SAME outer "-iΔt" factor,
                    #     we compute the partial as `bw·Σ c_jac·H_d / Δt`. Below
                    #     we instead absorb the cross-term into the integrand
                    #     directly (apply with weight (1/Δt) so the outer -iΔt
                    #     yields -i·∂_u H). When Δt>0 this is well-defined;
                    #     symbolically: ∂²(-iΔtH)/∂u_i∂Δt = -i·bw·c_jac·H_d, and
                    #     to express this as `-iΔt · X` we set X = bw·c_jac·H_d/Δt.
                    # (c) Both p_i, p_j = Δt: ∂²(-iΔtH)/∂Δt² = 0.
                    fill!(tmp2, 0)
                    any_t3 = false
                    if !is_dt_i
                        # m_i ≥ 1: contribution from (a) all m_j-mapped params,
                        # and (b) cross-term with the Δt parameter slot.
                        # (a) bw_i·bw_j·c_hess
                        # uₖ → m_j = j_mod, bw_j = onemτ, ode_idx = j_mod
                        # uₖ₊₁ → m_j = j_mod, bw_j = τ, ode_idx = u_dim + j_mod
                        for j = 1:u_dim
                            v_pj = p[v_p_offset+j]
                            v_pj == 0.0 && continue
                            bw_j = onemτ
                            for t_idx = 1:n_terms
                                d2c = drive_coeff_hess(drives[t_idx], u_interp, m_i, j)
                                d2c == 0.0 && continue
                                apply!(
                                    tmp2,
                                    drive_ops[t_idx],
                                    ψ′,
                                    bw_i * bw_j * d2c * v_pj,
                                    true,
                                )
                                any_t3 = true
                            end
                        end
                        for j = 1:u_dim
                            v_pj = p[v_p_offset+u_dim+j]
                            v_pj == 0.0 && continue
                            bw_j = τ
                            for t_idx = 1:n_terms
                                d2c = drive_coeff_hess(drives[t_idx], u_interp, m_i, j)
                                d2c == 0.0 && continue
                                apply!(
                                    tmp2,
                                    drive_ops[t_idx],
                                    ψ′,
                                    bw_i * bw_j * d2c * v_pj,
                                    true,
                                )
                                any_t3 = true
                            end
                        end
                        # (b) cross with Δt: ∂²(-iΔtH)/∂u_i∂Δt = -i·bw_i·c_jac·H_d.
                        # Outer "-iΔt" prefactor wants to scale `tmp2` to give
                        # this. We need apply weight = bw_i·c_jac/Δt so that
                        # -iΔt · (bw_i·c_jac/Δt) = -i·bw_i·c_jac. Stable since
                        # Δt > 0 in well-posed problems.
                        let v_p_dt = p[v_p_offset+2u_dim+1]
                            if v_p_dt != 0.0 && Δtₖ != 0.0
                                inv_dt = 1.0 / Δtₖ
                                for t_idx in control_to_drives[m_i]
                                    dc = drive_coeff_jac(drives[t_idx], u_interp, m_i)
                                    dc == 0.0 && continue
                                    apply!(
                                        tmp2,
                                        drive_ops[t_idx],
                                        ψ′,
                                        bw_i * dc * v_p_dt * inv_dt,
                                        true,
                                    )
                                    any_t3 = true
                                end
                            end
                        end
                    else  # is_dt_i: contribution only from (b) cross with u-params
                        # ∂²(-iΔtH)/∂Δt∂u_j = -i·bw_j·c_jac·H_d. Same scaling
                        # trick as above to keep the outer "-iΔt" prefactor.
                        if Δtₖ != 0.0
                            inv_dt = 1.0 / Δtₖ
                            for j = 1:u_dim
                                v_pj = p[v_p_offset+j]
                                v_pj == 0.0 && continue
                                bw_j = onemτ
                                for t_idx in control_to_drives[j]
                                    dc = drive_coeff_jac(drives[t_idx], u_interp, j)
                                    dc == 0.0 && continue
                                    apply!(
                                        tmp2,
                                        drive_ops[t_idx],
                                        ψ′,
                                        bw_j * dc * v_pj * inv_dt,
                                        true,
                                    )
                                    any_t3 = true
                                end
                            end
                            for j = 1:u_dim
                                v_pj = p[v_p_offset+u_dim+j]
                                v_pj == 0.0 && continue
                                bw_j = τ
                                for t_idx in control_to_drives[j]
                                    dc = drive_coeff_jac(drives[t_idx], u_interp, j)
                                    dc == 0.0 && continue
                                    apply!(
                                        tmp2,
                                        drive_ops[t_idx],
                                        ψ′,
                                        bw_j * dc * v_pj * inv_dt,
                                        true,
                                    )
                                    any_t3 = true
                                end
                            end
                        end
                    end

                    acc3 = zero(ComplexF64)
                    if any_t3
                        for k = 1:n
                            acc3 += conj(tmp2[k]) * λ[k]
                        end
                    end

                    # Per-slot prefactor:
                    # - u-component params: `-iΔt` (one Δt factor for `∂_u(-iΔtH) = -iΔt·∂_uH`)
                    # - Δt param: `-i`        (the Δt factor was already differentiated out:
                    #                         `∂_Δt(-iΔtH) = -iH`)
                    # The T3 contributions for Δt-slot already had a `1/Δt` scaling
                    # applied at apply!-time so the outer `-iΔt` factor produces the
                    # correct `-i·∂_uH` form for cross-Hessian. Now we shift to a
                    # `-i` outer factor for the Δt slot, and the inner `1/Δt` scaling
                    # cancels — re-scale T3 by `Δt` for is_dt_i=true to compensate.
                    if is_dt_i
                        # For Δt slot: T3's inner 1/Δt was tuned for outer -iΔt;
                        # outer is now -i, so multiply T3 contributions by Δt to
                        # restore the right magnitude.
                        dg[i] = -im * (acc1 + acc2 + Δtₖ * acc3)
                    else
                        dg[i] = -im * Δtₖ * (acc1 + acc2 + acc3)
                    end
                end

                return nothing
            end
        end

    else  # order == 3
        () -> begin
            u_interp = zeros(u_dim)
            tmp1 = Vector{ComplexF64}(undef, n)
            tmp2 = Vector{ComplexF64}(undef, n)
            # #353: closure scratch for the cubic Δt chain rule. `du_dΔt[m]` is
            # D_m = ∂u_m/∂Δt; `dc_du_dΔt[t]` is Σ_m ∂c_t/∂u_m · D_m, the drive-level
            # contraction that ∂A/∂Δt needs and that is reused by the δψ' forcing,
            # the ν forcing, both Δt-slot integrand terms and the (Δt, Δt) diagonal.
            # Preallocated per factory call, like every other buffer here (#342).
            du_dΔt = zeros(u_dim)
            dc_du_dΔt = zeros(n_terms)
            # #342: CONCRETE per-block buffers, not `@view`s of `x`/`dx`.
            #
            # This RHS is called ~35 times per per-knot backward solve, so
            # anything it allocates per call becomes a PER-KNOT term and breaks
            # the knot-flat contract (ADR-0009 decision 3). Measured on the
            # #337 cell (ketdim 2, 2 controls, cubic): the nine `@view`s below
            # cost 720 B per call — 24 832 B per knot, which is the whole of
            # this cell's residual growth once the per-knot `remake`/`solve`
            # was replaced by a preallocated integrator. `SubArray`s of the
            # solver-owned state escape into `apply!`, whose operator argument
            # comes out of a `Vector` of wrappers, so they cannot be elided;
            # copying the four state blocks in and the four rate blocks out
            # (8·n complex assignments, n = ketdim) is strictly cheaper and
            # allocation-free. The `g` block is written straight into `dx`.
            #
            # NOTHING about the second-order adjoint math changes here.
            ψb = Vector{ComplexF64}(undef, n)
            δψb = Vector{ComplexF64}(undef, n)
            λb = Vector{ComplexF64}(undef, n)
            νb = Vector{ComplexF64}(undef, n)
            dψb = Vector{ComplexF64}(undef, n)
            dδψb = Vector{ComplexF64}(undef, n)
            dλb = Vector{ComplexF64}(undef, n)
            dνb = Vector{ComplexF64}(undef, n)

            (dx, x, p, τ) -> begin
                Δtₖ = p[4u_dim+1]
                # tₖ  = p[4u_dim + 2]  # unused
                v_p_offset = 4u_dim + 2

                τ2 = τ * τ
                τ3 = τ2 * τ
                h00 = 2τ3 - 3τ2 + 1
                h10 = (τ3 - 2τ2 + τ) * Δtₖ
                h01 = -2τ3 + 3τ2
                h11 = (τ3 - τ2) * Δtₖ
                # Note: cubic du-parameter weights carry a Δt factor; we
                # treat (h10, h11) as the effective basis weights at runtime.
                # #353: ∂h10/∂Δt and ∂h11/∂Δt are the un-scaled Hermite tangent
                # bases — the `Δt`-derivative of the `du` basis weights.
                dh10_dΔt = τ3 - 2τ2 + τ
                dh11_dΔt = τ3 - τ2

                @inbounds for i = 1:u_dim
                    duₖᵢ = p[2u_dim+i]
                    duₖ₊₁ᵢ = p[3u_dim+i]
                    u_interp[i] =
                        h00 * p[i] + h10 * duₖᵢ + h01 * p[u_dim+i] + h11 * duₖ₊₁ᵢ
                    # D_i = ∂u_i/∂Δt, holding the knot parameters fixed.
                    du_dΔt[i] = dh10_dΔt * duₖᵢ + dh11_dΔt * duₖ₊₁ᵢ
                end

                # #353: Σ_m (∂c_t/∂u_m)·D_m per drive. This is the drive-level
                # weight of the Δt chain rule through u(τ; Δt); ∂A/∂Δt is
                # -i·H(u) - iΔt·Σ_t dc_du_dΔt[t]·H_t.
                @inbounds for t_idx = 1:n_terms
                    acc_dt = 0.0
                    for m in drive_controls[t_idx]
                        Dm = du_dΔt[m]
                        Dm == 0.0 && continue
                        acc_dt += drive_coeff_jac(drives[t_idx], u_interp, m) * Dm
                    end
                    dc_du_dΔt[t_idx] = acc_dt
                end

                ψ′ = ψb
                δψ′ = δψb
                λ = λb
                ν = νb
                @inbounds for i = 1:n
                    ψ′[i] = x[i]
                    δψ′[i] = x[n+i]
                    λ[i] = x[2n+i]
                    ν[i] = x[3n+i]
                end
                dψ′ = dψb
                dδψ′ = dδψb
                dλ = dλb
                dν = dνb
                # `dg[i]` writes go straight into `dx` at the `g` block offset.
                dg_offset = 4n

                # ── ψ'̇ = -iΔt H ψ' ──
                apply!(dψ′, drift_op, ψ′, -im * Δtₖ, false)
                @inbounds for t_idx = 1:n_terms
                    c = drive_coeff(drives[t_idx], u_interp)
                    apply!(dψ′, drive_ops[t_idx], ψ′, -im * Δtₖ * c, true)
                end

                # ── δψ'̇ = -iΔt H δψ' - iΔt M_v ψ' ──
                # Homogeneous part
                apply!(dδψ′, drift_op, δψ′, -im * Δtₖ, false)
                @inbounds for t_idx = 1:n_terms
                    c = drive_coeff(drives[t_idx], u_interp)
                    apply!(dδψ′, drive_ops[t_idx], δψ′, -im * Δtₖ * c, true)
                end
                # Forcing M_v ψ' : uₖ, uₖ₊₁, duₖ, duₖ₊₁ contributions
                _bws = (h00, h01, h10, h11)
                # #353: ∂bw/∂Δt per slot — zero for the (uₖ, uₖ₊₁) value weights,
                # the un-scaled Hermite tangent bases for the (duₖ, duₖ₊₁) ones.
                _dbws = (0.0, 0.0, dh10_dΔt, dh11_dΔt)
                _offsets = (0, u_dim, 2u_dim, 3u_dim)
                @inbounds for slot = 1:4
                    bw_slot = _bws[slot]
                    off = _offsets[slot]
                    for j = 1:u_dim
                        v_pj = p[v_p_offset+off+j]
                        v_pj == 0.0 && continue
                        for t_idx in control_to_drives[j]
                            dc = drive_coeff_jac(drives[t_idx], u_interp, j)
                            if dc != 0.0
                                apply!(
                                    dδψ′,
                                    drive_ops[t_idx],
                                    ψ′,
                                    -im * Δtₖ * bw_slot * dc * v_pj,
                                    true,
                                )
                            end
                        end
                    end
                end
                # Forcing M_v ψ' : Δt contribution.
                # #353: ∂A/∂Δt = -i·H(u) - iΔt·Σ_t (Σ_m ∂c_t/∂u_m·D_m)·H_t. The
                # second part is the cubic chain rule through the `du` basis
                # weights; `build_hvp_forward_ode` has always carried it (this is
                # the forward tangent's own Δt column) and this builder dropped it,
                # so the round-tripped δψ' diverged from the forward δψ.
                let v_pj = p[v_p_offset+4u_dim+1]
                    if v_pj != 0.0
                        apply!(dδψ′, drift_op, ψ′, -im * v_pj, true)
                        @inbounds for t_idx = 1:n_terms
                            c = drive_coeff(drives[t_idx], u_interp)
                            apply!(dδψ′, drive_ops[t_idx], ψ′, -im * c * v_pj, true)
                        end
                        @inbounds for t_idx = 1:n_terms
                            gdt = dc_du_dΔt[t_idx]
                            gdt == 0.0 && continue
                            apply!(dδψ′, drive_ops[t_idx], ψ′, -im * Δtₖ * gdt * v_pj, true)
                        end
                    end
                end

                # ── λ̇ = -iΔt H λ ──
                apply!(dλ, drift_op, λ, -im * Δtₖ, false)
                @inbounds for t_idx = 1:n_terms
                    c = drive_coeff(drives[t_idx], u_interp)
                    apply!(dλ, drive_ops[t_idx], λ, -im * Δtₖ * c, true)
                end

                # ── ν̇ = -iΔt H ν - iΔt M_v λ ──
                apply!(dν, drift_op, ν, -im * Δtₖ, false)
                @inbounds for t_idx = 1:n_terms
                    c = drive_coeff(drives[t_idx], u_interp)
                    apply!(dν, drive_ops[t_idx], ν, -im * Δtₖ * c, true)
                end
                @inbounds for slot = 1:4
                    bw_slot = _bws[slot]
                    off = _offsets[slot]
                    for j = 1:u_dim
                        v_pj = p[v_p_offset+off+j]
                        v_pj == 0.0 && continue
                        for t_idx in control_to_drives[j]
                            dc = drive_coeff_jac(drives[t_idx], u_interp, j)
                            if dc != 0.0
                                apply!(
                                    dν,
                                    drive_ops[t_idx],
                                    λ,
                                    -im * Δtₖ * bw_slot * dc * v_pj,
                                    true,
                                )
                            end
                        end
                    end
                end
                # #353: same Δt chain-rule part of ∂A/∂Δt on the ν forcing. This is
                # the term the (x, p) cross-Hessian block rides on — ν(0) is the
                # ψ̃-slot output — so dropping it also biased the state block.
                let v_pj = p[v_p_offset+4u_dim+1]
                    if v_pj != 0.0
                        apply!(dν, drift_op, λ, -im * v_pj, true)
                        @inbounds for t_idx = 1:n_terms
                            c = drive_coeff(drives[t_idx], u_interp)
                            apply!(dν, drive_ops[t_idx], λ, -im * c * v_pj, true)
                        end
                        @inbounds for t_idx = 1:n_terms
                            gdt = dc_du_dΔt[t_idx]
                            gdt == 0.0 && continue
                            apply!(dν, drive_ops[t_idx], λ, -im * Δtₖ * gdt * v_pj, true)
                        end
                    end
                end

                # ── ġ_i = three-term integrand ──
                # Param mapping for cubic:
                #   1..u_dim         → uₖ,    bw = h00
                #   u_dim+1..2u_dim  → uₖ₊₁,  bw = h01
                #   2u_dim+1..3u_dim → duₖ,   bw = h10
                #   3u_dim+1..4u_dim → duₖ₊₁, bw = h11
                #   4u_dim+1         → Δt
                #   4u_dim+2         → t
                # #342: THE LOOP BOUND IS ANNOTATED, AND THAT IS LOAD-BEARING.
                #
                # `n_params` is captured from the enclosing builder, where it is
                # assigned inside an `if order == …` chain; the closure sees it
                # as non-concrete (`code_warntype` reports `Union{}`), so
                # `1:n_params` made this whole loop dynamic. Measured cost on the
                # #337 cell: 352 B per RHS call → ~12 kB per knot at ~35 calls,
                # a per-knot term that breaks the knot-flat contract (ADR-0009
                # decision 3). With `::Int` it is 0 B. `is_dt_i` joins the
                # `local … ::T` declaration for the same reason — dropping the
                # declarations, or replacing the if/elseif chain with table
                # reads, both measured WORSE (3 680 B/call).
                @inbounds for i = 1:(n_params::Int)
                    # #353: `dbw_i` = ∂bw_i/∂Δt joins the `local … ::T` declaration
                    # for the same type-stability reason as the rest.
                    local m_i::Int, bw_i::Float64, dbw_i::Float64, is_dt_i::Bool
                    is_dt_i = false
                    dbw_i = 0.0
                    if i <= u_dim
                        m_i = i
                        bw_i = h00
                    elseif i <= 2u_dim
                        m_i = i - u_dim
                        bw_i = h01
                    elseif i <= 3u_dim
                        m_i = i - 2u_dim
                        bw_i = h10
                        dbw_i = dh10_dΔt
                    elseif i <= 4u_dim
                        m_i = i - 3u_dim
                        bw_i = h11
                        dbw_i = dh11_dΔt
                    elseif i == 4u_dim + 1
                        m_i = 0
                        bw_i = 0.0
                        is_dt_i = true
                    else
                        dx[dg_offset+i] = zero(ComplexF64)
                        continue
                    end

                    # T1: ⟨-iΔt ∂_{p_i}H δψ', λ⟩
                    fill!(tmp1, 0)
                    if is_dt_i
                        apply!(tmp1, drift_op, δψ′, one(ComplexF64), true)
                        for t_idx = 1:n_terms
                            c = drive_coeff(drives[t_idx], u_interp)
                            apply!(tmp1, drive_ops[t_idx], δψ′, ComplexF64(c), true)
                        end
                        # #353: ∂A/∂Δt also carries Δt·Σ_t (Σ_m ∂c_t/∂u_m·D_m)·H_t.
                        # The Δt-slot prefactor is `-i`, so the Δt factor is explicit.
                        for t_idx = 1:n_terms
                            gdt = dc_du_dΔt[t_idx]
                            gdt == 0.0 && continue
                            apply!(tmp1, drive_ops[t_idx], δψ′, ComplexF64(Δtₖ * gdt), true)
                        end
                    else
                        for t_idx in control_to_drives[m_i]
                            dc = drive_coeff_jac(drives[t_idx], u_interp, m_i)
                            dc == 0.0 && continue
                            apply!(tmp1, drive_ops[t_idx], δψ′, bw_i * dc, true)
                        end
                    end
                    acc1 = zero(ComplexF64)
                    for k = 1:n
                        acc1 += conj(tmp1[k]) * λ[k]
                    end

                    # T2: ⟨-iΔt ∂_{p_i}H ψ', ν⟩
                    fill!(tmp1, 0)
                    if is_dt_i
                        apply!(tmp1, drift_op, ψ′, one(ComplexF64), true)
                        for t_idx = 1:n_terms
                            c = drive_coeff(drives[t_idx], u_interp)
                            apply!(tmp1, drive_ops[t_idx], ψ′, ComplexF64(c), true)
                        end
                        # #353: the Δt chain-rule part of ∂A/∂Δt, as in T1.
                        for t_idx = 1:n_terms
                            gdt = dc_du_dΔt[t_idx]
                            gdt == 0.0 && continue
                            apply!(tmp1, drive_ops[t_idx], ψ′, ComplexF64(Δtₖ * gdt), true)
                        end
                    else
                        for t_idx in control_to_drives[m_i]
                            dc = drive_coeff_jac(drives[t_idx], u_interp, m_i)
                            dc == 0.0 && continue
                            apply!(tmp1, drive_ops[t_idx], ψ′, bw_i * dc, true)
                        end
                    end
                    acc2 = zero(ComplexF64)
                    for k = 1:n
                        acc2 += conj(tmp1[k]) * ν[k]
                    end

                    # T3: Σⱼ v_{p_j} ⟨-iΔt ∂²_{p_i p_j}H ψ', λ⟩
                    fill!(tmp2, 0)
                    any_t3 = false
                    if !is_dt_i
                        # u-component / u-component cross terms
                        for slot_j = 1:4
                            bw_jj = _bws[slot_j]
                            off_j = _offsets[slot_j]
                            for j = 1:u_dim
                                v_pj = p[v_p_offset+off_j+j]
                                v_pj == 0.0 && continue
                                for t_idx = 1:n_terms
                                    d2c = drive_coeff_hess(drives[t_idx], u_interp, m_i, j)
                                    d2c == 0.0 && continue
                                    apply!(
                                        tmp2,
                                        drive_ops[t_idx],
                                        ψ′,
                                        bw_i * bw_jj * d2c * v_pj,
                                        true,
                                    )
                                    any_t3 = true
                                end
                            end
                        end
                        # Cross with Δt. #353: ALL THREE paths of ∂²A/∂p_i∂Δt, not
                        # just (a). The outer prefactor on `tmp2` for a non-Δt slot
                        # is `-iΔt`, so each term's weight is its own coefficient
                        # divided by that:
                        #   (a) -i·bw_i·c_jac        → bw_i·c_jac/Δt
                        #   (c) -iΔt·(∂bw_i/∂Δt)·c_jac → (∂bw_i/∂Δt)·c_jac
                        #   (b) -iΔt·bw_i·Σ_m c_hess(m_i,m)·D_m → bw_i·Σ_m c_hess·D_m
                        # (c) is the term the `du` basis weights' explicit Δtₖ factor
                        # generates and (b) the one u's own Δt dependence generates;
                        # both were missing, which is why the `du` block was off by
                        # 4.5e-1 while the pure-control blocks were nearly right.
                        let v_p_dt = p[v_p_offset+4u_dim+1]
                            if v_p_dt != 0.0
                                a_c_scale = (Δtₖ == 0.0 ? 0.0 : bw_i / Δtₖ) + dbw_i
                                if a_c_scale != 0.0
                                    for t_idx in control_to_drives[m_i]
                                        dc = drive_coeff_jac(drives[t_idx], u_interp, m_i)
                                        dc == 0.0 && continue
                                        apply!(
                                            tmp2,
                                            drive_ops[t_idx],
                                            ψ′,
                                            a_c_scale * dc * v_p_dt,
                                            true,
                                        )
                                        any_t3 = true
                                    end
                                end
                                if bw_i != 0.0
                                    for t_idx = 1:n_terms
                                        s = 0.0
                                        for cm in drive_controls[t_idx]
                                            Dm = du_dΔt[cm]
                                            Dm == 0.0 && continue
                                            d2c = drive_coeff_hess(
                                                drives[t_idx],
                                                u_interp,
                                                m_i,
                                                cm,
                                            )
                                            d2c == 0.0 && continue
                                            s += d2c * Dm
                                        end
                                        s == 0.0 && continue
                                        apply!(
                                            tmp2,
                                            drive_ops[t_idx],
                                            ψ′,
                                            bw_i * s * v_p_dt,
                                            true,
                                        )
                                        any_t3 = true
                                    end
                                end
                            end
                        end
                    else
                        # is_dt_i. Outer prefactor on `tmp2` here is `-iΔt` too
                        # (`dg = -i(acc1 + acc2 + Δt·acc3)`), so the same weights as
                        # above apply, transposed onto p_j.
                        for slot_j = 1:4
                            bw_jj = _bws[slot_j]
                            dbw_jj = _dbws[slot_j]
                            off_j = _offsets[slot_j]
                            a_c_scale = (Δtₖ == 0.0 ? 0.0 : bw_jj / Δtₖ) + dbw_jj
                            for j = 1:u_dim
                                v_pj = p[v_p_offset+off_j+j]
                                v_pj == 0.0 && continue
                                if a_c_scale != 0.0
                                    for t_idx in control_to_drives[j]
                                        dc = drive_coeff_jac(drives[t_idx], u_interp, j)
                                        dc == 0.0 && continue
                                        apply!(
                                            tmp2,
                                            drive_ops[t_idx],
                                            ψ′,
                                            a_c_scale * dc * v_pj,
                                            true,
                                        )
                                        any_t3 = true
                                    end
                                end
                                if bw_jj != 0.0
                                    for t_idx = 1:n_terms
                                        s = 0.0
                                        for cm in drive_controls[t_idx]
                                            Dm = du_dΔt[cm]
                                            Dm == 0.0 && continue
                                            d2c = drive_coeff_hess(
                                                drives[t_idx],
                                                u_interp,
                                                j,
                                                cm,
                                            )
                                            d2c == 0.0 && continue
                                            s += d2c * Dm
                                        end
                                        s == 0.0 && continue
                                        apply!(
                                            tmp2,
                                            drive_ops[t_idx],
                                            ψ′,
                                            bw_jj * s * v_pj,
                                            true,
                                        )
                                        any_t3 = true
                                    end
                                end
                            end
                        end
                        # #353: the (Δt, Δt) DIAGONAL, which the cubic branch had
                        # nowhere. It vanishes for a linear spline (u carries no Δt
                        # dependence there) but not for cubic:
                        #   ∂²A/∂Δt² = -2i·Σ_t (Σ_m ∂c_t/∂u_m·D_m)·H_t
                        #              -iΔt·Σ_t (Σ_{m,m'} ∂²c_t/∂u_m∂u_{m'} D_m D_{m'})·H_t
                        # Against the `-iΔt` outer prefactor the weight is
                        # 2·dc_du_dΔt/Δt + Σ_{m,m'} c_hess·D_m·D_{m'}.
                        let v_p_dt = p[v_p_offset+4u_dim+1]
                            if v_p_dt != 0.0
                                inv_dt = Δtₖ == 0.0 ? 0.0 : 1.0 / Δtₖ
                                for t_idx = 1:n_terms
                                    w = 2.0 * dc_du_dΔt[t_idx] * inv_dt
                                    for cm in drive_controls[t_idx]
                                        Dm = du_dΔt[cm]
                                        Dm == 0.0 && continue
                                        for cn in drive_controls[t_idx]
                                            Dn = du_dΔt[cn]
                                            Dn == 0.0 && continue
                                            d2c = drive_coeff_hess(
                                                drives[t_idx],
                                                u_interp,
                                                cm,
                                                cn,
                                            )
                                            d2c == 0.0 && continue
                                            w += d2c * Dm * Dn
                                        end
                                    end
                                    w == 0.0 && continue
                                    apply!(tmp2, drive_ops[t_idx], ψ′, w * v_p_dt, true)
                                    any_t3 = true
                                end
                            end
                        end
                    end

                    acc3 = zero(ComplexF64)
                    if any_t3
                        for k = 1:n
                            acc3 += conj(tmp2[k]) * λ[k]
                        end
                    end

                    # Per-slot prefactor (see linear branch for explanation):
                    # u-comp: `-iΔt`; Δt-slot: `-i`. T3 1/Δt scaling already in.
                    if is_dt_i
                        dx[dg_offset+i] = -im * (acc1 + acc2 + Δtₖ * acc3)
                    else
                        dx[dg_offset+i] = -im * Δtₖ * (acc1 + acc2 + acc3)
                    end
                end

                # Copy the four rate blocks out of the concrete buffers (#342).
                @inbounds for i = 1:n
                    dx[i] = dψ′[i]
                    dx[n+i] = dδψ′[i]
                    dx[2n+i] = dλ[i]
                    dx[3n+i] = dν[i]
                end

                return nothing
            end
        end
    end

    return make_f!_2nd_adj, n_params
end

"""
    build_ket_sensitivity_ode(drift_op, drives, u_dim, ketdim, order, n_kets)

Build an ODE that propagates K kets and their parameter sensitivities as vectors.

Instead of propagating n×n sensitivity matrices Sⱼ = ∂Φ/∂pⱼ, this propagates
K n-dim ket sensitivity vectors sᵢⱼ(τ) = Sⱼ(τ)·ψᵢ(0), which satisfy:
    dsᵢⱼ/dτ = -iΔt·H·sᵢⱼ + forcing_ij(ψᵢ)

The forcing terms are identical to the matrix case but applied to ψᵢ instead of Φ.

State layout: [ψ₁, ..., ψ_K, s₁₁, ..., s₁J, s₂₁, ..., s_{K,J}]
Total: K·n·(1+J) complex values, vs n²·(1+J) for the matrix case.

Speedup: n/K when K < n (e.g., 3.75× for n=30, K=8).

Returns `(make_f!, n_params)`.
"""
function build_ket_sensitivity_ode(
    drift_op::AbstractDynamicsOperator,
    drives::AbstractVector{<:AbstractDrive},
    u_dim::Int,
    ketdim::Int,
    order::Int,
    n_kets::Int,
)
    n_terms = length(drives)

    if order == 1
        n_params = 2 * u_dim + 2  # [uₖ, uₖ₊₁, Δt, t]
    elseif order == 3
        n_params = 4 * u_dim + 2  # [uₖ, uₖ₊₁, duₖ, duₖ₊₁, Δt, t]
    else
        error("Unsupported spline order: $order")
    end

    if isempty(drives)
        return nothing, n_params
    end

    K = n_kets
    n = ketdim

    # Pre-build operator wrappers
    drive_ops = [_to_operator(d.H) for d in drives]

    # Control→drives mapping (same as build_sensitivity_ode)
    control_to_drives = [Int[] for _ = 1:u_dim]
    drives_active = [active_controls(d) for d in drives]
    for (t_idx, ac) in enumerate(drives_active)
        if isempty(ac)
            for j = 1:u_dim
                push!(control_to_drives[j], t_idx)
            end
        else
            for j in ac
                push!(control_to_drives[j], t_idx)
            end
        end
    end

    make_f!_ket_sens = if order == 1
        () -> begin
            u_interp = zeros(u_dim)
            jac_cache = zeros(u_dim, n_terms)
            (dx, x, p, τ) -> begin
                uₖ = @view p[1:u_dim]
                uₖ₊₁ = @view p[(u_dim+1):2u_dim]
                Δtₖ = p[end-1]

                onemτ = 1 - τ
                @inbounds for i = 1:u_dim
                    u_interp[i] = onemτ * uₖ[i] + τ * uₖ₊₁[i]
                end

                # Op-fusion: per-substep jac coefficient cache.
                # jac_cache[j, t_idx] = ∂c_t/∂u_j; non-active entries left at 0.
                # Mirrors build_ket_jvp_ode (line 1842).
                @inbounds for t_idx = 1:n_terms
                    ac = drives_active[t_idx]
                    iter = isempty(ac) ? (1:u_dim) : ac
                    for j = 1:u_dim
                        jac_cache[j, t_idx] = 0.0
                    end
                    for j in iter
                        jac_cache[j, t_idx] = drive_coeff_jac(drives[t_idx], u_interp, j)
                    end
                end

                # Apply -iΔt·H(u(τ)) to a vector
                @inline function apply_H_vec!(dy, y)
                    apply!(dy, drift_op, y, -im * Δtₖ, false)
                    @inbounds for t_idx = 1:n_terms
                        c = drive_coeff(drives[t_idx], u_interp)
                        apply!(dy, drive_ops[t_idx], y, -im * Δtₖ * c, true)
                    end
                end

                # --- Propagate K kets: dψᵢ/dτ = -iΔt·H·ψᵢ ---
                @inbounds for i = 1:K
                    off = (i - 1) * n
                    ψᵢ = @view x[(off+1):(off+n)]
                    dψᵢ = @view dx[(off+1):(off+n)]
                    apply_H_vec!(dψᵢ, ψᵢ)
                end

                # --- Ket sensitivities ---
                # Layout after kets: s_{i,j} at offset K*n + ((i-1)*n_params + (j-1))*n
                ket_block = K * n

                # Sensitivities for uₖ[j]: basis = (1-τ)
                @inbounds for j = 1:u_dim
                    for i = 1:K
                        off_s = ket_block + ((i - 1) * n_params + (j - 1)) * n
                        sᵢⱼ = @view x[(off_s+1):(off_s+n)]
                        dsᵢⱼ = @view dx[(off_s+1):(off_s+n)]
                        apply_H_vec!(dsᵢⱼ, sᵢⱼ)
                        # Forcing: (∂c/∂u_j)·(1-τ)·H_d·ψᵢ
                        off_ψ = (i - 1) * n
                        ψᵢ = @view x[(off_ψ+1):(off_ψ+n)]
                        for t_idx in control_to_drives[j]
                            dc = jac_cache[j, t_idx]
                            if dc != 0.0
                                apply!(dsᵢⱼ, drive_ops[t_idx], ψᵢ, -im * Δtₖ * onemτ * dc, true)
                            end
                        end
                    end
                end

                # Sensitivities for uₖ₊₁[j]: basis = τ
                @inbounds for j = 1:u_dim
                    for i = 1:K
                        j_param = u_dim + j
                        off_s = ket_block + ((i - 1) * n_params + (j_param - 1)) * n
                        sᵢⱼ = @view x[(off_s+1):(off_s+n)]
                        dsᵢⱼ = @view dx[(off_s+1):(off_s+n)]
                        apply_H_vec!(dsᵢⱼ, sᵢⱼ)
                        off_ψ = (i - 1) * n
                        ψᵢ = @view x[(off_ψ+1):(off_ψ+n)]
                        for t_idx in control_to_drives[j]
                            dc = jac_cache[j, t_idx]
                            if dc != 0.0
                                apply!(dsᵢⱼ, drive_ops[t_idx], ψᵢ, -im * Δtₖ * τ * dc, true)
                            end
                        end
                    end
                end

                # Sensitivity for Δt: forcing = -iH(u)·ψᵢ (no Δt factor)
                @inbounds for i = 1:K
                    j_param = 2 * u_dim + 1
                    off_s = ket_block + ((i - 1) * n_params + (j_param - 1)) * n
                    sᵢⱼ = @view x[(off_s+1):(off_s+n)]
                    dsᵢⱼ = @view dx[(off_s+1):(off_s+n)]
                    apply_H_vec!(dsᵢⱼ, sᵢⱼ)
                    off_ψ = (i - 1) * n
                    ψᵢ = @view x[(off_ψ+1):(off_ψ+n)]
                    apply!(dsᵢⱼ, drift_op, ψᵢ, -im, true)
                    for t_idx = 1:n_terms
                        c = drive_coeff(drives[t_idx], u_interp)
                        apply!(dsᵢⱼ, drive_ops[t_idx], ψᵢ, -im * c, true)
                    end
                end

                # Sensitivity for t (zero forcing for time-independent)
                @inbounds for i = 1:K
                    j_param = 2 * u_dim + 2
                    off_s = ket_block + ((i - 1) * n_params + (j_param - 1)) * n
                    sᵢⱼ = @view x[(off_s+1):(off_s+n)]
                    dsᵢⱼ = @view dx[(off_s+1):(off_s+n)]
                    apply_H_vec!(dsᵢⱼ, sᵢⱼ)
                end

                return nothing
            end
        end

    else  # order == 3
        () -> begin
            u_interp = zeros(u_dim)
            jac_cache = zeros(u_dim, n_terms)
            (dx, x, p, τ) -> begin
                uₖ = @view p[1:u_dim]
                uₖ₊₁ = @view p[(u_dim+1):2u_dim]
                duₖ = @view p[(2u_dim+1):3u_dim]
                duₖ₊₁ = @view p[(3u_dim+1):4u_dim]
                Δtₖ = p[end-1]

                τ2 = τ * τ
                τ3 = τ2 * τ
                h00 = 2τ3 - 3τ2 + 1
                h10 = (τ3 - 2τ2 + τ) * Δtₖ
                h01 = -2τ3 + 3τ2
                h11 = (τ3 - τ2) * Δtₖ
                dh10_dΔt = τ3 - 2τ2 + τ
                dh11_dΔt = τ3 - τ2

                @inbounds for i = 1:u_dim
                    u_interp[i] =
                        h00 * uₖ[i] + h10 * duₖ[i] + h01 * uₖ₊₁[i] + h11 * duₖ₊₁[i]
                end

                # Op-fusion: per-substep jac coefficient cache.
                # jac_cache[j, t_idx] = ∂c_t/∂u_j; non-active entries left at 0.
                # Mirrors build_ket_jvp_ode (line 1842).
                @inbounds for t_idx = 1:n_terms
                    ac = drives_active[t_idx]
                    iter = isempty(ac) ? (1:u_dim) : ac
                    for j = 1:u_dim
                        jac_cache[j, t_idx] = 0.0
                    end
                    for j in iter
                        jac_cache[j, t_idx] = drive_coeff_jac(drives[t_idx], u_interp, j)
                    end
                end

                @inline function apply_H_vec!(dy, y)
                    apply!(dy, drift_op, y, -im * Δtₖ, false)
                    @inbounds for t_idx = 1:n_terms
                        c = drive_coeff(drives[t_idx], u_interp)
                        apply!(dy, drive_ops[t_idx], y, -im * Δtₖ * c, true)
                    end
                end

                ket_block = K * n

                # Propagate K kets
                @inbounds for i = 1:K
                    off = (i - 1) * n
                    ψᵢ = @view x[(off+1):(off+n)]
                    dψᵢ = @view dx[(off+1):(off+n)]
                    apply_H_vec!(dψᵢ, ψᵢ)
                end

                # Sensitivities for uₖ[j]: cubic basis = h00
                @inbounds for j = 1:u_dim
                    for i = 1:K
                        off_s = ket_block + ((i - 1) * n_params + (j - 1)) * n
                        sᵢⱼ = @view x[(off_s+1):(off_s+n)]
                        dsᵢⱼ = @view dx[(off_s+1):(off_s+n)]
                        apply_H_vec!(dsᵢⱼ, sᵢⱼ)
                        off_ψ = (i - 1) * n
                        ψᵢ = @view x[(off_ψ+1):(off_ψ+n)]
                        for t_idx in control_to_drives[j]
                            dc = jac_cache[j, t_idx]
                            if dc != 0.0
                                apply!(dsᵢⱼ, drive_ops[t_idx], ψᵢ, -im * Δtₖ * h00 * dc, true)
                            end
                        end
                    end
                end

                # Sensitivities for uₖ₊₁[j]: basis = h01
                @inbounds for j = 1:u_dim
                    for i = 1:K
                        j_param = u_dim + j
                        off_s = ket_block + ((i - 1) * n_params + (j_param - 1)) * n
                        sᵢⱼ = @view x[(off_s+1):(off_s+n)]
                        dsᵢⱼ = @view dx[(off_s+1):(off_s+n)]
                        apply_H_vec!(dsᵢⱼ, sᵢⱼ)
                        off_ψ = (i - 1) * n
                        ψᵢ = @view x[(off_ψ+1):(off_ψ+n)]
                        for t_idx in control_to_drives[j]
                            dc = jac_cache[j, t_idx]
                            if dc != 0.0
                                apply!(dsᵢⱼ, drive_ops[t_idx], ψᵢ, -im * Δtₖ * h01 * dc, true)
                            end
                        end
                    end
                end

                # Sensitivities for duₖ[j]: basis = h10
                @inbounds for j = 1:u_dim
                    for i = 1:K
                        j_param = 2 * u_dim + j
                        off_s = ket_block + ((i - 1) * n_params + (j_param - 1)) * n
                        sᵢⱼ = @view x[(off_s+1):(off_s+n)]
                        dsᵢⱼ = @view dx[(off_s+1):(off_s+n)]
                        apply_H_vec!(dsᵢⱼ, sᵢⱼ)
                        off_ψ = (i - 1) * n
                        ψᵢ = @view x[(off_ψ+1):(off_ψ+n)]
                        for t_idx in control_to_drives[j]
                            dc = jac_cache[j, t_idx]
                            if dc != 0.0
                                apply!(dsᵢⱼ, drive_ops[t_idx], ψᵢ, -im * Δtₖ * h10 * dc, true)
                            end
                        end
                    end
                end

                # Sensitivities for duₖ₊₁[j]: basis = h11
                @inbounds for j = 1:u_dim
                    for i = 1:K
                        j_param = 3 * u_dim + j
                        off_s = ket_block + ((i - 1) * n_params + (j_param - 1)) * n
                        sᵢⱼ = @view x[(off_s+1):(off_s+n)]
                        dsᵢⱼ = @view dx[(off_s+1):(off_s+n)]
                        apply_H_vec!(dsᵢⱼ, sᵢⱼ)
                        off_ψ = (i - 1) * n
                        ψᵢ = @view x[(off_ψ+1):(off_ψ+n)]
                        for t_idx in control_to_drives[j]
                            dc = jac_cache[j, t_idx]
                            if dc != 0.0
                                apply!(dsᵢⱼ, drive_ops[t_idx], ψᵢ, -im * Δtₖ * h11 * dc, true)
                            end
                        end
                    end
                end

                # Δt sensitivity
                @inbounds for i = 1:K
                    j_param = 4 * u_dim + 1
                    off_s = ket_block + ((i - 1) * n_params + (j_param - 1)) * n
                    sᵢⱼ = @view x[(off_s+1):(off_s+n)]
                    dsᵢⱼ = @view dx[(off_s+1):(off_s+n)]
                    apply_H_vec!(dsᵢⱼ, sᵢⱼ)
                    off_ψ = (i - 1) * n
                    ψᵢ = @view x[(off_ψ+1):(off_ψ+n)]
                    # Forcing: -iH·ψ + -iΔt·Σ(∂c/∂Δt)·H_d·ψ
                    apply!(dsᵢⱼ, drift_op, ψᵢ, -im, true)
                    for t_idx = 1:n_terms
                        c = drive_coeff(drives[t_idx], u_interp)
                        apply!(dsᵢⱼ, drive_ops[t_idx], ψᵢ, -im * c, true)
                        # ∂u/∂Δt contribution via Hermite basis
                        for jj = 1:u_dim
                            dc = jac_cache[jj, t_idx]
                            if dc != 0.0
                                du_dDt = dh10_dΔt * duₖ[jj] + dh11_dΔt * duₖ₊₁[jj]
                                apply!(
                                    dsᵢⱼ,
                                    drive_ops[t_idx],
                                    ψᵢ,
                                    -im * Δtₖ * dc * du_dDt,
                                    true,
                                )
                            end
                        end
                    end
                end

                # t sensitivity (zero forcing)
                @inbounds for i = 1:K
                    j_param = 4 * u_dim + 2
                    off_s = ket_block + ((i - 1) * n_params + (j_param - 1)) * n
                    sᵢⱼ = @view x[(off_s+1):(off_s+n)]
                    dsᵢⱼ = @view dx[(off_s+1):(off_s+n)]
                    apply_H_vec!(dsᵢⱼ, sᵢⱼ)
                end

                return nothing
            end
        end
    end

    return make_f!_ket_sens, n_params
end

"""
    build_ket_sensitivity_problems(make_f!, n_params, ketdim, n_kets, p₀, N)

Build ODE problems for ket-level sensitivity. Initial conditions are set to
zero (kets are populated via `remake` at each knot point).
"""
function build_ket_sensitivity_problems(make_f!, n_params, ketdim, n_kets, p₀, N)
    K = n_kets
    n = ketdim
    state_dim = K * n * (1 + n_params)
    x₀ = zeros(ComplexF64, state_dim)
    probs = [ODEProblem(make_f!(), x₀, (0.0, 1.0), p₀) for _ = 1:N]
    return probs, state_dim
end

"""
    tri_idx(i, j, n) -> Int

Upper-triangle linear index for (i,j) with 1 ≤ i ≤ j ≤ n.
Maps to sequential 1-based index: (1,1)→1, (1,2)→2, ..., (1,n)→n, (2,2)→n+1, ...
"""
@inline function tri_idx(i::Int, j::Int, n::Int)
    return (i - 1) * n - ((i - 1) * (i - 2)) ÷ 2 + (j - i + 1)
end

"""
    build_second_order_sensitivity_ode(drift_op, drives, u_dim, ketdim, order)

Build the second-order sensitivity ODE for exact Hessian computation.

The augmented state layout is `[Φ; S₁...Sₙ; T₁₁, T₁₂, ..., Tₙₙ]` where:
- `Φ` is the propagator matrix (ketdim² complex elements)
- `Sⱼ` is the first-order sensitivity ∂Φ/∂pⱼ
- `T_{ij}` (upper triangle, i ≤ j) is the second-order sensitivity ∂²Φ/∂pᵢ∂pⱼ

The ODE for T_{ij} is:
    dT_{ij}/dτ = A·T_{ij} + (∂A/∂pᵢ)·Sⱼ + (∂A/∂pⱼ)·Sᵢ + (∂²A/∂pᵢ∂pⱼ)·Ψ

where A = -iΔt·H(u(τ)), and the derivatives of A involve `drive_coeff_jac` and
`drive_coeff_hess`.

Pairs involving the `t` parameter (time-independent systems) have zero forcing and
are skipped.

Returns `(make_f!, n_params, active_pairs, hess_state_dim)` where:
- `make_f!`: factory `() -> (dx, x, p, τ) -> nothing` producing thread-safe closures
- `n_params`: number of ODE parameters per knot point
- `active_pairs`: Vector{Tuple{Int,Int}} of (i,j) parameter index pairs being tracked
- `hess_state_dim`: total state dimension (complex elements)
"""
function build_second_order_sensitivity_ode(
    drift_op::AbstractDynamicsOperator,
    drives::AbstractVector{<:AbstractDrive},
    u_dim::Int,
    ketdim::Int,
    order::Int,
)
    Φ_dim = ketdim^2
    n_terms = length(drives)

    if order == 1
        n_params = 2 * u_dim + 2  # [uₖ, uₖ₊₁, Δt, t]
    elseif order == 3
        n_params = 4 * u_dim + 2  # [uₖ, uₖ₊₁, duₖ, duₖ₊₁, Δt, t]
    else
        error("Unsupported spline order: $order")
    end

    if isempty(drives)
        return nothing, n_params, Tuple{Int,Int}[], 0
    end

    drive_ops = [_to_operator(d.H) for d in drives]

    # Precompute control→drives mapping (same as first-order)
    control_to_drives = [Int[] for _ = 1:u_dim]
    drives_active = [active_controls(d) for d in drives]
    for (t_idx, ac) in enumerate(drives_active)
        if isempty(ac)
            for j = 1:u_dim
                push!(control_to_drives[j], t_idx)
            end
        else
            for j in ac
                push!(control_to_drives[j], t_idx)
            end
        end
    end

    # Precompute control pair → drives mapping for ∂²c/∂uᵢ∂uⱼ (Hessian sparsity)
    # hess_pair_to_drives[(ci, cj)] = [drive indices with potentially nonzero Hessian]
    hess_pair_to_drives = Dict{Tuple{Int,Int},Vector{Int}}()
    for ci = 1:u_dim, cj = ci:u_dim
        relevant = Int[]
        for (t_idx, ac) in enumerate(drives_active)
            if isempty(ac) || (ci ∈ ac && cj ∈ ac)
                push!(relevant, t_idx)
            end
        end
        if !isempty(relevant)
            hess_pair_to_drives[(ci, cj)] = relevant
        end
    end

    # Build active pairs list (skip t parameter which has zero forcing)
    t_param_idx = n_params  # t is always the last parameter
    active_pairs = Tuple{Int,Int}[]
    for i = 1:n_params, j = i:n_params
        if i == t_param_idx || j == t_param_idx
            continue  # skip t pairs (zero forcing for time-independent)
        end
        push!(active_pairs, (i, j))
    end
    n_active = length(active_pairs)

    # Build reverse lookup: (i,j) → index in active_pairs
    pair_to_idx = Dict{Tuple{Int,Int},Int}()
    for (k, (i, j)) in enumerate(active_pairs)
        pair_to_idx[(i, j)] = k
    end

    # State layout: [Φ; S₁...Sₙ; T_active₁...T_activeₘ]
    first_order_dim = Φ_dim * (1 + n_params)
    hess_state_dim = first_order_dim + n_active * Φ_dim

    # Helper: which control index does parameter p_idx map to? (0 = not a control)
    # And what is its spline basis function (returned as a function of τ, Δt)?
    # For linear spline:
    #   p 1..u_dim       → ctrl 1..u_dim (uₖ)
    #   p u_dim+1..2u_dim → ctrl 1..u_dim (uₖ₊₁)
    #   p 2u_dim+1       → Δt (special)
    #   p 2u_dim+2       → t (skipped)

    Δt_idx = n_params - 1  # Δt parameter index

    make_f! = if order == 1
        () -> begin
            u_interp = zeros(u_dim)

            (dx, x, p, τ) -> begin
                uₖ = @view p[1:u_dim]
                uₖ₊₁ = @view p[(u_dim+1):2u_dim]
                Δtₖ = p[end-1]

                onemτ = 1 - τ
                @inbounds for i = 1:u_dim
                    u_interp[i] = onemτ * uₖ[i] + τ * uₖ₊₁[i]
                end

                # === Apply H(u(τ)) to a matrix: drift + Σ coeff(d, u_τ) * d_op ===
                @inline function apply_H!(dM, M)
                    apply!(dM, drift_op, M, -im * Δtₖ, false)
                    @inbounds for t_idx = 1:n_terms
                        c = drive_coeff(drives[t_idx], u_interp)
                        apply!(dM, drive_ops[t_idx], M, -im * Δtₖ * c, true)
                    end
                end

                # === Apply ∂A/∂p_i forcing to matrix M (accumulate into dM) ===
                # For control params: ∂A/∂(uₖ[j]) = -iΔt · Σ_d dc_d · (1-τ) · H_d
                # For Δt: ∂A/∂Δt = -i · H(u)
                @inline function apply_dA_dp!(dM, M, p_idx)
                    if p_idx <= u_dim
                        # uₖ[j]: basis = (1-τ)
                        j = p_idx
                        for t_idx in control_to_drives[j]
                            dc = drive_coeff_jac(drives[t_idx], u_interp, j)
                            if dc != 0.0
                                apply!(dM, drive_ops[t_idx], M, -im * Δtₖ * onemτ * dc, true)
                            end
                        end
                    elseif p_idx <= 2 * u_dim
                        # uₖ₊₁[j]: basis = τ
                        j = p_idx - u_dim
                        for t_idx in control_to_drives[j]
                            dc = drive_coeff_jac(drives[t_idx], u_interp, j)
                            if dc != 0.0
                                apply!(dM, drive_ops[t_idx], M, -im * Δtₖ * τ * dc, true)
                            end
                        end
                    elseif p_idx == Δt_idx
                        # Δt: ∂A/∂Δt = -i·H(u)
                        apply!(dM, drift_op, M, -im, true)
                        @inbounds for t_idx = 1:n_terms
                            c = drive_coeff(drives[t_idx], u_interp)
                            apply!(dM, drive_ops[t_idx], M, -im * c, true)
                        end
                    end
                    # t parameter: zero forcing (do nothing)
                end

                # === Apply ∂²A/∂pᵢ∂pⱼ forcing to Ψ (accumulate into dT) ===
                @inline function apply_d2A_dp!(dT, Ψ, p_i, p_j)
                    # Get control indices and basis values for both parameters
                    ci, bi = if p_i <= u_dim
                        p_i, onemτ
                    elseif p_i <= 2 * u_dim
                        p_i - u_dim, τ
                    else
                        0, 0.0  # Δt or t
                    end
                    cj, bj = if p_j <= u_dim
                        p_j, onemτ
                    elseif p_j <= 2 * u_dim
                        p_j - u_dim, τ
                    else
                        0, 0.0
                    end

                    if ci > 0 && cj > 0
                        # Both are control params: ∂²A = -iΔt · Σ_d d²c · bᵢ · bⱼ · H_d
                        ci_lo, cj_lo = min(ci, cj), max(ci, cj)
                        pair_key = (ci_lo, cj_lo)
                        if haskey(hess_pair_to_drives, pair_key)
                            for t_idx in hess_pair_to_drives[pair_key]
                                d2c = drive_coeff_hess(drives[t_idx], u_interp, ci, cj)
                                if d2c != 0.0
                                    apply!(
                                        dT,
                                        drive_ops[t_idx],
                                        Ψ,
                                        -im * Δtₖ * bi * bj * d2c,
                                        true,
                                    )
                                end
                            end
                        end
                    elseif ci > 0 && p_j == Δt_idx
                        # One control, one Δt: ∂²A/∂(u[ci])∂Δt = -i · Σ_d dc · bᵢ · H_d
                        for t_idx in control_to_drives[ci]
                            dc = drive_coeff_jac(drives[t_idx], u_interp, ci)
                            if dc != 0.0
                                apply!(dT, drive_ops[t_idx], Ψ, -im * bi * dc, true)
                            end
                        end
                    elseif cj > 0 && p_i == Δt_idx
                        # Symmetric: Δt first, control second
                        for t_idx in control_to_drives[cj]
                            dc = drive_coeff_jac(drives[t_idx], u_interp, cj)
                            if dc != 0.0
                                apply!(dT, drive_ops[t_idx], Ψ, -im * bj * dc, true)
                            end
                        end
                    end
                    # Δt-Δt: ∂²A/∂Δt² = 0 for linear spline (no Δt dependence in u)
                    # t-anything: zero (skipped by active_pairs)
                end

                # ==================== Forward propagation ====================
                Φ_vec = @view x[1:Φ_dim]
                Φ_mat = reshape(Φ_vec, ketdim, ketdim)
                dΦ_vec = @view dx[1:Φ_dim]
                dΦ_mat = reshape(dΦ_vec, ketdim, ketdim)
                apply_H!(dΦ_mat, Φ_mat)

                # ==================== First-order sensitivities ====================
                @inbounds for j = 1:n_params
                    offset = Φ_dim + (j - 1) * Φ_dim
                    Sⱼ_mat = reshape(@view(x[(offset+1):(offset+Φ_dim)]), ketdim, ketdim)
                    dSⱼ_mat = reshape(@view(dx[(offset+1):(offset+Φ_dim)]), ketdim, ketdim)
                    apply_H!(dSⱼ_mat, Sⱼ_mat)
                    apply_dA_dp!(dSⱼ_mat, Φ_mat, j)
                end

                # ==================== Second-order sensitivities ====================
                @inbounds for (pair_k, (pi, pj)) in enumerate(active_pairs)
                    T_offset = first_order_dim + (pair_k - 1) * Φ_dim
                    Tij_mat = reshape(
                        @view(x[(T_offset+1):(T_offset+Φ_dim)]),
                        ketdim,
                        ketdim,
                    )
                    dTij_mat = reshape(
                        @view(dx[(T_offset+1):(T_offset+Φ_dim)]),
                        ketdim,
                        ketdim,
                    )

                    # Homogeneous term: A·T_{ij}
                    apply_H!(dTij_mat, Tij_mat)

                    # Cross terms: (∂A/∂pᵢ)·Sⱼ + (∂A/∂pⱼ)·Sᵢ
                    Sj_offset = Φ_dim + (pj - 1) * Φ_dim
                    Sⱼ_mat = reshape(
                        @view(x[(Sj_offset+1):(Sj_offset+Φ_dim)]),
                        ketdim,
                        ketdim,
                    )
                    apply_dA_dp!(dTij_mat, Sⱼ_mat, pi)

                    Si_offset = Φ_dim + (pi - 1) * Φ_dim
                    Sᵢ_mat = reshape(
                        @view(x[(Si_offset+1):(Si_offset+Φ_dim)]),
                        ketdim,
                        ketdim,
                    )
                    apply_dA_dp!(dTij_mat, Sᵢ_mat, pj)

                    # Second-order forcing: (∂²A/∂pᵢ∂pⱼ)·Ψ
                    apply_d2A_dp!(dTij_mat, Φ_mat, pi, pj)
                end

                return nothing
            end
        end

    else  # order == 3
        # Cubic spline second-order sensitivity ODE
        # Parameter layout: [uₖ, uₖ₊₁, duₖ, duₖ₊₁, Δt, t]
        () -> begin
            u_interp = zeros(u_dim)

            (dx, x, p, τ) -> begin
                uₖ = @view p[1:u_dim]
                uₖ₊₁ = @view p[(u_dim+1):2u_dim]
                duₖ = @view p[(2u_dim+1):3u_dim]
                duₖ₊₁ = @view p[(3u_dim+1):4u_dim]
                Δtₖ = p[end-1]

                τ2 = τ * τ
                τ3 = τ2 * τ
                h00 = 2τ3 - 3τ2 + 1
                h10 = (τ3 - 2τ2 + τ) * Δtₖ
                h01 = -2τ3 + 3τ2
                h11 = (τ3 - τ2) * Δtₖ
                dh10_dΔt = τ3 - 2τ2 + τ
                dh11_dΔt = τ3 - τ2

                @inbounds for i = 1:u_dim
                    u_interp[i] =
                        h00 * uₖ[i] + h10 * duₖ[i] + h01 * uₖ₊₁[i] + h11 * duₖ₊₁[i]
                end

                # Basis function for each parameter group
                # Group 0: uₖ    → h00
                # Group 1: uₖ₊₁  → h01
                # Group 2: duₖ   → h10
                # Group 3: duₖ₊₁ → h11
                @inline function param_basis(p_idx)
                    if p_idx <= u_dim
                        return h00
                    elseif p_idx <= 2 * u_dim
                        return h01
                    elseif p_idx <= 3 * u_dim
                        return h10
                    elseif p_idx <= 4 * u_dim
                        return h11
                    else
                        return 0.0
                    end
                end

                @inline function param_ctrl(p_idx)
                    if p_idx <= u_dim
                        return p_idx
                    elseif p_idx <= 2 * u_dim
                        return p_idx - u_dim
                    elseif p_idx <= 3 * u_dim
                        return p_idx - 2 * u_dim
                    elseif p_idx <= 4 * u_dim
                        return p_idx - 3 * u_dim
                    else
                        return 0
                    end
                end

                @inline function apply_H!(dM, M)
                    apply!(dM, drift_op, M, -im * Δtₖ, false)
                    @inbounds for t_idx = 1:n_terms
                        c = drive_coeff(drives[t_idx], u_interp)
                        apply!(dM, drive_ops[t_idx], M, -im * Δtₖ * c, true)
                    end
                end

                @inline function apply_dA_dp!(dM, M, p_idx)
                    ci = param_ctrl(p_idx)
                    if ci > 0
                        bi = param_basis(p_idx)
                        for t_idx in control_to_drives[ci]
                            dc = drive_coeff_jac(drives[t_idx], u_interp, ci)
                            if dc != 0.0
                                apply!(dM, drive_ops[t_idx], M, -im * Δtₖ * bi * dc, true)
                            end
                        end
                    elseif p_idx == Δt_idx
                        # ∂A/∂Δt = -i·H(u) + (-iΔt)·Σ_d ∂c/∂u_m · ∂u_m/∂Δt · H_d
                        apply!(dM, drift_op, M, -im, true)
                        @inbounds for t_idx = 1:n_terms
                            c = drive_coeff(drives[t_idx], u_interp)
                            apply!(dM, drive_ops[t_idx], M, -im * c, true)
                        end
                        # Chain rule: u depends on Δt through h10, h11
                        @inbounds for t_idx = 1:n_terms
                            op = drive_ops[t_idx]
                            ac = drives_active[t_idx]
                            iter = isempty(ac) ? (1:u_dim) : ac
                            total_du_dΔt = 0.0
                            for i in iter
                                dc = drive_coeff_jac(drives[t_idx], u_interp, i)
                                if dc != 0.0
                                    du_i_dΔt = dh10_dΔt * duₖ[i] + dh11_dΔt * duₖ₊₁[i]
                                    total_du_dΔt += dc * du_i_dΔt
                                end
                            end
                            if total_du_dΔt != 0.0
                                apply!(dM, op, M, -im * Δtₖ * total_du_dΔt, true)
                            end
                        end
                    end
                end

                @inline function apply_d2A_dp!(dT, Ψ, p_i, p_j)
                    ci = param_ctrl(p_i)
                    cj = param_ctrl(p_j)
                    bi = param_basis(p_i)
                    bj = param_basis(p_j)

                    if ci > 0 && cj > 0
                        ci_lo, cj_lo = min(ci, cj), max(ci, cj)
                        pair_key = (ci_lo, cj_lo)
                        if haskey(hess_pair_to_drives, pair_key)
                            for t_idx in hess_pair_to_drives[pair_key]
                                d2c = drive_coeff_hess(drives[t_idx], u_interp, ci, cj)
                                if d2c != 0.0
                                    apply!(
                                        dT,
                                        drive_ops[t_idx],
                                        Ψ,
                                        -im * Δtₖ * bi * bj * d2c,
                                        true,
                                    )
                                end
                            end
                        end
                    elseif ci > 0 && p_j == Δt_idx
                        # Term (a): -i · bi · dc/du_ci · H_d (from d(-iΔt)/dΔt = -i)
                        for t_idx in control_to_drives[ci]
                            dc = drive_coeff_jac(drives[t_idx], u_interp, ci)
                            if dc != 0.0
                                apply!(dT, drive_ops[t_idx], Ψ, -im * bi * dc, true)
                            end
                        end
                        # Term (b): chain rule for d²c through u-Δt coupling
                        for t_idx = 1:n_terms
                            ac = drives_active[t_idx]
                            iter = isempty(ac) ? (1:u_dim) : ac
                            for cm in iter
                                d2c = drive_coeff_hess(drives[t_idx], u_interp, ci, cm)
                                if d2c != 0.0
                                    du_m_dΔt = dh10_dΔt * duₖ[cm] + dh11_dΔt * duₖ₊₁[cm]
                                    if du_m_dΔt != 0.0
                                        apply!(
                                            dT,
                                            drive_ops[t_idx],
                                            Ψ,
                                            -im * Δtₖ * bi * du_m_dΔt * d2c,
                                            true,
                                        )
                                    end
                                end
                            end
                        end
                        # Term (c): basis function depends on Δt (h10, h11 ∝ Δt)
                        # ∂bi/∂Δt is nonzero only for du groups
                        dbi_dΔt = if p_i <= 2 * u_dim
                            0.0  # uₖ (h00) or uₖ₊₁ (h01) don't depend on Δt
                        elseif p_i <= 3 * u_dim
                            dh10_dΔt  # duₖ group
                        elseif p_i <= 4 * u_dim
                            dh11_dΔt  # duₖ₊₁ group
                        else
                            0.0
                        end
                        if dbi_dΔt != 0.0
                            for t_idx in control_to_drives[ci]
                                dc = drive_coeff_jac(drives[t_idx], u_interp, ci)
                                if dc != 0.0
                                    apply!(
                                        dT,
                                        drive_ops[t_idx],
                                        Ψ,
                                        -im * Δtₖ * dbi_dΔt * dc,
                                        true,
                                    )
                                end
                            end
                        end
                    elseif cj > 0 && p_i == Δt_idx
                        # Term (a): -i · bj · dc/du_cj · H_d
                        for t_idx in control_to_drives[cj]
                            dc = drive_coeff_jac(drives[t_idx], u_interp, cj)
                            if dc != 0.0
                                apply!(dT, drive_ops[t_idx], Ψ, -im * bj * dc, true)
                            end
                        end
                        # Term (b): chain rule for d²c through u-Δt coupling
                        for t_idx = 1:n_terms
                            ac = drives_active[t_idx]
                            iter = isempty(ac) ? (1:u_dim) : ac
                            for cm in iter
                                d2c = drive_coeff_hess(drives[t_idx], u_interp, cj, cm)
                                if d2c != 0.0
                                    du_m_dΔt = dh10_dΔt * duₖ[cm] + dh11_dΔt * duₖ₊₁[cm]
                                    if du_m_dΔt != 0.0
                                        apply!(
                                            dT,
                                            drive_ops[t_idx],
                                            Ψ,
                                            -im * Δtₖ * bj * du_m_dΔt * d2c,
                                            true,
                                        )
                                    end
                                end
                            end
                        end
                        # Term (c): basis function depends on Δt
                        dbj_dΔt = if p_j <= 2 * u_dim
                            0.0
                        elseif p_j <= 3 * u_dim
                            dh10_dΔt
                        elseif p_j <= 4 * u_dim
                            dh11_dΔt
                        else
                            0.0
                        end
                        if dbj_dΔt != 0.0
                            for t_idx in control_to_drives[cj]
                                dc = drive_coeff_jac(drives[t_idx], u_interp, cj)
                                if dc != 0.0
                                    apply!(
                                        dT,
                                        drive_ops[t_idx],
                                        Ψ,
                                        -im * Δtₖ * dbj_dΔt * dc,
                                        true,
                                    )
                                end
                            end
                        end
                    elseif p_i == Δt_idx && p_j == Δt_idx
                        # #353: the (Δt, Δt) DIAGONAL. It is identically zero for a
                        # linear spline — u carries no Δt dependence there — which is
                        # why "typically small; omit" survived. For cubic it is not
                        # small: with D_m = ∂u_m/∂Δt = h̃10·duₖ[m] + h̃11·duₖ₊₁[m],
                        #   ∂A/∂Δt   = -i·H(u) - iΔt·Σ_t (Σ_m ∂c_t/∂u_m·D_m)·H_t
                        #   ∂²A/∂Δt² = -2i·Σ_t (Σ_m ∂c_t/∂u_m·D_m)·H_t
                        #              -iΔt·Σ_t (Σ_{m,m'} ∂²c_t/∂u_m∂u_{m'} D_m D_{m'})·H_t
                        # The leading factor is 2, not 1: one copy from
                        # ∂/∂Δt of -i·H(u(Δt)) and one from ∂/∂Δt of the explicit Δt.
                        # Measured cost of omitting it on the Ket Tsit5 cubic cell:
                        # 1.5e-2 relative error in the Δt block of the assembled
                        # exact Hessian against a central FD of the residual gradient.
                        @inbounds for t_idx = 1:n_terms
                            ac = drives_active[t_idx]
                            iter = isempty(ac) ? (1:u_dim) : ac
                            g1 = 0.0
                            q = 0.0
                            for cm in iter
                                Dm = dh10_dΔt * duₖ[cm] + dh11_dΔt * duₖ₊₁[cm]
                                Dm == 0.0 && continue
                                g1 += drive_coeff_jac(drives[t_idx], u_interp, cm) * Dm
                                for cn in iter
                                    Dn = dh10_dΔt * duₖ[cn] + dh11_dΔt * duₖ₊₁[cn]
                                    Dn == 0.0 && continue
                                    d2c = drive_coeff_hess(drives[t_idx], u_interp, cm, cn)
                                    d2c == 0.0 && continue
                                    q += d2c * Dm * Dn
                                end
                            end
                            w = -2im * g1 - im * Δtₖ * q
                            w == 0.0 && continue
                            apply!(dT, drive_ops[t_idx], Ψ, w, true)
                        end
                    end
                    # t-anything: zero (skipped by active_pairs)
                end

                # ==================== Forward propagation ====================
                Φ_vec = @view x[1:Φ_dim]
                Φ_mat = reshape(Φ_vec, ketdim, ketdim)
                dΦ_vec = @view dx[1:Φ_dim]
                dΦ_mat = reshape(dΦ_vec, ketdim, ketdim)
                apply_H!(dΦ_mat, Φ_mat)

                # ==================== First-order sensitivities ====================
                @inbounds for j = 1:n_params
                    offset = Φ_dim + (j - 1) * Φ_dim
                    Sⱼ_mat = reshape(@view(x[(offset+1):(offset+Φ_dim)]), ketdim, ketdim)
                    dSⱼ_mat = reshape(@view(dx[(offset+1):(offset+Φ_dim)]), ketdim, ketdim)
                    apply_H!(dSⱼ_mat, Sⱼ_mat)
                    apply_dA_dp!(dSⱼ_mat, Φ_mat, j)
                end

                # ==================== Second-order sensitivities ====================
                @inbounds for (pair_k, (pi, pj)) in enumerate(active_pairs)
                    T_offset = first_order_dim + (pair_k - 1) * Φ_dim
                    Tij_mat = reshape(
                        @view(x[(T_offset+1):(T_offset+Φ_dim)]),
                        ketdim,
                        ketdim,
                    )
                    dTij_mat = reshape(
                        @view(dx[(T_offset+1):(T_offset+Φ_dim)]),
                        ketdim,
                        ketdim,
                    )

                    apply_H!(dTij_mat, Tij_mat)

                    Sj_offset = Φ_dim + (pj - 1) * Φ_dim
                    Sⱼ_mat = reshape(
                        @view(x[(Sj_offset+1):(Sj_offset+Φ_dim)]),
                        ketdim,
                        ketdim,
                    )
                    apply_dA_dp!(dTij_mat, Sⱼ_mat, pi)

                    Si_offset = Φ_dim + (pi - 1) * Φ_dim
                    Sᵢ_mat = reshape(
                        @view(x[(Si_offset+1):(Si_offset+Φ_dim)]),
                        ketdim,
                        ketdim,
                    )
                    apply_dA_dp!(dTij_mat, Sᵢ_mat, pj)

                    apply_d2A_dp!(dTij_mat, Φ_mat, pi, pj)
                end

                return nothing
            end
        end
    end

    return make_f!, n_params, active_pairs, hess_state_dim
end

"""
    build_forward_ode(drift_op, drives, u_dim, order)

Build the forward ODE right-hand side from drive terms.

Returns a **factory** `() -> (dx, x, p, τ) -> nothing`. Each call to the factory
produces an independent closure with its own workspace buffers, ensuring thread
safety when multiple ODE problems are solved in parallel via `Threads.@threads`.
Works with both vector (ket) and matrix (unitary/propagator) state `x`.
"""
function build_forward_ode(
    drift_op::AbstractDynamicsOperator,
    drives::AbstractVector{<:AbstractDrive},
    u_dim::Int,
    order::Int,
)
    n_terms = length(drives)
    drive_ops = [_to_operator(d.H) for d in drives]

    if order == 1
        return () -> begin
            u_interp = zeros(u_dim)
            (dx, x, p, τ) -> begin
                uₖ = @view p[1:u_dim]
                uₖ₊₁ = @view p[(u_dim+1):2u_dim]
                Δtₖ = p[end-1]

                onemτ = 1 - τ
                @inbounds for i = 1:u_dim
                    u_interp[i] = onemτ * uₖ[i] + τ * uₖ₊₁[i]
                end

                apply!(dx, drift_op, x, -im * Δtₖ, false)
                @inbounds for t_idx = 1:n_terms
                    c = drive_coeff(drives[t_idx], u_interp)
                    apply!(dx, drive_ops[t_idx], x, -im * Δtₖ * c, true)
                end
                return nothing
            end
        end

    else  # order == 3
        return () -> begin
            u_interp = zeros(u_dim)
            (dx, x, p, τ) -> begin
                uₖ = @view p[1:u_dim]
                uₖ₊₁ = @view p[(u_dim+1):2u_dim]
                duₖ = @view p[(2u_dim+1):3u_dim]
                duₖ₊₁ = @view p[(3u_dim+1):4u_dim]
                Δtₖ = p[end-1]

                τ2 = τ * τ
                τ3 = τ2 * τ
                h00 = 2τ3 - 3τ2 + 1
                h10 = (τ3 - 2τ2 + τ) * Δtₖ
                h01 = -2τ3 + 3τ2
                h11 = (τ3 - τ2) * Δtₖ

                @inbounds for i = 1:u_dim
                    u_interp[i] =
                        h00 * uₖ[i] + h10 * duₖ[i] + h01 * uₖ₊₁[i] + h11 * duₖ₊₁[i]
                end

                apply!(dx, drift_op, x, -im * Δtₖ, false)
                @inbounds for t_idx = 1:n_terms
                    c = drive_coeff(drives[t_idx], u_interp)
                    apply!(dx, drive_ops[t_idx], x, -im * Δtₖ * c, true)
                end
                return nothing
            end
        end
    end
end

"""
    build_hamiltonian_jacobian!(jac_buf, drift_op, drive_ops, drives, u_dim, ketdim, order)

Build a closure that materializes the dense Hamiltonian `H(τ; p)` into `jac_buf`
without scaling. Used to construct closed-form Jacobians for implicit ODE solvers
(e.g. Rodas5P) which require an analytic Jacobian when the state is `Vector{ComplexF64}`
— OrdinaryDiffEq's default `ForwardDiff.DerivativeConfig` rejects complex states.

The returned closure has signature `(τ, p) -> nothing` and writes `H(τ; p)` into
`jac_buf::Matrix{ComplexF64}` of size `(ketdim, ketdim)`. Caller multiplies by
`-im * Δt` to obtain the per-state Jacobian for the forward ODE
`dψ/dτ = -iΔt H(τ;p) ψ`. For the sensitivity ODE the same `H` populates each
diagonal block under the `(I_ketdim ⊗ H)` Kronecker structure of `vec(H·M)`.

The factory mirrors `build_forward_ode`: each call returns an independent closure
with its own `u_interp` workspace for thread safety.
"""
function build_hamiltonian_jacobian!(
    drift_op::AbstractDynamicsOperator,
    drives::AbstractVector{<:AbstractDrive},
    u_dim::Int,
    ketdim::Int,
    order::Int,
)
    n_terms = length(drives)
    drive_ops = [_to_operator(d.H) for d in drives]
    Iket = Matrix{ComplexF64}(I, ketdim, ketdim)

    if order == 1
        return () -> begin
            u_interp = zeros(u_dim)
            (H_buf::AbstractMatrix{ComplexF64}, p::AbstractVector, τ::Real) -> begin
                uₖ = @view p[1:u_dim]
                uₖ₊₁ = @view p[(u_dim+1):2u_dim]

                onemτ = 1 - τ
                @inbounds for i = 1:u_dim
                    u_interp[i] = onemτ * uₖ[i] + τ * uₖ₊₁[i]
                end

                # H = drift + Σ c_i(u) drive_op_i, applied to identity → dense matrix
                apply!(H_buf, drift_op, Iket, true, false)
                @inbounds for t_idx = 1:n_terms
                    c = drive_coeff(drives[t_idx], u_interp)
                    apply!(H_buf, drive_ops[t_idx], Iket, c, true)
                end
                return nothing
            end
        end

    else  # order == 3
        return () -> begin
            u_interp = zeros(u_dim)
            (H_buf::AbstractMatrix{ComplexF64}, p::AbstractVector, τ::Real) -> begin
                uₖ = @view p[1:u_dim]
                uₖ₊₁ = @view p[(u_dim+1):2u_dim]
                duₖ = @view p[(2u_dim+1):3u_dim]
                duₖ₊₁ = @view p[(3u_dim+1):4u_dim]
                Δtₖ = p[end-1]

                τ2 = τ * τ
                τ3 = τ2 * τ
                h00 = 2τ3 - 3τ2 + 1
                h10 = (τ3 - 2τ2 + τ) * Δtₖ
                h01 = -2τ3 + 3τ2
                h11 = (τ3 - τ2) * Δtₖ

                @inbounds for i = 1:u_dim
                    u_interp[i] =
                        h00 * uₖ[i] + h10 * duₖ[i] + h01 * uₖ₊₁[i] + h11 * duₖ₊₁[i]
                end

                apply!(H_buf, drift_op, Iket, true, false)
                @inbounds for t_idx = 1:n_terms
                    c = drive_coeff(drives[t_idx], u_interp)
                    apply!(H_buf, drive_ops[t_idx], Iket, c, true)
                end
                return nothing
            end
        end
    end
end

"""
    build_forward_ode_jacobian(drift_op, drives, u_dim, ketdim, order)

Build a Jacobian factory paired with `build_forward_ode`. Each factory call returns
an independent closure `(J, x, p, τ) -> nothing` that writes `J = -iΔt H(τ; p)`
(dense `ketdim × ketdim` complex matrix). For the linear ket ODE
`dψ/dτ = -iΔt H(τ;p) ψ`, this is the exact Jacobian w.r.t. ψ — independent of ψ
because the ODE is linear in the state.

Used to construct `ODEFunction(f!; jac = jac!)` for implicit solvers (Rodas5P)
that cannot autodiff complex states. Factory pattern (rather than direct closure)
mirrors `build_forward_ode` and ensures each per-knot ODE problem gets its own
workspace closure for thread safety.
"""
function build_forward_ode_jacobian(
    drift_op::AbstractDynamicsOperator,
    drives::AbstractVector{<:AbstractDrive},
    u_dim::Int,
    ketdim::Int,
    order::Int,
)
    make_H! = build_hamiltonian_jacobian!(drift_op, drives, u_dim, ketdim, order)

    return () -> begin
        H_buf = zeros(ComplexF64, ketdim, ketdim)
        H_writer! = make_H!()
        # ODE is linear in state, so the Jacobian doesn't depend on `_x`.
        (J::AbstractMatrix{ComplexF64}, _x::AbstractVector, p::AbstractVector, τ::Real) -> begin
            Δtₖ = p[end-1]
            H_writer!(H_buf, p, τ)
            # J = -i Δt H
            α = -im * Δtₖ
            @inbounds for j = 1:ketdim, i = 1:ketdim
                J[i, j] = α * H_buf[i, j]
            end
            return nothing
        end
    end
end

"""
    build_sensitivity_ode_jacobian(drift_op, drives, u_dim, ketdim, order, n_params)

Build a Jacobian factory for the sensitivity-augmented ODE. The augmented state
is `[Φ_vec; S_1_vec; ...; S_n_vec]` of length `Φ_dim · (1 + n_params)` where
`Φ_dim = ketdim^2`. The Jacobian w.r.t. state is block-diagonal: each block is
`(I_ketdim ⊗ (-iΔt H(τ;p)))`, the Kronecker form of `vec(H·M)`.

The forcing terms for sensitivity equations involve `(∂_p H) Φ` which depends only
on parameters `p` (not on the augmented state), so they do not contribute to the
state-Jacobian. Returns the same factory pattern as `build_forward_ode_jacobian`.

Note: for large `n_params`, the dense block-diagonal Jacobian is wasteful in memory
and Rodas5P factorization cost. A `jac_prototype` with sparse block-diagonal pattern
is the proper optimization but is deferred — Phase 0's goal is correctness, not the
full performance story for Rodas5P.
"""
function build_sensitivity_ode_jacobian(
    drift_op::AbstractDynamicsOperator,
    drives::AbstractVector{<:AbstractDrive},
    u_dim::Int,
    ketdim::Int,
    order::Int,
    n_params::Int,
)
    Φ_dim = ketdim^2
    n_blocks = 1 + n_params
    make_H! = build_hamiltonian_jacobian!(drift_op, drives, u_dim, ketdim, order)

    return () -> begin
        H_buf = zeros(ComplexF64, ketdim, ketdim)
        H_writer! = make_H!()
        # Sensitivity ODE is linear in augmented state, so Jacobian doesn't depend on `_x`.
        (J::AbstractMatrix{ComplexF64}, _x::AbstractVector, p::AbstractVector, τ::Real) -> begin
            Δtₖ = p[end-1]
            H_writer!(H_buf, p, τ)
            α = -im * Δtₖ
            fill!(J, zero(ComplexF64))
            # Each block is (I_ketdim ⊗ (-iΔt H)) of size (Φ_dim, Φ_dim).
            # vec(H · M)[i + (a-1)*ketdim] = Σ_b H[i,b] M[b,a],
            # so block[(i + (a-1)*ketdim), (b + (a-1)*ketdim)] = α * H[i,b].
            @inbounds for blk = 0:(n_blocks-1)
                base = blk * Φ_dim
                for a = 1:ketdim
                    col_base = base + (a - 1) * ketdim
                    row_base = base + (a - 1) * ketdim
                    for b = 1:ketdim, i = 1:ketdim
                        J[row_base+i, col_base+b] = α * H_buf[i, b]
                    end
                end
            end
            return nothing
        end
    end
end

"""
    build_sensitivity_problems(make_f!_sens, n_params, ketdim, ode_p₀, N)

Create ODE problems and preallocated state buffer for the sensitivity-augmented ODE.
The ODE state is complex-valued: Φ ∈ ℂ^{n×n}, Sⱼ ∈ ℂ^{n×n}.

# Arguments
- `make_f!_sens`: Closure factory from `build_sensitivity_ode`; called once per interval
- `n_params::Int`: Number of ODE parameters
- `ketdim::Int`: Ket space dimension
- `ode_p₀::AbstractVector`: Initial parameter vector `[p₀; Δt₀; t₀]`
- `N::Int`: Number of knot points (creates N-1 problems)

# Returns
`(sens_probs, sens_state)` — vector of ODEProblems and preallocated state buffer
"""
function build_sensitivity_problems(
    make_f!_sens,
    n_params::Int,
    ketdim::Int,
    ode_p₀::AbstractVector,
    N::Int,
)
    Φ_dim = ketdim^2
    sens_state_dim = Φ_dim * (1 + n_params)

    # Initial condition: Φ₀ = I (complex), all sensitivities = 0
    sens_x₀ = zeros(ComplexF64, sens_state_dim)
    sens_x₀[1:Φ_dim] = vec(Matrix{ComplexF64}(I, ketdim, ketdim))

    sens_probs = [ODEProblem(make_f!_sens(), sens_x₀, (0.0, 1.0), ode_p₀) for _ = 1:(N-1)]
    sens_state = zeros(ComplexF64, sens_state_dim)

    return sens_probs, sens_state
end

"""
    extract_sensitivity_solution!(result, final_state, Φ_dim, n_params)

Extract propagator and sensitivity matrices from the augmented ODE solution into PropagatorResult.

# Arguments
- `result::PropagatorResult`: Target for Φ (Φ_vec) and ∂Φ/∂p (S_mat)
- `final_state::AbstractVector`: Final augmented ODE solution `[Φ; S₁; ...; Sₙ]`
- `Φ_dim::Int`: Dimension of vectorized propagator
- `n_params::Int`: Number of parameters
"""
function extract_sensitivity_solution!(
    result::PropagatorResult,
    final_state::AbstractVector,
    Φ_dim::Int,
    n_params::Int,
)
    # Extract Φ
    copyto!(result.Φ_vec, 1, final_state, 1, Φ_dim)

    # Extract sensitivities ∂Φ/∂p
    @inbounds for j = 1:n_params
        offset = Φ_dim + (j - 1) * Φ_dim
        for i = 1:Φ_dim
            result.S_mat[i, j] = final_state[offset+i]
        end
    end

    return nothing
end

"""
    _compute_ode_jacobian_analytical!(𝒮::SplineIntegrator, pₖ, k)

Compute Jacobian using the sensitivity-augmented ODE.

Solves a single augmented ODE containing `[Φ; S₁; ...; Sₙ]` where `Sⱼ = ∂Φ/∂pⱼ`,
then extracts the results into `𝒮.prop_results[k]`.
"""
@views function _compute_ode_jacobian_analytical!(
    𝒮::SplineIntegrator,
    pₖ::AbstractVector,
    k::Int,
)
    Φ_dim = propagator_side_dim(𝒮)^2
    n_params = length(pₖ)

    sens_prob = remake(𝒮.sens_probs[k], p = pₖ)
    sol = if 𝒮.alg isa Rodas5PAlg
        solve(
            sens_prob,
            Rodas5P(autodiff = AutoFiniteDiff());
            abstol = 𝒮.tol,
            reltol = 𝒮.tol,
            save_everystep = false,
        )
    else
        solve(sens_prob, Tsit5(); abstol = 𝒮.tol, reltol = 𝒮.tol, save_everystep = false)
    end

    extract_sensitivity_solution!(𝒮.prop_results[k], sol.u[end], Φ_dim, n_params)

    return nothing
end

"""
    _compute_ode_hessian_analytical!(𝒮::SplineIntegrator, pₖ, k)

Compute exact Hessian using the second-order sensitivity ODE.

Solves the augmented ODE `[Φ; S₁...Sₙ; T₁₁, T₁₂, ..., Tₙₙ]` and extracts:
- Propagator and first-order sensitivities into `𝒮.prop_results[k]`
- Second-order sensitivities T_{ij} as reshaped complex matrices

Returns a vector of `Matrix{ComplexF64}` for the active pairs, where
`T_mats[pair_k]` corresponds to `𝒮.hess_active_pairs[pair_k]`.
"""
@views function _compute_ode_hessian_analytical!(
    𝒮::SplineIntegrator,
    pₖ::AbstractVector,
    k::Int,
)
    Φ_dim = propagator_side_dim(𝒮)^2
    n_params = 𝒮.hess_n_params
    active_pairs = 𝒮.hess_active_pairs
    n_active = length(active_pairs)
    first_order_dim = Φ_dim * (1 + n_params)

    hess_prob = remake(𝒮.hess_probs[k], p = pₖ)
    sol = solve(
        hess_prob,
        Tsit5();
        abstol = 𝒮.tol,
        reltol = 𝒮.tol,
        save_everystep = false,
        maxiters = 1_000_000,
    )
    final_state = sol.u[end]

    # Extract Φ and first-order sensitivities into prop_results (reuse existing helper)
    extract_sensitivity_solution!(𝒮.prop_results[k], final_state, Φ_dim, n_params)

    # Extract second-order sensitivity matrices T_{ij}
    pdim = propagator_side_dim(𝒮)
    T_mats = Vector{Matrix{ComplexF64}}(undef, n_active)
    @inbounds for pair_k = 1:n_active
        T_offset = first_order_dim + (pair_k - 1) * Φ_dim
        T_mats[pair_k] =
            reshape(copy(final_state[(T_offset+1):(T_offset+Φ_dim)]), pdim, pdim)
    end

    return T_mats
end
