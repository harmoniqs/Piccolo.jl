"""
    NonlinearSegmentConstraint{F} <: AbstractNonlinearConstraint

Constraint applied to trajectory segments (between adjacent knot points k and k+1).

Similar to `NonlinearKnotPointConstraint` but evaluates constraints on intervals between
knot points rather than at individual points. Useful for enforcing properties like:
- Cubic spline extrema remain bounded
- Monotonicity between knots
- Maximum curvature constraints
- Any property requiring access to both k and k+1

# Architecture

Follows the integrator pattern for accessing adjacent knot points. The constraint function
receives variables from BOTH knot points k and k+1, just like `BilinearIntegrator` does.

**Jacobian structure spans two knot points:**
```julia
∂g[segment_slice, slice(k, 1:2*traj.dim, traj.dim)] = jacobian_of_g_wrt_[zₖ, zₖ₊₁]
```

This enables efficient parallel evaluation (segments don't overlap in Jacobian) and
reuses the proven integrator architecture.

# Fields
- `g::F`: Constraint function `g(varsₖ, varsₖ₊₁, params) -> vector`
- `var_names::Vector{Symbol}`: Variable names to extract from trajectory
- `equality::Bool`: If true, g = 0; if false, g ≤ 0
- `segments::Vector{Int}`: Segment indices (segment k connects knots k and k+1)
- `params::Vector`: Per-segment parameters (e.g., bounds, timestep)
- `g_dim::Int`: Constraint output dimension per segment
- `var_dim::Int`: Combined dimension of variables (from a single knot point)
- `dim::Int`: Total constraint dimension (g_dim * length(segments))
- `∂g_full::SparseMatrixCSC{Float64, Int}`: Pre-allocated Jacobian
- `μ∂²g_full::SparseMatrixCSC{Float64, Int}`: Pre-allocated Hessian

# Constructor

```julia
NonlinearSegmentConstraint(
    g::Function,
    names::Union{Symbol, AbstractVector{Symbol}},
    traj::NamedTrajectory;
    equality::Bool=true,
    segments::AbstractVector{Int}=1:traj.N-1,
    params::AbstractVector=fill(nothing, length(segments)),
    jacobian_structure::Union{Nothing, SparseMatrixCSC}=nothing,
    hessian_structure::Union{Nothing, SparseMatrixCSC}=nothing,
)
```

# Arguments
- `g::Function`: Constraint function with signature:
  - `g(varsₖ, varsₖ₊₁, param) -> Vector` for constraints with params
  - `g(varsₖ, varsₖ₊₁) -> Vector` for constraints without params
  where `varsₖ` and `varsₖ₊₁` are extracted variables from knot points k and k+1
- `names`: Variable name(s) to extract. Single `:u` or multiple `[:u, :du]`
- `traj`: Trajectory defining the optimization problem

# Keyword Arguments
- `equality::Bool=true`: Constraint type (equality vs inequality)
- `segments::AbstractVector{Int}=1:traj.N-1`: Segment indices to constrain
- `params::AbstractVector=fill(nothing, length(segments))`: Per-segment parameters
- `jacobian_structure`: Custom sparsity (g_dim × 2*var_dim per segment block)
- `hessian_structure`: Custom sparsity (2*var_dim × 2*var_dim per segment block)

# Examples

## Cubic Spline Boundedness
```julia
# Ensure spline stays in [-1, 1] between all knots
function spline_bound_constraint(varsₖ, varsₖ₊₁, params)
    u_min, u_max, Δt = params
    u_k, du_k = varsₖ[1:end÷2], varsₖ[end÷2+1:end]
    u_kp1, du_kp1 = varsₖ₊₁[1:end÷2], varsₖ₊₁[end÷2+1:end]
    
    # Find extrema in [0, 1] for cubic Hermite spline
    extrema_vals = find_cubic_extrema(u_k, u_kp1, du_k, du_kp1, Δt)
    
    # Return inequality constraints: u_min ≤ u(τ) ≤ u_max
    return [u_min .- extrema_vals; extrema_vals .- u_max]
end

constraint = NonlinearSegmentConstraint(
    spline_bound_constraint,
    [:u, :du],
    traj;
    equality=false,
    params=[(u_min, u_max, traj.Δt[k]) for k in 1:traj.N-1]
)
```

## Monotonicity Constraint
```julia
# Ensure control is monotonically increasing
function monotonic_constraint(varsₖ, varsₖ₊₁, param)
    u_k = varsₖ[1]
    u_kp1 = varsₖ₊₁[1]
    return [u_kp1 - u_k]  # u_{k+1} ≥ u_k
end

constraint = NonlinearSegmentConstraint(
    monotonic_constraint,
    :u,
    traj;
    equality=false
)
```

## Sufficient Condition for Bounds (Conservative)
```julia
# Limit slope magnitude based on distance to bounds
function sufficient_bound_constraint(varsₖ, varsₖ₊₁, params)
    u_min, u_max, C, Δt = params
    u_k, du_k = varsₖ[1:end÷2], varsₖ[end÷2+1:end]
    
    gap_to_max = u_max .- u_k
    gap_to_min = u_k .- u_min
    max_slope = C .* min.(gap_to_max, gap_to_min) ./ Δt
    
    # Inequality: max_slope - |du_k| ≥ 0
    return max_slope .- abs.(du_k)
end

constraint = NonlinearSegmentConstraint(
    sufficient_bound_constraint,
    [:u, :du],
    traj;
    equality=false,
    params=[(u_min, u_max, 2.0, traj.Δt[k]) for k in 1:traj.N-1]
)
```
"""
struct NonlinearSegmentConstraint{F} <: AbstractNonlinearConstraint
    g::F
    var_names::Vector{Symbol}
    equality::Bool
    segments::Vector{Int}
    params::Vector
    g_dim::Int
    var_dim::Int
    dim::Int
    ∂g_full::SparseMatrixCSC{Float64,Int}
    μ∂²g_full::SparseMatrixCSC{Float64,Int}

    function NonlinearSegmentConstraint(
        g::Function,
        names::AbstractVector{Symbol},
        traj::NamedTrajectory,
        params::AbstractVector;
        equality::Bool = true,
        segments::AbstractVector{Int} = 1:(traj.N-1),
        jacobian_structure::Union{Nothing,SparseMatrixCSC} = nothing,
        hessian_structure::Union{Nothing,SparseMatrixCSC} = nothing,
    )
        @assert length(params) == length(segments) "params must have same length as segments"
        @assert all(1 .<= segments .< traj.N) "segment indices must be in [1, N-1]"

        # Get component indices for all variables
        var_comps = vcat([traj.components[name] for name in names]...)
        var_dim = length(var_comps)

        # Determine constraint dimension by evaluating with first parameter
        Z⃗ = vec(traj)
        varsₖ = Z⃗[slice(segments[1], var_comps, traj.dim)]
        varsₖ₊₁ = Z⃗[slice(segments[1]+1, var_comps, traj.dim)]

        # Try calling with or without params
        test_result = if isnothing(params[1])
            g(varsₖ, varsₖ₊₁)
        else
            g(varsₖ, varsₖ₊₁, params[1])
        end

        @assert test_result isa AbstractVector{<:Real} "Constraint function must return a vector"
        g_dim = length(test_result)

        dim = g_dim * length(segments)
        Z_dim = traj.dim * traj.N + traj.global_dim

        # Build full Jacobian structure
        # Each segment block is g_dim × (2*var_dim) spanning [zₖ, zₖ₊₁]
        if !isnothing(jacobian_structure)
            @assert size(jacobian_structure) == (g_dim, 2*var_dim) "jacobian_structure must be g_dim × 2*var_dim"
            # Replicate custom structure for all segments
            ∂g_rows = Int[]
            ∂g_cols = Int[]
            for (seg_idx, k) in enumerate(segments)
                # Map segment constraint index to full constraint vector
                constraint_offset = (seg_idx - 1) * g_dim

                for (row, col, _) in zip(findnz(jacobian_structure)...)
                    # Determine if column is in zₖ or zₖ₊₁
                    if col <= var_dim
                        # Column in zₖ
                        comp_idx = var_comps[col]
                        Z_idx = (k - 1) * traj.dim + comp_idx
                    else
                        # Column in zₖ₊₁
                        comp_idx = var_comps[col-var_dim]
                        Z_idx = k * traj.dim + comp_idx
                    end
                    push!(∂g_rows, constraint_offset + row)
                    push!(∂g_cols, Z_idx)
                end
            end
            ∂g_full = sparse(∂g_rows, ∂g_cols, ones(length(∂g_rows)), dim, Z_dim)
        else
            # Default: dense blocks (g_dim × 2*var_dim per segment)
            ∂g_rows = Int[]
            ∂g_cols = Int[]
            for (seg_idx, k) in enumerate(segments)
                constraint_offset = (seg_idx - 1) * g_dim
                for row = 1:g_dim
                    # Columns for zₖ
                    for comp_idx in var_comps
                        Z_idx = (k - 1) * traj.dim + comp_idx
                        push!(∂g_rows, constraint_offset + row)
                        push!(∂g_cols, Z_idx)
                    end
                    # Columns for zₖ₊₁
                    for comp_idx in var_comps
                        Z_idx = k * traj.dim + comp_idx
                        push!(∂g_rows, constraint_offset + row)
                        push!(∂g_cols, Z_idx)
                    end
                end
            end
            ∂g_full = sparse(∂g_rows, ∂g_cols, ones(length(∂g_rows)), dim, Z_dim)
        end

        # Build full Hessian structure
        # Each segment block couples (2*var_dim × 2*var_dim) variables
        if !isnothing(hessian_structure)
            @assert size(hessian_structure) == (2*var_dim, 2*var_dim) "hessian_structure must be 2*var_dim × 2*var_dim"
            # Replicate custom structure for all segments
            μ∂²g_rows = Int[]
            μ∂²g_cols = Int[]
            for k in segments
                # Map variables to trajectory indices
                Z_indices_k = [(k - 1) * traj.dim + comp for comp in var_comps]
                Z_indices_kp1 = [k * traj.dim + comp for comp in var_comps]
                Z_indices = vcat(Z_indices_k, Z_indices_kp1)

                for (row, col, _) in zip(findnz(hessian_structure)...)
                    push!(μ∂²g_rows, Z_indices[row])
                    push!(μ∂²g_cols, Z_indices[col])
                end
            end
            μ∂²g_full = sparse(μ∂²g_rows, μ∂²g_cols, ones(length(μ∂²g_rows)), Z_dim, Z_dim)
        else
            # Default: dense blocks (2*var_dim × 2*var_dim per segment)
            μ∂²g_rows = Int[]
            μ∂²g_cols = Int[]
            for k in segments
                Z_indices_k = [(k - 1) * traj.dim + comp for comp in var_comps]
                Z_indices_kp1 = [k * traj.dim + comp for comp in var_comps]
                Z_indices = vcat(Z_indices_k, Z_indices_kp1)

                for row_idx in Z_indices
                    for col_idx in Z_indices
                        push!(μ∂²g_rows, row_idx)
                        push!(μ∂²g_cols, col_idx)
                    end
                end
            end
            μ∂²g_full = sparse(μ∂²g_rows, μ∂²g_cols, ones(length(μ∂²g_rows)), Z_dim, Z_dim)
        end

        return new{typeof(g)}(
            g,
            names,
            equality,
            collect(segments),
            params,
            g_dim,
            var_dim,
            dim,
            ∂g_full,
            μ∂²g_full,
        )
    end
end

# Single variable name convenience constructor
function NonlinearSegmentConstraint(
    g::Function,
    name::Symbol,
    traj::NamedTrajectory,
    params::AbstractVector;
    kwargs...,
)
    return NonlinearSegmentConstraint(g, [name], traj, params; kwargs...)
end

# No-params convenience constructor
function NonlinearSegmentConstraint(
    g::Function,
    names::Union{Symbol,AbstractVector{Symbol}},
    traj::NamedTrajectory;
    segments::AbstractVector{Int} = 1:(traj.N-1),
    kwargs...,
)
    params = fill(nothing, length(segments))
    return NonlinearSegmentConstraint(
        g,
        names,
        traj,
        params;
        segments = segments,
        kwargs...,
    )
end

# ----------------------------------------------------------------------------- #
# Interface Implementation
# ----------------------------------------------------------------------------- #

"""
    evaluate!(values::AbstractVector, constraint::NonlinearSegmentConstraint, traj::NamedTrajectory)

Evaluate constraint on all segments in parallel.
"""
function evaluate!(
    values::AbstractVector,
    constraint::NonlinearSegmentConstraint,
    traj::NamedTrajectory,
)
    Z⃗ = vec(traj)
    var_comps = vcat([traj.components[name] for name in constraint.var_names]...)

    Threads.@threads for (seg_idx, k) in collect(enumerate(constraint.segments))
        # Extract variables from both knot points
        varsₖ = Z⃗[slice(k, var_comps, traj.dim)]
        varsₖ₊₁ = Z⃗[slice(k+1, var_comps, traj.dim)]

        # Evaluate constraint
        result = if isnothing(constraint.params[seg_idx])
            constraint.g(varsₖ, varsₖ₊₁)
        else
            constraint.g(varsₖ, varsₖ₊₁, constraint.params[seg_idx])
        end

        # Write to output
        offset = (seg_idx - 1) * constraint.g_dim
        values[(offset+1):(offset+constraint.g_dim)] .= result
    end

    return nothing
end

"""
    jacobian_structure(constraint::NonlinearSegmentConstraint)

Return pre-allocated Jacobian structure.
"""
function jacobian_structure(constraint::NonlinearSegmentConstraint)
    return constraint.∂g_full
end

"""
    jacobian!(constraint::NonlinearSegmentConstraint, traj::NamedTrajectory)

Compute Jacobian using ForwardDiff on each segment block.
"""
function jacobian!(constraint::NonlinearSegmentConstraint, traj::NamedTrajectory)
    Z⃗ = vec(traj)
    var_comps = vcat([traj.components[name] for name in constraint.var_names]...)
    ∂g = constraint.∂g_full

    # Zero out (preserve structure)
    ∂g.nzval .= 0.0

    Threads.@threads for (seg_idx, k) in collect(enumerate(constraint.segments))
        # Extract indices for this segment's block in full Jacobian
        constraint_slice = (seg_idx - 1) * constraint.g_dim .+ (1:constraint.g_dim)
        Z_slice_k = slice(k, var_comps, traj.dim)
        Z_slice_kp1 = slice(k+1, var_comps, traj.dim)
        Z_slice_segment = vcat(Z_slice_k, Z_slice_kp1)

        # Create local Jacobian block (g_dim × 2*var_dim)
        ∂g_local = zeros(constraint.g_dim, 2*constraint.var_dim)

        # Compute Jacobian via ForwardDiff
        ForwardDiff.jacobian!(
            ∂g_local,
            vars_concat -> begin
                varsₖ = vars_concat[1:constraint.var_dim]
                varsₖ₊₁ = vars_concat[(constraint.var_dim+1):end]
                if isnothing(constraint.params[seg_idx])
                    return constraint.g(varsₖ, varsₖ₊₁)
                else
                    return constraint.g(varsₖ, varsₖ₊₁, constraint.params[seg_idx])
                end
            end,
            Z⃗[Z_slice_segment],
        )

        # Map local block to full sparse matrix
        # This is safe for parallel writes because segments don't overlap
        for (local_row, constraint_row) in enumerate(constraint_slice)
            for local_col = 1:(2*constraint.var_dim)
                if local_col <= constraint.var_dim
                    Z_col = Z_slice_k[local_col]
                else
                    Z_col = Z_slice_kp1[local_col-constraint.var_dim]
                end

                # Find the position in sparse matrix
                for idx in nzrange(∂g, Z_col)
                    if ∂g.rowval[idx] == constraint_row
                        ∂g.nzval[idx] = ∂g_local[local_row, local_col]
                        break
                    end
                end
            end
        end
    end

    return nothing
end

"""
    hessian_structure(constraint::NonlinearSegmentConstraint)

Return pre-allocated Hessian structure.
"""
function hessian_structure(constraint::NonlinearSegmentConstraint)
    return constraint.μ∂²g_full
end

"""
    hessian_of_lagrangian(constraint::NonlinearSegmentConstraint, traj::NamedTrajectory, μ::AbstractVector)

Compute Hessian of Lagrangian using ForwardDiff on each segment block.
"""
function hessian_of_lagrangian(
    constraint::NonlinearSegmentConstraint,
    traj::NamedTrajectory,
    μ::AbstractVector,
)
    Z⃗ = vec(traj)
    var_comps = vcat([traj.components[name] for name in constraint.var_names]...)
    μ∂²g = constraint.μ∂²g_full

    # Zero out (preserve structure)
    μ∂²g.nzval .= 0.0

    for (seg_idx, k) in enumerate(constraint.segments)
        # Extract multipliers for this segment
        μ_seg = μ[((seg_idx-1)*constraint.g_dim+1):(seg_idx*constraint.g_dim)]

        # Skip if multipliers are zero
        if norm(μ_seg) < 1e-16
            continue
        end

        # Extract indices
        Z_slice_k = slice(k, var_comps, traj.dim)
        Z_slice_kp1 = slice(k+1, var_comps, traj.dim)
        Z_slice_segment = vcat(Z_slice_k, Z_slice_kp1)

        # Compute weighted Hessian via ForwardDiff
        μ∂²g_local = ForwardDiff.hessian(
            vars_concat -> begin
                varsₖ = vars_concat[1:constraint.var_dim]
                varsₖ₊₁ = vars_concat[(constraint.var_dim+1):end]
                g_result = if isnothing(constraint.params[seg_idx])
                    constraint.g(varsₖ, varsₖ₊₁)
                else
                    constraint.g(varsₖ, varsₖ₊₁, constraint.params[seg_idx])
                end
                return dot(μ_seg, g_result)
            end,
            Z⃗[Z_slice_segment],
        )

        # Map local block to full sparse matrix
        for local_row = 1:(2*constraint.var_dim)
            for local_col = 1:(2*constraint.var_dim)
                Z_row = Z_slice_segment[local_row]
                Z_col = Z_slice_segment[local_col]

                # Find the position in sparse matrix
                for idx in nzrange(μ∂²g, Z_col)
                    if μ∂²g.rowval[idx] == Z_row
                        μ∂²g.nzval[idx] += μ∂²g_local[local_row, local_col]
                        break
                    end
                end
            end
        end
    end

    return nothing
end

# Convenience methods for CommonInterface
function get_full_jacobian(constraint::NonlinearSegmentConstraint, traj::NamedTrajectory)
    return constraint.∂g_full
end

function get_full_hessian(constraint::NonlinearSegmentConstraint, traj::NamedTrajectory)
    return constraint.μ∂²g_full
end

function eval_jacobian(constraint::NonlinearSegmentConstraint, traj::NamedTrajectory)
    jacobian!(constraint, traj)
    return constraint.∂g_full
end

function eval_hessian_of_lagrangian(
    constraint::NonlinearSegmentConstraint,
    traj::NamedTrajectory,
    μ::AbstractVector,
)
    hessian_of_lagrangian(constraint, traj, μ)
    return constraint.μ∂²g_full
end

# ----------------------------------------------------------------------------- #
# Tests
# ----------------------------------------------------------------------------- #

@testitem "NonlinearSegmentConstraint - basic functionality" begin
    using NamedTrajectories
    using SparseArrays
    using LinearAlgebra
    using Piccolo: NonlinearSegmentConstraint
    using DirectTrajOpt.CommonInterface: evaluate!, jacobian_structure
    using DirectTrajOpt.Constraints: jacobian!
    using DirectTrajOpt.CommonInterface: hessian_structure, hessian_of_lagrangian
    using DirectTrajOpt.Constraints: get_full_jacobian, get_full_hessian

    # Create simple trajectory with u and du
    N = 10
    u_dim = 2
    traj = NamedTrajectory(
        (u = randn(u_dim, N), du = randn(u_dim, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
    )

    # Simple constraint: u_{k+1} - u_k ≥ 0 (monotonicity)
    function monotonic(varsₖ, varsₖ₊₁)
        u_k = varsₖ[1:2]
        u_kp1 = varsₖ₊₁[1:2]
        return u_kp1 - u_k  # Should be ≥ 0
    end

    constraint = NonlinearSegmentConstraint(monotonic, [:u, :du], traj; equality = false)

    # Test dimensions
    @test constraint.g_dim == u_dim
    @test constraint.var_dim == 2 * u_dim  # u and du
    @test constraint.dim == u_dim * (N - 1)
    @test length(constraint.segments) == N - 1

    # Test evaluate!
    values = zeros(constraint.dim)
    evaluate!(values, constraint, traj)
    @test length(values) == constraint.dim

    # Test Jacobian structure
    ∂g = jacobian_structure(constraint)
    @test size(∂g) == (constraint.dim, traj.dim * traj.N + traj.global_dim)
    @test nnz(∂g) > 0

    # Test Jacobian computation
    jacobian!(constraint, traj)
    ∂g_eval = get_full_jacobian(constraint, traj)
    @test ∂g_eval === constraint.∂g_full

    # Test Hessian structure
    μ∂²g = hessian_structure(constraint)
    @test size(μ∂²g) ==
          (traj.dim * traj.N + traj.global_dim, traj.dim * traj.N + traj.global_dim)

    # Test Hessian computation
    μ = randn(constraint.dim)
    hessian_of_lagrangian(constraint, traj, μ)
    μ∂²g_eval = get_full_hessian(constraint, traj)
    @test μ∂²g_eval === constraint.μ∂²g_full
end

@testitem "NonlinearSegmentConstraint - with parameters" begin
    using NamedTrajectories
    using LinearAlgebra
    using Piccolo: NonlinearSegmentConstraint
    using DirectTrajOpt.CommonInterface: evaluate!

    N = 5
    traj = NamedTrajectory(
        (u = randn(1, N), du = randn(1, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
    )

    # Constraint with time-varying target
    function target_constraint(varsₖ, varsₖ₊₁, target)
        u_kp1 = varsₖ₊₁[1]
        return [abs(u_kp1) - target]  # |u| ≤ target
    end

    # N-1 segments need N-1 targets
    targets = collect(range(1.0, 3.0, length = N-1))  # Different target per segment

    constraint = NonlinearSegmentConstraint(
        target_constraint,
        [:u, :du],
        traj,
        targets;
        equality = false,
    )

    @test length(constraint.params) == N - 1
    @test constraint.params == targets

    # Evaluate
    values = zeros(constraint.dim)
    evaluate!(values, constraint, traj)
    @test length(values) == N - 1
end

@testitem "NonlinearSegmentConstraint - subset of segments" begin
    using NamedTrajectories
    using Piccolo: NonlinearSegmentConstraint
    using DirectTrajOpt.CommonInterface: jacobian_structure

    N = 20
    traj = NamedTrajectory(
        (u = randn(1, N), du = randn(1, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
    )

    # Only constrain middle segments
    segments = 5:15

    function simple_constraint(varsₖ, varsₖ₊₁)
        return [varsₖ₊₁[1] - varsₖ[1]]
    end

    constraint =
        NonlinearSegmentConstraint(simple_constraint, [:u, :du], traj; segments = segments)

    @test length(constraint.segments) == length(segments)
    @test constraint.dim == length(segments)

    # Verify only specified segments are in Jacobian structure
    ∂g = jacobian_structure(constraint)
    @test size(∂g, 1) == length(segments)
end

@testitem "NonlinearSegmentConstraint - custom sparsity, single-name ctor, params Hessian, eval accessors" begin
    using NamedTrajectories
    using DirectTrajOpt
    using DirectTrajOpt.CommonInterface: evaluate!
    using Piccolo: NonlinearSegmentConstraint
    using TrajectoryIndexingUtils
    using LinearAlgebra
    using SparseArrays
    using ForwardDiff

    N = 5
    traj = NamedTrajectory(
        (u = randn(1, N), du = 0.1 * randn(1, N), Δt = fill(0.5, N));
        timestep = :Δt,
        controls = :u,
    )

    # Single-variable-name convenience constructor
    g_mono = (uₖ, uₖ₊₁) -> [uₖ₊₁[1] - uₖ[1]]
    c_mono = NonlinearSegmentConstraint(
        g_mono,
        :u,
        traj,
        [nothing for _ = 1:(N-1)];
        equality = false,
    )
    δ = zeros(c_mono.dim)
    evaluate!(δ, c_mono, traj)
    @test all(δ .≈ [traj.u[1, k+1] - traj.u[1, k] for k = 1:(N-1)])

    # Custom sparsity: g(vₖ, vₖ₊₁, p) = p·uₖ² + duₖ₊₁ has Jacobian
    # [2p·uₖ, 0, 0, 1] in local (uₖ, duₖ, uₖ₊₁, duₖ₊₁) order — declare it.
    g_sparse = (vₖ, vₖ₊₁, p) -> [p[1] * vₖ[1]^2 + vₖ₊₁[2]]
    js = sparse([1, 1], [1, 4], [1.0, 1.0], 1, 4)
    hs = sparse([1], [1], [1.0], 4, 4)
    params = [0.7 for _ = 1:(N-1)]
    c_sparse = NonlinearSegmentConstraint(
        g_sparse,
        [:u, :du],
        traj,
        params;
        equality = false,
        jacobian_structure = js,
        hessian_structure = hs,
    )

    # Declared patterns are replicated per segment
    @test nnz(c_sparse.∂g_full) == 2 * (N - 1)
    @test nnz(c_sparse.μ∂²g_full) == (N - 1)

    Z⃗ = vec(traj)
    U_SLICE(k) = slice(k, traj.components[:u], traj.dim)
    DU_SLICE(k) = slice(k, traj.components[:du], traj.dim)
    ĝ(Z) = vcat(
        [
            g_sparse(
                vcat(Z[U_SLICE(k)], Z[DU_SLICE(k)]),
                vcat(Z[U_SLICE(k+1)], Z[DU_SLICE(k+1)]),
                params[k],
            ) for k = 1:(N-1)
        ]...,
    )

    # eval_jacobian: the assembled Jacobian matches forward-mode AD on the
    # declared pattern (AD has exactly the declared support)
    ∂g_eval = DirectTrajOpt.CommonInterface.eval_jacobian(c_sparse, traj)
    ∂g_fd = ForwardDiff.jacobian(ĝ, Z⃗)
    @test ∂g_eval ≈ ∂g_fd

    # eval_hessian_of_lagrangian with a ZERO multiplier block on segment 1:
    # the zero-μ skip must leave segment 1's curvature out of the total
    μ = vcat(zeros(1), ones(N - 2))
    μ∂²g_eval = DirectTrajOpt.CommonInterface.eval_hessian_of_lagrangian(c_sparse, traj, μ)
    H_fd = ForwardDiff.hessian(Z -> dot(μ, ĝ(Z)), Z⃗)
    @test μ∂²g_eval[1:length(Z⃗), 1:length(Z⃗)] ≈ H_fd
    # the params-weighted uₖ² curvature is 2·p·μ per segment; segment 1's block
    # must be exactly zero (skipped), every other block exactly 2p
    z_dim = traj.dim
    @test iszero(μ∂²g_eval[1, 1])
    for k = 2:(N-1)
        @test μ∂²g_eval[(k-1)*z_dim+1, (k-1)*z_dim+1] ≈ 2 * params[k]
    end
end
