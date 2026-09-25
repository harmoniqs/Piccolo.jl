"""
    OptimizedNonlinearKnotPointConstraint{F} <: AbstractNonlinearConstraint
    Piccolo.Control.QuantumConstraints.SplineConstraints.NonlinearKnotPointConstraint{F} <: AbstractNonlinearConstraint

**Optimized** constraint implementation with pre-allocated full sparse matrices.

!!! note "Exported as OptimizedNonlinearKnotPointConstraint"
    Exported (from Piccolo, slice 3c #431; re-exported by Piccolissimo) as
    `OptimizedNonlinearKnotPointConstraint` to avoid shadowing
    `DirectTrajOpt.NonlinearKnotPointConstraint`. Use either:
    - `Piccolo.OptimizedNonlinearKnotPointConstraint` (exported, recommended)
    - `Piccolo.Control.QuantumConstraints.SplineConstraints.NonlinearKnotPointConstraint` (internal name)

# Optimization vs DirectTrajOpt Version

DirectTrajOpt's version stores compact per-knot-point blocks and assembles them on each call.
This Piccolissimo-originated version (Piccolo home since slice 3c, #431) pre-allocates the full sparse structure once and writes directly.

**Performance Benefits:**
- **`get_full_jacobian` / `get_full_hessian` are allocation-free** — they return the
  pre-allocated reference, with no `spzeros` and no `vcat`
- **10-100x less memory** (single sparse matrix vs N blocks + overhead)
- **Cache-friendly**: ForwardDiff still operates on small variable blocks

!!! warning "`jacobian!` is NOT zero-allocation"
    This docstring claimed zero allocation for `jacobian!` as well. Measurement
    contradicts it: `jacobian!` allocates a fresh ForwardDiff config plus a per-knot
    `vcat` on every call, and reads the trajectory through the `Any`-inferred
    `traj[t][name]` accessor — 133 kB at `N`=51 and 524 kB at `N`=201, i.e. **linear in
    the knot count** (`evaluate!`: 75.9 kB → 299 kB, ~98% of it from the accessor's
    boxing, not from the constraint algebra). Only `get_full_jacobian` /
    `get_full_hessian` are allocation-free. See `docs/adr/0010` and CONTEXT.md's
    "Optimized Constraints" correction; `HermiteSmoothAccelerationConstraint` shows the
    go-forward design (a `ConstraintStencilTable`, knot-flat rather than
    zero-allocation), and bringing this ForwardDiff-backed type to the same standard
    needs cached differentiation configs and is deliberately not attempted here.

# Fields
- `g::F`: Constraint function mapping (variables..., params) -> constraint values
- `var_names::Vector{Symbol}`: Names of trajectory variables the constraint depends on
- `equality::Bool`: If true, g(x) = 0; if false, g(x) ≤ 0
- `times::Vector{Int}`: Time indices where constraint is applied
- `params::Vector`: Parameters for each time index (e.g., time-varying targets)
- `g_dim::Int`: Dimension of constraint output at each knot point
- `var_dim::Int`: Combined dimension of all constrained variables
- `dim::Int`: Total constraint dimension (g_dim * length(times))
- `∂g_full::SparseMatrixCSC{Float64, Int}`: Pre-allocated full Jacobian (dim × Z_dim)
- `μ∂²g_full::SparseMatrixCSC{Float64, Int}`: Pre-allocated full Hessian (Z_dim × Z_dim)

# Constructor

The constructor signature is **identical** to DirectTrajOpt's version for drop-in compatibility:

```julia
OptimizedNonlinearKnotPointConstraint(
    g::Function,
    names::Union{Symbol, AbstractVector{Symbol}},
    traj::NamedTrajectory;
    equality::Bool=true,
    times::AbstractVector{Int}=1:traj.N,
    params::AbstractVector=fill(nothing, length(times)),
    jacobian_structure::Union{Nothing, SparseMatrixCSC}=nothing,
    hessian_structure::Union{Nothing, SparseMatrixCSC}=nothing,
)
```

# Arguments
- `g::Function`: Constraint function. For single variable: `g(x)`. For multiple: `g(x, u)` or `g(x_concat)`.
- `names`: Variable name(s) to constrain. Single: `:x`. Multiple: `[:x, :u]`.
- `traj::NamedTrajectory`: The trajectory on which the constraint is defined.

# Keyword Arguments
- `equality::Bool=true`: If true, enforces `g(x) = 0`. Otherwise `g(x) ≤ 0`.
- `times::AbstractVector{Int}=1:traj.N`: Time indices where constraint is enforced.
- `params::AbstractVector=fill(nothing, length(times))`: Per-knot-point parameters.
- `jacobian_structure::Union{Nothing, SparseMatrixCSC}=nothing`: Custom sparsity pattern (g_dim × var_dim per block).
- `hessian_structure::Union{Nothing, SparseMatrixCSC}=nothing`: Custom sparsity pattern (var_dim × var_dim per block).

# Usage Examples

```julia
using Piccolo

# Simple norm constraint
constraint = OptimizedNonlinearKnotPointConstraint(
    x -> [norm(x) - 1.0],
    :x, traj;
    equality=false
)

# Multi-variable constraint with separate arguments
constraint = OptimizedNonlinearKnotPointConstraint(
    (x, u) -> [x[1]^2 + u[1]^2 - 1.0],
    [:x, :u], traj;
    times=2:traj.N-1
)

# With custom sparsity (only first component matters)
∂g_structure = sparse([1], [1], [1.0], 1, 2)  # 1×2, only (1,1) non-zero
constraint = OptimizedNonlinearKnotPointConstraint(
    u -> [u[1]^2 - 1.0],
    :u, traj;
    jacobian_structure=∂g_structure
)
```

# Comparison with DirectTrajOpt

```julia
using DirectTrajOpt
using Piccolo

# DirectTrajOpt version (compact blocks, assembly on each call)
dto_constraint = NonlinearKnotPointConstraint(g, :u, traj)  # Uses DTO version

# The optimized version (pre-allocated full matrix, allocation-free `get_full_*`)
pic_constraint = Piccolo.Control.QuantumConstraints.SplineConstraints.NonlinearKnotPointConstraint(g, :u, traj)

# Both implement the same interface and work with DirectTrajOptProblem
prob = DirectTrajOptProblem(traj, objective, integrators; 
    constraints=[pic_constraint]  # Use optimized version
)
```
"""
struct NonlinearKnotPointConstraint{F} <: AbstractNonlinearConstraint
    g::F
    var_names::Vector{Symbol}
    equality::Bool
    times::Vector{Int}
    params::Vector
    g_dim::Int
    var_dim::Int
    dim::Int
    ∂g_full::SparseMatrixCSC{Float64,Int}
    μ∂²g_full::SparseMatrixCSC{Float64,Int}

    function NonlinearKnotPointConstraint(
        g::Function,
        names::AbstractVector{Symbol},
        traj::NamedTrajectory,
        params::AbstractVector;
        equality::Bool = true,
        times::AbstractVector{Int} = 1:traj.N,
        jacobian_structure::Union{Nothing,SparseMatrixCSC} = nothing,
        hessian_structure::Union{Nothing,SparseMatrixCSC} = nothing,
    )
        @assert length(params) == length(times) "params must have the same length as times"

        # Get component indices for all variables
        x_comps = vcat([traj.components[name] for name in names]...)
        var_dim = length(x_comps)

        # Determine constraint dimension by evaluating with first parameter
        Z⃗ = vec(traj)
        x_slice_test = slice(1, x_comps, traj.dim)
        test_result = g(Z⃗[x_slice_test], params[1])
        @assert test_result isa AbstractVector{<:Real} "Constraint function must return a vector"
        g_dim = length(test_result)

        dim = g_dim * length(times)
        Z_dim = traj.dim * traj.N + traj.global_dim

        # Build full Jacobian structure
        if !isnothing(jacobian_structure)
            @assert size(jacobian_structure) == (g_dim, var_dim) "jacobian_structure must be (g_dim=$g_dim × var_dim=$var_dim)"

            # Replicate block structure for all knot points
            ∂g_full = spzeros(dim, Z_dim)
            for (i, t) in enumerate(times)
                row_range = slice(i, g_dim)
                col_range = slice(t, x_comps, traj.dim)
                ∂g_full[row_range, col_range] = jacobian_structure
            end
        else
            # Default: assume all entries in each block may be non-zero
            ∂g_full = spzeros(dim, Z_dim)
            for (i, t) in enumerate(times)
                row_range = slice(i, g_dim)
                col_range = slice(t, x_comps, traj.dim)
                ∂g_full[row_range, col_range] .= 1.0  # Mark structure
            end
        end

        # Build full Hessian structure
        if !isnothing(hessian_structure)
            @assert size(hessian_structure) == (var_dim, var_dim) "hessian_structure must be (var_dim=$var_dim × var_dim=$var_dim)"

            # Block diagonal structure
            μ∂²g_full = spzeros(Z_dim, Z_dim)
            for (i, t) in enumerate(times)
                block_range = slice(t, x_comps, traj.dim)
                μ∂²g_full[block_range, block_range] = hessian_structure
            end
        else
            # Default: assume all entries in each block may be non-zero
            μ∂²g_full = spzeros(Z_dim, Z_dim)
            for (i, t) in enumerate(times)
                block_range = slice(t, x_comps, traj.dim)
                μ∂²g_full[block_range, block_range] .= 1.0  # Mark structure
            end
        end

        return new{typeof(g)}(
            g,
            names,
            equality,
            times,
            params,
            g_dim,
            var_dim,
            dim,
            ∂g_full,
            μ∂²g_full,
        )
    end
end

# Convenience constructor without params - creates wrapper that ignores param argument
function NonlinearKnotPointConstraint(
    g::Function,
    names::AbstractVector{Symbol},
    traj::NamedTrajectory;
    times::AbstractVector{Int} = 1:traj.N,
    kwargs...,
)
    num_vars = length(names)

    if num_vars == 1
        # Single variable: g(x) where x is the variable values
        params = [nothing for _ in times]
        g_param = (x, _) -> g(x)
        return NonlinearKnotPointConstraint(
            g_param,
            names,
            traj,
            params;
            times = times,
            kwargs...,
        )
    else
        # Multiple variables: determine if g expects separate arguments or concatenated

        # Get component ranges for splitting concatenated vector
        comp_ranges = Vector{UnitRange{Int}}(undef, num_vars)
        offset = 0
        for (i, name) in enumerate(names)
            comp_len = length(traj.components[name])
            comp_ranges[i] = (offset+1):(offset+comp_len)
            offset += comp_len
        end

        params = [nothing for _ in times]

        # Test if g accepts separate arguments
        Z⃗ = vec(traj)
        x_comps = vcat([traj.components[name] for name in names]...)
        x_slice = slice(1, x_comps, traj.dim)
        test_vec = Z⃗[x_slice]
        test_args = [test_vec[r] for r in comp_ranges]

        accepts_separate_args = false
        try
            result = g(test_args...)
            accepts_separate_args = result isa AbstractVector
        catch
            accepts_separate_args = false
        end

        if accepts_separate_args
            # Wrapper that splits concatenated vector into separate arguments
            g_param = function (x_concat, _)
                args = [x_concat[r] for r in comp_ranges]
                return g(args...)
            end
        else
            # g expects a single concatenated vector
            g_param = (x, _) -> g(x)
        end

        return NonlinearKnotPointConstraint(
            g_param,
            names,
            traj,
            params;
            times = times,
            kwargs...,
        )
    end
end

function NonlinearKnotPointConstraint(
    g::Function,
    name::Symbol,
    traj::NamedTrajectory;
    kwargs...,
)
    return NonlinearKnotPointConstraint(g, [name], traj; kwargs...)
end

# ----------------------------------------------------------------------------- #
# Core Evaluation Methods
# ----------------------------------------------------------------------------- #

"""
    evaluate!(values::AbstractVector, constraint::NonlinearKnotPointConstraint, traj::NamedTrajectory)

Evaluate the constraint, storing results in-place in `values`.
Part of the common interface with integrators.
"""
function CommonInterface.evaluate!(
    values::AbstractVector,
    constraint::NonlinearKnotPointConstraint,
    traj::NamedTrajectory,
)
    @views for (i, t) ∈ enumerate(constraint.times)
        x_vals = vcat([traj[t][name] for name in constraint.var_names]...)
        values[slice(i, constraint.g_dim)] = constraint.g(x_vals, constraint.params[i])
    end
    return nothing
end

# ----------------------------------------------------------------------------- #
# Optimized Jacobian Methods
# ----------------------------------------------------------------------------- #

"""
    jacobian_structure(constraint::NonlinearKnotPointConstraint, traj::NamedTrajectory)

Return the sparsity structure of the full Jacobian (already stored).
"""
function jacobian_structure(constraint::NonlinearKnotPointConstraint, traj::NamedTrajectory)
    return constraint.∂g_full
end

"""
    jacobian!(constraint::NonlinearKnotPointConstraint, traj::NamedTrajectory)

Compute Jacobian and store directly in pre-allocated `constraint.∂g_full`.

Writes directly into the sparse matrix's `nzval`, so no matrix is assembled — but this
call is **not** allocation-free: it builds a fresh ForwardDiff config and a per-knot
`vcat`, and reads the trajectory through the `Any`-inferred accessor, for a cost linear
in the knot count (133 kB at `N`=51, 524 kB at `N`=201). See the type's docstring.
"""
function DirectTrajOpt.Constraints.jacobian!(
    constraint::NonlinearKnotPointConstraint,
    traj::NamedTrajectory,
)
    # Zero out non-zero values
    fill!(constraint.∂g_full.nzval, 0.0)

    x_comps = vcat([traj.components[name] for name in constraint.var_names]...)

    @views for (i, t) ∈ enumerate(constraint.times)
        x_vals = vcat([traj[t][name] for name in constraint.var_names]...)
        row_range = slice(i, constraint.g_dim)
        col_range = slice(t, x_comps, traj.dim)

        # Compute Jacobian directly into pre-allocated sparse matrix block
        ForwardDiff.jacobian!(
            constraint.∂g_full[row_range, col_range],
            x -> constraint.g(x, constraint.params[i]),
            x_vals,
        )
    end
    return nothing
end

"""
    get_full_jacobian(constraint::NonlinearKnotPointConstraint, traj::NamedTrajectory)

Return the full Jacobian (already computed and stored).
**Allocation-free** - just returns a reference to the pre-allocated matrix. (This
method, and `get_full_hessian`, are the only allocation-free ones on this type.)
"""
function DirectTrajOpt.Constraints.get_full_jacobian(
    constraint::NonlinearKnotPointConstraint,
    traj::NamedTrajectory,
)
    return constraint.∂g_full
end

"""
    eval_jacobian(constraint::NonlinearKnotPointConstraint, traj::NamedTrajectory)

High-level method: compute and return the full Jacobian.
Calls `jacobian!` then returns the pre-allocated matrix.
"""
function eval_jacobian(constraint::NonlinearKnotPointConstraint, traj::NamedTrajectory)
    jacobian!(constraint, traj)
    return constraint.∂g_full
end

# ----------------------------------------------------------------------------- #
# Optimized Hessian Methods
# ----------------------------------------------------------------------------- #

"""
    hessian_structure(constraint::NonlinearKnotPointConstraint, traj::NamedTrajectory)

Return the sparsity structure of the full Hessian (already stored).
"""
function hessian_structure(constraint::NonlinearKnotPointConstraint, traj::NamedTrajectory)
    return constraint.μ∂²g_full
end

"""
    hessian_of_lagrangian!(constraint::NonlinearKnotPointConstraint, traj::NamedTrajectory, μ::AbstractVector)

Compute Hessian weighted by Lagrange multipliers, storing directly in `constraint.μ∂²g_full`.

Writes directly into the sparse matrix's `nzval`, so no matrix is assembled — but, like
`jacobian!`, this call is **not** allocation-free and its cost is linear in the knot
count. See the type's docstring.
"""
function DirectTrajOpt.Constraints.hessian_of_lagrangian!(
    constraint::NonlinearKnotPointConstraint,
    traj::NamedTrajectory,
    μ::AbstractVector,
)
    # Zero out non-zero values
    fill!(constraint.μ∂²g_full.nzval, 0.0)

    x_comps = vcat([traj.components[name] for name in constraint.var_names]...)

    @views for (i, t) ∈ enumerate(constraint.times)
        x_vals = vcat([traj[t][name] for name in constraint.var_names]...)
        μₖ = μ[slice(i, constraint.g_dim)]
        block_range = slice(t, x_comps, traj.dim)

        # Compute Hessian directly into pre-allocated sparse matrix block
        ForwardDiff.hessian!(
            constraint.μ∂²g_full[block_range, block_range],
            x -> μₖ' * constraint.g(x, constraint.params[i]),
            x_vals,
        )
    end
    return nothing
end

"""
    get_full_hessian(constraint::NonlinearKnotPointConstraint, traj::NamedTrajectory)

Return the full Hessian (already computed and stored).
**Allocation-free** - just returns a reference to the pre-allocated matrix.
"""
function DirectTrajOpt.Constraints.get_full_hessian(
    constraint::NonlinearKnotPointConstraint,
    traj::NamedTrajectory,
)
    return constraint.μ∂²g_full
end

"""
    eval_hessian_of_lagrangian(constraint::NonlinearKnotPointConstraint, traj::NamedTrajectory, μ::AbstractVector)

High-level method: compute and return the full Hessian of the Lagrangian.
Calls `hessian_of_lagrangian!` then returns the pre-allocated matrix.
"""
function eval_hessian_of_lagrangian(
    constraint::NonlinearKnotPointConstraint,
    traj::NamedTrajectory,
    μ::AbstractVector,
)
    hessian_of_lagrangian!(constraint, traj, μ)
    return constraint.μ∂²g_full
end

# ============================================================================= #
# Tests
# ============================================================================= #

@testitem "OptimizedNonlinearKnotPointConstraint is exported correctly" begin
    using Piccolo

    # Verify the exported name exists and is the optimized version (with full sparse matrices)
    @test isdefined(Piccolo, :OptimizedNonlinearKnotPointConstraint)

    # Verify it's the optimized version (has ∂g_full field) not DirectTrajOpt version (has ∂gs field)
    # We check this by looking at the fieldnames
    fields = fieldnames(Piccolo.OptimizedNonlinearKnotPointConstraint)
    @test :∂g_full ∈ fields  # the optimized version
    @test :μ∂²g_full ∈ fields  # the optimized version
    @test :∂gs ∉ fields  # NOT DirectTrajOpt's version
end

@testitem "NonlinearKnotPointConstraint - single variable" begin
    # NonlinearKnotPointConstraint is defined in this file, available in test scope
    # Import the interface functions from DirectTrajOpt where they're defined
    using DirectTrajOpt.CommonInterface: evaluate!
    import DirectTrajOpt.Constraints: get_full_jacobian, get_full_hessian
    import DirectTrajOpt.Constraints: jacobian!, hessian_of_lagrangian!
    using TrajectoryIndexingUtils
    using LinearAlgebra
    using NamedTrajectories
    using ForwardDiff

    # Create test trajectory
    N = 10
    traj = NamedTrajectory(
        (x = randn(2, N), u = randn(1, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
    )

    g(a) = [norm(a) - 1.0]
    times = 1:traj.N

    NLC = Piccolo.OptimizedNonlinearKnotPointConstraint(
        g,
        :u,
        traj;
        times = times,
        equality = false,
    )
    U_SLICE(k) = slice(k, traj.components[:u], traj.dim)

    ĝ(Z⃗) = vcat([g(Z⃗[U_SLICE(k)]) for k ∈ times]...)

    # Test evaluate!
    δ = zeros(NLC.dim)
    evaluate!(δ, NLC, traj)
    @test δ ≈ ĝ(vec(traj))

    # Test jacobian!
    jacobian!(NLC, traj)
    ∂g_full = get_full_jacobian(NLC, traj)
    ∂g_autodiff = ForwardDiff.jacobian(ĝ, vec(traj))

    @test ∂g_full[:, 1:(traj.dim*traj.N)] ≈ ∂g_autodiff

    # Test hessian_of_lagrangian
    μ = randn(length(times))
    hessian_of_lagrangian!(NLC, traj, μ)
    μ∂²g_full = get_full_hessian(NLC, traj)
    hessian_autodiff = ForwardDiff.hessian(Z -> μ'ĝ(Z), vec(traj))

    @test μ∂²g_full[1:(traj.dim*traj.N), 1:(traj.dim*traj.N)] ≈ hessian_autodiff
end

@testitem "NonlinearKnotPointConstraint - multiple variables" begin
    # NonlinearKnotPointConstraint is defined in this file, available in test scope
    # Import the interface functions from DirectTrajOpt where they're defined
    using DirectTrajOpt.CommonInterface: evaluate!
    import DirectTrajOpt.Constraints: get_full_jacobian, get_full_hessian
    import DirectTrajOpt.Constraints: jacobian!, hessian_of_lagrangian!
    using TrajectoryIndexingUtils
    using LinearAlgebra
    using NamedTrajectories
    using ForwardDiff

    N = 10
    traj = NamedTrajectory(
        (x = randn(2, N), u = randn(1, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
    )

    # Constraint with separate arguments
    g_separate(x, u) = [x[1]^2 + x[2]^2 - 1.0, u[1] - 0.5]
    times = 1:traj.N

    NLC = Piccolo.OptimizedNonlinearKnotPointConstraint(
        g_separate,
        [:x, :u],
        traj;
        times = times,
        equality = false,
    )

    X_SLICE(k) = slice(k, traj.components[:x], traj.dim)
    U_SLICE(k) = slice(k, traj.components[:u], traj.dim)
    ĝ(Z⃗) = vcat([g_separate(Z⃗[X_SLICE(k)], Z⃗[U_SLICE(k)]) for k ∈ times]...)

    # Test evaluate!
    δ = zeros(NLC.dim)
    evaluate!(δ, NLC, traj)
    @test δ ≈ ĝ(vec(traj))

    # Test jacobian!
    jacobian!(NLC, traj)
    ∂g_full = get_full_jacobian(NLC, traj)
    ∂g_autodiff = ForwardDiff.jacobian(ĝ, vec(traj))

    @test ∂g_full[:, 1:(traj.dim*traj.N)] ≈ ∂g_autodiff

    # Test hessian_of_lagrangian
    μ = randn(2 * traj.N)
    hessian_of_lagrangian!(NLC, traj, μ)
    μ∂²g_full = get_full_hessian(NLC, traj)
    hessian_autodiff = ForwardDiff.hessian(Z -> μ'ĝ(Z), vec(traj))

    @test μ∂²g_full[1:(traj.dim*traj.N), 1:(traj.dim*traj.N)] ≈ hessian_autodiff
end

@testitem "NonlinearKnotPointConstraint - custom sparsity" begin
    # NonlinearKnotPointConstraint is defined in this file, available in test scope
    # Import the interface functions from DirectTrajOpt where they're defined
    using DirectTrajOpt.CommonInterface: evaluate!
    import DirectTrajOpt.Constraints: get_full_jacobian
    import DirectTrajOpt.Constraints: jacobian!
    using TrajectoryIndexingUtils
    using LinearAlgebra
    using NamedTrajectories
    using ForwardDiff
    using SparseArrays

    N = 10
    traj = NamedTrajectory(
        (x = randn(2, N), u = randn(2, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
    )

    # Constraint that only depends on first component of u
    g(u) = [u[1]^2 - 1.0]

    g_dim = 1
    var_dim = 2

    # Custom sparsity: only first column non-zero
    ∂g_structure = spzeros(g_dim, var_dim)
    ∂g_structure[1, 1] = 1.0

    # Custom Hessian sparsity: only (1,1) entry
    μ∂²g_structure = spzeros(var_dim, var_dim)
    μ∂²g_structure[1, 1] = 1.0

    NLC = Piccolo.OptimizedNonlinearKnotPointConstraint(
        g,
        :u,
        traj;
        jacobian_structure = ∂g_structure,
        hessian_structure = μ∂²g_structure,
        equality = false,
    )

    # Verify sparsity structure was preserved
    @test nnz(NLC.∂g_full) == traj.N  # Only N non-zeros (one per knot point)
    @test nnz(NLC.μ∂²g_full) == traj.N  # Only N non-zeros (diagonal blocks)

    # Test computation still works
    U_SLICE(k) = slice(k, traj.components[:u], traj.dim)
    ĝ(Z⃗) = vcat([g(Z⃗[U_SLICE(k)]) for k ∈ 1:traj.N]...)

    δ = zeros(NLC.dim)
    evaluate!(δ, NLC, traj)
    @test δ ≈ ĝ(vec(traj))

    jacobian!(NLC, traj)
    ∂g_full = get_full_jacobian(NLC, traj)
    ∂g_autodiff = ForwardDiff.jacobian(ĝ, vec(traj))
    @test ∂g_full[:, 1:(traj.dim*traj.N)] ≈ ∂g_autodiff
end
