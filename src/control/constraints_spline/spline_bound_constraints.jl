# ----------------------------------------------------------------------------- #
# Approach 1: Analytical Extrema Constraints (Exact)
# ----------------------------------------------------------------------------- #

"""
    CubicSplineExtremaConstraint(
        traj::NamedTrajectory,
        u_bounds::Tuple{Real, Real};
        control_name::Symbol=:u,
        segments::AbstractVector{Int}=1:traj.N-1
    )

Create segment constraints ensuring cubic spline remains bounded via analytical extrema.

This is the **exact** approach: finds critical points where du/dτ = 0 and constrains
the spline value at those points.

# Arguments
- `traj`: Trajectory with cubic spline pulse (must have :u and :du components)
- `u_bounds`: (u_min, u_max) bounds for the control

# Keyword Arguments
- `control_name`: Name of control variable (default :u)
- `segments`: Which segments to constrain (default: all)

# Returns
`NonlinearSegmentConstraint` that can be added to a `DirectTrajOptProblem`.

# Example
```julia
using Piccolo

# Create spline pulse problem
qcp = SplinePulseProblem(qtraj, N; Q=100.0, R=1e-2)
traj = get_trajectory(qcp)

# Add bound constraints
u_min, u_max = -10.0, 10.0
constraint = CubicSplineExtremaConstraint(traj, (u_min, u_max))

# Recreate problem with constraint
prob = DirectTrajOptProblem(
    traj, 
    qcp.objective, 
    qcp.integrators;
    constraints=vcat(qcp.constraints, [constraint])
)
```

# Mathematical Details

For each segment k → k+1, the constraint finds all critical points τ* ∈ (0,1) where
du/dτ = 0 by solving a quadratic equation. It then evaluates u(τ*) at each critical
point and imposes:

```
u_min ≤ u(τ*) ≤ u_max  for all critical points
```

This ensures the spline remains bounded throughout each segment.

# Performance
- Per-segment evaluation: ~10-50 μs (quadratic solve + evaluations)
- Differentiable via ForwardDiff (implicit function theorem applies)
- Typically 0-2 critical points per segment
"""
function CubicSplineExtremaConstraint(
    traj::NamedTrajectory,
    u_bounds::Tuple{Real,Real};
    control_name::Symbol = :u,
    segments::AbstractVector{Int} = 1:(traj.N-1),
)
    u_min, u_max = u_bounds
    du_name = Symbol(:d, control_name)

    @assert control_name ∈ traj.names "Control $control_name not found in trajectory"
    @assert du_name ∈ traj.names "Derivative $du_name not found in trajectory"
    @assert traj.timestep ∈ traj.names "Timestep not found in trajectory"

    n_drives = traj.dims[control_name]

    # Constraint function: check extrema for each drive
    # Each drive gets 8 constraints: 2 endpoints + up to 2 critical points, each with 2 bounds
    # This gives a fixed dimension: 8 * n_drives per segment
    function extrema_constraint(varsₖ, varsₖ₊₁, params)
        u_min, u_max, Δt = params

        # Split variables
        n = length(varsₖ) ÷ 2
        u_k = varsₖ[1:n]
        du_k = varsₖ[(n+1):end]
        u_kp1 = varsₖ₊₁[1:n]
        du_kp1 = varsₖ₊₁[(n+1):end]

        # Generic type for ForwardDiff compatibility
        T = eltype(varsₖ)
        constraints = T[]

        for i = 1:n
            # Always check exactly 4 points: 2 endpoints + 2 potential critical points
            # This ensures fixed constraint dimension

            # Endpoints (always exist)
            u_vals = [u_k[i], u_kp1[i]]

            # Find critical points
            τ_crits = find_cubic_critical_points(u_k[i], u_kp1[i], du_k[i], du_kp1[i], Δt)

            # Evaluate at critical points (pad with endpoints if fewer than 2 critical points)
            while length(τ_crits) < 2
                push!(τ_crits, 0.5)  # Use midpoint as dummy (will be feasible if endpoints are)
            end

            for τ in τ_crits[1:2]  # Take at most 2 critical points
                u_τ = evaluate_hermite_spline(τ, u_k[i], u_kp1[i], du_k[i], du_kp1[i], Δt)
                push!(u_vals, u_τ)
            end

            # Create 8 constraints per drive: 4 points × 2 bounds each
            # Inequality constraints: u_min ≤ u(τ) ≤ u_max
            # Reformulated as: u_min - u(τ) ≤ 0 and u(τ) - u_max ≤ 0
            for u_val in u_vals
                push!(constraints, u_min - u_val)  # u_min ≤ u
                push!(constraints, u_val - u_max)  # u ≤ u_max
            end
        end

        return constraints
    end

    # Create parameters for each segment
    params = [(u_min, u_max, traj[k].timestep) for k in segments]

    return NonlinearSegmentConstraint(
        extrema_constraint,
        [control_name, du_name],
        traj,
        params;
        equality = false,  # Inequality constraints
        segments = segments,
    )
end

# ----------------------------------------------------------------------------- #
# Approach 2: Sufficient Conditions (Conservative)
# ----------------------------------------------------------------------------- #

"""
    CubicSplineSufficientBoundConstraint(
        traj::NamedTrajectory,
        u_bounds::Tuple{Real, Real};
        control_name::Symbol=:u,
        safety_factor::Real=2.0,
        segments::AbstractVector{Int}=1:traj.N-1
    )

Create segment constraints ensuring cubic spline remains bounded via sufficient conditions.

This is the **conservative** approach: limits the slope magnitude based on distance
to bounds. Guarantees boundedness but may be overly restrictive.

# Arguments
- `traj`: Trajectory with cubic spline pulse
- `u_bounds`: (u_min, u_max) bounds for the control

# Keyword Arguments
- `control_name`: Name of control variable (default :u)
- `safety_factor`: Multiplier for slope limit (default 2.0)
- `segments`: Which segments to constrain (default: all)

# Returns
`NonlinearSegmentConstraint` that can be added to problem.

# Example
```julia
# More conservative than extrema constraint, but simpler
constraint = CubicSplineSufficientBoundConstraint(
    traj, (-10.0, 10.0);
    safety_factor=1.5  # Lower = more conservative
)
```

# Mathematical Details

For each knot k, the constraint enforces:

```
|du_k| ≤ C * min(u_max - u_k, u_k - u_min) / Δt
```

where C is the safety factor. This ensures the tangent at each knot point cannot
"escape" the bounds within one time step.

**Pros:**
- Simple to evaluate (no critical point search)
- Always convex
- Faster than analytical extrema

**Cons:**
- Conservative (excludes some feasible solutions)
- May require tuning safety factor

# Safety Factor Guidelines
- `C = 1.0`: Very conservative, nearly monotonic
- `C = 2.0`: Moderate (recommended default)
- `C = 3.0`: Less conservative, allows more curvature
"""
function CubicSplineSufficientBoundConstraint(
    traj::NamedTrajectory,
    u_bounds::Tuple{Real,Real};
    control_name::Symbol = :u,
    safety_factor::Real = 2.0,
    segments::AbstractVector{Int} = 1:(traj.N-1),
)
    u_min, u_max = u_bounds
    du_name = Symbol(:d, control_name)

    @assert control_name ∈ traj.names "Control $control_name not found in trajectory"
    @assert du_name ∈ traj.names "Derivative $du_name not found in trajectory"
    @assert safety_factor > 0 "safety_factor must be positive"

    n_drives = traj.dims[control_name]

    # Constraint function: limit slope based on gap to bounds
    function sufficient_constraint(varsₖ, varsₖ₊₁, params)
        u_min, u_max, C, Δt = params

        # Split variables (only need varsₖ for this constraint)
        n = length(varsₖ) ÷ 2
        u_k = varsₖ[1:n]
        du_k = varsₖ[(n+1):end]

        # Generic type for ForwardDiff compatibility
        T = eltype(varsₖ)
        constraints = T[]

        for i = 1:n
            gap_to_max = u_max - u_k[i]
            gap_to_min = u_k[i] - u_min
            min_gap = min(gap_to_max, gap_to_min)

            # Maximum allowed slope magnitude
            max_slope = C * min_gap / Δt

            # Inequality: |du_k| − max_slope ≤ 0 (backend convention: feasible ⟺ residual ≤ 0)
            push!(constraints, abs(du_k[i]) - max_slope)
        end

        return constraints
    end

    # Create parameters for each segment
    params = [(u_min, u_max, safety_factor, traj[k].timestep) for k in segments]

    return NonlinearSegmentConstraint(
        sufficient_constraint,
        [control_name, du_name],
        traj,
        params;
        equality = false,  # Inequality constraints
        segments = segments,
    )
end

# ----------------------------------------------------------------------------- #
# Multi-drive convenience constructors
# ----------------------------------------------------------------------------- #

"""
    CubicSplineExtremaConstraint(
        traj::NamedTrajectory,
        u_bounds::Vector{Tuple{Real, Real}};
        control_name::Symbol=:u,
        segments::AbstractVector{Int}=1:traj.N-1
    )

Create extrema constraints with per-drive bounds.

# Example
```julia
# Different bounds for different drives
bounds = [(-5.0, 5.0), (-10.0, 10.0), (-8.0, 8.0)]
constraint = CubicSplineExtremaConstraint(traj, bounds)
```
"""
function CubicSplineExtremaConstraint(
    traj::NamedTrajectory,
    u_bounds::Vector{<:Tuple{Real,Real}};
    control_name::Symbol = :u,
    segments::AbstractVector{Int} = 1:(traj.N-1),
)
    du_name = Symbol(:d, control_name)
    n_drives = traj.dims[control_name]

    @assert length(u_bounds) == n_drives "Must provide bounds for each drive"

    # Fixed dimension: 8 constraints per drive (2 endpoints + 2 critical points, each with 2 bounds)
    function extrema_constraint_multi(varsₖ, varsₖ₊₁, params)
        u_bounds_vec, Δt = params

        n = length(varsₖ) ÷ 2
        u_k = varsₖ[1:n]
        du_k = varsₖ[(n+1):end]
        u_kp1 = varsₖ₊₁[1:n]
        du_kp1 = varsₖ₊₁[(n+1):end]

        # Generic type for ForwardDiff compatibility
        T = eltype(varsₖ)
        constraints = T[]

        for i = 1:n
            u_min, u_max = u_bounds_vec[i]

            # Always check exactly 4 points: 2 endpoints + 2 potential critical points
            u_vals = [u_k[i], u_kp1[i]]

            # Find critical points
            τ_crits = find_cubic_critical_points(u_k[i], u_kp1[i], du_k[i], du_kp1[i], Δt)

            # Pad with midpoint if fewer than 2 critical points
            while length(τ_crits) < 2
                push!(τ_crits, 0.5)
            end

            for τ in τ_crits[1:2]
                u_τ = evaluate_hermite_spline(τ, u_k[i], u_kp1[i], du_k[i], du_kp1[i], Δt)
                push!(u_vals, u_τ)
            end

            # 8 constraints per drive
            for u_val in u_vals
                push!(constraints, u_min - u_val)
                push!(constraints, u_val - u_max)
            end
        end

        return constraints
    end

    params = [(u_bounds, traj[k].timestep) for k in segments]

    return NonlinearSegmentConstraint(
        extrema_constraint_multi,
        [control_name, du_name],
        traj,
        params;
        equality = false,
        segments = segments,
    )
end

"""
    CubicSplineSufficientBoundConstraint(
        traj::NamedTrajectory,
        u_bounds::Vector{Tuple{Real, Real}};
        control_name::Symbol=:u,
        safety_factor::Real=2.0,
        segments::AbstractVector{Int}=1:traj.N-1
    )

Create sufficient bound constraints with per-drive bounds.
"""
function CubicSplineSufficientBoundConstraint(
    traj::NamedTrajectory,
    u_bounds::Vector{<:Tuple{Real,Real}};
    control_name::Symbol = :u,
    safety_factor::Real = 2.0,
    segments::AbstractVector{Int} = 1:(traj.N-1),
)
    du_name = Symbol(:d, control_name)
    n_drives = traj.dims[control_name]

    @assert length(u_bounds) == n_drives "Must provide bounds for each drive"

    function sufficient_constraint_multi(varsₖ, varsₖ₊₁, params)
        u_bounds_vec, C, Δt = params

        n = length(varsₖ) ÷ 2
        u_k = varsₖ[1:n]
        du_k = varsₖ[(n+1):end]

        # Generic type for ForwardDiff compatibility
        T = eltype(varsₖ)
        constraints = T[]

        for i = 1:n
            u_min, u_max = u_bounds_vec[i]
            gap_to_max = u_max - u_k[i]
            gap_to_min = u_k[i] - u_min
            min_gap = min(gap_to_max, gap_to_min)
            max_slope = C * min_gap / Δt

            # |du_k| − max_slope ≤ 0 (backend convention: feasible ⟺ residual ≤ 0)
            push!(constraints, abs(du_k[i]) - max_slope)
        end

        return constraints
    end

    params = [(u_bounds, safety_factor, traj[k].timestep) for k in segments]

    return NonlinearSegmentConstraint(
        sufficient_constraint_multi,
        [control_name, du_name],
        traj,
        params;
        equality = false,
        segments = segments,
    )
end

# ----------------------------------------------------------------------------- #
# Approach 3: Slope/Derivative Extrema Constraints
# ----------------------------------------------------------------------------- #

"""
    CubicSplineSlopeConstraint(
        traj::NamedTrajectory,
        du_bounds::Tuple{Real, Real};
        control_name::Symbol=:u,
        segments::AbstractVector{Int}=1:traj.N-1
    )

Create segment constraints ensuring cubic spline derivative (slope) remains bounded.

This finds critical points where d²u/dt² = 0 and constrains the slope at those points
and at segment endpoints.

# Arguments
- `traj`: Trajectory with cubic spline pulse (must have :u and :du components)
- `du_bounds`: (du_min, du_max) bounds for the derivative

# Keyword Arguments
- `control_name`: Name of control variable (default :u)
- `segments`: Which segments to constrain (default: all)

# Returns
`NonlinearSegmentConstraint` that can be added to a `DirectTrajOptProblem`.

# Example
```julia
# Constrain slope to ±5
du_bounds = (-5.0, 5.0)
slope_constraint = CubicSplineSlopeConstraint(traj, du_bounds)

prob = DirectTrajOptProblem(
    traj, obj, integrators;
    constraints=[slope_constraint]
)
```

# Note
Returns 4 constraints per drive per segment:
- 2 endpoints (k and k+1)
- Up to 1 interior critical point (padded to fixed size with midpoint)
- Each checked against both bounds (lower and upper)
"""
function CubicSplineSlopeConstraint(
    traj::NamedTrajectory,
    du_bounds::Tuple{Real,Real};
    control_name::Symbol = :u,
    segments::AbstractVector{Int} = 1:(traj.N-1),
)
    du_name = Symbol(:d, control_name)

    # Verify required components exist
    @assert control_name ∈ traj.names "Trajectory must have $control_name component"
    @assert du_name ∈ traj.names "Trajectory must have $du_name component"
    @assert :Δt ∈ traj.names "Trajectory must have Δt component"

    du_min, du_max = du_bounds

    n_drives = traj.dims[control_name]

    function slope_constraint(varsₖ, varsₖ₊₁, (du_bounds, Δt))
        # Split variables: assume [u..., du...] concatenated
        n = length(varsₖ) ÷ 2
        u_k = varsₖ[1:n]
        du_k = varsₖ[(n+1):end]
        u_kp1 = varsₖ₊₁[1:n]
        du_kp1 = varsₖ₊₁[(n+1):end]

        T = eltype(varsₖ)
        constraints = T[]

        for i = 1:n
            # Find critical point of derivative (where d²u/dt² = 0)
            τ_crits = find_slope_critical_points(u_k[i], u_kp1[i], du_k[i], du_kp1[i], Δt)

            # Collect evaluation points: endpoints + critical point
            τ_evals = T[0.0, 1.0]
            if !isempty(τ_crits)
                push!(τ_evals, τ_crits[1])
            else
                # Pad with midpoint for fixed dimension
                push!(τ_evals, 0.5)
            end

            # Evaluate slope at each point and constrain
            for τ in τ_evals
                slope =
                    evaluate_hermite_derivative(τ, u_k[i], u_kp1[i], du_k[i], du_kp1[i], Δt)

                # Two inequality constraints per evaluation point. Backend convention:
                # feasible ⟺ residual ≤ 0 (matches the AL's max(0, μ+ρc) penalty and the
                # bounds rows; the previous sign was inverted, making feasible slopes read
                # as violations of |du_bound| and crushing the pulse flat).
                push!(constraints, du_bounds[1] - slope)  # du_min − slope ≤ 0 ⟺ slope ≥ du_min
                push!(constraints, slope - du_bounds[2])  # slope − du_max ≤ 0 ⟺ slope ≤ du_max
            end
        end

        # Returns 3 eval points × 2 bounds × n_drives constraints
        return constraints
    end

    params = [(du_bounds, traj[k].timestep) for k in segments]

    return NonlinearSegmentConstraint(
        slope_constraint,
        [control_name, du_name],
        traj,
        params;
        equality = false,
        segments = segments,
    )
end

"""
    CubicSplineSlopeConstraint(
        traj::NamedTrajectory,
        du_bounds::Vector{Tuple{Real, Real}};
        control_names::Vector{Symbol}=[:u],
        segments::AbstractVector{Int}=1:traj.N-1
    )

Multi-drive version for constraining slopes of multiple control variables.

# Arguments
- `traj`: Trajectory with spline pulse
- `du_bounds`: Vector of (du_min, du_max) for each drive

# Example
```julia
# Different bounds for each drive
du_bounds = [(-5.0, 5.0), (-3.0, 3.0)]
control_names = [:u1, :u2]
slope_constraint = CubicSplineSlopeConstraint(traj, du_bounds; control_names)
```
"""
function CubicSplineSlopeConstraint(
    traj::NamedTrajectory,
    du_bounds::Vector{Tuple{T,T}};
    control_names::Vector{Symbol} = [:u],
    segments::AbstractVector{Int} = 1:(traj.N-1),
) where {T<:Real}

    du_names = [Symbol(:d, name) for name in control_names]

    # Verify all required components exist
    for (name, du_name) in zip(control_names, du_names)
        @assert name ∈ traj.names "Trajectory must have $name component"
        @assert du_name ∈ traj.names "Trajectory must have $du_name component"
    end
    @assert :Δt ∈ traj.names "Trajectory must have Δt component"

    # Calculate total number of drives across all control variables
    total_drives = sum(traj.dims[name] for name in control_names)

    function slope_constraint_multi(varsₖ, varsₖ₊₁, (du_bounds_vec, Δt))
        # Split variables: [all u's..., all du's...]  
        n = length(varsₖ) ÷ 2
        u_k = varsₖ[1:n]
        du_k = varsₖ[(n+1):end]
        u_kp1 = varsₖ₊₁[1:n]
        du_kp1 = varsₖ₊₁[(n+1):end]

        ElemType = eltype(varsₖ)
        constraints = ElemType[]

        # Iterate over all drives
        for i = 1:n
            # Find the corresponding bounds
            bounds = du_bounds_vec[i]

            τ_crits = find_slope_critical_points(u_k[i], u_kp1[i], du_k[i], du_kp1[i], Δt)

            τ_evals = ElemType[0.0, 1.0]
            if !isempty(τ_crits)
                push!(τ_evals, τ_crits[1])
            else
                push!(τ_evals, 0.5)
            end

            for τ in τ_evals
                slope =
                    evaluate_hermite_derivative(τ, u_k[i], u_kp1[i], du_k[i], du_kp1[i], Δt)
                # feasible ⟺ residual ≤ 0 (backend convention; see CubicSplineSlopeConstraint)
                push!(constraints, bounds[1] - slope)  # du_min − slope ≤ 0
                push!(constraints, slope - bounds[2])  # slope − du_max ≤ 0
            end
        end

        return constraints
    end

    params = [(du_bounds, traj[k].timestep) for k in segments]

    return NonlinearSegmentConstraint(
        slope_constraint_multi,
        vcat(control_names, du_names),
        traj,
        params;
        equality = false,
        segments = segments,
    )
end

# ----------------------------------------------------------------------------- #
# Tests
# ----------------------------------------------------------------------------- #

@testitem "Hermite spline math utilities" begin
    # Define functions locally for testing
    function hermite_basis_functions(τ)
        τ² = τ * τ
        τ³ = τ² * τ
        h00 = 2τ³ - 3τ² + 1
        h10 = τ³ - 2τ² + τ
        h01 = -2τ³ + 3τ²
        h11 = τ³ - τ²
        return (h00, h10, h01, h11)
    end

    function evaluate_hermite_spline(τ, u_k, u_kp1, du_k, du_kp1, Δt)
        h00, h10, h01, h11 = hermite_basis_functions(τ)
        return h00 * u_k + h10 * Δt * du_k + h01 * u_kp1 + h11 * Δt * du_kp1
    end

    # Test basis functions at key points
    h00_0, h10_0, h01_0, h11_0 = hermite_basis_functions(0.0)
    @test h00_0 ≈ 1.0
    @test h10_0 ≈ 0.0
    @test h01_0 ≈ 0.0
    @test h11_0 ≈ 0.0

    h00_1, h10_1, h01_1, h11_1 = hermite_basis_functions(1.0)
    @test h00_1 ≈ 0.0
    @test h10_1 ≈ 0.0
    @test h01_1 ≈ 1.0
    @test h11_1 ≈ 0.0

    # Test spline evaluation
    u_k, u_kp1 = 0.0, 1.0
    du_k, du_kp1 = 0.0, 0.0
    Δt = 1.0

    @test evaluate_hermite_spline(0.0, u_k, u_kp1, du_k, du_kp1, Δt) ≈ u_k
    @test evaluate_hermite_spline(1.0, u_k, u_kp1, du_k, du_kp1, Δt) ≈ u_kp1
    @test evaluate_hermite_spline(0.5, u_k, u_kp1, du_k, du_kp1, Δt) ≈ 0.5
end

@testitem "Find cubic critical points" begin
    # Define functions locally for testing
    function find_cubic_critical_points(u_k, u_kp1, du_k, du_kp1, Δt)
        a = 6*u_k + 3*Δt*du_k - 6*u_kp1 + 3*Δt*du_kp1
        b = -6*u_k - 4*Δt*du_k + 6*u_kp1 - 2*Δt*du_kp1
        c = Δt * du_k

        critical_points = Float64[]

        if abs(a) < 1e-12
            if abs(b) > 1e-12
                τ = -c / b
                if 0 < τ < 1
                    push!(critical_points, τ)
                end
            end
        else
            discriminant = b^2 - 4*a*c
            if discriminant >= 0
                sqrt_disc = sqrt(discriminant)
                τ1 = (-b + sqrt_disc) / (2*a)
                τ2 = (-b - sqrt_disc) / (2*a)
                if 0 < τ1 < 1
                    push!(critical_points, τ1)
                end
                if 0 < τ2 < 1
                    push!(critical_points, τ2)
                end
            end
        end

        return critical_points
    end

    function hermite_derivative_basis(τ)
        τ² = τ * τ
        dh00 = 6τ² - 6τ
        dh10 = 3τ² - 4τ + 1
        dh01 = -6τ² + 6τ
        dh11 = 3τ² - 2τ
        return (dh00, dh10, dh01, dh11)
    end

    # Linear spline (no overshoot) - should have no critical points
    u_k, u_kp1 = 0.0, 1.0
    du_k, du_kp1 = 1.0, 1.0
    Δt = 1.0

    τ_crits = find_cubic_critical_points(u_k, u_kp1, du_k, du_kp1, Δt)
    @test length(τ_crits) == 0  # Monotonic

    # Spline with overshoot - should have critical point
    u_k, u_kp1 = 0.0, 1.0
    du_k, du_kp1 = 3.0, -3.0  # Large opposite slopes → overshoot

    τ_crits = find_cubic_critical_points(u_k, u_kp1, du_k, du_kp1, Δt)
    @test length(τ_crits) >= 1  # Has at least one critical point
    @test all(0 .< τ_crits .< 1)  # In valid range

    # Verify critical point is actually an extremum
    for τ in τ_crits
        dh00, dh10, dh01, dh11 = hermite_derivative_basis(τ)
        du_dτ = dh00 * u_k + dh10 * Δt * du_k + dh01 * u_kp1 + dh11 * Δt * du_kp1
        @test abs(du_dτ) < 1e-10  # Derivative should be ~0 at critical point
    end
end

@testitem "CubicSplineExtremaConstraint - basic" begin
    using NamedTrajectories
    using LinearAlgebra
    using Piccolo: CubicSplineExtremaConstraint
    using Piccolo: NonlinearSegmentConstraint
    using DirectTrajOpt.CommonInterface: evaluate!

    N = 10
    T = 5.0
    times = collect(range(0, T, length = N))

    # Create trajectory with spline controls
    traj = NamedTrajectory(
        (u = 0.5 * randn(2, N), du = 0.1 * randn(2, N), Δt = fill(T/(N-1), N));
        timestep = :Δt,
        controls = :u,
    )

    # Apply bound constraints
    u_min, u_max = -1.0, 1.0
    constraint = CubicSplineExtremaConstraint(traj, (u_min, u_max))

    @test constraint isa NonlinearSegmentConstraint
    @test constraint.equality == false  # Inequality
    @test length(constraint.segments) == N - 1

    # Evaluate constraint
    values = zeros(constraint.dim)
    evaluate!(values, constraint, traj)

    # All constraints should be negative (feasible) or zero if trajectory is within bounds
    # (Initial random trajectory may not satisfy bounds)
    @test length(values) > 0
end

@testitem "CubicSplineExtremaConstraint - per drive bounds" begin
    using NamedTrajectories
    using Piccolo: CubicSplineExtremaConstraint
    using DirectTrajOpt.CommonInterface: evaluate!

    N = 5
    n_drives = 3
    traj = NamedTrajectory(
        (u = randn(n_drives, N), du = randn(n_drives, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
    )

    # Different bounds per drive
    bounds = [(-1.0, 1.0), (-2.0, 2.0), (-0.5, 0.5)]

    constraint = CubicSplineExtremaConstraint(traj, bounds)

    values = zeros(constraint.dim)
    evaluate!(values, constraint, traj)
    @test length(values) > 0
end

@testitem "CubicSplineSufficientBoundConstraint - basic" begin
    using NamedTrajectories
    using Piccolo: CubicSplineSufficientBoundConstraint
    using Piccolo: NonlinearSegmentConstraint
    using DirectTrajOpt.CommonInterface: evaluate!

    N = 10
    traj = NamedTrajectory(
        (u = randn(2, N), du = 0.1 * randn(2, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
    )

    u_min, u_max = -5.0, 5.0
    constraint =
        CubicSplineSufficientBoundConstraint(traj, (u_min, u_max); safety_factor = 2.0)

    @test constraint isa NonlinearSegmentConstraint
    @test constraint.equality == false

    # Evaluate
    values = zeros(constraint.dim)
    evaluate!(values, constraint, traj)
    @test length(values) == 2 * (N - 1)  # 2 drives × (N-1) segments
end

@testitem "CubicSplineSufficientBoundConstraint - safety factor effect" begin
    using NamedTrajectories
    using Piccolo: CubicSplineSufficientBoundConstraint
    using DirectTrajOpt.CommonInterface: evaluate!

    N = 5
    traj = NamedTrajectory(
        (u = zeros(1, N), du = ones(1, N), Δt = fill(1.0, N));
        timestep = :Δt,
        controls = :u,
    )

    u_min, u_max = -10.0, 10.0

    # With higher safety factor, constraint should be less restrictive
    constraint_conservative =
        CubicSplineSufficientBoundConstraint(traj, (u_min, u_max); safety_factor = 1.0)
    constraint_moderate =
        CubicSplineSufficientBoundConstraint(traj, (u_min, u_max); safety_factor = 3.0)

    values_conservative = zeros(constraint_conservative.dim)
    values_moderate = zeros(constraint_moderate.dim)

    evaluate!(values_conservative, constraint_conservative, traj)
    evaluate!(values_moderate, constraint_moderate, traj)

    # Higher safety factor → larger allowed slopes → SMALLER residuals under the
    # backend convention (feasible ⟺ residual = |du| − max_slope ≤ 0): a larger
    # max_slope subtracts more. (Expectation flipped with the 1eb1b84 sign fix.)
    @test all(values_moderate .<= values_conservative .+ 1e-10)
end

@testitem "Spline constraints - integration test" begin
    using NamedTrajectories
    using LinearAlgebra
    using SparseArrays
    using Piccolo: CubicSplineExtremaConstraint, CubicSplineSufficientBoundConstraint
    using DirectTrajOpt.CommonInterface: evaluate!
    using DirectTrajOpt.Constraints: jacobian!, get_full_jacobian

    # Create a trajectory that violates bounds
    N = 8
    traj = NamedTrajectory(
        (u = 5.0 * randn(1, N), du = 10.0 * randn(1, N), Δt = fill(0.5, N));
        timestep = :Δt,
        controls = :u,
        bounds = (u = (-3.0, 3.0), du = (-5.0, 5.0)),
    )

    # Apply extrema constraint
    extrema_con = CubicSplineExtremaConstraint(traj, (-3.0, 3.0))

    # Apply sufficient condition constraint
    sufficient_con =
        CubicSplineSufficientBoundConstraint(traj, (-3.0, 3.0); safety_factor = 2.0)

    # Evaluate both
    vals_extrema = zeros(extrema_con.dim)
    vals_sufficient = zeros(sufficient_con.dim)

    evaluate!(vals_extrema, extrema_con, traj)
    evaluate!(vals_sufficient, sufficient_con, traj)

    # Both should produce valid constraint vectors
    @test length(vals_extrema) > 0
    @test length(vals_sufficient) > 0
    @test all(isfinite.(vals_extrema))
    @test all(isfinite.(vals_sufficient))

    # Jacobian should be computable
    jacobian!(extrema_con, traj)
    jacobian!(sufficient_con, traj)

    ∂g_extrema = get_full_jacobian(extrema_con, traj)
    ∂g_sufficient = get_full_jacobian(sufficient_con, traj)

    @test nnz(∂g_extrema) > 0
    @test nnz(∂g_sufficient) > 0
    @test all(isfinite.(∂g_extrema.nzval))
    @test all(isfinite.(∂g_sufficient.nzval))
end

@testitem "Slope constraint - basic functionality" begin
    using NamedTrajectories
    using Piccolo: CubicSplineSlopeConstraint
    using DirectTrajOpt.CommonInterface: evaluate!

    # Create trajectory with large derivatives
    N = 6
    traj = NamedTrajectory(
        (u = randn(1, N), du = 10.0 * randn(1, N), Δt = fill(0.5, N));
        timestep = :Δt,
        controls = :u,
    )

    # Create slope constraint
    du_bounds = (-5.0, 5.0)
    slope_con = CubicSplineSlopeConstraint(traj, du_bounds)

    # Evaluate constraints
    g = zeros(slope_con.dim)
    evaluate!(g, slope_con, traj)

    # Should return 3 eval points × 2 bounds × 1 drive × (N-1) segments
    @test length(g) == 3 * 2 * 1 * (N-1)
    @test all(isfinite.(g))
end
