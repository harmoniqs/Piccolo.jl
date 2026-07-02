module QuantumConstraints

using ..QuantumObjectives
using ..QuantumObjectives: ket_fidelity_loss, unitary_fidelity_loss, coherent_ket_fidelity

using DirectTrajOpt
using LinearAlgebra
using NamedTrajectories
using TestItems
using ...Quantum

export FinalKetFidelityConstraint
export FinalKetFreePhaseConstraint
export FinalUnitaryFidelityConstraint
export FinalCoherentKetFidelityConstraint
export FinalDensityFidelityConstraint
export LeakageConstraint
export BoundStateL2Constraint

# ---------------------------------------------------------
#                        Kets
# ---------------------------------------------------------

function FinalKetFidelityConstraint(
    ψ_goal::AbstractVector{<:Complex{Float64}},
    ψ̃_name::Symbol,
    final_fidelity::Float64,
    traj::NamedTrajectory,
)
    terminal_constraint = ψ̃ -> [final_fidelity - ket_fidelity_loss(ψ̃, ψ_goal)]

    return NonlinearKnotPointConstraint(
        terminal_constraint,
        ψ̃_name,
        traj,
        equality = false,
        times = [traj.N],
    )
end

"""
    FinalKetFreePhaseConstraint(goal_fn, ψ̃_name, θ_names, final_fidelity, traj)

Free-phase version of `FinalKetFidelityConstraint` for single-ket trajectories.

`goal_fn(θ)` returns the phase-rotated goal ket. The constraint enforces
F(θ) = |⟨goal(θ)|ψ_final⟩|² ≥ final_fidelity, where θ are optimizable
per-subsystem phase variables stored in the trajectory's global_data.

Uses `NonlinearGlobalKnotPointConstraint` to include global phase variables
in the constraint evaluation.
"""
function FinalKetFreePhaseConstraint(
    goal_fn::Function,
    ψ̃_name::Symbol,
    θ_names::AbstractVector{Symbol},
    final_fidelity::Float64,
    traj::NamedTrajectory,
)
    d_state = traj.dims[ψ̃_name]

    function terminal_constraint(z)
        ψ̃ = z[1:d_state]
        θ = z[(d_state+1):end]
        phased_goal = goal_fn(θ)
        return [final_fidelity - ket_fidelity_loss(ψ̃, phased_goal)]
    end

    return NonlinearGlobalKnotPointConstraint(
        terminal_constraint,
        [ψ̃_name],
        collect(θ_names),
        traj,
        equality = false,
        times = [traj.N],
    )
end

# ---------------------------------------------------------
#                  Coherent Ket Fidelity
# ---------------------------------------------------------

"""
    FinalCoherentKetFidelityConstraint(ψ_goals, ψ̃_names, final_fidelity, traj)

Create a final fidelity constraint using coherent ket fidelity across multiple states.

Coherent fidelity: F = |1/n ∑ᵢ ⟨ψᵢ_goal|ψᵢ⟩|²

This constraint enforces that all state overlaps have aligned phases, which is 
essential when implementing a gate via multiple state transfers (e.g., MultiKetTrajectory).

# Arguments
- `ψ_goals::Vector{<:AbstractVector{<:Complex}}`: Target ket states
- `ψ̃_names::Vector{Symbol}`: Names of isomorphic state variables in trajectory
- `final_fidelity::Float64`: Minimum fidelity threshold (constraint: F ≥ final_fidelity)
- `traj::NamedTrajectory`: The trajectory

# Example
```julia
# For implementing X gate via |0⟩→|1⟩ and |1⟩→|0⟩
goals = [ComplexF64[0, 1], ComplexF64[1, 0]]
names = [:ψ̃1, :ψ̃2]
constraint = FinalCoherentKetFidelityConstraint(goals, names, 0.99, traj)
```
"""
function FinalCoherentKetFidelityConstraint(
    ψ_goals::Vector{<:AbstractVector{<:Complex}},
    ψ̃_names::Vector{Symbol},
    final_fidelity::Float64,
    traj::NamedTrajectory,
)
    n_states = length(ψ_goals)
    @assert length(ψ̃_names) == n_states "Number of names must match number of goals"

    # Convert goals to ComplexF64
    goals = [ComplexF64.(g) for g in ψ_goals]

    # Get component info for extracting states from concatenated vector
    state_dims = [traj.dims[name] for name in ψ̃_names]

    function terminal_constraint(z_terminal)
        # Extract each state from the concatenated vector
        ψ̃s = Vector{Vector{eltype(z_terminal)}}(undef, n_states)
        offset = 0
        for i = 1:n_states
            ψ̃s[i] = z_terminal[(offset+1):(offset+state_dims[i])]
            offset += state_dims[i]
        end

        # Constraint: final_fidelity - F_coherent ≤ 0
        return [final_fidelity - coherent_ket_fidelity(ψ̃s, goals)]
    end

    return NonlinearKnotPointConstraint(
        terminal_constraint,
        ψ̃_names,
        traj,
        equality = false,
        times = [traj.N],
    )
end

"""
    FinalCoherentKetFidelityConstraint(goals_fn, ψ̃_names, θ_names, final_fidelity, traj)

Free-phase version: `goals_fn(θ)` returns phase-adjusted goal kets.
Uses `NonlinearGlobalKnotPointConstraint` to include global phase variables.
"""
function FinalCoherentKetFidelityConstraint(
    goals_fn::Function,
    ψ̃_names::Vector{Symbol},
    θ_names::AbstractVector{Symbol},
    final_fidelity::Float64,
    traj::NamedTrajectory,
)
    n_states = length(ψ̃_names)
    state_dims = [traj.dims[name] for name in ψ̃_names]
    total_state_dim = sum(state_dims)

    function terminal_constraint(z)
        x = z[1:total_state_dim]
        θ = z[(total_state_dim+1):end]

        # Extract each state from the concatenated vector
        ψ̃s = Vector{Vector{eltype(x)}}(undef, n_states)
        offset = 0
        for i = 1:n_states
            ψ̃s[i] = x[(offset+1):(offset+state_dims[i])]
            offset += state_dims[i]
        end

        phased_goals = goals_fn(θ)
        return [final_fidelity - coherent_ket_fidelity(ψ̃s, phased_goals)]
    end

    return NonlinearGlobalKnotPointConstraint(
        terminal_constraint,
        ψ̃_names,
        θ_names,
        traj,
        equality = false,
        times = [traj.N],
    )
end

# ---------------------------------------------------------
#                        Unitaries
# ---------------------------------------------------------

function FinalUnitaryFidelityConstraint(
    U_goal::AbstractPiccoloOperator,
    Ũ⃗_name::Symbol,
    final_fidelity::Float64,
    traj::NamedTrajectory,
)
    terminal_constraint = Ũ⃗ -> [final_fidelity - unitary_fidelity_loss(Ũ⃗, U_goal)]

    return NonlinearKnotPointConstraint(
        terminal_constraint,
        Ũ⃗_name,
        traj,
        equality = false,
        times = [traj.N],
    )
end

function FinalUnitaryFidelityConstraint(
    U_goal::Function,
    Ũ⃗_name::Symbol,
    θ_names::AbstractVector{Symbol},
    final_fidelity::Float64,
    traj::NamedTrajectory,
)
    θ_dim = sum(traj.global_dims[n] for n in θ_names)
    function terminal_constraint(z)
        Ũ⃗, θ = z[1:(end-θ_dim)], z[(end-θ_dim+1):end]
        return [final_fidelity - unitary_fidelity_loss(Ũ⃗, U_goal(θ))]
    end

    return NonlinearGlobalKnotPointConstraint(
        terminal_constraint,
        Ũ⃗_name,
        θ_names,
        traj,
        equality = false,
        times = [traj.N],
    )
end

# ---------------------------------------------------------
#                     Density Matrices
# ---------------------------------------------------------

@doc raw"""
    FinalDensityFidelityConstraint(ρ_goal, ρ̃_name, final_fidelity, traj)

Enforce a minimum state-transfer fidelity on the final density matrix of a
`DensityTrajectory` knot-point variable.

The trajectory stores the density matrix in the compact real isomorphism
representation (see `density_to_compact_iso` / `compact_iso_to_density`): a
real vector of length ``n^2`` for a Hermitian ``n \times n`` density matrix.

The constraint enforces

```math
F(\rho_{\text{final}}) = \mathrm{Re}\,\mathrm{tr}(\rho_{\text{goal}}\, \rho_{\text{final}}) \geq F_{\text{threshold}},
```

which is the exact state-transfer fidelity when ``\rho_{\text{goal}}`` is pure.
``F`` is real-linear (but not complex-linear) in the compact-iso vector, so the
Jacobian is a constant row vector; `NonlinearKnotPointConstraint` recovers it
via automatic differentiation at constraint construction.

# Arguments
- `ρ_goal::AbstractMatrix{<:Number}`: Target density matrix (``n \times n``,
  Hermitian).
- `ρ̃_name::Symbol`: Name of the compact-iso density component in the
  trajectory (e.g. `:ρ⃗̃`).
- `final_fidelity::Float64`: Minimum fidelity threshold
  (constraint: ``F \geq F_{\text{threshold}}``).
- `traj::NamedTrajectory`: The trajectory carrying the density component.
"""
function FinalDensityFidelityConstraint(
    ρ_goal::AbstractMatrix{<:Number},
    ρ̃_name::Symbol,
    final_fidelity::Float64,
    traj::NamedTrajectory,
)
    ρ_goal_c = Matrix{ComplexF64}(ρ_goal)

    function terminal_constraint(ρ̃)
        ρ = compact_iso_to_density(ρ̃)
        F = real(tr(ρ_goal_c * ρ))
        return [final_fidelity - F]
    end

    return NonlinearKnotPointConstraint(
        terminal_constraint,
        ρ̃_name,
        traj,
        equality = false,
        times = [traj.N],
    )
end

# ---------------------------------------------------------
# Leakage Constraint
# ---------------------------------------------------------

"""
    LeakageConstraint(value, indices, name, traj::NamedTrajectory)

Construct a `KnotPointConstraint` that bounds leakage of `name` at the knot points specified by `times` at any `indices` that are outside the computational subspace.

"""
function LeakageConstraint(
    value::Float64,
    indices::AbstractVector{Int},
    name::Symbol,
    traj::NamedTrajectory;
    times = 1:(traj.N),
)
    leakage_constraint(x) = abs2.(x[indices]) .- value

    return NonlinearKnotPointConstraint(
        leakage_constraint,
        name,
        traj,
        equality = false,
        times = times,
    )
end

# ---------------------------------------------------------
# Bound State L2 Constraint
# ---------------------------------------------------------

"""
    _compact_iso_index_map(n::Int)

Precompute index arrays for the compact density isomorphism of an `n × n`
Hermitian matrix. Returns vectors (not Dicts) for closure-capture performance:
- `re_idx_pairs`: Re index for each off-diagonal pair (j<k), length n(n-1)/2
- `im_idx_pairs`: Im index for each off-diagonal pair, same length
- `diag_j_pairs`: diagonal index ρ_jj for each pair
- `diag_k_pairs`: diagonal index ρ_kk for each pair
"""
function _compact_iso_index_map(n::Int)
    # Re upper triangle indices (column-major: j <= k)
    re_map = Matrix{Int}(undef, n, n)
    idx = 0
    for k = 1:n, j = 1:k
        idx += 1
        re_map[j, k] = idx
    end

    # Im strict upper triangle indices (column-major: j < k)
    im_map = Matrix{Int}(undef, n, n)
    for k = 2:n, j = 1:(k-1)
        idx += 1
        im_map[j, k] = idx
    end

    # Build flat arrays for each off-diagonal pair
    n_pairs = n * (n - 1) ÷ 2
    re_idx_pairs = Vector{Int}(undef, n_pairs)
    im_idx_pairs = Vector{Int}(undef, n_pairs)
    diag_j_pairs = Vector{Int}(undef, n_pairs)
    diag_k_pairs = Vector{Int}(undef, n_pairs)
    p = 0
    for k = 2:n, j = 1:(k-1)
        p += 1
        re_idx_pairs[p] = re_map[j, k]
        im_idx_pairs[p] = im_map[j, k]
        diag_j_pairs[p] = re_map[j, j]
        diag_k_pairs[p] = re_map[k, k]
    end

    return re_idx_pairs, im_idx_pairs, diag_j_pairs, diag_k_pairs
end

"""
    BoundStateL2Constraint(name, traj, iso_layout; times=1:traj.N)

Constrain each complex component's magnitude via a layout-dependent nonlinear
inequality constraint.

`iso_layout` determines the Re/Im index pairing:
- `:block` — ket iso `ψ̃ = [Re(ψ); Im(ψ)]`, pairs `(k, n+k)` for `k=1:n`.
  Constraint: `Re² + Im² - 1 ≤ 0` per complex entry.
- `:interleaved_columns` — unitary iso `Ũ⃗`, per-column `[Re(col); Im(col)]`,
  pairs `(offset+j, offset+d+j)` within each `2d`-stride column block.
  Constraint: `Re² + Im² - 1 ≤ 0` per complex entry.
- `:compact_density` — compact density iso `ρ⃗̃`, with Re upper-triangle block
  followed by Im strict-upper-triangle block. Enforces the Cauchy-Schwarz
  bound `Re(ρ_jk)² + Im(ρ_jk)² - ρ_jj · ρ_kk ≤ 0` per off-diagonal pair.
"""
function BoundStateL2Constraint(
    name::Symbol,
    traj::NamedTrajectory,
    iso_layout::Symbol;
    times = 1:traj.N,
)
    dim = traj.dims[name]

    if iso_layout == :block
        n = dim ÷ 2
        iseven(dim) || throw(ArgumentError("block layout expects even dim; got $dim"))
        function block_constraint(x)
            re = @view x[1:n]
            im = @view x[(n+1):(2n)]
            return re .^ 2 .+ im .^ 2 .- 1.0
        end
        return NonlinearKnotPointConstraint(
            block_constraint,
            name,
            traj;
            equality = false,
            times = times,
        )
    elseif iso_layout == :interleaved_columns
        d = isqrt(dim ÷ 2)
        dim == 2 * d^2 ||
            throw(ArgumentError("interleaved_columns expects dim = 2d²; got $dim"))
        n_complex = d * d
        function interleaved_constraint(x)
            result = Vector{eltype(x)}(undef, n_complex)
            idx = 1
            for col = 0:(d-1)
                offset = col * 2d
                for row = 1:d
                    re = x[offset+row]
                    im = x[offset+d+row]
                    result[idx] = re^2 + im^2 - 1.0
                    idx += 1
                end
            end
            return result
        end
        return NonlinearKnotPointConstraint(
            interleaved_constraint,
            name,
            traj;
            equality = false,
            times = times,
        )
    elseif iso_layout == :compact_density
        n = isqrt(dim)
        dim == n^2 || throw(ArgumentError("compact_density expects dim = n²; got $dim"))
        re_idx, im_idx, dj_idx, dk_idx = _compact_iso_index_map(n)
        n_pairs = n * (n - 1) ÷ 2
        function density_constraint(x)
            result = Vector{eltype(x)}(undef, n_pairs)
            for p = 1:n_pairs
                re = x[re_idx[p]]
                im = x[im_idx[p]]
                ρ_jj = x[dj_idx[p]]
                ρ_kk = x[dk_idx[p]]
                result[p] = re^2 + im^2 - ρ_jj * ρ_kk
            end
            return result
        end
        return NonlinearKnotPointConstraint(
            density_constraint,
            name,
            traj;
            equality = false,
            times = times,
        )
    else
        throw(
            ArgumentError(
                "Unknown iso_layout :$iso_layout. " *
                "Expected :block, :interleaved_columns, or :compact_density.",
            ),
        )
    end
end

# ---------------------------------------------------------
#                         Tests
# ---------------------------------------------------------

@testitem "FinalDensityFidelityConstraint via MinimumTimeProblem" begin
    using DirectTrajOpt
    using LinearAlgebra
    using NamedTrajectories

    # 2-level open system, no dissipation, σx drive: |0⟩⟨0| → |1⟩⟨1|
    sys = OpenQuantumSystem(zeros(ComplexF64, 2, 2), [PAULIS.X], [1.0])

    ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρg = ComplexF64[0.0 0.0; 0.0 1.0]

    T = 10.0
    N = 50

    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    qtraj = DensityTrajectory(sys, pulse, ρ0, ρg)

    qcp_smooth = SmoothPulseProblem(qtraj, N; Q = 100.0, R = 1e-2, Δt_bounds = (0.01, 0.5))
    solve!(qcp_smooth; max_iter = 100, verbose = false, print_level = 1)

    # Convert to minimum-time — this is the path the dispatch stub used to block
    qcp_mintime = MinimumTimeProblem(qcp_smooth; final_fidelity = 0.95, D = 50.0)

    @test qcp_mintime isa QuantumControlProblem{<:DensityTrajectory}

    # Solve minimum-time problem
    solve!(qcp_mintime; max_iter = 100, verbose = false, print_level = 1)

    # Roll out with the optimized pulse and verify the fidelity constraint is respected
    traj = get_trajectory(qcp_mintime)
    ρ̃_final = traj[end][:ρ⃗̃]
    ρ_final = compact_iso_to_density(ρ̃_final)

    # Final density matrix should be Hermitian and trace-preserving
    @test ρ_final ≈ ρ_final' atol = 1e-6
    @test real(tr(ρ_final)) ≈ 1.0 atol = 1e-2

    # The constraint enforces F = Re tr(ρ_goal · ρ_final) ≥ 0.95 at the final knot
    fid = real(tr(ρg * ρ_final))
    @test fid ≥ 0.94  # small tolerance for solver feasibility slack
end

@testitem "FinalDensityFidelityConstraint direct construction" begin
    using DirectTrajOpt
    using LinearAlgebra

    # Smoke test the constructor independently of the minimum-time path.
    sys = OpenQuantumSystem(zeros(ComplexF64, 2, 2), [PAULIS.X], [1.0])

    ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρg = ComplexF64[0.5 0.5; 0.5 0.5]  # |+⟩⟨+|

    T = 2.0
    N = 10
    pulse = ZeroOrderPulse(0.1 * randn(1, N), collect(range(0.0, T, length = N)))
    qtraj = DensityTrajectory(sys, pulse, ρ0, ρg)

    qcp = SmoothPulseProblem(qtraj, N; Q = 10.0, R = 1e-2)
    traj = get_trajectory(qcp)

    constraint = FinalDensityFidelityConstraint(ρg, :ρ⃗̃, 0.9, traj)

    @test constraint isa DirectTrajOpt.AbstractNonlinearConstraint
    # Inequality constraint, evaluated only at final knot
    @test constraint.equality == false
    @test constraint.times == [traj.N]
    @test constraint.g_dim == 1
end

@testitem "BoundStateL2Constraint block layout" begin
    using NamedTrajectories
    using DirectTrajOpt

    N = 5
    # 4-dim iso-vec = 2 complex components (block layout: [Re; Im])
    traj = NamedTrajectory(
        (ψ̃ = rand(4, N), u = rand(1, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
    )

    con = BoundStateL2Constraint(:ψ̃, traj, :block)
    @test con isa DirectTrajOpt.AbstractNonlinearConstraint
    @test con.equality == false
    # 2 complex components → 2 inequality constraints
    @test con.g_dim == 2

    # Evaluate constraint at a known point: ψ̃ = [0.6, 0.8, 0.0, 0.0]
    # Complex components: z₁ = 0.6+0i, z₂ = 0.8+0i
    # |z₁|² - 1 = -0.64, |z₂|² - 1 = -0.36 (both satisfied)
    x = [0.6, 0.8, 0.0, 0.0]
    g = con.g(x, nothing)
    @test g ≈ [0.36 - 1.0, 0.64 - 1.0]

    # Point violating constraint: ψ̃ = [0.8, 0.0, 0.8, 0.0]
    # z₁ = 0.8+0.8i → |z₁|² = 1.28 > 1
    x2 = [0.8, 0.0, 0.8, 0.0]
    g2 = con.g(x2, nothing)
    @test g2[1] > 0  # violated
end

@testitem "BoundStateL2Constraint interleaved_columns layout" begin
    using NamedTrajectories
    using DirectTrajOpt

    N = 5
    # 2×2 unitary → 8-dim iso-vec (interleaved columns: [Re(col₁); Im(col₁); Re(col₂); Im(col₂)])
    traj = NamedTrajectory(
        (Ũ⃗ = rand(8, N), u = rand(1, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
    )

    con = BoundStateL2Constraint(:Ũ⃗, traj, :interleaved_columns)
    @test con isa DirectTrajOpt.AbstractNonlinearConstraint
    @test con.equality == false
    # 2×2 = 4 complex entries → 4 inequality constraints
    @test con.g_dim == 4

    # Identity matrix I₂: iso_vec = [1, 0, 0, 0, 0, 1, 0, 0]
    # col₀: Re=[1,0], Im=[0,0] → z₁=1+0i, z₂=0+0i → |z|²-1 = [0, -1]
    # col₁: Re=[0,1], Im=[0,0] → z₃=0+0i, z₄=1+0i → |z|²-1 = [-1, 0]
    x = [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    g = con.g(x, nothing)
    @test g ≈ [0.0, -1.0, -1.0, 0.0]
end

@testitem "BoundStateL2Constraint invalid layout" begin
    using NamedTrajectories
    using DirectTrajOpt

    N = 5
    traj = NamedTrajectory(
        (x = rand(4, N), u = rand(1, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
    )

    @test_throws ArgumentError BoundStateL2Constraint(:x, traj, :invalid)
end

@testitem "BoundStateL2Constraint compact_density layout" begin
    using NamedTrajectories
    using DirectTrajOpt

    N = 5
    # 2×2 density → compact iso dim = 4
    # Layout: [Re(ρ₁₁), Re(ρ₁₂), Re(ρ₂₂), Im(ρ₁₂)]
    traj = NamedTrajectory(
        (ρ⃗̃ = rand(4, N), u = rand(1, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
    )

    con = BoundStateL2Constraint(:ρ⃗̃, traj, :compact_density)
    @test con isa DirectTrajOpt.AbstractNonlinearConstraint
    @test con.equality == false
    # 2×2 has 1 off-diagonal pair → 1 Cauchy-Schwarz constraint
    @test con.g_dim == 1

    # Identity: ρ = I/2 → compact = [0.5, 0.0, 0.5, 0.0]
    # Cauchy-Schwarz: 0² + 0² - 0.5*0.5 = -0.25 ≤ 0 (satisfied)
    g = con.g([0.5, 0.0, 0.5, 0.0], nothing)
    @test g ≈ [-0.25]

    # Pure state |0⟩⟨0|: ρ = [1 0; 0 0] → compact = [1.0, 0.0, 0.0, 0.0]
    # Cauchy-Schwarz: 0² + 0² - 1.0*0.0 = 0.0 ≤ 0 (tight, satisfied)
    g2 = con.g([1.0, 0.0, 0.0, 0.0], nothing)
    @test g2 ≈ [0.0]

    # Violated: Re(ρ₁₂) = 0.8, ρ₁₁ = 0.5, ρ₂₂ = 0.5
    # Cauchy-Schwarz: 0.8² + 0² - 0.5*0.5 = 0.64 - 0.25 = 0.39 > 0
    g3 = con.g([0.5, 0.8, 0.5, 0.0], nothing)
    @test g3[1] > 0
end

@testitem "BoundStateL2Constraint compact_density 3x3" begin
    using NamedTrajectories
    using DirectTrajOpt

    N = 5
    # 3×3 density → compact iso dim = 9
    traj = NamedTrajectory(
        (ρ⃗̃ = rand(9, N), u = rand(1, N), Δt = fill(0.1, N));
        timestep = :Δt,
        controls = :u,
    )

    con = BoundStateL2Constraint(:ρ⃗̃, traj, :compact_density)
    @test con.equality == false
    # 3×3 has 3 off-diagonal pairs → 3 constraints
    @test con.g_dim == 3

    # Identity: ρ = I/3 → compact = [1/3, 0, 1/3, 0, 0, 1/3, 0, 0, 0]
    # All off-diagonal Re=Im=0, all diag=1/3
    # Cauchy-Schwarz: 0 - (1/3)*(1/3) = -1/9 for each pair
    x = [1/3, 0.0, 1/3, 0.0, 0.0, 1/3, 0.0, 0.0, 0.0]
    g = con.g(x, nothing)
    @test length(g) == 3
    @test all(g .≈ -1/9)
end

end
