# ============================================================================ #
# Extract Pulse from Optimized Controls
# ============================================================================ #

export extract_pulse

"""
    extract_pulse(qtraj::AbstractQuantumTrajectory, traj::NamedTrajectory)

Extract an optimized pulse from a NamedTrajectory.

This function extracts the control values from the optimized trajectory and creates
a new pulse object of the same type as the original pulse in `qtraj`.

The extraction process depends on the pulse type:
- `ZeroOrderPulse`, `LinearSplinePulse`: Extracts `u` (drive variable)
- `CubicSplinePulse`: Extracts both `u` and `du` (derivative variable)

# Arguments
- `qtraj`: Original quantum trajectory (provides pulse type and drive names)
- `traj`: Optimized NamedTrajectory with new control values

# Returns
A new pulse of the same type as `qtraj.pulse` with optimized control values.

# Example
```julia
# After optimization
solve!(prob)
new_pulse = extract_pulse(qtraj, prob.trajectory)
rollout!(qtraj, new_pulse)
```
"""
function extract_pulse end

# Dispatch on pulse type
function extract_pulse(
    qtraj::AbstractQuantumTrajectory{<:Union{ZeroOrderPulse,LinearSplinePulse}},
    traj::NamedTrajectory,
)
    times = collect(get_times(traj))
    u_name = drive_name(qtraj)
    u = Matrix(traj[u_name])
    return _rebuild_pulse(qtraj.pulse, u, times)
end

function extract_pulse(
    qtraj::AbstractQuantumTrajectory{<:CubicSplinePulse},
    traj::NamedTrajectory,
)
    times = collect(get_times(traj))
    u_name = drive_name(qtraj)
    du_name = Symbol(:d, u_name)
    u = Matrix(traj[u_name])
    du = Matrix(traj[du_name])
    return CubicSplinePulse(u, du, times; drive_name = u_name)
end

# MVP: B-spline trajectories were routed through the Hermite (u, du) representation
# at construction time (see `_get_control_data(::BSplinePulse, ...)`). After
# optimization the trajectory holds Hermite samples at the integration nodes.
# Recover a BSplinePulse via least-squares fit of control points to those samples
# (M unknowns vs N = M - k + 2 equations → underdetermined for k > 2; use pinv).
# The full Route-A path that has control points as decision variables directly
# is tracked in `amico/vault/specs/spec-20260527-173644-bspline-pulse.md`.
function extract_pulse(
    qtraj::AbstractQuantumTrajectory{<:BSplinePulse},
    traj::NamedTrajectory,
)
    times = collect(get_times(traj))
    u_name = drive_name(qtraj)
    du_name = Symbol(:d, u_name)
    u = Matrix(traj[u_name])
    du = haskey(traj.components, du_name) ? Matrix(traj[du_name]) : zeros(size(u))

    pulse_orig = qtraj.pulse
    new_cp = _fit_bspline_to_hermite(pulse_orig, u, du, times)
    τ = pulse_orig.basis.knot_vector
    return BSplinePulse(
        new_cp,
        [τ[1], τ[end]];
        order = get_order(pulse_orig),
        drive_name = u_name,
        initial_value = :free,  # already constrained by clamped basis; don't overwrite
        final_value = :free,
    )
end

# Fit M B-spline control points to the cubic Hermite curve implied by (u, du)
# samples by oversampling the Hermite curve at `Ndense` points and solving the
# overdetermined least-squares B * c = u_dense.
#
# This is the MVP path: the actual optimization happens in Hermite space, the
# B-spline acts as a smooth low-pass approximation of the optimized result so
# the user sees a recognizable smooth curve. The full Route-A path (control
# points as native decision variables) lives in
# `amico/vault/specs/spec-20260527-173644-bspline-pulse.md`.
function _fit_bspline_to_hermite(
    pulse::BSplinePulse{Order},
    u_values::AbstractMatrix,
    du_values::AbstractMatrix,
    times::AbstractVector;
    oversample::Int = 20,
) where {Order}
    M = pulse.basis.M
    n_d = size(u_values, 1)
    τ = pulse.basis.knot_vector

    # Build a Hermite curve from (u, du, times) and sample it densely
    hermite = CubicSplinePulse(u_values, du_values, collect(times))
    t_start, t_end = first(times), last(times)
    Ndense = oversample * length(times)
    t_dense = collect(range(t_start, t_end, length = Ndense))
    u_dense = hcat([hermite(t) for t in t_dense]...)

    # Basis matrix at dense times: B[i, j] = j-th basis at t_dense[i]
    B = zeros(Ndense, M)
    cp_indicator = zeros(1, M)
    for j = 1:M
        cp_indicator[1, j] = 1.0
        for (i, t) in enumerate(t_dense)
            t_clamped = clamp(t, τ[1], τ[end])
            B[i, j] = _deboor_basis_eval(cp_indicator, τ, Order, t_clamped)
        end
        cp_indicator[1, j] = 0.0
    end

    # Overdetermined least-squares: Ndense >> M
    C = zeros(n_d, M)
    for d = 1:n_d
        C[d, :] = B \ Vector(u_dense[d, :])
    end
    return C
end

# Wrapper that evaluates a single-drive B-spline with a given indicator control-point
# vector at a clamped t. Type-generic so ForwardDiff can differentiate through `t`.
function _deboor_basis_eval(
    cp_indicator::AbstractMatrix,
    knot_vector::AbstractVector,
    order::Int,
    t::Real,
)
    M = size(cp_indicator, 2)
    k = order
    Tout = promote_type(eltype(cp_indicator), eltype(knot_vector), typeof(t))
    j = clamp(searchsortedlast(knot_vector, t), k, M)
    p_work = Vector{Tout}(undef, k)
    @inbounds for r = 0:(k - 1)
        p_work[r + 1] = cp_indicator[1, j - r]
    end
    @inbounds for level = 1:(k - 1)
        for r = 0:(k - 1 - level)
            knot_low = knot_vector[j - r]
            knot_high = knot_vector[j + k - level - r]
            denom = knot_high - knot_low
            α = denom > zero(denom) ? Tout((t - knot_low) / denom) : zero(Tout)
            p_work[r + 1] = (one(Tout) - α) * p_work[r + 2] + α * p_work[r + 1]
        end
    end
    return p_work[1]
end

# SamplingTrajectory delegates to base_trajectory
function extract_pulse(qtraj::SamplingTrajectory, traj::NamedTrajectory)
    return extract_pulse(qtraj.base_trajectory, traj)
end

# Helper functions for pulse reconstruction
function _rebuild_pulse(p::ZeroOrderPulse, u::Matrix, times::Vector)
    return ZeroOrderPulse(u, times; drive_name = p.drive_name)
end

function _rebuild_pulse(p::LinearSplinePulse, u::Matrix, times::Vector)
    return LinearSplinePulse(u, times; drive_name = p.drive_name)
end

# ============================================================================ #
# Tests
# ============================================================================ #

@testitem "extract_pulse with ZeroOrderPulse - UnitaryTrajectory" begin
    using LinearAlgebra
    using NamedTrajectories

    system = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])

    T = 1.0
    X_gate = ComplexF64[0 1; 1 0]
    qtraj = UnitaryTrajectory(system, X_gate, T)

    N = 11
    traj = NamedTrajectory(qtraj, N)
    new_controls = 0.5 * randn(1, N)
    traj.u .= new_controls

    pulse = extract_pulse(qtraj, traj)

    @test pulse isa ZeroOrderPulse
    # Sample at the trajectory times to verify
    sampled = sample(pulse, collect(get_times(traj)))
    @test sampled ≈ new_controls
    @test pulse.drive_name == :u
    @test duration(pulse) ≈ T
end

@testitem "extract_pulse with LinearSplinePulse - KetTrajectory" begin
    using LinearAlgebra
    using NamedTrajectories

    system = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])

    T = 1.0
    times = collect(range(0.0, T, length = 11))
    u = 0.1 * randn(1, 11)
    pulse = LinearSplinePulse(u, times)

    ψ0 = ComplexF64[1.0, 0.0]
    ψg = ComplexF64[0.0, 1.0]
    qtraj = KetTrajectory(system, pulse, ψ0, ψg)

    traj = NamedTrajectory(qtraj, times)
    new_controls = 0.5 * randn(1, 11)
    traj.u .= new_controls

    pulse_new = extract_pulse(qtraj, traj)

    @test pulse_new isa LinearSplinePulse
    sampled = sample(pulse_new, times)
    @test sampled ≈ new_controls
    @test pulse_new.drive_name == :u
end

@testitem "extract_pulse with CubicSplinePulse - UnitaryTrajectory" begin
    using LinearAlgebra
    using NamedTrajectories

    system = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])

    T = 1.0
    times = collect(range(0.0, T, length = 11))
    u = 0.1 * randn(1, 11)
    du = zeros(1, 11)
    pulse = CubicSplinePulse(u, du, times)

    X_gate = ComplexF64[0 1; 1 0]
    qtraj = UnitaryTrajectory(system, pulse, X_gate)

    traj = NamedTrajectory(qtraj, times)
    new_u = 0.5 * randn(1, 11)
    new_du = 0.1 * randn(1, 11)
    traj.u .= new_u
    traj.du .= new_du

    pulse_new = extract_pulse(qtraj, traj)

    @test pulse_new isa CubicSplinePulse
    # Sample at trajectory times to verify controls were extracted
    sampled = sample(pulse_new, times)
    @test sampled ≈ new_u
    @test pulse_new.drive_name == :u
end

@testitem "extract_pulse with MultiKetTrajectory" begin
    using LinearAlgebra
    using NamedTrajectories

    system = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])

    T = 1.0
    initials = [ComplexF64[1.0, 0.0], ComplexF64[0.0, 1.0]]
    goals = [ComplexF64[0.0, 1.0], ComplexF64[1.0, 0.0]]
    weights = [0.6, 0.4]

    qtraj = MultiKetTrajectory(system, initials, goals, T; weights = weights)

    N = 11
    traj = NamedTrajectory(qtraj, N)
    new_controls = 0.5 * randn(1, N)
    traj.u .= new_controls

    pulse = extract_pulse(qtraj, traj)

    @test pulse isa ZeroOrderPulse
    sampled = sample(pulse, collect(get_times(traj)))
    @test sampled ≈ new_controls
end

@testitem "extract_pulse with SamplingTrajectory" begin
    using LinearAlgebra
    using NamedTrajectories

    system = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])

    T = 1.0
    X_gate = ComplexF64[0 1; 1 0]
    base_qtraj = UnitaryTrajectory(system, X_gate, T)

    # Create sampling trajectory with random samples
    samples = [system for _ = 1:3]
    qtraj = SamplingTrajectory(base_qtraj, samples)

    N = 11
    traj = NamedTrajectory(qtraj, N)
    new_controls = 0.5 * randn(1, N)
    traj.u .= new_controls

    pulse = extract_pulse(qtraj, traj)

    @test pulse isa ZeroOrderPulse
    sampled = sample(pulse, collect(get_times(traj)))
    @test sampled ≈ new_controls
end

@testitem "extract_pulse preserves drive_name" begin
    using LinearAlgebra
    using NamedTrajectories

    system = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])

    T = 1.0
    times = collect(range(0.0, T, length = 11))
    u = 0.1 * randn(1, 11)
    pulse = ZeroOrderPulse(u, times; drive_name = :a)

    X_gate = ComplexF64[0 1; 1 0]
    qtraj = UnitaryTrajectory(system, pulse, X_gate)

    traj = NamedTrajectory(qtraj, times)
    traj.a .= 0.5 * randn(1, 11)

    pulse_new = extract_pulse(qtraj, traj)

    @test pulse_new.drive_name == :a
    sampled = sample(pulse_new, times)
    @test sampled ≈ traj.a
end

@testitem "extract_pulse with multi-drive system" begin
    using LinearAlgebra
    using NamedTrajectories

    system = QuantumSystem(PAULIS.Z, [PAULIS.X, PAULIS.Y], [1.0, 1.0])

    T = 1.0
    X_gate = ComplexF64[0 1; 1 0]
    qtraj = UnitaryTrajectory(system, X_gate, T)

    N = 11
    traj = NamedTrajectory(qtraj, N)
    new_controls = randn(2, N)
    traj.u .= new_controls

    pulse = extract_pulse(qtraj, traj)

    @test pulse isa ZeroOrderPulse
    @test pulse.n_drives == 2
    sampled = sample(pulse, collect(get_times(traj)))
    @test sampled ≈ new_controls
end
