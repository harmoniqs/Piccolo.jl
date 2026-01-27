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
