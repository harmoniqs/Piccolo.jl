module QuantumIntegrators

using LinearAlgebra
using NamedTrajectories
using DirectTrajOpt
using ...Quantum
using ...Quantum:
    SamplingTrajectory, MultiKetTrajectory, state_name, state_names, drive_name
using SparseArrays
using TestItems

import DirectTrajOpt: BilinearIntegrator

# Import QuantumTrajectories types (will be loaded before this module)
using ...Quantum.QuantumTrajectories

const ⊗ = kron

# ----------------------------------------------------------------------------- #
# Default Integrators (require NamedTrajectory)
# ----------------------------------------------------------------------------- #

"""
    BilinearIntegrator(qtraj::UnitaryTrajectory, N::Int)

Create a BilinearIntegrator for unitary evolution.
"""
function BilinearIntegrator(qtraj::UnitaryTrajectory, N::Int)
    sys = get_system(qtraj)
    traj = NamedTrajectory(qtraj, N)
    Ĝ = u_ -> I(sys.levels) ⊗ sys.G(u_, 0.0)
    return BilinearIntegrator(Ĝ, state_name(qtraj), drive_name(qtraj), traj)
end

"""
    BilinearIntegrator(qtraj::KetTrajectory, N::Int)

Create a BilinearIntegrator for ket evolution.
"""
function BilinearIntegrator(qtraj::KetTrajectory, N::Int)
    sys = get_system(qtraj)
    traj = NamedTrajectory(qtraj, N)
    Ĝ = u_ -> sys.G(u_, 0.0)
    return BilinearIntegrator(Ĝ, state_name(qtraj), drive_name(qtraj), traj)
end

"""
    BilinearIntegrator(qtraj::DensityTrajectory, N::Int)

Create a BilinearIntegrator for density matrix evolution using the compact
Lindbladian generators (n² × n²) matching the compact density isomorphism.
"""
function BilinearIntegrator(qtraj::DensityTrajectory, N::Int)
    sys = get_system(qtraj)
    traj = NamedTrajectory(qtraj, N)

    # Build compact generator function: u -> 𝒢c_drift + Σ uᵢ 𝒢c_drives[i]
    𝒢c_drift, 𝒢c_drives = compact_lindbladian_generators(sys)
    if isempty(𝒢c_drives)
        𝒢c = u -> 𝒢c_drift
    else
        𝒢c = u -> 𝒢c_drift + sum(u .* 𝒢c_drives)
    end

    return BilinearIntegrator(𝒢c, state_name(qtraj), drive_name(qtraj), traj)
end

"""
    BilinearIntegrator(qtraj::MultiKetTrajectory, N::Int)

Create a vector of BilinearIntegrators for each ket in an MultiKetTrajectory.
"""
function BilinearIntegrator(qtraj::MultiKetTrajectory, N::Int)
    sys = get_system(qtraj)
    traj = NamedTrajectory(qtraj, N)
    Ĝ = u_ -> sys.G(u_, 0.0)
    control_sym = drive_name(qtraj)
    snames = state_names(qtraj)

    return [BilinearIntegrator(Ĝ, name, control_sym, traj) for name in snames]
end

# ----------------------------------------------------------------------------- #
# SamplingTrajectory Integrators
# ----------------------------------------------------------------------------- #

"""
    BilinearIntegrator(qtraj::SamplingTrajectory, N::Int)

Create a vector of BilinearIntegrators for each system in a SamplingTrajectory.

Each system in the sampling ensemble gets its own dynamics integrator, but they
all share the same control variables.

# Returns
- `Vector{BilinearIntegrator}`: One integrator per system in the ensemble
"""
function BilinearIntegrator(qtraj::SamplingTrajectory, N::Int)
    traj = NamedTrajectory(qtraj, N)
    snames = state_names(qtraj)
    control_sym = drive_name(qtraj)
    systems = qtraj.systems

    return [
        _sampling_integrator(qtraj.base_trajectory, sys, traj, name, control_sym) for
        (sys, name) in zip(systems, snames)
    ]
end

# Helper to create single integrator for sampling - dispatches on base trajectory type
function _sampling_integrator(
    base_qtraj::UnitaryTrajectory,
    sys::AbstractQuantumSystem,
    traj::NamedTrajectory,
    state_sym::Symbol,
    control_sym::Symbol,
)
    Ĝ = u_ -> I(sys.levels) ⊗ sys.G(u_, 0.0)
    return BilinearIntegrator(Ĝ, state_sym, control_sym, traj)
end

function _sampling_integrator(
    base_qtraj::KetTrajectory,
    sys::AbstractQuantumSystem,
    traj::NamedTrajectory,
    state_sym::Symbol,
    control_sym::Symbol,
)
    Ĝ = u_ -> sys.G(u_, 0.0)
    return BilinearIntegrator(Ĝ, state_sym, control_sym, traj)
end

function _sampling_integrator(
    base_qtraj::DensityTrajectory,
    sys::OpenQuantumSystem,
    traj::NamedTrajectory,
    state_sym::Symbol,
    control_sym::Symbol,
)
    return BilinearIntegrator(sys.𝒢, state_sym, control_sym, traj)
end

# ----------------------------------------------------------------------------- #
# Variational Integrators
# ----------------------------------------------------------------------------- #

function VariationalKetIntegrator(
    sys::VariationalQuantumSystem,
    traj::NamedTrajectory,
    ψ̃::Symbol,
    ψ̃_variations::AbstractVector{Symbol},
    u::Symbol;
    scale::Float64 = 1.0,
)
    var_ψ̃ = vcat(ψ̃, ψ̃_variations...)
    G = u_ -> Isomorphisms.var_G(sys.G(u_), [G(u_) / scale for G in sys.G_vars])
    return BilinearIntegrator(G, var_ψ̃, u, traj)
end

function VariationalUnitaryIntegrator(
    sys::VariationalQuantumSystem,
    traj::NamedTrajectory,
    Ũ⃗::Symbol,
    Ũ⃗_variations::AbstractVector{Symbol},
    u::Symbol;
    scales::AbstractVector{<:Float64} = fill(1.0, length(sys.G_vars)),
)
    var_Ũ⃗ = vcat(Ũ⃗, Ũ⃗_variations...)

    function Ĝ(u_)
        G0 = sys.G(u_)
        Gs = typeof(G0)[
            I(sys.levels) ⊗ G(u_) / scale for (scale, G) in zip(scales, sys.G_vars)
        ]
        return Isomorphisms.var_G(I(sys.levels) ⊗ G0, Gs)
    end
    return BilinearIntegrator(Ĝ, var_Ũ⃗, u, traj)
end

# ----------------------------------------------------------------------------- #
# Tests
# ----------------------------------------------------------------------------- #

@testitem "BilinearIntegrator dispatch on UnitaryTrajectory" begin
    using DirectTrajOpt
    using NamedTrajectories

    # Create system and pulse
    sys = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    N = 11
    times = collect(range(0, 1.0, length = N))
    controls = zeros(2, N)
    pulse = LinearSplinePulse(controls, times)

    # Create quantum trajectory
    qtraj = UnitaryTrajectory(sys, pulse, GATES[:H])
    traj = NamedTrajectory(qtraj, N)

    integrator = BilinearIntegrator(qtraj, N)

    @test integrator isa BilinearIntegrator
    test_integrator(integrator, traj; atol = 1e-3)
end

@testitem "BilinearIntegrator dispatch on KetTrajectory" begin
    using DirectTrajOpt
    using NamedTrajectories

    # Create system and pulse
    sys = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    N = 11
    times = collect(range(0, 1.0, length = N))
    controls = zeros(2, N)
    pulse = LinearSplinePulse(controls, times)

    # Create quantum trajectory
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)
    traj = NamedTrajectory(qtraj, N)

    integrator = BilinearIntegrator(qtraj, N)

    @test integrator isa BilinearIntegrator
    test_integrator(integrator, traj; atol = 1e-3)
end

@testitem "BilinearIntegrator dispatch on DensityTrajectory" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    # Create open system with dissipation (σ₋ decay operator)
    L = ComplexF64[0.0 0.1; 0.0 0.0]
    sys = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])

    ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρg = ComplexF64[0.0 0.0; 0.0 1.0]

    N = 11
    times = collect(range(0, 1.0, length = N))
    controls = zeros(1, N)
    pulse = ZeroOrderPulse(controls, times)

    qtraj = DensityTrajectory(sys, pulse, ρ0, ρg)
    traj = NamedTrajectory(qtraj, N)

    integrator = BilinearIntegrator(qtraj, N)

    @test integrator isa BilinearIntegrator

    # State dimension should be n² (compact iso)
    n = sys.levels
    @test integrator.x_dim == n^2

    test_integrator(integrator, traj; atol = 1e-3)
end

@testitem "BilinearIntegrator dispatch on SamplingTrajectory (Unitary)" begin
    using DirectTrajOpt
    using NamedTrajectories

    # Create systems with parameter variation
    sys1 = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    sys2 = QuantumSystem(1.1 * GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

    # Create pulse
    N = 11
    times = collect(range(0, 1.0, length = N))
    controls = zeros(2, N)
    pulse = LinearSplinePulse(controls, times)

    # Create base trajectory and sampling trajectory
    base_qtraj = UnitaryTrajectory(sys1, pulse, GATES[:H])
    sampling_qtraj = SamplingTrajectory(base_qtraj, [sys1, sys2])

    # Convert to NamedTrajectory
    expanded_traj = NamedTrajectory(sampling_qtraj, N)

    # Create integrators
    integrators = BilinearIntegrator(sampling_qtraj, N)

    @test integrators isa Vector{<:BilinearIntegrator}
    @test length(integrators) == 2

    # Test each integrator
    for integrator in integrators
        test_integrator(integrator, expanded_traj; atol = 1e-3)
    end
end

@testitem "BilinearIntegrator dispatch on SamplingTrajectory (Ket)" begin
    using DirectTrajOpt
    using NamedTrajectories

    # Create systems with parameter variation
    sys1 = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    sys2 = QuantumSystem(1.1 * GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]

    # Create pulse
    N = 11
    times = collect(range(0, 1.0, length = N))
    controls = zeros(2, N)
    pulse = LinearSplinePulse(controls, times)

    # Create base trajectory and sampling trajectory
    base_qtraj = KetTrajectory(sys1, pulse, ψ_init, ψ_goal)
    sampling_qtraj = SamplingTrajectory(base_qtraj, [sys1, sys2])

    # Convert to NamedTrajectory
    expanded_traj = NamedTrajectory(sampling_qtraj, N)

    # Create integrators
    integrators = BilinearIntegrator(sampling_qtraj, N)

    @test integrators isa Vector{<:BilinearIntegrator}
    @test length(integrators) == 2

    for integrator in integrators
        test_integrator(integrator, expanded_traj; atol = 1e-3)
    end
end

@testitem "BilinearIntegrator dispatch on MultiKetTrajectory" begin
    using DirectTrajOpt
    using NamedTrajectories

    # Shared system
    sys = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

    # Different initial/goal states
    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]

    # Create pulse
    N = 11
    times = collect(range(0, 1.0, length = N))
    controls = zeros(2, N)
    pulse = LinearSplinePulse(controls, times)

    # Create ensemble trajectory: |0⟩ → |1⟩ and |1⟩ → |0⟩
    qtraj = MultiKetTrajectory(sys, pulse, [ψ0, ψ1], [ψ1, ψ0])
    traj = NamedTrajectory(qtraj, N)

    # Create integrators
    integrators = BilinearIntegrator(qtraj, N)

    @test integrators isa Vector{<:BilinearIntegrator}
    @test length(integrators) == 2

    for integrator in integrators
        test_integrator(integrator, traj; atol = 1e-3)
    end
end

@testitem "BilinearIntegrator dispatch on time-dependent UnitaryTrajectory" begin
    using DirectTrajOpt
    using NamedTrajectories

    # Time-dependent Hamiltonian
    ω = 2π * 5.0
    H(u, t) = GATES[:Z] + u[1] * cos(ω * t) * GATES[:X] + u[2] * sin(ω * t) * GATES[:Y]

    T = 1.0
    N = 11
    sys = QuantumSystem(H, [1.0, 1.0])

    times = collect(range(0, T, length = N))
    controls = zeros(2, N)
    pulse = LinearSplinePulse(controls, times)

    qtraj = UnitaryTrajectory(sys, pulse, GATES[:X])
    traj = NamedTrajectory(qtraj, N)

    integrator = BilinearIntegrator(qtraj, N)

    @test integrator isa BilinearIntegrator

    # Test integrator derivatives
    test_integrator(integrator, traj; atol = 1e-2)
end

@testitem "BilinearIntegrator dispatch on time-dependent KetTrajectory" begin
    using DirectTrajOpt
    using NamedTrajectories

    # Time-dependent Hamiltonian
    ω = 2π * 5.0
    H(u, t) = GATES[:Z] + u[1] * cos(ω * t) * GATES[:X]

    T = 1.0
    N = 11
    sys = QuantumSystem(H, [1.0])

    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]

    times = collect(range(0, T, length = N))
    controls = zeros(1, N)
    pulse = LinearSplinePulse(controls, times)

    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)
    traj = NamedTrajectory(qtraj, N)

    integrator = BilinearIntegrator(qtraj, N)

    @test integrator isa BilinearIntegrator

    # Test integrator derivatives
    test_integrator(integrator, traj; atol = 1e-2)
end

@testitem "BilinearIntegrator dispatch on time-dependent MultiKetTrajectory" begin
    using DirectTrajOpt
    using NamedTrajectories

    # Time-dependent Hamiltonian
    ω = 2π * 5.0
    H(u, t) = GATES[:Z] + u[1] * cos(ω * t) * GATES[:X] + u[2] * sin(ω * t) * GATES[:Y]

    T = 1.0
    N = 11
    sys = QuantumSystem(H, [1.0, 1.0])

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]

    times = collect(range(0, T, length = N))
    controls = zeros(2, N)
    pulse = LinearSplinePulse(controls, times)

    qtraj = MultiKetTrajectory(sys, pulse, [ψ0, ψ1], [ψ1, ψ0])
    traj = NamedTrajectory(qtraj, N)

    integrators = BilinearIntegrator(qtraj, N)

    @test integrators isa Vector{<:BilinearIntegrator}
    @test length(integrators) == 2

    for integrator in integrators
        test_integrator(integrator, traj; atol = 1e-2)
    end
end

@testitem "BilinearIntegrator dispatch on time-dependent SamplingTrajectory (Unitary)" begin
    using DirectTrajOpt
    using NamedTrajectories

    # Time-dependent Hamiltonians with parameter variation
    ω = 2π * 5.0
    H1(u, t) = GATES[:Z] + u[1] * cos(ω * t) * GATES[:X]
    H2(u, t) = 1.1 * GATES[:Z] + u[1] * cos(ω * t) * GATES[:X]

    T = 1.0
    N = 11
    sys1 = QuantumSystem(H1, [1.0])
    sys2 = QuantumSystem(H2, [1.0])

    times = collect(range(0, T, length = N))
    controls = zeros(1, N)
    pulse = LinearSplinePulse(controls, times)

    base_qtraj = UnitaryTrajectory(sys1, pulse, GATES[:X])
    sampling_qtraj = SamplingTrajectory(base_qtraj, [sys1, sys2])

    traj = NamedTrajectory(sampling_qtraj, N)

    integrators = BilinearIntegrator(sampling_qtraj, N)

    @test integrators isa Vector{<:BilinearIntegrator}
    @test length(integrators) == 2

    for integrator in integrators
        test_integrator(integrator, traj; atol = 1e-2)
    end
end


end
