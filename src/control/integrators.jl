module QuantumIntegrators

using LinearAlgebra
using NamedTrajectories
using DirectTrajOpt
using ...Quantum
using ...Quantum:
    SamplingTrajectory,
    MultiKetTrajectory,
    MultiDensityTrajectory,
    state_name,
    state_names,
    sampling_member_states,
    drive_name
using SparseArrays
using TestItems

import DirectTrajOpt: BilinearIntegrator
import DirectTrajOpt: TimeDependentBilinearIntegrator

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
    if sys.time_dependent
        Ĝ = (u_, t) -> I(sys.levels) ⊗ sys.G(u_, t)
        return TimeDependentBilinearIntegrator(
            Ĝ,
            state_name(qtraj),
            drive_name(qtraj),
            :t,
            traj,
        )
    else
        Ĝ = u_ -> I(sys.levels) ⊗ sys.G(u_, 0.0)
        return BilinearIntegrator(Ĝ, state_name(qtraj), drive_name(qtraj), traj)
    end
end

"""
    BilinearIntegrator(qtraj::KetTrajectory, N::Int)

Create a BilinearIntegrator for ket evolution.
"""
function BilinearIntegrator(qtraj::KetTrajectory, N::Int)
    sys = get_system(qtraj)
    traj = NamedTrajectory(qtraj, N)
    if sys.time_dependent
        Ĝ = (u_, t) -> sys.G(u_, t)
        return TimeDependentBilinearIntegrator(
            Ĝ,
            state_name(qtraj),
            drive_name(qtraj),
            :t,
            traj,
        )
    else
        Ĝ = u_ -> sys.G(u_, 0.0)
        return BilinearIntegrator(Ĝ, state_name(qtraj), drive_name(qtraj), traj)
    end
end

"""
    BilinearIntegrator(qtraj::DensityTrajectory, N::Int)

Create a BilinearIntegrator for density matrix evolution using the compact
Lindbladian generators (n² × n²) matching the compact density isomorphism.
"""
function BilinearIntegrator(qtraj::DensityTrajectory, N::Int)
    sys = get_system(qtraj)
    traj = NamedTrajectory(qtraj, N)

    # Build compact generator via the shared helper — uses drive_coeff and
    # rate_coeff at call time so NonlinearDrive coefficients (including
    # global-reading ones) and NonlinearDissipator rates propagate correctly.
    # Pass sparse factors through; Julia's type promotion handles `Dual` u
    # during ForwardDiff jacobian eval and avoids `O(n^4)` densification.
    𝒢c_drift_ham, 𝒢c_drives, 𝒢c_dissipators = compact_lindbladian_parts(sys)
    𝒢c = compact_generator_closure(sys, 𝒢c_drift_ham, 𝒢c_drives, 𝒢c_dissipators)

    return BilinearIntegrator(𝒢c, state_name(qtraj), drive_name(qtraj), traj)
end

"""
    BilinearIntegrator(qtraj::MultiKetTrajectory, N::Int)

Create a vector of BilinearIntegrators for each ket in an MultiKetTrajectory.
"""
function BilinearIntegrator(qtraj::MultiKetTrajectory, N::Int)
    sys = get_system(qtraj)
    traj = NamedTrajectory(qtraj, N)
    control_sym = drive_name(qtraj)
    snames = state_names(qtraj)
    if sys.time_dependent
        Ĝ = (u_, t) -> sys.G(u_, t)
        return [
            TimeDependentBilinearIntegrator(Ĝ, name, control_sym, :t, traj) for
            name in snames
        ]
    else
        Ĝ = u_ -> sys.G(u_, 0.0)
        return [BilinearIntegrator(Ĝ, name, control_sym, traj) for name in snames]
    end
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
    control_sym = drive_name(qtraj)
    systems = qtraj.systems
    member_states = sampling_member_states(qtraj)

    per_member = map(zip(systems, member_states)) do (sys, states)
        _sampling_integrator(qtraj.base_trajectory, sys, traj, states, control_sym)
    end
    # Single-state bases return one integrator per member; multi-state bases
    # return one vector of integrators per member — flatten either way.
    return reduce(vcat, map(i -> i isa AbstractVector ? i : [i], per_member))
end

# Helper to create single integrator for sampling - dispatches on base trajectory type
function _sampling_integrator(
    base_qtraj::UnitaryTrajectory,
    sys::AbstractQuantumSystem,
    traj::NamedTrajectory,
    state_sym::Symbol,
    control_sym::Symbol,
)
    if sys.time_dependent
        Ĝ = (u_, t) -> I(sys.levels) ⊗ sys.G(u_, t)
        return TimeDependentBilinearIntegrator(Ĝ, state_sym, control_sym, :t, traj)
    else
        Ĝ = u_ -> I(sys.levels) ⊗ sys.G(u_, 0.0)
        return BilinearIntegrator(Ĝ, state_sym, control_sym, traj)
    end
end

function _sampling_integrator(
    base_qtraj::KetTrajectory,
    sys::AbstractQuantumSystem,
    traj::NamedTrajectory,
    state_sym::Symbol,
    control_sym::Symbol,
)
    if sys.time_dependent
        Ĝ = (u_, t) -> sys.G(u_, t)
        return TimeDependentBilinearIntegrator(Ĝ, state_sym, control_sym, :t, traj)
    else
        Ĝ = u_ -> sys.G(u_, 0.0)
        return BilinearIntegrator(Ĝ, state_sym, control_sym, traj)
    end
end

function _sampling_integrator(
    base_qtraj::MultiKetTrajectory,
    sys::AbstractQuantumSystem,
    traj::NamedTrajectory,
    state_syms::Vector{Symbol},
    control_sym::Symbol,
)
    # One integrator per ket sub-state of this member, all on the member's system
    if sys.time_dependent
        Ĝ = (u_, t) -> sys.G(u_, t)
        return [
            TimeDependentBilinearIntegrator(Ĝ, name, control_sym, :t, traj) for
            name in state_syms
        ]
    else
        Ĝ = u_ -> sys.G(u_, 0.0)
        return [BilinearIntegrator(Ĝ, name, control_sym, traj) for name in state_syms]
    end
end

function _sampling_integrator(
    base_qtraj::DensityTrajectory,
    sys::OpenQuantumSystem,
    traj::NamedTrajectory,
    state_sym::Symbol,
    control_sym::Symbol,
)
    # Compact Lindbladian (n² × n²) matching the compact density iso the sampling
    # trajectory carries — same construction as the non-sampling density integrator.
    𝒢c_drift_ham, 𝒢c_drives, 𝒢c_dissipators = compact_lindbladian_parts(sys)
    𝒢c = compact_generator_closure(sys, 𝒢c_drift_ham, 𝒢c_drives, 𝒢c_dissipators)
    return BilinearIntegrator(𝒢c, state_sym, control_sym, traj)
end

function _sampling_integrator(
    base_qtraj::MultiDensityTrajectory,
    sys::OpenQuantumSystem,
    traj::NamedTrajectory,
    state_syms::Vector{Symbol},
    control_sym::Symbol,
)
    # One compact-Lindbladian integrator per density sub-state of this member
    𝒢c_drift_ham, 𝒢c_drives, 𝒢c_dissipators = compact_lindbladian_parts(sys)
    𝒢c = compact_generator_closure(sys, 𝒢c_drift_ham, 𝒢c_drives, 𝒢c_dissipators)
    return [BilinearIntegrator(𝒢c, name, control_sym, traj) for name in state_syms]
end

# ----------------------------------------------------------------------------- #
# Variational Integrators
# ----------------------------------------------------------------------------- #

export VariationalKetIntegrator, VariationalUnitaryIntegrator

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

@testitem "VariationalKetIntegrator / VariationalUnitaryIntegrator construct and integrate (DTO ≥ 0.10.1)" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    # Requires the multi-state BilinearIntegrator constructor restored in
    # DirectTrajOpt 0.10.1 (harmoniqs/DirectTrajOpt.jl#139): these constructors
    # pass a stacked names vector vcat(ψ̃, ψ̃_variations...) and had thrown
    # MethodError on every call since the DTO c9fdeb7 refactor (issue #300).
    # Zero callers existed, so nothing fired until the coverage campaign.
    #
    # NOTE: the analytic-Jacobian gates of test_integrator are deliberately
    # NOT run here — the variational construction's ForwardDiff path disagrees
    # with finite differences (duals dropped through var_G / expv). Tracked as
    # #307; pin it there once fixed, not here.
    pkgversion(DirectTrajOpt) >= v"0.10.1" || return  # no-op on the 0.9 line

    H_drift = PAULIS[:Z] / 2
    H_drives = [PAULIS[:X], PAULIS[:Y]]
    H_vars = [PAULIS[:Z] / 2]  # drift-frequency uncertainty direction
    sys = VariationalQuantumSystem(H_drift, H_drives, H_vars, [1.0, 1.0])

    N = 5
    ψ̃₀ = [1.0; 0.0; 0.0; 0.0]  # iso(|0⟩) = vcat(real, imag)

    # Ket: state ψ̃ (4-dim) + one variational copy ψ̃_var (4-dim)
    ψ̃data = zeros(4, N)
    ψ̃data[:, 1] = ψ̃₀
    traj = NamedTrajectory(
        (ψ̃ = ψ̃data, ψ̃_var = zeros(4, N), u = zeros(2, N), Δt = fill(0.1, N));
        controls = (:u, :Δt),
        timestep = :Δt,
        initial = (ψ̃ = ψ̃₀, ψ̃_var = zeros(4)),
    )
    B = VariationalKetIntegrator(sys, traj, :ψ̃, [:ψ̃_var], :u)
    @test B isa BilinearIntegrator
    @test B.x_names == [:ψ̃, :ψ̃_var]
    @test B.x_dim == 8
    @test B.dim == 8 * (N - 1)
    δ = zeros(B.dim)
    evaluate!(δ, B, traj)
    @test norm(δ) > 0          # integrates without MethodError
    @test all(isfinite, δ)

    # The stacked dimension feeds get_nonlinear_constraints through x_names
    prob = DirectTrajOptProblem(traj, QuadraticRegularizer(:u, traj, 1.0), [B])
    @test length(Solvers.get_nonlinear_constraints(prob)) == B.dim

    # Unitary: Ũ⃗ = vcat(real.(vec(U)), imag.(vec(U))) (8-dim for 2 levels)
    # + one variational copy (8-dim). The generator is var_G(I ⊗ G, [I ⊗ G_var]).
    U₀ = Matrix{Float64}(I(2))
    Ũ⃗₀ = vcat(real.(vec(U₀)), imag.(vec(U₀)))
    Ũdata = zeros(8, N)
    Ũdata[:, 1] = Ũ⃗₀
    traj_U = NamedTrajectory(
        (Ũ⃗ = Ũdata, Ũ⃗_var = zeros(8, N), u = zeros(2, N), Δt = fill(0.1, N));
        controls = (:u, :Δt),
        timestep = :Δt,
        initial = (Ũ⃗ = Ũ⃗₀, Ũ⃗_var = zeros(8)),
    )
    B_U = VariationalUnitaryIntegrator(sys, traj_U, :Ũ⃗, [:Ũ⃗_var], :u)
    @test B_U isa BilinearIntegrator
    @test B_U.x_names == [:Ũ⃗, :Ũ⃗_var]
    @test B_U.x_dim == 16
    δ_U = zeros(B_U.dim)
    evaluate!(δ_U, B_U, traj_U)
    @test norm(δ_U) > 0
    @test all(isfinite, δ_U)
end

@testitem "BilinearIntegrator dispatch on UnitaryTrajectory" begin
    using DirectTrajOpt
    using NamedTrajectories

    # Create system and pulse
    sys = OpenQuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
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
    sys = OpenQuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
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
    sys1 = OpenQuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    sys2 = OpenQuantumSystem(1.1 * GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

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
    sys1 = OpenQuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    sys2 = OpenQuantumSystem(1.1 * GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

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

@testitem "BilinearIntegrator dispatch on SamplingTrajectory (Density)" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    # Open systems with dissipation and a drift variation
    L = ComplexF64[0.0 0.1; 0.0 0.0]
    sys1 = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    sys2 =
        OpenQuantumSystem(0.95 * PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])

    ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρg = ComplexF64[0.0 0.0; 0.0 1.0]

    N = 11
    times = collect(range(0, 1.0, length = N))
    pulse = LinearSplinePulse(zeros(1, N), times)

    base_qtraj = DensityTrajectory(sys1, pulse, ρ0, ρg)
    sampling_qtraj = SamplingTrajectory(base_qtraj, [sys1, sys2])

    expanded_traj = NamedTrajectory(sampling_qtraj, N)
    integrators = BilinearIntegrator(sampling_qtraj, N)

    @test integrators isa Vector{<:BilinearIntegrator}
    @test length(integrators) == 2

    # Compact Lindbladian: state dim is n² (compact iso), NOT 2n² (full vec iso)
    n = sys1.levels
    for integrator in integrators
        @test integrator.x_dim == n^2
        test_integrator(integrator, expanded_traj; atol = 1e-3)
    end

    # Parity: the member-1 sampling integrator must agree with the non-sampling
    # density integrator built on the same system (same compact generator).
    ref_integrator = BilinearIntegrator(base_qtraj, N)
    x = randn(n^2)
    x_next = randn(n^2)
    u = [0.3]
    Δt = 0.05
    @test integrators[1].f(x_next, x, u, Δt) ≈ ref_integrator.f(x_next, x, u, Δt)
end

@testitem "BilinearIntegrator dispatch on SamplingTrajectory (MultiKet)" begin
    using DirectTrajOpt
    using NamedTrajectories

    # Systems with drift variation
    sys1 = OpenQuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])
    sys2 = OpenQuantumSystem(1.1 * GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

    ψ0 = ComplexF64[1.0, 0.0]
    ψ1 = ComplexF64[0.0, 1.0]

    N = 11
    times = collect(range(0, 1.0, length = N))
    pulse = LinearSplinePulse(zeros(2, N), times)

    base_qtraj = MultiKetTrajectory(sys1, pulse, [ψ0, ψ1], [ψ1, ψ0])
    sampling_qtraj = SamplingTrajectory(base_qtraj, [sys1, sys2])

    expanded_traj = NamedTrajectory(sampling_qtraj, N)
    integrators = BilinearIntegrator(sampling_qtraj, N)

    # 2 members × 2 kets = 4 integrators, one per (member, ket) state component
    @test integrators isa Vector{<:BilinearIntegrator}
    @test length(integrators) == 4
    @test [i.x_name for i in integrators] == [:ψ̃1, :ψ̃2, :ψ̃3, :ψ̃4]

    for integrator in integrators
        test_integrator(integrator, expanded_traj; atol = 1e-3)
    end
end

@testitem "BilinearIntegrator dispatch on SamplingTrajectory (MultiDensity)" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    L = ComplexF64[0.0 0.1; 0.0 0.0]
    sys1 = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    sys2 =
        OpenQuantumSystem(0.95 * PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])

    ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρ1 = ComplexF64[0.0 0.0; 0.0 1.0]

    N = 11
    times = collect(range(0, 1.0, length = N))
    pulse = LinearSplinePulse(zeros(1, N), times)

    base_qtraj = MultiDensityTrajectory(sys1, pulse, [ρ0, ρ1], [ρ1, ρ0])
    sampling_qtraj = SamplingTrajectory(base_qtraj, [sys1, sys2])

    expanded_traj = NamedTrajectory(sampling_qtraj, N)
    integrators = BilinearIntegrator(sampling_qtraj, N)

    # 2 members × 2 densities = 4 compact-Lindbladian integrators
    @test integrators isa Vector{<:BilinearIntegrator}
    @test length(integrators) == 4
    @test [i.x_name for i in integrators] == [:ρ⃗̃1, :ρ⃗̃2, :ρ⃗̃3, :ρ⃗̃4]

    n = sys1.levels
    for integrator in integrators
        @test integrator.x_dim == n^2  # compact iso, not 2n²
        test_integrator(integrator, expanded_traj; atol = 1e-3)
    end
end

@testitem "BilinearIntegrator dispatch on MultiKetTrajectory" begin
    using DirectTrajOpt
    using NamedTrajectories

    # Shared system
    sys = OpenQuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]], [1.0, 1.0])

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
    sys = OpenQuantumSystem(H, [1.0, 1.0])

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
    sys = OpenQuantumSystem(H, [1.0])

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
    sys = OpenQuantumSystem(H, [1.0, 1.0])

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
    sys1 = OpenQuantumSystem(H1, [1.0])
    sys2 = OpenQuantumSystem(H2, [1.0])

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

@testitem "BilinearIntegrator dispatch on modulated UnitaryTrajectory" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    omega = 2pi * 2.0
    H_z = GATES[:Z]
    H_x = GATES[:X]

    # Build modulated system via Pair
    sys = OpenQuantumSystem(H_z, [H_x => t -> cos(omega * t)], [1.0])
    @test sys.time_dependent

    T = 1.0
    N = 11
    times = collect(range(0, T, length = N))
    controls = zeros(1, N)
    pulse = LinearSplinePulse(controls, times)

    qtraj = UnitaryTrajectory(sys, pulse, GATES[:X])
    integrator = BilinearIntegrator(qtraj, N)

    # Should dispatch to TimeDependentBilinearIntegrator
    @test integrator isa TimeDependentBilinearIntegrator

    traj = NamedTrajectory(qtraj, N)

    # Validate that evaluate! runs and Jacobian has correct size
    δ = zeros(integrator.dim)
    evaluate!(δ, integrator, traj)
    @test !all(iszero.(δ))

    ∂f = eval_jacobian(integrator, traj)
    @test size(∂f, 1) == integrator.dim
    @test size(∂f, 2) == traj.dim * traj.N + traj.global_dim
end

@testitem "BilinearIntegrator dispatch on modulated KetTrajectory" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    omega = 2pi * 2.0
    H_z = GATES[:Z]
    H_x = GATES[:X]

    sys = OpenQuantumSystem(H_z, [H_x => t -> cos(omega * t)], [1.0])

    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]

    T = 1.0
    N = 11
    times = collect(range(0, T, length = N))
    controls = zeros(1, N)
    pulse = LinearSplinePulse(controls, times)

    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)
    integrator = BilinearIntegrator(qtraj, N)

    @test integrator isa TimeDependentBilinearIntegrator

    traj = NamedTrajectory(qtraj, N)

    # Validate that evaluate! runs and Jacobian has correct size
    δ = zeros(integrator.dim)
    evaluate!(δ, integrator, traj)
    @test !all(iszero.(δ))

    ∂f = eval_jacobian(integrator, traj)
    @test size(∂f, 1) == integrator.dim
    @test size(∂f, 2) == traj.dim * traj.N + traj.global_dim
end

@testitem "BilinearIntegrator dispatch on modulated SamplingTrajectory (Unitary)" begin
    using DirectTrajOpt
    using NamedTrajectories
    using LinearAlgebra

    omega = 2pi * 2.0
    H_z = GATES[:Z]
    H_x = GATES[:X]

    sys1 = OpenQuantumSystem(H_z, [H_x => t -> cos(omega * t)], [1.0])
    sys2 = OpenQuantumSystem(1.1 * H_z, [H_x => t -> cos(omega * t)], [1.0])

    T = 1.0
    N = 11
    times = collect(range(0, T, length = N))
    controls = zeros(1, N)
    pulse = LinearSplinePulse(controls, times)

    base_qtraj = UnitaryTrajectory(sys1, pulse, GATES[:X])
    sampling_qtraj = SamplingTrajectory(base_qtraj, [sys1, sys2])

    traj = NamedTrajectory(sampling_qtraj, N)
    integrators = BilinearIntegrator(sampling_qtraj, N)

    @test integrators isa Vector
    @test length(integrators) == 2

    # Validate each integrator: evaluate! and Jacobian dimensions
    for integrator in integrators
        δ = zeros(integrator.dim)
        evaluate!(δ, integrator, traj)
        @test !all(iszero.(δ))

        ∂f = eval_jacobian(integrator, traj)
        @test size(∂f, 1) == integrator.dim
        @test size(∂f, 2) == traj.dim * traj.N + traj.global_dim
    end
end

# ----------------------------------------------------------------------------- #
# Exponential integrator family (moved from Piccolissimo — open-core slice 3a,
# harmoniqs/Piccolissimo.jl#429; see docs/moved-file-review.md rows 1-9 + the
# two sampling dispatch files that postdate the survey). The matrix-free
# MultiKet hook (`matrix_free_jacobian_op`) stays proprietary: the function is
# declared empty here, Piccolissimo's solvers define its methods.
# ----------------------------------------------------------------------------- #

include("integrators_exponential/_exponential_integrators.jl")

# ----------------------------------------------------------------------------- #
# SplineIntegrator struct + dense cells + shared interval-coefficient kernel
# (open-core slice 3b, harmoniqs/Piccolissimo.jl#430). The matrix-free
# cells/kernels, MF layout machinery and GPU variant stay in Piccolissimo as
# method extensions on the shared types below.
# ----------------------------------------------------------------------------- #

include("integrators_spline/_spline_integrators.jl")
using .SplineIntegrators
export AbstractSplineIntegrator
export SplineIntegrator
export du_name, ddu_name
export build_sensitivity_ode, build_sensitivity_problems, extract_sensitivity_solution!
export PropagatorResult, get_propagator, get_sensitivities, get_sensitivities_flat
export IntegrationAlgorithm, Tsit5Alg, MagnusGL4Alg, MagnusAdapt4Alg
export Tsit5Data, MagnusGL4Data, MagnusAdapt4Data
export ChebyshevAlg
export refresh_sensitivities!
export SplineType, LinearSpline, CubicSpline
export SplineIntervalCoeffs, interval_coeff!, interval_coeff_dir!, interval_vjp_scatter!
export interval_hvp_scatter!
export LindbladDuhamelTape, DensityLindbladData, compact_iso_hs_weights
# The submodule's second export block (Piccolo #326): rollout helpers + the
# sensitivity-kick family — moved in 3b but missed by the seam, leaving
# downstream `using Piccolo` callers (Piccolissimo's benchmark/correctness
# suite) with UndefVarError on names the package ships.
export gauss_legendre_01,
    compute_sensitivity_kick_first_order,
    compute_sensitivity_kick_exact,
    KnotPointPropagationData,
    setup_knot_point_propagation,
    ket_vjp,
    ket_jvp,
    unitary_rollout_trajectory
export ChebyshevData,
    Rodas5PAlg,
    Rodas5PData,
    SplineIntegrators,
    build_multiket_hvp_cache,
    canonical_hessian_knot_dim,
    expl_discretization_error,
    suggest_n_sub
using .ExponentialIntegrators
export AbstractExponentialIntegrator
export HermitianExponentialIntegrator
export NonHermitianExponentialIntegrator
export x_name, single_state_dim
export exp_eigen, exp_eigen!, exp_generator!
export dk_divided_difference!, dk_apply!, dk_first_order_derivative!
export DaleckiiKreinWorkspace, DK_DEGENERACY_RTOL
export dk_second_divided_difference,
    dk_second_order_apply!, dk_second_order_derivative!, dk_second_order_block!
export DaleckiiKreinSecondOrderWorkspace

end
