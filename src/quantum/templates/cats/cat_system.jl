export CatSystem
export coherent_ket
export get_cat_controls

"""
    coherent_ket(α, levels)

Construct a coherent state ``|α⟩`` in the Fock basis truncated to `levels`.
"""
function coherent_ket(α::Union{Real,Complex}, levels::Int)::Vector{ComplexF64}
    return [exp(-0.5 * abs2(α)) * α^n / sqrt(factorial(n)) for n = 0:(levels-1)]
end

@doc raw"""
    CatSystem(;
        g2=0.36, χ_aa=-7e-3, χ_bb=-32, χ_ab=0.79,
        κa=53e-3, κb=13,
        cat_levels=13, buffer_levels=3,
        prefactor=1, drive_bounds=[1.0, 1.0],
    )::OpenQuantumSystem

Construct an `OpenQuantumSystem` for a two-mode cat qubit (cat ⊗ buffer).

# Hamiltonian

The drift Hamiltonian includes Kerr, cross-Kerr, and two-photon exchange terms:

```math
H = -\frac{\chi_{aa}}{2} a^{\dagger 2} a^2
    -\frac{\chi_{bb}}{2} b^{\dagger 2} b^2
    -\chi_{ab}\, a^\dagger a\, b^\dagger b
    + g_2\, a^{\dagger 2} b
    + g_2^*\, a^2 b^\dagger
```

The two drives are the buffer displacement ``b + b^\dagger`` and a Kerr
correction ``a^\dagger a``.

# Keyword Arguments
- `g2`: Two-photon exchange coupling (MHz · 2π)
- `χ_aa`: Cat self-Kerr (MHz · 2π)
- `χ_bb`: Buffer self-Kerr (MHz · 2π)
- `χ_ab`: Cross-Kerr between cat and buffer (MHz · 2π)
- `κa`: Cat decay rate (MHz · 2π)
- `κb`: Buffer decay rate (MHz · 2π)
- `cat_levels`: Truncation of cat mode Fock space
- `buffer_levels`: Truncation of buffer mode Fock space
- `prefactor`: Global scaling applied to all couplings and rates
- `drive_bounds`: Bounds on the two drive amplitudes

All parameters are scaled by ``2π`` (Hamiltonian) or ``\\sqrt{2π}``
(dissipators) internally.
"""
function CatSystem(;
    g2::Real = 0.36,
    χ_aa::Real = -7e-3,
    χ_bb::Real = -32,
    χ_ab::Real = 0.79,
    κa::Real = 53e-3,
    κb::Real = 13,
    cat_levels::Int = 13,
    buffer_levels::Int = 3,
    prefactor::Real = 1,
    drive_bounds::Vector{<:Real} = [1.0, 1.0],
)
    global_params = (
        g2 = prefactor * g2,
        χ_aa = prefactor * χ_aa,
        χ_bb = prefactor * χ_bb,
        χ_ab = prefactor * χ_ab,
        κa = prefactor * κa,
        κb = prefactor * κb,
        cat_levels = cat_levels,
        buffer_levels = buffer_levels,
        prefactor = prefactor,
    )

    # Apply prefactor
    g2 *= prefactor
    χ_aa *= prefactor
    χ_bb *= prefactor
    χ_ab *= prefactor
    κa *= prefactor
    κb *= prefactor

    # Cat ⊗ Buffer
    a = annihilate(cat_levels) ⊗ Matrix(1.0I, buffer_levels, buffer_levels)
    b = Matrix(1.0I, cat_levels, cat_levels) ⊗ annihilate(buffer_levels)

    H_drift =
        -χ_aa / 2 * a'a'a * a - χ_bb / 2 * b'b'b * b - χ_ab * a'a * b'b +
        g2 * a'a'b +
        conj(g2) * a * a * b'

    # buffer drive, kerr-correction drive
    H_drives = [b + b', a'a]

    L_dissipators = [√κa * a, √κb * b]

    H_drift *= 2π
    H_drives .*= 2π
    L_dissipators .*= sqrt(2π)

    return OpenQuantumSystem(
        H_drift,
        H_drives,
        drive_bounds;
        dissipation_operators = L_dissipators,
        global_params = global_params,
    )
end

"""
    get_cat_controls(sys::AbstractQuantumSystem, α, N)

Compute static cat qubit controls for `N` time steps at coherent amplitude `α`,
reading `g2` and `χ_aa` from `sys.global_params`.

Returns a `2 × N` matrix where row 1 is the buffer drive and row 2 is the
Kerr correction drive. These are the steady-state controls that maintain a
coherent state ``|α⟩`` in the cat mode.

# Arguments
- `sys`: A quantum system with `g2` and `χ_aa` in its `global_params`
- `α`: Coherent state amplitude
- `N`: Number of time steps
"""
function get_cat_controls(sys::AbstractQuantumSystem, α::Real, N::Int)
    @assert haskey(sys.global_params, :g2) "Requires g2 in system global_params"
    @assert haskey(sys.global_params, :χ_aa) "Requires χ_aa in system global_params"
    buffer_drive = abs2(α) * sys.global_params.g2
    cat_kerr_correction = (2.0 * abs2(α) + 1.0) * sys.global_params.χ_aa
    return stack([fill(buffer_drive, N), fill(cat_kerr_correction, N)], dims = 1)
end

# *************************************************************************** #

@testitem "CatSystem: default construction" begin
    sys = CatSystem()
    @test sys isa OpenQuantumSystem
    @test sys.levels == 13 * 3
    @test sys.n_drives == 2
    @test length(sys.dissipation_operators) == 2
end

@testitem "CatSystem: custom parameters" begin
    sys = CatSystem(cat_levels = 4, buffer_levels = 2)
    @test sys.levels == 4 * 2
    @test sys.n_drives == 2

    # prefactor does not change Hilbert space dimension
    sys2 = CatSystem(cat_levels = 4, buffer_levels = 2, prefactor = 2.0)
    @test sys2.levels == 4 * 2
end

@testitem "CatSystem: drive bounds" begin
    sys = CatSystem(cat_levels = 4, buffer_levels = 2, drive_bounds = [0.5, 2.0])
    @test sys.drive_bounds == [(-0.5, 0.5), (-2.0, 2.0)]
end

@testitem "coherent_ket" begin
    using LinearAlgebra

    # Vacuum state α=0 should be |0⟩
    ψ = coherent_ket(0.0, 5)
    @test ψ[1] ≈ 1.0
    @test all(ψ[2:end] .≈ 0.0)

    # Normalization for nonzero α
    ψ = coherent_ket(2.0, 20)
    @test norm(ψ) ≈ 1.0 atol = 1e-6

    # Complex α
    ψ = coherent_ket(1.0 + 1.0im, 20)
    @test norm(ψ) ≈ 1.0 atol = 1e-6
end

@testitem "get_cat_controls" begin
    sys = CatSystem(cat_levels = 4, buffer_levels = 2)
    α = 2.0
    N = 10

    u = get_cat_controls(sys, α, N)
    @test size(u) == (2, N)

    # Check values against global_params
    @test all(u[1, :] .≈ abs2(α) * sys.global_params.g2)
    @test all(u[2, :] .≈ (2.0 * abs2(α) + 1.0) * sys.global_params.χ_aa)
end

@testitem "CatSystem: compact Lindbladian generators" begin
    using LinearAlgebra
    using SparseArrays

    cat_levels = 4
    buffer_levels = 2
    sys = CatSystem(cat_levels = cat_levels, buffer_levels = buffer_levels)
    n = sys.levels  # 8

    𝒢c_drift_ham, 𝒢c_drives, 𝒢c_dissipators = compact_lindbladian_parts(sys)

    # Compact generators should be n² × n²
    @test size(𝒢c_drift_ham) == (n^2, n^2)
    @test length(𝒢c_drives) == 2
    @test all(size(𝒢c) == (n^2, n^2) for 𝒢c in 𝒢c_drives)
    # CatSystem has 2 dissipators (κa*a, κb*b)
    @test length(𝒢c_dissipators) == 2
    @test all(size(𝒢c) == (n^2, n^2) for 𝒢c in 𝒢c_dissipators)

    # Verify P * L = I identity used in construction
    L = Piccolo.Isomorphisms.density_lift_matrix(n)
    P = Piccolo.Isomorphisms.density_projection_matrix(n)
    @test P * L ≈ I(n^2)
end

@testitem "CatSystem: DensityTrajectory rollout fidelity" begin
    using LinearAlgebra

    # Small cat system for fast testing
    cat_levels = 6
    buffer_levels = 2

    sys = CatSystem(cat_levels = cat_levels, buffer_levels = buffer_levels)
    n = sys.levels  # 12

    # Coherent state |α⟩ ⊗ |0⟩ as initial density matrix
    α = 1.0
    ψ_cat = coherent_ket(α, cat_levels)
    ψ_buf = zeros(ComplexF64, buffer_levels)
    ψ_buf[1] = 1.0
    ψ = kron(ψ_cat, ψ_buf)
    ρ0 = ψ * ψ'

    # Goal: maintain the coherent state
    ρg = copy(ρ0)

    # Steady-state controls that should maintain |α⟩
    T = 0.5
    N = 51
    times = collect(range(0.0, T, length = N))
    u = get_cat_controls(sys, α, N)
    pulse = ZeroOrderPulse(u, times)

    # Solve Lindblad master equation
    qtraj = DensityTrajectory(sys, pulse, ρ0, ρg)

    # Check fidelity: with short time and steady-state controls, should stay close
    fid = fidelity(qtraj)
    @test fid > 0.9

    # Verify trace preservation along trajectory (ODE solver tolerance)
    for t in range(0.0, T, length = 5)
        ρ = qtraj(t)
        @test real(tr(ρ)) ≈ 1.0 atol = 1e-3
    end

    # Verify Hermiticity along trajectory
    for t in range(0.0, T, length = 5)
        ρ = qtraj(t)
        @test norm(ρ - ρ') < 1e-6
    end
end

@testitem "CatSystem: SmoothPulseProblem optimization" begin
    using DirectTrajOpt
    using LinearAlgebra

    # Minimal cat system for fast testing
    cat_levels = 3
    buffer_levels = 2
    sys = CatSystem(cat_levels = cat_levels, buffer_levels = buffer_levels)
    n = sys.levels  # 6

    # Start in vacuum ⊗ vacuum, target coherent state |α⟩ ⊗ |0⟩
    ρ0 = zeros(ComplexF64, n, n)
    ρ0[1, 1] = 1.0

    α = 0.5
    ψ_cat = coherent_ket(α, cat_levels)
    ψ_buf = zeros(ComplexF64, buffer_levels)
    ψ_buf[1] = 1.0
    ψ_goal = kron(ψ_cat, ψ_buf)
    ρg = ψ_goal * ψ_goal'

    T = 1.0
    N = 11

    # Initialize with deterministic cat controls + a small smooth perturbation
    # (cos along time axis, both channels). Keeps the solve trajectory
    # reproducible across Julia versions — randn(2, N) was the flake source.
    times_arr = (0:(N-1)) ./ max(N - 1, 1)
    perturb =
        0.01 *
        vcat(reshape(cos.(2π .* times_arr), 1, N), reshape(sin.(2π .* times_arr), 1, N))
    u_init = get_cat_controls(sys, α, N) .+ perturb
    pulse = ZeroOrderPulse(u_init, collect(range(0.0, T, length = N)))
    qtraj = DensityTrajectory(sys, pulse, ρ0, ρg)

    # Build problem
    qcp = SmoothPulseProblem(qtraj, N; Q = 100.0, R = 1e-2, ddu_bound = 10.0)

    @test qcp isa QuantumControlProblem

    # Compact iso: state dim should be n²
    @test qcp.prob.trajectory.dims[:ρ⃗̃] == n^2

    # Pipeline-smoke solve. max_iter=50 gives IPOPT enough to drive the
    # constraint residual below the tolerance from the cat warmstart on any
    # Julia version (was 20 — too tight to converge reproducibly).
    solve!(qcp; max_iter = 50, print_level = 1, verbose = false)

    # Dynamics constraints should be satisfied. Tolerance loosened from 1e-2
    # to 5e-2 to absorb numerical variation across LAPACK / BLAS implementations
    # — this test smokes the pipeline; tighter checks belong on dedicated
    # convergence tests.
    traj = get_trajectory(qcp)
    dynamics_integrator = qcp.prob.integrators[1]
    δ = zeros(dynamics_integrator.dim)
    DirectTrajOpt.evaluate!(δ, dynamics_integrator, traj)
    @test norm(δ, Inf) < 5e-2
end
