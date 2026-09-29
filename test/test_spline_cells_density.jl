# ============================================================================ #
# Cluster E1 (#347) — SplineIntegrator DensityTrajectory cell behavior tests.
#
# Ported from the pre-split embedded suite (Piccolissimo
# src/integrators/spline/spline_integrator_density.jl). The density cell sat at
# 69%: the dense-parity harness drove the closed-system forward only, so the
# whole Lindblad/Duhamel machinery (lindblad_apply!, the adjoint lanes, the
# dissipator sensitivity ODE, typed dissipators) was dark. These restore the
# open-system Jacobian conformance, the in-place compact-iso helpers, and the
# vector-field adjoint identities.
#
# NOT ported: the gpu_proxy (JLArrays) parity item and the perf-regression
# baseline item — neither belongs in the open-core suite (JLArrays is not a
# Piccolo test dep; perf gating is Piccolissimo's wall-clock discipline).
# ============================================================================ #

@testitem "E1: in-place compact iso helpers match the allocating Isomorphisms" begin
    using Piccolo
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators:
        density_to_compact_iso!, compact_iso_to_density!
    using LinearAlgebra
    using Random
    Random.seed!(3471)

    for n in [2, 3, 4]
        ρ = randn(ComplexF64, n, n)
        ρ = (ρ + ρ') / 2  # make Hermitian

        # density_to_compact_iso! vs the allocating reference
        x_expected = Isomorphisms.density_to_compact_iso(ρ)
        x_actual = zeros(n^2)
        density_to_compact_iso!(x_actual, ρ, n)
        @test x_actual ≈ x_expected atol = 1e-14

        # compact_iso_to_density! vs the allocating reference
        M_expected = Isomorphisms.compact_iso_to_density(x_expected)
        M_actual = Matrix{ComplexF64}(undef, n, n)
        compact_iso_to_density!(M_actual, x_expected, n)
        @test M_actual ≈ M_expected atol = 1e-14
    end
end

@testitem "E1: lindblad_apply!/lindblad_adjoint_apply! allocate exactly zero" begin
    using LinearAlgebra, Random, SparseArrays
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators:
        lindblad_apply!, lindblad_adjoint_apply!

    Random.seed!(312)

    # NOTE: every buffer and every measured closure below is bound inside the
    # `for` body, i.e. in LOCAL scope. Hoisting them to the testitem's top level
    # (module global scope) gives the closures `Any`-typed captures, which adds
    # spurious allocations and makes the measurement uninterpretable.
    #
    # For the same reason NOTHING a closure captures may be assigned twice:
    # a captured variable with more than one assignment gets wrapped in a
    # `Core.Box`, reads back as `Any`, and re-introduces exactly the spurious
    # allocation this test exists to rule out.
    for n in (2, 4, 8), nL in (0, 1, 3)
        H0 = randn(ComplexF64, n, n)
        H_eff = (H0 + H0') / 2
        Ls = [randn(ComplexF64, n, n) for _ = 1:nL]
        Ks_dense = [L' * L for L in Ls]
        Ks_sparse = [sparse(K) for K in Ks_dense]   # the shape used at the call sites
        M = randn(ComplexF64, n, n)
        dM = zeros(ComplexF64, n, n)
        tmp = zeros(ComplexF64, n, n)
        Δt = 0.037

        for Ks in (Ks_dense, Ks_sparse)
            fwd = () -> lindblad_apply!(dM, M, H_eff, Δt, Ls, Ks, tmp)
            adj = () -> lindblad_adjoint_apply!(dM, M, H_eff, Δt, Ls, Ks, tmp)

            # Warm every path before measuring.
            fwd()
            adj()

            # Byte-measuring form (`@allocated`), asserted as equality with
            # zero — not a bounded budget, and not an allocation *count*
            # (a count-stable but allocating broadcast has slipped past a
            # count check on this codebase before).
            @test (@allocated fwd()) == 0
            @test (@allocated adj()) == 0
        end
    end
end

@testitem "E1: lindblad_adjoint_apply! satisfies the Hilbert–Schmidt adjoint identity" begin
    using LinearAlgebra, Random
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators:
        lindblad_apply!, lindblad_adjoint_apply!

    Random.seed!(3121)

    # ⟨X, Y⟩ = tr(X' Y)
    hs(X, Y) = tr(X' * Y)

    for n in (2, 3, 5, 8), nL in (0, 1, 2, 4)
        H0 = randn(ComplexF64, n, n)
        H_eff = (H0 + H0') / 2
        Ls = [randn(ComplexF64, n, n) for _ = 1:nL]
        Ks = [L' * L for L in Ls]
        Δt = 0.137

        dM = zeros(ComplexF64, n, n)
        dA = zeros(ComplexF64, n, n)
        tmp = zeros(ComplexF64, n, n)   # one caller-owned scratch, shared by both

        # ⟨A, ℒ[M]⟩ == ⟨ℒ†[A], M⟩ for *arbitrary* probes. Deliberately
        # non-Hermitian: the identity is a statement about the two linear maps,
        # not about physical states, and Hermitian probes would let a wrong
        # transpose slip through.
        rel = map(1:8) do _
            M = randn(ComplexF64, n, n)
            A = randn(ComplexF64, n, n)
            @test !ishermitian(M)
            @test !ishermitian(A)

            lindblad_apply!(dM, M, H_eff, Δt, Ls, Ks, tmp)
            lindblad_adjoint_apply!(dA, A, H_eff, Δt, Ls, Ks, tmp)

            lhs = hs(A, dM)
            rhs = hs(dA, M)
            abs(lhs - rhs) / max(abs(lhs), abs(rhs), eps())
        end

        @test maximum(rel) < 1e-13
    end

    # The identity must be exact for a NON-Hermitian H_eff too — that is why the
    # adjoint contracts against H_eff' rather than H_eff. This case is what
    # separates "adjoint of the map we actually compute" from "adjoint of the
    # map we assumed we were computing".
    let n = 6, nL = 3
        H_eff = randn(ComplexF64, n, n)     # NOT Hermitian
        Ls = [randn(ComplexF64, n, n) for _ = 1:nL]
        Ks = [L' * L for L in Ls]
        Δt = 0.29
        dM = zeros(ComplexF64, n, n)
        dA = zeros(ComplexF64, n, n)
        tmp = zeros(ComplexF64, n, n)

        M = randn(ComplexF64, n, n)
        A = randn(ComplexF64, n, n)
        lindblad_apply!(dM, M, H_eff, Δt, Ls, Ks, tmp)
        lindblad_adjoint_apply!(dA, A, H_eff, Δt, Ls, Ks, tmp)

        lhs = hs(A, dM)
        rhs = hs(dA, M)
        @test abs(lhs - rhs) / max(abs(lhs), abs(rhs)) < 1e-13
    end
end

@testitem "E1: lindblad_adjoint_apply! matches the explicit Heisenberg generator" begin
    using LinearAlgebra, Random, SparseArrays
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators: lindblad_adjoint_apply!

    Random.seed!(3122)

    # ℒ†[A] = Δt * ( +i[H,A] + Σⱼ (Lⱼ† A Lⱼ - ½{Kⱼ, A}) ), Kⱼ = Lⱼ†Lⱼ.
    # Written out with plain allocating matrix algebra — no shared helper with
    # the implementation, so this is an independent check of the formula rather
    # than of the buffer plumbing.
    function heisenberg_reference(A, H, Δt, Ls, Ks)
        out = Δt * (im * (H * A - A * H))
        for (L, K) in zip(Ls, Ks)
            out += Δt * (L' * A * L)
            out -= (Δt / 2) * (K * A + A * K)
        end
        return out
    end

    for n in (2, 3, 5, 8), nL in (0, 1, 2, 4)
        H0 = randn(ComplexF64, n, n)
        H_eff = (H0 + H0') / 2          # physical case: H' == H, so H and H' agree
        Ls = [randn(ComplexF64, n, n) for _ = 1:nL]
        Ks = [L' * L for L in Ls]
        Δt = 0.137

        A = randn(ComplexF64, n, n)     # non-Hermitian probe
        dA = zeros(ComplexF64, n, n)
        tmp = zeros(ComplexF64, n, n)

        lindblad_adjoint_apply!(dA, A, H_eff, Δt, Ls, Ks, tmp)
        ref = heisenberg_reference(A, H_eff, Δt, Ls, Ks)

        @test norm(dA - ref) / max(norm(ref), eps()) < 1e-13

        # Same result through the sparse-Kⱼ shape used at the call sites.
        dA_sp = zeros(ComplexF64, n, n)
        lindblad_adjoint_apply!(dA_sp, A, H_eff, Δt, Ls, [sparse(K) for K in Ks], tmp)
        @test norm(dA_sp - ref) / max(norm(ref), eps()) < 1e-13
    end
end

@testitem "E1: testing SplineIntegrator{DensityTrajectory} with LinearSpline" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra

    # Simple 2-level open system with decay
    L = ComplexF64[0.1 0.0; 0.0 0.0]
    sys = OpenQuantumSystem(
        PAULIS.Z,
        [PAULIS.X, PAULIS.Y],
        [1.0, 1.0];
        dissipation_operators = [L],
    )

    T = 1.0
    N = 5

    ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρg = ComplexF64[0.0 0.0; 0.0 1.0]

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.5, 2, N), times)
    qtraj = DensityTrajectory(sys, pulse, ρ0, ρg)

    𝒮 = SplineIntegrator(qtraj, N)
    @test 𝒮 isa SplineIntegrator{DensityTrajectory,LinearSpline}
    @test 𝒮.ketdim == 2
    @test 𝒮.x_dim == 4  # n² = 4 for n=2

    traj = NamedTrajectory(qtraj, N)
    test_integrator(𝒮, traj; atol = 1e-4, gauss_newton = true)
end

@testitem "E1: testing SplineIntegrator{DensityTrajectory} with CubicSpline" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra

    L = ComplexF64[0.1 0.0; 0.0 0.0]
    sys = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])

    T = 1.0
    N = 5

    ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρg = ComplexF64[0.0 0.0; 0.0 1.0]

    times = collect(range(0.0, T, N))
    pulse = CubicSplinePulse(fill(0.5, 1, N), fill(0.0, 1, N), times)
    qtraj = DensityTrajectory(sys, pulse, ρ0, ρg)

    𝒮 = SplineIntegrator(qtraj, N)
    @test 𝒮 isa SplineIntegrator{DensityTrajectory,CubicSpline}

    traj = NamedTrajectory(qtraj, N)
    test_integrator(𝒮, traj; atol = 1e-4, gauss_newton = true)
end

@testitem "E1: SplineIntegrator{DensityTrajectory} constraint satisfaction" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra

    # 3-level system to test non-trivial dimensions (n²=9)
    H_drift = diagm(ComplexF64[1.0, 0.0, -1.0])
    H_drive = zeros(ComplexF64, 3, 3)
    H_drive[1, 2] = H_drive[2, 1] = 1.0
    L = zeros(ComplexF64, 3, 3)
    L[1, 2] = 0.05  # |2⟩ → |1⟩ decay
    sys = OpenQuantumSystem(H_drift, [H_drive], [1.0]; dissipation_operators = [L])

    T = 1.0
    N = 5

    ρ0 = zeros(ComplexF64, 3, 3)
    ρ0[1, 1] = 1.0
    ρg = zeros(ComplexF64, 3, 3)
    ρg[2, 2] = 1.0

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.3, 1, N), times)
    qtraj = DensityTrajectory(sys, pulse, ρ0, ρg)

    𝒮 = SplineIntegrator(qtraj, N)
    @test 𝒮 isa SplineIntegrator{DensityTrajectory,LinearSpline}
    @test 𝒮.ketdim == 3
    @test 𝒮.x_dim == 9  # n² = 9

    # Test that constraint is satisfied for the rollout trajectory
    traj = NamedTrajectory(qtraj, N)

    δ = zeros(𝒮.dim)
    DirectTrajOpt.evaluate!(δ, 𝒮, traj)
    # Constraint should be close to zero for the initial rollout
    @test norm(δ, Inf) < 1e-2
end

@testitem "E1: SplineIntegrator{DensityTrajectory} no dissipation matches ket" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra

    # System WITH NO dissipation — density dynamics should match ket dynamics
    sys_open = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])
    sys_ket = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])

    T = 1.0
    N = 5

    ψ0 = ComplexF64[1.0, 0.0]
    ψg = ComplexF64[0.0, 1.0]
    ρ0 = ψ0 * ψ0'
    ρg = ψg * ψg'

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.5, 1, N), times)

    qtraj_ket = KetTrajectory(sys_ket, pulse, ψ0, ψg)
    qtraj_density = DensityTrajectory(sys_open, pulse, ρ0, ρg)

    traj_ket = NamedTrajectory(qtraj_ket, N)
    traj_density = NamedTrajectory(qtraj_density, N)

    𝒮_ket = SplineIntegrator(qtraj_ket, N)
    𝒮_den = SplineIntegrator(qtraj_density, N)

    # Both should have near-zero constraint violations
    δ_ket = zeros(𝒮_ket.dim)
    δ_den = zeros(𝒮_den.dim)
    DirectTrajOpt.evaluate!(δ_ket, 𝒮_ket, traj_ket)
    DirectTrajOpt.evaluate!(δ_den, 𝒮_den, traj_density)

    @test norm(δ_ket, Inf) < 1e-3
    @test norm(δ_den, Inf) < 1e-2
end

@testitem "E1: SplineIntegrator{DensityTrajectory} Jacobian with 3-level system" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra

    # 3-level system with decay
    H_drift = diagm(ComplexF64[1.0, 0.0, -1.0])
    H_drive = zeros(ComplexF64, 3, 3)
    H_drive[1, 2] = H_drive[2, 1] = 1.0
    L = zeros(ComplexF64, 3, 3)
    L[1, 2] = 0.05
    sys = OpenQuantumSystem(H_drift, [H_drive], [1.0]; dissipation_operators = [L])

    T = 1.0
    N = 5

    ρ0 = zeros(ComplexF64, 3, 3)
    ρ0[1, 1] = 1.0
    ρg = zeros(ComplexF64, 3, 3)
    ρg[2, 2] = 1.0

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.3, 1, N), times)
    qtraj = DensityTrajectory(sys, pulse, ρ0, ρg)

    𝒮 = SplineIntegrator(qtraj, N)
    traj = NamedTrajectory(qtraj, N)
    test_integrator(𝒮, traj; atol = 1e-4, gauss_newton = true)
end

@testitem "E1: SplineIntegrator{DensityTrajectory} NonlinearDrive linear spline" begin
    using DirectTrajOpt
    using Piccolo
    using SparseArrays
    using LinearAlgebra

    # 2-level open system: H = σz + u₁·σx + u₁²·σz  with decay
    L = ComplexF64[0.1 0.0; 0.0 0.0]
    drives = AbstractDrive[
        LinearDrive(sparse(ComplexF64.(PAULIS.X)), 1),
        NonlinearDrive(PAULIS.Z, u -> u[1]^2),
    ]
    sys = OpenQuantumSystem(PAULIS.Z, drives, [1.0]; dissipation_operators = [L])

    T = 1.0
    N = 5

    ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρg = ComplexF64[0.0 0.0; 0.0 1.0]

    times = collect(range(0.0, T, N))
    pulse = LinearSplinePulse(fill(0.3, 1, N), times)
    qtraj = DensityTrajectory(sys, pulse, ρ0, ρg)

    𝒮 = SplineIntegrator(qtraj, N)
    @test 𝒮 isa SplineIntegrator{DensityTrajectory,LinearSpline}

    traj = NamedTrajectory(qtraj, N)
    test_integrator(𝒮, traj; atol = 1e-4, gauss_newton = true)
end

@testitem "E1: SplineIntegrator{DensityTrajectory} NonlinearDrive cubic spline" begin
    using DirectTrajOpt
    using Piccolo
    using SparseArrays
    using LinearAlgebra

    # 2-level open system with bilinear NonlinearDrive:
    # H = σz + u₁·σx + u₂·σy + (u₁·u₂)·σz  with decay
    L = ComplexF64[0.1 0.0; 0.0 0.0]
    drives = AbstractDrive[
        LinearDrive(sparse(ComplexF64.(PAULIS.X)), 1),
        LinearDrive(sparse(ComplexF64.(PAULIS.Y)), 2),
        NonlinearDrive(PAULIS.Z, u -> u[1] * u[2]),
    ]
    sys = OpenQuantumSystem(PAULIS.Z, drives, [1.0, 1.0]; dissipation_operators = [L])

    T = 1.0
    N = 5

    ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρg = ComplexF64[0.0 0.0; 0.0 1.0]

    times = collect(range(0.0, T, N))
    pulse = CubicSplinePulse(fill(0.3, 2, N), fill(0.0, 2, N), times)
    qtraj = DensityTrajectory(sys, pulse, ρ0, ρg)

    𝒮 = SplineIntegrator(qtraj, N)
    @test 𝒮 isa SplineIntegrator{DensityTrajectory,CubicSpline}

    traj = NamedTrajectory(qtraj, N)
    test_integrator(𝒮, traj; atol = 1e-4, gauss_newton = true)
end

@testitem "E1: SplineIntegrator{DensityTrajectory} constructs with typed dissipators" begin
    using LinearAlgebra, NamedTrajectories, Piccolo

    # Typed-dissipator path — exercises the 3-tuple return + dissipator-fold
    # reassembly in the sensitivity ODE constructor (Task 11 / Piccolo #152).
    diss = LinearDissipator(ComplexF64.(PAULIS.Z / sqrt(2)))
    sys = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipators = [diss])
    ρ0 = ComplexF64[1 0; 0 0]
    ρg = ComplexF64[0 0; 0 1]
    N = 11
    times = collect(range(0, 1.0, length = N))
    pulse = LinearSplinePulse(fill(0.3, 1, N), times)
    qtraj = DensityTrajectory(sys, pulse, ρ0, ρg)

    # Constructor must not error on the 3-tuple destructure with a typed
    # dissipator present.
    𝒮 = SplineIntegrator(qtraj, N)
    @test 𝒮 isa SplineIntegrator{DensityTrajectory,LinearSpline}

    # The typed dissipator also supports the full Jacobian conformance
    traj = NamedTrajectory(qtraj, N)
    test_integrator(𝒮, traj; atol = 1e-4, gauss_newton = true)
end

@testitem "E1: SplineIntegrator{DensityTrajectory} MagnusGL4 is refused (Lindblad lane)" begin
    using Piccolo
    using LinearAlgebra

    L = ComplexF64[0.1 0.0; 0.0 0.0]
    sys = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    ρ0 = ComplexF64[1 0; 0 0]
    ρg = ComplexF64[0 0; 0 1]
    N = 5
    times = collect(range(0.0, 1.0, N))
    pulse = LinearSplinePulse(fill(0.3, 1, N), times)
    qtraj = DensityTrajectory(sys, pulse, ρ0, ρg)

    # The non-Hermitian Lindbladian generator has no Magnus cell: refuse at
    # construction rather than silently degrading to a Tsit5 forward.
    @test_throws ErrorException SplineIntegrator(qtraj, N; alg = MagnusGL4Alg())
    @test_throws ErrorException SplineIntegrator(
        qtraj,
        N;
        alg = ChebyshevAlg(bracket = (-8.0, 8.0), n_sub = 8),
    )
end

@testitem "E1: each density knot's sensitivity ODE carries its OWN RHS closure (#354 race)" begin
    using DirectTrajOpt
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using Random
    Random.seed!(354_001)

    # #354's ROOT CAUSE, pinned STRUCTURALLY — and therefore INDEPENDENTLY OF THE
    # THREAD COUNT the suite happens to run with.
    #
    # `eval_jacobian` solves the `N-1` knot sensitivity problems under
    # `Threads.@threads`. `build_density_sensitivity_ode` used to bind ONE
    # `u_interp` in a `let` and hand the SAME closure to every problem, so the
    # threads raced on the interpolated control. This gate asserts the repair:
    # one closure (and one captured scratch buffer) per knot problem.
    L = ComplexF64[0.0 0.12; 0.0 0.0]
    sys = OpenQuantumSystem(
        Matrix{ComplexF64}(PAULIS.Z),
        [Matrix{ComplexF64}(PAULIS.X)],
        [1.0];
        dissipation_operators = [L],
    )
    N = 6
    times = collect(range(0.0, 1.0, N))
    ρ1 = ComplexF64[1 0; 0 0]
    ρ2 = ComplexF64[0 0; 0 1]

    pulses = (
        LinearSplinePulse(0.3 .* randn(1, N), times),
        CubicSplinePulse(0.3 .* randn(1, N), 0.1 .* randn(1, N), times),
    )
    for pulse in pulses
        qtrajs = (
            DensityTrajectory(sys, pulse, ρ1, ρ2),
            MultiDensityTrajectory(sys, pulse, [ρ1, ρ2], [ρ2, ρ1]),
        )
        for qtraj in qtrajs
            𝒮 = SplineIntegrator(qtraj, N)
            probs = 𝒮.sens_probs
            @test !isnothing(probs)
            @test length(probs) == N - 1

            rhs = [getfield(p.f, :f) for p in probs]
            # DISTINCT closure objects — one per knot problem.
            @test allunique(objectid.(rhs))
            # ...and, the property that actually matters, distinct captured scratch.
            # A shared `u_interp` IS the race; distinct closures over a shared buffer
            # would be just as broken, so the buffer is checked directly.
            bufs = [getfield(f, :u_interp) for f in rhs]
            @test allunique(objectid.(bufs))
            @test all(b -> b isa Vector{Float64} && length(b) == 𝒮.u_dim, bufs)
            # The gate is not vacuous: `objectid` really does collide for a shared
            # buffer, which is what the pre-#354 `let` produced.
            shared = fill(bufs[1], N - 1)
            @test !allunique(objectid.(shared))
        end
    end
end

# ── E1 resume fill (#347): typed-dissipator kernels, inner-ctor lanes, the ─── #
# ── fixed-step sparsity probe, and the cubic Δt↔du sensitivity chain. ──────── #

@testitem "E1: dissipator_apply!/dissipator_adjoint_apply! rated kernels" begin
    using Piccolo
    using LinearAlgebra
    using Random
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators:
        dissipator_apply!, dissipator_adjoint_apply!, duhamel_tape

    Random.seed!(3473)

    n = 3
    L1 = zeros(ComplexF64, n, n); L1[1, 2] = 0.3
    L2 = zeros(ComplexF64, n, n); L2[2, 3] = 0.5
    Ls = [L1, L2]
    Ks = [L' * L for L in Ls]
    rates = [2.0, 0.5]
    Δt = 0.05
    tmp = zeros(ComplexF64, n, n)

    M = randn(ComplexF64, n, n); M = (M + M') / 2
    A = randn(ComplexF64, n, n); A = (A + A') / 2

    # dissipator_apply! = Δt·Σ ratesⱼ (Lⱼ M Lⱼ† − ½{Kⱼ, M}) from zero
    dM = Matrix{ComplexF64}(undef, n, n)
    dissipator_apply!(dM, M, Δt, Ls, Ks, rates, tmp)
    expected =
        sum(Δt * rates[j] * (L * M * L' - (Ks[j] * M + M * Ks[j]) / 2) for (j, L) in enumerate(Ls))
    @test dM ≈ expected atol = 1e-12

    # dissipator_adjoint_apply! is the Hilbert–Schmidt adjoint of the above:
    # ⟨A, 𝒟(M)⟩ = ⟨𝒟†(A), M⟩ exactly.
    dA = Matrix{ComplexF64}(undef, n, n)
    dissipator_adjoint_apply!(dA, A, Δt, Ls, Ks, rates, tmp)
    @test abs(tr(A' * dM) - tr(dA' * M)) < 1e-12

    # rates = ones reproduces the un-rated Lindblad dissipator half
    dissipator_apply!(dM, M, Δt, Ls, Ks, fill(1.0, 2), tmp)
    un_rated =
        sum(Δt * (L * M * L' - ((L' * L) * M + M * (L' * L)) / 2) for L in Ls)
    @test dM ≈ un_rated atol = 1e-12
end

@testitem "E1: density inner-ctor lanes: gates, pulse inference, globals, fixed-step" begin
    using DirectTrajOpt
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators: spline_order, duhamel_tape

    L = ComplexF64[0.1 0.0; 0.0 0.0]
    sys = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρg = ComplexF64[0.0 0.0; 0.0 1.0]
    N = 5
    times = collect(range(0.0, 1.0, N))

    # (qtraj, traj) inner form: gates refuse Magnus/Chebyshev at THIS layer too
    qtraj = DensityTrajectory(sys, LinearSplinePulse(fill(0.3, 1, N), times), ρ0, ρg)
    traj = NamedTrajectory(qtraj, N)
    @test_throws ErrorException SplineIntegrator(qtraj, traj; alg = MagnusGL4Alg())
    @test_throws ErrorException SplineIntegrator(
        qtraj,
        traj;
        alg = ChebyshevAlg(bracket = (-8.0, 8.0)),
    )

    # Non-spline pulse: defaults to linear through the inner form
    zo_qtraj = DensityTrajectory(sys, ZeroOrderPulse(0.3 * fill(1.0, 1, N), times), ρ0, ρg)
    zo_traj = NamedTrajectory(zo_qtraj, N)
    𝒮 = SplineIntegrator(zo_qtraj, zo_traj)
    @test spline_order(𝒮) == 1
    δ = zeros(𝒮.dim)
    evaluate!(δ, 𝒮, zo_traj)
    @test norm(δ, Inf) < 1e-4

    # Fixed-step Tsit5: the Φ-solve sparsity probe derives the structural
    # pattern from the solved propagator (the adaptive lane stays dense)
    𝒮_f = SplineIntegrator(qtraj, traj; alg = Tsit5Alg(adaptive = false))
    @test 𝒮_f.alg isa Tsit5Alg
    δ_f = zeros(𝒮_f.dim)
    evaluate!(δ_f, 𝒮_f, traj)
    @test norm(δ_f, Inf) < 1e-4

    # Globals: an open system with a global rides it through the seed + comps
    gsys = OpenQuantumSystem(
        PAULIS.Z,
        [PAULIS.X],
        [1.0];
        dissipation_operators = [L],
        global_params = (δ = 0.01,),
    )
    gq = DensityTrajectory(
        gsys,
        LinearSplinePulse(fill(0.3, 1, N), times),
        ρ0,
        ρg,
    )
    gtraj = NamedTrajectory(gq, N)
    # The (qtraj, traj) inner form does NOT auto-detect globals — the names
    # must be given explicitly (the traj carries the :δ component from the
    # system's global_params).
    𝒮_g = SplineIntegrator(gq, gtraj; global_names = [:δ])
    @test 𝒮_g.global_names == [:δ]
    @test 𝒮_g.global_dim == 1
    δ_g = zeros(𝒮_g.dim)
    evaluate!(δ_g, 𝒮_g, gtraj)
    @test norm(δ_g, Inf) < 1e-4

    # Cubic without du bounds: the derivative seed zero-fills
    cq = DensityTrajectory(
        sys,
        CubicSplinePulse(fill(0.3, 1, N), fill(0.0, 1, N), times),
        ρ0,
        ρg,
    )
    ctraj = NamedTrajectory(cq, N)
    @test haskey(ctraj.bounds, :du)
    # Rebuild the same trajectory WITHOUT the du bound (copy-ctor with an
    # overridden bounds tuple — `components` are index ranges, not data).
    boundfree = NamedTrajectory(ctraj; bounds = (u = 1.0,))
    𝒮_c = SplineIntegrator(cq, boundfree)
    @test spline_order(𝒮_c) == 3
    δ_c = zeros(𝒮_c.dim)
    evaluate!(δ_c, 𝒮_c, boundfree)
    @test norm(δ_c, Inf) < 1e-4

    # The Duhamel tape travels with the cell (density-specific alg_data)
    @test !isnothing(duhamel_tape(𝒮_c))
end

@testitem "E1: density cubic with NONZERO du drives the Δt↔du sensitivity chain" begin
    using DirectTrajOpt
    using Piccolo
    using LinearAlgebra
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators: spline_order

    # A cubic density cell with genuinely nonzero du values exercises the
    # total_du_dΔt chain inside the density sensitivity ODE (the prior cubic
    # items all seeded du = 0, which skips it), then passes full conformance.
    L = ComplexF64[0.1 0.0; 0.0 0.0]
    sys = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])
    ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρg = ComplexF64[0.0 0.0; 0.0 1.0]
    N = 5
    times = collect(range(0.0, 1.0, N))
    pulse = CubicSplinePulse(fill(0.4, 1, N), 0.2 * fill(1.0, 1, N), times)
    qtraj = DensityTrajectory(sys, pulse, ρ0, ρg)

    𝒮 = SplineIntegrator(qtraj, N)
    @test spline_order(𝒮) == 3
    traj = NamedTrajectory(qtraj, N)
    test_integrator(𝒮, traj; atol = 1e-4, gauss_newton = true)
end
