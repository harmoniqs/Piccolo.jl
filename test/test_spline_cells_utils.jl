# ============================================================================ #
# Cluster E1 (#347) — integrators_spline utility & algorithm-layer tests.
#
# Covers the family's 0% utility modules that no cell test can reach
# incidentally: the knot-point sensitivity kicks (knot_point_sensitivity.jl),
# the cached VJP/JVP machinery (vjp_jvp.jl), the cubic-Hermite interval
# coefficient kernel (spline_interval_coeffs.jl), the complex↔real boundary
# (complex_real_interface.jl), PropagatorResult's derived views
# (propagator_result.jl), the algorithm types' validation lanes
# (algorithms.jl), Tsit5Data's convenience constructors (alg_data.jl), the
# SplineType trait layer (spline_types.jl), and the family common module's
# unitary_rollout_trajectory (_spline_integrators.jl).
#
# Parity idiom: analytic directional derivatives against finite differences /
# ForwardDiff through the same forward maps — the same discipline the
# exponential-integrator family's testitems use.
# ============================================================================ #

# ── knot_point_sensitivity.jl ─────────────────────────────────────────────── #

@testitem "E1: gauss_legendre_01 nodes/weights are exact for polynomials" begin
    using Piccolo
    using LinearAlgebra

    for n in (2, 4, 7)
        nodes, weights = gauss_legendre_01(n)
        @test length(nodes) == n
        @test length(weights) == n
        @test all(0 .<= nodes .<= 1)
        # Partition of unity on [0, 1]
        @test sum(weights) ≈ 1 atol = 1e-12
        # Gauss–Legendre integrates degree < 2n polynomials EXACTLY on [0, 1]
        for m = 0:(2n-1)
            @test sum(weights .* nodes .^ m) ≈ 1 // (m + 1) atol = 1e-12
        end
        # And is NOT exact at degree 2n (the quadrature boundary — pins the
        # order, guarding against a lower-order rule passing the tests above)
        @test abs(sum(weights .* nodes .^ (2n)) - 1 / (2n + 1)) > 1e-12
    end

    # Known closed form: the 1-point rule is the midpoint
    n1, w1 = gauss_legendre_01(1)
    @test n1 ≈ [0.5] atol = 1e-12
    @test w1 ≈ [1.0] atol = 1e-12
end

@testitem "E1: compute_sensitivity_kick_exact matches a finite difference of the propagated state" begin
    using Piccolo
    using LinearAlgebra
    using Random
    Random.seed!(3472)

    n = 3
    H_drift = diagm(ComplexF64[0.0, 1.0, 2.0])
    H_drive = ComplexF64[0 1 0; 1 0 1; 0 1 0] / sqrt(2)
    Δt = 0.05
    ψ = normalize(randn(ComplexF64, n))
    u = 0.4

    H_of(u) = H_drift + u * H_drive
    # ∂(exp(-iΔtH)ψ)/∂u by central finite difference of the exact exponential
    ε = 1e-6
    prop(u) = exp(-im * Δt * H_of(u)) * ψ
    kick_fd = (prop(u + ε) - prop(u - ε)) / (2ε)

    # Exact quadrature kick: eigendecomposition of H at u
    λ, V = eigen(Hermitian(H_of(u)))
    kick = compute_sensitivity_kick_exact(H_drive, ψ, Δt, λ, V; n_quad = 8)
    @test norm(kick - kick_fd) / max(norm(kick_fd), eps()) < 1e-6

    # The quadrature order knob really tightens: n_quad=1 is a one-node rule,
    # still accurate for this mild generator but measurably coarser than 8
    kick_coarse = compute_sensitivity_kick_exact(H_drive, ψ, Δt, λ, V; n_quad = 1)
    @test norm(kick_coarse - kick_fd) > norm(kick - kick_fd)

    # First-order kick: -iΔt·H_drive·ψ_next — the leading-order approximation.
    # It agrees with the exact kick at O(Δt·‖H‖·kick) — loose but nonzero.
    ψ_next = prop(u)
    kick_fo = compute_sensitivity_kick_first_order(H_drive, ψ_next, Δt)
    @test norm(kick_fo - kick_fd) < 5 * Δt * norm(H_drive) * norm(kick_fd)
    # And the two approximations converge to each other as Δt → 0
    λs, Vs = eigen(Hermitian(H_of(u)))
    for Δt_small in (0.02, 0.01, 0.005)
        prop_s = exp(-im * Δt_small * H_of(u)) * ψ
        kick_e = compute_sensitivity_kick_exact(H_drive, ψ, Δt_small, λs, Vs; n_quad = 8)
        kick_o = compute_sensitivity_kick_first_order(H_drive, prop_s, Δt_small)
        @test norm(kick_o - kick_e) / max(norm(kick_e), eps()) < 10 * Δt_small
    end
end

# ── vjp_jvp.jl ────────────────────────────────────────────────────────────── #

@testitem "E1: setup_knot_point_propagation caches the exact forward rollout" begin
    using Piccolo
    using LinearAlgebra
    using Random
    Random.seed!(3473)

    n = 2
    sys = QuantumSystem(0.3 * ComplexF64.(PAULIS.Z), [ComplexF64.(PAULIS.X)], [1.0];)
    @test sys.n_drives == 1

    N = 6
    controls = 0.2 .* randn(1, N)
    Δts = fill(0.05, N - 1)
    ψ0 = ComplexF64[1.0, 0.0]

    data = setup_knot_point_propagation(sys, controls, Δts, ψ0)

    @test data isa KnotPointPropagationData
    @test length(data.propagators) == N - 1
    @test length(data.states) == N
    @test size(data.controls) == (1, N)
    @test data.controls ≈ controls
    @test data.Δts == Δts

    # states[1] is the seed; each states[j+1] = Φ_j * states[j] by construction
    @test data.states[1] ≈ ψ0
    for j = 1:(N-1)
        @test data.states[j+1] ≈ data.propagators[j] * data.states[j] atol = 1e-12
    end

    # The cached propagators are the exact matrix exponentials
    for j = 1:(N-1)
        H = Matrix{ComplexF64}(sys.H_drift) + controls[1, j] .* sys.H_drives[1].H
        @test data.propagators[j] ≈ exp(-im * Δts[j] * H) atol = 1e-10
    end

    # Every cached propagator is unitary (Lindblad-free chain stays on the group)
    for U in data.propagators
        @test norm(U'U - I) < 1e-12
    end
end

@testitem "E1: ket_vjp and ket_jvp agree with finite differences of the rollout" begin
    using Piccolo
    using LinearAlgebra
    using Random
    Random.seed!(3474)

    n = 2
    sys = QuantumSystem(0.3 * ComplexF64.(PAULIS.Z), [ComplexF64.(PAULIS.X)], [1.0];)
    N = 6
    controls = 0.2 .* randn(1, N)
    Δts = fill(0.05, N - 1)
    ψ0 = ComplexF64[1.0, 0.0]

    data = setup_knot_point_propagation(sys, controls, Δts, ψ0)

    # rollout(u) through the SAME per-knot eigendecomposition machinery
    function rollout(ctrls)
        ψ = ComplexF64.(ψ0)
        for j = 1:(N-1)
            H = Matrix{ComplexF64}(sys.H_drift) + ctrls[1, j] .* sys.H_drives[1].H
            ψ = exp(-im * Δts[j] * H) * ψ
        end
        return ψ
    end

    λ_seed = randn(ComplexF64, n)
    loss(ctrls) = real(dot(λ_seed, rollout(ctrls)))

    # ── VJP: grad[k, j] must be ∂Re(λᵀψ_N)/∂u_{k,j} ──────────────────────────
    grad = ket_vjp(data, λ_seed; use_exact_kick = true, n_quad = 8)
    @test size(grad) == (1, N - 1)
    ε = 1e-6
    for j = 1:(N-1)
        p = copy(controls);
        p[1, j] += ε
        m = copy(controls);
        m[1, j] -= ε
        fd = (loss(p) - loss(m)) / (2ε)
        @test grad[1, j] ≈ fd atol = 1e-6
    end

    # ── JVP: δψ_N in a control tangent direction ────────────────────────────
    # The tangent only spans the N-1 INTERVAL controls (the chain never reads
    # the last knot's control), so the FD perturbs only those columns.
    v = 0.5 .* randn(1, N - 1)
    δψ = ket_jvp(data, v; use_exact_kick = true, n_quad = 8)
    @test δψ isa Vector{ComplexF64}
    @test length(δψ) == n
    perturb(vs) = begin
        c = copy(controls)
        c[:, 1:(N-1)] .+= vs
        return c
    end
    δψ_fd = (rollout(perturb(+ε .* v)) - rollout(perturb(-ε .* v))) / (2ε)
    @test norm(δψ - δψ_fd) / max(norm(δψ_fd), eps()) < 1e-6

    # ── Adjoint identity through the SAME cached chain: ⟨J·v, λ⟩ = ⟨v, Jᵀλ⟩ ──
    lhs = real(dot(δψ, λ_seed))
    rhs = dot(vec(v), grad)
    @test isapprox(lhs, rhs; rtol = 1e-7, atol = 1e-9)

    # ── First-order-kick lanes are wired (use_exact_kick=false) and converge ─
    grad_fo = ket_vjp(data, λ_seed; use_exact_kick = false)
    @test size(grad_fo) == (1, N - 1)
    @test !iszero(grad_fo)
    δψ_fo = ket_jvp(data, v; use_exact_kick = false)
    @test !iszero(δψ_fo)
    # The first-order kick is the LEADING-ORDER approximation — its relative
    # error is O(Δt·‖H‖) (Δt=0.05, ‖H‖≈0.5 here), NOT quadrature-tight.
    @test norm(δψ_fo - δψ_fd) / max(norm(δψ_fd), eps()) < 0.05
end

# ── spline_interval_coeffs.jl ─────────────────────────────────────────────── #

@testitem "E1: SplineIntervalCoeffs forward/directional/VJP/HVP parity (incl. NonlinearDrive)" begin
    using Piccolo
    using SparseArrays
    using LinearAlgebra
    using Random
    using ForwardDiff
    Random.seed!(3475)

    u_dim = 2
    # A linear and a quadratic drive so drive_coeff_jac/hess are both live
    drives = [
        LinearDrive(sparse(ComplexF64.(PAULIS.X)), 1),
        NonlinearDrive(PAULIS.Z, u -> u[1]^2),
    ]

    p = [0.4, -0.2, 0.5, 0.3, 0.1, -0.05, 0.02, -0.01, 0.25, 0.0]
    @test length(p) == 4 * u_dim + 2
    sic = SplineIntervalCoeffs(drives, u_dim, p)

    # Constructor error lane: packed length must be 4·u_dim + 2
    @test_throws AssertionError SplineIntervalCoeffs(drives, u_dim, zeros(4 * u_dim + 1))

    τ = 0.3

    # The DOCUMENTED Hermite-basis spec (the file header), used as the analytic
    # reference the evaluator is asserted against — ForwardDiff-differentiable,
    # unlike the Float64-buffered evaluator itself.
    function u_at_ref(τr, pp)
        dΔt = pp[4*u_dim+1]
        τ2 = τr * τr
        τ3 = τ2 * τr
        h00 = 2τ3 - 3τ2 + 1
        h01 = -2τ3 + 3τ2
        h10 = (τ3 - 2τ2 + τr) * dΔt
        h11 = (τ3 - τ2) * dΔt
        return [
            h00 * pp[i] + h01 * pp[u_dim+i] + h10 * pp[2u_dim+i] + h11 * pp[3u_dim+i] for
            i = 1:u_dim
        ]
    end
    c_ref(l, pp) = drive_coeff(drives[l], u_at_ref(τ, pp))

    # ── forward: c_l(τ) = drive_coeff(d_l, u(τ)) — matches the documented basis
    c = zeros(2)
    interval_coeff!(c, sic, τ)
    @test c ≈ [c_ref(1, p), c_ref(2, p)] atol = 1e-14

    # ── directional: interval_coeff_dir! == ∇c·δp ─────────────────────────────
    δp = 0.3 .* randn(length(p))
    dc = zeros(2)
    interval_coeff_dir!(dc, sic, τ, δp)
    for l = 1:2
        fd = (c_ref(l, p .+ 1e-6 .* δp) - c_ref(l, p .- 1e-6 .* δp)) / (2e-6)
        @test dc[l] ≈ fd atol = 1e-7
    end

    # ── VJP scatter: one pairing scatters into ∇c_l (real part) ───────────────
    for l = 1:2
        g = zeros(length(p))
        interval_vjp_scatter!(g, sic, l, τ, 1.0 + 0.0im)
        ∇c_l = ForwardDiff.gradient(pp -> c_ref(l, pp), p)
        @test g ≈ ∇c_l atol = 1e-9

        # A complex pairing contributes ONLY its real part
        g_c = zeros(length(p))
        interval_vjp_scatter!(g_c, sic, l, τ, 2.0 + 3.0im)
        @test g_c ≈ 2 .* ∇c_l atol = 1e-9
    end

    # ── HVP scatter: the directional derivative of the VJP map ───────────────
    for l = 1:2
        ∇c_l = pp -> ForwardDiff.gradient(q -> c_ref(l, q), pp)
        hvp_ref = (∇c_l(p .+ 1e-6 .* δp) - ∇c_l(p .- 1e-6 .* δp)) / (2e-6)
        g = zeros(length(p))
        interval_hvp_scatter!(g, sic, l, τ, 1.0 + 0.0im, δp)
        @test g ≈ hvp_ref atol = 1e-6 rtol = 1e-5

        # val scales linearly (it is the SAME first-order pairing the VJP takes)
        g_v = zeros(length(p))
        interval_hvp_scatter!(g_v, sic, l, τ, -1.5 + 0.5im, δp)
        @test g_v ≈ (-1.5) .* hvp_ref atol = 1e-6 rtol = 1e-5
    end

    # ── buffer reuse: the evaluator is reusable across intervals in place ────
    c1 = zeros(2)
    interval_coeff!(c1, sic, 0.0)
    @test c1[1] ≈ p[1] atol = 1e-14   # τ=0 → the left knot's values
    @test c1[2] ≈ p[1]^2 atol = 1e-14
    interval_coeff!(c1, sic, 1.0)
    @test c1[1] ≈ p[1+u_dim] atol = 1e-14  # τ=1 → the right knot
    @test c1[2] ≈ p[1+u_dim]^2 atol = 1e-14
end

# ── complex_real_interface.jl ──────────────────────────────────────────────── #

@testitem "E1: complex↔real boundary utilities: round-trips and exact blocks" begin
    using Piccolo
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators:
        propagator_to_iso,
        complex_ket_to_iso,
        complex_ket_to_iso!,
        iso_to_complex_ket,
        iso_to_complex_ket!,
        iso_gather_ket!,
        iso_vec_to_complex_operator,
        complex_operator_to_iso_vec!,
        sensitivity_to_jac_col!,
        sensitivity_to_hess_col!,
        hessian_pp_contraction
    using LinearAlgebra
    using Random
    Random.seed!(3476)

    n = 2
    Φ = exp(-im * 0.7 * Matrix{ComplexF64}(PAULIS.Z))
    ψ = randn(ComplexF64, n)

    # ── propagator_to_iso: the 2n×2n block acts on iso kets like Φ on kets ──
    Φ̃ = propagator_to_iso(Φ)
    @test size(Φ̃) == (2n, 2n)
    @test Φ̃[1:n, 1:n] == real(Φ)
    @test Φ̃[1:n, (n+1):2n] == -imag(Φ)
    @test Φ̃[(n+1):2n, 1:n] == imag(Φ)
    @test Φ̃[(n+1):2n, (n+1):2n] == real(Φ)
    @test Φ̃ * complex_ket_to_iso(ψ) ≈ complex_ket_to_iso(Φ * ψ) atol = 1e-14

    # ── ket conversions: allocating and in-place agree; round-trip is identity
    @test complex_ket_to_iso(ψ) == [real(ψ); imag(ψ)]
    ψ̃_buf = zeros(2n)
    complex_ket_to_iso!(ψ̃_buf, ψ)
    @test ψ̃_buf == complex_ket_to_iso(ψ)
    @test iso_to_complex_ket(ψ̃_buf) ≈ ψ atol = 1e-15
    ψ_buf = zeros(ComplexF64, n)
    iso_to_complex_ket!(ψ_buf, ψ̃_buf)
    @test ψ_buf ≈ ψ atol = 1e-15

    # ── iso_gather_ket!: a blocked [Re; Im] slab gather ─────────────────────
    src = randn(5n)
    cols = [3, 7, 1, 5]  # length 2n slab-relative columns
    @test length(cols) == 2n
    ψ_g = zeros(ComplexF64, n)
    iso_gather_ket!(ψ_g, src, cols)
    @test ψ_g ≈ [src[cols[i]] + im * src[cols[n+i]] for i = 1:n] atol = 1e-15

    # ── operator vec conversions round-trip ────────────────────────────────
    Ũ⃗ = operator_to_iso_vec(Φ)
    @test iso_vec_to_complex_operator(Ũ⃗, n) ≈ Φ atol = 1e-14
    Ũ⃗2 = zeros(2n^2)
    complex_operator_to_iso_vec!(Ũ⃗2, Φ)
    @test Ũ⃗2 == Ũ⃗

    # ── Jacobian column helper: col = ket_to_iso(-S_j·ψ) ────────────────────
    S = randn(ComplexF64, n, n)
    col = zeros(2n)
    tmp = zeros(ComplexF64, n)
    sensitivity_to_jac_col!(col, S, ψ, tmp)
    @test col ≈ complex_ket_to_iso(-S * ψ) atol = 1e-14

    # ── Hessian cross-term column: col = ket_to_iso(-S_j'·μ_C) ─────────────
    μ_real = randn(2n)
    tmp1 = zeros(ComplexF64, n)
    tmp2 = zeros(ComplexF64, n)
    hcol = zeros(2n)
    sensitivity_to_hess_col!(hcol, S, μ_real, tmp1, tmp2)
    μ_C = iso_to_complex_ket(μ_real)
    @test hcol ≈ complex_ket_to_iso(-S' * μ_C) atol = 1e-13

    # ── (p,p) contraction: -Re(μ_C' · T · ψ_C) ──────────────────────────────
    T = randn(ComplexF64, n, n)
    ψ_real = randn(2n)
    val = hessian_pp_contraction(T, μ_real, ψ_real, tmp1, tmp2)
    @test val ≈ -real(dot(iso_to_complex_ket(μ_real), T * iso_to_complex_ket(ψ_real))) atol =
        1e-13
end

# ── propagator_result.jl ───────────────────────────────────────────────────── #

@testitem "E1: PropagatorResult views alias the flat storage; malformed inputs throw" begin
    using Piccolo
    using LinearAlgebra

    pdim = 3
    n_params = 5
    pr = PropagatorResult{ComplexF64}(pdim, n_params)

    @test length(pr.Φ_vec) == pdim^2
    @test size(pr.S_mat) == (pdim^2, n_params)

    # The derived views are built ONCE at construction (#358) and ALIAS the
    # flat fields: an in-place write to Φ_vec is visible through Φ_mat.
    pr.Φ_vec .= 1.0
    @test all(isone, pr.Φ_mat)
    pr.S_mat[1] = 2.0
    @test pr.S_3d[1, 1, 1] == 2.0

    # get_propagator / get_sensitivities return those pre-built views
    @test get_propagator(pr, pdim) === pr.Φ_mat
    @test get_sensitivities(pr, pdim) === pr.S_3d
    # get_sensitivities_flat is the NATIVE pdim²×n_params storage, unwrapped
    @test get_sensitivities_flat(pr) === pr.S_mat

    # Two-field construction derives the views; malformed inputs throw loudly
    Φ_vec = fill(1.0 + 0im, pdim^2)
    S_mat = zeros(ComplexF64, pdim^2, n_params)
    pr2 = PropagatorResult(Φ_vec, S_mat)
    @test size(pr2.Φ_mat) == (pdim, pdim)
    @test size(pr2.S_3d) == (pdim, pdim, n_params)

    # Φ_vec must be a perfect square
    @test_throws DimensionMismatch PropagatorResult(fill(1.0 + 0im, pdim^2 + 1), S_mat)
    # S_mat's first dimension must match Φ_vec
    @test_throws DimensionMismatch PropagatorResult(
        Φ_vec,
        zeros(ComplexF64, pdim^2 + 1, n_params),
    )
end

# ── algorithms.jl ──────────────────────────────────────────────────────────── #

@testitem "E1: IntegrationAlgorithm constructor contracts and validation lanes" begin
    using Piccolo
    using LinearAlgebra

    # Tsit5Alg kwarg shape
    a = Tsit5Alg(adaptive = false, tol = 1e-8, ode_h = 0.05)
    @test !a.adaptive
    @test a.tol == 1e-8
    @test a.ode_h == 0.05
    @test Tsit5Alg() isa Tsit5Alg && Tsit5Alg().adaptive

    # MagnusGL4Alg / MagnusAdapt4Alg / Rodas5PAlg kwarg shape
    @test MagnusGL4Alg(n_steps = 20).n_steps == 20
    @test MagnusAdapt4Alg(tol = 1e-10).tol == 1e-10
    @test Rodas5PAlg(tol = 1e-8).tol == 1e-8

    # ── ChebyshevAlg validation lanes (all ArgumentError) ──────────────────
    @test_throws ArgumentError ChebyshevAlg(n_sub = :bogus)
    @test_throws ArgumentError ChebyshevAlg(n_sub = 0)
    @test_throws ArgumentError ChebyshevAlg(sub_dt = -0.1)
    @test_throws ArgumentError ChebyshevAlg(bracket = (2.0, 1.0))
    @test_throws ArgumentError ChebyshevAlg(phase_budget = 0.0)

    # ── ChebyshevAlg normalization: Int → Float64, tuples → floats ─────────
    c = ChebyshevAlg(
        n_sub = 8,
        sub_dt = 2,
        bracket = (-8, 8),
        cheb_tol = 1e-12,
        dyn_tol = 1e-9,
        grad_tol = 1e-7,
        phase_budget = 3,
        tol = 1e-8,
    )
    @test c.n_sub == 8
    @test c.sub_dt == 2.0
    @test c.bracket == (-8.0, 8.0)
    @test c.cheb_tol == 1e-12
    @test c.dyn_tol == 1e-9
    @test c.grad_tol == 1e-7
    @test c.phase_budget == 3.0
    @test c.tol == 1e-8
    # :auto default and nothing-passthrough
    @test ChebyshevAlg().n_sub === :auto
    @test isnothing(ChebyshevAlg().sub_dt)
    @test isnothing(ChebyshevAlg().bracket)
end

# ── alg_data.jl ────────────────────────────────────────────────────────────── #

@testitem "E1: Tsit5Data convenience constructors wire jvp/vjp/hvp fields" begin
    using Piccolo
    using SciMLBase: ODEProblem
    using SparseArrays

    f!(dx, x, p, t) = (dx .= x; nothing)
    probs = [ODEProblem(f!, zeros(4), (0.0, 1.0), zeros(3)) for _ = 1:4]
    structure = sparse(ones(2, 2))

    # 2-field: everything optional is nothing
    d2 = Tsit5Data{ComplexF64}(probs, structure)
    @test length(d2.Φ_probs) == 4
    @test d2.Φ_structure === structure
    @test isnothing(d2.jvp_probs)
    @test isnothing(d2.vjp_probs)
    @test isnothing(d2.hvp_fwd_probs)
    @test isnothing(d2.hvp_bwd_probs)

    # 3-field: jvp_probs materialized, the rest nothing
    d3 = Tsit5Data{ComplexF64}(probs, structure, probs)
    @test !isnothing(d3.jvp_probs)
    @test length(d3.jvp_probs) == 4
    @test isnothing(d3.vjp_probs)

    # 4-field: vjp_probs materialized, hvp pair nothing
    d4 = Tsit5Data{ComplexF64}(probs, structure, probs, probs)
    @test !isnothing(d4.jvp_probs)
    @test !isnothing(d4.vjp_probs)
    @test isnothing(d4.hvp_fwd_probs)
    @test isnothing(d4.hvp_bwd_probs)
end

# ── spline_types.jl ────────────────────────────────────────────────────────── #

@testitem "E1: SplineType trait layer: pulse inference, orders, packed block layout" begin
    using Piccolo
    using LinearAlgebra
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators:
        _drift_matrix,
        spline_order,
        param_blocks,
        n_param_blocks,
        param_block_carries_globals,
        ControlValueBlock,
        ControlDerivBlock,
        ControlDeriv2Block

    # Trait inference from pulse types
    @test SplineType(LinearSplinePulse(zeros(1, 3), 0:0.5:1.0)) == LinearSpline()
    @test SplineType(CubicSplinePulse(zeros(1, 3), zeros(1, 3), 0:0.5:1.0)) == CubicSpline()
    @test LinearSpline() isa SplineType
    @test CubicSpline() isa SplineType

    # Orders
    @test spline_order(LinearSpline()) == 1
    @test spline_order(CubicSpline()) == 3
    @test spline_order(LinearSplinePulse(zeros(1, 3), 0:0.5:1.0)) == 1
    @test spline_order(CubicSplinePulse(zeros(1, 3), zeros(1, 3), 0:0.5:1.0)) == 3

    # Packed ODE-parameter block layouts (#338): linear = [uₖ, uₖ₊₁],
    # cubic = [uₖ, uₖ₊₁, duₖ, duₖ₊₁], each (role, knot_offset)
    pb_lin = param_blocks(LinearSpline())
    @test length(pb_lin) == 2
    @test all(r isa ControlValueBlock for (r, _) in pb_lin)
    @test [off for (_, off) in pb_lin] == [0, 1]
    @test n_param_blocks(LinearSpline()) == 2

    pb_cub = param_blocks(CubicSpline())
    @test length(pb_cub) == 4
    @test [r for (r, _) in pb_cub] == [
        ControlValueBlock(),
        ControlValueBlock(),
        ControlDerivBlock(),
        ControlDerivBlock(),
    ]
    @test [off for (_, off) in pb_cub] == [0, 1, 0, 1]
    @test n_param_blocks(CubicSpline()) == 4

    # Global slots: value blocks carry global VALUES, derivative blocks are zero
    @test param_block_carries_globals(ControlValueBlock())
    @test !param_block_carries_globals(ControlDerivBlock())
    @test !param_block_carries_globals(ControlDeriv2Block())

    # _drift_matrix: matrices pass through, operators materialize
    H = Matrix{ComplexF64}(PAULIS.Z)
    @test _drift_matrix(H) === H
    op = dynamics_operator(LinearDrive(H, 1))
    @test _drift_matrix(op) == H
end

# ── _spline_integrators.jl (family common module) ──────────────────────────── #

@testitem "E1: unitary_rollout_trajectory matches the analytic exponential" begin
    using Piccolo
    using NamedTrajectories
    using LinearAlgebra

    # Constant control: H(t) = 0.3·σz + 0.2·σx — the rollout IS exp(-iHT)
    a = 0.2
    u_fn = t -> [a]
    G(u, t) =
        Piccolo.Isomorphisms.G(0.3 * ComplexF64.(PAULIS.Z) + u[1] * ComplexF64.(PAULIS.X))
    T = 1.3
    samples = 40

    traj = unitary_rollout_trajectory(u_fn, G, T; samples = samples)

    @test traj isa NamedTrajectory
    @test traj.N == samples
    @test size(traj.Ũ⃗, 2) == samples
    @test size(traj.u, 2) == samples
    @test size(traj.u, 1) == 1
    @test all(traj.u .== a)

    U_exact = exp(-im * T * (0.3 * ComplexF64.(PAULIS.Z) + a * ComplexF64.(PAULIS.X)))
    Ũ⃗_roll = traj.Ũ⃗[:, end]
    @test norm(Ũ⃗_roll - operator_to_iso_vec(U_exact)) /
          norm(operator_to_iso_vec(U_exact)) < 1e-10

    # Δt is the uniform grid step (with the final step repeated)
    @test all(dt -> abs(dt - T / (samples - 1)) < 1e-12, diff(vec(collect(traj.t))))
    @test traj.Δt[end] ≈ T / (samples - 1)

    # The initial column is the identity operator's iso vec
    @test traj.Ũ⃗[:, 1] ≈ operator_to_iso_vec(1.0I(2)) atol = 1e-14

    # kwargs: U_goal lands in the goal; control_bounds land in bounds
    U_goal = ComplexF64[0 1; 1 0]
    traj2 = unitary_rollout_trajectory(
        u_fn,
        G,
        T;
        samples = 10,
        U_goal = U_goal,
        control_bounds = (-1.0, 1.0),
    )
    @test traj2.goal[:Ũ⃗] ≈ operator_to_iso_vec(U_goal)
    # NamedTrajectories normalizes Tuple{Real,Real} bounds to ([lo], [hi])
    @test traj2.bounds[:u][1][1] == -1.0
    @test traj2.bounds[:u][2][1] == 1.0
    # Without U_goal, no Ũ⃗ goal is attached
    @test !haskey(traj.goal, :Ũ⃗)
end

# ── E1 resume fill (#347): type-layer utilities, packed-layout traits, and ─── #
# ── the shared ODE-builder contract lanes. ────────────────────────────────── #

@testitem "E1: SplineType trait layer: orders, packed blocks, global carriage" begin
    using Piccolo
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators:
        SplineType, spline_order, param_blocks, n_param_blocks
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators:
        ControlValueBlock, ControlDerivBlock, ControlDeriv2Block
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators:
        param_block_carries_globals

    # The SplineType-value order methods — through the pulse trait chain and directly
    @test spline_order(LinearSpline()) == 1
    @test spline_order(CubicSpline()) == 3

    # Pulse-level inference routes through SplineType(pulse) → spline_order(S)
    times = collect(range(0.0, 1.0, 5))
    lp = LinearSplinePulse(fill(0.5, 1, 5), times)
    cp = CubicSplinePulse(fill(0.5, 1, 5), fill(0.0, 1, 5), times)
    @test spline_order(lp) == 1
    @test spline_order(cp) == 3
    @test SplineType(lp) isa LinearSpline
    @test SplineType(cp) isa CubicSpline

    # Packed block declarations: linear = [uₖ, uₖ₊₁], cubic adds du endpoints
    @test n_param_blocks(LinearSpline()) == 2
    @test n_param_blocks(CubicSpline()) == 4
    @test param_blocks(LinearSpline()) == ((ControlValueBlock(), 0), (ControlValueBlock(), 1))

    # Global carriage: only the control-VALUE blocks ride globals; derivative
    # roles are identically-zero slots.
    @test param_block_carries_globals(ControlValueBlock())
    @test !param_block_carries_globals(ControlDerivBlock())
    @test !param_block_carries_globals(ControlDeriv2Block())
end

@testitem "E1: packed-layout component helpers on live cells (ddu, knot dims)" begin
    using DirectTrajOpt
    using Piccolo
    using NamedTrajectories
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators:
        param_block_component,
        ddu_name,
        du_name,
        canonical_hessian_knot_dim,
        _spline_type,
        spline_order,
        ControlValueBlock,
        ControlDerivBlock,
        ControlDeriv2Block

    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    N = 5
    times = collect(range(0.0, 1.0, N))
    pulse = LinearSplinePulse(fill(0.3, 2, N), times)
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)
    𝒮 = SplineIntegrator(qtraj, N)

    # _spline_type reads the SplineType type-parameter straight off the cell
    @test _spline_type(𝒮) isa LinearSpline
    @test spline_order(𝒮) == 1

    # param_block_component resolves each role to the trajectory component name
    @test param_block_component(𝒮, ControlValueBlock()) == :u
    @test param_block_component(𝒮, ControlDerivBlock()) == du_name(𝒮)
    @test param_block_component(𝒮, ControlDeriv2Block()) == ddu_name(𝒮)
    @test du_name(𝒮) == Symbol("d", 𝒮.u_name)
    @test ddu_name(𝒮) == Symbol("dd", 𝒮.u_name)

    # canonical_hessian_knot_dim: [x, u, Δt, t] for linear (the canonical knot
    # counts ONE u slot), cubic adds the du slot
    @test canonical_hessian_knot_dim(𝒮) == 𝒮.x_dim + 𝒮.u_dim + 2
    𝒮c = SplineIntegrator(
        KetTrajectory(
            sys,
            CubicSplinePulse(fill(0.3, 2, N), fill(0.0, 2, N), times),
            ψ_init,
            ψ_goal,
        ),
        N,
    )
    @test _spline_type(𝒮c) isa CubicSpline
    @test spline_order(𝒮c) == 3
    @test canonical_hessian_knot_dim(𝒮c) == 𝒮c.x_dim + 2 * 𝒮c.u_dim + 2
end

@testitem "E1: canonical_block_hessian_structure covers order/global/ensemble lanes" begin
    using Piccolo
    using LinearAlgebra
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators:
        canonical_block_hessian_structure
    using SparseArrays

    # Single-state, linear: canonical knot = x + u + 2 (ONE u slot), total = 2 knots
    S1 = canonical_block_hessian_structure(4, 2, 1)
    @test size(S1) == (16, 16)
    @test issymmetric(S1)

    # Single-state, cubic: canonical knot = x + 2u + 2
    S3 = canonical_block_hessian_structure(4, 2, 3)
    @test size(S3) == (20, 20)

    # Globals extend the structure by global_dim on BOTH axes
    S1g = canonical_block_hessian_structure(4, 2, 1, 3)
    @test size(S1g) == (19, 19)
    S3g = canonical_block_hessian_structure(4, 2, 3, 1)
    @test size(S3g) == (21, 21)

    # Structure content, linear single-state (x = 1:4, knot dim 8):
    # p columns at knot k: uₖ 5:6, Δt 7, t 8; at k+1: uₖ₊₁ 13:14.
    # Knot-k state rows couple with ALL p columns (upper-triangle entries —
    # they survive the impl's `sparse(Symmetric(...))`). The knot-k+1 state
    # rows' knot-k-side couplings are filled but sit in the LOWER triangle, so
    # the Symmetric truncation drops them: k+1 states keep only their k+1-side
    # parameter columns. That stays a superset of the physical Hessian — the
    # constraint is linear in x_{k+1}, so no true x_{k+1}/p curvature exists.
    for r in 1:4
        for c in (5, 6, 7, 8, 13, 14)
            @test S1[r, c] == 1.0
        end
    end
    for r in 9:12
        for c in (13, 14)
            @test S1[r, c] == 1.0
        end
        for c in (5, 6, 7, 8)
            @test S1[r, c] == 0.0
        end
    end
    @test nnz(S1[1:4, 1:4]) == 0
    @test nnz(S1[1:4, 9:12]) == 0
    @test nnz(S1[5:8, 5:8]) == 0

    # globals couple to both knots' state rows (upper triangle on both sides:
    # global columns land after both knot blocks)
    for r in (1, 3, 9, 11)
        for c in 17:19
            @test S1g[r, c] == 1.0
        end
    end

    # Multi-state (ensemble) version: two kets of dim 4 share the p columns.
    # knot = (4 + 4) + 2 + 2 = 12; x rows 1:8; uₖ 9:10, Δt 11, t 12, uₖ₊₁ 21:22.
    Sm = canonical_block_hessian_structure([4, 4], 2, 1)
    @test size(Sm) == (24, 24)
    for r in 1:8
        for c in (9, 10, 11, 12, 21, 22)
            @test Sm[r, c] == 1.0
        end
    end
    # k+1 state rows (13:20) keep only their k+1-side parameter columns
    for r in 13:20
        for c in (21, 22)
            @test Sm[r, c] == 1.0
        end
        for c in (9, 10)
            @test Sm[r, c] == 0.0
        end
    end

    # Unsupported order: the impl has no else-branch, so knot_dim is unbound
    @test_throws Exception canonical_block_hessian_structure(4, 2, 2)
end

@testitem "E1: get_param_indices two-arg form and refresh_sensitivities! lanes" begin
    using DirectTrajOpt
    using Piccolo
    using NamedTrajectories
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators:
        get_param_indices, refresh_sensitivities!, _refresh_prop_results!

    sys = QuantumSystem(GATES.Z, [GATES.X, GATES.Y], [1.0, 1.0])
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    N = 5
    times = collect(range(0.0, 1.0, N))
    pulse = LinearSplinePulse(fill(0.3, 2, N), times)
    qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)
    traj = NamedTrajectory(qtraj, N)
    𝒮 = SplineIntegrator(qtraj, N)

    # Two-arg form: the per-knot parameter columns plus the Δt/t slots
    ctrl_indices, Δt_idx, t_idx, global_indices = get_param_indices(𝒮, traj)
    @test length(ctrl_indices) == 2 * 𝒮.u_dim  # [uₖ; uₖ₊₁] (linear)
    @test Δt_idx == 2 * 𝒮.u_dim + 1
    @test t_idx == 2 * 𝒮.u_dim + 2
    @test isempty(global_indices)

    # The per-knot form: traj columns [u at k; u at k+1; Δt; t] with matching
    # packed-vector ode indices — the Δt/t slots land at exactly the two-arg
    # form's Δt_idx/t_idx.
    cols_k1, oi_k1 = get_param_indices(𝒮, traj, 1)
    @test length(cols_k1) == 2 * 𝒮.u_dim + 2
    @test length(oi_k1) == 2 * 𝒮.u_dim + 2
    @test oi_k1[end-1:end] == [Δt_idx, t_idx]

    # _refresh_prop_results! runs the threaded compute_ode_jacobian! loop over
    # every knot; refresh_sensitivities! is the public seam over it. The
    # forward residual is unchanged by a sensitivity refresh (the refresh never
    # perturbs the forward state).
    δ = zeros(𝒮.dim)
    evaluate!(δ, 𝒮, traj)
    δ_before = copy(δ)

    _refresh_prop_results!(𝒮, traj, nothing)
    refresh_sensitivities!(𝒮, traj, nothing)

    δ_after = zeros(𝒮.dim)
    evaluate!(δ_after, 𝒮, traj)
    @test δ_after ≈ δ_before atol = 1e-12
    # After the refresh, every knot's prop_results carries the solved Φ
    for k = 1:(traj.N-1)
        @test !iszero(𝒮.prop_results[k].Φ_vec)
    end
end

@testitem "E1: shared ODE-builder contracts: drive-free systems and bad orders" begin
    using Piccolo
    using LinearAlgebra
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators:
        build_sensitivity_ode,
        build_ket_jvp_ode,
        build_hvp_forward_ode,
        build_second_order_adjoint_ode,
        build_ket_sensitivity_ode,
        build_second_order_sensitivity_ode,
        build_density_sensitivity_ode,
        tri_idx,
        _to_operator

    drift_op = _to_operator(Matrix{ComplexF64}(1.0I, 2, 2))

    # Drive-free systems: every builder returns (nothing, n_params) — the
    # documented "no analytic sensitivities without explicit drives" contract.
    n_params_1 = 2 * 1 + 2

    f_sens, np = build_sensitivity_ode(drift_op, AbstractDrive[], 1, 2, 1)
    @test isnothing(f_sens)
    @test np == n_params_1

    f_jvp, np = build_ket_jvp_ode(drift_op, AbstractDrive[], 1, 2, 1)
    @test isnothing(f_jvp)
    @test np == n_params_1

    f_hvp, np = build_hvp_forward_ode(drift_op, AbstractDrive[], 1, 2, 1)
    @test isnothing(f_hvp)
    @test np == n_params_1

    f_2adj, np = build_second_order_adjoint_ode(drift_op, AbstractDrive[], 1, 2, 1)
    @test isnothing(f_2adj)
    @test np == n_params_1

    f_ksens, np = build_ket_sensitivity_ode(drift_op, AbstractDrive[], 1, 2, 1, 2)
    @test isnothing(f_ksens)
    @test np == n_params_1

    f_hess, np, pairs, statedim =
        build_second_order_sensitivity_ode(drift_op, AbstractDrive[], 1, 2, 1)
    @test isnothing(f_hess)
    @test np == n_params_1
    @test isempty(pairs)
    @test statedim == 0

    # The density sensitivity builder shares the same contract (𝒢c real form)
    d_sens, np =
        build_density_sensitivity_ode(Matrix{Float64}(1.0I, 4, 4), Matrix{Float64}[], AbstractDrive[], Vector{Int}[], 1, 2, 1)
    @test isnothing(d_sens)
    @test np == n_params_1

    # Unsupported spline orders error loudly in every builder
    for order in (0, 2, 4)
        @test_throws ErrorException build_sensitivity_ode(drift_op, AbstractDrive[], 1, 2, order)
        @test_throws ErrorException build_ket_jvp_ode(drift_op, AbstractDrive[], 1, 2, order)
        @test_throws ErrorException build_hvp_forward_ode(drift_op, AbstractDrive[], 1, 2, order)
        @test_throws ErrorException build_second_order_adjoint_ode(drift_op, AbstractDrive[], 1, 2, order)
        @test_throws ErrorException build_ket_sensitivity_ode(drift_op, AbstractDrive[], 1, 2, order, 1)
        @test_throws ErrorException build_second_order_sensitivity_ode(drift_op, AbstractDrive[], 1, 2, order)
    end

    # tri_idx: the strictly-upper-triangular pair enumeration index — dense,
    # gap-free, in (i ≤ j) scan order.
    n = 4
    seen = Int[]
    for i = 1:n, j = i:n
        push!(seen, tri_idx(i, j, n))
    end
    @test seen == collect(1:(n * (n + 1) ÷ 2))
end

@testitem "E1: _solve_forward_tsit5 keeps non-Tsit5 algs on the adaptive solve" begin
    using Piccolo
    using SciMLBase: ODEProblem
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators: _solve_forward_tsit5

    # A unit-rate problem whose solution at t=1 is exp(1)·x₀
    prob = ODEProblem((dx, x, p, t) -> (dx[1] = x[1]; nothing), [1.0], (0.0, 1.0))

    # Non-Tsit5Alg dispatch (the Rodas5P-sensitivity fallback): adaptive solve
    sol_other = _solve_forward_tsit5(prob, MagnusGL4Alg(), 1e-8)
    @test length(sol_other.u) == 2          # saveat = 1.0 → [t=0, t=1]
    @test sol_other.u[end][1] ≈ exp(1.0) atol = 1e-6

    # Tsit5Alg fixed-step dispatch honors ode_h (#180)
    sol_fixed = _solve_forward_tsit5(prob, Tsit5Alg(adaptive = false, ode_h = 0.5), 1e-8)
    @test sol_fixed.u[end][1] ≈ exp(1.0) atol = 1e-3
end
