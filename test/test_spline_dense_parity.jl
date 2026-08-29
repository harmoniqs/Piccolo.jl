# ============================================================================ #
# Slice 3b — DENSE parity harness (OB-2, harmoniqs/Piccolissimo.jl#430, AC 1)
#
# The dense leg's population is frozen in test/parity_manifests/slice_3b_dense.toml.
#
# Oracle (3b pre-work, Piccolo #172 resolved by PR #318): for piecewise-constant
# controls, fixed-step MagnusGL4 is grid-EXACT up to roundoff on ANY grid — each
# sub-step sees a constant Hamiltonian, the Magnus commutators vanish (measured:
# d(coarse, n_save=11) = 5.4e-13 vs the n_save=10001 reference, while the adaptive
# residual is tolerance-driven and NOT the tighter of the two). So
# `UnitaryTrajectory(algorithm=MagnusGL4(), n_save=11)` on a PWC fixture IS the
# essentially-exact reference for every dense row.
#
# PWC-ALIGNMENT NOTE (measured while wiring this harness): a spline cell reads
# its control from the trajectory's knot values through ITS OWN spline basis —
# a LinearSpline cell interpolates piecewise-LINEARLY between adjacent knots —
# so the cell's control equals a zero-order hold only when the knot data is
# PWC-aligned (constant across the grid). The oracle rows therefore drive a
# CONSTANT control (the non-trivial PWC-aligned fixture: full drift+drive
# machinery, exact oracle); arbitrary knot values would compare a piecewise-
# linear control against a zero-order-hold reference and fail at O(Δt·Δu) —
# a fixture mismatch, not an integrator defect. Coverage under arbitrary knot
# data is provided by the CROSS-ALGORITHM CONSISTENCY block below: all dense
# algorithms share the cell's spline interpretation, so their mutual agreement
# at tight tolerance gates the dense leg on non-trivial controls.
#
# Piccolissimo-free by construction: every row constructs and runs WITHOUT the
# proprietary package (the Ket Magnus/Chebyshev cells are matrix-free and stay
# in Piccolissimo — a Piccolo-only manifest row would be unconstructible).
#
# RED→GREEN witness for the split: cannot pass before SplineIntegrator + its
# dense cells live in Piccolo.
# ============================================================================ #

@testitem "slice 3b dense parity: manifest rows against the GL4-on-PWC oracle" begin
    using LinearAlgebra
    using Random
    using TOML
    using DirectTrajOpt
    using NamedTrajectories
    using Piccolo
    using OrdinaryDiffEqLinear: MagnusGL4
    using Piccolo.Control.QuantumIntegrators.SplineIntegrators:
        SplineIntegrator, Tsit5Alg, MagnusGL4Alg, MagnusAdapt4Alg, Rodas5PAlg

    manifest_path = joinpath(@__DIR__, "parity_manifests", "slice_3b_dense.toml")
    @test isfile(manifest_path)
    manifest = TOML.parsefile(manifest_path)
    @test haskey(manifest, "rows") && !isempty(manifest["rows"])

    Random.seed!(430)
    T = 1.0
    N = 9
    times = collect(range(0.0, T, length = N))
    n_ctrl = 2

    sys = QuantumSystem(1.7 * PAULIS.Z, [PAULIS.X, PAULIS.Y], [1.0, 1.0])
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]

    alg_of = Dict(
        "Tsit5Alg" => () -> Tsit5Alg(adaptive = true, tol = 1e-12),
        "MagnusGL4Alg" => () -> MagnusGL4Alg(),
        "MagnusAdapt4Alg" => () -> MagnusAdapt4Alg(tol = 1e-12),
        "Rodas5PAlg" => () -> Rodas5PAlg(),
    )

    # Rollout THROUGH a cell: the call operator computes δ = xₖ₊₁ - f(xₖ), so
    # with the candidate xₖ₊₁ zeroed, f(xₖ) = -δ. Walk the knots forward.
    function cell_rollout!(traj, 𝒮, xname, x₀)
        comps = collect(traj.components[xname])
        for (i, c) in enumerate(comps)
            traj.data[c, 1] = x₀[i]
        end
        for k = 1:(traj.N-1)
            δ = zeros(𝒮.x_dim)
            for c in comps
                traj.data[c, k+1] = 0.0
            end
            𝒮(δ, traj[k], traj[k+1], k, nothing)
            fk = .-δ
            for (i, c) in enumerate(comps)
                traj.data[c, k+1] = fk[i]
            end
        end
        return traj.data[comps, traj.N]
    end

    # ── Oracle rows: PWC-aligned (constant) control, GL4 reference (#318) ──── #
    pwc_vals = fill(0.6, n_ctrl, N)
    pwc_pulse = ZeroOrderPulse(pwc_vals, times)

    U_ref = UnitaryTrajectory(
        sys,
        pwc_pulse,
        ComplexF64[0 1; 1 0];
        algorithm = MagnusGL4(),
        n_save = 11,
    ).solution.u[end]
    U_ref_iso = operator_to_iso_vec(U_ref)
    ψ_ref_iso = ket_to_iso(U_ref * ψ_init)

    for row in manifest["rows"]
        kind = row["trajectory"]
        alg_name = row["algorithm"]
        @test row["oracle"] == "gl4_pwc"
        if kind == "unitary"
            qtraj = UnitaryTrajectory(sys, pwc_pulse, ComplexF64[0 1; 1 0])
            traj = NamedTrajectory(qtraj, N)
            x₀ = traj.data[collect(traj.components[:Ũ⃗]), 1]
            𝒮 = SplineIntegrator(qtraj, N; alg = alg_of[alg_name]())
            xN = cell_rollout!(traj, 𝒮, :Ũ⃗, x₀)
            err = norm(xN - U_ref_iso) / norm(U_ref_iso)
            @test err < 1e-8
        elseif kind == "ket"
            qtraj = KetTrajectory(sys, pwc_pulse, ψ_init, ψ_goal)
            traj = NamedTrajectory(qtraj, N)
            x₀ = traj.data[collect(traj.components[:ψ̃]), 1]
            𝒮 = SplineIntegrator(qtraj, N; alg = alg_of[alg_name]())
            xN = cell_rollout!(traj, 𝒮, :ψ̃, x₀)
            err = norm(xN - ψ_ref_iso) / norm(ψ_ref_iso)
            @test err < 1e-8
        elseif kind == "multiket"
            kets = [ψ_init, ComplexF64[0.0, 1.0]]
            goals = [ψ_goal, ψ_init]
            qtraj = MultiKetTrajectory(sys, pwc_pulse, kets, goals)
            traj = NamedTrajectory(qtraj, N)
            𝒮 = SplineIntegrator(qtraj, N; alg = alg_of[alg_name]())
            # The MultiKet call operator fills δ for ALL ket blocks at once —
            # roll the whole state, then compare each ket's final column.
            all_comps = vcat([collect(traj.components[n]) for n in 𝒮.x_names]...)
            for k = 1:(N-1)
                δ = zeros(𝒮.x_dim)
                for c in all_comps
                    traj.data[c, k+1] = 0.0
                end
                𝒮(δ, traj[k], traj[k+1], k, nothing)
                for (i, c) in enumerate(all_comps)
                    traj.data[c, k+1] = -δ[i]
                end
            end
            for (j, name) in enumerate(𝒮.x_names)
                xN = traj.data[collect(traj.components[name]), N]
                ref = ket_to_iso(U_ref * kets[j])
                err = norm(xN - ref) / norm(ref)
                @test err < 1e-8
            end
        elseif kind == "density"
            # Closed-system density fixture: the oracle is the GL4 propagator
            # conjugating ρ₀ — exact for PWC controls (#318). Dissipator-
            # generator parity is covered by the density cells' own testitems
            # in Piccolissimo (unchanged, AC 2).
            osys = OpenQuantumSystem(1.7 * PAULIS.Z, [PAULIS.X, PAULIS.Y], [1.0, 1.0])
            ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]
            ρg = ComplexF64[0.0 0.0; 0.0 1.0]
            qtraj = DensityTrajectory(osys, pwc_pulse, ρ0, ρg)
            traj = NamedTrajectory(qtraj, N)
            x₀ = traj.data[collect(traj.components[:ρ⃗̃]), 1]
            𝒮 = SplineIntegrator(qtraj, N; alg = alg_of[alg_name]())
            xN = cell_rollout!(traj, 𝒮, :ρ⃗̃, x₀)
            ρ_ref_iso = density_to_compact_iso(U_ref * ρ0 * U_ref')
            err = norm(xN - ρ_ref_iso) / norm(ρ_ref_iso)
            @test err < 1e-8
        else
            error("fixture kind not wired in this harness: $kind")
        end
    end

    # ── Cross-algorithm consistency: arbitrary (random) knot data ─────────── #
    # All dense algorithms share the cell's spline interpretation of the knot
    # values, so on a non-PWC fixture their rollouts must STILL agree with
    # each other at tight tolerance (each is a different ODE integration of
    # the same piecewise-linear control problem).
    # (Unitary cells only: the Ket Magnus/Chebyshev cells are matrix-free and
    # stay in Piccolissimo — they are covered by the MF leg's manifest.)
    rand_vals = 0.6 .* randn(n_ctrl, N)
    rand_pulse = ZeroOrderPulse(rand_vals, times)
    rand_refs = Dict{String,Vector{Float64}}()
    for alg_name in ("Tsit5Alg", "MagnusGL4Alg", "MagnusAdapt4Alg")
        qtraj = UnitaryTrajectory(sys, rand_pulse, ComplexF64[0 1; 1 0])
        traj = NamedTrajectory(qtraj, N)
        x₀ = traj.data[collect(traj.components[:Ũ⃗]), 1]
        𝒮 = SplineIntegrator(qtraj, N; alg = alg_of[alg_name]())
        rand_refs[alg_name] = cell_rollout!(traj, 𝒮, :Ũ⃗, x₀)
    end
    ref_tsit = rand_refs["Tsit5Alg"]
    # Tolerance calibrated, not assumed: on a piecewise-linear control the
    # algorithms' own truncation errors disagree at ~1e-3 for Δt = 1/8
    # (measured 7.5e-4, seeded); a wiring/routing regression moves this to
    # O(1). The oracle rows above carry the essentially-exact bound.
    for alg_name in ("MagnusGL4Alg", "MagnusAdapt4Alg")
        err = norm(rand_refs[alg_name] - ref_tsit) / norm(ref_tsit)
        @test err < 2e-3
    end
end
