"""Calculate Weyl chamber coordinates"""
function c1c2c3(U)
    SySy = kron(PAULIS[:Y], PAULIS[:Y])

    if size(U) != (4, 4)
        throw(DimensionMismatch("Expected a 4x4 matrix"))
    end
    Ũ = SySy * transpose(U) * SySy
    ev = eigvals(U * Ũ / √complex(det(U)))
    two_S = angle.(ev) / π
    for i in eachindex(two_S)
        if two_S[i] <= -0.5
            two_S[i] += 2.0
        end
    end
    S = sort(two_S / 2.0, rev = true)
    n = Int(round(sum(S)))
    S -= vcat(ones(n), zeros(4 - n))
    S = circshift(S, -n)
    M = [1 1 0; 1 0 1; 0 1 1]
    c1, c2, c3 = M * S[1:3]
    if c3 < 0
        c1 = 1 - c1
        c3 = -c3
    end
    return Point3f(c1, c2, c3)
end


"""plot the weyl chamber and record animated trajectory"""
function plot_weyl_trajectory(traj::NamedTrajectory, output_mp4 = "weyl_trajectory.mp4")
    A1 = Point3f(1, 0, 0)
    A2 = Point3f(0.5, 0.5, 0)
    A3 = Point3f(0.5, 0.5, 0.5)
    O = Point3f(0, 0, 0)
    L = Point3f(0.5, 0, 0)
    M = Point3f(0.75, 0.25, 0)
    N = Point3f(0.75, 0.25, 0.25)
    P = Point3f(0.25, 0.25, 0.25)
    Q = Point3f(0.25, 0.25, 0)

    fig = Figure()
    ax = Axis3(
        fig[1, 1],
        xlabel = "c₁ / π",
        ylabel = "c₂ / π",
        zlabel = "c₃ / π",
        xgridvisible = false,
        ygridvisible = false,
        zgridvisible = false,
        azimuth = 1.7pi,
    )
    ax.xspinecolor_2 = :transparent
    ax.xspinecolor_3 = :transparent
    ax.xspinecolor_4 = :transparent
    ax.yspinecolor_2 = :transparent
    ax.yspinecolor_3 = :transparent
    ax.yspinecolor_4 = :transparent
    ax.zspinecolor_2 = :transparent
    ax.zspinecolor_3 = :transparent
    ax.zspinecolor_4 = :transparent

    # background - weyl edges
    lines!(ax, [O, A2], color = :black, linestyle = :dash)

    # background - PE edges
    lines!(ax, [M, L], color = :black, linestyle = :dash)
    lines!(ax, [Q, L], color = :black, linestyle = :dash)
    lines!(ax, [P, Q], color = :black, linestyle = :dash)
    lines!(ax, [P, A2], color = :black, linestyle = :dash)

    # scatter plots
    points_obs = Observable(Point3f[])
    scatter!(ax, points_obs, color = :red)

    # foreground - weyl edges
    lines!(ax, [O, A1], color = :black)
    lines!(ax, [A1, A2], color = :black)
    lines!(ax, [A2, A3], color = :black)
    lines!(ax, [A3, A1], color = :black)
    lines!(ax, [A3, O], color = :black)

    # foreground - PE edges
    lines!(ax, [L, N], color = :black)
    lines!(ax, [L, P], color = :black)
    lines!(ax, [N, P], color = :black)
    lines!(ax, [N, A2], color = :black)
    lines!(ax, [N, M], color = :black)

    traj_points =
        Point3f[c1c2c3(iso_vec_to_operator(traj[:Ũ⃗][:, k])) for k in axes(traj[:Ũ⃗], 2)]

    record(fig, output_mp4, eachindex(traj_points); framerate = 30) do i
        push!(points_obs[], traj_points[i])
        notify(points_obs)
    end

    return fig
end

# ============================================================================ #
# Tests
# ============================================================================ #

@testitem "c1c2c3 maps canonical two-qubit gates to their Weyl chamber classes" begin
    using Piccolo
    using LinearAlgebra

    c1c2c3 = Piccolo.Visualizations.QuantumObjectPlots.c1c2c3

    I4 = Matrix{ComplexF64}(I, 4, 4)
    SWAP = ComplexF64[1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]
    CNOT = ComplexF64[1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]
    iSWAP = ComplexF64[1 0 0 0; 0 0 im 0; 0 im 0 0; 0 0 0 1]
    CZ = ComplexF64[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 -1]

    # Corners of the Weyl chamber.
    @test collect(c1c2c3(I4)) ≈ [0.0, 0.0, 0.0] atol = 1e-6    # identity: O
    @test collect(c1c2c3(SWAP)) ≈ [0.5, 0.5, 0.5] atol = 1e-6  # SWAP: A3
    @test collect(c1c2c3(iSWAP)) ≈ [0.5, 0.5, 0.0] atol = 1e-6 # iSWAP: A1
    @test collect(c1c2c3(CNOT)) ≈ [0.5, 0.0, 0.0] atol = 1e-6  # CNOT: A2
    @test collect(c1c2c3(CZ)) ≈ [0.5, 0.0, 0.0] atol = 1e-6    # CZ is locally equivalent to CNOT

    # A product of single-qubit gates is locally equivalent to the identity.
    @test collect(c1c2c3(kron(PAULIS[:X], PAULIS[:X]))) ≈ [0.0, 0.0, 0.0] atol = 1e-6

    # Local gates conjugating CNOT leave the class invariant.
    L = kron(PAULIS[:X], PAULIS[:Y])
    R = kron(PAULIS[:Z], PAULIS[:X])
    @test collect(c1c2c3(L * CNOT * R)) ≈ collect(c1c2c3(CNOT)) atol = 1e-6

    # The global phase is divided out via the det normalization.
    @test collect(c1c2c3(exp(im * π / 7) * CNOT)) ≈ collect(c1c2c3(CNOT)) atol = 1e-6

    # Only 4×4 matrices are accepted.
    @test_throws DimensionMismatch c1c2c3(Matrix{ComplexF64}(I, 2, 2))
    @test_throws DimensionMismatch c1c2c3(Matrix{ComplexF64}(I, 8, 8))
end

@testitem "plot_weyl_trajectory renders headlessly and records an mp4" begin
    using Piccolo
    using CairoMakie
    using LinearAlgebra

    # A short trajectory: identity → CNOT family → SWAP.
    N = 5
    I4 = Matrix{ComplexF64}(I, 4, 4)
    SWAP = ComplexF64[1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]
    CNOT = ComplexF64[1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]
    Us = [I4, exp(-im * 0.2 * CNOT), exp(-im * 0.4 * CNOT), exp(-im * 0.7 * SWAP), SWAP]

    traj = NamedTrajectory(
        (Ũ⃗ = hcat(operator_to_iso_vec.(Us)...), u = zeros(1, N), Δt = fill(0.1, N));
        controls = :u,
        timestep = :Δt,
    )

    # Explicit output path.
    out = tempname() * ".mp4"
    fig = plot_weyl_trajectory(traj, out)
    @test fig isa Figure
    @test isfile(out)
    @test filesize(out) > 0
    rm(out; force = true)

    # Default filename lands in the working directory.
    fig2 = plot_weyl_trajectory(traj)
    @test fig2 isa Figure
    @test isfile("weyl_trajectory.mp4")
    rm("weyl_trajectory.mp4"; force = true)
end
