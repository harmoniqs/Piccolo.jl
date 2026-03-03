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
    S -= vcat(ones(n), zeros(4-n))
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
