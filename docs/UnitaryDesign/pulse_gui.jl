#!/usr/bin/env julia

using Piccolo
using GLMakie
using Observables
using LinearAlgebra

# Force the interactive window
GLMakie.activate!()

# 1. Setup Data & Constants
T, N = 10.0, 200
times = range(0, T, length=N)
U_goal = GATES[:X]
curr_idx = Observable(1)

# 2. UI Structure
fig = Figure(size = (1400, 950), backgroundcolor = :white)
Label(fig[0, :], "Piccolo.jl: Professional Trajectory Dashboard", fontsize = 26, font = :bold)

sg = SliderGrid(fig[1, 1:2],
    (label = "Ω (Amp)", range = 0.0:0.1:5.0, startvalue = 1.0),
    (label = "Δ (Freq)", range = 0.0:0.2:10.0, startvalue = 2.0),
    (label = "Chirp", range = -1.0:0.1:1.0, startvalue = 0.0),
    width = 700
)

# 3. Reactive Physics Engine (Fixed Syntax)
traj_obs = lift(sg.sliders[1].value, sg.sliders[2].value, sg.sliders[3].value) do A, f, χ
    # FIX: Removed the misplaced 'drive_bounds =' inside the array
    sys = QuantumSystem(0.5 * PAULIS[:Z], [PAULIS[:X], PAULIS[:Y]], [5.0, 5.0])
    
    # Construct Pulse
    phase = (f .* times .+ 0.5 .* χ .* times.^2) 
    u = A .* vcat(cos.(phase)', sin.(phase)')
    pulse = ZeroOrderPulse(u, collect(times))
    
    # rollout with n_save=N ensures data synchronization with 'times'
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)
    return rollout(qtraj, n_save=N) 
end

# 4. Visualizations
# Panel A: Pulse Fields
ax_pulse = Axis(fig[2, 1], title = "Control Fields (u_x, u_y)", xlabel = "Time")
u_matrix = lift(traj_obs) do t
    reduce(hcat, [t.pulse(ts) for ts in times])
end
lines!(ax_pulse, times, lift(m -> m[1,:], u_matrix), color = :cyan, label="u_x")
lines!(ax_pulse, times, lift(m -> m[2,:], u_matrix), color = :magenta, label="u_y")
vlines!(ax_pulse, lift(i -> [times[i]], curr_idx), color = :red, linestyle = :dash)
axislegend(ax_pulse)

# Panel B: Unitary Populations
ax_pop = Axis(fig[2, 2], title = "Unitary Populations (|U_{ij}|^2)", xlabel = "Time")
pop_data = lift(traj_obs) do t
    # Using .solution.u to access the raw unitary matrices from the ODE solution
    p1 = [abs2(U[1,1]) for U in t.solution.u]
    p2 = [abs2(U[2,1]) for U in t.solution.u]
    return (p1, p2)
end
lines!(ax_pop, times, lift(d -> d[1], pop_data), color = :blue, label="|0⟩ → |0⟩")
lines!(ax_pop, times, lift(d -> d[2], pop_data), color = :red, label="|0⟩ → |1⟩")
vlines!(ax_pop, lift(i -> [times[i]], curr_idx), color = :black, linestyle = :dash)
axislegend(ax_pop, position = :rc)

# Panel C: Bloch Sphere Trajectory
ax_bloch = Axis3(fig[3, 1], title = "Bloch Sphere Trajectory", aspect = :data)



points = lift(traj_obs) do t
    unitaries = t.solution.u
    # Project unitary evolution onto the Bloch Sphere starting from ground state |0>
    states = [U * [1, 0] for U in unitaries] 
    return [Point3f(real(s' * PAULIS[:X] * s), 
                    real(s' * PAULIS[:Y] * s), 
                    real(s' * PAULIS[:Z] * s)) for s in states]
end
lines!(ax_bloch, points, color = (:blue, 0.4))
scatter!(ax_bloch, lift((p, i) -> [p[i]], points, curr_idx), color = :red, markersize = 12)

# Panel D: Fidelity Metric
ax_fid = Axis(fig[3, 2], title = "Gate Fidelity (Target: X)")
barplot!(ax_fid, [1], lift(t -> [fidelity(t)], traj_obs), color = :orange)
ylims!(ax_fid, 0, 1.1)

# 5. Playback
play_btn = Button(fig[4, 1], label = "Run Playback", width = 150)
on(play_btn.clicks) do _
    @async for i in 1:N
        curr_idx[] = i
        sleep(0.02)
    end
end

# 6. Launch and Persist
display(fig)
println("Interactive dashboard running. Close the window or press Ctrl+C in terminal to stop.")

# Keeps the window open until the user manually closes it or kills the process
wait(Condition())
