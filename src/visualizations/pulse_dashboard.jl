using Piccolo
using GLMakie
using LinearAlgebra

# Force GLMakie to open a window
GLMakie.activate!()

# 1. Define the Quantum System
H_drift = zeros(ComplexF64, 2, 2)
H_control = [PAULIS[:X]] 
params = Float64[] 
system = QuantumSystem(H_drift, H_control, params)

# 2. Setup initial template Trajectory
ψ_init = ComplexF64[1.0, 0.0]
ψ_goal = ComplexF64[0.0, 1.0]
duration = 10.0
initial_pulse = GaussianPulse([1.5], 1.5, duration; center=5.0)
base_qtraj = KetTrajectory(system, initial_pulse, ψ_init, ψ_goal)

# 3. Setup GLMakie Observables with concrete types 
amp = Observable{Float64}(1.5)
sigma = Observable{Float64}(1.5)
center = Observable{Float64}(5.0)

# 4. Reactive Simulation Logic
traj_obs = lift(amp, sigma, center) do a, s, c
    new_pulse = GaussianPulse([a], s, duration; center=c)
    # Re-run the physics simulation
    return rollout(base_qtraj, new_pulse)
end

# 5. Create the Dashboard UI
fig = Figure(size = (1200, 800), fontsize = 20)

lsgrid = SliderGrid(
    fig[1, 1],
    (label = "Amplitude (A)", range = 0.1:0.1:5.0, startvalue = 1.5),
    (label = "Width (σ)", range = 0.5:0.1:4.0, startvalue = 1.5),
    (label = "Center (t₀)", range = 1.0:0.1:9.0, startvalue = 5.0),
    width = 350
)

connect!(amp, lsgrid.sliders[1].value)
connect!(sigma, lsgrid.sliders[2].value)
connect!(center, lsgrid.sliders[3].value)

# 6. Plotting the Live Data
ax_pulse = Axis(fig[1, 2], title="Control Pulse: u(t)", xlabel="Time", ylabel="Amplitude")
ax_pop = Axis(fig[2, 2], title="State Evolution: P(|1⟩)", xlabel="Time", ylabel="Probability")
ylims!(ax_pop, -0.05, 1.05)

# FIX: Access the solution field directly
# 't.solution' contains the states. We use 't.solution.t' for the time grid
# and abs2.(t.solution[2, :]) for the population of state |1⟩
times = lift(t -> t.solution.t, traj_obs)
u_vals = lift(t -> sample(t.pulse, t.solution.t)[1, :], traj_obs)
pop_vals = lift(t -> abs2.(t.solution[2, :]), traj_obs) 

lines!(ax_pulse, times, u_vals, color=:blue, linewidth=4)
lines!(ax_pop, times, pop_vals, color=:red, linewidth=4)

display(fig)

# Added functions for lifecycle management
println("Interactive dashboard running. Close the window or press Ctrl+C in terminal to stop.")

# Keeps the window open until the user manually closes it or kills the process
wait(Condition())