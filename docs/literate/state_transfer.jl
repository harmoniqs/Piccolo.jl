# # [State Transfer](@id state-transfer-tutorial)
#
# This tutorial covers state-to-state transfer using `KetTrajectory`.
# We'll prepare quantum states and implement gates via state mappings.
#
# ## Overview
#
# While `UnitaryTrajectory` optimizes for a full unitary gate,
# `KetTrajectory` optimizes for transferring a specific initial state
# to a target state. This is useful for:
# - State preparation
# - Gates where you only care about specific state mappings
# - Problems where full unitary tracking is expensive

using Piccolo
using CairoMakie
using Random
Random.seed!(123)

# ## Single State Transfer
#
# Let's prepare the |1⟩ state starting from |0⟩.

# ### Setup

## Create quantum system
H_drift = 0.5 * PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
sys = QuantumSystem(H_drift, H_drives, [1.0, 1.0])

## Time parameters
T, N = 10.0, 100
times = collect(range(0, T, length=N))

## Initial pulse
pulse = ZeroOrderPulse(0.1 * randn(2, N), times)

# ### Define State Transfer

## Initial state: |0⟩
ψ_init = ComplexF64[1.0, 0.0]

## Target state: |1⟩
ψ_goal = ComplexF64[0.0, 1.0]

## Create KetTrajectory
qtraj = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

# ### Solve

qcp = SmoothPulseProblem(qtraj, N; Q=100.0, R=1e-2)
solve!(qcp; max_iter=20, verbose=false, print_level=1)

println("State transfer fidelity: ", round(fidelity(qcp), digits=6))

# ### Visualize State Evolution

traj = get_trajectory(qcp)

## Convert isomorphic states to physical states
n_steps = size(traj[:ψ̃], 2)
populations = zeros(2, n_steps)

for k in 1:n_steps
    ψ = iso_to_ket(traj[:ψ̃][:, k])
    populations[:, k] = abs2.(ψ)
end

## Plot
fig = Figure(size=(800, 400))

ax1 = Axis(fig[1, 1], xlabel="Timestep", ylabel="Population", title="State Evolution")
lines!(ax1, 1:n_steps, populations[1, :], label="|0⟩", linewidth=2)
lines!(ax1, 1:n_steps, populations[2, :], label="|1⟩", linewidth=2)
axislegend(ax1, position=:rt)

fig

# ## Multi-State Transfer (Gate via States)
#
# We can define a gate by specifying how it maps multiple states.
# `MultiKetTrajectory` optimizes all mappings simultaneously with
# coherent phases.

# ### Define X Gate via State Mappings
#
# The X gate maps:
# - |0⟩ → |1⟩
# - |1⟩ → |0⟩

## Basis states
ψ0 = ComplexF64[1.0, 0.0]  # |0⟩
ψ1 = ComplexF64[0.0, 1.0]  # |1⟩

## Initial states and their targets
initial_states = [ψ0, ψ1]
goal_states = [ψ1, ψ0]  # X gate swaps them

## Create new pulse for this problem
pulse_multi = ZeroOrderPulse(0.1 * randn(2, N), times)

## Create MultiKetTrajectory
qtraj_multi = MultiKetTrajectory(sys, pulse_multi, initial_states, goal_states)

# ### Solve with Coherent Fidelity
#
# `MultiKetTrajectory` uses `CoherentKetInfidelityObjective`, which ensures
# not just that each state reaches its target, but that the relative phases
# between states are preserved correctly.

qcp_multi = SmoothPulseProblem(qtraj_multi, N; Q=100.0, R=1e-2)
solve!(qcp_multi; max_iter=20, verbose=false, print_level=1)

println("Multi-state gate fidelity: ", round(fidelity(qcp_multi), digits=6))

# ### Visualize Both State Evolutions

traj_multi = get_trajectory(qcp_multi)

## Extract populations for both initial states
## MultiKetTrajectory stores states as :ψ̃1, :ψ̃2, etc.
n_steps = size(traj_multi[:ψ̃1], 2)

pops1 = zeros(2, n_steps)  # Evolution from |0⟩
pops2 = zeros(2, n_steps)  # Evolution from |1⟩

for k in 1:n_steps
    ψ1_k = iso_to_ket(traj_multi[:ψ̃1][:, k])
    ψ2_k = iso_to_ket(traj_multi[:ψ̃2][:, k])
    pops1[:, k] = abs2.(ψ1_k)
    pops2[:, k] = abs2.(ψ2_k)
end

fig2 = Figure(size=(800, 400))

ax1 = Axis(fig2[1, 1], xlabel="Timestep", ylabel="Population", title="|0⟩ → |1⟩")
lines!(ax1, 1:n_steps, pops1[1, :], label="|0⟩", linewidth=2, color=:blue)
lines!(ax1, 1:n_steps, pops1[2, :], label="|1⟩", linewidth=2, color=:red)
axislegend(ax1, position=:rt)

ax2 = Axis(fig2[1, 2], xlabel="Timestep", ylabel="Population", title="|1⟩ → |0⟩")
lines!(ax2, 1:n_steps, pops2[1, :], label="|0⟩", linewidth=2, color=:blue)
lines!(ax2, 1:n_steps, pops2[2, :], label="|1⟩", linewidth=2, color=:red)
axislegend(ax2, position=:rt)

fig2

# ## When to Use Each Trajectory Type
#
# | Trajectory | Use Case | Pros | Cons |
# |------------|----------|------|------|
# | `UnitaryTrajectory` | Full gate synthesis | Complete gate, any input | Tracks d² elements |
# | `KetTrajectory` | Single state prep | Fast, simple | Single state only |
# | `MultiKetTrajectory` | Gate via state maps | Phase coherent | Need all relevant states |

# ## Superposition State Preparation
#
# Let's prepare a superposition state: |+⟩ = (|0⟩ + |1⟩)/√2

ψ_plus = ComplexF64[1, 1] / sqrt(2)

pulse_super = ZeroOrderPulse(0.1 * randn(2, N), times)
qtraj_super = KetTrajectory(sys, pulse_super, ψ0, ψ_plus)

qcp_super = SmoothPulseProblem(qtraj_super, N; Q=100.0, R=1e-2)
solve!(qcp_super; max_iter=20, verbose=false, print_level=1)

println("\n|+⟩ state preparation fidelity: ", round(fidelity(qcp_super), digits=6))

# Verify the final state:

traj_super = get_trajectory(qcp_super)
ψ_final = iso_to_ket(traj_super[:ψ̃][:, end])
println("Final state: ", round.(ψ_final, digits=3))
println("Target state: ", round.(ψ_plus, digits=3))

# ## Next Steps
#
# - [Multilevel Transmon](@ref multilevel-transmon-tutorial): Work with realistic superconducting qubits
# - [Robust Control](@ref robust-control-tutorial): Make pulses robust to parameter uncertainty
# - [Problem Templates](@ref problem-templates-overview): Full API documentation
