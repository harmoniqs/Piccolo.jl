# A/B gate for #309: R_bend default (1e-3) vs opt-out (0) on a cubic-spline gate solve.
# Acceptance: |F_default − F_off| < 0.005 AND bend_default < bend_off.
using Piccolo
using DirectTrajOpt
using NamedTrajectories
using LinearAlgebra
using Random
using Printf

Random.seed!(309)

# Transmon-ish 2-level X gate, cubic Hermite spline, small enough to solve fast.
ω = 2π * 4.0
H_drift = (ω / 2) * PAULIS[:Z]
H_drives = [(1 / 2) * PAULIS[:X]]
sys = QuantumSystem(H_drift, H_drives, [1.0])

T = 10.0
N = 41
times = collect(range(0.0, T, length = N))
amps = 0.05 * randn(1, N)          # cold-ish start
derivs = zeros(1, N)
pulse = CubicSplinePulse(amps, derivs, times)
qtraj = UnitaryTrajectory(sys, pulse, GATES[:X])

function solve_with(R_bend)
    Random.seed!(309)               # identical starts
    qcp = SplinePulseProblem(qtraj, N; Q = 100.0, R_bend = R_bend, integrator_type = :pwc)
    solve!(qcp; max_iter = 200, print_level = 0)
    rep = verify(qcp; verbose = false)   # independent rollout of the extracted pulse
    F = rep.F_rollout
    m = shape_metrics(extract_pulse(qcp.qtraj, qcp.prob.trajectory))
    return (F = F, bend = m.bend[1], traj = qcp.prob.trajectory)
end

@printf("solving R_bend = 0 (opt-out)…\n")
off = solve_with(0)
@printf("solving R_bend = 1e-3 (default)…\n")
on = solve_with(1e-3)

ΔF = on.F - off.F
bend_ratio = on.bend / off.bend

@printf("\nRESULTS\n")
@printf("F_off      = %.6f\n", off.F)
@printf("F_default  = %.6f\n", on.F)
@printf("ΔF         = %+.6f  (gate: |ΔF| < 0.005)\n", ΔF)
@printf("bend_off   = %.6e\n", off.bend)
@printf("bend_def   = %.6e\n", on.bend)
@printf("bend_ratio = %.4f   (gate: < 1.0)\n", bend_ratio)

open("ab_gate_results.txt", "w") do io
    println(io, "F_off = $(off.F)")
    println(io, "F_default = $(on.F)")
    println(io, "delta_F = $(ΔF)")
    println(io, "bend_off = $(off.bend)")
    println(io, "bend_default = $(on.bend)")
    println(io, "bend_ratio = $(bend_ratio)")
end
