# ```@copybutton
# literate/two_qubit_gate_validation.jl
# ```
#
# # [Two-Qubit Gate Validation](@id two-qubit-gate-validation)
#
# This tutorial walks through synthesizing a CNOT gate (Controlled NOT or Controlled X gate) on
# two coupled transmon qubits using three Pulse types with varying degrees of smoothness (numbers
# continuous derivatives).
#
# Furthermore, we validate the fidelities reported by Piccolo.jl when using
# those pulses by rolling out the same pulses in
# [QuantumToolbox.jl](https://github.com/qutip/QuantumToolbox.jl).
#
# ## Setup

using Piccolo
using QuantumToolbox
using CairoMakie
using LinearAlgebra
using Printf
using Random
Random.seed!(1234)  # For reproducibility

# !!! note "A naming clash"
#     Both Piccolo and QuantumToolbox export `fidelity`. With both packages in
#     scope we write `Piccolo.fidelity` for Piccolo's trajectory fidelity. The
#     gate fidelity from the QuantumToolbox rollout we compute by hand below, and
#     for the honest continuous-time Piccolo fidelity we use
#     `unitary_rollout_fidelity` (whose name does not clash).
#
# ## Step 1: Define the Two-Qubit System
#
# The system Hamiltonian for ``N`` transmon qubits with pairwise dipole coupling in
# the rotating frame is
# ```math
# H_{\textrm{drift}} = \sum_{i=1}^{N} \left[
#     (\omega_i - \omega_{\text{frame}})\, a_i^\dagger a_i
#     - \frac{\delta_i}{2}\, a_i^{\dagger 2} a_i^2
# \right]
# + \sum_{i < j} g_{ij} \left( a_i\, a_j^\dagger + a_i^\dagger\, a_j \right),
# ```
# and the control Hamiltonians are
# ```math
# H_{\textrm{control}}(t) = \sum_{i=1}^{N} \left[ u_{x,i}(t)(a_i + a_i^\dagger) + u_y(t)i(a_i^\dagger - a_i) \right]
# ```
# where ``g_{ij}`` determines the coupling strength between qubits ``i`` and ``j``.
#
# We could build the drift and control Hamiltonians ourselves and
# contruct the quantum system directly, but Piccolo.jl
# provides the `MultiTransmonSystem` to easily construct the system.

ωs = [4.0, 4.1]            # transmon frequencies (GHz)
δs = [0.2, 0.2]            # anharmonicities (GHz) — unused at 2 levels, kept for realism
g = 0.1                    # exchange coupling (GHz) — artificially large (see note)
gs = [0.0 g; g 0.0]

sys = MultiTransmonSystem(ωs, δs, gs; levels_per_transmon = 2, drive_bounds = 0.1)

# 
# !!! note "Accuracy of model"
#     We ignore two [best practices](@ref transmon-best-practices) to reduce the
#     computational cost of the gate synthesis in this tutorial.
#
#     1. We only model two levels of each qubit. This reduces the size of
#        the Hilbert space, which reduces the number of 
#        [decision variables](@ref decision-variables) in the NLP optimization.
#     2. The coupling between the qubits is artificially high, which
#        greatly reduces the pulse duration needed to synthesize the gate with high
#        fidelity.
#
#     When synthesizing a gate for real hardware, at least 3 levels should be
#     modeled for each qubit in order to suppress leakege to excited states
#     outside of the computational subspace. Additionally, ``g`` will by in the
#     1-10 MHz range, which will require a much longer gate duration.
#
# ## Step 2: Define the Gate
#
# Our target is the CNOT gate. Piccolo already defines this gate in `GATES[:CX]`.
# Because we only model two levels per qubit, we could use `GATES[:CX]` as our
# gate directly, but we will use an `EmbeddedOperator` to place the gate in
# the system's Hilbert space and record the indices of the computational
# subspace, so that the tutorial is still applicable when `levels_per_transmon` 
# ``\geq 2``.

U_goal = EmbeddedOperator(GATES[:CX], sys)

# The indices of computational subspace basis states ``|00\rangle, |01\rangle, |10\rangle,`` and ``|11\rangle`` are:

U_goal.subspace

# We set an initial gate duration of ``T = 50`` ns. The timestep ``\Delta t`` is
# left as a free optimization variable, so the optimizer may shorten the gate.
# Each warm-started stage below inherits the previous stage's optimized duration
# (and its optimized control parameters), so the three solves share a starting
# point rather than a fixed duration.

T = 50.0    # initial gate duration (ns)

# Every parameterization uses the same number of control parameters, so the only
# thing that changes across the three solves is the pulse *type*.

N_params = 100

# ## Step 3: Synthesize the gates
#
# We synthesize the gate using three different types of pulses:
#
# | Pulse type          | Problem template     | Continuity |
# |:--------------------|:---------------------|------------|
# | `ZeroOrderPulse`    | `SmoothPulseProblem` | ``C^{-1}`` |
# | `LinearSplinePulse` | `SplinePulseProblem` | ``C^{0}``  |
# | `CubicSplinePulse`  | `SplinePulseProblem` | ``C^{1}``  |
#
# The `ZeroOrderPulse` type implements piecewise-constant pulses. Although
# piecewise-constant pulses are unrealistic for hardware platforms which
# require smooth pulses, they are highly computationally efficient because the
# dynamics constraints in the gate synthesis optimization are [**exact**](@ref
# discretization) across each constant interval ``[t_{k}, t_{k+1}]``, and
# computed efficiently via Krylov subspace methods.
#
# For these reasons, we first optimize the `ZeroOrderPulse`. Then we convert the
# optimized pulse to the more realistic `LinearSplinePulse` and use that as a
# warm start for another pulse optimization. Finally, we convert the result of
# that optimization to a `CubicSplinePulse` and use that as a warm start for
# our final pulse optimization.
# 
#
# !!! note "Smoothness of `ZeroOrderPulse`"
#     `ZeroOrderPulse` is not smooth in the mathematical sense, but placing a
#     penalty on the magnitude of the amplitiude change between constant regions
#     can reduce the size of the discontinuities.
#
# ### Piecewise-Constant Pulse
#
# We start from a random piecewise-constant pulse with `N_params` amplitudes.
# Because the `ZeroOrderPulse` dynamics are exact on each constant interval, the
# number of amplitudes is also the number of optimization timesteps.

times_zoh = collect(range(0, T, length = N_params))

pulse_zoh = ZeroOrderPulse(0.02 * randn(sys.n_drives, N_params), times_zoh)
qtraj_zoh = UnitaryTrajectory(sys, pulse_zoh, U_goal)
qcp_zoh = SmoothPulseProblem(
    qtraj_zoh,
    N_params;
    Q = 100.0,
    R = 1e-2,
    ddu_bound = 1.0,
    piccolo_options = PiccoloOptions(timesteps_all_equal = true),
)

cached_solve!(qcp_zoh, "two_qubit_zoh"; max_iter = 4)

#-

Piccolo.fidelity(qcp_zoh), honest_fidelity(qcp_zoh, sys)

# We visualize the controls:

plot_pulse(qcp_zoh; title = "Optimized ZeroOrderPulse Controls")

# ### Linear Spline
#
#
#
# For the warm start we hand the optimized `ZeroOrderPulse` *parameters* straight
# to the `LinearSplinePulse`. Building it from the solved trajectory reuses the
# optimized control values as the spline knots, at the optimized knot times —
# there is no resampling, and the linear spline inherits the ZOH solve's
# duration.
zoh_traj = get_trajectory(qcp_zoh)
pulse_lin = LinearSplinePulse(zoh_traj)
qtraj_lin = UnitaryTrajectory(sys, pulse_lin, U_goal)

plot_pulse(qtraj_lin; title = "Warm Start LinearSplinePulse Controls")

# We now perform the optimization. Omitting the timestep count keeps the native
# knot grid, so the problem has the same `N_params` control parameters as the
# `ZeroOrderPulse` solve.

qcp_lin = SplinePulseProblem(
    qtraj_lin;
    Q = 100.0,
    piccolo_options = PiccoloOptions(timesteps_all_equal = true),
)

cached_solve!(qcp_lin, "two_qubit_linear_spline"; max_iter = 4)

#-

Piccolo.fidelity(qcp_lin), honest_fidelity(qcp_lin, sys)

# We visualize the pulse:

plot_pulse(qcp_lin; title = "Optimized LinearSplinePulse Controls")

# ### Cubic Spline
#
# Finally, a cubic Hermite spline — control values *and* tangents are optimized,
# giving the smoothest waveform. We seed the knot *values* from the optimized
# linear spline (again passed directly, at its optimized knot times) and start
# every Hermite *tangent* at zero. We deliberately do **not** reuse the linear
# spline's slopes: a linear spline is kinked at its knots, so its one-sided
# derivatives there are not a meaningful starting tangent. Zero tangents give the
# optimizer a clean, well-defined starting point.

lin_traj = get_trajectory(qcp_lin)
pulse_cub = CubicSplinePulse(lin_traj[:u], zero(lin_traj[:u]), get_times(lin_traj))
qtraj_cub = UnitaryTrajectory(sys, pulse_cub, U_goal)

plot_pulse(qtraj_cub; title = "Warm Start CubicSplinePulse Controls")

# We now perform the optimization, again on the native `N_params` knot grid.
qcp_cub = SplinePulseProblem(
    qtraj_cub;
    Q = 100.0,
    piccolo_options = PiccoloOptions(timesteps_all_equal = true),
)

cached_solve!(qcp_cub, "two_qubit_cubic_spline"; max_iter = 4)

#- 

Piccolo.fidelity(qcp_cub), honest_fidelity(qcp_cub, sys)

# We visualize the pulse:

plot_pulse(qcp_cub; title = "Optimized CubicSplinePulse Controls")

# ### Optimized Pulses
#
# Finally, we place the three optimized pulses side by side — one panel per
# parameterization. We reuse `plot_pulse!`, which draws a pulse onto an existing
# axis with rendering tailored to its type (stairs for the piecewise-constant
# `ZeroOrderPulse`, line segments for the kinked `LinearSplinePulse`, a smooth
# curve for the `CubicSplinePulse`). The y-axes are linked for a fair amplitude
# comparison.

let
    optimized = [
        ("ZeroOrderPulse", get_pulse(qcp_zoh.qtraj)),
        ("LinearSplinePulse", get_pulse(qcp_lin.qtraj)),
        ("CubicSplinePulse", get_pulse(qcp_cub.qtraj)),
    ]
    colors = Makie.wong_colors()[1:sys.n_drives]

    fig = Figure(size = (1200, 360))
    axes = map(enumerate(optimized)) do (col, (name, pulse))
        ax = Axis(
            fig[1, col];
            title = name,
            xlabel = "Time (μs)",
            ylabel = col == 1 ? "Amplitude (GHz)" : "",
        )
        plot_pulse!(ax, pulse; colors, show_knots = false)
        ax
    end
    linkyaxes!(axes...)

    entries = [LineElement(; color = colors[i], linewidth = 2) for i = 1:sys.n_drives]
    Legend(fig[1, length(optimized) + 1], entries, ["Drive $i" for i = 1:sys.n_drives])
    Label(fig[0, :], "Optimized pulses by parameterization"; font = :bold, fontsize = 18)
    fig
end

# ## Step 5: Validate Pulses
#
# Here is the validation machinery. Given the system and an optimized pulse, we
# reconstruct the implemented unitary with QuantumToolbox.
#
# `pulse_value(pulse, i, t)` samples drive `i` of any pulse at arbitrary time
# `t`. Every Piccolo pulse is callable — `pulse(t)` returns the full control
# vector — so this one helper works uniformly across `ZeroOrderPulse`,
# `LinearSplinePulse`, and `CubicSplinePulse`:

pulse_value(pulse, i, t) = pulse(t)[i]

# The rollout builds a time-dependent Hamiltonian as a `QobjEvo`, evolves each
# computational basis state through `sesolve`, and collects the final states
# into the implemented unitary `U_impl`. It also returns the per-state
# evolutions, which we use for the population plot later.

function rollout_qutip(sys, pulse, U_goal; n_save = 200)
    H_drift = sys.H_drift
    H_drives = get_drives(sys)

    ## A QobjEvo bundles a constant term with (operator, coefficient) pairs.
    ## Each coefficient is a function (params, t); we ignore params and read the
    ## pulse. The comprehension gives each closure its own `i`.
    H = QobjEvo((
        Qobj(H_drift),
        ((Qobj(H_drives[i]), (p, t) -> pulse_value(pulse, i, t)) for i = 1:sys.n_drives)...,
    ))

    tlist = range(0, duration(pulse), length = n_save)
    subspace = U_goal.subspace
    d = length(subspace)
    U_impl = zeros(ComplexF64, d, d)
    sols = Vector{Any}(undef, d)
    for (col, idx) in enumerate(subspace)
        ψ0 = Qobj(ComplexF64.(1:sys.levels .== idx))
        sol = sesolve(H, ψ0, tlist; progress_bar = Val(false))
        sols[col] = sol
        U_impl[:, col] = sol.states[end].data[subspace]
    end
    return U_impl, tlist, sols
end

# The gate fidelity is the standard phase-insensitive overlap on the
# computational subspace — the same convention Piccolo uses internally:
#
# ```math
# F = \frac{\left|\operatorname{tr}\!\left(U_\text{goal}^\dagger\, U_\text{impl}\right)\right|^2}{d^2}
# ```

function gate_fidelity_qutip(sys, pulse, U_goal; kwargs...)
    U_impl, _, _ = rollout_qutip(sys, pulse, U_goal; kwargs...)
    sub = U_goal.subspace
    U_target = U_goal.operator[sub, sub]
    d = length(sub)
    return abs2(tr(U_target' * U_impl)) / d^2
end

# ## Step XXX: Which Piccolo fidelity to trust
#
# Piccolo can report a gate fidelity two ways, and for spline pulses they
# *disagree* at the ``10^{-4}`` level — a subtlety worth understanding before we
# compare against QuantumToolbox.
#
# - **`Piccolo.fidelity(qcp)`** is the convenient default. After a solve it
#   reconstructs the pulse from the optimized trajectory and rolls it out. It is
#   spot-on for `ZeroOrderPulse`, but for splines it comes out slightly
#   *optimistic* (see the discussion at the end of the page).
# - **`unitary_rollout_fidelity(traj, sys; interpolation=...)`** rolls
#   the controls out through an adaptive ODE solve using the interpolation that
#   *matches the pulse type*. This is the honest continuous-time number — and the
#   one that agrees with QuantumToolbox.
#
# We pick the interpolation from the pulse type:

piccolo_interpolation(::ZeroOrderPulse)   = :constant
piccolo_interpolation(::LinearSplinePulse) = :linear
piccolo_interpolation(::CubicSplinePulse)  = :cubic

honest_fidelity(qcp, sys) = unitary_rollout_fidelity(
    get_trajectory(qcp),
    sys;
    interpolation = piccolo_interpolation(get_pulse(qcp.qtraj)),
)

# With these utilities, we can check the fidelity of the synthesized gates, as given by QuantumToolbox.

F_zoh_qutip = gate_fidelity_qutip(sys, get_pulse(qcp_zoh.qtraj), U_goal)
F_lin_qutip = gate_fidelity_qutip(sys, get_pulse(qcp_lin.qtraj), U_goal)
F_cub_qutip = gate_fidelity_qutip(sys, get_pulse(qcp_cub.qtraj), U_goal)

F_zoh_qutip, F_lin_qutip, F_cub_qutip

#-

function cache_stats(name)
    data_dir = joinpath(dirname(Base.active_project()), "data")
    path = Piccolo.DocsCache.find_cache_file(name, data_dir)
    log = isnothing(path) ? "" : Piccolo.DocsCache.load_cache(path).output
    m_iter = match(r"Number of Iterations[.]*:\s*(\d+)", log)
    m_secs = match(r"Total seconds in IPOPT\s*=\s*([\d.]+)", log)
    iters = isnothing(m_iter) ? -1 : parse(Int, m_iter.captures[1])
    secs = isnothing(m_secs) ? NaN : parse(Float64, m_secs.captures[1])
    return iters, secs
end

zoh_iters, zoh_secs = cache_stats("two_qubit_zoh")
lin_iters, lin_secs = cache_stats("two_qubit_linear_spline")
cub_iters, cub_secs = cache_stats("two_qubit_cubic_spline")

results = [
    ("ZeroOrderPulse", N_params, zoh_iters, zoh_secs, Piccolo.fidelity(qcp_zoh), F_zoh_qutip),
    ("LinearSplinePulse", N_params, lin_iters, lin_secs, Piccolo.fidelity(qcp_lin), F_lin_qutip),
    ("CubicSplinePulse", N_params, cub_iters, cub_secs, Piccolo.fidelity(qcp_cub), F_cub_qutip),
]

println(@sprintf("%-18s %5s %6s %9s %12s %12s %9s",
    "Pulse type", "N", "iters", "solve(s)", "F Piccolo", "F QuTiP", "|Δ|"))
for (name, n, it, s, fd,fq) in results
    println(@sprintf("%-18s %5d %6d %9.1f %12.7f %12.7f %9.2e",
        name, n, it, s, fd, fq, abs(fd - fq)))
end


# ## Step 8: Population Dynamics
#
# CNOT flips the target qubit when the control qubit is excited:
# ``|10\rangle \to |11\rangle``. We take the QuantumToolbox rollout of the
# ``|10\rangle`` input (reusing the cubic-spline solution) and plot how the
# computational-state populations evolve.

_, tlist, sols = rollout_qutip(sys, get_pulse(qcp_cub.qtraj), U_goal)

basis_labels = ["|00⟩", "|01⟩", "|10⟩", "|11⟩"]
sub = U_goal.subspace
sol_10 = sols[3]   # |10⟩ input
populations = [abs2.(getindex.(getfield.(sol_10.states, :data), j)) for j in sub]

fig_pop = Figure()
ax = Axis(
    fig_pop[1, 1];
    xlabel = "time (ns)",
    ylabel = "population",
    title = "CNOT: evolution of |10⟩ (cubic spline)",
)
for (j, label) in enumerate(basis_labels)
    lines!(ax, collect(tlist), populations[j]; label = label)
end
axislegend(ax; position = :rc)
fig_pop
