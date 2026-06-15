# ```@copybutton
# literate/two_qubit_gate_validation.jl
# ```
#
# # [Two-Qubit Gate Validation](@id two-qubit-gate-validation)
#
# This tutorial walks through synthesizing a CNOT gate (Controlled NOT or
# Controlled X gate) on two coupled transmon qubits using three Pulse types
# with varying degrees of smoothness (numbers of continuous derivatives).
#
# Furthermore, we validate the fidelities reported by Piccolo.jl when using
# those pulses by rolling out the same pulses in
# [QuantumToolbox.jl](https://github.com/qutip/QuantumToolbox.jl).
#
# ## Setup
# First we import the packages we will use throughout this tutorial.

using Piccolo
using QuantumToolbox
using CairoMakie
using LinearAlgebra
using Printf
using Random
Random.seed!(1234)  # For reproducibility

# !!! note "Naming Clash"
#     Both `Piccolo` and `QuantumToolbox` export `fidelity`. With both packages in
#     scope we write `Piccolo.fidelity` for Piccolo's trajectory fidelity.
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
# construct the quantum system directly, but we use Piccolo.jl's
# `MultiTransmonSystem` constructor to easily build the system.

ωs = [4.0, 4.1]            # transmon frequencies (GHz)
δs = [0.2, 0.2]            # anharmonicities (GHz) — unused at 2 levels_per_transmon
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
#     modeled for each qubit in order to suppress leakage to excited states
#     outside of the computational subspace. Additionally, ``g`` will be in the
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

# We set a gate duration of ``T = 10`` ns. The timestep size ``\Delta t`` is
# left as a free optimization variable, so the optimizer may adjust the total
# duration slightly. For each pulse type, we use the same number of control
# parameters (knots).

T = 10.0    # gate duration (ns)
N_params = 200

#
# !!! note "Number of Parameters"
#     With ``200`` parameters, the spacing of the knots is ``\approx 0.05`` ns.
#     This means the width of each constant section of a `ZeroOrderPulse` is
#     ``\approx 0.05`` ns, which may be unrealistic for the arbitrary waveform
#     generators that produce piecewise constant pulses in most quantum
#     hardware.  
#
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

cached_solve!(qcp_zoh, "two_qubit_zoh"; max_iter = 150)

#-

F_zoh_piccolo = Piccolo.fidelity(qcp_zoh)


# ### Linear Spline
#
# Next we optimize a piecewise-linear control pulse, whose control parameters
# are the pulse values at the edges of each linear interval.
#
# We use the results from the previous optimization to set up a warm start.
# `LinearSplinePulse` pulls the control parameters directly from
# `ZeroOrderPulse` to build a similar pulse shape (with the same duration,
# which has been optimized from the intial duration).

zoh_traj = get_trajectory(qcp_zoh)
pulse_lin = LinearSplinePulse(zoh_traj)
qtraj_lin = UnitaryTrajectory(sys, pulse_lin, U_goal)

# We now perform the optimization. 

N_timesteps_lin = N_params
qcp_lin = SplinePulseProblem(
    qtraj_lin,
    N_timesteps_lin;
    Q = 100.0,
    R_du = 0.1,
    piccolo_options = PiccoloOptions(timesteps_all_equal = true),
)

cached_solve!(qcp_lin, "two_qubit_linear_spline"; max_iter = 80)

# !!! note "Number of Timesteps"
#     The discretized dynamics constraints of the optimizer are not exact when
#     the control pulses are not piecewise-constant. Consequently, the number of
#     timesteps may need to be increased, depending the degree of physical
#     accuracy needed. If the discretized dynamics are not accurate, the
#     optimizer may maximize the fidelity for the *discretized* dynamics, while
#     the actual fidelity for the real, *continuous-time* dynamics is not as
#     good. 

#-

F_lin_piccolo = Piccolo.fidelity(qcp_lin)

# ### Cubic Spline
#
# Next we optimize a cubic spline control pulse, whose control parameters
# are control the pulse's value and the slope of the tangent line at several
# points in time.
#
# We use the results from the previous optimization to set up a warm start.
# From the `LinearSplinePulse`, we use the knot values as our pulse values, and
# set the initial slope of the tangent line at those points to be zero. 

lin_traj = get_trajectory(qcp_lin)
pulse_cub = CubicSplinePulse(lin_traj[:u], zero(lin_traj[:u]), get_times(lin_traj))
qtraj_cub = UnitaryTrajectory(sys, pulse_cub, U_goal)

# We now perform the optimization, again noting that the number of timesteps
# may need to be adjusted for greater physical accuracy.

N_timesteps_cub = N_params
qcp_cub = SplinePulseProblem(
    qtraj_cub,
    N_timesteps_cub;
    Q = 100.0,
    R_du = 0.1,
    piccolo_options = PiccoloOptions(timesteps_all_equal = true),
)

cached_solve!(qcp_cub, "two_qubit_cubic_spline"; max_iter = 80)

#- 

F_cub_piccolo = Piccolo.fidelity(qcp_cub)


# ### Visualizing the Optimized Pulses
# 
# Finally, we visualize the optimized pulses side by side. 

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
            xlabel = "Time (ns)",
            ylabel = col == 1 ? "Amplitude (GHz)" : "",
        )
        plot_pulse!(ax, pulse; colors, show_knots = false)
        ax
    end
    linkyaxes!(axes...)

    entries = [LineElement(; color = colors[i], linewidth = 2) for i = 1:sys.n_drives]
    Legend(fig[1, length(optimized)+1], entries, ["Drive $i" for i = 1:sys.n_drives])
    Label(fig[0, :], "Optimized Pulses"; font = :bold, fontsize = 18)
    fig
end

# ## Step 4: Validate Pulses
#
# Because the discretized dynamics used by the optimizer are not exact for the
# `LinearSplinePulse` and `CubicSplinePulse`,  it is important to validate that
# the reported fidelities are accurate. We can do this using
# [QuantumToolbox.jl](https://github.com/qutip/QuantumToolbox.jl).
#
#
# ### Helper Functions
#
# Piccolo ships a `rollout_with_qutip(system, pulse, ψ0)` helper in its
# `PiccoloQuantumToolboxExt` extension, which loads once both `QuantumToolbox`
# and a Makie backend are in scope. Given a a system, a pulse, and an initial
# condition, it allows us to compute solution to the system of ODEs
# (Schrodinger's equation or the Linlad master equation).
#
# The following function uses `rollout_with_qutip` to roll out each initial
# condition from the computational basis states and collect them into the
# columns of `U_impl` (and stores the solutions for later plotting).

function unitary_from_qutip(sys, pulse, U_goal; n_save = 200)
    subspace = U_goal.subspace
    d = length(subspace)
    U_impl = zeros(ComplexF64, d, d)
    sols = Vector{Any}(undef, d)
    for (col, idx) in enumerate(subspace)
        ψ0 = ComplexF64.(1:sys.levels .== idx)
        sol = rollout_with_qutip(sys, pulse, ψ0; n_save)
        sols[col] = sol
        U_impl[:, col] = sol.states[end].data[subspace]
    end
    return U_impl, sols
end

# We finally make a function for computing the fidelity using the
# QuantumToolbox rollout.
#
# ```math
# F = \frac{\left|\operatorname{tr}\!\left(U_\text{goal}^\dagger\, U_\text{impl}\right)\right|^2}{d^2}
# ```

function gate_fidelity_qutip(sys, pulse, U_goal; kwargs...)
    U_impl, _ = unitary_from_qutip(sys, pulse, U_goal; kwargs...)
    sub = U_goal.subspace
    U_target = U_goal.operator[sub, sub]
    d = length(sub)
    return abs2(tr(U_target' * U_impl)) / d^2
end


# ### Compute Fidelities with Quantum Toolbox
#
# With our helper functions in plase, we can compute the fidelities using
# QuantumToolbox.

F_zoh_qutip = gate_fidelity_qutip(sys, get_pulse(qcp_zoh.qtraj), U_goal)
F_lin_qutip = gate_fidelity_qutip(sys, get_pulse(qcp_lin.qtraj), U_goal)
F_cub_qutip = gate_fidelity_qutip(sys, get_pulse(qcp_cub.qtraj), U_goal)

F_zoh_qutip, F_lin_qutip, F_cub_qutip

# We need a helper function to get the number of iterations and wall time.
#
# !!! note "TODO"
#     I am currently grabbing the iterations and wall time by REGEXing the IPOPT
#     output generated by the cached solve. Is there an easier way to do this?
#     This is very hacky and user-unfriendly (the user should not have to deal
#     with caching), but I don't see a way to get the IPOPT iterations and
#     walltime from a `QuantumControlProblem`.

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
    ("ZeroOrderPulse", N_params, zoh_iters, zoh_secs, F_zoh_piccolo, F_zoh_qutip),
    ("LinearSplinePulse", N_params, lin_iters, lin_secs, F_lin_piccolo, F_lin_qutip),
    ("CubicSplinePulse", N_params, cub_iters, cub_secs, F_cub_piccolo, F_cub_qutip),
]

# We now summarize the results for each pulse type:

println(
    @sprintf(
        "%-18s %5s %6s %9s %11s %11s %10s",
        "Pulse type",
        "N",
        "iters",
        "solve(s)",
        "F Piccolo",
        "F QuTiP",
        "pic−q"
    )
)
for (name, n, it, s, fd, fq) in results
    println(
        @sprintf(
            "%-18s %5d %6d %9.1f %11.7f %11.7f %10.2e",
            name,
            n,
            it,
            s,
            fd,
            fq,
            fd - fq
        )
    )
end

for (name, _, _, _, fd, fq) in results
    @assert fd ≥ 0.999 "$name: fidelity $fd below 0.999 target"
    @assert abs(fd - fq) ≤ 1e-4 "$name: |F_Piccolo - F_QuTiP| = $(abs(fd - fq)) exceeds 1e-4"
end
println("All parameterizations reach ≥ 0.999 and agree with QuantumToolbox to ≤ 1e-4.")

# We see that all pulses synthesize the CNOT gate with ``\geq 99.9`` %
# fidelity, and that the reported fidelity agrees with the rollout fidelity
# computed with QuantumToolbox. Fidelity reported by Piccolo matches the
# rollout fidelity reported by QuantumToolbox matches more closely for
# `ZeroOrderPulse` than for the other pulses because the discretized dynamics
# used in the optimization are exact for `ZeroOrderPulse`.

# ## Step 6: Visualizing Population Dynamics
#
# Finally, we visualize the population dynamics for each pulse type, side by
# side. The CNOT gate flips the target qubit when the control qubit is excited:
# ``|10\rangle \to |11\rangle``. We plot the evolution of the computational-state
# populations starting from the initial condition ``|10\rangle``, using
# QuantumToolbox's rollout of the dynamics.

basis_labels = ["|00⟩", "|01⟩", "|10⟩", "|11⟩"]
sub = U_goal.subspace

let
    rollouts = [
        ("ZeroOrderPulse", get_pulse(qcp_zoh.qtraj)),
        ("LinearSplinePulse", get_pulse(qcp_lin.qtraj)),
        ("CubicSplinePulse", get_pulse(qcp_cub.qtraj)),
    ]
    colors = Makie.wong_colors()[1:length(sub)]

    fig = Figure(size = (1200, 360))
    axes = map(enumerate(rollouts)) do (col, (name, pulse))
        _, sols = unitary_from_qutip(sys, pulse, U_goal)
        sol_10 = sols[3]   # |10⟩ input
        tlist = sol_10.times
        populations = [abs2.(getindex.(getfield.(sol_10.states, :data), j)) for j in sub]
        ax = Axis(
            fig[1, col];
            title = name,
            xlabel = "Time (ns)",
            ylabel = col == 1 ? "Population" : "",
        )
        for (j, population) in enumerate(populations)
            lines!(ax, collect(tlist), population; color = colors[j])
        end
        ax
    end
    linkyaxes!(axes...)

    entries = [LineElement(; color = colors[j], linewidth = 2) for j = 1:length(sub)]
    Legend(fig[1, length(rollouts)+1], entries, basis_labels)
    Label(fig[0, :], "Population dynamics: evolution of |10⟩"; font = :bold, fontsize = 18)
    fig
end
