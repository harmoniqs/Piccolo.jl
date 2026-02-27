# # [Pulses](@id pulses-concept)
#
# Pulses parameterize how control amplitudes ``\boldsymbol{u}(t)`` vary over
# time.  The choice of pulse type determines both the NLP structure and the
# smoothness class of the resulting controls.
#
# ## Overview
#
# | Pulse Type | Continuity | Decision Variables | Use With |
# |------------|------------|-------------------|----------|
# | `ZeroOrderPulse` | ``C^{-1}`` (piecewise constant) | ``\boldsymbol{u}_k`` | `SmoothPulseProblem` |
# | `LinearSplinePulse` | ``C^0`` | ``\boldsymbol{u}_k`` (knot values) | `SplinePulseProblem` |
# | `CubicSplinePulse` | ``C^1`` | ``\boldsymbol{u}_k,\, \dot{\boldsymbol{u}}_k`` (values + tangents) | `SplinePulseProblem` |
# | `GaussianPulse` | ``C^\infty`` | ``A_i, \sigma_i, \mu_i`` (parametric) | Analytical |
# | `CompositePulse` | varies | union of sub-pulse variables | Various |
#
# ## ZeroOrderPulse
#
# Piecewise constant (zero-order hold) controls — the most common choice.
# On the interval ``[t_k, t_{k+1})``, the control is constant:
#
# ```math
# \boldsymbol{u}(t) = \boldsymbol{u}_k, \qquad t \in [t_k,\, t_{k+1})
# ```
#
# Smoothness is enforced indirectly via regularization of the discrete
# differences ``\Delta\boldsymbol{u}_k`` and ``\Delta^2\boldsymbol{u}_k``
# (see [Objectives](@ref objectives-concept)).
#
# ### Construction

using Piccolo

n_dr = 2
N = 100
T = 10.0

## Time points
times = collect(range(0, T, length = N))

## Control values: n_drives × N matrix
controls = 0.1 * randn(n_dr, N)

pulse_zop = ZeroOrderPulse(controls, times)

# ### Properties

## Evaluate at time t
u = pulse_zop(T / 2)
u

## Duration
duration(pulse_zop)

## Number of drives
n_drives(pulse_zop)

# ### Visualization
#
# ```
# Control Value
#     │    ┌──┐
#     │    │  │  ┌──┐
#     │────┘  │  │  └───
#     │       └──┘
#     └─────────────────── Time
# ```
#
# ### Use Case
#
# - **Primary use**: `SmoothPulseProblem`
# - **Characteristics**: Simple structure, smoothness via derivative regularization
# - **Best for**: Initial optimization, most quantum control problems
#
# 
# ## LinearSplinePulse
#
# Linear interpolation between knot values.  On ``[t_k, t_{k+1}]``:
#
# ```math
# \boldsymbol{u}(t) = \boldsymbol{u}_k + \frac{t - t_k}{t_{k+1} - t_k}\,(\boldsymbol{u}_{k+1} - \boldsymbol{u}_k)
# ```
#
# This gives ``C^0`` continuity (continuous values, discontinuous first
# derivative at knots).
#
# ### Construction

pulse_linear = LinearSplinePulse(controls, times)

# ### Properties
#
# - Continuous control values
# - Discontinuous first derivative (at knots)
# - Derivative = slope between knots

u_linear = pulse_linear(T / 2)
u_linear

# ### Visualization
#
# ```
# Control Value
#     │      /\
#     │     /  \    /
#     │    /    \  /
#     │───/      \/
#     └─────────────────── Time
# ```
#
# ## CubicSplinePulse
#
# Cubic Hermite spline interpolation with independent tangents
# ``\dot{\boldsymbol{u}}_k`` at each knot.  The Hermite basis gives ``C^1``
# continuity (continuous values and first derivatives).
#
# On ``[t_k, t_{k+1}]`` with ``s = (t - t_k) / (t_{k+1} - t_k) \in [0, 1]``
# and ``h = t_{k+1} - t_k``:
#
# ```math
# \boldsymbol{u}(t) = (2s^3 - 3s^2 + 1)\,\boldsymbol{u}_k
# + (s^3 - 2s^2 + s)\,h\,\dot{\boldsymbol{u}}_k
# + (-2s^3 + 3s^2)\,\boldsymbol{u}_{k+1}
# + (s^3 - s^2)\,h\,\dot{\boldsymbol{u}}_{k+1}
# ```
#
# Both ``\boldsymbol{u}_k`` and ``\dot{\boldsymbol{u}}_k`` are decision
# variables in `SplinePulseProblem`.
#
# ### Construction

tangents = zeros(n_dr, N)  # Initial tangents (slopes)
pulse_cubic = CubicSplinePulse(controls, tangents, times)

# ### Properties
#
# - Continuous control values AND first derivatives
# - Tangents are independent optimization variables
# - Smooth C¹ continuous curves

u_cubic = pulse_cubic(T / 2)
u_cubic

# ## GaussianPulse
#
# Parametric Gaussian envelope:
#
# ```math
# u_i(t) = A_i \exp\!\left(-\frac{(t - \mu_i)^2}{2\sigma_i^2}\right)
# ```
#
# Useful for analytical pulse design and warm-starting.
#
# ### Construction

amplitudes = [0.5, 0.3]
sigmas = [1.0, 1.5]
centers = [5.0, 5.0]

pulse_gauss = GaussianPulse(amplitudes, sigmas, centers, T)

u_gauss = pulse_gauss(5.0)
u_gauss

# ## CompositePulse
#
# Combine multiple pulses (e.g., different drives from different sources):

pulse1 = GaussianPulse([0.5], [0.5], [2.0], T)
pulse2 = GaussianPulse([0.3], [0.5], [8.0], T)

composite = CompositePulse([pulse1, pulse2])
n_drives(composite)

# ## Choosing a Pulse Type
#
# | Scenario | Recommended Pulse | Smoothness |
# |----------|-------------------|------------|
# | Starting fresh | `ZeroOrderPulse` + `SmoothPulseProblem` | ``C^{-1}`` (regularized) |
# | Refining a solution | `CubicSplinePulse` + `SplinePulseProblem` | ``C^1`` |
# | Hardware requires smooth pulses | `CubicSplinePulse` | ``C^1`` |
# | Simple continuous pulses | `LinearSplinePulse` | ``C^0`` |
# | Analytical pulse design | `GaussianPulse` | ``C^\infty`` |
#
# ## Converting Between Pulse Types
#
# ### ZeroOrderPulse → CubicSplinePulse

## Sample control values from the zero-order pulse
ctrl = hcat([pulse_zop(t) for t in times]...)

## Estimate tangents (finite differences)
tgts = similar(ctrl)
for k = 1:(N-1)
    tgts[:, k] = (ctrl[:, k+1] - ctrl[:, k]) / (times[k+1] - times[k])
end
tgts[:, N] = tgts[:, N-1]

## Create cubic spline
cubic_from_zop = CubicSplinePulse(ctrl, tgts, times)
duration(cubic_from_zop)

# ### Arbitrary → ZeroOrderPulse

## Sample any pulse type to create a ZeroOrderPulse
new_times = collect(range(0, T, length = 200))
new_ctrl = hcat([pulse_cubic(t) for t in new_times]...)
resampled = ZeroOrderPulse(new_ctrl, new_times)
length(new_times)

# ## Best Practices
#
# ### 1. Initialize Appropriately
#
# ```julia
# max_amp = 0.1 * maximum(drive_bounds)
# controls = max_amp * randn(n_drives, N)
# ```
#
# ### 2. Use Enough Time Points
#
# ```julia
# # Rule of thumb: ~10 points per characteristic time scale
# T = 10.0  # Total time
# τ = 1.0   # Shortest feature you want to capture
# N = ceil(Int, 10 * T / τ)
# ```
#
# ### 3. Save Optimized Pulses
#
# Pulses are the primary output of Piccolo.jl. Save them after every
# successful solve so you can warm-start future optimizations, share
# results, or deploy to hardware:
#
# ```julia
# using JLD2
# optimized_pulse = get_pulse(qcp.qtraj)
# jldsave("my_gate.jld2"; pulse=optimized_pulse, fidelity=fidelity(qcp))
# ```
#
# See [Saving and Loading Pulses](@ref saving-loading) for the full guide.
#
# ### 4. Start with ZeroOrderPulse
#
# Even if you need smooth pulses, optimize with `ZeroOrderPulse` first, then
# convert to `CubicSplinePulse` for refinement.
#
# ## See Also
#
# - [SmoothPulseProblem](@ref smooth-pulse) - Using `ZeroOrderPulse`
# - [SplinePulseProblem](@ref spline-pulse) - Using spline pulses
# - [Trajectories](@ref trajectories-concept) - Combining pulses with systems
