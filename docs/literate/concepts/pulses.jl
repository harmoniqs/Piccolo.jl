# # [Pulses](@id pulses-concept)
#
# Pulses in Piccolo.jl parameterize how control signals vary over time. The
# choice of pulse type affects both the optimization problem structure and the
# resulting control smoothness.
#
# ## Overview
#
# | Pulse Type | Description | Use With |
# |------------|-------------|----------|
# | `ZeroOrderPulse` | Piecewise constant | `SmoothPulseProblem` |
# | `LinearSplinePulse` | Linear interpolation | `SplinePulseProblem` |
# | `CubicSplinePulse` | Cubic Hermite splines | `SplinePulseProblem` |
# | `GaussianPulse` | Parametric Gaussian | Analytical evaluation |
# | `CompositePulse` | Combination of pulses | Various |
#
# ## ZeroOrderPulse
#
# Piecewise constant (zero-order hold) controls. The most common choice for
# initial optimization.
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
println("Control values at t=$(T/2): ", u)

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
# ## LinearSplinePulse
#
# Linear interpolation between control knots.
#
# ### Construction

pulse_linear = LinearSplinePulse(controls, times)

# ### Properties
#
# - Continuous control values
# - Discontinuous first derivative (at knots)
# - Derivative = slope between knots

u_linear = pulse_linear(T / 2)
println("Linear spline at t=$(T/2): ", u_linear)

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
# Cubic Hermite spline interpolation with independent tangents at each knot.
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
println("Cubic spline at t=$(T/2): ", u_cubic)

# ## GaussianPulse
#
# Parametric Gaussian envelope with analytical form.
#
# ### Construction

amplitudes = [0.5, 0.3]
sigmas = [1.0, 1.5]
centers = [5.0, 5.0]

pulse_gauss = GaussianPulse(amplitudes, sigmas, centers, T)

# ### Mathematical Form
#
# ```math
# u_i(t) = A_i \exp\left(-\frac{(t - \mu_i)^2}{2\sigma_i^2}\right)
# ```

u_gauss = pulse_gauss(5.0)
println("Gaussian at center: ", u_gauss)

# ## CompositePulse
#
# Combine multiple pulses.
#
# ### Construction

pulse1 = GaussianPulse([0.5], [0.5], [2.0], T)
pulse2 = GaussianPulse([0.3], [0.5], [8.0], T)

## Interleave the pulses (each pulse contributes different drives)
composite = CompositePulse([pulse1, pulse2])
println("Composite pulse drives: ", n_drives(composite))

# ## Choosing a Pulse Type
#
# ### Decision Guide
#
# ```
# Start
#   │
#   ▼
# Is this your first optimization attempt?
#   │
#   ├── Yes → ZeroOrderPulse + SmoothPulseProblem
#   │
#   └── No → Do you have a previous solution?
#               │
#               ├── Yes → CubicSplinePulse + SplinePulseProblem (warm-start)
#               │
#               └── No → Do you need smooth pulses?
#                           │
#                           ├── Yes → CubicSplinePulse
#                           │
#                           └── No → ZeroOrderPulse
# ```
#
# ### Practical Recommendations
#
# | Scenario | Recommended Pulse |
# |----------|-------------------|
# | Starting fresh | `ZeroOrderPulse` |
# | Refining a solution | `CubicSplinePulse` |
# | Hardware requires smooth pulses | `CubicSplinePulse` |
# | Simple continuous pulses | `LinearSplinePulse` |
# | Analytical pulse design | `GaussianPulse` |
#
# ## Converting Between Pulse Types
#
# ### ZeroOrderPulse → CubicSplinePulse

## Sample control values from the zero-order pulse
ctrl = hcat([pulse_zop(t) for t in times]...)

## Estimate tangents (finite differences)
tgts = similar(ctrl)
for k in 1:N-1
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
# # Scale by drive bounds
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
# ### 3. Start with ZeroOrderPulse
#
# Even if you need smooth pulses, optimize with `ZeroOrderPulse` first, then
# convert to `CubicSplinePulse` for refinement.
#
# ## See Also
#
# - [SmoothPulseProblem](@ref smooth-pulse) - Using `ZeroOrderPulse`
# - [SplinePulseProblem](@ref spline-pulse) - Using spline pulses
# - [Trajectories](@ref trajectories-concept) - Combining pulses with systems
