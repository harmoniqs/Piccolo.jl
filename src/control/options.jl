module Options

export PiccoloOptions

using ExponentialAction


# TODO: Add duration and symmetry options

"""
    PiccoloOptions

Options for the Piccolo quantum optimal control library.

# Fields
- `verbose::Bool = true`: Print verbose output
- `timesteps_all_equal::Bool = true`: Use equal timesteps
- `rollout_integrator::Function = expv`: Integrator to use for rollout
- `geodesic = true`: Use the geodesic to initialize the optimization.
- `zero_initial_and_final_derivative::Bool=false`: Zero the initial and final control pulse derivatives.
- `complex_control_norm_constraint_name::Union{Nothing, Symbol} = nothing`: Name of the complex control norm constraint.
- `complex_control_norm_constraint_radius::Float64 = 1.0`: Radius of the complex control norm constraint.
- `bound_state::Bool = false`: Bound the state variables <= 1.0.
- `leakage_constraint::Bool = false`: Suppress leakage with constraint and cost.
- `leakage_constraint_value::Float64 = 1e-2`: Value for the leakage constraint.
- `leakage_cost::Float64 = 1e-2`: Leakage suppression parameter.
"""
@kwdef mutable struct PiccoloOptions
    verbose::Bool = true
    timesteps_all_equal::Bool = true
    rollout_integrator::Function = expv
    geodesic::Bool = true
    zero_initial_and_final_derivative::Bool = false
    complex_control_norm_constraint_name::Union{Nothing,Symbol} = nothing
    complex_control_norm_constraint_radius::Float64 = 1.0
    bound_state::Bool = false
    leakage_constraint::Bool = false
    leakage_constraint_value::Float64 = 1e-2
    leakage_cost::Float64 = 1e-2
end

end
