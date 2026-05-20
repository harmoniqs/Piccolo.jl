module Options

export PiccoloOptions
export display_level
export DISPLAY_SILENT, DISPLAY_COMPACT, DISPLAY_STANDARD, DISPLAY_DETAILED

using ExponentialAction


# TODO: Add duration and symmetry options

# ============================================================================ #
# Display levels
# ============================================================================ #

const DISPLAY_SILENT = 0
const DISPLAY_COMPACT = 1
const DISPLAY_STANDARD = 2
const DISPLAY_DETAILED = 3

const _DISPLAY_LEVELS = (
    silent = DISPLAY_SILENT,
    compact = DISPLAY_COMPACT,
    standard = DISPLAY_STANDARD,
    detailed = DISPLAY_DETAILED,
)

"""
    display_level(s::Symbol) -> Int

Map a display-level symbol to its ordering integer. Comparison-friendly:
`display_level(opts.display) >= DISPLAY_STANDARD`.

| Level         | What you see                                                          |
|---------------|-----------------------------------------------------------------------|
| `:silent`     | Nothing during construction                                           |
| `:compact`    | One-line "constructing X" header per outer template call              |
| `:standard`   | Header + rich problem inspection (tree view, no plot) auto-printed    |
| `:detailed`   | Header + full inspection (tree view + sparsity + terminal pulse plot) |

The default for `PiccoloOptions` is `:standard`.
"""
function display_level(s::Symbol)
    haskey(_DISPLAY_LEVELS, s) || throw(
        ArgumentError(
            "Unknown display level :$s. Expected one of: $(keys(_DISPLAY_LEVELS))",
        ),
    )
    return _DISPLAY_LEVELS[s]
end

# ============================================================================ #
# PiccoloOptions
# ============================================================================ #

"""
    PiccoloOptions

Options for the Piccolo quantum optimal control library.

# Fields
- `display::Symbol = :standard`: Verbosity for construction-time output. One of
  `:silent`, `:compact`, `:standard`, `:detailed`. See [`display_level`](@ref).
- `timesteps_all_equal::Bool = false`: Constrain all timesteps to be equal. When
  `false` (default), each knot-point Δt_k is an independent variable bounded by
  `Δt_bounds`, letting the optimizer choose a non-uniform grid (concentrating
  samples where the pulse is changing fast). Set to `true` to force uniform
  spacing — useful when the target hardware has a fixed sample-rate.
- `rollout_integrator::Function = expv`: Integrator to use for rollout
- `geodesic = true`: Use the geodesic to initialize the optimization.
- `zero_initial_and_final_derivative::Bool=false`: Zero the initial and final control pulse derivatives.
- `complex_control_norm_constraint_name::Union{Nothing, Symbol} = nothing`: Name of the complex control norm constraint.
- `complex_control_norm_constraint_radius::Float64 = 1.0`: Radius of the complex control norm constraint.
- `bound_state::Bool = true`: Keep the default [-1, 1] box bounds on each component
  of the isomorphism state vector at every knot point. These bounds are set
  automatically during trajectory construction. Set to `false` to remove them,
  giving the optimizer full freedom on state variables.
- `bound_state_l2::Bool = false`: Add nonlinear constraints bounding each complex
  component's magnitude (Re² + Im²) ≤ 1 at every knot point. Tighter than `bound_state`
  but adds NLP complexity (Jacobian/Hessian entries). Not yet supported for
  DensityTrajectory.
- `leakage_constraint::Bool = false`: Suppress leakage with constraint and cost.
- `leakage_constraint_value::Float64 = 1e-2`: Value for the leakage constraint.
- `leakage_cost::Float64 = 1e-2`: Leakage suppression parameter.
"""
mutable struct PiccoloOptions
    display::Symbol
    timesteps_all_equal::Bool
    rollout_integrator::Function
    geodesic::Bool
    zero_initial_and_final_derivative::Bool
    complex_control_norm_constraint_name::Union{Nothing,Symbol}
    complex_control_norm_constraint_radius::Float64
    bound_state::Bool
    bound_state_l2::Bool
    leakage_constraint::Bool
    leakage_constraint_value::Float64
    leakage_cost::Float64
end

function PiccoloOptions(;
    display::Union{Symbol,Nothing} = nothing,
    verbose = nothing,
    timesteps_all_equal::Bool = false,
    rollout_integrator::Function = expv,
    geodesic::Bool = true,
    zero_initial_and_final_derivative::Bool = false,
    complex_control_norm_constraint_name::Union{Nothing,Symbol} = nothing,
    complex_control_norm_constraint_radius::Float64 = 1.0,
    bound_state::Bool = true,
    bound_state_l2::Bool = false,
    leakage_constraint::Bool = false,
    leakage_constraint_value::Float64 = 1e-2,
    leakage_cost::Float64 = 1e-2,
)
    if verbose !== nothing
        display === nothing || throw(
            ArgumentError(
                "PiccoloOptions: pass either `verbose` or `display`, not both " *
                "(got verbose=$verbose, display=:$display).",
            ),
        )
        Base.depwarn(
            "`PiccoloOptions(verbose=$verbose)` is deprecated; use " *
            "`display=:$(verbose ? :standard : :silent)` instead. " *
            "Equivalents: verbose=false → display=:silent, verbose=true → display=:standard.",
            :PiccoloOptions,
        )
        display = verbose ? :standard : :silent
    end
    display === nothing && (display = :standard)
    display_level(display)  # validate
    return PiccoloOptions(
        display,
        timesteps_all_equal,
        rollout_integrator,
        geodesic,
        zero_initial_and_final_derivative,
        complex_control_norm_constraint_name,
        complex_control_norm_constraint_radius,
        bound_state,
        bound_state_l2,
        leakage_constraint,
        leakage_constraint_value,
        leakage_cost,
    )
end

end
