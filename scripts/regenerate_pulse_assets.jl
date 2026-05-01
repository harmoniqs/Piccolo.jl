# Regenerate PNG demos used in the pulse-plotting PR description.
#
# Usage (from repo root):
#   julia --project=docs scripts/regenerate_pulse_assets.jl
#
# (`--project=docs` picks up CairoMakie, which is in the docs Project.toml. The
# package's own Project.toml only lists CairoMakie under `[targets].test`.)
#
# Writes PNGs to assets/ for every pulse type. Pass `--dark` to render with
# `theme_dark()` and write `assets/<name>_dark.png` variants too.

using Piccolo
using CairoMakie
using Random

const ASSETS = joinpath(@__DIR__, "..", "assets")
const DARK = "--dark" in ARGS

isdir(ASSETS) || mkpath(ASSETS)

function _save(fig, name)
    path = joinpath(ASSETS, name * ".png")
    save(path, fig)
    println("wrote $path")
end

function regenerate(; suffix = "")
    Random.seed!(0)
    T = 10.0
    N_demo = 8
    times = collect(range(0, T, length = N_demo))
    controls = 0.5 * randn(2, N_demo)

    # ZeroOrderPulse
    zop = ZeroOrderPulse(controls, times)
    _save(
        plot_pulse(zop; title = "ZeroOrderPulse", labels = ["Drive 1", "Drive 2"]),
        "zero_order" * suffix,
    )

    # LinearSplinePulse
    lsp = LinearSplinePulse(controls, times)
    _save(
        plot_pulse(lsp; title = "LinearSplinePulse", labels = ["Drive 1", "Drive 2"]),
        "linear_spline" * suffix,
    )

    # CubicSplinePulse — values only
    csp = CubicSplinePulse(controls, times)
    _save(
        plot_pulse(csp; title = "CubicSplinePulse", labels = ["Drive 1", "Drive 2"]),
        "cubic_spline" * suffix,
    )

    # CubicSplinePulse — with tangents
    Random.seed!(1)
    tangents = 0.3 * randn(2, N_demo)
    csp_t = CubicSplinePulse(controls, tangents, times)
    _save(
        plot_pulse(
            csp_t;
            title = "CubicSplinePulse (with tangents)",
            labels = ["Drive 1", "Drive 2"],
            show_tangents = true,
            tangent_scale = 0.05,
        ),
        "cubic_spline_tangents" * suffix,
    )

    # GaussianPulse
    gp = GaussianPulse([0.5, 0.3], [1.0, 1.5], [5.0, 5.0], T)
    _save(
        plot_pulse(gp; title = "GaussianPulse", labels = ["Drive 1", "Drive 2"]),
        "gaussian" * suffix,
    )

    # ErfPulse
    ep = ErfPulse([0.8], 2.0, T)
    _save(plot_pulse(ep; title = "ErfPulse", labels = ["Phase"]), "erf" * suffix)

    # CompositePulse — overlay layout
    amp = GaussianPulse([0.5], 1.5, T)
    phs = ErfPulse([0.8], 2.0, T)
    correction = CubicSplinePulse(
        [0.0 0.02 -0.03 0.05 -0.04 0.03 -0.01 0.0],
        collect(range(0, T, length = 8)),
    )
    cp = CompositePulse([amp, phs, correction], :concatenate)
    _save(
        plot_pulse(
            cp;
            layout = :overlay,
            title = "CompositePulse (Gaussian + Erf + Cubic)",
            labels = ["Amplitude", "Phase", "Correction"],
        ),
        "composite_overlay" * suffix,
    )

    # FunctionPulse
    Tf = 1.0
    fp = FunctionPulse(t -> [sin(π * t / Tf)^2, cos(π * t / Tf)^2], Tf, 2)
    _save(
        plot_pulse(fp; title = "FunctionPulse", labels = ["sin²", "cos²"]),
        "function_pulse" * suffix,
    )

    # Side-by-side comparison: same data, three interpolations
    cmp_t = collect(range(0, T, length = 7))
    cmp_u = [0.0 0.8 -0.5 1.2 0.3 -0.9 0.0]
    fig = Figure(size = (900, 700))
    ax1 = Axis(fig[1, 1]; title = "ZeroOrderPulse", ylabel = "Amplitude")
    plot_pulse!(ax1, ZeroOrderPulse(cmp_u, cmp_t))
    ax2 = Axis(fig[2, 1]; title = "LinearSplinePulse", ylabel = "Amplitude")
    plot_pulse!(ax2, LinearSplinePulse(cmp_u, cmp_t))
    ax3 = Axis(fig[3, 1]; title = "CubicSplinePulse", xlabel = "Time", ylabel = "Amplitude")
    plot_pulse!(ax3, CubicSplinePulse(cmp_u, cmp_t))
    Label(
        fig[0, 1],
        "Same Data, Three Interpolations";
        fontsize = 18,
        font = :bold,
        tellwidth = false,
    )
    _save(fig, "comparison" * suffix)
end

# Always render the default-theme set
set_theme!()
regenerate()

# Optionally render a dark-theme set with `_dark` suffix
if DARK
    set_theme!(theme_dark())
    regenerate(suffix = "_dark")
    set_theme!()
end
