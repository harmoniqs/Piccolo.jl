module LiveCallbacks

export LivePulsePlotCallback

"""
    LivePulsePlotCallback(
        qtraj::AbstractQuantumTrajectory,
        traj::NamedTrajectory;
        every::Int = 1,
        save_dir::Union{Nothing,String} = nothing,
        title_prefix::String = "iter",
    )

Solver-agnostic per-iteration callback for trajectory optimization. Each
time it fires it reconstructs the current pulse via `extract_pulse(qtraj, traj)`,
renders it with `plot_pulse(pulse; title = "\$title_prefix \$iter")`, and
(if `save_dir` is set) writes a PNG snapshot to disk.

The returned callback subtypes `DirectTrajOpt.AbstractIntermediateCallback`,
so it can be passed directly to any DTO solver backend:

```julia
cb = LivePulsePlotCallback(qtraj, qcp.prob.trajectory; save_dir = "out/")
solve!(qcp; options = MadNLPOptions(intermediate_callback = cb))
```

# Arguments
- `qtraj`: the original `AbstractQuantumTrajectory` (provides pulse type + drive names).
- `traj`: the `NamedTrajectory` being optimized in place; mutated by the callback each iter.

# Keyword Arguments
- `every`: redraw cadence in IPM iterations. `1` plots every iter.
- `save_dir`: if non-`nothing`, save `iter_NNN.png` per rendered iter.
- `title_prefix`: title-bar prefix, joined to the iter count.

# MadNLP note
The callback needs the full primal vector (including fixed boundary variables) to
reconstruct the trajectory faithfully. DTO sets `fixed_variable_treatment =
MadNLP.RelaxBound` automatically when it sees an `AbstractIntermediateCallback`
installed without an explicit value — you don't need to set it yourself.

Defined as a stub here; the implementation is provided by the `PiccoloMakieExt`
extension and is loaded when a Makie backend (`CairoMakie`, `GLMakie`, `WGLMakie`)
is available.
"""
function LivePulsePlotCallback end

end
