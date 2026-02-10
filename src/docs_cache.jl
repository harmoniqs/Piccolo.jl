module DocsCache

using JLD2
using NamedTrajectories: NamedTrajectory, load_traj
using ..Control: QuantumControlProblem, solve!, sync_trajectory!, get_trajectory

export cached_solve!

"""
    get_cache_hash(n=7)

Get the short git hash of the current HEAD commit.
Falls back to "unknown" if git is unavailable.
"""
function get_cache_hash(n::Int=7)
    try
        return String(strip(read(`git rev-parse --short=$n HEAD`, String)))
    catch
        return "unknown"
    end
end

"""
    find_cache_file(name, data_dir)

Find a cached `.jld2` file matching `{name}_*.jld2` in `data_dir`.
Returns the full path, or `nothing` if no cache exists.
"""
function find_cache_file(name::String, data_dir::String)
    isdir(data_dir) || return nothing
    for f in readdir(data_dir)
        if startswith(f, "$(name)_") && endswith(f, ".jld2")
            return joinpath(data_dir, f)
        end
    end
    return nothing
end

"""
    save_cache(traj, output, name, hash, data_dir)

Save a `NamedTrajectory` and captured solver output to a `.jld2` file.
"""
function save_cache(
    traj::NamedTrajectory,
    output::String,
    name::String,
    hash::String,
    data_dir::String,
)
    mkpath(data_dir)
    path = joinpath(data_dir, "$(name)_$(hash).jld2")
    jldsave(path; traj=traj, output=output, hash=hash)
    return path
end

"""
    load_cache(path)

Load a cached trajectory and solver output from a `.jld2` file.
Returns a `NamedTuple` with fields `traj`, `output`, and `hash`.
"""
function load_cache(path::String)
    data = load(path)
    return (
        traj=data["traj"]::NamedTrajectory,
        output=get(data, "output", "")::String,
        hash=get(data, "hash", "unknown")::String,
    )
end

"""
    clean_old_caches!(name, data_dir, current_hash)

Remove old cache files for `name` that don't match `current_hash`.
"""
function clean_old_caches!(name::String, data_dir::String, current_hash::String)
    isdir(data_dir) || return
    current_file = "$(name)_$(current_hash).jld2"
    for f in readdir(data_dir)
        if startswith(f, "$(name)_") && endswith(f, ".jld2") && f != current_file
            rm(joinpath(data_dir, f); force=true)
        end
    end
end

"""
    truncate_solver_output(output; head=5, tail=5)

Truncate long solver iteration tables, keeping the first `head` and last `tail`
iteration lines with a `⋮` marker in between. Non-iteration lines (preamble,
summary, etc.) are preserved. Returns the original string if there are fewer
than `head + tail` iterations.
"""
function truncate_solver_output(output::String; head::Int=5, tail::Int=5)
    lines = split(output, '\n')

    # Iteration data lines match: leading whitespace, number, then scientific notation
    iter_indices = findall(l -> occursin(r"^\s+\d+\s+[\d.e+-]", l), lines)

    n_iters = length(iter_indices)
    if n_iters <= head + tail
        return output
    end

    # Keep lines up to and including the last head iteration
    last_head = iter_indices[head]

    # First tail iteration line
    first_tail = iter_indices[n_iters - tail + 1]

    # Check for an iteration header line just before the tail section
    has_header = first_tail > 1 && occursin(r"^iter\s+objective", lines[first_tail - 1])

    result = lines[1:last_head]
    push!(result, "⋮")
    if has_header
        push!(result, lines[first_tail - 1])
    end
    append!(result, lines[first_tail:end])

    return join(result, '\n')
end

"""
    capture_output(f)

Execute `f()` while capturing all stdout output (including C-level I/O from
solvers like Ipopt). Returns the captured output as a `String`.
"""
function capture_output(f::Function)
    output_file = tempname()
    local result
    open(output_file, "w") do io
        redirect_stdio(; stdout=io) do
            result = f()
        end
    end
    output = read(output_file, String)
    rm(output_file; force=true)
    return output
end

_maybe_snip(output::String, snip::Bool) =
    snip ? truncate_solver_output(output) : output
_maybe_snip(output::String, snip::Tuple{Int,Int}) =
    truncate_solver_output(output; head=snip[1], tail=snip[2])

"""
    cached_solve!(qcp, name; data_dir, force, verbose, kwargs...)

Solve a `QuantumControlProblem` with transparent caching.

If a cached solution exists in `data_dir` for `name`, loads the trajectory and
replays the saved solver output. Otherwise, runs `solve!`, captures the output,
and saves both to a `.jld2` file named `{name}_{git_hash}.jld2`.

# Arguments
- `qcp::QuantumControlProblem`: The problem to solve.
- `name::String`: Unique cache identifier (e.g., `"multilevel_transmon"`).

# Keyword Arguments
- `data_dir::String`: Directory for cache files.
  Defaults to `{project_root}/data`.
- `force::Bool`: Force regeneration even if cache exists.
  Defaults to `false`, or `true` if `ENV["PICCOLO_REGENERATE_CACHE"] == "true"`.
- `verbose::Bool`: Print solver output and info messages. Default `true`.
- `snip::Union{Bool, Tuple{Int,Int}}`: Truncate long Ipopt iteration tables.
  `true` (default) keeps the first 5 and last 5 iterations.
  A tuple `(head, tail)` specifies custom counts. `false` disables truncation.
- `kwargs...`: Forwarded to `solve!` (e.g., `max_iter`, `print_level`).

# Example
```julia
qcp = SmoothPulseProblem(qtraj, N; Q=100.0)
cached_solve!(qcp, "my_gate"; max_iter=50)
fidelity(qcp)  # works as if solve! was called directly
```
"""
function cached_solve!(
    qcp::QuantumControlProblem,
    name::String;
    data_dir::String = joinpath(
        dirname(something(Base.active_project(), @__DIR__)), "data"
    ),
    force::Bool = get(ENV, "PICCOLO_REGENERATE_CACHE", "false") == "true",
    verbose::Bool = true,
    snip::Union{Bool, Tuple{Int,Int}} = true,
    kwargs...,
)
    mkpath(data_dir)
    cache_path = find_cache_file(name, data_dir)

    # Load from cache if available and not forcing regeneration
    if !isnothing(cache_path) && !force
        data = load_cache(cache_path)
        qcp.prob.trajectory = data.traj
        sync_trajectory!(qcp)
        verbose && print(_maybe_snip(data.output, snip))
        verbose && @info "Loaded cached trajectory from $(basename(cache_path))"
        return nothing
    end

    # Run solve! and capture output
    output = capture_output() do
        solve!(qcp; kwargs...)
    end
    verbose && print(_maybe_snip(output, snip))

    # Save cache
    hash = get_cache_hash()
    saved_path = save_cache(get_trajectory(qcp), output, name, hash, data_dir)
    clean_old_caches!(name, data_dir, hash)
    verbose && @info "Saved cache to $(basename(saved_path))"

    return nothing
end

end # module DocsCache
