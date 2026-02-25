# Wigner Marginals — Live Coding Cheat Sheet

## The story arc

```
Step 1  bare axes                        → show the layout system
Step 2  raw lines on subplots            → show something is happening
Step 3  correct marginals + animation    → fix the physics, fix the gotchas
Step 4  fill shading                     → band! + poly! polish
```

---

## What already exists (don't touch)

```julia
# QuantumToolbox creates this figure internally:
#   lyt = GridLayout(fig[1,1])   ← nested layout
#   lyt[1,1] = ax                ← heatmap axis
#   lyt[1,2] = Colorbar          ← only if colorbar=true
#
# Piccolo wraps it, adds:
#   fig[0,1] = Label("Timestep N")
#   fig.attributes[:hm] = hm     ← the heatmap plot object
#   fig.attributes[:ax] = ax

fig, ax, hm = QuantumToolbox.plot_wigner(state; library=Val(:Makie))
```

---

## Step 1 — Bare axes (layout only)

```julia
function QuantumToolbox.plot_wigner(traj, idx;
    show_marginals::Bool = false, kwargs...
)
    fig, ax, hm = QuantumToolbox.plot_wigner(state; library=Val(:Makie), kwargs...)

    fig.attributes[:show_marginals] = show_marginals

    if show_marginals
        # Recover the nested GridLayout that QuantumToolbox created
        lyt = ax.layoutobservables.gridcontent[].parent
        ncols = size(lyt)[2]   # 1 (no colorbar) or 2 (with colorbar)

        ax_top   = Axis(lyt[0, 1])          # row 0 = above heatmap
        ax_right = Axis(lyt[1, ncols + 1])  # next column = right of colorbar

        linkxaxes!(ax_top, ax)     # zoom stays in sync
        linkyaxes!(ax_right, ax)

        hidedecorations!(ax_top,   grid=false)
        hidedecorations!(ax_right, grid=false)

        rowsize!(lyt, 0,          Relative(0.2))
        colsize!(lyt, ncols + 1,  Relative(0.2))
    end
    return fig
end
```

---

## Step 2 — Raw lines (something visible, not yet correct)

```julia
# Naive approach: truncate the edge coords to match W, sum without integrating
W = hm[3][]
n, m = size(W)
xvec = hm[1][][1:n]   # ← edge coords, not centers (slightly wrong)
yvec = hm[2][][1:m]

lines!(ax_top,   xvec, vec(sum(W, dims=2)))   # raw sum, no Δp factor
lines!(ax_right, vec(sum(W, dims=1)), yvec)   # raw sum, no Δx factor
```

Lines show up. Two problems to fix next:
- **Coords are edges not centers** — Makie stores `n+1` edge points for an `n`-point grid
- **Values aren't normalized** — missing the `Δx`/`Δp` integration step

---

## Step 3 — Correct marginals + animation

**The edge→center fix:**

```julia
# hm[1][] has 201 pts, W = hm[3][] has shape 200×200  →  mismatch!
function _hm_centers(coords, n)
    length(coords) == n && return coords
    return (coords[1:end-1] .+ coords[2:end]) ./ 2   # midpoints of edges
end
```

**The marginal math:**

```
W  has shape  (nx, np)   where  W[ix, ip] = Wigner at (xvec[ix], yvec[ip])

x_marg = ∫W dp  ≈  sum(W, dims=2) * Δp   → length nx  → plotted vs xvec  (top)
p_marg = ∫W dx  ≈  sum(W, dims=1) * Δx   → length np  → plotted vs yvec  (right, horizontal)
```

**Replace the raw lines block:**

```julia
W    = hm[3][]
xvec = _hm_centers(hm[1][], size(W, 1))
yvec = _hm_centers(hm[2][], size(W, 2))
Δx   = (xvec[end] - xvec[1]) / (length(xvec) - 1)
Δp   = (yvec[end] - yvec[1]) / (length(yvec) - 1)

x_marg = vec(sum(W, dims=2)) .* Δp
p_marg = vec(sum(W, dims=1)) .* Δx

line_top   = lines!(ax_top,   xvec,   x_marg)
line_right = lines!(ax_right, p_marg, yvec)

# Store for animation
fig.attributes[:line_top]   = line_top
fig.attributes[:line_right] = line_right
fig.attributes[:xvec_marg]  = xvec
fig.attributes[:yvec_marg]  = yvec
```

**The animation gotcha — Lines stores `Vector{Point{2,Float64}}`, not separate x/y:**

```julia
# WRONG:  line[2][] = new_y_values      ← BoundsError, index [2] doesn't exist
# RIGHT:  line[1][] = Point2f.(x, y)    ← replaces the whole point vector
```

**Add to `plot_wigner!`:**

```julia
function Piccolo.plot_wigner!(fig, traj, idx)
    hm   = fig.attributes[:hm][]
    # Use _hm_centers here too — otherwise wigner() gets 201-pt edge vectors
    xvec = _hm_centers(hm[1][], size(hm[3][], 1))
    yvec = _hm_centers(hm[2][], size(hm[3][], 2))
    W    = transpose(wigner(state, xvec, yvec))
    hm[3][] = W
    # ... label update ...

    if fig.attributes[:show_marginals][]
        xvec_m = fig.attributes[:xvec_marg][]
        yvec_m = fig.attributes[:yvec_marg][]
        Δx = (xvec_m[end] - xvec_m[1]) / (length(xvec_m) - 1)
        Δp = (yvec_m[end] - yvec_m[1]) / (length(yvec_m) - 1)

        fig.attributes[:line_top][][1][]   = Point2f.(xvec_m, vec(sum(W, dims=2)) .* Δp)
        fig.attributes[:line_right][][1][] = Point2f.(vec(sum(W, dims=1)) .* Δx, yvec_m)
    end
end
```

`animate_wigner` needs **no changes** — `show_marginals` passes through `kwargs...` automatically.

---

## Step 4 — Fill shading

```julia
fill_color = (:slateblue, 0.4)
line_color = :slateblue

# Top: band fills between 0 and x_marg (vertical fill)
band!(ax_top, xvec, zeros(length(xvec)), x_marg; color=fill_color)
line_top = lines!(ax_top, xvec, x_marg; color=line_color)

# Right: closed polygon for horizontal fill
# Trace (p_marg, yvec) then close back along x=0
poly!(ax_right,
    vcat(Point2f.(p_marg, yvec), [Point2f(0, yvec[end]), Point2f(0, yvec[1])]);
    color=fill_color,
)
line_right = lines!(ax_right, p_marg, yvec; color=line_color)
```

Note: only the **lines** update each frame during animation. The fill is static per-frame
(band and poly would require full recreation to animate, which isn't worth it here).

---

## Quick reference card

| What | How |
|---|---|
| Get nested layout from axis | `ax.layoutobservables.gridcontent[].parent` |
| Add axis above | `Axis(lyt[0, 1])` |
| Add axis right of col N | `Axis(lyt[1, N+1])` |
| Link x-axes | `linkxaxes!(ax_top, ax)` |
| Link y-axes | `linkyaxes!(ax_right, ax)` |
| Size marginal rows/cols | `rowsize!(lyt, 0, Relative(0.2))` |
| Fix edge→center mismatch | `(edges[1:end-1] .+ edges[2:end]) ./ 2` |
| Update a Lines plot | `line[1][] = Point2f.(x, y)` |
| Vertical fill | `band!(ax, x, y_low, y_high)` |
| Horizontal fill | `poly!(ax, vcat(Point2f.(vals, yvec), [Point2f(0,y_end), Point2f(0,y_start)]))` |
