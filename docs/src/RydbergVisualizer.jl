using GLMakie
using Colors

"""
    visualize_rydberg_chain(N::Int, distance::Float64; C6=1.0)

Generates a 1D visualization of N atoms trapped in optical tweezers.
Interaction strengths are visualized as connective 'bonds' where 
thickness and opacity follow the 1/r^6 Rydberg blockade scaling.
"""
function visualize_rydberg_chain(N=7, distance=4.0; C6=1.0)
    # 1. Setup Data
    xs = collect(1:N) .* distance
    ys = zeros(N)
    
    fig = Figure(resolution = (1200, 600), backgroundcolor = :white)
    ax = Axis(fig[1, 1], title = "Neutral Atom Chain (Rydberg Interaction Visualization)",
         xlabel = "Position (μm)", ylabel = "", yticks = ([-1, 1], ["", ""]))
    
    hidedecorations!(ax, label = false, ticklabels = false, ticks = false)

    # 2. Calculate Interactions (C6 / r^6)
    # We visualize all-to-all interactions, though nearest-neighbor dominates.
    for i in 1:N, j in (i+1):N
        r = abs(xs[i] - xs[j])
        strength = C6 / (r^6)
        
        # Normalize strength for visual clarity (blockade radius representation)
        linewidth = clamp(strength * 5000, 0.5, 12.0)
        alpha = clamp(strength * 2000, 0.1, 0.8)
        
        lines!(ax, [xs[i], xs[j]], [0, 0], 
               color = (:black, alpha), 
               linewidth = linewidth,
               linestyle = :solid)
    end

    # 3. Atoms (Represented as Spheres)
    # State: 0.0 = Ground |g>, 1.0 = Rydberg |r>
    atom_states = Observable(zeros(N))
    
    # Map states to colors: Ground (Cyan) to Rydberg (Magenta) - Harmoniqs Palette
    colors = lift(atom_states) do states
        [interpolate_colors(cgrad([:cyan, :magenta]), s) for s in states]
    end
    
    # Glow effect for Rydberg atoms
    sizes = lift(atom_states) do states
        [20 + s * 15 for s in states]
    end

    scatter!(ax, xs, ys, 
             color = colors, 
             markersize = sizes, 
             strokewidth = 2, 
             strokecolor = :black)

    # 4. Animation logic: Simulate a Global π-pulse
    record(fig, "rydberg_chain_excitation.mp4", 1:120; framerate = 30) do frame
        # Sinusoidal excitation pattern
        phase = frame / 20
        atom_states[] = [0.5 * (1 + sin(phase - i*0.5)) for i in 1:N]
    end

    return fig
end

# Helper to map state value to color gradient
function interpolate_colors(gradient, val)
    return gradient[val]
end

# Generate the visualization
fig = visualize_rydberg_chain(9, 5.0)
display(fig)
