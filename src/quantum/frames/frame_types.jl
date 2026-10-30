# frame_types.jl — frame descriptors and drive/subsystem metadata.
# (Populated in Task A1.)

export AbstractFrame, RotatingFrame, LabFrame, FrameSpec, frame_frequencies

"""
    AbstractFrame

A reference frame for `to_lab_frame`/`to_rotating_frame`. A frame is defined by
its generator `H_f = Σ_i ω_i n_i` (per-subsystem frame frequencies).
"""
abstract type AbstractFrame end

"""
    RotatingFrame(ωs::Vector{Float64})
    RotatingFrame(ω::Real)

Rotating frame carrying one frame/carrier frequency per subsystem (single-carrier
= length 1). Angular units (rad/time) to match the systems it transforms.
"""
struct RotatingFrame <: AbstractFrame
    ωs::Vector{Float64}
end
RotatingFrame(ω::Real) = RotatingFrame(Float64[ω])

"""
    LabFrame()

The bare/lab frame — the ω = 0 frame (all frame frequencies zero).
"""
struct LabFrame <: AbstractFrame end

"""
    frame_frequencies(frame, n_subsystems) -> Vector{Float64}

Per-subsystem frame frequencies. `LabFrame` ⇒ zeros; `RotatingFrame` ⇒ its `ωs`
(broadcast to `n_subsystems` if length 1).
"""
frame_frequencies(::LabFrame, n::Int) = zeros(Float64, n)
function frame_frequencies(f::RotatingFrame, n::Int)
    length(f.ωs) == n && return f.ωs
    length(f.ωs) == 1 && return fill(f.ωs[1], n)
    error("RotatingFrame has $(length(f.ωs)) frequencies but system has $n subsystems")
end

"""
    FrameSpec(; number_ops, drive_map, drive_ops)

Subsystem + drive metadata the frame layer needs but a bare `QuantumSystem` does
not carry. **All three fields are the "input data" the layer reads and never
guesses** — both `to_lab_frame` and `to_rotating_frame` read operators from a
`FrameSpec`, never from a system's `H_drives` (a function-based lab system stores
none). Fields:
- `number_ops[i]` — number operator n_i of subsystem i, **lifted to the full
  Hilbert space** (the frame generator pieces `H_f = Σ ω_i n_i`).
- `drive_map[k] = (subsystem::Int, quadrature::Symbol, sign::Float64)` — physical
  drive index `k`, `quadrature ∈ (:x, :y, :real)`. `:x`/`:y` are the two members
  of an (u_x, u_y) pair combined into (Ω, φ); `:real` is a lone real drive to
  carrier-modulate directly.
- `drive_ops[k]` — the operator matrix physical control `k` multiplies in the
  rotating frame (e.g. `½(a+a†)` for `:x`, `½·i(a†−a)` for `:y`). The lab **field**
  operator for that drive is `2·drive_ops[k]` (undoes the ½), and the RWA inverse
  reuses `drive_ops[k]` directly — so the ½/sign convention lives entirely in this
  operator, pinned by the round-trip test.
"""
struct FrameSpec
    number_ops::Vector{Matrix{ComplexF64}}
    drive_map::Vector{Tuple{Int,Symbol,Float64}}
    drive_ops::Vector{Matrix{ComplexF64}}
end
function FrameSpec(; number_ops, drive_map, drive_ops)
    dmap = [(Int(s), Symbol(q), Float64(sg)) for (s, q, sg) in drive_map]
    length(drive_ops) == length(dmap) || error(
        "FrameSpec: drive_ops ($(length(drive_ops))) must match drive_map ($(length(dmap)))",
    )
    return FrameSpec(Matrix{ComplexF64}.(number_ops), dmap, Matrix{ComplexF64}.(drive_ops))
end

"""
    FrameSpec(sys::QuantumSystem, drive_map; number_ops) -> FrameSpec

Build a spec for a matrix-drive-built `QuantumSystem`, reading `drive_ops` from its
`H_drives` (`drive_matrix.(sys.H_drives)`, which exist on matrix/typed-drive
systems). `number_ops` must be supplied (a bare system has no subsystem structure).
"""
function FrameSpec(sys::QuantumSystem, drive_map; number_ops)
    ops = [Matrix{ComplexF64}(drive_matrix(d)) for d in sys.H_drives]
    return FrameSpec(number_ops = number_ops, drive_map = drive_map, drive_ops = ops)
end

"""
    FrameSpec(sys::CompositeQuantumSystem) -> FrameSpec

Fully auto-derive from a Piccolo-built multi-transmon composite: `n_i = a_i^† a_i`
lifted per subsystem; drive map = each subsystem's (x, y) quadrature pair in
drive-index order (coupling drives, if any, are skipped — no single-subsystem
carrier); `drive_ops` read from `sys.H_drives` at the matching indices. Uses
`subsystem_levels`/`subsystems`.
"""
function FrameSpec(sys::CompositeQuantumSystem)
    levels = sys.subsystem_levels
    nops = Matrix{ComplexF64}[]
    for i in eachindex(sys.subsystems)
        a = annihilate(levels[i])
        push!(nops, Matrix{ComplexF64}(lift_operator(a' * a, i, levels)))
    end
    all_ops = [Matrix{ComplexF64}(drive_matrix(d)) for d in sys.H_drives]
    dmap = Tuple{Int,Symbol,Float64}[]
    ops = Matrix{ComplexF64}[]
    n_sub_drives = sum(s.n_drives for s in sys.subsystems)
    n_coupling = sys.n_drives - n_sub_drives
    k = n_coupling
    for (i, s) in enumerate(sys.subsystems)
        s.n_drives == 0 && continue
        s.n_drives == 2 || error(
            "FrameSpec auto-derivation expects 2 quadrature drives per driven subsystem (got $(s.n_drives) on subsystem $i)",
        )
        push!(dmap, (i, :x, +1.0))
        push!(ops, all_ops[k+1])
        push!(dmap, (i, :y, +1.0))
        push!(ops, all_ops[k+2])
        k += 2
    end
    return FrameSpec(nops, dmap, ops)
end

@testitem "frame descriptors: RotatingFrame / LabFrame" begin
    using Piccolo
    rf = RotatingFrame([4.0, 4.1])
    @test rf.ωs == [4.0, 4.1]
    @test LabFrame() isa AbstractFrame            # exported by Frames
    @test RotatingFrame(4.0).ωs == [4.0]          # scalar convenience
    # LabFrame is the ω=0 rotating frame
    @test frame_frequencies(LabFrame(), 2) == [0.0, 0.0]
end

@testitem "FrameSpec: explicit, from-system, and from-composite" begin
    using Piccolo
    using LinearAlgebra
    levels = 3
    a = annihilate(levels)
    nop = Matrix(a' * a)
    Hx = 0.5 * Matrix(a + a')
    Hy = 0.5 * Matrix(im * (a' - a))
    spec = FrameSpec(
        number_ops = [nop],
        drive_map = [(1, :x, +1.0), (1, :y, +1.0)],
        drive_ops = [Hx, Hy],
    )
    @test length(spec.number_ops) == 1
    @test spec.drive_map[2] == (1, :y, +1.0)
    @test length(spec.drive_ops) == 2
    @test spec.drive_ops[1] ≈ Hx

    sys = QuantumSystem(nop, [Hx, Hy], [0.05, 0.05])
    sspec = FrameSpec(sys, [(1, :x, +1.0), (1, :y, +1.0)]; number_ops = [nop])
    @test length(sspec.drive_ops) == sys.n_drives
    @test sspec.drive_ops[2] ≈ Hy

    comp = MultiTransmonSystem([4.0, 4.1], [0.2, 0.2], [0.0 0.01; 0.01 0.0])
    dspec = FrameSpec(comp)
    @test length(dspec.number_ops) == 2
    @test all(op -> size(op, 1) == comp.levels, dspec.number_ops)
    @test length(dspec.drive_map) == comp.n_drives
    @test length(dspec.drive_ops) == comp.n_drives
end
