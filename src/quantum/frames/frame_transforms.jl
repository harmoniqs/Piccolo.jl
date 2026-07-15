# frame_transforms.jl — to_lab_frame / to_rotating_frame + drive reconstruction.
# (Populated in Tasks A2–A3.)

export to_lab_frame

# Group physical drive indices by subsystem into quadrature pairs / lone reals.
# Returns Vector of (subsystem, kind, ix, iy, sx, sy) where kind ∈ (:pair, :real).
function _grouped_drives(spec::FrameSpec, n_drives::Int)
    length(spec.drive_map) == n_drives ||
        error("FrameSpec drive_map has $(length(spec.drive_map)) entries; system has $n_drives drives")
    groups = Any[]
    pending = Dict{Int,Tuple{Int,Float64}}()   # subsystem => (x-index, x-sign) awaiting its :y
    for (k, (sub, q, sg)) in enumerate(spec.drive_map)
        if q === :real
            push!(groups, (sub, :real, k, 0, sg, 0.0))
        elseif q === :x
            haskey(pending, sub) && error("two :x drives on subsystem $sub without a :y")
            pending[sub] = (k, sg)
        elseif q === :y
            haskey(pending, sub) || error("drive $k is :y on subsystem $sub but no preceding :x — a single quadrature cannot form (Ω, φ); mark it :real to carrier-modulate directly")
            (ix, sx) = pending[sub]; delete!(pending, sub)
            push!(groups, (sub, :pair, ix, k, sx, sg))
        else
            error("unknown quadrature $q (expected :x, :y, :real)")
        end
    end
    isempty(pending) || error("unpaired :x drive(s) on subsystem(s) $(collect(keys(pending))) — a single quadrature cannot form (Ω, φ); mark :real to carrier-modulate directly")
    return groups
end

"""
    to_lab_frame(sys_rot::QuantumSystem, frame::RotatingFrame, spec::FrameSpec) -> QuantumSystem

Transform a rotating-frame system into the physically correct **lab frame**:
adds the frame generator `Σ_i ω_i n_i` back into the drift (bare oscillator) and
replaces each rotating-frame quadrature pair `(u_x, u_y)` on subsystem `i` with a
carrier-modulated real field `Ω_i cos(ω_{d,i} t + φ_i)(a_i + a_i^†)`,
`Ω = √(u_x²+u_y²)`, `φ = atan2(u_y, u_x)`. Returns a `time_dependent=true`
function-based `QuantumSystem` (rolled via the `Rollouts` ODE path, not
`SplineIntegrator`).
"""
function to_lab_frame(sys_rot::QuantumSystem, frame::RotatingFrame, spec::FrameSpec)
    levels = sys_rot.levels
    n_sub = length(spec.number_ops)
    ωs = frame_frequencies(frame, n_sub)

    Hf = sum(ωs[i] * spec.number_ops[i] for i in 1:n_sub; init = zeros(ComplexF64, levels, levels))
    H_rot_drift = Matrix{ComplexF64}(sys_rot.H(zeros(sys_rot.n_drives), 0.0))
    H_lab_drift = H_rot_drift + Hf

    ops = spec.drive_ops
    groups = _grouped_drives(spec, sys_rot.n_drives)
    ω_of_sub = Dict(i => ωs[i] for i in 1:n_sub)

    function H_lab(u, t)
        H = copy(H_lab_drift)
        for (sub, kind, ix, iy, sx, sy) in groups
            ωd = ω_of_sub[sub]
            if kind === :pair
                ux = sx * u[ix]; uy = sy * u[iy]
                Ω = sqrt(ux^2 + uy^2); φ = atan(uy, ux)
                H += Ω * cos(ωd * t + φ) * (2 * ops[ix])
            else  # :real — the op IS already the physical field; modulate directly (no 2×)
                H += (sx * u[ix]) * cos(ωd * t) * ops[ix]
            end
        end
        return H
    end

    return QuantumSystem(H_lab, sys_rot.drive_bounds; time_dependent = true,
                         global_params = sys_rot.global_params,
                         hermitian = sys_rot.hermitian)
end

export to_rotating_frame

"""
    to_rotating_frame(sys_lab::QuantumSystem, frame::RotatingFrame, spec::FrameSpec;
                      rwa::Bool=true) -> QuantumSystem

Inverse of `to_lab_frame`. Subtracts the frame generator `Σ ω_i n_i` from the
drift and, under `rwa=true`, reduces each carrier field back to its slow
quadrature pair `(u_x, u_y)` (dropping the counter-rotating e^{±2iωt} terms) — a
time-independent rotating-frame `QuantumSystem`. `rwa=false` retains explicit
time-dependence.
"""
function to_rotating_frame(sys_lab::QuantumSystem, frame::RotatingFrame, spec::FrameSpec;
                           rwa::Bool = true)
    levels = sys_lab.levels
    n_sub = length(spec.number_ops)
    ωs = frame_frequencies(frame, n_sub)
    Hf = sum(ωs[i] * spec.number_ops[i] for i in 1:n_sub; init = zeros(ComplexF64, levels, levels))
    H_lab_drift = Matrix{ComplexF64}(sys_lab.H(zeros(sys_lab.n_drives), 0.0))
    H_rot_drift = H_lab_drift - Hf

    ops = spec.drive_ops
    groups = _grouped_drives(spec, sys_lab.n_drives)

    if !rwa
        function H_full(u, t)
            H = copy(H_rot_drift)
            for (sub, kind, ix, iy, sx, sy) in groups
                ωd = ωs[sub]
                if kind === :pair
                    ux = sx * u[ix]; uy = sy * u[iy]
                    Ω = sqrt(ux^2 + uy^2); φ = atan(uy, ux)
                    H += Ω * cos(ωd * t + φ) * (2 * ops[ix])
                else
                    H += (sx * u[ix]) * cos(ωd * t) * ops[ix]
                end
            end
            return H
        end
        return QuantumSystem(H_full, sys_lab.drive_bounds; time_dependent = true,
                             global_params = sys_lab.global_params, hermitian = sys_lab.hermitian)
    end

    function H_slow(u, t)
        H = copy(H_rot_drift)
        for (sub, kind, ix, iy, sx, sy) in groups
            if kind === :pair
                H += (sx * u[ix]) * ops[ix]
                H += (sy * u[iy]) * ops[iy]
            else
                H += (sx * u[ix]) * ops[ix]
            end
        end
        return H
    end
    return QuantumSystem(H_slow, sys_lab.drive_bounds; time_dependent = false,
                         global_params = sys_lab.global_params, hermitian = sys_lab.hermitian)
end

@testitem "to_lab_frame: single transmon builds a carrier-modulated time-dependent system" begin
    using Piccolo
    using LinearAlgebra
    levels = 3
    a = annihilate(levels)
    nop = Matrix(a' * a)
    Hx = 0.5 * Matrix(a + a')            # rotating-frame quadratures
    Hy = 0.5 * Matrix(im * (a' - a))     #   (see substrate_transmon convention)
    α = -0.2; Δ = 0.0
    Hdrift_rot = Δ * nop + 0.5 * α * Matrix(a' * a' * a * a)
    sys_rot = QuantumSystem(Hdrift_rot, [Hx, Hy], [0.05, 0.05])

    ωd = 4.0
    spec = FrameSpec(number_ops = [nop], drive_map = [(1, :x, +1.0), (1, :y, +1.0)],
                     drive_ops = [Hx, Hy])
    sys_lab = to_lab_frame(sys_rot, RotatingFrame(ωd), spec)

    @test sys_lab isa QuantumSystem
    @test sys_lab.time_dependent
    @test sys_lab.n_drives == sys_rot.n_drives
    # lab drift adds back the frame generator ω n
    @test norm(Matrix(sys_lab.H([0.0, 0.0], 0.0)) - (Hdrift_rot + ωd * nop)) < 1e-10
    # the drive oscillates in absolute time (carrier present)
    u = [0.03, 0.01]
    H0 = Matrix(sys_lab.H(u, 0.0))
    Hq = Matrix(sys_lab.H(u, 2π / ωd / 4))     # quarter carrier period later
    @test norm(H0 - Hq) > 1e-6
    @test ishermitian(H0) && ishermitian(Hq)
end

@testitem "to_lab_frame: edge cases (no drives, unpaired quadrature)" begin
    using Piccolo
    using LinearAlgebra
    levels = 3
    a = annihilate(levels)
    nop = Matrix(a' * a)
    # (i) no drives ⇒ pure drift relabel (bare oscillator, no carrier)
    sys0 = QuantumSystem(0.0 * nop + 0.5 * (-0.2) * Matrix(a' * a' * a * a))
    lab0 = to_lab_frame(sys0, RotatingFrame(4.0),
        FrameSpec(number_ops = [nop], drive_map = [], drive_ops = []))
    @test lab0.n_drives == 0
    @test norm(Matrix(lab0.H(Float64[], 0.0)) - Matrix(sys0.H(Float64[], 0.0) + 4.0 * nop)) < 1e-10
    # (ii) unpaired quadrature (single :x with no matching :y) ⇒ clear error
    Hx = 0.5 * Matrix(a + a')
    sys1 = QuantumSystem(0.0 * nop, [Hx], [0.05])
    @test_throws ErrorException to_lab_frame(sys1, RotatingFrame(4.0),
        FrameSpec(number_ops = [nop], drive_map = [(1, :x, +1.0)], drive_ops = [Hx]))
    # (iii) lone :real drive ⇒ allowed, carrier-modulated directly
    lab1 = to_lab_frame(sys1, RotatingFrame(4.0),
        FrameSpec(number_ops = [nop], drive_map = [(1, :real, +1.0)], drive_ops = [Hx]))
    @test lab1.time_dependent && lab1.n_drives == 1
end

@testitem "frames: round-trip to_rotating_frame(to_lab_frame(s)) ≈ s (RWA)" begin
    using Piccolo
    using LinearAlgebra
    levels = 3
    a = annihilate(levels)
    nop = Matrix(a' * a)
    Hx = 0.5 * Matrix(a + a'); Hy = 0.5 * Matrix(im * (a' - a))
    α = -0.2; Δ = 0.01
    Hdrift = Δ * nop + 0.5 * α * Matrix(a' * a' * a * a)
    sys_rot = QuantumSystem(Hdrift, [Hx, Hy], [0.05, 0.05])
    spec = FrameSpec(number_ops = [nop], drive_map = [(1, :x, +1.0), (1, :y, +1.0)],
                     drive_ops = [Hx, Hy])
    frame = RotatingFrame(4.0)

    back = to_rotating_frame(to_lab_frame(sys_rot, frame, spec), frame, spec; rwa = true)

    for u in ([0.0, 0.0], [0.03, 0.0], [0.0, 0.02], [0.02, -0.03]), t in (0.0, 1.3, 7.1)
        @test norm(Matrix(back.H(u, t)) - Matrix(sys_rot.H(u, t))) < 1e-8
    end
end
