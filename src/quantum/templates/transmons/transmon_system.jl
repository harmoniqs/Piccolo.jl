export TransmonSystem
export TransmonDipoleCoupling
export MultiTransmonSystem
export QuantumSystemCoupling
export TransmonCavitySystem

@doc raw"""
    TransmonSystem(;
        ω::Float64=4.4153,  # GHz
        δ::Float64=0.17215, # GHz
        levels::Int=3,
        lab_frame::Bool=false,
        frame_ω::Float64=ω,
    ) -> QuantumSystem

Returns a `QuantumSystem` object for a transmon qubit, with the Hamiltonian

```math
H = \omega a^\dagger a - \frac{\delta}{2} a^\dagger a^\dagger a a
```

where `a` is the annihilation operator.

# Keyword Arguments
- `ω`: The frequency of the transmon, in GHz.
- `δ`: The anharmonicity of the transmon, in GHz.
- `levels`: The number of levels in the transmon.
- `lab_frame`: Whether to use the lab frame Hamiltonian, or an ω-rotating frame.
- `frame_ω`: The frequency of the rotating frame, in GHz.
- `mutiply_by_2π`: Whether to multiply the Hamiltonian by 2π, set to true by default because the frequency is in GHz.
- `lab_frame_type`: The type of lab frame Hamiltonian to use, one of (:duffing, :quartic, :cosine).
- `drives`: Whether to include drives in the Hamiltonian.
"""
function TransmonSystem(;
    drive_bounds::Vector{<:Union{Tuple{Float64,Float64},Float64}} = fill(1.0, 2),
    ω::Float64 = 4.0,  # GHz
    δ::Float64 = 0.2, # GHz
    levels::Int = 3,
    lab_frame::Bool = false,
    frame_ω::Float64 = lab_frame ? 0.0 : ω,
    multiply_by_2π::Bool = true,
    lab_frame_type::Symbol = :duffing,
    drives::Bool = true,
    kwargs...,
)
    @assert lab_frame_type ∈ (:duffing, :quartic, :cosine) "lab_frame_type must be one of (:duffing, :quartic, :cosine)"

    if lab_frame
        if frame_ω ≉ 0.0
            frame_ω = 0.0
        end
    end

    if frame_ω ≉ 0.0
        lab_frame = false
    end

    a = annihilate(levels)

    if lab_frame
        if lab_frame_type == :duffing
            H_drift = ω * a' * a - δ / 2 * a' * a' * a * a
        elseif lab_frame_type == :quartic
            ω₀ = ω + δ
            H_drift = ω₀ * a' * a - δ / 12 * (a + a')^4
        elseif lab_frame_type == :cosine
            ω₀ = ω + δ
            E_C = δ
            E_J = ω₀^2 / 8E_C
            n̂ = im / 2 * (E_J / 2E_C)^(1 / 4) * (a - a')
            φ̂ = (2E_C / E_J)^(1 / 4) * (a + a')
            H_drift = 4 * E_C * n̂^2 - E_J * cos(φ̂)
            # H_drift = 4 * E_C * n̂^2 - E_J * (I - φ̂^2 / 2 + φ̂^4 / 24)
        end
    else
        H_drift = (ω - frame_ω) * a' * a - δ / 2 * a' * a' * a * a
    end

    if drives
        H_drives = [a + a', 1.0im * (a - a')]
    else
        H_drives = Matrix{ComplexF64}[]
    end

    if multiply_by_2π
        H_drift *= 2π
        H_drives .*= 2π
    end

    return QuantumSystem(H_drift, H_drives, drive_bounds)
end

struct QuantumSystemCoupling
    op::Any
    g_ij::Any
    pair::Any
    subsystem_levels::Any
    coupling_type::Any
    params::Any
end

@doc raw"""
    TransmonDipoleCoupling(
        g_ij::Float64,
        pair::Tuple{Int, Int},
        subsystem_levels::Vector{Int};
        lab_frame::Bool=false,
    ) -> QuantumSystemCoupling

    TransmonDipoleCoupling(
        g_ij::Float64,
        pair::Tuple{Int, Int},
        sub_systems::Vector{QuantumSystem};
        kwargs...
    ) -> QuantumSystemCoupling

Returns a `QuantumSystemCoupling` object for a transmon qubit. In the lab frame, the Hamiltonian coupling term is

```math
H = g_{ij} (a_i + a_i^\dagger) (a_j + a_j^\dagger)
```

In the rotating frame, the Hamiltonian coupling term is

```math
H = g_{ij} (a_i a_j^\dagger + a_i^\dagger a_j)
```

where `a_i` is the annihilation operator for the `i`th transmon.

"""
function TransmonDipoleCoupling end

function TransmonDipoleCoupling(
    g_ij::Float64,
    pair::Tuple{Int,Int},
    subsystem_levels::Vector{Int};
    lab_frame::Bool = false,
    mulitply_by_2π::Bool = true,
)

    i, j = pair
    a_i = lift_operator(annihilate(subsystem_levels[i]), i, subsystem_levels)
    a_j = lift_operator(annihilate(subsystem_levels[j]), j, subsystem_levels)

    if lab_frame
        op = g_ij * (a_i + a_i') * (a_j + a_j')
    else
        op = g_ij * (a_i * a_j' + a_i' * a_j)
    end

    if mulitply_by_2π
        op *= 2π
    end

    params = Dict{Symbol,Any}(:lab_frame => lab_frame, :mulitply_by_2π => mulitply_by_2π)

    return QuantumSystemCoupling(
        op,
        g_ij,
        pair,
        subsystem_levels,
        TransmonDipoleCoupling,
        params,
    )
end

function TransmonDipoleCoupling(
    g_ij::Float64,
    pair::Tuple{Int,Int},
    sub_systems::AbstractVector{<:AbstractQuantumSystem};
    kwargs...,
)
    subsystem_levels = [sys.levels for sys ∈ sub_systems]
    return TransmonDipoleCoupling(g_ij, pair, subsystem_levels; kwargs...)
end

"""
    MultiTransmonSystem(
        ωs::AbstractVector{Float64},
        δs::AbstractVector{Float64},
        gs::AbstractMatrix{Float64};
        drive_bounds::Union{Float64, Vector{<:Union{Tuple{Float64, Float64}, Float64}}}=1.0,
        levels_per_transmon::Int = 3,
        subsystem_levels::AbstractVector{Int} = fill(levels_per_transmon, length(ωs)),
        lab_frame=false,
        subsystems::AbstractVector{Int} = 1:length(ωs),
        subsystem_drive_indices::AbstractVector{Int} = 1:length(ωs),
        kwargs...
    ) -> CompositeQuantumSystem

Returns a `CompositeQuantumSystem` object for a multi-transmon system.
"""
function MultiTransmonSystem(
    ωs::AbstractVector{Float64},
    δs::AbstractVector{Float64},
    gs::AbstractMatrix{Float64};
    drive_bounds::Union{Float64,Vector{<:Union{Tuple{Float64,Float64},Float64}}} = 1.0,
    levels_per_transmon::Int = 3,
    subsystem_levels::AbstractVector{Int} = fill(levels_per_transmon, length(ωs)),
    lab_frame = false,
    subsystems::AbstractVector{Int} = 1:length(ωs),
    subsystem_drive_indices::AbstractVector{Int} = 1:length(ωs),
    kwargs...,
)
    n_subsystems = length(ωs)
    @assert length(δs) == n_subsystems
    @assert size(gs) == (n_subsystems, n_subsystems)

    # Convert drive_bounds to vector if scalar
    if drive_bounds isa Float64
        drive_bounds_vec = fill(drive_bounds, 2)
    else
        drive_bounds_vec = drive_bounds
    end

    systems = QuantumSystem[]

    for (i, (ω, δ, levels)) ∈ enumerate(zip(ωs, δs, subsystem_levels))
        if i ∈ subsystems
            sysᵢ = TransmonSystem(
                drive_bounds = drive_bounds_vec,
                levels = levels,
                ω = ω,
                δ = δ,
                lab_frame = lab_frame,
                drives = i ∈ subsystem_drive_indices,
            )
            push!(systems, sysᵢ)
        end
    end

    couplings = QuantumSystemCoupling[]

    for local_i = 1:(length(systems)-1)
        for local_j = (local_i+1):length(systems)
            global_i = subsystems[local_i]
            global_j = subsystems[local_j]
            push!(
                couplings,
                TransmonDipoleCoupling(
                    gs[global_i, global_j],
                    (local_i, local_j),
                    [sys.levels for sys in systems];
                    lab_frame = lab_frame,
                ),
            )
        end
    end

    levels = prod([sys.levels for sys in systems])
    H_drift = sum(c -> c.op, couplings; init = zeros(ComplexF64, levels, levels))
    return CompositeQuantumSystem(H_drift, systems, drive_bounds_vec)
end

# *************************************************************************** #

@testitem "TransmonSystem: default and custom parameters" begin
    sys = TransmonSystem()
    @test sys isa QuantumSystem
    @test sys.levels == 3
    @test sys.n_drives == 2

    sys2 = TransmonSystem(
        ω = 5.0,
        δ = 0.3,
        levels = 4,
        lab_frame = true,
        frame_ω = 0.0,
        lab_frame_type = :duffing,
        drives = false,
    )
    @test sys2.levels == 4
    @test sys2.n_drives == 0
end

@testitem "TransmonSystem: lab_frame_type variations" begin
    sys_duffing = TransmonSystem(lab_frame = true, lab_frame_type = :duffing)
    sys_quartic = TransmonSystem(lab_frame = true, lab_frame_type = :quartic)
    sys_cosine = TransmonSystem(lab_frame = true, lab_frame_type = :cosine)
    @test sys_duffing isa QuantumSystem
    @test sys_quartic isa QuantumSystem
    @test sys_cosine isa QuantumSystem
end

@testitem "TransmonSystem: error on invalid lab_frame_type" begin
    @test_throws AssertionError TransmonSystem(lab_frame = true, lab_frame_type = :invalid)
end

@testitem "TransmonDipoleCoupling: both constructors and frames" begin
    levels = [3, 3]
    g = 0.01

    c1 = TransmonDipoleCoupling(g, (1, 2), levels, lab_frame = false)
    c2 = TransmonDipoleCoupling(g, (1, 2), levels, lab_frame = true)
    @test c1 isa QuantumSystemCoupling
    @test c2 isa QuantumSystemCoupling

    sys1 = TransmonSystem(levels = 3)
    sys2 = TransmonSystem(levels = 3)
    c3 = TransmonDipoleCoupling(g, (1, 2), [sys1, sys2], lab_frame = false)
    @test c3 isa QuantumSystemCoupling
end

@testitem "MultiTransmonSystem: minimal and custom" begin
    using LinearAlgebra: norm

    ωs = [4.0, 4.1]
    δs = [0.2, 0.21]
    gs = [0.0 0.01; 0.01 0.0]

    comp = MultiTransmonSystem(ωs, δs, gs)
    @test comp isa CompositeQuantumSystem
    @test length(comp.subsystems) == 2
    @test !iszero(comp.H(zeros(comp.n_drives), 0.0))

    comp2 = MultiTransmonSystem(
        ωs,
        δs,
        gs;
        levels_per_transmon = 4,
        subsystem_levels = [4, 4],
        subsystems = [1],
        subsystem_drive_indices = [1],
    )
    @test comp2 isa CompositeQuantumSystem
    @test length(comp2.subsystems) == 1
    @test !isapprox(norm(comp2.H(zeros(comp2.n_drives), 0.0)), 0.0; atol = 1e-12)
end

@testitem "MultiTransmonSystem: edge cases" begin
    ωs = [4.0, 4.1, 4.2]
    δs = [0.2, 0.21, 0.22]
    gs = [0.0 0.01 0.02; 0.01 0.0 0.03; 0.02 0.03 0.0]
    # Only a subset of subsystems
    comp =
        MultiTransmonSystem(ωs, δs, gs; subsystems = [1, 3], subsystem_drive_indices = [3])
    @test comp isa CompositeQuantumSystem
    @test length(comp.subsystems) == 2
    # Only one drive
    @test comp.subsystems[1].n_drives == 0
    @test comp.subsystems[2].n_drives == 2
end

# *************************************************************************** #

@doc raw"""
    TransmonCavitySystem(;
        qubit_levels::Int=4,
        cavity_levels::Int=12,
        χ::Float64=2π * 32.8e-6,    # Dispersive shift (GHz)
        χ′::Float64=2π * 1.5e-9,    # Higher-order dispersive shift (GHz)
        K_c::Float64=2π * 1e-9 / 2, # Cavity self-Kerr (GHz)
        K_q::Float64=2π * 193e-3 / 2, # Qubit self-Kerr (GHz)
        drive_bounds::Vector{<:Union{Tuple{Float64, Float64}, Float64}}=fill(1.0, 4),
        multiply_by_2π::Bool=false, # Already in GHz with 2π factors
    ) -> QuantumSystem

Returns a `QuantumSystem` for a transmon qubit dispersively coupled to a cavity mode.

This system models circuit QED architectures where a transmon (artificial atom) is coupled 
to a microwave resonator in the dispersive regime, enabling quantum state manipulation and 
readout.

# Hamiltonian

```math
H = \tilde{\Delta} \hat{b}^\dagger \hat{b} 
    - \chi \hat{a}^\dagger \hat{a} \hat{b}^\dagger \hat{b}
    - \chi' \hat{b}^{\dagger 2} \hat{b}^2 \hat{a}^\dagger \hat{a}
    - K_q \hat{a}^{\dagger 2} \hat{a}^2
    - K_c \hat{b}^{\dagger 2} \hat{b}^2
```

where:
- ``\hat{a}``, ``\hat{a}^\dagger`` are the transmon annihilation/creation operators
- ``\hat{b}``, ``\hat{b}^\dagger`` are the cavity annihilation/creation operators  
- ``\tilde{\Delta} = \chi/2`` is the shifted cavity frequency
- ``\chi`` is the dispersive shift (qubit-cavity coupling)
- ``\chi'`` is a higher-order dispersive correction
- ``K_q``, ``K_c`` are self-Kerr nonlinearities

The drives are:
1. ``\hat{a}^\dagger + \hat{a}`` - Real transmon drive
2. ``i(\hat{a}^\dagger - \hat{a})`` - Imaginary transmon drive
3. ``\hat{b}^\dagger + \hat{b}`` - Real cavity drive
4. ``i(\hat{b}^\dagger - \hat{b})`` - Imaginary cavity drive

# Keyword Arguments
- `qubit_levels`: Number of transmon Fock states (typically 3-5)
- `cavity_levels`: Number of cavity Fock states (typically 10-20)
- `χ`: Dispersive shift in GHz (with 2π). Typical: ~2π × 30-50 kHz
- `χ′`: Higher-order dispersive shift in GHz. Typically small (~2π × 1-2 Hz)
- `K_c`: Cavity self-Kerr in GHz. Typically ~2π × 1 Hz
- `K_q`: Qubit self-Kerr (anharmonicity/2) in GHz. Typical: ~2π × 100-200 MHz
- `drive_bounds`: Control bounds for [Ωᵣ(qubit), Ωᵢ(qubit), αᵣ(cavity), αᵢ(cavity)]
- `multiply_by_2π`: Whether to multiply by 2π (default false, assuming parameters already include it)

# Example
```julia
# Standard cQED parameters
sys = TransmonCavitySystem(
    qubit_levels=4,
    cavity_levels=15,
    χ=2π * 32.8e-6,   # 32.8 kHz dispersive shift
    K_q=2π * 193e-3/2, # ~193 MHz anharmonicity
)
```

# References
- Blais et al., "Circuit quantum electrodynamics," Rev. Mod. Phys. 93, 025005 (2021)
- Koch et al., "Charge-insensitive qubit design derived from Cooper pair box," 
  Phys. Rev. A 76, 042319 (2007)
"""
function TransmonCavitySystem(;
    qubit_levels::Int = 4,
    cavity_levels::Int = 12,
    χ::Float64 = 2π * 32.8e-6,
    χ′::Float64 = 2π * 1.5e-9,
    K_c::Float64 = 2π * 1e-9 / 2,
    K_q::Float64 = 2π * 193e-3 / 2,
    drive_bounds::Vector{<:Union{Tuple{Float64,Float64},Float64}} = fill(1.0, 4),
    multiply_by_2π::Bool = false,
)
    # Cavity detuning in dispersive frame
    Δ̃ = χ / 2

    # Build operators: qubit ⊗ cavity
    â = lift_operator(annihilate(qubit_levels), 1, [qubit_levels, cavity_levels])
    b̂ = lift_operator(annihilate(cavity_levels), 2, [qubit_levels, cavity_levels])

    # Drift Hamiltonian
    H_drift =
        Δ̃ * b̂' * b̂ - χ * â' * â * b̂' * b̂ - χ′ * b̂'^2 * b̂^2 * â' * â -
        K_q * â'^2 * â^2 - K_c * b̂'^2 * b̂^2

    # Drive operators
    H_drives = [
        â' + â,              # Real transmon drive
        1.0im * (â' - â),    # Imaginary transmon drive
        b̂' + b̂,              # Real cavity drive
        1.0im * (b̂' - b̂),    # Imaginary cavity drive
    ]

    if multiply_by_2π
        H_drift *= 2π
        H_drives = [2π * H for H in H_drives]
    end

    return QuantumSystem(H_drift, H_drives, drive_bounds)
end

# *************************************************************************** #

@testitem "TransmonCavitySystem: basic construction" begin

    sys = TransmonCavitySystem()
    @test sys isa QuantumSystem
    @test sys.levels == 4 * 12  # qubit_levels × cavity_levels
    @test sys.n_drives == 4     # 2 qubit drives + 2 cavity drives

    # Check Hamiltonian is Hermitian
    using LinearAlgebra: ishermitian
    @test ishermitian(sys.H_drift)
    @test all(d -> ishermitian(d.H), sys.H_drives)
end

@testitem "TransmonCavitySystem: custom parameters" begin

    sys = TransmonCavitySystem(
        qubit_levels = 3,
        cavity_levels = 10,
        χ = 2π * 50e-6,
        K_q = 2π * 200e-3 / 2,
        drive_bounds = [0.5, 0.5, 1.0, 1.0],
    )
    @test sys.levels == 3 * 10
    @test sys.n_drives == 4
end
