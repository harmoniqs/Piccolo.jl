export IonChainSystem
export MolmerSorensenCoupling

@doc raw"""
    IonChainSystem(;
        N_ions::Int=2,
        ion_levels::Int=2,
        N_modes::Int=1,
        mode_levels::Int=10,
        ωq::Union{Float64, Vector{Float64}}=1.0,
        ωm::Union{Float64, Vector{Float64}}=0.1,
        η::Union{Float64, Matrix{Float64}}=0.1,
        lab_frame::Bool=false,
        frame_ω::Float64=lab_frame ? 0.0 : (ωq isa Vector ? ωq[1] : ωq),
        multiply_by_2π::Bool=true,
        drive_bounds::Vector{<:Union{Tuple{Float64, Float64}, Float64}}=fill(1.0, 2*N_ions),
    ) -> QuantumSystem

Returns a `QuantumSystem` object for a chain of trapped ions coupled via motional modes.

The system consists of `N_ions` ions, each with `ion_levels` internal states, coupled to 
`N_modes` shared motional modes with `mode_levels` Fock states each.

# Hamiltonian

In the lab frame:
```math
H = \sum_{i=1}^{N_{\text{ions}}} \omega_{q,i} \sigma_i^+ \sigma_i^- 
    + \sum_{m=1}^{N_{\text{modes}}} \omega_{m} a_m^\dagger a_m
    + \sum_{i,m} \eta_{i,m} (\sigma_i^+ + \sigma_i^-)(a_m + a_m^\dagger)
    + \sum_i \Omega_{x,i}(t) \sigma_i^x + \sum_i \Omega_{y,i}(t) \sigma_i^y
```

In the rotating frame at frequency `frame_ω`:
```math
H = \sum_{i=1}^{N_{\text{ions}}} (\omega_{q,i} - \omega_{\text{frame}}) \sigma_i^+ \sigma_i^- 
    + \sum_{m=1}^{N_{\text{modes}}} \omega_{m} a_m^\dagger a_m
    + \sum_{i,m} \eta_{i,m} (\sigma_i^+ + \sigma_i^-)(a_m + a_m^\dagger)
    + \sum_i \Omega_{x,i}(t) \sigma_i^x + \sum_i \Omega_{y,i}(t) \sigma_i^y
```

where:
- ``\sigma_i^+, \sigma_i^-`` are raising/lowering operators for ion `i`
- ``\sigma_i^x, \sigma_i^y`` are Pauli operators for ion `i`
- ``a_m, a_m^\dagger`` are annihilation/creation operators for mode `m`
- ``\omega_{q,i}`` is the transition frequency of ion `i`
- ``\omega_m`` is the frequency of motional mode `m`
- ``\eta_{i,m}`` is the Lamb-Dicke parameter coupling ion `i` to mode `m`

# Keyword Arguments
- `N_ions`: Number of ions in the chain
- `ion_levels`: Number of internal levels per ion (typically 2 for qubit)
- `N_modes`: Number of motional modes to include
- `mode_levels`: Number of Fock states per motional mode
- `ωq`: Ion transition frequency (or frequencies). Scalar or vector of length `N_ions`. In GHz.
- `ωm`: Motional mode frequency (or frequencies). Scalar or vector of length `N_modes`. In GHz.
- `η`: Lamb-Dicke parameter(s). Scalar (uniform coupling), or `N_ions × N_modes` matrix.
- `lab_frame`: If true, use lab frame Hamiltonian. If false, use rotating frame.
- `frame_ω`: Rotating frame frequency in GHz. Defaults to first ion frequency.
- `multiply_by_2π`: Whether to multiply Hamiltonian by 2π (default true, since frequencies are in GHz).
- `drive_bounds`: Control bounds. Vector of length `2*N_ions` for [Ωx₁, Ωy₁, Ωx₂, Ωy₂, ...].

# Example
```julia
# Two ions, one motional mode, Mølmer-Sørensen setup
sys = IonChainSystem(
    N_ions=2,
    N_modes=1,
    ωq=1.0,      # 1 GHz qubit frequency
    ωm=0.1,      # 100 MHz mode frequency  
    η=0.1,       # Lamb-Dicke parameter
    mode_levels=5,
)
```

# References
- Sørensen, A. & Mølmer, K. "Quantum computation with ions in thermal motion." 
  Phys. Rev. Lett. 82, 1971 (1999).
- Sørensen, A. & Mølmer, K. "Entanglement and quantum computation with ions in thermal motion."
  Phys. Rev. A 62, 022311 (2000).
"""
function IonChainSystem(;
    N_ions::Int = 2,
    ion_levels::Int = 2,
    N_modes::Int = 1,
    mode_levels::Int = 10,
    ωq::Union{Float64,Vector{Float64}} = 1.0,
    ωm::Union{Float64,Vector{Float64}} = 0.1,
    η::Union{Float64,Matrix{Float64}} = 0.1,
    lab_frame::Bool = false,
    frame_ω::Float64 = lab_frame ? 0.0 : (ωq isa Vector ? ωq[1] : ωq),
    multiply_by_2π::Bool = true,
    drive_bounds::Vector{<:Union{Tuple{Float64,Float64},Float64}} = fill(1.0, 2 * N_ions),
)
    # Convert scalar parameters to vectors
    ωq_vec = ωq isa Vector ? ωq : fill(ωq, N_ions)
    ωm_vec = ωm isa Vector ? ωm : fill(ωm, N_modes)

    @assert length(ωq_vec) == N_ions "ωq must be scalar or vector of length N_ions"
    @assert length(ωm_vec) == N_modes "ωm must be scalar or vector of length N_modes"
    @assert length(drive_bounds) == 2 * N_ions "drive_bounds must have length 2*N_ions"

    # Convert η to matrix
    if η isa Float64
        η_mat = fill(η, N_ions, N_modes)
    else
        η_mat = η
        @assert size(η_mat) == (N_ions, N_modes) "η matrix must be N_ions × N_modes"
    end

    # Total Hilbert space dimension: ion_levels^N_ions × mode_levels^N_modes
    # Subsystem structure: [ion_1, ion_2, ..., ion_N, mode_1, mode_2, ..., mode_M]
    subsystem_levels = vcat(fill(ion_levels, N_ions), fill(mode_levels, N_modes))
    ion_dim = ion_levels^N_ions
    mode_dim = mode_levels^N_modes
    total_dim = ion_dim * mode_dim

    # Pauli operators for ion_levels=2 case
    if ion_levels == 2
        σ_plus = ComplexF64[0 1; 0 0]
        σ_minus = ComplexF64[0 0; 1 0]
        σ_x = ComplexF64[0 1; 1 0]
        σ_y = ComplexF64[0 -im; im 0]
    else
        # For multi-level ions, use general raising/lowering between |0⟩ and |1⟩
        σ_plus = zeros(ComplexF64, ion_levels, ion_levels)
        σ_plus[2, 1] = 1.0
        σ_minus = zeros(ComplexF64, ion_levels, ion_levels)
        σ_minus[1, 2] = 1.0
        σ_x = σ_plus + σ_minus
        σ_y = -1.0im * (σ_plus - σ_minus)
    end

    # Build drift Hamiltonian
    H_drift = zeros(ComplexF64, total_dim, total_dim)

    # Ion internal Hamiltonian: Σᵢ ωq,i σᵢ⁺ σᵢ⁻
    for i = 1:N_ions
        detuning = lab_frame ? ωq_vec[i] : (ωq_vec[i] - frame_ω)
        H_drift += detuning * lift_operator(σ_plus' * σ_plus, i, subsystem_levels)
    end

    # Motional mode Hamiltonian: Σₘ ωₘ aₘ⁺ aₘ
    for m = 1:N_modes
        a_m = annihilate(mode_levels)
        H_drift += ωm_vec[m] * lift_operator(a_m' * a_m, N_ions + m, subsystem_levels)
    end

    # Lamb-Dicke coupling: Σᵢ,ₘ ηᵢ,ₘ (σᵢ⁺ + σᵢ⁻)(aₘ + aₘ⁺)
    for i = 1:N_ions
        for m = 1:N_modes
            if abs(η_mat[i, m]) > 1e-12  # Skip negligible couplings
                σ_i_x = lift_operator(σ_x, i, subsystem_levels)
                a_m = annihilate(mode_levels)
                x_m = lift_operator(a_m + a_m', N_ions + m, subsystem_levels)
                H_drift += η_mat[i, m] * σ_i_x * x_m
            end
        end
    end

    # Drive operators: Ωx,i σᵢˣ and Ωy,i σᵢʸ for each ion
    H_drives = Matrix{ComplexF64}[]
    for i = 1:N_ions
        push!(H_drives, lift_operator(σ_x, i, subsystem_levels))  # X drive on ion i
        push!(H_drives, lift_operator(σ_y, i, subsystem_levels))  # Y drive on ion i
    end

    # Apply 2π factor if requested
    if multiply_by_2π
        H_drift *= 2π
        H_drives = [2π * H for H in H_drives]
    end

    return QuantumSystem(H_drift, H_drives, drive_bounds)
end

@doc raw"""
    MolmerSorensenCoupling(
        N_ions::Int,
        N_modes::Int,
        ion_levels::Int,
        mode_levels::Int,
        η::Union{Float64, Matrix{Float64}},
        ωm::Union{Float64, Vector{Float64}},
    ) -> Matrix{ComplexF64}

Returns the Mølmer-Sørensen coupling term for an ion chain, which mediates 
effective ion-ion interactions via the motional modes.

In the Lamb-Dicke regime with appropriate drive detunings, the effective Hamiltonian is:
```math
H_{\text{MS}} = \sum_{i<j} J_{ij} \sigma_i^x \sigma_j^x
```

where the coupling strength is:
```math
J_{ij} = \sum_m \frac{\eta_{i,m} \eta_{j,m} \Omega_i \Omega_j}{4 \Delta_m}
```

with ``\Delta_m`` being the detuning from mode ``m``.

# Arguments
- `N_ions`: Number of ions
- `N_modes`: Number of motional modes
- `ion_levels`: Internal levels per ion
- `mode_levels`: Fock states per mode
- `η`: Lamb-Dicke parameters (scalar or N_ions × N_modes matrix)
- `ωm`: Mode frequencies (scalar or vector)

# Returns
Matrix representing the σˣᵢ σˣⱼ interaction for use in building MS gates.
"""
function MolmerSorensenCoupling(
    N_ions::Int,
    N_modes::Int,
    ion_levels::Int,
    mode_levels::Int,
    η::Union{Float64,Matrix{Float64}},
    ωm::Union{Float64,Vector{Float64}},
)
    # Convert to matrix form
    η_mat = η isa Float64 ? fill(η, N_ions, N_modes) : η
    ωm_vec = ωm isa Vector ? ωm : fill(ωm, N_modes)

    # Subsystem structure: [ion_1, ion_2, ..., ion_N, mode_1, mode_2, ..., mode_M]
    subsystem_levels = vcat(fill(ion_levels, N_ions), fill(mode_levels, N_modes))
    ion_dim = ion_levels^N_ions
    mode_dim = mode_levels^N_modes
    total_dim = ion_dim * mode_dim

    # Pauli-X for two-level systems
    if ion_levels == 2
        σ_x = ComplexF64[0 1; 1 0]
    else
        σ_x = zeros(ComplexF64, ion_levels, ion_levels)
        σ_x[1, 2] = σ_x[2, 1] = 1.0
    end

    # Build MS interaction: Σᵢ<ⱼ σᵢˣ σⱼˣ
    H_MS = zeros(ComplexF64, total_dim, total_dim)
    for i = 1:(N_ions-1)
        for j = (i+1):N_ions
            σ_x_i = lift_operator(σ_x, i, subsystem_levels)
            σ_x_j = lift_operator(σ_x, j, subsystem_levels)
            H_MS += σ_x_i * σ_x_j
        end
    end

    return H_MS
end

# *************************************************************************** #

@testitem "IonChainSystem: basic construction" begin

    # Minimal system: 2 ions, 1 mode
    sys = IonChainSystem(N_ions = 2, N_modes = 1, mode_levels = 5)
    @test sys isa QuantumSystem
    @test sys.n_drives == 4  # 2 drives (X, Y) per ion
    @test sys.levels == 2^2 * 5  # 2 ions × 5 mode levels

    # Single ion with mode
    sys_single = IonChainSystem(N_ions = 1, N_modes = 1, mode_levels = 3)
    @test sys_single.levels == 2 * 3
    @test sys_single.n_drives == 2
end

@testitem "IonChainSystem: parameter variations" begin

    # Vector frequencies
    sys = IonChainSystem(
        N_ions = 3,
        N_modes = 2,
        ωq = [1.0, 1.01, 1.02],
        ωm = [0.1, 0.11],
        η = 0.05,
        mode_levels = 3,
    )
    @test sys isa QuantumSystem
    @test sys.n_drives == 6  # 3 ions × 2 drives each

    # Matrix Lamb-Dicke parameters
    η_mat = [0.1 0.05; 0.1 0.05]
    sys2 = IonChainSystem(N_ions = 2, N_modes = 2, η = η_mat, mode_levels = 3)
    @test sys2 isa QuantumSystem
end

@testitem "IonChainSystem: lab frame vs rotating frame" begin

    sys_lab = IonChainSystem(N_ions = 2, N_modes = 1, lab_frame = true, mode_levels = 3)
    sys_rot = IonChainSystem(
        N_ions = 2,
        N_modes = 1,
        lab_frame = false,
        frame_ω = 1.0,
        mode_levels = 3,
    )

    @test sys_lab isa QuantumSystem
    @test sys_rot isa QuantumSystem
    @test sys_lab.H_drift != sys_rot.H_drift  # Different frames
end

@testitem "IonChainSystem: multi-level ions" begin

    sys = IonChainSystem(N_ions = 2, ion_levels = 3, N_modes = 1, mode_levels = 3)
    @test sys.levels == 3^2 * 3  # 3-level ions
    @test sys isa QuantumSystem
end

@testitem "MolmerSorensenCoupling: basic" begin
    using LinearAlgebra: ishermitian

    H_MS = MolmerSorensenCoupling(2, 1, 2, 5, 0.1, 0.1)
    @test size(H_MS) == (2^2 * 5, 2^2 * 5)
    @test ishermitian(H_MS)

    # Three ions
    H_MS3 = MolmerSorensenCoupling(3, 1, 2, 3, 0.1, 0.1)
    @test size(H_MS3) == (2^3 * 3, 2^3 * 3)
end
