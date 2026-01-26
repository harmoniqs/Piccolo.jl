# ============================================================================ #
# UnitaryTrajectory
# ============================================================================ #

"""
    UnitaryTrajectory{P<:AbstractPulse, S<:ODESolution, G} <: AbstractQuantumTrajectory{P}

Trajectory for unitary gate synthesis. The ODE solution is computed at construction.

# Fields
- `system::QuantumSystem`: The quantum system
- `pulse::P`: The control pulse (stores drive_name)
- `initial::Matrix{ComplexF64}`: Initial unitary (default: identity)
- `goal::G`: Target unitary operator (AbstractPiccoloOperator or Matrix)
- `solution::S`: Pre-computed ODE solution

# Callable
`traj(t)` returns the unitary at time `t` by interpolating the solution.

# Conversion to NamedTrajectory
Use `NamedTrajectory(traj, N)` or `NamedTrajectory(traj, times)` for optimization.
"""
mutable struct UnitaryTrajectory{P<:AbstractPulse, S<:ODESolution, G} <: AbstractQuantumTrajectory{P}
    system::QuantumSystem
    pulse::P
    initial::Matrix{ComplexF64}
    goal::G
    solution::S
end

"""
    UnitaryTrajectory(system, pulse, goal; initial=I, algorithm=MagnusGL4())

Create a unitary trajectory by solving the Schrödinger equation.

# Arguments
- `system::QuantumSystem`: The quantum system
- `pulse::AbstractPulse`: The control pulse
- `goal`: Target unitary (Matrix or AbstractPiccoloOperator)

# Keyword Arguments
- `initial`: Initial unitary (default: identity matrix)
- `algorithm`: ODE solver algorithm (default: MagnusGL4())
"""
function UnitaryTrajectory(
    system::QuantumSystem,
    pulse::AbstractPulse,
    goal::G;
    initial::AbstractMatrix{<:Number}=Matrix{ComplexF64}(I, system.levels, system.levels),
    algorithm=MagnusGL4(),
) where G
    @assert n_drives(pulse) == system.n_drives "Pulse has $(n_drives(pulse)) drives, system has $(system.n_drives)"
    
    U0 = Matrix{ComplexF64}(initial)
    times = collect(range(0.0, duration(pulse), length=101))
    prob = UnitaryOperatorODEProblem(system, pulse, times; U0=U0)
    sol = solve(prob, algorithm; saveat=times)
    
    return UnitaryTrajectory{typeof(pulse), typeof(sol), G}(system, pulse, U0, goal, sol)
end

"""
    UnitaryTrajectory(system, goal, T::Real; drive_name=:u, algorithm=MagnusGL4())

Convenience constructor that creates a zero pulse of duration T.

# Arguments
- `system::QuantumSystem`: The quantum system
- `goal`: Target unitary (Matrix or AbstractPiccoloOperator)
- `T::Real`: Duration of the pulse

# Keyword Arguments
- `drive_name::Symbol`: Name of the drive variable (default: `:u`)
- `algorithm`: ODE solver algorithm (default: MagnusGL4())
"""
function UnitaryTrajectory(
    system::QuantumSystem,
    goal::G,
    T::Real;
    drive_name::Symbol=:u,
    algorithm=MagnusGL4(),
) where G
    times = [0.0, T]
    controls = zeros(system.n_drives, 2)
    pulse = ZeroOrderPulse(controls, times; drive_name)
    return UnitaryTrajectory(system, pulse, goal; algorithm)
end

# Callable: sample solution at any time
(traj::UnitaryTrajectory)(t::Real) = traj.solution(t)

# ============================================================================ #
# Tests
# ============================================================================ #

@testitem "UnitaryTrajectory construction" begin
    using LinearAlgebra
    
    # Simple 2-level system
    system = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])
    
    # Create with duration
    T = 1.0
    X_gate = ComplexF64[0 1; 1 0]
    qtraj = UnitaryTrajectory(system, X_gate, T)
    
    @test qtraj isa UnitaryTrajectory
    @test qtraj.system === system
    @test qtraj.goal === X_gate
    @test qtraj.initial ≈ Matrix{ComplexF64}(I, 2, 2)
    
    # Create with explicit pulse
    times = [0.0, 0.5, 1.0]
    controls = 0.1 * randn(1, 3)
    pulse = ZeroOrderPulse(controls, times)
    qtraj2 = UnitaryTrajectory(system, pulse, X_gate)
    
    @test qtraj2 isa UnitaryTrajectory
    @test duration(qtraj2) ≈ 1.0
end

@testitem "UnitaryTrajectory callable" begin
    using LinearAlgebra
    
    system = QuantumSystem([PAULIS.X], [1.0])
    
    T = 1.0
    X_gate = ComplexF64[0 1; 1 0]
    qtraj = UnitaryTrajectory(system, X_gate, T)
    
    # Test at initial time
    U0 = qtraj(0.0)
    @test U0 ≈ Matrix{ComplexF64}(I, 2, 2)
    
    # Test at intermediate time
    U_mid = qtraj(0.5)
    @test U_mid isa Matrix{ComplexF64}
    @test size(U_mid) == (2, 2)
    
    # Test at final time
    U_final = qtraj(T)
    @test U_final isa Matrix{ComplexF64}
end

@testitem "UnitaryTrajectory fidelity" begin
    using LinearAlgebra
    
    # System that naturally implements X gate
    σx = ComplexF64[0 1; 1 0]
    system = QuantumSystem([σx], [1.0])
    
    # Create pulse that implements X gate: exp(-i π/2 σx) = -i σx
    T = π / 2
    times = [0.0, T]
    controls = ones(1, 2)  # Constant amplitude 1
    pulse = ZeroOrderPulse(controls, times)
    
    X_gate = ComplexF64[0 1; 1 0]
    qtraj = UnitaryTrajectory(system, pulse, X_gate)
    
    # Fidelity should be high
    fid = fidelity(qtraj)
    @test fid > 0.99
end

@testitem "UnitaryTrajectory with EmbeddedOperator goal" begin
    using LinearAlgebra
    
    # 3-level system with embedded 2-level gate
    H_drift = diagm(ComplexF64[1.0, 0.0, -1.0])
    H_drive = zeros(ComplexF64, 3, 3)
    H_drive[1,2] = H_drive[2,1] = 1.0
    system = QuantumSystem(H_drift, [H_drive], [1.0])
    
    # Embedded X gate on levels 1,2
    X_embedded = EmbeddedOperator(:X, [1, 2], 3)
    
    T = 1.0
    qtraj = UnitaryTrajectory(system, X_embedded, T)
    
    @test qtraj.goal === X_embedded
    @test qtraj.goal isa EmbeddedOperator
    
    # Fidelity with subspace
    fid = fidelity(qtraj; subspace=[1, 2])
    @test fid isa Real
end
