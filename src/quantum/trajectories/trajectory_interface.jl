# ============================================================================ #
# AbstractQuantumTrajectory Interface
#
# Shared getters, name accessors, and display methods for all concrete
# subtypes of AbstractQuantumTrajectory.
#
# NOTE: Rollout and fidelity methods live in `rollouts_extensions.jl` — keep
# them there so users can find the user-facing entry points easily.
# ============================================================================ #

# ============================================================================ #
# Field Getters
# ============================================================================ #

"""
    get_system(qtraj)

Get the quantum system from a trajectory.
"""
get_system(qtraj::AbstractQuantumTrajectory) = qtraj.system

"""
    get_pulse(qtraj)

Get the control pulse from a trajectory.
"""
get_pulse(qtraj::AbstractQuantumTrajectory) = qtraj.pulse

"""
    get_initial(qtraj)

Get the initial state/operator from a trajectory.
"""
get_initial(qtraj::UnitaryTrajectory) = qtraj.initial
get_initial(qtraj::KetTrajectory) = qtraj.initial
get_initial(qtraj::MultiKetTrajectory) = qtraj.initials
get_initial(qtraj::DensityTrajectory) = qtraj.initial
get_initial(qtraj::MultiDensityTrajectory) = qtraj.initials

"""
    get_goal(qtraj)

Get the goal state/operator from a trajectory.
"""
get_goal(qtraj::UnitaryTrajectory) = qtraj.goal
get_goal(qtraj::KetTrajectory) = qtraj.goal
get_goal(qtraj::MultiKetTrajectory) = qtraj.goals
get_goal(qtraj::DensityTrajectory) = qtraj.goal
get_goal(qtraj::MultiDensityTrajectory) = qtraj.goals

"""
    get_solution(qtraj)

Get the ODE solution from a trajectory.
"""
get_solution(qtraj::AbstractQuantumTrajectory) = qtraj.solution

# ============================================================================ #
# Fixed Name Accessors (for NamedTrajectory conversion)
# ============================================================================ #

"""
    state_name(::AbstractQuantumTrajectory)

Get the fixed state variable name for a trajectory type.
- `UnitaryTrajectory` → `:Ũ⃗`
- `KetTrajectory` → `:ψ̃`
- `MultiKetTrajectory` → `:ψ̃` (with index appended: `:ψ̃1`, `:ψ̃2`, etc.)
- `DensityTrajectory` → `:ρ⃗̃`
"""
state_name(::UnitaryTrajectory) = :Ũ⃗
state_name(::KetTrajectory) = :ψ̃
state_name(::MultiKetTrajectory) = :ψ̃  # prefix for :ψ̃1, :ψ̃2, etc.
state_name(::DensityTrajectory) = :ρ⃗̃
state_name(::MultiDensityTrajectory) = :ρ⃗̃  # prefix for :ρ⃗̃1, :ρ⃗̃2, etc.

"""
    state_names(qtraj::MultiKetTrajectory)

Get all state names for an ensemble trajectory (`:ψ̃1`, `:ψ̃2`, etc.)
"""
function state_names(qtraj::MultiKetTrajectory)
    prefix = state_name(qtraj)
    return [Symbol(prefix, i) for i = 1:length(qtraj.initials)]
end

"""
    state_names(qtraj::MultiDensityTrajectory)

Get all state names for a multi-density trajectory (`:ρ⃗̃1`, `:ρ⃗̃2`, etc.)
"""
function state_names(qtraj::MultiDensityTrajectory)
    prefix = state_name(qtraj)
    return [Symbol(prefix, i) for i = 1:length(qtraj.initials)]
end

"""
    drive_name(qtraj::AbstractQuantumTrajectory)

Get the drive/control variable name from the trajectory's pulse.
"""
drive_name(qtraj::AbstractQuantumTrajectory) = drive_name(qtraj.pulse)

"""
    time_name(::AbstractQuantumTrajectory)

Get the time variable name (always `:t`).
"""
time_name(::AbstractQuantumTrajectory) = :t

"""
    timestep_name(::AbstractQuantumTrajectory)

Get the timestep variable name (always `:Δt`).
"""
timestep_name(::AbstractQuantumTrajectory) = :Δt

"""
    duration(qtraj)

Get the duration of a trajectory (from its pulse).
"""
duration(qtraj::AbstractQuantumTrajectory) = duration(qtraj.pulse)

# ============================================================================ #
# Display
# ============================================================================ #

"""
    Base.summary(io::IO, qtraj::AbstractQuantumTrajectory)

Compact one-line summary for quantum trajectories.
"""
function Base.summary(io::IO, qtraj::AbstractQuantumTrajectory)
    sys = get_system(qtraj)
    pulse = get_pulse(qtraj)
    nd = n_drives(pulse)
    drives_word = nd == 1 ? "drive" : "drives"
    print(
        io,
        nameof(typeof(qtraj)),
        "($(sys.levels)-level $(nameof(typeof(sys))), ",
        "$(nameof(typeof(pulse))) with $nd $drives_word, ",
        "T = ",
        duration(qtraj),
        ")",
    )
end

"""
    Base.show(io::IO, qtraj::AbstractQuantumTrajectory)

Single-line display for compact contexts (arrays, tuples).
"""
Base.show(io::IO, qtraj::AbstractQuantumTrajectory) = summary(io, qtraj)

"""
    Base.show(io::IO, ::MIME"text/plain", qtraj::AbstractQuantumTrajectory)

Richer plain-text display for REPL and notebooks.
"""
function Base.show(io::IO, ::MIME"text/plain", qtraj::AbstractQuantumTrajectory)
    sys = get_system(qtraj)
    pulse = get_pulse(qtraj)
    println(io, nameof(typeof(qtraj)))
    println(io, "  system: $(sys.levels)-level $(nameof(typeof(sys)))")
    println(io, "  drives: $(nameof(typeof(pulse))) with $(n_drives(pulse)) drives")
    print(io, "  duration: ", duration(qtraj))
end

# ============================================================================ #
# Tests
# ============================================================================ #

@testitem "Common interface - getters" begin
    using LinearAlgebra

    # UnitaryTrajectory
    system = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])

    X_gate = ComplexF64[0 1; 1 0]
    qtraj = UnitaryTrajectory(system, X_gate, 1.0)

    @test get_system(qtraj) === system
    @test get_pulse(qtraj) isa AbstractPulse
    @test get_initial(qtraj) ≈ Matrix{ComplexF64}(I, 2, 2)
    @test get_goal(qtraj) === X_gate
    @test duration(qtraj) ≈ 1.0

    # KetTrajectory
    ψ0 = ComplexF64[1.0, 0.0]
    ψg = ComplexF64[0.0, 1.0]
    qtraj_ket = KetTrajectory(system, ψ0, ψg, 1.0)

    @test get_system(qtraj_ket) === system
    @test get_initial(qtraj_ket) ≈ ψ0
    @test get_goal(qtraj_ket) ≈ ψg

    # MultiKetTrajectory
    initials = [ψ0, ψg]
    goals = [ψg, ψ0]
    qtraj_ens = MultiKetTrajectory(system, initials, goals, 1.0)

    @test get_system(qtraj_ens) === system
    @test get_initial(qtraj_ens) == qtraj_ens.initials
    @test get_goal(qtraj_ens) == qtraj_ens.goals
end

@testitem "Common interface - name accessors" begin
    using LinearAlgebra

    system = QuantumSystem([PAULIS.X], [1.0])

    # Test state_name for each trajectory type
    qtraj_u = UnitaryTrajectory(system, ComplexF64[0 1; 1 0], 1.0)
    @test state_name(qtraj_u) == :Ũ⃗

    qtraj_k = KetTrajectory(system, ComplexF64[1, 0], ComplexF64[0, 1], 1.0)
    @test state_name(qtraj_k) == :ψ̃

    qtraj_e = MultiKetTrajectory(system, [ComplexF64[1, 0]], [ComplexF64[0, 1]], 1.0)
    @test state_name(qtraj_e) == :ψ̃
    @test state_names(qtraj_e) == [:ψ̃1]

    # Test drive_name propagation
    times = [0.0, 1.0]
    pulse = ZeroOrderPulse(zeros(1, 2), times; drive_name = :control)
    qtraj_custom = UnitaryTrajectory(system, pulse, ComplexF64[0 1; 1 0])
    @test drive_name(qtraj_custom) == :control

    # Test time_name and timestep_name (always fixed)
    @test time_name(qtraj_u) == :t
    @test timestep_name(qtraj_u) == :Δt
end

@testitem "Interface - DensityTrajectory getters" begin
    using LinearAlgebra

    L = ComplexF64[0.1 0.0; 0.0 0.0]
    system = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])

    ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρg = ComplexF64[0.0 0.0; 0.0 1.0]
    qtraj = DensityTrajectory(system, ρ0, ρg, 1.0)

    @test get_system(qtraj) === system
    @test get_initial(qtraj) ≈ ρ0
    @test get_goal(qtraj) ≈ ρg
    @test state_name(qtraj) == :ρ⃗̃
    @test duration(qtraj) ≈ 1.0
end

@testitem "Display - show and summary" begin
    system = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])
    X_gate = ComplexF64[0 1; 1 0]
    qtraj = UnitaryTrajectory(system, X_gate, 1.0)

    # Compact summary (single line)
    buf = IOBuffer()
    summary(buf, qtraj)
    s = String(take!(buf))
    @test occursin("UnitaryTrajectory", s)
    @test occursin("QuantumSystem", s)
    @test occursin("drives", s)
    @test occursin("T =", s)

    # Plain-text (multi-line) display
    buf = IOBuffer()
    show(buf, MIME("text/plain"), qtraj)
    s = String(take!(buf))
    @test occursin("UnitaryTrajectory", s)
    @test occursin("system:", s)
    @test occursin("drives:", s)
    @test occursin("duration:", s)

    # Uses nameof(typeof(...)) — no module-qualified types
    @test !occursin("Piccolo.", s)
    @test !occursin("QuantumTrajectories.", s)
end
