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

#
# Two sets of names are exposed:
#
#   `state_name(qtraj)` — user-facing symbol addressing the native-form
#     rollout states (unitaries, kets, density matrices). Use for lifted
#     accessors: `qtraj[state_name(qtraj)]` returns rollout states.
#
#   `isomorphism_state_name(qtraj)` — internal symbol used by the isomorphism
#     representation that `NamedTrajectory` and integrators work with. Use
#     this when indexing into an optimization trajectory.
#

"""
    state_name(::AbstractQuantumTrajectory)

User-facing state variable name — addresses the native rollout states.
Use with `qtraj[state_name(qtraj)]` or the lifted property accessor to pull
the evolved states (unitaries, kets, density matrices).

- `UnitaryTrajectory` → `:U`
- `KetTrajectory` → `:ψ`
- `MultiKetTrajectory` → `:ψ` (prefix for `:ψ1`, `:ψ2`, etc.)
- `DensityTrajectory` → `:ρ`
- `MultiDensityTrajectory` → `:ρ` (prefix for `:ρ1`, `:ρ2`, etc.)
"""
state_name(::UnitaryTrajectory) = :U
state_name(::KetTrajectory) = :ψ
state_name(::MultiKetTrajectory) = :ψ
state_name(::DensityTrajectory) = :ρ
state_name(::MultiDensityTrajectory) = :ρ

"""
    isomorphism_state_name(::AbstractQuantumTrajectory)

Internal state variable name — addresses the isomorphism representation used
by `NamedTrajectory` and integrators during optimization.

- `UnitaryTrajectory` → `:Ũ⃗`
- `KetTrajectory` → `:ψ̃`
- `MultiKetTrajectory` → `:ψ̃` (prefix for `:ψ̃1`, `:ψ̃2`, etc.)
- `DensityTrajectory` → `:ρ⃗̃`
- `MultiDensityTrajectory` → `:ρ⃗̃` (prefix for `:ρ⃗̃1`, `:ρ⃗̃2`, etc.)
"""
isomorphism_state_name(::UnitaryTrajectory) = :Ũ⃗
isomorphism_state_name(::KetTrajectory) = :ψ̃
isomorphism_state_name(::MultiKetTrajectory) = :ψ̃
isomorphism_state_name(::DensityTrajectory) = :ρ⃗̃
isomorphism_state_name(::MultiDensityTrajectory) = :ρ⃗̃

"""
    state_names(qtraj::MultiKetTrajectory)

Per-component user-facing names for a multi-ket trajectory (`:ψ1`, `:ψ2`, etc.).
"""
function state_names(qtraj::MultiKetTrajectory)
    prefix = state_name(qtraj)
    return [Symbol(prefix, i) for i = 1:length(qtraj.initials)]
end

"""
    state_names(qtraj::MultiDensityTrajectory)

Per-component user-facing names for a multi-density trajectory (`:ρ1`, `:ρ2`, etc.).
"""
function state_names(qtraj::MultiDensityTrajectory)
    prefix = state_name(qtraj)
    return [Symbol(prefix, i) for i = 1:length(qtraj.initials)]
end

"""
    isomorphism_state_names(qtraj::MultiKetTrajectory)

Per-component isomorphism names for a multi-ket trajectory (`:ψ̃1`, `:ψ̃2`, etc.) —
use these to index a `NamedTrajectory`.
"""
function isomorphism_state_names(qtraj::MultiKetTrajectory)
    prefix = isomorphism_state_name(qtraj)
    return [Symbol(prefix, i) for i = 1:length(qtraj.initials)]
end

"""
    isomorphism_state_names(qtraj::MultiDensityTrajectory)

Per-component isomorphism names for a multi-density trajectory (`:ρ⃗̃1`, `:ρ⃗̃2`, etc.).
"""
function isomorphism_state_names(qtraj::MultiDensityTrajectory)
    prefix = isomorphism_state_name(qtraj)
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
# Lifted Accessors — treat a quantum trajectory like a NamedTrajectory
#
# `qtraj.t`       → ODE save times (from `qtraj.solution.t`)
# `qtraj.sol`     → alias for the underlying ODE solution (`qtraj.solution`)
# `qtraj.U|ψ|ρ`   → rollout states (`qtraj.solution.u`), keyed by the
#                   trajectory's `state_name`
# `qtraj[:sym]`   → same as getproperty, with Symbol indexing
#
# For multi-state trajectories (MultiKet, MultiDensity) the underlying ODE
# solution is an `EnsembleSolution`; lifted access to the whole ensemble is
# intentionally left out — use `qtraj[i]` to pull the i-th sub-solution.
# ============================================================================ #

function Base.getproperty(qtraj::UnitaryTrajectory, name::Symbol)
    name === :t && return getfield(qtraj, :solution).t
    name === :sol && return getfield(qtraj, :solution)
    name === :U && return getfield(qtraj, :solution).u
    return getfield(qtraj, name)
end

function Base.getproperty(qtraj::KetTrajectory, name::Symbol)
    name === :t && return getfield(qtraj, :solution).t
    name === :sol && return getfield(qtraj, :solution)
    name === :ψ && return getfield(qtraj, :solution).u
    return getfield(qtraj, name)
end

function Base.getproperty(qtraj::DensityTrajectory, name::Symbol)
    name === :t && return getfield(qtraj, :solution).t
    name === :sol && return getfield(qtraj, :solution)
    name === :ρ && return getfield(qtraj, :solution).u
    return getfield(qtraj, name)
end

# For MultiKet / MultiDensity: :t and :sol only. State access goes through
# the existing integer indexing (`qtraj[i]`) or explicit `state_names(qtraj)`.
function Base.getproperty(qtraj::MultiKetTrajectory, name::Symbol)
    name === :t && return getfield(qtraj, :solution).u[1].t
    name === :sol && return getfield(qtraj, :solution)
    return getfield(qtraj, name)
end

function Base.getproperty(qtraj::MultiDensityTrajectory, name::Symbol)
    name === :t && return getfield(qtraj, :solution).u[1].t
    name === :sol && return getfield(qtraj, :solution)
    return getfield(qtraj, name)
end

# SamplingTrajectory delegates to its base_trajectory for lifted access
function Base.getproperty(qtraj::SamplingTrajectory, name::Symbol)
    name === :t && return getproperty(getfield(qtraj, :base_trajectory), :t)
    name === :sol && return getproperty(getfield(qtraj, :base_trajectory), :sol)
    # Delegate any state_name access (:U, :ψ, :ρ) to the base trajectory
    base = getfield(qtraj, :base_trajectory)
    if name === state_name(base)
        return getproperty(base, name)
    end
    return getfield(qtraj, name)
end

# Extend propertynames so tab completion reveals the lifted names
Base.propertynames(::UnitaryTrajectory, private::Bool = false) =
    (:system, :pulse, :initial, :goal, :solution, :t, :sol, :U)
Base.propertynames(::KetTrajectory, private::Bool = false) =
    (:system, :pulse, :initial, :goal, :solution, :t, :sol, :ψ)
Base.propertynames(::DensityTrajectory, private::Bool = false) =
    (:system, :pulse, :initial, :goal, :solution, :t, :sol, :ρ)
Base.propertynames(::MultiKetTrajectory, private::Bool = false) =
    (:system, :pulse, :initials, :goals, :weights, :solution, :t, :sol)
Base.propertynames(::MultiDensityTrajectory, private::Bool = false) =
    (:system, :pulse, :initials, :goals, :weights, :solution, :t, :sol)

"""
    Base.getindex(qtraj::AbstractQuantumTrajectory, name::Symbol)

Symbol indexing is equivalent to property access — `qtraj[:t]` == `qtraj.t`,
`qtraj[state_name(qtraj)]` returns the rollout states, and
`qtraj[isomorphism_state_name(qtraj)]` returns the isomorphism form at save
times.
"""
function Base.getindex(qtraj::AbstractQuantumTrajectory, name::Symbol)
    if name === isomorphism_state_name(qtraj)
        return _iso_form(qtraj)
    end
    return getproperty(qtraj, name)
end

# Isomorphism form of the rollout states, materialized at save times.
# Used by `qtraj[isomorphism_state_name(qtraj)]` for users who want the
# iso representation without building a NamedTrajectory.
_iso_form(qtraj::UnitaryTrajectory) = [operator_to_iso_vec(U) for U in qtraj.solution.u]
_iso_form(qtraj::KetTrajectory) = [ket_to_iso(ψ) for ψ in qtraj.solution.u]
_iso_form(qtraj::DensityTrajectory) =
    [density_to_compact_iso(ρ) for ρ in qtraj.solution.u]
_iso_form(qtraj::MultiKetTrajectory) =
    [[ket_to_iso(ψ) for ψ in sol.u] for sol in qtraj.solution.u]
_iso_form(qtraj::MultiDensityTrajectory) =
    [[density_to_compact_iso(ρ) for ρ in sol.u] for sol in qtraj.solution.u]
_iso_form(qtraj::SamplingTrajectory) = _iso_form(qtraj.base_trajectory)

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
    print(
        io,
        nameof(typeof(qtraj)),
        "($(sys.levels)-level $(nameof(typeof(sys))), ",
        "$(nameof(typeof(pulse))) with $(n_drives(pulse)) drives, ",
        "T = ", duration(qtraj),
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

    # Test state_name (user-facing) and isomorphism_state_name (internal)
    qtraj_u = UnitaryTrajectory(system, ComplexF64[0 1; 1 0], 1.0)
    @test state_name(qtraj_u) == :U
    @test isomorphism_state_name(qtraj_u) == :Ũ⃗

    qtraj_k = KetTrajectory(system, ComplexF64[1, 0], ComplexF64[0, 1], 1.0)
    @test state_name(qtraj_k) == :ψ
    @test isomorphism_state_name(qtraj_k) == :ψ̃

    qtraj_e = MultiKetTrajectory(system, [ComplexF64[1, 0]], [ComplexF64[0, 1]], 1.0)
    @test state_name(qtraj_e) == :ψ
    @test isomorphism_state_name(qtraj_e) == :ψ̃
    @test state_names(qtraj_e) == [:ψ1]
    @test isomorphism_state_names(qtraj_e) == [:ψ̃1]

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
    @test state_name(qtraj) == :ρ
    @test isomorphism_state_name(qtraj) == :ρ⃗̃
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

@testitem "Lifted accessors - UnitaryTrajectory" begin
    using LinearAlgebra

    system = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])
    X_gate = ComplexF64[0 1; 1 0]
    qtraj = UnitaryTrajectory(system, X_gate, 1.0)

    # Lifted time
    @test qtraj.t === qtraj.solution.t
    @test qtraj[:t] === qtraj.solution.t
    @test length(qtraj.t) == 101

    # Lifted ODE solution (sol alias)
    @test qtraj.sol === qtraj.solution
    @test qtraj[:sol] === qtraj.solution

    # Lifted state access via state_name
    @test state_name(qtraj) == :U
    @test qtraj.U === qtraj.solution.u
    @test qtraj[:U] === qtraj.solution.u
    @test qtraj[state_name(qtraj)] === qtraj.solution.u

    # Isomorphism form on demand
    iso = qtraj[isomorphism_state_name(qtraj)]
    @test iso isa Vector
    @test length(iso) == length(qtraj.t)
    # Each entry is a real iso-vec of length 2 * levels^2
    @test length(iso[1]) == 2 * system.levels^2

    # Struct field access still works
    @test qtraj.system === system
    @test qtraj.pulse isa AbstractPulse
    @test qtraj.goal === X_gate
    @test qtraj.initial ≈ Matrix{ComplexF64}(I, 2, 2)

    # propertynames exposes the lifted names for tab completion
    pnames = propertynames(qtraj)
    @test :t ∈ pnames
    @test :sol ∈ pnames
    @test :U ∈ pnames
end

@testitem "Lifted accessors - KetTrajectory" begin
    using LinearAlgebra

    system = QuantumSystem(PAULIS.Z, [PAULIS.X], [1.0])
    ψ0 = ComplexF64[1.0, 0.0]
    ψg = ComplexF64[0.0, 1.0]
    qtraj = KetTrajectory(system, ψ0, ψg, 1.0)

    @test state_name(qtraj) == :ψ
    @test qtraj.ψ === qtraj.solution.u
    @test qtraj[:ψ] === qtraj.solution.u
    @test qtraj[state_name(qtraj)] === qtraj.solution.u
    @test qtraj.t === qtraj.solution.t
    @test qtraj.sol === qtraj.solution

    # Isomorphism form (real iso vector per save time)
    iso = qtraj[isomorphism_state_name(qtraj)]
    @test length(iso) == length(qtraj.t)
    @test length(iso[1]) == 2 * system.levels
end

@testitem "Lifted accessors - DensityTrajectory" begin
    using LinearAlgebra

    L = ComplexF64[0.1 0.0; 0.0 0.0]
    system = OpenQuantumSystem(PAULIS.Z, [PAULIS.X], [1.0]; dissipation_operators = [L])

    ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]
    ρg = ComplexF64[0.0 0.0; 0.0 1.0]
    qtraj = DensityTrajectory(system, ρ0, ρg, 1.0)

    @test state_name(qtraj) == :ρ
    @test qtraj.ρ === qtraj.solution.u
    @test qtraj[:ρ] === qtraj.solution.u
    @test qtraj.t === qtraj.solution.t
    @test qtraj.sol === qtraj.solution
end
