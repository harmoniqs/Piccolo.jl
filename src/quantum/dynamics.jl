module Rollouts

"""
Rollouts of quantum systems using SciML's DifferentialEquations.jl.

# Integration Algorithms

The default integrator is `MagnusAdapt4()` — a 4th-order adaptive Magnus method that:
- Preserves unitarity (Lie group structure)
- Adapts step size based on local error estimation
- Uses `abstol`/`reltol` instead of fixed point counts

`MagnusGL4()` is still available for fixed-step integration (pass `algorithm=MagnusGL4()`).
`Tsit5()` is used by default for `DensityTrajectory` (open systems).

# Two Ways to Check Fidelity

## 1. Fast fidelity from quantum trajectory (O(1) - recommended)
```julia
qtraj = UnitaryTrajectory(system, pulse, goal)
fid = fidelity(qtraj)  # Uses pre-computed ODE solution
```
**Use this for:** Post-optimization fidelity checks, analysis

## 2. Validate discrete controls with interpolation (O(solve))
```julia
traj = get_trajectory(qcp)  # NamedTrajectory with discrete controls
fid = rollout_fidelity(traj, system; interpolation=:cubic)
```
**Use this for:** Testing interpolation methods, validation against discrete trajectory

# Rolling Out New Pulses

```julia
# Roll out a new pulse through the system (creates new trajectory)
qtraj_new = rollout(qtraj, new_pulse)

# In-place update after optimization
pulse = extract_pulse(qtraj, get_trajectory(qcp))
rollout!(qtraj, pulse)
```

# Provided Functions

Domain-specific language for quantum system rollouts:
- `KetODEProblem`: Ket rollouts
- `UnitaryODEProblem`: Unitary rollouts
- `DensityODEProblem`: Density matrix rollouts (open systems)

SciML MatrixOperator versions for Lie group integrators (e.g., Magnus expansion):
- `KetOperatorODEProblem`
- `UnitaryOperatorODEProblem`

Fidelity and rollout methods:
- `fidelity(qtraj)`: Fast lookup from quantum trajectory
- `rollout(qtraj, pulse; kwargs...)`: Roll out a new pulse
- `rollout_fidelity(traj, sys; kwargs...)`: Validate discrete NamedTrajectory controls

"""

export fidelity
export unitary_fidelity
export rollout
export rollout!
export rollout_fidelity
export ket_rollout
export ket_rollout_fidelity
export unitary_rollout
export unitary_rollout_fidelity
export open_rollout
export open_rollout_fidelity
export update_global_params!

export KetODEProblem
export KetOperatorODEProblem
export UnitaryODEProblem
export UnitaryOperatorODEProblem
export DensityODEProblem

using LinearAlgebra
using NamedTrajectories
using DataInterpolations
using OrdinaryDiffEqLinear
using SymbolicIndexingInterface
const SII = SymbolicIndexingInterface
using TestItems

using ..Isomorphisms
using ..QuantumSystems

# ------------------------------------------------------------ #
# Rollout functions (stubs - extended in quantum_trajectories)
# ------------------------------------------------------------ #

"""
    rollout(qtraj, args...; kwargs...)

Roll out a quantum trajectory with new pulse or ODE parameters.
Extended in quantum_trajectories module for specific trajectory types.
"""
function rollout end

"""
    rollout!(qtraj, args...; kwargs...)

In-place rollout of quantum trajectory with new pulse or ODE parameters.
Extended in quantum_trajectories module for specific trajectory types.
"""
function rollout! end

"""
    update_global_params!(qtraj, traj::NamedTrajectory)

Update the global parameters in the quantum trajectory's system with the optimized
values from the NamedTrajectory after optimization. Handles immutable QuantumSystem
by reconstructing with updated global_params NamedTuple.
"""
function update_global_params!(qtraj, traj)
    # Check if trajectory has global components
    if !hasfield(typeof(traj), :global_components) || isempty(traj.global_components)
        return nothing
    end

    # Extract optimized global values preserving original key order from the system.
    # Using the original key order is critical: NamedTuples with different key orders are
    # different types in Julia, so reordering keys would cause type mismatches when rollout!
    # tries to assign the new solution back to the parametrically-typed trajectory.
    sys = qtraj.system
    original_keys = keys(sys.global_params)

    new_values = map(original_keys) do name
        indices = traj.global_components[name]
        if length(indices) == 1
            traj.global_data[indices[1]]
        else
            # Multi-dimensional globals (future support)
            traj.global_data[indices]
        end
    end

    new_global_params = NamedTuple{original_keys}(new_values)

    new_sys = _reconstruct_system(sys, new_global_params)

    # Update the quantum trajectory's system field (using internal method)
    _update_system!(qtraj, new_sys)

    return nothing
end

function _reconstruct_system(sys::QuantumSystem, new_global_params::NamedTuple)
    # Always use the inner constructor to preserve the exact H/G closure types.
    # H/G closures do not capture global_params directly (they are appended to u
    # by the integrator), so reusing sys.H and sys.G is always correct.
    return QuantumSystem(
        sys.H,
        sys.G,
        sys.H_drift,
        sys.H_drives,
        sys.drive_bounds,
        sys.n_drives,
        sys.levels,
        sys.time_dependent,
        new_global_params,
    )
end

function _reconstruct_system(sys::OpenQuantumSystem, new_global_params::NamedTuple)
    return OpenQuantumSystem(
        sys.H,
        sys.𝒢,
        sys.H_drift,
        sys.H_drives,
        sys.drive_bounds,
        sys.n_drives,
        sys.levels,
        sys.dissipation_operators,
        sys.time_dependent,
        new_global_params,
    )
end

"""
    _update_system!(qtraj, sys::QuantumSystem)

Internal method to update the system field in a quantum trajectory.
Extended in quantum_trajectories module for specific trajectory types.
"""
function _update_system! end

"""
    extract_globals(traj::NamedTrajectory, names::Vector{Symbol}=Symbol[])

Extract global variables from trajectory as a NamedTuple for easy access.
If names is empty, extracts all global variables.

# Example
```julia
traj = NamedTrajectory(...; global_data=[0.5, 1.0], global_components=(δ=1:1, Ω=2:2))
g = extract_globals(traj)  # (δ = 0.5, Ω = 1.0)
```
"""
function extract_globals(traj, names::Vector{Symbol} = Symbol[])
    # Check if trajectory has global components
    if !hasfield(typeof(traj), :global_components) || isempty(traj.global_components)
        return NamedTuple()
    end

    if isempty(names)
        names = collect(keys(traj.global_components))
    end

    global_dict = Dict{Symbol,Any}()
    for name in names
        indices = traj.global_components[name]
        if length(indices) == 1
            global_dict[name] = traj.global_data[indices[1]]
        else
            # Multi-dimensional globals - return vector
            global_dict[name] = traj.global_data[indices]
        end
    end

    return NamedTuple(global_dict)
end

# ------------------------------------------------------------ #
# Fidelity
# ------------------------------------------------------------ #

"""
    fidelity(ψ::AbstractVector{<:Number}, ψ_goal::AbstractVector{<:Number})

Calculate the fidelity between two quantum states `ψ` and `ψ_goal`.
"""
function fidelity(ψ::AbstractVector{<:Number}, ψ_goal::AbstractVector{<:Number})
    return abs2(ψ'ψ_goal)
end

"""
    fidelity(ρ::AbstractMatrix{<:Number}, ρ_goal::AbstractMatrix{<:Number})

Calculate the fidelity between two density matrices `ρ` and `ρ_goal`.
"""
function fidelity(ρ::AbstractMatrix{<:Number}, ρ_goal::AbstractMatrix{<:Number})
    return real(tr(ρ * ρ_goal))
end

"""
    unitary_fidelity(U::AbstractMatrix{<:Number}, U_goal::AbstractMatrix{<:Number})

Calculate the fidelity between unitary operators `U` and `U_goal` in the `subspace`.
"""
function unitary_fidelity(
    U::AbstractMatrix{<:Number},
    U_goal::AbstractMatrix{<:Number};
    subspace::AbstractVector{Int} = axes(U, 1),
)
    U = U[subspace, subspace]
    U_goal = U_goal[subspace, subspace]
    N = size(U, 1)
    return abs2(tr(U' * U_goal)) / N^2
end

# ------------------------------------------------------------ #
# DSL for Piccolo
# ------------------------------------------------------------ #

const SymbolIndex =
    Union{Int,AbstractVector{Int},CartesianIndex{N} where N,CartesianIndices{N} where N}

function _index(name::Symbol, n::Int)
    index = Dict{Symbol,SymbolIndex}()
    for i = 1:n
        index[Symbol(name, :_, i)] = i
    end
    index[name] = 1:n
    return index
end

function _index(name::Symbol, n1::Int, n2::Int)
    idx = Dict{Symbol,SymbolIndex}()
    for i = 1:n1, j = 1:n2
        idx[Symbol(name, :_, i, :_, j)] = CartesianIndex(i, j)
    end
    # block symbol: preserves matrix shape
    idx[name] = CartesianIndices((n1, n2))
    return idx
end

struct PiccoloRolloutSystem{T1<:SymbolIndex}
    state_index::Dict{Symbol,T1}
    t::Symbol
    defaults::Dict{Symbol,Float64}
end

function PiccoloRolloutSystem(
    state::Pair{Symbol,Int},
    timestep_name::Symbol = :t,
    defaults::Dict{Symbol,Float64} = Dict{Symbol,Float64}(),
)
    state_name, n_state = state
    state_index = _index(state_name, n_state)
    return PiccoloRolloutSystem(state_index, timestep_name, defaults)
end

function PiccoloRolloutSystem(
    state::Pair{Symbol,Tuple{Int,Int}},
    timestep_name::Symbol = :t,
    defaults::Dict{Symbol,Float64} = Dict{Symbol,Float64}(),
)
    state_name, (n1, n2) = state
    state_index = _index(state_name, n1, n2)
    return PiccoloRolloutSystem(state_index, timestep_name, defaults)
end

function _construct_operator(sys::AbstractQuantumSystem, u::F) where {F}
    A0 = zeros(ComplexF64, sys.levels, sys.levels)

    # Build u_vec function that appends globals to controls
    if length(sys.global_params) > 0
        global_vals = collect(values(sys.global_params))
        u_vec = t -> vcat(u(t), global_vals)
    else
        u_vec = u  # No globals, use controls directly
    end

    function update!(A, x, p, t)
        Ht = collect(sys.H(u_vec(t), t))
        @. A = -im * Ht
        return nothing
    end
    return SciMLOperators.MatrixOperator(A0; (update_func!) = update!)
end

function _construct_rhs(sys::AbstractQuantumSystem, u::F) where {F}
    # Build u_vec function that appends globals to controls
    if length(sys.global_params) > 0
        global_vals = collect(values(sys.global_params))
        u_vec = t -> vcat(u(t), global_vals)
    else
        u_vec = u  # No globals, use controls directly
    end

    function rhs!(dx, x, p, t)
        mul!(dx, sys.H(u_vec(t), t), x, -im, 0.0)
        return nothing
    end
    return rhs!
end

function _construct_rhs(sys::OpenQuantumSystem, u::F) where {F}
    Ls = sys.dissipation_operators
    Ks = map(L -> adjoint(L) * L, Ls)  # precompute L†L once
    tmp = similar(Matrix{ComplexF64}, (sys.levels, sys.levels))  # buffer

    rhs!(dρ, ρ, p, t) = begin
        Ht = sys.H(u(t), t)

        # dρ = -im*(Hρ - ρH)  (accumulate directly)
        mul!(dρ, Ht, ρ, -im, 0.0)   # dρ = -im*H*ρ
        mul!(dρ, ρ, Ht, im, 1.0)   # dρ +=  im*ρ*H

        # dρ += Σ [ LρL† - 1/2(Kρ + ρK) ]
        @inbounds for (L, K) in zip(Ls, Ks)
            mul!(tmp, L, ρ)
            mul!(dρ, tmp, adjoint(L), 1.0, 1.0)  # dρ += tmp*L†

            mul!(dρ, K, ρ, -0.5, 1.0)
            mul!(dρ, ρ, K, -0.5, 1.0)
        end

        return nothing
    end

    return rhs!
end

# ------------------------------------------
# Standard, sparse ODE integrators
# ------------------------------------------
# TODO: document solve kwarg defaults
# TODO: states must be vector (not sparse), but could infer eltype (NT eltype?)

function KetODEProblem(
    sys::AbstractQuantumSystem,
    u::F,
    ψ0::Vector{ComplexF64},
    times::AbstractVector{<:Real};
    state_name::Symbol = :ψ,
    control_name::Symbol = :u,
    kwargs...,
) where {F}
    rhs! = _construct_rhs(sys, u)
    sii_sys = PiccoloRolloutSystem(state_name => sys.levels)
    return ODEProblem(
        ODEFunction(rhs!; sys = sii_sys),
        ψ0,
        (0, times[end]);
        tstops = times,
        kwargs...,
    )
end

function UnitaryODEProblem(
    sys::AbstractQuantumSystem,
    u::F,
    times::AbstractVector{<:Real};
    U0::Matrix{ComplexF64} = Matrix{ComplexF64}(I, sys.levels, sys.levels),
    state_name::Symbol = :U,
    control_name::Symbol = :u,
    kwargs...,
) where {F}
    rhs! = _construct_rhs(sys, u)
    sii_sys = PiccoloRolloutSystem(state_name => (sys.levels, sys.levels))
    return ODEProblem(
        ODEFunction(rhs!; sys = sii_sys),
        U0,
        (0, times[end]);
        tstops = times,
        kwargs...,
    )
end

function DensityODEProblem(
    sys::OpenQuantumSystem,
    u::F,
    ρ0::Matrix{ComplexF64},
    times::AbstractVector{<:Real};
    state_name::Symbol = :ρ,
    control_name::Symbol = :u,
    kwargs...,
) where {F}
    n = sys.levels
    rhs! = _construct_rhs(sys, u)
    sii_sys = PiccoloRolloutSystem(state_name => (n, n))
    return ODEProblem(
        ODEFunction(rhs!; sys = sii_sys),
        ρ0,
        (0, times[end]);
        tstops = times,
        kwargs...,
    )
end

# ------------------------------------------
# Lie Group ODE solvers (e.g., Magnus)
# ------------------------------------------
# TODO: Operator integrator for Density

function KetOperatorODEProblem(
    sys::AbstractQuantumSystem,
    u::F,
    ψ0::Vector{ComplexF64},
    times::AbstractVector{<:Real};
    state_name::Symbol = :ψ,
    control_name::Symbol = :u,
    kwargs...,
) where {F}
    op! = _construct_operator(sys, u)
    sii_sys = PiccoloRolloutSystem(state_name => sys.levels)
    return ODEProblem(
        ODEFunction(op!; sys = sii_sys),
        ψ0,
        (0, times[end]);
        tstops = times,
        kwargs...,
    )
end

function UnitaryOperatorODEProblem(
    sys::AbstractQuantumSystem,
    u::F,
    times::AbstractVector{<:Real};
    U0::Matrix{ComplexF64} = Matrix{ComplexF64}(I, sys.levels, sys.levels),
    state_name::Symbol = :U,
    control_name::Symbol = :u,
    kwargs...,
) where {F}
    op! = _construct_operator(sys, u)
    sii_sys = PiccoloRolloutSystem(state_name => (sys.levels, sys.levels))
    return ODEProblem(
        ODEFunction(op!; sys = sii_sys),
        U0,
        (0, times[end]);
        tstops = times,
        kwargs...,
    )
end

# ------------------------------------------------------------ #
# Rollout fidelity methods
# ------------------------------------------------------------ #
# TODO: These can be extension methods for OrdinaryDiffEq
# TODO: Adapt these methods to use quantum trajectories (only _one_ rollout_fidelity method (remove unitary_rollout_fidelity), have ensemble trajectory for EnsembleProblem, etc.)

function rollout_fidelity(
    traj::NamedTrajectory,
    sys::AbstractQuantumSystem;
    state_name::Symbol = :ψ̃,
    control_name::Symbol = :u,
    algorithm = MagnusAdapt4(),
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
    interpolation::Symbol = :linear,  # :constant, :linear, or :cubic
)
    state_names = [n for n ∈ traj.names if startswith(string(n), string(state_name))]
    isempty(state_names) && error("Trajectory does not contain $(state_name).")

    # Select interpolation method for controls
    if interpolation == :constant
        u = ConstantInterpolation(traj, control_name)
    elseif interpolation == :linear
        u = LinearInterpolation(traj, control_name)
    elseif interpolation == :cubic
        u = CubicSplineInterpolation(traj, control_name)
    else
        error(
            "Unknown interpolation method: $(interpolation). Use :constant, :linear, or :cubic",
        )
    end
    times = get_times(traj)

    # Blank initial state
    tmp0 = zeros(ComplexF64, sys.levels)
    rollout = KetOperatorODEProblem(sys, u, tmp0, times, state_name = state_name)

    # Ensemble over initial states
    prob_func(prob, i, repeat) = remake(prob, u0 = iso_to_ket(traj.initial[state_names[i]]))
    ensemble_prob = EnsembleProblem(rollout, prob_func = prob_func)
    ensemble_sol = solve(
        ensemble_prob,
        algorithm;
        trajectories = length(state_names),
        saveat = [times[end]],
        abstol = abstol,
        reltol = reltol,
    )

    fids = map(zip(ensemble_sol, state_names)) do (sol, name)
        xf = sol[state_name][end]
        xg = iso_to_ket(traj.goal[name])
        fidelity(xf, xg)
    end
    return length(fids) == 1 ? fids[1] : fids
end

function unitary_rollout_fidelity(
    traj::NamedTrajectory,
    sys::AbstractQuantumSystem;
    state_name::Symbol = :Ũ⃗,
    control_name::Symbol = :u,
    algorithm = MagnusAdapt4(),
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
    interpolation::Symbol = :linear,  # :constant, :linear, or :cubic
)
    state_name ∉ traj.names && error("Trajectory does not contain $(state_name).")

    # Select interpolation method for controls
    if interpolation == :constant
        u = ConstantInterpolation(traj, control_name)
    elseif interpolation == :linear
        u = LinearInterpolation(traj, control_name)
    elseif interpolation == :cubic
        u = CubicSplineInterpolation(traj, control_name)
    else
        error(
            "Unknown interpolation method: $(interpolation). Use :constant, :linear, or :cubic",
        )
    end
    times = get_times(traj)

    x0 = iso_vec_to_operator(traj.initial[state_name])
    rollout = UnitaryOperatorODEProblem(sys, u, times, U0 = x0, state_name = state_name)
    sol = solve(rollout, algorithm; saveat = [times[end]], abstol = abstol, reltol = reltol)
    xf = sol[state_name][end]
    xg = iso_vec_to_operator(traj.goal[state_name])
    return unitary_fidelity(xf, xg)
end

function unitary_rollout(
    traj::NamedTrajectory,
    sys::AbstractQuantumSystem;
    state_name::Symbol = :Ũ⃗,
    control_name::Symbol = :u,
    algorithm = MagnusAdapt4(),
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
    interpolation::Symbol = :linear,  # :constant, :linear, or :cubic
)
    state_name ∉ traj.names && error("Trajectory does not contain $(state_name).")

    # Select interpolation method for controls
    if interpolation == :constant
        u = ConstantInterpolation(traj, control_name)
    elseif interpolation == :linear
        u = LinearInterpolation(traj, control_name)
    elseif interpolation == :cubic
        u = CubicSplineInterpolation(traj, control_name)
    else
        error(
            "Unknown interpolation method: $(interpolation). Use :constant, :linear, or :cubic",
        )
    end
    times = get_times(traj)

    x0 = iso_vec_to_operator(traj.initial[state_name])
    prob = UnitaryOperatorODEProblem(sys, u, times, U0 = x0, state_name = state_name)
    sol = solve(prob, algorithm; saveat = times, abstol = abstol, reltol = reltol)

    # Extract and convert to iso-vec trajectory
    Ũ⃗_traj = hcat([operator_to_iso_vec(sol[state_name][i]) for i = 1:length(times)]...)

    return Ũ⃗_traj
end

function ket_rollout_fidelity(
    traj::NamedTrajectory,
    sys::AbstractQuantumSystem;
    state_name::Symbol = :ψ̃,
    control_name::Symbol = :u,
    algorithm = MagnusAdapt4(),
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
    interpolation::Symbol = :linear,  # :constant, :linear, or :cubic
)
    return rollout_fidelity(
        traj,
        sys;
        state_name = state_name,
        control_name = control_name,
        algorithm = algorithm,
        abstol = abstol,
        reltol = reltol,
        interpolation = interpolation,
    )
end

function ket_rollout(
    traj::NamedTrajectory,
    sys::AbstractQuantumSystem;
    state_name::Symbol = :ψ̃,
    control_name::Symbol = :u,
    algorithm = MagnusAdapt4(),
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
    interpolation::Symbol = :linear,  # :constant, :linear, or :cubic
)
    state_name ∉ traj.names && error("Trajectory does not contain $(state_name).")

    # Select interpolation method for controls
    if interpolation == :constant
        u = ConstantInterpolation(traj, control_name)
    elseif interpolation == :linear
        u = LinearInterpolation(traj, control_name)
    elseif interpolation == :cubic
        u = CubicSplineInterpolation(traj, control_name)
    else
        error(
            "Unknown interpolation method: $(interpolation). Use :constant, :linear, or :cubic",
        )
    end
    times = get_times(traj)

    ψ0 = iso_to_ket(traj.initial[state_name])
    prob = KetOperatorODEProblem(sys, u, ψ0, times, state_name = state_name)
    sol = solve(prob, algorithm; saveat = times, abstol = abstol, reltol = reltol)

    # Extract and convert to iso-vec trajectory
    ψ̃_traj = hcat([ket_to_iso(sol[state_name][i]) for i = 1:length(times)]...)

    return ψ̃_traj
end

# ------------------------------------------------------------ #
# Minimal interface 
# https://docs.sciml.ai/SymbolicIndexingInterface/
# ------------------------------------------------------------ #

_name(sym::Symbol) = sym
_name(::Any) = nothing

SII.constant_structure(::PiccoloRolloutSystem) = true
SII.default_values(sys::PiccoloRolloutSystem) = sys.defaults

SII.is_time_dependent(sys::PiccoloRolloutSystem) = true
SII.is_independent_variable(sys::PiccoloRolloutSystem, sym) = _name(sym) === sys.t
SII.independent_variable_symbols(sys::PiccoloRolloutSystem) = [sys.t]

# solved variables (state)
SII.is_variable(sys::PiccoloRolloutSystem, sym) = haskey(sys.state_index, _name(sym))
SII.variable_index(sys::PiccoloRolloutSystem, sym) =
    get(sys.state_index, _name(sym), nothing)
SII.variable_symbols(sys::PiccoloRolloutSystem) = collect(keys(sys.state_index))

# parameters (none)
SII.is_parameter(::PiccoloRolloutSystem, _) = false
SII.parameter_index(::PiccoloRolloutSystem, _) = nothing
SII.parameter_symbols(::PiccoloRolloutSystem) = Symbol[]

SII.is_observed(sys::PiccoloRolloutSystem, sym) = false

# *************************************************************************** #
# TODO: Test rollout fidelity (after adpating to new interface)

@testitem "Test ket rollout symbolic interface" begin
    using OrdinaryDiffEqTsit5

    T, Δt = 1.0, 0.1
    sys = QuantumSystem([PAULIS.X, PAULIS.Y], [1.0, 1.0])
    ψ0 = ComplexF64[1, 0]
    u = t -> [t; 0.0]
    times = 0:Δt:T
    rollout = KetODEProblem(sys, u, ψ0, times)

    # test default
    sol1 = solve(rollout, Tsit5())
    @test sol1[:ψ] ≈ sol1.u

    # test solve kwargs
    sol2 = solve(rollout, Tsit5(), saveat = [times[end]])
    @test length(sol2[:ψ]) == 1
    @test length(sol2[:ψ][1]) == length(ψ0)

    # rename 
    rollout = KetODEProblem(sys, u, ψ0, times, state_name = :x)
    sol = solve(rollout, Tsit5())
    @test sol[:x] ≈ sol.u
end

@testitem "Test unitary rollout symbolic interface" begin
    using OrdinaryDiffEqLinear

    T, Δt = 1.0, 0.1
    sys = QuantumSystem([PAULIS.X, PAULIS.Y], [1.0, 1.0])
    u = t -> [t; 0.0]
    times = 0:Δt:T
    rollout = UnitaryOperatorODEProblem(sys, u, times)

    # test default
    sol1 = solve(rollout, MagnusGL4())
    @test sol1[:U] ≈ sol1.u

    # test solve kwargs
    sol2 = solve(rollout, MagnusGL4(), saveat = [times[end]])
    @test length(sol2[:U]) == 1
    @test size(sol2[:U][1]) == (sys.levels, sys.levels)

    # rename 
    rollout = UnitaryOperatorODEProblem(sys, u, times, state_name = :X)
    sol = solve(rollout, MagnusGL4())
    @test sol[:X] ≈ sol.u
end

@testitem "Test density rollout symbolic interface" begin
    using OrdinaryDiffEqTsit5

    T, Δt = 1.0, 0.1
    csys = QuantumSystem([PAULIS.X, PAULIS.Y], [1.0, 1.0])
    a = ComplexF64[0 1; 0 0]
    sys = OpenQuantumSystem(csys, dissipation_operators = [1e-3 * a])
    u = t -> [t; 0.0]
    times = 0:Δt:T

    ψ0 = ComplexF64[1, 0]
    ρ0 = ψ0 * ψ0'
    rollout = DensityODEProblem(sys, u, ρ0, times)

    # test default symbolic access
    sol1 = solve(rollout, Tsit5())
    @test sol1[:ρ] ≈ sol1.u

    # test solve kwargs
    sol2 = solve(rollout, Tsit5(), saveat = [times[end]])
    @test length(sol2[:ρ]) == 1
    @test size(sol2[:ρ][1]) == (sys.levels, sys.levels)

    # rename
    rollout = DensityODEProblem(sys, u, ρ0, times, state_name = :X)
    sol = solve(rollout, Tsit5())
    @test sol[:X] ≈ sol.u
end

@testitem "Rollout internal consistency (ket/unitary/density, closed system)" begin
    using OrdinaryDiffEqTsit5
    using OrdinaryDiffEqLinear

    T, Δt = 1.0, 0.1
    sys = QuantumSystem([PAULIS.X, PAULIS.Y], [1.0, 1.0])
    osys = OpenQuantumSystem(sys)

    u = t -> [t; 0.0]
    times = 0:Δt:T
    ψ0 = ComplexF64[1, 0]
    ρ0 = ψ0 * ψ0'

    ket_prob = KetODEProblem(sys, u, ψ0, times)
    U_prob = UnitaryOperatorODEProblem(sys, u, times)
    rho_prob = DensityODEProblem(osys, u, ρ0, times)

    # Save only final state so comparisons are well-defined
    kw = (dense = false, save_everystep = false, save_start = false, save_end = true)
    ket_sol = solve(ket_prob, Tsit5(); kw...)
    U_sol = solve(U_prob, MagnusGL4(); kw...)
    ρ_sol = solve(rho_prob, Tsit5(); kw...)

    ψT = ket_sol.u[end]
    UT = U_sol.u[end]
    ρT = ρ_sol.u[end]

    @test ψT ≈ UT * ψ0
    @test ρT ≈ ψT * ψT' atol = 1e-5
    @test ρT ≈ UT * ρ0 * UT' atol = 1e-5
end

@testitem "Rollouts with all Pulse types" begin
    using OrdinaryDiffEqTsit5
    using OrdinaryDiffEqLinear

    T, Δt = 1.0, 0.01
    sys = QuantumSystem([PAULIS.X, PAULIS.Y], [1.0, 1.0])
    osys = OpenQuantumSystem(sys)
    times = 0:Δt:T
    n_times = length(times)
    ψ0 = ComplexF64[1, 0]
    ρ0 = ψ0 * ψ0'

    # Generate test control values (smooth ramp)
    controls = [sin(π * t / T) for t in times]
    control_matrix = [controls zeros(n_times)]'  # 2 drives × n_times timesteps

    # Test all pulse types
    pulse_types = [
        ZeroOrderPulse(control_matrix, times),
        LinearSplinePulse(control_matrix, times),
        CubicSplinePulse(control_matrix, times),
    ]

    for pulse in pulse_types
        # Verify pulse is callable and returns correct shape
        @test length(pulse(0.0)) == 2
        @test pulse(0.0) ≈ [0.0, 0.0]

        # KetODEProblem
        ket_prob = KetODEProblem(sys, pulse, ψ0, times)
        ket_sol = solve(ket_prob, Tsit5(); saveat=times)
        @test length(ket_sol.u) == n_times
        @test length(ket_sol.u[end]) == 2  # 2-level system

        # UnitaryOperatorODEProblem (for MagnusGL4)
        U_prob = UnitaryOperatorODEProblem(sys, pulse, times)
        U_sol = solve(U_prob, MagnusGL4())
        @test length(U_sol.u) == n_times
        @test size(U_sol.u[end]) == (2, 2)

        # DensityODEProblem
        rho_prob = DensityODEProblem(osys, pulse, ρ0, times)
        rho_sol = solve(rho_prob, Tsit5(); saveat=times)
        @test length(rho_sol.u) == n_times
        @test size(rho_sol.u[end]) == (2, 2)

        # Check consistency: ψ_final should equal U_final * ψ0.
        #
        # Why Δt = 0.01 is required for ZeroOrderPulse:
        #   ZeroOrderPulse is piecewise constant — it jumps at each sample time
        #   t_k. The problem constructor sets `tstops = times`, which forces both
        #   solvers to evaluate at every t_k exactly.
        #
        #   Tsit5 is an FSAL (First Same As Last) method: its final RK stage in
        #   the step [t_k, t_{k+1}] is evaluated at t_{k+1} and reused as the
        #   first stage of the next step. Because t_{k+1} is a tstop (a pulse
        #   jump point), `pulse(t_{k+1})` returns the *next* interval's control
        #   value. This bleeds a small amount of the next Hamiltonian into the
        #   current step update, causing an O(Δt * Δu_k) error per step, where
        #   Δu_k is the control jump. Summing over all steps, the total error is
        #   O(Δt * TV(u)) where TV(u) is the total variation of the pulse — i.e.
        #   proportional to Δt and independent of the number of steps.
        #
        #   MagnusGL4, by contrast, evaluates its Gauss-Legendre quadrature
        #   nodes strictly inside (t_k, t_{k+1}), never at the boundary, so it
        #   always sees the correct constant value for the current interval.
        #
        #   With Δt = 0.1 the accumulated error is ~0.011, exceeding atol = 1e-2.
        #   With Δt = 0.01 the error shrinks to ~0.001, passing comfortably.
        #   This Δt-scaling confirms the O(Δt) theory above.
        ψT = ket_sol.u[end]
        UT = U_sol.u[end]
        @test ψT ≈ UT * ψ0 atol = 1e-2
    end
end

@testitem "Rollouts with GaussianPulse" begin
    using OrdinaryDiffEqTsit5
    using OrdinaryDiffEqLinear

    T = 1.0
    Δt = 0.1
    sys = QuantumSystem([PAULIS.X, PAULIS.Y], [1.0, 1.0])
    times = 0:Δt:T
    n_times = length(times)
    ψ0 = ComplexF64[1, 0]

    # Create GaussianPulse with 2 drives
    # Constructor: GaussianPulse(amplitudes, sigmas, centers, duration)
    amplitudes = [1.0, 0.5]
    sigmas = [T / 4, T / 4]
    centers = [T / 2, T / 2]
    pulse = GaussianPulse(amplitudes, sigmas, centers, T)

    # Verify pulse properties
    @test duration(pulse) == T
    @test n_drives(pulse) == 2
    @test length(pulse(T / 2)) == 2

    # Peak should be at t = center (T/2)
    @test pulse(T / 2)[1] ≈ 1.0 atol = 1e-10
    @test pulse(T / 2)[2] ≈ 0.5 atol = 1e-10

    # Should be symmetric around center
    @test pulse(0.25)[1] ≈ pulse(0.75)[1] atol = 1e-10

    # KetODEProblem
    ket_prob = KetODEProblem(sys, pulse, ψ0, times)
    ket_sol = solve(ket_prob, Tsit5(); saveat=times)
    @test length(ket_sol.u) == n_times

    # UnitaryOperatorODEProblem
    U_prob = UnitaryOperatorODEProblem(sys, pulse, times)
    U_sol = solve(U_prob, MagnusGL4())
    @test length(U_sol.u) == n_times

    # Check consistency
    # Note: different solvers (Tsit5 vs MagnusGL4) have different accuracy
    ψT = ket_sol.u[end]
    UT = U_sol.u[end]
    @test ψT ≈ UT * ψ0 atol = 1e-2
end

@testitem "Two ways to check fidelity" begin
    using OrdinaryDiffEqLinear
    using NamedTrajectories

    # Setup
    T = 1.0
    sys = QuantumSystem([PAULIS.X, PAULIS.Y], [1.0, 1.0])
    X_gate = ComplexF64[0 1; 1 0]

    # Method 1: Fast fidelity from quantum trajectory (O(1))
    pulse = ZeroOrderPulse([0.5 0.5; 0.1 0.1], [0.0, T])
    qtraj = UnitaryTrajectory(sys, pulse, X_gate)
    fid1 = fidelity(qtraj)  # Uses stored solution - FAST!
    @test fid1 isa Float64
    @test 0.0 <= fid1 <= 1.0

    # Method 2: Validate discrete controls (for NamedTrajectory)
    # This is useful when you have discrete trajectory and want to test interpolation
    I_matrix = ComplexF64[1 0; 0 1]
    traj = NamedTrajectory(
        (Ũ⃗ = randn(8, 11), u = randn(2, 11), Δt = fill(T / 10, 11));
        controls = :u,
        timestep = :Δt,
        initial = (Ũ⃗ = operator_to_iso_vec(I_matrix),),
        goal = (Ũ⃗ = operator_to_iso_vec(X_gate),),
    )

    # Test different interpolation methods (use unitary_rollout_fidelity for unitaries)
    fid_constant =
        unitary_rollout_fidelity(traj, sys; state_name = :Ũ⃗, interpolation = :constant)
    fid_linear =
        unitary_rollout_fidelity(traj, sys; state_name = :Ũ⃗, interpolation = :linear)

    @test fid_constant isa Float64
    @test fid_linear isa Float64
    @test 0.0 <= fid_constant <= 1.0
    @test 0.0 <= fid_linear <= 1.0
end

@testitem "rollout with new pulse" begin
    using OrdinaryDiffEqLinear: MagnusGL4

    # Setup
    T = 1.0
    sys = QuantumSystem([PAULIS.X, PAULIS.Y], [1.0, 1.0])
    X_gate = ComplexF64[0 1; 1 0]

    # Create initial trajectory
    pulse1 = ZeroOrderPulse([0.5 0.5; 0.1 0.1], [0.0, T])
    qtraj1 = UnitaryTrajectory(sys, pulse1, X_gate)
    fid1 = fidelity(qtraj1)

    # Roll out a new pulse
    pulse2 = ZeroOrderPulse([0.8 0.8; 0.2 0.2], [0.0, T])
    qtraj2 = rollout(qtraj1, pulse2)
    fid2 = fidelity(qtraj2)

    # Should have different fidelities (different pulses)
    @test fid2 != fid1
    @test qtraj2.pulse === pulse2
    @test qtraj2.system === qtraj1.system

    # Roll out with custom resolution
    qtraj3 = rollout(qtraj1, pulse2; n_save = 501)
    @test length(qtraj3.solution.u) == 501
end

@testitem "Global parameter updates" begin
    using LinearAlgebra
    using NamedTrajectories

    # Create a system with global parameters
    H_drives = [PAULIS[:X], PAULIS[:Y]]
    global_params = (δ = 0.5, Ω = 1.0)
    sys = QuantumSystem(H_drives, [1.0, 1.0]; global_params = global_params)

    # Create a unitary trajectory (2 drives × 2 timesteps)
    pulse = ZeroOrderPulse([0.5 0.3; 0.5 0.3], [0.0, 1.0])
    U_goal = PAULIS[:X]
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)

    # Verify initial global parameters
    @test qtraj.system.global_params.δ == 0.5
    @test qtraj.system.global_params.Ω == 1.0

    # Create a NamedTrajectory with different global values
    traj = NamedTrajectory(
        (u = rand(2, 10), Δt = fill(0.1, 10));
        timestep = :Δt,
        global_data = [0.8, 1.5],
        global_components = (δ = 1:1, Ω = 2:2),
    )

    # Update global parameters
    Rollouts.update_global_params!(qtraj, traj)

    # Verify updated values
    @test qtraj.system.global_params.δ == 0.8
    @test qtraj.system.global_params.Ω == 1.5

    # Verify system structure preserved
    @test qtraj.system.n_drives == 2
    @test qtraj.system.levels == 2
    @test length(qtraj.system.H_drives) == 2

    # Test with KetTrajectory
    ψ_init = ComplexF64[1.0, 0.0]
    ψ_goal = ComplexF64[0.0, 1.0]
    qtraj_ket = KetTrajectory(sys, pulse, ψ_init, ψ_goal)

    Rollouts.update_global_params!(qtraj_ket, traj)
    @test qtraj_ket.system.global_params.δ == 0.8
    @test qtraj_ket.system.global_params.Ω == 1.5
end

@testitem "extract_globals utility" begin
    using NamedTrajectories

    # Create trajectory with globals
    traj = NamedTrajectory(
        (u = rand(2, 10), Δt = fill(0.1, 10));
        timestep = :Δt,
        global_data = [0.8, 1.5, 2.0],
        global_components = (δ = 1:1, Ω = 2:2, α = 3:3),
    )

    # Extract all globals
    g_all = Rollouts.extract_globals(traj)
    @test g_all isa NamedTuple
    @test g_all.δ == 0.8
    @test g_all.Ω == 1.5
    @test g_all.α == 2.0

    # Extract specific globals
    g_partial = Rollouts.extract_globals(traj, [:δ, :Ω])
    @test g_partial isa NamedTuple
    @test g_partial.δ == 0.8
    @test g_partial.Ω == 1.5
    @test !haskey(g_partial, :α)

    # Test with trajectory without global components (edge case)
    traj_no_globals = NamedTrajectory((u = rand(2, 10), Δt = fill(0.1, 10)); timestep = :Δt)
    g_empty = Rollouts.extract_globals(traj_no_globals)
    @test g_empty isa NamedTuple
    @test isempty(g_empty)
end

@testitem "Multi-dimensional global parameters" begin
    using NamedTrajectories

    # Test extract_globals with multi-dimensional globals
    traj = NamedTrajectory(
        (u = rand(2, 10), Δt = fill(0.1, 10));
        timestep = :Δt,
        global_data = [0.8, 1.5, 2.0, 3.0],  # Two scalars and one 2D vector
        global_components = (δ = 1:1, Ω = 2:2, α = 3:4),
    )

    g = Rollouts.extract_globals(traj)
    @test g.δ == 0.8
    @test g.Ω == 1.5
    @test g.α isa Vector
    @test g.α == [2.0, 3.0]
end

@testitem "update_global_params! edge cases" begin
    using NamedTrajectories

    # Create a system with global parameters
    H_drives = [PAULIS[:X], PAULIS[:Y]]
    global_params = (δ = 0.5, Ω = 1.0)
    sys = QuantumSystem(H_drives, [1.0, 1.0]; global_params = global_params)
    pulse = ZeroOrderPulse([0.5 0.3; 0.5 0.3], [0.0, 1.0])
    U_goal = PAULIS[:X]
    qtraj = UnitaryTrajectory(sys, pulse, U_goal)

    # Test with trajectory without global components (should not error)
    traj_no_globals = NamedTrajectory((u = rand(2, 10), Δt = fill(0.1, 10)); timestep = :Δt)

    # Should return nothing without error
    result = Rollouts.update_global_params!(qtraj, traj_no_globals)
    @test result === nothing
    # Original global params should be unchanged
    @test qtraj.system.global_params.δ == 0.5
    @test qtraj.system.global_params.Ω == 1.0
end


end
