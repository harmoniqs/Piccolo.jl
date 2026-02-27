module QuantumSystemUtils

export is_reachable

using ..EmbeddedOperators
using ..QuantumObjectUtils
using ..QuantumSystems
using ..Gates

using LinearAlgebra
using SparseArrays
using TestItems


commutator(A::AbstractMatrix{<:Number}, B::AbstractMatrix{<:Number}) = A * B - B * A

is_hermitian(H::AbstractMatrix{<:Number}; atol = eps(Float32)) =
    all(isapprox.(H - H', 0.0, atol = atol))

# ----------------------------------------------------------------------------- #
# Linear independence
# ----------------------------------------------------------------------------- #

"""
    is_linearly_dependent(M::AbstractMatrix; eps=eps(Float32), verbose=true)

Check if the columns of the matrix `M` are linearly dependent.
"""
function is_linearly_dependent(M::AbstractMatrix; eps = eps(Float32), verbose = true)
    if size(M, 2) > size(M, 1)
        if verbose
            println(
                "Linearly dependent because columns > rows, $(size(M, 2)) > $(size(M, 1)).",
            )
        end
        return true
    end
    # QR decomposition has a zero R on diagonal if linearly dependent
    val = minimum(abs.(diag(qr(M).R)))
    return isapprox(val, 0.0, atol = eps)
end

is_linearly_dependent(basis::Vector{<:AbstractMatrix{<:Number}}; kwargs...) =
    is_linearly_dependent(stack(vec.(basis)); kwargs...)

function linearly_independent_indices(
    basis::Vector{<:AbstractMatrix{<:Number}};
    order = 1:length(basis),
    kwargs...,
)
    @assert issetequal(order, 1:length(basis)) "Order must enumerate entire basis."
    bᵢ = Int[]
    for i ∈ order
        if !is_linearly_dependent([basis[bᵢ]..., basis[i]]; kwargs...)
            push!(bᵢ, i)
        end
    end
    return bᵢ
end

function linearly_independent_subset(basis::Vector{<:AbstractMatrix}; kwargs...)
    bᵢ = linearly_independent_indices(basis; kwargs...)
    return deepcopy(basis[bᵢ])
end

function linearly_independent_subset!(basis::Vector{<:AbstractMatrix}; kwargs...)
    bᵢ = linearly_independent_indices(basis; kwargs...)
    deleteat!(basis, setdiff(1:length(basis), bᵢ))
    return nothing
end

# ----------------------------------------------------------------------------- #
# Operator algebra
# ----------------------------------------------------------------------------- #

traceless(M::AbstractMatrix) = M - tr(M) * I / size(M, 1)

function clean(M::AbstractMatrix; normalize::Bool = true, remove_trace::Bool = true)
    M_ = remove_trace ? traceless(M) : M
    return normalize && norm(M_) > 0 ? M_ / norm(M_) : M_
end

"""
    operator_algebra(generators; kwargs...)

Compute the Lie algebra basis for the given `generators`.

# Arguments
- `generators::Vector{<:AbstractMatrix}`: generators of the Lie algebra

# Keyword Arguments
- `return_layers::Bool=false`: return the Lie tree layers
- `normalize::Bool=false`: normalize the basis
- `verbose::Bool=false`: print information
- `remove_trace::Bool=true`: remove trace from generators
"""
function operator_algebra(
    generators::Vector{<:AbstractMatrix{T}};
    return_layers::Bool = false,
    remove_trace::Bool = true,
    normalize::Bool = true,
    verbose::Bool = false,
    atol::Float64 = 1e-10,
) where {T<:Number}
    N = size(first(generators), 1)
    max_dim = N^2 - 1

    # Initialize orthonormal basis of vectorized operators (Gram-Schmidt)
    q_basis = Matrix{T}(undef, N^2, 0)
    basis = Matrix{T}[]
    current_layer = Matrix{T}[]

    # Prepare initial generators
    for g in generators
        M = clean(g, normalize = normalize, remove_trace = remove_trace)
        v = vec(M)
        if size(q_basis, 2) == 0
            q_basis = reshape(v / norm(v), :, 1)
            push!(basis, M)
            push!(current_layer, M)
        else
            proj = q_basis * (q_basis' * v)
            residual = v - proj
            rnorm = norm(residual)
            if rnorm > atol
                q_basis = hcat(q_basis, residual / rnorm)
                push!(basis, M)
                push!(current_layer, M)
            end
        end
    end

    if verbose
        print("operator algebra depth = [1")
    end

    all_layers = return_layers ? [copy(current_layer)] : nothing
    layer_count = 1

    while !isempty(current_layer) && length(basis) < max_dim
        next_layer = Matrix{T}[]

        for A in current_layer, B in basis
            comm = commutator(A, B)

            if all(comm .≈ 0)
                continue
            end

            comm = is_hermitian(comm) ? comm : im * comm
            comm = clean(comm, normalize = normalize, remove_trace = remove_trace)

            v = vec(comm)
            proj = q_basis * (q_basis' * v)
            residual = v - proj
            rnorm = norm(residual)

            if rnorm > atol
                q_basis = hcat(q_basis, residual / rnorm)
                push!(basis, comm)
                push!(next_layer, comm)
            end
        end

        if isempty(next_layer)
            break
        end

        current_layer = next_layer
        layer_count += 1
        if verbose
            print(" $layer_count")
        end
        if return_layers
            push!(all_layers, copy(current_layer))
        end
    end

    if verbose
        println("]")
    end

    return return_layers ? (basis, all_layers) : basis
end

function fit_gen_to_basis(
    gen::AbstractMatrix{<:Number},
    basis::AbstractVector{<:AbstractMatrix{<:Number}},
)
    A = stack(vec.(basis))
    b = vec(gen)
    return A \ b
end

function is_in_span(
    gen::AbstractMatrix{<:Number},
    basis::AbstractVector{<:AbstractMatrix{<:Number}};
    subspace::AbstractVector{Int} = 1:size(gen, 1),
    atol = eps(Float32),
    return_effective_gen = false,
)
    g_basis = [deepcopy(b[subspace, subspace]) for b ∈ basis]
    linearly_independent_subset!(g_basis)
    # Traceless basis needs traceless fit
    x = fit_gen_to_basis(gen, g_basis)
    g_eff = sum(x .* g_basis)
    ε = norm(g_eff - gen, 2)
    if return_effective_gen
        return ε < atol, g_eff
    else
        return ε < atol
    end
end

# ----------------------------------------------------------------------------- #
# Reachability
# ----------------------------------------------------------------------------- #

"""
    is_reachable(gate, hamiltonians; kwargs...)

Check if the `gate` is reachable using the given `hamiltonians`.

# Arguments
- `gate::AbstractMatrix`: target gate
- `hamiltonians::AbstractVector{<:AbstractMatrix}`: generators of the Lie algebra

# Keyword Arguments
- `subspace::AbstractVector{<:Int}=1:size(gate, 1)`: subspace indices
- `compute_basis::Bool=true`: compute the basis or use the Hamiltonians directly
- `remove_trace::Bool=true`: remove trace from generators
- `verbose::Bool=true`: print information about the operator algebra
- `atol::Float32=eps(Float32)`: absolute tolerance

See also [`QuantumSystemUtils.operator_algebra`](@ref).
"""
function is_reachable(
    gate::AbstractMatrix{<:Number},
    hamiltonians::AbstractVector{<:AbstractMatrix{<:Number}};
    subspace::AbstractVector{Int} = 1:size(gate, 1),
    compute_basis = true,
    remove_trace = true,
    verbose = true,
    atol = eps(Float32),
)
    @assert size(gate, 1) == length(subspace) "Gate must be given in the subspace."
    generator = im * log(gate)

    if remove_trace
        generator = traceless(generator)
    end

    if compute_basis
        basis =
            operator_algebra(hamiltonians, remove_trace = remove_trace, verbose = verbose)
    else
        basis = hamiltonians
    end

    return is_in_span(generator, basis, subspace = subspace, atol = atol)
end

"""
    is_reachable(gate::AbstractMatrix{<:Number}, system::AbstractQuantumSystem; kwargs...)

Check if the `gate` is reachable using the given `system`.

# Keyword Arguments
- `use_drift::Bool=true`: include drift Hamiltonian in the generators
- `kwargs...`: keyword arguments for `is_reachable`
"""
function is_reachable(
    gate::AbstractMatrix{<:Number},
    system::AbstractQuantumSystem;
    use_drift::Bool = true,
    kwargs...,
)
    H_drift = get_drift(system)
    H_drives = get_drives(system)
    if use_drift && !all(H_drift .≈ 0)
        push!(H_drives, H_drift)
    end
    return is_reachable(gate, H_drives; kwargs...)
end

is_reachable(gate::EmbeddedOperator, args...; kwargs...) =
    is_reachable(unembed(gate), args...; subspace = gate.subspace, kwargs...)

# ****************************************************************************** #

@testitem "Lie algebra basis" begin
    using .QuantumSystemUtils: operator_algebra

    # Check 1 qubit with complete basis
    gen = operator_from_string.(["X", "Y"])
    basis = operator_algebra(gen, return_layers = false, verbose = false)
    @test length(basis) == size(first(gen), 1)^2 - 1

    # Check 1 qubit with complete basis and layers
    basis, layers = operator_algebra(gen, return_layers = true, verbose = false)
    @test length(basis) == size(first(gen), 1)^2 - 1

    # Check 1 qubit with subspace
    gen = operator_from_string.(["X"])
    basis = operator_algebra(gen, verbose = false)
    @test length(basis) == 1

    # Check 2 qubit with complete basis
    gen = operator_from_string.(["XX", "YY", "XI", "YI", "IY", "IX"])
    basis = operator_algebra(gen, verbose = false)
    @test length(basis) == size(first(gen), 1)^2 - 1

    # Check 2 qubit linearly dependent
    res = ["XX", "XI"]
    gen = operator_from_string.(["XX", "XX", "XI", "XI"])
    basis = operator_algebra(gen, verbose = false)
    @test length(basis) == length(res)

    # Check 2 qubit with linearly dependent basis
    res = ["XX", "YY", "XI", "IX", "ZY", "YZ", "ZZ"]
    gen = operator_from_string.(["XX", "YY", "XI", "XI", "IX"])
    basis = operator_algebra(gen, verbose = false)
    length(basis) == length(res)

    # Check 2 qubit with pair of 1-qubit subspaces
    gen = operator_from_string.(["XI", "YI", "IY", "IX"])
    basis = operator_algebra(gen, verbose = false)
    @test length(basis) == 2 * (2^2 - 1)
end

@testitem "Lie Algebra reachability single qubit" begin
    # Check 1 qubit with complete basis
    gen = [PAULIS[:X], PAULIS[:Y]]
    target = PAULIS[:Z]
    @test is_reachable(target, gen, compute_basis = true, verbose = false)

    # System
    sys = QuantumSystem(
        [PAULIS[:X], PAULIS[:Y], PAULIS[:Z]],
        [(-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0)],
    )
    target = PAULIS[:Z]
    @test is_reachable(target, sys, verbose = false)

    # System with drift
    sys = QuantumSystem(PAULIS[:Z], [PAULIS[:X]], [(-1.0, 1.0)])
    target = PAULIS[:Z]
    @test is_reachable(target, sys, verbose = false)
end

@testitem "Lie Algebra reachability two qubits" begin
    using LinearAlgebra
    ⊗ = kron

    # Check 2 qubit with complete basis
    XI = PAULIS[:X] ⊗ PAULIS[:I]
    IX = PAULIS[:I] ⊗ PAULIS[:X]
    YI = PAULIS[:Y] ⊗ PAULIS[:I]
    IY = PAULIS[:I] ⊗ PAULIS[:Y]
    XX = PAULIS[:X] ⊗ PAULIS[:X]
    YY = PAULIS[:Y] ⊗ PAULIS[:Y]
    ZI = PAULIS[:Z] ⊗ PAULIS[:I]
    IZ = PAULIS[:I] ⊗ PAULIS[:Z]
    ZZ = PAULIS[:Z] ⊗ PAULIS[:Z]

    complete_gen = [XX + YY, XI, YI, IX, IY]
    incomplete_gen = [XI, ZZ]
    r = [0, 1, 2, 3, 4]
    r /= norm(r)
    R2 = exp(-im * sum([θ * H for (H, θ) in zip(complete_gen, r)]))
    CZ = GATES[:CZ]
    CX = GATES[:CX]

    # Pass
    @test is_reachable(R2, complete_gen, verbose = false)
    @test is_reachable(CZ, complete_gen, verbose = false)
    @test is_reachable(CX, complete_gen, verbose = false)
    @test is_reachable(XI, complete_gen, verbose = false)

    # Mostly fail
    @test !is_reachable(R2, incomplete_gen, verbose = false)
    @test !is_reachable(CZ, incomplete_gen, verbose = false)
    @test !is_reachable(CX, incomplete_gen, verbose = false)
    @test is_reachable(XI, incomplete_gen, verbose = false)

    # QuantumSystems
    complete_gen_sys = QuantumSystem(complete_gen, fill((-1.0, 1.0), length(complete_gen)))
    incomplete_gen_sys =
        QuantumSystem(incomplete_gen, fill((-1.0, 1.0), length(incomplete_gen)))
    # Pass
    @test is_reachable(R2, complete_gen_sys, verbose = false)
    @test is_reachable(CZ, complete_gen_sys, verbose = false)
    @test is_reachable(CX, complete_gen_sys, verbose = true)
    @test is_reachable(XI, complete_gen_sys, verbose = false)

    # Mostly fail
    @test !is_reachable(R2, incomplete_gen_sys, verbose = false)
    @test !is_reachable(CZ, incomplete_gen_sys, verbose = false)
    @test !is_reachable(CX, incomplete_gen_sys, verbose = false)
    @test is_reachable(XI, incomplete_gen_sys, verbose = false)
end

@testitem "Lie Algebra embedded subspace reachability" begin
    # Check 1 qubit with complete basis
    gen = [PAULIS[:X], PAULIS[:Y]]
    target = EmbeddedOperator(PAULIS[:Z], 1:2, 4)
    @test is_reachable(target, gen, verbose = false)

    # System
    sys = QuantumSystem([GATES[:X], GATES[:Y]], [(-1.0, 1.0), (-1.0, 1.0)])
    @test is_reachable(target, sys, verbose = false)
end

end
