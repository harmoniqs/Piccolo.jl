"""
    KnotPointPropagationData

Cached forward-rollout data for efficient VJP / JVP computation.

Stores propagators, eigendecompositions, intermediate states, and system
information produced by [`setup_knot_point_propagation`](@ref).
"""
struct KnotPointPropagationData
    propagators::Vector{Matrix{ComplexF64}}
    eigenvalues::Vector{Vector{Float64}}
    eigenvectors::Vector{Matrix{ComplexF64}}
    states::Vector{Vector{ComplexF64}}
    H_drift::Matrix{ComplexF64}
    H_drives::Vector{Matrix{ComplexF64}}
    controls::Matrix{Float64}
    Δts::Vector{Float64}
end

"""
    setup_knot_point_propagation(sys, controls, Δts, ψ0)

Perform a forward rollout of the knot-point dynamics, computing and caching
all propagators via eigendecomposition at each interval.

# Arguments
- `sys`: a `QuantumSystem` (provides `H_drift`, `H_drives`, `n_drives`, `levels`)
- `controls::AbstractMatrix`: size `[n_drives × N]`, control amplitudes at each knot point
- `Δts::AbstractVector`: time-step durations (length ≥ N-1)
- `ψ0::AbstractVector`: initial ket state

Returns a [`KnotPointPropagationData`](@ref).
"""
function setup_knot_point_propagation(
    sys,
    controls::AbstractMatrix,
    Δts::AbstractVector,
    ψ0::AbstractVector,
)
    n_drives = sys.n_drives
    N = size(controls, 2)
    n_intervals = N - 1

    H_drift = Matrix{ComplexF64}(sys.H_drift)
    H_drives = [Matrix{ComplexF64}(d.H) for d in sys.H_drives]

    propagators = Vector{Matrix{ComplexF64}}(undef, n_intervals)
    eigenvalues = Vector{Vector{Float64}}(undef, n_intervals)
    eigenvectors = Vector{Matrix{ComplexF64}}(undef, n_intervals)
    states = Vector{Vector{ComplexF64}}(undef, N)

    states[1] = ComplexF64.(ψ0)

    for j = 1:n_intervals
        # Total Hamiltonian at knot point j
        H_total = copy(H_drift)
        for k = 1:n_drives
            H_total .+= controls[k, j] .* H_drives[k]
        end

        # Eigendecomposition
        λ, V = eigen(Hermitian(H_total))
        eigenvalues[j] = λ
        eigenvectors[j] = V

        # Propagator: U_j = exp(-iΔt_j H_total)
        Δt = Δts[j]
        propagators[j] = V * Diagonal(cis.(-Δt .* λ)) * V'

        # Forward-propagated state
        states[j+1] = propagators[j] * states[j]
    end

    return KnotPointPropagationData(
        propagators,
        eigenvalues,
        eigenvectors,
        states,
        H_drift,
        H_drives,
        Matrix{Float64}(controls),
        Vector{Float64}(Δts),
    )
end

"""
    ket_vjp(data, λ; use_exact_kick=true, n_quad=4)

Backward-pass vector-Jacobian product (VJP).

Propagates the covector `λ` backward through the knot-point chain, extracting
control gradients at each interval via the sensitivity kick.

Returns a `Matrix{Float64}` of size `[n_drives × (N-1)]` where entry `[k, j]`
is `∂L/∂u_{k,j}` when `λ` is the gradient of the loss w.r.t. the final state.
"""
function ket_vjp(
    data::KnotPointPropagationData,
    λ_seed::AbstractVector;
    use_exact_kick::Bool = true,
    n_quad::Int = 4,
)
    N = length(data.states)
    n_intervals = N - 1
    n_drives = length(data.H_drives)

    grad = zeros(Float64, n_drives, n_intervals)
    λ = ComplexF64.(λ_seed)

    for j = n_intervals:-1:1
        for k = 1:n_drives
            if use_exact_kick
                kick = compute_sensitivity_kick_exact(
                    data.H_drives[k],
                    data.states[j],
                    data.Δts[j],
                    data.eigenvalues[j],
                    data.eigenvectors[j];
                    n_quad = n_quad,
                )
            else
                kick = compute_sensitivity_kick_first_order(
                    data.H_drives[k],
                    data.states[j+1],
                    data.Δts[j],
                )
            end
            grad[k, j] = real(dot(λ, kick))
        end

        # Propagate covector backward: λ ← U_j† · λ
        λ = data.propagators[j]' * λ
    end

    return grad
end

"""
    ket_jvp(data, v; use_exact_kick=true, n_quad=4)

Forward-pass Jacobian-vector product (JVP).

Accumulates the tangent vector `δψ` forward through the knot-point chain:

    δψ_{j+1} = U_j · δψ_j + Σ_k v[k,j] · kick_{k,j}

# Arguments
- `data`: cached propagation data from [`setup_knot_point_propagation`](@ref)
- `v::AbstractMatrix`: tangent directions, size `[n_drives × (N-1)]`

Returns `δψ_N ∈ Vector{ComplexF64}`.
"""
function ket_jvp(
    data::KnotPointPropagationData,
    v::AbstractMatrix;
    use_exact_kick::Bool = true,
    n_quad::Int = 4,
)
    N = length(data.states)
    n_intervals = N - 1
    n_drives = length(data.H_drives)
    n = length(data.states[1])

    δψ = zeros(ComplexF64, n)

    for j = 1:n_intervals
        # δψ_{j+1} = U_j · δψ_j + Σ_k v[k,j] · kick_{k,j}
        δψ = data.propagators[j] * δψ

        for k = 1:n_drives
            if use_exact_kick
                kick = compute_sensitivity_kick_exact(
                    data.H_drives[k],
                    data.states[j],
                    data.Δts[j],
                    data.eigenvalues[j],
                    data.eigenvectors[j];
                    n_quad = n_quad,
                )
            else
                kick = compute_sensitivity_kick_first_order(
                    data.H_drives[k],
                    data.states[j+1],
                    data.Δts[j],
                )
            end
            δψ .+= v[k, j] .* kick
        end
    end

    return δψ
end

# ============================================================================ #
# Tests
# ============================================================================ #
