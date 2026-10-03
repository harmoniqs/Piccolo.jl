export NonHermitianExponentialIntegrator

"""
    NonHermitianExponentialIntegrator{T<:AbstractQuantumTrajectory} <: AbstractExponentialIntegrator

Exponential integrator for non-Hermitian generators — specifically, Lindbladian
superoperators acting on vectorized density matrices. Mirrors `HermitianExponentialIntegrator`
in sparsity/Jacobian/Hessian scaffolding; differs only in the forward `exp` call, which
uses `exp_generator!` (scaling-and-squaring) instead of `exp_eigen!` (eigendecomposition).

Type parameter `T`:
- `DensityTrajectory`: single-density evolution
- `MultiDensityTrajectory`: ensemble density evolution

Dropped relative to `HermitianExponentialIntegrator`:
- `H::Function` field (Liouvillian is not factored as -iH)
- The `exp_eigen!` buffer set; replaced with `G_buf`, `expG_buf` for `exp_generator!`.

# Fields
- `G::Function`: Generator function G(u) returning the (real, possibly non-Hermitian)
  superoperator in isomorphic form
- `x_names::Vector{Symbol}`: Names of state variables (single element for non-ensemble)
- `u_name::Symbol`: Name of control variable
- `x_dim::Int`: Total constraint dimension per knot point
- `u_dim::Int`: Dimension of control (without globals)
- `var_dim::Int`: Total variable dimension (states + controls + time)
- `dim::Int`: Total constraint dimension
- `statedim::Int`: Dimension of the state space; for `DensityTrajectory`, `statedim = n^2`
  (compact vectorized density matrix)
- `∂ℰs::Vector{SparseMatrixCSC{Float64, Int}}`: Preallocated Jacobian structures per knot point
- `μ∂²ℰs::Vector{SparseMatrixCSC{Float64, Int}}`: Preallocated Hessian structures per knot point
- `z_dim::Int`: Dimension of single knot point
- `var_comps::Vector{Int}`: Component indices for variables
- `Id::Union{Nothing, Matrix{Float64}}`: Preallocated identity matrix (trajectory-specific)
- `∂u∂Δt_bufs::Union{Nothing, Vector{Vector{Float64}}}`: Per-thread buffers for the
  u-Δt cross term (trajectory-specific; one entry per `Threads.maxthreadid()`)
- `∂²u_bufs::Union{Nothing, Vector{Matrix{Float64}}}`: Per-thread buffers for the
  u-u Hessian block (trajectory-specific; one entry per `Threads.maxthreadid()`)
- `G_bufs::Vector{Matrix{Float64}}`: Per-thread preallocated generator buffers, each size (statedim, statedim)
- `expG_bufs::Vector{Matrix{Float64}}`: Per-thread preallocated exp outputs, same size as `G_bufs[i]`
- `global_names::Vector{Symbol}`: Names of global (time-invariant) variables
- `global_dim::Int`: Dimension of global variables
- `gauss_newton::Bool`: Whether the integrator is configured for Gauss-Newton Hessian
  approximation

# Per-thread buffer set
`G_bufs` / `expG_bufs` are stored as `Vector{Matrix{Float64}}` with one entry per
OS thread (sized to `Threads.maxthreadid()` at construction). This enables
`Threads.@threads` in the per-knot Jacobian/Hessian loops without racing on
the `exp_generator!` scratch; concurrent knot points on different threads get
disjoint buffers by indexing `ℰ.G_bufs[Threads.threadid()]`.
"""
struct NonHermitianExponentialIntegrator{T<:AbstractQuantumTrajectory} <:
       AbstractExponentialIntegrator
    G::Function
    x_names::Vector{Symbol}
    u_name::Symbol
    x_dim::Int
    u_dim::Int
    var_dim::Int
    dim::Int
    statedim::Int                 # replaces `ketdim`: for DensityTrajectory, statedim = n^2 (compact iso)
    ∂ℰs::Vector{SparseMatrixCSC{Float64,Int}}
    μ∂²ℰs::Vector{SparseMatrixCSC{Float64,Int}}
    z_dim::Int
    var_comps::Vector{Int}
    # Type-specific caches
    Id::Union{Nothing,Matrix{Float64}}
    ∂u∂Δt_bufs::Union{Nothing,Vector{Vector{Float64}}}
    ∂²u_bufs::Union{Nothing,Vector{Matrix{Float64}}}
    # Per-thread preallocated work buffers for exp_generator! (one entry per thread,
    # sized to Threads.maxthreadid()).
    G_bufs::Vector{Matrix{Float64}}         # each size (statedim, statedim)
    expG_bufs::Vector{Matrix{Float64}}      # each size (statedim, statedim)
    # Global variable support
    global_names::Vector{Symbol}
    global_dim::Int
    # Gauss-Newton flag
    gauss_newton::Bool
end

# Convenience accessors for single state name (non-ensemble case)
x_name(ℰ::NonHermitianExponentialIntegrator) = ℰ.x_names[1]
single_state_dim(ℰ::NonHermitianExponentialIntegrator) = ℰ.x_dim ÷ length(ℰ.x_names)

# Implementations (constructor, evaluate!, jacobian!, hessian_of_lagrangian!)
# live in density / multidensity files added in Tasks 9 and 10.
