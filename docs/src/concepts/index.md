# [Concepts Overview](@id concepts-overview)

Piccolo.jl solves **quantum optimal control problems** via direct trajectory optimization. This page describes the mathematical problem and how Piccolo.jl's components map onto it.

## The Optimization Problem

Given a quantum system with Hamiltonian

```math
H(\boldsymbol{u}, t) = H_{\text{drift}} + \sum_{i=1}^{m} u_i(t)\, H_{\text{drive},i}
```

we seek piecewise-constant controls ``\boldsymbol{u}_1, \dots, \boldsymbol{u}_N`` that steer the system from an initial state ``x_1`` toward a goal, subject to hardware constraints. Piccolo.jl discretizes this into a finite-dimensional nonlinear program (NLP):

```math
\begin{aligned}
\min_{\boldsymbol{u},\, \Delta t} \quad & Q \cdot \ell(x_N,\, x_{\text{goal}}) \;+\; \sum_{k=1}^{N}\left( R_u \lVert \boldsymbol{u}_k \rVert^2 + R_{du} \lVert \Delta\boldsymbol{u}_k \rVert^2 + R_{ddu} \lVert \Delta^2\boldsymbol{u}_k \rVert^2 \right) \\[6pt]
\text{s.t.} \quad & x_{k+1} = \exp\!\bigl(\Delta t_k \cdot G(\boldsymbol{u}_k)\bigr)\, x_k, \qquad k = 1,\dots, N-1 \\
& x_1 = x_{\text{init}} \\
& \boldsymbol{u}_{\min} \;\leq\; \boldsymbol{u}_k \;\leq\; \boldsymbol{u}_{\max}
\end{aligned}
```

where:

- ``x_k`` is the quantum state at timestep ``k``, represented as a real vector via an [isomorphism](@ref isomorphisms-concept)
- ``G(\boldsymbol{u}) = G_{\text{drift}} + \sum_i u_i\, G_{\text{drive},i}`` is the generator of the dynamics (see below)
- ``\ell(x_N, x_{\text{goal}})`` is an infidelity measure at the final time
- ``\Delta \boldsymbol{u}_k`` and ``\Delta^2 \boldsymbol{u}_k`` are discrete first and second differences of the controls
- ``Q``, ``R_u``, ``R_{du}``, ``R_{ddu}`` are scalar weights

### Generators

The generator ``G`` depends on the type of evolution:

| Evolution | Generator | Equation |
|-----------|-----------|----------|
| **Closed** (Schrödinger) | ``G(\boldsymbol{u}) = -i\bigl(H_{\text{drift}} + \sum_i u_i H_i\bigr)`` | ``\dot{U} = G\, U`` or ``\dot{\psi} = G\, \psi`` |
| **Open** (Lindblad) | ``\mathcal{G}(\boldsymbol{u}) = \mathcal{L}_{\text{drift}} + \sum_i u_i\, \mathcal{L}_i`` | ``\dot{\rho} = \mathcal{G}\, \text{vec}(\rho)`` |

For open systems, the Lindbladian superoperator includes dissipation:

```math
\mathcal{L}[\rho] = -i[H, \rho] + \sum_k \left( L_k \rho L_k^\dagger - \tfrac{1}{2}\{L_k^\dagger L_k, \rho\} \right)
```

### Fidelity Metrics

The terminal cost ``\ell`` is ``1 - F`` for a trajectory-dependent fidelity ``F``:

| Trajectory | Fidelity ``F`` |
|------------|----------------|
| `UnitaryTrajectory` | ``\frac{1}{d^2} \lvert \operatorname{tr}(U_{\text{goal}}^\dagger\, U_N) \rvert^2`` |
| `KetTrajectory` | ``\lvert \langle \psi_{\text{goal}} \mid \psi_N \rangle \rvert^2`` |
| `DensityTrajectory` | ``\operatorname{tr}(\rho_{\text{goal}}\, \rho_N)`` |

### Discretization

The dynamics constraint ``x_{k+1} = \exp(\Delta t_k \cdot G(\boldsymbol{u}_k))\, x_k`` is an **exact matrix exponential** propagator for the piecewise-constant Hamiltonian on each interval ``[t_k, t_{k+1}]``. This preserves unitarity by construction and is computed efficiently via Krylov subspace methods (`ExponentialAction.jl`).

## How Piccolo.jl Maps to the NLP

```
  QuantumSystem          Pulse
  ┌──────────────┐    ┌──────────────┐
  │ H_drift      │    │ u(t), times  │
  │ H_drives     │    │ pulse type   │
  │ drive_bounds │    └──────┬───────┘
  └──────┬───────┘           │
         │                   │
         └─────────┬─────────┘
                   │
                   ▼
             Trajectory                    ← defines  x_init, x_goal, G(u)
          ┌──────────────┐
          │ system       │
          │ pulse        │
          │ goal         │
          └──────┬───────┘
                 │
                 ▼
       Problem Template                    ← assembles the NLP
     ┌───────────────────┐
     │ SmoothPulseProblem│
     │ SplinePulseProblem│
     │ MinimumTimeProblem│
     └──────┬────────────┘
            │
            ▼
     QuantumControlProblem                 ← the NLP
     ┌───────────────────┐
     │ objective  J(z)   │
     │ dynamics   G(u)   │
     │ constraints       │
     └──────┬────────────┘
            │
            ▼
        solve!(qcp)                        ← Ipopt (interior point)
            │
            ▼
     optimized pulse u*(t)
```

### Component Roles

| NLP element | Piccolo.jl component | Reference |
|-------------|---------------------|-----------|
| System Hamiltonian ``H(\boldsymbol{u})`` | [`QuantumSystem`](@ref systems-overview) | [Systems](@ref systems-overview) |
| Control parameterization ``\boldsymbol{u}(t)`` | [`Pulse`](@ref pulses-concept) | [Pulses](@ref pulses-concept) |
| State type, initial/goal, generator ``G`` | [`Trajectory`](@ref trajectories-concept) | [Trajectories](@ref trajectories-concept) |
| Infidelity ``\ell`` and regularization | [`Objective`](@ref objectives-concept) | [Objectives](@ref objectives-concept) |
| Bounds, fidelity/leakage constraints | [`Constraint`](@ref constraints-concept) | [Constraints](@ref constraints-concept) |
| Real vector representation ``x \in \mathbb{R}^n`` | [`Isomorphism`](@ref isomorphisms-concept) | [Isomorphisms](@ref isomorphisms-concept) |
| Subspace embeddings | [`Operators`](@ref operators-concept) | [Operators](@ref operators-concept) |

## Decision Variables

The full NLP decision vector ``z`` is a [`NamedTrajectory`](https://docs.harmoniqs.co/NamedTrajectories.jl/) containing, at each of ``N`` knot points:

| Variable | Symbol | Dimension | Description |
|----------|--------|-----------|-------------|
| State | ``x_k`` | ``n_x`` | Isomorphic quantum state |
| Controls | ``\boldsymbol{u}_k`` | ``m`` | Piecewise-constant amplitudes |
| 1st differences | ``\Delta\boldsymbol{u}_k`` | ``m`` | Control velocity |
| 2nd differences | ``\Delta^2\boldsymbol{u}_k`` | ``m`` | Control acceleration |
| Timestep | ``\Delta t_k`` | ``1`` | Interval duration (optionally free) |

The state dimension ``n_x`` depends on the trajectory type and the system dimension ``d``:

| Trajectory | State | ``n_x`` |
|------------|-------|---------|
| `UnitaryTrajectory` | ``\tilde{U} \in \mathbb{R}^{2d^2}`` | ``2d^2`` |
| `KetTrajectory` | ``\tilde{\psi} \in \mathbb{R}^{2d}`` | ``2d`` |
| `DensityTrajectory` | ``\tilde{\rho} \in \mathbb{R}^{d^2}`` | ``d^2`` (compact) |

## Workflow
A typical Piccolo.jl workflow follows these steps:

```julia
using Piccolo

# 1. Define the quantum system: H(u) = H_drift + u₁ H_x + u₂ H_y
sys = QuantumSystem(PAULIS[:Z], [PAULIS[:X], PAULIS[:Y]], [1.0, 1.0])

# 2. Initial control pulse
pulse = ZeroOrderPulse(0.1 * randn(2, 100), range(0, 10.0, length=100))

# 3. Trajectory = system + pulse + goal
qtraj = UnitaryTrajectory(sys, pulse, GATES[:X])

# 4. Assemble the NLP
qcp = SmoothPulseProblem(qtraj, 100; Q=100.0, R=1e-2)

# 5. Solve (Ipopt interior-point method)
solve!(qcp; max_iter=100)

# 6. Extract results
fidelity(qcp)
optimized_pulse = get_pulse(qcp.qtraj)
```

## Concept Pages

- [Quantum Systems](@ref systems-overview) — Hamiltonian structure, hardware constraints, and system templates
- [Trajectories](@ref trajectories-concept) — State types, goals, and generators
- [Pulses](@ref pulses-concept) — Control parameterizations, how controls vary in time
- [Objectives](@ref objectives-concept) — Fidelity and regularization
- [Constraints](@ref constraints-concept) — Bounds and equality constraints, leakage, fidelity
- [Operators](@ref operators-concept) — Subspace embeddings and lifted operators
- [Isomorphisms](@ref isomorphisms-concept) — Real vector representations

## Architecture

```
Piccolo.jl
├── Quantum          # Quantum mechanical building blocks
│   ├── Systems      # Hamiltonian representations
│   ├── Trajectories # Time evolution containers
│   ├── Pulses       # Control parameterizations
│   ├── Operators    # Embedded and lifted operators
│   └── Isomorphisms # Real vector representations
│
├── Control          # Optimal control framework
│   ├── Problems     # QuantumControlProblem wrapper
│   ├── Objectives   # Fidelity and regularization
│   ├── Constraints  # Bounds and equality constraints
│   └── Templates    # High-level problem constructors
│
└── Visualizations   # Plotting and analysis
    ├── Trajectories # State and control plots
    └── Populations  # Population dynamics
```

## Reexported Packages

| Package | Role |
|---------|------|
| [`DirectTrajOpt`](https://docs.harmoniqs.co/DirectTrajOpt.jl/) | NLP assembly and Ipopt interface |
| [`NamedTrajectories`](https://docs.harmoniqs.co/NamedTrajectories.jl/) | Decision variable storage |
| [`TrajectoryIndexingUtils`](https://docs.harmoniqs.co/TrajectoryIndexingUtils.jl/) | Trajectory slicing and indexing |

These are available when you `using Piccolo` without additional imports.
