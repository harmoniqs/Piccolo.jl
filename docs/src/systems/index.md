# [Quantum Systems](@id systems-overview)

Piccolo.jl models quantum control using the **bilinear control Hamiltonian**:

```math
H(\boldsymbol{u}, t) = H_{\text{drift}} + \sum_{i=1}^{m} u_i(t)\, H_{\text{drive},i}
```

where ``H_{\text{drift}}`` is the always-on Hamiltonian,
``\{H_{\text{drive},i}\}`` are the controllable interactions, and
``\boldsymbol{u}(t) \in \mathbb{R}^m`` are the control amplitudes.

For optimization, the Hamiltonian is converted to a **generator**:

```math
G(\boldsymbol{u}) = -i\!\left( H_{\text{drift}} + \sum_{i=1}^{m} u_i\, H_{\text{drive},i} \right)
```

so that ``\dot{\psi} = G\,\psi`` (ket) or ``\dot{U} = G\,U`` (unitary propagator).
For open systems, the Lindbladian superoperator ``\mathcal{G}(\boldsymbol{u})``
has the same bilinear structure (see [Isomorphisms](@ref isomorphisms-concept)).

## System Types

| Type | Dynamics | Use Case |
|------|----------|----------|
| `QuantumSystem` | ``\dot{U} = G\,U`` | Closed systems, gate synthesis |
| `OpenQuantumSystem` | ``\dot{\vec\rho} = \mathcal{G}\,\vec\rho`` | Dissipative systems (Lindblad) |
| `CompositeQuantumSystem` | ``\mathcal{H} = \mathcal{H}_1 \otimes \mathcal{H}_2`` | Multi-subsystem setups |

## Platform Templates

Each page below documents the physical Hamiltonian, constructor API, and a worked
optimization example.

| Platform | Template(s) | System Type | Key Physics |
|----------|-------------|-------------|-------------|
| [Transmon Qubits](@ref transmon-systems) | `TransmonSystem`, `MultiTransmonSystem`, `TransmonCavitySystem` | `QuantumSystem` / `CompositeQuantumSystem` | Anharmonicity, dipole coupling, dispersive readout |
| [Trapped Ions](@ref trapped-ion-systems) | `IonChainSystem`, `RadialMSGateSystem` | `QuantumSystem` | Motional modes, Lamb-Dicke coupling, Mølmer-Sørensen gates |
| [Rydberg Atoms](@ref rydberg-atom-systems) | `RydbergChainSystem` | `QuantumSystem` | Van der Waals blockade, global/local drives |
| [Cat Qubits](@ref cat-qubit-systems) | `CatSystem` | `OpenQuantumSystem` | Two-photon drive, Kerr nonlinearity, photon loss |
| [Silicon Spins](@ref silicon-spin-systems) | *(coming soon)* | `QuantumSystem` | Exchange-only qubits, detuning control |

## Construction

All templates produce a `QuantumSystem` (or subtype) that plugs directly into
[Trajectories](@ref trajectories-concept) and [Problem Templates](@ref problem-templates-overview):

```julia
using Piccolo

sys = TransmonSystem(levels=3, δ=0.2, drive_bounds=[0.2, 0.2])
pulse = ZeroOrderPulse(0.05 * randn(2, 100), collect(range(0, 10, 100)))
qtraj = UnitaryTrajectory(sys, pulse, GATES[:X])
qcp = SmoothPulseProblem(qtraj, 100; Q=100.0)
solve!(qcp)
```

## Custom Systems

For platforms not covered by a template, build directly from matrices:

```julia
H_drift = my_drift_matrix
H_drives = [my_drive_1, my_drive_2]
drive_bounds = [1.0, 1.0]

sys = QuantumSystem(H_drift, H_drives, drive_bounds)
```

All Hamiltonians must be Hermitian (``H = H^\dagger``); Piccolo.jl validates this
at construction.

## See Also

- [Trajectories](@ref trajectories-concept) — Combining systems with pulses and goals
- [Pulses](@ref pulses-concept) — Control parameterizations
- [Problem Templates](@ref problem-templates-overview) — Setting up optimization
- [Isomorphisms](@ref isomorphisms-concept) — How complex matrices become real generators
