"""
    Frames

System-agnostic frame transformations for `QuantumSystem`s. A frame is defined by
its generator ``H_f = \\sum_i \\omega_i n_i`` (per-subsystem frame frequencies).

- [`RotatingFrame`](@ref) / [`LabFrame`](@ref) — frame descriptors.
- [`FrameSpec`](@ref) — the subsystem/drive metadata (`number_ops`, `drive_map`,
  `drive_ops`) the transforms read; auto-derivable from a `CompositeQuantumSystem`
  or a matrix-drive `QuantumSystem`.
- [`to_lab_frame`](@ref) — rotating → **lab**: adds the frame generator back into
  the drift and reconstructs each quadrature pair `(u_x, u_y)` as a carrier-
  modulated real field ``\\Omega\\cos(\\omega_d t + \\phi)(a+a^\\dagger)``. Returns a
  `time_dependent=true` `QuantumSystem` rolled through the `Rollouts` **ODE** path
  (not `SplineIntegrator`).
- [`to_rotating_frame`](@ref) — the RWA inverse.

# Example
```julia
levels = 3; a = annihilate(levels)
Hx = 0.5 * Matrix(a + a'); Hy = 0.5 * Matrix(im * (a' - a))
sys_rot = QuantumSystem(0.5 * (-0.2) * Matrix(a' * a' * a * a), [Hx, Hy], [0.05, 0.05])
spec = FrameSpec(number_ops = [Matrix(a' * a)],
                 drive_map = [(1, :x, +1.0), (1, :y, +1.0)], drive_ops = [Hx, Hy])
sys_lab = to_lab_frame(sys_rot, RotatingFrame(4.0), spec)     # carrier-modulated device
# roll the lab-frame (time-dependent) system via the ODE path:
qtraj = UnitaryTrajectory(sys_lab, pulse, U_goal)             # default MagnusAdapt4()
```
"""
module Frames

using LinearAlgebra
using SparseArrays
using TestItems             # @testitem (no-op during normal package load)
using ..QuantumSystems       # QuantumSystem, CompositeQuantumSystem, AbstractDrive, LinearDrive
import ..QuantumSystems: drive_matrix   # materialize an AbstractDrive's operator matrix
using ..QuantumObjectUtils   # annihilate
using ..LiftedOperators      # lift_operator

include("frame_types.jl")
include("frame_transforms.jl")

# exports added incrementally in A1–A3:
export AbstractFrame, RotatingFrame, LabFrame, FrameSpec
export to_lab_frame
export to_rotating_frame

end
