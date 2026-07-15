module Frames

using LinearAlgebra
using SparseArrays
using ..QuantumSystems       # QuantumSystem, CompositeQuantumSystem, AbstractDrive, LinearDrive
import ..QuantumSystems: drive_matrix   # materialize an AbstractDrive's operator matrix
using ..QuantumObjectUtils   # annihilate
using ..LiftedOperators      # lift_operator

include("frame_types.jl")
include("frame_transforms.jl")

# exports added incrementally in A1–A3:
# export AbstractFrame, RotatingFrame, LabFrame, FrameSpec
# export to_lab_frame, to_rotating_frame

end
