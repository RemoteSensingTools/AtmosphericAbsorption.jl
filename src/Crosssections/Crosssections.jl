"""
    Crosssections

The compute core: `LineByLineModel`, the CPU pre-pass that prepares per-line
parameters, and the single KernelAbstractions kernel that sums `S · profile` over
lines. Specializes on the profile and CPF singleton types for CPU/CUDA/Metal.
"""
module Crosssections

using KernelAbstractions
using ..Architectures
using ..Constants
using ..LineShapes: AbstractLineProfile, Doppler, Lorentz, Voigt,
                    AbstractCPF, HumlicekWeideman32, evaluate
using ..PartitionFunctions: AbstractPartitionFunction, Q_ratio
using ..LineLists: LineDatabase

export AbstractCrossSectionModel, LineByLineModel, compute_cross_section

include("models.jl")
include("prepare.jl")
include("kernels.jl")

end # module
