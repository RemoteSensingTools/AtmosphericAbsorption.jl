"""
    Crosssections

The compute core, with two cross-section model types:
- `LineByLineModel` — the KernelAbstractions path: a CPU pre-pass prepares per-line
  parameters, then one kernel sums `S · profile` over lines, specialized on the profile
  and CPF singleton types for CPU/CUDA/Metal.
- `TabulatedCrossSection` — HITRAN `.xsc` panels interpolated on a CPU grid (linear in ν,
  linear in T at the nearest tabulated pressure), for molecules with no line list.
"""
module Crosssections

using KernelAbstractions
using Scratch: @get_scratch!
import Downloads, Serialization
using ..Architectures
using ..Constants
using ..LineShapes: AbstractLineProfile, Doppler, Lorentz, Voigt, SpeedDependentVoigt,
                    Rautian, SpeedDependentRautian, HartmannTran,
                    AbstractCPF, HumlicekWeideman32, evaluate
using ..PartitionFunctions: AbstractPartitionFunction, Q_ratio, pf_name
using ..LineLists: LineDatabase, molecules

export AbstractCrossSectionModel, LineByLineModel, compute_cross_section,
       TabulatedCrossSection, XscBand, read_xsc, load_xsc, fetch_hitran_xsc,
       InterpolationModel, build_interpolation_model,
       save_interpolation_model, load_interpolation_model

include("models.jl")
include("prepare.jl")
include("kernels.jl")
include("xsc.jl")
include("interpolation.jl")

end # module
