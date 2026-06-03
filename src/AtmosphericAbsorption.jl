"""
    AtmosphericAbsorption

Clean-slate, GPU-accelerated molecular absorption cross-sections for atmospheric
radiative transfer. Pluggable line-list sources (HITRAN, ExoMol, …), advanced line
shapes (Voigt → speed-dependent Voigt → Hartmann-Tran), line mixing, and continua.

Phase 0 surface: `Architectures`, `Constants`, `LineShapes`, `PartitionFunctions`.
"""
module AtmosphericAbsorption

include("Architectures/Architectures.jl")
include("Constants/Constants.jl")
include("LineShapes/LineShapes.jl")
include("PartitionFunctions/PartitionFunctions.jl")
include("LineLists/LineLists.jl")
include("Crosssections/Crosssections.jl")

using .Architectures
using .LineShapes
using .PartitionFunctions
using .LineLists
using .Crosssections

# Compute backends
export AbstractArchitecture, CPU, GPU, MetalGPU, default_architecture, array_type
# Line-shape math
export AbstractLineProfile, Doppler, Lorentz, Voigt,
       AbstractCPF, HumlicekWeideman32, ErfcxCPF, w, evaluate
# Partition functions
export AbstractPartitionFunction, TabulatedPF, Q, Q_ratio
# Line lists
export AbstractLineListPort, LineDatabase, SourceMetadata,
       load_lines, partition_function, source_metadata
# Cross-section compute core
export AbstractCrossSectionModel, LineByLineModel, compute_cross_section

end # module
