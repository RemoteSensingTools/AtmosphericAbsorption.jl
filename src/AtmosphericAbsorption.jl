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

using .Architectures
using .LineShapes
using .PartitionFunctions

# Compute backends
export AbstractArchitecture, CPU, GPU, MetalGPU, default_architecture, array_type
# Line-shape math
export AbstractLineProfile, Doppler, Lorentz, Voigt,
       AbstractCPF, HumlicekWeideman32, ErfcxCPF, w, evaluate
# Partition functions
export AbstractPartitionFunction, TabulatedPF, Q, Q_ratio

end # module
