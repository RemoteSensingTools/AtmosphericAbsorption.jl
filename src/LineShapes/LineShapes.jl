"""
    LineShapes

Pure, GPU-safe line-profile math: the complex probability function `w(z)` and the
area-normalized Doppler/Lorentz/Voigt profiles. No I/O, no allocations — every
function is FT-generic so the same code runs in Float32 (Metal) and Float64.
"""
module LineShapes

using SpecialFunctions: erfcx
using ..Constants

# `evaluate` is the kernel's internal per-line dispatch entry point — intentionally
# NOT exported; the public line-shape API is the bare profile functions below.
export AbstractCPF, HumlicekWeideman32, ErfcxCPF, w,
       AbstractLineProfile, Doppler, Lorentz, Voigt, SpeedDependentVoigt, HartmannTran,
       doppler, lorentz, voigt, pcqsdhc

include("cpf.jl")
include("profiles.jl")
include("pcqsdhc.jl")

end # module
