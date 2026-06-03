"""
    LineShapes

Pure, GPU-safe line-profile math: the complex probability function `w(z)` and the
area-normalized Doppler/Lorentz/Voigt profiles. No I/O, no allocations — every
function is FT-generic so the same code runs in Float32 (Metal) and Float64.
"""
module LineShapes

using SpecialFunctions: erfcx
using ..Constants

export AbstractCPF, HumlicekWeideman32, ErfcxCPF, w,
       AbstractLineProfile, Doppler, Lorentz, Voigt,
       doppler, lorentz, voigt, evaluate

include("cpf.jl")
include("profiles.jl")

end # module
