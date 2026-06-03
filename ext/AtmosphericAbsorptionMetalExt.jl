"""
    AtmosphericAbsorptionMetalExt

Loads Apple-Silicon Metal support when Metal.jl is present: binds the `MetalGPU`
architecture to a Metal backend and `MtlArray`. Metal is Float32-only, which the
FT-generic kernel already supports — build the model with `FT=Float32`.
"""
module AtmosphericAbsorptionMetalExt

using AtmosphericAbsorption.Architectures
using Metal

Architectures.devi(::MetalGPU) = Metal.MetalBackend()
Architectures.array_type(::MetalGPU) = MtlArray
Architectures.architecture(::MtlArray) = MetalGPU()
Architectures.synchronize_if_gpu(::MetalGPU) = Metal.synchronize()

function __init__()
    if Metal.functional()
        Architectures._has_metal[] = true
        Metal.allowscalar(false)
    end
end

end # module
