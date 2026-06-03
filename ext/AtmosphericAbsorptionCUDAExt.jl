"""
    AtmosphericAbsorptionCUDAExt

Loads NVIDIA CUDA support when CUDA.jl is present: binds the `GPU` architecture to a
CUDA backend and `CuArray`, so the cross-section kernel runs on the GPU unchanged.
"""
module AtmosphericAbsorptionCUDAExt

using AtmosphericAbsorption.Architectures
using CUDA

Architectures.devi(::GPU) = CUDA.CUDABackend(; always_inline = true)
Architectures.array_type(::GPU) = CuArray
Architectures.architecture(::CuArray) = GPU()
Architectures.synchronize_if_gpu(::GPU) = CUDA.synchronize()

function __init__()
    if CUDA.functional()
        Architectures._has_cuda[] = true
        CUDA.allowscalar(false)
    end
end

end # module
