"""
    Architectures

CPU/GPU compute-backend abstraction (vendored, trimmed from vSmartMOM/Oceananigans).
Concrete GPU bindings (`devi`, `array_type`) are injected by the CUDA/Metal package
extensions when those backends are loaded.
"""
module Architectures

using KernelAbstractions

export AbstractArchitecture, CPU, GPU, MetalGPU,
       devi, array_type, architecture, default_architecture,
       synchronize_if_gpu, has_cuda, has_metal

# Set to true by the CUDA/Metal extensions when the backend is functional.
const _has_cuda  = Ref(false)
const _has_metal = Ref(false)
has_cuda()  = _has_cuda[]
has_metal() = _has_metal[]

"""Abstract supertype for compute architectures."""
abstract type AbstractArchitecture end

"""Single-threaded / multithreaded CPU execution."""
struct CPU <: AbstractArchitecture end

"""NVIDIA CUDA GPU (bindings from the CUDA extension)."""
struct GPU <: AbstractArchitecture end

"""Apple-Silicon GPU via Metal (Float32 only; bindings from the Metal extension)."""
struct MetalGPU <: AbstractArchitecture end

@inline devi(::CPU) = KernelAbstractions.CPU()
# devi(::GPU)/devi(::MetalGPU) are injected by the backend extensions.

@inline array_type(::CPU) = Array
# array_type(::GPU)/array_type(::MetalGPU) are injected by the backend extensions.

@inline architecture(::Array) = CPU()

"""Return `GPU()` if CUDA is functional, `MetalGPU()` if Metal is, else `CPU()`."""
default_architecture() = has_cuda() ? GPU() : has_metal() ? MetalGPU() : CPU()

"""Synchronize the given backend after a kernel launch; a no-op on CPU. GPU/Metal
methods are provided by the backend extensions."""
@inline synchronize_if_gpu(::AbstractArchitecture) = nothing

end # module
