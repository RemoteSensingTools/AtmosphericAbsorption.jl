# GPU & precision

The cross-section kernel is written once with [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl). The *same* kernel runs on the CPU, on NVIDIA CUDA GPUs, and on Apple Silicon — you choose where by passing an `architecture` to `LineByLineModel`. Nothing else in your script changes.

```julia
architecture = CPU()           # default
architecture = GPU()           # NVIDIA CUDA
architecture = MetalGPU()      # Apple Silicon (Float32 only)
architecture = default_architecture()  # picks GPU if one is available, else CPU
```

## Running on a CUDA GPU

GPU support activates when you load the backend package. For NVIDIA, that is `using CUDA`:

```julia
using AtmosphericAbsorption
using CUDA   # activates the CUDA backend

port = HitranPort("CO2.par"; edition="HITRAN2016")
lines = load_lines(port; mol=2, iso=1, ν_min=6300, ν_max=6400)
pf    = partition_function(port, 2, 1)

model = LineByLineModel(lines, pf;
                        profile=Voigt(),
                        wing_cutoff=40.0,
                        architecture=GPU())

grid = collect(6300.0:0.01:6400.0)   # cm⁻¹
σ_dev = compute_cross_section(model, grid, 1013.25, 296.0)  # pressure in hPa, T in K
```

The result lives **on the device**. To bring it back to host memory for plotting or saving, wrap it with `Array(...)`:

```julia
σ = Array(σ_dev)   # copy back to the CPU; σ in cm²/molecule
```

## Apple Silicon (Metal)

On a Mac with Apple Silicon, load `Metal` and select `MetalGPU()`:

```julia
using AtmosphericAbsorption
using Metal   # activates the Metal backend

model = LineByLineModel(lines, pf;
                        profile=Voigt(),
                        architecture=MetalGPU())
```

`MetalGPU()` is **Float32 only** — Metal does not support double precision. Build a Float32 line list (see below) so the whole pipeline stays single-precision with no surprise promotion.

## Float32 vs Float64

The pipeline is precision-generic. The float type is fixed when you load the line list, via `FT` in `load_lines`:

```julia
lines32 = load_lines(port; mol=2, iso=1, ν_min=6300, ν_max=6400, FT=Float32)
```

A `Float32` line list gives a **Float32-clean** pipeline: the line parameters, the model, and the computed cross-section all stay `Float32`, with no accidental Float64 promotion anywhere along the way. This matters on Metal, where any stray `Float64` would either fail or silently fall back.

`Float64` (the default) and `Float32` agree to machine precision on real cases — the single-precision path is a true fast path, not an approximation. Pick `Float32` for speed and for Metal; keep `Float64` when you want the reference result.

```julia
# Float32 end-to-end on Apple Silicon
using Metal
lines32 = load_lines(port; mol=2, iso=1, ν_min=6300, ν_max=6400, FT=Float32)
model32 = LineByLineModel(lines32, pf; profile=Voigt(), architecture=MetalGPU())
σ32     = Array(compute_cross_section(model32, grid, 1013.25, 296.0))
```

## Speed

The GPU path is **55–460× faster than hapi2** (the HITRAN reference Python implementation), depending on the size of the wavenumber grid and line list — larger problems see the bigger speedups. The GPU result matches the CPU result to machine precision. See [Benchmarks](@ref) for the full comparison and methodology.
