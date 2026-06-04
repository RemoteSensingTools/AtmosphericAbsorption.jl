```@raw html
---
layout: home

hero:
  name: "AtmosphericAbsorption.jl"
  text: "GPU-accelerated molecular absorption"
  tagline: "Line-by-line absorption cross-sections for atmospheric radiative transfer — HITRAN & ExoMol, Voigt → Hartmann-Tran, line mixing, MT_CKD/CIA continuum, on CPU, CUDA and Metal. Validated against HAPI."
  image:
    src: /assets/logo.png
    alt: AtmosphericAbsorption logo
  actions:
    - theme: brand
      text: Getting started
      link: /getting_started
    - theme: alt
      text: Line shapes
      link: /line_shapes
    - theme: alt
      text: Benchmarks
      link: /benchmarks

features:
  - title: Many line shapes
    details: Doppler, Lorentz, Voigt, speed-dependent Voigt, and Hartmann-Tran — plus first-order (Rosenkranz) line mixing. Validated against HAPI to ~1e-6.
    link: /line_shapes
  - title: Pluggable data sources
    details: HITRAN (direct download, incl. authenticated non-Voigt parameters) and ExoMol behind one interface; provenance-tracked caching.
    link: /data_sources
  - title: One kernel, every device
    details: The same KernelAbstractions kernel runs on CPU, CUDA, and Metal, in Float32 and Float64 to machine precision. ~55–460× faster than hapi2 on a GPU.
    link: /gpu
---
```

## What it does

`AtmosphericAbsorption.jl` computes molecular absorption cross-sections σ(ν) [cm²/molecule]
from spectroscopic line lists, the engine that turns a HITRAN or ExoMol database into the
optical properties a radiative-transfer model needs.

```julia
using AtmosphericAbsorption

# download CO2 lines for the 1.6 µm band straight from hitran.org
db    = load_lines(HitranPort(; molecule = "CO2", numin = 6300, numax = 6400); mol = 2, iso = 1)
model = LineByLineModel(db, TIPS2021PF(); profile = Voigt(), wing_cutoff = 40.0)

grid  = collect(6300.0:0.01:6400.0)          # cm⁻¹
σ     = compute_cross_section(model, grid, 500.0, 250.0)   # pressure [hPa], temperature [K]
```

## Install

```julia
pkg> add AtmosphericAbsorption
```

Julia 1.10+. GPU support loads automatically when `CUDA.jl` (NVIDIA) or `Metal.jl`
(Apple Silicon) is available.

## Highlights

- **HITRAN & ExoMol** line lists behind a single `load_lines` interface, with on-demand
  downloads cached and version-tracked — plus tabulated HITRAN `.xsc` cross-sections for
  heavy molecules (CFCs, SF₆, …) that have no line list.
- **Voigt → speed-dependent Voigt → Rautian → Hartmann-Tran** line shapes and **line
  mixing**, all validated against [HAPI](https://hitran.org/hapi).
- **Authenticated non-Voigt parameters** (Hartmann-Tran / speed-dependent / line-mixing)
  from HITRANonline.
- **MT_CKD** water-vapor continuum and **CIA** collision-induced absorption.
- **CPU / CUDA / Metal** from one kernel; **Float32 & Float64** to machine precision.
- **Partition functions** (latest TIPS-2021, with TIPS-2017, and ExoMol) driving the
  line-strength temperature dependence.

This package is a clean-slate successor to the absorption module of
[vSmartMOM.jl](https://github.com/RemoteSensingTools/vSmartMOM.jl).
