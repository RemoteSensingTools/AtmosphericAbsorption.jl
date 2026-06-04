# Getting started

Welcome to **AtmosphericAbsorption.jl** — a GPU-accelerated package for computing molecular absorption cross-sections directly from spectroscopic line databases (HITRAN, ExoMol) and continuum models (MT-CKD, CIA). It validates against HAPI, the HITRAN reference implementation, to better than 5×10⁻³ on real CO₂/H₂O spectra, and runs the same code on CPU or GPU at single or double precision.

This page gets you from a fresh install to your first cross-section in a few minutes.

## Installation

The package requires **Julia 1.9 or newer**. From the Julia REPL, press `]` to enter the package manager and add it:

```julia
pkg> add AtmosphericAbsorption
```

Then load it like any other package:

```julia
using AtmosphericAbsorption
```

GPU acceleration is optional and loads automatically when you bring in the relevant backend. For NVIDIA cards add `using CUDA`; on Apple silicon add `using Metal`. Until then everything runs on the CPU, which is the default. You do not need a GPU to follow this guide.

## The mental model

Computing a cross-section is a three-step pipeline:

1. **Get a line list.** A *Port* points at a data source — a HITRAN `.par` file, a download from hitran.org, or an ExoMol line list. `load_lines(port)` reads it into a `LineDatabase`. Each port also gives you a matching `partition_function`.
2. **Build a model.** `LineByLineModel` bundles the line list, its partition function, a line-shape `profile` (Voigt by default), and numerical choices (wing cutoff, broadener VMR, architecture).
3. **Evaluate it.** `compute_cross_section(model, grid, pressure, temperature)` returns the cross-section sampled on your wavenumber grid.

```
data source (Port)  ──load_lines──▶  LineDatabase
                                          │
                    partition_function    │
                              └───────┐    │
                                   LineByLineModel
                                          │
                       compute_cross_section(grid, p, T)
                                          │
                                          ▼
                                σ  [cm²/molecule]
```

### Units (used everywhere)

| Quantity | Unit |
|---|---|
| Wavenumber grid, `ν0`, `numin`/`numax` | cm⁻¹ |
| Pressure | hPa |
| Temperature | K |
| Cross-section `σ` | cm²/molecule |
| `wing_cutoff` | cm⁻¹ |

## A first end-to-end example

Let's download CO₂ from hitran.org over 6300–6400 cm⁻¹ and compute its cross-section at a mid-troposphere state of **p = 500 hPa, T = 250 K**. The HITRAN download is public — no account or API key is needed.

```julia
using AtmosphericAbsorption

# 1. A Port that fetches CO2 lines from hitran.org for our spectral window
port = HitranPort(; molecule="CO2", numin=6300, numax=6400, edition="HITRAN2020")

# 2. Read the lines and grab the matching partition function
lines = load_lines(port)                         # -> LineDatabase{Float64}
pf    = partition_function(port, 2, 1)           # CO2 is HITRAN molecule 2, isotopologue 1

# 3. Build the line-by-line model (Voigt profile by default)
model = LineByLineModel(lines, pf; profile=Voigt(), wing_cutoff=40.0)

# 4. Evaluate on a 0.01 cm⁻¹ grid at p = 500 hPa, T = 250 K
grid = collect(6300.0:0.01:6400.0)               # cm⁻¹
σ    = compute_cross_section(model, grid, 500.0, 250.0)   # cm²/molecule
```

`σ` is a `Vector` the same length as `grid`, giving the absorption cross-section in **cm²/molecule** at each wavenumber. The line strengths are automatically scaled from their 296 K reference values to 250 K using the model's partition function, so you do not have to apply any temperature correction yourself.

<iframe title="CO2 cross-section" src="assets/plots/co2_crosssection.html" loading="lazy" style="width:100%;height:520px;border:1px solid var(--vp-c-divider);border-radius:8px;"></iframe>

The forest of peaks is the CO₂ line manifold in this band; each line is broadened by pressure (Lorentz wings) and temperature (Doppler core), which the Voigt profile combines. Lower the pressure and the lines sharpen toward their Doppler limit; raise it and they broaden and blend.

## Loading from a local `.par` file

If you already have a HITRAN `.par` file on disk — for example one you downloaded earlier — point a `HitranPort` at the file instead of requesting a download. Everything downstream is identical.

```julia
using AtmosphericAbsorption

# A Port backed by a local HITRAN .par file
port = HitranPort("CO2.par"; edition="HITRAN2016")

lines = load_lines(port; ν_min=6300, ν_max=6400, min_strength=0.0)
pf    = partition_function(port, 2, 1)

model = LineByLineModel(lines, pf; profile=Voigt(), wing_cutoff=40.0)

grid = collect(6300.0:0.01:6400.0)
σ    = compute_cross_section(model, grid, 500.0, 250.0)
```

The `load_lines` keywords let you trim the line list as you read it: `ν_min`/`ν_max` clip the spectral range, `min_strength` drops weak lines below a chosen intensity, and `mol`/`iso` restrict to a specific molecule or isotopologue. Passing `FT=Float32` produces a Float32 line list and a fully Float32 pipeline — handy for GPU runs and matched to the Float64 result to machine precision.

## What you get back

`compute_cross_section` returns a `Vector` of cross-sections in **cm²/molecule**, one value per grid point, living on the model's device. On the CPU that is an ordinary `Array`. If you build the model with `architecture=GPU()` (or `MetalGPU()` on Apple silicon), the result lives on the GPU — wrap it with `Array(σ)` to copy it back to the host for plotting or saving.

To turn a cross-section into an optical depth, multiply by the absorber number density and the path length yourself; the package deliberately returns the pure cross-section so you stay in control of the atmospheric bookkeeping.

## Where to next

- Swap `profile=Voigt()` for `Doppler()`, `Lorentz()`, `SpeedDependentVoigt()`, or `HartmannTran()` to explore different line-shape physics.
- Use `ExoMolPort` to pull line lists from exomol.com for molecules and temperature ranges beyond HITRAN's coverage.
- Add the MT-CKD water-vapor continuum (`load_mtckd`) or collision-induced absorption (`load_cia`) on top of the line-by-line cross-section.
- Run on a GPU by passing `architecture=GPU()` after `using CUDA` — the same pipeline, no code changes.
