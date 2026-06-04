# Examples

A cookbook of short, end-to-end recipes. Each example is self-contained and runs against the public API. All wavenumbers and grids are in cm⁻¹, pressures in hPa, temperatures in K, and cross-sections in cm²/molecule.

## 1. HITRAN CO₂ at 1.6 µm

Download the CO₂ lines for the 1.6 µm band directly from hitran.org (public, no API key), build a Voigt line-by-line model, and compute the cross-section on a fine grid. The `HitranPort` keyword form fetches the `.par` file; `load_lines` reads it into a `LineDatabase`, and `partition_function` supplies the iso-aware TIPS-2017 partition sums.

```julia
using AtmosphericAbsorption

# 6300–6400 cm⁻¹ is the 1.6 µm CO₂ band
port  = HitranPort(; molecule="CO2", numin=6300, numax=6400, edition="HITRAN2020")
lines = load_lines(port; mol=2, iso=1, ν_min=6300, ν_max=6400)
pf    = partition_function(port, 2, 1)        # TIPS2017PF

model = LineByLineModel(lines, pf; profile=Voigt(), wing_cutoff=40.0)

grid = collect(6320.0:0.01:6340.0)            # cm⁻¹
σ    = compute_cross_section(model, grid, 1013.25, 296.0)   # hPa, K -> cm²/molecule
```

## 2. ExoMol CO fundamental band

`ExoMolPort` downloads a line list from exomol.com. Pass the HITRAN molecule/isotopologue IDs so the lines can be cross-referenced; `partition_function` returns a `TabulatedPF` built from the ExoMol `.pf` file. Here we compute the CO fundamental band near 4.7 µm.

```julia
using AtmosphericAbsorption

port  = ExoMolPort("CO", "12C-16O", "Li2015"; hitran_mol=5, hitran_iso=1)
lines = load_lines(port; ν_min=2000, ν_max=2250)
pf    = partition_function(port, 5, 1)        # TabulatedPF from the ExoMol .pf

model = LineByLineModel(lines, pf; profile=Voigt(), wing_cutoff=40.0)

grid = collect(2000.0:0.01:2250.0)            # cm⁻¹
σ    = compute_cross_section(model, grid, 1013.25, 296.0)
```

## 3. HITRAN non-Voigt H₂O (speed-dependent Voigt)

The Hartmann-Tran / speed-dependent parameters (`γ2_air`, `δ2_air`, `η`, `νVC`) require an authenticated download. Provide your hitran.org API key with `activate_hitran!` (held in memory only, never written to disk), then `load_hitran_nonvoigt` fills the advanced columns of the `LineDatabase`. Select `SpeedDependentVoigt()` (or `HartmannTran()`) to use them.

```julia
using AtmosphericAbsorption

activate_hitran!("your-key")                  # or set ENV["HITRAN_API_KEY"]

db = load_hitran_nonvoigt("H2O"; numin=3700, numax=3850, min_strength=0.0)
pf = partition_function(HitranPort(; molecule="H2O", numin=3700, numax=3850), 1, 1)

model = LineByLineModel(db, pf; profile=SpeedDependentVoigt(), wing_cutoff=40.0)

grid = collect(3700.0:0.005:3850.0)           # cm⁻¹
σ    = compute_cross_section(model, grid, 1013.25, 296.0)
```

## 4. Run on the GPU

The pipeline runs on NVIDIA CUDA (`GPU()`), Apple Metal (`MetalGPU()`, Float32 only), or the CPU. GPU support activates once you `using CUDA` (or `using Metal`). Set `architecture=GPU()` on the model; the returned cross-section lives on the device, so wrap it with `Array(...)` to copy back to the host. For a fully Float32 pipeline, pass `FT=Float32` to `load_lines`.

```julia
using AtmosphericAbsorption
using CUDA                                     # loads the GPU backend

port  = HitranPort(; molecule="CO2", numin=6300, numax=6400, edition="HITRAN2020")
lines = load_lines(port; mol=2, iso=1, ν_min=6300, ν_max=6400, FT=Float32)
pf    = partition_function(port, 2, 1)

model = LineByLineModel(lines, pf; profile=Voigt(), wing_cutoff=40.0,
                        architecture=GPU())

grid  = collect(Float32, 6320.0:0.01:6340.0)   # cm⁻¹
σ_dev = compute_cross_section(model, grid, 1013.25f0, 296.0f0)
σ     = Array(σ_dev)                            # copy back to host
```

## 5. MT_CKD water-vapor continuum + O₂-O₂ CIA

Continuum and collision-induced absorption are returned as cross-section quantities — multiply by the appropriate number densities and path length yourself to get optical depth. The MT_CKD water-vapor model needs the temperature and the H₂O and dry-air partial pressures (hPa); the O₂-O₂ CIA returns a binary cross-section in cm⁵/molecule² as a function of temperature.

```julia
using AtmosphericAbsorption

grid = collect(2400.0:1.0:2600.0)             # cm⁻¹

# MT_CKD water-vapor continuum (bundled AER v4.2 reference)
tbl  = load_mtckd()
band = build_mtckd_band(tbl, grid)
σ_h2o = h2o_continuum(band, 296.0, 20.0, 993.25)   # T[K], p_h2o, p_dry [hPa]

# O2-O2 collision-induced absorption from a HITRAN .cia file
tcia  = load_cia("O2-O2_2024.cia", grid)
σ_o2  = cia_cross_section(tcia, 296.0)              # cm⁵/molecule²
```

See the [API reference](@ref) for the full list of functions, keywords, and types.
