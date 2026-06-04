# Examples

A cookbook of short, end-to-end recipes. Each example is self-contained and runs against the public API. All wavenumbers and grids are in cm⁻¹, pressures in hPa, temperatures in K, and cross-sections in cm²/molecule.

## 1. HITRAN CO₂ at 1.6 µm

Download the CO₂ lines for the 1.6 µm band directly from hitran.org (public, no API key), build a Voigt line-by-line model, and compute the cross-section on a fine grid. The `HitranPort` keyword form fetches the `.par` file; `load_lines` reads it into a `LineDatabase`, and `partition_function` supplies the iso-aware latest-edition (TIPS-2021) partition sums.

```julia
using AtmosphericAbsorption

# 6300–6400 cm⁻¹ is the 1.6 µm CO₂ band
port  = HitranPort(; molecule="CO2", numin=6300, numax=6400, edition="HITRAN2020")
lines = load_lines(port; mol=2, iso=1, ν_min=6300, ν_max=6400)
pf    = partition_function(port, 2, 1)        # TIPS2021PF (latest)

model = LineByLineModel(lines, pf; profile=Voigt(), wing_cutoff=40.0)

grid = collect(6320.0:0.01:6340.0)            # cm⁻¹
σ    = compute_cross_section(model, grid, 1013.25, 296.0)   # hPa, K -> cm²/molecule
```

## 2. ExoMol CO fundamental band (and a HITRAN cross-check)

`ExoMolPort` downloads a line list from exomol.com. Unlike HITRAN, ExoMol stores Einstein-A coefficients rather than precomputed line intensities, so `load_lines` derives each line strength on the fly from the Einstein-A coefficient and the ExoMol partition function. Pass the HITRAN molecule/isotopologue IDs so the species can be cross-referenced; `partition_function` returns a `TabulatedPF` built from the ExoMol `.pf` file. Here we compute the CO fundamental band near 4.7 µm.

```julia
using AtmosphericAbsorption

port  = ExoMolPort("CO", "12C-16O", "Li2015"; hitran_mol=5, hitran_iso=1)
lines = load_lines(port; ν_min=2000, ν_max=2250)
pf    = partition_function(port, 5, 1)        # TabulatedPF from the ExoMol .pf

model = LineByLineModel(lines, pf; profile=Voigt(), wing_cutoff=40.0)

grid = collect(2000.0:0.01:2250.0)            # cm⁻¹
σ    = compute_cross_section(model, grid, 1013.25, 296.0)
```

Because ExoMol and HITRAN encode line intensities completely differently, computing the same band from both is a strong validation of the intensity path. Load the matching HITRAN lines, build an identical model, and overlay the two cross-sections — every downstream call is the same; only the `Port` changes.

```julia
# Same band, same physics, from HITRAN instead of ExoMol
hport = HitranPort(; molecule="CO", numin=2000, numax=2250, edition="HITRAN2020")
hlines = load_lines(hport; mol=5, iso=1, ν_min=2000, ν_max=2250)
hpf    = partition_function(hport, 5, 1)       # TIPS2021PF (latest)
hmodel = LineByLineModel(hlines, hpf; profile=Voigt(), wing_cutoff=40.0)

σ_exomol = compute_cross_section(model,  grid, 1013.25, 296.0)
σ_hitran = compute_cross_section(hmodel, grid, 1013.25, 296.0)

# Compare matched-line intensities directly
for k in eachindex(lines.ν0)
    j = argmin(abs.(hlines.ν0 .- lines.ν0[k]))
    abs(hlines.ν0[j] - lines.ν0[k]) < 0.02 || continue
    relerr = abs(lines.S[k] - hlines.S[j]) / hlines.S[j]
    @info "CO line" ν=lines.ν0[k] S_exomol=lines.S[k] S_hitran=hlines.S[j] relerr
end
```

Line positions coincide and the ExoMol-derived intensities track HITRAN to **0.045 %**. The [Data sources](@ref) page shows the resulting cross-section overlay.

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

## 6. Tabulated cross-sections for a heavy molecule

Molecules like SF₆ or the CFCs have no line list — HITRAN ships them as measured `.xsc` cross-section panels. `load_xsc` reads a `.xsc` file into a `TabulatedCrossSection`, which answers the same `compute_cross_section` call as a line-by-line model, interpolating in wavenumber and temperature onto your grid.

```julia
using AtmosphericAbsorption

path  = fetch_hitran_xsc("SF6_..._N2.xsc")     # filename from hitran.org, or a local path
model = load_xsc(path)

grid = collect(940.0:0.05:960.0)               # cm⁻¹ (the 948 cm⁻¹ SF₆ ν₃ band)
σ    = compute_cross_section(model, grid, 1013.25, 296.0)   # hPa, K -> cm²/molecule
```

See the [API reference](@ref) for the full list of functions, keywords, and types.
