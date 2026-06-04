# Continuum Absorption

Discrete spectral lines are not the whole story. Two slowly varying "continuum" contributions matter for accurate atmospheric radiative transfer, and `AtmosphericAbsorption.jl` exposes both:

- **MT_CKD water-vapor continuum** — the self- and foreign-broadened water-vapor continuum (the far-wing and weak-interaction absorption that the line-by-line sum truncates away), using the bundled AER v4.2 reference data.
- **Collision-induced absorption (CIA)** — transient absorption by colliding pairs (e.g. O₂–O₂, N₂–N₂) that have no permanent dipole, read from HITRAN `.cia` files.

Both are returned as **cross-section quantities**. To get optical depth you multiply by the appropriate number densities and the path length yourself — that integration lives in your RT code, not here. The two continua use different cross-section conventions because they involve different numbers of molecules (see the units notes below).

## Water-vapor continuum (MT_CKD)

The workflow has three steps: load the bundled table, build a band interpolated onto your wavenumber grid, then evaluate at a given temperature and the relevant partial pressures.

```julia
using AtmosphericAbsorption

# Wavenumber grid in cm⁻¹
grid = collect(2350.0:0.01:2450.0)

# 1. Load the bundled AER MT_CKD v4.2 reference table
tbl = load_mtckd()

# 2. Interpolate the continuum onto your grid
band = build_mtckd_band(tbl, grid)

# 3. Evaluate the cross-section [cm²/molecule]
T      = 250.0    # temperature [K]
p_h2o  = 12.0     # water-vapor partial pressure [hPa]
p_dry  = 988.0    # dry-air partial pressure [hPa]
σ_h2o  = h2o_continuum(band, T, p_h2o, p_dry)   # Vector, cm²/molecule
```

`h2o_continuum` returns the water-vapor continuum cross-section in **cm²/molecule**, one value per grid point. It combines the self-continuum (scaling with the water-vapor partial pressure `p_h2o`) and the foreign-continuum (scaling with the dry-air partial pressure `p_dry`).

To convert to optical depth, multiply by the water-vapor number density and the path length over each layer in your RT code.

::: tip Temperature dependence
The water-vapor **self-continuum strengthens as temperature drops** — it has a strong negative temperature coefficient. Evaluating it at a single representative temperature for a thick, cold layer will misestimate the absorption; sample it per layer at the layer temperature.
:::

## Collision-induced absorption (CIA)

CIA is read directly from a HITRAN `.cia` file. As with MT_CKD, you load and interpolate onto your grid once, then evaluate per temperature.

```julia
using AtmosphericAbsorption

grid = collect(7000.0:0.1:8500.0)   # cm⁻¹

# Load and interpolate a HITRAN .cia file onto the grid
tcia = load_cia("O2-O2_2024.cia", grid)

# Binary cross-section at a given temperature [cm⁵/molecule²]
T   = 280.0
σ_cia = cia_cross_section(tcia, T)   # Vector, cm⁵/molecule²
```

`cia_cross_section` returns the **binary** cross-section σ_AB(ν, T) in **cm⁵/molecule²**. The extra factors of length and molecule count (relative to the cm²/molecule of a line-by-line cross-section) reflect that CIA scales with the product of *two* number densities. For an O₂–O₂ pair the layer optical depth is

```
τ(ν) = σ_AB(ν, T) · n_O₂² · L
```

where `n_O₂` is the O₂ number density [molecules/cm³] and `L` is the path length [cm]. For an unlike pair A–B, use `n_A · n_B` in place of `n_O₂²`.

::: tip Multi-band .cia files and temperature
HITRAN `.cia` files often contain several spectral bands, each tabulated on its own temperature grid. `load_cia` and `cia_cross_section` perform the temperature interpolation **per frequency**, so files spanning multiple bands (with differing temperature coverage) are handled correctly — each grid point is interpolated using the temperatures available for its band.
:::

## Units summary

| Quantity | Function | Units |
|----------|----------|-------|
| Wavenumber grid | — | cm⁻¹ |
| Temperature `T` | — | K |
| Partial pressures `p_h2o`, `p_dry` | — | hPa |
| MT_CKD water-vapor continuum | `h2o_continuum` | cm²/molecule |
| CIA binary cross-section | `cia_cross_section` | cm⁵/molecule² |

In both cases the package returns a cross-section, never an optical depth. Combine with number densities (× path length) in your RT code to build τ(ν).
