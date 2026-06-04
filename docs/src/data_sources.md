# Data sources

AtmosphericAbsorption.jl reads spectroscopic line lists from several providers — HITRAN and ExoMol — through a single abstraction called a **Port**. A Port knows how to locate, download, cache, and parse one kind of data source. Whatever the origin, `load_lines` returns a uniform `LineDatabase`, so the rest of the pipeline (partition functions, line shapes, cross-sections) is identical regardless of where the lines came from.

```julia
port  = HitranPort("CO2.par"; edition="HITRAN2016")
lines = load_lines(port; mol=:CO2, ν_min=6300, ν_max=6400)   # -> LineDatabase, partition attached
```

The returned `LineDatabase` already carries the right partition function (`lines.partition` — TIPS-2021 for HITRAN, the ExoMol `.pf` for ExoMol), so you build a model straight from it — no partition bookkeeping in the common path:

```julia
model = LineByLineModel(lines; profile=Voigt())     # uses lines.partition automatically
```

Downloads are written to a scratch cache and never re-fetched if already present. Alongside each cached artifact, a `.meta.toml` records its provenance — source, molecule, isotopologue, spectral window, edition/version — so a result can always be traced back to the exact data that produced it. `source_metadata(port, mol, iso)` returns this information programmatically.

## Species notation

Molecules and isotopologues are chosen with a readable, backend-agnostic notation on top of the raw HITRAN integer ids — the same names work for HITRAN and ExoMol:

```julia
load_lines(port; mol=:CO2)            # symbol  (≡ "CO2" ≡ 2)
load_lines(port; mol=:CO2, iso=:main) # principal isotopologue (≡ iso=1)
load_lines(port; mol=:ALL)            # every molecule in the file (≡ -1) — see what a band picks up
```

`:ALL` is the readable spelling of the `-1` "all" sentinel. After loading a multi-molecule band, `molecules(db)` lists what's present and `db[:CO2]` / `db[mask]` subsets it (preserving the partition function) so you can plot each species independently:

```julia
db = load_lines(HitranPort("scene.par"); mol=:ALL, ν_min=2000, ν_max=2400)
for m in molecules(db)                # e.g. [:CO2, :H2O, :N2O]
    σ = compute_cross_section(LineByLineModel(db[m]), grid, 1013.25, 296.0)
end
```

The mapping is **overridable** — adopt your own accepted codes without touching the package:

```julia
register_molecule!(:HDO, 1)                  # name an alias / a custom species
register_isotopologue!(:CO2, Symbol("636"), 2)
```

The full set of names is the [Species & isotopologues](isotopologues.md) table. `molecule_symbol(id)` and `isotopologue(mol, iso)` give the formula back for display (the model's `show` uses them).

## 1. HITRAN from a local `.par` file

If you already have a HITRAN `.par` file, wrap it in a `HitranPort` and load. The `edition` keyword records which HITRAN edition the file came from (it is provenance metadata, not a download trigger).

```julia
using AtmosphericAbsorption

port  = HitranPort("CO2.par"; edition="HITRAN2016")
lines = load_lines(port; mol=:CO2, iso=:main, ν_min=6300, ν_max=6400, min_strength=0.0)
```

`load_lines` accepts windowing and filtering keywords:

- `mol`, `iso` — select a molecule / isotopologue (a symbol `:CO2`/`:main`, a string, an integer id, or `:ALL` for all).
- `ν_min`, `ν_max` — spectral window in cm⁻¹.
- `min_strength` — drop lines weaker than this threshold (cm⁻¹/(molecule·cm⁻²)).
- `FT` — element type of the returned line list, `Float64` (default) or `Float32`.

## 2. HITRAN direct download

To fetch lines straight from hitran.org, construct a `HitranPort` with keyword arguments instead of a file path. The line list is downloaded once into the scratch cache (with a `.meta.toml` alongside) and reused thereafter; pass `force=true` to re-download. This endpoint is public and needs no API key.

```julia
using AtmosphericAbsorption

port  = HitranPort(; molecule="CO2", numin=6300, numax=6400, edition="HITRAN2020")
lines = load_lines(port)
pf    = partition_function(port, 2, 1)
```

If you only want the raw `.par` on disk, `fetch_hitran` downloads it and returns the file path:

```julia
path = fetch_hitran("CO2"; numin=6300, numax=6400, edition="HITRAN2020")
```

The default HITRAN download provides the standard Voigt line parameters (air/self half-widths, pressure shifts, and temperature exponents). It does **not** carry first-order line-mixing coefficients — those, along with the speed-dependent and Hartmann–Tran parameters, come from the authenticated endpoint below, so a basic `.par` load applies no line mixing.

## 3. HITRAN authenticated NON-Voigt parameters

HITRAN's advanced line-shape parameters — speed-dependent broadening, Hartmann–Tran, and full line-mixing — live behind an authenticated endpoint. You need a free API key from your profile page on hitran.org.

The key is held in memory only for the duration of your session. It is never written to disk and must never be committed to a repository. Provide it one of two ways:

```julia
using AtmosphericAbsorption

# Option A: activate explicitly (key kept in memory only, never stored)
activate_hitran!("your-key")

# Option B: set the environment variable before starting Julia
#   export HITRAN_API_KEY="your-key"
```

With the key active, `load_hitran_nonvoigt` fetches the advanced parameters and fills the corresponding columns of the `LineDatabase` — the speed-dependent / Hartmann–Tran terms (`γ2_air`, `δ2_air`, `η`, `νVC`), the first-order line-mixing coefficient `Y_LM` (from HITRAN's `Y_HT_air_296`, falling back to `Y_SDV_air_296`), and the self-broadening temperature exponent `n_self` and pressure shift `δ_self` (absent from the basic `.par` format). Line mixing then folds into every collisional profile automatically for the lines that carry it (e.g. CO₂ and O₂ Q-branches); the self parameters refine the width and shift at non-zero broadener `vmr` (`n_self` defaults to `n_air`, `δ_self` to 0 where HITRAN omits them):

```julia
db = load_hitran_nonvoigt("H2O"; numin=3700, numax=3850, min_strength=0.0, FT=Float64)
```

To actually use those parameters, choose a non-Voigt profile when building the model — `SpeedDependentVoigt()` for speed-dependent broadening, or `HartmannTran()` for the full Hartmann–Tran profile:

```julia
pf    = partition_function(HitranPort(; molecule="H2O", numin=3700, numax=3850), 1, 1)
model = LineByLineModel(db, pf; profile=HartmannTran(), wing_cutoff=40.0)

grid  = collect(3700.0:0.01:3850.0)      # cm⁻¹
σ     = compute_cross_section(model, grid, 1013.25, 296.0)   # p[hPa], T[K]
```

These advanced profiles reproduce HAPI's reference `pcqsdhc` evaluation to ~1e-6.

## 4. ExoMol

ExoMol line lists are loaded through an `ExoMolPort`. ExoMol stores transitions with Einstein-A coefficients rather than precomputed intensities, so line strengths are computed on the fly from the Einstein-A coefficient and the ExoMol partition function. The `hitran_mol` / `hitran_iso` keywords map the ExoMol species onto the corresponding HITRAN molecule and isotopologue IDs, so the downstream pipeline is unchanged.

```julia
using AtmosphericAbsorption

port  = ExoMolPort("CO", "12C-16O", "Li2015"; hitran_mol=5, hitran_iso=1)
lines = load_lines(port; ν_min=2000, ν_max=2200)
pf    = partition_function(port, 5, 1)    # TabulatedPF from the ExoMol .pf file
```

Pair the port with its own `partition_function` — for ExoMol this returns a `TabulatedPF` built from the ExoMol `.pf` file, which is the partition function used to derive the intensities. Computed ExoMol intensities agree with HITRAN to better than 0.05% for CO.

ExoMol `.trans` files can be enormous (often gigabytes). They are streamed and windowed during loading rather than read into memory all at once, so requesting a narrow `ν_min`/`ν_max` window keeps the footprint small. As with HITRAN, the download is cached in scratch with a `.meta.toml` recording the species, line-list name (`"Li2015"` above), and spectral window for full provenance.

```julia
model = LineByLineModel(lines, pf; profile=Voigt(), wing_cutoff=40.0)
grid  = collect(2000.0:0.01:2200.0)                     # cm⁻¹
σ     = compute_cross_section(model, grid, 1013.25, 296.0)   # cm²/molecule
```

### ExoMol ↔ HITRAN cross-check

Because the two providers store fundamentally different quantities — ExoMol gives Einstein-A coefficients and a partition function, HITRAN gives precomputed intensities — agreement between them is a strong end-to-end test of the intensity machinery. The figure below overlays the CO cross-section from ExoMol's *Li2015* line list (intensities derived from Einstein-A on the fly) and from HITRAN over a slice of the 4.7 µm fundamental band, at the same pressure and temperature. Line positions coincide and matched-line intensities agree to **0.045 %**; the residual difference in the wings is broadening-limited (the two providers carry slightly different air-broadening parameters), not an intensity error.

<iframe title="ExoMol vs HITRAN CO" src="assets/plots/exomol_co_xsec.html" loading="lazy" style="width:100%;height:520px;border:1px solid var(--vp-c-divider);border-radius:8px;"></iframe>

The recipe that produces this comparison is in [Examples](@ref) (§2); the figure itself is regenerated by `benchmark/gen_exomol_fig.jl`.

## 5. Tabulated cross-sections (HITRAN `.xsc`)

Many heavy or strongly-mixing molecules — CFCs, HFCs, SF₆, and assorted VOCs — have no usable line list. HITRAN distributes them instead as **tabulated absorption cross-sections** (`.xsc`): laboratory-measured σ(ν) panels at fixed temperature and pressure. These are loaded through a different model type, `TabulatedCrossSection`, but expose the same `compute_cross_section` call as the line-by-line path.

```julia
using AtmosphericAbsorption

# Download a HITRAN .xsc by filename (as listed on hitran.org), or point at a local file
path  = fetch_hitran_xsc("SF6_..._N2.xsc")        # or: path = "/data/SF6.xsc"
model = load_xsc(path)                             # -> TabulatedCrossSection

grid = collect(940.0:0.1:960.0)                    # cm⁻¹
σ    = compute_cross_section(model, grid, 1013.25, 296.0)   # p[hPa], T[K]
```

A `.xsc` file holds one or more (T, p) panels; `compute_cross_section` interpolates σ onto your grid — linearly in wavenumber, and linearly in temperature at the nearest tabulated pressure, clamping to the measured range. Wavenumbers outside every panel return zero. `read_xsc(path)` returns the raw `XscBand`s if you want the panels directly.

## Provenance and reproducibility

Every cached artifact carries a `.meta.toml` sidecar recording its source, molecule/isotopologue, spectral window, and edition or line-list version. This lets you confirm — months later, or on another machine — exactly which spectroscopic data underlies a cross-section, and is the recommended way to track data versions in reproducible workflows. Query it at runtime with `source_metadata(port, mol, iso)`.
