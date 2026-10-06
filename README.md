# AtmosphericAbsorption.jl

Clean-slate, GPU-accelerated molecular absorption cross-sections for atmospheric
radiative transfer — a standalone successor to vSmartMOM's `Absorption` module.

**Goals**
- Pluggable line-list sources behind one interface: HITRAN, ExoMol, … (the "Port" hierarchy).
- Advanced line shapes: Voigt → speed-dependent Voigt → Hartmann-Tran, plus line mixing.
- Water-vapor continuum (MT_CKD) and CIA.
- One KernelAbstractions compute core running on CPU, CUDA, and Metal.
- Type-stable and correct in both Float32 and Float64 to machine precision.
- Validated against [hapi2](https://github.com/hitranonline/hapi2) as a golden-file benchmark.

## Status

Implemented:

- **Core** — `Architectures`, TOML-defined `Constants`, columnar `LineDatabase`, the Port
  interface, and one KernelAbstractions compute core with CUDA/Metal extensions.
- **HITRAN** — `.par` parser and direct download (`HitranPort`, `load_hitran`), authenticated
  non-Voigt parameters (`load_hitran_nonvoigt`), and tabulated `.xsc` cross-sections. Voigt
  validated to <5×10⁻³ of HAPI on real CO2/H2O; ~10× faster than HAPI on CPU and ~1900× on an
  A100 GPU (see [benchmark/](benchmark/)).
- **ExoMol** — `ExoMolPort`, with line strengths derived from Einstein-A coefficients and the
  ExoMol partition function.
- **Line shapes** — Doppler, Lorentz, Voigt, speed-dependent Voigt, Rautian and Hartmann-Tran
  (HAPI's `pcqsdhc`, matched to ~1e-6), plus first-order (Rosenkranz) line mixing.
- **Partition functions** — TIPS-2021 (default), TIPS-2017, and tabulated ExoMol `.pf`.
- **Continuum** — MT_CKD water-vapor continuum and HITRAN CIA.
- **Lookup tables** — interpolated cross-section models, AER ABSCO tables (including GPU-native
  OCO-2 LUTs), and batched profile queries (`compute_cross_section_profile`).
- **vSmartMOM** — vSmartMOM.jl builds its absorption models with this package; its legacy
  `Absorption` module has not been retired yet.

See the [documentation](docs/src/index.md) for details.

## HITRAN API key (for non-Voigt parameters)

Standard HITRAN line data and ExoMol downloads need **no credentials**. Fetching HITRAN
**non-Voigt** parameters (Hartmann-Tran / speed-dependent / line-mixing) uses the
authenticated HITRANonline API and requires your own key from your profile at
<https://hitran.org>. Supply it in-memory (never stored on disk or in this repo):

```julia
activate_hitran!("your-key")          # or set the HITRAN_API_KEY environment variable
```

## Citing HITRAN and ExoMol

Line data come from the [HITRAN database](https://hitran.org). If you publish results computed
with HITRAN data through this package, **please cite the HITRAN edition you used**. Direct
downloads from hitran.org always return the current edition (HITRAN2024 at the time of writing).
The `edition` keyword is only a cache/provenance label. For local `.par` files, cite the edition
the file came from.

- **HITRAN2024:** Gordon IE, Rothman LS, Hargreaves RJ, et al. The HITRAN2024 molecular
  spectroscopic database. *J Quant Spectrosc Radiat Transf* 2026;353:109807.
  <https://doi.org/10.1016/j.jqsrt.2026.109807>
- **HITRAN2020:** Gordon IE, Rothman LS, Hargreaves RJ, et al. The HITRAN2020 molecular
  spectroscopic database. *J Quant Spectrosc Radiat Transf* 2022;277:107949.
  <https://doi.org/10.1016/j.jqsrt.2021.107949>
- **HITRAN2016:** Gordon IE, Rothman LS, Hill C, et al. The HITRAN2016 molecular spectroscopic
  database. *J Quant Spectrosc Radiat Transf* 2017;203:3–69.
  <https://doi.org/10.1016/j.jqsrt.2017.06.038>

If you use the speed-dependent / Hartmann–Tran profiles (based on HAPI's `pcqsdhc`) or
non-Voigt parameters from the HITRANonline API, please also cite **HAPI**:

- Kochanov RV, Gordon IE, Rothman LS, Wcisło P, Hill C, Wilzewski JS. HITRAN Application
  Programming Interface (HAPI): A comprehensive approach to working with spectroscopic data.
  *J Quant Spectrosc Radiat Transf* 2016;177:15–30. <https://doi.org/10.1016/j.jqsrt.2016.03.005>

If you use **ExoMol** line lists (`ExoMolPort`), please cite the ExoMol database and the paper
for each line list you used. That paper is listed on the line list's page at
<https://www.exomol.com>. ExoMol data are released under CC BY-SA 4.0.

- **ExoMol 2024 release:** Tennyson J, Yurchenko SN, Zhang J, et al. The 2024 release of the
  ExoMol database: Molecular line lists for exoplanet and other hot atmospheres.
  *J Quant Spectrosc Radiat Transf* 2024;326:109083. <https://doi.org/10.1016/j.jqsrt.2024.109083>
- **ExoMol database and data format:** Tennyson J, Yurchenko SN, Al-Refaie AF, et al. The ExoMol
  database: Molecular line lists for exoplanet and other hot atmospheres. *J Mol Spectrosc*
  2016;327:73–94. <https://doi.org/10.1016/j.jms.2016.05.002>
- **Line list, e.g. CO `Li2015`:** Li G, Gordon IE, Rothman LS, et al. Rovibrational line lists
  for nine isotopologues of the CO molecule in the X ¹Σ⁺ ground electronic state.
  *Astrophys J Suppl Ser* 2015;216:15. <https://doi.org/10.1088/0067-0049/216/1/15>

The full author lists are in [docs/src/citing.md](docs/src/citing.md).

## Develop

```julia
julia --project=. -e 'using Pkg; Pkg.test()'
```
