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

Early development, built in phases:

- **Phase 0 ✓** — line-shape math + core abstractions: `Architectures`, `Constants`
  (TOML-defined), `LineShapes` (complex probability function + Doppler/Lorentz/Voigt),
  `PartitionFunctions`.
- **Phase 1 ✓** — HITRAN Voigt parity + GPU: columnar `LineDatabase`, the Port
  interface, the `HitranPort` (.par parser + TIPS-2017), the unified KernelAbstractions
  compute core, and CUDA/Metal extensions. Validated to <5×10⁻³ of HAPI on real CO2/H2O;
  ~10× faster than HAPI on CPU and ~1900× on an A100 GPU (see [benchmark/](benchmark/)).
- Phase 2 — advanced shapes (HT/SDV); Phase 3 — ExoMol; Phase 4 — line mixing +
  continuum; Phase 5 — vSmartMOM integration.

## HITRAN API key (for non-Voigt parameters)

Standard HITRAN line data and ExoMol downloads need **no credentials**. Fetching HITRAN
**non-Voigt** parameters (Hartmann-Tran / speed-dependent / line-mixing) uses the
authenticated HITRANonline API and requires your own key from your profile at
<https://hitran.org>. Supply it in-memory (never stored on disk or in this repo):

```julia
activate_hitran!("your-key")          # or set the HITRAN_API_KEY environment variable
```

## Citing HITRAN

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

The full author lists are in [docs/src/citing.md](docs/src/citing.md).

## Develop

```julia
julia --project=. -e 'using Pkg; Pkg.test()'
```
