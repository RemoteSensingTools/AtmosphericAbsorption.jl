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

## Develop

```julia
julia --project=. -e 'using Pkg; Pkg.test()'
```
