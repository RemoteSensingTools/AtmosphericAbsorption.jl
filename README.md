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

- **Phase 0 (current)** — line-shape math + core abstractions, no I/O:
  `Architectures`, `Constants`, `LineShapes` (complex probability function +
  area-normalized Doppler/Lorentz/Voigt), `PartitionFunctions`.
- Phase 1 — HITRAN Voigt parity + GPU; Phase 2 — advanced shapes;
  Phase 3 — ExoMol; Phase 4 — line mixing + continuum; Phase 5 — vSmartMOM integration.

## Develop

```julia
julia --project=. -e 'using Pkg; Pkg.test()'
```
