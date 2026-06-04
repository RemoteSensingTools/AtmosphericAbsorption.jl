# Benchmarks

Speed of AtmosphericAbsorption.jl vs the HITRAN reference implementations on a realistic
workload: **4000 Voigt lines over a 400 cm⁻¹ band** (6000–6400 cm⁻¹, 40001 grid points),
swept over the wing **cutoff** — the dominant performance lever. hapi2 is run *windowed*
(the profile evaluated only within ±cutoff of each line, exactly as its
`absorptionCoefficient` does internally); a "no cutoff" run computes every line over the
whole grid.

## Run

```bash
# hapi2 numba (windowed) — needs numba + the standalone numba kernels in /tmp/hnb
PYTHONPATH=/tmp:/path/to/hapi python3 benchmark/hapi_reference/bench_hapi2.py
# Julia (CPU; GPU if CUDA is on the load path)
julia --project=benchmark benchmark/benchmarks.jl
```

## Results (NVIDIA A100, single CPU core), time in ms

| wing cutoff [cm⁻¹] | hapi2 (numba, windowed) | AtmosphericAbsorption CPU (f64) | AtmosphericAbsorption GPU (f64) | GPU speedup vs hapi2 |
|---|---|---|---|---|
| 2.5  | 114  | 177  | 2.1  | **55×** |
| 5    | 173  | 235  | 2.2  | **80×** |
| 10   | 286  | 349  | 2.4  | **120×** |
| 25   | 623  | 683  | 3.0  | **209×** |
| full (no cutoff) | 4691 | 4122 | 10.1 | **464×** |

vs the legacy **HAPI** (pure numpy, no JIT) the CPU path is ~3–8× faster (HAPI ≈ 4800 ms on
the Phase-1 CO₂ case where the Julia CPU path took ≈580 ms).

**Reading the table.** On the **CPU**, AtmosphericAbsorption is *comparable* to hapi2's
numba kernels (both are JIT-compiled to native code) — hapi2 edges ahead at tiny cutoffs
where Julia's per-line setup dominates, Julia pulls ahead as the window grows. The decisive
difference is the **GPU**: the *same* unified kernel runs unchanged on CUDA and is **55–460×
faster than hapi2** across every cutoff (and hapi2/HAPI have no GPU path at all). Float32
roughly halves the GPU time again. All engines compute the same physics — validated against
HAPI to ≤5×10⁻³.
