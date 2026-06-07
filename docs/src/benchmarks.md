# Benchmarks

Performance numbers reported honestly: where the CPU path is merely competitive, and where the GPU path is decisive. The point of this page is to let you predict runtime for your own workload, not to cherry-pick a flattering single number.

## Methodology

All timings below evaluate the **same cross-section** with `compute_cross_section`:

- **4000 Voigt lines** of CO₂ over a **400 cm⁻¹ band** (6000–6400 cm⁻¹),
- on a uniform grid of **40001 points** (0.01 cm⁻¹ spacing),
- at a single pressure/temperature point,
- swept over the **`wing_cutoff`** parameter — the half-width (in cm⁻¹) of the window around each line center within which the line profile is actually evaluated.

The hardware: one **NVIDIA A100** for the GPU rows, a **single CPU core** for the CPU and `hapi2` rows.

The reference is **`hapi2`** (the HITRAN Python API, numba-JIT backend), run **windowed** — its `absorptionCoefficient` evaluates each line's profile only within `±wing_cutoff` of the line center, exactly as the package does. The `full` row removes the window entirely: every one of the 4000 lines is evaluated over all 40001 grid points. That row is the apples-to-apples "no shortcuts" comparison, and the one where vectorized hardware matters most.

The package call being timed is the standard pipeline:

```julia
using AtmosphericAbsorption

port = HitranPort(; molecule="CO2", numin=6000, numax=6400, edition="HITRAN2020")
lines = load_lines(port; mol=2, iso=1, ν_min=6000, ν_max=6400)
pf    = partition_function(port, 2, 1)

grid = collect(range(6000.0, 6400.0; length=40001))   # cm⁻¹

# CPU f64
model_cpu = LineByLineModel(lines, pf; profile=Voigt(), wing_cutoff=10.0,
                            architecture=CPU())
σ_cpu = compute_cross_section(model_cpu, grid, 1013.25, 296.0)   # p [hPa], T [K]

# GPU f64 — identical kernel, different architecture
model_gpu = LineByLineModel(lines, pf; profile=Voigt(), wing_cutoff=10.0,
                            architecture=GPU())
σ_gpu = Array(compute_cross_section(model_gpu, grid, 1013.25, 296.0))
```

Only the `architecture` and `wing_cutoff` keywords change between rows.

## Results

| wing cutoff [cm⁻¹] | hapi2 (numba, windowed) | CPU f64 | GPU f64 | GPU vs hapi2 |
|---|---|---|---|---|
| 2.5 | 114 ms | 177 ms | 2.1 ms | 55× |
| 5 | 173 ms | 235 ms | 2.2 ms | 80× |
| 10 | 286 ms | 349 ms | 2.4 ms | 120× |
| 25 | 623 ms | 683 ms | 3.0 ms | 209× |
| full | 4691 ms | 4122 ms | 10.1 ms | 464× |

```@raw html
<iframe title="Benchmark" src="assets/plots/benchmark.html" loading="lazy" style="width:100%;height:520px;border:1px solid var(--vp-c-divider);border-radius:8px;"></iframe>
```

## How to read this

**On the CPU, the package is comparable to `hapi2`.** Both are JIT-compiled vectorized kernels, so neither has a structural advantage. At tiny wing cutoffs (2.5–5 cm⁻¹) `hapi2` edges ahead — there the per-line setup cost dominates and the absolute work is small. As the window grows, the package pulls ahead: at the `full` cutoff the CPU path (4122 ms) is faster than windowed `hapi2` (4691 ms), because the package's kernel amortizes setup better over the larger workload. If you have been using the **legacy HAPI** (pure NumPy, no numba), the CPU path here is roughly **3–8× faster** across this sweep.

So on a single core, treat the package as *a peer of `hapi2`*, not a leap — the real reason to adopt it is portability and the device behind it.

**On the GPU, the comparison changes character.** The exact same kernel — no rewrite, only `architecture=GPU()` — runs on CUDA and lands **55–460× faster than `hapi2`**, which has no GPU path at all. The GPU advantage widens as the window grows, because more lines × more points is exactly the regime that saturates the A100's throughput while leaving the CPU and `hapi2` to grind serially. At the `full` cutoff the GPU is still only **10.1 ms** against `hapi2`'s **4691 ms**.

**Float32 roughly halves GPU time again.** Loading the line list with `FT=Float32` produces a Float32-clean pipeline (no silent Float64 promotion), which on the A100 cuts memory traffic and runtime by about half with no measurable accuracy loss versus the Float64 path:

```julia
lines32 = load_lines(port; mol=2, iso=1, ν_min=6000, ν_max=6400, FT=Float32)
model32 = LineByLineModel(lines32, pf; profile=Voigt(), wing_cutoff=10.0,
                          architecture=GPU())
σ32 = Array(compute_cross_section(model32, grid, 1013.25, 296.0))
```

## Caveats

- These are **single-point** cross-sections (one pressure, one temperature). The GPU's per-call fixed overhead (kernel launch, host↔device transfer) is amortized far better when you batch many spectra; a single tiny grid is the *least* favorable case for the GPU and it still wins by 55×.
- `Array(...)` is included in the GPU timings — the result is copied back to the host. If you keep the result on the device for downstream work, subtract that transfer.
- CPU and `hapi2` rows are single-core. Both can be threaded; the comparison here isolates the per-core kernel, not a fully parallelized deployment.
- Numbers are for this specific band and line count on this specific hardware. Use the **shape** of the sweep — flat GPU, cutoff-linear CPU/`hapi2` — to extrapolate to your own grids rather than the absolute milliseconds.
