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

## Regenerating the reference data (goldens + docs figures)

The HAPI golden spectra and the precomputed docs-figure data are **not committed** — they
live in the lazy `reference_data` Pkg artifact (a GitHub-release tarball pulled on demand by
the test suite and the docs build; see [`Artifacts.toml`](../Artifacts.toml)). The generator
scripts here still read/write the canonical staging layout — `test/golden/` (`*.par`, `*.txt`,
`*.data`) and `docs/src/assets/*.txt` — so regenerating and republishing is:

**1. Materialize the current artifact into the staging dirs** (so scripts that *consume* a
golden `.par` — `gen_exomol_fig.jl`, the Python timing harnesses — find their inputs):

```bash
julia --project=. -e 'using Pkg.Artifacts; \
  d = ensure_artifact_installed("reference_data", "Artifacts.toml"); \
  mkpath("test/golden"); for f in readdir(joinpath(d,"golden")); cp(joinpath(d,"golden",f), joinpath("test/golden",f); force=true); end; \
  for f in readdir(joinpath(d,"figures")); cp(joinpath(d,"figures",f), joinpath("docs/src/assets",f); force=true); end'
```

**2. Regenerate** whichever files changed (all need network / a HITRAN key):
- HAPI goldens → `python3 benchmark/hapi_reference/generate_golden.py` and `generate_pcqsdhc_golden.py` (write `test/golden/`)
- docs figure data → `source ~/.bashrc; julia --project=. benchmark/gen_physics_figs.jl` and `julia --project=. benchmark/gen_exomol_fig.jl` (write `docs/src/assets/*.txt`)

**3. Repack, upload, rebind.** Build a fresh tree from the two staging dirs, archive it, upload
the tarball as a **new** `reference-data-vN` GitHub release asset, and point `Artifacts.toml` at it:

```julia
using Pkg.Artifacts
hash = create_artifact() do dir
    mkpath(joinpath(dir, "golden")); mkpath(joinpath(dir, "figures"))
    for f in readdir("test/golden"); cp(joinpath("test/golden", f), joinpath(dir, "golden", f)); end
    for f in ("co2_linemix_band.txt", "exomol_co_xsec.txt", "h2o_voigt_vs_ht.txt")
        cp(joinpath("docs/src/assets", f), joinpath(dir, "figures", f))
    end
end
sha = archive_artifact(hash, "reference_data.tar.gz")          # → upload this to the new release
url = "https://github.com/RemoteSensingTools/AtmosphericAbsorption.jl/releases/download/reference-data-vN/reference_data.tar.gz"
bind_artifact!("Artifacts.toml", "reference_data", hash; download_info = [(url, sha)], lazy = true, force = true)
```

A new release tag (don't overwrite an old asset) keeps existing checkouts reproducible.

**Reading the table.** On the **CPU**, AtmosphericAbsorption is *comparable* to hapi2's
numba kernels (both are JIT-compiled to native code) — hapi2 edges ahead at tiny cutoffs
where Julia's per-line setup dominates, Julia pulls ahead as the window grows. The decisive
difference is the **GPU**: the *same* unified kernel runs unchanged on CUDA and is **55–460×
faster than hapi2** across every cutoff (and hapi2/HAPI have no GPU path at all). Float32
roughly halves the GPU time again. All engines compute the same physics — validated against
HAPI to ≤5×10⁻³.
