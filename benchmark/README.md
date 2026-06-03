# Benchmarks

Speed of AtmosphericAbsorption.jl vs the HAPI reference, on the same real HITRAN case
(CO2 6300–6400 cm⁻¹, 4971 lines × 10001 grid points, Voigt, wing 40 cm⁻¹).

## Run

```bash
# Julia (CPU; also GPU if CUDA is on the load path)
julia --project=benchmark -e 'using Pkg; Pkg.instantiate()'
julia --project=benchmark benchmark/benchmarks.jl

# HAPI reference
PYTHONPATH=/path/to/hapi python3 benchmark/hapi_reference/time_hapi.py
```

## Results (NVIDIA A100, single CPU core)

| Engine | Precision | Time | Throughput | vs HAPI |
|--------|-----------|------|------------|---------|
| HAPI (numpy)         | f64 | 4799 ms | 0.010 Gline·pt/s | 1× |
| AtmosphericAbsorption CPU | f64 | 582 ms | 0.09 Gline·pt/s | **8.3×** |
| AtmosphericAbsorption CPU | f32 | 471 ms | 0.11 Gline·pt/s | **10×** |
| AtmosphericAbsorption GPU | f64 | 5.4 ms | 9.1 Gline·pt/s | **882×** |
| AtmosphericAbsorption GPU | f32 | 2.5 ms | 19.9 Gline·pt/s | **1925×** |

Same physics, validated to within 5×10⁻³ of HAPI (see `test/test_hitran_golden.jl`).
The CPU path is single-threaded here; the GPU path is the same unified kernel,
unchanged, dispatched to CUDA.
