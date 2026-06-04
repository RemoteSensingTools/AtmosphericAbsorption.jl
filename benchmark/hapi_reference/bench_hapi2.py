#!/usr/bin/env python3
"""Benchmark hapi2's numba Voigt profile (PROFILE_SDVOIGT with Γ2=0), WINDOWED by a wing
cutoff exactly as HAPI/hapi2's absorptionCoefficient does internally, over a range of
cutoffs on a large spectral band. Prints CSV: engine,cutoff_cm1,time_ms,Glinept_per_s.

Run: PYTHONPATH=/tmp:/home/cfranken/code/gitHub/hapi python3 bench_hapi2.py
(the numba kernels are copied to /tmp/hnb; see the docs build for setup).
"""
import sys, time
import numpy as np
sys.path.insert(0, "/tmp")
from hnb import sdv

NU_MIN, NU_MAX, STEP = 6000.0, 6400.0, 0.01
N = 4000                          # lines, uniformly spread over the band
GAMD, GAM0 = 0.008, 0.04          # Doppler / Lorentz HWHM [cm-1]
CUTOFFS = [2.5, 5.0, 10.0, 25.0, None]   # None = full grid (no cutoff)

grid = np.arange(NU_MIN, NU_MAX + 0.5 * STEP, STEP)
npts = len(grid)
nu0 = np.linspace(NU_MIN, NU_MAX, N)
Sw = np.full(N, 1e-23)
work = N * npts


def run(cutoff):
    x = np.zeros(npts)
    for j in range(N):
        if cutoff is None:
            x += Sw[j] * sdv.PROFILE_SDVOIGT(nu0[j], GAMD, GAM0, 0., 0., 0., grid)[0]
        else:
            lo = np.searchsorted(grid, nu0[j] - cutoff)
            hi = np.searchsorted(grid, nu0[j] + cutoff)
            if hi > lo:
                x[lo:hi] += Sw[j] * sdv.PROFILE_SDVOIGT(nu0[j], GAMD, GAM0, 0., 0., 0., grid[lo:hi])[0]
    return x


print("engine,cutoff_cm1,time_ms,Glinept_per_s")
for c in CUTOFFS:
    run(c); run(c)                 # numba JIT warmup
    ts = []
    for _ in range(3):
        s = time.perf_counter(); run(c); ts.append(time.perf_counter() - s)
    t = sorted(ts)[1]
    label = "full" if c is None else f"{c:g}"
    print("hapi2_numba,%s,%.1f,%.4f" % (label, 1e3 * t, work / t / 1e9))
