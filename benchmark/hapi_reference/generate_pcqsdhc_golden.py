#!/usr/bin/env python3
"""Golden values for the pCqSDHC line shape (HT/SDV) from HAPI's pcqsdhc, to validate
AtmosphericAbsorption.jl's port. Writes test/golden/pcqsdhc_<case>.txt with columns
ν, Re(LS), Im(LS). Run: PYTHONPATH=/home/cfranken/code/gitHub/hapi python3 generate_pcqsdhc_golden.py
"""
import os
import numpy as np
import hapi

HERE = os.path.dirname(os.path.abspath(__file__))
GOLDEN = os.path.join(HERE, "..", "..", "test", "golden")
os.makedirs(GOLDEN, exist_ok=True)

NU0, GAMD = 1000.0, 0.02
grid = np.arange(999.5, 1000.5 + 1e-9, 0.002)

# name: (Gam0, Gam2, Shift0, Shift2, anuVC, eta)
CASES = {
    "voigt":  (0.03, 0.0,   -0.01, 0.0,   0.0,  0.0),   # must reduce to Voigt
    "sdv":    (0.03, 0.006, -0.01, 0.002, 0.0,  0.0),   # speed-dependent Voigt
    "ht":     (0.03, 0.006, -0.01, 0.002, 0.012, 0.3),  # full Hartmann-Tran
}

for name, (Gam0, Gam2, Shift0, Shift2, anuVC, eta) in CASES.items():
    re, im = hapi.pcqsdhc(NU0, GAMD, Gam0, Gam2, Shift0, Shift2, anuVC, eta, grid)
    np.savetxt(os.path.join(GOLDEN, "pcqsdhc_%s.txt" % name),
               np.column_stack([grid, re, im]),
               header="nu Re(LS) Im(LS)  Gam0=%g Gam2=%g Shift0=%g Shift2=%g anuVC=%g eta=%g"
               % (Gam0, Gam2, Shift0, Shift2, anuVC, eta))
    print("wrote pcqsdhc_%s: max|Re|=%.4g" % (name, np.max(np.abs(re))))
