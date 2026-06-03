#!/usr/bin/env python3
"""Generate golden absorption cross-sections with HAPI (the HITRAN API) for
validating AtmosphericAbsorption.jl.

HAPI is the reference hapi2 wraps. For each case we feed HAPI the same HITRAN .par
lines (filtered to a window), the same wing cutoff, and HITRAN units (cm2/molecule),
then write (wavenumber, sigma) columns to test/golden/<name>.txt. The Julia tests load
the wavenumber grid from the golden file and compare their cross-section to sigma.

Run:  PYTHONPATH=/home/cfranken/code/gitHub/hapi python3 generate_golden.py
"""
import os
import json
import tempfile
import numpy as np
import hapi

HERE = os.path.dirname(os.path.abspath(__file__))
GOLDEN = os.path.join(HERE, "..", "..", "test", "golden")
os.makedirs(GOLDEN, exist_ok=True)

ART = os.path.expanduser("~/.julia/artifacts")
PAR = {
    2: os.path.join(ART, "bd1208b6552fe1d2c7419043680ba0b7486062a6",
                    "hitran_molec_id_2_CO2.par"),
    1: os.path.join(ART, "738b53c3ad6f795ebae1f9af6642f7b04efba4d6",
                    "hitran_molec_id_1_H2O.par"),
}

# name, molecule, isotopologue, nu_min, nu_max, step, p[hPa], T[K], wing[cm-1]
CASES = [
    ("co2_6300_6400_p500_T250", 2, 1, 6300.0, 6400.0, 0.01, 500.0, 250.0, 40.0),
    ("co2_6300_6400_p1013_T296", 2, 1, 6300.0, 6400.0, 0.01, 1013.25, 296.0, 40.0),
    ("h2o_7000_7100_p800_T280", 1, 1, 7000.0, 7100.0, 0.01, 800.0, 280.0, 40.0),
]

P_REF = 1013.25  # hPa


def isoid(ch):
    # HITRAN local_iso_id: 1-9 digits, '0' -> 10, then 'A' -> 11, 'B' -> 12, ...
    if ch.isdigit():
        return 10 if ch == "0" else int(ch)
    return 11 + (ord(ch.upper()) - ord("A"))


def write_table(db, name, par, mol, iso, nu_min, nu_max):
    """Filter .par to (mol, iso, window) and write HAPI .data + .header."""
    rows = []
    with open(par) as f:
        for ln in f:
            if len(ln) < 160:
                continue
            if int(ln[0:2]) != mol or isoid(ln[2]) != iso:
                continue
            nu = float(ln[3:15])
            if nu_min - 60.0 <= nu <= nu_max + 60.0:   # pad for the wing
                rows.append(ln.rstrip("\n"))
    with open(os.path.join(db, name + ".data"), "w") as f:
        f.write("\n".join(rows) + "\n")
    # Also ship the filtered .par subset so the Julia test is self-contained.
    with open(os.path.join(GOLDEN, name + ".par"), "w") as f:
        f.write("\n".join(rows) + "\n")
    header = dict(hapi.HITRAN_DEFAULT_HEADER)
    header["table_name"] = name
    header["number_of_rows"] = len(rows)
    with open(os.path.join(db, name + ".header"), "w") as f:
        json.dump(header, f)
    return len(rows)


def main():
    db = tempfile.mkdtemp(prefix="hapidb_")
    hapi.db_begin(db)
    for (name, mol, iso, nu_min, nu_max, step, p, T, wing) in CASES:
        n = write_table(db, name, PAR[mol], mol, iso, nu_min, nu_max)
        hapi.storage2cache(name)
        grid = np.arange(nu_min, nu_max + 0.5 * step, step)
        nu, coef = hapi.absorptionCoefficient_Voigt(
            Components=[(mol, iso)], SourceTables=[name],
            Environment={"p": p / P_REF, "T": T},
            OmegaGrid=grid, HITRAN_units=True, LineShift=True,
            WavenumberWing=wing, WavenumberWingHW=0.0,
            Diluent={"air": 1.0},
        )
        out = os.path.join(GOLDEN, name + ".txt")
        np.savetxt(out, np.column_stack([nu, coef]),
                   header="nu[cm-1] sigma[cm2/molecule]  (%d lines, p=%.2fhPa T=%.1fK wing=%.0f)"
                   % (n, p, T, wing))
        print("wrote %s: %d lines, %d points, max sigma=%.3e"
              % (name, n, len(nu), np.max(coef)))


if __name__ == "__main__":
    main()
