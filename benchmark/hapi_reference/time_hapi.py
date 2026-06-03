#!/usr/bin/env python3
"""Time HAPI's absorptionCoefficient_Voigt on the same CO2 case the Julia benchmark
uses, so the two throughputs are comparable. Times only the cross-section call
(line list already cached), median of a few runs.

Run:  PYTHONPATH=/home/cfranken/code/gitHub/hapi python3 time_hapi.py
"""
import os
import json
import time
import tempfile
import numpy as np
import hapi

HERE = os.path.dirname(os.path.abspath(__file__))
PAR = os.path.join(HERE, "..", "..", "test", "golden", "co2_6300_6400_p500_T250.par")
NU_MIN, NU_MAX, STEP, P, T, WING = 6300.0, 6400.0, 0.01, 500.0, 250.0, 40.0
P_REF = 1013.25


def main():
    db = tempfile.mkdtemp(prefix="hapibench_")
    hapi.db_begin(db)
    with open(PAR) as f:
        rows = [ln.rstrip("\n") for ln in f if len(ln) >= 160]
    with open(os.path.join(db, "co2.data"), "w") as f:
        f.write("\n".join(rows) + "\n")
    header = dict(hapi.HITRAN_DEFAULT_HEADER)
    header["table_name"] = "co2"
    header["number_of_rows"] = len(rows)
    with open(os.path.join(db, "co2.header"), "w") as f:
        json.dump(header, f)
    hapi.storage2cache("co2")

    grid = np.arange(NU_MIN, NU_MAX + 0.5 * STEP, STEP)
    npts = len(grid)
    work = len(rows) * npts

    def run():
        hapi.absorptionCoefficient_Voigt(
            Components=[(2, 1)], SourceTables=["co2"],
            Environment={"p": P / P_REF, "T": T}, OmegaGrid=grid,
            HITRAN_units=True, LineShift=True,
            WavenumberWing=WING, WavenumberWingHW=0.0, Diluent={"air": 1.0})

    run()  # warmup
    ts = []
    for _ in range(5):
        t0 = time.perf_counter()
        run()
        ts.append(time.perf_counter() - t0)
    med = sorted(ts)[len(ts) // 2]
    print("HAPI CPU: %.1f ms  (%.4f Gline*pt/s)  [%d lines x %d points]"
          % (1e3 * med, work / med / 1e9, len(rows), npts))


if __name__ == "__main__":
    main()
