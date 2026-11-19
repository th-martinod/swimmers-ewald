"""Parameter sweeps: python -m spinners.sweep configs/hex.toml"""
import argparse
import csv
import os
import time
try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10
    import tomli as tomllib
from multiprocessing import Pool

import numpy as np

from .backends import make_field
from .pair import observables

FIELDS = ["geometry", "lattice", "R", "h", "Q_over_h", "theta", "rate_conf", "omega_qdot", "A_conf",
          "A_qdot", "Vconf_par", "U_perp", "U_par", "U_perp_isolated", "err_conf", "err_qdot"]


def theta_grid(n, tmax):
    th = np.linspace(-tmax, tmax, n)
    th[np.argmin(np.abs(th))] = 0.0
    return th


def run_one(task):
    geometry, lattice, R, h, thetas, L_over_h = task
    field, kern = make_field(geometry, lattice, h, L_over_h=L_over_h)
    rows = []
    for th in thetas:
        o = observables(field, float(th), 1.0, R, kern)
        rows.append({"geometry": geometry, "lattice": lattice, "R": R, "h": h, "Q_over_h": 1.0 / h,
                     "theta": float(th), **o})
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("config")
    ap.add_argument("--out", default="results/data")
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    a = ap.parse_args(argv)
    with open(a.config, "rb") as fh:
        cfg = tomllib.load(fh)["sweep"]
    thetas = theta_grid(cfg.get("ntheta", 101), cfg.get("theta_max", np.pi / 4))
    L = tuple(cfg.get("L_over_h", (40, 60, 90, 120, 160)))
    tasks = [(g, cfg["lattice"], R, float(h), thetas, L)
             for g in cfg["geometries"] for R in cfg["R"] for h in cfg["h"]]
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    with Pool(a.workers) as pool:
        results = pool.map(run_one, tasks)
    for g in cfg["geometries"]:
        for R in cfg["R"]:
            fn = os.path.join(a.out, f"{cfg['lattice']}_{g}_R{R:g}.csv")
            with open(fn, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=FIELDS)
                w.writeheader()
                for res in results:
                    for r in res:
                        if r["geometry"] == g and r["R"] == R:
                            w.writerow(r)
            print("wrote", fn)
    print(f"done in {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
