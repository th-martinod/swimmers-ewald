"""Parameter sweeps: spinners-sweep configs/paper.toml [--out DIR] [--workers N]"""
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

FIELDS = ["geometry", "lattice", "R", "h", "Q_over_h", "theta", "T", "A_perp", "B_perp", "err_A", "err_B",
          "B_img", "gradp_img"]


def theta_grid(n, tmax):
    th = np.linspace(-tmax, tmax, n)
    th[np.argmin(np.abs(th))] = 0.0
    return th


def theta_union(ranges):
    """Union of symmetric grids [[tmax, n], ...], rounded to avoid near-duplicates, with theta = 0 exact."""
    th = np.concatenate([theta_grid(int(n), float(t)) for t, n in ranges])
    th = np.unique(np.round(th, 12))
    th[np.argmin(np.abs(th))] = 0.0
    return th


def run_one(task):
    geometry, lattice, R, h, thetas, L_over_h = task
    t0 = time.time()
    field = make_field(geometry, lattice, h, L_over_h=L_over_h)
    rows = []
    for th in thetas:
        o = observables(field, float(th), 1.0, R)
        rows.append({"geometry": geometry, "lattice": lattice, "R": R, "h": h, "Q_over_h": 1.0 / h,
                     "theta": float(th), **{k: o[k] for k in FIELDS[6:]}})
    return rows, time.time() - t0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("config")
    ap.add_argument("--out", default="results/data")
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    a = ap.parse_args(argv)
    with open(a.config, "rb") as fh:
        cfg = tomllib.load(fh)["sweep"]
    if "theta_ranges" in cfg:
        thetas = theta_union(cfg["theta_ranges"])
    else:
        thetas = theta_grid(cfg.get("ntheta", 101), cfg.get("theta_max", np.pi / 4))
    L = tuple(cfg.get("L_over_h", (40, 60, 90, 120, 160)))
    lattices = cfg.get("lattices", [cfg.get("lattice")])
    tasks = [(g, lat, R, float(h), thetas, L)
             for lat in lattices for g in cfg["geometries"] for R in cfg["R"] for h in cfg["h"]]
    # heaviest tasks first for better load balance
    tasks.sort(key=lambda t: (t[0] == "bulk3d", -t[3]))
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    print(f"{len(tasks)} tasks ({len(thetas)} theta values each), {a.workers} workers -> {a.out}", flush=True)
    results = []
    with Pool(a.workers) as pool:
        for i, (rows, dt) in enumerate(pool.imap_unordered(run_one, tasks, chunksize=1), 1):
            r = rows[0]
            el = time.time() - t0
            print(f"[{i:4d}/{len(tasks)}] {r['lattice']:3s} {r['geometry']:9s} R={r['R']:<4g} h={r['h']:<4g} "
                  f"{dt:6.1f} s | elapsed {el / 60:6.1f} min | ETA {el / i * (len(tasks) - i) / 60:6.1f} min",
                  flush=True)
            results.append(rows)
    for lat in lattices:
        for g in cfg["geometries"]:
            for R in cfg["R"]:
                fn = os.path.join(a.out, f"{lat}_{g}_R{R:g}.csv")
                rows = sorted((r for res in results for r in res
                               if r["lattice"] == lat and r["geometry"] == g and r["R"] == R),
                              key=lambda r: (r["h"], r["theta"]))
                with open(fn, "w", newline="") as fh:
                    w = csv.DictWriter(fh, fieldnames=FIELDS)
                    w.writeheader()
                    w.writerows(rows)
                print("wrote", fn, flush=True)
    print(f"done in {time.time() - t0:.1f} s", flush=True)


if __name__ == "__main__":
    main()
