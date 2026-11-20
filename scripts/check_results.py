"""Checks of the sweep results (results/data) against the acceptance criteria of the review.

Writes results/logs/checks.log with one PASS / FAIL / WARN / INFO line per check and the number that decides it.
Run: python3 scripts/check_results.py
"""
import csv
import os
import sys

import numpy as np

ROOT = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from spinners.backends import make_field  # noqa: E402
from spinners.pair import observables  # noqa: E402

DATA = os.path.join(ROOT, "results", "data")
LOG = os.path.join(ROOT, "results", "logs", "checks.log")
LATS = ["scl", "bcl", "fcl", "hex"]
GEOS = ["bulk3d", "monolayer"]
RS = [0.5, 0.25]
NROWS = 15 * 161

REF = {  # (geometry, lattice): (B_perp, A_perp) at h = 5, theta = 0.3, R = 0.5. None means ~0
    # B_perp with the common reference ball r_f = sqrt(2) h / 4 (A_perp does not depend on it)
    ("bulk3d", "scl"): (-2.1017e-3, None), ("bulk3d", "bcl"): (-1.8882e-3, -9.044e-4),
    ("bulk3d", "fcl"): (-9.662e-4, -2.132e-3), ("bulk3d", "hex"): (-2.4892e-3, -6.199e-4),
    ("monolayer", "scl"): (-1.9893e-3, None), ("monolayer", "bcl"): (-1.7127e-3, -6.200e-4),
    ("monolayer", "fcl"): (-1.0391e-3, -1.703e-3), ("monolayer", "hex"): (-2.3750e-3, -6.017e-4),
}
R_F_OVER_H = np.sqrt(2) / 4
SIX_FOLD = {"monolayer": 0.927, "bulk3d": 0.955}

lines = []
counts = {"PASS": 0, "FAIL": 0, "WARN": 0, "INFO": 0}


def report(status, name, detail):
    counts[status] += 1
    line = f"{status:4s}  {name}: {detail}"
    lines.append(line)
    print(line, flush=True)


def check(ok, name, detail):
    report("PASS" if ok else "FAIL", name, detail)


def load(lat, geo, R):
    with open(os.path.join(DATA, f"{lat}_{geo}_R{R:g}.csv")) as fh:
        rows = list(csv.DictReader(fh))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0] if k not in ("geometry", "lattice")}


def by_h(c):
    for h in np.unique(c["h"]):
        yield h, np.isclose(c["h"], h)


def main():
    # 1. files
    data = {}
    for lat in LATS:
        for geo in GEOS:
            for R in RS:
                name = f"{lat}_{geo}_R{R:g}.csv"
                if not os.path.exists(os.path.join(DATA, name)):
                    check(False, f"1 files {name}", "missing")
                    continue
                c = load(lat, geo, R)
                n = len(c["theta"])
                bad = sum(int(np.sum(~np.isfinite(v))) for v in c.values())
                check(n == NROWS and bad == 0, f"1 files {name}", f"rows = {n} (expected {NROWS}), non-finite = {bad}")
                data[lat, geo, R] = c
    report("INFO", "1 files count", f"{len(data)} of 16 CSV present")

    # 2. theta = 0
    worst = max(np.max(np.abs(c["A_perp"][c["theta"] == 0.0])) for c in data.values())
    nzero = min(int(np.sum(c["theta"] == 0.0)) for c in data.values())
    check(worst == 0.0 and nzero == 15, "2 theta = 0", f"max|A_perp| = {worst:.3e} (min rows with theta = 0 per file: {nzero})")

    # 3. SCL
    for geo in GEOS:
        m = np.max(np.abs(data["scl", geo, 0.5]["A_perp"]))
        check(m < 1e-13, f"3 SCL {geo} R=0.5 exact cancellation", f"max|A_perp| = {m:.3e} (< 1e-13)")
        m = np.max(np.abs(data["scl", geo, 0.25]["A_perp"]))
        check(m > 1e-10, f"3 SCL {geo} R=0.25 non-zero", f"max|A_perp| = {m:.3e} (> 1e-10)")

    # 4. HEX
    for geo in GEOS:
        for R in RS:
            c = data["hex", geo, R]
            th, A, B = c["theta"], c["A_perp"], c["B_perp"]
            at = np.abs(th)
            neg = (at > 0.02) & (at < 0.50)
            pos = (at > 0.55) & (at < np.pi / 4)
            tag = f"4 HEX {geo} R={R:g}"
            check(np.all(A[neg] < 0), f"{tag} A<0 on 0.02<|th|<0.50",
                  f"max A_perp = {np.max(A[neg]):.3e} over {neg.sum()} rows")
            check(np.all(A[pos] > 0), f"{tag} A>0 on 0.55<|th|<pi/4",
                  f"min A_perp = {np.min(A[pos]):.3e} over {pos.sum()} rows")
            worst_par, spreads = 0.0, []
            for h, m in by_h(c):
                t, a, b = th[m], A[m], B[m]
                key = {round(x, 10): y for x, y in zip(t, a)}
                diffs = [abs(y - key[round(-x, 10)]) for x, y in zip(t, a) if round(-x, 10) in key]
                worst_par = max(worst_par, max(diffs) / np.max(np.abs(a)))
                spreads.append((h, (b.max() - b.min()) / abs(b.mean())))
            check(worst_par <= 1e-6, f"{tag} A even in theta",
                  f"max_h max|A(th)-A(-th)|/max|A| = {worst_par:.3e} (<= 1e-6)")
            if R == 0.5:
                # square of forces: its coupling with the six-fold lattice cancels at low order
                worst = max(sp for _, sp in spreads)
                check(worst < 1e-3, f"{tag} B independent of theta",
                      f"max_h (max-min)/|mean| of B_perp = {worst:.3e} (< 1e-3)")
            else:
                # rectangle of forces: real oscillation of period pi/3 whose relative size decays as (Q/h)^2
                scaled = np.array([sp * h**2 for h, sp in spreads])
                check(scaled.max() / scaled.min() < 1.2, f"{tag} theta dependence of B decays as (Q/h)^2",
                      f"(max-min)/|mean| from {spreads[0][1]:.2e} (h={spreads[0][0]:g}) to {spreads[-1][1]:.2e} "
                      f"(h={spreads[-1][0]:g}), spread*h^2 varies by a factor {scaled.max() / scaled.min():.3f} (< 1.2)")

    # 5. sign of the theta-averaged B_perp (uniform grid of 80 angles in [-pi, pi), as in the figures)
    grid = np.round(np.linspace(-np.pi, np.pi, 81)[:-1], 12)
    for (lat, geo, R), c in data.items():
        sel = np.isin(np.round(c["theta"], 12), grid)
        means = np.array([np.mean(c["B_perp"][sel & m]) for _, m in by_h(c)])
        status = "PASS" if np.all(means < 0) else "WARN"
        report(status, f"5 B_perp mean sign {lat} {geo} R={R:g}",
               f"max over h of mean B_perp = {means.max():.4e} (must be < 0), {int(np.sum(means >= 0))} h with >= 0")

    # 6. point values
    for (geo, lat), (Bref, Aref) in REF.items():
        o = observables(make_field(geo, lat, 5.0), 0.3, 1.0, 0.5)
        eB = abs(o["B_perp"] - Bref) / abs(Bref)
        check(eB <= 2e-3, f"6 point {geo} {lat} B_perp", f"{o['B_perp']:.4e} vs {Bref:.4e}, rel. err {eB:.2e}")
        if Aref is None:
            check(abs(o["A_perp"]) < 1e-13, f"6 point {geo} {lat} A_perp", f"{o['A_perp']:.3e} vs ~0 (< 1e-13)")
        else:
            eA = abs(o["A_perp"] - Aref) / abs(Aref)
            check(eA <= 2e-3, f"6 point {geo} {lat} A_perp", f"{o['A_perp']:.4e} vs {Aref:.4e}, rel. err {eA:.2e}")

    # 7. six-fold amplitude (SI Section S3)
    h, R, th = 5.0, 0.5, 0.1
    K = 45 * R**2 / (64 * np.pi * h**5)
    for geo, ref in SIX_FOLD.items():
        T = observables(make_field(geo, "hex", h), th, 1.0, R)["T"]
        ratio = T / (126 * K * np.sin(6 * th))
        check(abs(ratio - ref) <= 2e-3, f"7 six-fold amplitude {geo}", f"T/(126 K sin 6th) = {ratio:.4f} vs {ref}")

    # 8. error columns, relative to the largest |B_perp| over theta at the same h (B_perp crosses zero in some
    # lattices, so a point-wise relative error is meaningless there)
    for (lat, geo, R), c in data.items():
        rel = max(np.max(c["err_B"][m]) / np.max(np.abs(c["B_perp"][m])) for _, m in by_h(c))
        if geo == "bulk3d":
            check(rel < 1e-4, f"8 err_B {lat} {geo} R={R:g}", f"max_h max err_B / max|B_perp| = {rel:.3e} (< 1e-4)")
        else:
            report("INFO", f"8 err_B {lat} {geo} R={R:g}",
                   f"max_h max err_B / max|B_perp| = {rel:.3e} (conservative bound, report only)")

    # 9. decomposition B_perp = -Q R/(10 pi r_f^3) + B_img + r_f^2 gradp_img / 10, with r_f = sqrt(2) h / 4
    for (lat, geo, R), c in data.items():
        rf = R_F_OVER_H * c["h"]
        rebuilt = -R / (10 * np.pi * rf**3) + c["B_img"] + rf**2 / 10 * c["gradp_img"]
        dev = np.max(np.abs(rebuilt - c["B_perp"])) / np.max(np.abs(c["B_perp"]))
        check(dev < 1e-10, f"9 B decomposition {lat} {geo} R={R:g}", f"max deviation / max|B_perp| = {dev:.2e} (< 1e-10)")

    summary = f"\nTOTAL: {counts['PASS']} PASS, {counts['FAIL']} FAIL, {counts['WARN']} WARN, {counts['INFO']} INFO"
    lines.append(summary)
    print(summary)
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    with open(LOG, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return 1 if counts["FAIL"] else 0


if __name__ == "__main__":
    sys.exit(main())
