"""Numbers behind every quantitative or comparative statement of the results section.

Reads results/data (spinners-sweep configs/paper.toml) and prints one block per claim of the article.
Run: python3 scripts/summary.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from paper_figures import load, on_grid, rows_by_qh, minima, NAME  # noqa: E402

LATS = ["scl", "bcl", "fcl", "hex"]
GEOS = {"bulk3d": "3D", "monolayer": "2D"}
QH_SHOW = (0.1, 0.2, 0.25)


def velocities(lat, geo, R):
    c = load(lat, geo, R)
    qh, rows = rows_by_qh(c, "B_perp", on_grid(c, 81, np.pi))
    return qh, rows


def mean_abs_A(lat, geo, R):
    c = load(lat, geo, R)
    qh, rows = rows_by_qh(c, "A_perp", on_grid(c, 101, np.pi / 4))
    return qh, np.array([np.mean(np.abs(a)) for _, a in rows])


def stable_window(th, a):
    """Interval (t_minus, t_plus) around theta = 0 where A_perp < 0, with linear interpolation of the zeros.
    t_plus = 0 if A_perp > 0 just after 0. Values at +-pi/4 mean 'stable up to the edge of the map'."""
    out = []
    for side in (+1, -1):
        m = side * th > 1e-12
        t, y = side * th[m], a[m]
        o = np.argsort(t)
        t, y = t[o], y[o]
        if y[0] >= 0:
            out.append(0.0)
            continue
        k = np.argmax(y >= 0) if np.any(y >= 0) else None
        if k is None:
            out.append(side * t[-1])
        else:
            z = t[k - 1] - y[k - 1] * (t[k] - t[k - 1]) / (y[k] - y[k - 1])
            out.append(side * z)
    return out[1], out[0]


def stability(lat, geo, R):
    c = load(lat, geo, R)
    sel = on_grid(c, 101, np.pi / 4)
    qh, rows = rows_by_qh(c, "A_perp", sel)
    return qh, [stable_window(th, a) for th, a in rows]


def main():
    print("== B_perp/h averaged over theta (units 1/Q), R = Q/2")
    for geo, tag in GEOS.items():
        for lat in LATS:
            qh, rows = velocities(lat, geo, 0.5)
            bb = {q: np.mean(b) for q, (_, b) in zip(qh, rows)}
            amp = {q: (b.max() - b.min()) / abs(np.mean(b)) for q, (_, b) in zip(qh, rows)}
            print(f"  {tag} {NAME[lat]}: " + "  ".join(f"Q/h={q:.3g}: {bb[q]:+.3e} (spread/|mean| {amp[q]:.2f})"
                                                      for q in qh if np.round(q, 3) in QH_SHOW))
    print("\n== 2D / 3D ratio of the theta-averaged velocity, R = Q/2")
    for lat in LATS:
        q2, r2 = velocities(lat, "monolayer", 0.5)
        q3, r3 = velocities(lat, "bulk3d", 0.5)
        rat = [np.mean(b2) / np.mean(b3) for (_, b2), (_, b3) in zip(r2, r3)]
        print(f"  {NAME[lat]}: min {min(rat):.2f}  max {max(rat):.2f}")
    print("\n== R = Q/4 vs R = Q/2: ratio of the theta-averaged velocity and of the mean |A_perp|")
    for geo, tag in GEOS.items():
        for lat in LATS:
            _, r5 = velocities(lat, geo, 0.5)
            _, r25 = velocities(lat, geo, 0.25)
            rv = np.array([np.mean(b25) / np.mean(b5) for (_, b25), (_, b5) in zip(r25, r5)])
            _, a5 = mean_abs_A(lat, geo, 0.5)
            _, a25 = mean_abs_A(lat, geo, 0.25)
            ok = a5 > 1e-12 * np.abs(np.mean(r5[-1][1]))
            ra = a25[ok] / a5[ok] if ok.any() else np.array([np.nan])
            print(f"  {tag} {NAME[lat]}: velocity {rv.min():.2f} to {rv.max():.2f}   "
                  f"|A| {np.nanmin(ra):.2f} to {np.nanmax(ra):.2f}" + ("" if ok.any() else "  (A = 0 at R = Q/2)"))
    print("\n== stable window around theta = 0, (theta_minus, theta_plus) in rad, R = Q/2 and R = Q/4")
    for R in (0.5, 0.25):
        for geo, tag in GEOS.items():
            for lat in ("bcl", "fcl", "hex"):
                qh, w = stability(lat, geo, R)
                txt = "  ".join(f"Q/h={q:.3g}: ({a:+.2f},{b:+.2f}) width {b - a:.2f}" for q, (a, b) in zip(qh, w)
                                if np.round(q, 3) in QH_SHOW + (0.025, 0.05))
                print(f"  R={R:g} {tag} {NAME[lat]}: {txt}")
    print("\n== velocity reachable inside the stable window vs extremes over theta, R = Q/2")
    for geo, tag in GEOS.items():
        for lat in ("bcl", "fcl", "hex"):
            c = load(lat, geo, 0.5)
            qh, w = stability(lat, geo, 0.5)
            qv, rows = velocities(lat, geo, 0.5)
            for q in (0.2, 0.25):
                i = int(np.argmin(np.abs(qh - q)))
                a, b = w[i]
                th, B = rows[int(np.argmin(np.abs(qv - q)))]
                m = (th >= a) & (th <= b)
                Bs = B[m] if m.any() else np.array([np.nan])
                big = B[np.argmax(np.abs(B))]
                print(f"  {tag} {NAME[lat]} Q/h={q}: stable-window B_perp/h in [{np.min(Bs):+.3e}, {np.max(Bs):+.3e}]  "
                      f"mean over theta {np.mean(B):+.3e}  largest |B| {big:+.3e}")
    print("\n== minima of A_perp/h (Fig. 7 and S6)")
    for R in (0.5, 0.25):
        for lat, geo in (("hex", "monolayer"), ("hex", "bulk3d"), ("bcl", "bulk3d"), ("bcl", "monolayer"),
                         ("fcl", "bulk3d"), ("fcl", "monolayer")):
            c = load(lat, geo, R)
            qh, mp = minima(c, +1)
            _, mm = minima(c, -1)
            sel = [i for i, q in enumerate(qh) if np.round(q, 3) in QH_SHOW]
            print(f"  R={R:g} {GEOS[geo]} {NAME[lat]}: " + "  ".join(
                f"Q/h={qh[i]:.3g}: ccw {mp[i]:+.2e} cw {mm[i]:+.2e}" for i in sel))
    print("\n== growth of the mean |A_perp|/h with Q/h, R = Q/2 (ratio Q/h = 0.25 vs 0.125)")
    for geo, tag in GEOS.items():
        for lat in ("bcl", "fcl", "hex"):
            qh, a = mean_abs_A(lat, geo, 0.5)
            i1, i2 = int(np.argmin(np.abs(qh - 0.125))), int(np.argmin(np.abs(qh - 0.25)))
            mono = bool(np.all(np.diff(a) > 0))
            print(f"  {tag} {NAME[lat]}: ratio {a[i2] / a[i1]:.1f} (slope in log-log {np.log(a[i2] / a[i1]) / np.log(2):.1f}),"
                  f" monotonic: {mono}")


if __name__ == "__main__":
    main()
