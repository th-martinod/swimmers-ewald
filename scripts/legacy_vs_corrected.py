"""Reproduce the diagnosis of base_sums.rkt with the numpy port (src/spinners/legacy.py).

Run:  PYTHONPATH=src python3 scripts/legacy_vs_corrected.py
"""
import numpy as np

from spinners import legacy
from spinners.backends import make_field
from spinners.pair import observables

pi = np.pi
print("1) Q rotated by 30 deg in the Fourier part (HEX, Q=1, h=5, alpha=1, theta=0.2)")
for dim, nF in (("2D", 60), ("3D", 22)):
    bug = legacy.conf_total(0.2, "hex", dim, 5.0, 1.0, nF, 4, bug=True)[0]
    ok = legacy.conf_total(0.2, "hex", dim, 5.0, 1.0, nF, 4, bug=False)[0]
    print(f"   {dim}: as coded {bug:+.4e}   same Q in both parts {ok:+.4e}")

print("2) Spurious prefactor: legacy(consistent) / corrected bulk = h^3/(4 pi^3)")
for h in (4.0, 5.0, 8.0):
    ok = legacy.conf_total(0.6, "hex", "3D", h, 1.0, int(4 * h) + 4, 4, bug=False)[0]
    fld, k = make_field("bulk3d", "hex", h)
    new = observables(fld, 0.6, 1.0, 0.5, k)["A_conf"]
    print(f"   h={h:g}: ratio {ok / new:.4f}   h^3/(4pi^3) {h**3 / (4 * pi**3):.4f}")

print("3) Truncation at the production n (HEX 2D, alpha=1, theta=0.2)")
for h, n in ((7, 10), (10, 15), (16, 20), (40, 50)):
    a = legacy.conf_total(0.2, "hex", "2D", float(h), 1.0, n, 2, bug=False)[0]
    ref = legacy.conf_total(0.2, "hex", "2D", float(h), 1.0, 4 * h, 2, bug=False)[0]
    print(f"   h={h:2d} n={n:2d}: {a:+.3e} vs converged {ref:+.3e}  (rel. error {abs(a - ref) / abs(ref):.1e})")

print("4) Lattice factor A in the velocity (3D, h=5): relative variation with theta")
ths = np.arange(-pi, pi, pi / 50)
nv = np.array([0, 0, 1.0])
for latt in ("hex", "scl"):
    for useA in (True, False):
        v = []
        for t in ths:
            q = np.array([np.cos(t), np.sin(t), 0])
            d = legacy.nz(np.cross(q, nv))
            v.append((legacy.rot_fourier(q, nv, 0.5, 5.0, 1.0, 18, latt, "3D", useA)
                      + legacy.rot_real(q, nv, 0.5, 5.0, 1.0, 3, latt, "3D", useA)) @ d)
        v = np.array(v)
        print(f"   {latt} {'with A ' if useA else 'without A'}: (max-min)/|mean| = {(v.max() - v.min()) / abs(v.mean()):.2e}")
