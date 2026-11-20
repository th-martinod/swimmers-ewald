"""Convergence studies of every lattice sum used for the figures.

Output: results/convergence/*.csv and results/figures/convergence/*.{eps,pdf} (+ previews)
  conv_ewald       3D lattices: relative error of T and B_perp vs the Gaussian cutoff of the Ewald sums,
                   for three splitting parameters xi (production: cut = 36)
  conv_xi          3D lattices: independence of T and B_perp on xi
  conv_direct      monolayers: partial sums vs cutoff radius L and extrapolation to L -> infinity
  conv_images      image sum of the pressure gradient used in the sub-volume average, 3D and monolayer
                   (3D: partial sum at L = 24h, monolayer: extrapolated with the other direct sums)
  conv_ball        sub-volume average: independent Fourier sum (ball form factor) vs the mean-value
                   evaluation used in production, as the number of wave vectors grows
  conv_production  every lattice, h and R of the sweep: production parameters vs a reference computed
                   with larger cutoffs

Point of evaluation: theta = 0.3, R = Q/4 unless stated (R = Q/4 so that T is non-zero in the SCL).
Run: PYTHONPATH=src python3 scripts/convergence.py
"""
import csv
import os
import sys
import time

import numpy as np
from matplotlib.ticker import NullFormatter

sys.path.insert(0, os.path.dirname(__file__))
from paper_figures import COLOR, MARK, NAME, LATS, INK, MUTED, plt, save  # noqa: E402  (shared style)

from spinners.backends import make_field
from spinners.direct import extrapolate
from spinners.ewald import Ewald3D
from spinners.lattices import basis, reciprocal
from spinners.pair import observables, sources, frame, eval_points, _scalars

OUTD = "results/convergence"
TH, R0, H0 = 0.3, 0.25, 5.0
PROD_L = (30, 45, 60, 90, 120)
REF_L = (60, 90, 120, 160, 240)
T0 = time.time()
SHADES = ["#2a4fd6", "#b58c00", "#cf2a2a"]
XI_MARK = ["o", "s", "^"]


def log(msg):
    print(f"[{(time.time() - T0) / 60:6.2f} min] {msg}", flush=True)


def write(name, rows):
    os.makedirs(OUTD, exist_ok=True)
    with open(os.path.join(OUTD, name + ".csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def bulk_obs(lat, h, th, R, **ew):
    f = make_field("bulk3d", lat, h)
    if ew:
        f.ewald = Ewald3D(f.A, **ew)
    return observables(f, th, 1.0, R)


def series(field, th, R):
    """Partial sums of T and B_perp for every cutoff radius of a direct (monolayer) field."""
    qh, e = frame(th)
    P, F = sources(th, 1.0, R)
    v = field.velocity_vs_L(eval_points(th, 1.0, R), P, F)
    gp = field.grad_p_vs_L(P, F)
    return _scalars(v, gp, th, qh, e, P, F, field.r_f)


# ------------------------------------------------------------------------------------------------
def log_xticks(ax, ticks):
    """Explicit major ticks on a log x axis and no minor labels (the default ones overlap)."""
    ax.set_xticks(ticks, [f"${t:g}$" for t in ticks])
    ax.xaxis.set_minor_formatter(NullFormatter())


def conv_ewald():
    rows = []
    cuts = np.arange(2, 42, 2)
    fig, axs = plt.subplots(2, 4, figsize=(9.5, 4.4), sharex=True, sharey=True)
    for i, lat in enumerate(LATS):
        V0 = abs(np.linalg.det(basis(lat, H0)))
        ref = bulk_obs(lat, H0, TH, R0, cut=60.0)
        for j, s in enumerate((0.25, 1.0, 4.0)):
            vals = [bulk_obs(lat, H0, TH, R0, xi=s * V0 ** (2 / 3), cut=float(c)) for c in cuts]
            for c, v in zip(cuts, vals):
                rows.append({"lattice": lat, "xi_over_V0^(2/3)": s, "cut": c, "T": v["T"], "B_perp": v["B_perp"],
                             "T_ref": ref["T"], "B_ref": ref["B_perp"]})
            for r, k in enumerate(("T", "B_perp")):
                err = np.maximum([abs(v[k] - ref[k]) / abs(ref[k]) for v in vals], 1e-17)
                axs[r, i].semilogy(cuts, err, color=SHADES[j], marker=XI_MARK[j], ms=3.5, mfc="white", lw=0.8,
                                   label=rf"$\xi={s:g}\,V_0^{{2/3}}$")
        log(f"conv_ewald {lat}")
        axs[0, i].set_title(f"{NAME[lat]} 3D", color=INK)
    for r, lab in enumerate((r"$T$", r"$B_\perp$")):
        axs[r, 0].set_ylabel(f"relative error of {lab}")
        for ax in axs[r]:
            ax.axvline(36, color=MUTED, lw=0.6, ls="--")
    for ax in axs[1]:
        ax.set_xlabel(r"cutoff  $\pi\xi k^2,\ \pi r^2/\xi \leq$ cut")
    axs[0, 0].legend(loc="lower left", fontsize=7)
    axs[0, 0].text(36.5, 1e-2, "production", color=MUTED, fontsize=7, rotation=90, va="bottom")
    fig.tight_layout()
    save(fig, "convergence", "conv_ewald")
    write("conv_ewald", rows)


def conv_xi():
    rows = []
    svals = np.geomspace(0.1, 10, 9)
    fig, axs = plt.subplots(1, 2, figsize=(6.4, 2.6), sharey=True)
    for lat in LATS:
        V0 = abs(np.linalg.det(basis(lat, H0)))
        vals = [bulk_obs(lat, H0, TH, R0, xi=s * V0 ** (2 / 3)) for s in svals]
        for s, v in zip(svals, vals):
            rows.append({"lattice": lat, "xi_over_V0^(2/3)": s, "T": v["T"], "B_perp": v["B_perp"]})
        for ax, k in zip(axs, ("T", "B_perp")):
            y = np.array([v[k] for v in vals])
            dev = np.abs(y / y[len(y) // 2] - 1)
            dev[len(y) // 2] = np.nan                                   # the reference point itself
            ax.loglog(svals, np.maximum(dev, 1e-17), color=COLOR[lat], marker=MARK[lat], ms=3.5, mfc="white",
                      lw=0.8, label=NAME[lat])
        log(f"conv_xi {lat}")
    for ax, lab in zip(axs, (r"$T$", r"$B_\perp$")):
        ax.set_xlabel(r"$\xi/V_0^{2/3}$")
        ax.set_title(f"{lab} (3D lattices)", color=INK)
        ax.set_ylim(1e-17, 1e-8)
    axs[0].set_ylabel(r"$|$value$/$value$(\xi=V_0^{2/3})-1|$")
    axs[0].legend(fontsize=7, ncol=2)
    fig.tight_layout()
    save(fig, "convergence", "conv_xi")
    write("conv_xi", rows)


def conv_direct():
    rows = []
    Ls = np.array([5, 7.5, 10, 15, 20, 30, 45, 60, 90, 120, 160, 240], float)
    fig, axs = plt.subplots(1, 2, figsize=(6.4, 2.8), sharey=True)
    for lat in LATS:
        f = make_field("monolayer", lat, H0, L_over_h=tuple(Ls))
        s = series(f, TH, R0)
        sel = np.isin(Ls, REF_L)
        for ax, k in zip(axs, ("T", "B_perp")):
            ref, _ = extrapolate(f.L[sel], s[k][sel][:, None])
            prod, _ = extrapolate(f.L[np.isin(Ls, PROD_L)], s[k][np.isin(Ls, PROD_L)][:, None])
            err = np.abs(s[k] - ref[0]) / abs(ref[0])
            ax.loglog(Ls, err, color=COLOR[lat], marker=MARK[lat], ms=3.5, mfc="white", lw=0.8, label=NAME[lat])
            ax.plot(Ls[-1] * 1.5, max(abs(prod[0] - ref[0]) / abs(ref[0]), 1e-12), ls="none", marker=MARK[lat],
                    ms=5, color=COLOR[lat])
            for L, val in zip(Ls, s[k]):
                rows.append({"lattice": lat, "quantity": k, "L_over_h": L, "partial_sum": val,
                             "extrapolated_ref": ref[0], "extrapolated_production": prod[0]})
        log(f"conv_direct {lat}")
    for ax, lab in zip(axs, (r"$T$", r"$B_\perp$")):
        ax.set_xlabel(r"cutoff radius $L/h$")
        ax.set_title(f"{lab} (monolayers)", color=INK)
        ax.axvspan(PROD_L[0], PROD_L[-1], color="#eeeeec", zorder=0, lw=0)
    axs[0].set_ylabel("relative distance to $L\\to\\infty$")
    axs[0].text(0.97, 0.95, "filled: production\nextrapolation", transform=axs[0].transAxes, ha="right",
                va="top", fontsize=7, color=MUTED)
    axs[1].legend(fontsize=7, loc="lower left")
    fig.tight_layout()
    save(fig, "convergence", "conv_direct")
    write("conv_direct", rows)


def conv_images():
    rows = []
    fig, axs = plt.subplots(1, 2, figsize=(6.4, 2.8))
    cases = (("bulk3d", np.array([2, 3, 4, 6, 8, 12, 16, 20, 24, 32, 40], float), axs[0]),
             ("monolayer", np.array([5, 10, 20, 30, 45, 60, 90, 120, 160, 240], float), axs[1]))
    for geo, Ls, ax in cases:
        for lat in LATS:
            f = make_field(geo, lat, H0, **({"gradp_L_3d": tuple(Ls)} if geo == "bulk3d" else {"L_over_h": tuple(Ls)}))
            _, e = frame(TH)
            P, F = sources(TH, 1.0, R0)
            a = f.r_f
            term = a**2 / 10 * f.grad_p_vs_L(P, F) @ e
            if geo == "bulk3d":     # absolutely convergent (d^-5): reference = largest cutoff, not plotted
                ref = term[-1:]
                Lp, tp = Ls[:-1], term[:-1]
            else:                   # part of the monolayer direct sums, extrapolated as in production
                ref, _ = extrapolate(f.L[-5:], term[-5:, None])
                Lp, tp = Ls, term
            B = observables(f, TH, 1.0, R0)["B_perp"]
            err = np.maximum(np.abs(tp - ref[0]) / abs(B), 1e-17)
            ax.loglog(Lp, err, color=COLOR[lat], marker=MARK[lat], ms=3.5, mfc="white", lw=0.8,
                      label=f"{NAME[lat]}  (term$/B_\\perp$ = {ref[0] / B:.2g})")
            for L, t in zip(Ls, term):
                rows.append({"geometry": geo, "lattice": lat, "L_over_h": L, "partial_sum": t,
                             "extrapolated": ref[0], "B_perp": B})
            log(f"conv_images {geo} {lat}")
        ax.set_xlabel(r"cutoff radius $L/h$")
        ax.set_title("3D lattices" if geo == "bulk3d" else "monolayers", color=INK)
        ax.legend(fontsize=6.5, loc="lower left")
        if geo == "bulk3d":
            log_xticks(ax, (2, 4, 8, 16, 32))
            ax.axvline(24, color=MUTED, lw=0.6, ls="--")
            ax.text(24.8, 3e-3, "production", color=MUTED, fontsize=7, rotation=90, va="bottom")
        else:
            ax.axvspan(PROD_L[0], PROD_L[-1], color="#eeeeec", zorder=0, lw=0)
    axs[0].set_ylabel(r"$|S_L - S_\infty|/|B_\perp|$")
    fig.tight_layout()
    save(fig, "convergence", "conv_images")
    write("conv_images", rows)


def _ball(k, a):
    x = 2 * np.pi * np.linalg.norm(k, axis=-1) * a
    return 3 * (np.sin(x) - x * np.cos(x)) / x**3


def _cube(n, d):
    m = np.arange(-n, n + 1)
    g = np.meshgrid(*([m] * d), indexing="ij")
    return np.stack([x.ravel() for x in g], 1).astype(float)


def _fourier(k, w, P, F):
    k2 = np.sum(k * k, axis=1)
    kh = k / np.sqrt(k2)[:, None]
    g = np.cos(2 * np.pi * k @ P.T) @ F
    g -= kh * np.sum(kh * g, axis=1)[:, None]
    return (w / k2) @ g


def fourier_ball(f, P, F, n, nq=240):
    a = f.r_f
    if f.geometry == "bulk3d":
        B, V0 = reciprocal(f.A)
        k = _cube(n, 3) @ B
        k = k[np.any(k != 0, axis=1)]
        return _fourier(k, _ball(k, a) / (4 * np.pi**2 * V0), P, F)
    a1, a2 = f.A[0], f.A[1]
    m = np.cross(a1, a2)
    A0 = np.linalg.norm(m)
    m /= A0
    G = np.array([a1, a2])
    kp = _cube(n, 2) @ np.linalg.solve(G @ G.T, G)
    x, wx = np.polynomial.legendre.leggauss(nq)
    phi, wphi = x * np.pi / 2, wx * np.pi / 2
    s = np.maximum(np.linalg.norm(kp, axis=1), 1 / np.linalg.norm(a1))
    kap = s[:, None] * np.tan(phi)[None]
    jac = s[:, None] / np.cos(phi)[None] ** 2
    k = (kp[:, None, :] + kap[..., None] * m).reshape(-1, 3)
    w = (jac * wphi[None]).reshape(-1) * _ball(k, a) / (4 * np.pi**2 * A0)
    return _fourier(k, w, P, F)


def conv_ball():
    rows = []
    fig, axs = plt.subplots(1, 2, figsize=(6.4, 2.8), sharey=True)
    for geo, ns, ax in (("bulk3d", (4, 6, 8, 12, 16, 24, 32, 48), axs[0]),
                        ("monolayer", (4, 6, 8, 12, 16, 24, 32), axs[1])):
        for lat in LATS:
            f = make_field(geo, lat, H0, L_over_h=REF_L)
            P, F = sources(TH, 1.0, R0)
            _, e = frame(TH)
            B = observables(f, TH, 1.0, R0)["B_perp"]
            vals = [fourier_ball(f, P, F, n) @ e for n in ns]
            ax.loglog(ns, np.abs(np.array(vals) - B) / abs(B), color=COLOR[lat], marker=MARK[lat], ms=3.5,
                      mfc="white", lw=0.8, label=NAME[lat])
            for n, v in zip(ns, vals):
                rows.append({"geometry": geo, "lattice": lat, "n": n, "B_fourier": v, "B_mean_value": B})
            log(f"conv_ball {geo} {lat}")
        log_xticks(ax, (4, 8, 16, 32))
        ax.set_xlabel(r"wave vectors with $|m_i|\leq n$")
        ax.set_title("3D lattices" if geo == "bulk3d" else "monolayers", color=INK)
    axs[0].set_ylabel(r"$|B_\perp^{\rm Fourier}-B_\perp|/|B_\perp|$")
    axs[1].legend(fontsize=7, loc="lower left")
    fig.tight_layout()
    save(fig, "convergence", "conv_ball")
    write("conv_ball", rows)


def conv_production():
    """Production parameters (sweep) against a reference with larger cutoffs, every lattice, h and R."""
    rows = []
    hs = (4, 5, 6, 7, 8, 9, 10, 12, 16, 20, 24, 28, 32, 36, 40)
    fig, axs = plt.subplots(1, 2, figsize=(6.4, 2.8), sharey=True)
    for ax, geo in zip(axs, ("bulk3d", "monolayer")):
        for lat in LATS:
            err = []
            for h in hs:
                if geo == "bulk3d":
                    fp = make_field(geo, lat, h)
                    fr = make_field(geo, lat, h, gradp_L_3d=(32, 40))
                    fr.ewald = Ewald3D(fr.A, cut=60.0)
                else:
                    fp = make_field(geo, lat, h, L_over_h=PROD_L)
                    fr = make_field(geo, lat, h, L_over_h=REF_L)
                e_h = 0.0
                for R in (0.5, 0.25):
                    p, r = observables(fp, TH, 1.0, R), observables(fr, TH, 1.0, R)
                    for k in ("A_perp", "B_perp"):
                        if not (lat == "scl" and R == 0.5 and k == "A_perp"):   # exact zero (round-off only)
                            rel = abs(p[k] - r[k]) / abs(r[k])
                            e_h = max(e_h, rel)
                            rows.append({"geometry": geo, "lattice": lat, "h": h, "R": R, "quantity": k,
                                         "production": p[k], "reference": r[k], "relative_error": rel})
                err.append(max(e_h, 1e-17))
            ax.semilogy(1 / np.array(hs), err, color=COLOR[lat], marker=MARK[lat], ms=3.5, mfc="white", lw=0.8,
                        label=NAME[lat])
            log(f"conv_production {geo} {lat}")
        ax.set_xlabel(r"$Q/h$")
        ax.set_title("3D lattices" if geo == "bulk3d" else "monolayers", color=INK)
    axs[0].set_ylabel(r"max relative error ($A_\perp$, $B_\perp$)")
    axs[1].legend(fontsize=7, ncol=2)
    fig.tight_layout()
    save(fig, "convergence", "conv_production")
    write("conv_production", rows)


if __name__ == "__main__":
    for job in (conv_xi, conv_ewald, conv_images, conv_direct, conv_ball, conv_production):
        job()
    log("done")
