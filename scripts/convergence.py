"""Convergence studies (corrected versions of the commented-out SI convergence figures).

Produces results/convergence/*.csv and results/figures/paper/conv_*.png:
  conv_ewald_<lattice>   relative error vs number of images n (index cube), several xi, for U_perp, U_par, rate_conf
  conv_xi                converged values vs xi (Ewald must not depend on the splitting parameter)
  conv_production        truncation of the sums with the production parameters of base_sums.rkt (alpha = 1)
  conv_direct            partial sums vs cutoff radius L for the monolayer and pure-2D direct sums

Run: PYTHONPATH=src python3 scripts/convergence.py
"""
import csv
import os
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from spinners.backends import make_field
from spinners.ewald import Ewald3D
from spinners.lattices import basis
from spinners.pair import observables

OUTD, OUTF = "results/convergence", "results/figures/paper"
INK, MUTED = "#1f1f1e", "#8a8984"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]  # validated categorical slots 1-4
MARK = ["o", "s", "^", "D"]
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
                     "savefig.dpi": 200, "savefig.bbox": "tight", "figure.facecolor": "#fcfcfb",
                     "axes.facecolor": "#fcfcfb",
                     # LaTeX look (Computer Modern) without needing a TeX install
                     "font.family": "serif", "font.serif": ["cmr10"], "mathtext.fontset": "cm",
                     "axes.formatter.use_mathtext": True, "axes.unicode_minus": False})
KEYS = [("U_perp", r"$U_\perp$"), ("U_par", r"$U_\parallel$"), ("rate_conf", r"rate$(V_{\mathrm{conf}})$")]
T0 = time.time()


def log(msg):
    """Progress message with elapsed wall time, flushed immediately."""
    print(f"[{(time.time() - T0) / 60:6.2f} min] {msg}", flush=True)


def ewald_obs(lattice, h, theta, R, **kw):
    ew = Ewald3D(basis(lattice, h), **kw)
    f = lambda X, P, F: (ew.velocity(X, P, F), np.zeros((len(X), 3)))
    return observables(f, theta, 1.0, R, "3d")


def write(name, rows):
    with open(os.path.join(OUTD, name + ".csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def conv_ewald(lattices=("scl", "bcl", "fcl", "hex"), h=5.0, theta=0.2, R=0.25):
    """R = Q/4 so that rate_conf is non-zero for every lattice, including SCL."""
    rows = []
    fig, axs = plt.subplots(len(lattices), 3, figsize=(11, 2.6 * len(lattices)), sharex=True)
    for i, lat in enumerate(lattices):
        V0 = abs(np.linalg.det(basis(lat, h)))
        log(f"conv_ewald  lattice {i + 1}/{len(lattices)} ({lat.upper()}): reference value")
        ref = ewald_obs(lat, h, theta, R, cut=60.0)
        for j, s in enumerate((0.25, 1.0, 4.0)):
            xi = s * V0 ** (2 / 3)
            ns = np.arange(1, 13)
            vals = []
            for n in ns:
                vals.append(ewald_obs(lat, h, theta, R, xi=xi, nk=n, nr=n))
                log(f"conv_ewald  lattice {i + 1}/{len(lattices)} ({lat.upper()})  xi {j + 1}/3 "
                    f"({s:g} V0^2/3)  n = {n:2d}/{ns[-1]}")
            for n, v in zip(ns, vals):
                rows.append({"lattice": lat, "h": h, "theta": theta, "R": R, "xi_over_V0^(2/3)": s, "n": n,
                             **{k: v[k] for k, _ in KEYS}, **{k + "_ref": ref[k] for k, _ in KEYS}})
            for c, (k, lab) in enumerate(KEYS):
                err = np.array([abs(v[k] - ref[k]) / abs(ref[k]) for v in vals])
                axs[i, c].semilogy(ns, np.maximum(err, 1e-16), color=SERIES[j], marker=MARK[j], ms=4, lw=1.5,
                                   label=rf"$\xi = {s:g}\,V_0^{{2/3}}$")
        for c, (k, lab) in enumerate(KEYS):
            axs[i, c].set_title(f"{lat.upper()} 3D: {lab}", fontsize=10, color=INK)
            axs[i, c].axhline(1e-12, color=MUTED, lw=0.6, ls=":")
        axs[i, 0].set_ylabel("relative error")
    for c in range(3):
        axs[-1, c].set_xlabel(r"$n$")
    fig.supxlabel(r"$n$: images with $|m_i| \leq n$ in real and reciprocal space", color=MUTED, fontsize=9)
    axs[0, 0].legend(frameon=False, fontsize=8)
    fig.savefig(os.path.join(OUTF, "conv_ewald.png"))
    plt.close(fig)
    write("conv_ewald", rows)


def conv_xi(lattices=("scl", "bcl", "fcl", "hex"), h=5.0, theta=0.2, R=0.25):
    rows = []
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.2))
    svals = np.geomspace(0.1, 10, 9)
    for i, lat in enumerate(lattices):
        V0 = abs(np.linalg.det(basis(lat, h)))
        vals = []
        for q, s in enumerate(svals):
            vals.append(ewald_obs(lat, h, theta, R, xi=s * V0 ** (2 / 3)))
            log(f"conv_xi     lattice {i + 1}/{len(lattices)} ({lat.upper()})  xi {q + 1}/{len(svals)} "
                f"({s:.3g} V0^2/3)")
        for s, v in zip(svals, vals):
            rows.append({"lattice": lat, "xi_over_V0^(2/3)": s, **{k: v[k] for k, _ in KEYS}})
        for c, (k, lab) in enumerate(KEYS):
            y = np.array([v[k] for v in vals])
            axs[c].semilogx(svals, y / y[len(y) // 2] - 1, color=SERIES[i], marker=MARK[i], ms=4, lw=1.5,
                            label=lat.upper())
    for c, (k, lab) in enumerate(KEYS):
        axs[c].set_title(lab, color=INK)
        axs[c].set_xlabel(r"$\xi / V_0^{2/3}$")
        axs[c].set_ylim(-2e-12, 2e-12)
    axs[0].set_ylabel(r"value / value($\xi=V_0^{2/3}$) - 1")
    axs[0].legend(frameon=False, fontsize=8)
    fig.suptitle(f"Ewald sums are independent of the splitting parameter (h = {h:g}, R = {R:g}Q)",
                 color=MUTED, fontsize=9, y=-0.02)
    fig.savefig(os.path.join(OUTF, "conv_xi.png"))
    plt.close(fig)
    write("conv_xi", rows)


def conv_production(theta=0.2, R=0.5):
    """alpha = xi = 1 fixed as in base_sums.rkt, n as in parametros_experimentos.csv (marked)."""
    prod_n = {7: 10, 10: 15, 16: 20}
    rows = []
    fig, ax = plt.subplots(figsize=(6, 3.8))
    for i, (h, n_prod) in enumerate(prod_n.items()):
        ref = ewald_obs("hex", float(h), theta, R)["rate_conf"]   # converged, balanced xi
        ns = np.unique(np.r_[np.arange(4, 2 * h + 1, max(1, h // 6)), n_prod, 2 * h, 3 * h])
        err = []
        for q, n in enumerate(ns):
            t = time.time()
            v = ewald_obs("hex", float(h), theta, R, xi=1.0, nk=int(n), nr=2)["rate_conf"]
            err.append(abs(v - ref) / abs(ref))
            log(f"conv_production  h = {h} ({i + 1}/{len(prod_n)})  n = {n:3d} ({q + 1}/{len(ns)}, "
                f"{(2 * n + 1) ** 3:,} k-vectors, {time.time() - t:5.1f} s)  rel. error = {err[-1]:.2e}"
                + ("   <- production n" if n == n_prod else ""))
            rows.append({"h": h, "Q_over_h": 1 / h, "n": int(n), "rate_conf": v, "converged": ref})
        ax.semilogy(ns / h, np.maximum(err, 1e-16), color=SERIES[i], marker=MARK[i], ms=4, lw=1.5,
                    label=f"$Q/h = {1 / h:.3f}$ ($h = {h}$)")
        e_prod = err[list(ns).index(n_prod)]
        ax.plot(n_prod / h, e_prod, marker="*", ms=14, color=SERIES[i], mec=INK, mew=0.8, ls="none")
    ax.set_xlabel(r"$n / h$")
    ax.set_ylabel(r"relative error of rate$(V_{\mathrm{conf}})$")
    ax.set_title(r"HEX 3D with the production setting $\alpha=\xi=1$ (stars: $n$ used in production)",
                 fontsize=9, color=INK)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(os.path.join(OUTF, "conv_production.png"))
    plt.close(fig)
    write("conv_production", rows)


def conv_direct(theta=0.2, R=0.5, h=5.0):
    Ls = (10, 15, 20, 30, 45, 60, 90, 120, 160, 200)
    rows = []
    fig, axs = plt.subplots(2, 2, figsize=(10, 6.2))
    for gi, geo in enumerate(("monolayer", "pure2d")):
        for i, lat in enumerate(("scl", "bcl", "fcl", "hex")):
            log(f"conv_direct geometry {gi + 1}/2 ({geo})  lattice {i + 1}/4 ({lat.upper()})")
            fld, k = make_field(geo, lat, h, L_over_h=Ls)
            from spinners.pair import frame, sources, eval_points, _scalars
            qh, e = frame(theta)
            P, F = sources(theta, 1.0, R)
            series = _scalars(fld.vs_L(eval_points(theta, 1.0, R), P, F), qh, e, 1.0)
            fin = observables(make_field(geo, lat, h)[0], theta, 1.0, R, k)
            for c, key in enumerate(("rate_conf", "U_perp")):
                y = series[key]
                if lat == "scl" and key == "rate_conf":
                    continue  # identically zero for R = Q/2
                rel = np.abs(y - fin[key]) / abs(fin[key])
                axs[gi, c].loglog(fld.L / fld.L[0] * Ls[0], np.maximum(rel, 1e-16), color=SERIES[i],
                                  marker=MARK[i], ms=4, lw=1.5, label=lat.upper())
                for L, val in zip(Ls, y):
                    rows.append({"geometry": geo, "lattice": lat, "L_over_h": L, "observable": key,
                                 "partial_sum": val, "extrapolated": fin[key]})
        axs[gi, 0].set_title(rf"{geo}: rate$(V_{{\mathrm{{conf}}}})$", color=INK, fontsize=10)
        axs[gi, 1].set_title(f"{geo}: $U_\\perp$", color=INK, fontsize=10)
        axs[gi, 0].set_ylabel(r"relative error  $|S_L - S_\infty| / |S_\infty|$")
    for c in range(2):
        axs[1, c].set_xlabel(r"cutoff radius $L / V_0^{1/3}$")
    axs[0, 0].legend(frameon=False, fontsize=8)
    fig.savefig(os.path.join(OUTF, "conv_direct.png"))
    plt.close(fig)
    write("conv_direct", rows)


if __name__ == "__main__":
    os.makedirs(OUTD, exist_ok=True)
    os.makedirs(OUTF, exist_ok=True)
    steps = [("conv_xi", conv_xi), ("conv_ewald", conv_ewald), ("conv_direct", conv_direct),
             ("conv_production", conv_production)]
    for i, (name, fn) in enumerate(steps, 1):
        log(f"=== step {i}/{len(steps)}: {name} ===")
        fn()
        log(f"=== step {i}/{len(steps)}: {name} done ===")
    log(f"all done. CSV in {OUTD}/, figures in {OUTF}/")
