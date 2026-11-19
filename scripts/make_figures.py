"""Figures from results/data/*.csv  ->  results/figures/*.png"""
import csv
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np

DATA, OUT = "results/data", "results/figures"
GEO = ["bulk3d", "monolayer", "pure2d"]
LABEL = {"bulk3d": "Bulk 3D (Ewald)", "monolayer": "Monolayer in 3D fluid", "pure2d": "Pure 2D (line forces)"}
COLOR = {"bulk3d": "#2a78d6", "monolayer": "#eb6834", "pure2d": "#1baf7a"}   # validated categorical slots 1-3
MARK = {"bulk3d": "o", "monolayer": "s", "pure2d": "^"}
STABLE, UNSTABLE, NEUTRAL, INK, MUTED = "#2a78d6", "#e34948", "#f0efec", "#1f1f1e", "#8a8984"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                     "savefig.dpi": 200, "savefig.bbox": "tight", "figure.facecolor": "#fcfcfb",
                     "axes.facecolor": "#fcfcfb"})


def load(lattice, geo, R):
    fn = os.path.join(DATA, f"{lattice}_{geo}_R{R:g}.csv")
    with open(fn) as fh:
        rows = list(csv.DictReader(fh))
    cols = {k: np.array([float(r[k]) if k not in ("geometry", "lattice") else 0 for r in rows]) for k in rows[0]}
    return cols


def grid(c, key):
    qh = np.unique(c["Q_over_h"]); th = np.unique(c["theta"])
    Z = np.full((len(qh), len(th)), np.nan)
    iq = np.searchsorted(qh, c["Q_over_h"]); it = np.searchsorted(th, c["theta"])
    Z[iq, it] = c[key]
    return th, qh, Z


def stability_maps(lattice="hex"):
    for R in (0.5, 0.25):
        fig, axs = plt.subplots(1, 3, figsize=(11, 3.4), sharey=True)
        for ax, g in zip(axs, GEO):
            th, qh, Z = grid(load(lattice, g, R), "A_conf")
            scale = np.nanmax(np.abs(Z), axis=1, keepdims=True)
            rel = Z / np.where(scale > 0, scale, 1)
            col = np.where(rel < -1e-9, -1, np.where(rel > 1e-9, 1, 0))
            ax.pcolormesh(th, qh, col, cmap=matplotlib.colors.ListedColormap([STABLE, NEUTRAL, UNSTABLE]),
                          vmin=-1.5, vmax=1.5, shading="nearest")
            for t in (-np.pi / 6, np.pi / 6):
                ax.axvline(t, color=INK, lw=0.8, ls=":")
            ax.set_title(LABEL[g], color=INK, fontsize=10)
            ax.set_xlabel(r"$\theta$ (rad)")
        axs[0].set_ylabel("Q/h")
        h1 = matplotlib.patches.Patch(color=STABLE, label=r"$A_\perp<0$ (pushes Q towards $a^{(1)}$)")
        h2 = matplotlib.patches.Patch(color=UNSTABLE, label=r"$A_\perp>0$ (pushes Q away)")
        fig.legend(handles=[h1, h2], loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.08))
        fig.suptitle(f"Corrected sign of the paper's criterion (V_conf), hexagonal lattice, R = {R:g}Q.  "
                     r"Dotted: $\theta=\pm\pi/6$", y=-0.04, fontsize=9, color=MUTED)
        fig.savefig(os.path.join(OUT, f"{lattice}_stability_maps_R{R:g}.png"))
        plt.close(fig)


def rates_vs_theta(lattice="hex", R=0.5):
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.4))
    for ax, qh0 in zip(axs, (0.25, 0.1)):
        for g in GEO:
            c = load(lattice, g, R)
            m = np.isclose(c["Q_over_h"], qh0)
            o = np.argsort(c["theta"][m])
            y = c["rate_conf"][m][o]
            ax.plot(c["theta"][m][o], y / np.abs(y).max(), color=COLOR[g], lw=2, label=LABEL[g])
        c = load(lattice, "bulk3d", R); m = np.isclose(c["Q_over_h"], qh0); o = np.argsort(c["theta"][m])
        sc = np.abs(c["rate_conf"][m]).max()
        ax.plot(c["theta"][m][o], c["omega_qdot"][m][o] / sc, color=INK, lw=1.2, ls="--",
                label=r"$\dot Q$ criterion (all geometries, $\equiv 0$)")
        ax.axhline(0, color=MUTED, lw=0.6)
        for t in (-np.pi / 6, np.pi / 6):
            ax.axvline(t, color=MUTED, lw=0.6, ls=":")
        ax.set_title(f"Q/h = {qh0:g}  (the three normalized curves coincide)", color=INK, fontsize=10)
        ax.set_xlabel(r"$\theta$ (rad)")
    axs[0].set_ylabel(r"$(\hat Q\times V_{conf})\cdot\hat n$ / max")
    axs[0].legend(frameon=False, fontsize=8, loc="lower left")
    fig.savefig(os.path.join(OUT, f"{lattice}_rates_vs_theta_R{R:g}.png"))
    plt.close(fig)


def amplitude_vs_qh(lattice="hex"):
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    for g in GEO:
        for R, ls in ((0.5, "-"), (0.25, "--")):
            c = load(lattice, g, R)
            qh = np.unique(c["Q_over_h"])
            amp = np.array([np.abs(c["rate_conf"][np.isclose(c["Q_over_h"], q)]).max() for q in qh])
            p = np.polyfit(np.log(qh[qh < 0.08]), np.log(amp[qh < 0.08]), 1)[0]
            ax.loglog(qh, amp, ls, color=COLOR[g], marker=MARK[g], ms=4, lw=1.5,
                      label=f"{LABEL[g]}, R={R:g}Q  (slope {p:.2f})")
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xticks([0.025, 0.05, 0.1, 0.25]); ax.set_xticklabels(["0.025", "0.05", "0.1", "0.25"])
    ax.set_xlabel("Q/h"); ax.set_ylabel(r"max$_\theta$ |rate(V$_{conf}$)|  [$\zeta V^S/(\mu Q)$]")
    ax.legend(frameon=False, fontsize=7)
    fig.savefig(os.path.join(OUT, f"{lattice}_vconf_amplitude.png"))
    plt.close(fig)


def velocity_vs_qh(lattice="hex"):
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.6))
    for g in GEO:
        for R, ls in ((0.5, "-"), (0.25, "--")):
            c = load(lattice, g, R)
            qh = np.unique(c["Q_over_h"])
            U = np.array([c["U_perp"][np.isclose(c["Q_over_h"], q)].mean() for q in qh])
            U0 = np.array([c["U_perp_isolated"][np.isclose(c["Q_over_h"], q)].mean() for q in qh])
            axs[0].plot(qh, -U, ls, color=COLOR[g], marker=MARK[g], ms=4, label=f"{LABEL[g]}, R={R:g}Q")
            axs[1].plot(qh, (U - U0) / U0, ls, color=COLOR[g], marker=MARK[g], ms=4)
    axs[0].set_ylabel(r"$-U_\perp$  [$\zeta V^S/(\mu Q)$]  (motion along $-\hat Q\times\hat n$)")
    axs[1].set_ylabel("lattice correction  (U - U$_{isolated}$)/U$_{isolated}$")
    for ax in axs:
        ax.set_xlabel("Q/h")
    fig.legend(*axs[0].get_legend_handles_labels(), loc="lower center", ncol=3, frameon=False,
               fontsize=8, bbox_to_anchor=(0.5, -0.22))
    fig.savefig(os.path.join(OUT, f"{lattice}_velocity_vs_qh.png"))
    plt.close(fig)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    lattice = sys.argv[1] if len(sys.argv) > 1 else "hex"
    jobs = [stability_maps, rates_vs_theta, amplitude_vs_qh, velocity_vs_qh]
    for i, fn in enumerate(jobs, 1):
        print(f"[{i}/{len(jobs)}] {fn.__name__}({lattice}) ...", end=" ", flush=True)
        fn(lattice)
        print("ok", flush=True)
    print("figures written to", OUT, flush=True)
