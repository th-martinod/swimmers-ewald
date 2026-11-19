"""Corrected versions of every data figure of the article and the SI.

Input : results/data_paper/{lattice}_{geometry}_R{R}.csv   (spinners-sweep configs/paper.toml)
Output: results/figures/paper/*.png  (see docs/FIGURE_MAP.md for the correspondence with the article)

Run: python3 scripts/paper_figures.py
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

DATA, OUT = "results/data_paper", "results/figures/paper"
LATS = ["scl", "bcl", "fcl", "hex"]
NAME = {"scl": "SCL", "bcl": "BCL", "fcl": "FCL", "hex": "HEX"}
COLOR = {"scl": "#2a78d6", "bcl": "#eb6834", "fcl": "#1baf7a", "hex": "#eda100"}  # validated, fixed order
MARK = {"scl": "s", "bcl": "^", "fcl": "o", "hex": "D"}
GEO = {"bulk3d": "3D", "monolayer": "monolayer (3D fluid)", "pure2d": "pure 2D (line forces)"}
STABLE, UNSTABLE, NEUTRAL, INK, MUTED = "#2a78d6", "#e34948", "#f0efec", "#1f1f1e", "#8a8984"
PI4 = np.pi / 4 + 1e-9
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "savefig.dpi": 200, "savefig.bbox": "tight", "figure.facecolor": "#fcfcfb",
                     "axes.facecolor": "#fcfcfb",
                     # LaTeX look (Computer Modern) without needing a TeX install
                     "font.family": "serif", "font.serif": ["cmr10"], "mathtext.fontset": "cm",
                     "axes.formatter.use_mathtext": True, "axes.unicode_minus": False})


def load(lat, geo, R):
    with open(os.path.join(DATA, f"{lat}_{geo}_R{R:g}.csv")) as fh:
        rows = list(csv.DictReader(fh))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0] if k not in ("geometry", "lattice")}


def by_qh(c, key, sel=None):
    qh = np.unique(c["Q_over_h"])
    m = np.ones_like(c["theta"], bool) if sel is None else sel
    return qh, [c[key][m & np.isclose(c["Q_over_h"], q)] for q in qh]


def full_circle(c):
    """Points of the [-pi, pi] grid (81 points) only, so averages are not biased by the denser map grid."""
    th = np.unique(c["theta"])
    grid = np.linspace(-np.pi, np.pi, 81)
    keep = np.isin(np.round(c["theta"], 12), np.round(grid, 12)) | (c["theta"] == 0)
    keep &= np.abs(c["theta"]) < np.pi - 1e-9  # drop theta = pi (duplicate of -pi)
    return keep


# ---------------------------------------------------------------- velocities (Figs. 2, 3, S1, S2, S7)
def velocity_figure(geo, R, fname, title):
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.6))
    for lat in LATS:
        c = load(lat, geo, R)
        circ = full_circle(c)
        qh, U = by_qh(c, "U_perp", circ)
        Ubar = np.array([u.mean() for u in U])
        U0 = c["U_perp_isolated"][0]
        kw = dict(color=COLOR[lat], marker=MARK[lat], ms=5, lw=1.4, mfc="none", label=NAME[lat])
        axs[0].plot(qh, Ubar, **kw)
        axs[1].plot(qh, 100 * (Ubar - U0) / abs(U0), **kw)
        for q, u in zip(qh, U):
            th = c["theta"][circ & np.isclose(c["Q_over_h"], q)]
            o = np.argsort(th)
            axs[2].plot(th[o], u[o] / u.mean(), color=COLOR[lat], lw=0.8, alpha=0.35 + 0.65 * q / qh.max())
    axs[0].axhline(U0, color=MUTED, lw=0.8, ls="--")
    axs[0].annotate("isolated pair", (0.03, U0), xytext=(0, -12), textcoords="offset points", color=MUTED, fontsize=8)
    axs[0].set_ylabel(r"$\bar U_\perp$  [$\zeta V^S/(\mu Q)$]")
    axs[1].set_ylabel(r"lattice correction  $100(\bar U_\perp-U_{iso})/|U_{iso}|$  (%)")
    axs[2].set_ylabel(r"$U_\perp/\bar U_\perp$")
    axs[0].set_xlabel(r"$Q/h$"); axs[1].set_xlabel(r"$Q/h$"); axs[2].set_xlabel(r"$\theta$ (rad)")
    axs[2].set_xlim(-np.pi, np.pi)
    axs[0].legend(frameon=False, fontsize=8)
    for lat in LATS:
        axs[2].plot([], [], color=COLOR[lat], lw=1.5, label=NAME[lat])
    axs[2].legend(frameon=False, fontsize=8, ncol=2)
    axs[2].set_title(r"all $Q/h$ superimposed (darker = larger $Q/h$)", fontsize=8, color=MUTED)
    fig.suptitle(title, color=INK, fontsize=10, y=1.02)
    fig.savefig(os.path.join(OUT, fname))
    plt.close(fig)


# ---------------------------------------------------------------- stability maps (Figs. 4, 5, S4, S5, S8, S9)
def stability_figure(geo, R, lats, fname, title):
    fig, axs = plt.subplots(1, len(lats), figsize=(3.6 * len(lats) + 0.6, 3.6), sharey=True, squeeze=False)
    axs = axs[0]
    for ax, lat, letter in zip(axs, lats, "ABCD"):
        c = load(lat, geo, R)
        m = np.abs(c["theta"]) <= PI4
        th = np.unique(c["theta"][m]); qh = np.unique(c["Q_over_h"])
        Z = np.full((len(qh), len(th)), np.nan)
        Z[np.searchsorted(qh, c["Q_over_h"][m]), np.searchsorted(th, c["theta"][m])] = c["A_conf"][m]
        scale = np.nanmax(np.abs(Z), axis=1, keepdims=True)
        rel = Z / np.where(scale > 0, scale, 1)
        zero = scale[:, 0] < 1e-10  # identically zero rows (SCL with R = Q/2)
        col = np.where(rel < -1e-6, -1, np.where(rel > 1e-6, 1, 0)).astype(float)
        col[zero] = 0
        ax.pcolormesh(th, qh, col, cmap=matplotlib.colors.ListedColormap([STABLE, NEUTRAL, UNSTABLE]),
                      vmin=-1.5, vmax=1.5, shading="nearest")
        if not zero.all():
            ax.contour(th, qh, np.where(zero[:, None], np.nan, rel), levels=np.linspace(-0.9, 0.9, 7),
                       colors=INK, linewidths=0.4, alpha=0.6)
        ax.set_title(f"{letter}) {NAME[lat]}", color=INK, fontsize=10, loc="left")
        ax.set_xlabel(r"$\theta$ (rad)")
        if zero.all():
            ax.text(0, qh.mean(), r"$A_\perp \equiv 0$" + "\n(neutral)", ha="center", va="center", color=INK)
    axs[0].set_ylabel(r"$Q/h$")
    h1 = matplotlib.patches.Patch(color=STABLE, label=r"$A_\perp<0$: pushes Q back towards $a^{(1)}$")
    h2 = matplotlib.patches.Patch(color=UNSTABLE, label=r"$A_\perp>0$: pushes Q away from $a^{(1)}$")
    fig.legend(handles=[h1, h2], loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.1))
    fig.suptitle(title + ".  Lines: level sets of $A_\\perp$ normalized at each $Q/h$.", color=MUTED,
                 fontsize=8, y=-0.04)
    fig.savefig(os.path.join(OUT, fname))
    plt.close(fig)


# ---------------------------------------------------------------- minima of A_perp (Figs. 7, S6, S10)
def minima_figure(R, fname, title, geos=("bulk3d", "monolayer")):
    """Same definition as the notebook of the article: min of A_perp over theta in (0, pi/4]
    (counterclockwise) and over [-pi/4, 0) (clockwise), times nothing else (A_perp/h already)."""
    fig, axs = plt.subplots(1, len(geos), figsize=(5.4 * len(geos), 3.8), sharey=False, squeeze=False)
    for ax, geo in zip(axs[0], geos):
        for lat in LATS:
            if lat == "scl" and R == 0.5:
                continue
            c = load(lat, geo, R)
            for side, ls, fill in ((1, "-", "full"), (-1, "--", "none")):
                m = (np.sign(c["theta"]) == side) & (np.abs(c["theta"]) <= PI4)
                qh, A = by_qh(c, "A_conf", m)
                ax.plot(qh, [a.min() for a in A], ls, color=COLOR[lat], marker=MARK[lat], ms=5, fillstyle=fill, lw=1.3,
                        label=f"{NAME[lat]} {'counterclockwise' if side > 0 else 'clockwise'}")
        ax.set_yscale("symlog", linthresh=1e-8)
        ax.axhline(0, color=MUTED, lw=0.6)
        ax.set_xlabel(r"$Q/h$"); ax.set_ylabel(r"min$_\theta\,A_\perp/h$")
        ax.set_title(GEO[geo], color=INK, fontsize=10)
    axs[0][0].legend(frameon=False, fontsize=7, ncol=1)
    fig.suptitle(title, color=MUTED, fontsize=8, y=-0.04)
    fig.savefig(os.path.join(OUT, fname))
    plt.close(fig)


# ---------------------------------------------------------------- |A_perp| averaged over theta (Fig. S3)
def amplitude_figure(R, fname, geos=("monolayer", "bulk3d", "pure2d")):
    fig, axs = plt.subplots(1, len(geos), figsize=(4.4 * len(geos), 3.6), sharey=True)
    for ax, geo in zip(axs, geos):
        for lat in LATS:
            c = load(lat, geo, R)
            m = np.abs(c["theta"]) <= PI4
            qh, A = by_qh(c, "A_conf", m)
            y = np.array([np.abs(a).mean() for a in A])
            ax.semilogy(qh, np.maximum(y, 1e-22), color=COLOR[lat], marker=MARK[lat], ms=5, mfc="none", lw=1.3,
                        label=NAME[lat])
        ax.set_title(GEO[geo], color=INK, fontsize=10)
        ax.set_xlabel(r"$Q/h$")
    axs[0].set_ylabel(r"$\overline{|A_\perp|}/h$  (average over $|\theta|\leq\pi/4$)")
    axs[0].legend(frameon=False, fontsize=8)
    fig.suptitle(f"$R = {R:g}\\,Q$. For SCL with $R = Q/2$ the values are at round-off level "
                 "(exact cancellation, SI Section S3).", color=MUTED, fontsize=8, y=-0.04)
    fig.savefig(os.path.join(OUT, fname))
    plt.close(fig)


# ---------------------------------------------------------------- converged U_par, U_perp (commented SI figures)
def converged_velocity_figure(fname):
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.6))
    for lat in ("scl", "hex"):
        c = load(lat, "bulk3d", 0.5)
        circ = full_circle(c)
        for key, ls in (("U_perp", "-"), ("U_par", "--")):
            qh, U = by_qh(c, key, circ)
            sub = r"\perp" if key == "U_perp" else r"\parallel"
            lab = NAME[lat] + " $U_{" + sub + "}$"
            axs[0].plot(qh, [u.mean() for u in U], ls, color=COLOR[lat], marker=MARK[lat], ms=5, mfc="none", label=lab)
            m = circ & np.isclose(c["Q_over_h"], 0.25)
            o = np.argsort(c["theta"][m])
            axs[1].plot(c["theta"][m][o], c[key][m][o], ls, color=COLOR[lat], lw=1.5, label=lab)
    axs[0].set_xlabel(r"$Q/h$"); axs[0].set_ylabel(r"$\theta$-average  [$\zeta V^S/(\mu Q)$]")
    axs[1].set_xlabel(r"$\theta$ (rad)"); axs[1].set_ylabel(r"value at $Q/h = 0.25$")
    axs[0].legend(frameon=False, fontsize=8, ncol=2)
    fig.suptitle("Converged translational velocity of the pair, 3D, $R = Q/2$", color=MUTED, fontsize=8, y=-0.04)
    fig.savefig(os.path.join(OUT, fname))
    plt.close(fig)


if __name__ == "__main__":
    import time
    os.makedirs(OUT, exist_ok=True)
    jobs = [
        (velocity_figure, ("bulk3d", 0.5, "fig2_velocity_3D.png", "Fig. 2 (corrected): pair velocity, 3D lattices, $R = Q/2$")),
        (velocity_figure, ("monolayer", 0.5, "fig3_velocity_monolayer.png", "Fig. 3 (corrected): pair velocity, monolayers in a 3D fluid, $R = Q/2$")),
        (stability_figure, ("bulk3d", 0.5, ["bcl", "fcl", "hex"], "fig4_stability_3D.png", "Fig. 4 (corrected): 3D, $R = Q/2$ (SCL: $A_\\perp = 0$ identically)")),
        (stability_figure, ("monolayer", 0.5, ["bcl", "fcl", "hex"], "fig5_stability_monolayer.png", "Fig. 5 (corrected): monolayers in a 3D fluid, $R = Q/2$ (SCL: $A_\\perp = 0$ identically)")),
        (minima_figure, (0.5, "fig7_minima.png", "Fig. 7 (corrected): $R = Q/2$. Filled: $\\theta > 0$, open: $\\theta < 0$")),
        (velocity_figure, ("bulk3d", 0.25, "figS1_velocity_3D_R025.png", "Fig. S1 (corrected): 3D, $R = Q/4$")),
        (velocity_figure, ("monolayer", 0.25, "figS2_velocity_monolayer_R025.png", "Fig. S2 (corrected): monolayers, $R = Q/4$")),
        (amplitude_figure, (0.5, "figS3_mean_abs_A.png")),
        (stability_figure, ("bulk3d", 0.25, LATS, "figS4_stability_3D_R025.png", "Fig. S4 (corrected): 3D, $R = Q/4$")),
        (stability_figure, ("monolayer", 0.25, LATS, "figS5_stability_monolayer_R025.png", "Fig. S5 (corrected): monolayers in a 3D fluid, $R = Q/4$")),
        (minima_figure, (0.25, "figS6_minima_R025.png", "Fig. S6 (corrected): $R = Q/4$. Filled: $\\theta > 0$, open: $\\theta < 0$")),
        (velocity_figure, ("pure2d", 0.5, "figS7_velocity_pure2d.png", "Extra SI: pair velocity, pure 2D, $R = Q/2$")),
        (stability_figure, ("pure2d", 0.5, ["bcl", "fcl", "hex"], "figS8_stability_pure2d.png", "Extra SI: pure 2D (line forces), $R = Q/2$")),
        (stability_figure, ("pure2d", 0.25, LATS, "figS9_stability_pure2d_R025.png", "Extra SI: pure 2D, $R = Q/4$")),
        (minima_figure, (0.5, "figS10_minima_pure2d.png", "Extra SI: pure 2D, $R = Q/2$"), {"geos": ("pure2d",)}),
        (converged_velocity_figure, ("figS11_converged_velocity_3D.png",)),
    ]
    t0 = time.time()
    for i, (fn, args, *kw) in enumerate(jobs, 1):
        name = next(a for a in args if isinstance(a, str) and a.endswith(".png"))
        print(f"[{i:2d}/{len(jobs)}] {name} ...", end=" ", flush=True)
        t = time.time()
        fn(*args, **(kw[0] if kw else {}))
        print(f"ok ({time.time() - t:.1f} s)", flush=True)
    print(f"{len(jobs)} figures written to {OUT} in {time.time() - t0:.1f} s", flush=True)
