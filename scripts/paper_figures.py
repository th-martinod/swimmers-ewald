"""Every data figure of the article and its SI, recomputed with the approved conventions.

Input : results/data/{lattice}_{geometry}_R{R}.csv     (spinners-sweep configs/paper.toml)
Output: results/figures/figs/*.{eps,pdf}    drop-in replacements for the files of the main text
        results/figures/figs_2/*.{eps,pdf}  drop-in replacements for the files of the SI
        results/figures/preview/*.png       previews
docs/FIGURE_MAP.md lists the correspondence with the figures of the manuscript.

The stability maps leave the top of the canvas empty because the manuscript adds the panel letter
with overpic at (57, 98). No transparency is used (EPS).

Run: python3 scripts/paper_figures.py
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

DATA, OUT = "results/data", "results/figures"
LATS = ["scl", "bcl", "fcl", "hex"]
NAME = {"scl": "SCL", "bcl": "BCL", "fcl": "FCL", "hex": "HEX"}
# hues and marker shapes of the original figures (red square, green triangle, blue circle, brown diamond),
# adjusted to pass the colour-vision-deficiency checks: teal instead of green, ochre instead of brown
COLOR = {"scl": "#cf2a2a", "bcl": "#109a86", "fcl": "#2a4fd6", "hex": "#b58c00"}
MARK = {"scl": "s", "bcl": "^", "fcl": "o", "hex": "D"}
STABLE, UNSTABLE, INK, MUTED = "#2a5bd7", "#d93a32", "#1f1f1e", "#5f5e5a"
PI4 = np.pi / 4 + 1e-9

plt.rcParams.update({
    "font.size": 9, "font.family": "serif", "font.serif": ["cmr10"], "mathtext.fontset": "cm",
    "axes.formatter.use_mathtext": True, "axes.unicode_minus": False,
    "axes.edgecolor": INK, "axes.linewidth": 0.6, "axes.labelcolor": INK,
    "xtick.color": INK, "ytick.color": INK, "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "legend.frameon": False, "legend.fontsize": 8, "legend.handletextpad": 0.3,
    "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
    "ps.fonttype": 42, "pdf.fonttype": 42, "savefig.dpi": 300,
})


def save(fig, sub, name, tight=True):
    os.makedirs(os.path.join(OUT, sub), exist_ok=True)
    os.makedirs(os.path.join(OUT, "preview"), exist_ok=True)
    kw = {"bbox_inches": "tight", "pad_inches": 0.03} if tight else {}
    for ext in ("eps", "pdf"):
        fig.savefig(os.path.join(OUT, sub, f"{name}.{ext}"), **kw)
    fig.savefig(os.path.join(OUT, "preview", f"{sub}__{name}.png"), dpi=150, **kw)
    plt.close(fig)
    print(f"  {sub}/{name}")


def load(lat, geo, R):
    with open(os.path.join(DATA, f"{lat}_{geo}_R{R:g}.csv")) as fh:
        rows = list(csv.DictReader(fh))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0] if k not in ("geometry", "lattice")}


def on_grid(c, n, tmax):
    """Rows of the uniform grid linspace(-tmax, tmax, n), without the duplicate theta = pi."""
    grid = np.round(np.linspace(-tmax, tmax, n), 12)
    keep = np.isin(np.round(c["theta"], 12), grid)
    if tmax == np.pi:
        keep &= c["theta"] < np.pi - 1e-9
    return keep


def rows_by_qh(c, key, sel):
    qh = np.unique(c["Q_over_h"])
    out = []
    for q in qh:
        m = sel & np.isclose(c["Q_over_h"], q)
        o = np.argsort(c["theta"][m])
        out.append((c["theta"][m][o], c[key][m][o]))
    return qh, out


def style_sci(ax):
    ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3), useMathText=True)


# ------------------------------------------------------------------ velocities (Figs. 2, 3, S1, S2)
def velocity_figures(geo, R, lats, sub, suffix):
    tag = "3D" if geo == "bulk3d" else "2D"
    figL, axL = plt.subplots(figsize=(3.3, 2.75))
    figR, axR = plt.subplots(figsize=(3.0, 2.75))
    for lat in lats:
        c = load(lat, geo, R)
        qh, rows = rows_by_qh(c, "B_perp", on_grid(c, 81, np.pi))
        Bbar = np.array([b.mean() for _, b in rows])
        axL.plot(qh, Bbar, ls="none", marker=MARK[lat], ms=4.5, mfc="white", mec=COLOR[lat], mew=0.9,
                 label=NAME[lat])
        for (th, b), bb in zip(rows, Bbar):
            axR.plot(th, b / bb, ls="none", marker=MARK[lat], ms=2.6, mfc="white", mec=COLOR[lat], mew=0.5)
        axR.plot([], [], ls="none", marker=MARK[lat], ms=4.5, mfc="white", mec=COLOR[lat], mew=0.9,
                 label=NAME[lat])
    axL.axhline(0, color=MUTED, lw=0.4)
    axL.set_xlabel(r"$Q/h$")
    axL.set_ylabel(r"$\overline{B}_\perp/h$")
    axL.set_xlim(0, 0.26)
    style_sci(axL)
    axL.legend(loc="best")
    axR.set_xlabel(r"$\theta$")
    axR.set_ylabel(r"$B_\perp/\overline{B}_\perp$")
    axR.set_xlim(-np.pi - 0.1, np.pi + 0.1)
    axR.set_xticks(np.arange(-3, 4))
    axR.axhline(1, color=MUTED, lw=0.4)
    axR.legend(loc="best", ncol=2 if len(lats) > 2 else 1)
    save(figL, sub, f"CM_av_all_{tag}{suffix}")
    save(figR, sub, f"CM_vs_ang_all_{tag}{suffix}")


# ------------------------------------------------------------------ stability maps (Figs. 4, 5, S4, S5)
def stability_map(lat, geo, R, sub, suffix):
    tag = "3D" if geo == "bulk3d" else "2D"
    c = load(lat, geo, R)
    sel = on_grid(c, 101, np.pi / 4) & (np.abs(c["theta"]) > 1e-12)   # theta = 0 is the trivial zero
    qh = np.unique(c["Q_over_h"])
    th = np.unique(c["theta"][sel])
    Z = np.empty((len(qh), len(th)))
    for i, q in enumerate(qh):
        m = sel & np.isclose(c["Q_over_h"], q)
        o = np.argsort(c["theta"][m])
        Z[i] = c["A_perp"][m][o]
    fig = plt.figure(figsize=(3.3, 3.2))
    ax = fig.add_axes([0.2, 0.14, 0.75, 0.76])      # top 8 % free for the overpic panel letter
    big = np.abs(Z).max() * 10
    ax.contourf(th, qh, Z, levels=[-big, 0, big], colors=[STABLE, UNSTABLE])
    ax.contour(th, qh, Z, levels=[0], colors=INK, linewidths=0.6)
    ax.set_xlim(-np.pi / 4, np.pi / 4)
    ax.set_ylim(qh.min(), qh.max())
    ax.set_xticks([-0.5, 0, 0.5])
    ax.set_xlabel(r"$\theta$")
    ax.set_ylabel(r"$Q/h$")
    for s in ax.spines.values():
        s.set_zorder(5)
    save(fig, sub, f"stability_{NAME[lat]}_{tag}{suffix}", tight=False)


# ------------------------------------------------------------------ minima of A_perp (Figs. 7, S6)
def minima(c, side):
    """min over theta > 0 (counterclockwise, side = +1) or theta < 0 (clockwise, side = -1) of A_perp,
    refined with a parabola through the grid minimum. NaN when there is no stable orientation."""
    sel = on_grid(c, 101, np.pi / 4) & (side * c["theta"] > 1e-12)
    qh, rows = rows_by_qh(c, "A_perp", sel)
    out = []
    for th, a in rows:
        j = int(np.argmin(a))
        v = a[j]
        if 0 < j < len(a) - 1:
            y0, y1, y2 = a[j - 1], a[j], a[j + 1]
            den = y0 - 2 * y1 + y2
            if den > 0:
                v = y1 - (y2 - y0) ** 2 / (8 * den)
        out.append(v if v < 0 else np.nan)
    return qh, np.array(out)


MIN_SERIES = [("hex", "monolayer", "HEX 2D", "#cf2a2a", "s"), ("hex", "bulk3d", "HEX 3D", "#2a4fd6", "o"),
              ("bcl", "bulk3d", "BCL 3D", "#109a86", "^")]


def minima_figure(R, sub, suffix, size):
    fig, ax = plt.subplots(figsize=size)
    for lat, geo, lab, col, mk in MIN_SERIES:
        c = load(lat, geo, R)
        for side, word, face, ms in ((+1, "counterclockwise", col, 5.5), (-1, "clockwise", "white", 4.0)):
            qh, m = minima(c, side)
            ax.plot(qh, m, ls="none", marker=mk, ms=ms, mfc=face, mec=col, mew=0.9, label=f"{lab} {word}")
    ax.axhline(0, color=MUTED, lw=0.4)
    ax.set_xlabel(r"$Q/h$")
    ax.set_ylabel(r"$\min A_\perp/h$")
    ax.set_xlim(0, 0.26)
    style_sci(ax)
    ax.legend(loc="lower left")
    save(fig, sub, f"conf_min_all{suffix}")


# ------------------------------------------------------------------ magnitude of A_perp (Fig. S3)
def magnitude_figure(geo, R=0.5):
    tag = "3D" if geo == "bulk3d" else "2D"
    fig, ax = plt.subplots(figsize=(3.2, 2.75))
    for lat in LATS:
        c = load(lat, geo, R)
        qh, rows = rows_by_qh(c, "A_perp", on_grid(c, 101, np.pi / 4))
        m = np.array([np.mean(np.abs(a)) for _, a in rows])
        m = np.maximum(m, 1e-22)                          # exact zeros of the SCL cancellation
        ax.semilogy(qh, m, ls="none", marker=MARK[lat], ms=4.5, mfc="white", mec=COLOR[lat], mew=0.9,
                    label=NAME[lat])
    ax.set_xlabel(r"$Q/h$")
    ax.set_ylabel(r"$\overline{|A_\perp|}/h$")
    ax.set_xlim(0, 0.26)
    ax.legend(loc="center right")
    save(fig, "figs_2", f"magstability_{tag}")


def main():
    print("main text (R = Q/2):")
    for geo in ("bulk3d", "monolayer"):
        velocity_figures(geo, 0.5, LATS, "figs", "")
        for lat in ("bcl", "fcl", "hex"):
            stability_map(lat, geo, 0.5, "figs", "")
    minima_figure(0.5, "figs", "", (6.2, 3.6))
    print("SI (R = Q/4 unless stated):")
    velocity_figures("bulk3d", 0.25, ["bcl", "hex"], "figs_2", "_R_p25")
    velocity_figures("monolayer", 0.25, ["hex"], "figs_2", "_R_p25")
    for geo in ("monolayer", "bulk3d"):
        magnitude_figure(geo)
    for lat in ("hex", "bcl"):
        stability_map(lat, "bulk3d", 0.25, "figs_2", "_R_p25")
    stability_map("hex", "monolayer", 0.25, "figs_2", "_R_p25")
    minima_figure(0.25, "figs_2", "_R_p25", (3.6, 3.0))


if __name__ == "__main__":
    main()
