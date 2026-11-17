#!/usr/bin/env python
"""
make_schematic_and_abstract.py
================================
Standalone script producing two schematic (non-data) figures for the
Rad4-DNA recognition ABFE paper:

  fig_thermodynamic_cycle.png  - the double-decoupling thermodynamic
                                  cycle used to define ddF^i = dG^i_bound
                                  - dG^i_docked for a structural motif i
                                  (illustrated with the partner bases,
                                  the largest contributor in Table 1).
  fig_graphical_abstract.png   - a TOC-style graphical abstract summarising
                                  docked -> bound lesion recognition and
                                  the ranked |ddF| of the four motifs.

Pure schematic/cartoon art (matplotlib patches only) -- no simulation
data is read or required. Okabe-Ito colorblind-safe palette, consistent
with the rest of the project:
    docked/reference   #0072B2 (blue)
    bound/productive   #E69F00 (orange)  -- also used for partner bases
    BHD2 hairpin       #009E73 (green)
    BHD3 hairpin       #D55E00 (vermillion) -- also used for the lesion
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, FancyArrowPatch, FancyBboxPatch
from matplotlib.path import Path
import matplotlib.patches as mpatches

# ---- palette -------------------------------------------------------------
C_DOCK = "#0072B2"
C_BOUND = "#E69F00"
C_BHD2 = "#009E73"
C_BHD3 = "#D55E00"
C_LESION = "#CC79A7"
C_GRAY = "#555555"
C_PROT = "#B0B7C6"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "mathtext.fontset": "dejavusans",
})


# ===========================================================================
# Reusable cartoon primitives
# ===========================================================================
def _backbone_xy(x0, x1, y, n=40, wobble=0.0, seed=0):
    x = np.linspace(x0, x1, n)
    rng = np.random.default_rng(seed)
    y_arr = np.full_like(x, y) + (wobble * np.sin((x - x0) * 2.4 + seed) * 0.15)
    return x, y_arr


def draw_protein_blob(ax, cx, cy, s, label="Rad4", engaged=False, alpha=1.0):
    """A simple irregular protein blob (rounded polygon) with a label."""
    n = 10
    rng = np.random.default_rng(7)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    radii = s * (1.0 + 0.18 * rng.standard_normal(n))
    xs = cx + radii * np.cos(angles) * (1.35 if engaged else 1.0)
    ys = cy + radii * np.sin(angles) * (0.85 if engaged else 1.0)
    xs = np.append(xs, xs[0]); ys = np.append(ys, ys[0])
    ax.fill(xs, ys, color=C_PROT, alpha=alpha * 0.9, zorder=1,
            edgecolor="#7A8296", linewidth=1.2)
    ax.text(cx, cy, label, ha="center", va="center", fontsize=9,
            style="italic", color="#3A3F4D", zorder=5, alpha=alpha)


def draw_hairpin(ax, base_xy, tip_xy, color, ghost=False, lw=3.0, label=None):
    """A small loop/hook representing a BHD2/BHD3 hairpin tip reaching
    from the protein body toward the duplex gap."""
    (x0, y0), (x1, y1) = base_xy, tip_xy
    mx, my = (x0 + x1) / 2, (y0 + y1) / 2 + 0.35 * abs(x1 - x0)
    verts = [(x0, y0), (mx, my), (x1, y1)]
    codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3]
    path = Path(verts, codes)
    alpha = 0.28 if ghost else 0.95
    ls = "--" if ghost else "-"
    patch = mpatches.PathPatch(path, facecolor="none", edgecolor=color,
                                lw=lw, alpha=alpha, linestyle=ls, zorder=4)
    ax.add_patch(patch)
    ax.plot([x1], [y1], "o", color=color, ms=6, alpha=alpha, zorder=4)
    if label:
        ax.text(x1, y1 + (0.22 if y1 > my else -0.28), label, color=color,
                fontsize=8, fontweight="bold", ha="center", alpha=min(1, alpha + 0.2))


def draw_duplex_closed(ax, ox, oy, s, ghost_motif=False):
    """Docked-state cartoon: an intact duplex, lesion + 2 partner bases
    intrahelical (drawn as coloured rungs), Rad4 not yet engaged."""
    x0, x1 = ox - 2.6 * s, ox + 2.6 * s
    yt = oy + 0.55 * s
    yb = oy - 0.55 * s
    xt, yt_arr = _backbone_xy(x0, x1, yt, wobble=s, seed=1)
    xb, yb_arr = _backbone_xy(x0, x1, yb, wobble=s, seed=2)
    ax.plot(xt, yt_arr, color=C_DOCK, lw=2.2, zorder=2)
    ax.plot(xb, yb_arr, color=C_DOCK, lw=2.2, zorder=2)

    rung_x = np.linspace(x0 + 0.35 * s, x1 - 0.35 * s, 9)
    special = {3: (C_LESION, "T64"), 4: (C_BOUND, "585"), 5: (C_BOUND, "584")}
    for i, rx in enumerate(rung_x):
        yt_i = yt + wob(rx, x0, s, 1)
        yb_i = yb + wob(rx, x0, s, 2)
        if i in special:
            color, tag = special[i]
            ghost = ghost_motif and (i in (4, 5))
            alpha = 0.25 if ghost else 1.0
            ls = "--" if ghost else "-"
            ax.plot([rx, rx], [yb_i, yt_i], color=color, lw=3.2, alpha=alpha,
                     linestyle=ls, zorder=3)
            ax.plot([rx], [(yt_i + yb_i) / 2], "o", color=color, ms=5,
                     alpha=alpha, zorder=3)
        else:
            ax.plot([rx, rx], [yb_i, yt_i], color="#9AA3B5", lw=1.6, zorder=2)

    draw_protein_blob(ax, ox - 0.3 * s, oy - 2.0 * s, 0.85 * s, engaged=False)
    ax.text(ox, oy + 1.15 * s, "DOCKED", ha="center", fontsize=11,
            fontweight="bold", color=C_DOCK)


def wob(x, x0, s, seed):
    return s * np.sin((x - x0) * 2.4 + seed) * 0.15


def draw_duplex_open(ax, ox, oy, s, ghost_motif=False):
    """Bound-state cartoon: duplex opened at the lesion window, BHD2/BHD3
    hairpins inserted from Rad4, lesion + partner bases flipped extra-
    helical into the protein pocket."""
    gap = 0.9 * s
    x0, x1 = ox - 2.6 * s, ox + 2.6 * s
    yt, yb = oy + 0.55 * s, oy - 0.55 * s

    for side in (-1, 1):
        xa = ox + side * gap
        xe = x0 if side < 0 else x1
        xs = np.linspace(xa, xe, 20)
        bow = 0.42 * s * (1 - abs(xs - xa) / abs(xe - xa))
        ax.plot(xs, yt + bow, color=C_BOUND, lw=2.2, zorder=2)
        ax.plot(xs, yb - bow, color=C_BOUND, lw=2.2, zorder=2)
        rung_x = np.linspace(xa + side * 0.3 * s, xe - side * 0.3 * s, 4)[
            (slice(None) if side < 0 else slice(None))
        ]
        for rx in rung_x:
            frac = 1 - abs(rx - xa) / abs(xe - xa)
            b = 0.42 * s * frac
            ax.plot([rx, rx], [yb - b, yt + b], color="#C9B27A", lw=1.6, zorder=2)

    draw_protein_blob(ax, ox, oy - 1.55 * s, 1.15 * s, engaged=True)

    ghost = ghost_motif
    alpha = 0.25 if ghost else 1.0
    ls = "--" if ghost else "-"
    ecol = "none" if ghost else None
    for (dx, color, tag) in ((-0.32 * s, C_BOUND, "585"), (0.0, C_LESION, "T64"),
                              (0.32 * s, C_BOUND, "584")):
        ex, ey = ox + dx, oy - 0.55 * s
        e = Ellipse((ex, ey), 0.45 * s, 0.28 * s, facecolor=color, alpha=alpha,
                    edgecolor=color if not ghost else color, lw=1.5,
                    linestyle=ls, zorder=4)
        ax.add_patch(e)

    draw_hairpin(ax, (ox - 0.75 * s, oy - 2.2 * s), (ox - 0.55 * s, oy - 0.05 * s),
                 C_BHD2, ghost=False, label="BHD2")
    draw_hairpin(ax, (ox + 0.75 * s, oy - 2.2 * s), (ox + 0.55 * s, oy - 0.05 * s),
                 C_BHD3, ghost=False, label="BHD3")

    ax.text(ox, oy + 1.15 * s, "BOUND", ha="center", fontsize=11,
            fontweight="bold", color=C_BOUND)


# ===========================================================================
# Deliverable 1 -- thermodynamic (double-decoupling) cycle
# ===========================================================================
def make_thermo_cycle(outpath):
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.set_xlim(-1, 11)
    ax.set_ylim(0.3, 8.9)
    ax.axis("off")

    s = 0.62
    TL, TR = (2.0, 7.3), (8.0, 7.3)
    BL, BR = (2.0, 2.3), (8.0, 2.3)

    draw_duplex_closed(ax, *TL, s, ghost_motif=False)
    draw_duplex_closed(ax, *TR, s, ghost_motif=True)
    draw_duplex_open(ax, *BL, s, ghost_motif=False)
    draw_duplex_open(ax, *BR, s, ghost_motif=True)

    def arrow(p0, p1, color, lw=2.4, style="-|>"):
        a = FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=18,
                             lw=lw, color=color, zorder=6,
                             connectionstyle="arc3,rad=0.0")
        ax.add_patch(a)

    arrow((TL[0] + 1.9 * s, TL[1]), (TR[0] - 1.9 * s, TR[1]), C_DOCK)
    ax.text((TL[0] + TR[0]) / 2, TL[1] + 0.55, r"$\Delta G^{i}_{\mathrm{docked}}$",
            ha="center", fontsize=12, color=C_DOCK, fontweight="bold")
    ax.text((TL[0] + TR[0]) / 2, TL[1] + 0.15, r"$\lambda:\,0\rightarrow1$",
            ha="center", fontsize=9, color=C_GRAY)

    arrow((BL[0] + 1.9 * s, BL[1]), (BR[0] - 1.9 * s, BR[1]), C_BOUND)
    ax.text((BL[0] + BR[0]) / 2, BL[1] + 0.55, r"$\Delta G^{i}_{\mathrm{bound}}$",
            ha="center", fontsize=12, color=C_BOUND, fontweight="bold")
    ax.text((BL[0] + BR[0]) / 2, BL[1] + 0.15, r"$\lambda:\,0\rightarrow1$",
            ha="center", fontsize=9, color=C_GRAY)

    arrow((TL[0], TL[1] - 1.75 * s), (BL[0], BL[1] + 1.75 * s), C_GRAY, lw=2.0,
          style="-|>")
    ax.text(TL[0] - 1.35, (TL[1] + BL[1]) / 2, "docked $\\rightarrow$ bound",
            rotation=90, ha="center", va="center", fontsize=10, color=C_GRAY)
    ax.text(TL[0] - 1.9, (TL[1] + BL[1]) / 2 - 1.1,
            "(real transition;\nnot directly\nsimulated here)",
            rotation=90, ha="center", va="center", fontsize=7.5, color=C_GRAY,
            style="italic")

    arrow((TR[0], TR[1] - 1.75 * s), (BR[0], BR[1] + 1.75 * s), C_GRAY, lw=1.6,
          style="-")
    ax.text(TR[0] + 1.05, (TR[1] + BR[1]) / 2, "decoupled\nreference\nstates",
            rotation=90, ha="center", va="center", fontsize=8.5, color=C_GRAY)

    ax.text(5.0, 4.75, r"$\Delta\Delta F^{i} \;=\; \Delta G^{i}_{\mathrm{bound}} "
                        r"- \Delta G^{i}_{\mathrm{docked}}$",
            ha="center", va="center", fontsize=13, color="#222222",
            fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                      edgecolor="#999999", lw=1.0))
    ax.text(5.0, 4.25, "motif $i$ shown here: partner bases   "
                       r"($|\Delta\Delta F| = +193.5$ kcal/mol)",
            ha="center", va="center", fontsize=9.5, color=C_GRAY)

    fig.tight_layout()
    fig.savefig(outpath, dpi=300, facecolor="white")
    plt.close(fig)


# ===========================================================================
# Deliverable 2 -- graphical abstract / TOC graphic
# ===========================================================================
def make_graphical_abstract(outpath):
    fig, ax = plt.subplots(figsize=(9, 5.2))
    ax.set_xlim(-1, 15)
    ax.set_ylim(-1, 9)
    ax.axis("off")

    s = 0.85
    L, R = (3.3, 4.2), (11.3, 4.2)
    draw_duplex_closed(ax, *L, s, ghost_motif=False)
    draw_duplex_open(ax, *R, s, ghost_motif=False)

    arr = FancyArrowPatch((L[0] + 2.6, L[1] + 0.3), (R[0] - 2.6, R[1] + 0.3),
                           connectionstyle="arc3,rad=-0.35", arrowstyle="-|>",
                           mutation_scale=30, lw=4.5, color="#333333", zorder=7)
    ax.add_patch(arr)
    ax.text((L[0] + R[0]) / 2, 7.55, "Rad4 lesion recognition",
            ha="center", fontsize=14, fontweight="bold", color="#222222")

    callouts = [
        ("Partner bases", 193.5, C_BOUND, (R[0] - 0.3, R[1] - 0.5), (14.2, 6.6), 15, "bold"),
        ("BHD3 hairpin", 73.5, C_BHD3, (R[0] + 0.55, R[1] + 1.1), (14.2, 5.2), 11, "bold"),
        ("BHD2 hairpin", 68.2, C_BHD2, (R[0] - 0.55, R[1] + 1.1), (14.2, 3.9), 11, "bold"),
        ("6-4PP lesion", 2.8, C_LESION, (R[0], R[1] - 0.5), (14.2, 2.6), 9, "normal"),
    ]
    for label, val, color, xy, xytext, fs, fw in callouts:
        ax.annotate(f"{label}\n" + r"$|\Delta\Delta F|=$" + f"{val:g} kcal/mol",
                    xy=xy, xytext=xytext, ha="left", va="center", fontsize=fs,
                    fontweight=fw, color=color,
                    arrowprops=dict(arrowstyle="-", color=color, lw=1.3,
                                     shrinkA=2, shrinkB=2, alpha=0.8))

    ax.text(7.3, 0.35,
            "Alchemical decoupling resolves recognition free energy by structural motif",
            ha="center", fontsize=10, style="italic", color="#444444")

    fig.tight_layout()
    fig.savefig(outpath, dpi=300, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    import os
    outdir = os.path.dirname(os.path.abspath(__file__))
    make_thermo_cycle(os.path.join(outdir, "fig_thermodynamic_cycle.png"))
    make_graphical_abstract(os.path.join(outdir, "fig_graphical_abstract.png"))
    print("wrote fig_thermodynamic_cycle.png and fig_graphical_abstract.png to", outdir)
