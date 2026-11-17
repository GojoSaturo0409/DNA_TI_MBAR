#!/usr/bin/env python3
"""TOC graphic for the revision: docked -> bound cartoon (original primitives) + signed motif ddF for 6-4PP and CPD."""
import os, sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import schematic_primitives as sp
A = "/media/rahul/EXTERNAL_USB/binding_affinity/TI/ti_rep1/state_2026-10-30/analysis"
DD = {"64pp": pd.read_csv(f"{A}/prod64v10/summary_64pp_ddG.csv").set_index("motif"),
      "cpd": pd.read_csv(f"{A}/cpd_v10/summary_cpd_ddG.csv").set_index("motif")}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.weight": "bold", "axes.labelweight": "bold", "mathtext.default": "bf"})
fig = plt.figure(figsize=(9.0, 4.8))
ax = fig.add_axes([0.0, 0.0, 0.52, 1.0]); ax.set_xlim(-1, 15); ax.set_ylim(-0.5, 9.5); ax.axis("off")
s = 0.85; L, R = (3.0, 4.2), (11.0, 4.2)
sp.draw_duplex_closed(ax, *L, s, ghost_motif=False); sp.draw_duplex_open(ax, *R, s, ghost_motif=False)
ax.add_patch(FancyArrowPatch((L[0] + 2.4, L[1] + 0.4), (R[0] - 2.4, R[1] + 0.4), connectionstyle="arc3,rad=-0.35",
                             arrowstyle="-|>", mutation_scale=28, lw=4.5, color="#333333", zorder=7))
ax.text(7.0, 8.4, "Rad4 lesion recognition", ha="center", fontsize=16, color="#222222")
ax.text(7.0, 0.3, "decouple each motif in both states", ha="center", fontsize=12.5, color="#333333")
ax.text(7.0, -0.35, r"$\Delta\Delta F = \Delta F_{\mathrm{bound}} - \Delta F_{\mathrm{docked}}$", ha="center", fontsize=13)
b = fig.add_axes([0.6, 0.2, 0.38, 0.66])
MOT = ["lesion", "partner", "bhd2", "bhd3"]; x = np.arange(4); w = 0.38; CL = {"64pp": "#0072B2", "cpd": "#009E73"}
for k, (s_, lab) in enumerate((("64pp", "6-4PP"), ("cpd", "CPD"))):
    v = [DD[s_].loc[m, "ddG"] for m in MOT]; e = [DD[s_].loc[m, "ddG_err"] for m in MOT]
    b.bar(x + (k - 0.5) * (w + 0.03), v, w, yerr=e, capsize=3, color=CL[s_], label=lab, error_kw=dict(lw=1.5))
b.axhline(0, color="black", lw=1.5)
b.set_xticks(x, ["lesion", "partners", "BHD2", "BHD3"], fontsize=12); plt.setp(b.get_xticklabels(), rotation=30, ha="right", rotation_mode="anchor")
b.set_ylabel(r"$\Delta\Delta F$ (kcal/mol)", fontsize=13); b.tick_params(labelsize=11.5, width=1.6)
for sp_ in ("top", "right"): b.spines[sp_].set_visible(False)
for sp_ in ("left", "bottom"): b.spines[sp_].set_linewidth(1.6)
b.set_ylim(-34, 22); b.legend(loc="upper left", fontsize=12.5, frameon=False)
b.annotate("CPD: BHD2\nheld in the\ndocked state", xy=(2 + 0.5 * (w + 0.03), -30), xytext=(2.45, -21), fontsize=10.5, color="#009E73", ha="left",
           arrowprops=dict(arrowstyle="-", color="#009E73", lw=1.3))
b.set_title("motif free energies", fontsize=14, pad=8)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures", "toc_graphic")
fig.savefig(out + ".png", dpi=300, facecolor="white"); fig.savefig(out + ".pdf", facecolor="white"); print("wrote", out)
