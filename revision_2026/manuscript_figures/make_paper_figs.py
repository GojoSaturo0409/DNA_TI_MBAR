#!/usr/bin/env python3
"""Manuscript + SI figures for the revision (17 Nov 2026), drawn at print size with bold text and thick lines.
Run from anywhere:  python make_paper_figs.py
Inputs: same sources as analysis/slides/make_slide_figs_v10.py (prod64v10, cpd_v10, unbiased MD of 6-4PP and CPD),
plus the per-leg MBAR overlap matrices (*_E_overlap.npy, *_V_{A,B_dense,C}_overlap.npy) and the old 12-window overlap
JSONs (analysis/old12_overlap_*.json).
Outputs: ../figures/*.pdf (main text) and ../figures_si/*.pdf (SI)."""
import json, os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import gaussian_kde

TI = "/media/rahul/EXTERNAL_USB/binding_affinity/TI"
ST = f"{TI}/ti_rep1/state_2026-10-30"; A = f"{ST}/analysis"; DATA = f"{A}/slides/data"
HERE = os.path.dirname(os.path.abspath(__file__))
OUTM = os.path.join(HERE, "..", "figures"); OUTS = os.path.join(HERE, "..", "figures_si")
os.makedirs(OUTM, exist_ok=True); os.makedirs(OUTS, exist_ok=True)
MD = {"64pp": f"{TI}/rahul_ti_unbiased/analysis_out", "cpd": f"{TI}/cpd_unbiased/analysis_mirror/analysis_out"}
DIH = {"64pp": f"{DATA}/64pp/dihedrals", "cpd": f"{TI}/cpd_unbiased/analysis_mirror/intermediates/dihedrals"}
MSDD = {"64pp": f"{DATA}/64pp/com_timeseries", "cpd": f"{MD['cpd']}/com_timeseries"}
LESD = {"64pp": f"{DATA}/64pp/lesion", "cpd": f"{DATA}/cpd/lesion"}
SRC = {"64pp": f"{A}/prod64v10", "cpd": f"{A}/cpd_v10"}

C = {("64pp", "docked"): "#0072B2", ("64pp", "pb"): "#E69F00", ("cpd", "docked"): "#009E73", ("cpd", "pb"): "#CC79A7"}
CL = {"64pp": "#0072B2", "cpd": "#009E73"}
LAB = {"64pp": "6-4PP", "cpd": "CPD"}; SLAB = {"docked": "docked", "pb": "bound"}
MOT = ["lesion", "partner", "bhd2", "bhd3"]
MN = {"lesion": "Lesion", "partner": "Partner\nbases", "bhd2": "BHD2", "bhd3": "BHD3"}
MN1 = {"lesion": "lesion", "partner": "partner bases", "bhd2": "BHD2 tip", "bhd3": "BHD3 tip"}

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10.5, "font.weight": "bold", "axes.labelweight": "bold",
    "axes.titleweight": "bold", "axes.labelsize": 11.5, "axes.titlesize": 11.5, "xtick.labelsize": 10,
    "ytick.labelsize": 10, "legend.fontsize": 10, "axes.linewidth": 1.6, "xtick.major.width": 1.6,
    "ytick.major.width": 1.6, "xtick.major.size": 4.5, "ytick.major.size": 4.5, "xtick.direction": "out",
    "ytick.direction": "out", "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
    "mathtext.default": "bf", "savefig.bbox": "tight", "savefig.dpi": 300, "pdf.fonttype": 42})
LW = 2.3
def letter(ax, s, x=-0.12, y=1.04): ax.text(x, y, s, transform=ax.transAxes, fontsize=14, fontweight="bold", va="bottom", ha="right")
def save(fig, out, n):
    fig.savefig(f"{out}/{n}.pdf"); fig.savefig(f"{out}/{n}.png", dpi=200); plt.close(fig); print("wrote", n)
def sgn(v, d=1): return f"{v:+.{d}f}".replace("-", "−")
def legend4(fig, y=1.0, ncol=4):
    hs = [plt.Line2D([], [], color=C[(s, st)], lw=4) for s in ("64pp", "cpd") for st in ("docked", "pb")]
    fig.legend(hs, [f"{LAB[s]} {SLAB[st]}" for s in ("64pp", "cpd") for st in ("docked", "pb")], loc="lower center",
               ncol=ncol, bbox_to_anchor=(0.5, y), handlelength=1.8, columnspacing=1.4, fontsize=10.5)

LG = {s: pd.DataFrame(json.load(open(f"{d}/summary_{s}_legs.json"))).set_index("leg") for s, d in SRC.items()}
DD = {s: pd.read_csv(f"{d}/summary_{s}_ddG.csv").set_index("motif") for s, d in SRC.items()}
J = {s: {f"{st}_{m}": json.load(open(f"{SRC[s]}/{s}_{st}_{m}_2stage.json")) for st in ("pb", "docked") for m in MOT} for s in SRC}

# ======================= MAIN: Figure 4, ddG and its electrostatic / vdW split =======================
x = np.arange(len(MOT)); w = 0.38
def bars(ax, get, fs=9.5, pad=2.6):
    for k, s in enumerate(("64pp", "cpd")):
        v = np.array([get(s, m)[0] for m in MOT]); e = np.array([get(s, m)[1] for m in MOT]); xs = x + (k - 0.5) * (w + 0.03)
        ax.bar(xs, v, w, color=CL[s], yerr=e, capsize=3.5, error_kw=dict(lw=1.6, capthick=1.6), label=LAB[s])
        for xx, vv, ee in zip(xs, v, e):
            ax.text(xx, vv + (ee + pad) * (1 if vv >= 0 else -1), sgn(vv), ha="center", va="bottom" if vv >= 0 else "top",
                    fontsize=fs, fontweight="bold", rotation=90)
    ax.axhline(0, color="black", lw=1.4); ax.set_xticks(x, ["lesion", "partners", "BHD2", "BHD3"])
def stage(s, m, stg):
    p, q = LG[s].loc[f"pb_{m}"], LG[s].loc[f"docked_{m}"]; return p[f"{stg}_MBAR"] - q[f"{stg}_MBAR"], np.hypot(p[f"{stg}_err"], q[f"{stg}_err"])
fig, axs = plt.subplots(1, 3, figsize=(7.1, 3.3), sharey=True, gridspec_kw=dict(wspace=0.1))
bars(axs[0], lambda s, m: (DD[s].loc[m, "ddG"], DD[s].loc[m, "ddG_err"])); axs[0].set_title(r"total $\Delta\Delta F$")
bars(axs[1], lambda s, m: stage(s, m, "E"), fs=9); axs[1].set_title("electrostatic stage")
bars(axs[2], lambda s, m: stage(s, m, "V"), fs=9); axs[2].set_title("vdW stage")
axs[0].set_ylabel(r"$\Delta\Delta F$ (kcal/mol)"); axs[0].set_ylim(-58, 38)
for ax, l in zip(axs, "abc"): letter(ax, l); ax.tick_params(axis="x", labelsize=10); plt.setp(ax.get_xticklabels(), rotation=35, ha="right", rotation_mode="anchor")
axs[2].legend(loc="lower right", fontsize=10)
save(fig, OUTM, "fig_ddF")

# ======================= MAIN: theta, distance CV and MSD distributions =======================
def load_dih(s, st): return np.vstack([np.loadtxt(f"{DIH[s]}/{st}_rep{r}.dat")[:, 1:] for r in (1, 2, 3, 4)])
def circ_mean(a): r = np.radians(a); return np.degrees(np.arctan2(np.sin(r).mean(), np.cos(r).mean()))
def ckde(a, grid):
    t = np.concatenate([a - 360, a, a + 360])
    bw = 1.06 * min(np.std(a), np.subtract(*np.percentile(a, [75, 25])) / 1.34) * len(a) ** -0.2
    return 3 * gaussian_kde(t, bw_method=bw / np.std(t, ddof=1))(grid)
grid = np.linspace(-180, 180, 721)
fig, axs = plt.subplots(1, 3, figsize=(7.1, 2.6), gridspec_kw=dict(wspace=0.32))
for ax, (ci, ttl), lt in zip(axs, ((0, "5′ partner (584)"), (1, "3′ partner (585)"), (2, "lesion (602)")), "abc"):
    for s in ("64pp", "cpd"):
        for st in ("docked", "pb"):
            a = load_dih(s, st)[:, ci]; y = ckde(a, grid); c = C[(s, st)]
            ax.fill_between(grid, y, color=c, alpha=0.25 if s == "64pp" else 0.15, lw=0)
            ax.plot(grid, y, color=c, lw=LW, ls="-" if s == "64pp" else (0, (4, 1.5))); ax.axvline(circ_mean(a), color=c, ls=":", lw=1.8)
    ax.set_xlim(-180, 180); ax.set_xticks([-180, -90, 0, 90, 180]); ax.set_ylim(bottom=0); ax.set_xlabel(r"$\theta$ (deg)")
    ax.set_title(ttl); letter(ax, lt)
axs[0].set_ylabel("density"); legend4(fig, y=1.0); save(fig, OUTM, "fig_theta_kde")

CB = {s: {st: pd.read_csv(f"{MD[s]}/{st}/all_replicates_combined.csv") for st in ("docked", "pb")} for s in MD}
CVS = [("d_BHD2_A", r"$d_{\mathrm{BHD2}}$ (Å)"), ("d_BHD3_A", r"$d_{\mathrm{BHD3}}$ (Å)"), ("d_lesion_A", r"$d_{\mathrm{lesion}}$ (Å)"),
       ("delta_3prime_flip_A", r"$\gamma$, 5$'$ partner 584 (Å)"), ("gamma_5prime_flip_A", r"$\delta$, 3$'$ partner 585 (Å)")]
fig, axs = plt.subplots(2, 3, figsize=(7.1, 4.6), gridspec_kw=dict(wspace=0.32, hspace=0.55)); axs = axs.ravel()
for ax, (col, xl), lt in zip(axs, CVS, "abcde"):
    allv = np.concatenate([CB[s][st][col].values for s in MD for st in ("docked", "pb")]); g = np.linspace(allv.min() - 1.5, allv.max() + 1.5, 500)
    for s in ("64pp", "cpd"):
        for st in ("docked", "pb"):
            v = CB[s][st][col].values; y = gaussian_kde(v)(g); c = C[(s, st)]
            ax.fill_between(g, y, color=c, alpha=0.25 if s == "64pp" else 0.15, lw=0)
            ax.plot(g, y, color=c, lw=LW, ls="-" if s == "64pp" else (0, (4, 1.5))); ax.axvline(v.mean(), color=c, ls=":", lw=1.8)
    ax.set_xlabel(xl); ax.set_ylim(bottom=0); letter(ax, lt)
    if lt in "ad": ax.set_ylabel("density")
axs[5].axis("off")
for s in ("64pp", "cpd"):
    for st in ("docked", "pb"): axs[5].plot([], [], color=C[(s, st)], lw=4, ls="-" if s == "64pp" else (0, (4, 1.5)), label=f"{LAB[s]} {SLAB[st]}")
axs[5].legend(loc="center", fontsize=11.5, handlelength=2.4, labelspacing=0.9)
save(fig, OUTM, "fig_distance_kde")

MSDF = [("BHD2_BHD3_bhd2", "BHD2 tip"), ("BHD2_BHD3_bhd3", "BHD3 tip"), ("lesion602", "lesion (602)"),
        ("partner585_partner584_partner584", "5′ partner (584)"), ("partner585_partner584_partner585", "3′ partner (585)")]
def msdf(s, key, st): return f"{LESD[s] if key == 'lesion602' else MSDD[s]}/{key}_{st}_msd_tau_perrep.dat"
fig, axs = plt.subplots(2, 3, figsize=(7.1, 4.6), gridspec_kw=dict(wspace=0.35, hspace=0.5)); axs = axs.ravel()
for ax, (key, ttl), lt in zip(axs, MSDF, "abcde"):
    for s in ("64pp", "cpd"):
        for st in ("docked", "pb"):
            d = np.loadtxt(msdf(s, key, st)); m = d[:, 1:].mean(1); se = d[:, 1:].std(1, ddof=1) / 2; c = C[(s, st)]
            ax.plot(d[:, 0], m, color=c, lw=LW, ls="-" if s == "64pp" else (0, (4, 1.5))); ax.fill_between(d[:, 0], m - se, m + se, color=c, alpha=0.2, lw=0)
    ax.set_xlim(0, 25); ax.set_ylim(bottom=0); ax.set_xlabel(r"$\tau$ (ns)"); ax.set_title(ttl); letter(ax, lt)
    if lt in "ad": ax.set_ylabel(r"MSD (Å$^2$)")
axs[5].axis("off")
for s in ("64pp", "cpd"):
    for st in ("docked", "pb"): axs[5].plot([], [], color=C[(s, st)], lw=4, ls="-" if s == "64pp" else (0, (4, 1.5)), label=f"{LAB[s]} {SLAB[st]}")
axs[5].legend(loc="center", fontsize=11.5, handlelength=2.4, labelspacing=0.9)
save(fig, OUTM, "fig_msd")

# ======================= SI: full MBAR overlap matrices of the new runs =======================
CM = LinearSegmentedColormap.from_list("ov", ["#FFFFFF", "#FFF3B0", "#7FCDBB", "#1D91C0", "#0C2C84"])
def vblock(s, leg):
    mats = [np.load(f"{SRC[s]}/{s}_{leg}_V_{seg}_overlap.npy") for seg in ("A", "B_dense", "C")]
    n = sum(m.shape[0] for m in mats); M = np.full((n, n), np.nan); i = 0; edges = []
    for m in mats: k = m.shape[0]; M[i:i + k, i:i + k] = m; i += k; edges.append(i)
    return M, edges[:-1]
for s in ("64pp", "cpd"):
    fig, axs = plt.subplots(4, 4, figsize=(7.2, 7.6), gridspec_kw=dict(wspace=0.28, hspace=0.42))
    for r, m in enumerate(MOT):
        for c0, st in ((0, "pb"), (2, "docked")):
            leg = f"{st}_{m}"; E = np.load(f"{SRC[s]}/{s}_{leg}_E_overlap.npy"); V, ed = vblock(s, leg)
            for c, M, nm, adj in ((c0, E, "E", J[s][leg]["E"]["MBAR"]["adj_min"]), (c0 + 1, V, "V", J[s][leg]["V"]["MBAR"]["adj_min"])):
                ax = axs[r, c]; im = ax.imshow(np.ma.masked_invalid(M), cmap=CM, vmin=0, vmax=0.5, interpolation="nearest")
                ax.set_facecolor("#E6E6E6")
                if nm == "V":
                    for e in ed: ax.axhline(e - 0.5, color="black", lw=0.8); ax.axvline(e - 0.5, color="black", lw=0.8)
                ax.set_xticks([]); ax.set_yticks([]); [sp.set_visible(True) for sp in ax.spines.values()]
                ax.set_title(f"{SLAB[st]} {MN1[m]}\nstage {nm}, min {adj:.2f}", fontsize=8.6, pad=3)
    cax = fig.add_axes([0.93, 0.25, 0.018, 0.5]); cb = fig.colorbar(im, cax=cax); cb.set_label("overlap $O_{jk}$", fontsize=10.5)
    cb.ax.tick_params(labelsize=9.5)
    fig.suptitle(f"{LAB[s]}: MBAR overlap matrices of all eight legs", fontsize=12, fontweight="bold", y=0.955)
    save(fig, OUTS, f"overlap_matrices_{s}")

# original (withdrawn) 12-window runs: adjacent overlaps from overlap_report.txt (= data of the original Figure S1)
import re
_t = open(f"{TI}/overlap_report.txt").read(); _p = re.split(r"─{10,}\n(64pp_\w+)\n─{10,}", _t); OLDADJ = {}
for name, body in zip(_p[1::2], _p[2::2]):
    leg = name.replace("64pp_", ""); leg = leg if leg.startswith("docked") else "pb_" + leg
    OLDADJ[leg] = np.array([float(v) for v in re.findall(r"λ\d+→λ\d+: ([\d.]+)", body)])
legs8 = [f"{st}_{m}" for m in MOT for st in ("docked", "pb")]
fig, ax = plt.subplots(figsize=(7.1, 3.0)); xr = np.arange(len(legs8))
for i, l in enumerate(legs8):
    o = OLDADJ.get(l)
    if o is not None:
        ax.scatter(np.full(len(o), i - 0.27), np.clip(o, 1e-4, None), s=10, color="#888888"); ax.scatter(i - 0.27, max(o.min(), 1e-4), s=60, color="#555555", marker="o", label="original 12 windows" if i == 0 else None)
    for k, s in enumerate(("64pp", "cpd")):
        for stg, mk, dx in (("E", "s", 0.0), ("V", "^", 0.12)):
            a = np.array(J[s][l][stg]["MBAR"]["adj"] if stg == "E" else sum((J[s][l]["V"]["MBAR"]["segments"][g]["adj"] for g in ("A", "B_dense", "C")), []))
            xx = i + dx + 0.24 * k - 0.05
            ax.scatter(np.full(len(a), xx), a, s=8, color=CL[s], alpha=0.5)
            ax.scatter(xx, a.min(), s=55, color=CL[s], marker=mk, edgecolor="black", lw=0.8,
                       label=f"{LAB[s]} new, stage {stg}" if i == 0 else None)
ax.set_yscale("log"); ax.axhline(0.03, color="#C0392B", ls="--", lw=1.6); ax.set_ylim(8e-5, 1.2)
ax.set_xticks(xr, [l.replace("pb_", "bound\n").replace("docked_", "docked\n").replace("partner", "partner").replace("bhd", "BHD") for l in legs8], fontsize=9)
ax.set_ylabel("adjacent-window overlap"); ax.legend(fontsize=8.5, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0))
save(fig, OUTS, "overlap_old_vs_new")

# ======================= SI: <dU/dlambda> profiles, all 16 legs =======================
def pw(s, leg, stg):
    d = {}
    for p in J[s][leg][stg]["per_window"]: d[round(p["lam"], 6)] = (p["mean"], p["sem"])
    l = np.array(sorted(d)); return l, np.array([d[k][0] for k in l]), np.array([d[k][1] for k in l])
for stg, nm in (("E", "dudl_E"), ("V", "dudl_V")):
    fig, axs = plt.subplots(2, 2, figsize=(7.1, 5.0), gridspec_kw=dict(wspace=0.32, hspace=0.45)); axs = axs.ravel()
    for ax, m, lt in zip(axs, MOT, "abcd"):
        for s in ("64pp", "cpd"):
            for st in ("docked", "pb"):
                l, mm, e = pw(s, f"{st}_{m}", stg)
                ax.plot(l, mm, color=C[(s, st)], lw=2.0, ls="-" if s == "64pp" else (0, (4, 1.5)), marker="o", ms=2.5)
                ax.fill_between(l, mm - e, mm + e, color=C[(s, st)], alpha=0.25, lw=0)
        ax.axhline(0, color="black", lw=1); ax.set_xlabel(r"$\lambda$" if stg == "E" else r"$\lambda'$"); ax.set_title(MN1[m]); letter(ax, lt)
        if lt in "ac": ax.set_ylabel(r"$\langle\partial U/\partial\lambda\rangle$ (kcal/mol)")
    legend4(fig, y=0.99); save(fig, OUTS, nm)

# forward / reverse convergence of every leg
fig, axs = plt.subplots(2, 4, figsize=(7.2, 4.0), sharey=True, gridspec_kw=dict(wspace=0.15, hspace=0.55)); fr = [25, 50, 75, 100]
for k, l in enumerate(legs8):
    ax = axs[k % 2, k // 2]
    for s, mk in (("64pp", "o"), ("cpd", "s")):
        f = np.array(LG[s].loc[l, "fwd"]); r = np.array(LG[s].loc[l, "rev"])
        ax.plot(fr, f - f[-1], color=CL[s], marker=mk, ms=4, lw=2, label=f"{LAB[s]} forward")
        ax.plot(fr, r - r[-1], color=CL[s], marker=mk, ms=4, lw=1.6, ls=":", mfc="white", label=f"{LAB[s]} reverse")
    ax.axhline(0, color="black", lw=1); ax.set_ylim(-5, 5); ax.set_title(l.replace("pb_", "bound ").replace("docked_", "docked ").replace("bhd", "BHD"), fontsize=9.5)
    ax.set_xticks([25, 50, 75, 100]); ax.tick_params(labelsize=8.5)
    if k % 2 == 1: ax.set_xlabel("% of data", fontsize=10)
for a in axs[:, 0]: a.set_ylabel(r"$\Delta F-\Delta F_{\mathrm{final}}$", fontsize=10)
h, lb = axs[0, 0].get_legend_handles_labels(); fig.legend(h, lb, loc="lower center", ncol=4, bbox_to_anchor=(0.5, 0.99), fontsize=9.5)
save(fig, OUTS, "fwd_rev_convergence")

# TI vs MBAR
fig, ax = plt.subplots(figsize=(3.6, 3.4))
for s in ("64pp", "cpd"):
    a = DD[s]; ax.errorbar(a["ddG_TI"], a["ddG"], xerr=a["ddG_TI_err"], yerr=a["ddG_err"], fmt="o", color=CL[s], ms=6, capsize=3, lw=1.6, label=LAB[s])
ax.plot([-40, 25], [-40, 25], color="black", lw=1.2, ls=":"); ax.set_xlabel(r"$\Delta\Delta F$, TI (kcal/mol)"); ax.set_ylabel(r"$\Delta\Delta F$, MBAR (kcal/mol)")
ax.legend(loc="upper left"); save(fig, OUTS, "ti_vs_mbar")

# ======================= SI: unbiased MD (SASA, per-replicate RMSD / MSD) =======================
SAS = [("BHD2_tip_sasa_A2", "BHD2 tip"), ("BHD3_tip_sasa_A2", "BHD3 tip"), ("lesion_whole_sasa_A2", "lesion (whole)"),
       ("lesion_ring1_sasa_A2", "lesion ring 1"), ("lesion_ring2_sasa_A2", "lesion ring 2"),
       ("partner_584_sasa_A2", "5′ partner (584)"), ("partner_585_sasa_A2", "3′ partner (585)")]
fig, axs = plt.subplots(3, 3, figsize=(7.1, 7.8), gridspec_kw=dict(wspace=0.38, hspace=0.75)); axs = axs.ravel()
for ax, (col, ttl), lt in zip(axs, SAS, "abcdefg"):
    allv = np.concatenate([CB[s][st][col].values for s in MD for st in ("docked", "pb")]); g = np.linspace(allv.min() - 10, allv.max() + 10, 500)
    for s in ("64pp", "cpd"):
        for st in ("docked", "pb"):
            v = CB[s][st][col].values; c = C[(s, st)]; y = 1e3 * gaussian_kde(v)(g)
            ax.fill_between(g, y, color=c, alpha=0.22 if s == "64pp" else 0.13, lw=0); ax.plot(g, y, color=c, lw=2.0, ls="-" if s == "64pp" else (0, (4, 1.5)))
            ax.axvline(v.mean(), color=c, ls=":", lw=1.6)
    ax.set_title(ttl, pad=6); ax.set_xlabel(r"SASA (Å$^2$)"); ax.set_ylim(bottom=0); letter(ax, lt)
    if lt in "adg": ax.set_ylabel(r"density ($10^{-3}$ Å$^{-2}$)", fontsize=10)
for a in axs[7:]: a.axis("off")
for s in ("64pp", "cpd"):
    for st in ("docked", "pb"): axs[7].plot([], [], color=C[(s, st)], lw=4, ls="-" if s == "64pp" else (0, (4, 1.5)), label=f"{LAB[s]} {SLAB[st]}")
axs[7].legend(loc="center left", fontsize=11, handlelength=2.4, labelspacing=0.8, bbox_to_anchor=(0.1, 0.5))
save(fig, OUTS, "sasa_kde")

fig, axs = plt.subplots(7, 4, figsize=(7.2, 9.2), sharex=True, gridspec_kw=dict(wspace=0.3, hspace=0.35))
for r, (col, ttl) in enumerate(SAS):
    for c, (s, st) in enumerate([("64pp", "docked"), ("64pp", "pb"), ("cpd", "docked"), ("cpd", "pb")]):
        ax = axs[r, c]; df = CB[s][st]
        for rep, al in zip(sorted(df.replicate.unique()), (1.0, 0.75, 0.55, 0.4)):
            q = df[df.replicate == rep]; ax.plot(q.time_ns, q[col], color=C[(s, st)], lw=1.0, alpha=al)
        ax.axhline(df[col].mean(), color="black", ls="--", lw=1.3); ax.tick_params(labelsize=8)
        if r == 0: ax.set_title(f"{LAB[s]} {SLAB[st]}", fontsize=10, color=C[(s, st)])
        if c == 0: ax.set_ylabel(ttl.replace(" (", "\n("), fontsize=8.5)
for ax in axs[-1]: ax.set_xlabel("time (ns)", fontsize=9.5)
fig.text(0.0, 0.5, r"SASA (Å$^2$)", rotation=90, va="center", ha="right", fontsize=11, fontweight="bold")
save(fig, OUTS, "sasa_timeseries")

RM = [("BHD2_BHD3_bhd2", "BHD2 tip"), ("BHD2_BHD3_bhd3", "BHD3 tip"),
      ("partner585_partner584_partner584", "5′ partner (584)"), ("partner585_partner584_partner585", "3′ partner (585)")]
for kind, files, ttlk, yl, xl in (("rmsd", RM, "rmsd_t", r"COM displacement (Å)", "time (ns)"),
                                 ("msd", MSDF, "msd_tau", r"MSD (Å$^2$)", r"$\tau$ (ns)")):
    fig, axs = plt.subplots(len(files), 4, figsize=(7.2, 1.55 * len(files) + 0.6), sharex=True, gridspec_kw=dict(wspace=0.3, hspace=0.4))
    for r, (key, ttl) in enumerate(files):
        for c, (s, st) in enumerate([("64pp", "docked"), ("64pp", "pb"), ("cpd", "docked"), ("cpd", "pb")]):
            ax = axs[r, c]
            f = (f"{LESD[s]}/{key}_{st}_msd_tau_perrep.dat" if key == "lesion602" else f"{MSDD[s]}/{key}_{st}_{ttlk}_perrep.dat")
            d = np.loadtxt(f)
            for i, al in zip(range(d.shape[1] - 1), (1.0, 0.75, 0.55, 0.4)): ax.plot(d[:, 0], d[:, i + 1], color=C[(s, st)], lw=1.3, alpha=al)
            ax.tick_params(labelsize=8); ax.set_ylim(bottom=0)
            if r == 0: ax.set_title(f"{LAB[s]} {SLAB[st]}", fontsize=10, color=C[(s, st)])
            if c == 0: ax.set_ylabel(ttl.replace(" (", "\n("), fontsize=8.5)
    for ax in axs[-1]: ax.set_xlabel(xl, fontsize=9.5)
    fig.text(0.0, 0.5, yl, rotation=90, va="center", ha="right", fontsize=11, fontweight="bold")
    save(fig, OUTS, f"{kind}_perreplicate")

# ======================= SI: time series of theta and distance CVs (all replicas) =======================
TSP = [("64pp", "docked"), ("64pp", "pb"), ("cpd", "docked"), ("cpd", "pb")]
def ts_fig(series, name, ylab_rows):
    fig, axs = plt.subplots(len(series), 4, figsize=(7.2, 1.55 * len(series) + 0.7), sharex=True, gridspec_kw=dict(wspace=0.3, hspace=0.42))
    axs = np.atleast_2d(axs)
    for r, (getter, ttl) in enumerate(series):
        for c, (s, st) in enumerate(TSP):
            ax = axs[r, c]
            for rep, al in zip((1, 2, 3, 4), (1.0, 0.75, 0.55, 0.4)):
                t, y = getter(s, st, rep); ax.plot(t, y, color=C[(s, st)], lw=1.1, alpha=al)
            ax.tick_params(labelsize=8)
            if r == 0: ax.set_title(f"{LAB[s]} {SLAB[st]}", fontsize=10, color=C[(s, st)])
            if c == 0: ax.set_ylabel(ttl, fontsize=8.5)
    for ax in axs[-1]: ax.set_xlabel("time (ns)", fontsize=9.5)
    save(fig, OUTS, name)
def dih_get(ci):
    def g(s, st, rep):
        d = np.loadtxt(f"{DIH[s]}/{st}_rep{rep}.dat"); t = np.arange(len(d)) * 50.0 / len(d)
        mu = circ_mean(load_dih(s, st)[:, ci]); y = (d[:, 1 + ci] - mu + 180) % 360 - 180 + mu   # unwrap around the state mean
        return t, y
    return g
ts_fig([(dih_get(0), "θ 5′ partner\n584 (deg)"), (dih_get(1), "θ 3′ partner\n585 (deg)"), (dih_get(2), "θ lesion\n(deg)")], "theta_timeseries", None)
def cv_get(col):
    def g(s, st, rep):
        df = CB[s][st]; q = df[df.replicate == sorted(df.replicate.unique())[rep - 1]]; return q.time_ns.values, q[col].values
    return g
ts_fig([(cv_get("d_BHD2_A"), "d$_{BHD2}$ (Å)"), (cv_get("d_BHD3_A"), "d$_{BHD3}$ (Å)"), (cv_get("d_lesion_A"), "d$_{lesion}$ (Å)"),
        (cv_get("delta_3prime_flip_A"), "γ, 5′ partner\n584 (Å)"), (cv_get("gamma_5prime_flip_A"), "δ, 3′ partner\n585 (Å)")], "distance_timeseries", None)

# ======================= SI: CPD vs 6-4PP discrimination by motif =======================
MU = f"{TI}/mech_unbiased"
E_ = pd.read_csv(f"{MU}/tables/discrim_effect_sizes.csv"); Cc = pd.read_csv(f"{MU}/tables/discrim_classification.csv")
DM = ["BHD2", "partner bases", "lesion", "BHD3", "DNA site"]
DCOL = {"BHD2": "#CC79A7", "partner bases": "#E69F00", "lesion": "#0072B2", "BHD3": "#009E73", "DNA site": "#888888"}
rng = np.random.default_rng(1)
fig, axs = plt.subplots(2, 2, figsize=(7.2, 6.0), gridspec_kw=dict(wspace=0.3, hspace=0.5))
for j, st in enumerate(("pb", "docked")):
    ax = axs[0, j]; e = E_[E_.state == st]
    for k, m in enumerate(DM):
        v = e[e.motif == m].d_replica.abs().values; ax.scatter(np.full(len(v), k) + rng.uniform(-0.12, 0.12, len(v)), v, s=30, color=DCOL[m], zorder=3)
        ax.plot([k - 0.28, k + 0.28], [np.median(v)] * 2, color="black", lw=2.2)
    ax.set_xticks(range(len(DM)), ["BHD2", "partners", "lesion", "BHD3", "DNA"], fontsize=9.5)
    plt.setp(ax.get_xticklabels(), rotation=35, ha="right", rotation_mode="anchor")
    ax.set_title(f"{SLAB[st]}: effect size per CV", fontsize=10.5); letter(ax, "ab"[j], x=-0.1)
    if j == 0: ax.set_ylabel("|effect size|, CPD vs 6-4PP")
    ax = axs[1, j]; c = Cc[Cc.state == st].set_index("motif"); xx = np.arange(len(DM))
    ax.bar(xx, [c.loc[m, "acc_all"] for m in DM], 0.6, color=[DCOL[m] for m in DM], edgecolor="black", lw=0.8)
    ax.scatter(xx, [c.loc[m, "null_max"] for m in DM], marker="_", s=420, color="black", lw=2.4, zorder=4, label="maximum of permutation null")
    ax.axhline(0.5, color="#666", ls=":", lw=1.5); ax.set_ylim(0.3, 1.05)
    ax.set_xticks(xx, ["BHD2", "partners", "lesion", "BHD3", "DNA"], fontsize=9.5)
    plt.setp(ax.get_xticklabels(), rotation=35, ha="right", rotation_mode="anchor")
    ax.set_title(f"{SLAB[st]}: classification accuracy", fontsize=10.5); letter(ax, "cd"[j], x=-0.1)
    if j == 0: ax.set_ylabel("held-out accuracy")
h_, l_ = axs[1, 0].get_legend_handles_labels(); fig.legend(h_, l_, loc="lower center", bbox_to_anchor=(0.5, -0.04), fontsize=9.5)
save(fig, OUTS, "discrimination")

# ======================= SI: interaction fingerprints of BHD2 / BHD3 =======================
def _reslabels(f):
    L = open(f).read().split("%FLAG RESIDUE_LABEL")[1].split("%FLAG")[0].splitlines()[2:]
    out = []
    for ln in L: out += [ln[i:i + 4].strip() for i in range(0, len(ln), 4) if ln[i:i + 4].strip()]
    return out
RN = {i + 1: n for i, n in enumerate(_reslabels(f"{TI}/rahul_ti_unbiased/pb/system.prmtop")[:620])}
ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "CYX": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "HIE": "H", "HID": "H",
       "HIP": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
def rnm(r): return f"{ONE[RN[r]]}{r}" if r <= 565 else ("lesion 602" if r == 602 else f"{RN[r].replace('D', '')[:1]}{r}")
PROBE = set(range(392, 407)) | set(range(475, 486)) | {584, 585, 602}
S = pd.read_csv(f"{MU}/network/tables/fingerprint_lesion_specificity.csv")
S = S[S.type.isin(["hbond", "saltbr", "wbridge", "hydroph"])]; S = S[~(S.partner.isin(PROBE) & (S.partner < S.probe))]
TL = {"hbond": "H-bond", "saltbr": "salt bridge", "wbridge": "water bridge", "hydroph": "hydrophobic"}
CLS = {"CPD-specific": "CPD only", "6-4PP-specific": "6-4PP only", "conserved": "both"}
for st in ("docked", "pb"):
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 5.4), gridspec_kw=dict(wspace=1.45))
    for ax, mot, lt in zip(axs, ("BHD2", "BHD3"), "ab"):
        a = S[(S.state == st) & (S.probe_motif == mot)].copy(); a["mx"] = a[["occ_64pp", "occ_cpd"]].max(1)
        spec = a[a.cls.isin(["CPD-specific", "6-4PP-specific"])].sort_values("diff", key=abs, ascending=False).head(10)
        cons = a[a.cls == "conserved"].sort_values("mx", ascending=False).head(14 - len(spec))
        b = pd.concat([spec, cons]).sort_values(["cls", "diff"]); y = np.arange(len(b))
        ax.barh(y + 0.2, b.occ_64pp, 0.4, color=C[("64pp", st)], label="6-4PP"); ax.barh(y - 0.2, b.occ_cpd, 0.4, color=C[("cpd", st)], label="CPD")
        ax.set_yticks(y, [f"{rnm(p)}–{rnm(q)}, {TL[t]}" for p, q, t in zip(b.probe, b.partner, b.type)], fontsize=8); ax.invert_yaxis()
        for yy, c_ in zip(y, b.cls): ax.text(1.03, yy, CLS[c_], va="center", fontsize=7.5)
        ax.set_xlim(0, 1.45); ax.set_xticks([0, 0.5, 1.0]); ax.set_xlabel("occupancy")
        n6, nc = (a.cls == "6-4PP-specific").sum(), (a.cls == "CPD-specific").sum()
        ax.set_title(f"{mot}: {n6} 6-4PP only, {nc} CPD only", fontsize=9.5, loc="right", pad=8); letter(ax, lt, x=-1.05, y=1.02)
    h_, l_ = axs[0].get_legend_handles_labels(); fig.legend(h_, l_, loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.97), fontsize=10)
    save(fig, OUTS, f"fingerprint_{st}")

# ======================= SI: per-residue RMSF of the hairpin tips, both lesions =======================
RR = {s_: pd.read_csv(f"{TI}/mech_unbiased/rmsf_res/rmsf_res_{s_}.csv") for s_ in ("64pp", "cpd")}
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.2), gridspec_kw=dict(wspace=0.3))
for ax, (lo, hi, ttl), lt in zip(axs, ((392, 406, "BHD2 tip"), (475, 485, "BHD3 tip")), "ab"):
    for s_ in ("64pp", "cpd"):
        for st in ("docked", "pb"):
            d = RR[s_][(RR[s_].state == st) & (RR[s_].resid.between(lo, hi))]
            g = d.groupby("resid").rmsf; m, e = g.mean(), g.std(ddof=1) / 2
            ax.errorbar(m.index, m.values, yerr=e.values, color=C[(s_, st)], lw=2.0, marker="o", ms=3.5, capsize=2,
                        ls="-" if s_ == "64pp" else (0, (4, 1.5)), label=f"{LAB[s_]} {SLAB[st]}")
    ax.set_xlabel("residue"); ax.set_title(ttl); ax.set_ylim(bottom=0); letter(ax, lt)
    ax.set_xticks(range(lo, hi + 1, 2))
axs[0].set_ylabel("RMSF (Å)"); legend4(fig, y=0.98)
save(fig, OUTS, "rmsf_per_residue")
