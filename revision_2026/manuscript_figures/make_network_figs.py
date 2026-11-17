#!/usr/bin/env python3
"""Hairpin clusters, residue interaction network, contact lifetimes and DCC (bold paper style).
Inputs: mech_unbiased/network/tables/{clusters,network_per_replica,contact_lifetimes,dcc_motif_blocks}.csv"""
import os, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
T = "/media/rahul/EXTERNAL_USB/binding_affinity/TI/mech_unbiased/network/tables"
H = os.path.dirname(os.path.abspath(__file__)); OM = f"{H}/../figures"; OS = f"{H}/../figures_si"
C = {("64pp", "docked"): "#0072B2", ("64pp", "pb"): "#E69F00", ("cpd", "docked"): "#009E73", ("cpd", "pb"): "#CC79A7"}
LAB = {"64pp": "6-4PP", "cpd": "CPD"}; SL = {"docked": "docked", "pb": "bound"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5, "font.weight": "bold", "axes.labelweight": "bold", "axes.titleweight": "bold",
    "axes.titlesize": 11, "axes.labelsize": 11, "xtick.labelsize": 9.5, "ytick.labelsize": 9.5, "axes.linewidth": 1.6, "xtick.major.width": 1.6,
    "ytick.major.width": 1.6, "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "savefig.bbox": "tight", "savefig.dpi": 300})
def letter(ax, s): ax.text(-0.14, 1.06, s, transform=ax.transAxes, fontsize=14, fontweight="bold", va="bottom", ha="right")
def save(fig, d, n): fig.savefig(f"{d}/{n}.pdf"); fig.savefig(f"{d}/{n}.png", dpi=200); plt.close(fig); print("wrote", n)
CL = pd.read_csv(f"{T}/clusters.csv"); NW = pd.read_csv(f"{T}/network_per_replica.csv"); LT = pd.read_csv(f"{T}/contact_lifetimes.csv"); DC = pd.read_csv(f"{T}/dcc_motif_blocks.csv")
def clusters(ax, hp, st, ttl):
    d = CL[(CL.hairpin == hp) & (CL.state == st)]; ks = sorted(d.cluster.unique()); w = 0.38
    for j, s in enumerate(("64pp", "cpd")):
        for k in ks:
            r = d[(d.cluster == k) & (d.lesion == s)].iloc[0]; x = k + (j - 0.5) * (w + 0.03)
            ax.bar(x, r["pop"], w, yerr=r.pop_sem, capsize=3, color=C[(s, st)], error_kw=dict(lw=1.5), label=LAB[s] if k == ks[0] else None)
            ax.scatter(np.full(4, x), r[["r1", "r2", "r3", "r4"]].values.astype(float), s=12, color="black", zorder=3)
    ax.set_xticks(ks, [f"C{k}" for k in ks]); ax.set_ylim(0, 1.05); ax.set_ylabel("population"); ax.set_title(f"{ttl} (silhouette {d.silhouette.iloc[0]:.2f})", fontsize=10)
    ax.legend(fontsize=9, loc="upper left", bbox_to_anchor=(0.0, 1.0))
def bars4(ax, vals, errs, groups, ylabel):
    w = 0.19
    for i, (s, st) in enumerate([("64pp", "docked"), ("64pp", "pb"), ("cpd", "docked"), ("cpd", "pb")]):
        x = np.arange(len(groups)) + (i - 1.5) * (w + 0.01)
        ax.bar(x, [vals[(s, st, g)] for g in groups], w, yerr=[errs[(s, st, g)] for g in groups], capsize=2.5, color=C[(s, st)], error_kw=dict(lw=1.3),
               label=f"{LAB[s]} {SL[st]}")
    ax.set_xticks(range(len(groups)), groups); ax.set_ylabel(ylabel)
# betweenness / lifetimes / dcc dicts
bv, be = {}, {}
for (sysn), g in NW.groupby("system"):
    s, st = sysn.split("_")
    for hp in ("BHD2", "BHD3"):
        v = g[f"bc_{hp}"].values; bv[(s, st, hp)] = v.mean(); be[(s, st, hp)] = v.std(ddof=1) / np.sqrt(len(v))
lv, le = {}, {}
for _, r in LT.iterrows():
    s, st = r.system.split("_"); lv[(s, st, r.hairpin)] = r.mean_lifetime_ns; le[(s, st, r.hairpin)] = r.lifetime_sem
dv, de = {}, {}
for _, r in DC.iterrows():
    s, st = r.system.split("_"); dv[(s, st, r.pair)] = r.dcc; de[(s, st, r.pair)] = r["sem"]

# ---- MAIN figure: distinct roles of BHD2 and BHD3
fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.6), gridspec_kw=dict(wspace=0.35, hspace=0.55)); axs = axs.ravel()
clusters(axs[0], "BHD2", "pb", "BHD2, bound"); axs[0].get_legend().remove(); letter(axs[0], "a")
bars4(axs[1], bv, be, ["BHD2", "BHD3"], "betweenness centrality"); axs[1].set_title("interaction network"); letter(axs[1], "b")
bars4(axs[2], lv, le, ["BHD2", "BHD3"], "contact lifetime (ns)"); axs[2].set_title("hairpin–DNA contacts"); letter(axs[2], "c")
pairs = ["BHD2–partner", "BHD2–lesion", "BHD3–partner", "BHD3–lesion"]
bars4(axs[3], dv, de, pairs, "DCC"); axs[3].set_xticklabels(["BHD2–\npartners", "BHD2–\nlesion", "BHD3–\npartners", "BHD3–\nlesion"], fontsize=8.5)
axs[3].axhline(0, color="black", lw=1.2); axs[3].set_title("correlated motion"); letter(axs[3], "d")
h, l = axs[1].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=4, bbox_to_anchor=(0.5, 0.99), fontsize=9.5)
save(fig, OM, "fig_network")

# ---- SI figure: all cluster panels
fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.6), gridspec_kw=dict(wspace=0.35, hspace=0.5)); axs = axs.ravel()
for ax, (hp, st), lt in zip(axs, [("BHD2", "pb"), ("BHD2", "docked"), ("BHD3", "pb"), ("BHD3", "docked")], "abcd"):
    clusters(ax, hp, st, f"{hp}, {SL[st]}"); letter(ax, lt)
save(fig, OS, "clusters_all")
# ---- SI figure: network hops and communities
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw=dict(wspace=0.35))
hv, he, cv, ce = {}, {}, {}, {}
for sysn, g in NW.groupby("system"):
    s, st = sysn.split("_")
    for nm, col in (("BHD2–lesion", "hops_BHD2_lesion"), ("BHD3–lesion", "hops_BHD3_lesion"), ("BHD2–BHD3", "hops_BHD2_BHD3")):
        v = g[col].values.astype(float); hv[(s, st, nm)] = v.mean(); he[(s, st, nm)] = v.std(ddof=1) / 2
    for nm, col in (("BHD2–lesion", "same_comm_BHD2_lesion"), ("BHD3–lesion", "same_comm_BHD3_lesion"), ("BHD2–BHD3", "same_comm_BHD2_BHD3")):
        v = g[col].values.astype(float); cv[(s, st, nm)] = v.mean(); ce[(s, st, nm)] = v.std(ddof=1) / 2
bars4(axs[0], hv, he, ["BHD2–lesion", "BHD3–lesion", "BHD2–BHD3"], "shortest-path hops"); letter(axs[0], "a")
bars4(axs[1], cv, ce, ["BHD2–lesion", "BHD3–lesion", "BHD2–BHD3"], "same community (fraction)"); letter(axs[1], "b")
for ax in axs: ax.tick_params(axis="x", labelsize=8.5)
h, l = axs[0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=4, bbox_to_anchor=(0.5, 0.98), fontsize=9.5)
save(fig, OS, "network_paths")
