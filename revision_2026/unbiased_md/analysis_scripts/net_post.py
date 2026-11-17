#!/usr/bin/env python
"""Post-processing of net_engine.py: (1) interaction fingerprints, (2) BHD2/BHD3 conformational clustering, (3) residue interaction
network, (4) time-resolved contacts, (5) dynamic cross-correlation (DCC).  Errors/robustness always over the 4 replicas.
Writes network/tables/*.csv, network/figs/*.png, network/RESULTS.md.  usage: net_post.py"""
import os, glob, itertools, warnings, numpy as np, pandas as pd, networkx as nx, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
warnings.filterwarnings("ignore")
MU = "/media/rahul/EXTERNAL_USB/binding_affinity/TI/mech_unbiased"; N = os.environ.get("NETDIR", f"{MU}/network"); T, F = f"{N}/tables", f"{N}/figs"
os.makedirs(T, exist_ok=True); os.makedirs(F, exist_ok=True)
SYS = [("64pp", "docked"), ("64pp", "pb"), ("cpd", "docked"), ("cpd", "pb")]
LAB = {"64pp": "6-4PP", "cpd": "CPD"}; SL = {"docked": "docked", "pb": "bound"}
C = {("64pp", "docked"): "#0072B2", ("64pp", "pb"): "#E69F00", ("cpd", "docked"): "#009E73", ("cpd", "pb"): "#CC79A7"}
Z = {(l, s, r): dict(np.load(f"{N}/raw/{l}_{s}_r{r}.npz")) for l, s in SYS for r in range(1, 5)}
P = Z[("64pp", "pb", 1)]["probe"]; S = Z[("64pp", "pb", 1)]["site"]
TYPES = ["hbond", "saltbr", "hydroph", "dna_bb", "dna_base", "wbridge", "contact"]
def motif(r): return "BHD2" if 392 <= r <= 406 else "BHD3" if 475 <= r <= 485 else "partner" if r in (584, 585) else "lesion" if r == 602 else "DNA" if r >= 566 else "protein"
REPORT = []
# ---------------- sanity: clustering frame == density-map frame (BHD2/BHD3 mean positions vs maps meanpos) ----------------
for l, s in SYS:
    mp = np.load(f"{MU}/maps/raw/{l}_{s}__unbiased.npz"); K = dict(zip(mp["keys"], mp["meanpos"]))
    for k in ("b2", "b3"):
        names = Z[(l, s, 1)][f"{k}_names"]; X = np.concatenate([Z[(l, s, r)][k] for r in range(1, 5)]).mean(0)
        d = np.sqrt(((X - np.array([K[n] for n in names])) ** 2).sum(1).mean())
        REPORT.append(f"sanity {l} {s} {k}: |mean(engine) - mean(density maps)| = {d:.3f} A (must be ~0)")
# ---------------- (1) fingerprints ----------------
rows = []
for l, s in SYS:
    for t in TYPES:
        occ = np.array([Z[(l, s, r)][f"fp_{t}"].mean(0) for r in range(1, 5)])   # rep x probe x site
        m = occ.mean(0)
        for i, j in zip(*np.where(occ.max(0) >= 0.1)):
            rows.append(dict(system=f"{l}_{s}", lesion=l, state=s, type=t, probe=int(P[i]), probe_motif=motif(P[i]), partner=int(S[j]),
                             partner_motif=motif(S[j]), occ=m[i, j], sem=occ[:, i, j].std(ddof=1) / 2, **{f"r{r + 1}": occ[r, i, j] for r in range(4)}))
FP = pd.DataFrame(rows); FP.round(3).to_csv(f"{T}/fingerprint_occupancy.csv", index=False)
# lesion specificity per state: CPD vs 6-4PP, Welch t over replica occupancies, |diff| >= 0.3
spec = []
for s in ("docked", "pb"):
    for t in TYPES:
        a = FP[(FP.state == s) & (FP.type == t)]
        keys = set(zip(a.probe, a.partner))
        for pr, pa in keys:
            v = {l: a[(a.lesion == l) & (a.probe == pr) & (a.partner == pa)][[f"r{k}" for k in range(1, 5)]].values for l in ("64pp", "cpd")}
            x = v["64pp"][0] if len(v["64pp"]) else np.zeros(4); y = v["cpd"][0] if len(v["cpd"]) else np.zeros(4)
            dlt = y.mean() - x.mean(); p = stats.ttest_ind(y, x, equal_var=False).pvalue if (x.std() + y.std()) > 0 else (0.0 if dlt != 0 else 1.0)
            cls = "conserved" if (x.mean() >= 0.3 and y.mean() >= 0.3 and abs(dlt) < 0.3) else ("CPD-specific" if dlt >= 0.3 and p < 0.05 else ("6-4PP-specific" if dlt <= -0.3 and p < 0.05 else "other"))
            spec.append(dict(state=s, type=t, probe=pr, probe_motif=motif(pr), partner=pa, partner_motif=motif(pa), occ_64pp=x.mean(), occ_cpd=y.mean(), diff=dlt, p=p, cls=cls))
SP = pd.DataFrame(spec); SP.round(3).to_csv(f"{T}/fingerprint_lesion_specificity.csv", index=False)
cnt = SP[SP.type != "contact"].groupby(["state", "probe_motif", "cls"]).size().unstack(fill_value=0)
REPORT.append("fingerprint classes (non-generic interaction types) by state and probe motif:\n" + cnt.to_string())
# figure: BHD2 + BHD3 polar/electrostatic interactions (hbond, saltbr, wbridge) with occupancy >= 0.3 somewhere
for s in ("pb", "docked"):
    a = SP[(SP.state == s) & SP.type.isin(["hbond", "saltbr", "wbridge", "hydroph"]) & SP.probe_motif.isin(["BHD2", "BHD3"]) & ((SP.occ_64pp >= 0.3) | (SP.occ_cpd >= 0.3))]
    a = a.sort_values(["probe", "type", "partner"])
    if not len(a): continue
    lab = [f"{r.probe}–{r.partner} {r.type}" for r in a.itertuples()]
    fig, ax = plt.subplots(figsize=(7.5, 0.26 * len(a) + 1.6)); y = np.arange(len(a))
    ax.barh(y + 0.2, a.occ_64pp, 0.4, color=C[("64pp", s)], label="6-4PP"); ax.barh(y - 0.2, a.occ_cpd, 0.4, color=C[("cpd", s)], label="CPD")
    for yy, c in zip(y, a.cls):
        if c != "other": ax.text(1.02, yy, c, va="center", fontsize=7.5, color="#333")
    ax.set_yticks(y, lab, fontsize=7.5); ax.invert_yaxis(); ax.set_xlim(0, 1.3); ax.set_xlabel("occupancy (mean of 4 replicas)")
    ax.set_title(f"{SL[s]}: BHD2/BHD3 polar, ionic, water-bridge and hydrophobic interactions", loc="left", fontsize=10.5)
    ax.legend(fontsize=8, loc="lower right"); ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(f"{F}/fingerprint_{s}.png", dpi=150); plt.close(fig)
# ---------------- (2) clustering of BHD2 / BHD3 (per state, both lesions pooled, common DNA-core frame) ----------------
cl_rows = []
for k, nm in (("b2", "BHD2"), ("b3", "BHD3")):
    for s in ("pb", "docked"):
        X = np.concatenate([Z[(l, s, r)][k].reshape(len(Z[(l, s, r)][k]), -1) for l in ("64pp", "cpd") for r in range(1, 5)])
        lab_l = np.concatenate([[l] * len(Z[(l, s, r)][k]) for l in ("64pp", "cpd") for r in range(1, 5)]); lab_r = np.concatenate([[r] * len(Z[(l, s, r)][k]) for l in ("64pp", "cpd") for r in range(1, 5)])
        best = None
        for nc in range(2, 7):
            km = KMeans(nc, n_init=20, random_state=0).fit(X); sc = silhouette_score(X, km.labels_, sample_size=min(800, len(X)), random_state=0)
            if best is None or sc > best[0]: best = (sc, nc, km)
        sc, nc, km = best; L = km.labels_
        ncoord = X.shape[1] // 3
        rms_between = np.array([[np.sqrt(((km.cluster_centers_[i] - km.cluster_centers_[j]) ** 2).reshape(-1, 3).sum(1).mean()) for j in range(nc)] for i in range(nc)])
        for c in range(nc):
            for l in ("64pp", "cpd"):
                pops = [np.mean(L[(lab_l == l) & (lab_r == r)] == c) for r in range(1, 5)]
                cl_rows.append(dict(hairpin=nm, state=s, k=nc, silhouette=sc, cluster=c, lesion=l, pop=np.mean(pops), pop_sem=np.std(pops, ddof=1) / 2,
                                    **{f"r{r}": pops[r - 1] for r in range(1, 5)}))
        REPORT.append(f"clustering {nm} {s}: k={nc} (silhouette {sc:.2f}); centroid separation (A RMSD) min {rms_between[rms_between > 0].min():.1f} max {rms_between.max():.1f}")
        # medoid frames -> pdb (heavy atoms of the hairpin only, site frame of the density maps)
        names = Z[("64pp", s, 1)][f"{k}_names"]
        for c in range(nc):
            idx = np.where(L == c)[0]; med = idx[np.argmin(((X[idx] - km.cluster_centers_[c]) ** 2).sum(1))]
            xyz = X[med].reshape(-1, 3)
            with open(f"{N}/figs/{nm}_{s}_cluster{c}_medoid.pdb", "w") as fh:
                for i_, (n_, x_) in enumerate(zip(names, xyz)):
                    r_, a_ = n_.split(":"); fh.write("ATOM  %5d %-4s RES A%4s    %8.3f%8.3f%8.3f  1.00  0.00\n" % (i_ + 1, a_ if len(a_) == 4 else " " + a_, r_, *x_))
CL = pd.DataFrame(cl_rows); CL.round(3).to_csv(f"{T}/clusters.csv", index=False)
fig, axs = plt.subplots(2, 2, figsize=(12, 7))
for ax, (nm, s) in zip(axs.flat, (("BHD2", "pb"), ("BHD2", "docked"), ("BHD3", "pb"), ("BHD3", "docked"))):
    d = CL[(CL.hairpin == nm) & (CL.state == s)]; ks = sorted(d.cluster.unique()); x = np.arange(len(ks))
    for j, l in enumerate(("64pp", "cpd")):
        e = d[d.lesion == l].set_index("cluster").loc[ks]
        ax.bar(x + (j - 0.5) * 0.38, e["pop"], 0.38, yerr=e.pop_sem, capsize=3, color=C[(l, s)], label=LAB[l])
        for jj, kk in enumerate(ks): ax.scatter(np.full(4, x[jj] + (j - 0.5) * 0.38), e.loc[kk, ["r1", "r2", "r3", "r4"]], color="black", s=9, zorder=3)
    ax.set_xticks(x, [f"C{k}" for k in ks]); ax.set_ylim(0, 1); ax.set_ylabel("population"); ax.spines[["top", "right"]].set_visible(False)
    ax.set_title(f"{nm} {SL[s]}: k={d.k.iloc[0]}, silhouette {d.silhouette.iloc[0]:.2f}", loc="left", fontsize=11); ax.legend(fontsize=9)
fig.suptitle("Hairpin conformational clusters (k-means on DNA-core-superposed heavy atoms; dots = replicas)", fontsize=12, x=0.01, ha="left")
fig.tight_layout(); fig.savefig(f"{F}/clusters.png", dpi=150); plt.close(fig)
# ---------------- (5) DCC (needed for network weights) ----------------
def kab(P_, Q_):
    pc, qc = P_.mean(0), Q_.mean(0); U, s_, Vt = np.linalg.svd((P_ - pc).T @ (Q_ - qc)); d = np.sign(np.linalg.det(Vt.T @ U.T)); R = Vt.T @ np.diag([1, 1, d]) @ U.T; return R, qc - R @ pc
def dcc(X):
    m = X.mean(0)
    for _ in range(3):
        X = np.array([x @ kab(x, m)[0].T + kab(x, m)[1] for x in X]); m = X.mean(0)
    D = X - m; c = np.einsum("fik,fjk->ij", D, D) / len(X); n_ = np.sqrt(np.diag(c)); return c / np.outer(n_, n_)
DCC = {(l, s, r): dcc(Z[(l, s, r)]["repx"].astype(float)) for l, s in SYS for r in range(1, 5)}
grp = {g: np.where([motif(r) == g for r in S])[0] for g in ("BHD2", "BHD3", "partner", "lesion")}
dr = []
for l, s in SYS:
    for a, b in (("BHD2", "partner"), ("BHD3", "partner"), ("BHD2", "BHD3"), ("BHD2", "lesion"), ("BHD3", "lesion"), ("partner", "lesion")):
        v = [DCC[(l, s, r)][np.ix_(grp[a], grp[b])].mean() for r in range(1, 5)]
        dr.append(dict(system=f"{l}_{s}", pair=f"{a}–{b}", dcc=np.mean(v), sem=np.std(v, ddof=1) / 2, **{f"r{r}": v[r - 1] for r in range(1, 5)}))
DR = pd.DataFrame(dr); DR.round(3).to_csv(f"{T}/dcc_motif_blocks.csv", index=False)
fig, ax = plt.subplots(figsize=(11, 4.2)); pairs = list(DR.pair.unique()); x = np.arange(len(pairs))
for k, (l, s) in enumerate(SYS):
    d = DR[DR.system == f"{l}_{s}"].set_index("pair").loc[pairs]; xs = x + (k - 1.5) * 0.2
    ax.bar(xs, d["dcc"], 0.2, yerr=d["sem"], capsize=2, color=C[(l, s)], label=f"{LAB[l]} {SL[s]}")
ax.axhline(0, color="black", lw=0.8); ax.set_xticks(x, pairs, fontsize=10); ax.set_ylabel("mean DCC (± SEM, 4 replicas)"); ax.legend(ncol=4, fontsize=9)
ax.spines[["top", "right"]].set_visible(False); fig.tight_layout(); fig.savefig(f"{F}/dcc_blocks.png", dpi=150); plt.close(fig)
# ---------------- (3) residue interaction network (contact occupancy >= 0.5, weight -log|DCC|; per replica) ----------------
nr = []
for l, s in SYS:
    for r in range(1, 5):
        z = Z[(l, s, r)]; occ = z["net"] / z["nfr"]; occ = np.maximum(occ, occ.T); D_ = DCC[(l, s, r)]
        G = nx.Graph()
        for i, j in zip(*np.where(np.triu(occ, 1) >= 0.5)):
            G.add_edge(int(S[i]), int(S[j]), w=float(-np.log(max(abs(D_[i, j]), 1e-3))))
        bc = nx.betweenness_centrality(G, weight="w", normalized=True)
        com = nx.algorithms.community.greedy_modularity_communities(G, weight=None); cid = {n_: k for k, cset in enumerate(com) for n_ in cset}
        def spath(a, b):
            best = (np.inf, None)
            for x_ in a:
                for y_ in b:
                    if x_ in G and y_ in G and nx.has_path(G, x_, y_):
                        L_ = nx.shortest_path_length(G, x_, y_, weight="w")
                        if L_ < best[0]: best = (L_, nx.shortest_path(G, x_, y_, weight="w"))
            return best
        row = dict(system=f"{l}_{s}", rep=r, nodes=G.number_of_nodes(), edges=G.number_of_edges(), n_comm=len(com),
                   bc_BHD2=np.mean([bc.get(x, 0) for x in range(392, 407)]), bc_BHD3=np.mean([bc.get(x, 0) for x in range(475, 486)]),
                   bc_partner=np.mean([bc.get(x, 0) for x in (584, 585)]), bc_lesion=bc.get(602, 0),
                   deg_BHD2=np.mean([G.degree(x) if x in G else 0 for x in range(392, 407)]), deg_BHD3=np.mean([G.degree(x) if x in G else 0 for x in range(475, 486)]),
                   same_comm_BHD2_lesion=float(any(cid.get(x) == cid.get(602) for x in range(392, 407))),
                   same_comm_BHD3_lesion=float(any(cid.get(x) == cid.get(602) for x in range(475, 486))),
                   same_comm_BHD2_BHD3=float(any(cid.get(x) == cid.get(y) for x in range(392, 407) for y in range(475, 486))))
        for a, b, nm in ((range(392, 407), [602], "BHD2_lesion"), (range(475, 486), [602], "BHD3_lesion"), (range(392, 407), range(475, 486), "BHD2_BHD3"),
                         (range(392, 407), [584, 585], "BHD2_partner"), (range(475, 486), [584, 585], "BHD3_partner")):
            L_, path = spath(list(a), list(b)); row[f"path_{nm}"] = L_; row[f"hops_{nm}"] = (len(path) - 1) if path else np.nan
            row[f"route_{nm}"] = "-".join(map(str, path)) if path else ""
        nr.append(row)
NR = pd.DataFrame(nr); NR.round(4).to_csv(f"{T}/network_per_replica.csv", index=False)
num = [c for c in NR.columns if c not in ("system", "rep") and not c.startswith("route")]
NS = NR.groupby("system")[num].agg(["mean", "sem"]); NS.round(3).to_csv(f"{T}/network_summary.csv")
fig, axs = plt.subplots(1, 3, figsize=(15, 4.2))
for ax, cols, ttl in ((axs[0], ["bc_BHD2", "bc_BHD3", "bc_partner", "bc_lesion"], "betweenness centrality (motif mean)"),
                      (axs[1], ["hops_BHD2_lesion", "hops_BHD3_lesion", "hops_BHD2_partner", "hops_BHD3_partner", "hops_BHD2_BHD3"], "shortest-path hops"),
                      (axs[2], ["same_comm_BHD2_lesion", "same_comm_BHD3_lesion", "same_comm_BHD2_BHD3"], "fraction of replicas in the same community")):
    x = np.arange(len(cols))
    for k, (l, s) in enumerate(SYS):
        d = NR[NR.system == f"{l}_{s}"]; ax.bar(x + (k - 1.5) * 0.2, d[cols].mean(), 0.2, yerr=d[cols].sem(), capsize=2, color=C[(l, s)], label=f"{LAB[l]} {SL[s]}")
    ax.set_xticks(x, [c.split("_", 1)[1].replace("_", "–") for c in cols], fontsize=8.5, rotation=20); ax.set_title(ttl, loc="left", fontsize=10.5); ax.spines[["top", "right"]].set_visible(False)
axs[0].legend(fontsize=8); fig.tight_layout(); fig.savefig(f"{F}/network.png", dpi=150); plt.close(fig)
# ---------------- (4) time-resolved contacts: persistence of BHD2 vs BHD3 DNA contacts ----------------
tr = []; lags = np.arange(0, 60)
fig, axs = plt.subplots(1, 2, figsize=(13, 4.2))
for ax, s in zip(axs, ("pb", "docked")):
    for l in ("64pp", "cpd"):
        for nm, lo, hi, ls in (("BHD2", 392, 406, "-"), ("BHD3", 475, 485, "--")):
            pi = [i for i, r in enumerate(P) if lo <= r <= hi]; dj = [j for j, r in enumerate(S) if r >= 566]
            acs, life, nsim = [], [], []
            for r in range(1, 5):
                fpc = Z[(l, s, r)]["fp_contact"][:, pi][:, :, dj].reshape(len(Z[(l, s, r)]["fp_contact"]), -1).astype(float)
                keep = fpc.mean(0) >= 0.1; x = fpc[:, keep]
                if x.shape[1] == 0: continue
                c0 = (x * x).mean(0); ac = [((x[:len(x) - k] * x[k:]).mean(0).mean() / c0.mean()) if k < len(x) else np.nan for k in lags]; acs.append(ac)
                runs = []
                for col in x.T:
                    d_ = np.diff(np.concatenate([[0], col, [0]])); st_, en_ = np.where(d_ == 1)[0], np.where(d_ == -1)[0]; runs += list((en_ - st_) * 0.4)
                life.append(np.mean(runs)); nsim.append(x.shape[1])
            if acs:
                ax.plot(lags * 0.4, np.nanmean(acs, 0), ls=ls, color=C[(l, s)], lw=2, label=f"{LAB[l]} {nm}")
                tr.append(dict(system=f"{l}_{s}", hairpin=nm, mean_lifetime_ns=np.mean(life), lifetime_sem=np.std(life, ddof=1) / np.sqrt(len(life)),
                               n_contacts=np.mean(nsim), persistence_10ns=np.nanmean([a[25] for a in acs])))
    ax.set_xlabel("lag (ns)"); ax.set_ylabel("contact persistence  ⟨c(t)c(t+τ)⟩/⟨c²⟩"); ax.set_title(f"{SL[s]}: hairpin–DNA residue contacts (solid BHD2, dashed BHD3)", loc="left", fontsize=10.5)
    ax.set_ylim(0, 1.02); ax.legend(fontsize=8.5); ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout(); fig.savefig(f"{F}/contact_persistence.png", dpi=150); plt.close(fig)
TR = pd.DataFrame(tr); TR.round(3).to_csv(f"{T}/contact_lifetimes.csv", index=False)
# ---------------- report ----------------
pd.set_option("display.width", 200)
REPORT += ["\nclusters:\n" + CL.round(2).to_string(index=False), "\nDCC motif blocks:\n" + DR[["system", "pair", "dcc", "sem"]].round(3).to_string(index=False),
           "\nnetwork (mean over replicas):\n" + NR.groupby("system")[["bc_BHD2", "bc_BHD3", "bc_partner", "bc_lesion", "hops_BHD2_lesion", "hops_BHD3_lesion",
                                                                       "hops_BHD2_BHD3", "same_comm_BHD2_lesion", "same_comm_BHD3_lesion", "same_comm_BHD2_BHD3"]].mean().round(3).to_string(),
           "\ncontact lifetimes:\n" + TR.round(2).to_string(index=False),
           "\nlesion-specific BHD2/BHD3 interactions (non-contact types):\n" + SP[SP.cls.isin(["CPD-specific", "6-4PP-specific"]) & SP.probe_motif.isin(["BHD2", "BHD3"]) & (SP.type != "contact")].sort_values(["state", "probe"]).round(2).to_string(index=False)]
open(f"{N}/RESULTS.md", "w").write("# Network / fingerprint / clustering / DCC results (auto-generated by net_post.py)\n\n```\n" + "\n".join(REPORT) + "\n```\n")
print("\n".join(REPORT[:12]))
