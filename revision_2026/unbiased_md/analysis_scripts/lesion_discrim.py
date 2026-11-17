#!/usr/bin/env python
"""Is BHD2 the lesion-sensitive element?  Same test for every motif (11 Nov 2026).
Data: allcv/unbiased (4 x 50 ns per lesion and state, 400 ps frames; heavy-atom CVs, theta branch-centred).
Each CV belongs to exactly ONE motif; CVs mixing motifs (par_b2, b2_b3, par_b3, par_prot) are excluded so no motif borrows another's signal.
 1. replica-level effect size per CV: d = (mean_CPD - mean_6-4PP) / pooled SD of REPLICA MEANS (n = 4 + 4); Welch t over replica means.
 2. leave-one-replica-out classification per motif (logistic regression on that motif's standardised CVs): train on 3+3 replicas,
    test on the held-out 6-4PP + CPD replica pair, all 16 pairings; balanced accuracy (chance 0.5). Feature count differs per motif,
    so the same test is repeated with the motif's single best CV and with a random 3-CV subset (mean over subsets).
 3. permutation null for (2): lesion labels shuffled at the REPLICA level (35 distinct relabelings) -> p value.
Bound (pb) is the primary state; docked reported separately (CPD docked has a global DNA bend from the 64pp-derived start, so DNA-level
features in docked are not interpreted as lesion effects).
Output: tables/discrim_*.csv, figs/discrim/*.png"""
import os, itertools, warnings, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
warnings.filterwarnings("ignore")
MU = "/media/rahul/EXTERNAL_USB/binding_affinity/TI/mech_unbiased"; OUT = f"{MU}/figs/discrim"; os.makedirs(OUT, exist_ok=True)
MOT = {"BHD2": ["b2tip_ax", "d_BHD2", "b2_dna", "b2_pbb", "b2_lbb"],
       "partner bases": ["p584_ax", "p585_ax", "th584", "th585", "gamma", "delta", "p584_p585_stack", "st584_583", "st585_586"],
       "lesion": ["les_ax", "thles", "d_lesion", "les_prot"],
       "BHD3": ["b3tip_ax", "d_BHD3", "b3_dna"],
       "DNA site": ["bend", "twist_site", "rise_site"]}
TH = ["th584", "th585", "thles"]
def cm(a): r = np.radians(a); return np.degrees(np.arctan2(np.sin(r).mean(), np.cos(r).mean()))
def load(st):
    D = {s: [pd.read_csv(f"{MU}/allcv/unbiased/{s}_{st}_r{r}.csv") for r in range(1, 5)] for s in ("64pp", "cpd")}
    for c in TH:   # common branch for both lesions (centred on the pooled circular mean of this state)
        mu = cm(np.concatenate([d[c].values for s in D for d in D[s]]))
        for s in D:
            for d in D[s]: d[c] = (d[c] - mu + 180) % 360 - 180 + mu
    return D
rng = np.random.default_rng(1)
rows_d, rows_c = [], []
for st in ("pb", "docked"):
    D = load(st)
    # 1. effect sizes
    for m, cvs in MOT.items():
        for cv in cvs:
            a = np.array([d[cv].mean() for d in D["64pp"]]); b = np.array([d[cv].mean() for d in D["cpd"]])
            sp = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2); t, p = stats.ttest_ind(b, a, equal_var=False)
            rows_d.append(dict(state=st, motif=m, cv=cv, mean_64pp=a.mean(), mean_cpd=b.mean(), diff=b.mean() - a.mean(),
                               d_replica=(b.mean() - a.mean()) / sp if sp > 0 else np.nan, p_welch=p))
    # 2./3. leave-one-replica-out classification
    def loro(cvs, lab_64=(0, 1, 2, 3), perm=None):
        X = {s: [d[cvs].values for d in D[s]] for s in D}
        reps = [("64pp", i) for i in range(4)] + [("cpd", i) for i in range(4)]
        y_of = {r: (1 if r[0] == "cpd" else 0) for r in reps}
        if perm is not None: y_of = {r: perm[k] for k, r in enumerate(reps)}
        pos = [r for r in reps if y_of[r] == 1]; neg = [r for r in reps if y_of[r] == 0]
        accs = []
        for tp, tn in itertools.product(pos, neg):
            tr = [r for r in reps if r not in (tp, tn)]
            Xtr = np.vstack([X[r[0]][r[1]] for r in tr]); ytr = np.concatenate([[y_of[r]] * len(X[r[0]][r[1]]) for r in tr])
            sc = StandardScaler().fit(Xtr); clf = LogisticRegression(max_iter=2000).fit(sc.transform(Xtr), ytr)
            a1 = (clf.predict(sc.transform(X[tp[0]][tp[1]])) == 1).mean(); a0 = (clf.predict(sc.transform(X[tn[0]][tn[1]])) == 0).mean()
            accs.append((a1 + a0) / 2)
        return float(np.mean(accs)), float(np.min(accs))
    perms = []
    for comb in itertools.combinations(range(8), 4):
        lab = np.zeros(8, int); lab[list(comb)] = 1
        if lab[:4].sum() in (0, 4): continue                          # skip the true labelling and its mirror
        if lab[0] == 1: continue                                     # count each split once (mirror images give the same accuracy)
        perms.append(lab)
    for m, cvs in MOT.items():
        acc, worst = loro(cvs)
        null = np.array([loro(cvs, perm=p_)[0] for p_ in perms]); pval = (1 + (null >= acc).sum()) / (1 + len(null))
        dd = pd.DataFrame([r for r in rows_d if r["state"] == st and r["motif"] == m])
        best = dd.loc[dd.d_replica.abs().idxmax(), "cv"]; acc1, _ = loro([best])
        sub = [loro(list(c))[0] for c in itertools.combinations(cvs, min(3, len(cvs)))] if len(cvs) > 3 else [acc]
        rows_c.append(dict(state=st, motif=m, n_cv=len(cvs), acc_all=acc, acc_worst_pair=worst, p_perm=pval, null_mean=null.mean(), null_max=null.max(),
                           best_cv=best, acc_best_cv=acc1, acc_3cv_mean=float(np.mean(sub)), max_abs_d=float(dd.d_replica.abs().max()),
                           n_cv_p05=int((dd.p_welch < 0.05).sum())))
        print(st, m, f"acc {acc:.2f} (worst {worst:.2f}) p_perm {pval:.3f} null max {null.max():.2f} | best CV {best} {acc1:.2f} | 3-CV {np.mean(sub):.2f}", flush=True)
E = pd.DataFrame(rows_d); Cc = pd.DataFrame(rows_c)
E.round(4).to_csv(f"{MU}/tables/discrim_effect_sizes.csv", index=False); Cc.round(4).to_csv(f"{MU}/tables/discrim_classification.csv", index=False)
# figure: per-motif |d| (replica level) and LORO accuracy, bound and docked
COL = {"BHD2": "#CC79A7", "partner bases": "#E69F00", "lesion": "#0072B2", "BHD3": "#009E73", "DNA site": "#888888"}
fig, axs = plt.subplots(2, 2, figsize=(14, 8.5))
for j, st in enumerate(("pb", "docked")):
    ax = axs[0, j]; e = E[E.state == st]
    for k, m in enumerate(MOT):
        v = e[e.motif == m].d_replica.abs().values; ax.scatter(np.full(len(v), k) + rng.uniform(-0.12, 0.12, len(v)), v, s=55, color=COL[m], zorder=3)
        ax.plot([k - 0.25, k + 0.25], [np.median(v)] * 2, color="black", lw=2)
    ax.set_xticks(range(len(MOT)), list(MOT), fontsize=11); ax.set_ylabel("|effect size| CPD vs 6-4PP\n(replica means, n = 4 + 4)")
    ax.set_title(f"{'bound' if st == 'pb' else 'docked'}: lesion effect on each motif's CVs", loc="left", fontsize=12.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax = axs[1, j]; c = Cc[Cc.state == st].set_index("motif")
    x = np.arange(len(MOT)); w = 0.26
    for k, (col, lab, a) in enumerate((("acc_all", "all motif CVs", 1.0), ("acc_3cv_mean", "3-CV subsets (mean)", 0.6), ("acc_best_cv", "best single CV", 0.35))):
        ax.bar(x + (k - 1) * w, [c.loc[m, col] for m in MOT], w, color=[COL[m] for m in MOT], alpha=a, label=lab, edgecolor="black", lw=0.5)
    ax.scatter(x - w, [c.loc[m, "null_max"] for m in MOT], marker="_", s=300, color="black", zorder=4, label="max of label-permutation null")
    ax.axhline(0.5, color="#888", ls=":"); ax.set_ylim(0.3, 1.02); ax.set_xticks(x, list(MOT), fontsize=11)
    ax.set_ylabel("held-out-replica balanced accuracy"); ax.set_title("can the motif alone tell CPD from 6-4PP?", loc="left", fontsize=12.5)
    ax.spines[["top", "right"]].set_visible(False)
axs[1, 0].legend(fontsize=9, loc="lower left")
fig.tight_layout(); fig.savefig(f"{OUT}/discrimination_by_motif.png", dpi=160); plt.close(fig)
pd.set_option("display.width", 220)
print(Cc.round(3).to_string(index=False))
print(E[E.p_welch < 0.05].sort_values(["state", "motif"]).round(3).to_string(index=False))
