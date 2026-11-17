#!/usr/bin/env python
"""Two-stage TI / MBAR analysis, v3 (SESSION_NOTES sec 47). Same estimators and rules as analyse_2stage_v2.py (sec 46), plus:
 * extended windows: wNN_s3_dvdl.txt (+2.5 ns, continuous from s2) and wNN_sp3_energies.txt are appended when present;
   equilibration = first EQ_PS (1 ns) of each window, so an extended window contributes 6.5 ns instead of 4 ns.
 * --dense: the vdW region REG = [0.788831, 0.879630] is taken from the fresh dense windows w55-w63 (lambdas_bhd2dense.csv);
   stage-V MBAR = MBAR(old windows, lambda' <= REG0) + MBAR(dense windows over REG) + MBAR(old windows, lambda' >= REG1).
   The old windows inside REG give an independent replicate of that segment (reported, checked).
Rules (fixed before the new data were seen): stage-E MBAR on all DV/DL data at 10 ps stride; errors from a paired moving-block
bootstrap (200 x 10 blocks per window, same blocks for TI and MBAR); stage V reported by MBAR, cross-checked by cubic-spline TI
over all windows; checks |TI - MBAR| <= 2 sd_boot(TI - MBAR), adjacent overlap >= 0.03, E-stage linearity < 0.01 kcal/mol.
Usage: analyse_2stage_v3.py <runs dir> <lambdas.csv> <outdir> --label 64pp_pb_bhd2 [--dense]"""
import os, json, argparse, warnings
import numpy as np, pandas as pd, pymbar
from pymbar import timeseries
from scipy.interpolate import CubicSpline
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
kT = 0.0019872041 * 300.0; beta = 1 / kT
ap = argparse.ArgumentParser(); ap.add_argument("runs"); ap.add_argument("lams"); ap.add_argument("out")
ap.add_argument("--label", default="leg"); ap.add_argument("--dense", action="store_true"); ap.add_argument("--nboot", type=int, default=200)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); L = pd.read_csv(a.lams)
EQ_PS = 1000.0; FRACS = (0.25, 0.5, 0.75, 1.0); ESTRIDE = 10; NBOOT, NBLK = a.nboot, 10
REG = (0.788831, 0.879630); DENSE = list(range(55, 64))
P = lambda f: os.path.join(a.runs, f)

def load_dvdl(w):
    parts = [pd.read_csv(P(f"w{w:02d}_dvdl.txt"), sep=r"\s+", header=None, names=["seg", "t", "d"])]
    if os.path.exists(P(f"w{w:02d}_s3_dvdl.txt")):
        parts.append(pd.read_csv(P(f"w{w:02d}_s3_dvdl.txt"), sep=r"\s+", header=None, names=["seg", "t", "d"]))
    X = pd.concat(parts); t = X.t.to_numpy() - X.t.iloc[0]; d = X.d.to_numpy()
    assert len(d) in (5000, 7500, 10000) and np.all(np.abs(np.diff(t) - 1.0) < 1e-6), f"w{w}: {len(d)} points / time not continuous"
    return t, d, len(parts) == 2

def load_sp(w, lams):
    """frames (time order: s1/s2 then s3) x lams; overflow -> +inf"""
    out = []
    for fn in (f"w{w:02d}_sp_energies.txt", f"w{w:02d}_sp3_energies.txt"):
        if not os.path.exists(P(fn)): continue
        rows = []
        for line in open(P(fn)):
            if line.startswith("#") or not line.strip(): continue
            c = line.split(); rows.append((float(c[0]), float(c[1]), np.inf if ("*" in line or c[2] == "inf") else float(c[2])))
        X = np.array(rows)
        for fr in np.unique(X[:, 0]):
            Y = X[X[:, 0] == fr]; row = []
            for l in lams:
                m = np.abs(Y[:, 1] - l) < 1e-6; assert m.sum() <= 1, f"w{w} frame {fr}: lambda {l} duplicated"; row.append(Y[m, 2][0] if m.sum() else np.nan)  # NaN = not computed (old windows at new dense lambdas)
            out.append(row)
    return np.array(out)

def ineff(x): return max(1.0, float(timeseries.statistical_inefficiency(x))) if len(x) > 10 and x.std() > 1e-10 else 1.0
def block_sem(x, nb=5): b = np.array([c.mean() for c in np.array_split(x, nb)]); return float(b.std(ddof=1) / np.sqrt(nb))
def sl(x, lo, hi): return x[int(lo * len(x)):int(hi * len(x))]
def rs(x, idx): b = np.array_split(np.arange(len(x)), NBLK); return x[np.concatenate([b[j] for j in idx])]
sd = lambda x: float(np.std(x, ddof=1))

def mbar_run(u_blocks, nf=0, nb=0, save=None):
    N = np.array([0] * nf + [b.shape[0] for b in u_blocks] + [0] * nb); u = np.concatenate([b.T for b in u_blocks], axis=1)
    u = u - np.min(np.where(np.isfinite(u), u, np.inf), axis=0, keepdims=True)
    m = pymbar.MBAR(u, N, solver_protocol="robust", maximum_iterations=100000); fe = m.compute_free_energy_differences()
    res = dict(dF=float(fe["Delta_f"][0, -1] * kT), err=float(fe["dDelta_f"][0, -1] * kT), profile=(fe["Delta_f"][0] * kT).tolist())
    if save:
        O = np.array(m.compute_overlap()["matrix"])[nf:len(N) - nb, nf:len(N) - nb]; np.save(os.path.join(a.out, f"{a.label}_{save}_overlap.npy"), O)
        fig, ax = plt.subplots(figsize=(6.5, 5.5)); im = ax.imshow(O, vmin=0, vmax=min(1, O.max()), cmap="viridis"); plt.colorbar(im)
        ax.set_title(f"{a.label}: MBAR overlap, {save}"); ax.set_xlabel("window (lambda order)"); ax.set_ylabel("window (lambda order)")
        fig.savefig(os.path.join(a.out, f"{a.label}_{save}_overlap.png"), dpi=150, bbox_inches="tight"); plt.close(fig)
        adj = [float(O[i, i + 1]) for i in range(O.shape[0] - 1)]; res.update(adj=adj, adj_min=min(adj), adj_median=float(np.median(adj)))
    return res

# ---------------- load ----------------
LAM = dict(zip(L.window, L["lambda"])); WGT = dict(zip(L.window, L.weight)); STG = dict(zip(L.window, L.stage))
wE = sorted([w for w in L.window if STG[w] == "E"], key=lambda w: LAM[w])
wOld = sorted([w for w in L.window if STG[w] == "V" and w <= 54], key=lambda w: LAM[w])
wDen = sorted([w for w in DENSE if w in LAM and os.path.exists(P(f"w{w:02d}_dvdl.txt"))], key=lambda w: LAM[w]) if a.dense else []
if a.dense: assert len(wDen) == int(os.environ.get("NDENSE", 9)), f"dense windows present: {wDen}"
D, EXT = {}, []
for w in wE + wOld + wDen:
    t, d, ext = load_dvdl(w); d = d[t >= EQ_PS]; D[w] = d
    if ext: EXT.append(w)
lamV_all = sorted(set(LAM[w] for w in wOld + wDen))
SPE = {w: load_sp(w, lamV_all) for w in wOld + wDen}
for w in wOld + wDen: assert np.all(np.isfinite(SPE[w][:, lamV_all.index(LAM[w])])), f"w{w}: non-finite energy at own lambda"
if a.dense:
    segA = [w for w in wOld if LAM[w] <= REG[0] + 1e-9]; segC = [w for w in wOld if LAM[w] >= REG[1] - 1e-9]
    segOld = [w for w in wOld if REG[0] - 1e-9 <= LAM[w] <= REG[1] + 1e-9]; SEGS = [("A", segA), ("B_dense", wDen), ("C", segC)]
else: SEGS = [("all", wOld)]

def stats(w):
    d = D[w]; g = ineff(d); h1, h2 = d[:len(d) // 2], d[len(d) // 2:]
    z = (h2.mean() - h1.mean()) / np.hypot(h1.std() / np.sqrt(len(h1) / ineff(h1)), h2.std() / np.sqrt(len(h2) / ineff(h2)))
    sem = max(d.std() / np.sqrt(len(d) / g), block_sem(d))
    return dict(w=int(w), lam=float(LAM[w]), weight=float(WGT[w]), mean=float(d.mean()), sem=float(sem), g=g, neff=len(d) / g,
                drift_z=float(z), ns_used=len(d) / 1000.0, extended=w in EXT)

# ---------------- estimators on (possibly resampled / time-sliced) data ----------------
def est_E(Dx):
    lam = np.array([LAM[w] for w in wE]); lam_all = np.concatenate([[0.0], lam, [0.5]])
    ti = float(sum(WGT[w] * Dx[w].mean() for w in wE))
    ub = [beta * (lam_all[None, :] - LAM[w]) * Dx[w][::ESTRIDE][:, None] for w in wE]
    return ti, mbar_run(ub, 1, 1)["dF"]

def seg_mbar(ws, Ex, save=None):
    cols = [lamV_all.index(LAM[w]) for w in ws]; U = [beta * Ex[w][:, cols] for w in ws]
    assert not any(np.isnan(u).any() for u in U), f"segment {ws}: energies missing at a segment lambda"
    return mbar_run(U, save=save)

def est_V(Dx, Ex, save=False):
    mb = 0.0; out = {}
    for name, ws in SEGS:
        r = seg_mbar(ws, Ex, save=(f"V_{name}" if save else None)); mb += r["dF"]; out[name] = r
    # spline TI over all windows used (duplicate lambdas at segment joins -> mean of the two windows)
    pts = {}
    for _, ws in SEGS:
        for w in ws: pts.setdefault(round(LAM[w], 6), []).append(Dx[w].mean())
    x = np.array(sorted(pts)); y = np.array([np.mean(pts[k]) for k in x])
    spl = float(CubicSpline(x, y).integrate(0.0, 1.0)); gl = float(sum(WGT[w] * Dx[w].mean() for w in wOld))
    reg = None
    if a.dense: reg = dict(old=seg_mbar(segOld, Ex)["dF"], dense=out["B_dense"]["dF"])
    return mb, spl, gl, out, reg

# ---------------- point estimates ----------------
rep = {"label": a.label, "eq_ps": EQ_PS, "extended_windows": EXT, "dense": a.dense, "segments": {k: [int(w) for w in v] for k, v in SEGS}}
checks = {}
E_TI, E_MB = est_E(D); lam = np.array([LAM[w] for w in wE])
mbn = mbar_run([beta * (lam[None, :] - LAM[w]) * D[w][::ESTRIDE][:, None] for w in wE], save="E")
devs = []
for w in wE:
    X = np.loadtxt(P(f"w{w:02d}_sp_energies.txt"), comments="#")
    for fr in np.unique(X[:, 0]):
        Y = X[X[:, 0] == fr]; Y = Y[np.argsort(Y[:, 1])]; j = int(np.argmin(np.abs(Y[:, 1] - LAM[w]))); e0, dv = Y[j, 2], Y[j, 3]
        devs += [abs((Y[k, 2] - e0) - (Y[k, 1] - Y[j, 1]) * dv) for k in range(len(Y)) if k != j and 0.0 <= Y[k, 1] <= 0.5]
V_MB, V_SPL, V_GL, VSEG, VREG = est_V(D, SPE, save=True)
V_adj = [x for r in VSEG.values() for x in r["adj"]]

# time convergence (fractions of each window's production part)
def frac(lo, hi): return {w: sl(D[w], lo, hi) for w in D}, {w: sl(SPE[w], lo, hi) for w in SPE}
conv = {k: [] for k in ("E_TI_fwd", "E_MB_fwd", "E_TI_rev", "E_MB_rev", "V_MB_fwd", "V_MB_rev", "V_SPL_fwd", "V_SPL_rev")}
for f in FRACS:
    for tag, (lo, hi) in (("fwd", (0, f)), ("rev", (1 - f, 1))):
        Dx, Ex = frac(lo, hi); ti, mb = est_E(Dx); vm, vs, _, _, _ = est_V(Dx, Ex)
        conv[f"E_TI_{tag}"].append(ti); conv[f"E_MB_{tag}"].append(mb); conv[f"V_MB_{tag}"].append(vm); conv[f"V_SPL_{tag}"].append(vs)

# ---------------- paired block bootstrap ----------------
rng = np.random.default_rng(20261106); BS = {k: [] for k in ("E_TI", "E_MB", "V_MB", "V_SPL", "V_GL", "REG_old", "REG_dense")}
for r in range(NBOOT):
    pick = {w: rng.integers(0, NBLK, NBLK) for w in D}
    Dx = {w: rs(D[w], pick[w]) for w in D}; Ex = {w: rs(SPE[w], pick[w]) for w in SPE}
    ti, mb = est_E(Dx); vm, vs, vg, _, rg = est_V(Dx, Ex)
    for k, v in (("E_TI", ti), ("E_MB", mb), ("V_MB", vm), ("V_SPL", vs), ("V_GL", vg)): BS[k].append(v)
    if rg: BS["REG_old"].append(rg["old"]); BS["REG_dense"].append(rg["dense"])
BS = {k: np.array(v) for k, v in BS.items()}; np.savez(os.path.join(a.out, f"{a.label}_bootstrap.npz"), **BS)

# ---------------- report ----------------
pwE = [stats(w) for w in wE]; pwV = [stats(w) for w in wOld + wDen]
rep["E"] = dict(TI=E_TI, TI_err_boot=sd(BS["E_TI"]), MBAR=dict(dF=E_MB, adj_min=mbn["adj_min"], adj_median=mbn["adj_median"], adj=mbn["adj"]),
                MBAR_err_boot=sd(BS["E_MB"]), diff_err_boot=sd(BS["E_TI"] - BS["E_MB"]), linearity_max_dev_kcal=float(max(devs)),
                TI_forward=conv["E_TI_fwd"], TI_reverse=conv["E_TI_rev"], MBAR_forward=conv["E_MB_fwd"], MBAR_reverse=conv["E_MB_rev"],
                MBAR_half1=conv["E_MB_fwd"][1], MBAR_half2=conv["E_MB_rev"][1], TI_half1=conv["E_TI_fwd"][1], TI_half2=conv["E_TI_rev"][1],
                n_windows_drift_z_gt3=sum(abs(p["drift_z"]) > 3 for p in pwE), per_window=pwE)
rep["V"] = dict(TI=V_GL, TI_GL_err_boot=sd(BS["V_GL"]), TI_spline_all=V_SPL, TI_spline_err_boot=sd(BS["V_SPL"]),
                MBAR=dict(dF=V_MB, adj_min=min(V_adj), adj_median=float(np.median(V_adj)), segments={k: dict(dF=v["dF"], adj=v["adj"]) for k, v in VSEG.items()}),
                MBAR_err_boot=sd(BS["V_MB"]), TIspl_minus_MBAR_err_boot=sd(BS["V_SPL"] - BS["V_MB"]),
                MBAR_forward=conv["V_MB_fwd"], MBAR_reverse=conv["V_MB_rev"], TIspl_forward=conv["V_SPL_fwd"], TIspl_reverse=conv["V_SPL_rev"],
                MBAR_half1=conv["V_MB_fwd"][1], MBAR_half2=conv["V_MB_rev"][1],
                n_windows_drift_z_gt3=sum(abs(p["drift_z"]) > 3 for p in pwV), per_window=pwV)
checks["E_linear"] = max(devs) < 0.01
checks["E_TI_vs_MBAR_boot"] = abs(E_TI - E_MB) <= 2 * rep["E"]["diff_err_boot"]
checks["E_overlap"] = mbn["adj_min"] >= 0.03
checks["V_TIspline_vs_MBAR_boot"] = abs(V_SPL - V_MB) <= 2 * rep["V"]["TIspl_minus_MBAR_err_boot"]
checks["V_overlap"] = min(V_adj) >= 0.03
if VREG:
    e = sd(BS["REG_old"] - BS["REG_dense"]); rep["V"]["region"] = dict(old=VREG["old"], dense=VREG["dense"], diff_err_boot=e,
                                                                        old_err_boot=sd(BS["REG_old"]), dense_err_boot=sd(BS["REG_dense"]))
    checks["V_region_old_vs_dense"] = abs(VREG["old"] - VREG["dense"]) <= 2 * e
tot = dict(recommended=E_MB + V_MB, recommended_err=sd(BS["E_MB"] + BS["V_MB"]),
           TI_E_GL_plus_V_spline=E_TI + V_SPL, TI_E_GL_plus_V_spline_err=sd(BS["E_TI"] + BS["V_SPL"]),
           MBAR_forward=[e + v for e, v in zip(conv["E_MB_fwd"], conv["V_MB_fwd"])], MBAR_reverse=[e + v for e, v in zip(conv["E_MB_rev"], conv["V_MB_rev"])],
           TI_forward=[e + v for e, v in zip(conv["E_TI_fwd"], conv["V_SPL_fwd"])], TI_reverse=[e + v for e, v in zip(conv["E_TI_rev"], conv["V_SPL_rev"])])
tot["MBAR_err_boot"] = tot["recommended_err"]; rep["total"] = tot; rep["checks"] = checks
json.dump(rep, open(os.path.join(a.out, f"{a.label}_2stage.json"), "w"), indent=1, default=float)
print(f"[{a.label}] extended windows: {EXT}  dense: {wDen}")
print(f"[{a.label}] E: TI {E_TI:.3f} +- {rep['E']['TI_err_boot']:.3f} | MBAR {E_MB:.3f} +- {rep['E']['MBAR_err_boot']:.3f} | diff {E_TI-E_MB:+.3f} +- {rep['E']['diff_err_boot']:.3f} | overlap min {mbn['adj_min']:.3f} | halves(MBAR) {rep['E']['MBAR_half1']:.2f}/{rep['E']['MBAR_half2']:.2f}")
print(f"[{a.label}] V: MBAR {V_MB:.3f} +- {rep['V']['MBAR_err_boot']:.3f} | spline TI {V_SPL:.3f} | diff {V_SPL-V_MB:+.3f} +- {rep['V']['TIspl_minus_MBAR_err_boot']:.3f} | GL TI {V_GL:.3f} | overlap min {min(V_adj):.3f} | halves {rep['V']['MBAR_half1']:.2f}/{rep['V']['MBAR_half2']:.2f}")
if VREG: print(f"[{a.label}] region {REG}: old windows {VREG['old']:.3f} +- {sd(BS['REG_old']):.3f} | dense {VREG['dense']:.3f} +- {sd(BS['REG_dense']):.3f}")
print(f"[{a.label}] TOTAL: MBAR {tot['recommended']:.3f} +- {tot['recommended_err']:.3f} | TI {tot['TI_E_GL_plus_V_spline']:.3f} +- {tot['TI_E_GL_plus_V_spline_err']:.3f}")
print(f"[{a.label}] CHECKS:", {k: ("PASS" if v else "FAIL") for k, v in checks.items()})
# figure: dU/dl and convergence
fig, ax = plt.subplots(2, 2, figsize=(12, 8.5))
for k, (st, pw) in enumerate((("E", pwE), ("V", pwV))):
    for sel, mk, lab in ((lambda p: p["w"] < 55 and p["weight"] > 0, "o-", "GL windows"), (lambda p: p["w"] < 55 and p["weight"] == 0, "s", "extra (MBAR-only)"),
                         (lambda p: p["w"] >= 55, "^", "dense resampling")):
        q = [p for p in pw if sel(p)]
        if q: ax[0, k].errorbar([p["lam"] for p in q], [p["mean"] for p in q], yerr=[p["sem"] for p in q], fmt=mk, ms=4, capsize=2, label=lab, mfc="none" if mk != "o-" else None)
    ax[0, k].set_xlabel("lambda" if st == "E" else "lambda' (vdW stage)"); ax[0, k].set_ylabel("<dU/dlambda> (kcal/mol)"); ax[0, k].set_title(f"{a.label}: stage {st}"); ax[0, k].legend(fontsize=8)
    R = rep[st]; F = [100 * f for f in FRACS]
    keys = (("TI_forward", "TI fwd"), ("TI_reverse", "TI rev"), ("MBAR_forward", "MBAR fwd"), ("MBAR_reverse", "MBAR rev")) if st == "E" else \
           (("TIspl_forward", "spline-TI fwd"), ("TIspl_reverse", "spline-TI rev"), ("MBAR_forward", "MBAR fwd"), ("MBAR_reverse", "MBAR rev"))
    for key, lab in keys: ax[1, k].plot(F, R[key], "o-" if "fwd" in lab else "s--", ms=4, label=lab)
    c = R["MBAR"]["dF"]; e = R["MBAR_err_boot"]; ax[1, k].axhspan(c - e, c + e, alpha=0.2, color="grey")
    ax[1, k].set_xlabel("% of production data used (from start: fwd, from end: rev)"); ax[1, k].set_ylabel("dG (kcal/mol)"); ax[1, k].legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(a.out, f"{a.label}_dUdl_convergence.png"), dpi=150); plt.close(fig)
