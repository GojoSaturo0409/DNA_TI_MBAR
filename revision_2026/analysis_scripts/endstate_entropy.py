#!/usr/bin/env python
"""End-state energies and solute-environment energy/entropy split of each decoupling free energy (reviewer R1.5/R1.6).

Identity (Ben-Naim; Yu & Karplus 1988; Gallicchio & Levy): for a coupling perturbation H1 = H0 + dV,
    G1 - G0 = <dV>_1 - T dS_uv,
where <dV>_1 is the perturbation energy averaged in the COUPLED ensemble and the environment-reorganisation energy cancels
exactly against part of the entropy. Hence only the coupled ensemble of each stage is needed (no decoupled end-state
sampling). For decoupling: dG_dec = dU_dec - T dS_dec with dU_dec = <U_decoupled - U_coupled>_coupled.

Two-stage path, applied per stage at the coupled end of that stage:
  stage E (decharge, lambda 0 -> 0.5, exactly linear): U(l) = U(0) + l*D, D = DV/DL per frame.
     dU_E = 0.5 <D>_{l=0}; frames of w01 (l=0.001111) reweighted exactly to l=0 with weights exp(+beta*l1*D).
  stage V (vdW, lambda' 0 -> 1, softcore): single-point energies of w26 (l'=0.002222) at l'=0.002222 and 0.997778,
     extrapolated to the end points with the exact DV/DL at those lambdas (first order, dl = 0.002222),
     reweighted from l'=0.002222 to 0 with exp(+beta*dl*DVDL(0.002222)).
  -T dS = dG - dU  (per stage and in total; path-dependent split, exact sum identity per stage).
Errors: moving-block bootstrap over time (10 blocks, 1000 replicates), combined in quadrature with the bootstrap error of dG.
Usage: endstate_entropy.py <runs root (runs64v3)> <analysis dir with 64pp_<leg>_2stage.json> <rocklin tsv> <out prefix>"""
import sys, os, json, numpy as np, pandas as pd
root, adir, rtsv, out = sys.argv[1:5]
kT = 0.0019872041 * 300.0; beta = 1 / kT; EQ = 1000.0; NB = 10; NREP = 1000
rng = np.random.default_rng(20261107)
LEGS = [f"{s}_{m}" for s in ("pb", "docked") for m in ("lesion", "partner", "bhd2", "bhd3")]
RC = pd.read_csv(rtsv, sep="\t"); SYSN = os.environ.get("SYSN", "64pp"); RC = RC[RC.system == SYSN].set_index("leg")["total_kcal"]
L1E, L1V = 0.001111, 0.002222; LVHI = 0.997778

def blocks_boot(fn, n):
    idx = np.array_split(np.arange(n), NB); vals = []
    for _ in range(NREP):
        pick = np.concatenate([idx[j] for j in rng.integers(0, NB, NB)]); vals.append(fn(pick))
    return float(np.std(vals, ddof=1))

def stage_E(d):
    X = []
    for f in ("w01_dvdl.txt", "w01_s3_dvdl.txt"):
        p = os.path.join(d, f)
        if os.path.exists(p): X.append(pd.read_csv(p, sep=r"\s+", header=None, names=["seg", "t", "d"]))
    X = pd.concat(X); t = X.t.to_numpy() - X.t.iloc[0]; D = X.d.to_numpy()[t >= EQ]
    def est(ix):
        x = D[ix]; a = beta * L1E * x; w = np.exp(a - a.max()); return 0.5 * float(np.sum(w * x) / np.sum(w))
    ix = np.arange(len(D)); return est(ix), blocks_boot(est, len(D)), float(0.5 * D.mean()), len(D)

def stage_V(d):
    rows = []
    for f in ("w26_sp_energies.txt", "w26_sp3_energies.txt"):
        p = os.path.join(d, f)
        if not os.path.exists(p): continue
        off = 0 if not rows else max(r[0] for r in rows) + 1
        for line in open(p):
            if line.startswith("#") or not line.strip(): continue
            c = line.split()
            if "*" in line or c[2] == "inf": continue
            rows.append((float(c[0]) + off, float(c[1]), float(c[2]), float(c[3])))
    X = pd.DataFrame(rows, columns=["fr", "lam", "E", "dv"])
    lo = X[np.isclose(X.lam, L1V)].set_index("fr"); hi = X[np.isclose(X.lam, LVHI)].set_index("fr")
    fr = lo.index.intersection(hi.index)
    U0 = lo.loc[fr, "E"].to_numpy() - L1V * lo.loc[fr, "dv"].to_numpy()          # extrapolate to l'=0
    U1 = hi.loc[fr, "E"].to_numpy() + (1 - LVHI) * hi.loc[fr, "dv"].to_numpy()   # extrapolate to l'=1
    raw = hi.loc[fr, "E"].to_numpy() - lo.loc[fr, "E"].to_numpy()                # without end-point extrapolation
    dU = U1 - U0; a = beta * L1V * lo.loc[fr, "dv"].to_numpy()
    def est(ix):
        w = np.exp(a[ix] - a[ix].max()); return float(np.sum(w * dU[ix]) / np.sum(w))
    ix = np.arange(len(dU)); return est(ix), blocks_boot(est, len(dU)), float(raw.mean()), len(dU)

R = []
for leg in LEGS:
    d = os.path.join(root, leg); J = json.load(open(os.path.join(adir, f"{SYSN}_{leg}_2stage.json")))
    uE, eE, uE_raw, nE = stage_E(d); uV, eV, uV_raw, nV = stage_V(d)
    gE, gV = J["E"]["MBAR"]["dF"], J["V"]["MBAR"]["dF"]; geE, geV = J["E"]["MBAR_err_boot"], J["V"]["MBAR_err_boot"]
    R.append(dict(leg=leg, state=leg.split("_")[0], motif=leg.split("_")[1],
                  dU_E=uE, dU_E_err=eE, dU_E_unweighted=uE_raw, nE=nE, dU_V=uV, dU_V_err=eV, dU_V_noextrap=uV_raw, nV=nV,
                  dG_E=gE, dG_E_err=geE, dG_V=gV, dG_V_err=geV,
                  mTdS_E=gE - uE, mTdS_E_err=float(np.hypot(eE, geE)), mTdS_V=gV - uV, mTdS_V_err=float(np.hypot(eV, geV)),
                  dU=uE + uV, dU_err=float(np.hypot(eE, eV)), dG=gE + gV, dG_err=J["total"]["MBAR_err_boot"],
                  rocklin=float(RC[leg])))
    r = R[-1]; r["mTdS"] = r["dG"] - r["dU"]; r["mTdS_err"] = float(np.hypot(r["dU_err"], r["dG_err"]))
D = pd.DataFrame(R); D.to_csv(f"{out}_legs.csv", index=False)
dd = []
for m in ("lesion", "partner", "bhd2", "bhd3"):
    p = D[(D.motif == m) & (D.state == "pb")].iloc[0]; q = D[(D.motif == m) & (D.state == "docked")].iloc[0]; h = np.hypot
    row = dict(motif=m)
    for k in ("dU_E", "dU_V", "dU", "dG_E", "dG_V", "mTdS_E", "mTdS_V", "mTdS"):
        row["dd" + k] = p[k] - q[k]; row["dd" + k + "_err"] = float(h(p[k + "_err"], q[k + "_err"]))
    row["ddG_corr"] = (p.dG - p.rocklin) - (q.dG - q.rocklin)
    dd.append(row)
DD = pd.DataFrame(dd); DD.to_csv(f"{out}_ddG.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
print("per leg (kcal/mol): dU = <U_dec - U_coup> in the coupled ensemble of each stage; -TdS = dG - dU")
print(D[["leg", "dU_E", "dU_E_err", "dU_E_unweighted", "dU_V", "dU_V_err", "dU_V_noextrap", "nV", "dG_E", "dG_V", "mTdS_E", "mTdS_V", "dU", "dG", "mTdS", "mTdS_err"]].round(2).to_string(index=False))
print("\npb - docked:"); print(DD.round(2).to_string(index=False))
