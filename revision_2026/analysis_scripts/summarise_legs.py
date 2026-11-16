#!/usr/bin/env python
"""Summary of the two-stage production analysis (analyse_2stage_v2.py JSONs) for one system:
per-leg dG_E, dG_V, totals (MBAR and TI), checks, convergence; Rocklin correction; ddG(m) = dG_dec,corr(pb) - dG_dec,corr(docked).
Sign conventions: METHODS_two_stage_TI_MBAR.md sec 6 (dG_dec = G(decoupled) - G(coupled); corrected = dG_dec - ddG_Rocklin).
Usage: summarise_legs.py <analysis dir> <system> <rocklin tsv> <out prefix>"""
import sys, json, math, numpy as np, pandas as pd
d, sysn, rtsv, out = sys.argv[1:5]
RC = pd.read_csv(rtsv, sep="\t"); RC = RC[RC.system == sysn].set_index("leg")["total_kcal"]
MOTIFS = ["lesion", "partner", "bhd2", "bhd3"]

def rnd(v, e):
    """editor's rule: error to 1 significant figure (2 if it starts with 1), value to the same decimal place"""
    if not np.isfinite(e) or e <= 0: return f"{v:.2f}"
    p = -int(math.floor(math.log10(e))); lead = int(e * 10 ** p + 1e-9)
    if lead == 1: p += 1
    p = max(p, 0); return f"{round(v, p):.{p}f} ± {round(e, p):.{p}f}"

rows = []
for st in ("pb", "docked"):
    for m in MOTIFS:
        leg = f"{st}_{m}"; J = json.load(open(f"{d}/{sysn}_{leg}_2stage.json")); E, V, T = J["E"], J["V"], J["total"]
        # convergence of the recommended estimate (MBAR E + V): halves and forward/reverse
        h1, h2 = E["MBAR_half1"] + V["MBAR_half1"], E["MBAR_half2"] + V["MBAR_half2"]
        sig = T["MBAR_err_boot"]; z_half = (h2 - h1) / (2 * sig)          # each half has ~sqrt(2)*sigma -> diff sd ~ 2 sigma
        rows.append(dict(leg=leg, motif=m, state=st,
            E_TI=E["TI"], E_MBAR=E["MBAR"]["dF"], E_err=E["MBAR_err_boot"], E_TIminusMBAR=E["TI"] - E["MBAR"]["dF"], E_diff_err=E["diff_err_boot"],
            E_overlap_min=E["MBAR"]["adj_min"], E_linear_maxdev=E["linearity_max_dev_kcal"],
            V_TI_GL=V["TI"], V_TI_spline=V["TI_spline_all"], V_MBAR=V["MBAR"]["dF"], V_err=V["MBAR_err_boot"],
            V_TIspl_minus_MBAR=V["TI_spline_all"] - V["MBAR"]["dF"], V_diff_err=V["TIspl_minus_MBAR_err_boot"], V_overlap_min=V["MBAR"]["adj_min"],
            dG_dec=T["recommended"], dG_dec_err=T["recommended_err"], dG_dec_TI=T["TI_E_GL_plus_V_spline"], dG_dec_TI_err=T["TI_E_GL_plus_V_spline_err"],
            rocklin=float(RC[leg]), dG_corr=T["recommended"] - float(RC[leg]),
            half1=h1, half2=h2, z_half=z_half, fwd=T["MBAR_forward"], rev=T["MBAR_reverse"],
            drift_windows=E["n_windows_drift_z_gt3"] + V["n_windows_drift_z_gt3"],
            checks=J["checks"]))
D = pd.DataFrame(rows); D.to_json(f"{out}_legs.json", orient="records", indent=1)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
print(f"===== {sysn}: per leg (kcal/mol; errors = paired block bootstrap, 200 x 10 blocks) =====")
print(D[["leg", "E_TI", "E_MBAR", "E_err", "E_TIminusMBAR", "E_diff_err", "E_overlap_min", "V_TI_GL", "V_TI_spline", "V_MBAR", "V_err",
         "V_TIspl_minus_MBAR", "V_diff_err", "V_overlap_min"]].round(3).to_string(index=False))
print("\n----- totals, Rocklin, convergence (MBAR E+V) -----")
print(D[["leg", "dG_dec", "dG_dec_err", "dG_dec_TI", "dG_dec_TI_err", "rocklin", "dG_corr", "half1", "half2", "z_half", "drift_windows"]].round(3).to_string(index=False))
print("\nforward (25/50/75/100% of 4 ns) | reverse:")
for r in rows: print(f"  {r['leg']:15s} " + " ".join(f"{x:8.2f}" for x in r["fwd"]) + "  | " + " ".join(f"{x:8.2f}" for x in r["rev"]))
print("\nchecks:")
for r in rows: print(f"  {r['leg']:15s} " + " ".join(f"{k}={'P' if v else 'F'}" for k, v in r["checks"].items()))
# state-to-state
dd = []
for m in MOTIFS:
    p = D[(D.motif == m) & (D.state == "pb")].iloc[0]; q = D[(D.motif == m) & (D.state == "docked")].iloc[0]
    v = p.dG_corr - q.dG_corr; e = math.hypot(p.dG_dec_err, q.dG_dec_err); vti = (p.dG_dec_TI - p.rocklin) - (q.dG_dec_TI - q.rocklin)
    dd.append(dict(motif=m, dG_corr_pb=p.dG_corr, dG_corr_docked=q.dG_corr, ddG=v, ddG_err=e, ddG_TI=vti,
                   ddG_TI_err=math.hypot(p.dG_dec_TI_err, q.dG_dec_TI_err), rocklin_diff=p.rocklin - q.rocklin, ddG_rounded=rnd(v, e),
                   contribution_c=-v, c_rounded=rnd(-v, e)))
DD = pd.DataFrame(dd); DD.to_csv(f"{out}_ddG.csv", index=False)
print(f"\n===== {sysn}: ddG(m) = dG_dec,corr(pb) - dG_dec,corr(docked)  [>0: motif interacts more favourably in pb]  c(m) = -ddG =====")
print(DD.round(3).to_string(index=False))
