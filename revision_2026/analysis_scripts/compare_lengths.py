#!/usr/bin/env python3
"""16 Nov 2026: per-leg and ddG comparison of three data sets (MBAR E+V, Rocklin-corrected, paired bootstrap SE):
  A = earlier 5 ns analysis (prod64v4 / cpd_v3; 54 windows except BHD2 legs, which already had the dense windows),
  B = first 5 ns of the final 63-window set (analysis/v5x63; isolates the effect of sampling length),
  C = 10 ns, 63 windows (prod64v10 / cpd_v10; the reported values).
B - A = effect of the schedule (extra dense windows), C - B = effect of sampling 5 -> 10 ns on the same schedule.
Run from state_2026-10-30/: python analysis/v5x63/compare_lengths.py  -> analysis/v5x63/compare_lengths.{csv,txt,tex}"""
import json, os, math, pandas as pd
A = "analysis"; RC = pd.read_csv("rocklin/rocklin_all_legs.tsv", sep="\t").set_index(["system", "leg"])["total_kcal"]
SETS = {"5ns_54win": {"64pp": "prod64v4", "cpd": "cpd_v3"}, "5ns_63win": {"64pp": "v5x63", "cpd": "v5x63"}, "10ns_63win": {"64pp": "prod64v10", "cpd": "cpd_v10"}}
MOT = ["lesion", "partner", "bhd2", "bhd3"]
def leg(s, d, l):
    f = f"{A}/{d}/{s}_{l}_2stage.json"
    if not os.path.exists(f): return None
    t = json.load(open(f))["total"]; return t["recommended"] - RC[(s, l)], t["recommended_err"]
rows, dd = [], []
for s in ("64pp", "cpd"):
    for m in MOT:
        r = dict(system=s, motif=m)
        for k, ds in SETS.items():
            p, q = leg(s, ds[s], f"pb_{m}"), leg(s, ds[s], f"docked_{m}")
            for st, v in (("pb", p), ("docked", q)):
                rows.append(dict(system=s, leg=f"{st}_{m}", set=k, dG_corr=v[0] if v else float("nan"), err=v[1] if v else float("nan")))
            r[k] = p[0] - q[0] if p and q else float("nan"); r[k + "_err"] = math.hypot(p[1], q[1]) if p and q else float("nan")
        r["sched_effect(B-A)"] = r["5ns_63win"] - r["5ns_54win"]; r["sampling_effect(C-B)"] = r["10ns_63win"] - r["5ns_63win"]
        dd.append(r)
L = pd.DataFrame(rows).pivot_table(index=["system", "leg"], columns="set", values="dG_corr")[list(SETS)]
L["C-B"] = L["10ns_63win"] - L["5ns_63win"]; L["B-A"] = L["5ns_63win"] - L["5ns_54win"]
D = pd.DataFrame(dd)
D.to_csv(f"{A}/v5x63/compare_lengths_ddG.csv", index=False); L.to_csv(f"{A}/v5x63/compare_lengths_legs.csv")
with open(f"{A}/v5x63/compare_lengths.txt", "w") as fh:
    fh.write(__doc__ + "\n\nddG (kcal/mol, bound - docked):\n" + D.round(2).to_string(index=False) + "\n\nper leg dG_dec,corr:\n" + L.round(2).to_string() + "\n")
print(open(f"{A}/v5x63/compare_lengths.txt").read())
# LaTeX rows (ddG): set | 54-win 5 ns | 63-win 5 ns | 63-win 10 ns | C-B
NM = {"64pp": r"\textcolor{dock}{\textbf{6-4PP}}", "cpd": r"\textcolor{cpdd}{\textbf{CPD}}"}; MN = {"lesion": "Lesion", "partner": "Partner", "bhd2": "BHD2", "bhd3": "BHD3"}
f = lambda v, e: f"{v:.1f} $\\pm$ {e:.1f}".replace("-", "$-$", 1) if v < 0 else f"{v:.1f} $\\pm$ {e:.1f}"
out = []
for s in ("64pp", "cpd"):
    for i, m in enumerate(MOT):
        r = D[(D.system == s) & (D.motif == m)].iloc[0]
        out.append(" & ".join([(r"\multirow{4}{*}{" + NM[s] + "}") if i == 0 else "", MN[m], f(r["5ns_54win"], r["5ns_54win_err"]), f(r["5ns_63win"], r["5ns_63win_err"]),
                               r"\textbf{" + f(r["10ns_63win"], r["10ns_63win_err"]) + "}", f"{r['sampling_effect(C-B)']:+.1f}".replace("-", "$-$")]) + r"\\")
    if s == "64pp": out.append(r"\midrule")
open(f"{A}/v5x63/compare_lengths_ddG.tex", "w").write("\n".join(out) + "\n")
