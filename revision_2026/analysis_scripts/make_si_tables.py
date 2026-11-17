#!/usr/bin/env python3
"""LaTeX tables for the revised SI (16 Nov 2026). Run from rev/:  python scripts/make_si_tables.py <state_dir>
<state_dir> = TI/ti_rep1/state_2026-10-30 (analysis/prod64v10, analysis/cpd_v10, analysis/endstate/v10, analysis/v5x63,
rocklin/rocklin_all_legs.tsv, prod2c/lambdas_bhd2dense.csv). Writes si_tables/*.tex (complete tabular environments).
Rounding (editor): SE to 1 significant figure (2 if leading digit 1), value to the same decimal place."""
import sys, os, json, math, pandas as pd
ST = sys.argv[1]; A = f"{ST}/analysis"; OUT = "si_tables"; os.makedirs(OUT, exist_ok=True)
MOT = ["lesion", "partner", "bhd2", "bhd3"]; MN = {"lesion": "lesion", "partner": "partner bases", "bhd2": "BHD2 tip", "bhd3": "BHD3 tip"}
SYS = {"64pp": ("prod64v10", "prod64v4", "6-4PP"), "cpd": ("cpd_v10", "cpd_v3", "CPD")}
RC = pd.read_csv(f"{ST}/rocklin/rocklin_all_legs.tsv", sep="\t").set_index(["system", "leg"])
def m(x): return x.replace("-", "$-$", 1) if x.startswith("-") else x
def rnd(v, e):
    p = -int(math.floor(math.log10(e))); lead = int(e * 10 ** p + 1e-9)
    if lead == 1: p += 1
    p = max(p, 0); return m(f"{v:.{p}f}") + f" $\\pm$ {e:.{p}f}"
f1 = lambda v, d=1: m(f"{v:.{d}f}"); fs = lambda v, d=1: m(f"{v:+.{d}f}")
LEGS = lambda s: pd.DataFrame(json.load(open(f"{A}/{SYS[s][0]}/summary_{s}_legs.json"))).set_index("leg")
DDG = lambda s: pd.read_csv(f"{A}/{SYS[s][0]}/summary_{s}_ddG.csv").set_index("motif")
ORDER = [f"{st}_{mm}" for mm in MOT for st in ("pb", "docked")]
lab = lambda l: ("bound " if l.startswith("pb") else "docked ") + MN[l.split("_", 1)[1]]
def write(name, head, rows, foot=r"\bottomrule\end{tabular}"):
    open(f"{OUT}/{name}.tex", "w").write(head + "\n" + "\n".join(rows) + "\n" + foot + "\n")

# S2: lambda schedule
L = pd.read_csv(f"{ST}/prod2c/lambdas_bhd2dense.csv")
rows = []
E, V = L[L.stage == "E"].reset_index(drop=True), L[L.stage == "V"].reset_index(drop=True)
for i in range(max(len(E), len(V))):
    e = E.iloc[i] if i < len(E) else None; v = V.iloc[i] if i < len(V) else None
    c = [str(int(e.window)) if e is not None else "", f"{e['lambda']:.6f}" if e is not None else "", f"{e.weight:.6f}" if e is not None else "",
         str(int(v.window)) if v is not None else "", f"{2*v['lambda'] if False else v['lambda']:.6f}" if v is not None else "",
         (f"{v.weight:.6f}" if v.weight > 0 else "MBAR only") if v is not None else ""]
    rows.append(" & ".join(c) + r"\\")
write("tab_lambda", r"\begin{tabular}{@{}rrr@{\qquad}rrr@{}}\toprule"+"\n"+r"\multicolumn{3}{c}{stage E (\texttt{clambda})} & \multicolumn{3}{c}{stage V ($\lambda'$)}\\\cmidrule(r){1-3}\cmidrule(l){4-6}"+"\n"+r"window & $\lambda$ & weight & window & $\lambda'$ & weight\\\midrule", rows)

# S3: Rocklin
rows = []
for s in SYS:
    for k, mm in enumerate(MOT):
        a, b = RC.loc[(s, f"pb_{mm}")], RC.loc[(s, f"docked_{mm}")]
        rows.append(" & ".join([SYS[s][2] if k == 0 else "", MN[mm], f1(b.total_kcal, 2), f1(a.total_kcal, 2), fs(a.total_kcal - b.total_kcal, 2)]) + r"\\")
    if s == "64pp": rows.append(r"\midrule")
write("tab_rocklin", r"\begin{tabular}{@{}llrrr@{}}\toprule"+"\n"+r"Lesion & Motif & docked & bound & bound $-$ docked\\\midrule", rows)

# S4/S5: per-leg stage table
for s in SYS:
    D = LEGS(s); rows = []
    for l in ORDER:
        r = D.loc[l]
        rows.append(" & ".join([lab(l), f1(r.E_TI), rnd(r.E_MBAR, r.E_err), f"{r.E_overlap_min:.3f}", f1(r.V_TI_spline), rnd(r.V_MBAR, r.V_err),
                                f"{r.V_overlap_min:.3f}", f1(r.rocklin, 2), rnd(r.dG_corr, r.dG_dec_err)]) + r"\\")
    write(f"tab_legs_{s}", r"\begin{tabular}{@{}lrrrrrrrr@{}}\toprule"+"\n"+r"& \multicolumn{3}{c}{stage E} & \multicolumn{3}{c}{stage V} & & \\\cmidrule(lr){2-4}\cmidrule(lr){5-7}"+"\n"+r"Leg & TI & MBAR & min ovl & TI & MBAR & min ovl & Rocklin & $\Delta G^{\mathrm{corr}}_{\mathrm{dec}}$\\\midrule", rows)

# S6: ddG TI vs MBAR
rows = []
for s in SYS:
    D, Lg = DDG(s), LEGS(s)
    for k, mm in enumerate(MOT):
        d = D.loc[mm]; a, b = Lg.loc[f"pb_{mm}"], Lg.loc[f"docked_{mm}"]
        se = math.sqrt(a.E_diff_err**2 + a.V_diff_err**2 + b.E_diff_err**2 + b.V_diff_err**2)
        rows.append(" & ".join([SYS[s][2] if k == 0 else "", MN[mm], r"\textbf{" + rnd(d.ddG, d.ddG_err) + "}", rnd(d.ddG_TI, d.ddG_TI_err),
                                fs(d.ddG_TI - d.ddG, 2) + f" $\\pm$ {se:.2f}", f"{abs(d.ddG_TI - d.ddG)/se:.1f}"]) + r"\\")
    if s == "64pp": rows.append(r"\midrule")
write("tab_ti_mbar", r"\begin{tabular}{@{}llrrrr@{}}\toprule"+"\n"+r"Lesion & Motif & $\Delta\Delta G$ (MBAR) & $\Delta\Delta G$ (TI) & TI $-$ MBAR & $|\Delta|/\sigma$\\\midrule", rows)

# S7: convergence per leg (both lesions)
rows = []
for s in SYS:
    D = LEGS(s)
    for k, l in enumerate(ORDER):
        r = D.loc[l]; t5 = json.load(open(f"{A}/{SYS[s][1]}/{s}_{l}_2stage.json"))["total"]["recommended"]
        rows.append(" & ".join([SYS[s][2] if k == 0 else "", lab(l), f1(t5), f1(r.dG_dec), fs(r.dG_dec - t5), f1(r.half1), f1(r.half2), fs(r.fwd[0] - r.rev[0])]) + r"\\")
    if s == "64pp": rows.append(r"\midrule")
write("tab_conv", r"\begin{tabular}{@{}llrrrrrr@{}}\toprule"+"\n"+r"Lesion & Leg & 5\,ns & 10\,ns & change & 1st half & 2nd half & fwd$-$rev (25\,\%)\\\midrule", rows)

# S8: 5 ns 54 win vs 5 ns 63 win vs 10 ns
C = pd.read_csv(f"{A}/v5x63/compare_lengths_ddG.csv"); rows = []
for s in SYS:
    for k, mm in enumerate(MOT):
        r = C[(C.system == s) & (C.motif == mm)].iloc[0]
        rows.append(" & ".join([SYS[s][2] if k == 0 else "", MN[mm], rnd(r["5ns_54win"], r["5ns_54win_err"]), rnd(r["5ns_63win"], r["5ns_63win_err"]),
                                r"\textbf{" + rnd(r["10ns_63win"], r["10ns_63win_err"]) + "}", fs(r["sched_effect(B-A)"]), fs(r["sampling_effect(C-B)"])]) + r"\\")
    if s == "64pp": rows.append(r"\midrule")
write("tab_lengths", r"\begin{tabular}{@{}llrrrrr@{}}\toprule"+"\n"+r"Lesion & Motif & A: 5\,ns, 54 win & B: 5\,ns, 63 win & C: 10\,ns, 63 win & B$-$A & C$-$B\\\midrule", rows)

# S9: lesion effect (CPD - 6-4PP) per state
L64, LC = LEGS("64pp"), LEGS("cpd"); rows = []
for mm in MOT:
    c = []
    for st in ("docked", "pb"):
        x, y = LC.loc[f"{st}_{mm}"], L64.loc[f"{st}_{mm}"]; c.append((x.dG_corr - y.dG_corr, math.hypot(x.dG_dec_err, y.dG_dec_err)))
    eff = c[1][0] - c[0][0]; ee = math.hypot(c[0][1], c[1][1])
    rows.append(" & ".join([MN[mm], rnd(*c[0]), rnd(*c[1]), rnd(eff, ee)]) + r"\\")
write("tab_lesion_effect", r"\begin{tabular}{@{}lrrr@{}}\toprule"+"\n"+r"Motif & $\Delta G^{\mathrm{corr}}_{\mathrm{dec}}$(CPD)$-$(6-4PP), docked & same, bound & effect on $\Delta\Delta G$\\\midrule", rows)

# S10: end-state per stage
rows = []
for s, fn in (("64pp", "endstate64_ddG.csv"), ("cpd", "endstate_cpd_ddG.csv")):
    E = pd.read_csv(f"{A}/endstate/v10/{fn}").set_index("motif")
    for k, mm in enumerate(MOT):
        r = E.loc[mm]
        rows.append(" & ".join([SYS[s][2] if k == 0 else "", MN[mm], rnd(r.dddG_E, r.dddG_E_err), rnd(r.dddU_E, r.dddU_E_err), rnd(r.ddmTdS_E, r.ddmTdS_E_err),
                                rnd(r.dddG_V, r.dddG_V_err), rnd(r.dddU_V, r.dddU_V_err), rnd(r.ddmTdS_V, r.ddmTdS_V_err)]) + r"\\")
    if s == "64pp": rows.append(r"\midrule")
write("tab_endstate", r"\begin{tabular}{@{}llrrrrrr@{}}\toprule"+"\n"+r"& & \multicolumn{3}{c}{stage E} & \multicolumn{3}{c}{stage V}\\\cmidrule(lr){3-5}\cmidrule(l){6-8}"+"\n"+r"Lesion & Motif & $\Delta\Delta G$ & $\Delta\Delta U$ & $-T\Delta\Delta S$ & $\Delta\Delta G$ & $\Delta\Delta U$ & $-T\Delta\Delta S$\\\midrule", rows)
print("wrote", sorted(os.listdir(OUT)))
