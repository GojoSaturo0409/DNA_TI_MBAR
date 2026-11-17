#!/usr/bin/env python3
"""Dynamics tables for the revision (main Table 2 and 3, SI Table S11), 6-4PP and CPD.
Inputs: mech_unbiased/welch/welch_{64pp,cpd}_perreplicate.csv (per-replica Welch tests, Holm over the 20 tests of each lesion).
B-factors: 6CFI (6-4PP) and 2QSG (CPD; lesion not resolved), heavy atoms, resolved residues, same selections as the MD motifs.
Rounding: uncertainty to 1 significant figure (2 if leading digit 1), value to the same decimal place."""
import math, os, pandas as pd
W = "/media/rahul/EXTERNAL_USB/binding_affinity/TI/mech_unbiased/welch"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "si_tables")
def dec(e):
    p = -int(math.floor(math.log10(e))); lead = int(e * 10 ** p + 1e-9)
    if lead == 1: p += 1
    return p
def pm(v, e):
    p = dec(e)
    if p <= 0:
        q = 10 ** (-p); return f"${int(round(v / q) * q)}\\pm{int(round(e / q) * q)}$".replace("$-", "$-")
    return f"${v:.{p}f}\\pm{e:.{p}f}$"
def pv(p):
    if p < 0.001: return "$<$0.001"
    if p >= 0.995: return "1.00"
    return f"{float(f'{p:.2g}'):g}"
ORDER = [("Lesion (602)", "lesion"), ("5' partner (584)", "5$'$ partner (584)"), ("3' partner (585)", "3$'$ partner (585)"),
         ("BHD2 tip", "BHD2 hairpin"), ("BHD3 tip", "BHD3 hairpin")]
MET = ["MSD (A^2, 25 ns)", "RMSF (A)", "RMSD (A)", "SASA (A^2)"]
BF = {"64pp": {"lesion": "214.1", "5$'$ partner (584)": "110.4", "3$'$ partner (585)": "125.6", "BHD2 hairpin": "140.7$^a$", "BHD3 hairpin": "129.8"},
      "cpd": {"lesion": "n.r.$^b$", "5$'$ partner (584)": "76.7", "3$'$ partner (585)": "72.9", "BHD2 hairpin": "71.6$^a$", "BHD3 hairpin": "73.2"}}
T = {s: pd.read_csv(f"{W}/welch_{s}_perreplicate.csv") for s in ("64pp", "cpd")}
def get(s, mot, met): r = T[s][(T[s].motif == mot) & (T[s].metric == met)].iloc[0]; return r

# ---- main Table 2: mean +- SEM, both lesions
rows = []
for s, lab in (("64pp", "6-4PP"), ("cpd", "CPD")):
    rows.append(r"\multicolumn{10}{@{}l}{\textbf{" + lab + r"}}\\")
    for mot, nm in ORDER:
        cells = []
        for met in MET:
            r = get(s, mot, met); cells += [pm(r.docked, r.docked_sem), pm(r.bound, r.bound_sem)]
        allc = [nm] + cells + [BF[s][nm]]
        if s == "cpd": allc = [r"\rev{" + c + "}" for c in allc]
        rows.append(" & ".join(allc) + r"\\")
    if s == "64pp": rows.append(r"\midrule")
open(f"{OUT}/tab_dyn_main.tex", "w").write(r"""\begin{tabular}{@{}lccccccccc@{}}
\toprule
& \multicolumn{2}{c}{\textbf{MSD (\AA\textsuperscript{2})}} & \multicolumn{2}{c}{\textbf{RMSF (\AA)}} & \multicolumn{2}{c}{\textbf{RMSD (\AA)}} & \multicolumn{2}{c}{\textbf{SASA (\AA\textsuperscript{2})}} & \textbf{B-factor (\AA\textsuperscript{2})}\\
\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-7}\cmidrule(lr){8-9}\cmidrule(l){10-10}
\textbf{Motif} & Docked & Bound & Docked & Bound & Docked & Bound & Docked & Bound & Bound\\
\midrule
""" + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n")

# ---- main Table 3: p (p_Holm), both lesions
rows = []
for mot, nm in ORDER:
    c = []
    for s in ("64pp", "cpd"):
        for met in MET:
            r = get(s, mot, met); x = f"{pv(r.p)} ({pv(r.p_holm)})"
            c.append(r"\textbf{" + x + "}" if r.p_holm < 0.05 else x)
    rows.append(f"{nm} & " + " & ".join(c) + r"\\")
open(f"{OUT}/tab_welch_main.tex", "w").write(r"""\begin{tabular}{@{}lllllllll@{}}
\toprule
 & \multicolumn{4}{c}{6-4PP} & \multicolumn{4}{c}{CPD}\\
\cmidrule(lr){2-5}\cmidrule(l){6-9}
Motif & MSD & RMSF & RMSD & SASA & MSD & RMSF & RMSD & SASA\\
\midrule
""" + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n")

# ---- SI Table S11: full, both lesions
rows = []
for s, lab in (("64pp", "6-4PP"), ("cpd", "CPD")):
    rows.append(r"\multicolumn{7}{@{}l}{\textbf{" + lab + r"}}\\")
    for mot, nm in ORDER:
        for k, met in enumerate(MET):
            r = get(s, mot, met)
            mn = met.replace("A^2", "\\AA$^2$").replace("(A)", "(\\AA)").replace(", 25 ns", "")
            first = f"{nm} (B {BF[s][nm]})" if k == 0 else ""
            rows.append(f"{first} & {mn} & {pm(r.docked, r.docked_sem)} & {pm(r.bound, r.bound_sem)} & ${r.t:+.1f}$ ({r.df:.1f}) & {pv(r.p)} & {pv(r.p_holm)}\\\\")
        rows.append(r"\midrule")
rows = rows[:-1]
open(f"{OUT}/tab_dyn_si.tex", "w").write(r"""\begin{tabular}{@{}llrrrrr@{}}
\toprule
Motif (B, \AA$^2$) & Metric & Docked & Bound & $t$ (df) & $p$ & $p_{\mathrm{Holm}}$\\\midrule
""" + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n")
print("wrote tab_dyn_main.tex, tab_welch_main.tex, tab_dyn_si.tex")
