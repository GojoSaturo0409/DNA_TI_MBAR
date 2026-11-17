#!/usr/bin/env python
"""17 Nov 2026: residue-level direct Coulomb interaction energy (kcal/mol, vacuum, no cutoff, solute only) between each
BHD2-tip residue (392-406) and every other solute residue (protein 1-565, DNA 566-620; excluding the tip itself),
per frame, averaged per replica. 2 lesions x 2 states x 4 replicas. Usage: elec_res.py <lesion> -> elec_res_<lesion>.csv"""
import sys, numpy as np, pandas as pd, MDAnalysis as mda, warnings; warnings.filterwarnings('ignore')
B = "/media/rahul/EXTERNAL_USB/binding_affinity/TI"; K = 332.0636
SYS = {"64pp": {"docked": (f"{B}/rahul_ti_unbiased/docked/system_docked_cleaned.prmtop", f"{B}/rahul_ti_unbiased/docked/dk_final"),
                "pb": (f"{B}/rahul_ti_unbiased/pb/system.prmtop", f"{B}/rahul_ti_unbiased/pb/pb_final")},
       "cpd": {"docked": (f"{B}/cpd_unbiased/analysis_mirror/docked/system_docked_cleaned.prmtop", f"{B}/cpd_unbiased/analysis_mirror/docked/dk_final"),
               "pb": (f"{B}/cpd_unbiased/analysis_mirror/pb/system.prmtop", f"{B}/cpd_unbiased/analysis_mirror/pb/pb_final")}}
s = sys.argv[1]; rows = []
for st, (top, pre) in SYS[s].items():
    for r in (1, 2, 3, 4):
        u = mda.Universe(top, f"{pre}_{r}.traj.nc", format='NCDF')
        tip = u.select_atoms('resid 392-406'); env = u.select_atoms('resid 1-620 and not resid 392-406')
        qt, qe = tip.charges, env.charges; rt, re_ = tip.resids, env.resids
        ut, ue = np.unique(rt), np.unique(re_)
        it = np.searchsorted(ut, rt); ie = np.searchsorted(ue, re_)
        acc = np.zeros((len(ut), len(ue))); n = 0
        for ts in u.trajectory:
            d = np.linalg.norm(tip.positions[:, None, :] - env.positions[None, :, :], axis=-1)
            E = K * qt[:, None] * qe[None, :] / d
            M = np.zeros((len(ut), len(qe))); np.add.at(M, it, E)
            P = np.zeros((len(ut), len(ue))); np.add.at(P.T, ie, M.T); acc += P; n += 1
        acc /= n
        for a, ra in enumerate(ut):
            for b, rb in enumerate(ue):
                if abs(acc[a, b]) > 0.05: rows.append(dict(lesion=s, state=st, rep=r, tip_res=ra, partner=rb, E=acc[a, b]))
        print(s, st, r, n, flush=True)
pd.DataFrame(rows).to_csv(f"elec_res_{s}.csv", index=False)
