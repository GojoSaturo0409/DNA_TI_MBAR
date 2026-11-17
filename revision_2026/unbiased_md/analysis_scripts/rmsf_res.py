#!/usr/bin/env python
"""17 Nov 2026: per-residue RMSF of BHD2 (392-406) and BHD3 (475-485) tips, all atoms, after fitting protein CA + DNA P
(same definition as the per-replica RMSF of SI Table S11). 4 systems x 4 replicas. -> rmsf_res_<sys>.csv"""
import sys, numpy as np, pandas as pd, MDAnalysis as mda, warnings; warnings.filterwarnings('ignore')
from MDAnalysis.analysis import align
B = "/media/rahul/EXTERNAL_USB/binding_affinity/TI"
SYS = {"64pp": {"docked": (f"{B}/rahul_ti_unbiased/docked/system_docked_cleaned.prmtop", f"{B}/rahul_ti_unbiased/docked/dk_final"),
                "pb": (f"{B}/rahul_ti_unbiased/pb/system.prmtop", f"{B}/rahul_ti_unbiased/pb/pb_final")},
       "cpd": {"docked": (f"{B}/cpd_unbiased/analysis_mirror/docked/system_docked_cleaned.prmtop", f"{B}/cpd_unbiased/analysis_mirror/docked/dk_final"),
               "pb": (f"{B}/cpd_unbiased/analysis_mirror/pb/system.prmtop", f"{B}/cpd_unbiased/analysis_mirror/pb/pb_final")}}
s = sys.argv[1]; rows = []
for st, (top, pre) in SYS[s].items():
    for r in (1, 2, 3, 4):
        u = mda.Universe(top, f"{pre}_{r}.traj.nc", format='NCDF'); ref = mda.Universe(top, f"{pre}_{r}.traj.nc", format='NCDF'); ref.trajectory[0]
        align.AlignTraj(u, ref, select='(resid 1-620) and (name CA or name P)', in_memory=True).run()
        ag = u.select_atoms('resid 392-406 or resid 475-485'); X = np.array([ag.positions.copy() for ts in u.trajectory])
        f = np.sqrt(((X - X.mean(0)) ** 2).sum(-1).mean(0))
        for rid in sorted(set(ag.resids)):
            rows.append(dict(state=st, rep=r, resid=rid, rmsf=f[ag.resids == rid].mean()))
        print(s, st, r, flush=True)
pd.DataFrame(rows).to_csv(f"{B}/mech_unbiased/rmsf_res/rmsf_res_{s}.csv", index=False)
