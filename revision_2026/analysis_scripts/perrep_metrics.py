#!/usr/bin/env python
"""2026-11-16 (R1.7): per-replicate RMSF / RMSD of each motif from the 4x50 ns unbiased runs (64pp: rahul_ti_unbiased; CPD: cpd_unbiased),
for several definitions, to find the one behind manuscript Table S2 and to run Welch tests on per-replicate values.
usage: perrep_metrics.py <docked top> <docked traj glob prefix> <pb top> <pb prefix> <out csv>"""
import sys, numpy as np, MDAnalysis as mda, pandas as pd, warnings; warnings.filterwarnings('ignore')
from MDAnalysis.analysis import align, rms
dtop, dpre, ptop, ppre, out = sys.argv[1:6]
MOT = {'lesion': 'resid 602', 'p585': 'resid 585', 'p584': 'resid 584', 'bhd2': 'resid 392-406', 'bhd3': 'resid 475-485'}
FITS = {'protdna_heavy': 'resid 1-620 and not name H*', 'protdna_ca_p': '(resid 1-620) and (name CA or name P)'}
rows = []
for st, top, pre in (('docked', dtop, dpre), ('bound', ptop, ppre)):
    for r in (1, 2, 3, 4):
        u = mda.Universe(top, f"{pre}_{r}.traj.nc", format='NCDF')
        for fn, fsel in FITS.items():
            ref = mda.Universe(top, f"{pre}_{r}.traj.nc", format='NCDF'); ref.trajectory[0]
            align.AlignTraj(u, ref, select=fsel, in_memory=True).run()
            for m, sel in MOT.items():
                for an, asel in (('heavy', f'{sel} and not name H*'), ('all', sel)):
                    ag = u.select_atoms(asel); X = np.array([ag.positions.copy() for ts in u.trajectory])
                    rmsf = np.sqrt(((X - X.mean(0)) ** 2).sum(-1).mean(0)).mean()
                    rmsd_t = np.sqrt(((X - X[0]) ** 2).sum(-1).mean(-1))
                    rows.append(dict(state=st, rep=r, fit=fn, atoms=an, motif=m, rmsf=rmsf, rmsd_mean=rmsd_t.mean(), rmsd_end=rmsd_t[-1], n=len(X)))
            u = mda.Universe(top, f"{pre}_{r}.traj.nc", format='NCDF')
        print(st, r, flush=True)
pd.DataFrame(rows).to_csv(out, index=False)
