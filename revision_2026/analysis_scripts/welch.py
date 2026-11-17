#!/usr/bin/env python
"""2026-11-16 R1.7: Welch t-tests bound vs docked on PER-REPLICATE values (n=4 per state), 64pp.
RMSF/RMSD: perrep_64pp.csv, fit protein CA + DNA P, all motif atoms, 125 frames (reproduces SI Table S2 within ~0.1 A).
MSD: tip/base COM MSD at tau = 25 ns (slides/data/64pp/*_msd_tau_perrep.dat; rebuilt pipeline). SASA: analysis_out/<state>/per_replicate_means.csv."""
import numpy as np, pandas as pd
from scipy import stats
B = "<PROJECT_ROOT>/TI"; D = f"{B}/ti_rep1/state_2026-10-30/analysis/slides/data/64pp"
pr = pd.read_csv("perrep_64pp.csv"); pr = pr[(pr.fit == "protdna_ca_p") & (pr.atoms == "all")]
def msd(f):
    X = np.loadtxt(f); i = np.argmin(np.abs(X[:, 0] - 25.0)); assert abs(X[i, 0] - 25) < 0.21; return X[i, 1:5]
MSDF = {"bhd2": "com_timeseries/BHD2_BHD3_bhd2_{s}_msd_tau_perrep.dat", "bhd3": "com_timeseries/BHD2_BHD3_bhd3_{s}_msd_tau_perrep.dat",
        "p585": "com_timeseries/partner585_partner584_partner585_{s}_msd_tau_perrep.dat", "p584": "com_timeseries/partner585_partner584_partner584_{s}_msd_tau_perrep.dat",
        "lesion": "lesion/lesion602_{s}_msd_tau_perrep.dat"}
SAS = {"bhd2": "BHD2_tip_sasa_A2", "bhd3": "BHD3_tip_sasa_A2", "p585": "partner_585_sasa_A2", "p584": "partner_584_sasa_A2", "lesion": "lesion_whole_sasa_A2"}
sa = {s: pd.read_csv(f"{B}/rahul_ti_unbiased/analysis_out/{s}/per_replicate_means.csv") for s in ("docked", "pb")}
NAME = {"bhd3": "BHD3 tip", "bhd2": "BHD2 tip", "lesion": "Lesion (602)", "p585": "3' partner (585)", "p584": "5' partner (584)"}
rows = []
for m in ("bhd3", "bhd2", "lesion", "p585", "p584"):
    vals = {"MSD (A^2, 25 ns)": (msd(f"{D}/{MSDF[m].format(s='docked')}"), msd(f"{D}/{MSDF[m].format(s='pb')}")),
            "RMSF (A)": (pr[(pr.motif == m) & (pr.state == "docked")].sort_values("rep").rmsf.values, pr[(pr.motif == m) & (pr.state == "bound")].sort_values("rep").rmsf.values),
            "RMSD (A)": (pr[(pr.motif == m) & (pr.state == "docked")].sort_values("rep").rmsd_mean.values, pr[(pr.motif == m) & (pr.state == "bound")].sort_values("rep").rmsd_mean.values),
            "SASA (A^2)": (sa["docked"][SAS[m]].values, sa["pb"][SAS[m]].values)}
    for k, (a, b) in vals.items():
        assert len(a) == 4 and len(b) == 4
        t, p = stats.ttest_ind(b, a, equal_var=False)
        va, vb = a.var(ddof=1) / 4, b.var(ddof=1) / 4; df = (va + vb) ** 2 / (va ** 2 / 3 + vb ** 2 / 3)
        rows.append(dict(motif=NAME[m], metric=k, docked=a.mean(), docked_sem=a.std(ddof=1) / 2, bound=b.mean(), bound_sem=b.std(ddof=1) / 2, t=t, df=df, p=p))
R = pd.DataFrame(rows); R["p_holm"] = np.minimum(1, np.maximum.accumulate((len(R) - np.arange(len(R))) * np.sort(R.p.values)))[np.argsort(np.argsort(R.p.values))]
R.to_csv("welch_64pp_perreplicate.csv", index=False); print(R.round(3).to_string())
