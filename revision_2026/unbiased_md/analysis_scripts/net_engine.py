#!/usr/bin/env python
"""Per-frame interaction data for fingerprints, residue networks, clustering, time-resolved contacts and DCC (11 Nov 2026).
One job = one unbiased replica (full system, with H and water; 125 frames, 400 ps).
Probes: BHD2 392-406, BHD3 475-485, partner bases 584-585, lesion 602.  Partners: the 209 'site' residues (maps/top/site_resids.txt).
Interaction types (per frame, per probe residue x site residue; distances with minimum image):
  contact   any heavy-atom pair <= 4.5 A
  hbond     polar H ... acceptor (O, or N without H) <= 2.5 A and angle D-H...A >= 120 deg
  saltbr    Lys NZ / Arg NE NH1 NH2 (cationic N)  to  Asp OD1 OD2 / Glu OE1 OE2 / DNA phosphate O (OP1 OP2; lesion O1 O2 O9 O10) <= 4.0 A
  hydroph   apolar carbon (C with no bonded N/O) pairs <= 4.0 A
  dna_bb    probe heavy atom to DNA sugar-phosphate heavy atom <= 4.5 A   (only for protein probes; DNA partner residues)
  dna_base  probe heavy atom to DNA base heavy atom <= 4.5 A
  wbridge   one water whose O is <= 3.5 A from a polar (N/O) atom of both residues
Also stored: residue-residue heavy-atom contact counts for all site pairs (network), representative atoms of every site residue
(CA; DNA C1'; lesion C5) for DCC, and BHD2 / BHD3 heavy atoms superposed on the per-state DNA fit core (maps/top/fit_core.json) for clustering.
usage: net_engine.py LESION STATE REP OUTDIR"""
import sys, os, json, numpy as np, MDAnalysis as mda, warnings
from MDAnalysis.lib.distances import capped_distance, calc_angles
warnings.filterwarnings("ignore")
LES, ST, REP, OUT = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]; os.makedirs(OUT, exist_ok=True)
NF = int(os.environ.get("NFRAMES", "0"))
TI = "/media/rahul/EXTERNAL_USB/binding_affinity/TI"; MU = f"{TI}/mech_unbiased"; BA = os.path.dirname(TI)
ROOT = {"64pp": f"{TI}/rahul_ti_unbiased", "cpd": f"{TI}/cpd_unbiased/analysis_mirror"}[LES]
TOP = f"{ROOT}/docked/system_docked_cleaned.prmtop" if ST == "docked" else f"{ROOT}/pb/system.prmtop"
TRJ = f"{ROOT}/docked/dk_final_{REP}.traj.nc" if ST == "docked" else f"{ROOT}/pb/pb_final_{REP}.traj.nc"
SITE = np.loadtxt(f"{MU}/maps/top/site_resids.txt", dtype=int); SIDX = {r: i for i, r in enumerate(SITE)}
PROBE = list(range(392, 407)) + list(range(475, 486)) + [584, 585, 602]; PIDX = {r: i for i, r in enumerate(PROBE)}
u = mda.Universe(TOP, TRJ, format="NCDF"); A = u.atoms
el = np.array([n[0] for n in A.names]); el[np.isin(A.names, ["Na+", "NA", "Cl-"])] = "X"
resid = A.resids; isH = el == "H"
site_mask = np.isin(resid, SITE); heavy = site_mask & ~isH
nb = [[] for _ in range(len(A))]
for b in u.bonds: i, j = b.indices; nb[i].append(j); nb[j].append(i)
polarH = np.array([i for i in np.where(site_mask & isH)[0] if any(el[j] in "NO" for j in nb[i])])
donor_of = {i: next(j for j in nb[i] if el[j] in "NO") for i in polarH}
acc = np.array([i for i in np.where(heavy & np.isin(el, ["O", "N"]))[0] if el[i] == "O" or not any(el[j] == "H" for j in nb[i])])
cat = np.where(site_mask & (((A.resnames == "LYS") & (A.names == "NZ")) | ((A.resnames == "ARG") & np.isin(A.names, ["NE", "NH1", "NH2"]))))[0]
ani = np.where(site_mask & (((A.resnames == "ASP") & np.isin(A.names, ["OD1", "OD2"])) | ((A.resnames == "GLU") & np.isin(A.names, ["OE1", "OE2"])) |
                            ((resid >= 566) & np.isin(A.names, ["OP1", "OP2"])) | ((resid == 602) & np.isin(A.names, ["O1", "O2", "O9", "O10"]))))[0]
apol = np.array([i for i in np.where(heavy & (el == "C"))[0] if not any(el[j] in "NO" for j in nb[i])])
SP_NAMES = {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'", "P1", "P2", "O1", "O2", "O9", "O10"}
dna_heavy = np.where(heavy & (resid >= 566))[0]; dna_sp = dna_heavy[np.isin(A.names[dna_heavy], list(SP_NAMES))]; dna_base = np.setdiff1d(dna_heavy, dna_sp)
probe_heavy = np.where(heavy & np.isin(resid, PROBE))[0]; prot_probe_heavy = probe_heavy[resid[probe_heavy] <= 565]
polar_site = np.where(heavy & np.isin(el, ["N", "O"]))[0]; polar_probe = polar_site[np.isin(resid[polar_site], PROBE)]
wat_O = u.select_atoms("resname WAT and name O").indices
hv = np.where(heavy)[0]
rep_idx = []
for r in SITE:
    nm = "CA" if r <= 565 else ("C5" if r == 602 else "C1'")
    rep_idx.append(u.select_atoms(f"resid {r} and name {nm}").indices[0])
rep_idx = np.array(rep_idx)
CORE = [(int(k.split(":")[0]), k.split(":")[1]) for k in json.load(open(f"{MU}/maps/top/fit_core.json"))[ST]["core_keys"]]
ref = mda.Universe(f"{TI}/rahul_ti_unbiased/{'pb/system.prmtop' if ST == 'pb' else 'docked/system_docked_cleaned.prmtop'}",
                   f"{BA}/{'system.inpcrd' if ST == 'pb' else 'system_docked_cleaned.inpcrd'}", format="INPCRD")
core_i = np.array([u.select_atoms(f"resid {r} and name {n}").indices[0] for r, n in CORE])
core_ref = np.array([ref.select_atoms(f"resid {r} and name {n}").positions[0] for r, n in CORE], float)
b2_i = np.where(heavy & (resid >= 392) & (resid <= 406))[0]; b3_i = np.where(heavy & (resid >= 475) & (resid <= 485))[0]
def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0); U, s, Vt = np.linalg.svd((P - pc).T @ (Q - qc)); d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T; return R, qc - R @ pc
TYPES = ["contact", "hbond", "saltbr", "hydroph", "dna_bb", "dna_base", "wbridge"]
nfr = len(u.trajectory) if NF == 0 else NF
FP = {t: np.zeros((nfr, len(PROBE), len(SITE)), np.uint8) for t in TYPES}
NET = np.zeros((len(SITE), len(SITE)), np.int32); REPX = np.zeros((nfr, len(SITE), 3), np.float32)
B2 = np.zeros((nfr, len(b2_i), 3), np.float32); B3 = np.zeros((nfr, len(b3_i), 3), np.float32)
def put(f, t, ia, ib):
    if len(ia) == 0: return
    ra, rb = resid[ia], resid[ib]; keep = ra != rb; ra, rb = ra[keep], rb[keep]
    for x, y in ((ra, rb), (rb, ra)):
        m = np.isin(x, PROBE)
        if m.any(): FP[t][f, [PIDX[v] for v in x[m]], [SIDX[v] for v in y[m]]] = 1
for f, ts in enumerate(u.trajectory[:nfr]):
    X = A.positions; box = ts.dimensions
    p = capped_distance(X[hv], X[hv], 4.5, box=box, return_distances=False); ia, ib = hv[p[:, 0]], hv[p[:, 1]]
    ra, rb = resid[ia], resid[ib]; k = ra < rb
    pairs = np.unique(np.stack([ra[k], rb[k]], 1), axis=0)
    NET[[SIDX[a] for a in pairs[:, 0]], [SIDX[b] for b in pairs[:, 1]]] += 1
    put(f, "contact", ia, ib)
    p = capped_distance(X[polarH], X[acc], 2.5, box=box, return_distances=False)
    if len(p):
        h = polarH[p[:, 0]]; a = acc[p[:, 1]]; d = np.array([donor_of[i] for i in h])
        ang = np.degrees(calc_angles(X[d], X[h], X[a], box=box)); ok = ang >= 120; put(f, "hbond", d[ok], a[ok])
    p = capped_distance(X[cat], X[ani], 4.0, box=box, return_distances=False)
    if len(p): put(f, "saltbr", cat[p[:, 0]], ani[p[:, 1]])
    p = capped_distance(X[apol], X[apol], 4.0, box=box, return_distances=False)
    if len(p): put(f, "hydroph", apol[p[:, 0]], apol[p[:, 1]])
    for t, grp in (("dna_bb", dna_sp), ("dna_base", dna_base)):
        p = capped_distance(X[prot_probe_heavy], X[grp], 4.5, box=box, return_distances=False)
        if len(p): put(f, t, prot_probe_heavy[p[:, 0]], grp[p[:, 1]])
    # water bridges: waters near probe polar atoms, then all site polar atoms near those waters
    p = capped_distance(X[polar_probe], X[wat_O], 3.5, box=box, return_distances=False)
    if len(p):
        W = np.unique(wat_O[p[:, 1]]); q = capped_distance(X[W], X[polar_site], 3.5, box=box, return_distances=False)
        if len(q):
            wr = {}
            for w, s in zip(W[q[:, 0]], polar_site[q[:, 1]]): wr.setdefault(w, set()).add(resid[s])
            ia2, ib2 = [], []
            for rs in wr.values():
                rs = sorted(rs)
                for i_ in range(len(rs)):
                    for j_ in range(i_ + 1, len(rs)):
                        if rs[i_] in PIDX or rs[j_] in PIDX:
                            FP["wbridge"][f, PIDX[rs[i_]] if rs[i_] in PIDX else PIDX[rs[j_]], SIDX[rs[j_]] if rs[i_] in PIDX else SIDX[rs[i_]]] = 1
                            if rs[i_] in PIDX and rs[j_] in PIDX: FP["wbridge"][f, PIDX[rs[j_]], SIDX[rs[i_]]] = 1
    REPX[f] = X[rep_idx]
    R, t = kabsch(X[core_i].astype(float), core_ref); B2[f] = X[b2_i] @ R.T + t; B3[f] = X[b3_i] @ R.T + t
np.savez_compressed(f"{OUT}/{LES}_{ST}_r{REP}.npz", **{f"fp_{t}": FP[t] for t in TYPES}, net=NET, nfr=nfr, repx=REPX, b2=B2, b3=B3,
                    probe=np.array(PROBE), site=SITE, b2_names=np.array([f"{resid[i]}:{A.names[i]}" for i in b2_i]),
                    b3_names=np.array([f"{resid[i]}:{A.names[i]}" for i in b3_i]))
print(LES, ST, REP, "frames", nfr, {t: int(FP[t].sum()) for t in TYPES})
