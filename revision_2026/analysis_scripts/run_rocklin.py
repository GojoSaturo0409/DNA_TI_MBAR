"""Rocklin correction (rocklinc + APBS) for one leg. Ions (Na+) are removed before building the PQRs (treated as solvent);
protein_netq = net charge of the solute without the mutated region; lig = mutated region (:602 here).
Usage: run_rocklin.py <prmtop> <rst7> <ligand resid selection> <outdir>"""
import sys, os, numpy as np, parmed as pmd, MDAnalysis as mda, rocklinc, warnings; warnings.filterwarnings("ignore")
top, crd, ligsel, out = sys.argv[1:5]; os.makedirs(out, exist_ok=True); os.chdir(out)
s = pmd.load_file(top, crd); box = list(s.box[:3])
radii = np.array([a.solvent_radius for a in s.atoms]); pos = np.array(s.coordinates)
u0 = mda.Universe(top, format="PARM7"); u0.load_new(pos[None, :, :], order="fac")
u0.add_TopologyAttr("radii", radii); u0.dimensions = list(box) + [90, 90, 90]
assert np.allclose(u0.atoms.charges, [a.charge for a in s.atoms]) and radii.min() > 0
u = mda.Merge(u0.select_atoms("not resname Na+ Cl- K+")); u.dimensions = list(box) + [90, 90, 90]
lig = u.select_atoms(ligsel); rest = u.select_atoms(f"not ({ligsel}) and not resname WAT")
qL, qP = round(float(lig.charges.sum()), 4), round(float(rest.charges.sum()), 4)
print(f"box {np.round(box,3)} A, ligand {ligsel}: {len(lig)} atoms q={qL}; rest of solute: {len(rest)} atoms q={qP}; waters {len(u.select_atoms('resname WAT').residues)}", flush=True)
rc = rocklinc.RocklinCorrection(box, qL, qP, 300.0, rocklinc.waters.TIP3P)
rc.make_APBS_input(u, ligsel, solvent_selection="resname WAT")
rc.run_APBS(apbs_exe="<ENV_DIR>/rocklinc/bin/apbs", apbs_in="apbs.in")
rc.read_APBS(); res = rc.compute(); rc.write("correction.txt")
print(open("correction.txt").read())
