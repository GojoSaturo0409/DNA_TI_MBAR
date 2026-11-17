#!/usr/bin/env python3
"""Mean crystallographic B-factor (heavy atoms, resolved atoms only, per chain) on exactly the simulated motifs.
6CFI: Rad4-6-4PP (bound 6-4PP complex); 2QSG: Rad4 on CPD-mismatch DNA (CPD unresolved).
Usage: python bfactors.py 6CFI.pdb 2QSG.pdb   (PDB files from https://www.rcsb.org)"""
import sys
SEL = [("lesion 602", "Y", 9, 9), ("BHD2 tip (Rad4 510-524)", "A", 510, 524), ("BHD3 tip (Rad4 593-603)", "A", 593, 603),
       ("5' partner 584 (W16)", "W", 16, 16), ("3' partner 585 (W17)", "W", 17, 17)]
for pdb in sys.argv[1:]:
    print(pdb)
    for nm, ch, a, b in SEL:
        v = [float(l[60:66]) for l in open(pdb) if l.startswith(("ATOM", "HETATM")) and l[21] == ch and a <= int(l[22:26]) <= b
             and (l[76:78].strip() or l[12:16].strip()[0]) != "H"]
        print(f"  {nm:28s} {sum(v)/len(v):7.1f}  ({len(v)} atoms)" if v else f"  {nm:28s}   not resolved")
