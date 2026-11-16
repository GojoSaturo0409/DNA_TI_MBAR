#!/usr/bin/env python3
"""2026-11-16 (user: 'calculate all using 63 windows using only 5 ns and compare'):
runs<sys>_5x63/<leg>: the 63-window set truncated to the FIRST 5 ns of every window.
 - extended windows: drop s3/sp3 files (keep s0-s2 = first 5 ns)
 - windows run 10 ns directly (dense w55-63 of non-BHD2 legs, 10000 DV/DL pts): keep first 5000 pts and first 51 sp frames
   (a 5 ns window has 51 sp frames at 80 ps; a direct 10 ns window 113)."""
import os, glob, shutil, sys
NFR = 51
def build(src_dirs, out):
    os.makedirs(out, exist_ok=True)
    for d in src_dirs:
        for f in glob.glob(f"{d}/w*"):
            b = os.path.basename(f)
            if "_s3" in b or "sp3" in b: continue
            shutil.copy2(f, f"{out}/{b}")
    for w in range(1, 64):
        dv = f"{out}/w{w:02d}_dvdl.txt"; L = [l for l in open(dv) if l.strip()]
        if len(L) == 10000:
            open(dv, "w").writelines(L[:5000]); sp = f"{out}/w{w:02d}_sp_energies.txt"
            if os.path.exists(sp):
                keep = [l for l in open(sp) if l.startswith("#") or (l.strip() and int(float(l.split()[0])) <= NFR)]
                open(sp, "w").writelines(keep)
        n = sum(1 for l in open(dv) if l.strip()); assert n == 5000, (out, w, n)
    print(out, "OK 63 x 5000")
legs = [f"{s}_{m}" for s in ("pb", "docked") for m in ("lesion", "partner", "bhd2", "bhd3")]
for l in legs: build([f"runs64_10/{l}"], f"runs64_5x63/{l}")
for l in legs:
    if os.path.isdir(f"runs_cpd10/{l}"): build([f"runs_cpd10/{l}"], f"runscpd_5x63/{l}")
    else: build([f"runs_cpd/{l}"] + glob.glob(f"dlcpd10/*/out/dense9/cpd/{l}"), f"runscpd_5x63/{l}")
