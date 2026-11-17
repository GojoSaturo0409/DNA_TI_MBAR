#!/usr/bin/env python3
"""2026-11-16: build runs_cpd10/<leg> = runs_cpd/<leg> (5 ns s0/sp, w01-54) + ext s3 files + dense9 direct-10ns w55-63,
gathered from dlcpd10/<account>/out (7 research + <account>). Asserts one source per file, 10000 pts per window."""
import sys, glob, os, shutil
legs = sys.argv[1:]
for leg in legs:
    od = f"runs_cpd10/{leg}"; os.makedirs(od, exist_ok=True)
    for f in glob.glob(f"runs_cpd/{leg}/w*"):
        if int(os.path.basename(f)[1:3]) <= 54: shutil.copy2(f, od)
    src = {}
    for f in glob.glob(f"dlcpd10/*/out/ext/cpd/{leg}/*") + glob.glob(f"dlcpd10/*/out/dense9/cpd/{leg}/*"):
        b = os.path.basename(f); assert b not in src, (b, src.get(b), f); src[b] = f
    for b, f in src.items(): shutil.copy2(f, od)
    for w in range(1, 64):
        W = f"w{w:02d}"; n = 0
        for g in (f"{od}/{W}_dvdl.txt", f"{od}/{W}_s3_dvdl.txt"):
            if os.path.exists(g): n += sum(1 for l in open(g) if l.strip())
        assert n == 10000, (leg, W, n)
    print(leg, "OK 63 x 10000", len(os.listdir(od)), "files")
