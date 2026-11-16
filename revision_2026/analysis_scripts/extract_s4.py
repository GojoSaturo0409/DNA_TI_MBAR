#!/usr/bin/env python3
"""(s4 copy of extract_s3.py, 2026-11-14) stdlib only. DV/DL (TI region 1) of extension segments: <root>/<sys>/<leg>/runs/wNN_s4.out.gz (with wNN_s4.done)
-> <out>/<sys>/<leg>/wNN_s4_dvdl.txt (seg time dvdl) + copies of wNN_sp4_energies.txt, wNN_s4.check. Same parsing as extract_dvdl.py.
Usage: extract_s3.py <root e.g. ~/prod2c/ext3> <outdir> <sys>"""
import sys, os, re, gzip, glob, shutil
root, out, sysn = sys.argv[1:4]
re_dv = re.compile(r"^ DV/DL\s*=\s*(-?[\d.]+)"); re_t = re.compile(r"TIME\(PS\)\s*=\s*([\d.]+)")
r1 = re.compile(r"\|\s*TI region\s+1\b"); r2 = re.compile(r"\|\s*TI region\s+2\b")
for done in sorted(glob.glob(f"{root}/{sysn}/*/runs/w*_s4.done")):
    runs = os.path.dirname(done); leg = runs.split("/")[-2]; w = os.path.basename(done)[:3]
    od = f"{out}/{sysn}/{leg}"; os.makedirs(od, exist_ok=True); of = f"{od}/{w}_s4_dvdl.txt"
    if os.path.exists(of): continue
    rows = []; t = None; ti1 = False
    for line in gzip.open(f"{runs}/{w}_s4.out.gz", "rt", errors="replace"):
        if "A V E R A G E S" in line: break
        if r1.search(line): ti1 = True; continue
        if r2.search(line): ti1 = False; continue
        if ti1:
            m = re_t.search(line)
            if m: t = m.group(1)
            m = re_dv.match(line)
            if m and t is not None: rows.append(f"s4 {t} {m.group(1)}"); t = None
    open(of + ".tmp", "w").write("\n".join(rows) + "\n"); os.rename(of + ".tmp", of)
    for f in (f"{runs}/{w}_sp4_energies.txt", f"{runs}/{w}_s4.check"):
        if os.path.exists(f): shutil.copy2(f, od)
    print(f"{sysn}/{leg}/{w}_s4 {len(rows)}", flush=True)
