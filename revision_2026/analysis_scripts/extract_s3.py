#!/usr/bin/env python3
"""stdlib only. DV/DL (TI region 1) of extension segments: <root>/<sys>/<leg>/runs/wNN_s3.out.gz (with wNN_s3.done)
-> <out>/<sys>/<leg>/wNN_s3_dvdl.txt (seg time dvdl) + copies of wNN_sp3_energies.txt, wNN_s3.check. Same parsing as extract_dvdl.py.
Usage: extract_s3.py <root e.g. ~/prod2c/ext3> <outdir> <sys>"""
import sys, os, re, gzip, glob, shutil
root, out, sysn = sys.argv[1:4]
re_dv = re.compile(r"^ DV/DL\s*=\s*(-?[\d.]+)"); re_t = re.compile(r"TIME\(PS\)\s*=\s*([\d.]+)")
r1 = re.compile(r"\|\s*TI region\s+1\b"); r2 = re.compile(r"\|\s*TI region\s+2\b")
for done in sorted(glob.glob(f"{root}/{sysn}/*/runs/w*_s3.done")):
    runs = os.path.dirname(done); leg = runs.split("/")[-2]; w = os.path.basename(done)[:3]
    od = f"{out}/{sysn}/{leg}"; os.makedirs(od, exist_ok=True); of = f"{od}/{w}_s3_dvdl.txt"
    if os.path.exists(of): continue
    rows = []; t = None; ti1 = False
    for line in gzip.open(f"{runs}/{w}_s3.out.gz", "rt", errors="replace"):
        if "A V E R A G E S" in line: break
        if r1.search(line): ti1 = True; continue
        if r2.search(line): ti1 = False; continue
        if ti1:
            m = re_t.search(line)
            if m: t = m.group(1)
            m = re_dv.match(line)
            if m and t is not None: rows.append(f"s3 {t} {m.group(1)}"); t = None
    open(of + ".tmp", "w").write("\n".join(rows) + "\n"); os.rename(of + ".tmp", of)
    for f in (f"{runs}/{w}_sp3_energies.txt", f"{runs}/{w}_s3.check"):
        if os.path.exists(f): shutil.copy2(f, od)
    print(f"{sysn}/{leg}/{w}_s3 {len(rows)}", flush=True)
