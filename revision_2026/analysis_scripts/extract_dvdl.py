#!/usr/bin/env python3
"""stdlib only (runs on ada login node). Same parsing as analyse_2stage.py dvdl_series():
per-step DV/DL of TI region 1 with TIME(PS), from wNN_s0/s1/s2.out.gz, up to 'A V E R A G E S'.
Writes <outdir>/<sys>/<leg>/wNN_dvdl.txt (seg time dvdl) only for windows with .done.
Usage: extract_dvdl.py <prod2c dir> <outdir> <sys>"""
import sys, os, re, gzip, glob
root, out, sysn = sys.argv[1:4]
re_dv = re.compile(r"^ DV/DL\s*=\s*(-?[\d.]+)"); re_t = re.compile(r"TIME\(PS\)\s*=\s*([\d.]+)")
r1 = re.compile(r"\|\s*TI region\s+1\b"); r2 = re.compile(r"\|\s*TI region\s+2\b")
for done in sorted(glob.glob(f"{root}/{sysn}/*/runs/w*.done")):
    runs = os.path.dirname(done); leg = runs.split("/")[-2]; w = os.path.basename(done)[:-5]
    od = f"{out}/{sysn}/{leg}"; os.makedirs(od, exist_ok=True); of = f"{od}/{w}_dvdl.txt"
    if os.path.exists(of): continue
    rows = []
    for seg in ("s0", "s1", "s2"):
        t = None; ti1 = False; n = 0
        for line in gzip.open(f"{runs}/{w}_{seg}.out.gz", "rt", errors="replace"):
            if "A V E R A G E S" in line: break
            if r1.search(line): ti1 = True; continue
            if r2.search(line): ti1 = False; continue
            if ti1:
                m = re_t.search(line)
                if m: t = m.group(1)
                m = re_dv.match(line)
                if m and t is not None: rows.append(f"{seg} {t} {m.group(1)}"); t = None; n += 1
        if n == 0: print(f"WARN {sysn}/{leg}/{w} {seg}: 0 DV/DL", flush=True)
    open(of + ".tmp", "w").write("\n".join(rows) + "\n"); os.rename(of + ".tmp", of)
    for f in glob.glob(f"{runs}/{w}_sp_energies.txt") + glob.glob(f"{runs}/{w}_s0.check"):
        os.system(f"cp -p {f} {od}/")
    print(f"{sysn}/{leg}/{w} {len(rows)}", flush=True)
