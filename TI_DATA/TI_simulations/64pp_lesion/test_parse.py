import re

_RE_RESULTS   = re.compile(r"4\.\s+RESULTS")
_RE_MBAR_HDR  = re.compile(r"MBAR Energy analysis:")
_RE_MBAR_LINE = re.compile(r"Energy at\s+([\d.]+)\s*=\s*([-\d.]+)")
_RE_DVDL_STEP = re.compile(r"DV/DL\s*=\s*([-\d.]+)")
_RE_STARS     = re.compile(r"\*")

TARGET_STEP = 2  # the Nth MBAR block / DV/DL value to inspect

filepath = "64pp_w1.out"  # change as needed

in_results  = False
in_mbar     = False
mbar_count  = 0
dvdl_count  = 0
current_block = {}

dvdl_target  = None
mbar_target  = None

with open(filepath) as fh:
    for line in fh:

        if not in_results:
            if _RE_RESULTS.search(line):
                in_results = True
            continue

        # DV/DL
        if not in_mbar:
            m = _RE_DVDL_STEP.search(line)
            if m and not _RE_STARS.search(m.group(1)):
                dvdl_count += 1
                if dvdl_count == TARGET_STEP:
                    dvdl_target = float(m.group(1))

        # MBAR block start
        if _RE_MBAR_HDR.search(line):
            in_mbar = True
            current_block = {}
            continue

        if in_mbar:
            if line.strip().startswith("---"):
                in_mbar = False
                if len(current_block) == 12:
                    mbar_count += 1
                    if mbar_count == TARGET_STEP:
                        mbar_target = dict(current_block)
                current_block = {}
                continue

            m = _RE_MBAR_LINE.search(line)
            if m and not _RE_STARS.search(m.group(2)):
                current_block[float(m.group(1))] = float(m.group(2))

# ── Report ──────────────────────────────────────────────────────────────────
print(f"File : {filepath}")
print(f"Target step : {TARGET_STEP}")
print()

print(f"DV/DL #{TARGET_STEP}:")
if dvdl_target is not None:
    print(f"  {dvdl_target}")
else:
    print(f"  [not found — only {dvdl_count} DV/DL values in file]")

print()
print(f"MBAR block #{TARGET_STEP}:")
if mbar_target is not None:
    for lam, val in sorted(mbar_target.items()):
        print(f"  Energy at {lam:.4f} = {val:.8f}")
else:
    print(f"  [not found — only {mbar_count} complete MBAR blocks in file]")
