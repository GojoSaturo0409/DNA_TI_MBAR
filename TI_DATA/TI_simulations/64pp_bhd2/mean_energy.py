"""
Mean MBAR Energy per Lambda – Log to Text File
===============================================
Reproduces the values shown in "Per-window Mean Sampled Energy vs λ" plot.
Writes lambda, mean energy, and std to mean_energy_per_lambda.txt.

Run from the folder containing your .out files:
    python mean_energy_per_lambda.py
"""

import os, re, glob, sys
import numpy as np
from datetime import datetime

# ═════════════════════════════════════════════════════════════════════════════
# SETTINGS
# ═════════════════════════════════════════════════════════════════════════════
LAMBDA_VALS = np.array([
    0.0092, 0.0479, 0.1150, 0.2063, 0.3161, 0.4374,
    0.5626, 0.6839, 0.7937, 0.8850, 0.9521, 0.9908
])
N_LAMBDA     = len(LAMBDA_VALS)
LAMBDA_TOL   = 1e-4
EQ_SKIP_FRAC = 0.20      # drop first 20% of frames (= first 2 ns)
OUT_FILE     = "mean_energy_per_lambda.txt"

SENTINEL = 999999999999999.0

# ═════════════════════════════════════════════════════════════════════════════
# REGEX
# ═════════════════════════════════════════════════════════════════════════════
_RE_RESULTS   = re.compile(r"4\.\s+RESULTS")
_RE_MBAR_HDR  = re.compile(r"MBAR Energy analysis:")
_RE_STARS     = re.compile(r"\*")
_RE_BASE_FILE = re.compile(r"^64pp_bhd2_w(\d+)\.out$")


def _match_lambda(lam):
    diffs = np.abs(LAMBDA_VALS - lam)
    idx   = int(np.argmin(diffs))
    return idx if diffs[idx] < LAMBDA_TOL else None


# ═════════════════════════════════════════════════════════════════════════════
# PARSE ONE FILE → list of MBAR energy blocks
# ═════════════════════════════════════════════════════════════════════════════
def parse_mbar_blocks(filepath):
    """
    Parse MBAR energy blocks from one AMBER .out file.

    *** values are replaced with SENTINEL (not skipped), so that all frames
    are retained and full N_LAMBDA blocks are preserved.  Frames containing
    sentinels are filtered out later when computing per-window statistics.
    """
    blocks        = []
    in_results    = False
    in_mbar       = False
    current_block = {}

    with open(filepath) as fh:
        for line in fh:
            if not in_results:
                if _RE_RESULTS.search(line):
                    in_results = True
                continue

            if _RE_MBAR_HDR.search(line):
                in_mbar, current_block = True, {}
                continue

            if in_mbar:
                if line.strip().startswith("---"):
                    in_mbar = False
                    if len(current_block) == N_LAMBDA:
                        blocks.append(
                            np.array([current_block[l] for l in LAMBDA_VALS])
                        )
                    current_block = {}
                    continue

                # Broader match that also captures starred value fields
                m_lam = re.search(r"Energy at\s+([\d.]+)\s*=\s*", line)
                if m_lam:
                    idx = _match_lambda(float(m_lam.group(1)))
                    if idx is not None:
                        val_field = line[m_lam.end():].strip()
                        if _RE_STARS.search(val_field):
                            # Replace *** with sentinel instead of skipping
                            current_block[LAMBDA_VALS[idx]] = SENTINEL
                        else:
                            try:
                                current_block[LAMBDA_VALS[idx]] = float(val_field.split()[0])
                            except (ValueError, IndexError):
                                current_block[LAMBDA_VALS[idx]] = SENTINEL

    return blocks


# ═════════════════════════════════════════════════════════════════════════════
# LOAD ALL WINDOWS
# ═════════════════════════════════════════════════════════════════════════════
def load_windows():
    files = sorted(
        f for f in glob.glob("64pp_bhd2_w*.out")
        if _RE_BASE_FILE.match(os.path.basename(f))
    )
    if not files:
        sys.exit("ERROR: No matching .out files found in the current directory.")

    windows = {}
    for f in files:
        wn     = int(_RE_BASE_FILE.match(os.path.basename(f)).group(1))
        blocks = parse_mbar_blocks(f)
        if not blocks:
            print(f"  [WARNING] No MBAR blocks in {f} – skipping.")
            continue
        n_raw  = len(blocks)
        skip   = int(np.ceil(n_raw * EQ_SKIP_FRAC))
        blocks = np.array(blocks[skip:])   # shape: (n_frames, N_LAMBDA)
        windows[wn] = {"mbar": blocks, "n_raw": n_raw}
        print(f"  w{wn:02d}: {n_raw} raw → {len(blocks)} frames kept")

    return windows


# ═════════════════════════════════════════════════════════════════════════════
# COMPUTE & WRITE
# ═════════════════════════════════════════════════════════════════════════════
def main():
    print("\nLoading windows …")
    windows = load_windows()

    if not windows:
        sys.exit("No data loaded.")

    # ── Build table: one row per lambda window ────────────────────────────────
    rows = []
    for i, wn in enumerate(sorted(windows.keys())):
        lam  = LAMBDA_VALS[i]
        diag = windows[wn]["mbar"][:, i]   # energy at window's own lambda

        # Filter out sentinel values before computing statistics
        valid = diag[diag < SENTINEL / 2]
        n_sentinel = len(diag) - len(valid)
        if n_sentinel:
            print(f"  [INFO] w{wn:02d}: {n_sentinel} sentinel frame(s) excluded "
                  f"from statistics ({len(valid)} valid frames remain).")

        rows.append({
            "window":     wn,
            "lambda":     lam,
            "n_frames":   len(valid),
            "n_sentinel": n_sentinel,
            "mean":       valid.mean() if valid.size else float("nan"),
            "std":        valid.std()  if valid.size else float("nan"),
            "sem":        valid.std() / np.sqrt(len(valid)) if valid.size else float("nan"),
        })

    # ── Print & write ─────────────────────────────────────────────────────────
    hdr = (f"{'Window':>6}  {'Lambda':>8}  {'N_valid':>8}  {'N_sentinel':>10}  "
           f"{'Mean (kcal/mol)':>16}  {'Std':>10}  {'SEM':>10}")
    sep = "-" * len(hdr)

    lines = [
        "=" * len(hdr),
        "  Per-window Mean Sampled MBAR Energy vs Lambda  (matches plot)",
        f"  Generated      : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"  Eq. skip       : first {EQ_SKIP_FRAC*100:.0f}% of frames (= {EQ_SKIP_FRAC*10:.1f} ns)",
        f"  Energy column  : diagonal – energy evaluated at the window's own lambda",
        f"  Note           : *** values replaced with sentinel; excluded from stats",
        "=" * len(hdr),
        "",
        hdr,
        sep,
    ]

    for r in rows:
        lines.append(
            f"  w{r['window']:02d}    {r['lambda']:>8.4f}  {r['n_frames']:>8d}  "
            f"{r['n_sentinel']:>10d}  {r['mean']:>+16.4f}  "
            f"{r['std']:>10.4f}  {r['sem']:>10.6f}"
        )

    lines += [sep, ""]

    output = "\n".join(lines)
    print("\n" + output)

    with open(OUT_FILE, "w") as fh:
        fh.write(output)

    print(f"Saved → {OUT_FILE}")


if __name__ == "__main__":
    main()