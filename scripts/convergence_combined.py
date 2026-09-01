"""
MBAR ΔF Convergence vs Simulation Time  –  Multi-Folder Automated Script
=========================================================================
Runs across all 8 simulation folders, writes per-folder and combined results
to a text report, and plots all systems together in one figure with a
distinct colour per system.

Folders expected (relative to this script's location):
  64pp_bhd2            64pp_docked_bhd2
  64pp_bhd3            64pp_docked_bhd3
  64pp_lesion          64pp_docked_lesion
  64pp_partner         64pp_docked_partner

Usage
-----
  python mbar_convergence_multi.py

Key parameters – edit the USER SETTINGS block below.
"""

import os, re, glob, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime

try:
    from alchemlyb.estimators import MBAR
except ImportError:
    sys.exit("alchemlyb is required.  pip install alchemlyb")

# ═════════════════════════════════════════════════════════════════════════════
# USER SETTINGS
# ═════════════════════════════════════════════════════════════════════════════
LAMBDA_VALS = np.array([
    0.0092, 0.0479, 0.1150, 0.2063, 0.3161, 0.4374,
    0.5626, 0.6839, 0.7937, 0.8850, 0.9521, 0.9908
])
N_LAMBDA      = len(LAMBDA_VALS)
LAMBDA_TOL    = 1e-4

EQ_SKIP_FRAC  = 0.20          # fraction discarded as equilibration
SIM_LENGTH_NS = 10.0          # total simulation length per window (ns)
CONV_TIME_NS  = [2, 4, 6, 8, 9, 10]   # time points at which MBAR is evaluated
N_WORKERS     = None          # None → os.cpu_count()

# Script location is taken as the base directory that contains the 8 folders
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

FOLDERS = [
    "64pp_bhd2",
    "64pp_bhd3",
    "64pp_docked_bhd2",
    "64pp_docked_bhd3",
    "64pp_docked_lesion",
    "64pp_docked_partner",
    "64pp_lesion",
    "64pp_partner",
]

# Maps folder name → actual file prefix when it differs from the folder name
FILE_PREFIX_OVERRIDE = {
    "64pp_lesion":        "64pp",
    "64pp_docked_lesion": "64pp_docked",
}

# Output files (written next to this script)
OUT_DIR      = BASE_DIR
REPORT_FILE  = os.path.join(OUT_DIR, "mbar_convergence_report.txt")
FIGURE_FILE  = os.path.join(OUT_DIR, "mbar_convergence_combined.png")

# Colour palette – 8 visually distinct colours
PALETTE = [
    "#e05c2a",   # burnt orange    – bhd2
    "#2a7de0",   # cobalt blue     – bhd3
    "#27ae60",   # emerald green   – docked_bhd2
    "#8e44ad",   # amethyst purple – docked_bhd3
    "#f39c12",   # amber           – docked_lesion
    "#16a085",   # teal            – docked_partner
    "#c0392b",   # crimson         – lesion
    "#2980b9",   # steel blue      – partner
]

LINE_STYLES = ["-", "--", "-.", ":", "-", "--", "-.", ":"]
MARKERS     = ["o", "s", "^", "D", "v", "P", "X", "*"]

# ═════════════════════════════════════════════════════════════════════════════
# REGEX
# ═════════════════════════════════════════════════════════════════════════════
_RE_RESULTS   = re.compile(r"4\.\s+RESULTS")
_RE_MBAR_HDR  = re.compile(r"MBAR Energy analysis:")
_RE_MBAR_LINE = re.compile(r"Energy at\s+([\d.]+)\s*=\s*([-\d.]+)")
_RE_STARS     = re.compile(r"\*")
_RE_BASE_FILE = re.compile(r"^64pp_[a-z_]+_w(\d+)\.out$")


# ═════════════════════════════════════════════════════════════════════════════
# HELPERS
# ═════════════════════════════════════════════════════════════════════════════

def _match_lambda(lam: float) -> int | None:
    diffs = np.abs(LAMBDA_VALS - lam)
    idx   = int(np.argmin(diffs))
    return idx if diffs[idx] < LAMBDA_TOL else None


def parse_file(filepath: str) -> list[np.ndarray]:
    """Parse a single .out window file and return a list of MBAR energy rows."""
    mbar_blocks: list[np.ndarray] = []
    in_results    = False
    in_mbar       = False
    current_block: dict[float, float] = {}

    with open(filepath, "r") as fh:
        for line in fh:
            if not in_results:
                if _RE_RESULTS.search(line):
                    in_results = True
                continue

            if _RE_MBAR_HDR.search(line):
                in_mbar       = True
                current_block = {}
                continue

            if in_mbar:
                if line.strip().startswith("---"):
                    in_mbar = False
                    if len(current_block) == N_LAMBDA:
                        mbar_blocks.append(
                            np.array([current_block[l] for l in LAMBDA_VALS])
                        )
                    current_block = {}
                    continue

                m = _RE_MBAR_LINE.search(line)
                if m:
                    lam_s, val_s = m.group(1), m.group(2)
                    if _RE_STARS.search(val_s):
                        continue
                    try:
                        idx = _match_lambda(float(lam_s))
                        if idx is not None:
                            current_block[LAMBDA_VALS[idx]] = float(val_s)
                    except ValueError:
                        pass

    return mbar_blocks


def load_windows(folder_path: str, folder_name: str) -> dict[int, np.ndarray]:
    """
    Load and pre-process all window .out files from a folder.
    Returns dict  {window_index: 2-D array of shape (n_productive_frames, N_LAMBDA)}.
    """
    prefix   = FILE_PREFIX_OVERRIDE.get(folder_name, folder_name)
    pattern  = os.path.join(folder_path, f"{prefix}_w*.out")
    all_out  = sorted(glob.glob(pattern))
    base_re  = re.compile(rf"^{re.escape(prefix)}_w(\d+)\.out$")
    files    = [f for f in all_out if base_re.match(os.path.basename(f))]

    if not files:
        print(f"  [WARNING] No window files found in '{folder_path}' "
              f"matching pattern '{prefix}_w*.out'")
        return {}

    print(f"  Found {len(files)} window file(s) in '{folder_name}'.")
    windows: dict[int, np.ndarray] = {}
    for f in files:
        wn     = int(base_re.match(os.path.basename(f)).group(1))
        blocks = parse_file(f)
        if not blocks:
            print(f"    [WARNING] No MBAR blocks in {f!r} – skipping.")
            continue
        n_raw  = len(blocks)
        skip   = int(np.ceil(n_raw * EQ_SKIP_FRAC))
        arr    = np.array(blocks[skip:])
        windows[wn] = arr

    return windows


# ═════════════════════════════════════════════════════════════════════════════
# MBAR
# ═════════════════════════════════════════════════════════════════════════════

def _build_u_nk(windows: dict[int, np.ndarray], frac: float) -> pd.DataFrame:
    sorted_wns  = sorted(windows.keys())
    lambda_cols = LAMBDA_VALS.tolist()
    parts: list[pd.DataFrame] = []

    for wi, wn in enumerate(sorted_wns):
        lam_owner = LAMBDA_VALS[wi]
        mbar_full = windows[wn]
        n_use     = max(1, int(np.ceil(frac * len(mbar_full))))
        sliced    = mbar_full[:n_use]
        index = pd.MultiIndex.from_arrays(
            [np.arange(n_use, dtype=float),
             np.full(n_use, lam_owner)],
            names=["time", "lambdas"],
        )
        parts.append(pd.DataFrame(sliced, index=index, columns=lambda_cols))

    return pd.concat(parts)


def _mbar_worker(args: tuple) -> tuple[float, float, float] | None:
    windows, frac, abs_t_ns = args
    try:
        u_nk = _build_u_nk(windows, frac)
        mbar = MBAR()
        mbar.fit(u_nk)
        dF  = float(mbar.delta_f_.iloc[0, -1])
        ddF = float(mbar.d_delta_f_.iloc[0, -1])
        return abs_t_ns, dF, ddF
    except Exception as exc:
        print(f"    [WARN] MBAR failed at frac={frac:.4f}: {exc}")
        return None


def convergence_scan(
    windows: dict[int, np.ndarray],
) -> list[tuple[float, float, float]]:
    productive_ns = SIM_LENGTH_NS * (1.0 - EQ_SKIP_FRAC)
    eq_end_ns     = SIM_LENGTH_NS * EQ_SKIP_FRAC

    job_args = []
    for t_ns in CONV_TIME_NS:
        if t_ns <= eq_end_ns:
            continue
        frac = (t_ns - eq_end_ns) / productive_ns
        job_args.append((windows, frac, float(t_ns)))

    raw_results: dict[float, tuple[float, float, float]] = {}
    with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
        futures = {pool.submit(_mbar_worker, arg): arg[2] for arg in job_args}
        for fut in as_completed(futures):
            r = fut.result()
            if r is not None:
                raw_results[r[0]] = r

    return [raw_results[k] for k in sorted(raw_results)]


# ═════════════════════════════════════════════════════════════════════════════
# REPORT WRITER
# ═════════════════════════════════════════════════════════════════════════════

def write_report(
    all_results: dict[str, list[tuple[float, float, float]]],
    fh,
) -> None:
    eq_end_ns = SIM_LENGTH_NS * EQ_SKIP_FRAC

    fh.write("=" * 70 + "\n")
    fh.write("MBAR ΔF Convergence Report  –  All Systems\n")
    fh.write(f"Generated : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    fh.write(f"Base dir  : {BASE_DIR}\n")
    fh.write(f"Sim length: {SIM_LENGTH_NS:.1f} ns   "
             f"Eq skip: {EQ_SKIP_FRAC*100:.0f}%  "
             f"({eq_end_ns:.1f} ns discarded)\n")
    fh.write(f"Lambda values ({N_LAMBDA}): "
             + ", ".join(f"{l:.4f}" for l in LAMBDA_VALS) + "\n")
    fh.write("=" * 70 + "\n\n")

    for folder in FOLDERS:
        fh.write(f"{'─'*60}\n")
        fh.write(f"System : {folder}\n")
        fh.write(f"{'─'*60}\n")

        results = all_results.get(folder)
        if not results:
            fh.write("  [NO DATA – folder missing or no parsable .out files]\n\n")
            continue

        fh.write(f"  {'Time (ns)':>10}  {'ΔF (kcal/mol)':>15}  "
                 f"{'±σ (kcal/mol)':>15}\n")
        fh.write(f"  {'-'*10}  {'-'*15}  {'-'*15}\n")
        for (t, dF, ddF) in results:
            fh.write(f"  {t:10.2f}  {dF:+15.6f}  {ddF:15.6f}\n")

        final_t, final_dF, final_ddF = results[-1]
        fh.write(f"\n  Final estimate ({final_t:.1f} ns): "
                 f"ΔF = {final_dF:+.4f} ± {final_ddF:.4f} kcal/mol\n\n")

    # Summary table
    fh.write("=" * 70 + "\n")
    fh.write("SUMMARY  –  Final ΔF Values\n")
    fh.write("=" * 70 + "\n")
    fh.write(f"  {'System':<25}  {'ΔF (kcal/mol)':>15}  {'±σ':>10}  {'at (ns)':>8}\n")
    fh.write(f"  {'-'*25}  {'-'*15}  {'-'*10}  {'-'*8}\n")
    for folder in FOLDERS:
        results = all_results.get(folder)
        if results:
            t, dF, ddF = results[-1]
            fh.write(f"  {folder:<25}  {dF:+15.6f}  {ddF:10.6f}  {t:8.2f}\n")
        else:
            fh.write(f"  {folder:<25}  {'N/A':>15}  {'N/A':>10}  {'N/A':>8}\n")
    fh.write("\n")


# ═════════════════════════════════════════════════════════════════════════════
# COMBINED PLOT
# ═════════════════════════════════════════════════════════════════════════════

def plot_combined(
    all_results: dict[str, list[tuple[float, float, float]]],
) -> None:
    eq_end_ns = SIM_LENGTH_NS * EQ_SKIP_FRAC

    fig, ax = plt.subplots(figsize=(13, 6))

    # Equilibration shading
    ax.axvspan(0, eq_end_ns,
               color="lightgrey", alpha=0.45, zorder=0,
               label=f"Equilibration ({eq_end_ns:.1f} ns discarded)")

    for i, folder in enumerate(FOLDERS):
        results = all_results.get(folder)
        if not results:
            continue

        color  = PALETTE[i % len(PALETTE)]
        ls     = LINE_STYLES[i % len(LINE_STYLES)]
        marker = MARKERS[i % len(MARKERS)]

        times = np.array([r[0] for r in results])
        dFs   = np.array([r[1] for r in results])
        ddFs  = np.array([r[2] for r in results])

        # Shaded ±σ band
        ax.fill_between(times, dFs - ddFs, dFs + ddFs,
                        color=color, alpha=0.12, zorder=1)

        # Error-bar line
        ax.errorbar(
            times, dFs, yerr=ddFs,
            fmt=f"{marker}{ls}",
            color=color,
            capsize=4,
            linewidth=1.8,
            markersize=6,
            markeredgecolor="white",
            markeredgewidth=0.7,
            label=folder,
            zorder=2,
        )

    ax.set_xlim(0, SIM_LENGTH_NS + 0.3)
    ax.set_xlabel("Simulation Time (ns)", fontsize=12)
    ax.set_ylabel("MBAR  ΔF  (kcal/mol)", fontsize=12)
    ax.set_title("MBAR ΔF Convergence vs Simulation Time  –  All Systems",
                 fontsize=13, fontweight="bold")

    # Legend outside the plot area
    ax.legend(
        fontsize=8.5,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        borderaxespad=0,
        frameon=True,
        framealpha=0.9,
    )

    ax.grid(True, alpha=0.30, linestyle=":")
    plt.tight_layout(rect=[0, 0, 0.82, 1])   # leave room for legend

    fig.savefig(FIGURE_FILE, dpi=150, bbox_inches="tight")
    print(f"\nCombined figure saved: {FIGURE_FILE}")
    plt.close(fig)


# ═════════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════════

def main():
    matplotlib.rcParams.update({
        "font.size": 11,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
    })

    all_results: dict[str, list[tuple[float, float, float]]] = {}

    with open(REPORT_FILE, "w") as report_fh:
        for folder in FOLDERS:
            folder_path = os.path.join(BASE_DIR, folder)
            print(f"\n{'═'*60}")
            print(f"  Processing: {folder}")
            print(f"{'═'*60}")

            if not os.path.isdir(folder_path):
                print(f"  [SKIP] Directory not found: {folder_path}")
                all_results[folder] = []
                continue

            windows = load_windows(folder_path, folder)

            if not windows:
                print(f"  [SKIP] No usable window data in '{folder}'.")
                all_results[folder] = []
                continue

            if len(windows) != N_LAMBDA:
                print(f"  [WARNING] Expected {N_LAMBDA} windows, "
                      f"found {len(windows)} in '{folder}'. "
                      f"Proceeding anyway.")

            print(f"\n  ── Convergence scan ──")
            results = convergence_scan(windows)

            for (t, dF, ddF) in results:
                print(f"    t = {t:.2f} ns  →  "
                      f"ΔF = {dF:+.4f} ± {ddF:.4f} kcal/mol")

            all_results[folder] = results

        # Write consolidated text report
        write_report(all_results, report_fh)

    print(f"\nReport saved : {REPORT_FILE}")

    # Combined plot
    plot_combined(all_results)

    print("\nDone.")


if __name__ == "__main__":
    main()