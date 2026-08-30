"""
MBAR ΔF Convergence vs Simulation Time  –  Standalone Script
=============================================================
Re-uses parsed window data to recompute MBAR ΔF at increasing fractions of
the productive trajectory and plots how the estimate evolves with time.

Usage
-----
  1. Place this script in the same directory as your 64pp_bhd2_w*.out files.
  2. Run:   python mbar_convergence_standalone.py

Key parameters (edit the USER SETTINGS block below)
----------------------------------------------------
  LAMBDA_VALS    – must match your simulation setup
  EQ_SKIP_FRAC   – fraction of each window discarded as equilibration
  SIM_LENGTH_NS  – total simulation length in ns
  CONV_TIME_NS   – absolute time points (ns) at which MBAR is evaluated
  N_WORKERS      – parallel worker processes (None = use all CPU cores)
"""

import os, re, glob, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed

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

EQ_SKIP_FRAC  = 0.20
SIM_LENGTH_NS = 10.0
CONV_TIME_NS  = [2, 4, 6, 8, 9, 10]  # ns at which MBAR is evaluated
N_WORKERS     = None                  # None → os.cpu_count()

OUT_DIR       = "."

COLOR_CONV    = "#e05c2a"
COLOR_EQ      = "lightgrey"

# ═════════════════════════════════════════════════════════════════════════════
# REGEX
# ═════════════════════════════════════════════════════════════════════════════
_RE_RESULTS   = re.compile(r"4\.\s+RESULTS")
_RE_MBAR_HDR  = re.compile(r"MBAR Energy analysis:")
_RE_MBAR_LINE = re.compile(r"Energy at\s+([\d.]+)\s*=\s*([-\d.]+)")
_RE_STARS     = re.compile(r"\*")
_RE_BASE_FILE = re.compile(r"^64pp_w(\d+)\.out$")


# ═════════════════════════════════════════════════════════════════════════════
# PARSERS
# ═════════════════════════════════════════════════════════════════════════════

def _match_lambda(lam: float) -> int | None:
    diffs = np.abs(LAMBDA_VALS - lam)
    idx   = int(np.argmin(diffs))
    return idx if diffs[idx] < LAMBDA_TOL else None


def parse_file(filepath: str) -> list[np.ndarray]:
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


def load_windows() -> dict[int, np.ndarray]:
    all_out = sorted(glob.glob("64pp_w*.out"))
    files   = [f for f in all_out if _RE_BASE_FILE.match(os.path.basename(f))]

    if not files:
        sys.exit("No window files found matching 64pp_w*.out")

    print(f"  Found {len(files)} window file(s).")
    windows: dict[int, np.ndarray] = {}
    for f in files:
        wn     = int(_RE_BASE_FILE.match(os.path.basename(f)).group(1))
        blocks = parse_file(f)
        if not blocks:
            print(f"  [WARNING] No MBAR blocks in {f!r} – skipping.")
            continue
        n_raw  = len(blocks)
        skip   = int(np.ceil(n_raw * EQ_SKIP_FRAC))
        arr    = np.array(blocks[skip:])
        windows[wn] = arr
        print(f"    w{wn:02d}: {n_raw} raw → {len(arr)} kept "
              f"(skipped {skip} eq. frames, {EQ_SKIP_FRAC*100:.0f}%)")
    return windows


# ═════════════════════════════════════════════════════════════════════════════
# MBAR HELPER
# ═════════════════════════════════════════════════════════════════════════════

def _build_u_nk(windows: dict[int, np.ndarray], frac: float) -> pd.DataFrame:
    sorted_wns = sorted(windows.keys())
    lambda_cols = LAMBDA_VALS.tolist()
    parts: list[pd.DataFrame] = []

    for wi, wn in enumerate(sorted_wns):
        lam_owner = LAMBDA_VALS[wi]
        mbar_full = windows[wn]
        n_use     = max(1, int(np.ceil(frac * len(mbar_full))))
        sliced    = mbar_full[:n_use]  # pre-slice the array once
        index = pd.MultiIndex.from_arrays(
            [np.arange(n_use, dtype=float),
             np.full(n_use, lam_owner)],
            names=["time", "lambdas"],
        )
        parts.append(pd.DataFrame(sliced, index=index, columns=lambda_cols))

    return pd.concat(parts)


def _mbar_worker(args: tuple) -> tuple[float, float, float] | None:
    """Top-level function so it can be pickled for multiprocessing."""
    windows, frac, abs_t_ns = args
    try:
        u_nk = _build_u_nk(windows, frac)
        mbar = MBAR()
        mbar.fit(u_nk)
        dF  = float(mbar.delta_f_.iloc[0, -1])
        ddF = float(mbar.d_delta_f_.iloc[0, -1])
        return abs_t_ns, dF, ddF
    except Exception as exc:
        print(f"    [WARN] MBAR failed at frac={frac:.2f}: {exc}")
        return None


# ═════════════════════════════════════════════════════════════════════════════
# CONVERGENCE SCAN + PLOT
# ═════════════════════════════════════════════════════════════════════════════

def compute_and_plot_convergence(
    windows: dict[int, np.ndarray],
) -> list[tuple[float, float, float]]:
    productive_ns = SIM_LENGTH_NS * (1.0 - EQ_SKIP_FRAC)
    eq_end_ns     = SIM_LENGTH_NS * EQ_SKIP_FRAC

    job_args = []
    for t_ns in CONV_TIME_NS:
        if t_ns <= eq_end_ns:
            print(f"  [SKIP] t={t_ns} ns is within equilibration – skipped.")
            continue
        frac = (t_ns - eq_end_ns) / productive_ns
        job_args.append((windows, frac, float(t_ns)))

    print(f"\n── Convergence scan  ({len(job_args)} points, parallel) ──")
    raw_results: dict[float, tuple[float, float, float]] = {}

    with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
        futures = {pool.submit(_mbar_worker, arg): arg[2] for arg in job_args}
        for fut in as_completed(futures):
            r = fut.result()
            if r is not None:
                abs_t_ns, dF, ddF = r
                raw_results[abs_t_ns] = (abs_t_ns, dF, ddF)
                print(f"  t = {abs_t_ns:.2f} ns  →  "
                      f"ΔF = {dF:+.4f} ± {ddF:.4f} kcal/mol")

    results = [raw_results[k] for k in sorted(raw_results)]

    if not results:
        print("[ERROR] No convergence data collected – check your window files.")
        return []

    times = np.array([r[0] for r in results])
    dFs   = np.array([r[1] for r in results])
    ddFs  = np.array([r[2] for r in results])

    final_dF  = dFs[-1]
    final_ddF = ddFs[-1]

    # ── Figure ────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 5))

    ax.axvspan(0, eq_end_ns,
               color=COLOR_EQ, alpha=0.50,
               label=f"Equilibration discarded (0–{eq_end_ns:.1f} ns)")

    ax.axhline(final_dF, color=COLOR_CONV, linestyle="--",
               linewidth=1.2, alpha=0.7)
    ax.axhspan(final_dF - final_ddF, final_dF + final_ddF,
               color=COLOR_CONV, alpha=0.10,
               label=(f"Final value: {final_dF:+.3f} ± {final_ddF:.3f} "
                      f"kcal/mol  (full {SIM_LENGTH_NS:.0f} ns)"))

    ax.errorbar(times, dFs, yerr=ddFs,
                fmt="o-", color=COLOR_CONV, capsize=5,
                linewidth=2.0, markersize=6,
                markeredgecolor="white", markeredgewidth=0.8,
                label="MBAR ΔF (running estimate)")

    ax.set_xlim(0, SIM_LENGTH_NS + 0.3)
    ax.set_xlabel("Simulation Time (ns)", fontsize=11)
    ax.set_ylabel("MBAR ΔF (kcal/mol)", fontsize=11)
    ax.legend(fontsize=9, loc="best")
    ax.grid(True, alpha=0.35)
    plt.tight_layout()

    fig_path = os.path.join(OUT_DIR, "mbar_convergence_vs_time.png")
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    print(f"Figure saved: {fig_path}")
    plt.close(fig)

    return results


# ═════════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════════

def main():
    matplotlib.rcParams.update({
        "font.size": 11, "axes.titlesize": 11, "axes.labelsize": 10,
    })
    print("\n── Loading windows ──")
    windows = load_windows()

    if len(windows) != N_LAMBDA:
        sys.exit(
            f"[ERROR] Expected {N_LAMBDA} windows, found {len(windows)}. "
            "Check file names and LAMBDA_VALS."
        )

    compute_and_plot_convergence(windows)


if __name__ == "__main__":
    main()