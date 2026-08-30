"""
AMBER Free Energy Analysis: MBAR + TI  (10 ns only)
=====================================================
Parses AMBER .out files for MBAR energy blocks and per-step DV/DL values,
computes binding free energy via alchemlyb MBAR, and produces convergence
plots for the 10 ns simulation only.

File conventions
----------------
  10 ns run  :  64pp_w1.out  ... 64pp_w12.out     (one per lambda window)

Key design choices
------------------
  * Strict regex anchor: 64pp_w1.out never matches 64pp_w1_1.out.
  * Equilibration skip: first EQ_SKIP_FRAC fraction of frames discarded
    (default 0.20 = first 2 ns of 10 ns).
  * dV/dλ read from per-step "DV/DL = " lines inside .out files.
    Lines containing '*' are skipped.
  * Final binding affinity + convergence diagnostics saved to .txt file.
"""

import os, re, glob, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")          # remove this line if you want interactive windows
import matplotlib.pyplot as plt
from datetime import datetime

# ── alchemlyb ────────────────────────────────────────────────────────────────
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
N_LAMBDA    = len(LAMBDA_VALS)
LAMBDA_TOL  = 1e-4

# First EQ_SKIP_FRAC fraction of frames is treated as equilibration
# and discarded.  0.20 → first 2 ns of a 10 ns run.
EQ_SKIP_FRAC = 0.20

OUT_DIR = "."   # where to write figures and the summary text file

COLOR = "steelblue"

# ═════════════════════════════════════════════════════════════════════════════
# REGEX
# ═════════════════════════════════════════════════════════════════════════════
_RE_RESULTS   = re.compile(r"4\.\s+RESULTS")
_RE_MBAR_HDR  = re.compile(r"MBAR Energy analysis:")
_RE_MBAR_LINE = re.compile(r"Energy at\s+([\d.]+)\s*=\s*([-\d.]+)")
_RE_DVDL_STEP = re.compile(r"DV/DL\s*=\s*([-\d.]+)")
_RE_STARS     = re.compile(r"\*")
_RE_BASE_FILE = re.compile(r"^64pp_partner_w(\d+)\.out$")   # e.g. 64pp_w7.out only

SENTINEL = 999999999999999.0
# ═════════════════════════════════════════════════════════════════════════════
# PARSERS
# ═════════════════════════════════════════════════════════════════════════════

def _match_lambda(lam: float) -> int | None:
    diffs = np.abs(LAMBDA_VALS - lam)
    idx   = int(np.argmin(diffs))
    return idx if diffs[idx] < LAMBDA_TOL else None


def parse_file(filepath: str) -> tuple[list[np.ndarray], np.ndarray]:
    """
    Parse one AMBER .out file.  Returns:
        mbar_blocks : list of 1-D float arrays (len N_LAMBDA), one per MD step
        dvdl_series : 1-D float array of per-step DV/DL values

    *** values in MBAR blocks are replaced with SENTINEL (not skipped),
    so that all frames are retained and full N_LAMBDA blocks are preserved.
    DV/DL lines with *** are still skipped.
    """
    mbar_blocks: list[np.ndarray] = []
    dvdl_vals:   list[float]      = []

    in_results    = False
    in_mbar       = False
    current_block: dict[float, float] = {}

    with open(filepath, "r") as fh:
        for line in fh:
            if not in_results:
                if _RE_RESULTS.search(line):
                    in_results = True
                continue

            # Per-step DV/DL (only when outside an MBAR block)
            # *** in DV/DL still skipped — unphysical gradient is unusable
            if not in_mbar:
                m = _RE_DVDL_STEP.search(line)
                if m and not _RE_STARS.search(m.group(1)):
                    try:
                        dvdl_vals.append(float(m.group(1)))
                    except ValueError:
                        pass

            # Start of MBAR block
            if _RE_MBAR_HDR.search(line):
                in_mbar       = True
                current_block = {}
                continue

            # Inside MBAR block
            if in_mbar:
                if line.strip().startswith("---"):
                    in_mbar = False
                    # Accept block if all N_LAMBDA entries present
                    # (including sentinels — missing entries still excluded)
                    if len(current_block) == N_LAMBDA:
                        mbar_blocks.append(
                            np.array([current_block[l] for l in LAMBDA_VALS])
                        )
                    current_block = {}
                    continue

                # Match "Energy at X = Y" or "Energy at X = ****"
                # Use a broader regex that also captures starred values
                m_lam = re.search(r"Energy at\s+([\d.]+)\s*=\s*", line)
                if m_lam:
                    lam_s = m_lam.group(1)
                    idx   = _match_lambda(float(lam_s))
                    if idx is not None:
                        # Check if the value field contains stars
                        val_field = line[m_lam.end():].strip()
                        if _RE_STARS.search(val_field):
                            # Replace *** with sentinel
                            current_block[LAMBDA_VALS[idx]] = SENTINEL
                        else:
                            try:
                                current_block[LAMBDA_VALS[idx]] = float(val_field.split()[0])
                            except (ValueError, IndexError):
                                current_block[LAMBDA_VALS[idx]] = SENTINEL

    return mbar_blocks, np.asarray(dvdl_vals)

# ═════════════════════════════════════════════════════════════════════════════
# WINDOW LOADING
# ═════════════════════════════════════════════════════════════════════════════

def load_windows(eq_skip_frac: float = EQ_SKIP_FRAC) -> dict[int, dict]:
    """
    Load 10 ns windows: files whose basename matches 64pp_w<N>.out exactly.
    The _1.out continuation files are excluded by the strict regex.

    First eq_skip_frac fraction of frames is dropped (equilibration).
    Returns dict {window_number: {"mbar": ndarray, "dvdl": ndarray, ...}}
    """
    all_out = sorted(glob.glob("64pp_partner_w*.out"))
    files   = [f for f in all_out if _RE_BASE_FILE.match(os.path.basename(f))]

    print(f"  Found {len(files)} window files.")
    windows: dict[int, dict] = {}
    for f in files:
        wn = int(_RE_BASE_FILE.match(os.path.basename(f)).group(1))
        blocks, dvdl = parse_file(f)
        if not blocks:
            print(f"  [WARNING] No MBAR blocks in {f!r} – skipping.")
            continue
        n_raw  = len(blocks)
        skip   = int(np.ceil(n_raw * eq_skip_frac))
        blocks = blocks[skip:]
        dvdl   = dvdl[skip:] if dvdl.size > skip else dvdl[0:0]
        windows[wn] = {
            "mbar":  np.array(blocks),
            "dvdl":  dvdl,
            "n_raw": n_raw,
            "file":  os.path.basename(f),
        }
        print(f"    w{wn:02d}: {n_raw} raw frames → {len(blocks)} kept "
              f"(dropped first {skip}, {eq_skip_frac*100:.0f}% eq.)")
    return windows


# ═════════════════════════════════════════════════════════════════════════════
# MBAR
# ═════════════════════════════════════════════════════════════════════════════

def build_u_nk(windows: dict[int, dict]) -> pd.DataFrame:
    sorted_wns = sorted(windows.keys())
    if len(sorted_wns) != N_LAMBDA:
        raise ValueError(
            f"Expected {N_LAMBDA} windows, got {len(sorted_wns)}: {sorted_wns}"
        )
    rows_idx:  list[tuple]      = []
    rows_data: list[np.ndarray] = []
    for wi, wn in enumerate(sorted_wns):
        lam_owner = LAMBDA_VALS[wi]
        for fi, block in enumerate(windows[wn]["mbar"]):
            rows_idx.append((float(fi), lam_owner))
            rows_data.append(block)
    index = pd.MultiIndex.from_tuples(rows_idx, names=["time", "lambdas"])
    return pd.DataFrame(rows_data, index=index, columns=LAMBDA_VALS.tolist())


def run_mbar(windows: dict[int, dict]) -> tuple[float, float] | None:
    if len(windows) != N_LAMBDA:
        print(f"  [SKIP] Need {N_LAMBDA} windows, have {len(windows)}.")
        return None
    try:
        u_nk = build_u_nk(windows)
        mbar = MBAR()
        mbar.fit(u_nk)
        dF  = float(mbar.delta_f_.iloc[0, -1])
        ddF = float(mbar.d_delta_f_.iloc[0, -1])
        return dF, ddF
    except Exception as exc:
        print(f"  [ERROR] MBAR failed: {exc}")
        return None


# ═════════════════════════════════════════════════════════════════════════════
# PLOT HELPERS
# ═════════════════════════════════════════════════════════════════════════════

def _savefig(fig: plt.Figure, name: str) -> None:
    path = os.path.join(OUT_DIR, name)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"  Saved: {path}")
    plt.close(fig)


def _rolling(arr: np.ndarray, frac: float = 0.05) -> np.ndarray:
    w = max(1, int(len(arr) * frac))
    return pd.Series(arr).rolling(w, min_periods=1).mean().values


def _panel_grid() -> tuple[plt.Figure, list[plt.Axes]]:
    fig, axes = plt.subplots(3, 4, figsize=(18, 12), sharex=False)
    return fig, axes.flatten().tolist()


# ─────────────────────────────────────────────────────────────────────────────
# Plot 1 – dV/dλ time series, one panel per window
# ─────────────────────────────────────────────────────────────────────────────
def plot_dvdl_per_window(windows: dict[int, dict]) -> None:
    fig, axes = _panel_grid()
    for i, wn in enumerate(sorted(windows.keys())):
        ax  = axes[i]
        lam = LAMBDA_VALS[i]
        d   = windows[wn]["dvdl"]
        if d.size:
            ax.plot(d, color=COLOR, alpha=0.4, linewidth=0.6)
            ax.plot(_rolling(d), color=COLOR, linewidth=2, label="10 ns")
        ax.set_title(f"λ = {lam:.4f}", fontsize=9)
        ax.set_xlabel("Frame", fontsize=8)
        ax.set_ylabel("dV/dλ (kcal/mol)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3)
    fig.suptitle("dV/dλ per Lambda Window  (raw + rolling mean) – 10 ns",
                 fontsize=12, y=1.01)
    plt.tight_layout()
    _savefig(fig, "dvdl_per_window.png")


# ─────────────────────────────────────────────────────────────────────────────
# Plot 2 – TI integrand  ⟨dV/dλ⟩ ± σ  vs λ
# ─────────────────────────────────────────────────────────────────────────────
def plot_dvdl_mean_vs_lambda(windows: dict[int, dict]) -> None:
    means, stds = [], []
    for wn in sorted(windows.keys()):
        d = windows[wn]["dvdl"]
        means.append(d.mean() if d.size else np.nan)
        stds.append(d.std()   if d.size else np.nan)
    means = np.array(means)
    stds  = np.array(stds)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.errorbar(LAMBDA_VALS, means, yerr=stds, fmt="o-", capsize=5,
                color=COLOR, linewidth=2, label="10 ns")
    ax.fill_between(LAMBDA_VALS, means - stds, means + stds,
                    alpha=0.15, color=COLOR)
    ax.set_xlabel("λ")
    ax.set_ylabel("⟨dV/dλ⟩ (kcal/mol)")
    ax.set_title("TI Integrand ⟨dV/dλ⟩ vs λ – 10 ns")
    ax.legend()
    ax.grid(True, alpha=0.4)
    plt.tight_layout()
    _savefig(fig, "dvdl_mean_vs_lambda.png")


# ─────────────────────────────────────────────────────────────────────────────
# Plot 3 – Cumulative mean dV/dλ vs frame (convergence check)
# ─────────────────────────────────────────────────────────────────────────────
def plot_cumulative_dvdl(windows: dict[int, dict]) -> None:
    fig, axes = _panel_grid()
    for i, wn in enumerate(sorted(windows.keys())):
        ax  = axes[i]
        lam = LAMBDA_VALS[i]
        d   = windows[wn]["dvdl"]
        if d.size:
            cum = np.cumsum(d) / np.arange(1, len(d) + 1)
            ax.plot(cum, color=COLOR, linewidth=1.5)
            ax.axhline(cum[-1], color="grey", linestyle="--",
                       linewidth=1, alpha=0.7, label=f"final = {cum[-1]:.2f}")
            ax.legend(fontsize=7)
        ax.set_title(f"λ = {lam:.4f}", fontsize=9)
        ax.set_xlabel("Frame", fontsize=8)
        ax.set_ylabel("Cumulative ⟨dV/dλ⟩", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3)
    fig.suptitle("Cumulative Mean dV/dλ per Window  (flat = converged) – 10 ns",
                 fontsize=12, y=1.01)
    plt.tight_layout()
    _savefig(fig, "cumulative_dvdl_per_window.png")


# ─────────────────────────────────────────────────────────────────────────────
# Plot 4 – MBAR diagonal energy vs frame, one panel per window
# ─────────────────────────────────────────────────────────────────────────────
def plot_mbar_energy_per_window(windows: dict[int, dict]) -> None:
    fig, axes = _panel_grid()
    for i, wn in enumerate(sorted(windows.keys())):
        ax   = axes[i]
        lam  = LAMBDA_VALS[i]
        mbar = windows[wn]["mbar"]
        if mbar.shape[0]:
            diag = mbar[:, i]
            ax.plot(diag, color=COLOR, alpha=0.35, linewidth=0.6)
            ax.plot(_rolling(diag), color=COLOR, linewidth=2)
        ax.set_title(f"λ = {lam:.4f}", fontsize=9)
        ax.set_xlabel("Frame", fontsize=8)
        ax.set_ylabel("Energy (kcal/mol)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3)
    fig.suptitle("MBAR Diagonal Energy per Window  (raw + rolling mean) – 10 ns",
                 fontsize=12, y=1.01)
    plt.tight_layout()
    _savefig(fig, "mbar_energy_per_window.png")


# ─────────────────────────────────────────────────────────────────────────────
# Plot 5 – Mean MBAR diagonal energy ± σ vs λ
# ─────────────────────────────────────────────────────────────────────────────
def plot_mean_energy_vs_lambda(windows: dict[int, dict]) -> None:
    means, stds = [], []
    for i, wn in enumerate(sorted(windows.keys())):
        d = windows[wn]["mbar"][:, i]
        means.append(d.mean())
        stds.append(d.std())
    means = np.array(means)
    stds  = np.array(stds)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.errorbar(LAMBDA_VALS, means, yerr=stds, fmt="o-", capsize=5,
                color=COLOR, linewidth=2, label="10 ns")
    ax.fill_between(LAMBDA_VALS, means - stds, means + stds,
                    alpha=0.15, color=COLOR)
    ax.set_xlabel("λ")
    ax.set_ylabel("Mean energy (kcal/mol)")
    ax.set_title("Per-window Mean Sampled Energy vs λ – 10 ns")
    ax.legend()
    ax.grid(True, alpha=0.4)
    plt.tight_layout()
    _savefig(fig, "mean_energy_vs_lambda.png")


# ─────────────────────────────────────────────────────────────────────────────
# Plot 6 – MBAR ΔF result (single bar with error)
# ─────────────────────────────────────────────────────────────────────────────
def plot_mbar_result(result: tuple[float, float] | None) -> None:
    if result is None:
        return
    dF, ddF = result
    fig, ax = plt.subplots(figsize=(4, 5))
    bar = ax.bar(["10 ns"], [dF], yerr=[ddF], capsize=12,
                 color=COLOR, edgecolor="black", linewidth=0.8, width=0.4)
    ax.text(0, dF + ddF + abs(dF) * 0.015,
            f"{dF:+.3f} ± {ddF:.3f}", ha="center", va="bottom", fontsize=11)
    ax.set_ylabel("ΔF (kcal/mol)")
    ax.set_title("MBAR Binding Free Energy – 10 ns")
    ax.grid(axis="y", alpha=0.4)
    plt.tight_layout()
    _savefig(fig, "mbar_binding_energy.png")


# ─────────────────────────────────────────────────────────────────────────────
# Plot 7 – Block analysis: running block mean of dV/dλ (all windows)
# ─────────────────────────────────────────────────────────────────────────────
def plot_block_analysis(windows: dict[int, dict]) -> None:
    """
    For each window: cumulative mean dV/dλ as a function of fraction of
    simulation used.  A flat curve means the window is converged.
    All 12 windows shown in one figure.
    """
    fig, axes = _panel_grid()
    fracs = np.linspace(0.05, 1.0, 50)

    for i, wn in enumerate(sorted(windows.keys())):
        ax  = axes[i]
        lam = LAMBDA_VALS[i]
        d   = windows[wn]["dvdl"]
        if not d.size:
            ax.set_title(f"λ = {lam:.4f}  (no data)", fontsize=9)
            continue

        ns_axis    = fracs * 10.0        # 10 ns total
        block_mean = []
        block_sem  = []
        for f in fracs:
            n = max(1, int(f * len(d)))
            s = d[:n]
            block_mean.append(s.mean())
            block_sem.append(s.std() / np.sqrt(len(s)))
        block_mean = np.array(block_mean)
        block_sem  = np.array(block_sem)

        ax.plot(ns_axis, block_mean, color=COLOR, linewidth=2)
        ax.fill_between(ns_axis,
                        block_mean - block_sem,
                        block_mean + block_sem,
                        alpha=0.2, color=COLOR)
        ax.axhline(block_mean[-1], color="grey", linestyle="--",
                   linewidth=1, alpha=0.8)
        ax.set_title(f"λ = {lam:.4f}", fontsize=9)
        ax.set_xlabel("Simulation used (ns)", fontsize=8)
        ax.set_ylabel("Mean dV/dλ (kcal/mol)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3)

    fig.suptitle("Block Analysis – Running Mean dV/dλ vs Simulation Length  "
                 "(flat tail = converged) – 10 ns",
                 fontsize=11, y=1.01)
    plt.tight_layout()
    _savefig(fig, "block_analysis_dvdl.png")


# ═════════════════════════════════════════════════════════════════════════════
# RESULTS TEXT FILE
# ═════════════════════════════════════════════════════════════════════════════

def write_summary(
    result:   tuple[float, float] | None,
    windows:  dict[int, dict],
    eq_skip:  float,
) -> None:
    path  = os.path.join(OUT_DIR, "binding_affinity_results.txt")
    lines = [
        "=" * 65,
        "  AMBER MBAR FREE ENERGY ANALYSIS  –  RESULTS SUMMARY",
        f"  Generated : {datetime.now().strftime('%Y-%m-%d  %H:%M:%S')}",
        "=" * 65,
        "",
        "SIMULATION SETTINGS",
        f"  Simulation length  : 10 ns",
        f"  Lambda windows     : {N_LAMBDA}",
        f"  Lambda values      : {list(np.round(LAMBDA_VALS, 4))}",
        f"  Equilibration skip : first {eq_skip*100:.0f}% of frames "
        f"(= {eq_skip*10:.1f} ns)",
        "",
        "FRAME COUNTS  (after equilibration removal)",
    ]
    for wn in sorted(windows.keys()):
        n  = windows[wn]["mbar"].shape[0]
        nr = windows[wn]["n_raw"]
        lam = LAMBDA_VALS[wn - 1] if wn <= N_LAMBDA else 0.0
        lines.append(f"  w{wn:02d}  λ={lam:.4f}:  {nr} raw → {n} frames used")

    lines += ["", "MBAR FREE ENERGY RESULT"]
    if result is not None:
        dF, ddF = result
        lines += [
            f"  ΔF (10 ns) = {dF:+.4f} ± {ddF:.4f} kcal/mol",
            "",
            "CONVERGENCE INDICATORS",
            "  (Inspect the plots for definitive assessment)",
        ]
        # Per-window dV/dλ drift: last-quarter mean vs overall mean
        lines.append("")
        lines.append("  Per-window dV/dλ stability (last 25% vs full mean):")
        all_stable = True
        for i, wn in enumerate(sorted(windows.keys())):
            d   = windows[wn]["dvdl"]
            lam = LAMBDA_VALS[i]
            if d.size < 4:
                lines.append(f"    w{wn:02d} λ={lam:.4f}: insufficient data")
                continue
            full_mean  = d.mean()
            tail_mean  = d[int(0.75 * len(d)):].mean()
            drift_pct  = abs(tail_mean - full_mean) / (abs(full_mean) + 1e-12) * 100
            stable     = drift_pct < 5.0   # < 5% drift = stable
            all_stable = all_stable and stable
            flag = "✓" if stable else "⚠ drift"
            lines.append(
                f"    w{wn:02d} λ={lam:.4f}:  full={full_mean:+.3f},  "
                f"tail={tail_mean:+.3f},  drift={drift_pct:.1f}%  {flag}"
            )

        verdict = "CONVERGED ✓" if all_stable else \
                  "SOME WINDOWS SHOW DRIFT – check plots carefully"
        lines += [
            "",
            f"  Overall verdict : {verdict}",
            "",
            "─" * 65,
            "  FINAL BINDING FREE ENERGY  (10 ns, after 2 ns eq. skip)",
            f"  ΔF = {dF:+.4f} ± {ddF:.4f} kcal/mol",
            "─" * 65,
        ]
    else:
        lines.append("  [No result – MBAR did not converge or windows missing]")

    lines.append("")
    text = "\n".join(lines)
    with open(path, "w") as fh:
        fh.write(text)
    print(f"\n{'='*65}")
    print(text)
    print(f"  Results also written to: {path}")


# ═════════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════════

def main():
    matplotlib.rcParams.update({
        "font.size": 11, "axes.titlesize": 11, "axes.labelsize": 10,
    })

    print("\n── Loading 10 ns windows ──")
    windows = load_windows(EQ_SKIP_FRAC)

    if len(windows) != N_LAMBDA:
        print(f"[WARNING] Only {len(windows)}/{N_LAMBDA} windows loaded.")

    # ── MBAR ─────────────────────────────────────────────────────────────────
    print("\n── MBAR Analysis ──")
    result = run_mbar(windows)
    if result:
        print(f"  ΔF (10 ns) = {result[0]:+.4f} ± {result[1]:.4f} kcal/mol")

    # ── Plots ─────────────────────────────────────────────────────────────────
    print("\n── Generating plots ──")
    plot_dvdl_per_window(windows)           # 1
    plot_dvdl_mean_vs_lambda(windows)       # 2
    plot_cumulative_dvdl(windows)           # 3
    plot_mbar_energy_per_window(windows)    # 4
    plot_mean_energy_vs_lambda(windows)     # 5
    plot_mbar_result(result)               # 6
    plot_block_analysis(windows)           # 7

    # ── Summary ───────────────────────────────────────────────────────────────
    write_summary(result, windows, EQ_SKIP_FRAC)


if __name__ == "__main__":
    main()