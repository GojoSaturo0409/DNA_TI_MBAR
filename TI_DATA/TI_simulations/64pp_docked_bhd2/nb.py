"""
extract_ti_energies.py  (v3 — corrected GQ-TI  +  NB-only MBAR)
================================================================
Extracts energy components from TI region 1 of each AMBER .out file.
Accepts ANY *.out file whose name ends in _w<N>.out (any prefix).

Two composite quantities are computed at every MD step:

  Case 1  – full non-bonded + 1-4 terms:
            1-4 NB  +  1-4 EEL  +  VDWAALS  +  EELEC

  Case 2  – long-range non-bonded only:
            VDWAALS  +  EELEC

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
TI FREE ENERGY — ENERGY DECOMPOSITION ANALYSIS (EDA)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
The goal is NOT to compute the total alchemical ΔF, but to decompose
the free energy contribution from non-bonded terms across the λ pathway.

For a linear Hamiltonian (AMBER default, klambda=1):
    V(x; λ) = (1-λ)·V₀(x)  +  λ·V₁(x)

The exact TI expression is:
    ΔF_NB = ∫₀¹ ⟨dV_NB/dλ⟩_λ dλ  =  ∫₀¹ ⟨V_NB,1(x) - V_NB,0(x)⟩_λ dλ

At each GL node λᵢ, the integrand ⟨dV_NB/dλ(λᵢ)⟩ is estimated by computing
dV_NB/dλ per frame via finite differences between adjacent window energies:

    dV_NB/dλ(xₙ) ≈ [V_NB(xₙ; λ_{i+1}) - V_NB(xₙ; λ_{i-1})]
                    / (λ_{i+1} - λ_{i-1})       (central, interior)

    dV_NB/dλ(xₙ) ≈ [V_NB(xₙ; λ₁) - V_NB(xₙ; λ₀)] / (λ₁ - λ₀)
                                                  (one-sided, boundary)

These per-frame derivatives are then averaged over the ensemble at each window
to give ⟨dV_NB/dλ(λᵢ)⟩, the TI integrand.  Gaussian quadrature gives:

    ΔF_NB = Σᵢ  wᵢ · ⟨dV_NB/dλ(λᵢ)⟩

A trapezoidal integration of ⟨dV_NB/dλ⟩ vs λ is also computed for comparison.

The per-window standard deviation of dV_NB/dλ(xₙ) is reported as the
uncertainty on the integrand (fluctuations of the derivative within each window).

NOTE: The raw ⟨V_NB(λ)⟩ energies are still plotted in ti_energies_mean_vs_lambda.png
for reference, but they are NOT the TI integrand.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
MBAR — NB-ONLY USING AMBER SOFTCORE FRAMEWORK
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
AMBER's softcore TI outputs per-frame NB energies in TI region 1.
For a linear Hamiltonian the full potential is:

    V(x; λ) = V_common(x) + (1-λ)·V_NB,0(x) + λ·V_NB,1(x)

where V_common contains all bonded and common-atom terms that are
IDENTICAL across all λ windows.  Because V_common is the same constant
for every frame at every λ, it cancels exactly in all MBAR free energy
differences:

    ΔF(k→k') = -kT ln [ Σₙ exp(-β·u_{k',n}) / Σₙ exp(-β·u_{k,n}) ]

    u_{k,n} = β·V(xₙ; λ_k)
            = β·[V_common(xₙ) + (1-λ_k)·V_NB,0(xₙ) + λ_k·V_NB,1(xₙ)]

The V_common term shifts all u_{k,n} by the same constant per frame and
cancels in differences → we only need the NB part.

For the softcore case (ifsc=1), AMBER's output in TI region 1 gives:
    V_NB(x; λ) = VDWAALS + EELEC + 1-4 NB + 1-4 EEL
where VDWAALS already contains the softcore-modified vdW energy between
SC and common atoms (Eq. 25.5/25.6 of AMBER manual).  The SC_VDW_DER
and SC_EEL_DER terms (the λ-derivative corrections from SC potentials)
are part of DV/DL, not of the energy itself, so they do not appear here.

The reduced potential used for MBAR is therefore:
    u_kn[k, n] = β · V_NB(xₙ; λ_k)

For configurations xₙ sampled at window k_src, we need V_NB evaluated
at all K target λ values.  Since we only have V_NB at the sampled λ,
we use the linear-Hamiltonian identity:

    V_NB(x; λ') = V_NB(x; λ) + (λ'-λ) · [V_NB,1(x) - V_NB,0(x)]
                = V_NB(x; λ) + (λ'-λ) · dV_NB/dλ(x)

where dV_NB/dλ(x) per frame is estimated by central finite differences
across neighbouring windows (same frame index n):

    dV_NB/dλ(xₙ) ≈ [V_NB(xₙ; λ_{k+1}) - V_NB(xₙ; λ_{k-1})]
                    / (λ_{k+1} - λ_{k-1})

This is the same approach AMBER uses internally when ifmbar=1 is set
(section 25.1.9 of AMBER manual) — evaluate energies at neighbouring
λ values to extrapolate the reduced potential.

The numerical mean-shift stabilisation (E_src - E_src.mean()) is applied
before multiplying by β; additive per-source-window constants cancel
exactly in MBAR's log-sum-exp differences.

Outputs
-------
  ti_energies_w<N>.csv           – per-window CSV
  ti_energies_all.csv            – concatenated table
  ti_energies_timeseries.png     – Case 1 & Case 2 energy timeseries
  ti_energies_mean_vs_lambda.png – ⟨Case 1⟩ and ⟨Case 2⟩ ± std vs λ
  nb_integrand.png               – TI integrand ⟨V_NB(λ)⟩ vs λ + ΔF bar
  nb_slope.png                   – d⟨V_NB⟩/dλ diagnostic plot
  nb_summary.csv                 – λ, ⟨V_NB⟩, std, GQ contrib, ΔF
  mbar_results.csv               – MBAR ΔF and uncertainty per window

Usage
-----
    python nb.py [data_dir]

Dependencies
------------
    numpy, pandas, matplotlib
    pymbar  (optional — MBAR section skipped if not installed)
        pip install pymbar
"""

import os, re, sys, glob, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ─── Compatibility shim: np.trapezoid was added in NumPy 2.0 ─────────────────
_trapz = getattr(np, "trapezoid", np.trapz)

# ─── 12-point Gauss-Legendre quadrature nodes (λ) and weights ────────────────
# Source: AMBER manual Table 25.1, 12-point row.
# Nodes must match the clambda values used in your AMBER runs exactly.
LAMBDA_VALS = np.array([
    0.00922, 0.04794, 0.11505, 0.20634, 0.31608, 0.43738,
    0.56262, 0.68392, 0.79366, 0.88495, 0.95206, 0.99078,
])

# Corresponding Gauss-Legendre weights (sum to 1.0 over [0,1]).
GAUSS_WEIGHTS = np.array([
    0.02359,   # λ = 0.00922
    0.05347,   # λ = 0.04794
    0.08004,   # λ = 0.11505
    0.10158,   # λ = 0.20634
    0.11675,   # λ = 0.31608
    0.12457,   # λ = 0.43738
    0.12457,   # λ = 0.56262
    0.11675,   # λ = 0.68392
    0.10158,   # λ = 0.79366
    0.08004,   # λ = 0.88495
    0.05347,   # λ = 0.95206
    0.02359,   # λ = 0.99078
])

assert len(LAMBDA_VALS) == len(GAUSS_WEIGHTS), "Mismatch between λ and weights!"
assert abs(GAUSS_WEIGHTS.sum() - 1.0) < 1e-6, \
    f"Gauss weights sum to {GAUSS_WEIGHTS.sum():.6f}, expected 1.0"

N_LAMBDA = len(LAMBDA_VALS)
KBT      = 0.5961  # kcal/mol at 298 K  (k_B × 298)

# ─── Regex patterns ───────────────────────────────────────────────────────────
_RE_TI1_START  = re.compile(r"\|\s*TI region\s+1\b")
_RE_TI2_START  = re.compile(r"\|\s*TI region\s+2\b")
_RE_NSTEP      = re.compile(r"NSTEP\s*=\s*(\d+)")
_RE_NB14_LINE  = re.compile(
    r"1-4 NB\s*=\s*([-\d.]+)\s+1-4 EEL\s*=\s*([-\d.]+)\s+VDWAALS\s*=\s*([-\d.]+)"
)
_RE_EELEC      = re.compile(r"EELEC\s*=\s*([-\d.]+)")
_RE_CLAMBDA    = re.compile(r"clambda\s*=\s*([\d.]+)")
_RE_STARS      = re.compile(r"\*")
_RE_FNAME      = re.compile(r"^.+_w(\d+)\.out$")

OUT_DIR = os.path.dirname(os.path.abspath(__file__))


# ═════════════════════════════════════════════════════════════════════════════
# PARSER
# ═════════════════════════════════════════════════════════════════════════════

def parse_file(filepath: str) -> tuple:
    steps      = []
    nb14_vals  = []
    eel14_vals = []
    vdw_vals   = []
    eelec_vals = []
    clambda      = np.nan
    in_ti1       = False
    current_step = None
    nb14_val = eel14_val = vdw_val = eelec_val = None

    def _commit():
        if current_step is not None and None not in (nb14_val, eel14_val, vdw_val, eelec_val):
            steps.append(current_step)
            nb14_vals.append(nb14_val)
            eel14_vals.append(eel14_val)
            vdw_vals.append(vdw_val)
            eelec_vals.append(eelec_val)

    with open(filepath, "r") as fh:
        for line in fh:
            if np.isnan(clambda):
                m = _RE_CLAMBDA.search(line)
                if m:
                    try:
                        clambda = float(m.group(1))
                    except ValueError:
                        pass

            if _RE_TI1_START.search(line):
                in_ti1 = True
                current_step = None
                nb14_val = eel14_val = vdw_val = eelec_val = None
                continue

            if in_ti1 and _RE_TI2_START.search(line):
                _commit()
                in_ti1 = False
                current_step = None
                nb14_val = eel14_val = vdw_val = eelec_val = None
                continue

            if not in_ti1:
                continue

            m = _RE_NSTEP.search(line)
            if m:
                _commit()
                current_step = int(m.group(1))
                nb14_val = eel14_val = vdw_val = eelec_val = None
                continue

            m = _RE_NB14_LINE.search(line)
            if m and not _RE_STARS.search(line):
                try:
                    nb14_val  = float(m.group(1))
                    eel14_val = float(m.group(2))
                    vdw_val   = float(m.group(3))
                except ValueError:
                    pass
                continue

            m = _RE_EELEC.search(line)
            if m and not _RE_STARS.search(m.group(1)):
                try:
                    eelec_val = float(m.group(1))
                except ValueError:
                    pass
                continue

    _commit()
    nb14_arr  = np.array(nb14_vals,  dtype=float)
    eel14_arr = np.array(eel14_vals, dtype=float)
    vdw_arr   = np.array(vdw_vals,   dtype=float)
    eelec_arr = np.array(eelec_vals, dtype=float)

    df = pd.DataFrame({
        "step":    steps,
        "nb14":    nb14_arr,
        "eel14":   eel14_arr,
        "vdwaals": vdw_arr,
        "eelec":   eelec_arr,
        "case1":   nb14_arr + eel14_arr + vdw_arr + eelec_arr,
        "case2":   vdw_arr  + eelec_arr,
    })
    return df, clambda


# ═════════════════════════════════════════════════════════════════════════════
# LOAD ALL WINDOWS
# ═════════════════════════════════════════════════════════════════════════════

def load_all_windows(data_dir: str) -> dict:
    all_out = sorted(glob.glob(os.path.join(data_dir, "*.out")))
    raw: list = []
    seen_wn: set = set()

    for fpath in all_out:
        bname = os.path.basename(fpath)
        m = _RE_FNAME.match(bname)
        if not m:
            print(f"  [SKIP] {bname}  (name does not match *_w<N>.out)")
            continue
        wn = int(m.group(1))
        if wn in seen_wn:
            print(f"  [SKIP] {bname}  (window {wn} already loaded)")
            continue
        seen_wn.add(wn)

        print(f"  Parsing {bname} ...", end=" ", flush=True)
        df, clambda_file = parse_file(fpath)

        if 1 <= wn <= len(LAMBDA_VALS):
            clambda = LAMBDA_VALS[wn - 1]
            gauss_w = GAUSS_WEIGHTS[wn - 1]
            if not np.isnan(clambda_file) and abs(clambda - clambda_file) > 1e-4:
                print(f"\n  [WARNING] w{wn}: file clambda={clambda_file:.5f}, "
                      f"using canonical λ={clambda:.5f}")
        else:
            clambda = clambda_file
            gauss_w = np.nan
            print(f"\n  [WARNING] w{wn} outside 1–{len(LAMBDA_VALS)}; "
                  f"Gaussian weight set to NaN")

        print(f"{len(df)} steps  (λ={clambda:.5f}, w_GL={gauss_w:.5f})")
        raw.append({"wn": wn, "df": df, "clambda": clambda,
                    "gauss_w": gauss_w, "file": bname})

    if not raw:
        return {}

    raw.sort(key=lambda d: d["wn"])
    windows = {}
    for rank, entry in enumerate(raw, start=1):
        windows[rank] = {
            "df":      entry["df"],
            "clambda": entry["clambda"],
            "gauss_w": entry["gauss_w"],
            "wn":      entry["wn"],
            "file":    entry["file"],
        }
    return windows


# ═════════════════════════════════════════════════════════════════════════════
# SAVE CSVs
# ═════════════════════════════════════════════════════════════════════════════

def save_csvs(windows: dict, out_dir: str) -> None:
    all_rows = []
    for rank in sorted(windows.keys()):
        entry = windows[rank]
        df    = entry["df"].copy()
        wn    = entry["wn"]
        lam   = entry["clambda"]
        csv_path = os.path.join(out_dir, f"ti_energies_w{wn:02d}.csv")
        df.to_csv(csv_path, index=False)
        print(f"  Saved: {csv_path}")
        df["window"] = wn
        df["lambda"] = lam
        all_rows.append(df)

    if all_rows:
        combined = pd.concat(all_rows, ignore_index=True)
        cols = ["window", "lambda", "step",
                "nb14", "eel14", "vdwaals", "eelec", "case1", "case2"]
        combined = combined[[c for c in cols if c in combined.columns]]
        all_path = os.path.join(out_dir, "ti_energies_all.csv")
        combined.to_csv(all_path, index=False)
        print(f"  Saved combined: {all_path}")


# ═════════════════════════════════════════════════════════════════════════════
# PLOTS — ENERGY TIMESERIES
# ═════════════════════════════════════════════════════════════════════════════

COLOR_C1  = "steelblue"
COLOR_C2  = "darkorange"
ALPHA_RAW = 0.25
LW_RAW    = 0.5
LW_ROLL   = 1.8
ROLL_FRAC = 0.05


def _rolling(arr: np.ndarray, frac: float = ROLL_FRAC) -> np.ndarray:
    w = max(1, int(len(arr) * frac))
    return pd.Series(arr).rolling(w, min_periods=1).mean().values


def _savefig(fig: plt.Figure, name: str) -> None:
    path = os.path.join(OUT_DIR, name)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"  Saved: {path}")
    plt.close(fig)


def plot_timeseries(windows: dict) -> None:
    n_win  = len(windows)
    n_cols = min(4, n_win)
    n_rows = (n_win + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(5.5 * n_cols, 4.5 * n_rows))
    axes_flat = np.array(axes).flatten()

    legend_handles = [
        plt.Line2D([0], [0], color=COLOR_C1, lw=2,
                   label="Case 1: 1-4 NB + 1-4 EEL + VDWAALS + EELEC"),
        plt.Line2D([0], [0], color=COLOR_C2, lw=2,
                   label="Case 2: VDWAALS + EELEC"),
    ]

    for i, rank in enumerate(sorted(windows.keys())):
        entry = windows[rank]
        df    = entry["df"]
        lam   = entry["clambda"]
        wn    = entry["wn"]
        x     = np.arange(len(df))
        ax    = axes_flat[i]

        c1 = df["case1"].values
        ax.plot(x, c1, color=COLOR_C1, alpha=ALPHA_RAW, linewidth=LW_RAW)
        ax.plot(x, _rolling(c1), color=COLOR_C1, linewidth=LW_ROLL)
        ax.set_ylabel("Case 1 (kcal/mol)", color=COLOR_C1, fontsize=7)
        ax.tick_params(axis="y", labelcolor=COLOR_C1, labelsize=6)

        c2  = df["case2"].values
        ax2 = ax.twinx()
        ax2.plot(x, c2, color=COLOR_C2, alpha=ALPHA_RAW, linewidth=LW_RAW)
        ax2.plot(x, _rolling(c2), color=COLOR_C2, linewidth=LW_ROLL)
        ax2.set_ylabel("Case 2 (kcal/mol)", color=COLOR_C2, fontsize=7)
        ax2.tick_params(axis="y", labelcolor=COLOR_C2, labelsize=6)

        ax.set_title(f"w{wn}  λ={lam:.5f}", fontsize=9)
        ax.set_xlabel("Frame", fontsize=7)
        ax.tick_params(axis="x", labelsize=6)
        ax.grid(True, alpha=0.22)

    for j in range(n_win, len(axes_flat)):
        axes_flat[j].set_visible(False)

    fig.legend(handles=legend_handles, loc="upper center",
               ncol=2, fontsize=10, bbox_to_anchor=(0.5, 1.03))
    fig.suptitle(
        "TI Region 1 – Energy Component Time Series per Lambda Window\n"
        "(transparent = raw,  solid = rolling mean)",
        fontsize=12, y=1.06,
    )
    plt.tight_layout()
    _savefig(fig, "ti_energies_timeseries.png")


def plot_mean_vs_lambda(windows: dict) -> None:
    lams  = []
    c1_mean, c1_std = [], []
    c2_mean, c2_std = [], []

    for rank in sorted(windows.keys()):
        entry = windows[rank]
        df    = entry["df"]
        lams.append(entry["clambda"])
        c1_mean.append(df["case1"].mean());  c1_std.append(df["case1"].std())
        c2_mean.append(df["case2"].mean());  c2_std.append(df["case2"].std())

    lams = np.array(lams)
    c1_mean = np.array(c1_mean);  c1_std = np.array(c1_std)
    c2_mean = np.array(c2_mean);  c2_std = np.array(c2_std)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 9), sharex=True)
    ax1.errorbar(lams, c1_mean, yerr=c1_std, fmt="o-", capsize=5,
                 color=COLOR_C1, linewidth=2)
    ax1.fill_between(lams, c1_mean - c1_std, c1_mean + c1_std,
                     alpha=0.15, color=COLOR_C1)
    ax1.set_ylabel("⟨Case 1⟩ (kcal/mol)", fontsize=11)
    ax1.set_title("Case 1:  ⟨1-4 NB + 1-4 EEL + VDWAALS + EELEC⟩  per λ window",
                  fontsize=11)
    ax1.grid(True, alpha=0.3)

    ax2.errorbar(lams, c2_mean, yerr=c2_std, fmt="s-", capsize=5,
                 color=COLOR_C2, linewidth=2)
    ax2.fill_between(lams, c2_mean - c2_std, c2_mean + c2_std,
                     alpha=0.15, color=COLOR_C2)
    ax2.set_xlabel("λ", fontsize=11)
    ax2.set_ylabel("⟨Case 2⟩ (kcal/mol)", fontsize=11)
    ax2.set_title("Case 2:  ⟨VDWAALS + EELEC⟩  per λ window", fontsize=11)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    _savefig(fig, "ti_energies_mean_vs_lambda.png")


# ═════════════════════════════════════════════════════════════════════════════
# TI — NB ENERGY DECOMPOSITION
# ═════════════════════════════════════════════════════════════════════════════

def compute_ti(windows: dict) -> dict:
    """
    Compute NB free energy contributions via Gaussian quadrature and
    trapezoidal integration.

    TI integrand
    ------------
    The correct TI integrand is ⟨dV_NB/dλ⟩ — the ensemble average of the
    alchemical derivative at each window.

        ΔF_NB = Σᵢ  wᵢ · ⟨dV_NB/dλ(λᵢ)⟩      (Gaussian quadrature)
        ΔF_NB = ∫ ⟨dV_NB/dλ(λ)⟩ dλ            (trapezoidal, for comparison)

    dV_NB/dλ is estimated per frame using finite differences of V_NB between
    adjacent windows (same frame index n).  Interior windows use central
    differences; boundary windows use one-sided differences:

        central:   dV/dλ(xₙ) ≈ [V(xₙ; λ_{i+1}) - V(xₙ; λ_{i-1})]
                                / (λ_{i+1} - λ_{i-1})
        boundary:  dV/dλ(xₙ) ≈ [V(xₙ; λ₁)  - V(xₙ; λ₀)]  / (λ₁  - λ₀)
                   dV/dλ(xₙ) ≈ [V(xₙ; λ_{K-1}) - V(xₙ; λ_{K-2})] / Δλ

    The per-window mean of dV/dλ is the integrand; its std is the uncertainty.

    The raw ⟨V_NB⟩ and σ(V_NB) are retained for plotting/diagnostics.
    """
    ranks    = sorted(windows.keys())
    lambdas  = np.array([windows[r]["clambda"] for r in ranks])
    gweights = np.array([windows[r]["gauss_w"] for r in ranks])
    K        = len(ranks)

    # ── Per-window energy arrays (truncate to shortest window) ────────────────
    series_c1 = [windows[r]["df"]["case1"].values for r in ranks]
    series_c2 = [windows[r]["df"]["case2"].values for r in ranks]
    min_frames = min(len(s) for s in series_c1)

    E_c1 = np.stack([s[:min_frames] for s in series_c1], axis=0)  # (K, N)
    E_c2 = np.stack([s[:min_frames] for s in series_c2], axis=0)  # (K, N)

    # ── Per-window mean and std of V_NB (for reference / diagnostics) ─────────
    c1_means = E_c1.mean(axis=1)
    c2_means = E_c2.mean(axis=1)
    c1_stds  = E_c1.std(axis=1)
    c2_stds  = E_c2.std(axis=1)

    # ── Per-frame dV_NB/dλ via finite differences of adjacent windows ─────────
    # dvdl_c1[i, n] = dV_NB/dλ for frame n at window i
    dvdl_c1 = np.empty_like(E_c1)
    dvdl_c2 = np.empty_like(E_c2)
    for i in range(K):
        if i == 0:
            dl = lambdas[1] - lambdas[0]
            dvdl_c1[i] = (E_c1[1] - E_c1[0]) / dl
            dvdl_c2[i] = (E_c2[1] - E_c2[0]) / dl
        elif i == K - 1:
            dl = lambdas[-1] - lambdas[-2]
            dvdl_c1[i] = (E_c1[-1] - E_c1[-2]) / dl
            dvdl_c2[i] = (E_c2[-1] - E_c2[-2]) / dl
        else:
            dl = lambdas[i+1] - lambdas[i-1]
            dvdl_c1[i] = (E_c1[i+1] - E_c1[i-1]) / dl
            dvdl_c2[i] = (E_c2[i+1] - E_c2[i-1]) / dl

    # ── TI integrand: ⟨dV/dλ⟩ and its std at each window ─────────────────────
    integrand_c1      = dvdl_c1.mean(axis=1)   # shape (K,)
    integrand_c2      = dvdl_c2.mean(axis=1)
    integrand_c1_stds = dvdl_c1.std(axis=1)
    integrand_c2_stds = dvdl_c2.std(axis=1)

    # ── Gaussian quadrature: ΔF_NB = Σᵢ wᵢ · ⟨dV_NB/dλ(λᵢ)⟩ ─────────────────
    dF_c1_gauss = float(np.dot(gweights, integrand_c1))
    dF_c2_gauss = float(np.dot(gweights, integrand_c2))

    # ── Trapezoidal integration of ⟨dV/dλ⟩ vs λ (for comparison) ─────────────
    dF_c1_trap = float(_trapz(integrand_c1, lambdas))
    dF_c2_trap = float(_trapz(integrand_c2, lambdas))

    # ── Cumulative ΔF across windows (for plotting convergence) ───────────────
    cum_c1_gauss = np.array([np.dot(gweights[:i+1], integrand_c1[:i+1])
                              for i in range(K)])
    cum_c2_gauss = np.array([np.dot(gweights[:i+1], integrand_c2[:i+1])
                              for i in range(K)])
    cum_c1_trap  = np.array([_trapz(integrand_c1[:i+1], lambdas[:i+1])
                              if i > 0 else 0.0
                              for i in range(K)])
    cum_c2_trap  = np.array([_trapz(integrand_c2[:i+1], lambdas[:i+1])
                              if i > 0 else 0.0
                              for i in range(K)])

    # ── slope_c1/c2: d⟨V_NB⟩/dλ from window means (kept for diagnostic plot) ─
    slope_c1 = integrand_c1   # same quantity — mean of per-frame dV/dλ
    slope_c2 = integrand_c2

    return dict(
        lambdas           = lambdas,
        gauss_weights     = gweights,
        # raw energy stats (diagnostics)
        c1_means          = c1_means,
        c2_means          = c2_means,
        c1_stds           = c1_stds,
        c2_stds           = c2_stds,
        # TI integrand: ⟨dV/dλ⟩ per window
        integrand_c1      = integrand_c1,
        integrand_c2      = integrand_c2,
        integrand_c1_stds = integrand_c1_stds,
        integrand_c2_stds = integrand_c2_stds,
        # integrated free energies
        dF_c1_gauss  = dF_c1_gauss,
        dF_c2_gauss  = dF_c2_gauss,
        dF_c1_trap   = dF_c1_trap,
        dF_c2_trap   = dF_c2_trap,
        cum_c1_gauss = cum_c1_gauss,
        cum_c2_gauss = cum_c2_gauss,
        cum_c1_trap  = cum_c1_trap,
        cum_c2_trap  = cum_c2_trap,
        slope_c1     = slope_c1,
        slope_c2     = slope_c2,
    )


# ═════════════════════════════════════════════════════════════════════════════
# MBAR — NB ONLY, SOFTCORE-CONSISTENT
# ═════════════════════════════════════════════════════════════════════════════

def compute_mbar(windows: dict, case: int) -> dict | None:
    """
    Run MBAR using NB-only energies from TI region 1.

    Physical basis
    --------------
    AMBER's alchemical Hamiltonian:

        V(x; λ) = V_common(x) + (1-λ)·V_NB,0(x) + λ·V_NB,1(x)

    V_common (bonded terms, common-atom interactions) is identical across
    all λ windows.  Because it cancels exactly in all MBAR free energy
    differences, we need only the NB part:

        V_NB(x; λ) = (1-λ)·V_NB,0(x) + λ·V_NB,1(x)

    For a configuration xₙ sampled at window k_src, the NB energy at any
    other target state k_tgt is obtained via the linear-Hamiltonian identity:

        V_NB(xₙ; λ_{k_tgt}) = V_NB(xₙ; λ_{k_src})
                              + (λ_{k_tgt} - λ_{k_src}) · dV_NB/dλ(xₙ)

    where dV_NB/dλ(xₙ) = V_NB,1(xₙ) - V_NB,0(xₙ) is the per-frame
    alchemical derivative.  We estimate this from sampled NB energies
    using central finite differences across neighbouring windows:

        dV_NB/dλ(xₙ) ≈ [V_NB(xₙ; λ_{k+1}) - V_NB(xₙ; λ_{k-1})]
                        / (λ_{k+1} - λ_{k-1})

    This is exactly what AMBER does when ifmbar=1 is set (AMBER manual
    section 25.1.9) — energies are evaluated at neighbouring λ values
    to construct the u_kn matrix for post-processing.

    For softcore potentials (ifsc=1), AMBER's VDWAALS in TI region 1
    already contains the softcore-modified vdW energy between SC and
    common atoms (Eqs. 25.5/25.6).  SC_VDW_DER / SC_EEL_DER are
    λ-derivative corrections that appear in DV/DL only — they are NOT
    part of the energy printed in TI region 1, so no special handling
    is needed here.

    The reduced potential used for MBAR:
        u_kn[k_tgt, n] = β · V_NB(xₙ; λ_{k_tgt})

    Numerical stability: subtract per-source-window mean before scaling
    by β.  These additive constants cancel in all MBAR log-sum-exp
    differences (only u differences matter in MBAR).

    Parameters
    ----------
    windows : dict   — loaded window data from load_all_windows()
    case    : int    — 1 (full NB + 1-4) or 2 (VDWAALS + EELEC only)

    Returns
    -------
    dict with lambdas, dF_mbar (kcal/mol, relative to λ=0 window), dF_std
    """
    try:
        from pymbar import MBAR
    except ImportError:
        print("  [SKIP] pymbar not installed.  pip install pymbar")
        return None

    col     = f"case{case}"
    ranks   = sorted(windows.keys())
    K       = len(ranks)
    lambdas = np.array([windows[r]["clambda"] for r in ranks])
    beta    = 1.0 / KBT

    # ── Collect NB energy timeseries, truncate to shortest window ─────────────
    # All windows must have the same frame count so we can index frame n
    # across windows k-1, k, k+1 for the finite-difference dV_NB/dλ.
    series_list  = [windows[r]["df"][col].values for r in ranks]
    min_frames   = min(len(s) for s in series_list)
    if min_frames < 2:
        print(f"  [SKIP] MBAR Case {case}: too few frames ({min_frames})")
        return None

    # E_matrix[k, n] = V_NB(xₙ; λ_k)  sampled at window k
    E_matrix = np.stack([s[:min_frames] for s in series_list], axis=0)  # (K, N)

    # ── Per-frame dV_NB/dλ via central finite differences ─────────────────────
    # dvdl[k, n] ≈ dV_NB/dλ for frame n at window k.
    # For the linear Hamiltonian this equals V_NB,1(xₙ) - V_NB,0(xₙ),
    # estimated here from neighbouring sampled windows.
    dvdl = np.empty_like(E_matrix)   # (K, N)
    for i in range(K):
        if i == 0:
            dl       = lambdas[1] - lambdas[0]
            dvdl[i]  = (E_matrix[1] - E_matrix[0]) / dl
        elif i == K - 1:
            dl       = lambdas[-1] - lambdas[-2]
            dvdl[i]  = (E_matrix[-1] - E_matrix[-2]) / dl
        else:
            dl       = lambdas[i+1] - lambdas[i-1]
            dvdl[i]  = (E_matrix[i+1] - E_matrix[i-1]) / dl

    # ── Build u_kn matrix  (K × K·N) ─────────────────────────────────────────
    # Samples from all K windows are concatenated column-wise.
    # For samples xₙ drawn at k_src, evaluated at k_tgt:
    #
    #   u_kn[k_tgt, col_n] = β · [E_NB(xₙ; λ_{k_src})
    #                              + (λ_{k_tgt} - λ_{k_src}) · dV_NB/dλ(xₙ)]
    #
    # Mean-shift per source window for numerical stability (cancels in MBAR).
    N_k     = np.full(K, min_frames, dtype=int)
    N_total = int(N_k.sum())
    u_kn    = np.zeros((K, N_total))

    n_start = 0
    for k_src in range(K):
        n_end    = n_start + min_frames
        E_src    = E_matrix[k_src]               # V_NB(x; λ_{k_src}), shape (N,)
        dv_src   = dvdl[k_src]                   # dV_NB/dλ(x) at k_src, shape (N,)
        E_shift  = E_src - E_src.mean()          # subtract mean for stability

        for k_tgt in range(K):
            dlam = lambdas[k_tgt] - lambdas[k_src]
            # V_NB(x; λ_{k_tgt}) ≈ E_shift + Δλ · dV_NB/dλ(x)  [+ constant]
            u_kn[k_tgt, n_start:n_end] = beta * (E_shift + dlam * dv_src)

        n_start = n_end

    # ── Run MBAR ──────────────────────────────────────────────────────────────
    print(f"  Running MBAR for Case {case}  (K={K}, N={N_total}, "
          f"frames/window={min_frames}) ...")
    try:
        import warnings as _warnings
        with _warnings.catch_warnings(record=True) as caught:
            _warnings.simplefilter("always")
            mbar = MBAR(u_kn, N_k, verbose=False)
        for w_ in caught:
            msg = str(w_.message).lower()
            if "converge" in msg or "tolerance" in msg:
                print(f"  [MBAR WARNING] {w_.message}")

        results = mbar.compute_free_energy_differences()
        Deltaf  = results["Delta_f"]
        dDeltaf = results["dDelta_f"]

        # dF_mbar[i] = ΔF(λ₀ → λᵢ) in kcal/mol
        dF_mbar = Deltaf[0,  :] / beta
        dF_std  = dDeltaf[0, :] / beta

        return dict(
            lambdas = lambdas,
            dF_mbar = dF_mbar,
            dF_std  = dF_std,
            case    = case,
        )

    except Exception as e:
        print(f"  [ERROR] MBAR failed for Case {case}: {e}")
        return None


# ═════════════════════════════════════════════════════════════════════════════
# PLOTS — TI INTEGRAND AND MBAR
# ═════════════════════════════════════════════════════════════════════════════

def plot_nb_integrand(ti: dict, mbar_c1: dict | None,
                      mbar_c2: dict | None) -> None:
    """
    Three-panel figure:
      Panel 1: TI integrand ⟨V_NB(λ)⟩ vs λ  (error bars = per-window std)
      Panel 2: Cumulative ΔF_NB vs λ
      Panel 3: Method comparison bar chart (GQ, Trapz, MBAR)
    """
    lams              = ti["lambdas"]
    c1_means          = ti["integrand_c1"]
    c2_means          = ti["integrand_c2"]
    c1_stds           = ti["integrand_c1_stds"]
    c2_stds           = ti["integrand_c2_stds"]
    gweights          = ti["gauss_weights"]
    dF_c1_gauss       = ti["dF_c1_gauss"]
    dF_c2_gauss       = ti["dF_c2_gauss"]
    dF_c1_trap        = ti["dF_c1_trap"]
    dF_c2_trap        = ti["dF_c2_trap"]
    cum_c1_gauss      = ti["cum_c1_gauss"]
    cum_c2_gauss      = ti["cum_c2_gauss"]
    cum_c1_trap       = ti["cum_c1_trap"]
    cum_c2_trap       = ti["cum_c2_trap"]

    has_mbar = (mbar_c1 is not None) or (mbar_c2 is not None)
    n_panels = 3 if has_mbar else 2
    fig, axes = plt.subplots(n_panels, 1,
                             figsize=(11, 5 * n_panels),
                             sharex=(n_panels < 3))

    # ── Panel 1: TI integrand ⟨V_NB(λ)⟩ vs λ ─────────────────────────────────
    ax1 = axes[0]
    ax1.errorbar(lams, c1_means, yerr=c1_stds, fmt="o-", capsize=5,
                 color=COLOR_C1, linewidth=2,
                 label=(f"Case 1  GQ ΔF={dF_c1_gauss:+.3f},"
                        f"  Trap ΔF={dF_c1_trap:+.3f} kcal/mol"))
    ax1.fill_between(lams, c1_means - c1_stds, c1_means + c1_stds,
                     alpha=0.15, color=COLOR_C1)
    ax1.errorbar(lams, c2_means, yerr=c2_stds, fmt="s--", capsize=5,
                 color=COLOR_C2, linewidth=2,
                 label=(f"Case 2  GQ ΔF={dF_c2_gauss:+.3f},"
                        f"  Trap ΔF={dF_c2_trap:+.3f} kcal/mol"))
    ax1.fill_between(lams, c2_means - c2_stds, c2_means + c2_stds,
                     alpha=0.15, color=COLOR_C2)
    ax1.axhline(0, color="black", linewidth=0.7, linestyle=":")
    ax1.set_ylabel("⟨dV_NB/dλ(λ)⟩  (kcal/mol)", fontsize=11)
    ax1.set_title(
        "TI Integrand: ⟨dV_NB/dλ(λᵢ)⟩ vs λ\n"
        "ΔF_NB = Σᵢ wᵢ·⟨dV_NB/dλ(λᵢ)⟩   (error bars = per-window σ of dV/dλ)",
        fontsize=11)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)
    if n_panels == 2:
        ax1.set_xlabel("λ", fontsize=11)

    # ── Panel 2: Cumulative ΔF ────────────────────────────────────────────────
    ax2 = axes[1]
    ax2.plot(lams, cum_c1_gauss, "o-",  color=COLOR_C1, linewidth=2,
             label="Case 1  Gaussian QD")
    ax2.plot(lams, cum_c1_trap,  "o--", color=COLOR_C1, linewidth=1.5,
             alpha=0.6, label="Case 1  Trapz")
    ax2.plot(lams, cum_c2_gauss, "s-",  color=COLOR_C2, linewidth=2,
             label="Case 2  Gaussian QD")
    ax2.plot(lams, cum_c2_trap,  "s--", color=COLOR_C2, linewidth=1.5,
             alpha=0.6, label="Case 2  Trapz")
    ax2.set_ylabel("Cumulative ΔF_NB  (kcal/mol)", fontsize=11)
    ax2.set_title("Cumulative NB Free Energy Contribution: GQ vs Trapezoidal",
                  fontsize=11)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)
    if not has_mbar:
        ax2.set_xlabel("λ", fontsize=11)

    # ── Panel 3: Method comparison bar chart ──────────────────────────────────
    if has_mbar:
        ax3 = axes[2]

        bar_rows = [
            ("TI (Gaussian QD)", dF_c1_gauss, 0.0, dF_c2_gauss, 0.0),
            ("TI (Trapz)",       dF_c1_trap,  0.0, dF_c2_trap,  0.0),
        ]
        if mbar_c1 or mbar_c2:
            c1_v = mbar_c1["dF_mbar"][-1] if mbar_c1 else float("nan")
            c1_e = mbar_c1["dF_std"][-1]  if mbar_c1 else 0.0
            c2_v = mbar_c2["dF_mbar"][-1] if mbar_c2 else float("nan")
            c2_e = mbar_c2["dF_std"][-1]  if mbar_c2 else 0.0
            bar_rows.append(("MBAR", c1_v, c1_e, c2_v, c2_e))

        methods = [r[0] for r in bar_rows]
        c1_vals = np.array([r[1] for r in bar_rows])
        c1_errs = np.array([r[2] for r in bar_rows])
        c2_vals = np.array([r[3] for r in bar_rows])
        c2_errs = np.array([r[4] for r in bar_rows])

        x_pos = np.arange(len(methods))
        w     = 0.35
        c1_plot = np.where(np.isfinite(c1_vals), c1_vals, 0.0)
        c2_plot = np.where(np.isfinite(c2_vals), c2_vals, 0.0)

        ax3.bar(x_pos - w/2, c1_plot, w, yerr=c1_errs, capsize=5,
                color=COLOR_C1, alpha=0.8, label="Case 1")
        ax3.bar(x_pos + w/2, c2_plot, w, yerr=c2_errs, capsize=5,
                color=COLOR_C2, alpha=0.8, label="Case 2")

        for j, (cv1, cv2) in enumerate(zip(c1_vals, c2_vals)):
            if not np.isfinite(cv1):
                ax3.text(x_pos[j] - w/2, 0, "n/a", ha="center",
                         va="bottom", fontsize=8, color=COLOR_C1)
            if not np.isfinite(cv2):
                ax3.text(x_pos[j] + w/2, 0, "n/a", ha="center",
                         va="bottom", fontsize=8, color=COLOR_C2)

        ax3.axhline(0, color="black", linewidth=0.8)
        ax3.set_xticks(x_pos)
        ax3.set_xticklabels(methods, fontsize=11)
        ax3.set_xlabel("λ", fontsize=11)
        ax3.set_ylabel("ΔF_NB  (kcal/mol)", fontsize=11)
        ax3.set_title("NB Free Energy Decomposition: TI (GQ) vs TI (Trapz) vs MBAR",
                      fontsize=11)
        ax3.legend(fontsize=10)
        ax3.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    _savefig(fig, "nb_integrand.png")


def plot_nb_slope(ti: dict) -> None:
    """
    Diagnostic plot: d⟨V_NB⟩/dλ vs λ.
    This is the slope of the NB energy across the alchemical pathway —
    useful to see where NB energy changes rapidly, but NOT the TI integrand.
    """
    lams     = ti["lambdas"]
    slope_c1 = ti["slope_c1"]
    slope_c2 = ti["slope_c2"]

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(lams, slope_c1, "o-", color=COLOR_C1, linewidth=2,
            label="Case 1: d⟨1-4 NB + 1-4 EEL + VDWAALS + EELEC⟩/dλ")
    ax.plot(lams, slope_c2, "s--", color=COLOR_C2, linewidth=2,
            label="Case 2: d⟨VDWAALS + EELEC⟩/dλ")
    ax.axhline(0, color="black", linewidth=0.7, linestyle=":")
    ax.set_xlabel("λ", fontsize=11)
    ax.set_ylabel("⟨dV_NB/dλ⟩  (kcal/mol)", fontsize=11)
    ax.set_title(
        "TI Integrand: ⟨dV_NB/dλ⟩ vs λ  (finite differences of adjacent windows)\n"
        "This IS the TI integrand — ΔF_NB = ∫ ⟨dV/dλ⟩ dλ",
        fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    _savefig(fig, "nb_slope.png")


# ═════════════════════════════════════════════════════════════════════════════
# SAVE SUMMARIES
# ═════════════════════════════════════════════════════════════════════════════

def save_nb_summary(ti: dict, windows: dict, out_dir: str,
                    mbar_c1: dict | None, mbar_c2: dict | None) -> None:
    ranks = sorted(windows.keys())
    rows  = []
    for i, rank in enumerate(ranks):
        wn  = windows[rank]["wn"]
        lam = ti["lambdas"][i]
        gw  = ti["gauss_weights"][i]
        row = {
            "window":              wn,
            "lambda":              lam,
            "gauss_weight":        gw,
            "c1_VNB_mean":         ti["c1_means"][i],
            "c1_VNB_std":          ti["c1_stds"][i],
            "c2_VNB_mean":         ti["c2_means"][i],
            "c2_VNB_std":          ti["c2_stds"][i],
            "c1_dvdl_mean":        ti["integrand_c1"][i],
            "c1_dvdl_std":         ti["integrand_c1_stds"][i],
            "c2_dvdl_mean":        ti["integrand_c2"][i],
            "c2_dvdl_std":         ti["integrand_c2_stds"][i],
            "gq_contrib_c1":       gw * ti["integrand_c1"][i],
            "gq_contrib_c2":       gw * ti["integrand_c2"][i],
        }
        if mbar_c1:
            row["mbar_dF_case1"]     = mbar_c1["dF_mbar"][i]
            row["mbar_dF_std_case1"] = mbar_c1["dF_std"][i]
        if mbar_c2:
            row["mbar_dF_case2"]     = mbar_c2["dF_mbar"][i]
            row["mbar_dF_std_case2"] = mbar_c2["dF_std"][i]
        rows.append(row)

    df   = pd.DataFrame(rows)
    path = os.path.join(out_dir, "nb_summary.csv")
    df.to_csv(path, index=False, float_format="%.6f")
    print(f"  Saved: {path}")

    # ── Pretty-print ─────────────────────────────────────────────────────────
    print("\n" + "═" * 110)
    print("  NB FREE ENERGY DECOMPOSITION SUMMARY  (TI integrand = ⟨dV/dλ⟩)")
    print("═" * 110)
    print(f"  {'λ':>8}  {'w_GL':>7}  {'⟨dV/dλ⟩C1':>12}  {'σ(dV/dλ)':>10}"
          f"  {'GQ·C1':>12}  {'⟨dV/dλ⟩C2':>12}  {'σ(dV/dλ)':>10}  {'GQ·C2':>12}")
    print("  " + "-" * 106)
    for _, row in df.iterrows():
        print(f"  {row['lambda']:8.5f}  {row['gauss_weight']:7.5f}"
              f"  {row['c1_dvdl_mean']:12.4f}  {row['c1_dvdl_std']:10.4f}"
              f"  {row['gq_contrib_c1']:12.4f}"
              f"  {row['c2_dvdl_mean']:12.4f}  {row['c2_dvdl_std']:10.4f}"
              f"  {row['gq_contrib_c2']:12.4f}")

    print(f"\n  ── TI — Gaussian Quadrature  (ΔF_NB = Σᵢ wᵢ·⟨dV_NB/dλ(λᵢ)⟩) ──")
    print(f"     ΔF Case 1 = {ti['dF_c1_gauss']:+.4f} kcal/mol")
    print(f"     ΔF Case 2 = {ti['dF_c2_gauss']:+.4f} kcal/mol")
    print(f"\n  ── TI — Trapezoidal  [comparison] ──")
    print(f"     ΔF Case 1 = {ti['dF_c1_trap']:+.4f} kcal/mol")
    print(f"     ΔF Case 2 = {ti['dF_c2_trap']:+.4f} kcal/mol")

    if mbar_c1:
        print(f"\n  ── MBAR (NB-only, linear Hamiltonian) ──")
        print(f"     ΔF Case 1 = {mbar_c1['dF_mbar'][-1]:+.4f}"
              f"  ±  {mbar_c1['dF_std'][-1]:.4f}  kcal/mol  (λ₀ → λ_max)")
    if mbar_c2:
        if not mbar_c1:
            print(f"\n  ── MBAR (NB-only, linear Hamiltonian) ──")
        print(f"     ΔF Case 2 = {mbar_c2['dF_mbar'][-1]:+.4f}"
              f"  ±  {mbar_c2['dF_std'][-1]:.4f}  kcal/mol  (λ₀ → λ_max)")
    print("═" * 100)


def save_mbar_results(mbar_c1: dict | None, mbar_c2: dict | None,
                      out_dir: str) -> None:
    ref = (mbar_c1 or mbar_c2)
    if ref is None:
        return
    lams = ref["lambdas"]
    rows = []
    for i in range(len(lams)):
        row = {"lambda": lams[i]}
        if mbar_c1:
            row["mbar_dF_case1"]     = mbar_c1["dF_mbar"][i]
            row["mbar_dF_std_case1"] = mbar_c1["dF_std"][i]
        if mbar_c2:
            row["mbar_dF_case2"]     = mbar_c2["dF_mbar"][i]
            row["mbar_dF_std_case2"] = mbar_c2["dF_std"][i]
        rows.append(row)
    df   = pd.DataFrame(rows)
    path = os.path.join(out_dir, "mbar_results.csv")
    df.to_csv(path, index=False, float_format="%.6f")
    print(f"  Saved: {path}")


# ═════════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════════

def main():
    data_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    data_dir = os.path.abspath(data_dir)
    print(f"\nData directory  : {data_dir}")
    print(f"Output directory: {OUT_DIR}")
    print(f"\nGauss-Legendre weights: {len(GAUSS_WEIGHTS)} points, "
          f"sum = {GAUSS_WEIGHTS.sum():.6f}  (should be 1.000000)\n")

    matplotlib.rcParams.update({"font.size": 10})

    print("── Parsing windows ──")
    windows = load_all_windows(data_dir)
    if not windows:
        sys.exit("No matching files found. Expected files named *_w<N>.out")
    print(f"\n  Loaded {len(windows)} window(s).")

    print("\n── Saving per-window CSVs ──")
    save_csvs(windows, OUT_DIR)

    print("\n── Plotting energy timeseries ──")
    plot_timeseries(windows)

    print("\n── Plotting ⟨V_NB⟩ vs λ ──")
    plot_mean_vs_lambda(windows)

    mbar_c1 = mbar_c2 = None

    if len(windows) < 3:
        print("\n[WARNING] Need ≥ 3 λ windows for finite-difference dV/dλ and MBAR.")
        print("  Skipping TI integration and MBAR.")
    else:
        print("\n── Computing TI (NB energy decomposition) ──")
        ti = compute_ti(windows)

        print("\n── Running MBAR (NB-only, linear Hamiltonian) ──")
        mbar_c1 = compute_mbar(windows, case=1)
        mbar_c2 = compute_mbar(windows, case=2)

        print("\n── Plotting TI integrand ⟨V_NB(λ)⟩ ──")
        plot_nb_integrand(ti, mbar_c1, mbar_c2)

        print("\n── Plotting diagnostic d⟨V_NB⟩/dλ ──")
        plot_nb_slope(ti)

        print("\n── Saving summaries ──")
        save_nb_summary(ti, windows, OUT_DIR, mbar_c1, mbar_c2)
        if mbar_c1 or mbar_c2:
            save_mbar_results(mbar_c1, mbar_c2, OUT_DIR)

    print("\nDone.")


if __name__ == "__main__":
    main()