# DNA_TI_MBAR

Data and code availability repository for:

**Thermodynamic Integration and MBAR Analysis of Rad4–DNA Recognition via Free Energy Decomposition**
Rahul Singal and Marimuthu Krishnan
Center for Computational Natural Sciences and Bioinformatics, IIIT-Hyderabad

This repository contains the input files, topologies, raw/processed alchemical free-energy
data, analysis scripts, and supplementary material needed to reproduce the thermodynamic
integration (TI) / Multistate Bennett Acceptance Ratio (MBAR) results reported in the paper.

## Study overview

The paper uses an alchemical TI/MBAR framework (AMBER24) to energetically decompose the
binding free energy of the Rad4–Rad23 nucleotide excision repair complex to a 6-4
photoproduct (6-4PP)-lesioned DNA duplex. Two end states are considered — a **docked**
encounter complex and the **productively bound (PB)** complex — and four structural motifs
are independently alchemically decoupled in each end state to obtain their differential
free-energy contribution, ΔΔF = ΔF<sub>bound</sub> − ΔF<sub>docked</sub>:

| Motif | Residues | ΔΔF (kcal/mol) |
|---|---|---|
| 6-4PP lesion | T64 (602) | 2.83 |
| Partner bases | 584, 585 | 193.48 |
| BHD2 β-hairpin | 392–406 | 68.23 |
| BHD3 β-hairpin | 475–485 | 73.54 |

Free energies were computed with thermodynamic integration across 12 λ windows (10 ns/window,
Gauss–Legendre quadrature) and cross-checked with MBAR (`alchemlyb` + `pymbar`), including
forward/backward convergence, overlap-matrix, and bootstrap-uncertainty diagnostics.

## Repository structure

```
DNA_TI_MBAR/
├── structures/
│   └── 64pp1.pdb                  Full solvated Rad4–DNA system (productively bound state,
│                                   6-4PP lesion + partner bases, protein chain X)
├── amber_inputs/
│   ├── min.in                     Energy minimization (λ = 0.5 hybrid state, example: BHD2 leg)
│   ├── heat.in                    Heating, 0 → 300 K over 0.6 ns (NVT)
│   ├── post_heat.in               NPT equilibration, 2.4 ns, TI/MBAR settings enabled
│   └── ti.in                      Production TI/MBAR run, 10 ns per λ window
├── scripts/
│   ├── ti_analysis.py             TI free-energy + convergence analysis for a single system
│   └── convergence_combined.py    MBAR ΔF-vs-time convergence, run across all 8 systems
└── TI_DATA/
    ├── topologies/                AMBER prmtop/inpcrd + lesion force-field parameters
    ├── TI_simulations/            Per-system raw dV/dλ, MBAR energies, and per-system plots
    │   ├── 64pp_bhd2/             BHD2 hairpin decoupled, productively bound state
    │   ├── 64pp_bhd3/             BHD3 hairpin decoupled, productively bound state
    │   ├── 64pp_lesion/           6-4PP lesion decoupled, productively bound state
    │   ├── 64pp_partner/          Partner bases decoupled, productively bound state
    │   ├── 64pp_docked_bhd2/      BHD2 hairpin decoupled, docked state
    │   ├── 64pp_docked_bhd3/      BHD3 hairpin decoupled, docked state
    │   ├── 64pp_docked_lesion/    6-4PP lesion decoupled, docked state
    │   └── 64pp_docked_partner/   Partner bases decoupled, docked state
    └── MBAR_FINAL_SUPPLEMENTARY/  Final MBAR convergence package (SI figures/tables/text)
```

### `structures/`

`64pp1.pdb` is the full atomistic starting structure of the productively bound Rad4–DNA
complex (28-bp duplex with the 6-4PP lesion at T64 and its two partner bases, plus the
Rad4 protein), used to build the AMBER topologies in `TI_DATA/topologies/`.

### `amber_inputs/`

Representative `pmemd.cuda` input files for one alchemical leg (BHD2 hairpin, residues
392–406; `timask1=':392-406'`), illustrating the minimization → heating → NPT equilibration →
production TI/MBAR protocol used for every motif/state combination. Key settings:

- Soft-core TI (`ifsc=1`, `gti_*` flags), 12 λ windows via `mbar_lambda` (Gauss–Legendre nodes)
- `ntc=1/ntf=1` (no bond constraints) in production, `dt = 1 fs`, `cut = 10 Å`, PME electrostatics
- ff14SB (protein) + parmbsc1 (DNA) force fields, 300 K / 1 atm (Langevin thermostat, weak-coupling barostat)
- Per-window production length: 10 ns; first 20% discarded as equilibration in all downstream analysis

The same templates were reused per motif by changing `timask1`/`scmask1` to the corresponding
residue range and `clambda` to each of the 12 quadrature nodes.

### `scripts/`

- `ti_analysis.py` — parses AMBER `.out` files for a single system (per-λ MBAR energy blocks
  and per-step dV/dλ), computes the TI free energy and MBAR estimate, and produces the
  per-window convergence plots found in each `TI_DATA/TI_simulations/<system>/` folder.
- `convergence_combined.py` — runs MBAR ΔF-vs-simulation-time convergence across all 8
  systems and produces the combined multi-system convergence figure.

Each `TI_DATA/TI_simulations/<system>/` folder also carries its own copy of these scripts
(`ti_analysis.py`, `convergence.py`, `mean_energy.py`, `nb.py`) with the system name/path
hard-coded, exactly as run to produce that folder's outputs — kept alongside the data for
full provenance.

### `TI_DATA/topologies/`

AMBER topology/coordinate files for the three system variants used across the study:

| File | Description |
|---|---|
| `system.prmtop` / `system.inpcrd` | Productively bound (PB) Rad4–DNA complex |
| `system_docked.prmtop` / `system_docked.inpcrd` | Docked encounter complex |
| `system_docked_cleaned.prmtop` / `system_docked_cleaned.inpcrd` | Docked complex, cleaned/re-solvated |
| `inp_2.frcmod` / `inp_2.prepi` | Force-field parameters and RESP charges for the 6-4PP lesion |

### `TI_DATA/TI_simulations/<system>/`

For each of the 8 motif × state alchemical legs (four motifs, docked and productively
bound states), each folder contains:

- `ti_energies_w01.csv` … `ti_energies_w12.csv`, `ti_energies_all.csv` — per-λ-window dV/dλ
  and MBAR energy time series (12 windows, 10 ns each)
- `dvdl_summary.csv`, `mbar_results.csv`, `nb_summary.csv` — per-window summary statistics
- `binding_affinity_results.txt`, `mbar_convergence_results.txt`, `mean_energy_per_lambda.txt` — final ΔF/ΔU values and diagnostics
- `*.png` — dV/dλ, MBAR, block-analysis, and convergence-vs-time plots for that system
- `ti_analysis.py`, `convergence.py`, `mean_energy.py`, `nb.py` — the analysis scripts as run for this system

### `TI_DATA/MBAR_FINAL_SUPPLEMENTARY/`

The final, cross-checked MBAR convergence package underlying the Supporting Information:

- `figures/` — Figures S1–S7 and the main-text convergence figures (PDF + 300 dpi PNG)
- `MBAR_FINAL_CONVERGENCE_TABLE.csv`, `FINAL_MBAR_CONVERGENCE_DECISION.csv` — full per-system diagnostic tables
- `Table_S1.csv` / `Table_S1.tex` — paper-ready summary table (final ΔF, bootstrap uncertainty, forward/backward z-score, median overlap, verdict) for all 8 systems
- `SI_methods.txt`, `SI_results.txt`, `SI_MBAR_METHODS.txt`, `SI_MBAR_RESULTS.txt` — publication-ready SI text
- `FINAL_VERDICT.txt`, `FINAL_CONVERGENCE_VERDICT.txt` — per-system convergence verdict and reasoning
- `analysis_parameters.json` — exact, machine-readable analysis parameters (λ schedule, kBT, decorrelation/convergence criteria, software versions)
- `README.txt` — detailed notes on the MBAR pipeline, including the system-specific handling of the `docked_bhd2` leg (uniform-stride decorrelation; see `analysis_parameters.json` for why)

Convergence summary (bootstrap-calibrated forward/backward z-test, z ≤ 2 → converged):

| System | Final MBAR ΔF (kcal/mol) | Bootstrap uncertainty | Verdict |
|---|---|---|---|
| BHD2 hairpin (bound) | −2157.19 | 13.77 | Converged |
| BHD2 hairpin (docked) | −3218.12 | 8.05 | Borderline |
| BHD3 hairpin (bound) | −3545.07 | 0.74 | Converged |
| BHD3 hairpin (docked) | −3269.92 | 1.71 | Borderline |
| Lesion (bound) | −928.95 | 1.62 | Converged |
| Lesion (docked) | −2979.26 | 0.80 | Converged |
| Partner bases (bound) | −3061.59 | 1.17 | Converged |
| Partner bases (docked) | −2772.80 | — | Converged |

6/8 legs converged; 2/8 (BHD2 and BHD3 hairpins, docked state) are borderline by the
bootstrap-calibrated z-test — see `FINAL_VERDICT.txt` for the full diagnostic reasoning.

## Software

- AMBER24 (`pmemd.cuda`) for all MD/TI simulations
- Python 3.12, `alchemlyb` 2.5.0, `pymbar` 4.2.0 for MBAR analysis
- Force fields: ff14SB (protein), parmbsc1 (DNA), RESP-derived charges for the 6-4PP lesion

## Not included

Raw AMBER trajectory (`.nc`) and per-frame `.out` output files (12 windows × 8 systems ×
10 ns) are not included due to size; the processed dV/dλ and MBAR energy time series in
`TI_DATA/TI_simulations/` and the summary data in `TI_DATA/MBAR_FINAL_SUPPLEMENTARY/` are
sufficient to reproduce every free-energy value and figure in the paper. Contact the
corresponding author for the raw trajectories.

## Citation

If you use this data or code, please cite:

> Singal, R.; Krishnan, M. Thermodynamic Integration and MBAR Analysis of Rad4–DNA
> Recognition via Free Energy Decomposition. (manuscript in preparation).

## Contact

Marimuthu Krishnan — m.krishnan@iiit.ac.in
Center for Computational Natural Sciences and Bioinformatics, IIIT-Hyderabad
