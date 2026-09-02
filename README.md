

Data and code availability repository for:

**Energetic and Dynamical Basis of Rad4
Recognition of UV-Induced DNA Damage**

Rahul Singal and Marimuthu Krishnan

Center for Computational Natural Sciences and Bioinformatics, IIIT-Hyderabad

This repository contains the input files, topologies, raw/processed alchemical free-energy
data, analysis scripts, and supplementary material needed to reproduce the thermodynamic
integration (TI) / Multistate Bennett Acceptance Ratio (MBAR) results reported in the paper.



## Repository structure

```
DNA_TI_MBAR/
├── structures/
│   ├── ds_non_solvated.pdb        Full Rad4–DNA system (non-solvated, docked state,
│   │                               6-4PP lesion + partner bases, protein chain X)
│   └── pb_non_solvated.pdb        Full Rad4–DNA system (non-solvated, productively bound state)
├── amber_inputs/
│   ├── min.in                     Energy minimization (λ = 0.5 hybrid state, example: BHD2 leg)
│   ├── heat.in                    Heating, 0 → 300 K over 0.6 ns (NVT)
│   ├── post_heat.in               NPT equilibration, 2.4 ns, TI/MBAR settings enabled
│   └── ti.in                      Production TI/MBAR run, 10 ns per λ window
└── TI_DATA/
    ├── topologies/                AMBER prmtop/inpcrd + lesion force-field parameters
    └── TI_simulations/            Per-system TI/MBAR energy data (CSV only)
        ├── 64pp_bhd2/             BHD2 hairpin decoupled, productively bound state
        ├── 64pp_bhd3/             BHD3 hairpin decoupled, productively bound state
        ├── 64pp_lesion/           6-4PP lesion decoupled, productively bound state
        ├── 64pp_partner/          Partner bases decoupled, productively bound state
        ├── 64pp_docked_bhd2/      BHD2 hairpin decoupled, docked state
        ├── 64pp_docked_bhd3/      BHD3 hairpin decoupled, docked state
        ├── 64pp_docked_lesion/    6-4PP lesion decoupled, docked state
        └── 64pp_docked_partner/   Partner bases decoupled, docked state
```

Each `TI_simulations/<system>/` folder contains `ti_energies_w01.csv` … `ti_energies_w12.csv`
(raw per-λ-window energies) and `ti_energies_all.csv` (combined), plus `mbar_results.csv`
and `dvdl_summary.csv`/`nb_summary.csv` where applicable — the CSV data needed to reproduce
the TI/MBAR free-energy and convergence analysis.


