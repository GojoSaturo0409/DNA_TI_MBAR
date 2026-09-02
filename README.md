# DNA_TI_MBAR

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


