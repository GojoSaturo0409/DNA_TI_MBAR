MBAR_FINAL_SUPPLEMENTARY -- README
====================================
Generated: 2026-08-29 20:53:18

CONTENTS
--------
  figures/                          -- Figures S1-S7 (PDF + PNG, 300 dpi)
  MBAR_FINAL_CONVERGENCE_TABLE.csv  -- full per-system diagnostic table
  Table_S1.csv / Table_S1.tex       -- compact paper-ready table
  SI_methods.txt                    -- publication-ready Methods section
  SI_results.txt                    -- publication-ready Results section
  analysis_parameters.json          -- exact parameters used (machine-readable)
  FINAL_CONVERGENCE_VERDICT.txt     -- per-system verdict + overall conclusion
  README.txt                        -- this file

INPUT FILES
-----------
Raw AMBER .out files, 12 lambda windows per system, in
  TI/64pp_<system>/64pp_<system>_w<N>.out   (N=1..12)
(64pp_lesion and 64pp_docked_lesion use "64pp"/"64pp_docked" prefixes,
matching every existing script in this repository.)

SCRIPTS USED (this analysis)
-----------------------------
  fixed_loader.py           -- sentinel-corrected AMBER MBAR-block parser
                                (overflow -> 1e15 kcal/mol sentinel, not dropped)
  mbar_pub.py / mbar_pub_robust.py
                             -- forward/backward MBAR (num=10) + overlap +
                                adjacent-pair diagnostics for 7 systems
                                (bhd2, bhd3, docked_bhd3, docked_lesion,
                                docked_partner, lesion, partner)
  fix_overlap.py             -- corrected overlap-matrix extraction
                                (pymbar 4.2's real method is
                                MBAR.compute_overlap(), not
                                compute_overlap_matrix(), which does not
                                exist in this pymbar version)
  docked_bhd2_audit/*.py     -- forensic audit that identified and
                                localized the docked_bhd2 anomaly
                                (see docked_bhd2_audit/docked_bhd2_mbar_audit_report.txt)
  dbhd2_validated.py         -- VALIDATED docked_bhd2 pipeline: uniform
                                stride (5/10/20/40), forward/backward at
                                stride=10, overlap, bootstrap, block
                                uncertainty
  build_all.py, make_figures.py, make_tables.py, write_text_deliverables.py
                             -- assembly of this final package

SOFTWARE VERSIONS
-----------------
  Python:     3.12 (conda environment "myenv")
  alchemlyb:  2.5.0
  pymbar:     4.2.0

PHYSICAL / SIMULATION PARAMETERS
----------------------------------
  Temperature:            300.0 K
  kBT:                    0.596161 kcal/mol
  Lambda windows:         12  (0.0092, 0.0479, 0.1150, 0.2063, 0.3161, 0.4374,
                                0.5626, 0.6839, 0.7937, 0.8850, 0.9521, 0.9908)
  Simulation length:      10 ns/window
  Equilibration discarded: first 20% (2 ns) of each window
  Raw frames/window:      10000
  Production frames/window: 8000

UNITS / REDUCED POTENTIAL
---------------------------
  Raw AMBER "Energy at lambda" block values are kcal/mol (verified: a
  step's own-state block value matches that step's separately-reported
  EPtot to printed precision). alchemlyb's MBAR.fit() requires input
  already divided by kBT (its docstring: "u_nk[n,k] is the reduced
  potential energy"). This analysis divides by kBT before fitting and
  multiplies the (correctly dimensionless) output by kBT for reporting.

CONVERGED SAMPLE HANDLING
---------------------------
  7 systems: per-window independent decorrelation (pymbar
  statistical_inefficiency + subsample_correlated_data on each window's
  own diagonal MBAR-block energy series).
  docked_bhd2: uniform stride=10 across all windows (validated pipeline;
  see analysis_parameters.json and docked_bhd2_audit/ for why per-window
  independent decorrelation was rejected for this specific system).

CONVERGENCE SETTINGS
----------------------
  Forward/backward: 10 fractions (10%, 20%, ..., 100%), each an
  independent full 12-window MBAR fit. Compared at the 90% fraction
  (first-90% vs last-90%) as a percentage of the total free energy.
  Threshold: <=0.5% CONVERGED, 0.5-1.5% BORDERLINE, >1.5% NOT CONVERGED.

UNCERTAINTY METHOD
---------------------
  Primary: analytic MBAR uncertainty (asymptotic covariance).
  Cross-checks (docked_bhd2): bootstrap (n=200 resamples) and 5-block
  time-series uncertainty. Analytic uncertainty is systematically
  tighter than both cross-checks (bootstrap ~11x larger for
  docked_bhd2) -- reported explicitly rather than used as if it were
  the complete picture of estimator uncertainty.

SUMMARY
-------
  6 of 8 systems: CONVERGED
  2 of 8 systems: BORDERLINE
  0 of 8 systems: NOT CONVERGED
  See FINAL_CONVERGENCE_VERDICT.txt for the full per-system breakdown.
