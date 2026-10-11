REGULATION ACTIVITY

available_data contains only CSV data for the three requested BRO figures and
the three requested IG figures.

available_data/BRO/ (updated to number-density data)
  structure_factors_c_h_rhoA0.8_rhoB0.7_h01.5_source_data.csv
    S_AA(q) and S_BB(q): six c curves at h=1 and seven h curves at c=0.
  delta_SBB_c_h_sweeps_rhoA0.8_rhoB0.7_h01.5_source.csv
    All displayed Delta S_B(q) points, S_BB and matched h=0 reference.
  sq_SAA_c_rhoA0.8_rhoB0.7_h01.5_source.csv
    Six displayed S_A(q) curves for the c sweep.
  sq_SBB_delta_c_h_combined_rhoA0.8_rhoB0.7_h01.5_source.csv
    Displayed panels A-D: 2862 raw A points, 5 B plateaus,
    3339 raw C points and 6 D plateaus. B x=S_AA0 from predictions.
  delta_SBB_BD_common_denominator_fit_predictions.csv
    Per-run plateau means, predictions, residuals and unique fit membership.
    included_in_fit marks uniquely counted runs; shared h=1,c=0 is counted
    in panel B only. c=-1 has no prediction guide; h=0 is not fitted.
  delta_SBB_BD_common_denominator_fit_summary.csv
    P=bar(D_B)(D_A+bar(D_B))=5.6072389091486485; 10 unique fit points.
  delta_SBB_BD_common_denominator_fit_curves.csv
    Exact continuous fit curves from panels B/D (300 points per panel).
    D curve domain h=0.2..1.2 and mean S_A(0) convention match the figure.

available_data/IG/ is organized into exactly three figure folders:
  sq_SA_vary_DA/
    sq_SA_vary_DA_source.csv: plotted S_A(q) curves, with D_A values.
  sq_SBB_delta_combined/
    sq_SBB_delta_combined_source.csv: displayed panels A-D.
    sq_delta_D_B_fit.csv and sq_delta_D_B_fit_predictions.csv:
    the 500-frame dataset's shared fit and pointwise predictions.
  sq_SBB_minus_h0_T_h/
    sq_SBB_minus_h0_T_h_styled_source.csv: all displayed Delta S_B(q)
    points and horizontal guide values, from the 302-frame dataset.
    sq_delta_D_B_fit.csv and sq_delta_D_B_fit_predictions.csv:
    separate D_A/h fits and predictions for this figure.
Do not exchange fit files between the two IG S_B figure folders.
All packaged CSVs use D_A directly: sweep=D_A and value=D_A.
D_A=T/2 has already been applied; do not divide value by two again.
Temperature columns have been removed; D_A_theory was renamed D_A.
Original filenames and run identifiers retain their provenance spelling.
Original project CSVs are unchanged; conversion applies to available_data.
No other IG figures or datasets are included.
No trajectories, raw-run duplicates, QA text or unrelated datasets remain in available_data.
Simulation source code files were not changed by this data update.

BRO data: number densities rho_A=0.8, rho_B=0.7; h0=1.5, R_A=3,
tau_ratio=1, epsilon_A=1, 100000 steps, 500 snapshots strictly after step 50000 (50100..100000),
and 477 q bins/run. Twelve distinct runs; the h=1,c=0 run is shared.
c=-1,-0.8,-0.6,-0.4,-0.2,0 at h=1; h=0,0.2,0.4,0.6,0.8,1,1.2 at c=0.
Delta S_B(q)=S_B(q;h,c)-S_B(q;h=0,c=0), on matched q grids.
The baseline is reusable across c because B decouples from A at h=0.
Plateaus use all four bins with q sigma < 0.1. No smoothing or refitting
was performed during packaging. Values are copied from the current fit.
Shared fit uses the five B points c>-1 and h=0.2..1.2, counting h=1,c=0
once. c=-1 remains in raw A and S_A/delta curves but is excluded from B/fit.
h=0 remains in raw C and delta curves but is excluded from D/fit.
The new h sweep ends at 1.2; h=1.4 and 1.6 are not part of this dataset.
Delta figure y limits are (-0.0144914814,0.01) for c and (-0.014,0.01) for h; the low-q dashed guides use the
packaged per-run predictions. The original CSV schemas are retained;
parameter-bearing filenames now identify the new density and h0 values.
Updated 2026-10-11 from the all-100000-step dataset.
Source directory: /home/peng0153/Yuhan/programs/modulated_BRO/new_plots/number_density/completed_combined_steps100000/rhoA0.8_rhoB0.7_h01.5
Sources: source_data.csv, plateaus_and_predictions.csv, fit_summary.csv,
fit_curves.csv, sq_delta_SB_c_h_rhoA0.8_rhoB0.7_h01.5_source.csv,
and sq_SA_c_rhoA0.8_rhoB0.7_h01.5_source.csv.
Original project data and IG package data remain unchanged.
Previous package archived at: /home/peng0153/Yuhan/programs/regulation_activity_backups/before_BRO_rhoA08_steps100000_20261011_231035
Current simulation code defaults h0 to 2; these data use explicit h0=1.5.
Both packaged simulation codes interpret rho as number density.

IG combined figure panel D omits the h=0 measured point; panel C retains its h=0 curve. The existing fit is unchanged.

IG combined figure panel A includes D_A=100 (stored value=100, sweep=D_A), with all 477 q bins.
