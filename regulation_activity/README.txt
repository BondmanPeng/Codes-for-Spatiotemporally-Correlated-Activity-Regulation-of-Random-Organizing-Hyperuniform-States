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

