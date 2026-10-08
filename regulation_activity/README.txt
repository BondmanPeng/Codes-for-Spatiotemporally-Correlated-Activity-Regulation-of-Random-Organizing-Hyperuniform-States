REGULATION ACTIVITY

available_data contains only CSV data for the three requested BRO figures and
the three requested IG figures.

available_data/BRO/
  structure_factors_c_h_rhoA0.8_source_data.csv
    Full S_AA(q) and S_BB(q) for all six c values and nine h values.
  delta_SBB_c_h_sweeps_rhoA0.8_source.csv
    All six c and nine h Delta S_B(q) curves, S_BB and the subtracted h=0 reference.
  sq_SAA_c_rhoA0.8_source.csv
    All displayed S_A(q) curves for the six c values.
  sq_SBB_delta_c_h_combined_rhoA0.8_source.csv
    All displayed structure-factor and plateau data in panels A-D.
  delta_SBB_BD_common_denominator_fit_predictions.csv
    Per-run plateau predictions, fit membership and residuals; also supplies
    the horizontal guides for the Delta S_B figure.
  delta_SBB_BD_common_denominator_fit_summary.csv
    Common fitted denominator and fit metrics.

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
