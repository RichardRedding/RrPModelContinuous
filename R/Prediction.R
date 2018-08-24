predict_numFutureEvents <- function(staticTerm, observedTerm, predictionTerm, v_r_paymentRate, v_r_settleRate){
  fn_imputeData_expected(inactiveTerm, observedTerm, v_r_paymentRate, v_r_settleRate)$censored *
    fn_truncatedExponential_exp(v_r_settleRate, 0, predictionTerm) * v_r_paymentRate
}