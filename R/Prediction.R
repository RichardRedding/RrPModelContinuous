predict_numFutureEvents <- function(staticTerm, observedTerm, predictionTerm, p, s){
  fn_imputeData_expected(staticTerm, observedTerm, p, s)$censored *
    fn_truncatedExponential_exp(s, 0, predictionTerm) * p
}