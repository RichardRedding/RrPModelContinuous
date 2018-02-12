fn_imputeData_simulated <- function(ls_data,v_r_observedTerm,v_r_paymentRate,v_r_settleRate){
  n_numClaims <- length(v_r_observedTerm)

  fn_truncatedExponential_sim <- function(n,rate,lower,upper){
    -log(1-runif(n)*(1-exp(-rate*(upper-lower))))/rate + lower
  }
  
  setkey(ls_data$paymentTimesObserved,claimIdx)
  
  v_r_timeLatestPayment   <- ls_data$paymentTimesObserved[,max(paymentTime),claimIdx][J(claimIdx = 1:n_numClaims)]$V1
  v_r_timeLatestPayment[is.na(v_r_timeLatestPayment)] <- 0
  v_r_probabilityCensored <- 1-1/(1+(v_r_paymentRate+v_r_settleRate)/(v_r_settleRate*(exp((v_r_observedTerm-v_r_timeLatestPayment)*(v_r_paymentRate+v_r_settleRate))-1)))
  v_r_indicatorCensored   <- rbinom(n_numClaims,1,v_r_probabilityCensored)
  v_r_timeOpen_observed   <- ifelse(v_r_indicatorCensored==1,
                                    v_r_observedTerm,
                                    fn_truncatedExponential_sim(n_numClaims,
                                                                v_r_paymentRate+v_r_settleRate,
                                                                v_r_timeLatestPayment,
                                                                v_r_observedTerm))
  list(v_b_currentlyOpen     = v_r_indicatorCensored,
       v_r_termOpen_observed = v_r_timeOpen_observed)
}

fn_imputeData_expected <- function(ls_data,v_r_observedTerm,v_r_paymentRate,v_r_settleRate){
  n_numClaims <- length(v_r_observedTerm)

  fn_truncatedExponential_exp <- function(rate,lower,upper){
    (1/rate - (upper - lower + 1/rate)*exp(-rate*(upper - lower)))/(1-exp(-rate*(upper - lower))) + lower
  }
  
  setkey(ls_data$paymentTimesObserved,claimIdx)
  
  v_r_timeLatestPayment   <- ls_data$paymentTimesObserved[,max(paymentTime),claimIdx][J(claimIdx = 1:n_numClaims)]$V1
  v_r_timeLatestPayment[is.na(v_r_timeLatestPayment)] <- 0
  v_r_probabilityCensored <- 1-1/(1+(v_r_paymentRate+v_r_settleRate)/(v_r_settleRate*(exp((v_r_observedTerm-v_r_timeLatestPayment)*(v_r_paymentRate+v_r_settleRate))-1)))
  v_r_timeOpen_observed   <- v_r_probabilityCensored*v_r_observedTerm + (1-v_r_probabilityCensored)*fn_truncatedExponential_exp(v_r_paymentRate+v_r_settleRate,
                                                                                                                                v_r_timeLatestPayment,
                                                                                                                                v_r_observedTerm)
  list(v_r_currentlyOpen_prob = v_r_probabilityCensored,
       v_r_termOpen_observed  = v_r_timeOpen_observed)
}
