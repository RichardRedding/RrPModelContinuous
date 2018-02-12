fn_pseudoData <- function(n_numClaims,v_r_observedTerm,v_r_paymentRate,v_r_settleRate){
  v_r_timeOpen                <- rexp(n_numClaims, v_r_settleRate)
  
  v_r_timeOpen_censored       <- pmin(v_r_timeOpen, v_r_observedTerm)
  
  v_n_numPayments_censored    <- rpois(n_numClaims, v_r_timeOpen_censored*v_r_paymentRate)
  
  v_n_numPayments    <- v_n_numPayments_censored + rpois(n_numClaims, v_r_paymentRate*(v_r_timeOpen-v_r_timeOpen_censored))
  
  dt <- data.table(claimIdx = unlist(lapply(1:n_numClaims,function(x) rep(x,v_n_numPayments[x]))))
  
  dt[,paymentIdx := seq(1,.N), claimIdx]
  
  dt[,observedInd := ifelse(paymentIdx <= v_n_numPayments_censored[claimIdx],1,0)]
  
  dt[observedInd == 1,paymentTime := runif(.N,0,v_r_timeOpen_censored[claimIdx])]
  
  dt[observedInd == 0,paymentTime := runif(.N,v_r_observedTerm[claimIdx],v_r_timeOpen[claimIdx])]
  
  list(paymentTimesObserved = dt[observedInd == 1,.(claimIdx,paymentIdx,paymentTime)],
       numPaymentsObserved = v_n_numPayments_censored,
       paymentTimes = dt[,.(claimIdx,paymentIdx,paymentTime)],
       numPayments = v_n_numPayments,
       openTime = v_r_timeOpen)
}