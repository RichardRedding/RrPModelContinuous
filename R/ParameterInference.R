fn_parameter_mle <- function(v_n_numPayments,v_r_timeOpen,v_r_currentlyOpen){
  c(p = sum(v_n_numPayments)/sum(v_r_timeOpen),
    s = sum(1-v_r_currentlyOpen)/sum(v_r_timeOpen))
}

fn_parameter_mleEM <- function(ls_data,v_r_observedTerm){
  n_numClaims <- length(v_r_observedTerm)
  v_r_parameter_mle <- c(p=.5,s=.5)
  eps <- Inf
  while(sum(abs(eps)) > .001){
    ls_dataAugmented_exp <- fn_imputeData_expected(ls_data,
                                                   v_r_observedTerm,
                                                   v_r_parameter_mle["p"],
                                                   v_r_parameter_mle["s"])
    eps <- fn_parameter_mle(ls_data$numPaymentsObserved,
                            ls_dataAugmented_exp$v_r_termOpen_observed,
                            ls_dataAugmented_exp$v_r_currentlyOpen_prob)-v_r_parameter_mle
    v_r_parameter_mle <- v_r_parameter_mle + eps
  }
  v_r_parameter_mle
}

fn_parameter_sim <- function(n_numSims,ls_data,v_r_observedTerm,v_r_paymentPriorParams,v_r_settlementPriorParams){
  n_numClaims <- length(v_r_observedTerm)
  v_r_parameter_init <- fn_parameter_mleEM(ls_data,v_r_observedTerm)
  v_r_parameter_sim  <- matrix(rep(v_r_parameter_init,each=n_numSims),n_numSims,2)
  ptm <- proc.time()["elapsed"]
  for(i in 1:(n_numSims-1)){
    ls_missingData                <- fn_imputeData_simulated(ls_data,
                                                             v_r_observedTerm,
                                                             v_r_parameter_sim[i,1],
                                                             v_r_parameter_sim[i,2])
    v_r_paymentPosteriorParams    <- v_r_paymentPriorParams    + c(sum(ls_data$numPaymentsObserved), sum(ls_missingData$v_r_termOpen_observed))
    v_r_settlementPosteriorParams <- v_r_settlementPriorParams + c(sum(1-ls_missingData$v_b_currentlyOpen), sum(ls_missingData$v_r_termOpen_observed))
    v_r_parameter_sim[i+1,]       <- c(rgamma(1,v_r_paymentPosteriorParams[1],v_r_paymentPosteriorParams[2]),
                                       rgamma(1,v_r_settlementPosteriorParams[1],v_r_settlementPosteriorParams[2]))
    if((i %% ceiling(n_numSims/min(20,n_numSims)))==0) {print(paste((n_numSims-i)*(proc.time()["elapsed"]-ptm)/i,"seconds remaining"))}
  }
  v_r_parameter_sim
}

f_llikelihood <- function(ls_data,v_r_observedTerm,v_r_parameters){
  setkey(ls_data$paymentTimesObserved,claimIdx)
  l <- ls_data$paymentTimesObserved[,max(paymentTime),claimIdx][J(claimIdx = 1:n_numClaims)]$V1
  l[is.na(l)] <- 0
  w <- sum(v_r_parameters)
  p <- v_r_parameters[1]
  u <- v_r_observedTerm
  n <- ls_data$numPaymentsObserved
  #( (exp(-w*l)-exp(-w*u))*s/w + exp(-w*u) )*(p^n)
  #(1-(1-exp(w*(u-l)))*s/w)*exp(-w*u)*(p^n)
  log(1-(1-exp(w*(u-l)))*(1-p/w)) -w*u + n*log(p)
}
