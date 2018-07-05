fn_parameter_mle <- function(v_n_numPayments,v_r_timeOpen,v_r_currentlyOpen){
  c(p = sum(v_n_numPayments)/sum(v_r_timeOpen),
    s = sum(1-v_r_currentlyOpen)/sum(v_r_timeOpen))
}

fn_parameter_mleEM <- function(staticTerm, observedTerm, numObservations){
  n_numClaims <- length(observedTerm)
  v_r_parameter_mle <- c(p=mean(pmax(0,numObservations-1))/mean(observedTerm-staticTerm),
                         s=1/mean(observedTerm-staticTerm))
  stopping_scale <- max(v_r_parameter_mle)
  eps <- Inf
  while(sum(abs(eps)) > stopping_scale/1000){
    ls_dataAugmented_exp <- fn_imputeData_expected(staticTerm,
                                                   observedTerm,
                                                   v_r_parameter_mle["p"],
                                                   v_r_parameter_mle["s"])
    eps <- fn_parameter_mle(numObservations,
                            ls_dataAugmented_exp$termOpen,
                            ls_dataAugmented_exp$censored)-v_r_parameter_mle
    v_r_parameter_mle <- v_r_parameter_mle + eps
  }
  v_r_parameter_mle
}

fn_parameter_sim <- function(n_numSims, staticTerm, observedTerm, numObservations, v_r_paymentPriorParams, v_r_settlementPriorParams){
  n_numClaims <- length(observedTerm)
  v_r_parameter_init <- fn_parameter_mleEM(ls_data,observedTerm)
  v_r_parameter_sim  <- matrix(rep(v_r_parameter_init,each=n_numSims),n_numSims,2)
  ptm <- proc.time()["elapsed"]
  for(i in 1:(n_numSims-1)){
    ls_missingData                <- fn_imputeData_simulated(staticTerm,
                                                             observedTerm,
                                                             v_r_parameter_sim[i,1],
                                                             v_r_parameter_sim[i,2])
    v_r_paymentPosteriorParams    <- v_r_paymentPriorParams    + c(sum(numObservations), sum(ls_missingData$v_r_termOpen_observed))
    v_r_settlementPosteriorParams <- v_r_settlementPriorParams + c(sum(1-ls_missingData$v_b_currentlyOpen), sum(ls_missingData$v_r_termOpen_observed))
    v_r_parameter_sim[i+1,]       <- c(rgamma(1,v_r_paymentPosteriorParams[1],v_r_paymentPosteriorParams[2]),
                                       rgamma(1,v_r_settlementPosteriorParams[1],v_r_settlementPosteriorParams[2]))
    if((i %% ceiling(n_numSims/min(20,n_numSims)))==0) {print(paste((n_numSims-i)*(proc.time()["elapsed"]-ptm)/i,"seconds remaining"))}
  }
  v_r_parameter_sim
}

f_llikelihood <- function(staticTerm, observedTerm, numObservations, v_r_parameters){
  w <- sum(v_r_parameters)
  p <- v_r_parameters[1]
  #( (exp(-w*l)-exp(-w*u))*s/w + exp(-w*u) )*(p^n)
  #(1-(1-exp(w*(u-l)))*s/w)*exp(-w*u)*(p^n)
  log(1 - (1 - exp(w * staticTerm)) * (1 - p / w)) - w*observedTerm + numObservations * log(p)
}
