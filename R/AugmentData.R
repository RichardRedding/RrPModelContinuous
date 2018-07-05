fn_truncatedExponential_sim <- function(n,rate,lower,upper){
  -log(1-runif(n)*(1-exp(-rate*(upper-lower))))/rate + lower
}

fn_imputeData_simulated <- function(staticTerm, observedTerm, observationRate, settlementRate){
  numClaims <- length(staticTerm)
  probabilityCensored <- 1-1/(1+(observationRate+settlementRate)/(settlementRate*(exp(staticTerm*(observationRate+settlementRate))-1)))
  censored <- rbinom(numClaims,1,v_r_probabilityCensored)
  termOpen   <- ifelse(censored==1,
                       observedTerm,
                       fn_truncatedExponential_sim(numClaims,
                                                   observationRate+settlementRate,
                                                   observedTerm-staticTerm,
                                                   observedTerm))
  list(censored = censored,
       termOpen = termOpen)
}

fn_truncatedExponential_exp <- function(rate,lower,upper){
  (1/rate - (upper - lower + 1/rate)*exp(-rate*(upper - lower)))/(1-exp(-rate*(upper - lower))) + lower
}

fn_imputeData_expected <- function(staticTerm, observedTerm, observationRate, settlementRate){
  probabilityCensored <- 1-1/(1+(observationRate+settlementRate)/(settlementRate*(exp(staticTerm*(observationRate+settlementRate))-1)))
  termOpen <- probabilityCensored*observedTerm + (1-probabilityCensored)*fn_truncatedExponential_exp(observationRate+settlementRate,
                                                                                                     observedTerm-staticTerm,
                                                                                                     observedTerm)
  list(censored = probabilityCensored,
       termOpen = termOpen)
}
