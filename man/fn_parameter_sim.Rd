\name{fn_parameter_sim}
\alias{fn_parameter_sim}
\title{fn_parameter_sim}
\usage{
fn_parameter_sim(n_numSims,ls_data,v_r_observedTerm,v_r_paymentPriorParams,v_r_settlementPriorParams)
}
\description{
Simulates from the posterior distribution of the payment and settlement rates. The prior for the payment and settlement rates are independent gamma distributions.
}
\arguments{
\item{ls_data}{A dataset in the form of a list}
\item{v_r_observedTerm}{A vector of reals representing the observed term of each claim}
\item{n_numSims}{The number of simulations to generate}
\item{v_r_paymentPriorParams}{The parameters of the prior distribution of the payment rate}
\item{v_r_settlementPriorParams}{The parameters of the prior distribution of the settlement rate}
}

