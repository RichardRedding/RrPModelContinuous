\name{fn_imputeData_simulated}
\alias{fn_imputeData_simulated}
\title{fn_imputeData_simulated}
\usage{
fn_imputeData_simulated(ls_data,v_r_observedTerm,v_r_paymentRate,v_r_settleRate)
}
\description{
Given an ls_data list, this function simulates if and when the claim closed during the observed term.
}
\arguments{
\item{ls_data}{A dataset in the form of a list}
\item{v_r_observedTerm}{A vector of reals representing the observed term of each claim}
\item{v_r_paymentRate}{The payment rate of each claim}
\item{v_r_settleRate}{The settlement rate of each claim}
}
