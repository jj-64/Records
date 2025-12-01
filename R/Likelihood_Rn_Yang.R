############## Exponential #################
#' Log-Likelihood for Record Process with Exponential Distribution under YNM Model
#'
#' Computes the log-likelihood of a record process where record values follow
#' an Exponential distribution, and record arrival probabilities follow the YNM-Nevzorov model.
#'
#' @param R A numeric vector of record values.
#' @param L A numeric vector of time indices at which the records occur.
#' @param T An integer indicating the total number of time points.
#' @param params A numeric vector of parameters: \eqn{\gamma} (YNM parameter), and 1/rate of the exponential distribution.
#'
#' @return A numeric value representing the log-likelihood.
#' @export
#'
#' @examples
#' Likelihood_Rn_Exp_YNM(R = c(1.2, 2.5, 3.8), L = c(2, 4, 6), T = 10, params = c(0.8, 2))
Likelihood_Rn_Exp_YNM <- function(R,L,T, params) {  ## vector of R, L, value of T, vector of (gamma, 1/rate)
  m=length(R)

  # Calculate Sum1
  s1 <- log(params[1]) * sum(L)

  # Loop to calculate Sum2
  s2a=0
  for (j in 1:m) {
    s2a[j]= log(dexp(R[j], rate = 1/params[2]))
  }
  s2=sum(s2a)

  # Loop to calculate Sum3
  s3a=0
  for (j in 1:m) {
    s3a[j]= + ((params[1]^L[j]) - 1) * log(pexp(R[j], rate = 1/params[2]))
  }
  s3=sum(s3a)

  # Loop to calculate Sum4
  s4a=0
  for (j in 1:(m - 1)) {
    s4a[j] = (((params[1]^(L[j] + 1)) - (params[1]^L[j + 1])) / (1 - params[1])) *
      log(pexp(R[j], rate = 1/params[2]))
  }
  s4=sum(s4a)

  # Calculate Sum5
  s5=0
  if(L[m]<T){
    s5 <- (((params[1]^(L[m] + 1)) - (params[1]^(T + 1))) / (1 - params[1])) *
      log(pexp(R[m], rate = 1/params[2]))
  }

  # Sum all parts to calculate `vrais`
  vrais <- s1+s2+s3+s4+s5

  # Return the result
  return(vrais)
}

#' Likelihood for Record Process with Exponential Distribution under YNM Model
#'
#' Computes the likelihood (not log) of a record process assuming
#' Exponential record values and YNM-type record frequency.
#'
#' @inheritParams Likelihood_Rn_Exp_YNM
#'
#' @return A numeric value representing the likelihood.
#' @export
#'
#' @examples
#' L_Rn_Exp_YNM(R = c(1.2, 2.5, 3.8), L = c(2, 4, 6), T = 10, params = c(0.8, 2))
L_Rn_Exp_YNM = function(R,L,T, params){

  m = length(R)

  p1 = params[1]^(sum(L))

  p2a= dexp(x=R,rate=1/params[2])
  p2b =pexp(q=R,rate=1/params[2])^(-1+params[1]^L)
  p2=prod(p2a*p2b)

  p3a=1
  for(i in 1:(m-1)){
    p3a[i]=(pexp(q=R[i],rate=1/params[2]))^((params[1]^(L[i]+1) - params[1]^(L[i+1]))*(1-params[1])^-1)
  }
  p3b=1
  if(L[m]<T){
    p3b = pexp(q=R[m],rate=1/params[2]) ^(  (params[1]^(L[m]+1) - params[1]^(T+1))/(1-params[1]) )
  }
  p3=prod(p3a)*prod(p3b)

  p = p1*p2*p3

  return(p)
}

############### Frechet #####################
#' Log-Likelihood for Record Process with Fréchet Distribution under YNM Model
#'
#' Calculates the log-likelihood assuming the record values follow a Fréchet distribution
#' with YNM-Nevzorov model record probabilities.
#'
#' @param R A numeric vector of record values.
#' @param L A numeric vector of time indices of the records.
#' @param T Integer, total observation time or maximum index.
#' @param params A numeric vector: \eqn{\gamma} (YNM parameter), scale, and shape of the Fréchet distribution.
#'
#' @return The log-likelihood as a numeric value, or \code{-Inf} if parameters are invalid.
#' @export
#'
#' @examples
#' Likelihood_Rn_Frechet_YNM(R = c(3.5, 2.1, 1.8), L = c(2, 4, 6), T = 10, params = c(0.9, 1, 2))
Likelihood_Rn_Frechet_YNM = function(R,L,T,params){
  m = length(R) #m= rec_counts(y)  ## number of records

  # Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or params[3] <= 0, we return -Inf)
  if (any(params <= 0)) {
    return(-Inf)  # Invalid parameters, return a large negative value
  }

  # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
  s1 = m *log(params[3] *params[2]^(-params[3]))
  if (is.nan(s1) || !is.finite(s1)) return(-Inf)

  s2 = -(params[3] +1) * sum(log(R))
  if (is.nan(s2) || !is.finite(s2)) return(-Inf)

  s3 = sum(L)* log (params[1])
  if (is.nan(s3) || !is.finite(s3)) return(-Inf)

  s4 =  -sum(  ((params[1]^L[-length(L)]) - (params[1]^L[-1]) ) * (params[2] * R[-length(R)] )^(-params[3])  /(1-params[1]))
  if (is.nan(s4) || !is.finite(s4)) return(-Inf)

  s5= 0
  if(L[m]<T){
    s5 = -(params[1]^L[m] - params[1]^(T+1)) * (params[2] * R[m])^(-params[3]) / (1 - params[1])
    if (is.nan(s5) || !is.finite(s5)) return(-Inf)
  }
  return(s1+s2+s3+s4+s5)
}

#################### Gumbel ###############

#' Log-Likelihood for Record Process with Gumbel Distribution under YNM Model
#'
#' Computes the log-likelihood for records modeled using the Gumbel distribution
#' with the YNM-Nevzorov mechanism for record arrival.
#'
#' @param R A numeric vector of record values.
#' @param L A numeric vector of time indices of the records.
#' @param T Integer indicating the total number of observations.
#' @param params A numeric vector: \eqn{\gamma} (YNM parameter), location \eqn{\mu}, and scale \eqn{\beta} of the Gumbel distribution.
#'
#' @return Log-likelihood value as a numeric scalar.
#' @export
#'
#' @examples
#' Likelihood_Rn_Gumbel_YNM(R = c(2.4, 3.0, 4.1), L = c(1, 3, 7), T = 10, params = c(0.8, 0, 1))
Likelihood_Rn_Gumbel_YNM <- function(R,L,T, params) {  ## vector of R, L, value of T, vector of gamma, location mu and scale beta
  m=length(R)

  # Calculate Sum1
  s1 <- log(params[1]) * sum(L)

  # Loop to calculate Sum2
  s2a=0
  for (j in 1:m) {
    s2a[j]= log(VGAM::dgumbel(R[j], loc = params[2], scale = params[3]))
  }
  s2=sum(s2a)

  # Loop to calculate Sum3
  s3a=0
  for (j in 1:m) {
    s3a[j]= + ((params[1]^L[j]) - 1) * log(VGAM::pgumbel(R[j], loc = params[2], scale = params[3]))
  }
  s3=sum(s3a)

  # Loop to calculate Sum4
  s4a=0
  for (j in 1:(m - 1)) {
    if((Ln[j]+1 <= Ln[j+1]-1) == TRUE){ ## we have non-records in between
    s4a[j] = (((params[1]^(L[j] + 1)) - (params[1]^L[j + 1])) / (1 - params[1])) *
      log(VGAM::dgumbel(R[j], location = params[2], scale = params[3]))
  }}
  s4=sum(s4a)

  # Calculate Sum5
  s5=0
  if(L[m]<T){
    s5 <- (((params[1]^(L[m] + 1)) - (params[1]^(T + 1))) / (1 - params[1])) *
      log(VGAM::pgumbel(R[m], loc = params[2], scale = params[3]))
  }

  # Sum all parts to calculate `vrais`
  vrais <- s1+s2+s3+s4+s5

  # Return the result
  return(vrais)
}

#' Likelihood for Record Process with Gumbel Distribution under YNM Model
#'
#' Computes the likelihood (non-log) assuming Gumbel-distributed record values
#' and record occurrence governed by the YNM process.
#'
#' @inheritParams Likelihood_Rn_Gumbel_YNM
#'
#' @return Numeric value representing the likelihood.
#' @export
#'
#' @examples
#' L_Rn_Gumbel_YNM(R = c(2.4, 3.0, 4.1), L = c(1, 3, 7), T = 10, params = c(0.8, 0, 1))
L_Rn_Gumbel_YNM = function(R,L,T, params){

  m = length(R)

  p1 = params[1]^(sum(L))

  p2a= dgumbel(x=R,loc=params[2],scale=params[3])
  p2b =pgumbel(q=R,loc=params[2],scale=params[3])^(-1+params[1]^L)
  p2=prod(p2a*p2b)

  p3a=0
  for(i in 1:(m-1)){
    p3a[i]=(pgumbel(q=R[i],loc=params[2],scale=params[3]))^((params[1]^(L[i]+1) - params[1]^(L[i+1]))*(1-params[1])^-1)
  }
  p3b=1
  if(L[m]<T){
    p3b = (pgumbel(q=R[m],loc=params[2],scale=params[3]))^((params[1]^(L[m]+1) - params[1]^T)*(1-params[1])^-1)
  }
  p3=prod(p3a)*prod(p3b)

  p = p1*p2*p3

  return(p)
}

############## Weibull ##################
#' Log-Likelihood for Record Process with Weibull Distribution under YNM Model
#'
#' Calculates the log-likelihood for a process of records with Weibull-distributed values
#' and record times based on the YNM-Nevzorov model.
#'
#' @param R A numeric vector of record values.
#' @param L A numeric vector of time indices of the records.
#' @param T Integer representing the total number of observations.
#' @param params A numeric vector of parameters: \eqn{\gamma}, shape, and scale of the Weibull distribution.
#'
#' @return A numeric scalar giving the log-likelihood value.
#' @export
#'
#' @examples
#' Likelihood_Rn_Weibull_YNM(R = c(1.5, 2.6, 3.7), L = c(2, 5, 8), T = 10, params = c(0.9, 1.5, 2))
Likelihood_Rn_Weibull_YNM <- function(R,L,T, params) {  ## vector of R, L, value of T, vector of (gamma, 1/rate)
  m=length(R)

  # Calculate Sum1
  s1 <- log(params[1]) * sum(L)

  # Loop to calculate Sum2
  s2a=0
  for (j in 1:m) {
    s2a[j]= log(dweibull(R[j], shape=params[2], scale=params[3]))
  }
  s2=sum(s2a)

  # Loop to calculate Sum3
  s3a=0
  for (j in 1:m) {
    s3a[j]= + ((params[1]^L[j]) - 1) * log(pweibull(R[j], shape=params[2], scale=params[3]))
  }
  s3=sum(s3a)

  # Loop to calculate Sum4
  s4a=0
  for (j in 1:(m - 1)) {
    s4a[j] = (((params[1]^(L[j] + 1)) - (params[1]^L[j + 1])) / (1 - params[1])) *
      log(pweibull(R[j], shape=params[2], scale=params[3]))
  }
  s4=sum(s4a)

  # Calculate Sum5
  s5=0
  if(L[m]<T){
    s5 <- (((params[1]^(L[m] + 1)) - (params[1]^(T + 1))) / (1 - params[1])) *
      log(pweibull(R[m], shape=params[2], scale=params[3]))
  }

  # Sum all parts to calculate `vrais`
  vrais <- s1+s2+s3+s4+s5

  # Return the result
  return(vrais)
}



