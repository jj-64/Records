
#' @title Log-Likelihood of First-Passage Times for DTRW
#' @description Computes the log-likelihood of observing a sequence of record times \eqn{L} in a discrete-time random walk (DTRW) framework.
#'
#' @details The function computes the log-likelihood of observing a sequence of record times \eqn{L} given a total time \eqn{T}. The likelihood is computed as:
#' \deqn{ \log \mathcal{L}(L|T) = \sum_{i=1}^{m-1} \log P_F(\Delta L_i) + \log P_S(T - L_m) }
#' where:
#' - \eqn{ \Delta L_i = L_{i+1} - L_i } represents the step differences,
#' - \eqn{ P_F(\cdot) } is the first-passage probability (\code{\link{FirstPass}}), and
#' - \eqn{ P_S(\cdot) } is the survival probability (for cases where \eqn{ L_m < T }). See (\code{\link{Survival}})
#'
#' @param L Vector of record times.
#' @param T Total observation time.
#' @return The log-likelihood value.
#' @export
#' @examples
#' Xt = DTRW_series_Norm(T=25,loc=0, sd=1)
#' Xt
#' [1]  0.0000000  1.4880661 -4.8147394 -5.2196358 -1.0585060  0.1213980 -1.9903581
#' [8] -0.2289111  0.1204328  1.0547639 -0.4275506  0.2533994  2.3335009  3.4109412
#' [15]  3.1307010  5.0952581  5.8665243  8.9529542  7.2892812  4.8499811  0.7474847
#' [22] -0.7717111 -2.3678676 -5.8535865 -4.9697589
#' Likelihood_Ln_DTRW(L = rec_times(Xt), T=25)
#' [1] -11.24239
Likelihood_Ln_DTRW = function(L, T) {
  p = 0
  m = length(L)
  delta = diff(L)  # Compute step differences

  ## First-passage probabilities
  for (i in 1:(m-1)) {
    p[i] = log(FirstPass(delta[i]))
  }

  ## Final survival probability
  q = 0
  if (L[m] < T) {
    q = log(Survival(T - L[m]))
  }

  return(sum(p) + q)
}


# Based on rec_values,rec_times ----------------------------------------------------------

##Log Likelihood
#' @title Log-Likelihood of Record Values in Normal DTRW
#' @description Computes the log-likelihood of observing record values \eqn{R} given record times \eqn{L} in a normal discrete-time random walk (DTRW).
#'
#' @details This function computes:
#' \deqn{ \log \mathcal{L}(R|L, T, \mu, \sigma^2) = \sum_{i=1}^{m-1} \log \phi(R_{i+1} - R_i | (L_{i+1}-L_i)\mu, \sqrt{(L_{i+1}-L_i)\sigma^2}) }
#' where:
#' - \eqn{ \phi(x|\mu, \sigma) } is the normal density function,
#' - \eqn{ \mu } and \eqn{ \sigma^2 } are the mean and variance parameters of the normal step distribution,
#' - \eqn{ R } is the sequence of record values, and
#' - \eqn{L} is the sequence of record times.
#'
#' Additionally, if \eqn{ L_m < T }, the likelihood includes the survival term:
#' \deqn{ \sum_{j=1}^{T-L_m} \log \Phi(0 | j\mu, \sqrt{j\sigma^2}) }
#' where \eqn{ \Phi(\cdot) } is the normal cumulative distribution function (CDF).
#'
#' @param R Vector of record values.
#' @param L Vector of record times.
#' @param T Total observation time.
#' @param params Vector of parameters: variance (\eqn{\sigma^2}) of the step distribution.
#' @return The log-likelihood value.
#' @export
#' @examples
#' Xt = DTRW_series_Norm(T=25,loc=0, sd=1)
#' Xt
#' [1]  0.10608755  2.10580288  4.07636870  2.44521434  1.77945837  1.47440228
#' [7]  1.25472403  3.35634849  5.13609509  3.39129849  3.57836155  3.95685562
#' [13]  2.58880190  3.35365594  3.37435163  3.70094121  3.45657870  2.09124483
#' [19]  0.65113456 -0.00403568 -1.21277245 -0.92442545 -0.71107076 -0.67390482
#' [25] -0.41385780
#' Likelihood_Rn_Norm_DTRW(R=rec_values(Xt), L= rec_times(Xt), T=25, params = c(1))
#' -9.654213
Likelihood_Rn_Norm_DTRW = function(R,L,T,params){ ##variance
  # Example: probability that S_1, S_2, ..., S_5 all < 0
  sigma_cov = function(k, sigma){  ## k should be less than 1000, sigma is variance
    x=matrix(0, nrow=k, ncol=k)
    for(i in 1:k){
      for(j in 1:k){
        x[i,j] = min(i,j)*sigma
      }
    }
    return(x)
  }

  m= length(R)  # Number of observed time points

  if (m < 2) {
    stop("The vector l must contain at least two indices.")
  }

  #mean = (L[i+1]-L[i])*params[1]
  s1b=0
  for(i in 1:(m-1)){
    s1b[i] = log(dnorm(R[i+1]-R[i], mean = 0, sd = sqrt((L[i+1]-L[i])*params) )  )
  }
  s1=sum(s1b)

  ## case where NT<T
  s2=0
  # if(L[m]<T){
  #   for(j in 1:(T-L[m])){
  #     s2 = s2 + log(pnorm(0, mean = (j)*params[1], sd = sqrt(j*params[2]) ) )
  #   }
  # }
  if (L[m] < T) {
    #s2=log(mvtnorm::pmvnorm(lower=-Inf, upper=rep(0,(T - L[m])), sigma = sigma_cov((T - L[m]),params))[1])
      s2 = -0.5*log(pi * (T - L[m]))
    }

  return(s1+s2)

}

# Likelihood_Rn_Norm_DTRW <- function(R,L,T,params) {
#   N_T <- length(L)  # Number of observed time points
#
#   if (N_T < 2) {
#     stop("The vector l must contain at least two indices.")
#   }
#
#   # Compute t = (l_{n+1} - l_n)
#   t_values <- diff(L)
#
#   # Compute the squared differences (x_{l_{n+1}} - x_{l_n})^2
#   increments_sq <- diff(R)^2
#
#   # Compute each term in the log-likelihood function
#   term1 <- - (N_T - 1) / 2 * log(2 * pi * params[2])
#   term2 <- - 0.5 * sum(log(t_values))
#   term3 <- - sum(increments_sq / (2 * t_values * params[2]))
#   term4 <- - (T - L[N_T]) * log(2)
#
#   # Compute the log-likelihood
#   logL <- term1 + term2 + term3 + term4
#
#   return(logL)
# }

## Likelihood L not Log
#' @title Likelihood of Record Values in Normal DTRW
#' @description Computes the likelihood of observing record values \eqn{ R } given record times \eqn{L} in a normal discrete-time random walk (DTRW).
#'
#' @details This function computes the likelihood as:
#' \deqn{ \mathcal{L}(R|L, T, \mu, \sigma^2) = \prod_{i=1}^{m-1} \phi(R_{i+1} - R_i | (L_{i+1}-L_i)\mu, \sqrt{(L_{i+1}-L_i)\sigma^2}) }
#' where:
#' - \eqn{ \phi(x|\mu, \sigma) } is the normal probability density function,
#' - \eqn{ \mu } and \eqn{ \sigma^2 } are the mean and variance of the step distribution,
#' - \eqn{ R } is the sequence of record values, and
#' - \eqn{L} is the sequence of record times.
#'
#' If \eqn{ L_m < T }, the likelihood is multiplied by the survival probability:
#' \deqn{ \prod_{j=1}^{T-L_m} \Phi(0 | j\mu, \sqrt{j\sigma^2}) }
#' where \eqn{ \Phi(\cdot) } is the normal cumulative distribution function (CDF).
#'
#'For the log-Likelihood, refer to \code{\link{Likelihood_Rn_Norm_DTRW}}.
#' @param R Vector of record values.
#' @param L Vector of record times.
#' @param T Total observation time.
#' @param params Vector of parameters: mean (\eqn{\mu}) and variance (\eqn{\sigma^2}) of the step distribution.
#' @return The likelihood value.
#' @examples
#' Xt = DTRW_series_Norm(T=25,loc=0, sd=1)
#' Xt
#' [1]  0.10608755  2.10580288  4.07636870  2.44521434  1.77945837  1.47440228
#' [7]  1.25472403  3.35634849  5.13609509  3.39129849  3.57836155  3.95685562
#' [13]  2.58880190  3.35365594  3.37435163  3.70094121  3.45657870  2.09124483
#' [19]  0.65113456 -0.00403568 -1.21277245 -0.92442545 -0.71107076 -0.67390482
#' [25] -0.41385780
#' Likelihood_Rn_Norm_DTRW(R=rec_values(Xt), L= rec_times(Xt), T=25, params = c(1))
#' -9.654213
#' L_Rn_Norm_DTRW(R=rec_values(Xt), L= rec_times(Xt), T=25, params = c(1))
#'  0.00006416188
#' log(L_Rn_Norm_DTRW (R=rec_values(Xt),L=rec_times(Xt),25,c(1)))
#' -9.654213
L_Rn_Norm_DTRW = function(R, L, T, params) {

  # probability that S_1, S_2, ..., S_5 all < 0
  sigma_cov = function(k, sigma){
    x=matrix(0, nrow=k, ncol=k)
    for(i in 1:k){
      for(j in 1:k){
        x[i,j] = min(i,j)*sigma
      }
    }
    return(x)
  }

  m = length(R)  # Number of observed record values
  s1b = 0

  #mean = (L[i+1] - L[i]) * params[1]
  ## Compute product of normal densities
  for (i in 1:(m-1)) {
    s1b[i] = dnorm(R[i+1] - R[i], mean = 0,
                   sd = sqrt((L[i+1] - L[i]) * params[1]))
  }
  s1 = prod(s1b)

  ## Case where last record time L[m] < T
  s2 = 1
  if (L[m] < T) {
      #s2=mvtnorm::pmvnorm(lower=-Inf, upper=rep(0,(T - L[m])), sigma = sigma_cov((T - L[m]),params[1]))[1]
      s2=1/sqrt(pi * (T - L[m]))
      }


  return(s1 * s2)
}


########### Estimate variance from DTRW Likelihood with 0 mean #########
Likelihood_Rn_Norm_DTRW_DerivativeSolve = function(R,L){
  RR=diff(R)
  t= diff(L)
  m= length(R)

  biased = sum(RR^2/t)/(m-1)
  correcting = sqrt(2/pi)
  const = 0.29795219028
  unbiased = (biased+ const^2 + const) /correcting
  #= biased / (m+1) sqrt(2*T/pi) #sqrt(pi)/2
  #unbiased = biased * pi/2
 return(unbiased)
}
