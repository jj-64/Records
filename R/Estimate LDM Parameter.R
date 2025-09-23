
############################### Based on moments, number of records #####################
#' Estimation of \eqn{\theta} in LDM Process Based on Number of Records
#'
#' Estimates the parameter \eqn{\theta} in the Linear Drift Model (LDM) process using
#' the number of records.
#'
#' @details
#' An LDM process is defined as:
#' \deqn{ X_t = Y_t + \theta t, }
#' where \eqn{Y_t} are independent and identically distributed (i.i.d.) random variables.
#' The properties are studied when \eqn{Y_t} follows a \eqn{Gumbel(\alpha, \beta)} distribution.
#'
#' The estimation formula is given by:
#' \deqn{ \hat{\theta} = -\log\left(1 - \frac{N_T}{T} \right) }
#' where \eqn{N_T} is the number of records observed in the process \eqn{X_t}, and \eqn{T} is the length of the process.
#'
#' This estimator is known to be biased, and a bias-corrected version is provided in
#' \code{\link{Estim_theta_NT_unbiased}}.
#'
#' @param X A numeric vector representing the time series data.
#'
#' @return A numeric estimate of the parameter \eqn{\theta}.
#'
#' @examples
#' Yt <- rnorm(25)
#' Xt <- Yt + 0.2 * (1:25)
#' rec_counts(Xt)  # Number of records
#' Estim_theta_NT(X = Xt)
#'
#' @export
Estim_theta_NT <- function(X) {
  -log(1 - rec_counts(X) / length(X))
}

#' Variance of the \eqn{\theta} Estimator Based on Number of Records
#'
#' Computes the variance of the estimator \code{\link{Estim_theta_NT}}.
#'
#' @details
#' The variance formula is given by:
#' \deqn{ V(\hat{\theta}) = \frac{1 - e^{-\theta / \text{scale}}}{e^{-\theta / \text{scale}}} }
#'
#' @param theta Numeric. The estimated parameter \eqn{\theta}.
#' @param scale Numeric. Default is 1. The scale parameter of the \eqn{Gumbel} underlying distribution.
#'
#' @return The variance of the estimator.
#'
#' @examples
#' Estim_theta_NT_Variance(theta = 0.5)
#'
#' @export
Estim_theta_NT_Variance = function(theta,scale=1){
  (1-exp(-theta/scale))/exp(-theta/scale)
}


################################## Unbiased Moments ###############################
### Based on moments, number of records - Bias  #theta_estm2B
#' Bias-Corrected Estimation of \eqn{\theta} in LDM Process based on number of records
#'
#' Corrects the bias in \code{\link{Estim_theta_NT}} to provide a more accurate estimate.
#'
#' @details
#' The correction is computed based on the expected bias of the estimator.
#'
#' @param X A numeric vector representing the time series data.
#' @param scale Numeric. Default is 1. Scale parameter for the \eqn{Gumbel} underlying distribution.
#'
#' @return A bias-corrected estimate of \eqn{\theta}.
#'
#' @examples
#' Yt <- rnorm(25)
#' Xt <- Yt + 0.2 * (1:25)
#' Estim_theta_NT_unbiased(X = Xt)
#'
#' @export
Estim_theta_NT_unbiased= function(X,scale=1){

  GT=function(theta,T,scale=1){
    sum(rec_rate_LDM(theta,1:T,scale))/T}

  G = function(theta,scale=1) {1-exp(-theta/scale)}

  u2=function(t,theta){(1-exp(-theta))^2/(1-exp(-theta*t))^2}

  bias = function(theta,T,scale=1){

    gT = GT(theta,T,scale)
    g = G(theta,scale)

    a=(gT-g)/(1-g)

    b= 1/(2*(1-g)^2)

    cc=0
    for(k in 1:T) {
      cc[k] = rec_rate_LDM(theta,k,scale)^2 + (gT-g)^2
    }

    d= gT/T - (sum(cc))/T^2

    return(a+b*d)
  }

  theta= Estim_theta_NT(X)
  theta_bias = theta - bias(theta=theta,T=length(X),scale)
  return(theta_bias)
}


#' Variance of the Unbiased Estimator \eqn{\vartheta(\theta)}
#'
#' Computes the variance \eqn{\vartheta(\theta)} for the unbiased estimator of \eqn{\theta} in the LDM process.
#' Variance of the Unbiased Estimator \eqn{\vartheta(\theta)}
#'
#' Computes the variance \eqn{\vartheta(\theta)} for the unbiased estimator of \eqn{\theta} in the LDM process.
#'
#' @details
#' The variance formula is given by:
#'
#' \deqn{\vartheta(\theta) = \lambda(\theta) \left( \frac{d H_T(\theta)}{d \theta} \right)^2}
#'
#' where:
#'
#' \deqn{\frac{d H_T(\theta)}{d \theta} = 1 - \frac{d G_T(\theta)}{d \theta}}
#'
#' \deqn{= 1 - \frac{\frac{d P(\theta)}{d \theta}}{(1 - P(\theta))^2} \left[ \frac{1}{T} \sum_{t=1}^{T} P_t(\theta) - P(\theta) \right]}
#'
#' \deqn{+ \frac{1}{(1 - P(\theta))} \left[ \frac{1}{T} \sum_{t=1}^{T} \frac{d P_t(\theta)}{d \theta} - \frac{d P(\theta)}{d \theta} \right]}
#'
#' \deqn{+ \frac{\frac{d P(\theta)}{d \theta}}{(1 - P(\theta))^3} \left[ \frac{1}{T^2} \sum_{t=1}^{T} P_t(\theta) (1 - P_t(\theta)) + \left( \frac{1}{T} \sum_{t=1}^{T} P_t(\theta) - P(\theta) \right)^2 \right]}
#'
#' \deqn{+ \frac{1}{2 (1 - P(\theta))^2} \left[ \frac{1}{T^2} \sum_{t=1}^{T} \frac{d P_t(\theta)}{d \theta} - \frac{1}{T^2} \sum_{t=1}^{T} \frac{d^2 P_t(\theta)}{d \theta^2} \right]}
#'
#' \deqn{+ 2 \left( \frac{1}{T} \sum_{t=1}^{T} P_t(\theta) - P(\theta) \right) \left( \frac{1}{T} \sum_{t=1}^{T} \frac{d P_t(\theta)}{d \theta} - \frac{d P(\theta)}{d \theta} \right)}
#'
#'where \eqn{P(\theta) = 1-e^{-\theta}}, \eqn{P_t(\theta) = \frac{1-e^{-\theta}}{1-e^{-\theta t}}}, and \eqn{\lambda(\theta)} is the variance of the estimated paramter with bias in \code{\link{Estim_theta_NT_Variance}}
#' @param theta The estimated parameter \eqn{\theta}.
#' @param T Length of the time series.
#' @param scale Scaling parameter (default = 1).
#'
#' @return The computed variance \eqn{\vartheta(\theta)}.
#'
#' @export
#' @examples
#' Estim_theta_NT_Unbiased_Variance(theta = 0.5, T = 25)
Estim_theta_NT_Unbiased_Variance <- function(T, theta, scale = 1) {

  # Compute P(theta) and P_t(theta)
  P_theta <- 1 - exp(-theta)
  dP_dtheta = exp(-theta)
  P_t_theta <- function(t, theta) (1 - exp(-theta)) / (1 - exp(-t * theta))
  d_P_t_theta <- function(t, theta) ((exp(t * theta) - t * exp(theta) + t - 1) * exp(t * theta - theta)) / (exp(t * theta) - 1)^2
  d_P2_t_theta <- function(t, theta) -((  (exp(2 * t * theta)) + (-t^2 * exp(theta) + t^2 + 2 * t - 2) * exp((t * theta)) - t^2 * exp(theta) + t^2 - 2 * t + 1) * exp((t * theta - theta))) / (exp((t * theta)) - 1)^3
  Sum_P_t_theta = 1/T * sum(P_t_theta(t=1:T, theta=theta))
  Sum_d_P_t_theta = 1/T * sum(d_P_t_theta(t=1:T, theta=theta))
  Sum_d_P2_t_theta = 1/T * sum(d_P2_t_theta(t=1:T, theta=theta))
  # Compute first derivative dP/dθ
  dP_dtheta <- exp(-theta)

  # Compute dH_T/dθ
  dH_T_dtheta <- 1 - (dP_dtheta / (1 - P_theta)^2) *( Sum_P_t_theta - P_theta)

  dH_T_dtheta <- dH_T_dtheta + (1 / (1 - P_theta)) *(Sum_d_P_t_theta - dP_dtheta)

  dH_T_dtheta <- dH_T_dtheta + (dP_dtheta / (1 - P_theta)^3) * (
    (1/T^2) * sum(P_t_theta(1:T, theta)*(1-P_t_theta(1:T, theta)))
    + (Sum_P_t_theta - P_theta)^2)

  dH_T_dtheta <- dH_T_dtheta + (1/(2*(1-P_theta)^2)) * (
    Sum_d_P_t_theta/T - Sum_d_P2_t_theta/T)

  dH_T_dtheta <- dH_T_dtheta + 2*(Sum_P_t_theta - P_theta)*(Sum_d_P_t_theta- dP_dtheta)

  # Define lambda(theta), here assumed as 1 for simplicity
  lambda_theta <- Estim_theta_NT_Variance(theta=theta, scale=scale)

  # Compute variance
  v_theta <- lambda_theta * (dH_T_dtheta)^2

  return(v_theta)
}


######################## MLE ########################################
### Likelihood estimator based on delta series  #theta_estm3
#' Maximum Likelihood Estimation of \eqn{\theta} Based on Delta Series
#'
#' Computes the maximum likelihood estimator (MLE) for \eqn{\theta} based on
#' the transformed delta series.
#'
#' @details
#' The estimator is obtained by maximizing the log-likelihood function:
#' \deqn{ \log L(\theta) = N \log(1 - z) + (T - N) \log z - \log(1 - z^T) - \sum \delta_t \log(1 - z^{t-1}) }
#' where \eqn{z = e^{-\theta}}, \eqn{N} is the number of records, \eqn{T} is the series length,
#' and \eqn{\delta_t} are delta-transformed values.
#'
#' @param X A numeric vector representing the time series data.
#' @param min Lower bound for \eqn{\theta} search space. Default is 0.0001.
#' @param max Upper bound for \eqn{\theta} search space. Default is 5.
#'
#' @return A numeric estimate of \eqn{\theta} obtained via maximum likelihood.
#'
#' @examples
#' Yt <- rnorm(25)
#' Xt <- Yt + 0.2 * (1:25)
#' Estim_theta_indicatorX = Xt)
#'
#' @export
Estim_theta_indicator = function(X,min=0.0001,max=5) {
  LogT_LDM = function(z,T,delta,N){
    v = function(t,z){ log(1-z^(t-1)) }
    x=0; for(i in 2:T) x[i] = v(t=i,z=z)
    N*log(1-z) + (T-N)*log(z) - log(1-z^T)-sum(delta*x)
  }

  d=is_rec(X)
  NT=rec_counts(X)
  T= length(X)

  z1 = seq(min,max,by=0.01)
  z=exp(-z1)
  max_likel = 0
  for(i in 1:length(z)) {
    max_likel[i] = LogT_LDM(z=z[i],T=T,delta=d,N=NT)}
  theta_hat3 = -log(z[which.max(max_likel)])
  return(theta_hat3)
}



#' Variance of Maximum Likelihood Estimator of \eqn{\theta}
#'
#' Computes the variance of the MLE estimator for \eqn{\theta}.
#'
#' @param T Integer. The length of the time series.
#' @param theta Numeric. The estimated parameter \eqn{\theta}.
#' @param scale Numeric. Default is 1. Scale parameter for the \eqn{Gumbel} underlying distribution.
#'
#' @return The variance of the estimator.
#'
#' @examples
#'Estim_theta_indicator_Variance(T = 25, theta = 0.5)
#'[1] 0.0316963
#'
#' @export
Estim_theta_indicator_Variance= function(T,theta,scale=1){  #V_theta_estm3
  theta=theta/scale

  a=exp(-2*theta)/(1-exp(-theta))^2
  b=sum(rec_rate_LDM(theta,1:T,scale=1))

  cc=T*exp(-T*theta)*(T+exp(-T*theta)-1)/(1-exp(-T*theta))^2
  d=function(t,theta){
    (t-1)*exp(-theta*(t-1))*(t-2+exp(-theta*(t-1)))*rec_rate_LDM(theta,t,scale=1)/(1-exp(-theta*(t-1)))^2
  }
  dd=sum(d(t=2:T,theta=theta))
  return((a*b+T-b-cc-dd)^-1)
}

