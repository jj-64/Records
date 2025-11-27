###########################Classical #########################
#' Compute Log-Likelihood for IID Frechet Model
#'
#' Computes the log-likelihood of independent and identically distributed (IID) Frechet variables based on full information of the series.
#'
#' @details
#' Uses the **VGAM** package to compute:
#' \deqn{llog L = \sum \log(f_{\text{Frechet}}(x; scale=1/A, shape=a))}
#' Note that the growth parameter \eqn{\gamma} is fixed to 1 in the IID case.
#'
#' @param T Integer. Length of the series.
#' @param x Numeric vector. Observed data.
#' @param params Numeric vector of parameters \eqn{( A, a )}:
#'  \itemize{
#'   \item \eqn{A}: 1/Scale parameter (1/scale).
#'   \item \eqn{a}: Shape parameter.
#'  }
#' @return Numeric value of the log-likelihood.
#' @export
Likelihood_Xt_Frechet_iid <- function(T, x, params) {
  sum(log(VGAM::dfrechet(x = x, 0, scale = 1 / params[1], shape = params[2])))

  # s1=T*log(params[1]*params[2])
  # s2 = -(params[2]+1)*sum(log(params[1]*x))
  # s3 = -sum((params[1]*x)^(-params[2]))
  # s1+s2+s3
}

#' Compute Log-Likelihood for IID Normal Model
#'
#' Computes the log-likelihood of independent and identically distributed (IID) normal variables.
#'
#' @details
#' Uses the **dnorm** function:
#' \deqn{log L = \sum \log(f_{\text{Normal}}(x; \mu, \sigma))}
#'Note that the growth parameter \eqn{\gamma} is fixed to 1 in the IID case.
#' @inheritParams Likelihood_Xt_Frechet_YNM
#' @return Numeric value of the log-likelihood.
#' @export
Likelihood_Xt_Norm_iid <- function(T, x, params) {
  sum(log(dnorm(x = x, mean = params[1], sd = params[2])))

  # s1=T*log(params[1]*params[2])
  # s2 = -(params[2]+1)*sum(log(params[1]*x))
  # s3 = -sum((params[1]*x)^(-params[2]))
  # s1+s2+s3
}

################### Likelihood DTRW Xt #######################

#' Log-Likelihood expression
#'
#'Compute Log-Likelihood for DTRW Model with Gaussian Noise
#' @param x Numeric vector. Observed data.
#' @param params Numeric vector of parameters \eqn{( A, a )}:
#'  \itemize{
#'   \item \eqn{A}: Mean of normal distribution assumed 0 for no drift models.
#'   \item \eqn{a}: variance of normal distribution
#'  }
#'
#' @returns numerical value
#' @export
#' @examples
#' Xt=DTRW_series_Norm(T=25,loc=0, sd=3)
#'> Xt
#'[1]  0.0000000  1.4880661 -4.8147394 -5.2196358 -1.0585060  0.1213980 -1.9903581
#'[8] -0.2289111  0.1204328  1.0547639 -0.4275506  0.2533994  2.3335009  3.4109412
#'[15]  3.1307010  5.0952581  5.8665243  8.9529542  7.2892812  4.8499811  0.7474847
#'[22] -0.7717111 -2.3678676 -5.8535865 -4.9697589
#' Likelihood_Xt_Norm_DTRW(x=Xt, params=c(0,3))
#' -61.32352
#' Likelihood_Xt_Norm_DTRW(x=Xt, params=c(0,1))
#' -89.51009
Likelihood_Xt_Norm_DTRW =  function(x,params){


  if (length(x) < 2) {
    stop("The time series must have at least two observations.")
  }
  logL=sum(log(dnorm(x=diff(x),mean=params[1],sd = sqrt(params[2]))))
  #T <- length(x) - 1  # Number of steps in the random walk
  #increments <- diff(x)  # Compute X_t - X_{t-1}

  #logL <- - (T / 2) * log(2 * pi * params[2]) - sum((increments-params[1])^2) / (2 * params[2])

  return(logL)

}

############## Likelihood Functions for YNM Process ##############

#' Compute Log-Likelihood for Frechet YNM Process
#'
#' Computes the log-likelihood for a series following a Frechet-distributed YNM process, based on full information.
#'
#' @details
#' The likelihood function is given by:
#' \deqn{log L(x | ( \gamma, A, a )) = T \log(a A^{-a}) - (a+1) \sum \log(x) + \frac{T(T+1)}{2} \log(\gamma) - \sum \gamma^t (A x)^{-a}}
#'
#' @param T Integer. Length of the series.
#' @param x Numeric vector. Observed data.
#' @param params Numeric vector of parameters \eqn{( \gamma, A, a )}:
#'  \itemize{
#'   \item \eqn{\gamma}: Growth parameter.
#'   \item \eqn{A}: 1/Scale parameter (1/scale).
#'   \item \eqn{a}: Shape parameter.
#'  }
#' @return Numeric value of the log-likelihood.
#' @export
#' @examples
#' Xt = YNM_series_Frechet(T=25, gamma=1.1, shape=1, scale=2)
#' Likelihood_Xt_Frechet_YNM(T=25, x=Xt, params = c(1.1,1/2,1))
#'
Likelihood_Xt_Frechet_YNM <- function(T, x, params) {
  s1 <- T * log(params[3] * params[2]^(-params[3]))
  s2 <- -(params[3] + 1) * sum(log(x))
  s3 <- (T * (T + 1) / 2) * log(params[1])
  s4 <- -sum(params[1]^(1:T) * (params[2] * x)^(-params[3]))
  return(s1 + s2 + s3 + s4)
}


############## Likelihood Functions for LDM  Process ##############
#' Compute Log-Likelihood for Frechet LDM Model
#'
#' Computes the log-likelihood of a series using a **Frechet-distributed** Linear Drift Model (LDM), based on full information
#'
#' @details
#' The transformation applied is:
#' \deqn{Y_t = X_t - \theta t}
#' where \eqn{Y_t > 0} since we are assuming it is coming from a positive Frechet distribution. If any \eqn{Y_t \leq 0}, the function returns a penalty.
#'
#' The likelihood function is:
#' \deqn{log L(X | ( \theta, A, a )) = T \log(a A^{-a}) - (a+1) \sum \log(X) - \sum (A X)^{-a}}
#'
#' @param T Integer. Length of the series.
#' @param x Numeric vector. Observed data.
#' @param params Numeric vector of parameters \eqn{( \theta, A, a )}:
#'  \itemize{
#'   \item \eqn{\theta}: slope parameter.
#'   \item \eqn{A}: 1/Scale parameter (1/scale).
#'   \item \eqn{a}: Shape parameter.
#'  }
#' @return Numeric value of the log-likelihood.
#' @export
#' @examples
#' Xt = LDM_series_Frechet(T=25, theta=0.2, shape=1, scale=2)
#' Likelihood_Xt_Frechet_LDM(T=25, x=Xt, params = c(0.2,1/2, 1) )
Likelihood_Xt_Frechet_LDM <- function(T, x, params) {
  X <- x - params[1] * (1:T)

  if (any(X <= 0)) return(-1000)

  A <- (sum(X^-params[3]) / T)^(1/params[3])

  s1 <- T * log(params[3] * A^(-params[3]))
  s2 <- -(params[3] + 1) * sum(log(X))
  s3 <- -sum((A * X)^(-params[3]))

  return(s1 + s2 + s3)
}





############## Functions to Maximize Log-Likelihood LDM ##############

#' Estimate Scale Parameter A in Frechet LDM Process
#'
#'Equation to solve the 1/Scale parameter in the case of Frechet distribution in LDM process.
#'Using the first derivative equations, \deqn{A = (\frac{\sum X^{-a}}{T})^{1/a}}
#' @param a Shape parameter.
#' @param x Numeric vector of observations.
#' @param theta Drift parameter.
#' @return Numeric value of scale parameter A.
eq_A <- function(a, x, theta) {
  T <- length(x)
  X <- x - theta * (1:T)
  return((sum(X^-a) / T)^(1/a))
}

#' Solve for Drift Parameter Theta (Left-Side Equation) in Frechet LDM process
#'
#' @inheritParams eq_A
#' @details
#' \eqn{Y_t = X_t - \theta t}
#' \deqn{\sum t \times Y_t^{-a} / \sum t \times Y_t^{-1}  }
#'
#' @return Numeric value for the equation.
eq_theta_Left <- function(a, x, theta) {
  T <- length(x)
  X <- x - theta * (1:T)
  return(sum((1:T) * X^-a) / sum((1:T) * X^-1))
}

#' Solve for Drift Parameter Theta (Right-Side Equation) in Frechet LDM process
#'
#' @inheritParams eq_A
#' @return Numeric value for the equation.
eq_theta_Right <- function(a, x, theta) {
  T <- length(x)
  X <- x - theta * (1:T)
  return((a + 1) * sum(X^-a) / (T * a))
}

#' Solve for Shape Parameter Alpha (Left-Side Equation) in Frechet LDM process
#'
#' @inheritParams eq_A
#' @return Numeric value for the equation.
eq_alpha_Left <- function(a, x, theta) {
  T <- length(x)
  X <- x - theta * (1:T)
  A <- eq_A(a, x, theta)
  return(sum((A * X)^(-a) * log(A * X)) + T * (1 - log(A)) / a)
}

#' Solve for Shape Parameter Alpha (Right-Side Equation) in Frechet LDM process
#'
#' @inheritParams eq_A
#' @return Numeric value for the equation.
eq_alpha_Right <- function(x, theta) {
  T <- length(x)
  return(sum(log(x - theta * (1:T))))
}





#' Solve the System of Equations for LDM
#'
#' Estimate the parameters of Log-Likelihood of LDM process with Frechet distribution and based on full dataset
#' \code{\link{Likelihood_Xt_Frechet_LDM}} using the first derivatives equations
#' @details Uses **nleqslv** to solve for parameters \eqn{(a, \theta)}.
#'
#' @param x Numeric vector of observations.
#' @param initial_guess Numeric vector. Initial guess for \eqn{(\theta, a)}.
#' @return Numeric vector of estimated parameters in order \eqn{\theta, A, a}
#' @export
#' @examples
#' x = LDM_series_Frechet(T=25,theta=0.5,scale=1,shape=2)
#' x
#'  [1]  1.608046  1.549613  2.788720  4.077933  3.886494  6.142666  4.709011
#'  [8]  5.274752  5.512561  7.164634  6.340624  8.823882  7.447867  7.925539
#'  [15]  8.550669  9.028977  9.253112 10.572839  9.988717 11.996890 12.049306
#'  [22] 12.226410 11.992932 12.747571 13.098460
#' Likelihood_Xt_Frechet_LDM_DerivativeSolve(x, c(1,0.2))
#'      theta         A         a
#'   0.4383729 0.5824738 1.8204038
Likelihood_Xt_Frechet_LDM_DerivativeSolve <- function(x, initial_guess) {  #Solve_Log
  system_of_eqs <- function(vars) {
    a <- vars[1]
    theta <- vars[2]

    X <- x - theta * (1:length(x))
    A <- eq_A(a, x, theta)

    f1 <- eq_theta_Left(a, x, theta) - eq_theta_Right(a, x, theta)
    f2 <- eq_alpha_Left(a, x, theta) - eq_alpha_Right(x, theta)

    return(c(f1, f2))
  }

  solution <- nleqslv(initial_guess, system_of_eqs)

  A <- eq_A(a = solution$x[1], x = x, theta = solution$x[2])

  return(c(theta=solution$x[2], A=A, a=solution$x[1]))
}


#' Solve Log-Likelihood for LDM Model with Frechet distribution
#'
#' Uses **nlminb** for numerical optimization of log-likelihood expression \code{\link{Likelihood_Xt_Frechet_LDM}} to estimate parameters
#'
#' @inheritParams eq_A
#' @return Estimated parameters \eqn{(\theta, A, a, L_{\max})}.
#' @export
#' @examples
#' x = LDM_series_Frechet(T=25,theta=0.5,scale=1,shape=2)
#' x
#'  [1]  1.608046  1.549613  2.788720  4.077933  3.886494  6.142666  4.709011
#'  [8]  5.274752  5.512561  7.164634  6.340624  8.823882  7.447867  7.925539
#'  [15]  8.550669  9.028977  9.253112 10.572839  9.988717 11.996890 12.049306
#'  [22] 12.226410 11.992932 12.747571 13.098460
#' Likelihood_Xt_Frechet_LDM_OptSolve(x)
#'theta           A           a        logL
#'0.4909235   0.9644116   2.5826999 -20.9373215
Likelihood_Xt_Frechet_LDM_OptSolve <- function(x) {
  T=length(x)
  objective_wrapper <- function(params) {
    -Likelihood_Xt_Frechet_LDM(T, x, params)
  }

  result <- nlminb(
    start = c(0.01, 0.01, 0.01),
    objective = objective_wrapper,
    lower = c(0.01, 0.01, 0.01)
  )

  X <- x - result$par[1] * (1:T)
  A <- (sum(X^-result$par[3]) / T)^(1/result$par[3])

  return(c(theta=result$par[1], A=A, a=result$par[3], logL=-result$objective))
}

#library(nleqslv)


