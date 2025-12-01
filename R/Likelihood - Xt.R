

## convenience wrapper for user
lik <- function(model, obs_type, dist, data, params) {
  logLik_records(model, obs_type, dist, data, params)
}
#lik("LDM","records","gumbel",data = list(rec_values = c(0.1,0.2,0.3,0.4), rec_times = c(1,4,7,8), time = 10), params = list(theta=0.3, loc=1, scale=1))


## Classical -----------------------------
#' Compute Log-Likelihood for IID Frechet Model
#'
#' Computes the log-likelihood of independent and identically distributed (IID)
#' variables based on full information of the series, and under different
#' underlying distribution.
#'
#' @details
#' For \code{Frechet} distribution, use the **VGAM** package to compute:
#' \deqn{llog L = \sum \log(f_{\text{Frechet}}(x; scale shape))}
#'
#' For \code{Normal} distribution, use the **dnorm** function:
#' \deqn{log L = \sum \log(f_{\text{Normal}}(x; \mu, \sigma))}
#' @param x Numeric vector. Observed data.
#' @param dist Character, distribution name. One of:
#'   "norm", "frechet"
#' @param params list of vector of parameter:
#'  \itemize{
#'   \item Frechet: scale and shape.
#'   \item Norm: mean and sd.
#'  }
#' @return Numeric value of the log-likelihood.
#' @export
#' @examples
#' Xt = VGAM::rfrechet(n = 100, shape = 1, scale = 2)
#' Likelihood_Xt_iid(x = Xt, dist = "frechet", params =  list(shape = 1, scale=2))
#' Likelihood_Xt_iid(x = Xt, dist = "norm", params =  list(mean = mean(X), sd= sd(X) ))
Likelihood_Xt_iid <- function(x, dist = c("frechet","norm"), params) {
  dist <- match.arg(dist)

  if(dist == "frechet"){
  return(sum(log(VGAM::dfrechet(x = x, 0, scale = params$scale, shape =params$shape))))
  }
  if(dist == "norm"){
    return(sum(log(dnorm(x = x, mean = params$mean, sd = params$sd ))))
  }
}

##  Likelihood DTRW Xt--------------------

#' Log-Likelihood expression
#'
#'Compute Log-Likelihood for DTRW Model with Gaussian, Cauchy or uniform noise
#' @param x Numeric vector. Observed data.
#' @param dist Character, distribution name. One of:
#'   "norm", "cauchy", "uniform"
#' @param params list of vector of parameter:
#'  \itemize{
#'   \item norm: mean and sd
#'   \item cauchy: loc and scale.
#'   \item unifrom: min and max
#'  }
#' @returns numerical value
#' @export
#' @examples
#' Xt = DTRW_series(T=25,dist = "norm", mean=0, sd=3)
#' Likelihood_Xt_DTRW(x=Xt, dist = "norm", params=list(mean =0, sd= 3))
#' Likelihood_Xt_DTRW(x=Xt, dist = "cauchy", params = list (loc =0, scale = 3))
#' Likelihood_Xt_DTRW(x=Xt, dist = "uniform", params = list (min =-1, max = 1))
Likelihood_Xt_DTRW =  function(x, dist = c("norm", "cauchy", "uniform"), params){
  if (length(x) < 2) {
    stop("The time series must have at least two observations.")
  }
  dist <- match.arg(dist)
  if(dist == "norm"){
 return(sum(log(dnorm(x=diff(x),mean=params$mean,sd = params$sd ))))
  }
  if(dist == "cauchy"){
    return(sum(log(dcauchy(x=diff(x),location=params$loc, scale = params$scale ))))
  }
  if(dist == "uniform"){
    return(sum(log(dunif(x=diff(x),min=params$min, max = params$max ))))
  }
}

## Likelihood Functions for YNM Process ------------------

#' Compute Log-Likelihood for Frechet YNM Process
#'
#' Computes the log-likelihood for a series following a Frechet-distributed YNM process, based on full information.
#'
#' @details
#' The likelihood function is given by:
#' \deqn{log L(x | ( \gamma, scale = A, shape =a )) = T \log(a A^{-a}) - (a+1) \sum \log(x) + \frac{T(T+1)}{2} \log(\gamma) - \sum \gamma^t (A x)^{-a}}
#'
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
#' Xt = YNM_series(T=25, dist = "frechet", gamma=1.1, shape=1, scale=2)
#' Likelihood_Xt_YNM(x=Xt, dist = "frechet", params = list(gamma = 1.1, shape= 1, scale=2))
Likelihood_Xt_YNM <- function (x, dist = "frechet", params) {
  dist <- match.arg(dist)

  n= length(x)
  if(dist == "frechet"){
  s1 <- n* log(params$scale * params$shape^(-params$scale) )
  s2 <- -(params$scale + 1) * sum(log(x))
  s3 <- (n* (n+ 1) / 2) * log(params$gamma)
  s4 <- -sum(params$gamma^(1:T) * (params$shape * x)^(-params$scale))
  Lik = s1 + s2 + s3 + s4
  }
  return(Lik)
}

## Likelihood Functions for LDM  Process --------------------
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
#' Xt = LDM_series(T=25, dist = "frechet", theta=0.2, shape=1, scale=2)
#' Likelihood_Xt_LDM(x=Xt, dist = "frechet", params = list(theta = 0.2, shape = 1/2, scale= 1) )
Likelihood_Xt_LDM <- function (x, dist = "frechet", params) {

  dist <- match.arg(dist)

  n = length(x)

  X <- x - params$theta * (1:n)

  if (any(X <= 0)) return(-1000)

  if(dist == "frechet"){
  A <- (sum(X^-params$scale) / n)^(1/params$scale)

  s1 <- n * log(params$scale * A^(-params$scale))
  s2 <- -(params$scale + 1) * sum(log(X))
  s3 <- -sum((A * X)^(-params$scale))
  Lik = s1 + s2 + s3
  }
  return(Lik)
}

## Functions to Maximize Log-Likelihood LDM -------------------

#' Estimate Scale Parameter in Frechet LDM Process
#'
#'Equation to solve the 1/Scale parameter in the case of Frechet distribution in LDM process.
#'Using the first derivative equations, \deqn{A = (\frac{\sum X^{-a}}{T})^{1/a}}
#' @param shape Shape parameter.
#' @param x Numeric vector of observations.
#' @param theta Drift parameter.
#' @return Numeric value of scale parameter A.
eq_A <- function(shape, x, theta) {
  n <- length(x)
  X <- x - theta * (1:n)
  return((sum(X^-shape) / n)^(1/shape))
}

#' Solve for Drift Parameter Theta (Left-Side Equation) in Frechet LDM process
#'
#' @inheritParams eq_A
#' @details
#' \eqn{Y_t = X_t - \theta t}
#' \deqn{\sum t \times Y_t^{-a} / \sum t \times Y_t^{-1}  }
#'
#' @return Numeric value for the equation.
eq_theta_Left <- function(shape, x, theta) {
  n <- length(x)
  X <- x - theta * (1:n)
  return(sum((1:n) * X^-shape) / sum((1:n) * X^-1))
}

#' Solve for Drift Parameter Theta (Right-Side Equation) in Frechet LDM process
#'
#' @inheritParams eq_A
#' @return Numeric value for the equation.
eq_theta_Right <- function(shape, x, theta) {
  n <- length(x)
  X <- x - theta * (1:n)
  return((shape + 1) * sum(X^-shape) / (n * shape))
}

#' Solve for Shape Parameter Alpha (Left-Side Equation) in Frechet LDM process
#'
#' @inheritParams eq_A
#' @return Numeric value for the equation.
eq_alpha_Left <- function(shape, x, theta) {
  n <- length(x)
  X <- x - theta * (1:n)
  A <- eq_A(shape, x, theta)
  return(sum((A * X)^(-shape) * log(A * X)) + n * (1 - log(A)) / shape)
}

#' Solve for Shape Parameter Alpha (Right-Side Equation) in Frechet LDM process
#'
#' @inheritParams eq_A
#' @return Numeric value for the equation.
eq_alpha_Right <- function(x, theta) {
  n <- length(x)
  return(sum(log(x - theta * (1:n))))
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
#' x = LDM_series(T=25,dist = "frechet", theta=0.5,scale=1,shape=2)
#' x
#'  [1]  1.608046  1.549613  2.788720  4.077933  3.886494  6.142666  4.709011
#'  [8]  5.274752  5.512561  7.164634  6.340624  8.823882  7.447867  7.925539
#'  [15]  8.550669  9.028977  9.253112 10.572839  9.988717 11.996890 12.049306
#'  [22] 12.226410 11.992932 12.747571 13.098460
#' Likelihood_Xt_Frechet_LDM_DerivativeSolve(x, c(1,0.2))
#'      theta         A         shape
#'   0.4383729 0.5824738 1.8204038
Likelihood_Xt_Frechet_LDM_DerivativeSolve <- function(x, initial_guess) {  #Solve_Log
  system_of_eqs <- function(vars) {
    shape <- vars[1]
    theta <- vars[2]

    X <- x - theta * (1:length(x))
    A <- eq_A(shape, x, theta)

    f1 <- eq_theta_Left(shape, x, theta) - eq_theta_Right(shape, x, theta)
    f2 <- eq_alpha_Left(shape, x, theta) - eq_alpha_Right(x, theta)

    return(c(f1, f2))
  }

  solution <- nleqslv(initial_guess, system_of_eqs)

  A <- eq_A(shape = solution$x[1], x = x, theta = solution$x[2])

  return(c(theta=solution$x[2], A=A, shape=solution$x[1]))
}


#' Solve Log-Likelihood for LDM Model with Frechet distribution
#'
#' Uses **nlminb** for numerical optimization of log-likelihood expression \code{\link{Likelihood_Xt_Frechet_LDM}} to estimate parameters
#'
#' @inheritParams eq_A
#' @return Estimated parameters \eqn{(\theta, A, a, L_{\max})}.
#' @export
#' @examples
#' x = LDM_series(T=25,dist = "frechet", theta=0.5,scale=1,shape=2)
#' x
#'  [1]  1.608046  1.549613  2.788720  4.077933  3.886494  6.142666  4.709011
#'  [8]  5.274752  5.512561  7.164634  6.340624  8.823882  7.447867  7.925539
#'  [15]  8.550669  9.028977  9.253112 10.572839  9.988717 11.996890 12.049306
#'  [22] 12.226410 11.992932 12.747571 13.098460
#' Likelihood_Xt_Frechet_LDM_OptSolve(x)
#'theta           A           shape        logL
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

  return(c(theta=result$par[1], A=A, shape=result$par[3], logL=-result$objective))
}

#library(nleqslv)


