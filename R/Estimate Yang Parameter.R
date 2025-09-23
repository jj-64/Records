################### Moments: Estimate gamma parameter ###############

#' Estimate Gamma Parameter Based on Moments and Number of Records
#'
#' Computes the estimator for the gamma parameter using the number of records in the dataset.
#'
#' @param X A numeric vector representing the dataset.
#' @details
#' The gamma parameter is estimated as:
#' \deqn{\hat{\gamma} = \frac{1}{1 - \frac{N}{T}}}
#' where:
#' - \eqn{N} is the number of records in the sequence.
#' - \eqn{T} is the length of the sequence.
#' This is a biased estimator, check #'\code{\link{Estim_gamma_NT_Unbiased}} for the corrected one.
#' @return The estimated gamma parameter.
#' @export
#' @examples
#' X <- c(1.2, 2.3, 3.1, 5.7, 7.9)
#' Estim_gamma_NT(X)
#' [1] 2.083333
Estim_gamma_NT <- function(X) {
  N <- rec_counts(X)
  T <- length(X)
  return(1 / (1 - N / T))
}

#' Variance of Gamma Estimator using number of records
#'
#' Computes the variance of the gamma estimator based on number of records.
#'\code{\link{Estim_gamma_NT}}.
#' @details
#' The variance of \eqn{\hat{\gamma}} is given by:
#' \deqn{\text{Var}(\hat{\gamma}) = \gamma^2 (\gamma - 1)}
#'
#' @param gamma The estimated gamma parameter (\eqn{\gamma \geq 1}).
#' @return The variance of the estimator.
#' @export
#' @examples
#' Estim_gamma_NT_Variance(0.5)
Estim_gamma_NT_Variance = function(gamma){
gamma^2 * (gamma-1)
}


###################### Unbiased based on Moments ##################
### Based on moments, number of records - Bias  #gama_hat2b
#' Unbiased Estimation of Gamma Parameter based on number of records
#'
#' Computes an unbiased estimator of the gamma parameter using a bias correction, correcting
#' \code{\link{Estim_gamma_NT}}
#'
#' @details
#' The unbiased estimator is given by:
#' \deqn{\hat{\gamma}^{(unbiased)} = \hat{\gamma} - \text{Bias}(\hat{\gamma})}
#' where the bias is calculated using:
#' \deqn{\text{Bias}(\hat{\gamma}) = \frac{h_T - h}{(1-h)^2} + \frac{1}{(1-h)^3} \left[ \frac{h_T}{T} - \frac{1}{T^2} \sum_{k=1}^{T} P_T^2 + (h_T - h)^2 \right]}
#'
#' @param X A numeric vector representing the sequence.
#' @return The unbiased estimate of gamma.
#' @export
#' @examples
#' #' Estim_gamma_NT_unbiased(c(3, 1, 4, 1, 5, 9))
Estim_gamma_NT_unbiased = function(X){ ## compute the second estimator*

  HT = function(gamma,T){
    sum(rec_rate_Yang(gamma,1:T))/T
  }

  H = function(gamma){
    1-1/gamma
  }

  ### Bias estimator function
  bias = function(gamma,T) { ## the bias of estimator 2 hat
    hT = HT(gamma,T)
    h = H(gamma)

    a = (hT-h)/(1-h)^2

    b= 1/(1-h)^3

    cc=0
    for(k in 1:T) {
      cc[k] = rec_rate_Yang(gamma,k)^2 + (hT-h)^2
    }

    d=hT/T - (sum(cc))/T^2

    return(a+b*d)
  }

  g2 = Estim_gamma_NT(X) ## estimator 2
  b = bias(gamma=g2,T=length(X))     ## bias
  return(g2-b)
}


####################### MLE delta ##############################
### Likelihood estimator based on indicator series
#' Maximum Likelihood Estimation of Gamma using indicators Series
#'
#' Estimates the gamma parameter using the maximum likelihood method based on indicator series.
#'
#' @details
#' This method finds the gamma parameter by maximizing the log-likelihood function:
#' \deqn{\log ell(v = 1/\gamma) = -N \times log(1-v) + (T-N)log(v) -  log(1-v^T) - \sum_{t=2}^{T} \delta_t log(1-v^{t-1})}
#' where:
#' - \eqn{N} is the number of records.
#' - \eqn{T} is the length of the sequence.
#' - \eqn{\delta_t} represents the indicator series.
#' - \eqn{v} = \eqn{t/\gamma}
#'
#' @param X A numeric vector representing the sequence.
#' @param min Minimum value for search range (default = 1).
#' @param max Maximum value for search range (default = 5).
#' @return Estimated gamma parameter.
#' @export
#' @examples
#' Estim_gamma_indicatorc(3, 1, 4, 1, 5, 9))
Estim_gamma_indicator <- function(X, min = 1, max = 5, step = 0.001) {
  # X: input series
  # min, max: search range for gamma
  # step: resolution of search (smaller -> more precise but slower)

  # --- Precompute useful values ---
  T <- length(X)
  delta <- is_rec(X)        # indicator of records
  N <- rec_counts(X)        # number of records

  # Vectorized version of derivative of log-likelihood wrt v
  dLog <- function(v) {
    k <- 2:T
    x <- delta[k] * (k - 1) * v^(k - 2) / (1 - v^(k - 1))
    -N / (1 - v) + (T - N) / v + (T * v^(T - 1)) / (1 - v^T) + sum(x)
  }

  # Search grid for gamma
  gamma_grid <- seq(min, max, by = step)
  v_grid <- 1 / gamma_grid

  # Evaluate derivative across grid
  dvals <- vapply(v_grid, dLog, numeric(1))

  # Pick gamma where derivative is closest to zero
  Est_gamma <- gamma_grid[which.min(abs(dvals))]

  return(Est_gamma)
}


#' Variance of the Gamma Estimator using Indicator Series
#'
#' Computes the variance for the gamma estimator obtained through the indicator
#' series using \code{\link{Estim_gamma_indicator}}.
#'
#' @details
#' Two methods are available:
#' \itemize{
#'   \item \code{approximate = FALSE} (default): Uses the exact Fisher information formula:
#'     \deqn{I_T(\gamma) = \frac{1}{\gamma^2 (\gamma - 1)^2} \sum_{t=1}^{T} P_t(\gamma)
#'     + \frac{1}{\gamma^2} (T - \sum_{t=1}^{T} P_t(\gamma))
#'     - \frac{T (1+\gamma^T (T-1))}{\gamma^2 (\gamma^T - 1)^2}
#'     - \sum_{i=2}^{T} \frac{(i-1)(1+(i-2)\gamma^{i-1}) P_i(\gamma)}{\gamma^2 (\gamma^{i-1} -1)^2}}
#'
#'   \item \code{approximate = TRUE}: Uses the simplified asymptotic variance:
#'     \deqn{\text{Var}(\hat{\gamma}) \approx \frac{1-v}{T v^3}, \quad v=1/\gamma}
#' }
#'
#' @param T Integer. Length of the sequence.
#' @param gamma Numeric. The estimated gamma parameter.
#' @param approximate Logical. If \code{TRUE}, compute the simplified variance approximation.
#' @return Numeric. Variance of the gamma estimator.
#' @export
#'
#' @examples
#' Estim_gamma_indicator_Variance(T=100, gamma=1.2)
#' [1] 0.003272339
#' Estim_gamma_indicator_Variance(T=100, 1.2, approximate = TRUE)
#' [1] 0.00288
#' Estim_gamma_indicator_Variance(T=300, gamma=1.2)
#' [1] 0.0009999638
#' Estim_gamma_indicator_Variance(T=300, 1.2, approximate = TRUE)
#' [1] 0.00096
Estim_gamma_indicator_Variance <- function(T, gamma, approximate = FALSE) {
  if (approximate) {
    # --- Approximate variance formula ---
    v <- 1 / gamma
    Inv = (1-v)/v^3

    return(Inv/ T)
  }

  # --- Exact Fisher information ---
  ent <- ENT_Yang(T=T, gamma = gamma )

  # Components of Fisher information
  a <- (1 / (gamma^2 * (gamma - 1)^2)) * ent
  b <- (1 / gamma^2) * (T - ent)
  c <- T * (1 + gamma^T * (T - 1)) / (gamma^2 * (gamma^T - 1)^2)

  i <- 2:T
  d <- (i - 1) * (1 + (i - 2) * gamma^(i - 1)) * rec_rate_Yang(gamma, i) /
    (gamma^2 * (gamma^(i - 1) - 1)^2)

  I <- a + b - c - sum(d)

  # Return variance as inverse of Fisher info
  return(1 / I)
}

# Wrapper: estimate + variance together
#' MLE Estimator of Gamma with Variance (indicator series)
#'
#' Computes the MLE of gamma using indicator series, and its asymptotic variance
#' based on Fisher information.
#'
#' @inheritParams Estim_gamma_indicator
#' @return A list with:
#'   \item{gamma_hat}{Estimated gamma parameter.}
#'   \item{var}{Estimated variance of gamma_hat.}
#' @export
#' @examples
#' res <- Estim_gamma_indicator_full(rnorm(30))
#' res$gamma_hat
#' res$var
Estim_gamma_indicator_full <- function(X, min = 1, max = 5, step = 0.001, approximate= FALSE) {
  # Step 1: Estimate gamma
  gamma_hat <- Estim_gamma_indicator(X, min = min, max = max, step = step)

  # Step 2: Variance based on Fisher information
  T <- length(X)
  var_hat <- Estim_gamma_indicator_Variance(T, gamma_hat, approximate = approximate)

  # Output
  return(list(
    estim = gamma_hat,
    var = var_hat
  ))
}
