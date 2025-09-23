################## Generic function ###############
#' Generate a DTRW Process
#'
#' Simulates a discrete-time random Walk process of length T under different underlying distributions.
#'
#' For DTRW without a drift, the distribution of the errors shall be continuous and symmetric.
#'
#' @param T Integer, length of the series.
#' @param dist Character, distribution name. One of:
#'  "norm", "cauchy, "uniform"
#' @param ... Additional parameters specific to the chosen distribution:
#'   \describe{
#'     \item{norm}{`loc`, `sd`}
#'     \item{cauchy}{`loc`, `scale`}
#'     \item{uniform}{`min`, `max`}
#'   }
#'
#' @return A numeric vector of length T, the simulated DTRW process.
#' @examples
#' DTRW_series(10,  dist = "cauchy", loc = 0, scale = 1)
#' # [1] -0.6905644  2.1214308  2.7874249  4.1135190  3.3054300  2.6729198  2.5556620  1.6442579 15.6275793 15.4525462
#'
#' DTRW_series(100,  dist = "uniform", min = -1, max = 1)
#' DTRW_series(100,  dist = "norm", loc = 0, sd = 1)
#' @export
DTRW_series <- function(T, dist = c("norm", "cauchy", "uniform"), ...) {
  dist <- match.arg(dist)
  args <- list(...)
  X <- numeric(T)

  for (i in 1:T) {

    e <- switch(
      dist,

      norm = {
        loc <- args$loc %||% 0
        sd  <- args$sd  %||% 1
        if (sd <= 0) stop("Enter positive value for sd")
        rnorm(T, loc, sd)

      },

      cauchy = {
        loc <- args$loc %||% 0
        scale  <- args$scale  %||% 1
        if (scale <= 0) stop("Enter positive value for scale")
        rcauchy(T, loc, scale)
      },

      uniform = {
        min <- args$min %||% -1
        max <- args$max %||% 1
        runif(T, min, max)
      }
    )
  }

  X = cumsum(e)
  X <- X[is.finite(X)] # drop inf/nan if any
  return(X)
}

############## Generate DTRW Series ##############

#' Generate a Discrete-Time Random Walk (DTRW) Series with Cauchy Noise
#'
#' Generates a DTRW series where the increments follow a Cauchy distribution.
#'
#' @details
#' The series is generated recursively as:
#' \deqn{x_t = x_{t-1} + e_t, \quad e_t \sim \text{Cauchy}(\text{loc}, \text{scale})}
#' where:
#' - \eqn{T} is the length of the series.
#' - \eqn{e_t} are i.i.d. random variables from a Cauchy distribution.
#' - \eqn{\text{loc}} is the location parameter.
#' - \eqn{\text{scale}} is the scale parameter.
#'
#' @param T Integer. The length of the series.
#' @param loc Numeric. The location parameter of the Cauchy distribution (default = 0).
#' @param scale Positive numeric. The scale parameter of the Cauchy distribution (default = 1).
#' @return A numeric vector representing the DTRW series.
#' @export
#' @examples
#' DTRW_series_Cauchy(100, loc = 0, scale = 1)
DTRW_series_Cauchy <- function(T, loc = 0, scale = 1) {
  if (scale <= 0) stop("Enter a positive value for scale")
  e <- rcauchy(T, location = loc, scale = scale) ## Generate increments
  x <- numeric(T)

  for (i in 2:T) {
    x[i] <- x[i - 1] + e[i]
  }

  return(x)
}

#' Generate a Discrete-Time Random Walk (DTRW) Series with Normal Noise
#'
#' Generates a DTRW series where the increments follow a Normal distribution.
#'
#' @details
#' The series is generated recursively as:
#' \deqn{x_t = x_{t-1} + e_t, \quad e_t \sim \mathcal{N}(\text{loc}, \text{sd})}
#' where:
#' - \eqn{T} is the length of the series.
#' - \eqn{e_t} are i.i.d. random variables from a Normal distribution.
#' - \eqn{\text{loc}} is the mean of the Normal distribution.
#' - \eqn{\text{sd}} is the standard deviation of the Normal distribution.
#'
#' @param T Integer. The length of the series.
#' @param loc Numeric. The mean (location) parameter of the normal distribution.
#' @param sd Positive numeric. The standard deviation of the normal distribution.
#' @return A numeric vector representing the DTRW series.
#' @export
#' @examples
#' DTRW_series_Norm(100, loc = 0, sd = 1)
DTRW_series_Norm <- function(T, loc, sd) {
  if (sd <= 0) stop("Enter a positive value for standard deviation")
  # e <- rnorm(T, mean = loc, sd = sd) ## Generate increments
  # x <- numeric(T)

  x <- cumsum(rnorm(T, mean = loc, sd = sd))
  #x <- x - x[1] # Normalize the series so it starts at 0
  return(x)
}

