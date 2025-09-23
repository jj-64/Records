################## Generic function ###############
#' Generate a Yang Process
#'
#' Simulates a Yang process of length T under different underlying distributions.
#' The Yang process is defined by a parameter gamma > 0 controlling the block maxima.
#'
#' @param T Integer, length of the series.
#' @param gamma Numeric > 0, Yang process parameter.
#' @param dist Character, distribution name. One of:
#'   "beta", "gumbel", "weibull", "frechet", "exp", "pareto", "norm".
#' @param ... Additional parameters specific to the chosen distribution:
#'   \describe{
#'     \item{beta}{`shape1`, `shape2`}
#'     \item{gumbel}{`loc`, `scale`}
#'     \item{weibull}{`shape`, `scale`}
#'     \item{frechet}{`shape`, `scale`}
#'     \item{exp}{`rate`}
#'     \item{pareto}{`scale`, `shape`}
#'     \item{norm}{`loc`, `sd`}
#'   }
#'
#' @return A numeric vector of length T, the simulated Yang process.
#' @examples
#' Yang_series(100, gamma = 1.5, dist = "gumbel", loc = 0, scale = 1)
#' Yang_series(100, gamma = 2, dist = "beta", shape1 = 2, shape2 = 5)
#' Yang_series(100, gamma = 1.2, dist = "norm", loc = 0, sd = 1)
#' @export
Yang_series <- function(T, gamma, dist = c("beta", "gumbel", "weibull", "frechet", "exp", "pareto", "norm"), ...) {
  dist <- match.arg(dist)
  args <- list(...)
  X <- numeric(T)

  for (i in 1:T) {
    m <- gamma^i
    u <- runif(1)

    X[i] <- switch(
      dist,

      beta = {
        shape1 <- args$shape1 %||% 1
        shape2 <- args$shape2 %||% 1
        if (shape1 <= 0 | shape2 <= 0) stop("Enter positive values for shape1 and shape2")
        max(rbeta(m, shape1, shape2))
      },

      gumbel = {
        loc   <- args$loc   %||% 0
        scale <- args$scale %||% 1
        if (scale <= 0) stop("Enter a positive value for scale")
        loc - scale * log(-log(u^(1/m)))
      },

      weibull = {
        shape <- args$shape %||% 1
        scale <- args$scale %||% 1
        if (scale <= 0 | shape <= 0) stop("Enter positive values for shape and scale")
        (-scale^shape * log(1 - u^(1/m)))^(1/shape)
      },

      frechet = {
        shape <- args$shape %||% 2
        scale <- args$scale %||% 1
        if (scale <= 0 | shape <= 0) stop("Enter positive values for scale and shape")
        scale * (-(1/m) * log(u))^(-1/shape)
      },

      exp = {
        rate <- args$rate %||% 1
        if (rate <= 0) stop("Enter positive value for rate")
        -log(1 - u^(1/m)) / rate
      },

      pareto = {
        scale <- args$scale %||% 0.5
        shape <- args$shape %||% 4
        if (scale <= 0 | shape <= 0) stop("Enter positive values for scale and shape")
        scale * (1 - u^(1/m))^(-1/shape)
      },

      norm = {
        loc <- args$loc %||% 0
        sd  <- args$sd  %||% 1
        if (sd <= 0) stop("Enter positive value for sd")
        qnorm(u^(1/m), mean = loc, sd = sd)
      }
    )
  }

  X <- X[is.finite(X)] # drop inf/nan if any
  return(X)
}



############## Generate Yang Series ##############

#' Generate a Yang Series with Beta Noise
#'
#' Generates a Yang series where the noise follows a Beta distribution.
#'
#' @details
#' The series is generated as:
#' \deqn{X_k = \max(Y_k), \quad Y_k \sim \text{Beta}(\text{shape1}, \text{shape2}), \quad \text{length}(Y_k) = \gamma^k}
#' where:
#' - \eqn{\gamma} controls the sample size growth.
#' - \eqn{X_k} represents the maximum value of the Beta-distributed sample.
#'
#' @param T Integer. The length of the series.
#' @param gamma Positive numeric. Growth parameter controlling sample size \eqn{\gamma \geq 1}.
#' @param shape1 Positive numeric. The first shape parameter of the Beta distribution.
#' @param shape2 Positive numeric. The second shape parameter of the Beta distribution.
#' @return A numeric vector representing the Yang series.
#' @export
#' @examples
#' Yang_series_Beta(100, gamma = 1.1, shape1 = 2, shape2 = 5)
Yang_series_Beta <- function(T, gamma, shape1 = 1, shape2 = 1) {
  if (shape1 <= 0 | shape2 <= 0) stop("Enter positive values for shape parameters")

  X <- numeric(T)
  for (k in 1:T) {
    y <- rbeta(gamma^k, shape1, shape2)
    X[k] <- max(y)
  }
  return(X)
}

#' Generate a Yang Series with Gumbel Noise
#'
#' @details
#' \deqn{X_i = \text{loc} - \text{scale} \cdot \log(-\log(U^{1/\gamma^i}))}
#' where \eqn{U \sim \text{Uniform}(0,1)}.
#'
#' @inheritParams Yang_series_Beta
#' @param loc Numeric. The location parameter of the Gumbel distribution.
#' @param scale Positive numeric. The scale parameter of the Gumbel distribution.
#' @return A numeric vector representing the Yang series.
#' @export
Yang_series_Gumbel <- function(T, gamma, loc = 0, scale = 1) {
  if (scale <= 0) stop("Enter a positive value for scale")

  y <- numeric(T)
  for (i in 1:T) {
    y[i] <- loc - scale * log(-log(runif(1)^(1/gamma^i)))
  }
  return(y)
}

#' Generate a Yang Series with Weibull Noise
#'
#' @details
#' \deqn{X_i = (-\text{scale}^\text{shape} \cdot \log(1 - U^{1/\gamma^i}))^{1/\text{shape}}}
#' where \eqn{U \sim \text{Uniform}(0,1)}.
#'
#' @inheritParams Yang_series_Beta
#' @param scale Positive numeric. The scale parameter of the Weibull distribution.
#' @param shape Positive numeric. The shape parameter of the Weibull distribution.
#' @return A numeric vector representing the Yang series.
#' @export
Yang_series_Weibull <- function(T, gamma, scale = 1, shape = 1) {
  if (scale <= 0 | shape <= 0) stop("Enter positive values for scale and shape")

  X <- numeric(T)
  for (i in 1:T) {
    X[i] <- (-scale^shape * log(1 - runif(1)^(1/gamma^i)))^(1/shape)
  }
  return(X)
}

#' Generate a Yang Series with Frechet Noise
#'
#' @details
#' \deqn{X_i = \text{scale} \cdot \left(-\frac{1}{\gamma^i} \log(U) \right)^{-1/\text{shape}}}
#' where \eqn{U \sim \text{Uniform}(0,1)}.
#'
#' @inheritParams Yang_series_Beta
#' @param scale Positive numeric. The scale parameter of the Frechet distribution.
#' @param shape Positive numeric. The shape parameter of the Frechet distribution.
#' @return A numeric vector representing the Yang series.
#' @export
Yang_series_Frechet <- function(T, gamma, shape, scale) {
  if (scale <= 0 | shape <= 0) stop("Enter positive values for scale and shape")

  y <- numeric(T)
  for (i in 1:T) {
    y[i] <- scale * (-(1/gamma^i) * log(runif(1)))^(-1/shape)
  }
  return(y)
}

#' Generate a Yang Series with Exponential Noise
#'
#' @details
#' \deqn{X_i = -\frac{\log(1 - U^{1/\gamma^i})}{\text{rate}}}
#' where \eqn{U \sim \text{Uniform}(0,1)}.
#'
#' @inheritParams Yang_series_Beta
#' @param rate Positive numeric. The rate parameter of the Exponential distribution.
#' @return A numeric vector representing the Yang series.
#' @export
Yang_series_Exp <- function(T, gamma, rate = 1) {
  if (rate <= 0) stop("Enter a positive value for rate")

  y <- numeric(T)
  for (i in 1:T) {
    y[i] <- -log(1 - runif(1)^(1/gamma^i)) / rate
  }
  return(y)
}

#' Generate a Yang Series with Pareto Noise
#'
#' @details
#' \deqn{X_i = \text{scale} \cdot (1 - U^{1/\gamma^i})^{-1/\text{shape}}}
#' where \eqn{U \sim \text{Uniform}(0,1)}.
#'
#' @inheritParams Yang_series_Beta
#' @inheritParams Yang_series_Weibull
#' @return A numeric vector representing the Yang series.
#' @export
Yang_series_Pareto <- function(T, gamma, scale = 0.5, shape = 4) {
  if (scale <= 0 | shape <= 0) stop("Enter positive values for scale and shape")

  X <- numeric(T)
  for (i in 1:T) {
    X[i] <- scale * (1 - runif(1)^(1/gamma^i))^(-1/shape)
  }
  return(X)
}

#' Generate a Yang Series with Normal Noise
#'
#' @details
#' The series is generated as:
#' \deqn{X_i = \max(Y_i), \quad Y_i \sim \mathcal{N}(\text{loc}, \text{sd}), \quad \text{length}(Y_i) = \gamma^i}
#'
#' @inheritParams Yang_series_Beta
#' @param loc Numeric. The mean (location) of the normal distribution.
#' @param sd Positive numeric. The standard deviation of the normal distribution.
#' @return A numeric vector representing the Yang series.
#' @export
Yang_series_Norm <- function(T, gamma, loc = 0, sd = 1) {
if (sd <= 0) stop("Enter a positive value for standard deviation")

# Pre-allocate
y <- numeric(T)

for (i in 1:T) {
  m <- gamma^i
  u <- runif(1)                # uniform random variable
  y[i] <- qnorm(u^(1/m), mean = loc, sd = sd)
}
y= y[is.finite(y)]

return(y)
}

