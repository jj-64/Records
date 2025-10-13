############## Generic function ##################
#' Generate a Linear Drift Model (LDM) Series with different Noise distribution
#'
#' Generates an LDM series where the noise follows a user defined distribution.
#'
#' @details
#' The series is generated as:
#' \deqn{x_t = \theta t + y_t, \quad y_t \sim \text{F}(\text{.}, \beta)}
#' where:
#' - \eqn{\theta} is the linear drift.
#' - \eqn{y_t} are i.i.d. distributed random variables.
#' - \eqn{\beta} are specific parameters of the distribution F of \eqn{Y} .
#'
#' @param T Integer. The length of the series.
#' @param theta Numeric. The linear drift coefficient \eqn{\theta > 0}.
#' @param dist Character, distribution name. One of:
#'   "beta", "gumbel", "weibull", "frechet", "exp", "pareto", "norm", "exp", "pareto", "uniform".
#' @param ... Additional parameters specific to the chosen distribution:
#'   \describe{
#'     \item{beta}{`shape1`, `shape2`}
#'     \item{gumbel}{`loc`, `scale`}
#'     \item{weibull}{`shape`, `scale`}
#'     \item{frechet}{`shape`, `scale`}
#'     \item{norm}{`mean`, `sd`}
#'     \item{exp}{`rate`}
#'     \item{pareto}{`scale`, `shape`}
#'     \item{unifom}{`min`, `max`}
#'   }
#' @return A numeric vector representing the LDM series.
#' @export
#' @examples
#' #Gumbel with loc, scale
#'LDM_series(100, theta = 0.1, dist = "gumbel", loc = 0, scale = 2)
#'
#' #Weibull
#'LDM_series(100, theta = 0.05, dist = "weibull", shape = 2, scale = 1)
#'
#' #Beta
#'LDM_series(100, theta = 0.2, dist = "beta", shape1 = 2, shape2 = 5)
#'
#' #Normal
#'LDM_series(100, theta = 0.1, dist = "norm", mean = 0, sd = 1)
LDM_series <- function(T, theta, dist = c("beta", "gumbel", "weibull", "frechet", "norm", "exp", "pareto", "uniform"), ...) {
  dist <- match.arg(dist)   # enforce valid choice
  args <- list(...)

  # Simulate noise y according to chosen distribution
  y <- switch(
    dist,

    beta = {
      shape1 <- args$shape1 %||% 1
      shape2 <- args$shape2 %||% 1
      if (shape1 <= 0 | shape2 <= 0) stop("Enter positive values for shape1 and shape2")
      rbeta(T, shape1, shape2)
    },

    gumbel = {
      loc   <- args$loc   %||% 0
      scale <- args$scale %||% 1
      if (scale <= 0) stop("Enter a positive value for scale")
      VGAM::rgumbel(T, location = loc, scale = scale)
    },

    weibull = {
      shape <- args$shape %||% 1
      scale <- args$scale %||% 1
      if (scale <= 0 | shape <= 0) stop("Enter positive values for shape and scale")
      rweibull(T, shape = shape, scale = scale)
    },

    frechet = {
      scale <- args$scale %||% 1
      shape <- args$shape %||% 2
      if (scale <= 0 | shape <= 0) stop("Enter positive values for scale and shape")
      VGAM::rfrechet(T, loc = 0, scale = scale, shape = shape)
    },

    norm = {
      mean <- args$mean %||% 0
      sd   <- args$sd   %||% 1
      if (sd <= 0) stop("Enter positive value for sd")
      rnorm(T, mean = mean, sd = sd)
    },

    exp = {
      rate <- args$rate %||% 1
      if (rate <= 0) stop("Enter positive value for rate")
      rexp(T, rate=rate)
    },

    pareto = {
      scale <- args$scale %||% 1
      shape <- args$shape %||% 2
      if (scale <= 0 | shape <= 0) stop("Enter positive values for scale and shape")
      rpareto(T * 1.5, scale, shape)
    },

    uniform ={
      min <- args$min %||% -1
      max <- args$max %||% 1
      runif(T, min, max)
    }

  )

  # Add deterministic linear drift
  x <- theta * (1:T) + y
  return(x)
}


############## Generate LDM series ##############

#' Generate a Linear Drift Model (LDM) Series with Beta Noise
#'
#' Generates an LDM series where the noise follows a Beta distribution.
#'
#' @details
#' The series is generated as:
#' \deqn{x_t = \theta t + y_t, \quad y_t \sim \text{Beta}(\text{shape1}, \text{shape2})}
#' where:
#' - \eqn{\theta} is the linear drift.
#' - \eqn{y_t} are i.i.d. Beta-distributed random variables.
#'
#' @param T Integer. The length of the series.
#' @param theta Numeric. The linear drift coefficient \eqn{\theta > 0}.
#' @param shape1 Positive numeric. The first shape parameter of the Beta distribution.
#' @param shape2 Positive numeric. The second shape parameter of the Beta distribution.
#' @return A numeric vector representing the LDM series.
#' @examples
#' LDM_series_Beta(100, theta = 0.5, shape1 = 2, shape2 = 5)
LDM_series_Beta <- function(T, theta, shape1 = 1, shape2 = 1) {
  if (shape1 <= 0 | shape2 <= 0) stop("Enter positive values for shape parameters")
  y <- rbeta(T, shape1, shape2)
  x <- theta * (1:T) + y
  return(x)
}


#' Generate an LDM Series with Gumbel Noise
#'
#' @details
#' \deqn{x_t = \theta t + y_t, \quad y_t \sim \text{Gumbel}(\text{loc}, \text{scale})}
#'
#' @inheritParams LDM_series_Beta
#' @param loc Numeric. The location parameter of the Gumbel distribution.
#' @param scale Positive numeric. The scale parameter of the Gumbel distribution.
#' @return A numeric vector representing the LDM series.
LDM_series_Gumbel <- function(T, theta, loc = 0, scale = 1) {
  if (scale <= 0) stop("Enter a positive value for scale")
  y <- VGAM::rgumbel(T, loc, scale)
  x <- theta * (1:T) + y
  return(x)
}

#' Generate an LDM Series with Weibull Noise
#'
#' @details
#' \deqn{x_t = \theta t + y_t, \quad y_t \sim \text{Weibull}(\text{scale}, \text{shape})}
#'
#' @inheritParams LDM_series_Beta
#' @param scale Positive numeric. The scale parameter of the Weibull distribution.
#' @param shape Positive numeric. The shape parameter of the Weibull distribution.
#' @return A numeric vector representing the LDM series.
LDM_series_Weibull <- function(T, theta, shape = 1,scale = 1) {
  if (scale <= 0 | shape <= 0) stop("Enter positive values for scale and shape")
  y <- rweibull(T, shape=shape, scale=scale)
  x <- theta * (1:T) + y
  return(x)
}

#' Generate an LDM Series with Frechet Noise
#'
#' @details
#' \deqn{x_t = \theta t + y_t, \quad y_t \sim \text{Frechet}(\text{scale}, \text{shape})}
#'
#' @inheritParams LDM_series_Beta
#' @param scale Positive numeric. The scale parameter of the Frechet distribution.
#' @param shape Positive numeric. The shape parameter of the Frechet distribution.
#' @return A numeric vector representing the LDM series.
LDM_series_Frechet <- function(T, theta, scale = 1, shape = 2) {
  if (scale <= 0 | shape <= 0) stop("Enter positive values for scale and shape")
  y <- VGAM::rfrechet(T, loc = 0, scale, shape)
  x <- theta * (1:T) + y
  return(x)
}

#' Generate an LDM Series with Exponential Noise
#'
#' @details
#' \deqn{x_t = \theta t + y_t, \quad y_t \sim \text{Exp}(\text{rate})}
#'
#' @inheritParams LDM_series_Beta
#' @param rate Positive numeric. The rate parameter of the Exponential distribution.
#' @return A numeric vector representing the LDM series.
LDM_series_Exp <- function(T, theta, rate = 1) {
  if (rate <= 0) stop("Enter a positive value for scale:1/ rate")
  y <- rexp(T, rate=rate)
  x <- theta * (1:T) + y
  return(x)
}

#' Generate an LDM Series with Pareto Noise
#'
#' @details
#' \deqn{x_t = \theta t + y_t, \quad y_t \sim \text{Pareto}(\text{scale}, \text{shape})}
#'
#' The top 10% of values are filtered out to reduce extreme values.
#'
#' @inheritParams LDM_series_Beta
#' @inheritParams LDM_series_Weibull
#' @return A numeric vector representing the LDM series.
LDM_series_Pareto <- function(T, theta, scale = 0.5, shape = 4) {
  if (scale <= 0 | shape <= 0) stop("Enter positive values for scale and shape")
  y <- rpareto(T * 1.5, scale, shape)
  y <- y[y <= quantile(y, 0.9)] # Remove top 10% of extreme values
  x <- theta * (1:T) + y[1:T]
  return(x)
}

#' Generate an LDM Series with Normal Noise
#'
#' @details
#' \deqn{x_t = \theta t + y_t, \quad y_t \sim \mathcal{N}(\text{loc}, \text{sd})}
#'
#' @inheritParams LDM_series_Beta
#' @param loc Numeric. The mean (location) of the normal distribution.
#' @param sd Positive numeric. The standard deviation of the normal distribution.
#' @return A numeric vector representing the LDM series.
LDM_series_Norm <- function(T, theta, loc = 0, sd = 1) {
  if (sd <= 0) stop("Enter a positive value for standard deviation")
  y <- rnorm(T, loc, sd)
  x <- theta * (1:T) + y
  return(x)
}

#' Generate an LDM Series with Uniform Noise
#'
#' @details
#' \deqn{x_t = \theta t + y_t, \quad y_t \sim \text{Uniform}(\text{min}, \text{max})}
#'
#' @inheritParams LDM_series_Beta
#' @param min Numeric. The lower bound of the uniform distribution.
#' @param max Numeric. The upper bound of the uniform distribution.
#' @return A numeric vector representing the LDM series.
LDM_series_Unif <- function(T, theta, min = -1, max = 1) {
  y <- runif(T, min, max)
  x <- theta * (1:T) + y
  return(x)
}




