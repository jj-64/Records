#' LDM expected records (mix of closed form and simulation)
#' @param T Integer
#' @param theta Numeric parameter
#' @param dist Character distribution; one of several supported
#' @param n_sim Number of simulation replicates for simulation-based estimate
#' @param ... Additional args passed to LDM_series or rec_rate_LDM when needed
#' @export
record_mean_LDM <- function(T, theta, dist = c(\"beta\", \"gumbel\", \"weibull\", \"frechet\", \"norm\", \"exp\", \"pareto\", \"uniform\"), n_sim = 1000, ...) {
  dist <- match.arg(dist)
  args <- list(...)
  if (dist == \"gumbel\") {
    s <- rec_rate_LDM(t = 1:T, theta = theta, location = args$location, scale = args$scale)
    return(sum(s))
  } else {
    recs <- numeric(n_sim)
    for (i in seq_len(n_sim)) {
      X <- LDM_series(T = T, theta = theta, dist = dist, ...)
      recs[i] <- rec_counts(X)
    }
    mean(recs)
  }
}

#' LDM variance
#' @export
record_var_LDM <- function(T, theta, dist = c(\"beta\", \"gumbel\", \"weibull\", \"frechet\", \"norm\", \"exp\", \"pareto\", \"uniform\"), n_sim = 1000, ...) {
  dist <- match.arg(dist)
  args <- list(...)
  if (dist == \"gumbel\") {
    s <- rec_rate_LDM(t = 1:T, theta = theta, location = args$location, scale = args$scale)
    return(sum(s * (1 - s)))
  } else {
    recs <- numeric(n_sim)
    for (i in seq_len(n_sim)) {
      X <- LDM_series(T = T, theta = theta, dist = dist, ...)
      recs[i] <- rec_counts(X)
    }
    var(recs)
  }
}

#' LDM exact distribution of record counts (requires Stirling-like objects/functions)
#' @export
record_dist_LDM <- function(m, T, theta, scale = 1, s = NULL) {
  if (is.null(s)) s <- Stirling_2nd_LDM(T = T, theta = theta, scale = scale)
  p <- prod(u_t_LDM(t = 1:T, theta = theta, scale = scale))
  exp(-theta * T / scale) * s[T, m] / p
}
