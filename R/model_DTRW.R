#' DTRW expected number of records
#' @param T Integer
#' @param approx Logical
#' @export
record_mean_DTRW <- function(T, approx = FALSE) {
  if (approx) {
    return(sqrt(4 * T / pi))
  } else {
    return((2 * T + 1) * choose(2 * T, T) * 2^(-2 * T))
  }
}

#' DTRW variance
#' @export
record_var_DTRW <- function(T, approx = FALSE) {
  if (approx) {
    return(2 * (1 - 2 / pi) * T)
  } else {
    m <- record_mean_DTRW(T, approx = FALSE)
    v <- 2 * T + 2 - m - m^2
    return(v)
  }
}

#' DTRW distribution (exact or approx)
#' @export
record_dist_DTRW <- function(m, T, approx = FALSE) {
  if (approx) {
    return(exp(-m^2 / (4 * T)) / sqrt(pi * T))
  } else {
    return(choose(2 * T - m + 1, T) * 2^(-2 * T + m - 1))
  }
}
