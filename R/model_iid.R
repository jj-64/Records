#' Expected number of records (IID)
#'
#' @param T Integer.
#' @param approx Logical; if TRUE use asymptotic approximation.
#' @return Numeric expected number of records.
#' @export
record_mean_iid <- function(T, approx = FALSE) {
  if (approx) {
    return(log(T) + 0.57721566490153)
  } else {
    return(sum(1 / (1:T)))
  }
}

#' Variance of records (IID)
#' @param T Integer.
#' @export
record_var_iid <- function(T) {
  sum(1 / (1:T)) - sum(1 / (1:T)^2)
}

#' Exact distribution of number of records (IID)
#' @param m Integer.
#' @param T Integer.
#' @param s Optional precomputed Stirling matrix (first-kind absolute values).
#' @export
record_dist_iid <- function(m, T, s = NULL) {
  if (is.null(s)) {
    if (!requireNamespace(\"partitions\", quietly = TRUE)) {
      # Try simple fallback: compute unsigned Stirling numbers of the first kind for small T
      stop(\"Please install 'partitions' or provide 's' (Stirling matrix) for exact IID distribution.\")
    }
    # If partitions is available, user can precompute; here we assume 's' provided.
  }
  s / factorial(T)
}
