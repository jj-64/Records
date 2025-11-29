#' Yule--Simon / YNM expected number of records
#' @param T Integer
#' @param gamma Numeric
#' @export
record_mean_YNM <- function(T, gamma) {
  s <- numeric(T)
  for (k in seq_len(T)) s[k] <- rec_rate_YNM(gamma, k)
  sum(s)
}

#' YNM variance
#' @export
record_var_YNM <- function(T, gamma) {
  s <- numeric(T); s2 <- numeric(T)
  for (k in seq_len(T)) {
    s[k] <- rec_rate_YNM(gamma, k)
    s2[k] <- s[k]^2
  }
  sum(s) - sum(s2)
}

#' YNM exact distribution
#' @export
record_dist_YNM <- function(m, T, gamma, s = NULL) {
  if (is.null(s)) s <- Stirling_2nd_YNM(T = T, gamma = gamma)
  p <- prod(u_t_YNM(t = 1:T, gamma = gamma))
  s[T, m] / ((gamma^T) * p)
}
