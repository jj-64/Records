#' Estimate theta parameter in LDM model
#'
#' Estimates the drift parameter theta in the Linear Drift Model (LDM)
#' using different methods: "NT", "NT_unbiased", or "MLE".
#'
#' @param X Numeric vector of the LDM process.
#' @param method Character string, one of "NT", "NT_unbiased", or "MLE".
#' @param model Character string, one of "LDM", "YNM"
#' @param variance logical to return variance or nor (default = TRUE)
#' @param scale Scale parameter of the underlying Gumbel distribution (default = 1).
#' @param min,max Search bounds for MLE method.
#' @return Numeric estimate of theta and variance (if TRUE)
#' @export
Estimate_model_param <- function(X, method = c("NT", "NT_unbiased", "MLE_indicator"), model="LDM", variance = TRUE,
                                 scale = 1, min = 0.0001, max = 5) {
  method <- match.arg(method)
  model <- match.arg(model)

  if (model == "LDM"){
    if (method == "NT") {
      est <- Estimate_LDM_NT(X, variance = variance, scale=scale)

    } else if (method == "NT_unbiased") {
      est = Estimate_LDM_NT_unbiased(X, variance = variance, scale=scale)

    } else if (method == "MLE_indicator") {
      d <- is_rec(X)
      NT <- rec_counts(X)
      T <- length(X)
      LogT_LDM <- function(z, T, delta, N) {
        v <- function(t, z) log(1 - z^(t - 1))
        x <- sapply(2:T, v, z = z)
        N * log(1 - z) + (T - N) * log(z) - log(1 - z^T) - sum(delta * x)
      }
      z1 <- seq(min, max, by = 0.01)
      z <- exp(-z1)
      likel <- sapply(z, LogT_LDM, T = T, delta = d, N = NT)
      theta <- -log(z[which.max(likel)])
    }}
  return(theta)
}

