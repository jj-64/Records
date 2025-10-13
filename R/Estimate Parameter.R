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
#' @examples
#' X = c(0.428,1.311,2.023,2.882,2.096,-0.197,1.339,  1.748,1.418, 0.711, 1.999,3.598, 3.308,
#' 3.942,2.025,3.282,4.043, 0.492, 4.639, 1.408, 3.525, 5.398,  3.719, 3.741, 4.729)
#'
#' Estimate_model_param(X, method="NT", model = "LDM", scale=1)
Estimate_model_param <- function(X, method = c("NT", "NT_unbiased", "MLE_indicator"), model="LDM", variance = TRUE, ...) {
  method <- match.arg(method)
  model <- match.arg(model)
  args <- list(...)

  if (model == "LDM"){
    if (method == "NT") {
      est <- Estimate_LDM_NT(X, variance = variance, scale=args$scale)

    } else if (method == "NT_unbiased") {
      est = Estimate_LDM_NT_unbiased(X, variance = variance, scale=args$scale)

    } else if (method == "MLE_indicator") {
      est = Estimate_LDM_MLE_indicator(X, variance = variance, scale=args$scale, min = args$min, max=args$max, step=args$step)
    }
  } else if(model == "YNM"){
    if (method == "NT") {
      est <- Estimate_YNM_NT(X, variance = variance)

    } else if (method == "NT_unbiased") {
      est = Estimate_LDM_NT_unbiased(X, variance = variance)

    } else if (method == "MLE_indicator") {
      est = Estimate_YNM_MLE_indicator(X, variance = variance, approximate=args$approximate, min = args$min, max=args$max, step=args$step)
    }
  }
  return(est)
}

