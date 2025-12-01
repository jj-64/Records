#' Estimate theta parameter in LDM model
#'
#' Estimates the drift parameter theta in the Linear Drift Model (LDM)
#' using different methods: "moments", "moments_unbias", or "MLE".
#'
#' @param X Numeric vector of the LDM process.
#' @param method Character string, one of "moments", "moments_unbias", or "MLE_indicator".
#' @param bias Logical, if biased estimator when moments estimation (default = TRUE)
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
#' estimate_model_param(X, method="moments", bias = TRUE, model = "LDM", scale=1)
#' #  $param
#' # [1] 0.4462871
#'
#' #  $variance
#' # [1] 0.5625
#'
#'> estimate_model_param(X, method="moments", bias = FALSE, model = "LDM", scale=1)
#' #  $param
#' #  [1] 0.3601306
# '
#' #  $variance
#' #  [1] 0.3116952
#'
#' estimate_model_param(X, method="moments", model = "YNM")
#'
#' # $param
#' # 1.5625
#'
#' # $variance
#' # [1] 1.373291
#'
#' estimate_model_param(X, method="moments", bias = FALSE, model = "YNM")
#' # $param
#' # [1] 1.388625
#'
#' # $variance
#' # [1] 0.7497066
#'
#' estimate_model_param(X, method="MLE_indicator", model = "YNM")
#' # $param
#' # [1] 1.388625
#'
#' # $variance
#' # [1] 0.7497066
estimate_model_param <- function(X, method = c("moments","MLE_indicator"), bias = TRUE ,model=c("LDM","YNM"), variance = TRUE, obs_type = "records", ...) {
  method <- match.arg(method)
  model <- match.arg(model)
  args <- list(...)

  if(obs_type == "records"){
  if (model == "LDM"){
    if (method == "moments" && bias == TRUE) {
      est <- Estimate_LDM_NT(X, variance = variance, scale=args$scale)

    } else if (method == "moments" && bias == FALSE) {
      est = Estimate_LDM_NT_unbiased(X, variance = variance, scale=args$scale)

    } else if (method == "MLE_indicator") {
      est = Estimate_LDM_MLE_indicator(X, variance = variance, scale=args$scale, min = args$min, max=args$max, step=args$step)

    }
  } else if(model == "YNM"){
    if (method == "moments" && bias == TRUE) {
      est <- Estimate_YNM_NT(X, variance = variance)

    } else if (method == "moments" && bias == FALSE) {
      est = Estimate_YNM_NT_unbiased(X, variance = variance)

    } else if (method == "MLE_indicator") {
      est = Estimate_YNM_MLE_indicator(X, variance = variance, approximate=args$approximate, min = args$min, max=args$max, step=args$step)
    }
  }
  }
  return(est)
}

