# helper for defaults
`%||%` <- function(a, b) if (!is.null(a)) a else b

########################################################################
##To compute the bounds/ CI intervals
#' bounds: Confidence-Interval
#'
#'Compute the confidence interval bounds
#' @param value the value or the mean of which we are computing the CI
#' @param z number for z-score from normal distribution
#' @param variance number for the variance (not standard deviation)
#'
#' @returns vector of two: lower and upper bounds of C.I.
#' @export
#'
#' @examples bounds(5, z=qnorm(1-0.05/2), variance = 2)
#' 2.228192 7.771808
bounds = function(value,z,variance){
  upper = value+ z*sqrt(variance)
  lower = value- z*sqrt(variance)
  return(c(lower,upper))
}
