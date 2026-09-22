# Based on number of records --------
#' Arnold's Classical i.i.d. Test Based on Record Counts
#'
#' Performs a record-based goodness-of-fit test for the classical
#' independent and identically distributed (i.i.d.) model using the
#' asymptotic distribution of the number of upper records.
#'
#' @details
#'
#' Let
#'
#' \deqn{N_T}
#'
#' denote the observed number of upper records in a sample of size
#' \eqn{T}.
#'
#' For an i.i.d. sequence with a continuous distribution, the expected
#' number of records satisfies
#'
#' \deqn{
#' E(N_T) \approx \log(T)
#' }
#'
#' and
#'
#' \deqn{
#' Var(N_T) \approx \log(T).
#' }
#'
#' Consequently,
#'
#' \deqn{
#' Z_T
#' =
#' \frac{N_T-\log(T)}
#' {\sqrt{\log(T)}}
#' }
#'
#' converges in distribution to a standard normal random variable as
#' \eqn{T \to \infty}.
#'
#' Large positive or negative values of \eqn{Z_T} indicate departures
#' from the classical i.i.d. record process.
#'
#' The test therefore evaluates:
#'
#' \deqn{
#' H_0:
#' \text{observations are i.i.d.}
#' }
#'
#' against
#'
#' \deqn{
#' H_A:
#' \text{observations are not i.i.d.}
#' }.
#'
#' The decision rule:
#' \itemize{
#'   \item Reject H0 (no) if \eqn{p\_value < \alpha}.
#'   \item Otherwise, fail to reject H0 (classical i.i.d.).
#' }
#'
#' @param X Numeric vector of observations.
#'
#' @param alpha Significance level.
#' Default is \code{0.05}.
#'
#' @return A list containing:
#'
#' \describe{
#'
#' \item{rec_count}{
#' Observed number of records.
#' }
#'
#' \item{statistic}{
#' Standardized Arnold test statistic.
#' }
#'
#' \item{p_value}{
#' Two-sided asymptotic p-value.
#' }
#'
#' \item{decision}{
#' \code{"iid"} if the null hypothesis is not rejected and
#' \code{"no"} otherwise.
#' }
#'
#' }
#'
#'
#' @references
#' Arnold, B. C., Balakrishnan, N. and Nagaraja, H. N. (1998).
#' Records. Wiley.
#'
#' Nevzorov, V. B. (2001).
#' Records: Mathematical Theory.
#' American Mathematical Society.
#'
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#' test_iid_rec_count(x)
#'
#' # $rec_count
#' # [1] 7
#'
#' # $statistic
#' # [1] 1.115968
#'
#' # $p_value
#' # [1] 0.2644358
#'
#' # $decision
#' # [1] "iid"
#'
#' @export
test_iid_rec_count <- function(X, alpha = 0.05) {

  T <- length(X)

  n_records <- rec_count(X)

  z <- (n_records - log(T)) / sqrt(log(T))

  p_value <- 2 * pnorm(-abs(z))

  decision <- ifelse(
    p_value > alpha,
    "iid",
    "no"
  )

  list(
    rec_count = n_records,
    statistic = z,
    p_value = p_value,
    decision = decision
  )
}

## Ljung-Box --------------

#' Ljung-Box Test for Serial Independence
#'
#' Performs the Ljung-Box (or Box-Pierce) portmanteau test to assess
#' whether a time series exhibits statistically significant
#' autocorrelation up to a specified lag.
#'
#' @details
#' The test evaluates the joint null hypothesis that all
#' autocorrelations up to lag \code{lags} are equal to zero:
#'
#' \deqn{
#' H_0:
#' \rho_1=\rho_2=\cdots=\rho_h=0
#' }
#'
#' where \eqn{h} denotes the maximum lag included in the test.
#'
#' Failure to reject the null hypothesis indicates that the series is
#' consistent with serial independence (white noise) up to the tested
#' lag.
#'
#' Note that absence of autocorrelation does not imply that the data are
#' independent and identically distributed (i.i.d.). Nonlinear
#' dependence and heteroskedasticity may still be present even when the
#' Ljung-Box test is not significant.
#'
#' Two portmanteau statistics are available:
#'
#' \itemize{
#'   \item \code{"Ljung-Box"}:
#'   Finite-sample corrected statistic (recommended).
#'
#'   \item \code{"Box-Pierce"}:
#'   Classical large-sample approximation.
#' }
#'
#' @param X Numeric vector representing a time series.
#'
#' @param lags Integer.
#' Number of autocorrelation lags included in the test.
#' Default is \code{10}.
#'
#' @param alpha Significance level.
#' Default is \code{0.05}.
#'
#' @param type Character string specifying the portmanteau statistic.
#' One of:
#' \code{"Ljung-Box"} or \code{"Box-Pierce"}.
#'
#' @return A list containing:
#'
#' \describe{
#'
#' \item{method}{
#' Test statistic used.
#' }
#'
#' \item{statistic}{
#' Observed portmanteau statistic.
#' }
#'
#' \item{p_value}{
#' Test p-value.
#' }
#'
#' \item{decision}{
#' \code{"Independent"} if the null hypothesis is not rejected,
#' otherwise \code{"NO"}.
#' }
#'
#' }
#'
#' @references
#' Box, G. E. P. and Pierce, D. A. (1970).
#' Distribution of Residual Autocorrelations in Autoregressive
#' Integrated Moving Average Time Series Models.
#' Journal of the American Statistical Association,
#' 65(332), 1509-1526.
#'
#' Ljung, G. M. and Box, G. E. P. (1978).
#' On a Measure of Lack of Fit in Time Series Models.
#' Biometrika, 65(2), 297-303.
#'
#' @examples
#' set.seed(123)
#' X <- rnorm(100)
#'
#' test_iid_serial_independence(X)
#'
#' @export
test_iid_serial_independence <- function(
    X,
    lags = 10,
    alpha = 0.05,
    type = c("Ljung-Box", "Box-Pierce")
    ) {

  type <- match.arg(type)

  if(is.null(lags)) lags = floor(log(length(X)))

  bt <- Box.test(
    X,
    lag = lags,
    type = type
  )

  decision <- ifelse(
    bt$p.value > alpha,
    "iid",
    "no"
  )

  list(
    method = type,
    statistic = unname(bt$statistic),
    p_value = bt$p.value,
    decision = decision
  )
}

# test_iid_serial_independence <- function(X, lags = 10, alpha = 0.05, type = "Ljung-Box") {
#   bt <- Box.test(X, lag = lags, type = type)
#
#   decision <- ifelse(bt$p.value > alpha, "iid", "no")
#
#   return(list(
#     stat = bt$statistic,
#     p_value = bt$p.value,
#     decision = decision
#   ))
# }
