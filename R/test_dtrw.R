
############## Perform DTRW test ##############
#' Composite Test of Random-Walk Increment Assumptions
#'
#' Performs a composite hypothesis test to assess whether a given time series
#' follows a Discrete-Time Random Walk (DTRW) process. The test combines two
#' component tests — independence, and symmetry — using
#' multiple-comparison correction methods (Bonferroni, Holm-Bonferoni, Holm-Sidák, or chi-square).
#'
#' @details
#' The function applies the following tests:
#' \enumerate{
  # \item \strong{Stationarity:} Augmented Dickey–Fuller (ADF) test
  #       using \code{tseries::adf.test}.
  #       Null hypothesis (\eqn{H_0}): series is non-stationary.
  #       Desired outcome: fail to reject \eqn{H_0} (\eqn{p > \alpha}).
#'
#'   \item \strong{Independence:} Ljung–Box test (\code{Box.test}) on the increments.
#'         Null hypothesis (\eqn{H_0}): increments are independently distributed.
#'         Desired outcome: fail to reject \eqn{H_0}.
#'
#'   \item \strong{Symmetry:} Wilcoxon signed-rank test (\code{wilcox.test}) for zero median increments.
#'         Null hypothesis (\eqn{H_0}): increments are symmetric about zero.
#'         Desired outcome: fail to reject \eqn{H_0}.
#' }
#'
#' The overall decision is based on combining p-values using one of the following methods:
#' \itemize{
#'   \item \code{"Bonf"} — Bonferroni correction (default)
#'   \item \code{"Holm"} — Holm–Bonferroni sequential correction
#'   \item \code{"Sidak"} — Holm-Sidák correction
#'   \item \code{"Chisq"} — Fisher’s chi-squared combination test
#' }
#'
#' The decision rule is "dtrw" if all component tests fail to reject their nulls
#' under the chosen adjustment method; otherwise "no".
#'
#' @param X Numeric vector of observations representing the process path.
#' @param alpha Numeric, significance level (default = 0.05).
#' @param method Character, p-value combination method: one of
#'   \code{"Bonf"}, \code{"Holm"}, \code{"Sidak"}, or \code{"Chisq"} (default = "Bonf").
#'
#' @return A list with the following elements:
#' \item{method}{Method used to combine tests.}
#' \item{alternative_ljung_box} {Alternative hypthesis for Ljung-Box independence test.}
#' \item{p_value_ljung_box}{P-value from Ljung–Box independence test.}
#' \item{alternative_wilcoxon} {Alternative hypthesis for Wilcoxon symmetry test.}
#' \item{p_value_wilcoxon}{P-value from Wilcoxon symmetry test.}
#' \item{decision}{Overall decision: "dtrw" if consistent with random walk assumptions, "no" otherwise.}
#'
#' @references
#' Ljung, G. M. and Box, G. E. P. (1978). “On a Measure of Lack of Fit in Time Series Models.”
#' \emph{Biometrika}, 65(2), 297–303.
#'
#' Holm, S. (1979). “A Simple Sequentially Rejective Multiple Test Procedure.”
#' \emph{Scandinavian Journal of Statistics}, 6(2), 65–70.
#'
#' Fisher, R. A. (1932). \emph{Statistical Methods for Research Workers.}
#' @export
#' @examples
#' set.seed(123)
#' X <- cumsum(rnorm(100))  # Simulated random walk
#' test_dtrw_increment(X, alpha = 0.05, method = "Bonf")
test_dtrw_increment = function(X,alpha=0.05, method="Bonf"){

  X=X-X[1]
  increments <- diff(X)
  ## Test1: Perform ADF test H1: stationary
    #adf_test <- tseries::adf.test(X)
    #p_value1 = adf_test$p.value ## we want H0: series is not stationary, so fail to reject H0, so p_value>alpha

  ## Test2:  Perform Ljung-Box test on increments: H0:independently distributed (no autocorrelation)
    ljung_box_test <- Box.test(increments, type = "Ljung-Box")
    p_value2 = ljung_box_test$p.value  ## we want H0: fail to reject so p_value>p

  ##Test3: Wilcoxon signed-rank test for symmetry H0: symmetric
    wilcoxon_test = wilcox.test(increments, mu =0, alternative = "two.sided", exact = FALSE) #median(increments) # DTRW without drift
    p_value3 = wilcoxon_test$p.value

    m <- 2 # number of tests

  ## Decision logic
    if (method == "Bonf") {
      #dec <- ifelse(p_value2 > alpha / 2 & p_value3 > alpha / 2, "dtrw", "no")
      dec <- ifelse(all(c(p_value2,p_value3) > alpha/m), "dtrw","no")

    } else if (method == "max") {
      dec <- ifelse(all(c(p_value2,p_value3) > alpha), "dtrw","no")

    }
    else if (method == "Holm") {
      pp <- sort(c(p_value2, p_value3))
      #Holm <- pp > c(1 - (1 - alpha)^(1/2), 1 - (1 - alpha)^(1/2), alpha)
      #dec <- ifelse(sum(Holm) == 2, "dtrw", "no") ## Holm-Bonf is true true true
      Holm <- (pp[1] > alpha/2) & (pp[2] > alpha)
      dec <- ifelse(Holm,"dtrw","no")

    } else if (method == "Sidak") {
      # pp <- c( p_value2, p_value3)
      # Sidak <- pp > (1 - (1 - alpha)^(1/2))
      # dec <- ifelse(sum(Sidak) == 2, "dtrw", "no") ## Holm-Sidak is true true true
      critic <- 1-(1-alpha)^(1/m)
      dec <- ifelse( all(c(p_value2,p_value3) > critic), "dtrw","no" )

    } else if (method == "Chisq") {
      stat <- -2 * ( log(p_value2) + log(p_value3)) #log(p_value1)
      dec <- ifelse(stat <= qchisq(1 - alpha, df = 2* m), "dtrw", "no") #high typr I error

    } else {
      stop("method must be one of 'max' ,'Bonf', 'Holm', 'Sidak', or 'Chisq'")
    }

    return(list(
      method = method,
      alpha = alpha,
      alternative_ljung_box = "correlated, not independent",
      p_value_ljung_box = p_value2,
      alternative_wilcoxon = "NotSymmetric",
      p_value_wilcoxon = p_value3,
      decision = dec

    ))
}

#' Composite Test of Random-Walk Increment Assumptions
#'
#' Evaluates whether the increments of a stochastic process are
#' consistent with the fundamental assumptions of a discrete-time
#' random walk (DTRW) without drift. Specifically, the function
#' assesses whether the increment process:
#'
#' \enumerate{
#'   \item is independently distributed;
#'   \item has zero drift (median increment equal to zero);
#'   \item is symmetrically distributed around zero.
#' }
#'
#' The test is implemented as a composite procedure combining
#' the following hypothesis tests:
#'
#' \itemize{
#'   \item \strong{Independence:}
#'         Ljung-Box test applied to the increment series.
#'
#'         Null hypothesis:
#'         \deqn{H_0: \rho_k = 0}
#'         for all tested lags.
#'
#'         Failure to reject indicates no evidence of serial dependence.
#'
#'   \item \strong{Zero Drift:}
#'         Wilcoxon signed-rank test applied to the increment series.
#'
#'         Null hypothesis:
#'         \deqn{H_0: \mathrm{median}(\Delta X_t)=0}
#'
#'         Failure to reject indicates no evidence of systematic upward
#'         or downward drift.
#'
#'   \item \strong{Symmetry:}
#'         Miao-Gel-Gastwirth symmetry test
#'         (\code{lawstat::symmetry.test}).
#'
#'         Null hypothesis:
#'         \deqn{H_0: F(x)=1-F(-x)}
#'
#'         Failure to reject indicates consistency with a symmetric
#'         increment distribution.
#' }
#'
#' Individual p-values are combined using a family-wise error control
#' procedure. A process is classified as \code{"dtrw"} only when
#' all component tests fail to reject their respective null hypotheses
#' after adjustment.
#'
#' Available combination methods:
#'
#' \itemize{
#'   \item \code{"Bonf"}:
#'         Bonferroni adjustment.
#'
#'   \item \code{"Holm"}:
#'         Holm-Bonferroni sequential adjustment.
#'
#'   \item \code{"Sidak"}:
#'         Sidak adjustment.
#' }
#'
#' The Holm procedure is recommended because it controls the
#' family-wise error rate while providing greater power than
#' Bonferroni.
#'
#' @param X Numeric vector representing the observed process path.
#'
#' @param alpha Numeric significance level.
#' Default is \code{0.05}.
#'
#' @param method Character string specifying the multiple-testing
#' adjustment method. Must be one of:
#' \code{"Bonf"}, \code{"Holm"}, or \code{"Sidak"}.
#'
#' @param lag Integer lag used in the Ljung-Box test.
#' Default is \code{floor(log(length(X)))}.
#'
#' @return A list containing:
#'
#' \describe{
#'   \item{method}{
#'     Multiple-testing adjustment method used.
#'   }
#'
#'   \item{alpha}{
#'     Significance level.
#'   }
#'
#'   \item{p_value_independence}{
#'     Ljung-Box p-value.
#'   }
#'
#'   \item{p_value_zero_drift}{
#'     Wilcoxon signed-rank p-value.
#'   }
#'
#'   \item{p_value_symmetry}{
#'     Miao-Gel-Gastwirth symmetry-test p-value.
#'   }
#'
#'   \item{decision}{
#'     Overall classification:
#'     \code{"dtrw"} or \code{"no"}.
#'   }
#' }
#'
#' @references
#' Ljung, G. M. and Box, G. E. P. (1978).
#' On a Measure of Lack of Fit in Time Series Models.
#' Biometrika, 65(2), 297-303.
#'
#' Holm, S. (1979).
#' A Simple Sequentially Rejective Multiple Test Procedure.
#' Scandinavian Journal of Statistics, 6(2), 65-70.
#'
#' Miao, W., Gel, Y., and Gastwirth, J. (2006).
#' A New Test of Symmetry About an Unknown Median.
#' Random Structures and Algorithms, 29(1), 81-97.
#'
#' @export
#'
#' @examples
#' set.seed(123)
#'
#' X <- cumsum(rnorm(1000))
#'
#' test_dtrw_assumptions(
#'   X,
#'   alpha = 0.05,
#'   method = "Holm"
#' )
test_dtrw_assumptions <- function(
    X,
    alpha = 0.05,
    method = c("max","Holm", "Bonf", "Sidak"),
    lag = floor(log(length(X)))
) {

  method <- match.arg(method)

  if (!is.numeric(X))
    stop("'X' must be numeric.")

  if (length(X) < 20)
    stop("'X' must contain at least 20 observations.")

  increments <- diff(X)

  ## Independence
  lb_test <- Box.test(
    increments,
    lag = lag,
    type = "Ljung-Box"
  )

  p_indep <- lb_test$p.value

  ## Zero drift
  drift_test <- wilcox.test(
    increments,
    mu = 0,
    alternative = "two.sided",
    exact = FALSE
  )

  p_drift <- drift_test$p.value

  ## Symmetry
  sym_test <- lawstat::symmetry.test(
    increments,
    option = "MGG"
  )

  p_sym <- sym_test$p.value

  pvals <- c(
    independence = p_indep,
    zero_drift   = p_drift,
    symmetry     = p_sym
  )

  m <- length(pvals)

  if (method == "Bonf") {

    pass <- all(pvals > alpha / m)

  } else if (method == "Holm") {

    pp <- sort(pvals)

    pass <- all(
      pp > alpha / (m:1)
    )

  } else if (method == "Sidak") {

    crit <- 1 - (1 - alpha)^(1 / m)

    pass <- all(
      pvals > crit
    )

  }

  decision <- ifelse(pass, "dtrw", "no")

  return(
    list(
      method = method,
      alpha = alpha,
      p_value_independence = p_indep,
      p_value_zero_drift = p_drift,
      p_value_symmetry = p_sym,
      decision = decision
    )
  )
}
# test_dtrw_indep2 = function(X,alpha=0.05){
#
#   ## Test2:  Perform Ljung-Box test on increments: H0:independently distributed (no autocorrelation)
#   ljung_box_test <- Box.test(diff(X), type = "Ljung-Box")
#   p_value2 = ljung_box_test$p.value  ## we want H0: fail to reject so p_value>p
#
#   ##Test3: Wilcoxon signed-rank test for symmetry H0: symmetric
#   wilcoxon_test = wilcox.test(diff(X), mu =0, alternative = "two.sided", exact = FALSE) #median(increments) # DTRW without drift
#   p_value3 = wilcoxon_test$p.value
#
#   dec = ifelse(p_value2 < alpha/2 | p_value3< alpha/2 , "no", "dtrw")
#   return(list(
#     aleternative = "correlated",
#     p_value = p_value2,
#     decision = dec
#   ))
# }
### based on exact distribution of number of records #######

#' Test for Discrete-Time Random Walk (DTRW) Using Record Counts
#'
#'Tests whether the observed number of records in a time series is
#' consistent with the distribution expected under a discrete-time
#' random walk (DTRW).
#'
#' The test is based on the universal record statistics of random walks
#' with independent and symmetric increments. The observed record count
#' is compared either to exact finite-sample critical bounds or to an
#' asymptotic normal approximation.
#'
#' @details
#'
#' Let
#'
#' \deqn{N_T}
#'
#' denote the number of upper records observed in a series of length
#' \eqn{T} (obtained via \code{rec_count()}).
#'
#' Under the DTRW hypothesis, the distribution of \eqn{N_T} depends
#' only on the random-walk structure and not on the specific increment
#' distribution, provided increments are independent and symmetric.
#'
#'
#'For very large T, under the DTRW null hypothesis, the expected number of records and its variance
#' have asymptotic forms:
#' \deqn{
#'   E[N_T] \approx c_1 \sqrt{T}, \qquad Var[N_T] \approx c_2 T,
#' }
#' where constants \eqn{c_1} and \eqn{c_2} are model-specific.
#' A standardized test statistic is computed as:
#' \deqn{
#'   Z_T = \frac{N_T - E[N_T]}{\sqrt{Var[N_T]}}.
#' }
#'
#' Two testing approaches are available:
#'
#' \strong{Exact test}
#'
#' When \code{approximate = FALSE}, the observed record count is
#' compared to the theoretical acceptance region returned by
#' \code{\link{rec_count_bounds}}.
#'
#' The DTRW hypothesis is not rejected if:
#'
#' \deqn{
#' L_\alpha \le N_T \le U_\alpha
#' }
#'
#' where \eqn{L_\alpha} and \eqn{U_\alpha} denote the lower and upper
#' critical record counts.
#'
#'\strong{Asymptotic test}
#'
#' For large \eqn{T}, the record count admits the approximation
#'
#' \deqn{
#' E(N_T) \approx c_1 \sqrt{T}
#' }
#'
#' and
#'
#' \deqn{
#' Var(N_T) \approx c_2 T.
#' }
#'
#' A standardized statistic is computed as
#'
#' \deqn{
#' Z_T=
#' \frac{N_T-E(N_T)}
#' {\sqrt{Var(N_T)}}.
#' }
#'
#' Under the DTRW hypothesis,
#' \eqn{Z_T} is approximately standard normal.
#'
#' Decision rule:
#' \itemize{
#'   \item For two-sided tests (\code{one.sided = FALSE}), accept DTRW if
#'         the observed statistic lies between the lower and upper quantiles.
#'   \item For one-sided tests (\code{one.sided = TRUE}), accept if below the upper quantile.
#' }
#'
#' @param X Numeric vector representing the observed series.
#' @param alpha Numeric, significance level (default = 0.05).
#' @param approximate Logical, if \code{TRUE} use the asymptotic normal approximation
#'   (default = \code{FALSE} for the exact quantile test).
#' @param one.sided Logical.
#' If \code{TRUE}, performs an upper-tail test against an unusually
#' large number of records (one-sided test).
#' Default is \code{FALSE} for two-sided.
#'
#' @return A list containing:
#'
#' \describe{
#'
#' \item{rec_count}{
#' Observed number of records.
#' }
#'
#' \item{stat}{
#' Standardized test statistic when
#' \code{approximate = TRUE}.
#' }
#'
#' \item{critical_bounds}{
#' Lower and upper record-count bounds when
#' \code{approximate = FALSE}.
#' }
#'
#' \item{p_value}{
#' Asymptotic p-value when
#' \code{approximate = TRUE}.
#' }
#'
#' \item{decision}{
#' \code{"dtrw"} if the null hypothesis is not rejected and
#' \code{"no"} otherwise.
#' }
#'
#' }
#'
#' @references
#' Majumdar, S. N. and Ziff, R. M. (2008).
#' Universal Record Statistics of Random Walks and Lévy Flights.
#' Physical Review Letters, 101(5).
#'
#' Sparre Andersen, E. (1954).
#' On the Fluctuations of Sums of Random Variables.
#' Mathematica Scandinavica, 2, 195-223.
#'
#' @seealso
#' \code{\link{rec_count_bounds}},
#' \code{\link{rec_count}},
#' \code{\link{rec_count_mean_DTRW}},
#' \code{\link{rec_count_var_DTRW}}
#'
#' @examples
#' set.seed(123)
#' X <- cumsum(rnorm(100))
#'
#' # Exact quantile-based test
#' test_dtrw_rec_count(X, alpha = 0.05, approximate = FALSE)
#'
#' # Asymptotic normal approximation
#' test_dtrw_rec_count(X, alpha = 0.05, approximate = TRUE)
#'
#' @export
test_dtrw_rec_count <- function(X, alpha = 0.05, approximate = FALSE, one.sided = FALSE) {

  obs <- rec_count(X)

  if (!approximate) { # exact

    bounds <- rec_count_bounds(alpha = alpha, model = "dtrw", T = length(X), approximate = approximate)
    z_theo = c("lower_bound" = bounds$lower_bound, "upper_bounds" = bounds$upper_bound)

    if (!one.sided) {
      decision <- ifelse(obs <= z_theo[2] & obs >= z_theo[1], "dtrw", "no")
    } else {
      decision <- ifelse(obs <= z_theo[2], "dtrw", "no")
    }

    #return(list(stat = obs, stat_theo = z_theo, decision = decision))

  } else {

    z_obs <- (obs - rec_count_stats("dtrw", stat="mean", T=length(X), approximate=TRUE) )/ sqrt(rec_count_stats("dtrw", stat="var", T=length(X), approximate=TRUE))

    if (!one.sided) {
      decision <- ifelse(abs(z_obs) <= qnorm(1 - alpha / 2), "dtrw", "no")
      p_value <- 2 * pnorm(-abs(z_obs))
    } else {
      decision <- ifelse(z_obs <= qnorm(1 - alpha), "dtrw", "no")
      p_value <- pnorm(abs(z_obs))
    }

    #return(list(stat = z_obs, p_value = p_value, decision = decision))
  }

  return(
    list(
    method = ifelse(approximate, "asymptotic", "exact"),
    rec_count = obs,
    stat = ifelse(approximate, z_obs, NA),
    critical_bounds = ifelse(approximate, NA, list(z_theo)),
    p_value = ifelse(approximate, p_value, NA),
    decision = decision
  ))
}

######################

# Test_DTRW_ENT <- function(X, alpha= 0.05) {
#   # X: numeric time X (one realization)
#
#   T <- length(X)
#
#   # Step 1: identify records
#   records <- record_times(X)  # times of records
#   N_T <- length(records)                      # number of records
#
#   # Step 2: Fit log-log regression across subsamples
#   # Split into blocks to see growth of R_t vs t
#   block_sizes <- floor(seq(T/4, T, length.out = 5))  # subsample lengths
#   R_block <- sapply(block_sizes, function(t) {
#     length(record_values(X[1:t]))  ## number of records by blocks
#   })
#
#   df <- data.frame(logT = log(block_sizes), logR = log(R_block))
#
#   fit <- lm(logR ~ logT, data = df)
#
#   beta_hat <- coef(fit)[2]
#   se_beta  <- summary(fit)$coefficients[2,2]
#
#   # Wald test: H0: beta = 0.5
#   z <- (beta_hat - 0.5) / se_beta
#   pval <- 2 * (1 - pnorm(abs(z)))
#
#   return(list(
#     regression = summary(fit),
#     beta_hat = beta_hat,
#     se = se_beta,
#     z = z,
#     pval = pval,
#     decision = ifelse(pval < alpha, "no", "dtrw")
#   ))
# }

#######################

## FOUR TESTS FOR THE RANDOM WALK HYPOTHESIS - Handa Test 4 - Detection rte of 10% at 5%
# Test_DTRW_LSE<- function(X, alpha= 0.05) {
#   T <- length(X)
#
#   y1=X[-1]
#   y2=X[-length(X)]
#
#   num = sum(y1*y2)
#
#   denom = sum(X^2)
#
#   alpha_hat = num/denom
#
#   T_p = T * (alpha_hat -1 )/sqrt(2)
#
#   #p_val = pnorm(T_p)
#
#   return(list(
#     z = T_p,
#     decision = ifelse(T_p <= -5.79 | T_p > 0.922, "no", "dtrw")
#   ))
# }

## FOUR TESTS FOR THE RANDOM WALK HYPOTHESIS - Handa Test 8 - Detection rte of 10% at 5%
# Test_DTRW_tratio <- function(X, alpha= 0.05) {
#   T <- length(X)
#
#   y1=X[-1]  ##yt
#   y2=X[-length(X)]  ##yt-1
#
#   alpha_hat = sum(y1*y2)/sum(X^2)
#
#   denom = y1-alpha_hat*y2
#
#   T_p = T * sqrt(sum(X^2)) * (alpha_hat -1 )/sum(denom^2)
#
#   #p_val = pnorm(T_p)
#
#   return(list(
#     z = T_p,
#     "lb" = -5.79,
#     "ub"= 0.922,
#     decision = ifelse(T_p <= -5.79 | T_p > 0.922, "no", "dtrw") ## Null is DTRW
#   ))
# }
