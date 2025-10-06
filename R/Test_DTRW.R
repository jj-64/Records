
############## Perform DTRW test ##############
#' Test for Discrete-Time Random Walk (DTRW) Structure
#'
#' Performs a composite hypothesis test to assess whether a given time series
#' follows a Discrete-Time Random Walk (DTRW) process. The test combines three
#' component tests — stationarity, independence, and symmetry — using
#' multiple-comparison correction methods (Bonferroni, Holm-Bonferoni, Holm-Sidák, or chi-square).
#'
#' @details
#' The function applies the following tests:
#' \enumerate{
#'   \item \strong{Stationarity:} Augmented Dickey–Fuller (ADF) test
#'         using \code{tseries::adf.test}.
#'         Null hypothesis (\eqn{H_0}): series is non-stationary.
#'         Desired outcome: fail to reject \eqn{H_0} (\eqn{p > \alpha}).
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
#' The decision rule is `"DTRW"` if all component tests fail to reject their nulls
#' under the chosen adjustment method; otherwise `"NO"`.
#'
#' @param X Numeric vector of observations representing the process path.
#' @param alpha Numeric, significance level (default = 0.05).
#' @param method Character, p-value combination method: one of
#'   \code{"Bonf"}, \code{"Holm"}, \code{"Sidak"}, or \code{"Chisq"} (default = "Bonf").
#'
#' @return A list with the following elements:
#' \item{method}{Method used to combine tests.}
#' \item{p_valueStationary}{P-value from ADF stationarity test.}
#' \item{p_valueIndep}{P-value from Ljung–Box independence test.}
#' \item{p_valueSymm}{P-value from Wilcoxon symmetry test.}
#' \item{decision}{Overall decision: `"DTRW"` if consistent with random walk assumptions, `"NO"` otherwise.}
#'
#' @references
#' Ljung, G. M. and Box, G. E. P. (1978). “On a Measure of Lack of Fit in Time Series Models.”
#' \emph{Biometrika}, 65(2), 297–303.
#'
#' Holm, S. (1979). “A Simple Sequentially Rejective Multiple Test Procedure.”
#' \emph{Scandinavian Journal of Statistics}, 6(2), 65–70.
#'
#' Fisher, R. A. (1932). \emph{Statistical Methods for Research Workers.}
#'
#' @examples
#' set.seed(123)
#' X <- cumsum(rnorm(100))  # Simulated random walk
#' Test_DTRW(X, alpha = 0.05, method = "Bonf")
#'
#' @export
Test_DTRW = function(X,alpha=0.05, method="Bonf"){

  X=X-X[1]
  increments <- diff(X)
  ## Test1: Perform ADF test H1: stationary
    adf_test <- tseries::adf.test(X)
    p_value1 = adf_test$p.value ## we want H0: series is not stationary, so fail to reject H0, so p_value>alpha

  ## Test2:  Perform Ljung-Box test on increments: H0:independently distributed (no autocorrelation)
    ljung_box_test <- Box.test(increments, type = "Ljung-Box")
    p_value2 = ljung_box_test$p.value  ## we want H0: fail to reject so p_value>p

  ##Test3: Wilcoxon signed-rank test for symmetry H0: symmetric
    wilcoxon_test = wilcox.test(increments, mu =0, alternative = "two.sided", exact = FALSE) #median(increments) # DTRW without drift
    p_value3 = wilcoxon_test$p.value


  ## Decision logic
    if (method == "Bonf") {
      dec <- ifelse(p_value1 > alpha / 3 & p_value2 > alpha / 3 & p_value3 > alpha / 3, "DTRW", "NO")

    } else if (method == "Holm") {
      pp <- sort(c(p_value1, p_value2, p_value3))
      Holm <- pp > c(1 - (1 - alpha)^(1/3), 1 - (1 - alpha)^(1/2), alpha)
      dec <- ifelse(sum(Holm) == 3, "DTRW", "NO") ## Holm-Bonf is true true true

    } else if (method == "Sidak") {
      pp <- c(p_value1, p_value2, p_value3)
      Sidak <- pp > (1 - (1 - alpha)^(1/3))
      dec <- ifelse(sum(Sidak) == 3, "DTRW", "NO") ## Holm-Sidak is true true true

    } else if (method == "Chisq") {
      stat <- -2 * (log(p_value1) + log(p_value2) + log(p_value3))
      dec <- ifelse(stat <= qchisq(1 - alpha, df = 6), "DTRW", "NO") #high typr I error

    } else {
      stop("method must be one of 'Bonf', 'Holm', 'Sidak', or 'Chisq'")
    }

    return(list(
      method = method,
      p_valueStationary = p_value1,
      p_valueIndep = p_value2,
      p_valueSymm = p_value3,
      decision = dec
    ))
}


############## based on exact distribution of number of records ##############
Quantile_DTRW=function(alpha=0.05, T){
  Prob = 0
  for(i in 1:T){
    Prob[i]= NT_DTRW(m=i, T=T)
  }
  CDF = cumsum(Prob)  ## cumulative distribution
  #plot(x=1:T, CDF)

  #return( c( which.min(abs(CDF-(alpha/2))) ,which.min(abs(CDF-(1-alpha/2)))))
  return( c( which.min(abs(CDF-(alpha))) ,which.min(abs(CDF-(1-alpha)))))
}

Test_DTRW_NT =function(X,alpha=0.05){
  z_theo = Quantile_DTRW(alpha=alpha, T=length(X))

  obs = rec_counts(X)
  #decision
  decision=ifelse(obs<=z_theo[2] & obs >= z_theo[1] , "DTRW","NO")

  return(list("stat"=obs,"stat_theo"=z_theo,"decision"=decision))
}

######################

Test_DTRW_ENT <- function(X, alpha= 0.05) {
  # X: numeric time X (one realization)

  T <- length(X)

  # Step 1: identify records
  records <- record_times(X)  # times of records
  N_T <- length(records)                      # number of records

  # Step 2: Fit log-log regression across subsamples
  # Split into blocks to see growth of R_t vs t
  block_sizes <- floor(seq(T/4, T, length.out = 5))  # subsample lengths
  R_block <- sapply(block_sizes, function(t) {
    length(record_values(X[1:t]))  ## number of records by blocks
  })

  df <- data.frame(logT = log(block_sizes), logR = log(R_block))

  fit <- lm(logR ~ logT, data = df)

  beta_hat <- coef(fit)[2]
  se_beta  <- summary(fit)$coefficients[2,2]

  # Wald test: H0: beta = 0.5
  z <- (beta_hat - 0.5) / se_beta
  pval <- 2 * (1 - pnorm(abs(z)))

  return(list(
    regression = summary(fit),
    beta_hat = beta_hat,
    se = se_beta,
    z = z,
    pval = pval,
    decision = ifelse(pval < alpha, "NO", "DTRW")
  ))
}

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
#     decision = ifelse(T_p <= -5.79 | T_p > 0.922, "NO", "DTRW")
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
#     decision = ifelse(T_p <= -5.79 | T_p > 0.922, "NO", "DTRW") ## Null is DTRW
#   ))
# }
