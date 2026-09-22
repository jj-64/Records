## --------  LDM_Regression Hubert ------------------
####Robust Trend Test for the Linear Drift Model
# #'
# #' Tests whether the observed sequence \code{X} follows the Linear Drift Model (LDM)
# #' by fitting a robust linear regression of \code{X} on time and testing
# #' if the slope (drift parameter) is significantly different from zero.
# #'
# #' @param X Numeric vector of observations.
# #' @param alpha Numeric, significance level (default = 0.05).
# #' @param RSq Numeric, minimum adjusted R-squared required to accept the LDM
# #'   hypothesis (default = 0.8).
# #' @param obs_type String. "all" if data provided is the whole series \eqn{X_t} or
# #' "records" if the underlying series is \eqn{R_n}. In this case, the parameter
# #' record_times must be provided.
# #' @param record_times Numeric vector of the occurence times of records. (Default is NA).
# #' Forced in case "obs_type" = "records"
# #' @return A list with components:
# #' \describe{
# #'   \item{stat}{Estimated slope coefficient (drift parameter).}
# #'   \item{p_value}{p_value for testing \eqn{H_0: \theta = 0}.}
# #'   \item{RS}{Adjusted R-squared of the robust regression.}
# #'   and adjusted R-squared exceeds the threshold, \code{"no"} otherwise.}
# #'   \item{MSE}{Mean squared error of the fitted model (normalized by \eqn{\sum X^2}).}
# #'   \item{STD}{Estimated standard error of the slope coefficient.}
# #'   \item{decision}{Decision: \code{"ldm"} if the null hypothesis is rejected
# #' }
# #'
# #' @details
# #' The test fits a robust regression model of the form:
# #' \deqn{ X_t = \theta (t - \bar{t}) + \varepsilon_t }
# #' where time \eqn{t} and observations \eqn{X_t} are centralized.
# #' The decision rule is:
# #' \enumerate{
# #'   \item Reject \eqn{H_0: \theta = 0} if the p_value is below \code{alpha}.
# #'   \item Accept LDM only if the adjusted R-squared exceeds \code{RSq}.
# #' }
# #'
# #' @examples
# #' set.seed(123)
# #' t <- 1:50
# #' X <- 0.3 * t + rnorm(50, sd = 5)  # Linear drift with noise
# #' test_ldm_trend(X, alpha = 0.05, RSq = 0.7)
# #' @export


#' Robust Trend Test for the Linear Drift Model
#'
#' Tests whether an observed time series exhibits a significant
#' linear drift using robust regression. The procedure estimates
#' the drift parameter and evaluates whether it differs
#' significantly from zero.
#'
#' @details
#' The Linear Drift Model (LDM) assumes that observations follow
#'
#' \deqn{
#' X_t = \theta t + \varepsilon_t,
#' }
#'
#' where \eqn{\theta} denotes a deterministic drift parameter and
#' \eqn{\varepsilon_t} is a random error process.
#'
#' A robust MM-estimator (\code{robustbase::lmrob}) is used to fit
#' the model and reduce the influence of outliers.
#'
#' The null and alternative hypotheses are:
#'
#' \deqn{
#' H_0 : \theta = 0
#' }
#'
#' versus
#'
#' \deqn{
#' H_A : \theta \neq 0.
#' }
#'
#' A Wald-type statistic is computed as
#'
#' \deqn{
#' Z =
#' \frac{\hat{\theta}}
#' {\mathrm{SE}(\hat{\theta})}.
#' }
#'
#' The p-value is obtained from the corresponding t-distribution.
#'
#' In addition to testing significance, the function reports the
#' normalized mean squared error (NMSE):
#'
#' \deqn{
#' \mathrm{NMSE}
#' =
#' \frac{\sum_t (X_t-\hat X_t)^2}
#' {\sum_t (X_t-\bar X)^2}.
#' }
#'
#' NMSE measures the proportion of variability not explained by the
#' estimated linear drift. Smaller values indicate a better fit to
#' the Linear Drift Model.
#'
#' \enumerate{
#' \item {NMSE < 0.25}{ : strong support}
#' \item {NMSE < 0.50}{ : moderate support}
#' \item {NMSE > 0.75}{ : weak practical support}
#' }
#' @param X Numeric vector of observations.
#'
#' @param alpha Significance level.
#' Default is \code{0.05}.
#'
#' @param obs_type Character string.
#' Either \code{"all"} when \code{X} contains the complete time
#' series, or \code{"records"} when \code{X} contains only record
#' values.
#'
#' @param record_times Numeric vector containing the occurrence times
#' of the records when \code{obs_type = "records"}.
#'
#' @return A list containing:
#'
#' \describe{
#'
#' \item{theta_hat}{
#' Estimated drift parameter.
#' }
#'
#' \item{std_error}{
#' Standard error of the drift estimate.
#' }
#'
#' \item{statistic}{
#' Wald t-statistic for testing \eqn{\theta=0}.
#' }
#'
#' \item{p_value}{
#' Two-sided p-value.
#' }
#'
#' \item{adj_r_squared}{
#' Adjusted coefficient of determination.
#' }
#'
#' \item{nmse}{
#' Normalized mean squared error.
#' Values near zero indicate a strong linear drift component.
#' }
#'
#' \item{decision}{
#' \code{"LDM"} if the null hypothesis is rejected at level
#' \code{alpha}; otherwise \code{"NO"}.
#' }
#'
#' }
#'
#' @examples
#' set.seed(123)
#' t <- 1:100
#' X <- 0.2 * t + rnorm(100)
#'
#' # If we assume we have the whole series
#' test_ldm_trend(X)
#'
#' # $theta_hat
#' #  [1] 0.1897423
#'
#' #  $std_error
#' #  [1] 0.01246424
#'
#' #  $statistic
#' #  [1] 15.22293
#'
#' #  $p_value
#' #  [1] 5.320148e-20
#'
#' #  $conf_int
#' #  [1] 0.1646812 0.2148033
#'
#' #  $adj_r_squared
#' #  [1] 0.8538915
#'
#' #  $nmse
#' #  [1] 0.1362929
#'
#' #  $decision
#' #  [1] "ldm"
#'
#' # If we assume we have record values and record times only
#' test_ldm_trend(rec_values(X), obs_type = "records", record_times = rec_times(X))
#'
#' #' $theta_hat
#' # [1] 0.2156731
#'
#' # $std_error
#' # [1] 0.008436433
#'
#' # $statistic
#' # [ 1] 25.56449
#'
#' # $p_value
#' # [1] 8.774396e-14
#'
#' # $conf_int
#' # [1] 0.1976913 0.2336550
#'
#' # $adj_r_squared
#' # [1] 0.9752673
#'
#' # $nmse
#' # [1] 0.0206316
#'
#' # $decision
#' # [1] "ldm"
#' #'
#' @export
test_ldm_trend <- function(
    X,
    alpha = 0.05,
    obs_type = c("all", "records"),
    record_times = NULL
  ) {

  obs_type <- match.arg(obs_type)

  ##Force to provide record times in case only records are used
  if (obs_type == "records") {

    if (is.null(record_times))
      stop("'record_times' must be supplied.")

    if (length(record_times) != length(X))
      stop("'record_times' and 'X' must have the same length.")

    t <- record_times

  } else {

    t <- seq_along(X)

  }

  # # Centralize data
  # t <- t - mean(t)
  # X <- X - mean(X)
  dat <- data.frame(
    X = X,
    t = t
  )

  fit <- robustbase::lmrob(
    X ~ t,
    data = dat,
    method = "MM"
  )

  # Fit robust regression without intercept
  #train <- data.frame(X = X, t = t, rec= ifelse(is_rec(X) ==1, rec_times(X), 0))
  #robust_model <- robustbase::lmrob(X ~ t - 1, data = train, method = "MM")

  # Extract coefficient and standard error
  sm <- summary(fit)

  theta_hat <- coef(fit)[["t"]]

  se_theta <- sm$coefficients["t", "Std. Error"]

  t_stat <- theta_hat / se_theta

  p_value <- 2 * pt(
    -abs(t_stat),
    df = fit$df.residual
  )

  ci <- theta_hat +
    c(-1, 1) *
    qt(1 - alpha/2, fit$df.residual) *
    se_theta

  rsq <- sm$adj.r.squared

  nmse <- sum(residuals(fit)^2) /
    sum((X - mean(X))^2)

  decision <- ifelse(
    p_value < alpha,
    "ldm",
    "no"
  )

  return(list(
    theta_hat = theta_hat,
    std_error = se_theta,
    statistic = t_stat,
    p_value = p_value,
    conf_int = ci,
    adj_r_squared = rsq,
    nmse = nmse,
    decision = decision
  ))

  # coef_summary <- coef(summary(robust_model))
  # coefficient <- robust_model$coefficients
  # std_error <- coef_summary[1, 2]

  # # Compute p_value for H0: theta = 0
  # t_stat <- coef_summary[1, 1] / std_error
  # p_value <- 1 - pt(t_stat, df = robust_model$df.residual)
  #
  # # Compute fit quality
  # MSE <- sum(robust_model$residuals)^2 / sum(X^2)
  # RS <- summary(robust_model)$adj.r.squared

  # Decision rule
  #decision <- ifelse(p_value < alpha & RS >= RSq, "ldm", "no")
  # decision <- ifelse(p_value < alpha & RS >= RSq, "ldm", "no")
  #
  # return(list(
  #   stat = coefficient,
  #   p_value = p_value,
  #   RS = RS,
  #   MSE = MSE,
  #   STD = std_error,
  #   decision = decision
  # ))
}

## ------- Test based on observed number of records ------------------ ################

# #' Quantile Function for LDM Distribution of Record Numbers in case of Gumbel underlying distribution
# #'
# #' Computes the lower and upper quantile bounds for the distribution
# #' of the number of records under the Linear Drift Model (LDM) under Gumble (µ, σ).
# #'
# #' @param T Integer, sample size (number of observations).
# #' @param theta Numeric, drift parameter \eqn{\theta > 0}.
# #' @param scale Numeric, scale parameter σ (default = 1).
# #' @param alpha Numeric, significance level for the two-sided interval (default = 0.05).
# #'
# #' @return A numeric vector of length 2 giving the lower and upper quantile indices
# #' corresponding to probabilities \eqn{\alpha/2} and \eqn{1-\alpha/2}.
# #'
# #' @details
# #' The function computes the distribution of the number of records under the LDM,
# #' using the Stirling numbers of the second kind (\code{Stirling_2nd_LDM}) and
# #' the probability mass function \code{rec_count_dist_LDM}.
# #' The cumulative distribution function (CDF) is then compared to the
# #' desired quantile levels.
# #'
# #' @seealso \code{\link{test_ldm_rec_count}} for hypothesis testing of the number of records.
# #' @export
# #' @examples
# #' Quantile_LDM(T = 20, theta = 0.5, alpha = 0.05)
# #' # Quantile_LDM <- function(T, theta, scale = 1, alpha = 0.05) {#   S <- Stirling_2nd_LDM(T, theta, scale)##   # Vectorized probability computation#   Prob <- vapply(1:T, function(i) {#     rec_count_dist_LDM(m = i, T = T, theta = theta, scale = scale, s = S)#   }, numeric(1))##   # Cumulative distribution#   CDF <- cumsum(Prob)##   # Return indices closest to alpha/2 and 1 - alpha/2#   return(c(#     which.min(abs(CDF - alpha / 2)),#     which.min(abs(CDF - (1 - alpha / 2)))#   ))#

#' Hypothesis Test for Number of Records under LDM and Gumbel underlying distribution
#'
#' Record-Count Test for the Linear Drift Model (LDM)
#'
#' Tests whether an observed time series is consistent with a
#' Linear Drift Model (LDM) using a two-stage procedure based on
#' drift estimation and record statistics.
#'
#' @details
#' The Linear Drift Model (LDM) assumes that observations can be
#' represented as
#'
#' \deqn{
#' X_t = Y_t + \theta t,
#' }
#'
#' where \eqn{Y_t} are independent observations drawn from a
#' Gumbel distribution and \eqn{\theta} denotes a linear drift
#' parameter.
#'
#' The procedure consists of two stages:
#'
#' \enumerate{
#'
#' \item Drift detection.
#'
#' The drift parameter \eqn{\theta} is estimated using
#' \code{estimate_LDM_mle_indicator()}.
#'
#' A Wald-type statistic is then computed:
#'
#' \deqn{
#' Z_\theta=
#' \frac{\hat\theta}
#' {\sqrt{Var(\hat\theta)}}.
#' }
#'
#' If the estimated drift is not significantly different from zero,
#' the LDM hypothesis is rejected.
#'
#' \item Record-count consistency.
#'
#' Conditional on detecting a significant drift, the observed number
#' of records is compared with the theoretical acceptance region of
#' the LDM obtained from
#' \code{\link{rec_count_bounds}}.
#'
#' The LDM hypothesis is accepted only if the observed record count
#' falls inside the corresponding acceptance interval.
#'
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
#' \item{record_count}{
#' Observed number of records.
#' }
#'
#' \item{theta_hat}{
#' Estimated drift parameter.
#' }
#'
#' \item{var_theta_hat}{
#' Estimated variance of the drift estimator.
#' }
#'
#' \item{z_theo}{
#' Wald statistic for testing
#' \eqn{\theta=0}.
#' }
#'
#' \item{critical_bounds}{
#' Lower and upper record-count acceptance bounds.
#' }
#'
#' \item{decision}{
#' \code{"LDM"} if the series is consistent with the Linear Drift
#' Model and \code{"NO"} otherwise.
#' }
#'
#' }
#'
#' @seealso
#' \code{\link{estimate_LDM_mle_indicator}},
#' \code{\link{rec_count_bounds}}
#'
#' @examples
#' set.seed(123)
#' X <- rnorm(50, mean = 0.2 * (1:50))  # Linear drift
#'
#' test_ldm_rec_count(X, alpha = 0.05)
#' @export
test_ldm_rec_count <-  function(X, alpha = 0.05) {

    obs <- rec_count(X)

    estimate <- estimate_LDM_mle_indicator(
      X = X,
      variance = TRUE,
      scale = 1,
      min = 0.0001,
      max = 5,
      step = 0.01
    )

    theta_hat <- estimate$param
    var_theta_hat <- estimate$variance

    z_theo <- theta_hat / sqrt(var_theta_hat)

    drift_detected <-
      abs(z_theo) > qnorm(1 - alpha/2)


    # Decision rule
    if (!drift_detected) {

      decision <- "no"

      bounds <- NULL

    } else {

      bounds <- rec_count_bounds(
        T = length(X),
        model = "ldm",
        theta = theta_hat,
        scale = 1,
        alpha = alpha
      )

      decision <- ifelse(
        obs >= bounds$lower_bound &&
          obs <= bounds$upper_bound,
        "ldm",
        "no"
      )
    }

    list(
      record_count = obs,
      theta_hat = theta_hat,
      var_theta_hat = var_theta_hat,
      z_theo = z_theo,
      critical_bounds = c("lower_bound" = bounds$lower_bound, "upper_bound" =bounds$upper_bound ),
      decision = decision
    )
  }

## ----------- Two-stage based on Şen Slope  ---------------- ##########

#' Sequential Drift Detection Test for the Linear Drift Model
#'
#' Implements a two-stage sequential testing procedure (based on Şen-type tests)
#' to detect the presence of a linear drift in a time series.
#'
#' The procedure partitions the data into blocks and checks whether block means
#' and slopes are consistent with a monotonic non-zero linear drift. If stage 1 accepts, stage 2
#' refines the decision by testing slope significance across partitions.
#' In other words, the procedure first evaluates whether the series is compatible
#' with a constant linear drift and subsequently tests whether the
#' drift is significantly different from zero.
#'
#' @details
#'
#' The series is partitioned into three consecutive blocks of equal
#' size.
#'
#' Stage 1 evaluates the linearity of the trend using
#'
#' \deqn{
#' Z = 2\bar X_2-\bar X_1-\bar X_3.
#' }
#'
#' Under a Linear Drift Model,
#'
#' \deqn{
#' E(Z)=0,
#' }
#'
#' since the middle block mean should lie approximately halfway
#' between the first and third block means.
#'
#' Stage 2 estimates the drift separately within adjacent blocks and
#' computes standardized slope statistics
#'
#' \deqn{
#' z_A=\hat\theta_1/\hat\sigma_\theta
#' }
#'
#' and
#'
#' \deqn{
#' z_B=\hat\theta_2/\hat\sigma_\theta.
#' }
#'
#' The LDM hypothesis is supported only when:
#'
#' \enumerate{
#'   \item Stage 1 does not reject linearity;
#'   \item both slope statistics are significant;
#'   \item both slopes have the same sign.
#' }
#'
#' This combination provides evidence for a persistent positive or
#' negative linear drift over the entire observation period.
#'
#' @param X Numeric vector of observations.
#' @param time Numeric vector of the record times. (default = NA means we are using all observed X series)
#' @param alpha Numeric, significance level (default = 0.05).
#' @param pooled Logical, if \code{TRUE}, pooled variance across blocks is used
#'   in stage 1 for the noise variance estimation. Otherwise, variance is estimated directly
#'   from residuals (default = \code{FALSE}).
#' @return A list with results from both stages:
#' \describe{
#'   \item{block_size}{Size of each partitioned block.}
#'   \item{block_means}{Block means.}
#'   \item{block_slopes}{Estimated block-to-block slopes.}
#'   \item{var_noise}{Estimated variance of the noise process.}
#'   \item{Z}{Test statistic from stage 1.}
#'   \item{Var_Z}{Estimated variance of \code{Z}.}
#'   \item{stage1_statistic}{Standardized stage 1 test statistic.}
#'   \item{stage1_p_value}{p_value of the stage 1 test.}
#'   \item{stage1_decision}{Decision of stage 1: \code{"Sameslope"} or \code{"no"}.}
#'   \item{stage2_z1}{Stage 2 standardized slope for first block comparison.}
#'   \item{stage2_z2}{Stage 2 standardized slope for third block comparison.}
#'   \item {drift_direction}{sign of the average value of both slopes.}
#'   \item{drift_sd}{Estimated standard deviation of slope differences (stage 2).}
#'   \item{decision}{Final decision: \code{"ldm"} or \code{"no"}.}
#' }
#'
#' @details
#' Stage 1 partitions the series into 3 blocks of size \eqn{m = \lfloor n/3 \rfloor}.
#' It computes block means and slopes and tests whether the central block mean
#' aligns with the average of the first and third block means:
#' \deqn{ Z = \frac{2 \bar{X}_2 - \bar{X}_1 - \bar{X}_3}{m} }
#' with variance estimated either via pooling (\code{pooled=TRUE}) or directly.
#'
#' Stage 2 tests whether the slopes between blocks are significantly different
#' from zero:
#' \deqn{ z_A = \frac{s_1}{\hat{\sigma}_s}, \quad z_B = \frac{s_2}{\hat{\sigma}_s} }
#' where \eqn{\hat{\sigma}_s = \sqrt{2 \hat{\sigma}_y / m^3}} and \eqn{Y} is the underlying distribution of the LDM process.
#'
#' The procedure stops at stage 1 if no drift is detected.
#'
#' @examples
#' set.seed(123)
#' t <- 1:60
#' X <- 0.1 * t + rnorm(60, sd = 1)
#' test_ldm_sequential(X, alpha = 0.05)
#'
#' # $block_size
#' #  [1] 16
#'
#' #  $block_means
#' #  [1] 1.681402 4.984581 8.185253
#'
#' #  $block_slopes
#' #  [1] 0.2064487 0.2000420
#'
#' #  $var_noise
#' #  [1] 1.238635
#'
#' #  $Z
#' #  [1] 0.1025064
#'
#' #  $Var_Z
#' #  [1] 0.02903052
#'
#' #  $stage1_statistic
#' #  [1] 0.6016217
#'
#' #  $stage1_p_value
#' #  [1] 0.273713
#'
#' #  $stage1_decision
#' #  [1] "Sameslope"
#'
#' #  $stage2_z1
#' #  [1] 8.394703
#'
#' #  $stage2_z2
#' #  [1] 8.134193
#'
#' #  $drift_direction
#' #  [1] 1
#'
#' #  $drift_sd
#' #  [1] 0.02459273
#'
#' #  $decision
#' #  [1] "ldm"
#'
#' test_ldm_sequential(rnorm(100), alpha = 0.05)
#'
#' # $block_size
#' #  [1] 33
#'
#' #  $block_means
#' #  [1] 0.08574701 0.11501847 0.10558325
#'
#' #  $block_slopes
#' #  [1]  0.0008870138 -0.0002859157
#'
#' #  $var_noise
#' #  [1] 1.161635
#'
#' #  $Z
#' #  [1] 0.03870667
#'
#' #  $Var_Z
#' #  [1] 0.006400193
#'
#' #  $stage1_statistic
#' #  [1] 0.4838261
#'
#' #  $stage1_p_value
#' #  [1] 0.3142546
#'
#' #  $stage1_decision
#' #  [1] "Sameslope"
#'
#' #  $stage2_z1
#' #  [1] 0.1103193
#'
#' #  $stage2_z2
#' #  [1] -0.03555978
#'
#' #  $drift_direction
#' #  [1] 1
#'
#' #  $drift_sd
#' #  [1] 0.008040424
#'
#' #  $decision
#' #  [1] "no"
#'
#' @export
test_ldm_sequential <- function(X, time = NA,alpha = 0.05) {
  if(length(X)<4) {return(decision = "no")}

  # Stage 1
  res1 <- test_ldm_sequential_linearity_stage1(X, time=time, alpha = 2 * alpha, pooled = FALSE)  # one-sided

  if (res1$stage1_decision == "no") {
    #names(res1)[which(names(res1) == 'stage1_decision')] = "decision"
    res1$decision = "no"
    return(res1)
  }

  # Stage 2
  res2 <- test_ldm_sequential_drift_detection_stage2(res1, alpha = alpha)

  return(c(res1, res2))
}

#' @rdname test_ldm_sequential
#' @param pooled boolean (Default = FALSE), if variance is assumed pooled
#' @export
test_ldm_sequential_linearity_stage1 <- function(X, time=NA, alpha = 0.05, pooled = FALSE) {

  if(is.na(time[1])) {time = 1:length(X)}
  if (length(X) != length(time)) stop("X and time must have the same length")

  #n <- length(X)
  n_blocks <- 3

  # Sort by time
  ord <- order(time)
  X <- X[ord]
  time <- time[ord]

  block_means <- numeric(n_blocks)
  blocks <- vector("list", n_blocks)

  # Partition data into equal blocks
  m <- floor(length(X) / n_blocks)
  for (i in 1:n_blocks) {
    start <- (i - 1) * m + 1
    end <- ifelse( i< n_blocks, i * m, length(X) )
    blocks[[i]] <- X[start:end]
    block_means[i] <- mean(blocks[[i]])}


  # Slopes between consecutive block means
  slopes <- numeric(n_blocks)
  for (i in 1:(n_blocks - 1)) {
    slopes[i + 1] <- (block_means[i + 1] - block_means[i]) / m
  }

  # Residual estimates (noise)
  x_hat <- list(
    blocks[[1]] - slopes[2] * seq_len(m),
    blocks[[2]] - slopes[2] * ((m + 1):(2 * m)),
    blocks[[3]] - slopes[3] * (((2 * m) + 1):(length(X)))
  )

  # Stage 1 statistic (equality of slopes)

  z <- (2 * block_means[2] - block_means[1] - block_means[3])

  # Variance of Z
  if (pooled) {
    s_sq <- sapply(x_hat, var)
    var_noise <- (m - 1) * sum(s_sq) / (3 * m - 3)
  } else {
    var_noise <- var(unlist(x_hat))
  }

  # Standardized test statistic
  var_z <- (6 * var_noise) / m^2
  T_p <- z / sqrt(var_z)

  # p_value
  p_value <- if (max(time) >= 30) {
    #1 - pnorm(abs(T_p))
    2 * pnorm(-abs(T_p))
  } else {
    #1 - pt(q = abs(T_p), df = length(X) - 3)
    2 * pt(-abs(T_p), df = length(X) - 3)
  }

  # Decision: same slope or reject
  decision <- ifelse(p_value > alpha / 2, "Sameslope", "no")

  return(list(
    block_size = m,
    #block_times = block_times,
    block_means = block_means,
    block_slopes = slopes[-1],
    var_noise = var_noise,
    Z = z,
    Var_Z = var_z,
    stage1_statistic = T_p,
    stage1_p_value = p_value,
    stage1_decision = decision
  ))
}

#' @rdname test_ldm_sequential
#' @param stage1 output of \code{\link{test_ldm_sequential_linearity_stage1}}
#' @export
test_ldm_sequential_drift_detection_stage2 <- function(stage1, alpha = 0.05) {
  m <- stage1$block_size
  s1 <- stage1$block_slopes[1]
  s2 <- stage1$block_slopes[2]
  var_noise <- stage1$var_noise

  # Standard deviation of slope estimator
  sigma_s <- sqrt(2 * var_noise / m^3)

  # Test statistics
  z_A <- s1 / sigma_s
  z_B <- s2 / sigma_s

  # Decision
  # decision <- ifelse(abs(z_A) <= qnorm(1 - alpha / 2) | abs(z_B) <= qnorm(1 - alpha / 2),
  #               "no", "ldm")
  crit =  qnorm(1 - alpha / 2)

  same_sign <- sign(s1) == sign(s2)

  significant <- abs(z_A )> crit & abs(z_B) > crit

  decision <- ifelse(
    significant & same_sign,
    "ldm",
    "no"
  )

  return(list(stage2_z1 = z_A,
              stage2_z2 = z_B,
              drift_direction = sign(mean(c(s1,s2))),
              drift_sd = sigma_s,
              decision = decision))
}


##### ------------------  Based on Sen-Slope ------------------ #############

Test_LDM_Sen = function(X, alpha=0.05){

  n= length(X)

  m <- floor(n/2)
  y1 <- X[1:m]
  y2 <- X[(m+1):(2*m)]
  if (length(y1) != length(y2)) y2 = y2[1:length(y1)]

  # block means
  y1_bar <- mean(y1)
  y2_bar <- mean(y2)

  # slope and intercept
  s <- (y2_bar - y1_bar) / m
  a <- mean(X) - s * mean(1:n)

  ## Fit
  fit <- a + (1:n) * s
  noise <- X - fit

  # correlation:order-stat correlation
  rho <- cor(sort(y1), sort(y2))     # paired correlation between blocks: closer to one, then no trend bcz on 1:1 line

  ## Noise estimation
  sigma_s = sqrt(8 / n^3 * var(X)) #* sqrt((1 - abs(rho)))
  #sigma_s_theo = 8* var(x_hat)* (2)/(n^3)
  #CL_low_H0 = qnorm(alpha/2)* sigma_s  #-1.959 for alpha =0.05
  #CL_up_H0 = qnorm(1-alpha/2)* sigma_s
  #"CL_H0" = c(CL_low_H0, CL_up_H0)
  #decision= ifelse(s >= CL_up_H0 | s<= CL_low_H0, "ldm","no")
  z=s/sigma_s
  decision = ifelse(abs(z) >= qnorm(1-alpha/2), "ldm", "no")

   return(list( "stat" = z,"theta_hat"=s,"theta_SD"=sigma_s,"decision"=decision))
  #abline(a=a, b=s)
  #plot(sort(y2),sort(y1))
  #abline(a=0, b=1)
            }

# test_ldm_sequential_linearity_stage1 <- function(X, alpha=0.05, pooled = FALSE) {
#   n <- length(X)
#   n_blocks <- 3    #floor(length(X) / m)
#   m = floor(n/n_blocks)
#
#   rho_all = numeric(n_blocks)
#   block_means <- numeric(n_blocks)
#
#   slopes <- 0
#   var_slopes = 0
#   blocks= list()
#   x_hat = list()
#
#   var_slope = function(n,m, sigma2, correlation) {2/3 * sigma2/m^3 * (1-correlation) }
#
#   for (i in 1:n_blocks) {
#     start <- (i - 1) * m + 1
#     end <- i * m
#     block <- X[start:end]
#
#     # Store block mean
#     block_means[i] <- mean(block)
#     blocks[[i]] = block
#   }
#
#   ## slopes
#   for (i in 1:(n_blocks-1)){
#     slopes[i+1] = (block_means[i+1] - block_means[i])/m
#   }
#
#   rho_all[1] =  cor(sort(blocks[[1]]), sort(blocks[[2]])  )  # 1 and 2
#   rho_all[2] =  cor(sort(blocks[[2]]), sort(blocks[[3]])  )  # 2 and 3
#   rho_all[3] =  cor(sort(blocks[[1]]), sort(blocks[[3]])  )  # 1 and 3
#
#   ## Estimation of x_hat, the noise
#   x_hat[[1]] = blocks[[1]] - slopes[2] * c(1:m)
#   x_hat[[2]] = blocks[[2]] - slopes[2] * c((m+1):(2*m))
#   x_hat[[3]] = blocks[[3]] - slopes[3] * c(((2*m)+1):(3*m))
#
#
#   # Results
#   w <- data.frame(
#     block = 1:n_blocks,
#     mean_value = block_means,
#     slopes = slopes,
#     var_slopes = var_slopes
#   )
#
#   # 3. Compute z
#   z = (2*w[2,2]-w[1,2]-w[3,2])/m
#
#   # 4. Estimate variance of z, xhat and V(x)
#
#   if(pooled){
#     # Option A: pooled variance across all three groups (assuming equal variances)
#       s1_sq <- var(x_hat[[1]]) #var(blocks[[1]])
#       s2_sq <- var(x_hat[[2]]) #var(blocks[[2]])
#       s3_sq <- var(x_hat[[3]])# var(blocks[[3]])
#
#       var_noise <- (m-1) * (s1_sq + s2_sq + s3_sq) / (3*m - 3)  ## variance of X
#       var_z <- (6 *var_noise) / m^3
#   # } else if (unpooled) {
#   #   # Option B: group-wise (if you don’t assume equal variances)
#   #     var_z <- (4*s2_sq + s1_sq + s3_sq) / m^3
#   } else{
#     var_noise = var(unlist(x_hat))
#     var_z =  6 * var_noise/m^3 #2 * (var(x_hat)/m^3) * (3- 2* rho_all[1] - 2*rho_all[2] + rho_all[3])  #
#   }
#
#   ## estimated variance of slope
#   #var_slopes[2] = var_slope(n=n, m=m, sigma2 = var_noise, correlation = rho_all[1])
#   #var_slopes[3] = var_slope(n=n, m=m, sigma2 = var_noise, correlation = rho_all[2])
#
#   #5.Compute test statistic
#   T_p <- z / sqrt(var_z)
#
#   #p_value <- 2 * pt(-abs(T_p), df = n - 3)
#   if(n >=30){
#   p_value = 1-pnorm(abs(T_p))
#   }else{
#   p_value = 1-pt(q=abs(T_p), df = n-3)
#   }
#
#   #decision= ifelse((T_p <= qnorm(alpha,0,1) | T_p >= qnorm(1-alpha,0,1)), "no", "SameSlope")
#   decision = ifelse (p_value > alpha/2 , "Sameslope", "no")
#   #decision= ifelse( z >= qnorm(alpha/2,0,1)*sqrt(var_z) &  z <=  qnorm(1-alpha/2,0,1)*sqrt(var_z), "ldm", "no")
#   return(list("block_size"= m,"means"= w$mean_value, "slopes"= w$slopes[-1],"var_noise" = var_noise,"Z"=z, "Var_Z" = var_z, "stat"=  T_p,"p_value" =  p_value, "decision" = decision))
# }
#
# test_ldm_sequential_drift_detection_stage2 = function(stage1, alpha=0.05){
#
#   ## as if testing Y1-bar = y3_bar or if slope from 1st and third partition is zero
#   m = stage1$block_size
#
#   ## slopes
#   s1 = stage1$slopes[1]
#   s2 = stage1$slopes[2]
#   #s_hat = mean(s1, s2)
#
#   ## var (s-hat)
#   var_noise = stage1$var_noise
#   sigma_s = sqrt(2*var_noise / m^3)
#
#   ## Test
#   z_A=  s1/sigma_s
#   z_B=  s2/sigma_s
#   decision = ifelse(abs(z_A) <= qnorm(1-alpha/2) | abs(z_B) <= qnorm(1-alpha/2),  "no", "ldm")
#
#   return(list( "z_A" = z_A, "z_B" = z_B, "drift_sd"=sigma_s,"decision"=decision))
# }
#
# test_ldm_sequential <- function(X, alpha = 0.05) {
#   t= seq_along(X)
#   m = floor(length(X)/3)
#
#   # Step 1: Run test_ldm_sequential_linearity_stage1
#   res1 <- test_ldm_sequential_linearity_stage1(X, alpha = 2*alpha, pooled = FALSE)  ## one sided test
#
#   if (res1$decision == "no") {
#     # Stop immediately
#     return(res1)
#   }
#
#   # Step 2: Only run if stage1 says "ldm"
#   n= floor(length(X)/2)
#   # res2_A <- Test_LDM_Sen(X[1:n], alpha = alpha)
#   # res2_B <- Test_LDM_Sen(X[(n +1):length(X)], alpha = alpha)
#   # res2 = ifelse(res2_A$decision == "no" |  res2_B$decision == "no", "no", "ldm")
#   res2 = test_ldm_sequential_drift_detection_stage2(res1, alpha=alpha)
#
#   return( c(res1,res2 ))
# }

