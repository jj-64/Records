# | Feature block          | Detects              | Useful under          |
#   | ---------------------- | -------------------- | --------------------- |
#   | Record rates (fwd/bwd) | Drift direction      | RW vs drift           |
#   | Log-record slope       | Record law deviation | Non-i.i.d.            |
#   | Gap CV / max gap       | Burstiness           | Long memory           |
#   | Entropy                | Regularity           | Deterministic trend   |
#   | Spectral peak          | Cycles               | Periodic alternatives |
#   | Rolling instability    | Local breaks         | Regime switching      |
#   | Min/max symmetry       | Skewed innovation    | Non-Gaussian noise    |

FEATURE_VERSION <- "v1.0-paper3b"
FEATURE_DATE <- "2025-12-15"

# ----- 0. Rolling instability & fast slope -----
#' Rolling Instability
#'
#' A generic function for applying a series of functions to rolling margins
#' of a vector. The functions are mean, variance, and sd-to-mean ratio.
#'
#' @param x the data to be used (representing a series of observations).
#' @param width numeric value, an integer specifying the window width (in numbers of observations).
#' @return A named vector the results of the rolling functions.
#' @export
#' @examples
#' \dontrun{
#' rolling_instability(x=c(1,2,3,4,5) , width = 2)
#' #     mean_instab    var_instab mean_sd_ratio
#' #      1.290994      0.000000      1.825742
#' }
rolling_instability <- function(x, width = 20) {
  n <- length(x)
  if (n < 2 * width) return(c(NA, NA, NA))

  rolls <- zoo::rollapply(
    x, width = width,
    FUN = function(x) c(mean(x), sd(x)),
    by.column = FALSE, align = "right"
  )

  mean_sd_ratio <- sd(rolls[,1]) / mean(rolls[,2])
  var_instab <- sd(rolls[,2])
  mean_instab <- sd(rolls[,1])

  c(mean_instab = mean_instab, var_instab = var_instab, mean_sd_ratio = mean_sd_ratio)
}

fast_slope <- function(s) {
  n <- length(s)
  t <- seq_len(n)
  return(cov(s, t) / var(t))
}

# ----- 1. Simulate time series under multiple record models --------------------

## ---- Helper function: Generate many series and keep labels, returns a numeric vector.

generate_series <- function(generator, ## function:the function generating the series
                            series_args=list(), ## arguments of the generator function other than "T" and the "param_name" we are simulating
                            n_arg="T",
                            T_val = 200) {
  ## --- Helper: Generate series under true model H0
  args <- series_args
  args[[n_arg]] <- T_val   # could be "T" or "n"

  X <- do.call(generator, args)

  return(X)
}


## ---- Helper function: Generate many series and keep labels
#n_per_model: number of series generated for each model to have a balanced database
generate_series_multiple <- function(
    n_per_model = 50,
    T_vals = c(100, 200, 500) ,
    normalized = TRUE) {
  all_series <- list()
  labels <- c()
  labels_m <- c()
  series_id <- c()
  Ts <- c()  # store T for each series

  for (T_val in T_vals) {

    ## ---- DTRW
    i = 1
    while(i <= n_per_model) {
      s <- generate_series(
        DTRW_series,
        series_args = list(dist = "norm", mean = 0, sd = 1),
        T_val = T_val
      )
      if(length(rec_gaps(s)) <2 ) next;
      all_series[[length(all_series) + 1]] <- if(normalized) {s/max(s) } else {s}
      labels <- c(labels, "DTRW")
      series_id <- c(series_id, paste0("DTRW_T",T_val,"_",i))
      Ts <- c(Ts, T_val)
      labels_m = c(labels_m, "norm")
      i = i + 1
    }

    i = 1
    while(i <= n_per_model) {
      s <- generate_series(
        DTRW_series,
        series_args = list(dist = "cauchy", mean = 0, scale = 1),
        T_val = T_val
      )
      if(length(rec_gaps(s)) <2 ) next;
      all_series[[length(all_series) + 1]] <- if(normalized) {s/max(s) } else {s}
      labels <- c(labels, "DTRW")
      series_id <- c(series_id, paste0("DTRW_T",T_val,"_",i))
      Ts <- c(Ts, T_val)
      labels_m = c(labels_m, "cauchy")
      i = i + 1
    }

    i = 1
    while(i <= n_per_model) {
      s <- generate_series(
        DTRW_series,
        series_args = list(dist = "uniform", min = -1, scale = 1),
        T_val = T_val
      )
      if(length(rec_gaps(s)) <2 ) next;
      all_series[[length(all_series) + 1]] <- if(normalized) {s/max(s) } else {s}
      labels <- c(labels, "DTRW")
      series_id <- c(series_id, paste0("DTRW_T",T_val,"_",i))
      Ts <- c(Ts, T_val)
      labels_m = c(labels_m, "uniform")
      i = i + 1
    }

    ## ---- LDM
        ## Frechet
    i = 1
    while(i <= n_per_model) {
      s <- generate_series(
        LDM_series,
        series_args = list(theta = runif(1,0.02,0.1),
                           dist = "frechet", shape=5, scale=1),
        T_val = T_val
      )
      if(length(rec_gaps(s)) <2 ) next;
      all_series[[length(all_series) + 1]] <- if(normalized) {s/max(s) } else {s}
      labels <- c(labels, "LDM")
      series_id <- c(series_id, paste0("LDM_T",T_val,"_",i))
      Ts <- c(Ts, T_val)
      labels_m = c(labels_m, "frechet")
      i = i + 1
    }
        ## Weibull
    i = 1
    while(i <= n_per_model) {
      s <- generate_series(
        LDM_series,
        series_args = list(theta = runif(1,0.02,0.1),
                           dist = "weibull", shape = 2, scale=1),
        T_val = T_val
      )
      if(length(rec_gaps(s)) <2 ) next;
      all_series[[length(all_series) + 1]] <- if(normalized) {s/max(s) } else {s}
      labels <- c(labels, "LDM")
      series_id <- c(series_id, paste0("LDM_T",T_val,"_",i))
      Ts <- c(Ts, T_val)
      labels_m = c(labels_m, "weibull")
      i = i + 1
    }
        ## Normal
    i = 1
    while(i <= n_per_model) {
      s <- generate_series(
        LDM_series,
        series_args = list(theta = runif(1,0.02,0.1),
                           dist = "norm", mean =0 , sd =1),
        T_val = T_val
      )
      if(length(rec_gaps(s)) <2 ) next;
      all_series[[length(all_series) + 1]] <- if(normalized) {s/max(s) } else {s}
      labels <- c(labels, "LDM")
      series_id <- c(series_id, paste0("LDM_T",T_val,"_",i))
      Ts <- c(Ts, T_val)
      labels_m = c(labels_m, "norm")
      i = i + 1
    }
    ## ---- YNM
        ## Frechet
    i=1
    while(i <= n_per_model) {
      s <- generate_series(
        YNM_series,
        series_args = list(gamma = runif(1,1.1,1.2),
                           dist = "frechet", shape=3, scale=1),
        T_val = T_val
      )
      if(length(rec_gaps(s)) <2 ) next;
      all_series[[length(all_series) + 1]] <- if(normalized) {s/max(s) } else {s}
      labels <- c(labels, "YNM")
      series_id <- c(series_id, paste0("YNM_T",T_val,"_",i))
      Ts <- c(Ts, T_val)
      labels_m = c(labels_m, "frechet")
      i=i+1
    }
        ## Weibull
    i=1
    while(i <= n_per_model) {
      s <- generate_series(
        YNM_series,
        series_args = list(gamma = runif(1,1.2,1.5),
                           dist = "weibull", shape= 1/2, scale =0.1),
        T_val = T_val
      )
      if(length(rec_gaps(s)) <2 ) next;
      all_series[[length(all_series) + 1]] <- if(normalized) {s/max(s) } else {s}
      labels <- c(labels, "YNM")
      series_id <- c(series_id, paste0("YNM_T",T_val,"_",i))
      Ts <- c(Ts, T_val)
      labels_m = c(labels_m, "weibull")
      i=i+1
    }
        ## Pareto truncated
    i=1
    while(i <= n_per_model) {
      s <- generate_series(
        YNM_series,
        series_args = list(gamma = runif(1,1.2,1.5),
                           dist = "pareto_trunc", shape= 2.8, scale =1, xmax = 1000),
        T_val = T_val
      )
      if(length(rec_gaps(s)) <2 ) next;
      all_series[[length(all_series) + 1]] <- if(normalized) {s/max(s) } else {s}
      labels <- c(labels, "YNM")
      series_id <- c(series_id, paste0("YNM_T",T_val,"_",i))
      Ts <- c(Ts, T_val)
      labels_m = c(labels_m, "pareto")
      i=i+1
    }

    ## ---- iid (Classical)
    i=1
    while(i <= n_per_model) {
      s <- VGAM::rfrechet(T_val, shape = 4, scale=1)
      if(length(rec_gaps(s)) <2 ) next;
      all_series[[length(all_series) + 1]] <- if(normalized) {s/max(s) } else {s}
      labels <- c(labels, "Classical")
      series_id <- c(series_id, paste0("Classical_T",T_val,"_",i))
      Ts <- c(Ts, T_val)
      labels_m = c(labels_m, "frechet")
      i=i+1
    }

    i=1
    while(i <= n_per_model) {
      s <- VGAM::rgumbel(T_val, 0, 1)
      if(length(rec_gaps(s)) <2 ) next;
      all_series[[length(all_series) + 1]] <- if(normalized) {s/max(s) } else {s}
      labels <- c(labels, "Classical")
      series_id <- c(series_id, paste0("Classical_T",T_val,"_",i))
      Ts <- c(Ts, T_val)
      labels_m = c(labels_m, "gumbel")
      i=i+1
    }

    i=1
    while(i <= n_per_model) {
      s <- rweibull(T_val, shape = 2, scale=1)
      if(length(rec_gaps(s)) <2 ) next;
      all_series[[length(all_series) + 1]] <- if(normalized) {s/max(s) } else {s}
      labels <- c(labels, "Classical")
      series_id <- c(series_id, paste0("Classical_T",T_val,"_",i))
      Ts <- c(Ts, T_val)
      labels_m = c(labels_m, "weibull")
      i=i+1
    }

    ## end
  }

  return(list(series = all_series, labels = labels, labels_m = labels_m, series_id = series_id, T_vals = Ts))
}


# Helper: `%||%` operator - return left if not null else right
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ----- 2A. Feature extraction per series :Custom ---------------------------------------

# We'll compute a rich feature vector for each series.
# The feature extraction function returns a named list.

extract_custom_features <- function(series) {

  ## --- Preliminaries --------------------------------------------------------
  s <- as.numeric(series)
  n <- length(s)

  if (n < 10) stop("Series too short for stable feature extraction")

  ## Record helpers exist:
  ## is_rec(), rec_count(), rec_times(), rec_values(), rec_gaps()

  ## --- 1. Basic distributional statistics ----------------------------------
  ave  <- mean(s)
  med  <- median(s)
  std  <- sd(s)
  iqrv <- IQR(s)
  minv <- min(s)
  maxv <- max(s)
  rng  <- maxv - minv
  cv   <- ifelse(ave != 0, std / ave, NA)

  skew <- if (n > 2) moments::skewness(s) else NA
  kurt <- if (n > 3) moments::kurtosis(s) else NA

  ## --- 2. Trend features ----------------------------------------------------
  t <- seq_len(n)
  lmfit <- lm(s ~ t)

  slope        <- coef(lmfit)[2]  #fast_slope(s) #
  slope_pval   <- summary(lmfit)$coefficients[2, 4]
  slope_R2     <- summary(lmfit)$r.squared
  kendall_tau  <- cor(t, s, method = "kendall")

  ## --- 3. Increment / difference features ----------------------------------
  diffs <- diff(s)
  mean_diff <- mean(diffs)
  sd_diff   <- sd(diffs)
  mad_diff  <- mean(abs(diffs))

  diff_skew <- if (length(diffs) > 2)
    moments::skewness(diffs) else NA

  signs <- sign(diffs)
  sign_change_rate <- if (length(signs) > 2)
    mean(signs[-1] * signs[-length(signs)] < 0) else NA

  ## --- 4. Record counts and asymptotics ------------------------------------
  rec_nb <- rec_count(s)
  rec_rate <- rec_nb / n

  R_cum <- cumsum(is_rec(s))
  model_log <- lm(R_cum ~ log(t))

  beta_log      <- coef(model_log)[2]
  beta_log_pval <- summary(model_log)$coefficients[2, 4]
  beta_log_R2   <- summary(model_log)$r.squared

  ## --- 5. Forward / backward record rates ----------------------------------
  rec_times_f <- rec_times(s)
  rec_times_b <- rec_times(rev(s))

  rec_high_rate      <- length(rec_times_f) / n
  rec_high_back_rate <- length(rec_times_b) / n

  rec_rate_ratio <- ifelse(rec_high_back_rate > 0,
                           rec_high_rate / rec_high_back_rate, NA)

  ## --- 6. Record timing and span -------------------------------------------
  if (length(rec_times_f) >= 2) {
    span_norm <- (tail(rec_times_f, 1) - head(rec_times_f, 1)) / n
  } else {
    span_norm <- NA
  }

  frac_rec_first_half  <- mean(rec_times_f <= n / 2)
  frac_rec_last_quart  <- mean(rec_times_f > 3 * n / 4)

  ## --- 7. Inter-record gap statistics --------------------------------------
  rec_gap_f <- if (length(rec_times_f) >= 2) rec_gaps(s) else NA
  rec_gap_b <- if (length(rec_times_b) >= 2) rec_gaps(rev(s)) else NA

  rec_gap_med_f <- median(rec_gap_f, na.rm = TRUE)
  rec_gap_sd_f  <- sd(rec_gap_f, na.rm = TRUE)
  rec_gap_cv_f  <- ifelse(rec_gap_med_f > 0,
                          rec_gap_sd_f / mean(rec_gap_f, na.rm = TRUE), NA)

  rec_gap_med_b <- if (length(rec_times_b) >= 2) median(rec_gap_b, na.rm = TRUE) else n

  max_gap_over_n <- if (!all(is.na(rec_gap_f)))
    max(rec_gap_f, na.rm = TRUE) / n else NA

  ## --- 8. Record entropy ----------------------------------------------------
  p_rec <- mean(is_rec(s))
  # entropy <- if (p_rec > 0 && p_rec < 1)
  #   - (p_rec * log2(p_rec) + (1 - p_rec) * log2(1 - p_rec)) else 0

  entropy <- ifelse(p_rec > 0 & p_rec < 1, -p_rec*log(p_rec) - (1-p_rec)*log(1-p_rec), 0)

  ## --- 9. Spectral features -------------------------------------------------
  spec <- stats::spec.pgram(s, plot = FALSE, taper = 0)

  if (length(spec$spec) > 0) {
    dom_idx   <- which.max(spec$spec)
    dom_freq  <- 1/spec$freq[dom_idx]
    dom_power <- spec$spec[dom_idx]
  } else {
    dom_freq <- dom_power <- NA
  }

  ## --- 10. Crossing, extremes, dependence ---------------------------------
  cross_mean <- sum(diff(s > ave) != 0)

  extreme_2sd <- mean(abs(s - ave) > 2 * std)
  extreme_3sd <- mean(abs(s - ave) > 3 * std)

  #acf1 <- tryCatch(acf(s, plot = FALSE)$acf[2], error = function(e) NA)

  ## Stationarity proxies
  ndiff_needed <- tryCatch(forecast::ndiffs(s), error = function(e) NA)

  ljung_p <- tryCatch(
    Box.test(s, lag = 10, type = "Ljung-Box")$p.value,
    error = function(e) NA
  )

  ## --- 11. Local extrema ----------------------------------------------------
  diff_sign <- diff(sign(diffs))
  local_maxima <- mean(diff_sign == -2, na.rm = TRUE)
  local_minima <- mean(diff_sign ==  2, na.rm = TRUE)

  ## --- 12. Low (minimum) records via symmetry -----------------------------------
  s_low <- -s

  rec_times_low_f <- rec_times(s_low)
  rec_times_low_b <- rec_times(rev(s_low))

  rec_low_rate      <- length(rec_times_low_f) / n
  rec_low_back_rate <- length(rec_times_low_b) / n

  rec_low_rate_ratio <- ifelse(rec_low_back_rate > 0,
                               rec_low_rate / rec_low_back_rate, NA)

  rec_gap_low <- if (length(rec_times_low_f) >= 2) rec_gaps(s_low) else NA
  rec_gap_low_med <- if (length(rec_times_low_f) >= 2) median(rec_gap_low, na.rm = TRUE) else n

  ## --- 13. Rolling instability ---
  roll_feats <- rolling_instability(s)

  rolling_mean_instab = roll_feats[1]
  rolling_var_instab  = roll_feats[2]
  rolling_mean_sd_ratio = roll_feats[3]

  ## --- Feature vector -------------------------------------------------------
  features <- c(
    ave = ave, std = std, cv = cv,
    median = med, iqr = iqrv,
    skewness = skew, kurtosis = kurt,
    min = minv,
    #max = maxv,
    range = rng,

    slope = ifelse(slope_pval < 0.05, slope, 0),
    slope_R2 = ifelse(slope_pval < 0.05, slope_R2, 0),
    kendall_tau = kendall_tau,

    diff_mean = mean_diff,
    diff_sd = sd_diff,
    diff_mad = mad_diff,
    diff_skew = diff_skew,
    sign_change_rate = sign_change_rate,

    rec_rate = rec_rate,
    rec_high_rate = rec_high_rate,
    rec_high_back_rate = rec_high_back_rate,
    rec_rate_ratio = rec_rate_ratio,

    rec_gap_median = rec_gap_med_f,
    rec_gap_sd = rec_gap_sd_f,
    rec_gap_cv = rec_gap_cv_f,
    rec_gap_back_median = rec_gap_med_b,
    max_gap_over_n = max_gap_over_n,

    frac_rec_first_half = frac_rec_first_half,
    frac_rec_last_quarter = frac_rec_last_quart,
    span_norm = span_norm,

    entropy_shann = entropy,

    dom_period = dom_freq,
    dom_power = dom_power,

    #cross_mean = cross_mean, # we have crossing_points in ts
    extreme_2sd = extreme_2sd,
    #extreme_3sd = extreme_3sd,

    #acf1 = acf1,
    ndiff_needed = ndiff_needed,
    ljung_pvalue = ljung_p,

    local_minima = local_minima,
    local_maxima = local_maxima,

    rec_low_rate = rec_low_rate,
    rec_low_back_rate = rec_low_back_rate,
    rec_low_rate_ratio = rec_low_rate_ratio,
    rec_low_gap_median = rec_gap_low_med,

    rolling_mean_instab = rolling_mean_instab,
    rolling_var_instab  = rolling_var_instab,
    rolling_mean_sd_ratio = rolling_mean_sd_ratio

  )

  attr(features, "version") <- FEATURE_VERSION
  attr(features, "date") <- FEATURE_DATE
  return(features)
}


# ----- 2B. Feature extraction per series :tsfeatures ---------------------------------------

## acf_features: autocorrelation function of the series, the differenced series,
  ## and the twice-differenced series. It produces a vector comprising the first autocorrelation coefficient
  ## in each case, and the sum of squares of the first 10 autocorrelation coefficients
  ## "x_acf1, x_acf10, diff1_acf1, diff1_acf10, diff2_acf1, diff2_acf10"
## We compute the partial autocorrelation function of the series, the differenced series, and the second-order
  ## differenced series. Then pacf_features produces a vector comprising the sum of squares of the first 5
  ## partial autocorrelation coefficients in each case.
## arch_stat: Computes a statistic based on the Lagrange Multiplier (LM) test of Engle (1982)
  ## for autoregressive conditional heteroscedasticity (ARCH). The statistic returned is the R2
  ## value of an autoregressive model of order specified as lags applied to x2. "ARCH.LM"
## entropy: The spectral entropy is the Shannon entropy
## crossing_points defined as the number of times a time series crosses the median line.
## flat_spots are computed by dividing the sample space of a time series into ten equal-sized
  ##intervals, and computing the maximum run length within any single interval.
## The heterogeneity features measure the heterogeneity of the time series. First, we pre-whiten
  ## the time series to remove the mean, trend, and autoregressive (AR) information. Then we fit a GARCH(1,1)
  ## model to the pre-whitened time series, xt to measure for autoregressive conditional heteroskedasticity (ARCH) effects.
  ## The residuals from this model, zt are also measured for ARCH effects using a second GARCH(1,1)
  ## arch_acf is the sum of squares of the first 12 autocorrelations of {x2t}
  ## garch_acf is the sum of squares of the first 12 autocorrelations of {z2t}
  ## arch_r2 is the R2 value of an AR model applied to {x2t}
  ## garch_r2 is the R2 value of an AR model applied to {z2t}
## holt_parameters Estimate the smoothing parameter for the level-alpha and the smoothing parameter
  ## for the trend-beta of Holt’s linear trend method
## We use a measure of the long-term memory of a time series (hurst), computed as 0.5
  ## plus the maximum likelihood estimate of the fractional differencing order d
  ## given by Haslett & Raftery (1989). We add 0.5 to make it consistent with the Hurst coefficient.
## Stability and lumpiness are two time series features based on tiled (non-overlapping)
  ## windows. Means or variances are produced for all tiled windows. stability is the variance
  ## of the means, while lumpiness is the variance of the variances.
## The nonlinearity coefficient is computed using a modification of the statistic used in Teräsvirta’s nonlinearity test.
  ## Teräsvirta’s test uses a statistic X2=Tlog(SSE1/SSE0) where SSE1 and SSE0 are the sum of squared residuals
  ## from a nonlinear and linear autoregression respectively. This is non-ergodic, so instead, we define it as 10X2/T
  ## which will converge to a value indicating the extent of nonlinearity as T→∞
  ## This takes large values when the series is linear, and values around 0 when the series is nonlinear.
## station_features: std1st_der returns the standard deviation of the first derivative of the time series.
## spreadrandomlocal_meantaul_50: 100 time-series segments of length l are selected at random from the time series
    ## and the mean of the first zero-crossings of the autocorrelation function in each segment
extract_tsfeatures <- function(x) {
  tsf <- tsfeatures::tsfeatures(
    x,
    features = c("acf_features","pacf_features",
                 "arch_stat",
                 "crossing_points","entropy","flat_spots",
                 "heterogeneity", "hurst",
                 "holt_parameters", "lumpiness", "stability",
                 "nonlinearity", "station_features",
                 "stl_features"
                  )
  )
  return(as.list(tsf[1, ]))
}

# ----- 2C. Feature extraction per series : LogLik ---------------------------------------

## Extract Likelihood features
extract_LogLik_features <- function(series) {
  series <- as.numeric(series)
  n <- length(series)
  if (n < 10) stop("Series too short for stable features")
  if (length(rec_gaps(series)) <2) warning("No records found")

  ## theta hat
  time <- seq_len(length(series))
  lmfit <- lm(series ~ time)
  theta_hat = as.numeric(coef(lmfit)[2])

  ## gamma hat
  gamma_hat = estimate_model_param(series, method="mle_indicator", model = "YNM", min= 1, max=5, step = 0.001, approximate = FALSE)$param

  ## record data
  data_rec = data.frame(rec_values = rec_values(series),
                        rec_times = rec_times(series),
                        time = length(series),
                        theta = theta_hat,
                        gamma = gamma_hat)

  ## measures
  mean_all = mean(series)
  var_all = var(series)
  shape_all = abs(moments::skewness(series))
  mean_rec = mean(rec_values(series))
  var_rec =  var(rec_values(series))
  shape_rec = abs(moments::skewness(rec_values(series)))

  ## Classical -  all
  logLik_all_iid_gumbel = logLik_records(model = "iid", obs_type = "all",
                 dist = "gumbel", data = series,
                 params = c(loc = mean_all, scale=var_all))

  logLik_all_iid_norm = logLik_records(model = "iid", obs_type = "all",
                                         dist = "norm", data = series,
                                         params = c(mean = mean_all, sd=sqrt(var_all)) )

  logLik_all_iid_frechet = logLik_records(model = "iid", obs_type = "all",
                                       dist = "frechet", data = series,
                                       params = c(shape = shape_all, scale=var_all))

  logLik_all_iid_weibull = logLik_records(model = "iid", obs_type = "all",
                                          dist = "weibull", data = series,
                                          params = c(shape = shape_all, scale=var_all))

  ## Classical - Rn
  logLik_rec_iid_norm = logLik_records(model = "iid", obs_type = "records",
                                   dist = "norm", data = data_rec,
                                   params = c(mean = mean_rec, sd= sqrt(var_rec)) )

  logLik_rec_iid_gumbel = logLik_records(model = "iid", obs_type = "records",
                                       dist = "gumbel", data = data_rec,
                                       params = c(loc = mean_rec, scale=var_rec))

  logLik_rec_iid_frechet = logLik_records(model = "iid", obs_type = "records",
                                         dist = "frechet", data = data_rec,
                                         params = c(shape=shape_rec, scale=var_rec))

  logLik_rec_iid_weibull = logLik_records(model = "iid", obs_type = "records",
                                         dist = "weibull", data = data_rec,
                                         params =  c(shape=shape_rec, scale=var_rec))

  ## DTRW - all
  logLik_all_DTRW_norm = logLik_records(model = "DTRW", obs_type = "all",
                                          dist = "norm", data = series,
                                          params = c(mean = mean_all, sd= sqrt(var_all)))

  logLik_all_DTRW_cauchy = logLik_records(model = "DTRW", obs_type = "all",
                                        dist = "cauchy", data = series,
                                        params = c(loc = mean_all, scale=var_all))

  ## DTRW - rec
  logLik_rec_DTRW_norm = logLik_records(model = "DTRW", obs_type = "records",
                                        dist = "norm", data = data_rec,
                                        params = c(mean = mean_rec, sd= sqrt(var_all) ) )

  logLik_rec_DTRW_cauchy = logLik_records(model = "DTRW", obs_type = "records",
                                          dist = "cauchy", data = data_rec,
                                          params = c(loc = mean_rec, scale=var_all ))

  ## LDM - Xt
  logLik_all_LDM_norm = logLik_records(model = "LDM", obs_type = "all",
                                       dist = "norm", data = series,
                                       params = c(theta = theta_hat, mean = mean_all, sd=sqrt(var_all)))

  logLik_all_LDM_gumbel = logLik_records(model = "LDM", obs_type = "all",
                                         dist = "gumbel", data = series,
                                         params = c(theta = theta_hat, loc = mean_all, scale=var_all))

  logLik_all_LDM_frechet = logLik_records(model = "LDM", obs_type = "all",
                                          dist = "frechet", data = series,
                                          params = c(theta = theta_hat, shape=shape_all, scale=var_all))

  logLik_all_LDM_weibull = logLik_records(model = "LDM", obs_type = "all",
                                          dist = "weibull", data = series,
                                          params = c(theta = theta_hat, shape=shape_all, scale=var_all))

  ## LDM - rec
  logLik_rec_LDM_norm = logLik_records(model = "LDM", obs_type = "records",
                                       dist = "norm", data = data_rec,
                                       params = c(theta = theta_hat, mean = mean_rec, sd= sqrt(var_rec)))

  logLik_rec_LDM_gumbel = logLik_records(model = "LDM", obs_type = "records",
                                         dist = "gumbel", data = data_rec,
                                         params = c(theta = theta_hat, loc = mean_rec, scale=var_rec))

  logLik_rec_LDM_frechet = logLik_records(model = "LDM", obs_type = "records",
                                          dist = "frechet", data = data_rec,
                                          params = c(theta = theta_hat, shape=shape_rec, scale=var_rec))

  logLik_rec_LDM_weibull = logLik_records(model = "LDM", obs_type = "records",
                                          dist = "weibull", data = data_rec,
                                          params = c(theta = theta_hat, shape = shape_rec, scale=var_rec))

  ## YNM - Xt
  logLik_all_YNM_norm = logLik_records(model = "YNM", obs_type = "all",
                                       dist = "norm", data = series,
                                       params = c(gamma = gamma_hat, mean = mean_all, sd=sqrt(var_all)))

  logLik_all_YNM_gumbel = logLik_records(model = "YNM", obs_type = "all",
                                         dist = "gumbel", data = series,
                                         params = c(gamma = gamma_hat, loc = mean_all, scale=var_all))

  logLik_all_YNM_frechet = logLik_records(model = "YNM", obs_type = "all",
                                          dist = "frechet", data = series,
                                          params = c(gamma = gamma_hat, shape=shape_all, scale=var_all))

  logLik_all_YNM_weibull = logLik_records(model = "YNM", obs_type = "all",
                                          dist = "weibull", data = series,
                                          params = c(gamma = gamma_hat, shape =shape_all, scale=var_all))

  ## YNM - rec
  logLik_rec_YNM_norm = logLik_records(model = "YNM", obs_type = "records",
                                       dist = "norm", data = data_rec,
                                       params = c(gamma = gamma_hat, mean = mean_rec, sd= sqrt(var_rec)) )

  logLik_rec_YNM_gumbel = logLik_records(model = "YNM", obs_type = "records",
                                         dist = "gumbel", data = data_rec,
                                         params = c(gamma = gamma_hat, loc = mean_rec, scale=var_rec))

  logLik_rec_YNM_frechet = logLik_records(model = "YNM", obs_type = "records",
                                          dist = "frechet", data = data_rec,
                                          params = c(gamma = gamma_hat, shape=shape_rec, scale=var_rec))

  logLik_rec_YNM_weibull = logLik_records(model = "YNM", obs_type = "records",
                                          dist = "weibull", data = data_rec,
                                          params = c(gamma = gamma_hat, shape=shape_rec, scale=var_rec))
  Log_values = c(
    logLik_all_iid_gumbel =   logLik_all_iid_gumbel,
    logLik_all_iid_norm =   logLik_all_iid_norm,
    logLik_all_iid_frechet = logLik_all_iid_frechet,
    logLik_all_iid_weibull = logLik_all_iid_weibull,

    logLik_rec_iid_gumbel =   logLik_rec_iid_gumbel,
    logLik_rec_iid_norm =   logLik_rec_iid_norm,
    logLik_rec_iid_frechet = logLik_rec_iid_frechet,
    logLik_rec_iid_weibull = logLik_rec_iid_weibull,

    logLik_all_DTRW_norm = logLik_all_DTRW_norm,
    logLik_all_DTRW_cauchy = logLik_all_DTRW_cauchy,

    logLik_rec_DTRW_norm = logLik_rec_DTRW_norm,
    logLik_rec_DTRW_cauchy = logLik_rec_DTRW_cauchy,

    logLik_all_LDM_norm = logLik_all_LDM_norm,
    logLik_all_LDM_gumbel = logLik_all_LDM_gumbel,
    logLik_all_LDM_frechet = logLik_all_LDM_frechet,
    logLik_all_LDM_weibull = logLik_all_LDM_weibull,

    logLik_rec_LDM_norm = logLik_rec_LDM_norm,
    logLik_rec_LDM_gumbel = logLik_rec_LDM_gumbel,
    logLik_rec_LDM_frechet = logLik_rec_LDM_frechet,
    logLik_rec_LDM_weibull = logLik_rec_LDM_weibull,

    logLik_all_YNM_norm = logLik_all_YNM_norm,
    logLik_all_YNM_gumbel = logLik_all_YNM_gumbel,
    logLik_all_YNM_frechet = logLik_all_YNM_frechet,
    logLik_all_YNM_weibull = logLik_all_YNM_weibull,

    logLik_rec_YNM_norm = logLik_rec_YNM_norm,
    logLik_rec_YNM_gumbel = logLik_rec_YNM_gumbel,
    logLik_rec_YNM_frechet = logLik_rec_YNM_frechet,
    logLik_rec_YNM_weibull = logLik_rec_YNM_weibull
  )

  return(Log_values )
  }

# ----- 2D. Feature extraction per series : All ---------------------------------------

extract_all_features <- function(x) {
  # print("Extract custom features ...")
  a = extract_custom_features(x)

  # message("Extract time series features ...")
  b =  extract_tsfeatures(x)

  # message("Extract LogLik features ...")
  cc = extract_LogLik_features(x)

  Max_logLik = substr(names(which.max(cc)[1]),start = 12, stop = 25)

  idx_rec <- grepl("^logLik_rec_", names(cc))
  idx_all <- grepl("^logLik_all_", names(cc))

  Max_logLik_rec <- substr(names(cc[idx_rec][which.max(cc[idx_rec])]), start=12, stop = 25)
  Max_logLik_all <- substr(names(cc[idx_all][which.max(cc[idx_all])]), start=12, stop = 25)

  return(c(a, b, cc, Max_logLik = Max_logLik,
           Max_logLik_rec = Max_logLik_rec,
           Max_logLik_all = Max_logLik_all))
}

# ----- 3. Build labeled feature_matrix of features ---------------------------------------

## --- Helper function
generate_create_feature_dataset <- function(n_per_model = 10,
                                   T_vals = c(100, 200, 500), normalized = TRUE) {

  message("Generating series...")
  data <- generate_series_multiple(n_per_model, T_vals, normalized = normalized)

  n <- length(data$series)
  feature_list <- vector("list", n)

  message("Extracting features...")
  for (i in seq_len(n)) {

    x <- data$series[[i]]

    # extract all your features
    f <- extract_all_features(x)

    # add labels and ID
    f$label <- data$labels[i]
    f$series_id <- data$series_id[i]
    f$label_m = data$labels_m[i]

    # add T as a feature
    f$T_length <- data$T_vals[i]

    feature_list[[i]] <- f
  }

  features_df <- dplyr::bind_rows(feature_list)
  features_df$label = as.factor(features_df$label)
  message("Done.")
  return(list(feature = features_df, data = data ))
}

## --- Helper function
create_feature_dataset <- function(data) {

  n <- length(data$series)
  feature_list <- vector("list", n)

  message("Extracting features...")
  for (i in seq_len(n)) {

    x <- data$series[[i]]

    # extract all your features
    f <- extract_all_features(x)

    # add labels and ID
    f$label <- data$labels[i]
    f$series_id <- data$series_id[i]
    f$label_m = data$labels_m[i]

    # add T as a feature
    f$T_length <- data$T_vals[i]

    feature_list[[i]] <- f
  }

  features_df <- dplyr::bind_rows(feature_list)
  features_df$label = as.factor(features_df$label)
  message("Done.")
  return(feature = features_df)
}

# ----- 6. Train different classification methods --------------------------------

# We'll create a single function that takes a data frame (predictors + label) and trains:
# - multinomial logistic (glmnet, multinom)
# - random forest
# - xgboost
# - svmRadial (from caret/kernlab)
# - knn
# - naiveBayes
# - simple neural net (nnet)
#
# We'll use caret train with a consistent resampling scheme (repeated CV) and return models and results.

train_and_compare <- function(df, label_col = "label", id_col = "series_id",
                              seed = 42, do_parallel = FALSE) {
  set.seed(seed)

  ## Prepare
  df <- df %>% as.data.frame()
  df = na.omit(df)
  rownames(df) <- df[[id_col]]
  y <- as.factor(df[[label_col]])
  X_df <- df %>% dplyr::select(-one_of(label_col, id_col))

  ## partition into train/test (we'll use 80/20 stratified)
  train_index <- createDataPartition(y, p = 0.8, list = FALSE)
  train_data <- X_df[train_index, ]
  train_label <- y[train_index]
  test_data <- X_df[-train_index, ]
  test_label <- y[-train_index]

  # caret control
  trctrl <- trainControl(method = "repeatedcv",
                         number = 5,
                         repeats = 2,
                         classProbs = TRUE,
                         summaryFunction = multiClassSummary,
                         savePredictions = "final",
                         verboseIter = FALSE)

  models_list <- list()
  results <- list()

  # 1) multinom (multinomial logistic using nnet::multinom)
  cat("Training multinom...\n")
  m_multinom <- train(x = train_data, y = train_label,
                      method = "multinom",
                      trControl = trctrl,
                      trace = FALSE)
  models_list$multinom <- m_multinom

  # 2) glmnet (multinomial)
  cat("Training glmnet (multinomial)...\n")
  tunegrid_glmnet <- expand.grid(alpha = c(0, 0.5, 1), lambda = 10^seq(-3, 1, length = 10))
  m_glmnet <- train(x = train_data, y = train_label,
                    method = "glmnet",
                    tuneGrid = tunegrid_glmnet,
                    trControl = trctrl,
                    family = "multinomial")
  models_list$glmnet <- m_glmnet

  # 3) random forest
  #set.seed(seed)
  cat("Training randomForest...\n")
  m_rf <- train(x = train_data, y = train_label,
                method = "rf",
                trControl = trctrl,
                importance = TRUE,
                tuneLength = 5)
  models_list$rf <- m_rf

  # 4) xgboost (multiclass)
  # caret's xgbTree supports multiclass with label encoded as numeric starting at 0
  #set.seed(seed)
  cat("Training xgboost (caret xgbTree)...\n")
  # m_xgb <- train(x = train_data, y = train_label,
  #                method = "xgbTree",
  #                trControl = trctrl,
  #                tuneLength = 4)
  # models_list$xgb <- m_xgb

  # 5) SVM radial
  #set.seed(seed)
  cat("Training SVM radial...\n")
  m_svm <- train(x = train_data, y = train_label,
                 method = "svmRadial",
                 trControl = trctrl,
                 tuneLength = 4)
  models_list$svm <- m_svm

  # 6) kNN
  #set.seed(seed)
  cat("Training kNN...\n")
  m_knn <- train(x = train_data, y = train_label,
                 method = "knn",
                 trControl = trctrl,
                 tuneLength = 5)
  models_list$knn <- m_knn

  # 7) naive Bayes
  #set.seed(seed)
  cat("Training naiveBayes...\n")
  m_nb <- train(x = train_data, y = train_label,
                method = "naive_bayes",
                trControl = trctrl,
                tuneLength = 3)
  models_list$nb <- m_nb

  # 8) neural net (nnet)
  #set.seed(seed)
  cat("Training nnet...\n")
  m_nnet <- train(x = train_data, y = train_label,
                  method = "nnet",
                  trControl = trctrl,
                  tuneLength = 4,
                  trace = FALSE)
  models_list$nnet <- m_nnet

  # Evaluate on test set for each model
  evaluate_model <- function(model, test_data, test_label) {
    preds <- predict(model, newdata = test_data)
    probs <- tryCatch(predict(model, newdata = test_data, type = "prob"), error = function(e) NULL)
    cm <- confusionMatrix(preds, test_label)
    # multiclass AUC: average of one-vs-rest AUC if probs available
    auc_avg <- NA
    if (!is.null(probs)) {
      # compute multiclass AUC via micro-average: average of pairwise or one-vs-rest AUC
      labels <- levels(test_label)
      aucs <- c()
      for (lab in labels) {
        true_bin <- ifelse(test_label == lab, 1, 0)
        roc_obj <- tryCatch(pROC::roc(true_bin, probs[, lab], quiet = TRUE), error = function(e) NULL)
        if (!is.null(roc_obj)) {
          aucs <- c(aucs, pROC::auc(roc_obj))
        }
      }
      if (length(aucs) > 0) auc_avg <- mean(aucs)
    }
    list(confusion = cm, auc = auc_avg, predictions = preds, probs = probs)
  }

  for (mname in names(models_list)) {
    cat("Evaluating", mname, "on test data...\n")
    results[[mname]] <- evaluate_model(models_list[[mname]], test_data, test_label)
  }

  return(list(models = models_list, results = results, train_index = train_index,
       train_data = train_data, train_label = train_label,
       test_data = test_data, test_label = test_label))
}
