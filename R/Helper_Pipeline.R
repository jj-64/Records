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
rolling_instability <- function(x, width = 10) {
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

#' Rolling Instability for records
#'
#' A generic function for applying a series of functions to rolling margins
#' of a vector. The functions are mean, variance, and sd-to-mean ratio.
#'
#' @param x the data to be used (representing a series of observations).
#' @param width numeric value, NA by default. an integer specifying the window width (in numbers of observations).
#' @return A named vector the results of the rolling functions.
#' @export
#' @examples
#' \dontrun{
#' rolling_records(x=rnorm(100) , width = 2)}
rolling_records <- function(x, width = NA) {
  n <- length(x)
  if(!is.na(width)){
      if (n < 2 * width) return(c(NA, NA, NA))
  } else { width = n/2
  }

  rolls <- zoo::rollapply(
    x, width = width,
    FUN = function(x) c(rec_count(x), max(rec_gaps(x)), tail(rec_times(x),1) ),
    by.column = FALSE, align = "right"
  )
  ## sub-vectors with only one trivial record will have a record gap equal to width
  rolls[!is.finite(rolls[,2]),2] = width

  rec_rate <- mean(rolls[,1] / width)
  gaps <- median(rolls[,2])
  age_last_rec <- (width - tail(mean(rolls[,3]), 1) ) / width

  c(rec_rate_instab = rec_rate,
    max_rec_gaps_instab = gaps,
    age_last_rec_insta = age_last_rec )
}

fast_slope <- function(s) {
  n <- length(s)
  t_vec <- seq_len(n)
  return(cov(s, t_vec) / var(t_vec))
}

# Extreme Value Theory helper functions
fit_gev_block_maxima <- function(x, block_size = "sqrt") {
  # Fit GEV to block maxima
  if(block_size == "sqrt") block_size <- floor(sqrt(length(x)))

  n_blocks <- floor(length(x) / block_size)
  if(n_blocks < 5) return(list(shape = NA, scale = NA, loc = NA, se_shape = NA))

  block_maxima <- sapply(1:n_blocks, function(i) {
    max(x[((i-1)*block_size + 1):(i*block_size)])
  })

  tryCatch({
    fit <- fExtremes::gevFit(block_maxima, type = "pwm")
    params <- fit@fit$par.ests
    std_err <- fit@fit$par.ses
    list(shape = params["xi"], scale = params["beta"],
         loc = params["mu"], se_shape = std_err["xi"])
  }, error = function(e) {
    list(shape = NA, scale = NA, loc = NA, se_shape = NA)
  })
}

fit_gpd <- function(x, threshold_quantile = 0.85) {
  # Fit GPD to peaks over threshold
  threshold <- quantile(x, threshold_quantile, na.rm = TRUE)
  exceedances <- x[x > threshold] - threshold

  if(length(exceedances) < 3) return(list(shape = NA, scale = NA, threshold = threshold, exceedance_rate = NA))

  tryCatch({
    fit <- fExtremes::gpdFit(x, u = threshold, type = "pwm")
    params <- fit@fit$par.ests
    list(shape = params["xi"], scale = params["beta"],
         threshold = threshold, exceedance_rate = length(exceedances)/length(x))
  }, error = function(e) {
    list(shape = NA, scale = NA, threshold = threshold, exceedance_rate = NA)
  })
}

hill_estimator <- function(x, k = NULL) {
  # Hill estimator for tail index
  if(is.null(k)) k <- floor(0.1 * length(x))
  sorted_x <- sort(x, decreasing = TRUE)

  if(k < 2 || k >= length(x)) return(NA)

  log_ratios <- log(sorted_x[1:k] / sorted_x[k+1])
  alpha <- 1 / mean(log_ratios)
  alpha
}

extremal_index <- function(x, threshold_quantile = 0.95, r = 1) {
  # Estimate extremal index (clustering of extremes)
  threshold <- quantile(x, threshold_quantile, na.rm = TRUE)
  exceedances <- which(x > threshold)

  if(length(exceedances) < 2) return(NA)

  # Runs estimator
  clusters <- 1
  for(i in 2:length(exceedances)) {
    if(exceedances[i] - exceedances[i-1] > r) {
      clusters <- clusters + 1
    }
  }

  theta <- clusters / length(exceedances)
  theta
}

## accelaration index
#Acceleration (Record Arrival Rate): In a Linear Drift model, the probability of
# a record stays relatively high because the trend keeps pushing the values up.
#In Yang-Nevzorov, the probability $P(R_n) = \frac{\alpha}{\alpha + n - 1}$
#decays according to a very specific power-law-like curve.
compute_acceleration <- function(x) {
  n <- length(x)
  records <- as.numeric(x == cummax(x))
  cum_records <- cumsum(records)
  time_indices <- 1:n

  # Regress cumulative records against log(time)
  # The slope represents the 'intensity' or 'acceleration'
  fit <- lm(cum_records ~ log(time_indices))
  return(as.numeric(coef(fit)[2])) # Return the slope
}

## Detrending and Residual Records
#This is your "silver bullet." If you remove the linear trend from a Linear Drift series,
#the residuals become i.i.d. noise (the Classical Model). However, removing a
# linear trend from a Yang-Nevzorov series won't "fix" it because its structure
# isn't linear—it's probabilistic.
compute_residual_records <- function(x) {
  time <- 1:length(x)
  # Extract residuals from a linear model
  resids <- stats::resid(lm(x ~ time))

  # Count how many records occur in the 'cleaned' data
  num_resid_records <- sum(resids == cummax(resids))
  return(num_resid_records)
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

    ## ---- DTRW ------------
    i = 1
    while(i <= n_per_model) {
      s <- generate_series(
        dtrw_series,
        series_args = list(dist = "norm", mean = 0, sd = 1),
        T_val = T_val
      )
      if(length(rec_gaps(s)) <2 ) next;
      all_series[[length(all_series) + 1]] <- if(normalized) {s/max(s) } else {s}
      labels <- c(labels, "dtrw")
      series_id <- c(series_id, paste0("DTRW_T",T_val,"_",i))
      Ts <- c(Ts, T_val)
      labels_m = c(labels_m, "norm")
      i = i + 1
    }

    i = 1
    while(i <= n_per_model) {
      s <- generate_series(
        dtrw_series,
        series_args = list(dist = "norm", mean = 0, scale = 1),
        T_val = T_val
      )
      if(length(rec_gaps(s)) <2 ) next;
      all_series[[length(all_series) + 1]] <- if(normalized) {s/max(s) } else {s}
      labels <- c(labels, "dtrw")
      series_id <- c(series_id, paste0("DTRW_T",T_val,"_",i))
      Ts <- c(Ts, T_val)
      labels_m = c(labels_m, "cauchy")
      i = i + 1
    }

    i = 1
    while(i <= n_per_model) {
      s <- generate_series(
        dtrw_series,
        series_args = list(dist = "uniform", min = -1, scale = 1),
        T_val = T_val
      )
      if(length(rec_gaps(s)) <2 ) next;
      all_series[[length(all_series) + 1]] <- if(normalized) {s/max(s) } else {s}
      labels <- c(labels, "dtrw")
      series_id <- c(series_id, paste0("DTRW_T",T_val,"_",i))
      Ts <- c(Ts, T_val)
      labels_m = c(labels_m, "uniform")
      i = i + 1
    }

    ## ---- LDM ----------
        ## Frechet
    i = 1
    while(i <= n_per_model) {
      s <- generate_series(
        ldm_series,
        series_args = list(theta = runif(1,0.02,0.15),
                           dist = "frechet", shape=5, scale=1),
        T_val = T_val
      )
      if(length(rec_gaps(s)) <2 ) next;
      all_series[[length(all_series) + 1]] <- if(normalized) {s/max(s) } else {s}
      labels <- c(labels, "ldm")
      series_id <- c(series_id, paste0("LDM_T",T_val,"_",i))
      Ts <- c(Ts, T_val)
      labels_m = c(labels_m, "frechet")
      i = i + 1
    }
        ## Weibull
    i = 1
    while(i <= n_per_model) {
      s <- generate_series(
        ldm_series,
        series_args = list(theta = runif(1,0.02,0.15),
                           dist = "weibull", shape = 2, scale=1),
        T_val = T_val
      )
      if(length(rec_gaps(s)) <2 ) next;
      all_series[[length(all_series) + 1]] <- if(normalized) {s/max(s) } else {s}
      labels <- c(labels, "ldm")
      series_id <- c(series_id, paste0("LDM_T",T_val,"_",i))
      Ts <- c(Ts, T_val)
      labels_m = c(labels_m, "weibull")
      i = i + 1
    }
        ## Gumbel
    i = 1
    while(i <= n_per_model) {
      s <- generate_series(
        ldm_series,
        series_args = list(theta = runif(1,0.09,0.2),
                           dist = "gumbel", loc =0 , scale =1),
        T_val = T_val
      )
      if(length(rec_gaps(s)) <2 ) next;
      all_series[[length(all_series) + 1]] <- if(normalized) {s/max(s) } else {s}
      labels <- c(labels, "ldm")
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
        ynm_series,
        series_args = list(gamma = runif(1,1.2,1.4),
                           dist = "frechet", shape=5, scale=0.1),
        T_val = T_val
      )
      if(length(rec_gaps(s)) <2 ) next;
      all_series[[length(all_series) + 1]] <- if(normalized) {s/max(s) } else {s}
      labels <- c(labels, "ynm")
      series_id <- c(series_id, paste0("YNM_T",T_val,"_",i))
      Ts <- c(Ts, T_val)
      labels_m = c(labels_m, "frechet")
      i=i+1
    }
        ## Weibull
    i=1
    while(i <= n_per_model) {
      s <- generate_series(
        ynm_series,
        series_args = list(gamma = runif(1,1.4,1.7),
                           dist = "weibull", shape= 1/2, scale =0.1),
        T_val = T_val
      )
      if(length(rec_gaps(s)) <2 ) next;
      all_series[[length(all_series) + 1]] <- if(normalized) {s/max(s) } else {s}
      labels <- c(labels, "ynm")
      series_id <- c(series_id, paste0("YNM_T",T_val,"_",i))
      Ts <- c(Ts, T_val)
      labels_m = c(labels_m, "weibull")
      i=i+1
    }
        ## Pareto truncated
    i=1
    while(i <= n_per_model) {
      s <- generate_series(
        ynm_series,
        series_args = list(gamma = runif(1,1.3,2),
                           dist = "pareto", shape= 10, scale =1),
        T_val = T_val
      )
      if(length(rec_gaps(s)) <2 ) next;
      all_series[[length(all_series) + 1]] <- if(normalized) {s/max(s) } else {s}
      labels <- c(labels, "ynm")
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
      labels <- c(labels, "iid")
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
      labels <- c(labels, "iid")
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
      labels <- c(labels, "iid")
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

# ----- 2A. Feature extraction per series :Rec ---------------------------------

# We'll compute a rich feature vector for each series.
# The feature extraction function returns a named list.

extract_record_features <- function(series) {

  ## --- Preliminaries --------------------------------------------------------
  s <- as.numeric(series)
  s_n <- as.numeric (s/max(s))
  n <- length(s)

  if (n < 10) stop("Series too short for stable feature extraction")

  ## Record helpers exist:
  ## is_rec(), rec_count(), rec_times(), rec_values(), rec_gaps()

  ## --- 0. Pre-eliminary statistics ----------------------------------
  t_vec <- seq_len(n)
  rec_times_vec <- rec_times(s)
  rec_times_b <- rec_times(rev(s))
  rec_vals <- rec_values(s)
  rec_seq <- is_rec(s)
  rec_nb <- rec_count(s)
  diffs = diff(s)

  ## --- 1. Basic distributional statistics ----------------------------------
  ave  <- mean(s_n)
  med  <- median(s_n)
  std  <- sd(s_n)
  iqrv <- IQR(s_n)
  minv <- min(s)
  #maxv <- max(s)
  #rng  <- maxv - minv
  cv   <- ifelse(ave != 0, std / ave, NA)

  skew <- if (n > 2) moments::skewness(s) else NA
  kurt <- if (n > 3) moments::kurtosis(s) else NA

  ## --- 2. Trend features ----------------------------------------------------
  lmfit <- lm(s_n ~ t_vec)

  slope        <- coef(lmfit)[2]  #fast_slope(s) #
  slope_pval   <- summary(lmfit)$coefficients[2, 4]
  slope_R2     <- summary(lmfit)$r.squared
  #kendall_tau  <- cor(t_vec, s, method = "kendall")

  ## --- 3. Increment / difference features ----------------------------------
  diff_mean <- mean(diffs)
  diff_sd   <- sd(diffs)
  diff_mad  <- mean(abs(diffs))

  diff_skew <- if (length(diffs) > 2)
    moments::skewness(diffs) else NA

  signs <- sign(diffs)
  sign_change_rate <- if (length(signs) > 2)
    mean(signs[-1] * signs[-length(signs)] < 0) else NA

  ## --- 4. Record counts and asymptotics ------------------------------------

    # Estimate decay exponent (α in R_n ~ n^α)
    # Fit power law: log(R) ~ α*log(t)
  R_cum <- cumsum(rec_seq)
  model_log <- lm(log(R_cum[R_cum > 0]) ~ log(t_vec[R_cum > 0])) #lm(R_cum ~ log(t))

  beta_log      <- coef(model_log)[2]
  # beta_log_pval <- summary(model_log)$coefficients[2, 4]
  # beta_log_R2   <- summary(model_log)$r.squared

  ## --- 5. Forward / backward record rates ----------------------------------

  rec_rate <- rec_nb / n
  rec_rate_b <- length(rec_times_b) / n

  rec_rate_ratio <- ifelse(rec_rate_b > 0,
                          rec_rate / rec_rate_b, NA)

  #intensity = sum(rec_count(s)) / log(n)
  ## --- 6. Record timing and span -------------------------------------------
  age_last_rec <- (n-tail(rec_times_vec, 1)) / n

  # maximum jump
  max_jump = max(diff(cummax(s)), na.rm = TRUE)

  # Length of longest record streak
  longest_record_streak <- max(c(rec_gaps(s), (n-rec_times(s))[rec_nb]))

  frac_rec_first_half  <- mean(rec_times_vec <= n / 2)
  frac_rec_last_quart  <- mean(rec_times_vec > 3 * n / 4)

  # Record ratio
  frac_rec_last_first  <- if(rec_nb >2 ) {
            tail(rec_vals,1)/head(rec_vals,2)[2]
  } else { tail(rec_vals,1)/head(rec_vals,1)}

  #mean_frac_rec_last_first <- mean(rec_vals[-1]/ rec_vals[-rec_nb])

  ## --- 7. Inter-record gap statistics --------------------------------------
  mean_inter_time = mean(rec_times_vec, na.rm = TRUE)/n

  rec_gap <- if (length(rec_times_vec) >= 2) rec_gaps(s) else NA
  rec_gap_b <- if (length(rec_times_b) >= 2) rec_gaps(rev(s)) else NA

  rec_gap_mean = mean(rec_gap, na.rm = TRUE)
  rec_gap_med <- median(rec_gap, na.rm = TRUE)
  rec_gap_sd  <- sd(rec_gap, na.rm = TRUE)
  rec_gap_cv  <- ifelse(rec_gap_med> 0,
                          rec_gap_sd / rec_gap_mean , NA)

  rec_gap_med_b <- if (length(rec_times_b) >= 2) median(rec_gap_b, na.rm = TRUE) else n

  max_gap_over_n <- if (!all(is.na(rec_gap))) max(rec_gap, na.rm = TRUE) / n else NA

  # Random observation process proxy (using gaps between records)
  # if(length(rec_gap) > 2 && !all(is.na(rec_gap))) {
  #   # Kolmogorov-Smirnov test for exponentiality
  #   if(rec_gap_mean > 0 && length(rec_gap) > 5) {
  #     ks_test <- ks.test(rec_gap, "pexp", rate = 1/rec_gap_mean)
  #     rec_gap_exp_ks <- ks_test$statistic
  #   }
  # } else{
  #   rec_gap_exp_ks <- NA
  # }

  ## --- 8. Record entropy ----------------------------------------------------
  # entropy <- if (rec_rate > 0 && rec_rate < 1)
  #   - (rec_rate * log2(rec_rate) + (1 - rec_rate) * log2(1 - rec_rate)) else 0

  # entropy <- ifelse(rec_rate > 0 & rec_rate < 1, -rec_rate*log(rec_rate) - (1-rec_rate)*log(1-rec_rate), 0)

  ## --- 9. Spectral features -------------------------------------------------
  # spec <- stats::spec.pgram(s, plot = FALSE, taper = 0)
  #
  # if (length(spec$spec) > 0) {
  #   dom_idx   <- which.max(spec$spec)
  #   dom_freq  <- 1/spec$freq[dom_idx]
  #   dom_power <- spec$spec[dom_idx]
  # } else {
  #   dom_freq <- dom_power <- NA
  # }

  ## --- 10. Crossing, extremes, dependence ---------------------------------
  arith_ave = mean(s)
  cross_mean <- sum(diff(s > arith_ave) != 0)

  extreme_2sd <- mean(abs(s_n - mean(s_n)) > 2 * std)
  #extreme_3sd <- mean(abs(s - arith_ave) > 3 * std)

  #acf1 <- tryCatch(stats::acf(s, plot = FALSE)$acf[2], error = function(e) NA)

  ## Stationarity proxies
  ndiff_needed <- tryCatch(forecast::ndiffs(s), error = function(e) NA)

  # ljung_p <- tryCatch(
  #   Box.test(s, lag = 10, type = "Ljung-Box")$p.value,
  #   error = function(e) NA
  # )

  ## --- 11. Local extrema ----------------------------------------------------
  diff_sign <- diff(sign(diffs))
  #local_maxima <- mean(diff_sign == -2, na.rm = TRUE)
  #local_minima <- mean(diff_sign ==  2, na.rm = TRUE)

  ## --- 12. Low (minimum) records via symmetry -----------------------------------
  # s_low <- -s
  #
  # rec_times_low_f <- rec_times(s_low)
  # rec_times_low_b <- rec_times(rev(s_low))
  #
  # rec_low_rate      <- length(rec_times_low_f) / n
  # rec_low_back_rate <- length(rec_times_low_b) / n
  #
  # rec_low_rate_ratio <- ifelse(rec_low_back_rate > 0,
  #                              rec_low_rate / rec_low_back_rate, NA)
  #
  # rec_gap_low <- if (length(rec_times_low_f) >= 2) rec_gaps(s_low) else NA
  # rec_gap_low_med <- if (length(rec_times_low_f) >= 2) median(rec_gap_low, na.rm = TRUE) else n

  ## --- 13. Rolling instability ---
  roll_feats <- rolling_instability(s)

  rolling_mean_instab = roll_feats[1]
  rolling_var_instab  = roll_feats[2]
  rolling_mean_sd_ratio = roll_feats[3]

  rolling_rec <- rolling_records(s)

  rec_rate_instab = rolling_rec[1]
  max_rec_gaps_instab = rolling_rec[2]
  age_last_rec_instab = rolling_rec[3]

  ## --- 14. YNM vs LDM ----------------------
  residual_records = compute_residual_records(s)
  acceleration = compute_acceleration(s)

  ## --- 14. Sequence Pattern Features ----

  # Record clusters (records within k observations)
  # k_cluster <- floor(sqrt(n))
  # if(length(rec_times_vec) > 1) {
  #   rec_clusters <- 1
  #   for(i in 2:length(rec_times_vec)) {
  #     if(rec_times_vec[i] - rec_times_vec[i-1] > k_cluster) {
  #       rec_clusters <- rec_clusters + 1
  #     }
  #   }
  #   record_clusters <- rec_clusters
  # } else {
  #   record_clusters <- 1
  # }

  ## Transition probabilities
  transitions <- table(paste0(rec_seq[-n], "->", rec_seq[-1]))
  p_nonrec_to_rec <- transitions["0->1"] / sum(transitions["0->0"], transitions["0->1"], na.rm = TRUE)
  p_rec_to_rec <- transitions["1->1"] / sum(transitions["1->0"], transitions["1->1"], na.rm = TRUE)

  ## New
  slope = mean(diff((s-min(s))/(max(s)-min(s))))/var(diff((s-min(s))/(max(s)-min(s))))
  acf_diff1 <- acf(diff(s), plot=FALSE)$acf[2]
  vr2 =  vrtest::Auto.VR(s)$stat
  curv  <- sd(diff(diff(s)))

  ## --- Feature vector -------------------------------------------------------
  features <- c(
    ave = ave,
    std = std,
    cv = cv,
    median = med,
    iqr = iqrv,
    skewness = skew,
    kurtosis = kurt,
    #min = minv,
    #range = rng,

    slope = slope,# ifelse(slope_pval < 0.05 & slope_R2 >= 0.8, slope, 0),
    #slope_R2 = ifelse(slope_pval < 0.05, slope_R2, 0),
    #kendall_tau = kendall_tau,

    diff_mean = diff_mean,
    diff_sd = diff_sd,
    diff_mad = diff_mad,
    diff_skew = diff_skew,
    sign_change_rate = sign_change_rate,

    rec_rate = rec_rate,
    rec_rate_b = rec_rate_b,
    rec_rate_ratio = rec_rate_ratio,
    #intensity = intensity,
    convexity = mean(diff(diff(s))),

    age_last_rec = age_last_rec,
    max_jump = max_jump,
    longest_record_streak = longest_record_streak,
    frac_rec_first_half = frac_rec_first_half,
    frac_rec_last_quarter = frac_rec_last_quart,
    frac_rec_last_first = frac_rec_last_first,
    #mean_frac_rec_last_first = mean_frac_rec_last_first,

    mean_inter_time = mean_inter_time ,
    # rec_gap_mean = rec_gap_mean,
    rec_gap_median = rec_gap_med,
    rec_gap_sd = rec_gap_sd,
    rec_gap_cv = rec_gap_cv,
    rec_gap_back_median = rec_gap_med_b,
    max_gap_over_n = max_gap_over_n,
    #rec_gap_exp_ks = rec_gap_exp_ks,


#     entropy_shann = entropy,
#     dom_period = dom_freq,
#     dom_power = dom_power,

    #cross_mean = cross_mean, # we have crossing_points in ts
    extreme_2sd = extreme_2sd,
    #extreme_3sd = extreme_3sd,

    #acf1 = acf1,
    ndiff_needed = ndiff_needed,
    #ljung_pvalue = ljung_p,

    #local_minima = local_minima,
    #local_maxima = local_maxima,

    # rec_low_rate = rec_low_rate,
    # rec_low_back_rate = rec_low_back_rate,
    # rec_low_rate_ratio = rec_low_rate_ratio,
    # rec_low_gap_median = rec_gap_low_med,

    rolling_mean_instab = as.numeric(rolling_mean_instab),
    rolling_var_instab  = as.numeric(rolling_var_instab),
    rolling_mean_sd_ratio = as.numeric(rolling_mean_sd_ratio),

    rec_rate_instab = as.numeric(rec_rate_instab),
    max_rec_gaps_instab = as.numeric(max_rec_gaps_instab),
    age_last_rec_instab = as.numeric(age_last_rec_instab),

    # rec_clusters = rec_clusters,
    p_nonrec_to_rec = as.numeric(p_nonrec_to_rec),
    p_rec_to_rec = as.numeric(p_rec_to_rec),

    ## LDM vs YNM
    residual_records = residual_records,
    acceleration = acceleration,

    vr2 = vr2,
    curv = curv,
    acf_diff1 = acf_diff1
  )

  attr(features, "version") <- FEATURE_VERSION
  attr(features, "date") <- FEATURE_DATE
  return(features)
}

# ----- 2B. Feature extraction per series : EVT---------------------------------

extract_EVT_features <- function(series) {

  ## --- Preliminaries --------------------------------------------------------
  s <- as.numeric(series)
  n <- length(s)

  if (n < 10) stop("Series too short for stable feature extraction")

  ## Record helpers exist:
  ## is_rec(), rec_count(), rec_times(), rec_values(), rec_gaps()

  ## --- 0. Pre-eliminary statistics ----------------------------------
  t_vec <- seq_len(n)
  rec_vals <- rec_values(s)
  rec_nb <- rec_count(s)
  rec_rate <- rec_nb/n
  threshold <- quantile(s, 0.9, na.rm = TRUE)
  exceed_series <- as.numeric(s > threshold)

  ## --- 15. EXTREME VALUE THEORY FEATURES ------------------------------------

  # 3.1 GEV Parameter Estimates
  gev_fit <- fit_gev_block_maxima(s)
  gev_shape <- gev_fit$shape
  gev_scale <- gev_fit$scale
  gev_loc <- gev_fit$loc
  #gev_shape_se <- gev_fit$se_shape

  # 3.2 GPD Parameter Estimates
  gpd_fit <- fit_gpd(s)
  gpd_shape <- gpd_fit$shape
  gpd_scale <- gpd_fit$scale
  gpd_threshold <- gpd_fit$threshold
  gpd_exceedance_rate <- gpd_fit$exceedance_rate

  # Mean excess function slope
  threshold_seq <- quantile(s, probs = seq(0.7, 0.95, by = 0.05), na.rm = TRUE)
  mean_excess <- sapply(threshold_seq, function(u) {
    exceed <- s[s > u] - u
    if(length(exceed) > 0) mean(exceed) else NA
  })

  if(sum(!is.na(mean_excess)) > 2) {
    me_fit <- lm(mean_excess ~ threshold_seq)
    mean_excess_slope <- coef(me_fit)[2]  ##
  } else {mean_excess_slope = NA }

  ### --- 16. Tail Index Features ------------
  hill_tail_index <- hill_estimator(s)  ##
  pickands_tail_index <- tryCatch({ ##
    Dowd::PickandsEstimator(s, tail.size = floor(0.1*n))
  }, error = function(e) NA)

  # Tail heaviness ratio
  q99 <- quantile(s, 0.99, na.rm = TRUE)
  q95 <- quantile(s, 0.95, na.rm = TRUE)
  #tail_heaviness_ratio <- mean(s > q99, na.rm = TRUE) / mean(s > q95, na.rm = TRUE) ##

  ## ---- 17. Extreme Dependence Features --------------
  # extremal_index <- extremal_index(s) ##

  ## Autocorrelation of exceedances
  if(sum(exceed_series) > 3) { ##
    exceedance_acf1 <- stats::acf(exceed_series, plot = FALSE, na.action = stats::na.pass)$acf[2,1,1]
  } else { exceedance_acf1 <- NA}

  ## Return level estimates
  if(!is.na(gev_fit$shape) && !is.na(gev_fit$scale) && !is.na(gev_fit$loc)) {
    # 100-year return level (assuming 1 observation per time unit)
    T <- 100 ##
    return_level_100 <- gev_fit$loc + gev_fit$scale/gev_fit$shape * ((-log(1-1/T))^(-gev_fit$shape) - 1)
  } else { return_level_100 <- NA}

  # 3.5 Extreme Value Mixture Features
  # Probability of being in tail (using mixture threshold)
  #tail_probability_estimate <- mean(s > threshold, na.rm = TRUE) ##

  ## --- 18. CROSS-DOMAIN FEATURES (Record Theory + EVT) ----------------------

  # 4.1 Record-EVT Relationship Features
  if(!is.na(gpd_exceedance_rate) | !is.null(gpd_exceedance_rate)) {
    records_to_exceedances_ratio <- rec_rate/ gpd_exceedance_rate
  }

  # Difference between record value distribution and GEV fit
  # if( rec_nb >= 5 && !is.na(gev_fit$loc)) {
  #   # KS test between record values and GEV distribution
  #   pgev = fExtremes::pgev
  #   ks_gev <- ks.test(rec_vals, "pgev", xi = gev_fit$shape, mu = gev_fit$loc, beta = gev_fit$scale)
  #   record_gev_ks <- ks_gev$statistic ##
  # } else {record_gev_ks <- NA}
  #
  # # 4.3 Distributional Comparison Features
  # # KS statistic between record values and full dataset
  # if(length(rec_vals) > 5) {
  #   ks_full <- ks.test(rec_vals, s)
  #   record_full_ks <- ks_full$statistic ##
  # } else {record_full_ks <- NA}


  ## --- Feature vector -------------------------------------------------------
  features <- c(
    "gev_shape" = as.numeric(gev_shape),
    #"gev_scale" = as.numeric(gev_scale),
    #"gev_loc" = as.numeric(gev_loc),

    "gpd_shape" = as.numeric(gpd_shape),
    #"gpd_scale" = as.numeric(gpd_scale),
    #"gpd_threshold" = as.numeric(gpd_threshold),
    #"gpd_exceedance_rate" = as.numeric(gpd_exceedance_rate),

    "mean_excess_slope" = as.numeric(mean_excess_slope),

    "hill_tail_index"  = hill_tail_index ,
    "pickands_tail_index" =  pickands_tail_index,
    #"tail_heaviness_ratio" = tail_heaviness_ratio,

    #"exceedance_acf1" = as.numeric(exceedance_acf1),
    "return_level_100" =  as.numeric(return_level_100),
    #"tail_probability_estimate" = as.numeric(tail_probability_estimate),

    "records_to_exceedances_ratio" =  as.numeric(records_to_exceedances_ratio)
    #"record_gev_ks" = as.numeric(record_gev_ks),
    #"record_full_ks" =  as.numeric(record_full_ks)
  )

  attr(features, "version") <- FEATURE_VERSION
  attr(features, "date") <- FEATURE_DATE
  return(features)
}

# ----- 2C. Feature extraction per series :tsfeatures ---------------------------------------

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
  results = tsf[1, ]             ##ACF & PACF Features
  results = results %>% dplyr::select("x_acf1", # Lag-1 autocorrelation of original series
                               "x_acf10", # Sum of squared autocorrelations up to lag 10 (original series)
                               #"diff1_acf1", # Lag-1 autocorrelation of first-differenced series
                               #"diff1_acf10", # Sum of squared autocorrelations up to lag 10 (1st diff)
                               #"diff2_acf1", # Lag-1 autocorrelation of second-differenced series
                               #"diff2_acf10", #Sum of squared autocorrelations up to lag 10 (2nd diff)
                               #"x_pacf5",# Sum of absolute partial autocorrelations (PACF) up to lag 5
                               #"diff1x_pacf5" ,# Same as above but for 1st-differenced series
                               #"diff2x_pacf5",# Same for 2nd-differe

                               ## ARCH/GARCH Features
                               #"ARCH.LM", #Test statistic from ARCH LM test; detects conditional heteroskedasticity
                               #"arch_acf", "arch_r2", # Autocorrelations and R² from residuals of ARCH model
                               #"garch_acf", "garch_r2", #Same, but using GARCH model

                               ## Structural & Complexity Features
                               "crossing_points", # Number of times the series crosses its median
                               "entropy",##Spectral entropy; measures predictability/complexity
                               "flat_spots" , ##Number of flat segments in the series (constant value)
                               "hurst", ## Hurst exponent; >0.5 indicates long-term memory/persistence
                               "nonlinearity", ## Test statistic for nonlinearity (Teräsvirta test)

                               ## Trend/Seasonality (STL Features)
                               "trend", ## Strength of trend (from STL decomposition)
                               "spike", ##Measure of spikiness in the series
                               "linearity", ##Degree to which trend is linear
                               "curvature", ## Amount of curvature (non-linearity) in the trend
                               "e_acf1", "e_acf10", ## ACF at lag 1 and sum of ACFs up to lag 10 of remainder (after STL)

                               ##Heterogeneity & Local Variation
                               "stability", ##Variance change across time windows
                               "lumpiness", ##Variance across non-overlapping windows

                               ##Others
                               "alpha","beta", ##Parameters from Holt’s linear exponential smoothing (level & trend)
                               #"std1st_der" ## Std dev of first derivative (local change)
  )

  return(as.list(results))
}

# ----- 2D. Feature extraction per series : LogLik ---------------------------------------

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
  gamma_hat = estimate_model_param(series, method="mle_indicator", model = "ynm", min= 1.01, max=5, step = 0.001, approximate = FALSE)$param

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
  logLik_all_DTRW_norm = logLik_records(model = "dtrw", obs_type = "all",
                                          dist = "norm", data = series,
                                          params = c(mean = mean_all, sd= sqrt(var_all)))

  logLik_all_DTRW_cauchy = logLik_records(model = "dtrw", obs_type = "all",
                                        dist = "cauchy", data = series,
                                        params = c(loc = mean_all, scale=var_all))

  ## DTRW - rec
  logLik_rec_DTRW_norm = logLik_records(model = "dtrw", obs_type = "records",
                                        dist = "norm", data = data_rec,
                                        params = c(mean = mean_rec, sd= sqrt(var_all) ) )

  logLik_rec_DTRW_cauchy = logLik_records(model = "dtrw", obs_type = "records",
                                          dist = "cauchy", data = data_rec,
                                          params = c(loc = mean_rec, scale=var_all ))

  ## LDM - Xt
  logLik_all_LDM_norm = logLik_records(model = "ldm", obs_type = "all",
                                       dist = "norm", data = series,
                                       params = c(theta = theta_hat, mean = mean_all, sd=sqrt(var_all)))

  logLik_all_LDM_gumbel = logLik_records(model = "ldm", obs_type = "all",
                                         dist = "gumbel", data = series,
                                         params = c(theta = theta_hat, loc = mean_all, scale=var_all))

  logLik_all_LDM_frechet = logLik_records(model = "ldm", obs_type = "all",
                                          dist = "frechet", data = series,
                                          params = c(theta = theta_hat, shape=shape_all, scale=var_all))

  logLik_all_LDM_weibull = logLik_records(model = "ldm", obs_type = "all",
                                          dist = "weibull", data = series,
                                          params = c(theta = theta_hat, shape=shape_all, scale=var_all))

  ## LDM - rec
  logLik_rec_LDM_norm = logLik_records(model = "ldm", obs_type = "records",
                                       dist = "norm", data = data_rec,
                                       params = c(theta = theta_hat, mean = mean_rec, sd= sqrt(var_rec)))

  logLik_rec_LDM_gumbel = logLik_records(model = "ldm", obs_type = "records",
                                         dist = "gumbel", data = data_rec,
                                         params = c(theta = theta_hat, loc = mean_rec, scale=var_rec))

  logLik_rec_LDM_frechet = logLik_records(model = "ldm", obs_type = "records",
                                          dist = "frechet", data = data_rec,
                                          params = c(theta = theta_hat, shape=shape_rec, scale=var_rec))

  logLik_rec_LDM_weibull = logLik_records(model = "ldm", obs_type = "records",
                                          dist = "weibull", data = data_rec,
                                          params = c(theta = theta_hat, shape = shape_rec, scale=var_rec))

  ## YNM - Xt
  logLik_all_YNM_norm = logLik_records(model = "ynm", obs_type = "all",
                                       dist = "norm", data = series,
                                       params = c(gamma = gamma_hat, mean = mean_all, sd=sqrt(var_all)))

  logLik_all_YNM_gumbel = logLik_records(model = "ynm", obs_type = "all",
                                         dist = "gumbel", data = series,
                                         params = c(gamma = gamma_hat, loc = mean_all, scale=var_all))

  logLik_all_YNM_frechet = logLik_records(model = "ynm", obs_type = "all",
                                          dist = "frechet", data = series,
                                          params = c(gamma = gamma_hat, shape=shape_all, scale=var_all))

  logLik_all_YNM_weibull = logLik_records(model = "ynm", obs_type = "all",
                                          dist = "weibull", data = series,
                                          params = c(gamma = gamma_hat, shape =shape_all, scale=var_all))

  ## YNM - rec
  logLik_rec_YNM_norm = logLik_records(model = "ynm", obs_type = "records",
                                       dist = "norm", data = data_rec,
                                       params = c(gamma = gamma_hat, mean = mean_rec, sd= sqrt(var_rec)) )

  logLik_rec_YNM_gumbel = logLik_records(model = "ynm", obs_type = "records",
                                         dist = "gumbel", data = data_rec,
                                         params = c(gamma = gamma_hat, loc = mean_rec, scale=var_rec))

  logLik_rec_YNM_frechet = logLik_records(model = "ynm", obs_type = "records",
                                          dist = "frechet", data = data_rec,
                                          params = c(gamma = gamma_hat, shape=shape_rec, scale=var_rec))

  logLik_rec_YNM_weibull = logLik_records(model = "ynm", obs_type = "records",
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

# ----- 2E. Feature extraction per series : anomaly

extract_anomaly_features = function(series){
  s = as.numeric(series)

  ## mad_outlier_rate
        # High values → frequent abrupt deviations.
  # median_s <- median(s, na.rm = TRUE)
  # mad_s <- mad(s, constant = 1, na.rm = TRUE)
  #
  # robust_z <- (s - median_s) / mad_s
  #
  # mad_outlier_rate <- mean(abs(robust_z) > 3, na.rm = TRUE)

  ##Rolling Local Anomaly Score (Adaptive Z-score)
      # Detect local anomalies relative to rolling window behavior.

  # window_size <- floor(length(s) * 0.1)
  #
  # roll_mean <- zoo::rollapply(s, window_size, mean, fill = NA, align = "right")
  # roll_sd   <- zoo::rollapply(s, window_size, sd, fill = NA, align = "right")
  #
  # roll_z <- (s - roll_mean) / roll_sd
  #
  # rolling_outlier_rate <- mean(abs(roll_z) > 3, na.rm = TRUE)
  # rolling_max_z <- max(abs(roll_z), na.rm = TRUE)

  ## Change-Point Instability (Structural Break Count)
  # Anomalies often correspond to regime shifts.
  #
  # Use mean-shift detection via cumulative sum (CUSUM proxy).
  # s_centered <- s - mean(s, na.rm = TRUE)
  # cusum <- cumsum(s_centered) / sd(s, na.rm = TRUE)
  #
  # mean_shift_score <- max(abs(cusum), na.rm = TRUE) / length(s)

  # Spectral Residual Anomaly Score
  # Unexpected spikes relative to dominant frequency structure.
  #
  # Remove dominant spectral component and measure residual spikes.
  spec <- stats::spec.pgram(s, plot = FALSE)

  if(length(spec$spec) > 0) {
    dominant_power <- max(spec$spec)
    mean_power <- mean(spec$spec)
    spectral_spike_ratio <- dominant_power / mean_power
  } else {
    spectral_spike_ratio <- NA
  }
  # Isolation Forest Anomaly Score (Model-Based)
        # Model-free anomaly detection using isolation trees.
  # iso_model <- isotree::isolation.forest(matrix(s, ncol = 1), ntrees = 100)
  #
  # iso_scores <- predict(iso_model, matrix(s, ncol = 1), type = "score")
  #
  # isolation_anomaly_score <- mean(iso_scores, na.rm = TRUE)

  ##Extreme Jump Score (Derivative-Based)
  #Captures sudden jumps:

  # diff_s <- diff(s)
  # jump_threshold <- 3 * sd(diff_s, na.rm = TRUE)
  #
  # extreme_jump_rate <- mean(abs(diff_s) > jump_threshold, na.rm = TRUE)

  ## --- Feature vector -------------------------------------------------------
  features <- c(
    #mad_outlier_rate = mad_outlier_rate,
    #rolling_outlier_rate = rolling_outlier_rate,
    #rolling_max_z = rolling_max_z,
    # mean_shift_score= mean_shift_score,
    spectral_spike_ratio = spectral_spike_ratio,
    #isolation_anomaly_score = isolation_anomaly_score,
    #extreme_jump_rate = extreme_jump_rate
    )
  return(features)

}

# ----- 2F. Feature extraction per series : All ---------------------------------------

extract_all_features <- function(x) {
  # print("Extract custom features ...")
  rec = extract_record_features(x)

  evt = extract_EVT_features(x)

  # message("Extract time series features ...")
  tss =  extract_tsfeatures(x)

  # message("Extract Anomaly series features ...")
  #anom =  extract_anomaly_features(x)

  # message("Extract LogLik features ...")
  #loG = extract_LogLik_features(x)

  #Max_logLik = substr(names(which.max(loG)[1]),start = 12, stop = 25)

  #idx_rec <- grepl("^logLik_rec_", names(loG))
  #idx_all <- grepl("^logLik_all_", names(loG))

  #Max_logLik_rec <- substr(names(loG[idx_rec][which.max(loG[idx_rec])]), start=12, stop = 25)
  #Max_logLik_all <- substr(names(loG[idx_all][which.max(loG[idx_all])]), start=12, stop = 25)

  return(c(rec, evt = evt, tss
           # loG,
           # Max_logLik = Max_logLik,
           # Max_logLik_rec = Max_logLik_rec,
           # Max_logLik_all = Max_logLik_all
           ))
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
