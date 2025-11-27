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
generate_series_multiple_T <- function(
    n_per_model = 50,
    T_vals = c(100, 200, 500) ) {
  all_series <- list()
  labels <- c()
  series_id <- c()
  Ts <- c()  # store T for each series

  for (T_val in T_vals) {

    ## ---- DTRW
    for (i in 1:n_per_model) {
      s <- generate_series(
        DTRW_series,
        series_args = list(dist = "norm", mean = 0, sd = 1),
        T_val = T_val
      )
      all_series[[length(all_series) + 1]] <- s
      labels <- c(labels, "DTRW")
      series_id <- c(series_id, paste0("DTRW_T",T_val,"_",i))
      Ts <- c(Ts, T_val)
    }

    ## ---- LDM
    for (i in 1:n_per_model) {
      s <- generate_series(
        LDM_series,
        series_args = list(theta = runif(1,0.1,0.3),
                           dist = "frechet", shape=1, scale=2),
        T_val = T_val
      )
      all_series[[length(all_series) + 1]] <- s
      labels <- c(labels, "LDM")
      series_id <- c(series_id, paste0("LDM_T",T_val,"_",i))
      Ts <- c(Ts, T_val)
    }

    ## ---- YNM
    for (i in 1:n_per_model) {
      s <- generate_series(
        YNM_series,
        series_args = list(gamma = runif(1,1.1,1.5),
                           dist = "frechet", shape=5, scale=5),
        T_val = T_val
      )
      all_series[[length(all_series) + 1]] <- s
      labels <- c(labels, "YNM")
      series_id <- c(series_id, paste0("YNM_T",T_val,"_",i))
      Ts <- c(Ts, T_val)
    }

    ## ---- iid (Classical)
    for (i in 1:n_per_model) {
      s <- rnorm(T_val)
      all_series[[length(all_series) + 1]] <- s
      labels <- c(labels, "Classical")
      series_id <- c(series_id, paste0("Classical_T",T_val,"_",i))
      Ts <- c(Ts, T_val)
    }
  }

  return(list(series = all_series, labels = labels, series_id = series_id, T_vals = Ts))
}


# Helper: `%||%` operator - return left if not null else right
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ----- 2. Feature extraction per series ---------------------------------------

# We'll compute a rich feature vector for each series.
# The feature extraction function returns a named list.

extract_custom_features <- function(series) {
  s <- as.numeric(series)
  n <- length(s)
  if (n < 10) stop("Series too short for stable features")

  # Basic stats
  ave = mean(s)
  cv <- sd(s)/mean(s)
  std <- sd(s)
  med <- median(s)
  iqrv <- IQR(s)
  skew <- ifelse(length(na.omit(s))>2, moments::skewness(s), NA)
  kurt <- ifelse(length(na.omit(s))>3, moments::kurtosis(s), NA) # excess? moments::kurtosis returns fisher or pearson? okay

  # extremes
  minv <- min(s)
  maxv <- max(s)
  rng <- maxv - minv

  # record info
  rec_high_idx <- rec_times(s)   # Foreword record highs
  rec_high_values <- rec_values(s)
  rec_high_idx_back <- rec_times(rev(s))  # Backward record highs
  rec_high_values_back <- rec_values(rev(s))
  #rec_low_idx <- rec_times(-s)   # Foreword record lows
  #rec_low_values <- rec_values(-s)
  #rec_low_idx_back <- rec_times(rev(-s))  # Backward record lows
  #rec_low_values_back <- rec_values(rev(-s))

  rec_high_rate <- length(rec_high_idx)/n
  rec_high_back_rate <- length(rec_high_idx_back)/n
  #rec_low_rate <- length(rec_low_idx)/n
  #rec_low_back_rate <- length(rec_low_idx_back)/n


  # inter-record intervals (for highs and lows, forward and backward)
  rec_high_gap <- if (length(rec_high_idx) >= 2) rec_gaps(series) else NA
  #rec_low_gap <-  if (length(rec_low_idx) >= 2) rec_gaps(-series) else NA
  rec_high_gap_back <- if (length(rec_high_idx_back) >= 2) rec_gaps(rev(series)) else NA
  #rec_low_gap_back <-  if (length(rec_low_idx_back) >= 2) rec_gaps(-rev(series)) else NA

  rec_high_gap_median <- if (!all(is.na(rec_high_gap))) median(rec_high_gap, na.rm = TRUE) else NA
  #rec_low_gap_median <-  if (!all(is.na(rec_low_gap))) median(rec_low_gap, na.rm = TRUE) else NA
  rec_high_gap_back_median <- if (!all(is.na(rec_high_gap_back))) median(rec_high_gap_back, na.rm = TRUE) else NA
  #rec_low_gap_back_median <-  if (!all(is.na(rec_low_gap_back))) median(rec_low_gap_back, na.rm = TRUE) else NA

  rec_high_gap_sd <- if (!all(is.na(rec_high_gap))) sd(rec_high_gap, na.rm = TRUE) else NA
  #rec_low_gap_sd <-  if (!all(is.na(rec_low_gap))) sd(rec_low_gap, na.rm = TRUE) else NA
  rec_high_gap_back_sd <- if (!all(is.na(rec_high_gap_back))) sd(rec_high_gap_back, na.rm = TRUE) else NA
  #rec_low_gap_back_sd <-  if (!all(is.na(rec_low_gap_back))) sd(rec_low_gap_back, na.rm = TRUE) else NA

  # Entropy of record indicator
  p <- mean(is_rec(s))
  entropy <- if (p > 0 && p < 1) {
    - (p * log2(p) + (1 - p) * log2(1 - p))
  } else {
    0
  }

  # trend: linear slope and p-value
  time <- seq_len(n)
  lmfit <- lm(s ~ time)
  slope <- as.numeric(coef(lmfit)[2])
  slope_pval <- summary(lmfit)$coefficients[2, 4]

  # Mann-Kendall proxy: Kendall correlation between time and series
  kendall_tau <- cor(time, s, method = "kendall")

  # autocorrelation
  acf_vals <- acf(s, lag.max = 5, plot = FALSE)$acf
  acf1 <- ifelse(length(acf_vals) >= 2, acf_vals[2], NA)
  acf2 <- ifelse(length(acf_vals) >= 3, acf_vals[3], NA)


  # spectral: dominant frequency via periodogram
  spec <- stats::spec.pgram(s, plot = FALSE, taper = 0)
  if (length(spec$spec) > 0) {
    dom_idx <- which.max(spec$spec)
    dom_freq <- spec$freq[dom_idx]
    dom_power <- spec$spec[dom_idx]
  } else {
    dom_freq <- NA; dom_power <- NA
  }

  # crossing and peak counts
  cross_mean <- sum(diff(s > ave) != 0)  # number of times series crosses its mean
  prop_pos <- mean(s > mean(rec_values(s)))

  extreme_2sd <- sum(abs(s - ave) > 2 * std) / n
  extreme_3sd <- sum(abs(s - ave) > 3 * std) / n

  # variability of first differences
  diffs <- diff(s)
  diff_mean <- mean(diffs)
  diff_sd <- sd(diffs)
  diff_skew <- ifelse(length(diffs) > 2, moments::skewness(diffs), NA)

  # adf-like: use ndiffs from forecast to estimate number of differences needed (proxy for non-stationarity)
  ndiff_needed <- tryCatch(forecast::ndiffs(s), error = function(e) NA)
  ljung <- tryCatch(Box.test(s, lag = 10, type = "Ljung-Box"), error = function(e) list(statistic = NA, p.value = NA))

  # proportion of positive first differences
  prop_inc <- mean(diffs > 0)

  # Turning points
  diff_sign <- diff(sign(diff(s)))
  turning_points <- sum(diff_sign != 0, na.rm = TRUE)

  # Local extrema
  local_maxima <- length(which(diff_sign == -2))/n
  local_minima <- length(which(diff_sign == 2))/n

  # feature vector
  features <- c(
    ave = ave,
    std = std,
    median = med,
    iqr = iqrv,
    skewness = skew,
    kurtosis = kurt,
    min = minv,
    max = maxv,
    range = rng,
    rec_high_rate = rec_high_rate,
    #rec_low_rate = rec_low_rate,
    rec_high_back_rate = rec_high_back_rate,
    #rec_low_back_rate = rec_low_back_rate,

    rec_high_gap_median = rec_high_gap_median,
    #rec_low_gap_median = rec_low_gap_median,
    rec_high_gap_back_median = rec_high_gap_back_median,
    #rec_low_gap_back_median = rec_high_gap_back_median,

    rec_high_gap_sd = rec_high_gap_sd,
    #rec_low_gap_sd = rec_low_gap_sd,
    rec_high_gap_back_sd = rec_high_gap_back_sd,
    #rec_low_gap_back_sd = rec_high_gap_back_sd,
    entropy_c = entropy,

    slope = slope,
    slope_p = slope_pval,
    kendall_tau = kendall_tau,

    acf1 = acf1,
    acf2 = acf2,
    #pacf1 = pacf1,

    dom_freq = dom_freq,
    dom_power = dom_power,

    cross_mean = cross_mean,
    prop_pos = prop_pos,
    extreme_2sd = extreme_2sd,
    extreme_3sd = extreme_3sd,

    diff_mean = diff_mean,
    diff_sd = diff_sd,
    diff_skew = diff_skew,

    ndiff_needed = ndiff_needed,
    ljung = ljung$p.value,

    prop_inc = prop_inc,

    turning_points = turning_points,
    local_minima = local_minima,
    local_maxima = local_maxima
  )

  # make numeric and named
  #as.numeric(features) %>% setNames(names(features))
  return(features)
}

## extract tsfeatures
extract_tsfeatures <- function(x) {
  tsf <- tsfeatures::tsfeatures(
    x,
    features = c("mean", "var", "acf_features", "entropy",
                 "lumpiness", "stability", "heterogeneity",
                 "flat_spots", "hurst")
  )
  return(as.list(tsf[1, ]))
}

## extract TSEntropies
extract_ts_entropies <- function(x) {
  list(
    perm_entropy = TSEntropies::ApEn(x),
    samp_entropy = TSEntropies::SampEn(x)
  )
}

## extract all features
extract_all_features <- function(x) {
  c(
    extract_custom_features(x),
    extract_tsfeatures(x),
    extract_ts_entropies(x)
  )
}

# ----- 3. Build labeled feature_matrix of features ---------------------------------------

## --- Helper function
create_feature_dataset <- function(n_per_model = 50,
                                   T_vals = c(100, 200, 500)) {

  message("Generating multiple series...")

  data <- generate_series_multiple_T(n_per_model, T_vals)

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

    # add T as a feature
    f$T_length <- data$T_vals[i]

    feature_list[[i]] <- f
  }

  features_df <- dplyr::bind_rows(feature_list)
  return(features_df)
}
