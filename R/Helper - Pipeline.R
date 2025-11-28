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
      i = i + 1
    }

    ## ---- LDM
    i = 1
    while(i <= n_per_model) {
      s <- generate_series(
        LDM_series,
        series_args = list(theta = runif(1,0.1,0.3),
                           dist = "frechet", shape=1, scale=2),
        T_val = T_val
      )
      if(length(rec_gaps(s)) <2 ) next;
      all_series[[length(all_series) + 1]] <- if(normalized) {s/max(s) } else {s}
      labels <- c(labels, "LDM")
      series_id <- c(series_id, paste0("LDM_T",T_val,"_",i))
      Ts <- c(Ts, T_val)
      i = i + 1
    }

    ## ---- YNM
    i=1
    while(i <= n_per_model) {
      s <- generate_series(
        YNM_series,
        series_args = list(gamma = runif(1,1.1,1.5),
                           dist = "frechet", shape=5, scale=5),
        T_val = T_val
      )
      if(length(rec_gaps(s)) <2 ) next;
      all_series[[length(all_series) + 1]] <- if(normalized) {s/max(s) } else {s}
      labels <- c(labels, "YNM")
      series_id <- c(series_id, paste0("YNM_T",T_val,"_",i))
      Ts <- c(Ts, T_val)
      i=i+1
    }

    ## ---- iid (Classical)
    i=1
    while(i <= n_per_model) {
      s <- rnorm(T_val)
      if(length(rec_gaps(s)) <2 ) next;
      all_series[[length(all_series) + 1]] <- if(normalized) {s/max(s) } else {s}
      labels <- c(labels, "Classical")
      series_id <- c(series_id, paste0("Classical_T",T_val,"_",i))
      Ts <- c(Ts, T_val)
      i=i+1
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
  if (length(rec_gaps(s)) <2) warning("No records found")
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
  rec_high_gap <- if (length(rec_high_idx) >= 2) rec_gaps(series) else 0
  #rec_low_gap <-  if (length(rec_low_idx) >= 2) rec_gaps(-series) else NA
  rec_high_gap_back <- if (length(rec_high_idx_back) >= 2) rec_gaps(rev(series)) else 0
  #rec_low_gap_back <-  if (length(rec_low_idx_back) >= 2) rec_gaps(-rev(series)) else NA

  rec_high_gap_median <- if (!all(is.na(rec_high_gap))) median(rec_high_gap, na.rm = TRUE) else 0
  #rec_low_gap_median <-  if (!all(is.na(rec_low_gap))) median(rec_low_gap, na.rm = TRUE) else NA
  rec_high_gap_back_median <- if (!all(is.na(rec_high_gap_back))) median(rec_high_gap_back, na.rm = TRUE) else 0
  #rec_low_gap_back_median <-  if (!all(is.na(rec_low_gap_back))) median(rec_low_gap_back, na.rm = TRUE) else NA

  rec_high_gap_sd <- if (!all(is.na(rec_high_gap))) sd(rec_high_gap, na.rm = TRUE) else 0
  #rec_low_gap_sd <-  if (!all(is.na(rec_low_gap))) sd(rec_low_gap, na.rm = TRUE) else NA
  #rec_high_gap_back_sd <- if (!all(is.na(rec_high_gap_back))) sd(rec_high_gap_back, na.rm = TRUE) else 0
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
    #rec_high_gap_back_sd = rec_high_gap_back_sd,
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


## extract all features
extract_all_features <- function(x) {
  c(
    extract_custom_features(x),
    extract_tsfeatures(x)
  )
}

# ----- 3. Build labeled feature_matrix of features ---------------------------------------

## --- Helper function
create_feature_dataset <- function(n_per_model = 50,
                                   T_vals = c(100, 200, 500)) {

  message("Generating series...")
  data <- generate_series_multiple(n_per_model, T_vals, normalized = TRUE)

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
  message("Done.")
  return(features_df)
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
