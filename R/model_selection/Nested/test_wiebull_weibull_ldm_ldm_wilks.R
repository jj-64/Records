###############################################################################
# PART A: TYPE-I ERROR OF THE WILKS LIKELIHOOD RATIO TEST
#
# H0 : Exponential LDM
# H1 : Weibull LDM
#
# Test statistic:
#
#     LR = 2 ( l_H1 - l_H0 )
#
# Under Wilks theorem and assuming H0 is nested in H1:
#
#     LR ~ Chi-square(df = p1 - p0)
#
###############################################################################

set.seed(12345)

###############################################################################
# DATA GENERATING MECHANISMS --------
###############################################################################

## Generate a series under H0
series_H0 <- function(T_val, par_H0){

  # Exponential Linear Drift Model
  ldm_series(
    T      = T_val,
    dist   = "weibull",
    theta  = par_H0[["trend"]],
    shape  = 1,
    scale  = par_H0[["scale"]]
  )
}

## Generate a series under H1
series_H1 <- function(T_val, trend, par_H1){

  # Weibull Linear Drift Model
  ldm_series(
    T      = T_val,
    dist   = "weibull",
    theta  = trend,
    shape  = par_H1[["shape"]],
    scale  = par_H1[["scale"]]
  )
}

###############################################################################
# LIKELIHOOD FUNCTIONS --------
###############################################################################

Likelihood_under_H0 <- function(data_rec, params){

  logLik_fun_rec <-
    loglik_registry[["LDM"]][["records"]][["weibull"]]

  logLik_fun_rec(
    data   = data_rec,
    params = params
  )
}

Likelihood_under_H1 <- function(data_rec, params){

  logLik_fun_rec <-
    loglik_registry[["LDM"]][["records"]][["weibull"]]

  logLik_fun_rec(
    data   = data_rec,
    params = params
  )
}

###############################################################################
# Part A: SIMULATION SETTINGS ------
###############################################################################

T_val      <- 100
alpha      <- 0.05
simulation <- 100

trend_values <- seq(0.2, 0.50, by = 0.01)

## True parameters under H0
true_params <- list(
  trend = NA,
  shape  = 1,
  scale = 2
)


final <- matrix(0, nrow = length(trend_values), ncol = 4)
colnames(final) = c("trend", "Type_I_Error", "Power", "Nb_record")

kk <- 1

###############################################################################
# LOOP OVER TRUE DRIFT VALUES
###############################################################################

for(trend_val in trend_values){

  cat("\n====================================================\n")
  cat("True drift =", trend_val, "\n")
  cat("====================================================\n")

  true_params["trend"] = trend_val

  results <- data.frame(
    theta_H0 = rep(NA, simulation),
    scale_H0  = NA,
    logLik_H0 = NA,

    theta_H1 = NA,
    shape_H1 = NA,
    scale_H1 = NA,
    logLik_H1 = NA,

    LR       = NA,
    p_value  = NA,
    nRecords = NA
  )

  trial <- 1

  ###########################################################################
  # MONTE-CARLO LOOP
  ###########################################################################

  while(trial <= simulation){

    print( paste( "Trend value = ", trend_val, ", Trial = ", trial))
    #__________________________________________________________________________
    # STEP 1: GENERATE DATA UNDER H0
    #__________________________________________________________________________

    xt <- series_H0(
      T      = T_val,
      par_H0 = true_params
    )

    R <- rec_values(xt)
    L <- rec_times(xt)

    ## Reject degenerate samples with only one record
    while(length(R) <= 1){

      xt <- series_H0(
        T      = T_val,
        par_H0 = true_params
      )

      R <- rec_values(xt)
      L <- rec_times(xt)
    }

    m = length(R)

    data_rec <- list(
      rec_values = R,
      rec_times  = L,
      time       = T_val
    )

    #__________________________________________________________________________
    # STEP 2: FIT H0 (EXPONENTIAL LDM)
    #__________________________________________________________________________

    lb_0 <- c(
      theta = 0.01,
      shape = 1,
      scale = 0.01
    )

    ub_0 <- c(
      theta = min(R / L),
      shape = 1,
      scale = 10
    )

    negLogLik_H0 <- function(params) {
      -Likelihood_under_H0(data_rec, params)
    }

    ## Generate multiple feasible starting values

    n_starts <- 10
    starts_H0 <- vector("list", n_starts)
    for (s in seq_len(n_starts)) {

      starts_H0[[s]] <- c(
        theta = runif(
          1,
          lb_0[["theta"]],
          ub_0[["theta"]]
        ),

        shape = lb_0[["shape"]],

        scale = runif(
          1,
          lb_0[["scale"]],
          ub_0[["scale"]]
        )
      )
    }

    ## Add a deterministic starting value
    ## Useful because shape is fixed at its boundary.

    starts_H0[[1]] <- c(
      theta = lb_0["theta"],
      shape = lb_0[["shape"]],
      scale = lb_0["scale"]
    )

    ## Run optimization from each starting value
    fits_H0 <- vector("list",length(starts_H0) )

    for (s in seq_along(starts_H0)) {
      fits_H0[[s]] <- tryCatch(
        nlminb(
          start = as.list(starts_H0[[s]]),
          objective = negLogLik_H0,
          lower = as.list(lb_0),
          upper = as.list(ub_0),
          control = list(
            eval.max = 1000,
            iter.max = 500,
            rel.tol = 1e-8
          )
        ),
        error = function(e) NULL
      )
    }

    ## Keep finite optimization results
    valid_fits <- vapply(
      fits_H0,
      function(fit) {
        !is.null(fit) &&
          is.finite(fit$objective) &&
          all(is.finite(fit$par))

      },
      logical(1)
    )

    fits_H0_valid <- fits_H0[valid_fits]

    ## Check whether at least one solution was obtained
    if (length(fits_H0_valid) == 0) {
      next
    }

    ## Select the solution with the smallest negative log-likelihood
    objectives_H0 <- vapply(
      fits_H0_valid,
      function(fit) fit$objective,
      numeric(1)
    )

    best_H0 <- which.min(objectives_H0)
    fit_H0 <- fits_H0_valid[[best_H0]]

    ## Final validity check
    if (
      !is.finite(fit_H0$objective) ||
      any(!is.finite(fit_H0$par))
    ) {
      next
    }

    ## Store estimated parameters and log-likelihood
    param_H0 <- c(
      fit_H0$par,
      logLik = -fit_H0$objective
    )

    #__________________________________________________________________________
    # STEP 3: FIT H1 (WEIBULL LDM)
    #__________________________________________________________________________

    lb_1 <- c(
      theta = 0.01,
      shape = 0.1,
      scale = 0.1
    )

    ub_1 <- c(
      theta = min(R / L),
      shape = 10,
      scale = 10
    )

    negLogLik_H1 <- function(params) {
      -Likelihood_under_H1(data_rec, params)
    }

    ## Generate multiple feasible starting values

    n_starts <- 10
    starts_H1 <- vector("list", n_starts)

    for (s in seq_len(n_starts)) {
      starts_H1[[s]] <- c(
        theta = runif(
          1,
          lb_1["theta"],
          ub_1["theta"]
        ),

        shape = runif(
          1,
          lb_1["shape"],
          ub_1["shape"]
        ),

        scale = runif(
          1,
          lb_1["scale"],
          ub_1["scale"]
        )
      )
    }

    ## Add a deterministic starting value
    ## Useful because shape is fixed at its boundary.

    starts_H1[[1]] <- c(
      theta = lb_1["theta"],
      shape = lb_1["shape"],
      scale = lb_1["scale"]
    )

    ## Run optimization from each starting value
    fits_H1 <- vector(
      "list",
      length(starts_H1)
    )

    for (s in seq_along(starts_H1)) {
      fits_H1[[s]] <- tryCatch(
        nlminb(
          start = starts_H1[[s]],
          objective = negLogLik_H1,
          lower = lb_1,
          upper = ub_1,
          control = list(
            eval.max = 1000,
            iter.max = 500,
            rel.tol = 1e-8
          )
        ),
        error = function(e) NULL
      )
    }

    ## Keep finite optimization results
    valid_fits <- vapply(
      fits_H1,
      function(fit) {
        !is.null(fit) &&
          is.finite(fit$objective) &&
          all(is.finite(fit$par))

      },
      logical(1)
    )

    fits_H1_valid <- fits_H1[valid_fits]

    ## Check whether at least one solution was obtained
    if (length(fits_H1_valid) == 0) {
      next
    }

    ## Select the solution with the smallest negative log-likelihood
    objectives_H1 <- vapply(
      fits_H1_valid,
      function(fit) fit$objective,
      numeric(1)
    )

    best_H1 <- which.min(objectives_H1)
    fit_H1 <- fits_H1_valid[[best_H1]]

    ## Final validity check
    if (
      !is.finite(fit_H1$objective) ||
      any(!is.finite(fit_H1$par))
    ) {
      next
    }

    ## Store estimated parameters and log-likelihood
    param_H1 <- c(
      fit_H1$par,
      logLik = -fit_H1$objective
    )

    #__________________________________________________________________________
    # STEP 4: WILKS LIKELIHOOD RATIO TEST
    #__________________________________________________________________________

    logLik_H0 <- param_H0["logLik"]
    logLik_H1 <- param_H1["logLik"]

    LR <- 2 * (logLik_H1 - logLik_H0)

    ## Numerical safeguard
    LR <- max(0, LR)
    #if(LR <0) next;
    df_wilks <- 1   # Weibull adds one free parameter (shape)

    p_value <- 1 - pchisq(
      q  = LR,
      df = df_wilks
    )

    #__________________________________________________________________________
    # STEP 5: STORE RESULTS
    #__________________________________________________________________________

    results[trial, ] <- c(
      param_H0["theta"],
      param_H0["scale"],
      logLik_H0,

      param_H1["theta"],
      param_H1["shape"],
      param_H1["scale"],
      logLik_H1,

      LR,
      p_value,
      m
    )

    trial <- trial + 1
  }

  #__________________________________________________________________________
  # ESTIMATED TYPE-I ERROR
  #__________________________________________________________________________

  rejection_rate <- mean(
    results$p_value < alpha,
    na.rm = TRUE
  )

  results$reject_H0 = ifelse( results$p_value < alpha, 1, 0)

  cat("Estimated Type-I Error =", rejection_rate, "\n")

  ###############################################################################
  # PART B: EMPIRICAL POWER OF THE WILKS TEST
  #
  # Data generated under H1 (Weibull-LDM)
  #
  # H0 : Exponential LDM
  # H1 : Weibull LDM
  #
  ###############################################################################

  results_power <- data.frame(
    theta_H0  = rep(NA, simulation),
    rate_H0   = NA,
    logLik_H0 = NA,

    theta_H1  = NA,
    shape_H1  = NA,
    scale_H1  = NA,
    logLik_H1 = NA,

    LR        = NA,
    reject    = NA,
    nRecords  = NA
  )

  trial <- 1

  ###############################################################################
  # MONTE-CARLO LOOP
  ###############################################################################

  while(trial <= simulation){

    #_________________________________________________________________________
    # STEP 1: GENERATE DATA UNDER H1
    #_________________________________________________________________________

    xt <- series_H1(
      T      = T_val,
      trend  = true_params[["trend"]],
      par_H1 = c(
        shape = 1.20,
        scale = true_params[["rate"]]
      )
    )

    R <- rec_values(xt)
    L <- rec_times(xt)

    while(length(R) <= 1){

      xt <- series_H1(
        T      = T_val,
        trend  = true_params[["trend"]],
        par_H1 = c(
          shape = 1.20,
          scale = true_params[["rate"]]
        )
      )

      R <- rec_values(xt)
      L <- rec_times(xt)
    }

    m <- length(R)

    data_rec <- list(
      rec_values = R,
      rec_times  = L,
      time       = T_val
    )

    #_________________________________________________________________________
    # STEP 2: FIT H0 (EXPONENTIAL LDM)
    #_________________________________________________________________________

    lb_0 <- c(
      theta = 0.01,
      rate  = 0.01
    )

    ub_0 <- c(
      theta = min(R / L),
      rate  = 10
    )

    x0_0 <- c(
      theta = 0.10,
      rate  = 0.10
    )

    negLogLik_H0 <- function(params){
      -Likelihood_under_H0(data_rec, params)
    }

    fit_H0 <- NULL

    for(iterations in 1:50){

      fit_H0 <- tryCatch(
        nlminb(
          start     = x0_0,
          objective = negLogLik_H0,
          lower     = lb_0,
          upper     = ub_0
        ),
        error = function(e) NULL
      )

      if(!is.null(fit_H0) &&
         is.finite(fit_H0$objective) &&
         fit_H0$convergence == 0){
        break
      }

      x0_0 <- x0_0 + c(0.10, 0.10)
    }

    if(is.null(fit_H0) || fit_H0$convergence != 0){
      next
    }

    param_H0 <- c(
      fit_H0$par,
      logLik = -fit_H0$objective
    )

    #_________________________________________________________________________
    # STEP 3: FIT H1 (WEIBULL LDM)
    #_________________________________________________________________________

    lb_1 <- c(
      theta = 0.01,
      shape = 0.01,
      scale = 0.01
    )

    ub_1 <- c(
      theta = min(R / L),
      shape = 5,
      scale = 5
    )

    x0_1 <- c(
      theta = param_H0[["theta"]],
      shape = 1.20,
      scale = param_H0[["rate"]]
    )

    negLogLik_H1 <- function(params){
      -Likelihood_under_H1(data_rec, params)
    }

    fit_H1 <- NULL

    for(iterations in 1:70){

      fit_H1 <- tryCatch(
        nlminb(
          start     = x0_1,
          objective = negLogLik_H1,
          lower     = lb_1,
          upper     = ub_1
        ),
        error = function(e) NULL
      )

      if(!is.null(fit_H1) &&
         is.finite(fit_H1$objective) &&
         fit_H1$convergence == 0){
        break
      }

      x0_1 <- x0_1 + c(0.010, 0.10, 0.10)
    }

    if(is.null(fit_H1) || fit_H1$convergence != 0){
      next
    }

    param_H1 <- c(
      fit_H1$par,
      logLik = -fit_H1$objective
    )

    #_________________________________________________________________________
    # STEP 4: WILKS LIKELIHOOD RATIO TEST
    #_________________________________________________________________________

    logLik_H0 <- param_H0["logLik"]
    logLik_H1 <- param_H1["logLik"]

    LR <- 2 * (logLik_H1 - logLik_H0)

    ## Numerical safeguard
    LR <- max(0, LR)
    if(LR < 0) next;

    critical_value <- qchisq(
      p  = 1 - alpha,
      df = 1
    )

    reject <- as.integer(
      LR >= critical_value
    )

    #_________________________________________________________________________
    # STEP 5: STORE RESULTS
    #_________________________________________________________________________

    results_power[trial, ] <- c(
      param_H0["theta"],
      param_H0["rate"],
      logLik_H0,

      param_H1["theta"],
      param_H1["shape"],
      param_H1["scale"],
      logLik_H1,

      LR,
      reject,
      m
    )

    trial <- trial + 1
  }

  ###############################################################################
  # TYPE I ERROR ---------
  ###############################################################################

  type1_error <- mean(
    results$p_value < alpha,
    na.rm = TRUE
  )

  ###############################################################################
  # EMPIRICAL POWER
  ###############################################################################

  empirical_power <- mean(
    results_power$reject,
    na.rm = TRUE
  )

  ###############################################################################
  # REPORT
  ###############################################################################

  cat(
    "\nType I Error =",
    round(100 * type1_error, 2),
    "%\n"
  )

  cat(
    "Power =",
    round(100 * empirical_power, 2),
    "%\n"
  )

  ###############################################################################
  # SAVE FINAL RESULTS
  ###############################################################################

  final[kk, ] <- c(
    trend_val,
    100 * type1_error,
    100 * empirical_power,
    mean(results$nRecords, na.rm = TRUE)
  )

  kk <- kk + 1
}

###############################################################################
# Plot
###############################################################################
library(ggplot2)
library(scales)  # for percentage formatting if needed

# Ensure your data is a dataframe
final <- as.data.frame(final)

# Rescale Nb_record to match the scale of the left y-axis
scale_factor <- max(final$Power) / max(final$Nb_record)

ggplot(final, aes(x = trend)) +
  geom_line(aes(y = Type_I_Error, color = "Type I Error"), size = 1) +
  geom_line(aes(y = Power, color = "Power"), size = 1) +
  geom_line(aes(y = Nb_record * scale_factor, color = "Nb_record"), linetype = "dashed", size = 1) +

  # Left y-axis
  scale_y_continuous(
    name = "Test Power | Type I Error (%)",

    # Right y-axis
    sec.axis = sec_axis(~ . / scale_factor, name = "Number of Records")
  ) +

  scale_color_manual(values = c("Type I Error" = "lightblue", "Power" = "darkblue", "Nb_record" = "#2E8B57")) +

  labs(x = "Theta", color = "Metric") +
  theme_classic() +
  theme(
    axis.title.y.right = element_text(color = "#2E8B57"),
    axis.title.y.left = element_text(color = "black"),
    legend.position = "top"
  )
