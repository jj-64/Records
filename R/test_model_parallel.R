########---------Parallel Testing ------------- #############
#' Parallel Testing of Record Process Models
#'
#' Applies all 10 possible Tests (if applicable) of four record-based tests (D, L, C, Y)
#' on the same series `X` independently.
#'
#' @param X Numeric vector. The observed series.
#' @param obs_type String. "all" or "records" depending on the series provided.
#' @param alpha Numeric. Significance level for all tests (default = 0.05).
#' @param lag Integer. Lag parameter for Box-Jenkins test (default = 10).
#' @param warmup Integer. Warm-up parameter for YNM test (default = 2).
#' @param approximate Logical, if \code{TRUE} use the asymptotic normal approximation
#'   (default = \code{FALSE} for the exact quantile test).
#' @param one.sided Logical, if \code{TRUE} perform a one-sided test
#'   (default = \code{FALSE} for two-sided).
#' @param method Character, p-value combination method: one of
#'   \code{"Bonf"}, \code{"Holm"}, \code{"Sidak"}, or \code{"Chisq"} (default = "Bonf").
#' @param K Optional. Number of partitions in \link{test_ynm_chisq} test.
#'    If given, force exactly K partitions using quantiles.
#' @param estimate_gamma Logical. If to estimate \eqn{\gamma} through minimizing \eqn{\chi^2} (Default = TRUE)
#' @param gamma Numeric. Optional. Force if estimated_gamma = FALSE.
#' @param RSq Numeric, minimum adjusted R-squared required to accept the LDM
#'   hypothesis (default = 0.8).
#' @param record_times (Default = NA) record times in case obs_type = "records"
#' @return A list containing all:
#' \describe{
#'   \item{results}{List of all Test outcomes}
#'   \item{decision}{Final decision of each test}
#' }
#' @export
#' @examples
#' tests = Test_Parallel(X=rnorm(50), obs_type = "all", alpha = 0.05)
#' tests$decision
#' # iid_NT.decision        DTRW_NT.decision
#' # "iid"                  "dtrw"
#' # YNM_NT.decision       YNM_Pearson.decision
#' #           "no"                      NA
#' # YNM_Geom.decision         LDM_NT.decision
#' #            "no"                    "no"
#' # LDM_Sequential.decision     LDM_Regression.decision
#' #     "no"                           "no"
#' # iid_Box.decision     DTRW_Indep.decision
#' #  "iid"             "no"
test_model_parallel <- function(
    X,
    obs_type = c("all", "records"),
    record_times = NULL,
    alpha = 0.05,
    lag = 10,
    warmup = 2,
    approximate = FALSE,
    one.sided = FALSE,
    method = "Holm",
    K = NULL,
    gamma = NULL
) {

  obs_type <- match.arg(obs_type)

  results <- list()

  estimate_gamma = ifelse(is.null(gamma) || is.na(gamma), TRUE, FALSE)

  ## ------------------------------------------------------
  ## Record-statistic tests
  ## ------------------------------------------------------
  if (obs_type == "records") {

    if (is.null(record_times))
      stop("'record_times' must be supplied.")

    if (length(record_times) != length(X))
      stop("'record_times' and 'X' must have the same length.")

    }

  results$iid_rec_count <- test_iid_rec_count(
    X = X,
    alpha = alpha
  )

  results$dtrw_rec_count <- test_dtrw_rec_count(
    X = X,
    alpha = alpha,
    approximate = approximate,
    one.sided = one.sided
  )

  results$ynm_rec_count <- test_ynm_rec_count(
    X = X,
    gamma = gamma,
    alpha = alpha
  )

  results$ynm_chisq <- test_ynm_chisq(
    X = X,
    gamma = gamma,
    K = K,
    estimated = estimate_gamma,
    alpha = alpha
  )

  results$ynm_rec_gaps <- test_ynm_rec_gaps(
    X = X,
    alpha = alpha,
    K = K,
    warmup = warmup,
    obs_type = obs_type,
    record_times = record_times
  )

  results$ldm_reC_count <- test_ldm_rec_count(
    X = X,
    alpha = alpha
  )

  results$ldm_sequential <- test_ldm_sequential(
    X = X,
    alpha = alpha
  )

  results$ldm_tred <- test_ldm_trend(
    X = X,
    alpha = alpha
  )

  ## ------------------------------------------------------
  ## Full-series tests
  ## ------------------------------------------------------

  if (obs_type == "all") {

    results$IID_Independence <-
      test_serial_independence(
        X = X,
        alpha = alpha,
        lags = lag
      )

    results$DTRW_Increment <-
      test_dtrw_increment(
        X = X,
        alpha = alpha,
        method = method
      )
  }

  ## ------------------------------------------------------
  ## Extract decisions
  ## ------------------------------------------------------

  decision_table <- data.frame(
    test = names(results),
    decision = sapply(
      results,
      function(x) {

        if (!is.null(x$decision))
          toupper(as.character(x$decision))
        else
          NA_character_

      }
    ),
    stringsAsFactors = FALSE
  )

  ## ------------------------------------------------------
  ## Model vote count
  ## ------------------------------------------------------

  valid_decisions <- decision_table$decision[
    !is.na(decision_table$decision)
  ]

  model_scores <- as.data.frame(
    table(valid_decisions)
  )

  names(model_scores) <- c(
    "model",
    "count"
  )

  model_scores$proportion <-
    model_scores$count /
    sum(model_scores$count)

  model_scores <-
    model_scores[
      order(
        model_scores$count,
        decreasing = TRUE
      ),
    ]

  ## ------------------------------------------------------
  ## Final recommendation
  ## ------------------------------------------------------

  best_model <- model_scores$model[1]

  agreement <-
    model_scores$proportion[1]

  confidence <-
    agreement -
    ifelse(
      nrow(model_scores) > 1,
      model_scores$proportion[2],
      0
    )

  final_decision <- ifelse(
    agreement >= 0.60,
    best_model,
    "UNCERTAIN"
  )

  ## ------------------------------------------------------
  ## Test summary
  ## ------------------------------------------------------

  summary <- list(

    final_decision = final_decision,

    leading_model = best_model,

    agreement = agreement,

    confidence = confidence,

    model_scores = model_scores

  )

  return(

    list(

      summary = summary,

      decisions = decision_table,

      results = results

    )

  )

}

## older version ------------
Test_Parallel <- function(X, obs_type = c("all","records") , record_times = NA,
                          alpha = 0.05, lag = 10, warmup = 2, approximate = FALSE,
                          one.sided = FALSE, method="Bonf",
                          K= NULL, estimate_gamma = TRUE, gamma = NULL, RSq = 0.8) {
  results = list()
  obs_type <- match.arg(obs_type)

  if (obs_type == "all" | obs_type == "records"){
    ## iid_NT
    results$"iid_NT" = test_iid_rec_count(X=X, alpha= alpha)

    ## DTRW_NT
    results$"DTRW_NT" = test_dtrw_rec_count(X=X, alpha = alpha, approximate = approximate, one.sided = one.sided)

    ## YNM_NT
    results$"YNM_NT" = test_ynm_rec_count(X= X, gamma = NA, alpha = alpha)

    ## YNM_Pearson
    results$"YNM_Pearson" = test_ynm_chisq(X=X, Partition = NA, gamma = NULL, K=K, estimated = estimate_gamma, alpha = alpha)

    ##YNM_Geom
    results$"YNM_Geom" = test_ynm_record_gap(X = X, alpha=alpha, K=K, warmup=warmup, record_times= record_times)

    ##  LDM_NT
    results$"LDM_NT" = test_ldm_rec_count(X = X, alpha = alpha)

    ##LDM_Sequential
    results$"LDM_Sequential" = test_ldm_sequential(X=X, alpha = alpha)

    ##LDM_Regression
    results$"LDM_Regression" = test_ldm_trend(X=X, alpha = alpha, RSq = RSq)
  }

  if (obs_type == "all"){
    ## iid_Box
    results$"iid_Box" = test_iid_serial_independence(X=X, alpha= alpha, lags = lag)

    ## DTRW_Indep
    results$"DTRW_Indep" = test_dtrw_increment(X = X, alpha= alpha, method=method)
  }

  # Filter elements whose names end with ".decision"
  unlisted = unlist(results)
  decision_items <- unlisted[grepl("\\.decision$", names(unlisted))]

  return(list(results = results, decision= decision_items))
}

