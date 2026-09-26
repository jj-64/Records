#########------------------Monte Carlo Simulation ---------------######################
#' Monte Carlo Simulation of Record Model Identification via Test Permutations
#'
#' Runs multiple simulations of time series generated under a specified null model (e.g. DTRW, LDM, iid, YNM),
#' applies all 24 possible permutations of four record-based tests (DTRW, LDM, Classical, YNM),
#' and records which test (if any) is accepted by each permutation.
#'
#' @param n_sim Number of simulated series (default = 1000)
#' @param T_value Length of each time series (default = 50)
#' @param generator String. The name of series generator
#' @param series_args List. Arguments of the generator of the series, the
#' argument "T" is not included
#' @param n_arg Character. default is "T" as some generator functions take "n"
#' @param H0 Character. The true generating process: "dtrw", "ldm",
#' "iid", or "ynm"
#' @param alpha Numeric. Significance level for all tests (default = 0.05)
#' @param lag Integer. Lag parameter for Box-Jenkins test (default = 10)
#' @param warmup Integer. Warm-up parameter for YNM test (default = 2)
#' @param K Integer. minimum number of observation in bin for ynm test (default is NULL)
#' @param print logical default is FALSE, summary is not printed
#' @param obs_type String. "all" if data provided is the whole series \eqn{X_t} or
#' "records" if the underlying series is \eqn{R_n}. In this case, the parameter
#' record_times must be provided.
#' Forced in case "obs_type" = "records"
#' @param approximate Logical, if \code{TRUE} use the asymptotic normal approximation for DTRW \eqn{N_T}-test
#'   (default = \code{FALSE} for the exact quantile test).
#' @param one.sided Logical, if \code{TRUE} perform a one-sided test for DTRW \eqn{N_T}-test
#'   (default = \code{FALSE} for two-sided).
#' @return A list containing:
#' \describe{
#'   \item{results_all}{Data frame of all simulation × permutation outcomes}
#'   \item{summary_total}{Frequency of accepted decisions across all runs}
#'   \item{perm_dec_table}{Contingency Table of waht each permutation returns}
#'   \item{accuracy_by_perm}{How often each permutation recovered the true H0}
#' }
#' @examples
#' \dontrun{
#' simulation_permutation_analysis(n_sim = 1000, T_value = 50, generator = dtrw_series,
#' series_args = list(dist="cauchy",location=0, scale=1),
#' H0 = "dtrw",
#' obs_type = "all")
#' }
#' ===== Overall Accuracy =====
#' #' #' # [1] 0.957
#' #' #' #
#' # ===== Model Frequencies =====
#' #   majority_decision count proportion
#' # 1              dtrw   957      0.957
#' # 2               ldm     1      0.001
#' # 3                no    33      0.033
#' # 4               ynm     9      0.009
#' #
#' # ===== Stability Summary =====
#' #   mean_stability median_stability mean_entropy
#' # 1      0.9044998                1    0.1334375
#' @export
simulation_permutation_analysis <- function(
    n_sim = 1000,
    T_value = 50,
    generator,
    series_args = list(),
    n_arg = "T",
    H0 = c("dtrw", "ldm", "iid", "ynm"),
    alpha = 0.05,
    lag = 10,
    warmup = 2,
    obs_type = "all",
    approximate = FALSE,
    one.sided = FALSE,
    verbose = TRUE
) {

  H0 <- match.arg(H0)

  args <- series_args
  args[[n_arg]] <- T_value

  all_perm_results <- list()

  simulation_summary <- data.frame(
    sim_id = seq_len(n_sim),
    majority_decision = NA_character_,
    stability = NA_real_,
    #confidence = NA_real_,
    entropy = NA_real_,
    correct = NA,
    stringsAsFactors = FALSE
  )

  pb <- txtProgressBar(
    min = 0,
    max = n_sim,
    style = 3
  )

  for(i in seq_len(n_sim)) {

    X <- do.call(generator, args)

    result <- test_model_permutations(
      X = X,
      alpha = alpha,
      lag = lag,
      warmup = warmup,
      approximate = approximate,
      one.sided = one.sided,
      obs_type = obs_type,
      verbose = FALSE
    )

    sim_res <- result$permutation_results

    sim_res$sim_id <- i

    all_perm_results[[i]] <- sim_res

    simulation_summary$majority_decision[i] <-
      as.character(result$majority_decision)

    simulation_summary$stability[i] <-
      result$stability

    simulation_summary$entropy[i] <-
      result$entropy

    simulation_summary$correct[i] <-
      result$majority_decision == H0

    setTxtProgressBar(pb, i)

  }

  close(pb)

  all_perm_results <-
    dplyr::bind_rows(all_perm_results)

  ## --------------------------
  ## GLOBAL PERFORMANCE
  ## --------------------------

  model_frequency <-
    simulation_summary |>
    dplyr::count(
      majority_decision,
      name = "count"
    )

  model_frequency$proportion <-
    model_frequency$count / n_sim

  overall_accuracy <-
    mean(
      simulation_summary$correct,
      na.rm = TRUE
    )

  ## --------------------------
  ## CONFUSION MATRIX
  ## --------------------------

  confusion_matrix <-
    table(
      Truth = rep(H0, n_sim),
      Predicted =
        simulation_summary$majority_decision
    )

  ## --------------------------
  ## PERMUTATION STABILITY
  ## --------------------------

  permutation_stability <-
    all_perm_results |>
    dplyr::count(
      permutation,
      model
    ) |>
    tidyr::pivot_wider(
      names_from = model,
      values_from = n,
      values_fill = 0
    )

  ## --------------------------
  ## STABILITY METRICS
  ## --------------------------

  stability_summary <- data.frame(

    mean_stability =
      mean(
        simulation_summary$stability,
        na.rm = TRUE
      ),

    median_stability =
      median(
        simulation_summary$stability,
        na.rm = TRUE
      ),

    # mean_confidence =
    #   mean(
    #     simulation_summary$confidence,
    #     na.rm = TRUE
    #   ),

    mean_entropy =
      mean(
        simulation_summary$entropy,
        na.rm = TRUE
      )
  )

  if(verbose) {

    cat("\n===== Overall Accuracy =====\n")
    print(round(overall_accuracy,4))

    cat("\n===== Model Frequencies =====\n")
    print(model_frequency)

    cat("\n===== Stability Summary =====\n")
    print(stability_summary)

  }

  list(

    truth = H0,

    overall_accuracy =
      overall_accuracy,

    model_frequency =
      model_frequency,

    stability_summary =
      stability_summary,

    confusion_matrix =
      confusion_matrix,

    simulation_summary =
      simulation_summary,

    permutation_stability =
      permutation_stability,

    permutation_results =
      all_perm_results
  )
}

## Older Version -----------
#
# # sim_results <- Simulation_Permutation_Analysis(n_sim=2,
# # T=50,generator = DTRW_series, series_args =list(dist="cauchy",
# #  loc=0, scale=1),H0 = "dtrw", obs_type = "all)
#
# ### 75% of the permutations trees return "dtrw" and 25% return "ynm".
# ### On average, one simulation will return the following:
# #  summary_total
# #  decision Freq
# #  1     DTRW 0.75
# #  2      YNM 0.25
#
#
#  ### Contingency Table of what each permutation returns
# #  perm_dec_table
# #  Permutation DTRW YNM
# #  1         CDLY  1.0 0.0
# #  2         CDYL  1.0 0.0
# #  3         CLDY  1.0 0.0
# #  4         CLYD  0.5 0.5
# #  5         CYDL  0.5 0.5
# #  6         CYLD  1.0 0.0
# #  7         DCLY  1.0 0.0
# #  8         DCYL  1.0 0.0
# #  9         DLCY  1.0 0.0
# #  10        DLYC  1.0 0.0
# #  11        DYCL  0.5 0.5
# #  12        DYLC  0.5 0.5
# #  13        LCDY  0.5 0.5
# #  14        LCYD  0.5 0.5
# #  15        LDCY  0.5 0.5
# #  16        LDYC  1.0 0.0
# #  17        LYCD  1.0 0.0
# #  18        LYDC  0.5 0.5
# #  19        YCDL  0.5 0.5
# #  20        YCLD  0.5 0.5
# #  21        YDCL  0.5 0.5
# #  22        YDLC  0.5 0.5
# #  23        YLCD  1.0 0.0
# #  24        YLDC  1.0 0.0
#
#  ## Under HO: "dtrw", those permutations having accuracy 1.0 returned
#  # "dtrw" across all simulations
# #     Permutation Success_Rate
# #     1         CDLY          1.0
# #     2         CDYL          1.0
# #     3         CLDY          1.0
# #     4         CLYD          0.5
# #     5         CYDL          0.5
# #     6         CYLD          1.0
# #     7         DCLY          1.0
# #     8         DCYL          1.0
# #     9         DLCY          1.0
# #     10        DLYC          1.0
# #     11        DYCL          0.5
# #     12        DYLC          0.5
# #     13        LCDY          0.5
# #     14        LCYD          0.5
# #     15        LDCY          0.5
# #     16        LDYC          1.0
# #     17        LYCD          1.0
# #     18        LYDC          0.5
# #     19        YCDL          0.5
# #     20        YCLD          0.5
# #     21        YDCL          0.5
# #     22        YDLC          0.5
# #     23        YLCD          1.0
# #     24        YLDC          1.0
# Simulation_Permutation_Analysis <- function(
#     n_sim = 1000,
#     T = 50,
#     generator, ## function:the function generating the series
#     series_args=list(), ## arguments of the generator function other than "T" and the "param_name" we are simulating
#     n_arg="T",
#     H0 = c("dtrw", "ldm", "iid", "ynm"),
#     alpha = 0.05,
#     lag = 10,
#     warmup = 2,
#     print = FALSE,
#     obs_type="all",
#     approximate = FALSE,
#     one.sided = FALSE
# ) {
#   H0 <- match.arg(H0)
#   if (!requireNamespace("combinat", quietly = TRUE)) {
#     stop("Please install 'combinat' package first.")
#   }
#
#   ## --- Helper: Generate series under true model H0
#   args <- series_args
#   args[[n_arg]] <- T   # could be "T" or "n"
#
#   ## --- Get all 24 permutations
#   test_ids <- c("D", "L", "C", "Y")
#   perms <- sapply(combinat::permn(test_ids), paste, collapse = "")
#
#   ## --- Storage
#   all_results <- expand.grid(
#     sim_id = 1:n_sim,
#     Permutation = perms,
#     stringsAsFactors = FALSE
#   )
#   all_results$decision <- NA_character_
#
#   ## --- Simulation loop
#   pb <- txtProgressBar(min = 0, max = n_sim, style = 3)
#   for (i in seq_len(n_sim)) {
#     X <- do.call(generator, args)
#     record_times = rec_times(X)
#     perm_result <-  test_model_permutations(X, alpha = alpha, lag = lag, warmup = warmup, K=K,
#                                             print= print, obs_type = obs_type, record_times = record_times,
#                                             approximate = approxiamte, one.sided = one.sided)
#
#     # store results for this simulation
#     all_results$decision[all_results$sim_id == i] <- (perm_result$decision)$decision
#     summary_per_sim =  perm_result$summary
#     setTxtProgressBar(pb, i)
#   }
#   close(pb)
#
#   ############## ------------------------ SUMARIES --------------------
#   # 1) global counts of Decisions (exclude NA)
#   summary_total <- all_results %>%
#     filter(!is.na(decision)) %>%
#     count(decision, name = "Freq") %>%
#     arrange(desc(Freq))
#   summary_total = as.data.frame(summary_total)
#   summary_total$Freq = summary_total$Freq/(n_sim*length(perms))
#
#
#   # 2) contingency table: how many times each permutation returned each decision
#   perm_dec_table <- all_results %>%
#     count(Permutation, decision) %>%
#     pivot_wider(names_from = decision, values_from = n, values_fill = 0)
#   perm_dec_table = as.data.frame( perm_dec_table)
#   perm_dec_table[,-1] =  perm_dec_table[,-1]/(n_sim)
#
#   ## 3) --- Accuracy by permutation (fraction of runs where decision matches true H0)
#   accuracy_by_perm <- aggregate(
#     I(all_results$decision == H0) ~ Permutation, data = all_results, FUN = mean
#   )
#   names(accuracy_by_perm)[2] <- "Success_Rate"
#
#
#   if(print){
#     message("\n===== Overall Summary =====")
#     print(summary_total)
#     message("\n===== Contingency Table =====")
#     print( perm_dec_table)
#     message("\n===== Accuracy by Permutation =====")
#     print(accuracy_by_perm[order(-accuracy_by_perm$Success_Rate), ])
#
#   }
#
#   return(list(
#     results_all = all_results,
#     global_accuracy = summary_total,
#     perm_stability =  perm_dec_table,
#     accuracy_by_perm = accuracy_by_perm
#   ))
# }
