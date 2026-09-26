
####---------------------Sequential Permutation Testing-----------------------------###

#' Sequential Permutation Testing of Record Process Models
#'
#' Applies all 24 possible permutations of four record-based tests (D, L, C, Y)
#' on the same series `X`. Each permutation runs the tests sequentially until
#' one fails to reject (returns something other than "no" or NA). The name of
#' that test is recorded as the output of the permutation.
#'
#' @param X Numeric vector. The observed series.
#' @param alpha Numeric. Significance level for all tests (default = 0.05).
#' @param lag Integer. Lag parameter for Box-Jenkins test (default = 10).
#' @param warmup Integer. Warm-up parameter for YNM test (default = 2).
#' @param obs_type String. "all" if data provided is the whole series \eqn{X_t} or
#' "records" if the underlying series is \eqn{R_n}. In this case, the parameter
#' record_times must be provided.
#' @param record_times Numeric vector of the occurence times of records. (Default is NA).
#' Forced in case "obs_type" = "records"
#' @param approximate Logical, if \code{TRUE} use the asymptotic normal approximation for DTRW \eqn{N_T}-test
#'   (default = \code{FALSE} for the exact quantile test). Forced when obs_type = "records".
#' @param one.sided Logical, if \code{TRUE} perform a one-sided test for DTRW \eqn{N_T}-test
#'   (default = \code{FALSE} for two-sided). Forced when obs_type = "records".
#' @param print logical default is FALSE, summary is not printed
#' @return A list containing:
#'
#' 1 - decision: A data frame with 24 rows:
#' \describe{
#'   \item{majority_decision}{The most winning test according to all permutations}
#'   \item{stability}{The proportion for the most winning test among all permutations}
#'   \item{entropy}{entropy}
#'   \item{permutation_results}{dataframe of all permutaion results}
#'   \item{summary}{dataframe of summay count for each test}
#' }
#'
#' 2-  summary: A Summary of accepted models across permutations
#' @export
#' @examples
#' \dontrun{
#' result <- test_model_permutations(X = rnorm(50),
#' obs_type = "all",
#' verbose = FALSE,
#' lag = 10,
#' alpha = 0.05)
#' }
#' # result
#' # $majority_decision
#' # [1] iid
#' # Levels: iid
#' #
#' # $stability
#' # [1] 1
#' #
#' # $entropy
#' # [1] 0
#' #
#' # $permutation_results
#' # permutation accepted_test accepted_after model
#' # 1         DLCY             C              3   iid
#' # 2         DLYC             C              4   iid
#' # 3         DYLC             C              4   iid
#' # 4         YDLC             C              4   iid
#' # 5         YDCL             C              3   iid
#' # 6         DYCL             C              3   iid
#' # 7         DCYL             C              2   iid
#' # 8         DCLY             C              2   iid
#' # 9         CDLY             C              1   iid
#' # 10        CDYL             C              1   iid
#' # 11        CYDL             C              1   iid
#' # 12        YCDL             C              2   iid
#' # 13        YCLD             C              2   iid
#' # 14        CYLD             C              1   iid
#' # 15        CLYD             C              1   iid
#' # 16        CLDY             C              1   iid
#' # 17        LCDY             C              2   iid
#' # 18        LCYD             C              2   iid
#' # 19        LYCD             C              3   iid
#' # 20        YLCD             C              3   iid
#' # 21        YLDC             C              4   iid
#' # 22        LYDC             C              4   iid
#' # 23        LDYC             C              4   iid
#' # 24        LDCY             C              3   iid
#' #
#' # $summary
#' # model count proportion
#' # 1   iid    24          1
test_model_permutations <- function(
    X,
    alpha = 0.05,
    lag = 10,
    warmup = 2,
    K = NULL,
    verbose = TRUE,
    obs_type = c("all", "records"),
    record_times = NULL,
    approximate = FALSE,
    one.sided = FALSE
) {
  obs_type <- match.arg(obs_type)

  # Ensure combinat is available
  if (!requireNamespace("combinat", quietly = TRUE)) {
    stop("Please install the 'combinat' package with install.packages('combinat').")
  }

  ## ---------- Helper: run one test safely ----------
  ## Test dispatcher
  if (obs_type == "all") {

    test_functions <- list(

      D = function()
        test_dtrw_increment(
          X,
          alpha = alpha
        )$decision,

      L = function()
        test_ldm_sequential(
          X,
          alpha = alpha
        )$decision,

      C = function()
        test_iid_serial_independence(
          X,
          alpha = alpha,
          lags = lag
        )$decision,

      Y = function()
        test_ynm_rec_gap(
          X,
          alpha = alpha,
          warmup = warmup,
          K = K
        )$decision
    )

  } else {

    test_functions <- list(

      D = function()
        test_dtrw_rec_count(
          X,
          alpha = alpha,
          approximate = approximate,
          one.sided = one.sided
        )$decision,

      L = function()
        test_ldm_sequential(
          X,
          alpha = alpha
        )$decision,

      C = function()
        test_iid_rec_count(
          X,
          alpha = alpha
        )$decision,

      Y = function()
        test_ynm_rec_count(
          X,
          alpha = alpha
        )$decision
    )
  }

  ## ---------- Generate 24 permutations ----------

  test_ids <- names(test_functions)

  perms <- sapply(
    combinat::permn(test_ids),
    paste,
    collapse = ""
  )

  ## ---------- Prepare result storage ----------
  results <- data.frame(
    permutation = perms,
    accepted_test = NA_character_,
    accepted_after = NA_integer_,
    model = NA_character_,
    stringsAsFactors = FALSE
  )

  ## ---------- Test Acceptance Matrix
  # pass_matrix <- matrix(
  #   FALSE,
  #   nrow = length(perms),
  #   ncol = 4,
  #   dimnames = list(NULL, c("dtrw","ldm","iid","ynm"))
  # )
  ## ---------- Run all permutations ----------

  for (i in seq_along(perms)) {

    order_tests <- strsplit(
      perms[i],
      ""
    )[[1]]

    accepted_model <- "no"
    accepted_test <- NA
    n_checked <- length(order_tests)

    for (j in seq_along(order_tests)) {

      test_id <- order_tests[j]

      decision <- tryCatch(
        test_functions[[test_id]](),
        error = function(e) NA
      )

      if (!is.na(decision) &&
          tolower(decision) != "no") {

        accepted_model <- tolower(decision)
        accepted_test <- test_id
        n_checked <- j

        break
      }
    }

    results$accepted_test[i] <- accepted_test
    results$accepted_after[i] <- n_checked
    results$model[i] <- accepted_model
    # pass_matrix[i,"dtrw"] <- dtrw_accept
    # pass_matrix[i,"ldm"] <- ldm_accept
    # pass_matrix[i,"iid"] <- dtrw_accept
    # pass_matrix[i,"ynm"] <- ldm_accept
  }

  ## ---------- Optional: summary of how many times each model accepted ----------
  summary <- as.data.frame(table(results$model) )

  names(summary) <- c(
    "model",
    "count"
  )

  summary$proportion <-
    round(
      summary$count /
        sum(summary$count),
      4
    )

  summary <-
    summary[
      order(
        summary$count,
        decreasing = TRUE
      ),
    ]

  if (verbose) {

    message(
      "Model selection frequencies across all permutations:"
    )

    print(summary)
  }

  ## ------------------------
  ## Majority decision
  ## ------------------------

  majority_decision <- summary$model[1]

  stability <- summary$proportion[1]

  ## Entropy of decisions
  entropy <- -sum(
    summary$proportion *
      log(summary$proportion),
    na.rm = TRUE
  )

  return(
    list(
      majority_decision = majority_decision,
      stability = stability,
      entropy = entropy,
      permutation_results = results,
      summary = summary
      #acceptance_matrix = colMeans(pass_matrix)
    )
  )

}

######---------Old Function ------------###############
# Test_Permutation_OLD = function(x, sig=0.05){
#
#   ## Tree: dataframe of nrow=simulations, returning the decision of each tree[column] at each simulation[row]
#   tree = as.data.frame(matrix(0,ncol=24,nrow=1))
#   colnames(tree) = c("CLYD","CLDY","CYLD","CYDL","CDYL","CDLY",
#                      "LCYD","LCDY","LYCD","LYDC","LDYC","LDCY",
#                      "YCLD","YCDL","YLCD","YLDC","YDLC","YDCL",
#                      "DCLY","DCYL","DYCL","DYLC","DLYC","DLCY"
#   )
#
#   ## Result for the average: what the whole sequence of trees is voting for (simulation[row] and model[col])
#   result=data.frame(matrix(0,nrow=1,ncol=5))
#   colnames(result)=c("None","iid","dtrw","ldm","YN")
#
#   ############################ Perform simulation ################################
#
#   tree[1,1] = DT_CLYD(X=x, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL)
#   tree[1,2] = DT_CLDY(X=x, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL)
#   tree[1,3] = DT_CYLD(X=x, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL)
#   tree[1,4] = DT_CYDL(X=x, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL)
#   tree[1,5] = DT_CDYL(X=x, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL)
#   tree[1,6] = DT_CDLY(X=x, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL)
#
#   tree[1,7] = DT_LCYD(X=x, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL)
#   tree[1,8] = DT_LCDY(X=x, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL)
#   tree[1,9] = DT_LYCD(X=x, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL)
#   tree[1,10] = DT_LYDC(X=x, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL)
#   tree[1,11] = DT_LDYC(X=x, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL)
#   tree[1,12] = DT_LDCY(X=x, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL)
#
#   tree[1,13] = DT_YCLD(X=x, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL)
#   tree[1,14] = DT_YCDL(X=x, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL)
#   tree[1,15] = DT_YLCD(X=x, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL)
#   tree[1,16] = DT_YLDC(X=x, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL)
#   tree[1,17] = DT_YDLC(X=x, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL)
#   tree[1,18] = DT_YDCL(X=x, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL)
#
#   tree[1,19] = DT_DCLY(X=x, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL)
#   tree[1,20] = DT_DCYL(X=x, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL)
#   tree[1,21] = DT_DYCL(X=x, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL)
#   tree[1,22] = DT_DYLC(X=x, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL)
#   tree[1,23] = DT_DLYC(X=x, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL)
#   tree[1,24] = DT_DLCY(X=x, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL)
#
#   ## Summary
#   ## Results
#   result[1,"None"] = sum(tree[1,]==0)
#   result[1,"iid"] = sum(tree[1,]==1)
#   result[1,"dtrw"] = sum(tree[1,]==2)
#   result[1,"ldm"] = sum(tree[1,]==3)
#   result[1,"YN"] = sum(tree[1,]==4)
#
#   return(list("Tree"=tree,"summary" = result))
# }


##########--------------- Sequential Limited obs_type --------------########

# Sequential Permutation Testing of Record Process Models
#
# Applies all 24 possible permutations of four record-based tests (D, L, C, Y)
# on limited information of records and their occurence times.
# Each permutation runs the tests sequentially until
# one fails to reject (returns something other than "no" or NA). The name of
# that test is recorded as the output of the permutation.
#
# @param X Numeric vector. The observed series.
# @param alpha Numeric. Significance level for all tests (default = 0.05).
# @param print logical default is FALSE, summary is not printed
# @return A list containing:
#
# 1 - decision: A data frame with 24 rows:
# \describe{
#   \item{Permutation}{Order of tests applied (e.g., "DLCY")}
#   \item{Accepted_Test}{The first test that did not reject ("no")}
#   \item{decision}{Either the accepting test name or "All_Rejected"}
# }
#
#
# 2-  summary: A Summary of accepted models across permutations
# @export
# @examples
# \dontrun{
# result <- Test_Permutation_limited(X = rnorm(50), print = FALSE)
#
# print(result$decision)
# #  Permutation Accepted_Test  decision
# # 1         DLCY             C Classical
# # 2         DLYC             Y       YNM
# # 3         DYLC             Y       YNM
# # 4         YDLC             Y       YNM
# # 5         YDCL             Y       YNM
# # 6         DYCL             Y       YNM
# # 7         DCYL             C Classical
# # 8         DCLY             C Classical
# # 9         CDLY             C Classical
# # 10        CDYL             C Classical
# # 11        CYDL             C Classical
# # 12        YCDL             Y       YNM
# # 13        YCLD             Y       YNM
# # 14        CYLD             C Classical
# # 15        CLYD             C Classical
# # 16        CLDY             C Classical
# # 17        LCDY             C Classical
# # 18        LCYD             C Classical
# # 19        LYCD             Y       YNM
# # 20        YLCD             Y       YNM
# # 21        YLDC             Y       YNM
# # 22        LYDC             Y       YNM
# # 23        LDYC             Y       YNM
# # 24        LDCY             C Classical
#
# print(result$summary)
# #       Var1 Freq
# #        1 Classical   12
# #        2       YNM   12
# }
# Test_Permutation_limited <- function(X, alpha = 0.05, print = TRUE) {
#   # Ensure combinat is available
#   if (!requireNamespace("combinat", quietly = TRUE)) {
#     stop("Please install the 'combinat' package with install.packages('combinat').")
#   }
#
#   ## ---------- Helper: run one test safely ----------
#   run_test <- function(test_id, X, alpha) {
#     res <- switch(test_id,
#                   D = tryCatch(test_dtrw_rec_count(X, alpha = alpha)$decision, error = function(e) NA),
#                   L = tryCatch(test_ldm_sequential(X, alpha = alpha)$decision, error = function(e) NA),
#                   C = tryCatch(test_iid_rec_count(X, alpha = alpha)$decision, error = function(e) NA),
#                   Y = tryCatch(test_ynm_rec_count(X, alpha = alpha)$decision, error = function(e) NA),
#                   NA)
#     if (is.null(res) || length(res) == 0) res <- NA
#     return(res)
#   }
#
#   ## ---------- Generate 24 permutations ----------
#   test_ids <- c("D", "L", "C", "Y")
#   perms_list <- combinat::permn(test_ids)
#   perms <- sapply(perms_list, paste, collapse = "")
#
#   ## ---------- Prepare result storage ----------
#   results <- data.frame(
#     Permutation = perms,
#     Accepted_Test = NA_character_,
#     decision = NA_character_,
#     stringsAsFactors = FALSE
#   )
#
#   ## ---------- Run all permutations ----------
#   for (i in seq_along(perms)) {
#     order <- strsplit(perms[i], "")[[1]]
#     accepted <- NA
#     decision <- "All_Rejected"
#
#     for (test in order) {
#       dec <- run_test(test, X=X, alpha=alpha)
#       if (!is.na(dec) && dec != "no") {
#         accepted <- test
#         decision <- dec
#         break
#       }
#     }
#
#     results$Accepted_Test[i] <- accepted
#     results$decision[i] <- decision
#   }
#
#   ## ---------- Optional: summary of how many times each model accepted ----------
#   summary <- as.data.frame(table(results$decision))
#   summary <- summary[order(summary$Freq, decreasing = TRUE), ]
#   colnames(summary) = c("decision", "Freq")
#
#
#   if (print) {
#     message("Summary of accepted models across permutations:")
#     print(summary)
#   }
#   return(list(decision = results, summary = summary))
# }
