#' Sequential Permutation Testing of Record Process Models
#'
#' Applies all 24 possible permutations of four record-based tests (D, L, C, Y)
#' on the same series `X`. Each permutation runs the tests sequentially until
#' one fails to reject (returns something other than "NO" or NA). The name of
#' that test is recorded as the output of the permutation.
#'
#' @param X Numeric vector. The observed series.
#' @param alpha Numeric. Significance level for all tests (default = 0.05).
#' @param lag Integer. Lag parameter for Box-Jenkins test (default = 10).
#' @param warmup Integer. Warm-up parameter for YNM test (default = 2).
#' @return A list containing:
#'
#' 1 - decision: A data frame with 24 rows:
#' \describe{
#'   \item{Permutation}{Order of tests applied (e.g., "DLCY")}
#'   \item{Accepted_Test}{The first test that did not reject ("NO")}
#'   \item{Decision}{Either the accepting test name or "All_Rejected"}
#' }
#'
#'
#' 2-  summary: A Summary of accepted models across permutations
#' @export
#' @examples
#' \dontrun{
#' result <- Test_Permutation(X = rnorm(50), print = FALSE)
#'
#' print(result$decision)
#' #  Permutation Accepted_Test  Decision
#' # 1         DLCY             C Classical
#' # 2         DLYC             Y       YNM
#' # 3         DYLC             Y       YNM
#' # 4         YDLC             Y       YNM
#' # 5         YDCL             Y       YNM
#' # 6         DYCL             Y       YNM
#' # 7         DCYL             C Classical
#' # 8         DCLY             C Classical
#' # 9         CDLY             C Classical
#' # 10        CDYL             C Classical
#' # 11        CYDL             C Classical
#' # 12        YCDL             Y       YNM
#' # 13        YCLD             Y       YNM
#' # 14        CYLD             C Classical
#' # 15        CLYD             C Classical
#' # 16        CLDY             C Classical
#' # 17        LCDY             C Classical
#' # 18        LCYD             C Classical
#' # 19        LYCD             Y       YNM
#' # 20        YLCD             Y       YNM
#' # 21        YLDC             Y       YNM
#' # 22        LYDC             Y       YNM
#' # 23        LDYC             Y       YNM
#' # 24        LDCY             C Classical
#'
#' print(result$summary)
#' #       Var1 Freq
#' #        1 Classical   12
#' #        2       YNM   12
#' }
Test_Permutation <- function(X, alpha = 0.05, lag = 10, warmup = 2, print = TRUE) {
  # Ensure combinat is available
  if (!requireNamespace("combinat", quietly = TRUE)) {
    stop("Please install the 'combinat' package with install.packages('combinat').")
  }

  ## ---------- Helper: run one test safely ----------
  run_test <- function(test_id, X, alpha) {
    res <- switch(test_id,
                  D = tryCatch(Test_DTRW_Indep(X, alpha = alpha)$decision, error = function(e) NA),
                  L = tryCatch(Test_LDM_Sequential(X, alpha = alpha)$decision, error = function(e) NA),
                  C = tryCatch(Test_iid_BoxJenkins(X, alpha = alpha, lag = lag)$decision, error = function(e) NA),
                  Y = tryCatch(Test_YNM_Geom(X, alpha = alpha, warmup = warmup)$decision, error = function(e) NA),
                  NA)
    if (is.null(res) || length(res) == 0) res <- NA
    return(res)
  }

  ## ---------- Generate 24 permutations ----------
  test_ids <- c("D", "L", "C", "Y")
  perms_list <- combinat::permn(test_ids)
  perms <- sapply(perms_list, paste, collapse = "")

  ## ---------- Prepare result storage ----------
  results <- data.frame(
    Permutation = perms,
    Accepted_Test = NA_character_,
    Decision = NA_character_,
    stringsAsFactors = FALSE
  )

  ## ---------- Run all permutations ----------
  for (i in seq_along(perms)) {
    order <- strsplit(perms[i], "")[[1]]
    accepted <- NA
    decision <- "All_Rejected"

    for (test in order) {
      dec <- run_test(test, X=X, alpha=alpha)
      if (!is.na(dec) && dec != "NO") {
        accepted <- test
        decision <- dec
        break
      }
    }

    results$Accepted_Test[i] <- accepted
    results$Decision[i] <- decision
  }

  ## ---------- Optional: summary of how many times each model accepted ----------
  summary <- as.data.frame(table(results$Decision))
  summary <- summary[order(summary$Freq, decreasing = TRUE), ]
 colnames(summary) = c("decision", "Freq")


  if (print) {
    message("Summary of accepted models across permutations:")
    print(summary)
        }
  return(list(decision = results, summary = summary))
}



#############################################
Forest = function(x, sig=0.05){

## Tree: dataframe of nrow=simulations, returning the decision of each tree[column] at each simulation[row]
tree = as.data.frame(matrix(0,ncol=24,nrow=1))
colnames(tree) = c("CLYD","CLDY","CYLD","CYDL","CDYL","CDLY",
                   "LCYD","LCDY","LYCD","LYDC","LDYC","LDCY",
                   "YCLD","YCDL","YLCD","YLDC","YDLC","YDCL",
                   "DCLY","DCYL","DYCL","DYLC","DLYC","DLCY"
                      )

## Result for the average: what the whole sequence of trees is voting for (simulation[row] and model[col])
result=data.frame(matrix(0,nrow=1,ncol=5))
colnames(result)=c("None","Classical","DTRW","LDM","YN")

############################ Perform simulation ################################

  tree[1,1] = DT_CLYD(X=x, p =sig)
  tree[1,2] = DT_CLDY(X=x, p =sig)
  tree[1,3] = DT_CYLD(X=x, p =sig)
  tree[1,4] = DT_CYDL(X=x, p =sig)
  tree[1,5] = DT_CDYL(X=x, p =sig)
  tree[1,6] = DT_CDLY(X=x, p =sig)

  tree[1,7] = DT_LCYD(X=x, p =sig)
  tree[1,8] = DT_LCDY(X=x, p =sig)
  tree[1,9] = DT_LYCD(X=x, p =sig)
  tree[1,10] = DT_LYDC(X=x, p =sig)
  tree[1,11] = DT_LDYC(X=x, p =sig)
  tree[1,12] = DT_LDCY(X=x, p =sig)

  tree[1,13] = DT_YCLD(X=x, p =sig)
  tree[1,14] = DT_YCDL(X=x, p =sig)
  tree[1,15] = DT_YLCD(X=x, p =sig)
  tree[1,16] = DT_YLDC(X=x, p =sig)
  tree[1,17] = DT_YDLC(X=x, p =sig)
  tree[1,18] = DT_YDCL(X=x, p =sig)

  tree[1,19] = DT_DCLY(X=x, p =sig)
  tree[1,20] = DT_DCYL(X=x, p =sig)
  tree[1,21] = DT_DYCL(X=x, p =sig)
  tree[1,22] = DT_DYLC(X=x, p =sig)
  tree[1,23] = DT_DLYC(X=x, p =sig)
  tree[1,24] = DT_DLCY(X=x, p =sig)

## Summary
  ## Results
  result[1,"None"] = sum(tree[1,]==0)
  result[1,"Classical"] = sum(tree[1,]==1)
  result[1,"DTRW"] = sum(tree[1,]==2)
  result[1,"LDM"] = sum(tree[1,]==3)
  result[1,"YN"] = sum(tree[1,]==4)

  return(list("Tree"=tree,"summary" = result))
}

##########################################################
#' Monte Carlo Simulation of Record Model Identification via Test Permutations
#'
#' Runs multiple simulations of time series generated under a specified null model (e.g. DTRW, LDM, iid, YNM),
#' applies all 24 possible permutations of four record-based tests (DTRW, LDM, Classical, YNM),
#' and records which test (if any) is accepted by each permutation.
#'
#' @param n_sim Number of simulated series (default = 1000)
#' @param T Length of each time series (default = 50)
#' @param generator String. The name of series generator
#' @param series_args List. Arguments of the generator of the series, the
#' argument "T" is not included
#' @param n_arg Character. default is "T" as some generator functions take "n"
#' @param H0 Character. The true generating process: "DTRW", "LDM",
#' "Classical", or "YNM"
#' @param alpha Numeric. Significance level for all tests (default = 0.05)
#' @param lag Integer. Lag parameter for Box-Jenkins test (default = 10)
#' @param warmup Integer. Warm-up parameter for YNM test (default = 2)
#' @param print logical default is FALSE, summary is not printed
#' @return A list containing:
#' \describe{
#'   \item{results_all}{Data frame of all simulation × permutation outcomes}
#'   \item{summary_total}{Frequency of accepted decisions across all runs}
#'   \item{perm_dec_table}{Contingency Table of waht each permutation returns}
#'   \item{accuracy_by_perm}{How often each permutation recovered the true H0}
#' }
#' @export
#' @examples
#' \dontrun{
#' sim_results <- Simulation_Permutation_Analysis(n_sim=2, T=50,
#' generator = DTRW_series, series_args = list(dist="cauchy",loc=0, scale=1),
#' H0 = "DTRW")
#'
#'
#' ### 75% of the permutations trees return "DTRW" and 25% return "YNM".
#' ### On average, one simulation will return the following:
#' #  summary_total
#' #  Decision Freq
#' #  1     DTRW 0.75
#' #  2      YNM 0.25
#'
#'
#'  ### Contingency Table of what each permutation returns
#' #  perm_dec_table
#' #  Permutation DTRW YNM
#' #  1         CDLY  1.0 0.0
#' #  2         CDYL  1.0 0.0
#' #  3         CLDY  1.0 0.0
#' #  4         CLYD  0.5 0.5
#' #  5         CYDL  0.5 0.5
#' #  6         CYLD  1.0 0.0
#' #  7         DCLY  1.0 0.0
#' #  8         DCYL  1.0 0.0
#' #  9         DLCY  1.0 0.0
#' #  10        DLYC  1.0 0.0
#' #  11        DYCL  0.5 0.5
#' #  12        DYLC  0.5 0.5
#' #  13        LCDY  0.5 0.5
#' #  14        LCYD  0.5 0.5
#' #  15        LDCY  0.5 0.5
#' #  16        LDYC  1.0 0.0
#' #  17        LYCD  1.0 0.0
#' #  18        LYDC  0.5 0.5
#' #  19        YCDL  0.5 0.5
#' #  20        YCLD  0.5 0.5
#' #  21        YDCL  0.5 0.5
#' #  22        YDLC  0.5 0.5
#' #  23        YLCD  1.0 0.0
#' #  24        YLDC  1.0 0.0
#'
#'
#'  ## Under HO: "DTRW", those permutations having accuracy 1.0 returned "DTRW" across all simulations
#' #     Permutation Success_Rate
#' #     1         CDLY          1.0
#' #     2         CDYL          1.0
#' #     3         CLDY          1.0
#' #     4         CLYD          0.5
#' #     5         CYDL          0.5
#' #     6         CYLD          1.0
#' #     7         DCLY          1.0
#' #     8         DCYL          1.0
#' #     9         DLCY          1.0
#' #     10        DLYC          1.0
#' #     11        DYCL          0.5
#' #     12        DYLC          0.5
#' #     13        LCDY          0.5
#' #     14        LCYD          0.5
#' #     15        LDCY          0.5
#' #     16        LDYC          1.0
#' #     17        LYCD          1.0
#' #     18        LYDC          0.5
#' #     19        YCDL          0.5
#' #     20        YCLD          0.5
#' #     21        YDCL          0.5
#' #     22        YDLC          0.5
#' #     23        YLCD          1.0
#' #     24        YLDC          1.0
#' }
Simulation_Permutation_Analysis <- function(
    n_sim = 1000,
    T = 50,
    generator, ## function:the function generating the series
    series_args=list(), ## arguments of the generator function other than "T" and the "param_name" we are simulating
    n_arg="T",
    H0 = c("DTRW", "LDM", "Classical", "YNM"),
    alpha = 0.05,
    lag = 10,
    warmup = 2,
    print = FALSE
) {
  H0 <- match.arg(H0)
  if (!requireNamespace("combinat", quietly = TRUE)) {
    stop("Please install 'combinat' package first.")
  }

  ## --- Helper: Generate series under true model H0
  args <- series_args
  args[[n_arg]] <- T   # could be "T" or "n"

  ## --- Get all 24 permutations
  test_ids <- c("D", "L", "C", "Y")
  perms <- sapply(combinat::permn(test_ids), paste, collapse = "")

  ## --- Storage
  all_results <- expand.grid(
    sim_id = 1:n_sim,
    Permutation = perms,
    stringsAsFactors = FALSE
  )
  all_results$Decision <- NA_character_

  ## --- Simulation loop
  pb <- txtProgressBar(min = 0, max = n_sim, style = 3)
  for (i in seq_len(n_sim)) {
    X <- do.call(generator, args)
    perm_result <-  Test_Permutation(X, alpha = alpha, lag = lag, warmup = warmup, print= print)

    # store results for this simulation
    all_results$Decision[all_results$sim_id == i] <- (perm_result$decision)$Decision
    summary_per_sim =  perm_result$summary
    setTxtProgressBar(pb, i)
  }
  close(pb)

  ############## ------------------------ SUMARIES --------------------#################
  # 1) global counts of Decisions (exclude NA)
  summary_total <- all_results %>%
    filter(!is.na(Decision)) %>%
    count(Decision, name = "Freq") %>%
    arrange(desc(Freq))
  summary_total = as.data.frame(summary_total)
  summary_total$Freq = summary_total$Freq/(n_sim*length(perms))


  # 2) contingency table: how many times each permutation returned each Decision
  perm_dec_table <- all_results %>%
    count(Permutation, Decision) %>%
    pivot_wider(names_from = Decision, values_from = n, values_fill = 0)
  perm_dec_table = as.data.frame( perm_dec_table)
  perm_dec_table[,-1] =  perm_dec_table[,-1]/(n_sim)

  ## 3) --- Accuracy by permutation (fraction of runs where Decision matches true H0)
  accuracy_by_perm <- aggregate(
    I(all_results$Decision == H0) ~ Permutation, data = all_results, FUN = mean
  )
  names(accuracy_by_perm)[2] <- "Success_Rate"


  if(print){
    message("\n===== Overall Summary =====")
    print(summary_total)
    message("\n===== Contingency Table =====")
    print( perm_dec_table)
    message("\n===== Accuracy by Permutation =====")
    print(accuracy_by_perm[order(-accuracy_by_perm$Success_Rate), ])

  }

  return(list(
    results_all = all_results,
    summary_total = summary_total,
    perm_dec_table =  perm_dec_table,
    accuracy_by_perm = accuracy_by_perm
  ))
}


