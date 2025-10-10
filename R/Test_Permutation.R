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
#' @return A data frame with 24 rows:
#' \describe{
#'   \item{Permutation}{Order of tests applied (e.g., "DLCY")}
#'   \item{Accepted_Test}{The first test that did not reject ("NO")}
#'   \item{Decision}{Either the accepting test name or "All_Rejected"}
#' }
#' A Summary of accepted models across permutations
#' @export
#' @examples
#' \dontrun{
#' result <- Test_Permutation(X = rnorm(50))$decision
#' print(result)
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


  if (print) {
    message("Summary of accepted models across permutations:")
    print(summary)
        }
  return(list(decision = results, Summary = summary))
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
