
############## Perform DTRW test ##############

Test_DTRW = function(X,alpha=0.05, method="Bonf"){

  X=X-X[1]
  increments <- diff(X)
  ## Test1: Perform ADF test H1: stationary
    adf_test <- tseries::adf.test(X)
    p_value1 = adf_test$p.value ## we want H0: series is not stationary, so fail to reject H0, so p_value>alpha

  ## Test2:  Perform Ljung-Box test on increments: H0:independently distributed (no autocorrelation)
    ljung_box_test <- Box.test(increments, type = "Ljung-Box")
    p_value2 = ljung_box_test$p.value  ## we want H0: fail to reject so p_value>p

  ##Test3: Wilcoxon signed-rank test for symmetry H0: symmetric
    wilcoxon_test = wilcox.test(increments, mu =0) #median(increments) # DTRW without drift
    p_value3 = wilcoxon_test$p.value

    #p_value = c(p_value1,p_value2,p_value3)

  #decision
    # decision based on Bonferroni
    if (method == "Bonf"){
      dec = ifelse(p_value1> alpha/3 & p_value2 >alpha/3 & p_value3>alpha/3, "DTRW","NO")
      #decision Based on Holm
    } else if (method== "Holm"){
      pp = sort(c(p_value1,p_value2, p_value3))
      Holm = pp > c( 1-(1-alpha)^(1/3),  1-(1-alpha)^(1/2),  alpha)
      dec = ifelse( sum(Holm) == 3, "DTRW","NO") ## Holm is true true true
      #decision Based on Chisquare
    } else if (method == "Chisq"){
      dec = ifelse(-2*(log(p_value1)+ log(p_value2) + log(p_value3) ) >= qchisq(alpha, 6), "DTRW", "NO") #high typr I error
    } else{
      stop("method must be in Bonf, Holm, Chisq")
    }

  list("p_valueStationary" = p_value1,"p_valueIndep" = p_value2, "p_valueSymm" = p_value3, "dec"=dec)
}


Test_DTRW_multiple = function(X,alpha=0.05, method = "Bonf"){

  X=X-X[1]
  ## Calculate the increments / errors e = Xi-X[i-1]
  increments <- diff(X)

  ##  Test 1: identical distribution (stationarity of increments)
      # mid <- floor(length(X)/2)
      # p_value1 <- ks.test(increments[1:mid], increments[(mid+1):length(X)])$p.value  ## we want H0: stationary
      adf_test <- tseries::adf.test(X)
  p_value1 = adf_test$p.value ## we want H0: series is not stationary, so fail to reject H0, so p_value>alpha

  ## Test 2: Independence
      # p_runs <- tseries::runs.test(as.factor(increments>0))$p.value  #small p → evidence against i.i.d. increments (serial dependence).
      ## Perform Ljung-Box test on increments: H0:independently distributed (no autocorrelation)
      ljung_box_test <- Box.test(increments, type = "Ljung-Box")
  p_value2 = ljung_box_test$p.value  ## we want H0: fail to reject so p_value>p

  ##Test3: symmetry H0: symmetric
      # sign test (exact, robust) sign test only checks median=0.
      # p_sign <- binom.test(sum(increments>0), sum(increments!=0), p = 0.5)$p.value
      # Wilcoxon signed-rank (two-sided) # Wilcoxon assumes continuity (no ties) and symmetric null;
      wilcoxon_test <- wilcox.test(increments, mu = median(increments), alternative = "two.sided", exact = FALSE)
  p_value3 = wilcoxon_test$p.value   ## we want H0: symmetric


    # decision based on Bonferroni
  if (method == "Bonf"){
        dec = ifelse(p_value1> alpha/3 & p_value2 >alpha/3 & p_value3>alpha/3, "DTRW","NO")
    #decision Based on Holm
  } else if (method== "Holm"){
    pp = sort(c(p_value1,p_value2, p_value3))
    Holm = pp > c( 1-(1-alpha)^(1/3),  1-(1-alpha)^(1/2),  alpha)
    dec = ifelse( sum(Holm) == 3, "DTRW","NO") ## Holm is true true true
    #decision Based on Chisquare
  } else if (method == "Chisq"){
    dec = ifelse(-2*(log(p_value1)+ log(p_value2) + log(p_value3) ) >= qchisq(alpha, 6), "DTRW", "NO") #high typr I error
  } else{
    stop("method must be in Bonf, Holm, Chisq")
  }
  list("p_valueStationary" = p_value1,"p_valueIndep" = p_value2,
       "p_valueSymm" = p_value3,"dec"=dec)
}


############## based on exact distribution of number of records ##############
Quantile_DTRW=function(alpha=0.05, T){
  Prob = 0
  for(i in 1:T){
    Prob[i]= NT_DTRW(m=i, T=T)
  }
  CDF = cumsum(Prob)  ## cumulative distribution
  #plot(x=1:T, CDF)

  #return( c( which.min(abs(CDF-(alpha/2))) ,which.min(abs(CDF-(1-alpha/2)))))
  return( c( which.min(abs(CDF-(alpha))) ,which.min(abs(CDF-(1-alpha)))))
}

Test_DTRW_NT =function(X,alpha=0.05){
  z_theo = Quantile_DTRW(alpha=alpha, T=length(X))

  obs = rec_counts(X)
  #decision
  dec=ifelse(obs<=z_theo[2] & obs >= z_theo[1] , "DTRW","NO")

  return(list("stat"=obs,"stat_theo"=z_theo,"dec"=dec))
}

######################

Test_DTRW_ENT <- function(X, alpha= 0.05) {
  # X: numeric time X (one realization)

  T <- length(X)

  # Step 1: identify records
  records <- record_times(X)  # times of records
  N_T <- length(records)                      # number of records

  # Step 2: Fit log-log regression across subsamples
  # Split into blocks to see growth of R_t vs t
  block_sizes <- floor(seq(T/4, T, length.out = 5))  # subsample lengths
  R_block <- sapply(block_sizes, function(t) {
    length(record_values(X[1:t]))  ## number of records by blocks
  })

  df <- data.frame(logT = log(block_sizes), logR = log(R_block))

  fit <- lm(logR ~ logT, data = df)

  beta_hat <- coef(fit)[2]
  se_beta  <- summary(fit)$coefficients[2,2]

  # Wald test: H0: beta = 0.5
  z <- (beta_hat - 0.5) / se_beta
  pval <- 2 * (1 - pnorm(abs(z)))

  return(list(
    regression = summary(fit),
    beta_hat = beta_hat,
    se = se_beta,
    z = z,
    pval = pval,
    dec = ifelse(pval < alpha, "NO", "DTRW")
  ))
}

#######################

## FOUR TESTS FOR THE RANDOM WALK HYPOTHESIS - Handa Test 4 - Detection rte of 10% at 5%
Test_DTRW_LSE<- function(X, alpha= 0.05) {
  T <- length(X)

  y1=X[-1]
  y2=X[-length(X)]

  num = sum(y1*y2)

  denom = sum(X^2)

  alpha_hat = num/denom

  T_p = T * (alpha_hat -1 )/sqrt(2)

  #p_val = pnorm(T_p)

  return(list(
    z = T_p,
    dec = ifelse(T_p <= -5.79 | T_p > 0.922, "NO", "DTRW")
  ))
}

## FOUR TESTS FOR THE RANDOM WALK HYPOTHESIS - Handa Test 8 - Detection rte of 10% at 5%
Test_DTRW_tratio <- function(X, alpha= 0.05) {
  T <- length(X)

  y1=X[-1]  ##yt
  y2=X[-length(X)]  ##yt-1

  alpha_hat = sum(y1*y2)/sum(X^2)

  denom = y1-alpha_hat*y2

  T_p = T * sqrt(sum(X^2)) * (alpha_hat -1 )/sum(denom^2)

  #p_val = pnorm(T_p)

  return(list(
    z = T_p,
    "lb" = -5.79,
    "ub"= 0.922,
    dec = ifelse(T_p <= -5.79 | T_p > 0.922, "NO", "DTRW") ## Null is DTRW
  ))
}
