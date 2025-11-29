
##################Variance estimators - Xt - Frechet - YNM ##########################
#' Variance of LogL Gamma Estimator for full Frechet YANG-Nevzorov Process
#'
#' Computes the variance of the gamma estimator using the Fisher information matrix.
#' The estimator is computed by maximizing the log-likelihood function using the full dataset information.
#'
#' @details
#' The Fisher information is given by:
#' \deqn{ I(\gamma) = -\frac{T (T + 1)}{2 \gamma^2} - \sum_{t=1}^{T} t (t-1) \gamma^{t-2} (A y)^{-a} }
#'
#' The variance is then computed as:
#' \deqn{ V(\gamma) = -\frac{1}{I(\gamma)} }
#'
#' @param x Numeric vector of observations.
#' @param params Numeric vector of three values:  estimated gamma parameter \eqn{ >1}; estimated 1/scale parameter \eqn{ >0}; estimated shape parameter \eqn{ >0}.
#' @return Numeric, variance of gamma estimator.
#' @examples
#' Xt = YNM_series_Frechet(T=25,gamma=1.25, shape=2, scale=1)
#'> Xt
#'[1]  1.150213  1.967684  2.815744  5.612472  1.894868  1.545855  2.086920
#'[8]  2.786223  5.900306  3.548544  6.620712  2.774885  6.421587  7.112956
#'[15]  2.554307  3.329944 16.307641 13.773209 19.188507  6.757204  7.669945
#'[22] 10.595793 11.145220 35.001555 11.652975
#'  # Set parameter bounds: Lower and Upper bounds for gAa[1], gAa[2], and gAa[3]
#' lower_bounds <- c(1.001, 0.001, 0.001)
#' upper_bounds <- c(5,20,20)  # Example: Bound gAa[1] by min(R/L)
#' start_values <- c(1.001, 0.001, 0.001)
#'
#' Log_L= function(gAa){ x=Xt; T=T; return(-Likelihood_Xt_Frechet_YNM(T=T, x=Xt, gAa)) }
#'
#' ## Run optimization using optim with L-BFGS-B method
#' MLE_C <- optim(par = start_values,           # Starting values for gAa
#'               fn = Log_L,        # The likelihood function
#'               method = "L-BFGS-B",
#'               lower = lower_bounds,         # Lower bounds
#'               upper = upper_bounds)$par
#'MLE_C
#'1.2312049 0.7976267 2.2797543  ## Gamma, A, a
#'V= V_Xt_Frechet_YNM_gamma(x=Xt, g=MLE_C[1], A= MLE_C[2], a= MLE_C[3])
#'V
#'[1] 0.0002850398
#' #standard deviation
#' sqrt(V)
#' [1] 0.01688312]
#' #Confidence interval
#' bounds(value = MLE_C[1], z = 1.96 , variance = V)
#' [1] 1.198114 1.264296
V_Xt_Frechet_YNM_gamma <- function(x, params) {
  g = params[1]
  A = params[2]
  a = params[3]
  T <- length(x)
  sum_val <- sum((1:T) * (0:(T-1)) * (g^((1:T)-2)) * (A * x)^(-a))
  fisher <- (-T * (T + 1) / (2 * g^2)) - sum_val
  return(-1 / fisher)
}

#' Variance of 1/Scale LogL Estimator for full Frechet YANG-Nevzorov Process
#'
#' Computes the variance of the A; 1/scale estimator using the Fisher information matrix.
#' The estimator is computed by maximizing the log-likelihood function using the full dataset information.
#'
#' @details
#' The Fisher information is given by:
#' \deqn{ I(A) = \frac{a T}{A^2} - a (a+1) A^{-a-2} \sum_{t=1}^{T} \gamma^t y^{-a} }
#'
#' @inheritParams V_Xt_Frechet_YNM_gamma
#' @return Numeric, variance of gamma estimator.
#' @examples
#' Xt = YNM_series_Frechet(T=25,gamma=1.25, shape=2, scale=1)
#'> Xt
#'[1]  1.150213  1.967684  2.815744  5.612472  1.894868  1.545855  2.086920
#'[8]  2.786223  5.900306  3.548544  6.620712  2.774885  6.421587  7.112956
#'[15]  2.554307  3.329944 16.307641 13.773209 19.188507  6.757204  7.669945
#'[22] 10.595793 11.145220 35.001555 11.652975
#'  # Set parameter bounds: Lower and Upper bounds for gAa[1], gAa[2], and gAa[3]
#' lower_bounds <- c(1.001, 0.001, 0.001)
#' upper_bounds <- c(5,20,20)  # Example: Bound gAa[1] by min(R/L)
#' start_values <- c(1.001, 0.001, 0.001)
#'
#' Log_L= function(gAa){ x=Xt; T=T; return(-Likelihood_Xt_Frechet_YNM(T=T, x=Xt, gAa)) }
#'
#' ## Run optimization using optim with L-BFGS-B method
#' MLE_C <- optim(par = start_values,           # Starting values for gAa
#'               fn = Log_L,        # The likelihood function
#'               method = "L-BFGS-B",
#'               lower = lower_bounds,         # Lower bounds
#'               upper = upper_bounds)$par
#'MLE_C
#'1.2312049 0.7976267 2.2797543  ## Gamma, A, a
#'V= V_Xt_Frechet_YNM_scale(x=Xt, g=MLE_C[1], A= MLE_C[2], a= MLE_C[3])
#'V
#'[1] 0.004896407
#' #standard deviation
#' sqrt(V)
#' [1] 0.06997433
#' #Confidence interval
#' bounds(value = MLE_C[2], z = 1.96 , variance = V)
#' [1] 0.6604770 0.9347764
V_Xt_Frechet_YNM_scale = function(x,params){
  g = params[1]
  A = params[2]
  a = params[3]

  T=length(x)
  sum_val = sum( (g^(1:T))*(x^-a)  )
  fisher =  (a*T/(A^2)) -  a*(a+1) * A^(-a-2)*  sum_val
  return(-1/fisher)
}


#' Variance of shape LogL Estimator for full Frechet YANG-Nevzorov Process
#'
#' Computes the variance of the a:shape estimator using the Fisher information matrix.
#' The estimator is computed by maximizing the log-likelihood function using the full dataset information.
#'
#' @details
#' The Fisher information is given by:
#' \deqn{ I(a) = -\frac{T}{a^2} - \sum_{t=1}^{T} \gamma^t (\log(A y))^2 (A y)^{-a} }
#'
#' @inheritParams V_Xt_Frechet_YNM_gamma
#' @return Numeric, variance of gamma estimator.
#' @examples
#' Xt = YNM_series_Frechet(T=25,gamma=1.25, shape=2, scale=1)
#'> Xt
#'[1]  1.150213  1.967684  2.815744  5.612472  1.894868  1.545855  2.086920
#'[8]  2.786223  5.900306  3.548544  6.620712  2.774885  6.421587  7.112956
#'[15]  2.554307  3.329944 16.307641 13.773209 19.188507  6.757204  7.669945
#'[22] 10.595793 11.145220 35.001555 11.652975
#'  # Set parameter bounds: Lower and Upper bounds for gAa[1], gAa[2], and gAa[3]
#' lower_bounds <- c(1.001, 0.001, 0.001)
#' upper_bounds <- c(5,20,20)  # Example: Bound gAa[1] by min(R/L)
#' start_values <- c(1.001, 0.001, 0.001)
#'
#' Log_L= function(gAa){ x=Xt; T=T; return(-Likelihood_Xt_Frechet_YNM(T=T, x=Xt, gAa)) }
#'
#' ## Run optimization using optim with L-BFGS-B method
#' MLE_C <- optim(par = start_values,           # Starting values for gAa
#'               fn = Log_L,        # The likelihood function
#'               method = "L-BFGS-B",
#'               lower = lower_bounds,         # Lower bounds
#'               upper = upper_bounds)$par
#'MLE_C
#'1.2312049 0.7976267 2.2797543  ## Gamma, A, a
#'V= V_Xt_Frechet_YNM_shape(x=Xt, g=MLE_C[1], A= MLE_C[2], a= MLE_C[3])
#'V
#'[1] 0.02437784
#' #standard deviation
#' sqrt(V)
#' [1]  0.156134
#' #Confidence interval
#' bounds(value = MLE_C[3], z = 1.96 , variance = V)
#' [1] 1.973732 2.585777
V_Xt_Frechet_YNM_shape = function(x,params){
  g = params[1]
  A = params[2]
  a = params[3]
  T=length(x)
  sum_val = sum(  (g^(1:T))  *  (log(A *x))^2    * (A*x)^(-a)  )
  fisher =  (-T / (a^2)) - sum_val
  return(-1/fisher)
}


##################Variance estimators - Xt - Frechet - LDM ##########################

#' Variance of Theta LogL Estimator for Frechet LDM using full dataset
#'
#' @details
#' Estimators are obtained by maximizing the Log-Likelihood based on full dataset \eqn{X_t} \code{\link{Likelihood_Xt_Frechet_LDM}}
#' The Fisher information is given by:
#' \deqn{ I(\theta) = A^{-a} a (-a-1) \sum_{t=1}^{T} t^2 x^{-a-2} + (a+1) \sum_{t=1}^{T} \frac{t^2}{x^2} }
#'
#' @param x Numeric vector of observations.
#' @param params Numeric vector of three values:  drift parameter \eqn{\theta > 0};  scale parameter \eqn{A > 0};  shape parameter \eqn{a > 0}.
#' @return Numeric, variance of theta estimator.
#' @examples
#' Xt = LDM_series(T=25,dist = "frechet", theta=0.2, shape=2, scale=1)
#'Xt
#'[1] 1.252804 1.465630 2.622688 2.061684 2.723450 4.932948 2.717134 2.368716
#'[9] 2.487995 2.687706 3.656084 4.580848 4.215754 4.284083 4.110642 3.851272
#'[17] 4.031203 5.637333 5.372408 4.897131 5.071756 5.137329 5.698935 5.624777
#'[25] 9.041731
#'  # Set parameter bounds: Lower and Upper bounds for gAa[1], gAa[2], and gAa[3]
#' lower_bounds <- c(0.001, 0.001, 0.001)
#' upper_bounds <- c(5,20,20)  # Example: Bound gAa[1] by min(R/L)
#' start_values <- c(0.001, 0.001, 0.001)
#'
#' Log_L= function(gAa){ x=Xt; T=T; return(-Likelihood_Xt_Frechet_LDM(T=T, x=Xt, gAa)) }
#'
#' ## Run optimization using optim with L-BFGS-B method
#' MLE_C <- optim(par = start_values,           # Starting values for gAa
#'               fn = Log_L,        # The likelihood function
#'               method = "L-BFGS-B",
#'               lower = lower_bounds,         # Lower bounds
#'               upper = upper_bounds)$par
#' MLE_C
#' [1] 0.1915248 0.0010000 2.8231442  #theta, A, a
#' MLE_C[2] = eq_scale(a=MLE_C[3],Xt,theta=MLE_C[1])  ## fix for A estimation
#' MLE_C
#' 0.1915248 0.9054421 2.8231442
#' ## Or alternatively
#' Likelihood_Xt_Frechet_LDM_OptSolve(x=Xt)
#'       theta           A           a        logL
#'    0.1915428   0.9056486   2.8225928 -21.5907410
#'  ##Or more accurately using first-order derivatives
#'  Likelihood_Xt_Frechet_LDM_OptSolve(x=Xt)
#'  theta           A           a        logL
#'  0.1915428   0.9056486   2.8225928 -21.5907410
#'V= V_Xt_Frechet_LDM_theta(x=Xt, g=MLE_C[1], A= MLE_C[2], a= MLE_C[3])
#'V
#'[1] 0.00001897734
#' #Confidence interval
#' bounds(value = MLE_C[1], z = 1.96 , variance = V)
#' [1]  0.1829864 0.2000631
V_Xt_Frechet_LDM_theta =  function(x,params) { ## x vector and theta, A and alpha
  g = params[1]
  A = params[2]
  a = params[3]

  T= length(x)
  y = x-g*(1:T) ## frechet iid
  s1 = A^(-a)*a*(-a-1) * sum((1:T)^2 * y^(-a-2))
  s2 = (a+1)*sum(((1:T)/y)^2)
  fisher= s1+s2
  return(-1/fisher)
}

#'  Variance of A:1/Scale LogL Estimator for Frechet LDM using full dataset
#'
#' @details
#' Estimators are obtained by maximizing the Log-Likelihood based on full dataset \eqn{X_t} \code{\link{Likelihood_Xt_Frechet_LDM}}
#' The Fisher information is given by:
#' \deqn{ I(A) = - (a+1) a A^{-a-2} \sum_{t=1}^{T} x^{-a} + \frac{a T}{A^2} }
#'
#' @inheritParams V_Xt_Frechet_LDM_theta
#' @return Numeric, variance of A:1/Scale estimator.
#' @examples
#' Xt = LDM_series(T=25,dist = "frechet", theta=0.2, shape=2, scale=1)
#'Xt
#'[1] 1.252804 1.465630 2.622688 2.061684 2.723450 4.932948 2.717134 2.368716
#'[9] 2.487995 2.687706 3.656084 4.580848 4.215754 4.284083 4.110642 3.851272
#'[17] 4.031203 5.637333 5.372408 4.897131 5.071756 5.137329 5.698935 5.624777
#'[25] 9.041731
#'  # Set parameter bounds: Lower and Upper bounds for gAa[1], gAa[2], and gAa[3]
#' lower_bounds <- c(0.001, 0.001, 0.001)
#' upper_bounds <- c(5,20,20)  # Example: Bound gAa[1] by min(R/L)
#' start_values <- c(0.001, 0.001, 0.001)
#'
#' Log_L= function(gAa){ x=Xt; T=T; return(-Likelihood_Xt_Frechet_LDM(T=T, x=Xt, gAa)) }
#'
#' ## Run optimization using optim with L-BFGS-B method
#' MLE_C <- optim(par = start_values,           # Starting values for gAa
#'               fn = Log_L,        # The likelihood function
#'               method = "L-BFGS-B",
#'               lower = lower_bounds,         # Lower bounds
#'               upper = upper_bounds)$par
#' MLE_C
#' [1] 0.1915248 0.0010000 2.8231442  #theta, A, a
#' ### Estimate A correctly
#' MLE_C[2] = eq_scale(a=MLE_C[3],Xt,theta=MLE_C[1])
#' MLE_C
#' 0.1915248 0.9054421 2.8231442
#'
#' ## Alternatively, we can obtain accurate estimates, using first-order derivatives
#' Likelihood_Xt_Frechet_LDM_OptSolve(x=Xt)
#'       theta           A           a        logL
#'    0.1915428   0.9056486   2.8225928 -21.5907410
#'
#'V= V_Xt_Frechet_LDM_scale(x=Xt, g=MLE_C[1], A= MLE_C[2], a= MLE_C[3])
#'V
#'[1] 0.004114482
#'
#' #Confidence interval
#' bounds(value = MLE_C[2], z = 1.96 , variance = V)
#' [1]   0.7797194 1.0311647
V_Xt_Frechet_LDM_scale = function(x,params){
  g = params[1]
  A = params[2]
  a = params[3]

  T=length(x)
  y = x-g*(1:T)
  s1 = -(a+1)*a*A^(-a-2)*sum(y^-a)
  s2=a*T/A^2
  fisher= s1+s2
  return(-1/fisher)
}

#' Variance of Shape LogL Parameter a for Frechet LDM using full dataset
#'
#' @details
#' Estimators are obtained by maximizing the Log-Likelihood based on full dataset \eqn{X_t} \code{\link{Likelihood_Xt_Frechet_LDM}}
#' The Fisher information is given by:
#' \deqn{ I(a) = - A^{-a} \sum_{t=1}^{T} x^{-a} (\log(A x))^2 - \frac{T}{a^2} }
#'
#' @inheritParams V_Xt_Frechet_LDM_theta
#' @return Numeric, variance of a shape estimator.
#' @examples
#' Xt = LDM_series(T=25,dist = "frechet", theta=0.2, shape=2, scale=1)
#'Xt
#'[1] 1.252804 1.465630 2.622688 2.061684 2.723450 4.932948 2.717134 2.368716
#'[9] 2.487995 2.687706 3.656084 4.580848 4.215754 4.284083 4.110642 3.851272
#'[17] 4.031203 5.637333 5.372408 4.897131 5.071756 5.137329 5.698935 5.624777
#'[25] 9.041731
#'  # Set parameter bounds: Lower and Upper bounds for gAa[1], gAa[2], and gAa[3]
#' lower_bounds <- c(0.001, 0.001, 0.001)
#' upper_bounds <- c(5,20,20)  # Example: Bound gAa[1] by min(R/L)
#' start_values <- c(0.001, 0.001, 0.001)
#'
#' Log_L= function(gAa){ x=Xt; T=T; return(-Likelihood_Xt_Frechet_LDM(T=T, x=Xt, gAa)) }
#'
#' ## Run optimization using optim with L-BFGS-B method
#' MLE_C <- optim(par = start_values,           # Starting values for gAa
#'               fn = Log_L,        # The likelihood function
#'               method = "L-BFGS-B",
#'               lower = lower_bounds,         # Lower bounds
#'               upper = upper_bounds)$par
#' MLE_C
#' [1] 0.1915248 0.0010000 2.8231442  #theta, A, a
#' MLE_C[2] = eq_scale(a=MLE_C[3],Xt,theta=MLE_C[1])  ## fix for A estimation
#' MLE_C
#' 0.1915248 0.9054421 2.8231442
#' ## Or alternatively
#' Likelihood_Xt_Frechet_LDM_OptSolve(x=Xt)
#'       theta           A           a        logL
#'    0.1915428   0.9056486   2.8225928 -21.5907410
#'  ##Or more accurately using first-order derivatives
#'  Likelihood_Xt_Frechet_LDM_OptSolve(x=Xt)
#'  theta           A           a        logL
#'  0.1915428   0.9056486   2.8225928 -21.5907410
#' V= V_Xt_Frechet_LDM_shape(x=Xt, g=MLE_C[1], A= MLE_C[2], a= MLE_C[3])
#' V
#' [1] 0.1863364
#' #Confidence interval
#' bounds(value = MLE_C[3], z = 1.96 , variance = V)
#' [1]1.977077 3.669211
V_Xt_Frechet_LDM_shape = function(x,params){
  g = params[1]
  A = params[2]
  a = params[3]

  T=length(x)
  y = x-g*(1:T)
  p1 = y^(-a)*(log(A*y))^2
  s1 = -A^(-a)*sum(p1)
  s2 = -T/(a^2)
  fisher =  s1+s2
  return(-1/fisher)
}

##################Variance estimators - Xt - Frechet - iid ##########################

## Variance of 1/scale A
#' Variance of 1/Scale Parameter A for IID Frechet Process full dataset
#'
#' @details
#' Estimators are obtained by maximizing the Log-Likelihood based on full dataset \eqn{X_t} \code{\link{Likelihood_Xt_Frechet_iid}}
#' The Fisher information is given by:
#' \deqn{ I(A) = \frac{aT}{A^2} - \frac{(a^2 + a)}{A^2} \sum_{t=1}^{T} (x_t A)^{-a} }
#'
#' @param x Numeric vector of observations \eqn{x>0}.
#' @param params Numeric vector of two values: 1/scale parameter of Frechet distribution \eqn{A >0}; shape parameter of Frechet distribution \eqn{a>0}.
#' @return Numeric, variance of A estimator.
#' @examples
#' Xt = VGAM::rfrechet(n=25,shape=2, scale=1)
#' Xt
#' [1] 1.4654507 1.0621522 1.4127990 0.4982294 6.4634729 0.4731131 0.9917193
#' [8] 2.5978689 1.1208717 0.7365184 0.8867496 3.6458321 0.8661168 0.7883508
#' [15] 2.7997818 3.8250220 2.5340879 2.3363193 1.4448465 1.8548649 3.0123339
#' [22] 0.7486512 1.2412334 1.1041981 1.4358524
#'  # Set parameter bounds: Lower and Upper bounds for gAa[1], gAa[2], and gAa[3]
#' lower_bounds <- c( 0.1, 0.1)
#' upper_bounds <- c(20,20)
#' start_values <- c( 0.1, 0.1)
#'
#' Log_L= function(gAa){ x=Xt; T=T; return(-Likelihood_Xt_Frechet_iid(T=T, x=Xt, gAa)) }
#'
#' ## Run optimization using optim with L-BFGS-B method
#' MLE_C <- optim(par = start_values,           # Starting values for gAa
#'               fn = Log_L,        # The likelihood function
#'               method = "L-BFGS-B",
#'               lower = lower_bounds,         # Lower bounds
#'               upper = upper_bounds)$par
#' MLE_C
#' [1]  0.9397696 1.7604020   # A, a
#' V= V_Xt_Frechet_iid_scale(x=Xt, A= MLE_C[1], a= MLE_C[2])
#' V
#' [1] 0.01139935
#' #Confidence interval
#' bounds(value = MLE_C[1], z = 1.96 , variance = V)
#' [1] 0.7305049 1.1490343
V_Xt_Frechet_iid_scale = function(x,params){
  A = params[1]
  a = params[2]

  T=length(x)

  s1 = (a*T)/A^2

  s2 = -((a^2+a)/A^2) * sum((x*A)^(-a))

  var_A <- -1 / (s1 + s2)

  return(var_A)
}


## Variance of shape a
#' Variance of Shape Parameter a for IID Frechet Process using full dataset
#'
#' @details
#' Estimators are obtained by maximizing the Log-Likelihood based on full dataset \eqn{X_t} \code{\link{Likelihood_Xt_Frechet_iid}}
#' The Fisher information is given by:
#' \deqn{ I(a) = -\frac{T}{a^2} - \sum_{t=1}^{T} \frac{[\log(A x_t)]^2}{(A x_t)^a} }
#'
#' @inheritParams V_Xt_Frechet_iid_scale
#' @return Numeric, variance of a estimator.
#' @examples
#' Xt = VGAM::rfrechet(n=25,shape=2, scale=1)
#' Xt
#' [1] 1.4654507 1.0621522 1.4127990 0.4982294 6.4634729 0.4731131 0.9917193
#' [8] 2.5978689 1.1208717 0.7365184 0.8867496 3.6458321 0.8661168 0.7883508
#' [15] 2.7997818 3.8250220 2.5340879 2.3363193 1.4448465 1.8548649 3.0123339
#' [22] 0.7486512 1.2412334 1.1041981 1.4358524
#'  # Set parameter bounds: Lower and Upper bounds for gAa[1], gAa[2], and gAa[3]
#' lower_bounds <- c( 0.1, 0.1)
#' upper_bounds <- c(20,20)
#' start_values <- c( 0.1, 0.1)
#'
#' Log_L= function(gAa){ x=Xt; T=T; return(-Likelihood_Xt_Frechet_iid(T=T, x=Xt, gAa)) }
#'
#' ## Run optimization using optim with L-BFGS-B method
#' MLE_C <- optim(par = start_values,           # Starting values for gAa
#'               fn = Log_L,        # The likelihood function
#'               method = "L-BFGS-B",
#'               lower = lower_bounds,         # Lower bounds
#'               upper = upper_bounds)$par
#' MLE_C
#' [1]  0.9397696 1.7604020   # A, a
#' V= V_Xt_Frechet_iid_shape(x=Xt, A= MLE_C[1], a= MLE_C[2])
#' V
#' [1] 0.06489838
#' #Confidence interval
#' bounds(value = MLE_C[2], z = 1.96 , variance = V)
#' [1]  1.261089 2.259715
V_Xt_Frechet_iid_shape = function(x,params){
  A = params[1]
  a = params[2]

  T=length(x)

  s1=-T/(a^2)

  s2 = - sum( log(A*x)^2/(A*x)^a )

  var_a = -1 / (s1+s2)

  return(var_a)
}
