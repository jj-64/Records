Max_Likelihood_Xt = function(X,Likelihood_function,
                             lower_bounds=c(1.001, 0.001, 0.001),
                             upper_bounds = c(5,20,20),  # Example: Bound params[1] by min(R/L)
                             start_values = c(1.001, 0.001, 0.001)
                             ){
  if(length(rec_values(X))<=1) {return("Impossible to solve since we have only one trivial record")}  ## only one record, ignore

  success <- FALSE

  ##Reverse to minimize
  Likelihood = function(params){ x=X; T=T; return(-Likelihood_function(T=T, x=X, params)) }

 MLE_C= nlminb (start_values,Likelihood,  lower=lower_bounds, upper = upper_bounds)

  # Check if the result is successful (i.e., no error occurred)
  if (inherits(MLE_C$par, "try-error")){
    MLE_C <- maxLik(  Likelihood_function, start=start_values)
  } else if(inherits(MLE_C$par, "try-error")){
    MLE_C=   optim(par = start_values,           # Starting values for params
                   fn = Likelihood,        # The likelihood function
                   method = "L-BFGS-B",
                   lower = lower_bounds,         # Lower bounds
                   upper = upper_bounds)
  }

  return(MLE_C)
}

####### to be fixed
Max_Likelihood_Rn = function(X,Likelihood_function,T,
                             lower_bounds=c(1.001, 0.001, 0.001),
                             upper_bounds = c(5,20,20),  # Example: Bound params[1] by min(R/L)
                             start_values = c(1.001, 0.001, 0.001)
                                      ){
  if(length(rec_values(X))<=1) {return("Impossible to solve since we have only one trivial record")}  ## only one record, ignore

  success <- FALSE

  ##Reverse to minimize
  Likelihood = function(params){ y=y; T=T; return(-Likelihood_function(R=rec_values(X), L=rec_times(X),T=T, X=X, params)) }

  MLE_C <- optim(par = start_values,           # Starting values for params
                 fn = Likelihood,        # The likelihood function
                 method = "L-BFGS-B",
                 lower = lower_bounds,         # Lower bounds
                 upper = upper_bounds)

  result=MLE_C$par
  # Check if the result is successful (i.e., no error occurred)
  if (inherits(result, "try-error")){
    MLE_C <- maxLik(  Likelihood_function, start=start_values)
  }

  return(MLE_C)
}
#### Example

# Max_Likelihood_Xt(X=y, Likelihood_function = Likelihood_Xt_Frechet_Yang)
# Max_Likelihood_Xt(X=y, Likelihood_function = Likelihood_Xt_Frechet_LDM,
#                                             lower_bounds=c(0.001, 0.001, 0.001),
#                                             upper_bounds = c(5,20,20),  # Example: Bound params[1] by min(R/L)
#                                             start_values = c(0.001, 0.001, 0.001))
#
# Max_Likelihood_Xt(X=y, Likelihood_function = Likelihood_Xt_Frechet_iid,
#                   lower_bounds=c(0.1, 0.1),
#                    upper_bounds = c(20,20),  # Example: Bound params[1] by min(R/L)
#                   start_values = c(0.1, 0.1))
