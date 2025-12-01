logLik_records <- function(model, obs_type, dist, data, params) {
  # lookup
  f <- loglik_registry[[model]][[obs_type]][[dist]]

  if (is.null(f))
    stop("Likelihood expression not registered for this (model, obs_type, dist).")

  f(data, params)
}

## Helper ---------------------------------
# empty container
loglik_registry <- new.env(parent = emptyenv())

## Function to register likelihood expressions later
register_loglik <- function(model, obs_type, dist, fun) {
  if (!exists(model, envir = loglik_registry))
    loglik_registry[[model]] <- list()

  if (!obs_type %in% names(loglik_registry[[model]]))
    loglik_registry[[model]][[obs_type]] <- list()

  loglik_registry[[model]][[obs_type]][[dist]] <- fun
}

## Classical, Xt -------------------------

  ##Normal
register_loglik( model = "iid", obs_type = "all", dist = "norm",
  fun = function(data, params) {
    if(!is.numeric(data)) stop("data should be a numerical vector")
    if( all(c("mean", "sd") %in% names(params)) == FALSE ) stop("parameters mean, sd should be present in a list.")
    mean <- as.numeric(params['mean'])
    sd <- as.numeric(params["sd"])
    sum(dnorm(data, mean =mean , sd = sd, log = TRUE))
  }
)
  ## Frechet
register_loglik( "iid", "all", "frechet",
  fun = function(data, params) {
    if(!is.numeric(data)) stop("data should be a numerical vector")
    if( all(c("shape", "scale") %in% names(params)) == FALSE ) stop("parameters shape, scale should be present in a list.")

    sum(VGAM::dfrechet(x = data, location = 0, shape = as.numeric(params["shape"]), scale = params["scale"] ,log = TRUE))
  }
)

  ## Gumbel
register_loglik( "iid", "all", "gumbel",
                 fun = function(data, params) {
                   if(!is.numeric(data)) stop("data should be a numerical vector")
                   if( all(c("loc", "scale") %in% names(params)) == FALSE ) stop("parameters loc, scale should be present in a list.")

                   sum(VGAM::dgumbel(x = data, location = as.numeric(params["loc"]), scale = as.numeric(params["scale"]) ,log = TRUE))
                 }
)
  ## Weibull
register_loglik( "iid", "all", "weibull",
                 fun = function(data, params) {
                   if(!is.numeric(data)) stop("data should be a numerical vector")
                   if( all(c("shape", "scale") %in% names(params)) == FALSE ) stop("parameters shape, scale should be present in a list.")

                   sum(log(dweibull(x = data, shape = as.numeric(params["shape"]), scale = as.numeric(params["scale"]) )))
                 }
)
## Classical, Rn ------------------------

register_loglik("iid", "records", "norm",
  fun = function(data, params) {
    if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

    Rn <- data$rec_values
    Ln <- data$rec_times
    n  <- data$time[1]

    if( all(c("mean", "sd") %in% names(params)) == FALSE ) stop("parameters mean, sd should be present in a list.")

    mean <- as.numeric(params['mean'])
    sd <- as.numeric(params["sd"])
    m = length(Rn) #m= rec_counts(y)

    ## sum of log(f_rn)
    s1 <- sum(log(dnorm(Rn[-1], mean = mean, sd = sd)))
    #Compute s2 using vectorized approach
    intervals <- diff(Ln)-1
    s2b <- intervals * log(pnorm(Rn[-m], mean = mean, sd = sd))
    s2 <- sum(s2b)
    ##Compute s3 only if needed
    s3 <- if (Ln[m] < n) {
      (n - Ln[m]) * log(pnorm(Rn[m], mean = mean, sd = sd))
    } else {
      0
    }

    #Return total log-likelihood
    return(s1 + s2 + s3)
  }
)

register_loglik("iid", "records", "frechet",
                fun = function(data, params) {
                  if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                  Rn <- data$rec_values
                  Ln <- data$rec_times
                  n  <- data$time[1]

                  if( all(c("shape", "scale") %in% names(params)) == FALSE ) stop("parameters shape, scale should be present in a list.")

                  scale = as.numeric(params["scale"])
                  shape = as.numeric(params["shape"])
                  m = length(Rn) #m= rec_counts(y)

                  ## sum of log(f_rn)
                  s1 <- sum(log(VGAM::dfrechet(Rn[-1], shape, scale)))
                  #Compute s2 using vectorized approach
                  intervals <- diff(Ln)-1
                  s2b <- intervals * log(VGAM::pfrechet(Rn[-m], shape, scale))
                  s2 <- sum(s2b)
                  ##Compute s3 only if needed
                  s3 <- if (Ln[m] < n) {
                    (n - Ln[m]) * log(VGAM::pfrechet(Rn[m], shape, scale))
                  } else {
                    0
                  }

                  #Return total log-likelihood
                  return(s1 + s2 + s3)
                }
)

register_loglik("iid", "records", "gumbel",
                fun = function(data, params) {
                  if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                  Rn <- data$rec_values
                  Ln <- data$rec_times
                  n  <- data$time[1]

                  if( all(c("loc", "scale") %in% names(params)) == FALSE ) stop("parameters loc, scale should be present in a list.")

                  scale = as.numeric(params["scale"])
                  loc = as.numeric(params["loc"])
                  m = length(Rn) #m= rec_counts(y)

                  ## sum of log(f_rn)
                  s1 <- sum(log(VGAM::dgumbel(Rn[-1], loc, scale)))
                  #Compute s2 using vectorized approach
                  intervals <- diff(Ln)-1
                  s2b <- intervals * log(VGAM::pgumbel(Rn[-m], loc, scale))
                  s2 <- sum(s2b)
                  ##Compute s3 only if needed
                  s3 <- if (Ln[m] < n) {
                    (n - Ln[m]) * log(VGAM::pfrechetgumbel(Rn[m], loc, scale))
                  } else {
                    0
                  }

                  #Return total log-likelihood
                  return(s1 + s2 + s3)
                }
)

register_loglik("iid", "records", "weibull",
                fun = function(data, params) {
                  if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                  Rn <- data$rec_values
                  Ln <- data$rec_times
                  n  <- data$time[1]

                  if( all(c("shape", "scale") %in% names(params)) == FALSE ) stop("parameters shape, scale should be present in a list.")

                  scale = as.numeric(params["scale"])
                  shape = as.numeric(params["shape"])
                  m = length(Rn) #m= rec_counts(y)

                  ## sum of log(f_rn)
                  s1 <- sum(log(dweibull(Rn[-1], shape, scale)))
                  #Compute s2 using vectorized approach
                  intervals <- diff(Ln)-1
                  s2b <- intervals * log(pweibull(Rn[-m], shape, scale))
                  s2 <- sum(s2b)
                  ##Compute s3 only if needed
                  s3 <- if (Ln[m] < n) {
                    (n - Ln[m]) * log(pweibull(Rn[m], shape, scale))
                  } else {
                    0
                  }

                  #Return total log-likelihood
                  return(s1 + s2 + s3)
                }
)


## DTRW, Xt ---------------------------------

  ## Normal
register_loglik("DTRW", "all", "norm",
  fun = function(data, params) {
    if(!is.numeric(data)) stop("data should be a numerical vector")

    if( all(c("mean", "sd") %in% names(params)) == FALSE ) stop("parameters mean, sd should be present in a list.")

    mean <- as.numeric(params['mean'])
    sd <- as.numeric(params["sd"])

    sum(log(dnorm(x = diff(data), mean =mean , sd = sd)))
  }
)

  ##cauchy
register_loglik("DTRW", "all", "cauchy",
                fun = function(data, params) {
                  if(!is.numeric(data)) stop("data should be a numerical vector")
                  if( all(c("loc", "scale") %in% names(params)) == FALSE ) stop("parameters loc, scale should be present in a list.")

                  sum(log(dnorm(x = diff(data), location =params["loc"] , scale = as.numeric(params["scale"]))))
                }
)

## DTRW, Rn -------------------------------------
  ## Normal
register_loglik( "DTRW", "records", "norm",
  fun = function(data, params) {
    if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

    Rn <- data$rec_values
    Ln <- data$rec_times
    n <- data$time[1]
    if( all(c( "scale") %in% names(params)) == FALSE ) stop("parameters scale should be present in a list.")

    loc = 0
    sd =as.numeric(params["sd"])

    # Example: probability that S_1, S_2, ..., S_5 all < 0
    # sigma_cov = function(k, sigma){  ## k should be less than 1000, sigma is variance
    #   x=matrix(0, nrow=k, ncol=k)
    #   for(i in 1:k){
    #     for(j in 1:k){
    #       x[i,j] = min(i,j)*sigma
    #     }
    #   }
    #   return(x)
    # }

    m= length(Rn)  # Number of observed time points

    if (m < 2) {
      stop("The vector l must contain at least two indices.")
    }

    s1b = log(dnorm(diff(Rn), mean = 0, sd = sqrt(diff(Ln)*sd) )  )

    s1=sum(s1b)

    ## case where NT<T
    s2=0
    if (Ln[m] < n) {
      #s2=log(mvtnorm::pmvnorm(lower=-Inf, upper=rep(0,(n - Ln[m])), sigma = sigma_cov((n - Ln[m]),params))[1])
      s2 = -0.5*log(pi * (n - Ln[m]))
    }

    return(s1+s2)
  }
)

  ## cauchy
register_loglik( "DTRW", "records", "cauchy",
                 fun = function(data, params) {
                   if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                   Rn <- data$rec_values
                   Ln <- data$rec_times
                   n <- data$time[1]

                   if( all(c("scale") %in% names(params)) == FALSE ) stop("parameters  scale should be present in a list.")

                   loc = 0
                   scale = as.numeric(params["scale"])
                   m= length(Rn)  # Number of observed time points

                   if (m < 2) {
                     stop("The vector l must contain at least two indices.")
                   }

                   s1b = log(dcauchy(diff(Rn), location = 0, scale = sqrt(diff(Ln)*scale) )  )
                   s1=sum(s1b)

                   ## case where NT<T
                   s2=0
                   if (Ln[m] < n) { s2 = -0.5*log(pi * (n - Ln[m])) }

                  ## Return
                   return(s1+s2)
                 }
)

## YNM, Xt -----------------------------------------

  ##Frechet
register_loglik( "YNM", "all", "frechet",
                 fun = function(data, params) {
                   x=data
                   if(!is.numeric(x)) stop("data should be a numerical vector")
                   n= length(as.numeric(x))
                   if( all(c("gamma", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters gamma, shape, scale should be present.")

                   gamma <- as.numeric(params["gamma"])
                   shape <- as.numeric(params["shape"])
                   scale <- as.numeric(params["scale"])

                   ## Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or shape <= 0, we return -Inf)
                   if (shape <=0 | scale <= 0) { ##
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## pdf
                   pdf=function(x,par) {VGAM::dfrechet(x, shape=par$shape, scale=par$scale)}
                   cdf=function(x,par) {VGAM::pfrechet(x, shape=par$shape, scale=par$scale)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum(log(gamma^(1:n)))
                   s2 = sum(log(cdf(x, par=list(shape= shape, scale= scale))^(-1+gamma^(1:n))))
                   s3 = sum(log(pdf(x, par=list(shape= shape, scale= scale))))

                   # s1 <- n* log(shape * scale^(-shape) )
                   # s2 <- -(shape + 1) * sum(log(x))
                   # s3 <- (n* (n+ 1) / 2) * log(gamma)
                   # s4 <- -sum( gamma^(1:n) * (scale * x)^(-shape) )
                   return(s1+s2+s3)
                 }
)

  ##Gumbel
register_loglik( "YNM", "all", "gumbel",
                 fun = function(data, params) {
                   x=data
                   if(!is.numeric(x)) stop("data should be a numerical vector")

                   n= length(as.numeric(x))
                   if( all(c("gamma", "loc", "scale") %in% names(params)) == FALSE ) stop("parameters gamma, loc, scale should be present.")
                   gamma <- as.numeric(params["gamma"])
                   loc <- as.numeric(params["loc"])
                   scale <- as.numeric(params["scale"])

                   ## Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or loc <= 0, we return -Inf)
                   if (scale <= 0) { ##
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## pdf
                   pdf=function(x,par) {VGAM::dgumbel(x, loc=par$loc, scale=par$scale)}
                   cdf=function(x,par) {VGAM::pgumbel(x, loc=par$loc, scale=par$scale)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum(log(gamma^(1:n)))
                   s2 = sum(log(cdf(x, par=list(loc= loc, scale= scale))^(-1+gamma^(1:n))))
                   s3 = sum(log(pdf(x, par=list(loc= loc, scale= scale))))

                   return(s1+s2+s3)
                 }
)

  ##norm
register_loglik( "YNM", "all", "norm",
                 fun = function(data, params) {
                   x=data
                   if(!is.numeric(x)) stop("data should be a numerical vector")

                   n= length(as.numeric(x))
                   if( all(c("gamma", "mean", "sd") %in% names(params)) == FALSE ) stop("parameters gamma, mean, sd should be present.")
                   gamma <- as.numeric(params["gamma"])
                   mean <- as.numeric(params["mean"])
                   sd <- as.numeric(params["sd"])

                   ## Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or mean <= 0, we return -Inf)
                   if (sd <= 0) { ##
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## pdf
                   pdf=function(x,par) {dnorm(x, mean=par$mean, sd=par$sd)}
                   cdf=function(x,par) {pnorm(x, mean=par$mean, sd=par$sd)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum(log(gamma^(1:n)))
                   s2 = sum(log(cdf(x, par=list(mean= mean, sd= sd))^(-1+gamma^(1:n))))
                   s3 = sum(log(pdf(x, par=list(mean= mean, sd= sd))))

                   return(s1+s2+s3)
                 }
)

##weibull
register_loglik( "YNM", "all", "weibull",
                 fun = function(data, params) {
                   x=data
                   if(!is.numeric(x)) stop("data should be a numerical vector")

                   n= length(as.numeric(x))
                   if( all(c("gamma", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters gamma, scale, shape should be present.")
                   gamma <- as.numeric(params["gamma"])
                   shape <- as.numeric(params["shape"])
                   scale <- as.numeric(params["scale"])

                   ## Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or shape <= 0, we return -Inf)
                   if (shape <=0 | scale <= 0) { ##
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## pdf
                   pdf=function(x,par) {dweibull(x, shape=par$shape, scale=par$scale)}
                   cdf=function(x,par) {pweibull(x, shape=par$shape, scale=par$scale)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum(log(gamma^(1:n)))
                   s2 = sum(log(cdf(x, par=list(shape= shape, scale= scale))^(-1+gamma^(1:n))))
                   s3 = sum(log(pdf(x, par=list(shape= shape, scale= scale))))

                   return(s1+s2+s3)
                 }
)
## YNM, Rn -----------------------------------------

  ## Frechet
register_loglik( "YNM", "records", "frechet",
  fun = function(data, params) {
    if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

    Rn <- data$rec_values
    Ln <- data$rec_times
    n <- data$time[1]

    if( all(c("gamma", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters gamma, shape, scale should be present.")

    gamma <- as.numeric(params["gamma"])
    shape <- as.numeric(params["shape"])
    scale <- as.numeric(params["scale"])
    m = length(Rn) #m= rec_counts(y)  ## number of records

    # Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or shape <= 0, we return -Inf)
    if (scale <= 0 | shape<=0) {
      return(-Inf)  # Invalid parameters, return a large negative value
    }

    ## pdf
    pdf=function(x,par) {VGAM::dfrechet(x=x, shape=par$shape, scale=par$scale)}
    ## cdf
    cdf=function(x,par) {VGAM::pfrechet(q=x, shape= par$shape, scale=par$scale)}


    # your exact YNM record-pair likelihood
    s1 = log(gamma^(sum(Ln)))
    s2 = sum(log(pdf(Rn, par = list(shape=shape, scale=scale)) * cdf(Rn, par = list(shape=shape, scale=scale))^(-1+gamma^Ln)))
    s3b = 0
    for( i in 1:(m-1)  ){
      if((Ln[i]+1 <= Ln[i+1]) == TRUE){ ## we have non-records in between
        s3b[i]= sum(log(cdf(Rn[i], par = list(shape=shape, scale=scale))^( (gamma^(Ln[i]+1) - gamma^(Ln[i+1]))/(1-gamma)) ) )
      }
    }
    s3 = sum(na.omit(s3b))

    s4=0
    if (Ln[m] < n) {
      s4 = sum(log(cdf(Rn[m], par = list(shape=shape, scale=scale))^((gamma^(Ln[m]+1) - gamma^(n+1))/(1-gamma) ) ))
      }
     return(s1+s2+s3+s4)
    }
)

  ## Gumbel
register_loglik( "YNM", "records", "gumbel",
                 fun = function(data, params) {
                   if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                   Rn <- data$rec_values
                   Ln <- data$rec_times
                   n = data$time[1]
                   if( all(c("gamma", "loc", "scale") %in% names(params)) == FALSE ) stop("parameters gamma, loc, scale should be present.")

                   gamma <- as.numeric(params["gamma"])
                   loc <- as.numeric(params["loc"])
                   scale <- as.numeric(params["scale"])
                   m = length(Rn) #m= rec_counts(y)  ## number of records

                   # Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or shape <= 0, we return -Inf)
                   if (scale <= 0) {
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## pdf
                   pdf=function(x,par) {VGAM::dgumbel(x=x, location=par$loc, scale=par$scale)}
                   ## cdf
                   cdf=function(x,par) {VGAM::pgumbel(q=x, location=par$loc, scale=par$scale)}


                   # your exact YNM record-pair likelihood
                   s1 = log(gamma^(sum(Ln)))

                   s2 = sum(log(pdf(Rn, par = list(loc = loc, scale=scale)) * cdf(Rn, par = list(loc=loc, scale=scale))^(-1+gamma^Ln)))

                   s3b = 0
                   for( i in 1:(m-1)  ){
                     if((Ln[i]+1 <= Ln[i+1]-1) == TRUE){ ## we have non-records in between
                       s3b[i]= sum(log(cdf(Rn[i], par = list(loc = loc, scale=scale))^( (gamma^(Ln[i]+1) - gamma^(Ln[i+1]))/(1-gamma)) ) )
                     }
                   }
                   s3 = sum(na.omit(s3b))

                   s4=0
                   if (Ln[m] < n) {
                     s4 = sum(log(cdf(Rn[m], par = list(loc = loc, scale=scale))^( (gamma^(Ln[m]+1) - gamma^(n+1))/(1-gamma) ) ))
                   }
                   return(s1+s2+s3+s4)
                 }
)

  ## Norm
register_loglik( "YNM", "records", "norm",
                 fun = function(data, params) {
                   if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                   Rn <- data$rec_values
                   Ln <- data$rec_times
                   n = data$time[1]
                   if( all(c("gamma", "mean", "sd") %in% names(params)) == FALSE ) stop("parameters gamma, mean, sd should be present.")

                   gamma <- as.numeric(params["gamma"])
                   mean <- as.numeric(params["mean"])
                   sd <- as.numeric(params["sd"])
                   m = length(Rn) #m= rec_counts(y)  ## number of records

                   # Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or shape <= 0, we return -Inf)
                   if (sd <= 0) {
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## pdf
                   pdf=function(x,par) {dnorm(x=x, mean=par$mean, sd=par$sd)}
                   ## cdf
                   cdf=function(x,par) {pnorm(q=x, mean=par$mean, sd=par$sd)}


                   # your exact YNM record-pair likelihood
                   s1 = log(gamma^(sum(Ln)))
                   s2 = sum(log(pdf(Rn, par = list(mean = mean, sd=sd)) * cdf(Rn, par = list(mean=mean, sd=sd))^(-1+gamma^Ln)))
                   s3b = 0
                   for( i in 1:(m-1)  ){
                     if((Ln[i]+1 <= Ln[i+1]-1) == TRUE){ ## we have non-records in between
                       s3b[i]= sum(log(cdf(Rn[i], par = list(mean = mean, sd=sd))^( (gamma^(Ln[i]+1) - gamma^(Ln[i+1]-1))/(1-gamma)) ) )
                     }
                   }
                   s3 = sum(na.omit(s3b))

                   s4=0
                   if (Ln[m] < n) {
                     s4 = sum(log(cdf(Rn[m], par = list(mean = mean, sd=sd))^( (gamma^(Ln[m]+1) - gamma^(n+1))/(1-gamma) ) ))
                   }
                   return(s1+s2+s3+s4)
                 }
)

  ## weibull
register_loglik( "YNM", "records", "weibull",
                 fun = function(data, params) {
                   if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                   Rn <- data$rec_values
                   Ln <- data$rec_times
                   n <- data$time[1]
                   if( all(c("gamma", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters gamma, shape, scale should be present.")

                   gamma <- as.numeric(params["gamma"])
                   shape <- as.numeric(params["shape"])
                   scale <- as.numeric(params["scale"])
                   m = length(Rn) #m= rec_counts(y)  ## number of records

                   # Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or shape <= 0, we return -Inf)
                   if (scale <= 0 | shape<=0) {
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## pdf
                   pdf=function(x,par) {dweibull(x=x, shape=par$shape, scale=par$scale)}
                   ## cdf
                   cdf=function(x,par) {pweibull(q=x, shape= par$shape, scale=par$scale)}


                   # your exact YNM record-pair likelihood
                   s1 = log(gamma^(sum(Ln)))
                   s2 = sum(log(pdf(Rn, par = list(shape=shape, scale=scale)) * cdf(Rn, par = list(shape=shape, scale=scale))^(-1+gamma^Ln)))
                   s3b = 0
                   for( i in 1:(m-1)  ){
                     if((Ln[i]+1 <= Ln[i+1]) == TRUE){ ## we have non-records in between
                       s3b[i]= sum(log(cdf(Rn[i], par = list(shape=shape, scale=scale))^( (gamma^(Ln[i]+1) - gamma^(Ln[i+1]))/(1-gamma)) ) )
                     }
                   }
                   s3 = sum(na.omit(s3b))

                   s4=0
                   if (Ln[m] < n) {
                     s4 = sum(log(cdf(Rn[m], par = list(shape=shape, scale=scale))^((gamma^(Ln[m]+1) - gamma^(n+1))/(1-gamma) ) ))
                   }
                   return(s1+s2+s3+s4)
                 }
)
## LDM, Xt -----------------------------------------

  ##Norm
register_loglik( "LDM", "all", "norm",
                 fun = function(data, params) {
                   y=data
                   if(!is.numeric(y)) stop("data should be a numerical vector")

                   if( all(c("theta", "mean", "sd") %in% names(params)) == FALSE ) stop("parameterstheta, mean, sd should be present.")

                   theta <- as.numeric(params["theta"])
                   mean <- as.numeric(as.numeric(params["mean"]))
                   sd <- as.numeric(params["sd"])

                   ## Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or shape <= 0, we return -Inf)
                   if ( sd <= 0) { ##
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## transform
                   x = y - theta * (1:length(y)) ## y-theta*t

                   ## pdf
                   pdf=function(x,par) {dnorm(x=x, mean=par$mean, sd=par$sd)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum(log(pdf(x=x, par=list(mean = mean, sd = sd))))

                   return(s1)
                 }
)
  ##Frechet
register_loglik( "LDM", "all", "frechet",
                 fun = function(data, params) {
                   y=data
                   if(!is.numeric(y)) stop("data should be a numerical vector")

                   if( all(c("theta", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters theta, shape, scale should be present.")

                   theta <- as.numeric(params["theta"])
                   shape <- as.numeric(params["shape"])
                   scale <- as.numeric(params["scale"])

                   ## Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or shape <= 0, we return -Inf)
                   if (scale <= 0) { ##
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## transform
                   x = y - theta * (1:length(y)) ## y-theta*t

                   ## pdf
                   pdf=function(x,par) {VGAM::dfrechet(x=x, shape=par$shape, scale=par$scale)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum(log(pdf(x=x, par=list(shape= shape, scale= scale))))

                   # scale = 1/params["scale"]
                   # if (any(x <= 0)) return(-1000)
                   # A <- (sum(x^-scale) / n)^(1/scale)
                   # s1 <- n * log(scale * A^(-scale))
                   # s2 <- -(scale + 1) * sum(log(x))
                   # s3 <- -sum((A * x)^(-scale))
                   # return (s1 + s2 + s3)
                   return(s1)
                 }
)

  ##gumbel
register_loglik( "LDM", "all", "gumbel",
                 fun = function(data, params) {
                   y=data
                   if(!is.numeric(y)) stop("data should be a numerical vector")

                   if( all(c("theta", "loc", "scale") %in% names(params)) == FALSE ) stop("parameterstheta, loc, scale should be present.")
                   theta <- as.numeric(params["theta"])
                   loc <- as.numeric(params["loc"])
                   scale <- as.numeric(params["scale"])

                   ## Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or shape <= 0, we return -Inf)
                   if (scale <= 0) { ##
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## transform
                   x = y - theta * (1:length(y)) ## y-theta*t

                   ## pdf
                   pdf=function(x,par) {VGAM::dgumbel(x=x, location=par$loc, scale=par$scale)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum(log(pdf(x=x, par=list(loc = loc, scale= scale))))
                   return(s1)
                 }
)

## LDM, Rn -----------------------------------------

    ##Gumbel
register_loglik( "LDM", "records", "gumbel",
                 fun = function(data, params) {

                   if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                   Rn <- data$rec_values
                   Ln <- data$rec_times
                   n <- data$time[1]

                   if( all(c("theta", "loc", "scale") %in% names(params)) == FALSE ) stop("parameterstheta, loc, scale should be present.")
                   theta <- as.numeric(params["theta"])
                   loc <- as.numeric(params["loc"])
                   scale <- as.numeric(params["scale"])
                   m = length(Rn)

                   ## Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or shape <= 0, we return -Inf)
                   if ( scale <= 0) { ##
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## transform
                   x = Rn - theta * Ln ## Rn-theta*ln

                   ## pdf
                   pdf=function(x,par) {VGAM::dgumbel(x=x, loc=par$loc, scale=par$scale)}
                   ## cdf
                   cdf=function(x,par) {VGAM::pgumbel(q=x, loc=par$loc, scale=par$scale)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum(log(pdf(x=x, par=list(loc = loc, scale = scale))))
                   if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                   s2b=0
                   for( i in 1:(m-1)  ){
                     if((Ln[i]+1 <= Ln[i+1]-1) == TRUE){ ## we have non-records in between
                       s2b[i]=sum(log(cdf( Rn[i]-theta * (Ln[i]+1):(Ln[i+1]-1), par=list(loc = loc, scale = scale) )  ))
                     }
                   }
                   s2 = sum(na.omit(s2b))
                   if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                   ## case where last record is not last observation
                   s3=0
                   if( (Ln[m] < n) == TRUE  ) {
                     s3a=log(cdf( Rn[m]-theta * (Ln[m]+1):n, par=list(loc = loc, scale = scale) )  )
                     s3 = sum(s3a[is.finite(s3a)])
                   }
                   if (is.nan(s3) || !is.finite(s3)) return(-Inf)

                   return(s1+s2+s3)
                 }
)
    ##normal
register_loglik( "LDM", "records", "norm",
                 fun = function(data, params) {
                   if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                   Rn <- data$rec_values
                   Ln <- data$rec_times
                   n <- data$time[1]

                   if( all(c("theta", "mean", "sd") %in% names(params)) == FALSE ) stop("parameterstheta, mean, sd should be present.")
                   theta <- as.numeric(params["theta"])
                   mean <- as.numeric(params["mean"])
                   sd <- as.numeric(params["sd"])
                   m = length(Rn)

                   ## Check for invalid parameter values
                   if (sd <= 0) {
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## transform
                   x = Rn - theta * Ln ## Rn-theta*ln

                   ## pdf
                   pdf=function(x,par) {dnorm(x=x, mean=par$mean, sd=par$sd)}
                   ## cdf
                   cdf=function(x,par) {pnorm(q=x, mean=par$mean, sd=par$sd)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum(log(pdf(x=x, par=list(mean = mean, sd = sd))))
                   if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                   s2b=0
                   for( i in 1:(m-1)  ){
                     if((Ln[i]+1 <= Ln[i+1]-1) == TRUE){ ## we have non-records in between
                       s2b[i]=sum(log(cdf( Rn[i]-theta * (Ln[i]+1):(Ln[i+1]-1), par=list(mean = mean, sd = sd) )  ))
                     }
                   }
                   s2 = sum(na.omit(s2b))
                   if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                   ## case where last record is not last observation
                   s3=0
                   if( (Ln[m] < n) == TRUE  ) {
                     s3a= log(cdf( Rn[m]-theta * (Ln[m]+1):n, par=list(mean = mean, sd = sd) )  )
                     s3 = sum(s3a[is.finite(s3a)])
                     }
                   if (is.nan(s3) || !is.finite(s3)) return(-Inf)

                   return(s1+s2+s3)
                 }
)

##frechet
register_loglik( "LDM", "records", "frechet",
                 fun = function(data, params) {
                   if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                   Rn <- data$rec_values
                   Ln <- data$rec_times
                   n <- data$time[1]

                   if( all(c("theta", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters theta, shape, scale should be present.")
                   theta <- as.numeric(params["theta"])
                   shape <- as.numeric(params["shape"])
                   scale <- as.numeric(params["scale"])
                   m = length(Rn)

                   ## Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or shape <= 0, we return -Inf)
                   if ( shape <0 |scale <= 0) { ##
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## transform
                   x = Rn - theta * Ln ## Rn-theta*ln

                   ## pdf
                   pdf=function(x,par) {VGAM::dfrechet(x=x, shape=par$shape, scale=par$scale)}
                   ## cdf
                   cdf=function(x,par) {VGAM::pfrechet(q=x, shape=par$shape, scale=par$scale)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum(log(pdf(x=x, par=list(shape = shape, scale = scale))))
                   if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                   s2b=0
                   for( i in 1:(m-1)  ){
                     if((Ln[i]+1 <= Ln[i+1]-1) == TRUE){ ## we have non-records in between
                       s2b[i]=sum(log(cdf( Rn[i]-theta * (Ln[i]+1):(Ln[i+1]-1), par=list(shape = shape, scale = scale) )  ))
                     }
                   }
                   s2 = sum(na.omit(s2b))
                   if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                   ## case where last record is not last observation
                   s3=0
                   if( (Ln[m] < n) == TRUE  ) {
                     s3a=log(cdf( Rn[m]-theta * (Ln[m]+1):n, par=list(shape = shape, scale = scale) )  )
                     s3 = sum(s3a[is.finite(s3a)])
                   }
                   if (is.nan(s3) || !is.finite(s3)) return(-Inf)

                   return(s1+s2+s3)
                 }
)
