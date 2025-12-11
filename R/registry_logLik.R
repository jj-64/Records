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
    location = ifelse("location"  %in% names(params), as.numeric(params["location"]) ,(min(data)-1e-6))

    sum(VGAM::dfrechet(x = data, location = location, shape = as.numeric(params["shape"]), scale = params["scale"] ,log = TRUE))
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

                   return(sum(dweibull(x = data, shape = as.numeric(params["shape"]), scale = as.numeric(params["scale"]), log= TRUE )))
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
    s1 <- sum(dnorm(Rn[-1], mean = mean, sd = sd, log = TRUE))
    if (is.nan(s1) || !is.finite(s1)) return(-Inf)
    #Compute s2 using vectorized approach
    intervals <- diff(Ln)-1
    s2b <- intervals * pnorm(Rn[-m], mean = mean, sd = sd, log = TRUE)
    s2 <- sum(s2b)
    if (is.nan(s2) || !is.finite(s2)) return(-Inf)

    ##Compute s3 only if needed
    s3 <- if (Ln[m] < n) {
      (n - Ln[m]) * pnorm(Rn[m], mean = mean, sd = sd, log = TRUE)
    } else {
      0
    }
    if (is.nan(s3) || !is.finite(s3)) return(-Inf)

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
                  location = ifelse("location"  %in% names(params), as.numeric(params["location"]) ,(min(Rn)-1e-6))
                  m = length(Rn) #m= rec_counts(y)

                  ## sum of log(f_rn)
                  s1 <- sum(VGAM::dfrechet(Rn[-1],  location = location, shape = shape, scale= scale, log = TRUE))
                  if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                  #Compute s2 using vectorized approach
                  intervals <- diff(Ln)-1
                  s2b <- intervals * VGAM::pfrechet(Rn[-m], location = location, shape = shape, scale= scale, log = TRUE)
                  s2 <- sum(s2b)
                  if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                  ##Compute s3 only if needed
                  s3 <- if (Ln[m] < n) {
                    (n - Ln[m]) * VGAM::pfrechet(Rn[m], location = location,  shape = shape, scale= scale, log = TRUE)
                  } else {
                    0
                  }
                  if (is.nan(s3) || !is.finite(s3)) return(-Inf)

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
                  s1 <- sum(VGAM::dgumbel(Rn[-1], loc, scale, log = TRUE))
                  if (is.nan(s1) || !is.finite(s1)) return(-Inf)
                  #Compute s2 using vectorized approach
                  intervals <- diff(Ln)-1
                  s2b <- intervals * VGAM::pgumbel(Rn[-m], loc, scale, log = TRUE)
                  s2 <- sum(s2b)
                  if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                  ##Compute s3 only if needed
                  s3 <- if (Ln[m] < n) {
                    (n - Ln[m]) * VGAM::pgumbel(Rn[m], loc, scale, log = TRUE)
                  } else {
                    0
                  }
                  if (is.nan(s3) || !is.finite(s3)) return(-Inf)

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
                  s1 <- sum(dweibull(Rn[-1], shape, scale, log = TRUE))
                  if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                  #Compute s2 using vectorized approach
                  intervals <- diff(Ln)-1
                  s2b <- intervals * pweibull(Rn[-m], shape, scale, log = TRUE)
                  s2 <- sum(s2b)
                  if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                  ##Compute s3 only if needed
                  s3 <- if (Ln[m] < n) {
                    (n - Ln[m]) * pweibull(Rn[m], shape, scale, log = TRUE)
                  } else {
                    0
                  }
                  if (is.nan(s3) || !is.finite(s3)) return(-Inf)

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

    sum(dnorm(x = diff(data), mean =mean , sd = sd, log = TRUE))
  }
)

  ##cauchy
register_loglik("DTRW", "all", "cauchy",
                fun = function(data, params) {
                  if(!is.numeric(data)) stop("data should be a numerical vector")
                  if( all(c("loc", "scale") %in% names(params)) == FALSE ) stop("parameters loc, scale should be present in a list.")

                  sum(dcauchy(x = diff(data), location =as.numeric(params["loc"]) , scale = as.numeric(params["scale"]) , log = TRUE))
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
    if( all(c( "sd") %in% names(params)) == FALSE ) stop("parameters sd should be present in a list.")

    mean = 0
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

    s1b = dnorm(diff(Rn), mean = 0, sd = sqrt(diff(Ln)*sd), log = TRUE )
    s1=sum(s1b)
    if (is.nan(s1) || !is.finite(s1)) return(-Inf)

    ## case where NT<T
    s2=0
    if (Ln[m] < n) {
      #s2=log(mvtnorm::pmvnorm(lower=-Inf, upper=rep(0,(n - Ln[m])), sigma = sigma_cov((n - Ln[m]),params))[1])
      s2 = -0.5*log(pi * (n - Ln[m]))
    }
    if (is.nan(s2) || !is.finite(s2)) return(-Inf)

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

                   s1b = dcauchy(diff(Rn), location = 0, scale = sqrt(diff(Ln)*scale), log = TRUE)
                   s1=sum(s1b)
                   if (is.nan(s1) || !is.finite(s1)) return(-Inf)
                   ## case where NT<T
                   s2=0
                   if (Ln[m] < n) { s2 = -0.5*log(pi * (n - Ln[m])) }
                   if (is.nan(s2) || !is.finite(s2)) return(-Inf)

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
                   location = ifelse("location"  %in% names(params), as.numeric(params["location"]) ,(min(x)-1e-6))

                   ## Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or shape <= 0, we return -Inf)
                   if (shape <=0 | scale <= 0) { ##
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## pdf
                   pdf=function(x,par) {VGAM::dfrechet(x, location = par$location, shape=par$shape, scale=par$scale, log = TRUE)}
                   cdf=function(x,par) {VGAM::pfrechet(x, location = par$location, shape=par$shape, scale=par$scale, log = TRUE)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 =  sum( (1:n) * log(gamma))
                   if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                   s2 = sum( (-1+gamma^(1:n))  * cdf(x, par=list(location = location, shape= shape, scale= scale))  )
                   if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                   s3 = sum(pdf(x, par=list(location = location, shape= shape, scale= scale)))
                   if (is.nan(s3) || !is.finite(s3)) return(-Inf)
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
                   tvec = seq_len(n)

                   ## Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or loc <= 0, we return -Inf)
                   if (scale <= 0 || gamma <= 0) { ##
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## pdf
                   pdf=function(x,par) {VGAM::dgumbel(x, location=par$loc, scale=par$scale, log = TRUE)}
                   cdf=function(x,par) {VGAM::pgumbel(x, location=par$loc, scale=par$scale, log = TRUE)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum( tvec * log(gamma) )
                   if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                   s2 = sum( (gamma^tvec - 1) * cdf(x, par=list(loc= loc, scale= scale) ) )
                   if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                   s3 = sum(pdf(x, par=list(loc= loc, scale= scale)))
                   if (is.nan(s3) || !is.finite(s3)) return(-Inf)

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
                   pdf=function(x,par) {dnorm(x, mean=par$mean, sd=par$sd, log = TRUE)}
                   cdf=function(x,par) {pnorm(x, mean=par$mean, sd=par$sd, log = TRUE)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum( (1:n) * log(gamma))
                   if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                   s2 = sum((-1+gamma^(1:n)) * cdf(x, par=list(mean= mean, sd= sd)) )
                   if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                   s3 = sum( pdf(x, par=list(mean= mean, sd= sd)) )
                   if (is.nan(s3) || !is.finite(s3)) return(-Inf)


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
                   if (shape <=0 || scale <= 0) { ##
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## pdf
                   pdf=function(x,par) {dweibull(x, shape=par$shape, scale=par$scale, log = TRUE)}
                   cdf=function(x,par) {pweibull(x, shape=par$shape, scale=par$scale, log = TRUE)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   tvec = seq_len(n)
                   s1 = sum( tvec * log(gamma) )
                   if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                   s2 = sum( (gamma^tvec - 1) * cdf(x, par=list(shape= shape, scale= scale) ) )
                   if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                   s3 = sum(pdf(x, par=list(shape= shape, scale= scale)))
                   if (is.nan(s3) || !is.finite(s3)) return(-Inf)

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
    location = ifelse("location"  %in% names(params), as.numeric(params["location"]) ,(min(Rn)-1e-6))

    # Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or shape <= 0, we return -Inf)
    if (scale <= 0 || shape<=0 || gamma <0) {
      return(-Inf)  # Invalid parameters, return a large negative value
    }

    ## pdf
    pdf=function(x,par) {VGAM::dfrechet(x=x, location = par$location, shape=par$shape, scale=par$scale, log = TRUE)}
    ## cdf
    cdf=function(x,par) {VGAM::pfrechet(q=x, location = par$location, shape= par$shape, scale=par$scale, log = TRUE)}


    # your exact YNM record-pair likelihood
    s1 = sum(Ln) *log(gamma)
    if (is.nan(s1) || !is.finite(s1)) return(-Inf)

    s2 = sum( pdf(Rn, par = list(location = location, shape=shape, scale=scale)) + (-1+gamma^Ln) * cdf(Rn, par = list(location = location, shape=shape, scale=scale)))
    if (is.nan(s2) || !is.finite(s2)) return(-Inf)
    s3b = 0
    for( i in 1:(m-1)  ){
      if((Ln[i]+1 <= Ln[i+1]) == TRUE){ ## we have non-records in between
        s3b[i]= sum( ( (gamma^(Ln[i]+1) - gamma^(Ln[i+1]))/(1-gamma)) * cdf(Rn[i], par = list(location = location, shape=shape, scale=scale)) )
      }
    }
    s3 = sum(na.omit(s3b))
    if (is.nan(s3) || !is.finite(s3)) return(-Inf)

    s4=0
    if (Ln[m] < n) {
      s4 = sum( ((gamma^(Ln[m]+1) - gamma^(n+1))/(1-gamma) ) * cdf(Rn[m], par = list(location = location, shape=shape, scale=scale) ))
    }
    if (is.nan(s4) || !is.finite(s4)) return(-Inf)

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
                   pdf=function(x,par) {VGAM::dgumbel(x=x, location=par$loc, scale=par$scale, log = TRUE)}
                   ## cdf
                   cdf=function(x,par) {VGAM::pgumbel(q=x, location=par$loc, scale=par$scale, log = TRUE)}


                   # your exact YNM record-pair likelihood
                   s1 = sum(Ln) * log(gamma)
                   if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                   s2 = sum(pdf(Rn, par = list(loc = loc, scale=scale)) + (-1+gamma^Ln) * cdf(Rn, par = list(loc=loc, scale=scale)))
                   if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                   s3b = 0
                   for( i in 1:(m-1)  ){
                     if((Ln[i]+1 <= Ln[i+1]-1) == TRUE){ ## we have non-records in between
                       s3b[i]= sum(( (gamma^(Ln[i]+1) - gamma^(Ln[i+1]))/(1-gamma)) * cdf(Rn[i], par = list(loc = loc, scale=scale)) )
                     }
                   }
                   s3 = sum(na.omit(s3b))
                   if (is.nan(s3) || !is.finite(s3)) return(-Inf)

                   s4=0
                   if (Ln[m] < n) {
                     s4 = sum( ( (gamma^(Ln[m]+1) - gamma^(n+1))/(1-gamma) ) * cdf(Rn[m], par = list(loc = loc, scale=scale)) )
                   }
                   if (is.nan(s4) || !is.finite(s4)) return(-Inf)

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
                   pdf=function(x,par) {dnorm(x=x, mean=par$mean, sd=par$sd, log = TRUE)}
                   ## cdf
                   cdf=function(x,par) {pnorm(q=x, mean=par$mean, sd=par$sd, log = TRUE)}


                   # your exact YNM record-pair likelihood
                   s1 = sum(Ln) * log(gamma)
                   if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                   s2 = sum(pdf(Rn, par = list(mean = mean, sd=sd)) + (-1+gamma^Ln) *  cdf(Rn, par = list(mean=mean, sd=sd)))
                   if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                   s3b = 0
                   for( i in 1:(m-1)  ){
                     if((Ln[i]+1 <= Ln[i+1]-1) == TRUE){ ## we have non-records in between
                       s3b[i]= sum( ( (gamma^(Ln[i]+1) - gamma^(Ln[i+1]-1))/(1-gamma)) * cdf(Rn[i], par = list(mean = mean, sd=sd)) )
                     }
                   }
                   s3 = sum(na.omit(s3b))
                   if (is.nan(s3) || !is.finite(s3)) return(-Inf)

                   s4=0
                   if (Ln[m] < n) {
                     s4 = sum(( (gamma^(Ln[m]+1) - gamma^(n+1))/(1-gamma) ) * cdf(Rn[m], par = list(mean = mean, sd=sd)) )
                   }
                   if (is.nan(s4) || !is.finite(s4)) return(-Inf)

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
                   pdf=function(x,par) {dweibull(x=x, shape=par$shape, scale=par$scale, log = TRUE)}
                   ## cdf
                   cdf=function(x,par) {pweibull(q=x, shape= par$shape, scale=par$scale, log = TRUE)}


                   # your exact YNM record-pair likelihood
                   s1 = log(gamma^(sum(Ln)))
                   if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                   s2 = sum( pdf(Rn, par = list(shape=shape, scale=scale)) + (-1 + gamma^Ln) * cdf(Rn, par = list(shape=shape, scale=scale)) )
                   if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                   s3b = 0
                   for( i in 1:(m-1)  ){
                     if((Ln[i]+1 <= Ln[i+1]) == TRUE){ ## we have non-records in between
                       s3b[i]= sum( ( (gamma^(Ln[i]+1) - gamma^(Ln[i+1]))/(1-gamma)) * cdf(Rn[i], par = list(shape=shape, scale=scale) ) )
                     }
                   }
                   s3 = sum(na.omit(s3b))
                   if (is.nan(s3) || !is.finite(s3)) return(-Inf)

                   s4=0
                   if (Ln[m] < n) {
                     s4 = sum( ((gamma^(Ln[m]+1) - gamma^(n+1))/(1-gamma) ) * cdf(Rn[m], par = list(shape=shape, scale=scale)) )
                   }
                   if (is.nan(s4) || !is.finite(s4)) return(-Inf)

                   return(s1+s2+s3+s4)
                 }
)
## LDM, Xt -----------------------------------------

  ##Norm
register_loglik( "LDM", "all", "norm",
                 fun = function(data, params) {
                   y=data
                   if(!is.numeric(y)) stop("data should be a numerical vector")

                   if( all(c("theta", "mean", "sd") %in% names(params)) == FALSE ) stop("parameters theta, mean, sd should be present.")

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
                   pdf=function(x,par) {dnorm(x=x, mean=par$mean, sd=par$sd, log = TRUE)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum(pdf(x=x, par=list(mean = mean, sd = sd)))

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
                   location = ifelse("location"  %in% names(params), as.numeric(params["location"]) ,(min(x)-1e-6))

                   ## pdf
                   pdf=function(x,par) {VGAM::dfrechet(x=x, location = par$location, shape=par$shape, scale=par$scale, log = TRUE)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum(pdf(x=x, par=list(location = location, shape= shape, scale= scale)))

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
                   pdf=function(x,par) {VGAM::dgumbel(x=x, location=par$loc, scale=par$scale ,log = TRUE)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum(pdf(x=x, par=list(loc = loc, scale= scale)))
                   return(s1)
                 }
)

##weibull
register_loglik( "LDM", "all", "weibull",
                 fun = function(data, params) {
                   y=data
                   if(!is.numeric(y)) stop("data should be a numerical vector")

                   if( all(c("theta", "shape", "scale") %in% names(params)) == FALSE ) stop("parameterstheta, shape, scale should be present.")
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
                   pdf=function(x,par) {dweibull(x=x, shape=par$shape, scale=par$scale, log = TRUE)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum(pdf(x=x, par=list(shape = shape, scale= scale)))
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
                   pdf=function(x,par) {VGAM::dgumbel(x=x, loc=par$loc, scale=par$scale , log = TRUE)}
                   ## cdf
                   cdf=function(x,par) {VGAM::pgumbel(q=x, loc=par$loc, scale=par$scale, log = TRUE)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum(pdf(x=x, par=list(loc = loc, scale = scale)))
                   if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                   s2b=0
                   for( i in 1:(m-1)  ){
                     if((Ln[i]+1 <= Ln[i+1]-1) == TRUE){ ## we have non-records in between
                       s2b[i]=sum(cdf( Rn[i]-theta * (Ln[i]+1):(Ln[i+1]-1), par=list(loc = loc, scale = scale) )  )
                     }
                   }
                   s2 = sum(na.omit(s2b))
                   if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                   ## case where last record is not last observation
                   s3=0
                   if( (Ln[m] < n) == TRUE  ) {
                     s3a=cdf( Rn[m]-theta * (Ln[m]+1):n, par=list(loc = loc, scale = scale) )
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

                   if( all(c("theta", "mean", "sd") %in% names(params)) == FALSE ) stop("parameters theta, mean, sd should be present.")
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
                   pdf=function(x,par) {dnorm(x=x, mean=par$mean, sd=par$sd, log = TRUE)}
                   ## cdf
                   cdf=function(x,par) {pnorm(q=x, mean=par$mean, sd=par$sd, log = TRUE)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum(pdf(x=x, par=list(mean = mean, sd = sd)))
                   if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                   s2b=0
                   for( i in 1:(m-1)  ){
                     if((Ln[i]+1 <= Ln[i+1]-1) == TRUE){ ## we have non-records in between
                       s2b[i]=sum(cdf( Rn[i]-theta * (Ln[i]+1):(Ln[i+1]-1), par=list(mean = mean, sd = sd) )  )
                     }
                   }
                   s2 = sum(na.omit(s2b))
                   if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                   ## case where last record is not last observation
                   s3=0
                   if( (Ln[m] < n) == TRUE  ) {
                     s3a= cdf( Rn[m]-theta * (Ln[m]+1):n, par=list(mean = mean, sd = sd) )
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
                   location = ifelse("location"  %in% names(params), as.numeric(params["location"]) ,(min(x)-1e-6))

                   ## pdf
                   pdf=function(x,par) {VGAM::dfrechet(x=x, location = par$location, shape=par$shape, scale=par$scale, log = TRUE)}

                   ## cdf
                   cdf=function(x,par) {VGAM::pfrechet(q=x, location = par$location, shape=par$shape, scale=par$scale, log = TRUE)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum((pdf(x=x, par=list(location = location, shape = shape, scale = scale))))
                   if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                   s2b=0
                   for( i in 1:(m-1)  ){
                     if((Ln[i]+1 <= Ln[i+1]-1) == TRUE){ ## we have non-records in between
                       s2b[i]=sum(  cdf( Rn[i]- theta * (Ln[i]+1):(Ln[i+1]-1), par=list(location= location, shape = shape, scale = scale) )  )
                     }
                   }
                   s2b = s2b[is.finite(s2b)]
                   s2 = sum(na.omit(s2b))

                   if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                   ## case where last record is not last observation
                   s3=0
                   if( (Ln[m] < n) == TRUE  ) {
                     s3a=(cdf( Rn[m]-theta * (Ln[m]+1):n, par=list(location = location , shape = shape, scale = scale) )  )
                     s3 = sum(s3a[is.finite(s3a)])
                   }
                   if (is.nan(s3) || !is.finite(s3)) return(-Inf)

                   return(s1+s2+s3)
                 }
)

##Weibull
register_loglik( "LDM", "records", "weibull",
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
                   pdf=function(x,par) {dweibull(x=x,  shape=par$shape, scale=par$scale, log = TRUE)}

                   ## cdf
                   cdf=function(x,par) {pweibull(q=x,  shape=par$shape, scale=par$scale, log = TRUE)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum((pdf(x=x, par=list( shape = shape, scale = scale))))
                   if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                   s2b=0
                   for( i in 1:(m-1)  ){
                     if((Ln[i]+1 <= Ln[i+1]-1) == TRUE){ ## we have non-records in between
                       s2b[i]=sum(  cdf( Rn[i]- theta * (Ln[i]+1):(Ln[i+1]-1), par=list( shape = shape, scale = scale) )  )
                     }
                   }
                   s2b = s2b[is.finite(s2b)]
                   s2 = sum(na.omit(s2b))

                   if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                   ## case where last record is not last observation
                   s3=0
                   if( (Ln[m] < n) == TRUE  ) {
                     s3a=(cdf( Rn[m]-theta * (Ln[m]+1):n, par=list( shape = shape, scale = scale) )  )
                     s3 = sum(s3a[is.finite(s3a)])
                   }
                   if (is.nan(s3) || !is.finite(s3)) return(-Inf)

                   return(s1+s2+s3)
                 }
)
