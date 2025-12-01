var_logLik_records <- function(model, obs_type, dist, param_name) {
  # lookup
  f <- var_registry[[model]][[obs_type]][[dist]][[param_name]]

  if (is.null(f))
    stop("variance expression not registered for this (model, obs_type, dist, param_name).")

  #f(data, params)
  return(f)
}

## Helper ---------------------------------
# empty container
var_registry <- new.env(parent = emptyenv())

## Function to register likelihood expressions later
register_var <- function(model, obs_type, dist, param_name, fun) {
  if (!exists(model, envir = var_registry))
    var_registry[[model]] <- list()

  if (!obs_type %in% names(var_registry[[model]]))
    var_registry[[model]][[obs_type]] <- list()

  var_registry[[model]][[obs_type]][[dist]][[param_name]] <- fun
}

## Classical, Xt ---------------

register_var( model = "iid", obs_type = "all", dist = "frechet", param_name = "shape",
              fun = function(data, params) {

                if(!is.numeric(data)) stop("data should be a numerical vector")
                n <- length(data)

                if( all(c("shape", "scale") %in% names(params)) == FALSE ) stop("parameters shape, scale should be present.")

                shape <- as.numeric(params$shape)
                scale <- 1/as.numeric(params$scale)

                s1=-n/(shape^2)

                s2 = - sum( log(scale*data)^2/(scale*data)^shape )

                var_a = -1 / (s1+s2)

                return(var_a)


              }
)

register_var( model = "iid", obs_type = "all", dist = "frechet", param_name = "scale",
              fun = function(data, params) {

                if(!is.numeric(data)) stop("data should be a numerical vector")
                n <- length(data)

                if( all(c("shape", "scale") %in% names(params)) == FALSE ) stop("parameters shape, scale should be present.")

                shape <- as.numeric(params$shape)
                scale <- 1/as.numeric(params$scale)

                s1 = (shape*n)/scale^2

                s2 = -((shape^2+shape)/scale^2) * sum((data*scale)^(-shape))

                var_A <- -1 / (s1 + s2)

              }
)
## Classical, Rn ---------------

register_var( model = "iid", obs_type = "records", dist = "frechet", param_name = "scale",
              fun = function(data, params) {
                if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                Rn <- data$rec_values
                Ln <- data$rec_times
                T <- data$time

                if( all(c("shape", "scale") %in% names(params)) == FALSE ) stop("parameters shape, scale should be present.")

                shape <- as.numeric(params$shape)
                scale <- 1/as.numeric(params$scale)

                m = length(Rn)

                s1=-(m-1)/(scale^2)

                s2 = -((shape^2+shape)/scale^2) * sum((Rn[-1]*scale)^(-shape))

                s3 = (shape+1)*(m-1)/scale^2

                s4a=0
                for (i in 1:(m - 1)) {
                  s4a[i] = (Ln[i+1]-Ln[i]-1)*(scale*Rn[i])^(-shape)
                }
                s4 = -(shape^2+shape)*sum(s4a)/scale^2

                ## case where last record is not last observation
                s5=0
                if( (Ln[m] < T) == TRUE) {
                  s5 = -(shape*(shape+1))*(T-Ln[m])*(scale*Rn[m])^(-shape)/scale^2
                }

                var_A <- -1 / (s1 + s2 + s3 + s4 + s5)

                return(var_A)
              }
)

register_var( model = "iid", obs_type = "records", dist = "frechet", param_name = "shape",
              fun = function(data, params) {
                if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                Rn <- data$rec_values
                Ln <- data$rec_times
                T <- data$time

                if( all(c("shape", "scale") %in% names(params)) == FALSE ) stop("parameters shape, scale should be present.")

                shape <- as.numeric(params$shape)
                scale <- 1/as.numeric(params$scale)

                m = length(Rn)

                s1 = -(m - 1) / shape^2
                s2 = -sum((log(scale * Rn[-1]))^2 / (scale * Rn[-1])^shape)

                s3a = 0
                for (i in 1:(m - 1)) {
                  s3a[i] = (Ln[i + 1] - Ln[i] - 1) * (scale * Rn[i])^(-shape) * (log(scale * Rn[i]))^2
                }
                s3 = -sum(s3a)

                ## Case where last record is not the last observation
                s4 = 0
                if (Ln[m] < T) {
                  s4 = -(T - Ln[m]) * (scale * Rn[m])^(-shape) * (log(scale * Rn[m]))^2
                }

                var_a = -1 / (s1 + s2 + s3 + s4)

                return(var_a)
              }
)

## DTRW, Xt -------------
  ##norm
register_var( model = "iid", obs_type = "all", dist = "norm", param_name = "sd",
              fun = function(data, params) {
                if(!is.numeric(data)) stop("data should be a numerical vector")

                if( all(c("sd") %in% names(params)) == FALSE ) stop("parameter sd should be present.")

                n = length(data)
                sd <- 1/as.numeric(params$sd)

                m = length(Rn)
                return( 2 * sd^2 / (1:n) )
              }
)
## DTRW, Rn -------------
  ##norm
register_var( model = "iid", obs_type = "records", dist = "norm", param_name = "sd",
              fun = function(data, params) {
                if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                Rn <- data$rec_values
                Ln <- data$rec_times
                T <- data$time

                if( all(c("sd") %in% names(params)) == FALSE ) stop("parameter sd should be present.")

                mean = 0
                sd <- 1/as.numeric(params$sd)

                m = length(Rn)
                return(sd^2/(2*sqrt(m)))
              }
)
## LDM, Xt --------------

  ##frechet
register_var( model = "LDM", obs_type = "all", dist = "frechet", param_name = "theta",
              fun = function(data, params) {
                if(!is.numeric(data)) stop("data should be a numerical vector")
                n = length(data)

                if( all(c("theta", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters theta, shape, scale should be present.")

                theta <- as.numeric(params$theta)
                shape <- as.numeric(params$shape)
                scale <- 1/as.numeric(params$scale)

                y = data-theta*(1:n) ## frechet iid
                s1 = scale^(-shape)*shape*(-shape-1) * sum((1:n)^2 * y^(-shape-2))
                s2 = (shape+1)*sum(((1:n)/y)^2)
                fisher= s1+s2
                return(-1/fisher)

              }
)

register_var( model = "LDM", obs_type = "all", dist = "frechet", param_name = "scale",
              fun = function(data, params) {
                if(!is.numeric(data)) stop("data should be a numerical vector")
                n = length(data)

                if( all(c("theta", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters theta, shape, scale should be present.")

                theta <- as.numeric(params$theta)
                shape <- as.numeric(params$shape)
                scale <- 1/as.numeric(params$scale)

                y = data-theta*(1:n) ## frechet iid
                s1 = -(shape+1)*shape*scale^(-shape-2)*sum(y^-shape)
                s2=shape*n/scale^2
                fisher= s1+s2
                return(-1/fisher)

              }
)

register_var( model = "LDM", obs_type = "all", dist = "frechet", param_name = "shape",
              fun = function(data, params) {
                if(!is.numeric(data)) stop("data should be a numerical vector")
                n = length(data)

                if( all(c("theta", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters theta, shape, scale should be present.")

                theta <- as.numeric(params$theta)
                shape <- as.numeric(params$shape)
                scale <- 1/as.numeric(params$scale)

                y = data-theta*(1:n) ## frechet iid
                p1 = y^(-shape)*(log(scale*y))^2
                s1 = -scale^(-shape)*sum(p1)
                s2 = -n/(shape^2)
                fisher =  s1+s2
                return(-1/fisher)

              }
)

## LDM, Rn --------------

  ##frechet
register_var( model = "LDM", obs_type = "records", dist = "frechet", param_name = "theta",
                   fun = function(data, params) {
                     if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                     Rn <- data$rec_values
                     Ln <- data$rec_times
                     T <- data$time

                     if( all(c("theta", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters theta, shape, scale should be present.")

                     theta <- as.numeric(params$theta)
                     shape <- as.numeric(params$shape)
                     scale <- 1/as.numeric(params$scale)
                     m = length(Rn)

                     x= Rn - theta * Ln

                     ## 1) theta Variance
                     s1 =  (shape + 1) * sum(Ln^2 / x^2)
                     s2 = sum( Ln^2 * x^(-shape - 2))
                     s3a=0
                     for (i in 1:(m - 1)) {
                       # Inner summation from (l[i] + 1) to (l[i+1] - 1)
                       for (t in (Ln[i] + 1):(Ln[i + 1] - 1)) {
                         if( ((Ln[i]+1) <=(Ln[i+1]-1) ) == TRUE ){
                           s3a <- s3a + t^2 * (Rn[i] - theta * t)^(-shape - 2) }}
                     }

                        ## case where last record is not last observation
                     s3b=0
                     if( (Ln[m] < T) == TRUE) {
                       for (i in (Ln[m]+1):T) {
                         s3b[i] <- (i^2) * ((Rn[m] - theta * i)^(-shape - 2))
                       } }
                     s3=s3a+sum(na.omit(s3b))
                     s4 = -shape * (shape + 1) * scale^(-shape) * (s3+s2)

                     var_theta = -1 / (s1 + s4)

                   return(var_theta)
                   }
)

register_var( model = "LDM", obs_type = "records", dist = "frechet", param_name = "shape",
              fun = function(data, params) {
                if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                Rn <- data$rec_values
                Ln <- data$rec_times
                T <- data$time

                if( all(c("theta", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters theta, shape, scale should be present.")

                theta <- as.numeric(params$theta)
                shape <- as.numeric(params$shape)
                scale <- 1/as.numeric(params$scale)
                m = length(Rn)

                x= Rn - theta * Ln

                ## variance of shape
                s1_h = - m / (shape^2)

                s2_h= sum((log(scale*x))^2 * x^(-shape))

                s3a_h=0
                for (i in 1:(m - 1)) {
                  # Inner summation from (l[i] + 1) to (l[i+1] - 1)
                  for (t in (Ln[i] + 1):(Ln[i + 1] - 1)) {
                    if( ((Ln[i]+1) <=(Ln[i+1]-1) ) == TRUE ){
                      term <- scale * (Rn[i] - theta * t)
                      s3a_h <- s3a_h + (term^(-shape)) * (log(term)^2)  }}
                }

                ## case where last record is not last observation
                s3b_h=0
                if( (Ln[m] < T) == TRUE) {
                  for (i in (Ln[m]+1):T) {
                    term <- scale * (Rn[m] - theta * t)
                    s3b_h[i] <-(term^(-shape)) * (log(term)^2)
                  } }

                s3_h=s3a_h+sum(na.omit(s3b_h))

                s4_h <- - (scale^(-shape)) * (s2_h + s3_h)

                var_shape <- -1 / (s1_h+s4_h)

                return(var_shape )
              }
)

register_var( model = "LDM", obs_type = "records", dist = "frechet", param_name = "scale",
              fun = function(data, params) {
                if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                Rn <- data$rec_values
                Ln <- data$rec_times
                T <- data$time

                if( all(c("theta", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters theta, shape, scale should be present.")
                theta <- as.numeric(params$theta)
                shape <- as.numeric(params$shape)
                scale <- 1/as.numeric(params$scale)
                m = length(Rn)

                x= Rn - theta * Ln


                ## 1/Scale variance

                s1_s = m * shape / (scale^2)
                s2_s = sum(x^(-shape))

                s3a_s=0
                for (i in 1:(m - 1)) {
                  # Inner summation from (l[i] + 1) to (l[i+1] - 1)
                  for (t in (Ln[i] + 1):(Ln[i + 1] - 1)) {
                    if( ((Ln[i]+1) <=(Ln[i+1]-1) ) == TRUE ){
                      s3a_s <- s3a_s + (Rn[i] - theta * t)^(-shape ) }}
                }

                ## case where last record is not last observation
                s3b_s=0
                if( (Ln[m] < T) == TRUE) {
                  for (i in (Ln[m]+1):T) {
                    s3b_s[i] <- ((Rn[m] - theta * i)^(-shape))
                  } }

                s3_s=s3a_s+sum(na.omit(s3b_s))

                s4_s <- -shape * (shape + 1) * (scale^(-shape - 2)) * (s2_s + s3_s)

                var_scale <- -1 / (s1_s + s4_s)

                return(var_scale)
              }
)

## YNM, Xt -----------------

  ## frechet
register_var( model = "YNM", obs_type = "all", dist = "frechet", param_name = "gamma",
              fun = function(data, params) {
                if(!is.numeric(data)) stop("data should be a numerical vector")
                n = length(data)
                if( all(c("gamma", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters gamma, shape, scale should be present.")

                gamma <- as.numeric(params$gamma)
                shape <- as.numeric(params$shape)
                scale <- 1/as.numeric(params$scale)

                sum_val <- sum((1:n) * (0:(n-1)) * (gamma^((1:n)-2)) * (scale * data)^(-shape))
                fisher <- (-n * (n + 1) / (2 * gamma^2)) - sum_val
                return(-1 / fisher)

              }
)

register_var( model = "YNM", obs_type = "all", dist = "frechet", param_name = "scale",
              fun = function(data, params) {
                if(!is.numeric(data)) stop("data should be a numerical vector")
                n = length(data)
                if( all(c("gamma", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters gamma, shape, scale should be present.")

                gamma <- as.numeric(params$gamma)
                shape <- as.numeric(params$shape)
                scale <- 1/as.numeric(params$scale)

                sum_val = sum( (gamma^(1:n))*(data^-shape)  )
                fisher =  (shape*n/(scale^2)) -  shape*(shape+1) * scale^(-shape-2)*  sum_val
                return(-1 / fisher)

              }
)

register_var( model = "YNM", obs_type = "all", dist = "frechet", param_name = "shape",
              fun = function(data, params) {
                if(!is.numeric(data)) stop("data should be a numerical vector")
                n = length(data)
                if( all(c("gamma", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters gamma, shape, scale should be present.")

                gamma <- as.numeric(params$gamma)
                shape <- as.numeric(params$shape)
                scale <- 1/as.numeric(params$scale)

                sum_val = sum(  (gamma^(1:n))  *  (log(scale *data))^2    * (scale*data)^(-shape)  )
                fisher =  (-n / (shape^2)) - sum_val
                return(-1 / fisher)

              }
)


## YNM, Rn -----------------

## frechet
register_var( model = "YNM", obs_type = "records", dist = "frechet", param_name = "gamma",
              fun = function(data, params) {
                if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                Rn <- data$rec_values
                Ln <- data$rec_times
                T <- data$time

                if( all(c("gamma", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters gamma, shape, scale should be present.")

                v_Rn_Exp_YNM_gamma <- function(a, b, x) {
                  # Terms involving b
                  b_term <- x^b * ((b^2 - 3*b + 2) * x^2 + (4*b - 2*b^2) * x + b^2 - b)

                  # Terms involving a
                  a_term <- x^a * ((-a^2 + 3*a - 2) * x^2 + (2*a^2 - 4*a) * x - a^2 + a)

                  # Denominator
                  denominator <- (x - 1)^3 * x^2

                  # Final result
                  result <- (b_term + a_term) / denominator

                  return(result)
                }

                gamma <- as.numeric(params$gamma)
                shape <- as.numeric(params$shape)
                scale <- 1/as.numeric(params$scale)

                m = length(Rn)

                term1 = -sum(Ln)/gamma^2

                term2 = numeric(m)
                for(i in 1:(m-1)){
                  term2[i]=(scale*Rn[i])^(-shape) * v_Rn_Exp_YNM_gamma(a=Ln[i],b=Ln[i+1],x=gamma)
                }
                if(Ln[m]<T){
                  term2[m] = (scale*Rn[m])^(-shape)* v_Rn_Exp_YNM_gamma(a=Ln[m],b=(T+1),x=gamma)}

                fisher= term1-sum(term2)
                return(-1/fisher)
              }
)

register_var( model = "YNM", obs_type = "records", dist = "frechet", param_name = "scale",
              fun = function(data, params) {
                if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                Rn <- data$rec_values
                Ln <- data$rec_times
                T <- data$time

                if( all(c("gamma", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters gamma, shape, scale should be present.")

                gamma <- as.numeric(params$gamma)
                shape <- as.numeric(params$shape)
                scale <- 1/as.numeric(params$scale)

                m = length(Rn)

                term1=shape*m/scale^2

                term2 = numeric(m)
                for(i in 1:(m-1)){
                  term2[i] = (gamma^Ln[i]-gamma^Ln[i+1])*Rn[i]^(-shape)/(1-gamma)
                }

                if(Ln[m]<T){
                  term2[m] = (gamma^Ln[m]-gamma^(T+1))*Rn[m]^(-shape)/(1-gamma) }

                term2=term2*(shape^2+shape) * scale^(-shape-2)

                fisher=term1-sum(term2)

                return(-1/fisher)
              }
)

register_var( model = "YNM", obs_type = "records", dist = "frechet", param_name = "shape",
              fun = function(data, params) {
                if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                Rn <- data$rec_values
                Ln <- data$rec_times
                T <- data$time

                if( all(c("gamma", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters gamma, shape, scale should be present.")

                gamma <- as.numeric(params$gamma)
                shape <- as.numeric(params$shape)
                scale <- 1/as.numeric(params$scale)

                m = length(Rn)

                term1= -m/shape^2

                term2=numeric(m)
                for(i in 1:(m-1)){
                  term2[i] = (gamma^Ln[i]-gamma^Ln[i+1])/(1-gamma)  * (log(scale*Rn[i]))^2 * (scale*Rn[i])^(-shape)
                }

                if(Ln[m]<T){
                  term2[m] = (gamma^Ln[m]-gamma^(T+1))/(1-gamma) * (log(scale*Rn[m]))^2 * (scale*Rn[m])^(-shape) }

                fisher = term1 - sum(term2)

                return(-1/fisher)
              }
)


