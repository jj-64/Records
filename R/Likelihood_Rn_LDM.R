## Likelihood function LDM
Likelihood_Rn_Frechet_LDM = function(R,L,T,params){
  m = length(R) #m= rec_counts(y)  ## number of records
  x = R - params[1]* L ## Rn-theta*ln

  # Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or params[3] <= 0, we return -Inf)
  if (any(is.na(params) | is.nan(params) | params <= 0)) { ##
    return(-Inf)  # Invalid parameters, return a large negative value
  }

  # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
  s1 = sum(log(VGAM::dfrechet(x, loc=0, scale=1/params[2],shape=params[3])))
  if (is.nan(s1) || !is.finite(s1)) return(-Inf)

  s2b=0
  for( i in 1:(m-1)  ){

    if((L[i]+1 <= L[i+1]-1) == TRUE){ ## we have non-records in between
      s2b[i]=sum(log(VGAM::pfrechet( R[i]-params[1] * (L[i]+1):(L[i+1]-1), loc=0, scale=1/params[2],shape=params[3] )  ))
    }
  }
  s2 = sum(na.omit(s2b))
  if (is.nan(s2) || !is.finite(s2)) return(-Inf)

  ## case where last record is not last observation
  s3=0
  if( (L[m] < T) == TRUE  ) {
    s3=sum(log(VGAM::pfrechet( R[m]-params[1] * (L[m]+1):T, loc=0, scale=1/params[2],shape=params[3] )  ))
  }
  if (is.nan(s3) || !is.finite(s3)) return(-Inf)

  vrais= s1+s2+s3

  return(vrais)
}

################### Gumbel ##################
## Likelihood function LDM
Likelihood_Rn_Gumbel_LDM = function(R,L,T,params){
  m = length(R) #m= rec_counts(y)  ## number of records
  x = R - params[1]* L ## Rn-theta*ln

  ## pdf
  pdf=function(x,par) {VGAM::dgumbel(x=x, loc=par[1], scale=par[2])}
  ## cdf
  cdf=function(x,par) {VGAM::pgumbel(q=x, loc=par[1], scale=par[2])}

  # Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or params[3] <= 0, we return -Inf)
  # if (any(params <= 0)) {
  #   return(-Inf)  ## Invalid parameters, return a large negative value
  # }

  # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
  s1 = sum(log(pdf(x=x, par=params[2:3])))
  if (is.nan(s1) || !is.finite(s1)) return(-Inf)

  s2b=0
  for( i in 1:(m-1)  ){

    if((L[i]+1 <= L[i+1]-1) == TRUE){ ## we have non-records in between
      s2b[i]=sum(log(cdf( R[i]-params[1] * (L[i]+1):(L[i+1]-1), par=params[2:3] )  ))
    }
  }

  s2 = sum(na.omit(s2b))
  if (is.nan(s2) || !is.finite(s2)) return(-Inf)

  ## case where last record is not last observation
  s3=0
  if( (L[m] < T) == TRUE  ) {
    s3=sum(log(cdf( R[m]-params[1] * (L[m]+1):T, par=params[2:3] )  ))
  }
  if (is.nan(s3) || !is.finite(s3)) return(-Inf)

  vrais= s1+s2+s3

  return(vrais)
}

################### Weibull ##################
Likelihood_Rn_Weibull_LDM = function(R,L,T,params){
  m = length(R) #m= rec_counts(y)  ## number of records
  x = R - params[1]* L ## Rn-theta*ln

  ## pdf
  pdf=function(x,par) {dweibull(x=x, shape=par[1],scale=par[2])}
  ## cdf
  cdf=function(x,par) {pweibull(q=x,  shape=par[1],scale=par[2])}

  # Check for invalid parameter values
  # if (any(params <= 0)) {
  #                   return(-Inf)  # Invalid parameters, return a large negative value
  #                             }

  # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
  s1 = sum(log(pdf(x=x, par=params[2:3])))
  if (is.nan(s1) || !is.finite(s1)) return(-Inf)

  s2b=0
  for( i in 1:(m-1)  ){

    if(  (L[i]+1 <= L[i+1]-1) == TRUE){ ## we have non-records in between
      s2b[i]=sum(log(cdf( R[i]-params[1] * (L[i]+1):(L[i+1]-1), par=params[2:3] )  ))
    }
  }

  s2 = sum(na.omit(s2b))
  if (is.nan(s2) || !is.finite(s2)) return(-Inf)

  ## case where last record is not last observation
  s3=0
  if( (L[m] < T) == TRUE  ) {
    s3=sum(log(cdf( x=R[m]-params[1] * (L[m]+1):T, par=params[2:3] )  ))
  }
  if (is.nan(s3) || !is.finite(s3)) return(-Inf)

  vrais= s1+s2+s3

  return(vrais)
}

################### Exp #########################

Likelihood_Rn_Exp_LDM = function(R,L,T,params){  ## Theta and scale= 1/rate
  m = length(R) #m= rec_counts(y)  ## number of records
  x = R - params[1]* L ## Rn-theta*ln

  ## pdf
  pdf=function(x,par) {dexp(x=x, rate=1/par[1])}  ## rate is 1/scale in weibull
  ## cdf
  cdf=function(x,par) {pexp(q=x, rate=1/par[1])}

  # Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or params[3] <= 0, we return -Inf)
  # if (any(params <= 0)) {
  #   return(-Inf)  # Invalid parameters, return a large negative value
  # }

  # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
  s1 = sum(log(pdf(x=x, par=params[2])))
  if (is.nan(s1) || !is.finite(s1)) return(-Inf)

  s2b=0
  for( i in 1:(m-1)  ){

    if(  (L[i]+1 <= (L[i+1]-1)) == TRUE){ ## we have non-records in between
      s2b[i]=sum(log(cdf( R[i]-params[1] * (L[i]+1):(L[i+1]-1), par=params[2] )  ))
    }
  }

  s2 = sum(na.omit(s2b))
  if (is.nan(s2) || !is.finite(s2)) return(-Inf)

  ## case where last record is not last observation
  s3=0
  if( (L[m] < T) == TRUE  ) {
    s3=sum(log(cdf( x=R[m]-params[1] * (L[m]+1):T, par=params[2] )  ))
  }
  if (is.nan(s3) || !is.finite(s3)) return(-Inf)

  vrais= s1+s2+s3

  return(vrais)
}





