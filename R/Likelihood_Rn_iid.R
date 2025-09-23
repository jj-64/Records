# Define the function in R

Likelihood_Rn_Frechet_iid= function(R,L,T,params){  ## Frechet parametrs 1/scale, shape
  m = length(R) #m= rec_counts(y)  ## number of records

  ## sum of log(f_rn) expect first record
  s1 = sum(log( VGAM::dfrechet(R[-1],loc=0, scale=1/params[1], shape=params[2])  ) )

  ## sum of log(F_rn) expect last record
  s2b=0
  for(i in 1:(m-1)){
    s2b[i] = (L[i+1]-L[i]-1 )*log( VGAM::pfrechet(R[i],loc=0, scale=1/params[1], shape=params[2]) )
  }
  s2=sum(s2b)

  ## Case where last record is not last observation
  s3=0
  if(L[m]<T){
    s3=(T-L[m])*log(VGAM::pfrechet(R[m],loc=0, scale=1/params[1], shape=params[2]) )
  }

  return(s1+s2+s3)
}

## detailed
Likelihood_Rn_Frechet2_iid = function(R,L,T,params){  ## Frechet parametrs 1/scale, shape
  m = length(R) #m= rec_counts(y)  ## number of records

  ## sum of log(f_rn) expect first record
  s1a = (m-1) * log(params[1]*params[2])
  s1b = -sum((params[1]*R[-1])^(-params[2]))
  s1c = (-1-params[2])*sum(log(params[1]*R[-1]))
  s1 = s1a+s1b+s1c

  ## sum of log(F_rn) expect last record
  s2b=0
  for(i in 1:(m-1)){
    s2b[i] = (L[i+1]-L[i]-1 )* (params[1]*R[i] )^(-params[2])
  }
  s2=-sum(s2b)

  ## Case where last record is not last observation
  s3=0
  if(L[m]<T){
    s3=-(T-L[m])*(params[1]*R[m] )^(-params[2])
  }

  return(s1+s2+s3)
}

Likelihood_Rn_Norm_iid = function(R,L,T,params){  ## mean, sd
  m = length(R) #m= rec_counts(y)

  ## sum of log(f_rn)
  s1 = sum(log( dnorm(R[-1],mean = params[1], sd=params[2])  ) )


  s2b=0
  for(i in 1:(m-1)){
    s2b[i] = (L[i+1]-L[i]-1 )*log( pnorm(R[i],mean = params[1], sd=params[2]) )
  }
  s2=sum(s2b)

  s3=0
  if(L[m]<T){
    s3=(T-L[m])*log( pnorm(R[m],mean = params[1], sd=params[2]) )
  }

  return(s1+s2+s3)
}

Likelihood_Rn_Exp_iid = function(R,L,T,params){  ## mean, sd
  m = length(R) #m= rec_counts(y)

  ## sum of log(f_rn)
  s1 = sum(log( dexp(R[-1],rate=params)  ) )


  s2b=0
  for(i in 1:(m-1)){
    s2b[i] = (L[i+1]-L[i]-1 )*log( pexp(R[i],rate=params) )
  }
  s2=sum(s2b)

  s3=0
  if(L[m]<T){
    s3=(T-L[m])*log( pexp(R[m],rate=params) )
  }

  return(s1+s2+s3)
}

Likelihood_Rn_Gumbel_iid = function(R,L,T,params){  ## mean, sd
  m = length(R) #m= rec_counts(y)

  ## sum of log(f_rn)
  s1 = sum(log( VGAM::dgumbel(R[-1],loc=params[1], scale=params[2])  ) )


  s2b=0
  for(i in 1:(m-1)){
    s2b[i] = (L[i+1]-L[i]-1 )*log(VGAM::pgumbel(R[i],loc=params[1], scale=params[2]) )
  }
  s2=sum(s2b)

  s3=0
  if(L[m]<T){
    s3=(T-L[m])*log(VGAM::pgumbel(R[m],loc=params[1], scale=params[2]) )
  }

  return(s1+s2+s3)
}
