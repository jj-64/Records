############################# Variance estimators - Rn - Norm - DTRW ############################
## Variance of mean A

V_Rn_Norm_DTRW_shape = function(R,L,T,params){  ## A=mean, a= sigma, sd
  A = params[1]
  a= params[2]
  m=length(R)
  t = diff(L)

  s1= -1/(a^2) * sum(t)

  s2b=0
  # for(j in 1:(T-L[m])){
  #   k= -sqrt(j)/a
  #   f=dnorm(0,mean=j*A,sd=sqrt(j)*a)
  #   FF=pnorm(0,mean=j*A,sd=sqrt(j)*a)
  #   r=f/FF
  #
  #     s2b[j]=k*r *(-j*A/a^2) - (k*r)^2
  # }
  #
  s2 = sum(s2b)


  var_A <- -1 / (s1 + s2)

  return(var_A)
}

## Variance of variance alpha or a  ## revisit
V_Rn_Norm_DTRW_scale2 = function(R,L,T,params){
  A = params[1]
  a= params[2]


  m=length(R)
  RR=diff(R)
  t = diff(L)

  #s1 = (m-1)/a^2

  #s2 = -sum(3*(RR-t*A)^2/(t*a^4))

  s1 = (m-1)/(2*a^2)

  s2 = - sum( (RR)^2/t) / (a^3)

  s3a=0
  for(j in 1:(T-L[m])){

    f=dnorm(0,mean=j*A,sd=sqrt(j*a))
    FF=pnorm(0,mean=j*A,sd=sqrt(j*a))
    r=f/FF

    s3a[j] = (sqrt(j)*A)*(r/a^3) * (-1+ ((-j*A)^2/(j*a)^2) + (r*(-j*A)/(sqrt(j)*a)) )
  }

  s3 = sum(s3a)


  var_a = -1 / (s1+s2 + s3)

  return(var_a)
}

## Variance of the variance estimator sigma2
V_Rn_Norm_DTRW_scale = function(R,L,T,params){
  a= params[2]
  m=length(R)
  #RR=diff(R)
  #t = diff(L)

  #s1 = (m-1)/(2*a^2)
  #s2 = - sum( (RR)^2/t) / (a^3)
  #var_a = -1 / (s1+s2)
  #var_a = 2*a^2/(m)
  var_a = a^2/(2*sqrt(m))
  return(var_a)
}


