##################Variance estimators - Rn - Frechet - YNM ##########################
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

## Variance of gamma
V_Rn_Frechet_YNM_gamma =  function(R,L,T,params){  ## vector of Rn, Ln, values of Time, gamma, A and a of Frechet parameters
  m = length(R) ## NT = m
  g = params[1]
  A = params[2]
  a = params[3]

  term1 = -sum(L)/g^2

  term2 = numeric(m)
  for(i in 1:(m-1)){
    term2[i]=(A*R[i])^(-a) * v_Rn_Exp_YNM_gamma(a=L[i],b=L[i+1],x=g)
  }
  if(L[m]<T){
    term2[m] = (A*R[m])^(-a)* v_Rn_Exp_YNM_gamma(a=L[m],b=(T+1),x=g)}

  fisher= term1-sum(term2)
  return(-1/fisher)
}

## Variance of 1/scale A
V_Rn_Frechet_YNM_scale = function(R,L,T,params){
  g = params[1]
  A = params[2]
  a = params[3]

  m=length(R)
  term1=a*m/A^2

  term2 = numeric(m)
  for(i in 1:(m-1)){
    term2[i] = (g^L[i]-g^L[i+1])*R[i]^(-a)/(1-g)
  }

  if(L[m]<T){
    term2[m] = (g^L[m]-g^(T+1))*R[m]^(-a)/(1-g) }

  term2=term2*(a^2+a) * A^(-a-2)

  fisher=term1-sum(term2)

  return(-1/fisher)
}

## Variance of shape a
V_Rn_Frechet_YNM_shape = function(R,L,T,params){
  m=length(R)
  g = params[1]
  A = params[2]
  a = params[3]

  term1= -m/a^2

  term2=numeric(m)
  for(i in 1:(m-1)){
    term2[i] = (g^L[i]-g^L[i+1])/(1-g)  * (log(A*R[i]))^2 * (A*R[i])^(-a)
  }

  if(L[m]<T){
    term2[m] = (g^L[m]-g^(T+1))/(1-g) * (log(A*R[m]))^2 * (A*R[m])^(-a) }

  fisher = term1 - sum(term2)

  return(-1/fisher)

}
##################Variance estimators - Rn - Frechet - LDM ##########################

## Variance of theta
V_Rn_Frechet_LDM_theta =  function(R,L,T, params) { ## y vector and theta, A and alpha
  theta = params[1]
  A = params[2]
  a = params[3]

  m = length(R)  ## number of records

  x= R - theta * L

  s1 =  (a + 1) * sum(L^2 / x^2)

  s2 = sum( L^2 * x^(-a - 2))

  s3a=0
  for (i in 1:(m - 1)) {
    # Inner summation from (l[i] + 1) to (l[i+1] - 1)
    for (t in (L[i] + 1):(L[i + 1] - 1)) {
      if( ((L[i]+1) <=(L[i+1]-1) ) == TRUE ){
        s3a <- s3a + t^2 * (R[i] - theta * t)^(-a - 2) }}
  }

  ## case where last record is not last observation
  s3b=0
  if( (L[m] < T) == TRUE) {
    for (i in (L[m]+1):T) {
      s3b[i] <- (i^2) * ((R[m] - theta * i)^(-a - 2))
    } }

  s3=s3a+sum(na.omit(s3b))

  s4 = -a * (a + 1) * A^(-a) * (s3+s2)

  var_theta = -1 / (s1 + s4)

  return(var_theta)
}

## Variance of 1/scale A

V_Rn_Frechet_LDM_scale = function(R,L,T,params){
  theta = params[1]
  A = params[2]
  a = params[3]

  m=length(R)

  s1 = m * a / (A^2)

  x=R - theta * L
  s2 = sum(x^(-a))

  s3a=0
  for (i in 1:(m - 1)) {
    # Inner summation from (l[i] + 1) to (l[i+1] - 1)
    for (t in (L[i] + 1):(L[i + 1] - 1)) {
      if( ((L[i]+1) <=(L[i+1]-1) ) == TRUE ){
        s3a <- s3a + (R[i] - theta * t)^(-a ) }}
  }

  ## case where last record is not last observation
  s3b=0
  if( (L[m] < T) == TRUE) {
    for (i in (L[m]+1):T) {
      s3b[i] <- ((R[m] - theta * i)^(-a))
    } }

  s3=s3a+sum(na.omit(s3b))

  s4 <- -a * (a + 1) * (A^(-a - 2)) * (s2 + s3)

  var_A <- -1 / (s1 + s4)

  return(var_A)
}

## Variance of shape a
V_Rn_Frechet_LDM_shape = function(R,L,T,params){
  m=length(R)
  theta = params[1]
  A = params[2]
  a = params[3]

  x=R-theta*L

  s1 = - m / (a^2)

  s2= sum((log(A*x))^2 * x^(-a))

  s3a=0
  for (i in 1:(m - 1)) {
    # Inner summation from (l[i] + 1) to (l[i+1] - 1)
    for (t in (L[i] + 1):(L[i + 1] - 1)) {
      if( ((L[i]+1) <=(L[i+1]-1) ) == TRUE ){
        term <- A * (R[i] - theta * t)
        s3a <- s3a + (term^(-a)) * (log(term)^2)  }}
  }

  ## case where last record is not last observation
  s3b=0
  if( (L[m] < T) == TRUE) {
    for (i in (L[m]+1):T) {
      term <- A * (R[m] - theta * t)
      s3b[i] <-(term^(-a)) * (log(term)^2)
    } }

  s3=s3a+sum(na.omit(s3b))

  s4 <- - (A^(-a)) * (s2 + s3)

  var_a <- -1 / (s1+s4)

  return(var_a)
}

##################Variance estimators - Rn - Frechet - iid ##########################


## Variance of 1/scale A
#' @title Variance Estimator for Scale Parameter (A) in Fréchet IID Model
#' @description Computes the variance of the scale parameter \eqn{ A } in an independent and identically distributed (IID) Fréchet model for record values.
#'
#' @details The variance is computed as:
#' \deqn{ V(A) = - I(A)^{-1}}
#' where:
#' \deqn{I(A) = s_1 + s_2 + s_3 + s_4 + s_5}
#' \deqn{ s_1 = -\frac{(N_T-1)}{A^2} }
#' \deqn{ s_2 = -\frac{(\alpha^2 + a)}{A^2} \sum_{n=2}^{N_T} (x_{L_n} A)^{-\alpha} }
#' \deqn{ s_3 = \frac{(\alpha+1)(N_T-1)}{A^2} }
#' \deqn{ s_4 = -\frac{(\alpha^2 + \alpha)}{A^2} \sum_{n=1}^{N_T-1} (L_{n+1} - L_n - 1) (A x_{L_n})^{-\alpha} }
#' \deqn{ s_5 = -\frac{(\alpha(\alpha+1))(T-L_{N_T})}{A^2} (A x_{L_{N_T}})^{-a}, \quad \text{if } L_{N_T} < T }
#'
#' @param R Vector of record values.
#' @param L Vector of record times.
#' @param T Total observation time.
#' @param params Numeric vector of two values: Scale parameter (\eqn{1/A} > 0) and Shape parameter (\eqn{\alpha} > 0).
#' @return The variance estimate of \eqn{ A }.
#' @examples
#' Xt = VGAM::rfrechet(n=25,location=0, scale=1, shape=2)
#' R=record_values(Xt)
#' L = record_times(Xt)
#' Xt
#' [1] 1.8981377 1.0773402 1.3607312 1.5417757 0.8582019 1.0779739 1.5744627
#' [8] 2.2239265 2.0618936 0.8036526 0.4404408 0.5775581 1.2024731 1.2722281
#' [15] 1.6804531 4.1700017 1.2460309 0.7485379 2.7474130 3.9600508 2.5000833
#' [22] 2.7688841 2.8051514 1.0738091 0.7334323
#' R
#' [1] 1.898138 2.223927 4.170002
#' L
#' [1]  1  8 16
#' V_Rn_Frechet_iid_scale(R=R,L=L,T=25,A=1,a=2)
#' [1] 0.05222625
#' V_Rn_Frechet_iid_shape(R=R,L=L,T=25,A=1,a=2)
#' [1] 0.3695643
V_Rn_Frechet_iid_scale = function(R,L,T,params){
  A = params[1]
  a = params[2]

  m=length(R)

  s1=-(m-1)/(A^2)

  s2 = -((a^2+a)/A^2) * sum((R[-1]*A)^(-a))

  s3 = (a+1)*(m-1)/A^2

  s4a=0
  for (i in 1:(m - 1)) {
    s4a[i] = (L[i+1]-L[i]-1)*(A*R[i])^(-a)
  }
  s4 = -(a^2+a)*sum(s4a)/A^2

  ## case where last record is not last observation
  s5=0
  if( (L[m] < T) == TRUE) {
    s5 = -(a*(a+1))*(T-L[m])*(A*R[m])^(-a)/A^2
  }

  var_A <- -1 / (s1 + s2 + s3 + s4 + s5)

  return(var_A)
}

# Variance of shape a
#' @title Variance Estimator for Shape Parameter (a) in Fréchet IID Model
#' @description Computes the variance of the shape parameter \eqn{ \alpha } in an IID Fréchet model for record values.
#'
#' @details The variance is computed as:
#' \deqn{ V(\alpha) = - I(\alpha)^{-1}}
#' where:
#' \deqn{I(\alpha) = s_1 + s_2 + s_3 + s_4 + s_5}
#' \deqn{ s_1 = -\frac{(N_T-1)}{\alpha^2} }
#' \deqn{ s_2 = - \sum_{n=2}^{N_T} \frac{(\log(A x_{L_n}))^2}{(A x_{L_n})^\alpha} }
#' \deqn{ s_3 = - \sum_{n=1}^{N_T-1} (L_{n+1} - L_n - 1) (A x_{L_n})^{-\alpha} (\log(A x_{L_n}))^2 }
#' \deqn{ s_4 = - (T - L_{N_T}) (A x_{L_{N_T}})^{-\alpha} (\log(A x_{L_{N_T}}))^2, \quad \text{if } L{N_T} < T }
#'
#' @param R Numeric vector of record values.
#' @param L Numeric vector of record times.
#' @param T Numeric value representing total observation time.
#' @param params Numeric vector of two values: Scale parameter (\eqn{1/A} > 0) and Shape parameter (\eqn{\alpha} > 0).
#' @return The variance estimate of \eqn{ \alpha }.
V_Rn_Frechet_iid_shape = function(R, L, T, params) {
  A = params[1]
  a = params[2]

  m = length(R)

  s1 = -(m - 1) / a^2
  s2 = -sum((log(A * R[-1]))^2 / (A * R[-1])^a)

  s3a = 0
  for (i in 1:(m - 1)) {
    s3a[i] = (L[i + 1] - L[i] - 1) * (A * R[i])^(-a) * (log(A * R[i]))^2
  }
  s3 = -sum(s3a)

  ## Case where last record is not the last observation
  s4 = 0
  if (L[m] < T) {
    s4 = -(T - L[m]) * (A * R[m])^(-a) * (log(A * R[m]))^2
  }

  var_a = -1 / (s1 + s2 + s3 + s4)

  return(var_a)
}

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



