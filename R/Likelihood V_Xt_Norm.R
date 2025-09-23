## Variance of the variance estimator sigma2
V_Xt_Norm_DTRW_a = function(x,params){
  a= params
  T=length(x)
  XX=diff(x)

  #s1 = T/(2*a^2)

  #s2 = - sum( (XX)^2) / (a^3)

  #var_a = -1 / (s1+s2)
  var_a = 2*a^2/(T)  ## 2*sigma^2/T
  return(var_a)
}

