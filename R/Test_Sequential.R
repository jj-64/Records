## Score 0: nothing
## Score 1: Classical
## Score 2: DTRW
## Score 3: LDM
## Score 4: YANG

######## 1- CLYD ########
DT_CLYD = function(X,p=0.05){
  score=0

  ##Classical
  dec=Test_iid(X, p=p)$dec
  if(dec!="NO"){     return(score=1)}

  ## LDM
  if(dec=="NO"){ dec=Test_LDM_Regression(X, p=p)$dec }
  if(dec!="NO"){     return(score=3)}

  ## YANG
  if(dec=="NO"){
    gamma_hat = Estim_gamma_indicator(X=X,min=1,max=5)
    part= partition(X)
    if(length(part$j)==1){dec="NO"     }else{     dec=Test_Yang_Pearson(X=X,gamma=gamma_hat,estimated=1,p=p)$dec}
  }
  if(dec!="NO"){    return(score=4)}

  ##DTRW
  if(dec=="NO"){
    dec=Test_DTRW_bonf(X, p=p)$dec
  }
  if(dec!="NO"){     return(score=2)}

  ## No model
  if(dec=="NO"){ return(score)}
}

######## 2 - CLDY ########
DT_CLDY = function(X,p=0.05){
  score=0

  ##Classical
  dec=Test_iid(X, p=p)$dec
  if(dec!="NO"){    return( score=1)}

  ## LDM
  if(dec=="NO"){
    dec=Test_LDM_Regression(X, p=p)$dec }
  if(dec!="NO"){    return( score=3)}

  ##DTRW
  if(dec=="NO"){
    dec=Test_DTRW_bonf(X, p=p)$dec
  }
  if(dec!="NO"){    return( score=2)}


  ## YANG
  if(dec=="NO"){
    gamma_hat = Estim_gamma_indicator(X=X,min=1,max=5)
    part= partition(X)
    if(length(part$j)==1){dec="NO"     }else{    dec=Test_Yang_Pearson(X=X,gamma=gamma_hat,estimated=1,p=p)$dec}
  }
  if(dec!="NO"){    return( score=4)}

  ## No model
  if(dec=="NO"){ return(score)}
}

######## 3 - CYLD ########
DT_CYLD = function(X,p=0.05){
  score=0

  ##Classical
  dec=Test_iid(X, p=p)$dec
  if(dec!="NO"){    return( score=1)}

  ## YANG
  if(dec=="NO"){
    gamma_hat = Estim_gamma_indicator(X=X,min=1,max=5)
    part= partition(X)
    if(length(part$j)==1){dec="NO"     }else{     dec=Test_Yang_Pearson(X=X,gamma=gamma_hat,estimated=1,p=p)$dec}
  }
  if(dec!="NO"){    return( score=4)}

  ## LDM
  if(dec=="NO"){
    dec=Test_LDM_Regression(X, p=p)$dec
  }
  if(dec!="NO"){    return( score=3)}

  ##DTRW
  if(dec=="NO"){
    dec=Test_DTRW_bonf(X, p=p)$dec
  }
  if(dec!="NO"){    return( score=2)}

  ## No model
  if(dec=="NO"){ return(score)}
}

######## 4 - CYDL ########
DT_CYDL = function(X,p=0.05){

  score=0

  ##Classical
  dec=Test_iid(X, p=p)$dec
  if(dec!="NO"){    return( score=1)}

  ## YANG
  if(dec=="NO"){
    gamma_hat = Estim_gamma_indicator(X=X,min=1,max=5)
    part= partition(rec_gaps(X))
    if(length(part$j)==1){dec="NO"     }else{     dec=Test_Yang_Pearson(X=X,gamma=gamma_hat,estimated=1,p=p)$dec}
  }
  if(dec!="NO"){    return( score=4)}

  ##DTRW
  if(dec=="NO"){
    dec=Test_DTRW_bonf(X, p=p)$dec
  }
  if(dec!="NO"){    return( score=2)}

  ## LDM
  if(dec=="NO"){
    dec=Test_LDM_Regression(X, p=p)$dec
  }
  if(dec!="NO"){    return( score=3)}

  ## No model
  if(dec=="NO"){ return(score)}


}

######## 5 - CDYL ########
DT_CDYL = function(X,p=0.05){

  score=0

  ##Classical
  dec=Test_iid(X, p=p)$dec
  if(dec!="NO"){    return( score=1)}

  ##DTRW
  if(dec=="NO"){
    dec=Test_DTRW_bonf(X, p=p)$dec
  }
  if(dec!="NO"){    return( score=2)}

  ## YANG
  if(dec=="NO"){
    gamma_hat = Estim_gamma_indicator(X=X,min=1,max=5)
    part= partition(X)
    if(length(part$j)==1){dec="NO"     }else{     dec=Test_Yang_Pearson(X=X,gamma=gamma_hat,estimated=1,p=p)$dec}
  }
  if(dec!="NO"){    return( score=4)}

  ## LDM
  if(dec=="NO"){
    dec=Test_LDM_Regression(X, p=p)$dec
  }
  if(dec!="NO"){    return( score=3)}

  ## No model
  if(dec=="NO"){ return(score)}
}

######## 6 - CDLY ########
DT_CDLY = function(X,p=0.05){

  score=0

  ##Classical
  dec=Test_iid(X, p=p)$dec
  if(dec!="NO"){    return( score=1)}

  ##DTRW
  if(dec=="NO"){
    dec=Test_DTRW_bonf(X, p=p)$dec
  }
  if(dec!="NO"){    return( score=2)}

  ## LDM
  if(dec=="NO"){
    dec=Test_LDM_Regression(X, p=p)$dec
  }
  if(dec!="NO"){    return( score=3)}

  ## YANG
  if(dec=="NO"){
    gamma_hat = Estim_gamma_indicator(X=X,min=1,max=5)
    part= partition(X)
    if(length(part$j)==1){dec="NO"     }else{     dec=Test_Yang_Pearson(X=X,gamma=gamma_hat,estimated=1,p=p)$dec}
  }
  if(dec!="NO"){    return( score=4)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 7- LCYD ########
DT_LCYD = function(X,p=0.05){

  score=0

  ## LDM
  dec=Test_LDM_Regression(X, p=p)$dec
  if(dec!="NO"){    return( score=3)}

  ##Classical
  if(dec=="NO"){  dec=Test_iid(X, p=p)$dec}
  if(dec!="NO"){    return( score=1)}

  ## YANG
  if(dec=="NO"){
    gamma_hat = Estim_gamma_indicator(X=X,min=1,max=5)
    part= partition(X)
    if(length(part$j)==1){dec="NO"     }else{     dec=Test_Yang_Pearson(X=X,gamma=gamma_hat,estimated=1,p=p)$dec}
  }
  if(dec!="NO"){    return( score=4)}

  ##DTRW
  if(dec=="NO"){
    dec=Test_DTRW_bonf(X, p=p)$dec
  }
  if(dec!="NO"){    return( score=2)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 8- LCDY ########
DT_LCDY = function(X,p=0.05){

  score=0

  ## LDM
  dec=Test_LDM_Regression(X, p=p)$dec
  if(dec!="NO"){    return( score=3)}

  ##Classical
  if(dec=="NO"){  dec=Test_iid(X, p=p)$dec}
  if(dec!="NO"){    return( score=1)}

  ##DTRW
  if(dec=="NO"){
    dec=Test_DTRW_bonf(X, p=p)$dec
  }
  if(dec!="NO"){    return( score=2)}

  ## YANG
  if(dec=="NO"){
    gamma_hat = Estim_gamma_indicator(X=X,min=1,max=5)
    part= partition(X)
    if(length(part$j)==1){dec="NO"     }else{    dec=Test_Yang_Pearson(X=X,gamma=gamma_hat,estimated=1,p=p)$dec}
  }
  if(dec!="NO"){    return( score=4)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 9- LYCD ########
DT_LYCD = function(X,p=0.05){
  score=0

  ## LDM
  dec=Test_LDM_Regression(X, p=p)$dec
  if(dec!="NO"){    return( score=3)}

  ## YANG
  if(dec=="NO"){
    gamma_hat = Estim_gamma_indicator(X=X,min=1,max=5)
    part= partition(X)
    if(length(part$j)==1){dec="NO"
    }else{dec=Test_Yang_Pearson(X=X,gamma=gamma_hat,estimated=1,p=p)$dec}}
  if(dec!="NO"){    return( score=4)}

  ##Classical
  if(dec=="NO"){  dec=Test_iid(X, p=p)$dec}
  if(dec!="NO"){    return( score=1)}

  ##DTRW
  if(dec=="NO"){
    dec=Test_DTRW_bonf(X, p=p)$dec
  }
  if(dec!="NO"){    return( score=2)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 10- LYDC ########
DT_LYDC = function(X,p=0.05){

  score=0

  ## LDM
  dec=Test_LDM_Regression(X, p=p)$dec
  if(dec!="NO"){    return( score=3)}

  ## YANG
  if(dec=="NO"){
    gamma_hat = Estim_gamma_indicator(X=X,min=1,max=5)
    part= partition(X)
    if(length(part$j)==1){dec="NO"
    }else{
      dec=Test_Yang_Pearson(X=X,gamma=gamma_hat,estimated=1,p=p)$dec}
  }
  if(dec!="NO"){    return( score=4)}

  ##DTRW
  if(dec=="NO"){
    dec=Test_DTRW_bonf(X, p=p)$dec
  }
  if(dec!="NO"){    return( score=2)}

  ##Classical
  if(dec=="NO"){  dec=Test_iid(X, p=p)$dec}
  if(dec!="NO"){    return( score=1)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 11- LDYC ########
DT_LDYC = function(X,p=0.05){
  score=0

  ## LDM
  dec=Test_LDM_Regression(X, p=p)$dec
  if(dec!="NO"){    return( score=3)}

  ##DTRW
  if(dec=="NO"){  dec=Test_DTRW_bonf(X, p=p)$dec
  }
  if(dec!="NO"){    return( score=2)}

  ## YANG
  if(dec=="NO"){
    gamma_hat = Estim_gamma_indicator(X=X,min=1,max=5)
    part= partition(X)
    if(length(part$j)==1){dec="NO"
    }else{
      dec=Test_Yang_Pearson(X=X,gamma=gamma_hat,estimated=1,p=p)$dec}
  }
  if(dec!="NO"){    return( score=4)}

  ##Classical
  if(dec=="NO"){ dec=Test_iid(X, p=p)$dec}
  if(dec!="NO"){    return( score=1)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 12- LDCY ########
DT_LDCY = function(X,p=0.05){
  score=0

  ## LDM
  dec=Test_LDM_Regression(X, p=p)$dec
  if(dec!="NO"){    return( score=3)}

  ##DTRW
  if(dec=="NO"){
    dec=Test_DTRW_bonf(X, p=p)$dec
  }
  if(dec!="NO"){    return( score=2)}

  ##Classical
  if(dec=="NO"){dec=Test_iid(X, p=p)$dec}
  if(dec!="NO"){    return( score=1)}

  ## YANG
  if(dec=="NO"){
    gamma_hat = Estim_gamma_indicator(X=X,min=1,max=5)
    part= partition(X)
    if(length(part$j)==1){dec="NO"     }else{     dec=Test_Yang_Pearson(X=X,gamma=gamma_hat,estimated=1,p=p)$dec}
  }
  if(dec!="NO"){    return( score=4)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 13- YCLD ########
DT_YCLD = function(X,p=0.05){
  score=0

  ## YANG
  gamma_hat = Estim_gamma_indicator(X=X,min=1,max=5)
  part= partition(X)
  if(length(part$j)==1){dec="NO"     }else{     dec=Test_Yang_Pearson(X=X,gamma=gamma_hat,estimated=1,p=p)$dec}
  if(dec!="NO"){    return( score=4)}

  ##Classical
  if(dec=="NO"){dec=Test_iid(X, p=p)$dec}
  if(dec!="NO"){    return( score=1)}

  ## LDM
  if(dec=="NO"){dec=Test_LDM_Regression(X, p=p)$dec}
  if(dec!="NO"){    return( score=3)}

  ##DTRW
  if(dec=="NO"){dec=Test_DTRW_bonf(X, p=p)$dec}
  if(dec!="NO"){    return( score=2)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 14- YCDL ########
DT_YCDL = function(X,p=0.05){
  score=0

  ## YANG
  gamma_hat = Estim_gamma_indicator(X=X,min=1,max=5)
  part= partition(X)
  if(length(part$j)==1){dec="NO"     }else{    dec=Test_Yang_Pearson(X=X,gamma=gamma_hat,estimated=1,p=p)$dec}
  if(dec!="NO"){    return( score=4)}

  ##Classical
  if(dec=="NO"){dec=Test_iid(X, p=p)$dec}
  if(dec!="NO"){    return( score=1)}

  ##DTRW
  if(dec=="NO"){dec=Test_DTRW_bonf(X, p=p)$dec }
  if(dec!="NO"){    return( score=2)}

  ## LDM
  if(dec=="NO"){dec=Test_LDM_Regression(X, p=p)$dec}
  if(dec!="NO"){    return( score=3)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 15- YLCD ########
DT_YLCD = function(X,p=0.05){
  score=0

  ## YANG
  gamma_hat = Estim_gamma_indicator(X=X,min=1,max=5)
  part= partition(X)
  if(length(part$j)==1){dec="NO"     }else{     dec=Test_Yang_Pearson(X=X,gamma=gamma_hat,estimated=1,p=p)$dec}
  if(dec!="NO"){    return( score=4)}

  ## LDM
  dec=Test_LDM_Regression(X, p=p)$dec
  if(dec!="NO"){    return( score=3)}

  ##Classical
  if(dec=="NO"){  dec=Test_iid(X, p=p)$dec}
  if(dec!="NO"){    return( score=1)}

  ##DTRW
  if(dec=="NO"){ dec=Test_DTRW_bonf(X, p=p)$dec}
  if(dec!="NO"){    return( score=2)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 16- YLDC ########
DT_YLDC = function(X,p=0.05){
  score=0

  ## YANG
  gamma_hat = Estim_gamma_indicator(X=X,min=1,max=5)
  part= partition(X)
  if(length(part$j)==1){dec="NO"     }else{     dec=Test_Yang_Pearson(X=X,gamma=gamma_hat,estimated=1,p=p)$dec}
  if(dec!="NO"){    return( score=4)}

  ## LDM
  if(dec=="NO"){dec=Test_LDM_Regression(X, p=p)$dec}
  if(dec!="NO"){    return( score=3)}

  ##DTRW
  if(dec=="NO"){
    dec=Test_DTRW_bonf(X, p=p)$dec
  }
  if(dec!="NO"){    return( score=2)}

  ##Classical
  if(dec=="NO"){dec=Test_iid(X, p=p)$dec}
  if(dec!="NO"){    return( score=1)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 17- YDLC ########
DT_YDLC = function(X,p=0.05){
  score=0

  ## YANG
  gamma_hat = Estim_gamma_indicator(X=X,min=1,max=5)
  part= partition(X)
  if(length(part$j)==1){dec="NO"     }else{     dec=Test_Yang_Pearson(X=X,gamma=gamma_hat,estimated=1,p=p)$dec}
  if(dec!="NO"){    return( score=4)}

  ##DTRW
  if(dec=="NO"){dec=Test_DTRW_bonf(X, p=p)$dec}
  if(dec!="NO"){    return( score=2)}

  ## LDM
  if(dec=="NO"){dec=Test_LDM_Regression(X, p=p)$dec}
  if(dec!="NO"){    return( score=3)}

  ##Classical
  if(dec=="NO"){dec=Test_iid(X, p=p)$dec}
  if(dec!="NO"){    return( score=1)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 18- YDCL ########
DT_YDCL = function(X,p=0.05){
  score=0

  ## YANG
  gamma_hat = Estim_gamma_indicator(X=X,min=1,max=5)
  part= partition(X)
  if(length(part$j)==1){dec="NO"     }else{     dec=Test_Yang_Pearson(X=X,gamma=gamma_hat,estimated=1,p=p)$dec}
  if(dec!="NO"){    return( score=4)}

  ##DTRW
  if(dec=="NO"){
    dec=Test_DTRW_bonf(X, p=p)$dec
  }
  if(dec!="NO"){    return( score=2)}

  ##Classical
  if(dec=="NO"){dec=Test_iid(X, p=p)$dec}
  if(dec!="NO"){    return( score=1)}

  ## LDM
  if(dec=="NO"){dec=Test_LDM_Regression(X, p=p)$dec}
  if(dec!="NO"){    return( score=3)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 19- DCLY ########
DT_DCLY = function(X,p=0.05){
  score=0

  ##DTRW
  dec=Test_DTRW_bonf(X, p=p)$dec
  if(dec!="NO"){    return( score=2)}

  ##Classical
  if(dec=="NO"){dec=Test_iid(X, p=p)$dec}
  if(dec!="NO"){    return( score=1)}

  ## LDM
  if(dec=="NO"){dec=Test_LDM_Regression(X, p=p)$dec}
  if(dec!="NO"){    return( score=3)}

  ## YANG
  if(dec=="NO"){
    gamma_hat = Estim_gamma_indicator(X=X,min=1,max=5)
    part= partition(X)
    if(length(part$j)==1){dec="NO"     }else{     dec=Test_Yang_Pearson(X=X,gamma=gamma_hat,estimated=1,p=p)$dec}}

  if(dec!="NO"){    return( score=4)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 20- DCYL ########
DT_DCYL = function(X,p=0.05){
  score=0

  ##DTRW
  dec=Test_DTRW_bonf(X, p=p)$dec
  if(dec!="NO"){    return( score=2)}

  ##Classical
  if(dec=="NO"){dec=Test_iid(X, p=p)$dec}
  if(dec!="NO"){    return( score=1)}

  ## YANG
  if(dec=="NO"){
    gamma_hat = Estim_gamma_indicator(X=X,min=1,max=5)
    part= partition(X)
    if(length(part$j)==1){dec="NO"     }else{     dec=Test_Yang_Pearson(X=X,gamma=gamma_hat,estimated=1,p=p)$dec}}
  if(dec!="NO"){    return( score=4)}

  ## LDM
  if(dec=="NO"){ dec=Test_LDM_Regression(X, p=p)$dec}
  if(dec!="NO"){    return( score=3)}


  ## No model
  if(dec=="NO"){ return(score)}

}

######## 21- DYCL ########
DT_DYCL = function(X,p=0.05){
  score=0

  ##DTRW
  dec=Test_DTRW_bonf(X, p=p)$dec
  if(dec!="NO"){    return( score=2)}

  ## YANG
  if(dec=="NO"){
    gamma_hat = Estim_gamma_indicator(X=X,min=1,max=5)
    part= partition(X)
    if(length(part$j)==1){dec="NO"     }else{    dec=Test_Yang_Pearson(X=X,gamma=gamma_hat,estimated=1,p=p)$dec}}
  if(dec!="NO"){    return( score=4)}

  ##Classical
  if(dec=="NO"){dec=Test_iid(X, p=p)$dec}
  if(dec!="NO"){    return( score=1)}

  ## LDM
  if(dec=="NO"){ dec=Test_LDM_Regression(X, p=p)$dec}
  if(dec!="NO"){    return( score=3)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 22- DYLC ########
DT_DYLC = function(X,p=0.05){
  score=0

  ##DTRW
  dec=Test_DTRW_bonf(X, p=p)$dec
  if(dec!="NO"){    return( score=2)}

  ## YANG
  if(dec=="NO"){
    gamma_hat = Estim_gamma_indicator(X=X,min=1,max=5)
    part= partition(X)
    if(length(part$j)==1){dec="NO"     }else{ dec=Test_Yang_Pearson(X=X,gamma=gamma_hat,estimated=1,p=p)$dec}}
  if(dec!="NO"){    return( score=4)}

  ## LDM
  if(dec=="NO"){
    dec=Test_LDM_Regression(X, p=p)$dec}
  if(dec!="NO"){    return( score=3)}

  ##Classical
  if(dec=="NO"){dec=Test_iid(X, p=p)$dec}
  if(dec!="NO"){    return( score=1)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 23- DLYC ########
DT_DLYC = function(X,p=0.05){
  score=0

  ##DTRW
  dec=Test_DTRW_bonf(X, p=p)$dec
  if(dec!="NO"){    return( score=2)}

  ## LDM
  if(dec=="NO"){ dec=Test_LDM_Regression(X, p=p)$dec}
  if(dec!="NO"){    return( score=3)}

  ## YANG
  if(dec=="NO"){
    gamma_hat = Estim_gamma_indicator(X=X,min=1,max=5)
    part= partition(X)
    if(length(part$j)==1){dec="NO"     }else{dec=Test_Yang_Pearson(X=X,gamma=gamma_hat,estimated=1,p=p)$dec}}
  if(dec!="NO"){    return( score=4)}

  ##Classical
  if(dec=="NO"){dec=Test_iid(X, p=p)$dec}
  if(dec!="NO"){    return( score=1)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 24- DLCY ########
DT_DLCY = function(X,p=0.05){
  score=0

  ##DTRW
  dec=Test_DTRW_bonf(X, p=p)$dec
  if(dec!="NO"){    return( score=2)}

  ## LDM
  if(dec=="NO"){ dec=Test_LDM_Regression(X, p=p)$dec}
  if(dec!="NO"){    return( score=3)}

  ##Classical
  if(dec=="NO"){ dec=Test_iid(X, p=p)$dec}
  if(dec!="NO"){    return( score=1)}

  ## YANG
  if(dec=="NO"){
    gamma_hat = Estim_gamma_indicator(X=X,min=1,max=5)
    part= partition(X)
    if(length(part$j)==1){dec="NO"     }else{dec=Test_Yang_Pearson(X=X,gamma=gamma_hat,estimated=1,p=p)$dec}}
  if(dec!="NO"){    return( score=4)}

  ## No model
  if(dec=="NO"){ return(score)}

}

