
# Load necessary packages
library(stats)
devtools::load_all("~/Records")
## change functions to generate H0 and H1
## Change Likelihood functions H0 and H1
## change T, simulations and alpha
## change gama_values
## change H0 Parameters "par_H0"
## change bounds and starting values lb_0, ub_0 and x0_0
## change par_H1_test
################################# functions H0| H1 #################################

## Generate series H0
series_H0 = function(T,gama,par_H0){  ## Gumbel
  Yang_series_Gumbel(T=T, gamma=gama, loc = par_H0[1], scale = par_H0[2])
}

## Generate series H1
series_H1=function(T,par_H1){  ## Normal Gumbel

  VGAM::rgumbel(n=T, loc = par_H1[1], scale = par_H1[2])

}

## Likelihood function H0
Likelihood_under_H0 = function(R,L,T,gAa){
  Likelihood_Rn_Gumbel_Yang(R=R,L=L,T=T,gAa=gAa)
}

## Likelihood function H1
Likelihood_under_H1 = function(R,L,T,gAa){
  Likelihood_Rn_Gumbel_iid(R=R,L=L,T=T,Aa=gAa)
}

################################# Parameters ###################################

##### Parameters
T <- 50  ## Time, series length
alpha <- 0.05
simulation=100  ## Simulations numbers

########## step 1: Generate a series under H0
true_gAa = c(1.2,0,1) ## true parameters were are generating mean and variance

xt = series_H0(T=T, gama= true_gAa[1], par_H0 = c(true_gAa[2], true_gAa[3]))

## Records and indicators
L = Ln(xt)
R = Rn(xt)
m <- length(L)

########## step 2: Fit the parameters using Likelihood of H0

lb_0 <- c(1.0001, -5, 0.01)  ## Lower bounds for gamma, location, scale
ub_0 <- c(5, 5, 5)             ## Upper bounds for gamma, location, scale
x0_0 <- c(1.01, 0.01,0.01)      ## Initial values

## Define the likelihood function
Likelihood_H0 <- function(gAa) { R = R; L = L ; T = T; return(-Likelihood_under_H0(R = R, L = L, T = T, gAa))}

## Repeat optimization until successful or after a maximum number of retries
max_retries <- 50      # Define a maximum number of retries to prevent infinite loops
retry_count <- 0

repeat {
  ## MLE try and error
  MLE_0 <- tryCatch({
    nlminb (x0_0,Likelihood_H0,  lower=lb_0, upper = ub_0)
  }, error = function(e) {
    message("Error in iteration of H0 Trial ",retry_count,": ", e$message)
    return(NULL)                    # Return NULL to indicate failure
  })

  # Check if optimization was successful
  if (!is.null(MLE_0) && is.finite(MLE_0$objective)) {

    param_H0_start = c(MLE_0$par, -MLE_0$objective)
    break                            # Exit the repeat loop if optimization is successful
  } else {
    # If optimization failed, increment x0_1 slightly and retry
    x0_0 <- x0_0 + c(0,0.1,0.1)              # Increment starting values by a small amount
    retry_count <- retry_count + 1

    if (retry_count >= max_retries) {
      message("Max retries reached for iteration ")
      break                        # Exit after max retries
    }
  }
}

param_H0_start  ## Gamma, mean, scale, LogL

##########  step 3: Fit the parameters using Likelihood of H1
lb_1 <- c( -5, 0.01)   ## Lower bounds for gamma, location, scale
ub_1 <- c( 100, 100)             ## Upper bounds for gamma, location, scale
x0_1 <- c( 0.1, 1)      ## Initial values

## Define the likelihood function
Likelihood_H1 <- function(gAa) {R = R; L = L ; T = T; return(-Likelihood_under_H1(R = R, L = L, T = T, gAa))}

## Repeat optimization until successful or after a maximum number of retries
max_retries <- 50      # Define a maximum number of retries to prevent infinite loops
retry_count <- 0

repeat {
  ## MLE try and error
  MLE_1 <- tryCatch({
    nlminb (x0_1,Likelihood_H1,  lower=lb_1, upper = ub_1)
  }, error = function(e) {
    message("Error in iteration of H1 Trial ",retry_count,": ", e$message)
    return(NULL)                    # Return NULL to indicate failure
  })

  # Check if optimization was successful
  if (!is.null(MLE_1) && is.finite(MLE_1$objective)) {

    param_H1_start = c(1,MLE_1$par, -MLE_1$objective)
    break                            # Exit the repeat loop if optimization is successful
  } else {
    # If optimization failed, increment x0_1 slightly and retry
    x0_1 <- x0_1 + c(0.1,0.1)              # Increment starting values by a small amount
    retry_count <- retry_count + 1

    if (retry_count >= max_retries) {
      message("Max retries reached for iteration ")
      break                        # Exit after max retries
    }
  }
}

## parameters, LogL
param_H1_start

##########  View and compare
par(mfrow=c(2,2))
my_bar = barplot(c(true_gAa[1],param_H0_start[1], param_H1_start[1]), width = 0.1, names.arg = c("True", "H0","H1"), col=c("red", "pink", "lightblue"), main= "Gamma")
text(my_bar, rep(mean(c(true_gAa[1],param_H0_start[1],param_H1_start[1])),3), round(c(true_gAa[1],param_H0_start[1],param_H1_start[1]),3))

my_bar = barplot(c(true_gAa[2],param_H0_start[2], param_H1_start[2]), width = 0.1, names.arg = c("True", "H0","H1"), col=c("red", "pink", "lightblue"), main= "Param1")
text(my_bar, rep(mean(c(true_gAa[2],param_H0_start[2], param_H1_start[2])),3), round(c(true_gAa[2],param_H0_start[2], param_H1_start[2]),3))

my_bar = barplot(c(true_gAa[3],param_H0_start[3], param_H1_start[3]), width = 0.1, names.arg = c("True", "H0","H1"), col=c("red", "pink", "lightblue"), main= "Param2")
text(my_bar, rep(mean(c(true_gAa[3],param_H0_start[3], param_H1_start[3])),3), round(c(true_gAa[3],param_H0_start[3], param_H1_start[3]),3))

param_H0_start
param_H1_start

plot(y=xt, x=c(1:T), xlim=c(1,T), ylim=c(min(xt),max(xt)), xlab = "Time", ylab="Xt")
par(new=TRUE)
plot(y=Rn(xt), x=Ln(xt),xlim=c(1,T),ylim=c(min(xt),max(xt)), col="red", xlab="",ylab="")
##################### Go to the code of hypothesis testing and put the parameters there






################################# Data Frames ###################################
# Initialize the variables and start loop over gamma values
gama <- param_H0_start[1] ## the one estimated for the series above
CI = as.data.frame(matrix(0,ncol=4,nrow=length(gama)))
colnames(CI)=c("Gamma","Lower", "Upper","Power")
CI[,1]=gama


## H0 Parameters
par_H0 = param_H0_start[2:3] ## estimated for the series above
par_H1 = param_H1_start[2:3] ## Dummy iid

## create matrixes to store MaxLik results
param_H0 <- as.data.frame(matrix(NA, ncol = 3 + length(par_H0), nrow = simulation*length(gama)))
param_H1 <- as.data.frame(matrix(NA, ncol = 3 + length(par_H1), nrow = simulation*length(gama)))

colnames(param_H0)=c("Gamma","GammaHat",paste0("Param",1:length(par_H0)),"MaxLik")
colnames(param_H1)=c("Gamma","GammaHat",paste0("Param",1:length(par_H1)),"MaxLik")

ignore_sim_H0 = numeric()
ignore_sim_H1 = numeric()

################################# Part A: The Critical Region of Lambda #################################
for (d in 1) {

  ## create matrixes to store records
  L_x <- matrix(NA, nrow = T, ncol = simulation)
  R_x <- matrix(NA, nrow = T, ncol = simulation)
  nbr_records <- numeric(simulation)

  #######################  1- Generate simulated records from series under H0 #######################
  sim=1
  while(sim <= simulation) {
    x <- series_H0(T=T,gama=gama[d],par_H0=c(par_H0[1], sqrt(par_H0[2]))) ##under H0: DTRW

    ## Records and indicators: consider cases of more than 2 records
    L = Ln(x)
    R = Rn(x)
    m <- length(L)
    if(m>2){
      L_x[1:m, sim] <- Ln(x)
      R_x[1:m, sim] <- Rn(x)

      nbr_records[sim] <- m
      sim=sim+1}
  }

  ####################### 2 - Perform H0 likelihood maximization #######################

  for (sim in 1:simulation) {
    lb_0 <- c( 1.00001,-5,0.001)  ## Lower bounds for location, scale
    ub_0 <- c(5,5,5)             ## Upper bounds for gamma, location, scale
    x0_0 <- c(1.001, 0,1)      ## Initial values

    Likelihood_H0 = function(gAa){ R=R_x[1:nbr_records[sim], sim]; L=L_x[1:nbr_records[sim], sim]; T=T; return(-Likelihood_under_H0(R=R,L=L,T=T,gAa)) }

    ## Repeat optimization until successful or after a maximum number of retries
    max_retries <- 50      # Define a maximum number of retries to prevent infinite loops
    retry_count <- 0

    repeat {
      ## MLE try and error
      MLE_0 <- tryCatch({
        nlminb (x0_0,Likelihood_H0,  lower=lb_0, upper = ub_0)
      }, error = function(e) {
        message("Error in iteration of H0 nb ", sim, " of gamma ", gama[d], "Trial ",retry_count,": ", e$message)
        return(NULL)                    # Return NULL to indicate failure
      })

      # Check if optimization was successful
      if (!is.null(MLE_0) && is.finite(MLE_0$objective)) {
        # Store results if optimization succeeded
        param_H0[sim + (simulation * (d - 1)), ] <- c(gama[d], MLE_0$par, -MLE_0$objective)
        break                            # Exit the repeat loop if optimization is successful
      } else {
        # If optimization failed, increment x0_1 slightly and retry
        x0_0 <- x0_0 + c(0,0.2,0.2)              # Increment starting values by a small amount
        retry_count <- retry_count + 1

        if (retry_count >= max_retries) {
          message("Max retries reached for iteration ", sim, " of gamma ", gama[d])
          ignore_sim_H0 <- c(ignore_sim_H0, sim + (simulation * (d - 1)))
          break                        # Exit after max retries
        }
      }
    }
  }

  ####################### 3 - Perform H1 likelihood maximization #######################

  for (sim in 1:simulation) {
    lb_1 <- c(-5,0.01)  ## Lower bounds for location, scale
    ub_1 <- c(5,5)             ## Upper bounds for gamma, location, scale
    x0_1 <- c( 0,1)      ## Initial values

    ## Define the likelihood function
    Likelihood_H1 <- function(gAa) { R <- R_x[1:nbr_records[sim], sim]; L <- L_x[1:nbr_records[sim], sim]; T <- T; return(-Likelihood_under_H1(R = R, L = L, T = T, gAa))}

    ## Repeat optimization until successful or after a maximum number of retries
    max_retries <- 100      # Define a maximum number of retries to prevent infinite loops
    retry_count <- 0

    repeat {
      ## MLE try and error
      MLE_1 <- tryCatch({
        # optim(par = x0_1,
        #       fn = Likelihood_Yang2_H1,
        #       method = "L-BFGS-B",
        #       lower = lb_1,
        #       upper = ub_1)
        nlminb (x0_1,Likelihood_H1,  lower=lb_1, upper = ub_1)
      }, error = function(e) {
        message("Error in iteration of H1 nb ", sim, "Trial ",retry_count,": ", e$message)
        return(NULL)                    # Return NULL to indicate failure
      })

      # Check if optimization was successful
      if (!is.null(MLE_1) && is.finite(MLE_1$objective)) {
        # Store results if optimization succeeded
        param_H1[sim + (simulation * (d - 1)), ] <- c(gama[d],1, MLE_1$par, -MLE_1$objective)
        break                            # Exit the repeat loop if optimization is successful
      } else {
        # If optimization failed, increment x0_1 slightly and retry
        x0_1 <- x0_1 + c(0.2,0.2)              # Increment starting values by a small amount
        retry_count <- retry_count + 1

        if (retry_count >= max_retries) {
          message("Max retries reached for iteration ", sim, " of gamma ", gama[d])
          ignore_sim_H1 <- c(ignore_sim_H1, sim + (simulation * (d - 1)))
          break                        # Exit after max retries
        }
      }
    }
  }




}

####################### 4 - Wilk's Lamda and CI #######################
### Ignore simulations
ignore_sim = c(ignore_sim_H1,ignore_sim_H0)
if(length(ignore_sim)>0){
  param_H1=param_H1[-ignore_sim,]
  param_H0=param_H0[-ignore_sim,]
}


for(d in 1:length(gama)){
  maxlogL_H1 = param_H1[param_H1$Gamma==gama[d],"MaxLik"]
  maxlogL_H0 = param_H0[param_H0$Gamma==gama[d],"MaxLik"]

  Rapport <- maxlogL_H1 - maxlogL_H0  ## Log L1 - Log L0

  Rapport_sorted <- sort(Rapport[is.finite(Rapport)])
  ####################### 4 - confidence bounds #######################

  CI[d,2]= quantile(Rapport_sorted,alpha/2) ##LOwer
  CI[d,3]= quantile(Rapport_sorted,1-alpha/2) ##Upper
}

#######################################################################################################
#################### Part B: Perform  power calculations ####################

## parameters under H0 and H1
par_H0_test = param_H0_start[2:3]##Dummy
par_H1_test = param_H1_start[2:3] ## the ones estimated from above

## Store values
param_H0_test <- as.data.frame(matrix(NA, ncol = 3 + length(par_H0_test), nrow = simulation*length(gama)))
param_H1_test <- as.data.frame(matrix(NA, ncol = 3 + length(par_H1_test), nrow = simulation*length(gama)))

colnames(param_H0_test)=c("Gamma","GammaHat",paste0("Param",1:length(par_H0_test)),"MaxLik")
colnames(param_H1_test)=c("Gamma","GammaHat",paste0("Param",1:length(par_H1_test)),"MaxLik")


for(d in 1) {

  ## Initializing Ltest and Rtest with the first record
  L_x_test <- matrix(NA, nrow = T, ncol = simulation)
  R_x_test <- matrix(NA, nrow = T, ncol = simulation)

  maxlogL_H1_test <- numeric(simulation)
  maxlogL_H0_test <- numeric(simulation)

  nbr_records_test <- numeric(simulation)

  ############ 1 - Generating data under H1    ############

  sim=1
  while (sim <= simulation) {

    x <- series_H1(T=T, c(par_H1_test[1],sqrt(par_H1_test[2]))  )

    ## Records and indicators
    L = Ln(x)
    R = Rn(x)
    m <- length(L)

    if(m>2){
      L_x_test[1:m, sim] <- L
      R_x_test[1:m, sim] <- R

      nbr_records_test[sim] <- m
      sim=sim+1}
  }

  ######################## 2 - H0 likelihood maximization ########################

  for (sim in 1:simulation) {
    lb_0 <- c( 1.0001,-5,0.01)  ## gamma, location, scale
    ub_0 <- c(5,5,5)
    x0_0 <- c(1.1,0.1,0.1)

    Likelihood_H0 = function(gAa){ R=R_x_test[1:nbr_records_test[sim], sim]; L <- L_x_test[1:nbr_records_test[sim], sim]; T=T; return(-Likelihood_under_H0(R=R,L=L,T=T,gAa)) }

    ## Repeat optimization until successful or after a maximum number of retries
    max_retries <- 50      # Define a maximum number of retries to prevent infinite loops
    retry_count <- 0

    repeat {
      ## MLE try and error
      MLE_0 <- tryCatch({
        # optim(par = x0_0,fn = Likelihood_Yang2_H0,method = "L-BFGS-B",lower = lb_0,upper = ub_0)
        nlminb (x0_0,Likelihood_H0,  lower=lb_0, upper = ub_0)
      },  error = function(e) {message("Error in iteration of H0 nb ", sim, " of gamma ", gama[d], "Trial ",retry_count,": ", e$message); return(NULL)}                    # Return NULL to indicate failure
      )

      # Check if optimization was successful
      if (!is.null(MLE_0) && is.finite(MLE_0$objective)) {
        param_H0_test[sim + (simulation * (d - 1)), ] <- c(gama[d], MLE_0$par, -MLE_0$objective)
        break                            # Exit the repeat loop if optimization is successful
      } else {
        # If optimization failed, increment x0_1 slightly and retry
        x0_0 <- x0_0 + c(0,0.1,0.1)              # Increment starting values by a small amount
        retry_count <- retry_count + 1

        if (retry_count >= max_retries) {
          message("Max retries reached for iteration ", sim, " of gamma ", gama[d])
          ignore_sim_H0 <- c(ignore_sim_H0, sim + (simulation * (d - 1)))
          break                        # Exit after max retries
        }
      }

    }
  } ## end of sim of H0

  ########################  3 - H1 Likelihood Maximization  ########################

  for (sim in 1:simulation) {
    lb_1 <- c(-5, 0.01)  ## Lower bounds for gamma, location, scale
    ub_1 <- c(5,5)             ## Upper bounds for gamma, location, scale
    x0_1 <- c( 0.01,0.01)      ## Initial values

    ## Define the likelihood function
    Likelihood_H1 <- function(gAa) { R <- R_x_test[1:nbr_records_test[sim], sim]; L <- L_x_test[1:nbr_records_test[sim], sim]; T = T; return(-Likelihood_under_H1(R = R, L = L, T = T, gAa))}

    ## Repeat optimization until successful or after a maximum number of retries
    max_retries <- 50      # Define a maximum number of retries to prevent infinite loops
    retry_count <- 0

    repeat {
      ## MLE try and error
      MLE_1 <- tryCatch({
        # optim(par = x0_1,fn = Likelihood_Yang2_H1,method = "L-BFGS-B",lower = lb_1,upper = ub_1)
        nlminb (x0_1,Likelihood_H1,  lower=lb_1, upper = ub_1)
      }, error = function(e) { message("Error in iteration of H1 nb ", sim, " of gamma ", gama[d], "of trial", retry_count, ": ", e$message)
        return(NULL)                    # Return NULL to indicate failure
      })

      # Check if optimization was successful
      if (!is.null(MLE_1) && is.finite(MLE_1$objective)) {
        param_H1_test[sim + (simulation * (d - 1)), ] <-  c(gama[d],1,MLE_1$par,-MLE_1$objective)
        break                                 } else { ## If optimization failed, increment x0_1 slightly and retry
          x0_1 <- x0_1 + c(0.1,0.1)
          retry_count <- retry_count + 1

          if (retry_count >= max_retries) {
            message("Max retries reached for iteration ", sim, " of gamma ", gama[d])
            ignore_sim_H1 <- c(ignore_sim_H1, sim + (simulation * (d - 1)))
            break                        # Exit after max retries
          }
        }
    }
  }


}

######################## Calculate likelihood ratios if valid estimates are obtained ########################
ignore_sim_test = c(which(is.na(param_H0_test)),which(is.na(param_H1_test)))
if(length(ignore_sim_test)>0){
  param_H0_test = param_H0_test[-ignore_sim_test,]
  param_H1_test = param_H1_test[-ignore_sim_test,]
}

###########Compute Rapport Lambda
for(d in 1){

  maxlogL_H1_test = param_H1_test[param_H1_test$Gamma==gama[d],"MaxLik"]
  maxlogL_H0_test = param_H0_test[param_H0_test$Gamma==gama[d],"MaxLik"]

  Rapport_test <- maxlogL_H1_test - maxlogL_H0_test  ## Log L1 - Log L0

  ########################Test power calculation based on confidence interval bounds ########################

  CI[d,"Power"] = sum(Rapport_test > CI[d,"Upper"] | Rapport_test < CI[d,"Lower"])*100/length(maxlogL_H0_test)
}

###################################################################################################

################### Bias ################
Bias_H0 = as.data.frame(matrix(0, nrow=length(gama), ncol = 2+ 2*length(par_H0)))  ## 2 for the gamma + *2 for parameters
colnames(Bias_H0) = c("Gamma", "Bias",paste(c("Param_", "Bias_"),rep(1:length(par_H0),each=2)))
Bias_H1_test = Bias_H0

Estim_H1 = as.data.frame(matrix(0, nrow=length(gama), ncol = 2+length(par_H1)))  ## 2 for the gamma + each for parameters
colnames(Estim_H1) = c("Gamma","Bias",paste("Param_",1:length(par_H1)) )
Estim_H0_test = Estim_H1

## Bias of H0
for (d in 1:length(gama)){
  temp=subset(param_H0,Gamma==gama[d])
  Bias_H0[d,1] = gama[d]
  Bias_H0[d,2] = mean(temp[,2]-temp[,1])
  Bias_H0[d,3] = par_H0[1]
  Bias_H0[d,4] = mean(temp[,3]-par_H0[1])
  if(length(par_H0)>1){
    for(p in 2:length(par_H0)){
      Bias_H0[d,2*p+1] = par_H0[p]
      Bias_H0[d,2*p+2] = mean(temp[,p+1]-par_H0[p])
    }
  }
}

## average estimation under H1
for (d in 1:length(gama)){
  temp=subset(param_H1,Gamma==gama[d])
  Estim_H1[d,1] = 0## gama[d]
  Estim_H1[d,2] = 0 ## mean(temp[,2]-temp[,1])
  Estim_H1[d,3] = mean(temp[,3])  ## first parameter
  if(length(par_H1)>1){
    for(p in 2:length(par_H1)){
      Estim_H1[d,2+p] = mean(temp[,p+2])
    }
  }
}

## Test: bias of H1
for (d in 1:length(gama)){
  temp=subset(param_H1_test,Gamma==gama[d])
  Bias_H1_test[d,1] = 0## gama[d]
  Bias_H1_test[d,2] = 0 ## mean(temp[,2]-temp[,1])
  Bias_H1_test[d,3] = par_H1_test[1]
  Bias_H1_test[d,4] = mean(temp[,3]-par_H1_test[1])
  if(length(par_H1_test)>1){
    for(p in 2:length(par_H1_test)){
      Bias_H1_test[d,2*p+1] = par_H1_test[p]
      Bias_H1_test[d,2*p+2] = mean(temp[,p+1]-par_H1_test[p])
    }
  }
}

## Test: average estimation under H0
for (d in 1:length(gama)){
  temp=subset(param_H0_test,Gamma==gama[d])
  Estim_H0_test[d,1] = gama[d]
  Estim_H0_test[d,2] = mean(temp[,2]-temp[,1])
  Estim_H0_test[d,3] = mean(temp[,3])  ## first parameter
  if(length(par_H0_test)>1){
    for(p in 2:length(par_H0_test)){
      Estim_H0_test[d,2+p] = mean(temp[,p+2])
    }
  }
}

################################### Show Results ####################################
CI

cat("your input parameters to generate one series are\n", true_gAa)
cat("The estimated parameters of the input series under H0:\n", param_H0_start)
cat("The estimated parameters of the input series under H1:\n", param_H1_start)


## Dataset of estimated parameters under H0 of a H0 series
param_H0
param_H1 ## Dataset of estimated parameters under H1 of a H0 series

## Dataset of estimated parameters under H0 of a H1 series
param_H0_test
param_H1_test ## Dataset of estimated parameters under H1 of a H1 series

Bias_H0 ##average estimated values under H0 of the original series
Bias_H1_test ##average estimated values under H1 of the original series

Estim_H1  ## estimated values of H0 generated series under H1
Estim_H0_test ## estimated values of H1 generated series under H0
