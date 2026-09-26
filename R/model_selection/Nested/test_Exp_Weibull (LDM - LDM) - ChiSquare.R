
# Load necessary packages
library(stats)
#devtools::load_all("~/Records")
## Note that Exp(scale=x) is equivalant to Weibull (shape = 1, scale =x)
## change functions to generate H0 and H1
## Change Likelihood functions H0 and H1
## change T, simulations and alpha
## change trend_values
## change H0 Parameters "par_H0"
## change bounds and starting values lb_0, ub_0 and x0_0
## change par_H1_test
################################# functions H0| H1 #################################

## Generate series H0
series_H0 = function(T,trend,par_H0){  ##
  ldm_series(T=T, dist = "exp",theta=trend, rate = par_H0[1])  ## rate is 1/scale
}

## Generate series H1
series_H1=function(T,trend,par_H1){  ## Normal Gumbel

  ldm_series(T=T,dist="weibull", theta=trend, scale = par_H1[2], shape = par_H1[1])
}

## Likelihood function H0
Likelihood_under_H0 = function(data_rec, params){
  # Likelihood_Rn_Exp_LDM(R=R,L=L,T=T,params=params)
  logLik_fun_rec = loglik_registry[["LDM"]][["records"]][["exp"]]
  return(logLik_fun_rec(data = data_rec, params = params))
}

## Likelihood function H1
Likelihood_under_H1 = function(data_rec, params){
  logLik_fun_rec = loglik_registry[["LDM"]][["records"]][["weibull"]]
  return(logLik_fun_rec(data = data_rec, params = params))
}

################################# Parameters ###################################

##### Parameters
T_val <- 100  ## Time, series length
alpha <- 0.05
simulation=100  ## Simulations numbers

final = matrix(0, nrow=11, ncol=4)
kk=1

for(trend_val in seq(0.1, 0.2, by=0.05)){

    ########## step 1: Generate a series under H0

true_params = c(trend = trend_val, rate = 0.5) ## true parameters were are generating scale: "1/rate"

#results = as.data.frame(matrix(0, ncol=10, nrow=simulation))

results <- data.frame(
  theta_H0 = rep(NA, simulation),
  rate_H0  = NA,
  logLik_H0 = NA,

  theta_H1 = NA,
  shape_H1 = NA,
  scale_H1 = NA,
  logLik_H1 = NA,

  LR       = NA,
  p_value  = NA,
  nRecords = NA
)

############################## Part A: Type I Error ##########################
trial=1

while (trial <= simulation){

  xt = series_H0(T=T_val, trend= true_params["trend"], par_H0 =true_params["rate"])

  ## Records and indicators
  R = rec_values(xt)
  L = rec_times(xt)

  while(length(R)<=1) {
    xt = series_H0(T = T_val, trend= true_params["trend"], par_H0 =true_params["rate"])
    R = rec_values(xt)
    L = rec_times(xt)
  }  ## only one record, ignore

  m <- length(L)

  data_rec = list(rec_values = R, rec_times = L, time = T_val)

  ########## step 2: Fit the parameters using Likelihood of H0

  lb_0 <- list(theta = 0.01, rate = 0.01)  ## Lower bounds
  ub_0 <- list(theta = min(R/L), rate = 10) ## Upper bounds
  x0_0 <- list(theta = 0.01, rate = 0.01)      ## Initial values

  ## Define the likelihood function
  Likelihood_H0 <- function(params) { data_rec= data_rec ; return(-Likelihood_under_H0(data_rec, params))}

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
      x0_0 <- x0_0 + c(0.01,0.1)              # Increment starting values by a small amount
      retry_count <- retry_count + 1

      if (retry_count >= max_retries) {
        message("Max retries reached for iteration ")
        break                        # Exit after max retries
      }
    }
  }

  AIC_0 = 2*param_H0_start[length(param_H0_start)]-2*(length(param_H0_start)-2)
  param_H0_start  ## Gamma, rate, LogL
  AIC_0

  ##########  step 3: Fit the parameters using Likelihood of H1
  lb_1 <- list(theta = 0.001, shape = 0.01, scale = param_H0_start["rate"])   ## Lower bounds for gamma, location, scale
  ub_1 <- list(theta = min(R/L), shape = 5, scale =param_H0_start["rate"])             ## Upper bounds for gamma, location, scale
  x0_1 <- list(theta = 0.1, shape = 0.01, scale = param_H0_start["rate"])      ## Initial values

  ## Define the likelihood function
  Likelihood_H1 <- function(params) {data_rec = data_rec; return(-Likelihood_under_H1(data_rec, params))}

  ## Repeat optimization until successful or after a maximum number of retries
  max_retries <- 70      # Define a maximum number of retries to prevent infinite loops
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

      param_H1_start = c(MLE_1$par, -MLE_1$objective)
      break                            # Exit the repeat loop if optimization is successful
    } else {
      # If optimization failed, increment x0_1 slightly and retry
      x0_1 <- x0_1 + c(0.01,0.1,0)              # Increment starting values by a small amount
      retry_count <- retry_count + 1

      if (retry_count >= max_retries) {
        message("Max retries reached for iteration ")
        break                        # Exit after max retries
      }
    }
  }

  ## gamma, parameters, LogL
  AIC_1 = 2*param_H1_start[length(param_H1_start)]-2*(length(param_H1_start)-2)
  param_H1_start
  AIC_1


  ################################# Data Frames ###################################
  # Initialize the variables and start loop over gamma values
  wilks = 2*(param_H1_start[length(param_H1_start)]- param_H0_start[length(param_H0_start)])
  if(wilks <=0) {next } else{
  p_value = 1-pchisq(q=wilks, df= 2-1)
  results[trial,] = as.numeric(c(param_H0_start,
                                        param_H1_start,
                                        wilks,
                                        p_value,
                                        m))
  trial=trial+1
  rm(param_H0_start); rm( param_H1_start)}
}

apply(results,2,mean)

#################### Part B: Perform  power calculations ####################

results_power <- data.frame(
  theta_H0  = rep(NA, simulation),
  rate_H0   = NA,
  logLik_H0 = NA,

  theta_H1  = NA,
  shape_H1  = NA,
  scale_H1  = NA,
  logLik_H1 = NA,

  LR        = NA,
  reject    = NA,
  nRecords  = NA
)

trial=1
while(trial <=simulation) {

      x <- series_H1(T=T_val,
                     trend=true_params["trend"],
                     par_H1=c(1.2,true_params["rate"]))  ## shape scale

      ## Records and indicators
      L = rec_times(x)
      R = rec_values(x)
      m <- length(L)

      while(m<2 ){
        x <- series_H1(T=T_val,trend=as.numeric(results[trial,4]), par_H1=c(1.2,true_params["rate"]) )

        ## Records and indicators
        L = rec_times(x)
        R = rec_values(x)
        m <- length(L)
      }


    ######################## 2 - H0 likelihood maximization ########################

      lb_0 <- c(theta =  0.01, rate = 0.01)  ## gamma,
      ub_0 <- c(theta = min(R/L)-0.01, rate = 10)
      x0_0 <- c(theta = 0.1,rate = 0.1)

      Likelihood_H0 = function(params){-Likelihood_under_H0(data_rec = data_rec,params) }

      ## Repeat optimization until successful or after a maximum number of retries
      max_retries <- 50      # Define a maximum number of retries to prevent infinite loops
      retry_count <- 0

      repeat {
        ## MLE try and error
        MLE_0 <- tryCatch({
          # optim(par = x0_0,fn = Likelihood_Yang2_H0,method = "L-BFGS-B",lower = lb_0,upper = ub_0)
          nlminb (x0_0,Likelihood_H0,  lower=lb_0, upper = ub_0)
        }   )

        # Check if optimization was successful
        if (!is.null(MLE_0) && is.finite(MLE_0$objective)) {
          param_H0_test=c(MLE_0$par, -MLE_0$objective)
          break                            # Exit the repeat loop if optimization is successful
        } else {
          # If optimization failed, increment x0_1 slightly and retry
          x0_0 <- x0_0 + c(0.1,0.1)              # Increment starting values by a small amount
          retry_count <- retry_count + 1

          if (retry_count >= max_retries) {break}}


      }
      if (retry_count >= max_retries) {next;}
      param_H0_test
    ########################  3 - H1 Likelihood Maximization  ########################

      lb_1 <- c( theta = 0.1, shape = 0.01, scale = param_H0_test[["rate"]])  ## gamma, shape, scale
      ub_1 <- c(theta = min(R/L)-0.001,shape = 5, scale = param_H0_test[["rate"]])
      x0_1 <- c(theta = 0.1, shape=0.1, scale = param_H0_test[["rate"]])

      ## Define the likelihood function
      Likelihood_H1 <- function(params) { - Likelihood_under_H1(data_rec = data_rec, params)}

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

          param_H1_test = c(MLE_1$par, -MLE_1$objective)
          break                            # Exit the repeat loop if optimization is successful
        } else {
          # If optimization failed, increment x0_1 slightly and retry
          x0_1 <- x0_1 + c(0.01,0.1,0.1)              # Increment starting values by a small amount
          retry_count <- retry_count + 1

          if (retry_count >= max_retries) {
            message("Max retries reached for iteration ")
            break                        # Exit after max retries
          }
        }
      }
      #param_H1_test ## Gamma, shape and scale
      if (retry_count >= max_retries) {next;}
      param_H1_test

      wilks_test = 2*(param_H1_test[length(param_H1_test)]-param_H0_test[length(param_H0_test)])

      if(wilks_test<0) {next} else{

    power= ifelse(wilks_test>= qchisq(0.95, df=1), 1, 0)

    results_power[trial,]= c(param_H0_test, param_H1_test, wilks_test, power,m)
    trial=trial+1
    rm(param_H0_test); rm(param_H1_test)}
}

apply(results_power,2,mean)
################################ Results #####################################

cat("Type I Error:",nrow(results[results$sign<alpha,])*100/nrow(results),"%")  ## Should be 5%
cat("Power:",nrow(results_power[results_power$sign==1,])*100/nrow(results_power),"%")  ## Should be 5%


final[kk,] = c(trend_val, nrow(results[results$sign<alpha,])*100/nrow(results), nrow(results[results_power$sign==1,])*100/nrow(results_power), mean(results[,"m"]) )
kk=kk+1
}
## for Weibull (1.2, 0.5)
## 5% and Power 58%
## 0% and 50%
## 1% and 53#
## 3% and 54%
## 0 and 54
## 2% and 47%
## 3 and 49%
## 0 and 56%
## 1.8 and 51.4%
## 1.6 and 49.4%

## at 0.1, 0.5  we get 2% and 50%
## at 0.11, 0.5 we get 4.4% and 56%
## at 0.12, 0.5  we get 10% and 60%
## at 0.13, 0.5  we get 13% and 66%
## at 0.14, 0.5  we get %20 and 75%
## at 0.15, 0.5  we get 23% and 77%  m=15

colnames(final) = c("theta", "Type I Error", "Power", "Nb_record")

library(ggplot2)
library(scales)  # for percentage formatting if needed

# Ensure your data is a dataframe
final <- as.data.frame(final)

# Rescale Nb_record to match the scale of the left y-axis
scale_factor <- max(final$Power) / max(final$Nb_record)


ggplot(final, aes(x = theta)) +
  geom_line(aes(y = `Type I Error`, color = "Type I Error"), size = 1) +
  geom_line(aes(y = Power, color = "Power"), size = 1) +
  geom_line(aes(y = Nb_record * scale_factor, color = "Nb_record"), linetype = "dashed", size = 1) +

  # Left y-axis
  scale_y_continuous(
    name = "Test Power | Type I Error (%)",

    # Right y-axis
    sec.axis = sec_axis(~ . / scale_factor, name = "Number of Records")
  ) +

  scale_color_manual(values = c("Type I Error" = "lightblue", "Power" = "darkblue", "Nb_record" = "#2E8B57")) +

  labs(x = "Theta", color = "Metric") +
  theme_classic() +
  theme(
    axis.title.y.right = element_text(color = "#2E8B57"),
    axis.title.y.left = element_text(color = "black"),
    legend.position = "top"
  )
