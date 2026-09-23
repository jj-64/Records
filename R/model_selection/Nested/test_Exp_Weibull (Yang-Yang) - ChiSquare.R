
# Load necessary packages
library(stats)
devtools::load_all("~/Records")
## Note that Exp(scale=x) is equivalant to Weibull (shape = 1, scale =x)
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
  Yang_series_Exp(T=T, gamma=gama, rate = par_H0[1])

}

## Generate series H1
series_H1=function(T,gama,par_H1){  ## Normal Gumbel

  Yang_series_Weibull(T=T, gamma=gama, scale = par_H1[2], shape = par_H1[1])
}

## Likelihood function H0
Likelihood_under_H0 = function(R,L,T,gAa){
  Likelihood_Rn_Exp_Yang(R=R,L=L,T=T,gAa=gAa)
}

## Likelihood function H1
Likelihood_under_H1 = function(R,L,T,gAa){
  Likelihood_Rn_Weibull_Yang(R=R,L=L,T=T,gAa=gAa)
}

################################# Parameters ###################################

##### Parameters
T <- 75  ## Time, series length
alpha <- 0.05
simulation=100  ## Simulations numbers

final = matrix(0, nrow=6, ncol=4)
kk=1
for(thetas in seq(1.1, 1.6, by=0.1)){
########## step 1: Generate a series under H0
true_gAa = c(thetas,0.5) ## true parameters were are generating scale: "1/rate"
summary = as.data.frame(matrix(0, ncol=10, nrow=simulation))
############################## Part A: Type I Error ##########################
trial=1
while (trial <= simulation){

  xt = series_H0(T=T, gama= true_gAa[1], par_H0 =1/true_gAa[2])

  ## Records and indicators
  L = Ln(xt)
  R = Rn(xt)
  m <- length(L)
  while(m<2){
    xt = series_H0(T=T, gama= true_gAa[1], par_H0 =1/true_gAa[2])

    ## Records and indicators
    L = Ln(xt)
    R = Rn(xt)
    m <- length(L)
  }
  ########## step 2: Fit the parameters using Likelihood of H0

  lb_0 <- c(1.0001, 0.01)  ## Lower bounds for gamma,  scale
  ub_0 <- c(5, 10)             ## Upper bounds for gamma, scale
  x0_0 <- c(1.01, 0.01)      ## Initial values

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
      x0_0 <- x0_0 + c(0.1,0.1)              # Increment starting values by a small amount
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
  lb_1 <- c( 1.001, 0.01,0.01)#param_H0_start[2])   ## Lower bounds for gamma, location, scale
  ub_1 <- c( 5,5,5)#param_H0_start[2])             ## Upper bounds for gamma, location, scale
  x0_1 <- c( 0.1, 0.01,0.01)#param_H0_start[2])      ## Initial values

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

      param_H1_start = c(MLE_1$par, -MLE_1$objective)
      break                            # Exit the repeat loop if optimization is successful
    } else {
      # If optimization failed, increment x0_1 slightly and retry
      x0_1 <- x0_1 + c(0.1,0.1,0.1)              # Increment starting values by a small amount
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
  if(wilks <=0) next;
  p_value = 1-pchisq(q=wilks, df= 2-1)
  summary[trial,] = as.numeric(c(param_H0_start,param_H1_start,wilks, p_value,m))
  trial=trial+1
}

colnames(summary) = c("Gamma_0", "H_0[1]", "LogL_0" ,"Gamma_1","H_1[1]", "H_1[2]","LogL_1","Wilks", "sign","m")


#################### Part B: Perform  power calculations ####################

summary_test =  as.data.frame(matrix(0, ncol=10, nrow=simulation))
trial=1
while(trial <=simulation) {

      x <- series_H1(T=T,gama=as.numeric(summary[trial,4]), par_H1=c(0.8,true_gAa[2]))  ## shape scale

      ## Records and indicators
      L = Ln(x)
      R = Rn(x)
      m <- length(L)

      while(m<2 ){
        x <- series_H1(T=T,gama=as.numeric(summary[trial,4]), par_H1=as.numeric(summary[trial,5:6]  ))

        ## Records and indicators
        L = Ln(x)
        R = Rn(x)
        m <- length(L)
      }


    ######################## 2 - H0 likelihood maximization ########################

      lb_0 <- c( 1.0001,0.01)  ## gamma,
      ub_0 <- c(5,10)
      x0_0 <- c(1.1,0.1)

      Likelihood_H0 = function(gAa){ R=R; L=L ; T=T; return(-Likelihood_under_H0(R=R,L=L,T=T,gAa)) }

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

          if (retry_count >= max_retries) {break}
        }

      }
      if (retry_count >= max_retries) {next;}
      param_H0_test
    ########################  3 - H1 Likelihood Maximization  ########################

      lb_1 <- c( 1.01,0.01,0.01)  ## gamma, shape, scale
      ub_1 <- c(5,10,10)
      x0_1 <- c(2,1,1)

      ## Define the likelihood function
      Likelihood_H1 <- function(gAa) { R=R; L=L; T = T; return(-Likelihood_under_H1(R = R, L = L, T = T, gAa))}

      ## Repeat optimization until successful or after a maximum number of retries
      max_retries <- 50      # Define a maximum number of retries to prevent infinite loops
      retry_count <- 0

      repeat {
        ## MLE try and error
        MLE_1 <- tryCatch({   nlminb (x0_1,Likelihood_H1,  lower=lb_1, upper = ub_1) }, error = function(e) {return(NULL)})

        # Check if optimization was successful
        if (!is.null(MLE_1) && is.finite(MLE_1$objective)) {
          param_H1_test = c(MLE_1$par,-MLE_1$objective)
          break                                 } else { ## If optimization failed, increment x0_1 slightly and retry
            x0_1 <- x0_1 + c(0.1,0.1,0.1)
            retry_count <- retry_count + 1

            if (retry_count >= max_retries) {break}
          }
      }

      if (retry_count >= max_retries) {next;}
      wilks_test = 2*(param_H1_test[length(param_H1_test)]-param_H0_test[length(param_H0_test)])
      if(wilks_test<0) {next} else{
        power= ifelse(wilks_test>= qchisq(0.95, df=1), 1, 0)
        summary_test[trial,]= c(param_H0_test, param_H1_test, wilks_test, power,m)
        trial=trial+1
        rm(param_H0_test); rm(param_H1_test)}
}

colnames(summary_test) = c("Gamma_0", "H_0[1]", "LogL_0" ,"Gamma_1","H_1[1]", "H_1[2]","LogL_1","Wilks", "sign", "m")

################################ Results #####################################

cat("Type I Error:",nrow(summary[summary$sign<=alpha,])*100/nrow(summary),"%")  ## Should be 5%
cat("Power:",nrow(summary[summary_test$sign==1,])*100/nrow(summary_test),"%")  ## Should be 5%

final[kk,] = c(thetas, nrow(summary[summary$sign<alpha,])*100/nrow(summary), nrow(summary[summary_test$sign==1,])*100/nrow(summary_test),  mean(summary[,"m"]))
kk=kk+1
}
## Gamma 1.2, scale:2, shape =0.8
#Type I Error: 7.7 %
#Power: 87.6 %

## Gamma 1.2, scale:1, shape =0.8
#Type I Error: 7 %
#Power: 23%

## Gamma 1.2, scale:3, shape =0.8
#Type I Error: 32 %
#Power: 94%

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
