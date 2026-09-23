library(stats)
library(GA)
library(VGAM)
library(pracma)
library(maxLik)
#library("xlsx")
library(openxlsx)
options(scipen = 6)

Plot_auto = function(variable_name, matrice = sum_matrices, file_name = "plot_bias.png"){

  # Extract values from sum_matrices
  gamma_vals <- as.numeric(names(matrice))
  gamma_bias <- sapply(matrice, function(x) x[variable_name, "gamma"])
  A_bias <- sapply(matrice, function(x) x[variable_name, "scale"])
  shape_bias <- sapply(matrice, function(x) x[variable_name, "shape"])

  # Save the plot as an image with a unique name
  png(file_name, width = 800, height = 600, res = 150)

  # Create the plot
  ylim_range <- range(c(gamma_bias, A_bias, shape_bias), na.rm = TRUE)
  plot(gamma_vals, gamma_bias, type = "o", col = "plum", xlab = "Gamma", ylab = variable_name, ylim = ylim_range)
  lines(gamma_vals, A_bias, type = "o", col = "cyan4")
  lines(gamma_vals, shape_bias, type = "o", col = "lightgreen")

  # Add a legend (uncomment if needed)
  # legend("topright", legend = c("Power", "1/Scale", "Shape"), col = c("plum", "cyan4", "lightgreen"), lty = 1, pch = 1, cex = 0.7)

  dev.off()  # Close the graphics device

  return(file_name)
}

Plot_auto_var = function(matrice = sum_matrices, file_name = "plot_variance.png"){

  # Extract empirical and theoretical variances for each parameter
  Emp_gamma_var <- sapply(matrice, function(x) x["AVG_Emp_std", "gamma"])
  Emp_A_var <- sapply(matrice, function(x) x["AVG_Emp_std", "scale"])
  Emp_shape_var <- sapply(matrice, function(x) x["AVG_Emp_std", "shape"])

  Theo_gamma_var <- sapply(matrice, function(x) x["AVG_Theo_std", "gamma"])
  Theo_A_var <- sapply(matrice, function(x) x["AVG_Theo_std", "scale"])
  Theo_shape_var <- sapply(matrice, function(x) x["AVG_Theo_std", "shape"])

  # Save the plot as an image with a unique name
  png(file_name, width = 800, height = 600, res = 150)

  # Automatically determine the y-axis limits
  ylim_range <- range(c(Emp_gamma_var, Emp_A_var, Emp_shape_var,
                        Theo_gamma_var, Theo_A_var, Theo_shape_var), na.rm = TRUE)
  #X-axis
  gamma_vals <- as.numeric(names(matrice))

  # Plot Empirical vs Theoretical Variance for gamma
  plot(gamma_vals, Emp_gamma_var, type = "o", col = "plum", xlab = "Gamma", ylab = "Std",
       main = "Empirical vs Theoretical Std", ylim = ylim_range)
  lines(gamma_vals, Theo_gamma_var, type = "o", col = "plum", lty = 2)

  # Add Empirical and Theoretical Variance for scale
  lines(gamma_vals, Emp_A_var, type = "o", col = "cyan4")
  lines(gamma_vals, Theo_A_var, type = "o", col = "cyan4", lty = 2)

  # Add Empirical and Theoretical Variance for shape
  lines(gamma_vals, Emp_shape_var, type = "o", col = "lightgreen")
  lines(gamma_vals, Theo_shape_var, type = "o", col = "lightgreen", lty = 2)

  dev.off()  # Close the graphics device

  return(file_name)
}

list_to_df = function(sum_matrices, var="T"){
  df <- data.frame(matrix(unlist(sum_matrices), nrow=length(sum_matrices), byrow=TRUE))
  rownames(df) = paste0(var,"=",names(sum_matrices))
  colnames(df) =
    with(expand.grid(rownames(sum_matrices[[1]]), colnames(sum_matrices[[1]])), paste0(Var1,"_",Var2))
  return(df)
}

###########################   Parameters ###########################
simulation=1000

## Save Results in workbook
save = FALSE
save_plot = FALSE
save_details = FALSE

## Frechet parameters
scale_obs = 2
scale_obs_inv=1/scale_obs  ## 1/scale
shape_obs=2  ## shape

## Series length
T_values=seq(25,200,by =50 )

# Define the sequence of theta values
gamma =1

## significance level
sign = 0.05

## Likelihood expression
info = "records"

## file name if saving
file_name = paste0("Simulation_iid_Rn_overT_frechet_shape=",shape_obs,"_scale_inv=", scale_obs, ".xlsx")
###########################  Simulation  ###########################
# Create a new workbook
wb <- createWorkbook()

# Create a list to store sum_matrix for each gamma
sum_matrices <- list()

# Loop over each gamma value
for (T_val in T_values) {

  ## Reset the matrices at the start of each iteration
  params <- matrix(0, nrow = simulation, ncol = 4)
  colnames(params) <- c("gamma", "scale", "shape", "N_T")

  LogL <- matrix(0, nrow = simulation, ncol = 1)

  Emp_var <- matrix(0, nrow = simulation, ncol = 3)
  colnames(Emp_var) <- c("gamma", "scale", "shape")

  Theo_var <- matrix(0, nrow = simulation, ncol = 3)
  colnames(Theo_var) <- c("gamma", "scale", "shape")

  sum_matrix = matrix(0,nrow=5,ncol=4)
  rownames(sum_matrix)= c("AVG_param", "AVG_bias","AVG_Emp_std", "AVG_Theo_std", "Coverage_proba")
  colnames(sum_matrix) = c("gamma","scale","shape","N_T")

  conf_int <- matrix(0, nrow = simulation, ncol = 6)
  colnames(conf_int) <- c("g_L", "g_U", "A_L", "A_U", "a_L", "a_U")

  ############################ Simulation Code #############################

  no=1:simulation

  while(length(no)>=1){ #while (sim <= simulation){
    for(sim in no){

      ## generate yang series of Frechet underlying distr
      y=VGAM::rfrechet(n=T_val,shape=shape_obs,scale = scale_obs)
      R=rec_values(y)
      L = rec_times(y)

      while(length(R)<=2) {
        y=VGAM::rfrechet(n=T_val,shape=shape_obs,scale= scale_obs)
        R=rec_values(y)
        L = rec_times(y)
        }  ## only one record, ignore

      if(info == "all"){
        logLik_fun_rec <- loglik_registry[["iid"]][["all"]][["frechet_inv"]]
        MLE_C = estimate_model_mle(logLik_fun_rec, obs_type = "all", data = y,
                           lower_bounds=c(location=0, shape = 0.01, scale=0.01),
                           upper_bounds = c(location=0, shape = 10, scale=10),
                           start_values = c(location=0, shape = 0.01, scale=0.01)
        )
        MLE_C = c(MLE_C$par, "objective" = - MLE_C$objective )
        #MLE_C = c(-max_est, max_est, - logLik_fun_rec(data = y, params = c(max = max_est)) )

      } else{
        logLik_fun_rec <- loglik_registry[["iid"]][["records"]][["frechet_inv"]]
        ## Records Data
        data_rec = list(rec_values = R, rec_times = L, time = T_val)
        data_rec = data.frame(rec_values = R, rec_times = L, time = T_val)
        MLE_C = estimate_model_mle(logLik_fun_rec, obs_type = "records", data = data_rec,
                                   lower_bounds=c(location=0, shape = 0.01, scale=0.01),
                                   upper_bounds = c(location=0, shape = 10, scale=10),
                                   start_values = c(location=0, shape = 0.01, scale=0.01)
                                   )
        MLE_C = c(MLE_C$par, "objective" = - MLE_C$objective ) ## return location, shape, scale, objective value
        #MLE_C
        }

      ## store estimated parameters
      params[sim,1] <- 1  ## Gamma always 1 in iid
      params[sim,2] <- MLE_C[["scale"]]
      params[sim,3] <- MLE_C[["shape"]]
      params[sim,"N_T"] <- rec_count(y)
      LogL[sim,1] = MLE_C[["objective"]]

      ########################## Variance estimation ##################

      ### Empirical variance
      Emp_var[sim,1] = 0
      if(info == "all"){
        # var_fun = var_logLik_records("iid", "all","frechet_inv_scale", "scale")
        # Emp_var[sim,2] = var_fun(data = y, params=list(shape = MLE_C["shape"], scale = MLE_C["scale"] )) ## var of inverted scale
        # var_fun = var_logLik_records("iid", "all","frechet_inv_scale", "shape")
        # Emp_var[sim,3] = var_fun(data = y, params=list(shape = MLE_C["shape"], scale = MLE_C["scale"] ))
        var_fun = var_logLik_records("iid", "all","frechet_inv_scale", "all")
        emp_var = var_fun(data = y, params=list(shape = MLE_C["shape"], scale = MLE_C["scale"] ))
        Emp_var[sim,"scale"] = emp_var[["scale"]]
        Emp_var[sim,"shape"] = emp_var[["shape"]]

        } else{
        # var_fun = var_logLik_records("iid", "records","frechet_inv", "scale")
        # Emp_var[sim,2] = var_fun(data = data_rec, params=list(shape = MLE_C["shape"], scale = MLE_C["scale"] )) ## var of inverted scale
        # var_fun = var_logLik_records("iid", "records","frechet", "shape")
        # Emp_var[sim,3] = var_fun(data = data_rec, params=list(shape = MLE_C["shape"], scale = MLE_C["scale"] ))
          var_fun = var_logLik_records("iid", "records","frechet_inv_scale", "all")
          emp_var = var_fun(data = data_rec, params=list(shape = MLE_C["shape"], scale = MLE_C["scale"] ))
          Emp_var[sim,"scale"] = emp_var[["scale"]]
          Emp_var[sim,"shape"] = emp_var[["shape"]]

          }

      ### Theoretical variance
      Theo_var[sim,1] =0
      if(info == "all"){
      #   var_fun = var_logLik_records("iid", "all","frechet", "scale")
      #   Theo_var[sim,2] = var_fun(data = y, params=list(shape = shape_obs, scale = scale_obs_inv )) ## var of inverted scale
      #   var_fun = var_logLik_records("iid", "all","frechet", "shape")
      #   Theo_var[sim,3] = var_fun(data = y, params=list(shape = shape_obs, scale = scale_obs_inv))
        var_fun = var_logLik_records("iid", "all","frechet_inv_scale", "all")
        theo_var = var_fun(data = y, params=list(shape = shape_obs, scale = scale_obs_inv ))
        Theo_var[sim,"scale"] = theo_var[["scale"]]
        Theo_var[sim,"shape"] = theo_var[["shape"]]

      } else{

        # var_fun = var_logLik_records("iid", "records","frechet", "scale")
        # Theo_var[sim,2] = var_fun(data = data_rec, params=list(shape = shape_obs, scale = scale_obs )) ## var of inverted scale
        # var_fun = var_logLik_records("iid", "records","frechet", "shape")
        # Theo_var[sim,3] = var_fun(data = data_rec, params=list(shape = shape_obs, scale = scale_obs))
        var_fun = var_logLik_records("iid", "records","frechet_inv_scale", "all")
        theo_var = var_fun(data = data_rec, params=list(shape = shape_obs, scale = scale_obs_inv ))
        Theo_var[sim,"scale"] = theo_var[["scale"]]
        Theo_var[sim,"shape"] = theo_var[["shape"]]
        }

    }

    no =
      which(
        Theo_var[,1] < 0 | Theo_var[,2] < 0 | Theo_var[,3] < 0 |
          Emp_var[,1]  < 0 | Emp_var[,2]  < 0 | Emp_var[,3]  < 0
      )

    # print message
    cat("Number of simulations to repeat due to negative variance:", length(no), "\n")
  }


  ##########################################################################

  ## Compute the summary matrix
  sum_matrix[1, ] <- round(colMeans(params), 3)  # Average of parameters
  sum_matrix[2, ] <- round(colMeans(params) - c(gamma, scale_obs_inv, shape_obs, 0), 5)  # Average bias
  sum_matrix[3, ] <- round(c(sqrt(colMeans(Emp_var)), 0), 5)  # Average Empirical variance
  sum_matrix[4, ] <- round(c(sqrt(colMeans(Theo_var)), 0), 5)  # Average asymptotic variance

  ## Compute coverage probability
  z <- qnorm(1 - (sign / 2), 0, 1)

  power_g <- 0
  power_A <- 0
  power_a <- 0

  for (sim in 1:simulation) {
    ## CI gamma
    bound_g = bounds(  params[sim,"gamma"] ,z=z,Theo_var[sim,"gamma"] )

    if (gamma < bound_g[2] && gamma > bound_g[1]) {
      power_g <- power_g + 1 }

    ## CI scale
    bound_A = bounds(  params[sim,"scale"] ,z=z,Theo_var[sim,"scale"] )

    if (scale_obs_inv < bound_A[2] && scale_obs_inv > bound_A[1]) {
      power_A <- power_A + 1 }

    ##CI a
    bound_a = bounds(  params[sim,"shape"] ,z=z,Theo_var[sim,"shape"] )

    if (shape_obs < bound_a[2] && shape_obs > bound_a[1]) {
      power_a <- power_a + 1 }

    ## fill into matrix
    conf_int[sim,]= c(bound_g,bound_A,bound_a)
  }

  CP_gamma <- power_g / simulation
  CP_A <- power_A / simulation
  CP_shape <- power_a / simulation
  sum_matrix[5, ] <- c(CP_gamma, CP_A, CP_shape, 0)

  ## Store the sum_matrix for this gamma
  sum_matrices[[as.character(T_val)]] <- sum_matrix

  ## Saving
  if(save_details){
  # Add each matrix to a new sheet
  # addWorksheet(wb, paste0("Sum_", T_val))
  # writeData(wb, sheet = paste0("Sum_", T_val), sum_matrix,rowNames = TRUE)
  #
  # addWorksheet(wb, paste0("TheoVar_", T_val))
  # writeData(wb, sheet = paste0("TheoVar_", T_val), Theo_var)
  #
  # addWorksheet(wb, paste0("EmpVar_", T_val))
  # writeData(wb, sheet = paste0("EmpVar_", T_val), Emp_var)
  #
  # addWorksheet(wb, paste0("params_", T_val))
  # writeData(wb, sheet = paste0("params_", T_val), params)
  #
  # addWorksheet(wb, paste0("LogL_", T_val))
  # writeData(wb, sheet = paste0("LogL_", T_val), LogL)
    # combine into one dataframe
    deatiled_results <- data.frame(
      #gamma_theo  = Theo_var[, "gamma"],
      scale_var_theo  = Theo_var[, "scale"],
      shape_var_theo  = Theo_var[, "shape"],

      #gamma_emp   = Emp_var[, "gamma"],
      scale_var_emp   = Emp_var[, "scale"],
      shape_var_emp   = Emp_var[, "shape"],

      #gamma_hat   = params[, "gamma"],
      scale_hat   = params[, "scale"],
      shape_hat   = params[, "shape"],
      N_T         = params[, "N_T"],

      LogLik      = LogL[, 1]
    )
    addWorksheet(wb, paste0("T_val_", T_val))
    writeData(wb, paste0("T_val_", T_val), deatiled_results)
  }
}

## Check the results
print(sum_matrices)

#Add to EXCEL
if(save) {
addWorksheet(wb, c("All Par"))
writeData(wb, c("All Par"), list_to_df(sum_matrices, var="T"), rowNames = TRUE , colNames = TRUE)
}

# Plots -------------------------------------------------------------------

if(save_plot){
# Bias plot vs Gamma
p = Plot_auto(variable_name = "AVG_bias", file_name = "plot_bias.png")
addWorksheet(wb, "Bias vs. gamma")
insertImage(wb, "Bias vs. gamma", p, startRow = 2, startCol = 2, width = 6, height = 4)

# Coverage Probability plot vs Gamma
p2 = Plot_auto(variable_name = "Coverage_proba", file_name = "plot_coverage.png")
addWorksheet(wb, "Cov vs. gamma")
insertImage(wb, "Cov vs. gamma", p2, startRow = 2, startCol = 2, width = 6, height = 4)

# Variance plot
p3 = Plot_auto_var(sum_matrices, file_name = "plot_variance.png")
addWorksheet(wb, "Var vs. gamma")
insertImage(wb, "Var vs. gamma", p3, startRow = 2, startCol = 2, width = 6, height = 4)
}

# Save the workbook -------------------------------------------------------
if(save){
  ## Save the workbook
  saveWorkbook(wb, paste0("data/param_est/iid/frechet/records/",file_name), overwrite = TRUE)
}
