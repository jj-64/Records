library(stats)
library(GA)
library(VGAM)
library(pracma)
library(maxLik)
#library("xlsx")
library(openxlsx)
options(scipen = 6)

Plot_auto = function(variable_name, matrice = sum_matrices, file_name = "plot_bias.png", xlab = "Series length (T)"){

  # Extract values from sum_matrices
  gamma_vals <- as.numeric(names(matrice))
  gamma_bias <- sapply(matrice, function(x) x[variable_name, "gammaHat"])
  A_bias <- sapply(matrice, function(x) x[variable_name, "minHat"])
  alpha_bias <- sapply(matrice, function(x) x[variable_name, "maxHat"])

  # Save the plot as an image with a unique name
  png(file_name, width = 800, height = 600, res = 150)

  # Create the plot
  ylim_range <- range(c(gamma_bias, A_bias, alpha_bias), na.rm = TRUE)
  plot(gamma_vals, alpha_bias, type = "o", col = "plum", xlab = xlab , ylab = variable_name, ylim = ylim_range)
  #lines(gamma_vals, A_bias, type = "o", col = "cyan4")
  #lines(gamma_vals, alpha_bias, type = "o", col = "lightgreen")

  # Add a legend (uncomment if needed)
  # legend("topright", legend = c("Power", "1/Scale", "Shape"), col = c("plum", "cyan4", "lightgreen"), lty = 1, pch = 1, cex = 0.7)

  dev.off()  # Close the graphics device

  return(file_name)
}

Plot_auto_var = function(matrice = sum_matrices, file_name = "plot_variance.png"){

  # Extract empirical and theoretical variances for each parameter
  Emp_gamma_var <- sapply(matrice, function(x) x["AVG_Emp_std", "gammaHat"])
  Emp_A_var <- sapply(matrice, function(x) x["AVG_Emp_std", "minHat"])
  Emp_alpha_var <- sapply(matrice, function(x) x["AVG_Emp_std", "maxHat"])

  Theo_gamma_var <- sapply(matrice, function(x) x["AVG_Theo_std", "gammaHat"])
  Theo_A_var <- sapply(matrice, function(x) x["AVG_Theo_std", "minHat"])
  Theo_alpha_var <- sapply(matrice, function(x) x["AVG_Theo_std", "maxHat"])

  # Save the plot as an image with a unique name
  png(file_name, width = 800, height = 600, res = 150)

  # Automatically determine the y-axis limits
  ylim_range <- range(c(Emp_gamma_var, Emp_A_var, Emp_alpha_var,
                        Theo_gamma_var, Theo_A_var, Theo_alpha_var), na.rm = TRUE)
  #X-axis
  gamma_vals <- as.numeric(names(matrice))

  # Plot Empirical vs Theoretical Variance for gammaHat
  plot(gamma_vals, Emp_gamma_var, type = "o", col = "plum", xlab = "Gamma", ylab = "Std",
       main = "Empirical vs Theoretical Std", ylim = ylim_range)
  lines(gamma_vals, Theo_gamma_var, type = "o", col = "plum", lty = 2)

  # Add Empirical and Theoretical Variance for min
  lines(gamma_vals, Emp_A_var, type = "o", col = "cyan4")
  lines(gamma_vals, Theo_A_var, type = "o", col = "cyan4", lty = 2)

  # Add Empirical and Theoretical Variance for maxHat
  lines(gamma_vals, Emp_alpha_var, type = "o", col = "lightgreen")
  lines(gamma_vals, Theo_alpha_var, type = "o", col = "lightgreen", lty = 2)

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
save = TRUE
save_plot = TRUE
save_details = FALSE

## Series length
T_values=seq(25,200,by =1 )

## maximum (parameter of unfiorm)
max_obs = 3

# Define the sequence of theta values
gamma =1

## significance level
sign = 0.05

## Likelihood expression
info = "all"

## biased?
biased = TRUE ## TRUE: no bias correction; FALSE: bias correction

## file name if saving
file_name = paste0("Simulation_DTRW_Xt_overT_unif_max=",max_obs,"_biased.xlsx")
###########################  Simulation  ###########################
# Create a new workbook
wb <- openxlsx::createWorkbook()

# Create a list to store sum_matrix for each gamma
sum_matrices <- list()

# Loop over each gamma value
for (T_val in T_values) {

  ## Reset the matrices at the start of each iteration
  params <- matrix(0, nrow = simulation, ncol = 4)
  colnames(params) <- c("gammaHat", "minHat", "maxHat", "N_T")

  LogL <- matrix(0, nrow = simulation, ncol = 1)

  Emp_var <- matrix(0, nrow = simulation, ncol = 3)
  colnames(Emp_var) <- c("gamma", "min", "max")

  Theo_var <- matrix(0, nrow = simulation, ncol = 3)
  colnames(Theo_var) <- c("gamma", "min", "max")

  sum_matrix = matrix(0,nrow=5,ncol=4)
  rownames(sum_matrix)= c("AVG_param", "AVG_bias","AVG_Emp_std", "AVG_Theo_std", "Coverage_proba")
  colnames(sum_matrix) = c("gammaHat","minHat","maxHat","N_T")

  conf_int <- matrix(0, nrow = simulation, ncol = 6)
  colnames(conf_int) <- c("g_L", "g_U", "min_L", "min_U", "max_L", "max_U")

  ############################ Simulation Code #############################

  no=1:simulation

  while(length(no)>=1){ #while (sim <= simulation){
    for(sim in no){

      ## generate yang series of Frechet underlying distr
      y=dtrw_series(T =  T_val, dist = "unif", min = - max_obs, max = max_obs )
      R=rec_values(y)
      L = rec_times(y)

      while(length(R)<=1) {
        y=dtrw_series(T=T_val, dist = "unif", min = -max_obs, max = max_obs )
        R=rec_values(y)
        L = rec_times(y)
      }  ## only one record, ignore

      lower_bounds <- c(0.1)
      upper_bounds <- c(10)  # Example: Bound params[1] by min(R/L)
      start_values <- c(0.1)

      ## Likelihood optimizer
      if(info == "all"){
        logLik_fun_rec <- loglik_registry[["DTRW"]][["all"]][["unif"]]
        # MLE_C = estimate_model_mle(logLik_fun_rec, obs_type = "all", data = y,
        #                            lower_bounds=c(max = 0.01),
        #                            upper_bounds = c(max = 50),
        #                            start_values = c(max = sd((y)) ) )
        # MLE_C = c(-MLE_C$par, MLE_C$par, - MLE_C$objective ) ## return mean, max, objective value
        max_est = max(abs(diff(y)))
        if(biased == FALSE) {max_est = (T_val)/(T_val-1) * max_est}
        MLE_C = c(-max_est, max_est, - logLik_fun_rec(data = y, params = c(max = max_est)) )

      } else{

        logLik_fun_rec <- loglik_registry[["DTRW"]][["records"]][["unif"]]
        ## Records Data
        data_rec = list(rec_values = R, rec_times = L, time = T_val)
        data_rec = data.frame(rec_values = R, rec_times = L, time = T_val)
        # MLE_C = estimate_model_mle(logLik_fun_rec, obs_type = "records", data = data_rec,
        #                            lower_bounds=c(max = 0.01),
        #                            upper_bounds = c(max = 50),
        #                            start_values = c(max = sd((y)) ) )
        # MLE_C = c(-MLE_C$par, MLE_C$par, - MLE_C$objective ) ## return mean, max, objective value
        max_est = max(abs(diff(R)))
        if(biased == FALSE) {max_est = max_est * (T_val)/(T_val-1)} #bias correction
        MLE_C = c(-max_est, max_est, - logLik_fun_rec(data = data_rec, params = c(max = max_est)) )
      }

      ## store estimated parameters
      params[sim,1] <- 1  ## Gamma always 1 in iid
      params[sim,2] <- MLE_C[1]
      params[sim,3] <- MLE_C[2]
      params[sim,"N_T"] <- rec_count(y)
      LogL[sim,1] = MLE_C[3]

      ########################## Variance estimation ##################
      ### Empirical variance
      Emp_var[sim,1] = 0
      Emp_var[sim,2] = 0
      if(info == "all"){
        var_fun = var_logLik_records("dtrw", "all","unif", "max")
        Emp_var[sim,3] = var_fun(data = y, params=list("max" = max_est ), biased = biased)
      } else{
        var_fun = var_logLik_records("dtrw", "records","unif", "max")
        Emp_var[sim,3] = var_fun(data = data_rec, params=list("max" = max_est ), biased = biased)
      }

      ### Theoretical variance
      Theo_var[sim,1] =0
      Theo_var[sim,2] = 0
      if(info == "all"){
        var_fun = var_logLik_records("dtrw", "all","unif", "max")
        Theo_var[sim,3] = var_fun(data = y, params=list("max" = max_obs ), biased = biased)
      } else{
        var_fun = var_logLik_records("dtrw", "records","unif", "max")
        Theo_var[sim,3] = var_fun(data = data_rec, params=list("max" = max_obs ), biased = biased)
      }

      no=which(Theo_var[,1]<0 |Theo_var[,2]<0 |Theo_var[,3]<0  |Emp_var[,1]<0 |Emp_var[,2]<0 |Emp_var[,3]<0 | params[,1] <1)
    }

    ##########################################################################

    ## Compute the summary matrix
    sum_matrix[1, ] <- round(colMeans(params), 3)  # Average of parameters
    sum_matrix[2, ] <- round(colMeans(params) - c(1,-max_obs, max_obs, 0), 5)  # Average bias
    sum_matrix[3, ] <- round(c(sqrt(colMeans(Emp_var)), 0), 5)  # Average Empirical variance
    sum_matrix[4, ] <- round(c(sqrt(colMeans(Theo_var)), 0), 5)  # Average asymptotic variance

    ## Compute coverage probability
    z <- qnorm(1 - (sign / 2), 0, 1)

    power_sd <- 0

    for (sim in 1:simulation) {

      ##CI a
      bound_sd = bounds(  params[sim,"maxHat"] ,z=z,Theo_var[sim,"max"] )

      if (max_obs < bound_sd[2] && max_obs > bound_sd[1]) {
        power_sd <- power_sd + 1 }

      ## fill into matrix
      conf_int[sim,]= c(rep(0,4),bound_sd)
    }

    CP_sd <- power_sd / simulation
    sum_matrix[5, ] <- c(0,0, CP_sd, 0)

    ## Store the sum_matrix for this gamma
    sum_matrices[[as.character(T_val)]] <- sum_matrix

    ## Saving
    if(save_details){
      # Add each matrix to a new sheet
      addWorksheet(wb, paste0("Sum_", T_val))
      writeData(wb, sheet = paste0("Sum_", T_val), sum_matrix,rowNames = TRUE)

      addWorksheet(wb, paste0("TheoVar_", T_val))
      writeData(wb, sheet = paste0("TheoVar_", T_val), Theo_var)

      addWorksheet(wb, paste0("EmpVar_", T_val))
      writeData(wb, sheet = paste0("EmpVar_", T_val), Emp_var)

      addWorksheet(wb, paste0("params_", T_val))
      writeData(wb, sheet = paste0("params_", T_val), params)

      addWorksheet(wb, paste0("LogL_", T_val))
      writeData(wb, sheet = paste0("LogL_", T_val), LogL)
    }
  }
}

## Check the results
print(sum_matrices)

## Add to EXCEL
if(save) {
  addWorksheet(wb, c("All Par"))
  var_to_save = list_to_df(sum_matrices, var="T_val")[, grepl("maxHat", colnames(list_to_df(sum_matrices, var="T_val")))]
  var_to_save = cbind(var_to_save, AVG_param_N_T = list_to_df(sum_matrices, var="T_val")[,"AVG_param_N_T"] )
  writeData(wb, c("All Par"), var_to_save, rowNames = TRUE , colNames = TRUE)
}

# Plots -------------------------------------------------------------------

if(save_plot){
  # Bias plot vs Gamma
  p = Plot_auto(variable_name = "AVG_bias", file_name = "plot_bias.png", matrice = sum_matrices)
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
  saveWorkbook(wb, paste0("data/param_est/DTRW/unif/all/",file_name), overwrite = TRUE)
}

