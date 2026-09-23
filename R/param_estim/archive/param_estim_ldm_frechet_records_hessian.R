library(stats)
library(GA)
library(pracma)
library(maxLik)
#library("xlsx")
library(VGAM)
options(scipen = 6)

Plot_auto = function(variable_name, matrice = sum_matrices, file_name = "plot_bias.png"){

  # Extract values from sum_matrices
  gamma_vals <- as.numeric(names(matrice))
  gamma_bias <- sapply(matrice, function(x) x[variable_name, "thetaHat"])
  A_bias <- sapply(matrice, function(x) x[variable_name, "scaleHat"])
  alpha_bias <- sapply(matrice, function(x) x[variable_name, "shapeHat"])

  # Save the plot as an image with a unique name
  png(file_name, width = 800, height = 600, res = 150)

  # Create the plot
  ylim_range <- range(c(gamma_bias, A_bias, alpha_bias), na.rm = TRUE)
  plot(gamma_vals, gamma_bias, type = "o", col = "plum", xlab = "Theta", ylab = variable_name, ylim = ylim_range)
  lines(gamma_vals, A_bias, type = "o", col = "cyan4")
  lines(gamma_vals, alpha_bias, type = "o", col = "lightgreen")

  # Add a legend (uncomment if needed)
  # legend("topright", legend = c("Power", "1/Scale", "Shape"), col = c("plum", "cyan4", "lightgreen"), lty = 1, pch = 1, cex = 0.7)

  dev.off()  # Close the graphics device

  return(file_name)
}

Plot_auto_var = function(matrice = sum_matrices, file_name = "plot_variance.png"){

  # Extract empirical and theoretical variances for each parameter
  Emp_gamma_var <- sapply(matrice, function(x) x["AVG_Emp_std", "thetaHat"])
  Emp_A_var <- sapply(matrice, function(x) x["AVG_Emp_std", "scaleHat"])
  Emp_alpha_var <- sapply(matrice, function(x) x["AVG_Emp_std", "shapeHat"])

  Theo_gamma_var <- sapply(matrice, function(x) x["AVG_Theo_std", "thetaHat"])
  Theo_A_var <- sapply(matrice, function(x) x["AVG_Theo_std", "scaleHat"])
  Theo_alpha_var <- sapply(matrice, function(x) x["AVG_Theo_std", "shapeHat"])

  # Save the plot as an image with a unique name
  png(file_name, width = 800, height = 600, res = 150)

  # Automatically determine the y-axis limits
  ylim_range <- range(c(Emp_gamma_var, Emp_A_var, Emp_alpha_var,
                        Theo_gamma_var, Theo_A_var, Theo_alpha_var), na.rm = TRUE)
  #X-axis
  gamma_vals <- as.numeric(names(matrice))

  # Plot Empirical vs Theoretical Variance for gammaHat
  plot(gamma_vals, Emp_gamma_var, type = "o", col = "plum", xlab = "Theta", ylab = "Std",
       main = "Empirical vs Theoretical Std", ylim = ylim_range)
  lines(gamma_vals, Theo_gamma_var, type = "o", col = "plum", lty = 2)

  # Add Empirical and Theoretical Variance for scaleHat
  lines(gamma_vals, Emp_A_var, type = "o", col = "cyan4")
  lines(gamma_vals, Theo_A_var, type = "o", col = "cyan4", lty = 2)

  # Add Empirical and Theoretical Variance for shapeHat
  lines(gamma_vals, Emp_alpha_var, type = "o", col = "lightgreen")
  lines(gamma_vals, Theo_alpha_var, type = "o", col = "lightgreen", lty = 2)

  dev.off()  # Close the graphics device

  return(file_name)
}

list_to_df= function(sum_matrices, var="T"){
  df <- data.frame(matrix(unlist(sum_matrices), nrow=length(sum_matrices), byrow=TRUE))
  rownames(df) = paste0(var,"=",names(sum_matrices))
  colnames(df) =
    with(expand.grid(rownames(sum_matrices[[1]]), colnames(sum_matrices[[1]])), paste0(Var1,"_",Var2))
  return(df)
}

###########################   Parameters and Dataframes ###########################
simulation=100

## Save Results in workbook
save = FALSE
save_plot = FALSE
save_details = FALSE

## Frechet parameters
scale_obs = 1
scale_inv_obs=1/scale_obs
shape_obs=2  ## shape

## Series length
#T_val=100
T_values <- 50# c(50, 100, 200)

# Define the sequence of theta values
gamma_values <- 0.25# seq(0.05, 0.2, by = 0.05)

## significance level
sign = 0.05

## Likelihood expression
info = "records"

## file name if saving
file_name = paste0("Simulation_ldm_Rn_overT_frechet_shape=",shape_obs,"_scale=", scale_obs, ".xlsx")

###########################  Simulation  ###########################

# Load the openxlsx package for writing Excel files
library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

# Create a list to store sum_matrix for each gamma
sum_matrices <- list()


results_T <- list()

for (T_val in T_values) {

  sum_matrices <- list()

  # Loop over each gamma value
  for (gamma in gamma_values) {

  ## Reset the matrices at the start of each iteration
  params <- matrix(0, nrow = simulation, ncol = 4)
  colnames(params) <- c("gammaHat", "scaleHat", "shapeHat", "N_T")

  LogL <- matrix(0, nrow = simulation, ncol = 1)

  Emp_var <- matrix(0, nrow = simulation, ncol = 3)
  colnames(Emp_var) <- c("gamma", "scale", "shape")

  Theo_var <- matrix(0, nrow = simulation, ncol = 3)
  colnames(Theo_var) <- c("gamma", "scale", "shape")

  sum_matrix =  matrix(0,nrow=5,ncol=4)
  rownames(sum_matrix)= c("AVG_param", "AVG_bias","AVG_Emp_std", "AVG_Theo_std", "Coverage_proba")
  colnames(sum_matrix) = c("thetaHat","scaleHat","shapeHat","N_T")

  conf_int <- matrix(0, nrow = simulation, ncol = 6)
  colnames(conf_int) <- c("g_L", "g_U", "A_L", "A_U", "a_L", "a_U")

  ############################ Simulation Code #############################

  no=1:simulation

  while(length(no)>=1){ #while (sim <= simulation){
    for(sim in no){

      ## generate LDM series of Frechet underlying distr
      y= LDM_series(T = T_val,theta = gamma, dist = "frechet",
                    shape = shape_obs, scale=scale_obs)
      R=rec_values(y)
      L = rec_times(y)

      while(length(R)<=1) {
        y= LDM_series(T = T_val,theta = gamma, dist = "frechet",
                      shape = shape_obs, scale=scale_obs)
        R=rec_values(y)
        L = rec_times(y)
      }  ## only one record, ignore

      if(info == "all"){
        logLik_fun_rec <- loglik_registry[["LDM"]][["all"]][["frechet"]]
        MLE_C = estimate_model_mle(logLik_fun_rec, obs_type = "all", data = y,
                                   lower_bounds=c(theta=0.01, shape = 0.01, scale=0.01),
                                   upper_bounds = c(theta=10, shape = 10, scale=10),
                                   start_values = c(theta=0.01, shape = 0.01, scale=0.01)
        )
        MLE_C = c(MLE_C$par, "objective" = - MLE_C$objective )
        #MLE_C = c(-max_est, max_est, - logLik_fun_rec(data = y, params = c(max = max_est)) )

      } else{
        ## Records Data
        logLik_fun_rec <- loglik_registry[["LDM"]][["records"]][["frechet_inv_scale_explicit"]]

        data_rec = list(rec_values = R, rec_times = L, time = T_val)
        data_rec = data.frame(rec_values = R, rec_times = L, time = T_val)
        MLE_C = estimate_model_mle(logLik_fun_rec, obs_type = "records", data = data_rec,
                                   lower_bounds=c(theta=0.01, shape = 0.1, scale=0.05),
                                   upper_bounds = c(theta= min(R/L)-1e-10, shape = 10, scale=10),
                                   start_values = c(
                                     theta = 0.01,
                                     shape = 0.01,
                                     scale = sd(R - median(R / L) * L)
                                   )
        )
        MLE_C = c(MLE_C$par, "objective" = - MLE_C$objective ) ## return location, shape, scale, objective value
      MLE_C
        }

      ## store estimated parameters
      params[sim,"gammaHat"] <- MLE_C["theta"]
      params[sim,"scaleHat"] <- MLE_C["scale"]
      params[sim,"shapeHat"] <- MLE_C["shape"]
      params[sim,"N_T"] <- rec_count(y)
      LogL[sim,1] = MLE_C["objective"]

      ## ____________________ Variance estimation (HESSIAN) ____________________

      # --- Empirical variance (at MLE) ---
      par_hat <- c(
        MLE_C["theta"],
        MLE_C["shape"],
        MLE_C["scale"]
      )

      var_emp <- compute_var_hessian(logLik_fun_rec, data_rec, par_hat,eps = 1e-6)

      Emp_var[sim, "gamma"] <- var_emp[1]
      Emp_var[sim, "shape"] <- var_emp[3]
      Emp_var[sim, "scale"] <- var_emp[2]


      # --- Theoretical variance (at TRUE params) ---
      par_true <- c(
        theta = gamma,
        shape = shape_obs,
        scale = scale_obs
      )

      var_theo <- compute_var_hessian(logLik_fun_rec, data_rec, par_true)

      Theo_var[sim, "gamma"] <- var_theo[1]
      Theo_var[sim, "shape"] <- var_theo[3]
      Theo_var[sim, "scale"] <- var_theo[2]
    }

    no = which(
      rowSums(is.na(Emp_var)) > 0 |
        rowSums(is.na(Theo_var)) > 0 |
        rowSums(Emp_var <= 0) > 0 |
        rowSums(Theo_var <= 0) > 0
    )
  }

  ##___________________________________________________________

  ## Compute the summary matrix
  sum_matrix[1, ] <- round(colMeans(params), 3)  # Average of parameters
  sum_matrix[2, ] <- round(colMeans(params) - c(gamma, scale_obs, shape_obs, 0), 5)  # Average bias
  sum_matrix[3, ] <- round(c(sqrt(colMeans(Emp_var)), 0), 5)  # Average Empirical variance
  sum_matrix[4, ] <- round(c(sqrt(colMeans(Theo_var)), 0), 5)  # Average asymptotic variance

  ## Compute coverage probability
  z <- qnorm(1 - (sign / 2), 0, 1)

  power_g <- 0
  power_A <- 0
  power_a <- 0

  for (sim in 1:simulation) {
    ## CI gamma
    bound_g = bounds(  params[sim,"gammaHat"] ,z=z,Theo_var[sim,"gamma"] )

    if (gamma < bound_g[2] && gamma > bound_g[1]) {
      power_g <- power_g + 1 }

    ## CI A
    bound_A = bounds(  params[sim,"scaleHat"] ,z=z,Theo_var[sim,"scale"] )

    if (scale_obs < bound_A[2] && scale_obs > bound_A[1]) {
      power_A <- power_A + 1 }

    ##CI a
    bound_a = bounds(  params[sim,"shapeHat"] ,z=z,Theo_var[sim,"shape"] )

    if (shape_obs < bound_a[2] && shape_obs > bound_a[1]) {
      power_a <- power_a + 1 }

    ## fill into matrix
    conf_int[sim,]= c(bound_g,bound_A,bound_a)
  }

  CP_gamma <- power_g / simulation
  CP_A <- power_A / simulation
  CP_alpha <- power_a / simulation
  sum_matrix[5, ] <- c(CP_gamma, CP_A, CP_alpha, 0)

  ## Store the sum_matrix for this gamma
  sum_matrices[[as.character(gamma)]] <- sum_matrix


  }

  results_T[[paste0("T=", T_val)]] <- sum_matrices
}

  ## Saving
  if(save_details){
  # Add each matrix to a new sheet
  addWorksheet(wb, paste0("Sum_", gamma))
  writeData(wb, sheet = paste0("Sum_", gamma), sum_matrix)

  addWorksheet(wb, paste0("TheoVar_", gamma))
  writeData(wb, sheet = paste0("TheoVar_", gamma), Theo_var)

  addWorksheet(wb, paste0("EmpVar_", gamma))
  writeData(wb, sheet = paste0("EmpVar_", gamma), Emp_var)

  addWorksheet(wb, paste0("params_", gamma))
  writeData(wb, sheet = paste0("params_", gamma), params)

  addWorksheet(wb, paste0("LogL_", gamma))
  writeData(wb, sheet = paste0("LogL_", gamma), LogL)
}

# Check the results
print(sum_matrices)

#Add to EXCEL
if(save){
addWorksheet(wb, c("All Par"))
writeData(wb, c("All Par"), list_to_df(sum_matrices, var="theta"), rowNames = TRUE , colNames = TRUE)
}
# Plots -------------------------------------------------------------------

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

## Save the workbook--------------------
if(save){
saveWorkbook(wb, paste0("Simulation_LDM_Rn_T",T,"_A=",A,",a=",shape,".xlsx"), overwrite = TRUE)
}
