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
  gamma_bias <- sapply(matrice, function(x) x[variable_name, "gamma"])
  A_bias <- sapply(matrice, function(x) x[variable_name, "mean"])
  alpha_bias <- sapply(matrice, function(x) x[variable_name, "sd"])

  # Save the plot as an image with a unique name
  png(file_name, width = 800, height = 600, res = 150)

  # Create the plot
  ylim_range <- range(c(gamma_bias, A_bias, alpha_bias), na.rm = TRUE)
  plot(gamma_vals, alpha_bias, type = "o", col = "plum", xlab = xlab , ylab = variable_name, ylim = ylim_range)
  #lines(gamma_vals, A_bias, type = "o", col = "cyan4")
  #lines(gamma_vals, alpha_bias, type = "o", col = "lightgreen")

  # Add a legend (uncomment if needed)
  # legend("topright", legend = c("Power", "1/sd_obs", "Shape"), col = c("plum", "cyan4", "lightgreen"), lty = 1, pch = 1, cex = 0.7)

  dev.off()  # Close the graphics device

  return(file_name)
}

Plot_auto_var = function(matrice = sum_matrices, file_name = "plot_variance.png"){

  # Extract empirical and theoretical variances for each parameter
  Emp_gamma_var <- sapply(matrice, function(x) x["AVG_Emp_var", "gamma"])
  Emp_A_var <- sapply(matrice, function(x) x["AVG_Emp_var", "mean"])
  Emp_alpha_var <- sapply(matrice, function(x) x["AVG_Emp_var", "sd"])

  Theo_gamma_var <- sapply(matrice, function(x) x["AVG_Theo_var", "gamma"])
  Theo_A_var <- sapply(matrice, function(x) x["AVG_Theo_var", "mean"])
  Theo_alpha_var <- sapply(matrice, function(x) x["AVG_Theo_var", "sd"])

  # Save the plot as an image with a unique name
  png(file_name, width = 800, height = 600, res = 150)

  # Automatically determine the y-axis limits
  ylim_range <- range(c(Emp_gamma_var, Emp_A_var, Emp_alpha_var,
                        Theo_gamma_var, Theo_A_var, Theo_alpha_var), na.rm = TRUE)
  #X-axis
  gamma_vals <- as.numeric(names(matrice))

  # Plot Empirical vs Theoretical Variance for gamma
  plot(gamma_vals, Emp_gamma_var, type = "o", col = "plum", xlab = "Gamma", ylab = "var",
       main = "Empirical vs Theoretical var", ylim = ylim_range)
  lines(gamma_vals, Theo_gamma_var, type = "o", col = "plum", lty = 2)

  # Add Empirical and Theoretical Variance for mean
  lines(gamma_vals, Emp_A_var, type = "o", col = "cyan4")
  lines(gamma_vals, Theo_A_var, type = "o", col = "cyan4", lty = 2)

  # Add Empirical and Theoretical Variance for sd
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
# How many to simulate in one boot
simulation=1000

# number of bootstraps
B= 1000

## Save Results in workbook
save = TRUE
save_plot = FALSE
save_details = TRUE

## Series length
T_values=  c(50, 100) #seq(50,300,by =10 )

## Normal sd parameter
mean_obs = 0
sd_obs = 3
obs_vals = c(1, mean_obs, sd_obs)

# Define the sequence of theta values
gamma =1

## significance level
sign = 0.05

## Likelihood expression
info = "all"

## file name if saving
file_name = paste0("Simulation_dtrw_Xt_Rn_overT_norm_sd=",sd_obs,".xlsx")
out_path = "data/param_est/dtrw/norm/boot/"
###########################  Simulation  ###########################
# Create a new workbook
wb <- openxlsx::createWorkbook()

# Create a list to store sum_matrix
sum_matrices <- list()
sum_matrices_record <- list()

# Create a list to save detailed results
detailed_results <- list()
detailed_results_record <- list()
boot_results <- list()
boot_results_record <- list()

# Loop over each gamma value
for (T_val in T_values) {

  ## Reset the matrices at the start of each iteration
  params <- matrix(0, nrow = simulation * B, ncol = length(obs_vals)+1)
  colnames(params) <- c("gamma", "mean", "sd", "N_T")
  rownames(params) = rep(c(1:B), each = simulation)

  LogL <- matrix(0, nrow = simulation * B, ncol = 1)
  rownames(LogL) = rep(c(1:B), each = simulation)

  Emp_avg <- matrix(0, nrow = B, ncol = ncol(params))
  colnames(Emp_avg) <- colnames(params)

  Emp_var <- matrix(0, nrow = B, ncol = ncol(params))
  colnames(Emp_var) <- colnames(params)

  obs_q <- matrix(0, nrow = B, ncol =  ncol(params) * 3)
  colnames(obs_q) <- c(paste0(rep(c(paste0("q_",sign), paste0("q_",1-sign)), ncol(params)), "_", rep(colnames(params), each=2)),
  paste0("CP_", colnames(params))
                       )
  sum_matrix = matrix(0,nrow= 6,ncol= ncol(params))
  rownames(sum_matrix)= c("AVG_param", "AVG_bias","AVG_Emp_var", "Coverage_proba", "bottom_q", "upper_q")
  colnames(sum_matrix) = colnames(params)

  ## Reset the matrices at the start of each iteration
  params_record <- matrix(0, nrow = simulation * B, ncol =  ncol(params) )
  colnames(params_record) <- colnames(params)
  rownames(params_record) = rep(c(1:B), each = simulation)

  LogL_record <- matrix(0, nrow = simulation * B, ncol = 1)
  rownames(LogL_record) = rep(c(1:B), each = simulation)

  Emp_avg_record <- matrix(0, nrow = B, ncol = ncol(params))
  colnames(Emp_avg_record) <- colnames(params)

  Emp_var_record <- matrix(0, nrow = B, ncol =  ncol(params))
  colnames(Emp_var_record) <- colnames(params)

  obs_q_record <- matrix(0, nrow = B, ncol =  ncol(params) * 3)
  colnames(obs_q_record) <- c(paste0(rep(c(paste0("q_",sign), paste0("q_",1-sign)), ncol(params)), "_", rep(colnames(params), each=2)),
                              paste0("CP_", colnames(params))
  )

  sum_matrix_record = matrix(0,nrow=6,ncol= ncol(params))
  rownames(sum_matrix_record)= c("AVG_param", "AVG_bias","AVG_Emp_var", "Coverage_proba", "bottom_q", "upper_q")
  colnames(sum_matrix_record) = colnames(params)

  ############################ Simulation Code #############################

  for(boot in 0:(B-1)){

    for(sim in 1:simulation){

      ## generate yang series of Frechet underlying distr
      y = DTRW_series(T = T_val, dist = "norm", mean = 0, sd= sd_obs )
      R = rec_values(y)
      L = rec_times(y)

      while(length(R)<=1) {
        y=DTRW_series(T = T_val, dist = "norm", mean = 0, sd= sd_obs )
        R=rec_values(y)
        L = rec_times(y)
      }  ## only one record, ignore

      lower_bounds <- c(0.1)
      upper_bounds <- c(10)  # Example: Bound params[1] by min(R/L)
      start_values <- c(0.1)

      ## Likelihood optimizer
      ## All
        logLik_fun <- loglik_registry[["DTRW"]][["all"]][["norm"]]
        # MLE_C = estimate_model_mle(logLik_fun_rec, obs_type = "all", data = y,
        #                            lower_bounds=c(mean=0, sd = 0.01),
        #                            upper_bounds = c(mean=0, sd = 20),
        #                            start_values = c(mean = 0, sd = 11) )
        # MLE_C = c(MLE_C$par,MLE_C$objective ) ## return mean, sd, objective value
        scale_est = sum(diff(y)^2)/(length(y)-1)
        sd_est = sqrt(scale_est)
        MLE_C = c(mean = 0, sd = sd_est , objective = logLik_fun(data = y, params = c(mean =0, sd= sd_est)) )

        ## Records
        logLik_fun_rec <- loglik_registry[["DTRW"]][["records"]][["norm"]]
        data_rec = list(rec_values = R, rec_times = L, time = T_val)
        data_rec = data.frame(rec_values = R, rec_times = L, time = T_val)
        ## Optimize
        # MLE_C = estimate_model_mle(logLik_fun_rec, obs_type = "records", data = data_rec,
        #                            lower_bounds=c(sd = 0.1),
        #                            upper_bounds = c(sd = 20),
        #                            start_values = c(sd = 0.1) )
        # MLE_C = c(0,MLE_C$par, MLE_C$objective )
        scale_est_record = sum(diff(R)^2/diff(L))/(-1 + length(R))
        sd_est_record = sqrt(scale_est_record)
        MLE_C_record = c(mean = 0, sd = sd_est_record , objective = logLik_fun_rec(data_rec, params = c(sd= sd_est_record)) )

      ## store estimated parameters
      params[sim + boot * simulation,1] <- 1  ## Gamma always 1 in iid
      params[sim + boot * simulation,"mean"] <- MLE_C[["mean"]]
      params[sim + boot * simulation,"sd"] <- MLE_C["sd"]
      params[sim + boot * simulation,"N_T"] <- rec_count(y)
      LogL[sim + boot * simulation,1] = MLE_C[["objective"]]

      params_record[sim + boot * simulation,1] <- 1  ## Gamma always 1 in iid
      params_record[sim + boot * simulation, "mean"] <- MLE_C_record[["mean"]]
      params_record[sim + boot * simulation, "sd"] <- MLE_C_record[["sd"]]
      params_record[sim + boot * simulation,"N_T"] <- rec_count(y)
      LogL_record[sim + boot * simulation,1] = MLE_C_record[["objective"]]
    }

    # Empirical average and variance per boot
    temp_params = params[(1:simulation)+simulation*boot,]
    Emp_avg[boot+1, ] = apply( temp_params ,2,FUN = mean)
    Emp_avg_record[boot+1, ] = apply( params_record ,2,FUN = mean)
    Emp_var[boot+1, ] = apply( temp_params ,2,FUN = var)
    Emp_var_record[boot+1, ] = apply( params_record ,2,FUN = var)

    #observed quantiles : Checks if obs is >= (5%) AND <= (95%) across parameters except N_T
    temp_q = apply( temp_params ,2,FUN = function(x) return(c(quantile(x, sign), quantile(x, 1-sign)))  )
    temp_CP <- c(as.integer(obs_vals >= temp_q[1, 1:length(obs_vals)] & obs_vals <= temp_q[2, 1:length(obs_vals)]), 1)
    obs_q[boot+1, ] <- c(as.vector(temp_q), temp_CP)

    temp_q_record = apply( temp_params ,2,FUN = function(x) return(c(quantile(x, sign), quantile(x, 1-sign)))  )
    temp_CP_record <- c(as.integer(obs_vals >= temp_q_record[1, 1:length(obs_vals)] & obs_vals <= temp_q_record[2, 1:length(obs_vals)]), 1)
    obs_q_record[boot+1, ] <- c(as.vector(temp_q_record), temp_CP_record)

    }

    ## Compute the summary matrix per T
        ## All data
    sum_matrix[1, ] <- round(colMeans(params), 3)  # Average of parameters
    sum_matrix[2, ] <- round(colMeans(params) - c( 1, mean_obs, sd_obs, 0), 5)  # Average bias
    sum_matrix[3, ] <- round(c(colMeans(Emp_var)), 5)  # Average Empirical variance
    sum_matrix[4, ] <- round(c(colMeans( obs_q %>% as.data.frame() %>% dplyr::select( contains( c("CP") )  )  )), 5)  # Average Empirical variance
    sum_matrix[5, ] <- round(c(colMeans( obs_q %>% as.data.frame() %>% dplyr::select( contains( c(paste0("q_",sign)) )  )  )), 5)  # bottom q
    sum_matrix[6, ] <- round(c(colMeans( obs_q %>% as.data.frame() %>% dplyr::select( contains( c(paste0("q_",1 - sign)) )  )  )), 5)  # upper q

        ##record
    sum_matrix_record[1, ] <- round(colMeans(params_record), 3)  # Average of parameters
    sum_matrix_record[2, ] <- round(colMeans(params_record) - c( 1, mean_obs, sd_obs, 0), 5)  # Average bias
    sum_matrix_record[3, ] <- round(c(colMeans(Emp_var_record)), 5)  # Average Empirical variance
    sum_matrix_record[4, ] <- round(c(colMeans( obs_q_record %>% as.data.frame() %>% dplyr::select( contains( c("CP") )  )  )), 5)  # Average Empirical variance
    sum_matrix_record[5, ] <- round(c(colMeans( obs_q_record %>% as.data.frame() %>% dplyr::select( contains( c(paste0("q_",sign)) )  )  )), 5)  # bottom q
    sum_matrix_record[6, ] <- round(c(colMeans( obs_q_record %>% as.data.frame() %>% dplyr::select( contains( c(paste0("q_",1 - sign)) )  )  )), 5)  # upper q

    ## Store the sum_matrix for this gamma
    sum_matrices[[as.character(T_val)]] <- sum_matrix
    sum_matrices_record[[as.character(T_val)]] <- sum_matrix_record

    ## Saving
    if(save_details){
      # Add each matrix to a new sheet
      # 1. Convert to data frames and add distinct prefixes to column names
      df_q     <- as.data.frame(obs_q)
      df_avg   <- as.data.frame(Emp_avg) %>% rename_with(~paste0("avg_", .))
      df_emp    <- as.data.frame(Emp_var)  %>% rename_with(~paste0("Emp_", .))

      df_params <- as.data.frame(params)   %>% rename_with(~paste0("Param_", .))
      df_logl   <- as.data.frame(LogL)     %>% rename_with(~"LogL")

      # 2. Combine them all horizontally
      combined_data <- cbind(df_params, df_logl)
      combined_data <- combined_data %>%
        dplyr::select(-contains(c("gamma", "mean"))) %>%
        mutate(T_value = T_val, B = rep(1:B, each = simulation)) %>% # 4. Add the T_val column at the beginning (or end) of the dataframe
        relocate(B, T_value) # Moves T_value to the first column

      combined_data_q <- cbind(df_q, df_avg, df_emp) %>%
        dplyr::select(-contains(c("gamma", "mean"))) %>%
        mutate(T_value = T_val, B = 1:B) %>% # 4. Add the T_val column at the beginning (or end) of the dataframe
        relocate(B, T_value) # Moves T_value to the first column


      df_q_record   <- as.data.frame(obs_q_record)
      df_avg_record   <- as.data.frame(Emp_avg_record) %>% rename_with(~paste0("avg_", .))
      df_emp_record    <- as.data.frame(Emp_var_record)  %>% rename_with(~paste0("Emp_", .))
      df_params_record <- as.data.frame(params_record)   %>% rename_with(~paste0("Param_", .))
      df_logl_record   <- as.data.frame(LogL_record)     %>% rename_with(~"LogL")

      # 2. Combine them all horizontally
      combined_data_record <- cbind(df_params_record, df_logl_record) %>%
        dplyr::select(-contains(c("gamma", "mean"))) %>%
        mutate(T_value = T_val, B = rep(1:B, each = simulation)) %>% # 4. Add the T_val column at the beginning (or end) of the dataframe
        relocate(B, T_value) # Moves T_value to the first column

      combined_data_q_record <- cbind(df_q_record, df_avg_record, df_emp_record) %>%
        dplyr::select(-contains(c("gamma", "mean"))) %>%
        mutate(T_value = T_val, B = 1:B) %>% # 4. Add the T_val column at the beginning (or end) of the dataframe
        relocate(B, T_value) # Moves T_value to the first column

       # 5. Save this iteration's data frame into the list using T_val as the key
      detailed_results[[paste0("T_", T_val)]] <- combined_data
      detailed_results_record[[paste0("T_", T_val)]] <- combined_data_record
      boot_results[[paste0("T_", T_val)]] <- combined_data_q
      boot_results_record[[paste0("T_", T_val)]] <- combined_data_q_record

    }

    }

## Check the results
# print(sum_matrices)

if(save_details){
  final_combined_data <- bind_rows(detailed_results)
  addWorksheet(wb, "All_Results") # 3. Save to a single sheet in Excel
  writeData(wb, sheet = "All_Results", final_combined_data, rowNames = FALSE, colNames = TRUE)

  final_combined_data_record <- bind_rows(detailed_results_record)
  addWorksheet(wb, "All_Results_record") # 3. Save to a single sheet in Excel
  writeData(wb, sheet = "All_Results_record", final_combined_data_record, rowNames = FALSE, colNames = TRUE)

  boot_combined_data <- bind_rows(boot_results)
  addWorksheet(wb, "boot_Results") # 3. Save to a single sheet in Excel
  writeData(wb, sheet = "boot_Results", boot_combined_data, rowNames = FALSE, colNames = TRUE)

  boot_combined_data_record <- bind_rows(boot_results_record)
  addWorksheet(wb, "boot_Results_record") # 3. Save to a single sheet in Excel
  writeData(wb, sheet = "boot_Results_record", boot_combined_data_record, rowNames = FALSE, colNames = TRUE)


  ## summary table by Number of records
  summary_table_NT <- final_combined_data %>%
    dplyr::group_by(Param_N_T) %>%
    dplyr::summarise(
      count_sim = n(),
      Proba_N_T = n() / nrow(final_combined_data),
      across(everything(), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::relocate(count_sim, Proba_N_T, .after = last_col())
  addWorksheet(wb, "summary_NT")
  writeData(wb, sheet = "summary_NT", summary_table_NT, rowNames = FALSE, colNames = TRUE)

  summary_table_NT_record <- final_combined_data_record %>%
    dplyr::group_by(Param_N_T) %>%
    dplyr::summarise(
      count_sim = n(),
      Proba_N_T = n() / nrow(final_combined_data_record),
      across(everything(), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::relocate(count_sim, Proba_N_T, .after = last_col())
  addWorksheet(wb, "summary_NT_record")
  writeData(wb, sheet = "summary_NT_record", summary_table_NT_record, rowNames = FALSE, colNames = TRUE)
}

#Add to EXCEL
if(save) {
  # Remove any column containing "gamma"
  All_Par <- list_to_df(sum_matrices, var="T")
  All_Par <- All_Par %>%
    dplyr::select(-contains(c("gamma", "mean")))
  addWorksheet(wb, c("All Par"))
  writeData(wb, c("All Par"), All_Par, rowNames = TRUE , colNames = TRUE)

  All_Par_record <- list_to_df(sum_matrices_record, var="T")
  All_Par_record <- All_Par_record %>%
    dplyr::select(-contains(c("gamma", "mean")))
  addWorksheet(wb, c("All Par_record"))
  writeData(wb, c("All Par_record"), All_Par_record, rowNames = TRUE , colNames = TRUE)
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
  saveWorkbook(wb, paste0(out_path,file_name), overwrite = TRUE)
  print(paste0("file saved to ",out_path,file_name ))
}

