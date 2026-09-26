library("xlsx")
library("ggplot2")
library(openxlsx)
devtools::load_all(".")

# library(devtools)
# install_github("jj-64/Records")
# library(Records)

n_sim <- 1000
T <- seq(40, 100, by = 10)
save = TRUE
save_path ="C:/Users/User/OneDrive - Lebanese University/Lebanese University/PhD/Publication 3 - Decision Tree test/Power of Test/"

# -------------------------------
# Generic Simulation Function
 -------------------------------

simulate_model <- function(param_values, ## vector of values of the parameter that we are simulating upon
                           param_name, ## string: the parameter we are simulation upon "gamma", "theta", "scale"...
                           T, ## integer: length of the simulated series
                           n_sim,  ## integer: number of simulations
                           generator, ## function:the function generating the series
                           series_args=list(), ## arguments of the generator function other than "T" and the "param_name" we are simulating
                           test_fun, ## function: test function
                           test_args=list(),
                           n_arg="T"){ # could be "T" or "n"

  results <- as.data.frame(matrix(0, nrow=length(param_values), ncol=length(T)+2))
  results[,1] <- param_values
  colnames(results) <- c(param_name, paste0("T_", T), "average")

  for (k in seq_along(T)) {
    for (j in seq_along(param_values)) {
      dec <- rep(NA, n_sim)
      for (i in 1:n_sim) {
        # Build generator args
        args <- series_args
        args[[param_name]] <- param_values[j]
        args[[n_arg]] <- T[k]   # could be "T" or "n"

        x <- do.call(generator, args)

        # Apply test function
        test_call <- c(list(X=x), test_args)
        dec[i] <- do.call(test_fun, test_call)$decision
      }
      valid <- na.omit(dec)
      results[j, k+1] <- mean(valid=="NO") * 100
    }
  }

  results[,"average"] <- rowMeans(results[,2:(length(T)+1)])
  return(results)
}

# -------------------------------
# Plotting Function
# -------------------------------

plot_results <- function(df, param_name, title, ylab_name = "Power of test (1-ß, %)", xlab_name = NULL, ymin=0, ymax=100 ) {
  if (is.null(xlab_name)) xlab_name <- param_name

  # detect T_ columns
  T_cols <- grep("^T_", names(df), value = TRUE)

  # reshape to long format
  df_long <- reshape2::melt(df,
                            id.vars = c(param_name, "average"),
                            measure.vars = T_cols,
                            variable.name = "T",
                            value.name = "Power")


  # clean legend labels (remove "T_")
  df_long$T <- as.numeric(gsub("^T_", "", df_long$T))

  # nice color scale (one color per T)
  n_T <- length(T_cols)
  #colors <- scales::hue_pal()(n_T)   # dynamic palette
  #colors <- RColorBrewer::brewer.pal(min(8, n_T), "Dark2")
  #colors <- rev(viridisLite::viridis(n_T))


  p <- ggplot(df_long, aes(x = .data[[param_name]], y = Power, color = T, group = T)) +
    geom_line(linewidth = 1.1, alpha = 0.8) +
    geom_line(aes(y = average), df_long, color = "black", linewidth = 1.2, linetype = "dashed") +
    #scale_color_manual(values = colors, name = "Sample size") +
    scale_color_viridis_c(
      option = "D", direction = -1, name = "Sample size (T)",
      breaks = seq(40, 100, by = 10),   # show fewer ticks
      labels = seq(40, 100, by = 10)
    ) +
    ylab(ylab_name) + xlab(xlab_name) +
    ggtitle(title) +
    scale_x_continuous(n.breaks = 10)+
    scale_y_continuous(n.breaks = 10, limits = c(ymin, ymax))+
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5, size=16),
      axis.title = element_text(face = "bold"),
      #legend.title = element_text(face = "bold"),
      panel.grid.major = element_line(color = "grey95"),
      panel.grid.minor = element_blank(),
      # plot.background = element_rect(fill = "#f7f7f7", color = NA),
      # panel.background = element_rect(fill = "#f7f7f7", color = NA)
    )

  return(p)
}

save_plot <- function(path, filename, plot){
  ggsave(path = path, filename = filename, plot=plot, width = 7, height = 5, dpi = 600)
}

# -------------------------------
# Excel Writer
# -------------------------------
save_results <- function(df, file, sheet) {
  xlsx::write.xlsx(df, file, sheetName=sheet, append=TRUE, row.names=FALSE)
}

save_results_with_plot <- function(df, file = "results.xlsx", sheet ,p) {
  # 1. Create workbook
  wb <- createWorkbook()

  # 2. Add worksheet
  addWorksheet(wb, sheet)

  # 3. Write dataframe
  writeData(wb, sheet, df, startRow = 1, startCol = 1)

  # 4. Generate plot
  #p <- plot_results(df, param_name, title)

  # 5. Save plot temporarily as image
  img_file <- tempfile(fileext = ".png")
  ggsave(img_file, p, width = 7, height = 5, dpi = 600)

  # 6. Insert image into worksheet (e.g. below the dataframe)
  insertImage(wb, sheet, img_file,
              startRow = nrow(df) + 3, startCol = 1,
              width = 7, height = 5)

  # 7. Save Excel file
  saveWorkbook(wb, file, overwrite = TRUE)
}

######################## 1️⃣ Part 1 -  Classical Model 1️⃣ ##############################
# -------------------------------
# Run H0: Classical vs H1: Yang
# -------------------------------
gamma <- seq(1.01, 1.4, by=0.01)
m_c_y <- simulate_model(param_values = gamma,
  T = T,
  n_sim = n_sim,
  generator = YNM_series,
  param_name = "gamma", # varying param
  n_arg = "T",             # custom generator expects T=
  test_fun = Test_iid_BoxJenkins,#Test_iid_NT,
  series_args = list(dist = "gumbel", loc=0, scale=1)
)
plot_results(m_c_y, param_name="gamma", title="Classical vs Yang - Gumbel", xlab_name = "γ")
if(save == TRUE) {save_results(m_c_y, paste0(save_path,"/Classical BoxJenkins.xlsx"), "YNM_Gumbel_0_1")}

# -------------------------------
# Run H0: Classical vs H1: LDM
# -------------------------------
theta_vals <- seq(0.01, 0.5, by=0.05)
m_c_L <- simulate_model(
  param_values = theta_vals,
  T = T,
  n_sim = n_sim,
  generator = LDM_series,
  param_name = "theta",
  n_arg = "T",
  test_fun = Test_iid_BoxJenkins,#Test_iid_NT,
  series_args = list(dist="frechet",shape=5, scale=5)
)
plot_results(m_c_L, param_name="theta",  title = "Classical vs LDM - Frechet", xlab_name = "Θ")
if(save == TRUE) {save_results(m_c_L, paste0(save_path,"/Classical BoxJenkins.xlsx"), "LDM_Frechet_5_1")}

# -------------------------------
# Run H0: Classical vs H1: DTRW
# -------------------------------
scale_vals <- seq(1, 5, by=0.5)
m_c_R <- simulate_model(
  param_values = scale_vals,
  T = T,
  n_sim = n_sim,
  generator = DTRW_series,
  param_name = "scale",  ## sd for norm, scale for cauchy
  n_arg = "T",
  test_fun = Test_iid_BoxJenkins,#Test_iid_NT,
  series_args = list(dist="cauchy",loc=0)
)
plot_results(m_c_R, param_name="scale", title="Classical vs DTRW - Cauchy", xlab_name = "σ")
if(save == TRUE) {save_results(m_c_R, paste0(save_path,"/Classical BoxJenkins.xlsx"), "DTRW_Cauchy")}
# v1=NT_DTRW(0:10, 10)*100
# v=NA
# for(i in 0:(length(v1)-1)) v[i+1] = 100*NT_iid(i, (length(v1)-1))
# plot(0:(length(v1)-1), y=v1, type = "l", xlab = "X", ylim=c(0,30))
# lines(0:(length(v1)-1), y=v, type = "l", col = "red")

# -------------------------------
# Detection Rate: Classical vs Classical
# -------------------------------
scale_vals <- seq(1, 5, by=0.5)
m_c_c <- simulate_model(
  param_values = scale_vals,
  T = T,
  n_sim = n_sim,
  generator = rnorm,
  param_name = "sd",  ## sd for norm
  n_arg = "n",
  test_fun = Test_iid_BoxJenkins,#Test_iid_NT,
  series_args = list(mean=0)
)
plot_results(m_c_c, param_name="sd", title="Classical Detection Rate", xlab_name = "σ", ymax=25)
if(save == TRUE) {save_results(m_c_c, paste0(save_path,"/Classical BoxJenkins.xlsx"), "Detection")}

############################ 2️⃣ PART 2 : LDM 2️⃣  ####################################################
# ----------------------------------
# Run: H0: LDM vs H1: Yang
# ----------------------------------
gamma <- c(1.01,seq(1.05, 1.4, by=0.05))
m_L_y <- simulate_model(
  param_values = gamma,
  T = T,
  n_sim = n_sim,
  generator = YNM_series,
  param_name = "gamma",
  n_arg = "T",
  test_fun = Test_LDM_Sequential,
  series_args = list(dist="norm",loc=0, scale=1)
)
plot_results(m_L_y, "gamma", "LDM vs Yang-Nevzorov - Weibull", xlab_name="Gamma (γ)", ymax=100)
if(save == TRUE) {save_results(m_L_y, paste0(save_path,"/LDM_Sequential.xlsx"), "YNM_Weibull_5_1")}

# ----------------------------------
# Run H0: LDM vs H1: Classical
# ----------------------------------
b <- sqrt(seq(1, 5,1))
m_L_c <- simulate_model(
  param_values = b,
  T = T,
  n_sim = n_sim,
  generator = rnorm, #VGAM::rgumbel,
  param_name = "sd",   # param goes into scale=
  n_arg = "n",           # rgumbel expects n=
  test_fun =  Test_LDM_Sequential,
  series_args = list(mean=0)
)
plot_results(m_L_c, "sd", title="LDM vs Classical", xlab_name="scale parameter for normal")
if(save == TRUE) {save_results(m_L_c, paste0(save_path,"/LDM_Sequential.xlsx"), "Classical_Norm")}

# ----------------------------------
# H0: LDM vs H1: DTRW
# ----------------------------------
scale_vals <- seq(1, 5, 1)
m_L_R <- simulate_model(
  param_values = scale_vals,
  T = T,
  n_sim = n_sim,
  generator = DTRW_series,
  param_name =  "scale",      # sd for DTRW and scale for Cauchy
  n_arg = "T",
  test_fun = Test_LDM_Sequential,
  series_args = list(dist="cauchy",loc=0)
              )
plot_results(m_L_R, "scale", title= "LDM vs DTRW - Cauchy", xlab_name="Scale (σ²)")
if(save == TRUE) {save_results(m_L_R, paste0(save_path,"/LDM_Sequential.xlsx"), "DTRW_Norm")}

# ----------------------------------
# H0: LDM vs H1: LDM (should be Low)
# ----------------------------------
theta_vals <-c(0.01,seq(0.05,0.3,0.05))
m_L_L <- simulate_model(
  param_values = theta_vals,
  T = T,
  n_sim = n_sim,
  generator =LDM_series,
  param_name = "theta",      # sd for DTRW and scale for Cauchy
  n_arg = "T",
  test_fun = Test_LDM_Regression,
  series_args = list(dist= "frechet", shape=5, scale=1)
)
plot_results(m_L_L, "theta", "LDM vs LDM", xlab_name=" Theta (Θ) ", ymax= 25)
if(save == TRUE) {save_results(m_L_L, paste0(save_path,"/LDM_Sequential.xlsx"), "Detection")}

############################ 3️⃣ PART 3: DTRW 3️⃣ ################################################

# ----------------------------------
# H0: DTRW vs H1: Yang
# ----------------------------------
gamma <- seq(1.01, 1.4, by=0.1)
m_R_y <- simulate_model(
  param_values = gamma,
  T = T,
  n_sim = n_sim,
  generator = YNM_series,
  param_name = "gamma",
  n_arg = "T",
  test_fun = Test_DTRW_Indep, #Test_DTRW_bonf,
  series_args = list(dist="frechet",shape=5, scale=1)
)
plot_results(m_R_y, "gamma", "DTRW vs YNM - Gumbel", xlab = "Gamma (γ)")
if(save == TRUE) {save_results(m_R_y, paste0(save_path,"/DTRW_Indep.xlsx"), "YNM_Frechet_5_1")}

# ----------------------------------
# H0: DTRW vs H1: Classical
# ----------------------------------
b <- seq(1, 5, 1)
m_R_c <- simulate_model(
  param_values = b,
  T = T,
  n_sim = n_sim,
  generator = rnorm,
  param_name= "sd",
  n_arg = "n",
  test_fun = Test_DTRW_Indep,
  #test_args = list(alpha=0.05, method= "Bonf"),
  series_args = list(mean=0)
)
plot_results(m_R_c, "sd", "DTRW vs Normal i.i.d", xlab= "σ")
if(save == TRUE) {save_results(m_R_c,paste0(save_path,"/DTRW_Indep.xlsx"), "Classical_Norm")}

# ----------------------------------
# DTRW vs LDM
# ----------------------------------
theta_vals <- seq(0.01, 0.1, by=0.02)
m_R_L <- simulate_model(
  param_values = theta_vals,
  T = T,
  n_sim = n_sim,
  generator = LDM_series,
  param_name = "theta",
  n_arg = "T",
  test_fun = Test_DTRW_Indep,
  #test_args = list(method="Bonf"),
  series_args = list(dist="weibull", shape=1, scale=1)
)
plot_results(m_R_L, "theta", "DTRW vs LDM - Normal", xlab="theta (Θ)")
if(save == TRUE) {save_results(m_R_L, paste0(save_path,"/DTRW_Indep.xlsx"), "LDM_Weibull_2_1")}

# ----------------------------------
# DTRW Detection Rate: should be low (5%)
# ----------------------------------
scale_vals <- seq(1, 2, by=1)
m_R_R <- simulate_model(
  param_values = scale_vals,
  T = T,
  n_sim = n_sim,
  generator = DTRW_series,
  param_name = "sd",
  n_arg = "T",
  test_fun = Test_DTRW_Indep,
  #test_args = list(method="Bonf"),
  series_args = list(dist="norm",loc=0)
)
plot_results(m_R_R, "sd", "Detetction", ylab_name = "Type I Error (%)",xlab_name="Scale (σ²)", ymax = 25)
if(save == TRUE) {save_results(m_R_R, paste0(save_path,"/DTRW_Indep.xlsx"), "Detection")}

########################## 4️⃣ Part 4 - YANG 4️⃣ #######################################
# ----------------------------------
#  Yang vs Classiacl
# ----------------------------------
b <- seq(1, 2, 1)
m_y_c <- simulate_model(
  param_values = b,
  T = T,
  n_sim = n_sim,
  generator = rnorm,       # built-in uniform
  param_name = "sd",
  n_arg = "n",
  test_fun = Test_YNM_Pearson,
  series_args = list(mean=0),
  #test_args = list(K=NULL, warmup=2)#list(alpha=0.05, Partition=NA, gamma=1, estimated=1)
)
p = plot_results(m_y_c, "sd", "Yang vs Classical", xlab_name = "Scale (σ²)")
if(save == TRUE) {save_results(m_y_c, paste0(save_path,"/YNM_Pearson.xlsx"), "Classical_Norm")
save_plot(path =paste0(save_path, "Figures"), filename = "YNM_Pearson vs Classical_Norm.png" , plot = p)}

# ----------------------------------
# Yang vs DTRW (Normal increments)
# ----------------------------------
scale <- seq(1, 2, 1)
m_y_R <- simulate_model(
  param_values = scale,
  T = T,
  n_sim = n_sim,
  generator = DTRW_series,
  param_name = "sd",
  n_arg = "T",
  test_fun = Test_YNM_Pearson,
  series_args = list(dist="norm",loc=0),
  #test_args = list(K=NULL, warmup=NULL) #list(alpha=0.05, Partition=NA,gamma=1, estimated=1)
)
p=plot_results(m_y_R, "sd", "", xlab_name = "Scale (σ²)")

if(save == TRUE) {
  save_results(m_y_R, paste0(save_path,"/YNM_Pearson.xlsx"), "DTRW_Norm")
  save_plot(path =paste0(save_path, "Figures"), filename = "YNM_Pearson vs DTRW_Norm.png" , plot = p)
}
# ----------------------------------
# Yang vs LDM (Frechet)
# ----------------------------------
theta_vals <- c(0.01,seq(0.05, 0.2, 0.05))
m_y_L <- simulate_model(
  param_values = theta_vals,
  T = T,
  n_sim = n_sim,
  generator = LDM_series,
  param_name = "theta",
  n_arg = "T",
  test_fun = Test_YNM_Geom,
  series_args = list(dist="frechet",shape=5, scale=1),
  #test_args = list(K= 4)
)
p= plot_results(m_y_L, "theta", "", xlab = "Theta (Θ)")
if(save == TRUE) {save_results(m_y_L, paste0(save_path,"/YNM_Pearson.xlsx"), "LDM-Frechet_5_1")
  save_plot(path =paste0(save_path, "Figures"), filename = "YNM_Pearson vs LDM_Frechet_5_1.png" , plot = p)
}

for(i in 1:1000){
X= LDM_series(100,0.3,"frechet", shape=5, scale=1)
rr[i] = Test_YNM_Pearson(X)$decision}
table(rr)
# ----------------------------------
#  Yang vs Yang (Detection Rate)
# ----------------------------------
gamma <- seq(1.01, 1.5, by=0.05)
m_y_y <- simulate_model(param_values = gamma,
                        T = T,
                        n_sim = n_sim,
                        generator = YNM_series,
                        param_name = "gamma", # varying param
                        n_arg = "T",             # custom generator expects T=
                        test_fun = Test_YNM_Pearson,
                        series_args = list(dist = "gumbel", loc=0, scale=1),
                        test_args = list(K=4)#list(alpha=0.05, Partition=NA,gamma=1, estimated=1)
                        )
p = plot_results(m_y_y, "gamma", "", xlab="Gamma (γ)", ymax= 100)
if(save == TRUE) {save_results(m_y_y, paste0(save_path,"/YNM_Pearson.xlsx"), "Detection")
  save_plot(path =paste0(save_path, "Figures"), filename = "YNM_Pearson_typeI.png" , plot = p) }

