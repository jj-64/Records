library("ggplot2")
library(openxlsx)
devtools::load_all(".")

# library(devtools)
# install_github("jj-64/Records")
# library(Records)

n_sim <- 1000
T <- seq(40, 100, by = 10)
alpha = 0.05
save = TRUE
save_path ="inst/extdata/test_power_two_by_two/"
# ______________________________________
# Generic Simulation Function ----------
# ______________________________________

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
      results[j, k+1] <- mean(valid=="no") * 100
    }
  }

  results[,"average"] <- rowMeans(results[,2:(length(T)+1)])

  results = results %>% mutate(
    se_average  =  sqrt(average * (100 - average) / n_sim),

  average_LCL = pmax(
    0,
    average - 1.96 * se_average
  ),

  average_UCL= pmin(
    100,
    average + 1.96 * se_average
  )
  )
  return(results)
}

# ______________________________________
# Plotting Function ------------
# ______________________________________

plot_results <- function(
    df,
    param_name,
    title = NULL,
    ylab_name = "Power of test (%)",
    xlab_name = NULL,
    ymin = 0,
    ymax = 100
) {

  if (is.null(xlab_name))
    xlab_name <- param_name

  ## _____________________________________________________
  ## Detect sample-size columns
  ## _____________________________________________________

  T_cols <- grep(
    "^T_[0-9]+$",
    names(df),
    value = TRUE
  )

  ## _____________________________________________________
  ## Long format
  ## _____________________________________________________

  df_long <- reshape2::melt(
    df,
    id.vars = c(param_name, "average"),
    measure.vars = T_cols,
    variable.name = "T",
    value.name = "Power"
  )

  df_long$T <- factor(
    gsub("^T_", "", df_long$T),
    levels = sort(
      unique(
        as.numeric(
          gsub("^T_", "", df_long$T)
        )
      )
    )
  )

  ## _____________________________________________________
  ## Labels at right endpoint
  ## _____________________________________________________

  label_df <- df_long |>
    dplyr::group_by(T) |>
    dplyr::slice(round(seq(
      0.2 * dplyr::n(),
      0.5 * dplyr::n(),
      length.out = 1
    ))) |>
    dplyr::ungroup()

  ## _____________________________________________________
  ## Linetypes
  ## _____________________________________________________

  linetypes <- c(
    "solid",
    "longdash",
    "dashed",
    "dotdash",
    "twodash",
    "dotted"
  )

  linetypes <- rep(
    linetypes,
    length.out = length(levels(df_long$T))
  )

  names(linetypes) <- levels(df_long$T)

  ## _____________________________________________________
  ## Base plot
  ## _____________________________________________________

  p <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(
      x = .data[[param_name]],
      y = Power,
      group = T,
      linetype = T
    )
  )

  ## _____________________________________________________
  ## Average confidence band
  ## _____________________________________________________

  if (all(c(
    "average_LCL",
    "average_UCL"
  ) %in% names(df))) {

    p <- p +
      ggplot2::geom_ribbon(
        data = df,
        ggplot2::aes(
          x = .data[[param_name]],
          ymin = average_LCL,
          ymax = average_UCL
        ),
        inherit.aes = FALSE,
        fill = "grey80",
        alpha = 0.4
      )
  }

  ## _____________________________________________________
  ## Curves
  ## _____________________________________________________

  p <- p +

    ggplot2::geom_line(
      colour = "grey35",
      linewidth = 0.8
    ) +

    ## Average curve
    ggplot2::geom_line(
      ggplot2::aes(y = average),
      colour = "black",
      linewidth = 1.6
    ) +

    ## Direct labels
    ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(
        label = T
      ),
      direction = "y",
      hjust = 0,
      nudge_x =
        0.03 *
        diff(
          range(df_long[[param_name]])
        ),
      size = 3,
      segment.color = "grey60",
      segment.size = 0.25,
      box.padding = 0.15,
      point.padding = 0,
      min.segment.length = 0,
      show.legend = FALSE
    ) +

    ggplot2::scale_linetype_manual(
      values = linetypes
    ) +

    ggplot2::scale_x_continuous(
      n.breaks = 8,
      expand = ggplot2::expansion(
        mult = c(0.02, 0.15)
      )
    ) +

    ggplot2::scale_y_continuous(
      limits = c(ymin, ymax),
      breaks = seq(
        ymin,
        ymax,
        by = 10
      )
    ) +

    ggplot2::labs(
      title = title,
      x = xlab_name,
      y = ylab_name
    ) +

    ggplot2::theme_classic(
      base_size = 13
    ) +

    ggplot2::theme(

      legend.position = "none",

      plot.title =
        ggplot2::element_text(
          hjust = 0.5,
          face = "bold",
          size = 14
        ),

      axis.title =
        ggplot2::element_text(
          face = "bold",
          size = 13
        ),

      axis.text =
        ggplot2::element_text(
          colour = "black",
          size = 11
        ),

      panel.border =
        ggplot2::element_rect(
          colour = "black",
          fill = NA,
          linewidth = 0.5
        )
    )

  return(p)
}

save_plot <- function(path, filename, plot){
  ggsave(path = path, filename = filename, plot=plot, width = 7, height = 5, dpi = 600)
}

# ______________________________________
# Excel Writer ---------------
# ______________________________________
# save_results <- function(df, file, sheet) {
#   openxlsx::write.xlsx(df, file, sheetName=sheet, append=TRUE, row.names=FALSE)
# }

save_results <- function(
    df,
    file,
    sheet,
    overwrite_sheet = TRUE
) {

  if (file.exists(file)) {

    wb <- openxlsx::loadWorkbook(file)

  } else {

    wb <- openxlsx::createWorkbook()

  }

  ## Remove sheet if already exists

  if (overwrite_sheet &&
      sheet %in% names(wb)) {

    openxlsx::removeWorksheet(
      wb,
      sheet
    )

  }

  openxlsx::addWorksheet(
    wb,
    sheet
  )

  openxlsx::writeData(
    wb,
    sheet,
    df,
    withFilter = TRUE
  )

  ## Header style

  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    fgFill = "#D9EAD3",
    halign = "center",
    border = "Bottom"
  )

  openxlsx::addStyle(
    wb,
    sheet,
    style = header_style,
    rows = 1,
    cols = 1:ncol(df),
    gridExpand = TRUE
  )

  openxlsx::freezePane(
    wb,
    sheet,
    firstRow = TRUE
  )

  openxlsx::setColWidths(
    wb,
    sheet,
    cols = 1:ncol(df),
    widths = "auto"
  )

  openxlsx::saveWorkbook(
    wb,
    file,
    overwrite = TRUE
  )

}

save_results_with_plot <- function(
    df,
    file = "results.xlsx",
    sheet,
    p,
    figure_width = 7,
    figure_height = 5,
    dpi = 600
) {

  ## ---------------------------
  ## Workbook
  ## ---------------------------

  if (file.exists(file)) {

    wb <- openxlsx::loadWorkbook(file)

  } else {

    wb <- openxlsx::createWorkbook()

  }

  ## Remove sheet if exists

  if (sheet %in% names(wb)) {

    openxlsx::removeWorksheet(
      wb,
      sheet
    )

  }

  openxlsx::addWorksheet(
    wb,
    sheet
  )

  ## ---------------------------
  ## Metadata
  ## ---------------------------

  title_style <- openxlsx::createStyle(
    textDecoration = "bold",
    fontSize = 14
  )

  openxlsx::writeData(
    wb,
    sheet,
    paste(
      "Generated:",
      Sys.time()
    ),
    startRow = 1,
    startCol = 1
  )

  ## ---------------------------
  ## Data table
  ## ---------------------------

  openxlsx::writeData(
    wb,
    sheet,
    df,
    startRow = 3,
    startCol = 1,
    withFilter = TRUE
  )

  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    fgFill = "#D9EAD3",
    border = "Bottom"
  )

  openxlsx::addStyle(
    wb,
    sheet,
    style = header_style,
    rows = 3,
    cols = 1:ncol(df),
    gridExpand = TRUE
  )

  openxlsx::freezePane(
    wb,
    sheet,
    firstActiveRow = 4
  )

  openxlsx::setColWidths(
    wb,
    sheet,
    cols = 1:ncol(df),
    widths = "auto"
  )

  ## ---------------------------
  ## Save plot
  ## ---------------------------

  img_file <- tempfile(
    fileext = ".png"
  )

  ggplot2::ggsave(
    filename = img_file,
    plot = p,
    width = figure_width,
    height = figure_height,
    dpi = dpi,
    bg = "white"
  )

  ## Position plot below data

  plot_row <-
    nrow(df) + 8

  openxlsx::insertImage(
    wb,
    sheet,
    file = img_file,
    startRow = plot_row,
    startCol = 1,
    width = figure_width,
    height = figure_height,
    units = "in"
  )

  ## ---------------------------
  ## Save workbook
  ## ---------------------------

  openxlsx::saveWorkbook(
    wb,
    file,
    overwrite = TRUE
  )

}

############################ 2️⃣ PART 2 : ldm 2️⃣  ####################################################
# ______________________________________
# Run: H0: ldm vs H1: Yang
# ______________________________________
gamma <- c(1.05,seq(1.05, 1.3, by=0.05))
m_L_y <- simulate_model(
  param_values = gamma,
  T = T,
  n_sim = n_sim,
  generator = ynm_series,
  param_name = "gamma",
  n_arg = "T",
  test_fun = test_ldm_sequential,
  series_args = list(dist="norm",location=0, scale=1),
  test_args = list(time = NA)
)
p = plot_results(m_L_y, "gamma", "ldm vs ynm - Weibull", xlab_name="Gamma (γ)", ymax=100)
p
if(save == TRUE) {
  #save_results(m_c_y, paste0(save_path,"/test_iid_serial_independence.xlsx"), "ynm_gumbel_0_1")
  save_results_with_plot( m_L_y, file = paste0(save_path,"/test_ldm_sequential_plot.xlsx"),
                          sheet = "ynm_norm_0_1",
                          p = p,
                          figure_width = 7,
                          figure_height = 5,
                          dpi = 600
  )
}
# ______________________________________
# Run H0: ldm vs H1: Classical
# ______________________________________
b <- sqrt(seq(1, 2,1))
m_L_c <- simulate_model(
  param_values = b,
  T = T,
  n_sim = n_sim,
  generator = rnorm, #VGAM::rgumbel,
  param_name = "sd",   # param goes into scale=
  n_arg = "n",           # rgumbel expects n=
  test_fun =  test_ldm_sequential,
  series_args = list(mean=0),
  test_args = list(time = NA)
)
p = plot_results(m_L_c, "sd", title="ldm vs Classical", xlab_name="scale parameter for normal")
p
if(save == TRUE) {
  #save_results(m_c_y, paste0(save_path,"/test_iid_serial_independence.xlsx"), "ynm_gumbel_0_1")
  save_results_with_plot( m_L_c, file = paste0(save_path,"/test_ldm_sequential_plot.xlsx"),
                          sheet = "iid_norm_0",
                          p = p,
                          figure_width = 7,
                          figure_height = 5,
                          dpi = 600
  )
}
# ______________________________________
# H0: ldm vs H1: dtrw
# ______________________________________
scale_vals <- seq(1, 2, 1)
m_L_R <- simulate_model(
  param_values = scale_vals,
  T = T,
  n_sim = n_sim,
  generator = dtrw_series,
  param_name =  "scale",      # sd for dtrw and scale for Cauchy
  n_arg = "T",
  test_fun =  test_ldm_sequential,
  series_args = list(dist="cauchy",location=0),
  test_args = list(time = NA)
              )
p = plot_results(m_L_R, "scale", title= "ldm vs dtrw - Cauchy", xlab_name="Scale (σ²)")
p
if(save == TRUE) {
  #save_results(m_c_y, paste0(save_path,"/test_iid_serial_independence.xlsx"), "ynm_gumbel_0_1")
  save_results_with_plot( m_L_R, file = paste0(save_path,"/test_ldm_sequential_plot.xlsx"),
                          sheet = "dtrw_cauchy_0",
                          p = p,
                          figure_width = 7,
                          figure_height = 5,
                          dpi = 600
  )
}
# ______________________________________
# H0: ldm vs H1: ldm (should be Low)
# ______________________________________
theta_vals <- seq(0.2,0.5,0.1)
m_L_L <- simulate_model(
  param_values = theta_vals,
  T = T,
  n_sim = n_sim,
  generator =ldm_series,
  param_name = "theta",      # sd for dtrw and scale for Cauchy
  n_arg = "T",
  test_fun = test_ldm_sequential,
  series_args = list(dist= "norm", mean = 0, sd = 1),
  test_args = list(time = NA)
)

p = plot_results(m_L_L, "theta", "ldm vs ldm", xlab_name=" Theta (Θ) ", ymax= 100)
p
if(save == TRUE) {
  #save_results(m_c_y, paste0(save_path,"/test_iid_serial_independence.xlsx"), "ynm_gumbel_0_1")
  save_results_with_plot( m_L_L, file = paste0(save_path,"/test_ldm_sequential_plot.xlsx"),
                          sheet = "detection",
                          p = p,
                          figure_width = 7,
                          figure_height = 5,
                          dpi = 600
  )
}
