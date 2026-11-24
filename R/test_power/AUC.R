library(readxl)
library(pracma)
library(dplyr)
library(purrr)
library(stringr)

# 📂 1. Folder containing all Excel files
folder_path <- "inst/extdata/test_power_two_by_two"

# 📜 2. Get all Excel file paths
excel_files <- list.files(folder_path, pattern = "\\.xlsx$", full.names = TRUE)

# 🧮 3. Function to compute AUCs for one data frame
compute_auc <- function(x,y) {
  ## Apply to all other numeric columns
  auc_raw <- sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
  #auc_raw <- sum(diff(x)*rollmean(y,2))
  auc_norm <- auc_raw / ((max(x) - min(x)) * 100)  # since y max = 100
  return(auc_norm)
}

# 🧩 4. Loop through each file and sheet
all_auc_results <- list()

for (file in excel_files) {

  file_name <- tools::file_path_sans_ext(basename(file))
  sheet_names <- readxl::excel_sheets(file)

  for (sheet in sheet_names) {
    message("Processing: ", file_name, " | Sheet: ", sheet)

    # Read sheet
    df <- as.data.frame(suppressMessages(readxl::read_excel(file, sheet = sheet, skip = 2)))

    # Compute AUCs
    auc_values <- tryCatch({
      apply(df[, -1], 2, compute_auc, x = df[,1])
    }, error = function(e) {
      warning(paste("Error in", file_name, sheet, ":", e$message))
      return(NULL)
    })

    if (!is.null(auc_values)) {
      key <- paste(file_name, sheet, sep = ".")
      all_auc_results[[key]] <- auc_values
    }
  }

}

# 🧾 5. Convert to tidy data frame (optional)
auc_df <- bind_rows(lapply(names(all_auc_results), function(k) {
  data.frame(
    File_Sheet = k,
    Variable = names(all_auc_results[[k]]),
    AUC = as.numeric(all_auc_results[[k]]),
    stringsAsFactors = FALSE
  )
}), .id = NULL)

# ✅ View results
head(auc_df)

# ✅  Extract Average Power and Type I Error
auc_Power_avg = subset(auc_df, Variable == "average")
auc_typeI = auc_df[grepl(pattern = "\\.detection", x=auc_df$File_Sheet),]
auc_typeI_avg = subset(auc_typeI, Variable == "average")


auc_Power_avg
auc_typeI_avg

# ✅  Extract Average Power  for each test
auc_Power_avg[grepl(pattern = "iid_serial_independence", x=auc_Power_avg$File_Sheet) ,]
auc_Power_avg[grepl(pattern = "iid_rec_count", x=auc_Power_avg$File_Sheet) ,]

auc_Power_avg[grepl(pattern = "ynm_rec_count", x=auc_Power_avg$File_Sheet) ,]
auc_Power_avg[grepl(pattern = "ynm_chisq", x=auc_Power_avg$File_Sheet) ,]
auc_Power_avg[grepl(pattern = "ynm_rec_gap", x=auc_Power_avg$File_Sheet) ,]

auc_Power_avg[grepl(pattern = "ldm_rec_count", x=auc_Power_avg$File_Sheet) ,]
auc_Power_avg[grepl(pattern = "ldm_trend", x=auc_Power_avg$File_Sheet) ,]
auc_Power_avg[grepl(pattern = "ldm_sequential", x=auc_Power_avg$File_Sheet) ,]

auc_Power_avg[grepl(pattern = "dtrw_rec_count", x=auc_Power_avg$File_Sheet) ,]
auc_Power_avg[grepl(pattern = "dtrw_increment", x=auc_Power_avg$File_Sheet) ,]

## Save

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

save_results( auc_df, file = paste0("inst/extdata/AUC/AUC.xlsx"),sheet = "AUC")
save_results( auc_Power_avg, file = paste0("inst/extdata/AUC/AUC.xlsx"),sheet = "auc_Power_avg")
save_results( auc_typeI, file = paste0("inst/extdata/AUC/AUC.xlsx"),sheet = "auc_typeI")
save_results( auc_typeI_avg, file = paste0("inst/extdata/AUC/AUC.xlsx"),sheet = "auc_typeI_avg")
