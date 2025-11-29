record_registry <- new.env()

register_record_model <- function(model, mean_fun, var_fun, required_args = NULL) {
  record_registry[[model]] <- list(
    mean_fun = mean_fun,
    var_fun = var_fun,
    required_args = required_args
  )
}

register_record_model(
  "iid",
  mean_fun = rec_count_mean_iid,
  var_fun  = rec_count_var_iid
)

register_record_model(
  "DTRW",
  mean_fun = rec_count_mean_DTRW,
  var_fun  = rec_count_var_DTRW
)

register_record_model(
  "LDM",
  mean_fun = rec_count_mean_LDM,
  var_fun  = rec_count_var_LDM,
  required_args = c("theta")
)

register_record_model(
  "YNM",
  mean_fun = rec_count_mean_YNM,
  var_fun  = rec_count_var_YNM,
  required_args = c("gamma")
)
