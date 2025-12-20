
## =====Quantile-level regime differentiation (very important for records) ===========

## This captures features that differ only in tails.
feature_only2 = feature_only
feature_only2$regime = regime
quantile_diff <- feature_only2 %>%
  group_by(regime) %>%
  summarise(across(
    .cols = where(is.numeric),
    .fns = list(
      q05 = ~ quantile(.x, 0.05, na.rm = TRUE),
      q50 = ~ quantile(.x, 0.50, na.rm = TRUE),
      q95 = ~ quantile(.x, 0.95, na.rm = TRUE)
    ),
    .names = "{.col}_{.fn}"
  ))
rm(feature_only2)

## Measure range across regimes
quantile_range <- quantile_diff %>%
  pivot_longer(-regime) %>%
  group_by(name) %>%
  summarise(
    range_across_regimes = max(value) - min(value)
  ) %>%
  separate(name, into = c("feature", "quantile"), sep = "_q")

# quantile_range = quantile_range %>% pivot_wider(id_cols = "feature",
#                                names_from = "quantile",
#                                values_from = "range_across_regimes"
#                                )


View(quantile_range)

## recommended to remove
quantile_range[which(quantile_range$range_across_regimes <0.05), ]

# Interpretation
# Large q95 separation → upper-record sensitivity
# Large q05 separation → lower-record sensitivity



##ptional: Save each plot
#walk2(plots, seq_along(plots), ~ ggsave(paste0("density_plot_", .y, ".png"), .x, width = 12, height = 8))


## =====9. Decision rule for feature removal (explicit)=======

# I recommend retaining features that satisfy at least one:
#   η² > 0.02
# Kruskal–Wallis p < 0.01
# Quantile-range (q95 or q05) above 75th percentile
# Theoretical justification (record asymmetry, renewal structure)
# Everything else is defensibly removable.

saveRDS(feature_matrix_reduced, "data/feature_matrix_reduced.rds")
