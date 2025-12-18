
## ============================================================
## 4. Distribution visualization
## ============================================================

X_long <- feature_only %>%
  mutate(id = row_number()) %>%
  pivot_longer(-id, names_to = "feature", values_to = "value")

##Get unique features and split into groups of, say, 16
feature_list <- unique(X_long$feature)
feature_chunks <- split(feature_list, ceiling(seq_along(feature_list) / 8))

## Histograms
plots <- map(feature_chunks, ~ {
  ggplot(filter(X_long, feature %in% .x), aes(x = value, fill = regime)) +
    geom_histogram(bins = 50, fill = "grey70", color = "black") +
    facet_wrap(~ feature, scales = "free", ncol = 2) +
    theme_minimal()
})

# Interpretation
# Outliers dominating boxplots → consider robust scaling
# Long right tails → logs or ranks
# Multimodality → regime sensitivity (good candidates)

## ============================================================
## 6. Scale diagnostics
## ============================================================

scale_summary <- tibble(
  feature = names(feature_only),
  sd = apply(feature_only, 2, sd, na.rm = TRUE),
  IQR = apply(feature_only, 2, IQR, na.rm = TRUE)
) %>% arrange(desc(sd))

View(scale_summary)


