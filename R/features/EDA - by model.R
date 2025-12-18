## ============================================================
## 0. Setup
## ============================================================
set.seed(12345)

required_packages <- c(
  "tidyverse",
  "moments",      # skewness, kurtosis
  "psych",        # describe()
  "corrplot",     # correlation visualization
  "GGally",       # ggpairs
  "matrixStats",  # fast row/column stats
  "FactoMineR","factoextra","effectsize",  ## Dimensionality / structure
  "rstatix",    # robust stats
  "reshape2",
  "dplyr",
  "ggplot2",
  "purrr"
)

install_if_missing(required_packages)

## Load libraries
load_package(required_packages = required_packages)

options (scipen= 99999)

## ============================================================
## Group-wise summaries
## ============================================================

group_desc <- feature_only %>%
  group_by(feature_matrix$label) %>%
  summarise(across(
    .cols = where(is.numeric),
    .fns = list(
      mean = ~ mean(.x, na.rm = TRUE),
      sd   = ~ sd(.x, na.rm = TRUE),
      med  = ~ median(.x, na.rm = TRUE),
      IQR  = ~ IQR(.x, na.rm = TRUE)
    ),
    .names = "{.col}_{.fn}"
  ))


group_desc = as.data.frame(t(group_desc))
colnames(group_desc) = c("Classical", "DTRW", "LDM", "YNM")
#group_desc$feature = rownames(group_desc)
group_desc = group_desc[-1,]
group_desc[,1:4] = apply(group_desc[,1:4], 2, as.numeric)
group_desc$variability = apply(group_desc[,1:4], 1, FUN = function(x) abs(sd(x)/mean(x)))
View(group_desc)

## recommended to remove
group_desc[which(group_desc$variability ==0 ), ]

group_desc[which(group_desc$variability <0.05 ), ]
# Interpretation
# Identical summaries across regimes → drop candidate
# Large IQR shifts → tail-sensitive feature

## ============================================================
## Effect size screening (preferred over p-values)
## ============================================================

eta2_summary <- lapply(
  names(feature_only)[sapply(feature_only, is.numeric)],
  function(f) {
    fit <- aov(feature_only[[f]] ~ regime)
    data.frame(
      feature = f,
      eta2 = round(effectsize::eta_squared(fit, partial = FALSE)$Eta2,3)
    )
  }
) %>% bind_rows()

eta2_summary <- eta2_summary %>%
  arrange(desc(eta2))

## recommended to remove
print(eta2_summary[which(eta2_summary$eta2 <0.05 | is.na(eta2_summary$eta2)),])

# η² < 0.01 → negligible
# 0.01–0.06 → small
# 0.06 → meaningful

## ============================================================
## Non-Parameteric Seperation: Kruskal
## ============================================================
kw_summary <- lapply(
  names(feature_only)[sapply(feature_only, is.numeric)],
  function(f) {
    kt <- kruskal.test(feature_only[[f]] ~ regime)
    data.frame(
      feature = f,
      p_value = round(kt$p.value,3)
    )
  }
) %>% bind_rows() %>%
  arrange(p_value)

## recommended to remove
print(kw_summary[which(kw_summary$p_value >0.05 | is.na(kw_summary$p_value)) , ])

## ============================================================
## Quantile-level regime differentiation (very important for records)
## ============================================================
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

## ============================================================
## Regime-conditioned plots
## ============================================================

X_long <- cbind(feature_only,regime) %>%
  pivot_longer(
    cols = where(is.numeric),
    names_to = "feature",
    values_to = "value"
  )

##Get unique features and split into groups of, say, 16
feature_list <- unique(X_long$feature)
feature_chunks <- split(feature_list, ceiling(seq_along(feature_list) / 8))

## Density overlays
plots <- map(feature_chunks, ~ {
  ggplot(filter(X_long, feature %in% .x), aes(x = value, fill = regime)) +
    geom_density(alpha = 0.3) +
    facet_wrap(~ feature, scales = "free", ncol = 2) +
    theme_minimal()
})

## Boxplots by regime
plotsb <- map(feature_chunks, ~ {
  ggplot(filter(X_long, feature %in% .x), aes(x = regime, y = value)) +
    geom_boxplot(outlier.alpha = 0.2) +
    facet_wrap(~ feature, scales = "free", ncol = 2) +
    theme_minimal()
})

## display the first plot
p=3
plots[[p]]
plotsb[[p]]

##ptional: Save each plot
#walk2(plots, seq_along(plots), ~ ggsave(paste0("density_plot_", .y, ".png"), .x, width = 12, height = 8))




## ============================================================
## 9. Decision rule for feature removal (explicit)
## ============================================================

# I recommend retaining features that satisfy at least one:
#   η² > 0.02
# Kruskal–Wallis p < 0.01
# Quantile-range (q95 or q05) above 75th percentile
# Theoretical justification (record asymmetry, renewal structure)
# Everything else is defensibly removable.

## ============================================================
## 11. Output: filtered feature set
## ============================================================
    ##eta2 removal
remove_features = eta2_summary %>%
  filter(eta2 < 0.02)%>%
  pull(feature)

remove_features = c(remove_features, kw_summary %>%
  filter(p_value > 0.05) %>%
  pull(feature)
)

remove_features = c(remove_features, quantile_range %>%
                      filter(range_across_regimes == 0) %>%
                      pull(feature)
)


# selected_features <- eta2_summary %>%
#   filter(eta2 > 0.06) %>%
#   pull(feature)

selected_features <- setdiff(names(feature_matrix), unique(remove_features))

feature_matrix_reduced <- feature_matrix %>%
  select(all_of(selected_features))

saveRDS(feature_matrix_reduced, "data/feature_matrix_reduced.rds")
