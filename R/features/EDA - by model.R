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

## ============================================================
## 1. Input validation
## ============================================================

load("data/feature_matrix.rda")

## get feature matrix (numerics only)
feature_names <- setdiff(names(feature_matrix), c("Max_logLik","series_id", "label", "series"))
feature_only <- feature_matrix[, colnames(feature_matrix) %in% (feature_names)]
# feature_only <- feature_matrix[, sapply(feature_matrix, is.numeric)]

##Convert to numeric if needed
feature_only <- feature_only %>% mutate_if(is.factor, as.character) %>% mutate_all(as.numeric)

regime = feature_matrix$label

#### 4.2  Replace -Inf with -1e5 ( only fix LogLik values)
feature_only[grep("^logLik_", names(feature_only))] <-
  lapply(feature_only[grep("^logLik_", names(feature_only))], function(x) {
    x[x < -1e5] <- -1e5#Inf
    x[!is.finite(x)] <- -1e5
    x
  })

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

print(group_desc)
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
      eta2 = round(eta_squared(fit, partial = FALSE)$Eta2,3)
    )
  }
) %>% bind_rows()

eta2_summary <- eta2_summary %>%
  arrange(desc(eta2))

print(eta2_summary[eta2_summary$eta2 <0.06,])

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

print(kw_summary[kw_summary$p_value >0.05, ])

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

print(quantile_range)

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
feature_chunks <- split(feature_list, ceiling(seq_along(feature_list) / 16))

## Boxplots by regime
ggplot(X_long[X_long$feature == "cv",], aes(x = regime, y = value)) +
  geom_boxplot(outlier.alpha = 0.2) +
  facet_wrap(~ feature, scales = "free", ncol = 4) +
  theme_minimal()


## Density overlays
plots <- map(feature_chunks, ~ {
  ggplot(filter(X_long, feature %in% .x), aes(x = value, fill = regime)) +
    geom_density(alpha = 0.3) +
    facet_wrap(~ feature, scales = "free", ncol = 4) +
    theme_minimal()
})

## display the first plot
plots[[1]]


##ptional: Save each plot
#walk2(plots, seq_along(plots), ~ ggsave(paste0("density_plot_", .y, ".png"), .x, width = 12, height = 8))


## ============================================================
## 9. Supervised PCA (sanity check)
## ============================================================

pca <- PCA(
  feature_only,
  scale.unit = TRUE,
  graph = FALSE
)

fviz_pca_ind(pca,
             habillage = regime,
             addEllipses = TRUE,
             ellipse.level = 0.95)
## ============================================================
## 10. Decision rule for feature removal (explicit)
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
selected_features <- eta2_summary %>%
  filter(eta2 > 0.02) %>%
  pull(feature)

X_reduced <- X %>%
  select(all_of(selected_features), regime)

saveRDS(X_reduced, "features_regime_discriminative.rds")
