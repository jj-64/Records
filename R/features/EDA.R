## ============================================================
## 0. Setup
## ============================================================

## Core tidy tools
library(tidyverse)

## Statistics / diagnostics
library(moments)      # skewness, kurtosis
library(psych)        # describe()
library(corrplot)     # correlation visualization
library(GGally)       # ggpairs
library(matrixStats)  # fast row/column stats

## Dimensionality / structure
library(FactoMineR)
library(factoextra)

## Optional (comment out if not installed)
# library(rstatix)    # robust stats
# library(reshape2)

## ============================================================
## 1. Input validation
## ============================================================

load("data/feature_matrix.rda")

    ## get feature matrix (numerics only)
feature_names <- setdiff(names(feature_matrix), c("Max_logLik","series_id", "label","label_m", "series"))
feature_only <- feature_matrix[, colnames(feature_matrix) %in% (feature_names)]
# feature_only <- feature_matrix[, sapply(feature_matrix, is.numeric)]

    ##Convert to numeric if needed
feature_only <- feature_only %>% mutate_if(is.factor, as.character) %>% mutate_all(as.numeric)

    ## Assume feature_only is your object
stopifnot(is.matrix(feature_only) || is.data.frame(feature_only))
feature_only <- as.data.frame(feature_only)

    ## Enforce numeric-only
non_numeric <- names(feature_only)[!sapply(feature_only, is.numeric)]
if (length(non_numeric) > 0) {
  stop("Non-numeric features detected: ", paste(non_numeric, collapse = ", "))
}

## Basic dimensions
n_obs  <- nrow(feature_only)
n_feat <- ncol(feature_only)

cat("Observations:", n_obs, "\n")
cat("Features:", n_feat, "\n")

## ============================================================
## 2. Missingness and degeneracy
## ============================================================

## Missing values per feature
na_summary <- tibble(
  feature = names(feature_only),
  n_na    = colSums(is.na(feature_only)),
  pct_na  = colMeans(is.na(feature_only))
) %>% arrange(desc(pct_na))

print(na_summary %>% filter(n_na > 0))

## Constant or near-constant features
variance_summary <- tibble(
  feature = names(feature_only),
  variance = apply(feature_only, 2, var, na.rm = TRUE),
  unique_vals = apply(feature_only, 2, function(z) length(unique(z)))
) %>% arrange(variance)

print(variance_summary %>% filter(variance == 0))

## Flag unusable features
degenerate_features <- variance_summary %>%
  filter(variance == 0 | unique_vals <= 2)

if (nrow(degenerate_features) > 0) {
  warning("Degenerate or quasi-degenerate features detected")
  print(degenerate_features)
}
      # Zero variance → remove
      # Two unique values → possible indicator, verify meaning
      # High NA → consider imputation or exclusion

## Remove degenerate features
feature_only = feature_only[, !(colnames(feature_only) %in%  degenerate_features$feature)]

## ============================================================
## 3. Univariate descriptive statistics
## ============================================================

desc_stats <- psych::describe(feature_only)

View(desc_stats)
# What to look for
# Heavy skewness → log / rank transforms
# Extreme kurtosis → tail-driven statistics (records, extremes)
# Features with very different scales → standardization mandatory

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
## 5. Correlation analysis
## ============================================================

## Pairwise correlation (Spearman recommended)
cor_mat <- cor(feature_only, use = "pairwise.complete.obs", method = "spearman")

## Visual inspection
corrplot(cor_mat,
         method = "color",
         type = "upper",
         tl.cex = 0.6,
         order = "hclust")

## Identify highly correlated pairs
high_corr <- which(abs(cor_mat) > 0.95 & abs(cor_mat) < 1, arr.ind = TRUE)

high_corr_pairs <- tibble(
  feature_1 = rownames(cor_mat)[high_corr[,1]],
  feature_2 = colnames(cor_mat)[high_corr[,2]],
  rho = cor_mat[high_corr]
) %>% filter(feature_1 != feature_2)

View(high_corr_pairs)
# Decision rule
#
# |ρ| > 0.95 → keep one representative
# Prefer:
#   More interpretable
# Lower variance inflation
# Better downstream SHAP signal

## ============================================================
## 6. Scale diagnostics
## ============================================================

scale_summary <- tibble(
  feature = names(feature_only),
  sd = apply(feature_only, 2, sd, na.rm = TRUE),
  IQR = apply(feature_only, 2, IQR, na.rm = TRUE)
) %>% arrange(desc(sd))

View(scale_summary)

## Z-score normalization preview
X_scaled <- scale(feature_only)

summary(X_scaled)
# Rule
#
# Always classify on scaled features
# Keep raw scale only for interpretability and plots

## ============================================================
## 7. PCA
## ============================================================

pca <- prcomp(X_scaled, center = TRUE, scale. = TRUE)

## Variance explained Scree plot
fviz_eig(pca, addlabels = TRUE)

## Feature loadings
loadings <- pca$rotation %>%
  as.data.frame() %>%
  rownames_to_column("feature")

print(head(loadings))

## Contribution to first PCs
fviz_pca_var(pca, col.var = "contrib")

# Interpretation
# Few PCs explaining most variance → redundancy
# Record / extreme features often dominate first PCs
# Later PCs may encode regime-specific structure

## ============================================================
## 8. Unsupervised structure
## ============================================================

## Distance matrix
d <- dist(X_scaled)

## Hierarchical clustering
hc <- hclust(d, method = "ward.D2")

plot(hc, labels = FALSE, main = "Feature-only clustering")

## Optional k-means preview
set.seed(123)
fviz_nbclust(t(X_scaled), kmeans, method = "wss")
km <- kmeans(t(X_scaled), centers = 5, nstart = 20)
km$tot.withinss

fviz_cluster(list(data = t(X_scaled), cluster = km$cluster))
clusters = data.frame( features = names(feature_only), cluster = km$cluster)

# Purpose
# Check whether features contain structure before labels
# If clusters align with known regimes → strong signal
# If not → classification may rely on subtle margins

## ============================================================
## 9. Outliers
## ============================================================

## Row-wise Mahalanobis distance
md <- mahalanobis(t(X_scaled),
                  center = colMeans(t(X_scaled)),
                  cov = cov(t(X_scaled)) )

## Flag extreme observations
cutoff <- qchisq(0.99, df = n_feat)
outliers <- which(md > cutoff)

cat("Outlier rows:", length(outliers), "\n")

## Visual
plot(md, type = "h")
abline(h = cutoff, col = "red")

## ============================================================
## 10. Outputs
## ============================================================

eda_output <- list(
  desc_stats = desc_stats,
  na_summary = na_summary,
  variance_summary = variance_summary,
  correlation_matrix = cor_mat,
  pca = pca,
  outliers = outliers
)

saveRDS(eda_output, "feature_EDA_summary.rds")
