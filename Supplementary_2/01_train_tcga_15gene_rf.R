# TCGA preprocessing, training, and cross-validation of 15-gene random forest classifier
#
# Input:
# - dds_*_matched_group.rds files
#
# Output:
# - tcga_15gene_expression.rds
# - rf_15gene_tcga_model.rds
# - rf_15gene_tcga_cv.rds
# - tcga_sample_summary.csv
# - rf_15gene_variable_importance.csv
# - rf_15gene_model_summary.csv
#
# Notes:
# - Edit 'tcga_dir' and 'output_dir' below before running
# - Ensembl IDs are mapped to gene symbols using org.Hs.eg.db
# - Duplicate symbols are collapsed by mean expression

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(randomForest)
  library(caret)
  library(pROC)
  library(ggplot2)
})

# ------------------------------------------------------------------
# Input/output paths
# ------------------------------------------------------------------

tcga_dir   <- "path-to-tcga-dds-files"
output_dir <- "path-to-output-folder"

if (!dir.exists(tcga_dir)) {
  stop("Please set 'tcga_dir' to the folder containing dds_*_matched_group.rds files.")
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ------------------------------------------------------------------
# Define 15-gene panel
# ------------------------------------------------------------------

ttll_genes <- c(
  "TTLL1", "TTLL2", "TTLL4", "TTLL5", "TTLL6", "TTLL7", "TTLL9", "TTLL11", "TTLL13"
)

ccp_genes <- c(
  "AGBL1", "AGBL2", "AGBL3", "AGBL4", "AGBL5", "AGTPBP1"
)

feature_genes <- unique(c(ttll_genes, ccp_genes))

# ------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------

aggregate_dup_symbols <- function(mat, symbols, mode = "mean") {
  ok <- !is.na(symbols) & symbols != ""
  mat <- mat[ok, , drop = FALSE]
  symbols <- symbols[ok]

  if (mode == "sum") {
    out <- rowsum(mat, group = symbols)
  } else {
    sp <- split(seq_len(nrow(mat)), symbols)
    out <- vapply(
      sp,
      function(ix) colMeans(mat[ix, , drop = FALSE]),
      numeric(ncol(mat))
    )
    out <- t(out)
    rownames(out) <- names(sp)
  }

  out
}

plot_cm <- function(pred, truth, title, high_col = "#00B2EE") {
  pred <- factor(pred, levels = levels(truth))
  cm <- confusionMatrix(pred, truth)
  cm_df <- as.data.frame(cm$table)
  colnames(cm_df) <- c("Prediction", "Reference", "Freq")

  ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Freq)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Freq), size = 6, fontface = "bold") +
    scale_fill_gradient(low = "white", high = high_col) +
    labs(title = title, x = "Actual", y = "Predicted") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

# ==================================================================
# 1. PREPROCESS TCGA DATA
# ==================================================================

dds_files <- list.files(
  tcga_dir,
  pattern = "dds_.*_matched_group\\.rds$",
  full.names = TRUE
)

if (length(dds_files) == 0) {
  stop("No TCGA DESeq2 files found in 'tcga_dir'.")
}

tcga_data <- list()

for (file in dds_files) {
  project <- sub("dds_(.*)_matched_group\\.rds", "\\1", basename(file))
  message("Processing ", project)

  dds <- readRDS(file)
  vsd <- vst(dds, blind = TRUE)
  expr_mat <- assay(vsd)

  ens_ids <- sub("\\..*$", "", rownames(expr_mat))
  symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = ens_ids,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )

  expr_sym <- aggregate_dup_symbols(expr_mat, symbols, mode = "mean")

  matched_genes <- intersect(feature_genes, rownames(expr_sym))
  if (length(matched_genes) == 0) next

  expr_df <- t(expr_sym[matched_genes, , drop = FALSE]) %>% as.data.frame()
  expr_df$sample <- rownames(expr_df)

  cond <- tryCatch({
    as.character(colData(dds)[expr_df$sample, "condition"])
  }, error = function(e) NULL)

  if (is.null(cond)) next

  expr_df$condition <- cond
  expr_df$project <- project
  tcga_data[[length(tcga_data) + 1]] <- expr_df
}

tcga_df <- bind_rows(tcga_data) %>%
  dplyr::select(any_of(feature_genes), condition, sample, project) %>%
  filter(if_all(all_of(feature_genes), ~ !is.na(.))) %>%
  mutate(
    condition_bin = factor(
      ifelse(condition == "Primary Tumor", "Tumor", "Normal"),
      levels = c("Normal", "Tumor")
    )
  )

if (nrow(tcga_df) == 0) {
  stop("No TCGA samples available after preprocessing.")
}

saveRDS(tcga_df, file.path(output_dir, "tcga_15gene_expression.rds"))

# ==================================================================
# 2. TRAIN THE 15-GENE MODEL
# ==================================================================

x_tcga <- tcga_df[, feature_genes, drop = FALSE]
y_tcga <- tcga_df$condition_bin

set.seed(123)
rf_15 <- randomForest(
  x = x_tcga,
  y = y_tcga,
  ntree = 500,
  importance = TRUE
)

importance_15 <- importance(rf_15, type = 1) %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  arrange(desc(MeanDecreaseAccuracy))

print(importance_15)

saveRDS(rf_15, file.path(output_dir, "rf_15gene_tcga_model.rds"))

# ==================================================================
# 3. CROSS-VALIDATION OF THE 15-GENE MODEL
# ==================================================================

ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

set.seed(123)
rf_15_cv <- train(
  x = x_tcga,
  y = y_tcga,
  method = "rf",
  metric = "ROC",
  trControl = ctrl
)

print(rf_15_cv)
print(rf_15_cv$results)

saveRDS(rf_15_cv, file.path(output_dir, "rf_15gene_tcga_cv.rds"))

cv_preds <- rf_15_cv$pred %>%
  filter(mtry == rf_15_cv$bestTune$mtry)

roc_tcga <- roc(cv_preds$obs, cv_preds$Tumor, quiet = TRUE)
auc_tcga <- as.numeric(auc(roc_tcga))

roc_tcga_df <- data.frame(
  fpr = 1 - roc_tcga$specificities,
  tpr = roc_tcga$sensitivities
)

ggplot(roc_tcga_df, aes(x = fpr, y = tpr)) +
  geom_line(linewidth = 1.4, color = "#2C7FB8") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
  labs(
    title = "ROC Curve: 15-Gene Model Cross-Validation in TCGA",
    subtitle = sprintf("AUC = %.3f", auc_tcga),
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

p_cm_cv <- plot_cm(
  pred = cv_preds$pred,
  truth = cv_preds$obs,
  title = "Confusion Matrix: 15-Gene Model Cross-Validation in TCGA"
)

print(p_cm_cv)

# ==================================================================
# 4. SAVE SUMMARY TABLES
# ==================================================================

tcga_sample_summary <- tcga_df %>%
  count(project, condition) %>%
  tidyr::pivot_wider(
    names_from = condition,
    values_from = n,
    values_fill = 0
  ) %>%
  arrange(project)

model_summary <- tibble(
  n_samples = nrow(tcga_df),
  n_tumor = sum(tcga_df$condition_bin == "Tumor"),
  n_normal = sum(tcga_df$condition_bin == "Normal"),
  n_genes = length(feature_genes),
  cv_auc = auc_tcga,
  best_mtry = rf_15_cv$bestTune$mtry
)

write_csv(
  tcga_sample_summary,
  file.path(output_dir, "tcga_sample_summary.csv")
)

write_csv(
  importance_15,
  file.path(output_dir, "rf_15gene_variable_importance.csv")
)

write_csv(
  model_summary,
  file.path(output_dir, "rf_15gene_model_summary.csv")
)

message("Saved:")
message(" - ", file.path(output_dir, "tcga_15gene_expression.rds"))
message(" - ", file.path(output_dir, "rf_15gene_tcga_model.rds"))
message(" - ", file.path(output_dir, "rf_15gene_tcga_cv.rds"))
message(" - ", file.path(output_dir, "tcga_sample_summary.csv"))
message(" - ", file.path(output_dir, "rf_15gene_variable_importance.csv"))
message(" - ", file.path(output_dir, "rf_15gene_model_summary.csv"))
