# TCGA training and cross-validation of 9-gene random forest classifier
#
# This script:
# 1. reads matched TCGA tumour-normal DESeq2 objects
# 2. extracts VST-normalized expression for a fixed 9-gene panel
# 3. trains a random forest classifier on TCGA
# 4. evaluates the model by 5-fold cross-validation
# 5. saves the trained model and cross-validation object
#
# Input:
# - dds_*_matched_group.rds files
#
# Output:
# - rf_9gene_tcga_model.rds
# - rf_9gene_tcga_cv.rds
# - tcga_9gene_expression.rds
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
# Fixed 9-gene panel
# ------------------------------------------------------------------

genes_9 <- c(
  "TTLL4", "TTLL5", "TTLL6", "TTLL7", "TTLL11",
  "AGBL1", "AGBL2", "AGBL3", "AGBL5"
)

# ------------------------------------------------------------------
# Helper function
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

# ------------------------------------------------------------------
# Load and process TCGA data
# ------------------------------------------------------------------

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

  matched_genes <- intersect(genes_9, rownames(expr_sym))
  if (length(matched_genes) == 0) next

  expr_df <- t(expr_sym[matched_genes, , drop = FALSE]) %>% as.data.frame()

  for (g in setdiff(genes_9, colnames(expr_df))) {
    expr_df[[g]] <- NA_real_
  }
  expr_df <- expr_df[, genes_9, drop = FALSE]

  for (g in colnames(expr_df)) {
    if (anyNA(expr_df[[g]])) {
      med_g <- median(expr_df[[g]], na.rm = TRUE)
      expr_df[[g]][is.na(expr_df[[g]])] <- if (is.finite(med_g)) med_g else 0
    }
  }

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
  dplyr::select(all_of(genes_9), condition, sample, project) %>%
  filter(if_all(all_of(genes_9), ~ !is.na(.))) %>%
  mutate(
    condition_bin = factor(
      ifelse(condition == "Primary Tumor", "Tumor", "Normal"),
      levels = c("Normal", "Tumor")
    )
  )

if (nrow(tcga_df) == 0) {
  stop("No TCGA samples available after preprocessing.")
}

saveRDS(tcga_df, file.path(output_dir, "tcga_9gene_expression.rds"))

# ------------------------------------------------------------------
# Train 9-gene random forest
# ------------------------------------------------------------------

x_tcga <- tcga_df[, genes_9, drop = FALSE]
y_tcga <- tcga_df$condition_bin

set.seed(123)
rf_9 <- randomForest(
  x = x_tcga,
  y = y_tcga,
  ntree = 500,
  importance = TRUE
)

saveRDS(rf_9, file.path(output_dir, "rf_9gene_tcga_model.rds"))

# ------------------------------------------------------------------
# 5-fold cross-validation
# ------------------------------------------------------------------

ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

set.seed(123)
rf_9_cv <- train(
  x = x_tcga,
  y = y_tcga,
  method = "rf",
  metric = "ROC",
  trControl = ctrl
)

saveRDS(rf_9_cv, file.path(output_dir, "rf_9gene_tcga_cv.rds"))

# ------------------------------------------------------------------
# Optional console summary
# ------------------------------------------------------------------

cv_preds <- rf_9_cv$pred %>%
  filter(mtry == rf_9_cv$bestTune$mtry)

roc_tcga <- roc(cv_preds$obs, cv_preds$Tumor, quiet = TRUE)
auc_tcga <- as.numeric(auc(roc_tcga))

message("Saved:")
message(" - ", file.path(output_dir, "rf_9gene_tcga_model.rds"))
message(" - ", file.path(output_dir, "rf_9gene_tcga_cv.rds"))
message(" - ", file.path(output_dir, "tcga_9gene_expression.rds"))
message(sprintf("TCGA 5-fold CV AUC: %.3f", auc_tcga))
