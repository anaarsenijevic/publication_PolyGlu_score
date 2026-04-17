# Permutation-label sanity check for the TCGA 9-gene random forest classifier
#
# This script:
# 1. reads matched TCGA tumour-normal DESeq2 objects
# 2. extracts VST-transformed expression for the fixed 9-gene panel
# 3. maps Ensembl IDs to HGNC symbols using org.Hs.eg.db
# 4. trains a 9-gene random forest classifier using 5-fold cross-validation
# 5. computes the cross-validated AUC for the real labels
# 6. permutes the labels repeatedly and recomputes cross-validated AUC
# 7. builds a null distribution of AUC values
# 8. reports an empirical p-value and saves the null results
#
# Input:
# - dds_*_matched_group.rds files
#
# Output:
# - perm_null_auc.csv
# - perm_null_auc_hist.png
#
# Notes:
# - Edit 'input_dir' and 'output_dir' below before running
# - Ensembl IDs are mapped to HGNC symbols using org.Hs.eg.db
# - Duplicate gene symbols are collapsed by mean expression
# - The classifier uses the fixed 9-gene panel:
#   TTLL4, TTLL5, TTLL6, TTLL7, TTLL11, AGBL1, AGBL3, AGBL4, AGBL5
# - The same cross-validation folds are reused for the real and permuted runs

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(caret)
  library(randomForest)
  library(pROC)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# ------------------------------------------------------------------
# Input/output paths
# ------------------------------------------------------------------

input_dir <- "path-to-your-dds-files"
output_dir <- "path-to-output-folder"

if (!dir.exists(input_dir)) {
  stop("Please set 'input_dir' to the folder containing dds_*_matched_group.rds files.")
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------

set.seed(123)

feature_genes <- c(
  "TTLL4", "TTLL5", "TTLL6", "TTLL7", "TTLL11",
  "AGBL1", "AGBL3", "AGBL4", "AGBL5"
)

label_levels <- c("Tumor", "Normal")
n_perm <- 100

# ------------------------------------------------------------------
# Helper: collapse duplicate gene symbols
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
# Load matched TCGA DESeq2 objects and extract 9-gene expression
# ------------------------------------------------------------------

message("Loading matched TCGA DESeq2 objects...")

dds_files <- list.files(
  input_dir,
  pattern = "dds_.*_matched_group\\.rds$",
  full.names = TRUE
)

if (length(dds_files) == 0) {
  stop("No TCGA DESeq2 files found in 'input_dir'.")
}

all_data <- list()

for (file in dds_files) {
  project <- sub("dds_(.*)_matched_group\\.rds", "\\1", basename(file))
  dds <- readRDS(file)

  vsd <- DESeq2::vst(dds, blind = TRUE)
  norm_counts <- assay(vsd)

  ensembl_ids <- sub("\\..*$", "", rownames(norm_counts))
  symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = ensembl_ids,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )

  norm_counts_sym <- aggregate_dup_symbols(norm_counts, symbols, mode = "mean")

  matched_genes <- intersect(feature_genes, rownames(norm_counts_sym))
  if (!length(matched_genes)) next

  expr <- t(norm_counts_sym[matched_genes, , drop = FALSE]) %>%
    as.data.frame()

  for (g in setdiff(feature_genes, colnames(expr))) {
    expr[[g]] <- NA_real_
  }

  expr <- expr[, feature_genes, drop = FALSE]
  expr$sample <- rownames(expr)

  cond <- tryCatch(
    SummarizedExperiment::colData(dds)[expr$sample, "condition"],
    error = function(e) NULL
  )

  if (is.null(cond)) next

  expr$condition <- as.character(cond)
  expr$project <- project
  all_data[[length(all_data) + 1]] <- expr
}

tcga_df <- dplyr::bind_rows(all_data)

tcga_clean <- tcga_df %>%
  dplyr::select(dplyr::all_of(feature_genes), condition, sample, project) %>%
  dplyr::filter(dplyr::if_all(dplyr::all_of(feature_genes), ~ !is.na(.))) %>%
  dplyr::mutate(
    condition_bin = factor(
      ifelse(condition == "Primary Tumor", "Tumor", "Normal"),
      levels = label_levels
    )
  )

if (nrow(tcga_clean) == 0) {
  stop("No TCGA samples available after preprocessing.")
}

# ------------------------------------------------------------------
# Training matrix and labels
# ------------------------------------------------------------------

x_train <- tcga_clean %>%
  dplyr::select(dplyr::all_of(feature_genes)) %>%
  as.data.frame()

y_train <- factor(
  droplevels(tcga_clean$condition_bin),
  levels = label_levels
)

stopifnot(all(levels(y_train) == c("Tumor", "Normal")))

# ------------------------------------------------------------------
# Fixed 5-fold cross-validation folds
# ------------------------------------------------------------------

folds <- caret::createFolds(y_train, k = 5, returnTrain = TRUE)

ctrl_cv <- caret::trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = caret::twoClassSummary,
  savePredictions = "final",
  index = folds,
  allowParallel = TRUE
)

# ------------------------------------------------------------------
# Helper to compute CV AUC from caret predictions
# ------------------------------------------------------------------

auc_from_pred <- function(pred_df, mtry_val) {
  pred_best <- pred_df %>%
    dplyr::filter(mtry == mtry_val)

  roc_obj <- pROC::roc(
    response = pred_best$obs,
    predictor = pred_best[["Tumor"]],
    levels = c("Normal", "Tumor"),
    quiet = TRUE
  )

  as.numeric(pROC::auc(roc_obj))
}

# ------------------------------------------------------------------
# Real-label model
# ------------------------------------------------------------------

set.seed(123)

rf_cv_real <- caret::train(
  x = x_train,
  y = y_train,
  method = "rf",
  trControl = ctrl_cv,
  metric = "ROC",
  tuneLength = 5,
  ntree = 1000,
  importance = TRUE,
  nodesize = 5
)

best_mtry <- rf_cv_real$bestTune$mtry
auc_real_cv <- auc_from_pred(rf_cv_real$pred, best_mtry)

cat("CV AUC (real labels):", round(auc_real_cv, 3), "\n")

# ------------------------------------------------------------------
# One shuffled-label run (quick check)
# ------------------------------------------------------------------

y_perm1 <- factor(sample(y_train), levels = levels(y_train))

set.seed(123)

rf_cv_perm1 <- caret::train(
  x = x_train,
  y = y_perm1,
  method = "rf",
  trControl = ctrl_cv,
  metric = "ROC",
  tuneGrid = data.frame(mtry = best_mtry),
  ntree = 1000,
  importance = FALSE,
  nodesize = 5
)

auc_perm1 <- auc_from_pred(rf_cv_perm1$pred, best_mtry)
cat("CV AUC (1 shuffled labels):", round(auc_perm1, 3), "\n")

# ------------------------------------------------------------------
# Null distribution from repeated label permutations
# ------------------------------------------------------------------

perm_aucs <- numeric(n_perm)

pb <- txtProgressBar(min = 0, max = n_perm, style = 3)

for (b in seq_len(n_perm)) {
  y_perm <- factor(sample(y_train), levels = levels(y_train))

  rf_cv_perm <- caret::train(
    x = x_train,
    y = y_perm,
    method = "rf",
    trControl = ctrl_cv,
    metric = "ROC",
    tuneGrid = data.frame(mtry = best_mtry),
    ntree = 1000,
    importance = FALSE,
    nodesize = 5
  )

  perm_aucs[b] <- auc_from_pred(rf_cv_perm$pred, best_mtry)
  setTxtProgressBar(pb, b)
}

close(pb)

# ------------------------------------------------------------------
# Empirical p-value
# ------------------------------------------------------------------

p_emp <- (1 + sum(perm_aucs >= auc_real_cv)) / (n_perm + 1)

cat(
  sprintf(
    "Real AUC = %.3f | Null mean = %.3f (sd = %.3f) | Empirical p ≈ %.4g\n",
    auc_real_cv, mean(perm_aucs), sd(perm_aucs), p_emp
  )
)

# ------------------------------------------------------------------
# Save null results
# ------------------------------------------------------------------

readr::write_csv(
  tibble::tibble(
    perm_id = seq_len(n_perm),
    auc = perm_aucs,
    auc_real = auc_real_cv,
    p_emp = p_emp
  ),
  file.path(output_dir, "perm_null_auc.csv")
)

# ------------------------------------------------------------------
# Plot histogram of null distribution
# ------------------------------------------------------------------

png(
  file.path(output_dir, "perm_null_auc_hist.png"),
  width = 1800,
  height = 1200,
  res = 300
)

hist(
  perm_aucs,
  breaks = seq(0, 1, by = 0.05),
  col = "lightgray",
  border = "white",
  xlim = c(0, 1),
  xaxt = "n",
  main = sprintf("Null distribution of CV AUC from permuted models (n = %d)", n_perm),
  xlab = "Cross-validated AUC (permuted labels)"
)

axis(
  1,
  at = seq(0, 1, by = 0.1),
  labels = sprintf("%.1f", seq(0, 1, by = 0.1))
)

abline(v = auc_real_cv, col = "red", lwd = 2)

mtext(
  sprintf("Observed = %.3f | p = %.4f", auc_real_cv, p_emp),
  side = 3
)

dev.off()
