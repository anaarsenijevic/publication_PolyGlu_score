# External validation of 9-gene random forest models in Korean CRC cohort
#
# This script:
# 1. reads matched TCGA tumour-normal DESeq2 objects
# 2. extracts VST-transformed expression for the fixed 9-gene panel
# 3. trains two random forest classifiers:
#    - a pan-cancer TCGA model
#    - a CRC-only TCGA model (TCGA-COAD and TCGA-READ)
# 4. evaluates both models by 5-fold cross-validation
# 5. reads the external Korean CRC count matrix
# 6. builds a VST-transformed external expression matrix for the same 9-gene panel
# 7. applies both models to the external cohort
# 8. plots ROC curves and confusion matrices
# 9. plots PCA of the external CRC cohort using the 9-gene panel
#
# Input:
# - dds_*_matched_group.rds files
# - CMCBSN_expectedcount_342.txt
#
# Notes:
# - Edit the paths below before running
# - Download CMCBSN_expectedcount_342.txt from:
#   https://zenodo.org/records/8333650
# - Place CMCBSN_expectedcount_342.txt in 'crc_input_dir'
# - The external file is expected to have gene symbols in the first column
# - Sample labels are inferred from IDs containing '-T-' or '-N-'

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(biomaRt)
  library(randomForest)
  library(caret)
  library(pROC)
  library(ggplot2)
})

# ------------------------------------------------------------------
# Input paths
# ------------------------------------------------------------------

tcga_dir <- "path-to-your-dds-files"
crc_input_dir <- "path-to-korean-crc-folder"

crc_file <- file.path(crc_input_dir, "CMCBSN_expectedcount_342.txt")

if (!dir.exists(tcga_dir)) {
  stop("Please set 'tcga_dir' to the folder containing dds_*_matched_group.rds files.")
}

if (!dir.exists(crc_input_dir)) {
  stop("Please set 'crc_input_dir' to the folder containing CMCBSN_expectedcount_342.txt.")
}

if (!file.exists(crc_file)) {
  stop("Could not find: ", crc_file)
}

# ------------------------------------------------------------------
# 9-gene panel
# ------------------------------------------------------------------

feature_genes <- c(
  "TTLL4", "TTLL5", "TTLL6", "TTLL7", "TTLL11",
  "AGBL1", "AGBL3", "AGBL4", "AGBL5"
)

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

get_proba <- function(model_or_train, newdata) {
  p <- predict(model_or_train, newdata = newdata, type = "prob")

  if (is.data.frame(p)) {
    if ("Tumor" %in% names(p)) return(p$Tumor)
    return(p[[ncol(p)]])
  }

  if (is.matrix(p)) {
    if (is.null(colnames(p))) {
      colnames(p) <- c("Normal", "Tumor")[seq_len(ncol(p))]
    }
    if ("Tumor" %in% colnames(p)) return(p[, "Tumor"])
    return(p[, ncol(p)])
  }

  as.numeric(p)
}

best_thr_bacc <- function(truth, p_tumor) {
  r <- pROC::roc(truth, p_tumor, quiet = TRUE, levels = c("Normal", "Tumor"))
  cd <- as.data.frame(
    pROC::coords(
      r,
      x = r$thresholds,
      ret = c("threshold", "sensitivity", "specificity"),
      transpose = FALSE
    )
  )
  cd$bacc <- (cd$sensitivity + cd$specificity) / 2
  cd$threshold[which.max(cd$bacc)]
}

roc_df_from <- function(roc_obj, label) {
  data.frame(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities,
    model = label
  )
}

plot_roc_pair <- function(roc_cv, roc_ext, title_txt, sub_fmt) {
  ggplot(
    bind_rows(
      roc_df_from(roc_cv, "TCGA CV"),
      roc_df_from(roc_ext, "External")
    ),
    aes(x = fpr, y = tpr, color = model)
  ) +
    geom_line(linewidth = 1.4) +
    geom_abline(linetype = "dashed", color = "gray70") +
    scale_color_manual(values = c("TCGA CV" = "coral1", "External" = "gray70")) +
    labs(
      title = title_txt,
      subtitle = sprintf(sub_fmt, pROC::auc(roc_cv), pROC::auc(roc_ext)),
      x = "False Positive Rate",
      y = "True Positive Rate",
      color = "Dataset"
    ) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom"
    )
}

cm_plot <- function(pred, truth, title_txt) {
  tb <- table(Prediction = pred, Reference = truth)
  df <- as.data.frame(tb)
  colnames(df) <- c("Prediction", "Reference", "Freq")

  ggplot(df, aes(x = Reference, y = Prediction, fill = Freq)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Freq), size = 5, fontface = "bold") +
    scale_fill_gradient(low = "white", high = "slategray4") +
    labs(title = title_txt, x = "True Label", y = "Predicted Label") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

# ------------------------------------------------------------------
# Gene symbol mapping for TCGA rows
# ------------------------------------------------------------------

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_map <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  mart = ensembl
) %>%
  filter(external_gene_name != "") %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

# ------------------------------------------------------------------
# Load and process TCGA matched cohorts
# ------------------------------------------------------------------

dds_files <- list.files(
  tcga_dir,
  pattern = "dds_.*_matched_group\\.rds$",
  full.names = TRUE
)

all_data <- list()

for (file in dds_files) {
  project <- sub("dds_(.*)_matched_group\\.rds", "\\1", basename(file))
  dds <- readRDS(file)

  ensembl_ids <- sub("\\..*$", "", rownames(dds))
  symbols <- gene_map$external_gene_name[match(ensembl_ids, gene_map$ensembl_gene_id)]

  vsd <- vst(dds, blind = TRUE)
  norm_counts <- assay(vsd)
  rownames(norm_counts) <- symbols

  matched_genes <- intersect(feature_genes, rownames(norm_counts))
  if (!length(matched_genes)) next

  expr <- t(norm_counts[matched_genes, , drop = FALSE]) %>% as.data.frame()

  for (g in setdiff(feature_genes, colnames(expr))) {
    expr[[g]] <- NA_real_
  }

  expr <- expr[, feature_genes, drop = FALSE]
  expr$sample <- rownames(expr)

  cond <- tryCatch(
    colData(dds)[expr$sample, "condition"],
    error = function(e) NULL
  )

  if (is.null(cond)) next

  expr$condition <- as.character(cond)
  expr$project <- project
  all_data[[length(all_data) + 1]] <- expr
}

tcga_df <- bind_rows(all_data)

tcga_clean <- tcga_df %>%
  dplyr::select(all_of(feature_genes), condition, sample, project) %>%
  filter(if_all(all_of(feature_genes), ~ !is.na(.))) %>%
  mutate(
    condition_bin = factor(
      ifelse(condition == "Primary Tumor", "Tumor", "Normal"),
      levels = c("Normal", "Tumor")
    )
  )

# ------------------------------------------------------------------
# Train models: pan-cancer and CRC-only
# ------------------------------------------------------------------

ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = TRUE
)

# Pan-cancer model
x_pan <- tcga_clean[, feature_genes, drop = FALSE]
y_pan <- tcga_clean$condition_bin

set.seed(123)
rf_cv_pan <- train(
  x = x_pan,
  y = y_pan,
  method = "rf",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5
)

# CRC-only model
tcga_crc <- tcga_clean %>%
  filter(project %in% c("TCGA-COAD", "TCGA-READ"))

x_crc <- tcga_crc[, feature_genes, drop = FALSE]
y_crc <- tcga_crc$condition_bin

set.seed(123)
rf_cv_crc <- train(
  x = x_crc,
  y = y_crc,
  method = "rf",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5
)

roc_pan_cv <- with(
  rf_cv_pan$pred %>% filter(mtry == rf_cv_pan$bestTune$mtry),
  pROC::roc(obs, Tumor, quiet = TRUE)
)

roc_crc_cv <- with(
  rf_cv_crc$pred %>% filter(mtry == rf_cv_crc$bestTune$mtry),
  pROC::roc(obs, Tumor, quiet = TRUE)
)

# ------------------------------------------------------------------
# External Korean CRC cohort
# ------------------------------------------------------------------

raw_tab <- read.delim(crc_file, check.names = FALSE)

gene_symbols_infile <- raw_tab[[1]]
counts_mat <- as.data.frame(raw_tab[, -1, drop = FALSE])
rownames(counts_mat) <- gene_symbols_infile
counts_mat[] <- lapply(counts_mat, function(x) as.numeric(as.character(x)))

# collapse duplicate symbols by sum of counts
counts_sym <- as.data.frame(
  rowsum(
    as.matrix(counts_mat),
    group = rownames(counts_mat),
    reorder = FALSE,
    na.rm = TRUE
  )
)

sids <- colnames(counts_sym)
tissue <- ifelse(grepl("-T-", sids), "T", ifelse(grepl("-N-", sids), "N", NA))

sample_info <- data.frame(
  Sample = sids,
  Project = "CRC",
  Tissue = tissue,
  row.names = sids,
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(Tissue))

counts_sym <- counts_sym[, rownames(sample_info), drop = FALSE]

dds_ext <- DESeqDataSetFromMatrix(
  countData = round(counts_sym),
  colData = sample_info,
  design = ~ Tissue
)

vsd_ext <- vst(dds_ext, blind = TRUE)
vst_expr <- assay(vsd_ext)

available_genes <- intersect(feature_genes, rownames(vst_expr))
missing_genes <- setdiff(feature_genes, available_genes)

if (length(missing_genes) > 0) {
  message("Missing genes in external set: ", paste(missing_genes, collapse = ", "))
}

X_ext <- t(vst_expr[available_genes, , drop = FALSE]) %>% as.data.frame()

for (g in missing_genes) {
  X_ext[[g]] <- NA_real_
}

X_ext <- X_ext[, feature_genes, drop = FALSE]

for (g in colnames(X_ext)) {
  if (anyNA(X_ext[[g]])) {
    med_g <- median(X_ext[[g]], na.rm = TRUE)
    X_ext[[g]][is.na(X_ext[[g]])] <- if (is.finite(med_g)) med_g else 0
  }
}

truth_ext <- factor(
  ifelse(sample_info$Tissue == "T", "Tumor", "Normal"),
  levels = c("Normal", "Tumor")
)

# ------------------------------------------------------------------
# Predictions on external CRC cohort
# ------------------------------------------------------------------

p_pan_ext <- get_proba(rf_cv_pan, X_ext)
p_crc_ext <- get_proba(rf_cv_crc, X_ext)

roc_pan_ext <- pROC::roc(truth_ext, p_pan_ext, quiet = TRUE)
roc_crc_ext <- pROC::roc(truth_ext, p_crc_ext, quiet = TRUE)

# ------------------------------------------------------------------
# Thresholds: 0.5 and optimal balanced accuracy
# ------------------------------------------------------------------

thr_pan_opt <- best_thr_bacc(truth_ext, p_pan_ext)
thr_crc_opt <- best_thr_bacc(truth_ext, p_crc_ext)

pred_pan_05 <- factor(
  ifelse(p_pan_ext >= 0.5, "Tumor", "Normal"),
  levels = c("Normal", "Tumor")
)

pred_pan_opt <- factor(
  ifelse(p_pan_ext >= thr_pan_opt, "Tumor", "Normal"),
  levels = c("Normal", "Tumor")
)

pred_crc_05 <- factor(
  ifelse(p_crc_ext >= 0.5, "Tumor", "Normal"),
  levels = c("Normal", "Tumor")
)

pred_crc_opt <- factor(
  ifelse(p_crc_ext >= thr_crc_opt, "Tumor", "Normal"),
  levels = c("Normal", "Tumor")
)

# ------------------------------------------------------------------
# Console confusion matrices
# ------------------------------------------------------------------

cat("\n=== Pan model — threshold 0.5 ===\n")
print(confusionMatrix(pred_pan_05, truth_ext, positive = "Tumor"))

cat("\n=== Pan model — optimal threshold ===\n")
print(confusionMatrix(pred_pan_opt, truth_ext, positive = "Tumor"))

cat("\n=== CRC model — threshold 0.5 ===\n")
print(confusionMatrix(pred_crc_05, truth_ext, positive = "Tumor"))

cat("\n=== CRC model — optimal threshold ===\n")
print(confusionMatrix(pred_crc_opt, truth_ext, positive = "Tumor"))

# ------------------------------------------------------------------
# ROC plots
# ------------------------------------------------------------------

pROC_pan <- plot_roc_pair(
  roc_cv = roc_pan_cv,
  roc_ext = roc_pan_ext,
  title_txt = "ROC Curve: Pancancer/CRC",
  sub_fmt = "AUC - TCGA: %.3f | Test Set: %.3f"
)

pROC_crc <- plot_roc_pair(
  roc_cv = roc_crc_cv,
  roc_ext = roc_crc_ext,
  title_txt = "ROC Curve: CRC/CRC",
  sub_fmt = "AUC - TCGA: %.3f | Test Set: %.3f"
)

print(pROC_pan)
print(pROC_crc)

# ------------------------------------------------------------------
# Confusion matrix plots
# ------------------------------------------------------------------

pCM_pan_05 <- cm_plot(
  pred_pan_05,
  truth_ext,
  "Confusion Matrix — Pancancer model tested on Korean CRC (threshold 0.5)"
)

pCM_pan_opt <- cm_plot(
  pred_pan_opt,
  truth_ext,
  "Confusion Matrix — Pancancer model tested on Korean CRC (optimal threshold)"
)

pCM_crc_05 <- cm_plot(
  pred_crc_05,
  truth_ext,
  "Confusion Matrix — CRC model tested on Korean CRC (threshold 0.5)"
)

pCM_crc_opt <- cm_plot(
  pred_crc_opt,
  truth_ext,
  "Confusion Matrix — CRC model tested on Korean CRC (optimal threshold)"
)

print(pCM_pan_05)
print(pCM_pan_opt)
print(pCM_crc_05)
print(pCM_crc_opt)

# ------------------------------------------------------------------
# PCA plot of Korean CRC cohort using the pan-cancer 9-gene model
# ------------------------------------------------------------------

pca_res <- prcomp(X_ext, center = TRUE, scale. = TRUE)

pca_df <- as.data.frame(pca_res$x) %>%
  rownames_to_column("Sample") %>%
  mutate(
    RF_Pred = pred_pan_05,
    Tissue = truth_ext
  )

p_pca_crc <- ggplot(pca_df, aes(x = PC1, y = PC2, color = RF_Pred, shape = Tissue)) +
  geom_point(size = 4, alpha = 0.9) +
  scale_color_manual(values = c("Tumor" = "#FF6A6A", "Normal" = "lightskyblue")) +
  labs(
    title = "PCA: Korean CRC Cohort — 9-gene panel",
    subtitle = "Color = pan-cancer RF prediction, Shape = true label",
    x = paste0("PC1 (", round(100 * summary(pca_res)$importance[2, 1], 1), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca_res)$importance[2, 2], 1), "%)")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

print(p_pca_crc)
