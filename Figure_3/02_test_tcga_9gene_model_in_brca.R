# External validation of TCGA-trained 9-gene random forest classifier in BRCA patient samples
#
# This script:
# 1. reads Salmon quantifications from the external BRCA patient cohort
# 2. builds a 9-gene VST-normalized expression matrix
# 3. loads the TCGA-trained 9-gene random forest model and CV object
# 4. plots per-project TCGA PCA using the 9-gene panel
# 5. predicts tumour vs normal labels in the external BRCA dataset
# 6. plots ROC curve and confusion matrix
#
# Input:
# - Salmon quant.sf files for the patient cohort
# - Gencode GTF annotation file
# - rf_9gene_tcga_model.rds
# - rf_9gene_tcga_cv.rds
# - tcga_9gene_expression.rds
#
# Notes:
# - Run the TCGA training script first
# - Edit the paths below before running
# - Duplicate symbols are collapsed by mean expression

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(tidyverse)
  library(AnnotationDbi)
  library(GenomicFeatures)
  library(org.Hs.eg.db)
  library(randomForest)
  library(caret)
  library(pROC)
  library(ggplot2)
  library(patchwork)
})

# ------------------------------------------------------------------
# Input paths
# ------------------------------------------------------------------

annotation_dir <- "path-to-annotation-folder"
salmon_dir     <- "path-to-salmon-folder"
model_dir      <- "path-to-output-folder"

gtf_file    <- file.path(annotation_dir, "annotations", "gencode.v29.annotation.gtf")
model_file  <- file.path(model_dir, "rf_9gene_tcga_model.rds")
cv_file     <- file.path(model_dir, "rf_9gene_tcga_cv.rds")
tcga_file   <- file.path(model_dir, "tcga_9gene_expression.rds")

if (!dir.exists(annotation_dir)) {
  stop("Please set 'annotation_dir' to the folder containing the GTF annotation.")
}

if (!dir.exists(salmon_dir)) {
  stop("Please set 'salmon_dir' to the folder containing the Salmon quantification folders.")
}

if (!file.exists(gtf_file)) {
  stop("GTF file not found: ", gtf_file)
}

if (!file.exists(model_file)) {
  stop("Missing trained model file: ", model_file)
}

if (!file.exists(cv_file)) {
  stop("Missing cross-validation file: ", cv_file)
}

if (!file.exists(tcga_file)) {
  stop("Missing TCGA expression file: ", tcga_file)
}

# ------------------------------------------------------------------
# Fixed 9-gene panel
# ------------------------------------------------------------------

genes_9 <- c(
  "TTLL4", "TTLL5", "TTLL6", "TTLL7", "TTLL11",
  "AGBL1", "AGBL2", "AGBL3", "AGBL5"
)

# ------------------------------------------------------------------
# Sample metadata
# ------------------------------------------------------------------

samples_patient <- tibble::tribble(
  ~sample, ~tissue, ~group, ~path,
  "1NAT",  "NAT", "non-invasive", "salmon_1NAT/quant.sf",
  "2PT",   "PT",  "non-invasive", "salmon_2PT/quant.sf",
  "3NAT",  "NAT", "non-invasive", "salmon_3NAT/quant.sf",
  "4PT",   "PT",  "non-invasive", "salmon_4PT/quant.sf",
  "5NAT",  "NAT", "metastatic",   "salmon_5NAT/quant.sf",
  "6PT",   "PT",  "metastatic",   "salmon_6PT/quant.sf",
  "7NAT",  "NAT", "metastatic",   "salmon_7NAT/quant.sf",
  "8PT",   "PT",  "metastatic",   "salmon_8PT/quant.sf",
  "9NAT",  "NAT", "metastatic",   "salmon_9NAT/quant.sf",
  "10PT",  "PT",  "metastatic",   "salmon_10PT/quant.sf",
  "11NAT", "NAT", "TNBC",         "salmon_11NAT/quant.sf",
  "12PT",  "PT",  "TNBC",         "salmon_12PT/quant.sf",
  "13NAT", "NAT", "TNBC",         "salmon_13NAT/quant.sf",
  "14PT",  "PT",  "TNBC",         "salmon_14PT/quant.sf",
  "15NAT", "NAT", "TNBC",         "salmon_15NAT/quant.sf",
  "16PT",  "PT",  "TNBC",         "salmon_16PT/quant.sf",
  "17NAT", "NAT", "non-invasive", "salmon_17NAT/quant.sf",
  "18PT",  "PT",  "non-invasive", "salmon_18PT/quant.sf"
)

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

make_pca_plot <- function(expr_df, meta_df, color_var, title_txt, color_values) {
  sample_ids <- intersect(rownames(expr_df), meta_df$sample)

  expr_df <- expr_df[sample_ids, , drop = FALSE]
  meta_df <- meta_df %>%
    filter(sample %in% sample_ids) %>%
    slice(match(sample_ids, sample))

  keep <- apply(expr_df, 2, var, na.rm = TRUE) > 0
  expr_df <- expr_df[, keep, drop = FALSE]

  if (nrow(expr_df) < 2 || ncol(expr_df) < 2) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = "Not enough data for PCA", size = 5) +
        theme_void() +
        labs(title = title_txt)
    )
  }

  pca_res <- prcomp(expr_df, center = TRUE, scale. = TRUE)
  var_exp <- 100 * (pca_res$sdev^2 / sum(pca_res$sdev^2))

  pca_df <- as.data.frame(pca_res$x[, 1:2, drop = FALSE]) %>%
    rownames_to_column("sample") %>%
    left_join(meta_df, by = "sample")

  ggplot(pca_df, aes(x = PC1, y = PC2, color = .data[[color_var]])) +
    geom_point(size = 4, alpha = 1) +
    scale_color_manual(values = color_values) +
    labs(
      title = title_txt,
      x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
      y = paste0("PC2 (", round(var_exp[2], 1), "%)")
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "bottom"
    )
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

# ------------------------------------------------------------------
# Load saved TCGA objects
# ------------------------------------------------------------------

rf_9 <- readRDS(model_file)
rf_9_cv <- readRDS(cv_file)
tcga_df <- readRDS(tcga_file)

# ------------------------------------------------------------------
# Per-project TCGA PCA using the 9-gene panel
# ------------------------------------------------------------------

tcga_expr_9 <- tcga_df[, genes_9, drop = FALSE]
rownames(tcga_expr_9) <- tcga_df$sample

tcga_meta <- tcga_df %>%
  dplyr::select(sample, project, condition, condition_bin)

condition_colors <- c("Normal" = "lightskyblue", "Tumor" = "#FF6A6A")

projects <- sort(unique(tcga_meta$project))

project_pca_list <- lapply(projects, function(prj) {
  meta_sub <- tcga_meta %>% filter(project == prj)
  expr_sub <- tcga_expr_9[meta_sub$sample, , drop = FALSE]

  make_pca_plot(
    expr_df = expr_sub,
    meta_df = meta_sub,
    color_var = "condition_bin",
    title_txt = prj,
    color_values = condition_colors
  )
})

names(project_pca_list) <- projects

combined_project_pca <- wrap_plots(project_pca_list) +
  plot_annotation(
    title = "TCGA PCA by Project (9-Gene Model)",
    subtitle = "Tumour vs normal separation using the 9-gene model"
  )

print(combined_project_pca)

# ------------------------------------------------------------------
# Process BRCA patient cohort
# ------------------------------------------------------------------

txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")

tx2gene <- AnnotationDbi::select(
  txdb,
  keys = keys(txdb, keytype = "TXNAME"),
  columns = "GENEID",
  keytype = "TXNAME"
) %>%
  mutate(
    TXNAME = sub("\\..*", "", TXNAME),
    GENEID = sub("\\..*", "", GENEID)
  ) %>%
  distinct()

files <- file.path(salmon_dir, samples_patient$path)
names(files) <- samples_patient$sample

missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0) {
  stop("The following Salmon files were not found:\n", paste(missing_files, collapse = "\n"))
}

txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE
)

sample_table <- samples_patient %>%
  dplyr::select(sample, tissue, group) %>%
  column_to_rownames("sample")

dds_patient <- DESeqDataSetFromTximport(
  txi = txi,
  colData = sample_table,
  design = ~ tissue
)

vsd_patient <- vst(dds_patient, blind = TRUE)
vst_expr <- assay(vsd_patient)

ens_ids <- sub("\\..*", "", rownames(vst_expr))
symbols <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = ens_ids,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

vst_expr_sym <- aggregate_dup_symbols(vst_expr, symbols, mode = "mean")

available_genes <- intersect(genes_9, rownames(vst_expr_sym))
patient_expr_df <- t(vst_expr_sym[available_genes, , drop = FALSE]) %>%
  as.data.frame()

for (g in setdiff(genes_9, colnames(patient_expr_df))) {
  patient_expr_df[[g]] <- NA_real_
}
patient_expr_df <- patient_expr_df[, genes_9, drop = FALSE]

for (g in colnames(patient_expr_df)) {
  if (anyNA(patient_expr_df[[g]])) {
    med_g <- median(patient_expr_df[[g]], na.rm = TRUE)
    patient_expr_df[[g]][is.na(patient_expr_df[[g]])] <- if (is.finite(med_g)) med_g else 0
  }
}

truth_labels <- factor(
  ifelse(samples_patient$tissue == "PT", "Tumor", "Normal"),
  levels = c("Normal", "Tumor")
)

# ------------------------------------------------------------------
# ROC curve: TCGA CV vs BRCA patient cohort
# ------------------------------------------------------------------

cv_preds <- rf_9_cv$pred %>%
  filter(mtry == rf_9_cv$bestTune$mtry)

roc_tcga <- roc(cv_preds$obs, cv_preds$Tumor, quiet = TRUE)

patient_probs <- predict(rf_9, newdata = patient_expr_df, type = "prob")
patient_preds <- predict(rf_9, newdata = patient_expr_df)

roc_patient <- roc(truth_labels, patient_probs[, "Tumor"], quiet = TRUE)

roc_df <- bind_rows(
  data.frame(
    fpr = 1 - roc_tcga$specificities,
    tpr = roc_tcga$sensitivities,
    dataset = "TCGA CV (9-gene)"
  ),
  data.frame(
    fpr = 1 - roc_patient$specificities,
    tpr = roc_patient$sensitivities,
    dataset = "Patient (9-gene)"
  )
)

auc_tcga <- as.numeric(auc(roc_tcga))
auc_patient <- as.numeric(auc(roc_patient))

p_roc <- ggplot(roc_df, aes(x = fpr, y = tpr, color = dataset)) +
  geom_line(linewidth = 1.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
  scale_color_manual(values = c("TCGA CV (9-gene)" = "coral1", "Patient (9-gene)" = "#27408B")) +
  labs(
    title = "ROC Curve: 9-Gene Model",
    subtitle = sprintf("AUC — TCGA: %.3f | Patient: %.3f", auc_tcga, auc_patient),
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

print(p_roc)

# ------------------------------------------------------------------
# Confusion matrix: BRCA patient cohort
# ------------------------------------------------------------------

p_cm_9 <- plot_cm(
  pred = patient_preds,
  truth = truth_labels,
  title = "Confusion Matrix — 9-Gene Model",
  high_col = "#00B2EE"
)

print(p_cm_9)
