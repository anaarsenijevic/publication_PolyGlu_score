# External validation of TCGA-trained 15-gene random forest classifier in BRCA patient samples
#
#
#
# Input:
# - Salmon quant.sf files for the patient cohort
# - Gencode GTF annotation file
# - rf_15gene_tcga_model.rds
#
# Notes:
# - Run the TCGA training script first
# - Edit the paths below before running
# - Ensembl IDs are mapped to gene symbols using org.Hs.eg.db
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
})

# ------------------------------------------------------------------
# Input paths
# ------------------------------------------------------------------

annotation_dir <- "path-to-annotation-folder"
salmon_dir     <- "path-to-salmon-folder"
model_dir      <- "path-to-output-folder"

gtf_file   <- file.path(annotation_dir, "annotations", "gencode.v29.annotation.gtf")
model_file <- file.path(model_dir, "rf_15gene_tcga_model.rds")

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
# 1. PREPROCESS BRCA PATIENT DATA
# ==================================================================

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
patient_mat <- assay(vsd_patient)

ens_ids_patient <- sub("\\..*", "", rownames(patient_mat))
symbols_patient <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = ens_ids_patient,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

patient_sym <- aggregate_dup_symbols(patient_mat, symbols_patient, mode = "mean")

available_genes <- intersect(feature_genes, rownames(patient_sym))
missing_genes <- setdiff(feature_genes, available_genes)

patient_expr_df <- t(patient_sym[available_genes, , drop = FALSE]) %>%
  as.data.frame()

for (g in missing_genes) {
  patient_expr_df[[g]] <- NA_real_
}

patient_expr_df <- patient_expr_df[, feature_genes, drop = FALSE]

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

# ==================================================================
# 2. LOAD MODEL AND TEST IN BRCA
# ==================================================================

rf_15 <- readRDS(model_file)

patient_probs <- predict(rf_15, newdata = patient_expr_df, type = "prob")
patient_preds <- predict(rf_15, newdata = patient_expr_df)

roc_patient <- roc(truth_labels, patient_probs[, "Tumor"], quiet = TRUE)
auc_patient <- as.numeric(auc(roc_patient))

roc_patient_df <- data.frame(
  fpr = 1 - roc_patient$specificities,
  tpr = roc_patient$sensitivities
)

ggplot(roc_patient_df, aes(x = fpr, y = tpr)) +
  geom_line(linewidth = 1.4, color = "#8B2500") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
  labs(
    title = "ROC Curve: 15-Gene TCGA Model in External BRCA Dataset",
    subtitle = sprintf("AUC = %.3f", auc_patient),
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

p_cm_patient <- plot_cm(
  pred = patient_preds,
  truth = truth_labels,
  title = "Confusion Matrix: 15-Gene TCGA Model in External BRCA Dataset"
)

print(p_cm_patient)
