# External validation of pancancer and lung-trained 9-gene random forest models
# in GSE87410 lung RNA-seq samples
#
# This script:
# 1. reads matched TCGA tumour-normal DESeq2 objects
# 2. extracts VST-transformed expression for the fixed 9-gene panel
# 3. trains two random forest classifiers:
#    - a pancancer TCGA model
#    - a lung-specific TCGA model (TCGA-LUAD and TCGA-LUSC)
# 4. evaluates both models by 5-fold cross-validation
# 5. reads the external GSE87410 count matrix
# 6. builds a VST-transformed external expression matrix for the same 9-gene panel
# 7. applies both models to the external cohort
# 8. plots:
#    - ROC curves (pancancer/lung and lung/lung)
#    - confusion matrices at threshold 0.5
#    - PCA of the external lung cohort colored by each model's prediction
#
# Input:
# - dds_*_matched_group.rds files
# - GSE87410_raw_counts_GRCh38.p13_NCBI.tsv.gz
# - GSM_sample_map.csv
#
# Notes:
# - Edit the paths below before running
# - Ensembl IDs are mapped using biomaRt for TCGA data
# - External gene IDs are converted from Entrez ID to gene symbol using org.Hs.eg.db

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(biomaRt)
  library(randomForest)
  library(caret)
  library(pROC)
  library(ggplot2)
  library(patchwork)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# ------------------------------------------------------------------
# Input paths
# ------------------------------------------------------------------

tcga_dir <- "path-to-your-dds-files"
gse_dir  <- "path-to-GSE87410-folder"

counts_file <- file.path(gse_dir, "GSE87410_raw_counts_GRCh38.p13_NCBI.tsv.gz")
map_file    <- file.path(gse_dir, "GSM_sample_map.csv")

if (!dir.exists(tcga_dir)) {
  stop("Please set 'tcga_dir' to the folder containing dds_*_matched_group.rds files.")
}

if (!dir.exists(gse_dir)) {
  stop("Please set 'gse_dir' to the folder containing the GSE87410 files.")
}

if (!file.exists(counts_file)) {
  stop("Could not find: ", counts_file)
}

if (!file.exists(map_file)) {
  stop("Could not find: ", map_file)
}

# ------------------------------------------------------------------
# Fixed 9-gene panel
# ------------------------------------------------------------------

feature_genes <- c(
  "TTLL4", "TTLL5", "TTLL6", "TTLL7", "TTLL11",
  "AGBL1", "AGBL3", "AGBL4", "AGBL5"
)

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

cm_plot <- function(pred, truth, title_txt) {
  tb <- table(Prediction = pred, Reference = truth)
  df <- as.data.frame(tb)
  colnames(df) <- c("Prediction", "Reference", "Freq")

  ggplot(df, aes(x = Reference, y = Prediction, fill = Freq)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Freq), size = 5, fontface = "bold") +
    scale_fill_gradient(low = "white", high = "aquamarine4") +
    labs(title = title_txt, x = "Actual Label", y = "Predicted Label") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

roc_df_from <- function(roc_obj, label) {
  data.frame(
    fpr   = 1 - roc_obj$specificities,
    tpr   = roc_obj$sensitivities,
    model = label
  )
}

plot_roc_pair <- function(roc_cv, roc_ext, title_txt, sub_fmt) {
  ggplot(dplyr::bind_rows(
    roc_df_from(roc_cv, "TCGA CV"),
    roc_df_from(roc_ext, "External")
  ), aes(x = fpr, y = tpr, color = model)) +
    geom_line(linewidth = 1.3) +
    geom_abline(linetype = "dashed", color = "gray70") +
    scale_color_manual(values = c("TCGA CV" = "coral1", "External" = "slateblue4")) +
    labs(
      title = title_txt,
      subtitle = sprintf(sub_fmt, pROC::auc(roc_cv), pROC::auc(roc_ext)),
      x = "False Positive Rate",
      y = "True Positive Rate",
      color = NULL
    ) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom"
    )
}

# ------------------------------------------------------------------
# TCGA: read dds_* .rds -> VST -> map symbols -> build wide panel matrix
# ------------------------------------------------------------------

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_map <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  mart = ensembl
) %>%
  dplyr::filter(external_gene_name != "") %>%
  dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)

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

  matched <- intersect(feature_genes, rownames(norm_counts))
  if (!length(matched)) next

  expr <- t(norm_counts[matched, , drop = FALSE]) %>% as.data.frame()

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
  expr$project   <- project
  all_data[[length(all_data) + 1]] <- expr
}

tcga_df <- bind_rows(all_data)

tcga_clean <- tcga_df %>%
  dplyr::select(all_of(feature_genes), condition, sample, project) %>%
  dplyr::filter(if_all(all_of(feature_genes), ~ !is.na(.))) %>%
  dplyr::mutate(
    condition_bin = factor(
      ifelse(condition == "Primary Tumor", "Tumor", "Normal"),
      levels = c("Normal", "Tumor")
    )
  )

# ------------------------------------------------------------------
# Train BOTH models (same 5-fold CV, 9-gene panel)
# ------------------------------------------------------------------

ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = TRUE
)

# Pan-cancer
set.seed(123)
rf_cv_pan <- train(
  x = tcga_clean[, feature_genes, drop = FALSE],
  y = tcga_clean$condition_bin,
  method = "rf",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5
)

pred_pan_cv <- rf_cv_pan$pred
if (!is.null(rf_cv_pan$bestTune)) {
  for (nm in names(rf_cv_pan$bestTune)) {
    pred_pan_cv <- pred_pan_cv[pred_pan_cv[[nm]] == rf_cv_pan$bestTune[[nm]], , drop = FALSE]
  }
}
roc_tcga_pan <- pROC::roc(
  pred_pan_cv$obs,
  pred_pan_cv$Tumor,
  levels = c("Normal", "Tumor"),
  quiet = TRUE
)

cat(sprintf("\nRF Pan CV AUC: %.3f\n", pROC::auc(roc_tcga_pan)))

# Lung-only (LUAD + LUSC)
tcga_lung <- tcga_clean %>%
  dplyr::filter(project %in% c("TCGA-LUAD", "TCGA-LUSC"))

set.seed(123)
rf_cv_lung <- train(
  x = tcga_lung[, feature_genes, drop = FALSE],
  y = tcga_lung$condition_bin,
  method = "rf",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5
)

pred_lung_cv <- rf_cv_lung$pred
if (!is.null(rf_cv_lung$bestTune)) {
  for (nm in names(rf_cv_lung$bestTune)) {
    pred_lung_cv <- pred_lung_cv[pred_lung_cv[[nm]] == rf_cv_lung$bestTune[[nm]], , drop = FALSE]
  }
}
roc_tcga_lung <- pROC::roc(
  pred_lung_cv$obs,
  pred_lung_cv$Tumor,
  levels = c("Normal", "Tumor"),
  quiet = TRUE
)

cat(sprintf("RF Lung-only CV AUC: %.3f\n", pROC::auc(roc_tcga_lung)))

# ------------------------------------------------------------------
# External test: GSE87410 -> VST -> Entrez -> SYMBOL
# ------------------------------------------------------------------

raw_counts <- readr::read_tsv(counts_file)
gene_ids   <- raw_counts[[1]]
counts_mat <- as.data.frame(raw_counts[, -1])
rownames(counts_mat) <- gene_ids

gsm_map <- readr::read_csv(map_file, col_names = c("GSM", "Label")) %>%
  tidyr::separate(Label, into = c("Project", "Patient", "Tissue"), sep = "_", remove = FALSE) %>%
  dplyr::mutate(Tissue = ifelse(Tissue %in% c("T", "N"), Tissue, NA)) %>%
  dplyr::filter(!is.na(Tissue))

valid_samples <- intersect(colnames(counts_mat), gsm_map$GSM)
counts_mat <- counts_mat[, valid_samples, drop = FALSE]
counts_mat[] <- lapply(counts_mat, as.numeric)

sample_info <- gsm_map %>%
  dplyr::filter(GSM %in% valid_samples) %>%
  dplyr::select(GSM, Project, Tissue) %>%
  tibble::column_to_rownames("GSM")

dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData = sample_info,
  design = ~ Project + Tissue
)

vst_expr <- assay(vst(dds, blind = TRUE))

# Entrez -> SYMBOL
entrez_ids <- rownames(vst_expr)
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = entrez_ids,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

gene_symbols_clean <- ifelse(is.na(gene_symbols) | gene_symbols == "", entrez_ids, gene_symbols)
rownames(vst_expr) <- gene_symbols_clean
vst_expr <- vst_expr[rownames(vst_expr) != "", , drop = FALSE]

# Keep LUAD/LUSC samples only
keep <- rownames(sample_info)[sample_info$Project %in% c("LUAD", "LUSC")]
vst_expr <- vst_expr[, keep, drop = FALSE]
sample_info <- sample_info[keep, , drop = FALSE]

available <- intersect(feature_genes, rownames(vst_expr))
missing <- setdiff(feature_genes, available)
if (length(missing)) {
  stop("Missing feature genes: ", paste(missing, collapse = ", "))
}

patient_expr_df <- t(vst_expr[feature_genes, , drop = FALSE]) %>%
  as.data.frame()

truth <- factor(
  ifelse(sample_info$Tissue == "T", "Tumor", "Normal"),
  levels = c("Normal", "Tumor")
)

# ------------------------------------------------------------------
# Predictions (Pan + Lung models) at standard threshold 0.5
# ------------------------------------------------------------------

# Pan
probs_pan <- predict(rf_cv_pan, newdata = patient_expr_df, type = "prob")[, "Tumor"]
pred_pan_05 <- factor(
  ifelse(probs_pan >= 0.5, "Tumor", "Normal"),
  levels = c("Normal", "Tumor")
)
roc_ext_pan <- pROC::roc(
  truth,
  probs_pan,
  levels = c("Normal", "Tumor"),
  quiet = TRUE
)

# Lung-only
probs_lung <- predict(rf_cv_lung, newdata = patient_expr_df, type = "prob")[, "Tumor"]
pred_lung_05 <- factor(
  ifelse(probs_lung >= 0.5, "Tumor", "Normal"),
  levels = c("Normal", "Tumor")
)
roc_ext_lung <- pROC::roc(
  truth,
  probs_lung,
  levels = c("Normal", "Tumor"),
  quiet = TRUE
)

# ------------------------------------------------------------------
# Quick numeric evaluation
# ------------------------------------------------------------------

cat("\n=== External GSE87410 — Pan model ===\n")
cat(sprintf("AUC = %.3f\n", pROC::auc(roc_ext_pan)))
print(table(Predicted = pred_pan_05, Actual = truth))

cat("\n=== External GSE87410 — Lung-only model ===\n")
cat(sprintf("AUC = %.3f\n", pROC::auc(roc_ext_lung)))
print(table(Predicted = pred_lung_05, Actual = truth))

eval_df <- tibble(
  Sample = rownames(patient_expr_df),
  Project = sample_info$Project[rownames(patient_expr_df)],
  Truth = truth,
  PredPan = pred_pan_05,
  PredLung = pred_lung_05
) %>%
  group_by(Project) %>%
  summarise(
    Acc_Pan = mean(PredPan == Truth),
    Acc_Lung = mean(PredLung == Truth),
    N = n(),
    .groups = "drop"
  )

print(eval_df)

# ------------------------------------------------------------------
# Plots: CMs and paired ROCs
# ------------------------------------------------------------------

# Confusion matrices (0.5 only)
p_cm_pan_05 <- cm_plot(
  pred_pan_05,
  truth,
  "Confusion Matrix: Pan-Trained Model on Lung Test (thr = 0.5)"
)

p_cm_lung_05 <- cm_plot(
  pred_lung_05,
  truth,
  "Confusion Matrix: Lung-Trained Model on Lung Test (thr = 0.5)"
)

# ROC comparisons (CV vs External)
p_roc_pan <- plot_roc_pair(
  roc_cv = roc_tcga_pan,
  roc_ext = roc_ext_pan,
  title_txt = "ROC — Pan: TCGA CV vs GSE87410",
  sub_fmt = "AUC (TCGA Pan CV) = %.3f | AUC (Test) = %.3f"
)

p_roc_lung <- plot_roc_pair(
  roc_cv = roc_tcga_lung,
  roc_ext = roc_ext_lung,
  title_txt = "ROC Curve: Lung/Lung",
  sub_fmt = "AUC: TCGA Lung = %.3f | Test Set = %.3f"
)

print(p_roc_pan)
print(p_roc_lung)
print(p_cm_pan_05)
print(p_cm_lung_05)

# ------------------------------------------------------------------
# PCA (colored by each model's prediction)
# ------------------------------------------------------------------

pca_input <- patient_expr_df[, apply(patient_expr_df, 2, var) > 0, drop = FALSE]
pca_res <- prcomp(pca_input, center = TRUE, scale. = TRUE)

pca_df <- as.data.frame(pca_res$x) %>%
  tibble::rownames_to_column("Sample") %>%
  dplyr::left_join(
    tibble(
      Sample = rownames(patient_expr_df),
      PredPan = pred_pan_05,
      PredLung = pred_lung_05,
      Truth = truth
    ),
    by = "Sample"
  )

g1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = PredPan, shape = Truth)) +
  geom_point(size = 4, alpha = 0.9) +
  scale_color_manual(values = c("Tumor" = "#FF6A6A", "Normal" = "lightskyblue")) +
  labs(
    title = "PCA of Samples — Pan RF Prediction",
    x = paste0("PC1 (", round(100 * summary(pca_res)$importance[2, 1], 1), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca_res)$importance[2, 2], 1), "%)")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

g2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = PredLung, shape = Truth)) +
  geom_point(size = 4, alpha = 0.9) +
  scale_color_manual(values = c("Tumor" = "#FF6A6A", "Normal" = "lightskyblue")) +
  labs(
    title = "PCA — colored by Lung-only model",
    subtitle = "Shape = ground truth",
    x = paste0("PC1 (", round(100 * summary(pca_res)$importance[2, 1], 1), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca_res)$importance[2, 2], 1), "%)")
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

print(g1 + g2)
