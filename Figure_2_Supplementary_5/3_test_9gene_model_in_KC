# External validation of pancancer and kidney-trained 9-gene random forest models
# in CPTAC kidney RNA-seq samples
#
# This script:
# 1. reads matched TCGA tumour-normal DESeq2 objects
# 2. extracts VST-transformed expression for the fixed 9-gene panel
# 3. trains two random forest classifiers:
#    - a pancancer TCGA model
#    - a kidney-specific TCGA model (TCGA-KIRC, TCGA-KICH, TCGA-KIRP)
# 4. evaluates both models by 5-fold cross-validation
# 5. reads downloaded CPTAC kidney STAR-count files
# 6. builds a VST-transformed external expression matrix for the same 9-gene panel
# 7. applies both models to the external cohort
# 8. plots:
#    - ROC curves (pancancer/kidney and kidney/kidney)
#    - confusion matrices at threshold 0.5 and optimal threshold
#    - PCA of the external kidney cohort using the 9-gene panel
#
# Input:
# - dds_*_matched_group.rds files
# - downloaded CPTAC kidney STAR-count files in GDC_kidney_filtered/
#
# Notes:
# - Download the CPTAC kidney STAR-count files first using TCGAbiolinks
# - Edit the paths below before running
# - Ensembl IDs are mapped using biomaRt for TCGA data

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(biomaRt)
  library(randomForest)
  library(caret)
  library(pROC)
  library(ggplot2)
  library(TCGAbiolinks)
})

# ------------------------------------------------------------------
# Input paths
# ------------------------------------------------------------------

tcga_dir <- "path-to-your-dds-files"
cptac_dir <- "path-to-folder-containing-GDC_kidney_filtered"

cptac_counts_dir <- file.path(cptac_dir, "GDC_kidney_filtered")

if (!dir.exists(tcga_dir)) {
  stop("Please set 'tcga_dir' to the folder containing dds_*_matched_group.rds files.")
}

if (!dir.exists(cptac_counts_dir)) {
  stop("Please set 'cptac_dir' so that GDC_kidney_filtered/ exists inside it.")
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

get_proba <- function(model_or_train, newdata) {
  p <- predict(model_or_train, newdata = newdata, type = "prob")

  if (is.data.frame(p)) {
    if ("Tumor" %in% names(p)) return(p$Tumor)
    return(p[[ncol(p)]])
  } else if (is.matrix(p)) {
    colnames(p) <- if (is.null(colnames(p))) c("Normal", "Tumor")[seq_len(ncol(p))] else colnames(p)
    if ("Tumor" %in% colnames(p)) return(p[, "Tumor"])
    return(p[, ncol(p)])
  } else {
    return(as.numeric(p))
  }
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
# Load and process TCGA matched cohorts
# ------------------------------------------------------------------

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_map <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  mart = ensembl
) %>%
  filter(external_gene_name != "") %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

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
# Train pancancer and kidney models
# ------------------------------------------------------------------

ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = TRUE
)

# Pancancer model
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

# Kidney-specific model
kid_data <- tcga_clean %>%
  filter(project %in% c("TCGA-KIRC", "TCGA-KICH", "TCGA-KIRP"))

x_kid <- kid_data[, feature_genes, drop = FALSE]
y_kid <- kid_data$condition_bin

set.seed(123)
rf_cv_kid <- train(
  x = x_kid,
  y = y_kid,
  method = "rf",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5
)

cv_preds_pan <- rf_cv_pan$pred %>%
  filter(mtry == rf_cv_pan$bestTune$mtry)

roc_pan_cv <- roc(cv_preds_pan$obs, cv_preds_pan$Tumor, quiet = TRUE)

cv_preds_kid <- rf_cv_kid$pred %>%
  filter(mtry == rf_cv_kid$bestTune$mtry)

roc_kid_cv <- roc(cv_preds_kid$obs, cv_preds_kid$Tumor, quiet = TRUE)

# ------------------------------------------------------------------
# Load CPTAC kidney counts
# ------------------------------------------------------------------

count_files <- list.files(
  cptac_counts_dir,
  pattern = "\\.rna_seq\\.augmented_star_gene_counts\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(count_files) == 0) {
  stop("No CPTAC STAR-count files were found in GDC_kidney_filtered/.")
}

expr_list <- list()

for (f in count_files) {
  dat <- tryCatch({
    read_tsv(f, skip = 1, show_col_types = FALSE) %>%
      dplyr::select(gene_id, gene_name, count = unstranded) %>%
      mutate(
        gene_id = sub("\\..*", "", gene_id),
        count = as.numeric(count)
      ) %>%
      filter(!is.na(gene_name), !is.na(count)) %>%
      group_by(gene_name) %>%
      summarise(count = sum(count), .groups = "drop")
  }, error = function(e) {
    message("Skipping file: ", f, " due to: ", e$message)
    NULL
  })

  if (!is.null(dat)) {
    sample_id <- basename(dirname(f))
    colnames(dat)[colnames(dat) == "count"] <- sample_id
    expr_list[[sample_id]] <- dat
  }
}

expr_merged <- reduce(expr_list, full_join, by = "gene_name") %>%
  column_to_rownames("gene_name")

# ------------------------------------------------------------------
# Build kidney metadata from downloaded files
# ------------------------------------------------------------------

clinical_data <- GDCquery_clinic(project = "CPTAC-3", type = "clinical")
colnames(clinical_data) <- make.names(colnames(clinical_data), unique = TRUE)

kidney_cases <- clinical_data %>%
  filter(primary_diagnosis %in% c(
    "Renal cell carcinoma, NOS",
    "Papillary renal cell carcinoma",
    "Renal cell carcinoma, chromophobe type"
  )) %>%
  pull(submitter_id) %>%
  unique()

query_expr <- GDCquery(
  project = "CPTAC-3",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

metadata_expr <- getResults(
  query_expr,
  cols = c("file_id", "cases.submitter_id", "sample_type", "file_name")
)

metadata_kidney <- metadata_expr %>%
  filter(cases.submitter_id %in% kidney_cases) %>%
  filter(sample_type %in% c("Primary Tumor", "Solid Tissue Normal")) %>%
  dplyr::select(sample = file_id, condition = sample_type) %>%
  mutate(
    condition = factor(
      ifelse(grepl("Normal", condition), "Normal", "Tumor"),
      levels = c("Normal", "Tumor")
    )
  ) %>%
  distinct(sample, .keep_all = TRUE) %>%
  column_to_rownames("sample")

expr_filtered <- expr_merged[, rownames(metadata_kidney), drop = FALSE]

dds_kidney <- DESeqDataSetFromMatrix(
  countData = round(expr_filtered),
  colData = metadata_kidney,
  design = ~ condition
)

vsd_kidney <- vst(dds_kidney, blind = TRUE)
vst_counts_kidney <- assay(vsd_kidney)

available_genes <- intersect(feature_genes, rownames(vst_counts_kidney))
missing_genes <- setdiff(feature_genes, available_genes)

if (length(missing_genes) > 0) {
  message("Missing genes in kidney test set: ", paste(missing_genes, collapse = ", "))
}

test_expr <- t(vst_counts_kidney[available_genes, , drop = FALSE]) %>%
  as.data.frame()

for (g in missing_genes) {
  test_expr[[g]] <- NA_real_
}

test_expr <- test_expr[, feature_genes, drop = FALSE]

for (g in colnames(test_expr)) {
  if (anyNA(test_expr[[g]])) {
    med_g <- median(test_expr[[g]], na.rm = TRUE)
    test_expr[[g]][is.na(test_expr[[g]])] <- if (is.finite(med_g)) med_g else 0
  }
}

truth_labels <- metadata_kidney$condition

# ------------------------------------------------------------------
# Predictions on kidney test set
# ------------------------------------------------------------------

rf_pred_pan <- predict(rf_cv_pan, newdata = test_expr)
rf_probs_pan <- predict(rf_cv_pan, newdata = test_expr, type = "prob")

rf_pred_kidney <- predict(rf_cv_kid, newdata = test_expr)
rf_probs_kidney <- predict(rf_cv_kid, newdata = test_expr, type = "prob")

roc_pan_test <- roc(truth_labels, rf_probs_pan[, "Tumor"], quiet = TRUE)
roc_kid_test <- roc(truth_labels, rf_probs_kidney[, "Tumor"], quiet = TRUE)

# ------------------------------------------------------------------
# ROC plots
# ------------------------------------------------------------------

p_roc_pan <- ggplot(bind_rows(
  data.frame(
    fpr = 1 - roc_pan_cv$specificities,
    tpr = roc_pan_cv$sensitivities,
    model = "TCGA Pan CV"
  ),
  data.frame(
    fpr = 1 - roc_pan_test$specificities,
    tpr = roc_pan_test$sensitivities,
    model = "Kidney Test"
  )
), aes(x = fpr, y = tpr, color = model)) +
  geom_line(linewidth = 1.4) +
  geom_abline(linetype = "dashed", color = "gray70") +
  scale_color_manual(values = c("TCGA Pan CV" = "coral1", "Kidney Test" = "#27408B")) +
  labs(
    title = "ROC: Pan-Trained Model",
    subtitle = sprintf("AUC TCGA: %.3f | Test: %.3f", auc(roc_pan_cv), auc(roc_pan_test)),
    x = "False Positive Rate",
    y = "True Positive Rate",
    color = "Dataset"
  ) +
  theme_minimal(base_size = 14)

p_roc_kid <- ggplot(bind_rows(
  data.frame(
    fpr = 1 - roc_kid_cv$specificities,
    tpr = roc_kid_cv$sensitivities,
    model = "TCGA Kidney CV"
  ),
  data.frame(
    fpr = 1 - roc_kid_test$specificities,
    tpr = roc_kid_test$sensitivities,
    model = "Kidney Test"
  )
), aes(x = fpr, y = tpr, color = model)) +
  geom_line(linewidth = 1.4) +
  geom_abline(linetype = "dashed", color = "gray70") +
  scale_color_manual(values = c("TCGA Kidney CV" = "coral1", "Kidney Test" = "#27408B")) +
  labs(
    title = "ROC: Kidney-Trained Model",
    subtitle = sprintf("AUC TCGA: %.3f | Test: %.3f", auc(roc_kid_cv), auc(roc_kid_test)),
    x = "False Positive Rate",
    y = "True Positive Rate",
    color = "Dataset"
  ) +
  theme_minimal(base_size = 14)

print(p_roc_pan)
print(p_roc_kid)

# ------------------------------------------------------------------
# Confusion matrices at threshold 0.5
# ------------------------------------------------------------------

p_cm_pan <- cm_plot(
  rf_pred_pan,
  truth_labels,
  "Confusion Matrix: Pan-Trained Model on Kidney Test"
)

p_cm_kid <- cm_plot(
  rf_pred_kidney,
  truth_labels,
  "Confusion Matrix: Kidney-Trained Model on Kidney Test"
)

print(p_cm_pan)
print(p_cm_kid)

# ------------------------------------------------------------------
# Confusion matrices at optimal threshold
# ------------------------------------------------------------------

thr_pan_opt <- best_thr_bacc(truth_labels, rf_probs_pan[, "Tumor"])
pred_pan_opt <- factor(
  ifelse(rf_probs_pan[, "Tumor"] >= thr_pan_opt, "Tumor", "Normal"),
  levels = c("Normal", "Tumor")
)

thr_kid_opt <- best_thr_bacc(truth_labels, rf_probs_kidney[, "Tumor"])
pred_kid_opt <- factor(
  ifelse(rf_probs_kidney[, "Tumor"] >= thr_kid_opt, "Tumor", "Normal"),
  levels = c("Normal", "Tumor")
)

p_cm_pan_opt <- cm_plot(
  pred_pan_opt,
  truth_labels,
  sprintf("Confusion Matrix: Pan-Trained Model (thr = %.3f)", thr_pan_opt)
)

p_cm_kid_opt <- cm_plot(
  pred_kid_opt,
  truth_labels,
  sprintf("Confusion Matrix: Kidney-Trained Model (thr = %.3f)", thr_kid_opt)
)

print(p_cm_pan_opt)
print(p_cm_kid_opt)

# ------------------------------------------------------------------
# PCA of kidney cohort colored by kidney-trained RF prediction
# ------------------------------------------------------------------

expr_kidney_full <- t(vst_counts_kidney) %>%
  as.data.frame()

common_samples <- intersect(rownames(expr_kidney_full), rownames(metadata_kidney))
expr_kidney_full <- expr_kidney_full[common_samples, , drop = FALSE]
metadata_kidney_pca <- metadata_kidney[common_samples, , drop = FALSE]

expr_kidney_filtered <- expr_kidney_full[, apply(expr_kidney_full, 2, var) > 0, drop = FALSE]

pca_kidney <- prcomp(expr_kidney_filtered, center = TRUE, scale. = TRUE)

names(rf_pred_kidney) <- rownames(test_expr)

pca_df <- as.data.frame(pca_kidney$x) %>%
  rownames_to_column("sample") %>%
  left_join(
    metadata_kidney_pca %>%
      rownames_to_column("sample"),
    by = "sample"
  ) %>%
  mutate(
    rf_prediction = rf_pred_kidney[sample],
    tissue = condition
  )

p_pca_kidney <- ggplot(pca_df, aes(x = PC1, y = PC2, color = rf_prediction, shape = tissue)) +
  geom_point(size = 4, alpha = 0.9) +
  scale_color_manual(values = c("Tumor" = "#FF6A6A", "Normal" = "lightskyblue")) +
  labs(
    title = "PCA of Kidney Samples Colored by Kidney RF Prediction",
    x = paste0("PC1 (", round(100 * summary(pca_kidney)$importance[2, 1], 1), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca_kidney)$importance[2, 2], 1), "%)")
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

print(p_pca_kidney)
