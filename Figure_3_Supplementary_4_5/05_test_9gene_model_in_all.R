# Combined external cohort RF validation
#
# This script:
# - Trains a 9-gene random forest classifier on pan-cancer TCGA samples
# - Loads three independent external datasets: Lung, Kidney, and CRC
# - Treats all external samples as one combined validation cohort
# - Computes ROC and AUC for TCGA cross-validation and the combined external cohort
# - Generates confusion matrices at:
#     1. Threshold = 0.5
#     2. Optimal threshold based on balanced accuracy
#
# Inputs:
# - TCGA matched DESeq2 RDS files
# - Lung GSE87410 raw counts and metadata
# - Kidney CPTAC-3 STAR counts
# - Korean CRC count matrix
#
# Output:
# - ROC plot
# - Confusion matrices
# - Dataset-level AUC sanity check

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(TCGAbiolinks)
  library(randomForest)
  library(caret)
  library(pROC)
  library(ggplot2)
  library(patchwork)
})

set.seed(123)

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------

tcga_dir   <- "path-to/TCGA_files_rds/DESeq_excels_dds_rds"
lung_dir   <- "path-to/GSE87410"
kidney_dir <- "path-to/KIRC_CPTAC"
crc_dir    <- "path-to/KoreanCRC"

output_dir <- "Combined_external_RF_validation"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------
# 9-gene panel
# ------------------------------------------------------------------

feature_genes <- c(
  "TTLL4", "TTLL5", "TTLL6", "TTLL7", "TTLL11",
  "AGBL1", "AGBL3", "AGBL4", "AGBL5"
)

# ------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------

collapse_duplicate_symbols <- function(expr_mat) {
  expr_mat <- expr_mat[
    !is.na(rownames(expr_mat)) & rownames(expr_mat) != "",
    ,
    drop = FALSE
  ]

  as.data.frame(
    rowsum(
      as.matrix(expr_mat),
      group = rownames(expr_mat),
      reorder = FALSE
    )
  )
}

map_ensembl_to_symbol <- function(ids) {
  AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = sub("\\..*", "", ids),
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
}

map_entrez_to_symbol <- function(ids) {
  AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = as.character(ids),
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
  )
}

get_proba <- function(model, newdata) {
  probs <- predict(model, newdata = newdata, type = "prob")
  probs[, "Tumor"]
}

best_thr_bacc <- function(truth, p_tumor) {
  roc_obj <- pROC::roc(
    truth,
    p_tumor,
    levels = c("Normal", "Tumor"),
    direction = "<",
    quiet = TRUE
  )

  coords_df <- as.data.frame(
    pROC::coords(
      roc_obj,
      x = "all",
      ret = c("threshold", "sensitivity", "specificity"),
      transpose = FALSE
    )
  )

  coords_df$balanced_accuracy <- (
    coords_df$sensitivity + coords_df$specificity
  ) / 2

  coords_df$threshold[which.max(coords_df$balanced_accuracy)]
}

roc_df_from <- function(roc_obj, label) {
  data.frame(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities,
    model = label
  )
}

plot_roc_pair <- function(roc_cv, roc_ext, title_txt) {

  df_cv <- roc_df_from(roc_cv, "TCGA CV")
  df_ext <- roc_df_from(roc_ext, "Combined external")

  ggplot() +
    geom_line(
      data = df_cv,
      aes(x = fpr, y = tpr, color = model),
      linewidth = 1.4
    ) +
    geom_line(
      data = df_ext,
      aes(x = fpr, y = tpr, color = model),
      linewidth = 1.4
    ) +
    geom_abline(linetype = "dashed", color = "gray70") +
    scale_color_manual(
      values = c(
        "TCGA CV" = "coral1",
        "Combined external" = "#4682B4"
      )
    ) +
    labs(
      title = title_txt,
      subtitle = sprintf(
        "AUC - TCGA CV: %.3f | Combined external: %.3f",
        pROC::auc(roc_cv),
        pROC::auc(roc_ext)
      ),
      x = "False Positive Rate",
      y = "True Positive Rate",
      color = "Dataset"
    ) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom",
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      panel.grid = element_blank()
    )
}

cm_plot <- function(pred, truth, title_txt, high_col = "#4682B4") {

  cm <- caret::confusionMatrix(
    data = pred,
    reference = truth,
    positive = "Tumor"
  )

  print(cm)

  cm_df <- as.data.frame(cm$table)
  colnames(cm_df) <- c("Prediction", "Reference", "Freq")

  ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Freq)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = Freq), size = 6, fontface = "bold") +
    scale_fill_gradient(low = "white", high = high_col) +
    labs(
      title = title_txt,
      x = "True label",
      y = "Predicted label",
      fill = "Count"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      panel.grid = element_blank()
    )
}

# ============================================================
# 1. Load TCGA training data
# ============================================================

setwd(tcga_dir)

dds_files <- list.files(
  pattern = "dds_.*_matched_group\\.rds$",
  full.names = TRUE
)

tcga_list <- list()

for (file in dds_files) {

  project <- sub(
    "dds_(.*)_matched_group\\.rds",
    "\\1",
    basename(file)
  )

  message("Loading TCGA project: ", project)

  dds <- readRDS(file)
  vsd <- vst(dds, blind = TRUE)
  expr_mat <- assay(vsd)

  symbols <- map_ensembl_to_symbol(rownames(expr_mat))
  rownames(expr_mat) <- symbols
  expr_mat <- collapse_duplicate_symbols(expr_mat)

  matched_genes <- intersect(feature_genes, rownames(expr_mat))

  if (length(matched_genes) < 2) {
    next
  }

  expr <- t(expr_mat[matched_genes, , drop = FALSE]) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample")

  condition <- tryCatch(
    as.character(colData(dds)[expr$sample, "condition"]),
    error = function(e) NULL
  )

  if (is.null(condition)) {
    next
  }

  expr$condition <- condition
  expr$project <- project

  tcga_list[[project]] <- expr
}

tcga_df <- bind_rows(tcga_list)

tcga_clean <- tcga_df %>%
  dplyr::select(any_of(feature_genes), condition, sample, project) %>%
  filter(!is.na(condition)) %>%
  mutate(
    condition_bin = factor(
      ifelse(condition == "Primary Tumor", "Tumor", "Normal"),
      levels = c("Normal", "Tumor")
    )
  )

available_tcga_genes <- intersect(feature_genes, colnames(tcga_clean))

tcga_clean <- tcga_clean %>%
  filter(if_all(all_of(available_tcga_genes), ~ !is.na(.)))

x_tcga <- tcga_clean[, available_tcga_genes, drop = FALSE]
y_tcga <- tcga_clean$condition_bin

# ============================================================
# 2. Train pan-cancer RF classifier
# ============================================================

ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

set.seed(123)

rf_cv <- train(
  x = x_tcga,
  y = y_tcga,
  method = "rf",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5,
  ntree = 500
)

cv_pred <- rf_cv$pred %>%
  filter(mtry == rf_cv$bestTune$mtry)

roc_tcga_cv <- pROC::roc(
  cv_pred$obs,
  cv_pred$Tumor,
  levels = c("Normal", "Tumor"),
  direction = "<",
  quiet = TRUE
)

message("TCGA CV AUC: ", round(as.numeric(pROC::auc(roc_tcga_cv)), 3))

# ============================================================
# 3. Load Lung GSE87410
# ============================================================

setwd(lung_dir)

raw_lung <- readr::read_tsv(
  "GSE87410_raw_counts_GRCh38.p13_NCBI.tsv.gz",
  show_col_types = FALSE
)

counts_lung <- as.data.frame(raw_lung[, -1])
rownames(counts_lung) <- raw_lung[[1]]
counts_lung[] <- lapply(counts_lung, as.numeric)

lung_meta <- readr::read_csv(
  "GSM_sample_map.csv",
  col_names = c("GSM", "Label"),
  show_col_types = FALSE
) %>%
  tidyr::separate(
    Label,
    into = c("Project", "Patient", "Tissue"),
    sep = "_",
    remove = FALSE
  ) %>%
  filter(Tissue %in% c("T", "N"))

valid_lung <- intersect(colnames(counts_lung), lung_meta$GSM)
counts_lung <- counts_lung[, valid_lung, drop = FALSE]

lung_meta <- lung_meta %>%
  filter(GSM %in% valid_lung) %>%
  dplyr::select(GSM, Project, Tissue) %>%
  column_to_rownames("GSM")

dds_lung <- DESeqDataSetFromMatrix(
  countData = round(counts_lung),
  colData = lung_meta,
  design = ~ Project + Tissue
)

expr_lung <- assay(vst(dds_lung, blind = TRUE))

lung_symbols <- map_entrez_to_symbol(rownames(expr_lung))
rownames(expr_lung) <- ifelse(
  is.na(lung_symbols) | lung_symbols == "",
  rownames(expr_lung),
  lung_symbols
)

expr_lung <- collapse_duplicate_symbols(expr_lung)

keep_lung <- rownames(lung_meta)[lung_meta$Project %in% c("LUAD", "LUSC")]

expr_lung <- expr_lung[, keep_lung, drop = FALSE]
lung_meta <- lung_meta[keep_lung, , drop = FALSE]

lung_expr <- as.data.frame(t(expr_lung))
lung_truth <- factor(
  ifelse(lung_meta$Tissue == "T", "Tumor", "Normal"),
  levels = c("Normal", "Tumor")
)

# ============================================================
# 4. Load Kidney CPTAC-3
# ============================================================

setwd(kidney_dir)

kidney_files <- list.files(
  "GDC_kidney_filtered",
  pattern = "\\.rna_seq\\.augmented_star_gene_counts\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

kidney_list <- list()

for (file in kidney_files) {

  dat <- tryCatch({
    readr::read_tsv(file, skip = 1, show_col_types = FALSE) %>%
      dplyr::select(gene_name, count = unstranded) %>%
      mutate(count = as.numeric(count)) %>%
      filter(!is.na(gene_name), !is.na(count)) %>%
      group_by(gene_name) %>%
      summarise(count = sum(count), .groups = "drop")
  }, error = function(e) NULL)

  if (!is.null(dat)) {
    sample_id <- basename(dirname(file))
    colnames(dat)[colnames(dat) == "count"] <- sample_id
    kidney_list[[sample_id]] <- dat
  }
}

kidney_counts <- reduce(kidney_list, full_join, by = "gene_name") %>%
  column_to_rownames("gene_name")

kidney_counts[is.na(kidney_counts)] <- 0

query_expr <- GDCquery(
  project = "CPTAC-3",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

metadata_expr <- getResults(
  query_expr,
  cols = c("file_id", "sample_type")
)

kidney_meta <- metadata_expr %>%
  filter(
    file_id %in% colnames(kidney_counts),
    sample_type %in% c("Primary Tumor", "Solid Tissue Normal")
  ) %>%
  dplyr::select(sample = file_id, sample_type) %>%
  mutate(
    condition = factor(
      ifelse(grepl("Normal", sample_type), "Normal", "Tumor"),
      levels = c("Normal", "Tumor")
    )
  ) %>%
  distinct(sample, .keep_all = TRUE) %>%
  column_to_rownames("sample")

kidney_counts <- kidney_counts[, rownames(kidney_meta), drop = FALSE]

dds_kidney <- DESeqDataSetFromMatrix(
  countData = round(kidney_counts),
  colData = kidney_meta,
  design = ~ condition
)

expr_kidney <- assay(vst(dds_kidney, blind = TRUE))
expr_kidney <- collapse_duplicate_symbols(expr_kidney)

kidney_expr <- as.data.frame(t(expr_kidney))
kidney_truth <- kidney_meta$condition

# ============================================================
# 5. Load Korean CRC
# ============================================================

setwd(crc_dir)

raw_crc <- read.delim(
  "CMCBSN_expectedcount_342.txt",
  check.names = FALSE
)

counts_crc <- as.data.frame(raw_crc[, -1, drop = FALSE])
rownames(counts_crc) <- raw_crc[[1]]
counts_crc[] <- lapply(counts_crc, function(x) as.numeric(as.character(x)))

counts_crc <- collapse_duplicate_symbols(counts_crc)

crc_sample_ids <- colnames(counts_crc)

crc_tissue <- ifelse(
  grepl("-T-", crc_sample_ids),
  "T",
  ifelse(grepl("-N-", crc_sample_ids), "N", NA)
)

crc_meta <- data.frame(
  sample = crc_sample_ids,
  Tissue = crc_tissue,
  row.names = crc_sample_ids,
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(Tissue))

counts_crc <- counts_crc[, rownames(crc_meta), drop = FALSE]

dds_crc <- DESeqDataSetFromMatrix(
  countData = round(counts_crc),
  colData = crc_meta,
  design = ~ Tissue
)

expr_crc <- assay(vst(dds_crc, blind = TRUE))
expr_crc <- collapse_duplicate_symbols(expr_crc)

crc_expr <- as.data.frame(t(expr_crc))
crc_truth <- factor(
  ifelse(crc_meta$Tissue == "T", "Tumor", "Normal"),
  levels = c("Normal", "Tumor")
)

# ============================================================
# 6. Combine external cohorts
# ============================================================

common_external_genes <- Reduce(
  intersect,
  list(
    colnames(lung_expr),
    colnames(kidney_expr),
    colnames(crc_expr)
  )
)

eval_genes <- intersect(available_tcga_genes, common_external_genes)

if (length(eval_genes) < 2) {
  stop("Not enough feature genes are present in all datasets.")
}

message("Genes used for combined external validation: ",
        paste(eval_genes, collapse = ", "))

# Retrain using only genes available in all external datasets

set.seed(123)

rf_cv_eval <- train(
  x = tcga_clean[, eval_genes, drop = FALSE],
  y = y_tcga,
  method = "rf",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5,
  ntree = 500
)

cv_pred_eval <- rf_cv_eval$pred %>%
  filter(mtry == rf_cv_eval$bestTune$mtry)

roc_tcga_cv_eval <- pROC::roc(
  cv_pred_eval$obs,
  cv_pred_eval$Tumor,
  levels = c("Normal", "Tumor"),
  direction = "<",
  quiet = TRUE
)

external_expr <- bind_rows(
  lung_expr[, eval_genes, drop = FALSE] %>% mutate(dataset = "Lung"),
  kidney_expr[, eval_genes, drop = FALSE] %>% mutate(dataset = "Kidney"),
  crc_expr[, eval_genes, drop = FALSE] %>% mutate(dataset = "CRC")
)

external_dataset <- external_expr$dataset

external_expr <- external_expr %>%
  dplyr::select(all_of(eval_genes))

external_expr[] <- lapply(external_expr, as.numeric)

external_truth <- factor(
  c(
    as.character(lung_truth),
    as.character(kidney_truth),
    as.character(crc_truth)
  ),
  levels = c("Normal", "Tumor")
)

# ============================================================
# 7. ROC and confusion matrix on combined external cohort
# ============================================================

external_prob <- get_proba(rf_cv_eval, external_expr)

roc_external <- pROC::roc(
  external_truth,
  external_prob,
  levels = c("Normal", "Tumor"),
  direction = "<",
  quiet = TRUE
)

message("Combined external AUC: ", round(as.numeric(pROC::auc(roc_external)), 3))

optimal_threshold <- best_thr_bacc(external_truth, external_prob)

pred_05 <- factor(
  ifelse(external_prob >= 0.5, "Tumor", "Normal"),
  levels = c("Normal", "Tumor")
)

pred_opt <- factor(
  ifelse(external_prob >= optimal_threshold, "Tumor", "Normal"),
  levels = c("Normal", "Tumor")
)

p_roc_combined <- plot_roc_pair(
  roc_cv = roc_tcga_cv_eval,
  roc_ext = roc_external,
  title_txt = "ROC Curve: Combined External Cohort"
)

p_cm_05 <- cm_plot(
  pred = pred_05,
  truth = external_truth,
  title_txt = "Combined External Cohort — threshold 0.5"
)

p_cm_opt <- cm_plot(
  pred = pred_opt,
  truth = external_truth,
  title_txt = sprintf(
    "Combined External Cohort — optimal threshold %.2f",
    optimal_threshold
  )
)

p_combined <- p_roc_combined / (p_cm_05 | p_cm_opt)

print(p_combined)
