#!/usr/bin/env Rscript

# RF tumor probability across fibroblast transformation stages
#
# This script:
# 1. trains a pan-cancer random forest classifier on matched TCGA tumour-normal cohorts
# 2. uses the fixed 9-gene TTLL/AGBL panel
# 3. applies the trained model to Danielsson et al. fibroblast transformation data
# 4. predicts RF tumor probability across Primary, hTERT, SV40, and Ras stages
# 5. plots RF tumor probability by transformation stage
# 6. performs a Spearman trend test across ordered stages
#
# Input:
# - dds_*_matched_group.rds files
# - annotated_transcriptomes_Danielsson2013.csv (download from doi: 10.1073/pnas.1216436110)
#
# Notes:
# - Edit all input paths below before running
# - TCGA training data are VST-transformed
# - Danielsson expression values are FPKM
# - Danielsson values are scaled to the TCGA training distribution per gene
# - Missing panel genes in Danielsson are imputed at the TCGA training mean
# - The fixed 9-gene panel is:
#   TTLL4, TTLL5, TTLL6, TTLL7, TTLL11, AGBL1, AGBL3, AGBL4, AGBL5

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(biomaRt)
  library(randomForest)
  library(caret)
  library(pROC)
  library(ggplot2)
})

set.seed(123)

# ------------------------------------------------------------------
# Input paths
# ------------------------------------------------------------------

tcga_dir <- "path-to-tcga-dds-files"
danielsson_file <- "path-to-danielsson-folder/annotated_transcriptomes_Danielsson2013.csv"
out_dir <- "path-to-output-folder"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(tcga_dir)) {
  stop("Please set 'tcga_dir' to the folder containing dds_*_matched_group.rds files.")
}

if (!file.exists(danielsson_file)) {
  stop("Could not find: ", danielsson_file)
}

# ------------------------------------------------------------------
# Fixed 9-gene panel
# ------------------------------------------------------------------

feature_genes <- c(
  "TTLL4", "TTLL5", "TTLL6", "TTLL7", "TTLL11",
  "AGBL1", "AGBL3", "AGBL4", "AGBL5"
)

# ------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------

aggregate_dup_symbols <- function(mat, symbols) {
  keep <- !is.na(symbols) & symbols != ""
  mat <- mat[keep, , drop = FALSE]
  symbols <- symbols[keep]

  split_idx <- split(seq_len(nrow(mat)), symbols)

  out <- vapply(
    split_idx,
    function(ix) colMeans(mat[ix, , drop = FALSE]),
    numeric(ncol(mat))
  )

  out <- t(out)
  rownames(out) <- names(split_idx)
  out
}

# ------------------------------------------------------------------
# Ensembl to gene symbol mapping
# ------------------------------------------------------------------

message("Fetching Ensembl to gene symbol mapping via biomaRt...")

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_map <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  mart = ensembl
) %>%
  filter(external_gene_name != "") %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

# ------------------------------------------------------------------
# Load TCGA matched cohorts and train pan-cancer RF
# ------------------------------------------------------------------

message("Loading matched TCGA DESeq2 objects...")

dds_files <- list.files(
  tcga_dir,
  pattern = "dds_.*_matched_group\\.rds$",
  full.names = TRUE
)

if (length(dds_files) == 0) {
  stop("No matched TCGA DESeq2 files found in 'tcga_dir'.")
}

tcga_list <- list()

for (file in dds_files) {
  project <- sub("dds_(.*)_matched_group\\.rds$", "\\1", basename(file))
  message("  training: ", project)

  dds <- readRDS(file)
  vsd <- vst(dds, blind = TRUE)
  mat <- assay(vsd)

  ens <- sub("\\..*$", "", rownames(mat))
  sym <- gene_map$external_gene_name[match(ens, gene_map$ensembl_gene_id)]
  mat_sym <- aggregate_dup_symbols(mat, sym)

  matched <- intersect(feature_genes, rownames(mat_sym))
  if (!length(matched)) next

  expr <- t(mat_sym[matched, , drop = FALSE]) %>%
    as.data.frame()

  expr$sample <- rownames(expr)
  expr$project <- project

  cond <- tryCatch(
    colData(dds)[expr$sample, "condition"],
    error = function(e) NULL
  )

  if (is.null(cond)) next

  expr$condition <- as.character(cond)

  tcga_list[[length(tcga_list) + 1]] <- expr
}

tcga_df <- bind_rows(tcga_list)

if (nrow(tcga_df) == 0) {
  stop("No training data could be assembled from the matched TCGA cohorts.")
}

tcga_df <- tcga_df %>%
  mutate(
    condition_bin = factor(
      ifelse(condition == "Primary Tumor", "Tumor", "Normal"),
      levels = c("Normal", "Tumor")
    )
  ) %>%
  dplyr::select(sample, project, condition, all_of(feature_genes), condition_bin) %>%
  filter(if_all(all_of(feature_genes), ~ !is.na(.)))

x_tcga <- tcga_df %>%
  dplyr::select(all_of(feature_genes))

y_tcga <- tcga_df$condition_bin

ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

rf_cv <- train(
  x = x_tcga,
  y = y_tcga,
  method = "rf",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5
)

best_mtry <- rf_cv$bestTune$mtry

final_rf <- randomForest(
  x = x_tcga,
  y = y_tcga,
  mtry = best_mtry,
  ntree = 500,
  importance = TRUE
)

tcga_preds <- rf_cv$pred %>%
  filter(mtry == best_mtry)

roc_tcga <- pROC::roc(
  tcga_preds$obs,
  tcga_preds$Tumor,
  quiet = TRUE
)

cat(sprintf(
  "\nTCGA cross-validation AUC = %.3f\n",
  as.numeric(pROC::auc(roc_tcga))
))

saveRDS(
  final_rf,
  file.path(out_dir, "rf_pan_tcga_9gene_RAWVST.rds")
)

# ------------------------------------------------------------------
# Load Danielsson fibroblast transformation data
# ------------------------------------------------------------------

message("Loading Danielsson fibroblast transformation data...")

dan <- readr::read_csv(
  danielsson_file,
  show_col_types = FALSE
)

expr_mat <- dan %>%
  dplyr::select(
    hgnc_symbol,
    tidyselect::starts_with("FPKM.")
  ) %>%
  filter(!is.na(hgnc_symbol), hgnc_symbol != "") %>%
  distinct(hgnc_symbol, .keep_all = TRUE) %>%
  tibble::column_to_rownames("hgnc_symbol")

genes_train <- colnames(x_tcga)
genes_dan <- intersect(genes_train, rownames(expr_mat))

if (length(genes_dan) == 0) {
  stop("None of the training genes are present in the Danielsson dataset.")
}

if (length(genes_dan) < length(genes_train)) {
  warning(
    "Danielsson dataset contains ",
    length(genes_dan),
    "/",
    length(genes_train),
    " training genes: ",
    paste(genes_dan, collapse = ", ")
  )
}

# ------------------------------------------------------------------
# Prepare Danielsson prediction matrix
# ------------------------------------------------------------------

x_dan <- t(expr_mat[genes_dan, , drop = FALSE]) %>%
  as.data.frame()

x_dan$sample <- rownames(x_dan)

x_dan <- x_dan %>%
  mutate(
    stage = case_when(
      grepl("primary", sample, ignore.case = TRUE) ~ "Primary",
      grepl("hTERT", sample, ignore.case = TRUE) ~ "hTERT",
      grepl("SV40", sample, ignore.case = TRUE) ~ "SV40",
      grepl("Ras", sample, ignore.case = TRUE) ~ "Ras",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(stage))

train_means <- sapply(
  genes_train,
  function(g) mean(x_tcga[[g]], na.rm = TRUE)
)

train_sds <- sapply(
  genes_train,
  function(g) sd(x_tcga[[g]], na.rm = TRUE)
)

train_sds[train_sds == 0 | is.na(train_sds)] <- 1

for (g in genes_dan) {
  x_dan[[g]] <- (x_dan[[g]] - train_means[g]) / train_sds[g]
}

missing_in_new <- setdiff(genes_train, colnames(x_dan))

if (length(missing_in_new) > 0) {
  for (g in missing_in_new) {
    x_dan[[g]] <- 0
  }
}

extra_in_new <- setdiff(
  colnames(x_dan),
  c("sample", "stage", genes_train)
)

message(
  "Missing genes imputed at training mean: ",
  ifelse(length(missing_in_new) == 0, "none", paste(missing_in_new, collapse = ", "))
)

message(
  "Extra columns ignored: ",
  ifelse(length(extra_in_new) == 0, "none", paste(extra_in_new, collapse = ", "))
)

x_pred <- x_dan[, c("sample", "stage", genes_train), drop = FALSE]

x_pred_mat <- x_pred[, genes_train, drop = FALSE]
x_pred_mat[] <- lapply(
  x_pred_mat,
  function(x) as.numeric(as.character(x))
)

# ------------------------------------------------------------------
# Predict RF tumor probability in Danielsson samples
# ------------------------------------------------------------------

x_pred$rf_prob_tumor <- predict(
  final_rf,
  newdata = x_pred_mat,
  type = "prob"
)[, "Tumor"]

stage_order <- c(
  "Primary" = 1,
  "hTERT" = 2,
  "SV40" = 3,
  "Ras" = 4
)

x_pred$stage <- factor(
  x_pred$stage,
  levels = names(stage_order)
)

x_pred$stage_idx <- stage_order[as.character(x_pred$stage)]

write.csv(
  x_pred,
  file.path(out_dir, "danielsson_rf_tumor_probability.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Spearman trend test
# ------------------------------------------------------------------

trend_test <- cor.test(
  x_pred$stage_idx,
  x_pred$rf_prob_tumor,
  method = "spearman",
  exact = FALSE
)

print(trend_test)

trend_result <- tibble(
  method = "Spearman correlation",
  rho = unname(trend_test$estimate),
  p = trend_test$p.value
)

write.csv(
  trend_result,
  file.path(out_dir, "danielsson_rf_tumor_probability_spearman.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Plot RF tumor probability across transformation stages
# ------------------------------------------------------------------

p_danielsson <- ggplot(
  x_pred,
  aes(x = stage, y = rf_prob_tumor)
) +
  geom_jitter(
    width = 0.15,
    size = 3,
    alpha = 0.9
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 5,
    color = "black"
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    color = "black"
  ) +
  labs(
    title = "RF tumor probability across fibroblast transformation stages",
    subtitle = paste0(
      "Missing genes imputed at training mean: ",
      ifelse(
        length(missing_in_new) == 0,
        "none",
        paste(missing_in_new, collapse = ", ")
      )
    ),
    x = "Transformation stage",
    y = "RF predicted tumor probability"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    plot.margin = ggplot2::margin(15, 20, 15, 15)
  )

print(p_danielsson)

ggsave(
  file.path(out_dir, "danielsson_rf_tumor_probability_by_stage.pdf"),
  p_danielsson,
  width = 6.5,
  height = 4.5
)

ggsave(
  file.path(out_dir, "danielsson_rf_tumor_probability_by_stage.png"),
  p_danielsson,
  width = 6.5,
  height = 4.5,
  dpi = 300
)

message("\nFinished. Outputs written to:\n", out_dir)
