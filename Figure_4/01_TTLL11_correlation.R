#!/usr/bin/env Rscript

# TCGA pan-cancer correlation between TTLL11 expression and RF tumor score
#
# This script:
# 1. reads matched TCGA tumour-normal DESeq2 objects
# 2. extracts VST-normalized expression for a fixed 9-gene TTLL/AGBL panel
# 3. trains a pan-cancer random forest classifier
# 4. predicts RF tumor probability for TCGA samples
# 5. plots the association between TTLL11 expression and RF tumor score:
#    - across all tumour samples
#    - separately by TCGA cancer type
#
# Input:
# - dds_*_matched_group.rds files
#
# Notes:
# - Edit all input paths below before running
# - The fixed 9-gene panel is:
#   TTLL4, TTLL5, TTLL6, TTLL7, TTLL11, AGBL1, AGBL3, AGBL4, AGBL5
# - Expression values are VST-normalized
# - No z-scoring is applied

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(biomaRt)
  library(caret)
  library(randomForest)
  library(pROC)
  library(ggplot2)
  library(ggpubr)
})

set.seed(123)

# ------------------------------------------------------------------
# Input paths
# ------------------------------------------------------------------

tcga_dir <- "path-to-tcga-dds-files"
out_dir <- "path-to-output-folder"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(tcga_dir)) {
  stop("Please set 'tcga_dir' to the folder containing dds_*_matched_group.rds files.")
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

aggregate_dup_symbols <- function(mat, symbols, mode = "mean") {
  keep <- !is.na(symbols) & symbols != ""
  mat <- mat[keep, , drop = FALSE]
  symbols <- symbols[keep]

  if (mode == "sum") {
    rowsum(mat, group = symbols)
  } else {
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
# Load matched TCGA cohorts
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
  message("  processing: ", project)

  dds <- readRDS(file)
  vsd <- vst(dds, blind = TRUE)
  mat <- assay(vsd)

  ens <- sub("\\..*$", "", rownames(mat))
  sym <- gene_map$external_gene_name[match(ens, gene_map$ensembl_gene_id)]
  mat_sym <- aggregate_dup_symbols(mat, sym, mode = "mean")

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
  stop("No TCGA expression data could be assembled.")
}

tcga_clean <- tcga_df %>%
  dplyr::select(all_of(feature_genes), sample, project, condition) %>%
  filter(if_all(all_of(feature_genes), ~ !is.na(.))) %>%
  mutate(
    condition_bin = factor(
      ifelse(condition == "Primary Tumor", "Tumor", "Normal"),
      levels = c("Normal", "Tumor")
    )
  )

# ------------------------------------------------------------------
# Train pan-cancer random forest classifier
# ------------------------------------------------------------------

message("Training pan-cancer random forest model...")

x_pan <- tcga_clean %>%
  dplyr::select(all_of(feature_genes))

y_pan <- tcga_clean$condition_bin

ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = TRUE
)

rf_pan <- caret::train(
  x = x_pan,
  y = y_pan,
  method = "rf",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5
)

cv_pan <- rf_pan$pred %>%
  filter(mtry == rf_pan$bestTune$mtry)

roc_pan <- pROC::roc(
  cv_pan$obs,
  cv_pan$Tumor,
  levels = c("Normal", "Tumor"),
  quiet = TRUE
)

cat(sprintf(
  "\nPan-cancer RF CV AUC = %.3f\n",
  as.numeric(pROC::auc(roc_pan))
))

saveRDS(
  rf_pan,
  file.path(out_dir, "rf_pan_tcga_9gene_RAWVST.rds")
)

# ------------------------------------------------------------------
# Predict RF tumor probability for all TCGA samples
# ------------------------------------------------------------------

tcga_scores <- tcga_clean %>%
  mutate(
    rf_pan_prob = as.numeric(
      predict(
        rf_pan,
        newdata = dplyr::select(., all_of(feature_genes)),
        type = "prob"
      )[,"Tumor"]
    ),
    TTLL11_expr = TTLL11
  )

write.csv(
  tcga_scores,
  file.path(out_dir, "tcga_rf_scores_with_TTLL11.csv"),
  row.names = FALSE
)

tcga_tumor <- tcga_scores %>%
  filter(condition_bin == "Tumor")

# ------------------------------------------------------------------
# Overall tumour-only correlation plot
# ------------------------------------------------------------------

p_tumor <- ggplot(
  tcga_tumor,
  aes(x = TTLL11_expr, y = rf_pan_prob)
) +
  geom_point(alpha = 0.5, size = 1.8, color = "firebrick") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1.1) +
  stat_cor(
    method = "spearman",
    aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
    label.x.npc = "left",
    label.y.npc = "top",
    color = "black",
    size = 5
  ) +
  labs(
    title = "TCGA tumours: TTLL11 expression vs RF tumor score",
    x = "TTLL11 expression (VST)",
    y = "RF tumor probability"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold"),
    plot.margin = ggplot2::margin(15, 20, 15, 15)
  )

print(p_tumor)

ggsave(
  file.path(out_dir, "TCGA_tumours_TTLL11_vs_RF_score_overall.pdf"),
  p_tumor,
  width = 6.5,
  height = 5
)

ggsave(
  file.path(out_dir, "TCGA_tumours_TTLL11_vs_RF_score_overall.png"),
  p_tumor,
  width = 6.5,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------------
# Per-cancer tumour-only correlation plot
# ------------------------------------------------------------------

p_tumor_project <- ggplot(
  tcga_tumor,
  aes(x = TTLL11_expr, y = rf_pan_prob)
) +
  geom_point(alpha = 0.5, size = 1.4, color = "firebrick") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  stat_cor(
    method = "spearman",
    aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
    label.x.npc = "left",
    label.y.npc = "top",
    color = "black",
    size = 3
  ) +
  facet_wrap(~ project, scales = "free") +
  labs(
    title = "TCGA tumours by cancer type: TTLL11 expression vs RF tumor score",
    x = "TTLL11 expression (VST)",
    y = "RF tumor probability"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold", size = 9),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold"),
    plot.margin = ggplot2::margin(15, 20, 15, 15)
  )

print(p_tumor_project)

ggsave(
  file.path(out_dir, "TCGA_tumours_TTLL11_vs_RF_score_by_cancer.pdf"),
  p_tumor_project,
  width = 13,
  height = 10
)

ggsave(
  file.path(out_dir, "TCGA_tumours_TTLL11_vs_RF_score_by_cancer.png"),
  p_tumor_project,
  width = 13,
  height = 10,
  dpi = 300
)

message("\nFinished. Outputs written to:\n", out_dir)
