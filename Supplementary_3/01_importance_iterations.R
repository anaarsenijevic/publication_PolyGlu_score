# TCGA random forest feature elimination and breakpoint analysis
#
# This script:
# 1. Loads matched TCGA DESeq2 datasets
# 2. Extracts VST expression for TTLL/CCP genes
# 3. Performs backward feature elimination using random forest importance
# 4. Tracks cross-validated AUC at each step
# 5. Plots:
#    - Feature importance across iterations (faceted)
#    - TCGA AUC vs gene number with segmented regression breakpoint
#
# Input:
# - dds_*_matched_group.rds files
#
# Notes:
# - Set 'input_dir' before running
# - Requires 'condition' column in colData
# - Tumour defined as "Primary Tumor"

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(randomForest)
  library(caret)
  library(segmented)
  library(tibble)
})

# ------------------------------------------------------------------
# Input directory
# ------------------------------------------------------------------

input_dir <- "path-to-your-dds-files"

if (!dir.exists(input_dir)) {
  stop("Set 'input_dir' to folder with dds_*_matched_group.rds files.")
}

# ------------------------------------------------------------------
# Gene panel
# ------------------------------------------------------------------

ttll_genes <- c("TTLL1","TTLL2","TTLL4","TTLL5","TTLL6","TTLL7","TTLL9","TTLL11","TTLL13")
ccp_genes  <- c("AGBL1","AGBL2","AGBL3","AGBL4","AGBL5","AGTPBP1")
feature_genes <- unique(c(ttll_genes, ccp_genes))

gene_labels <- c(
  TTLL1="TTLL1", TTLL2="TTLL2", TTLL4="TTLL4", TTLL5="TTLL5",
  TTLL6="TTLL6", TTLL7="TTLL7", TTLL9="TTLL9", TTLL11="TTLL11",
  TTLL13="TTLL13",
  AGTPBP1="CCP1", AGBL2="CCP2", AGBL3="CCP3",
  AGBL1="CCP4", AGBL5="CCP5", AGBL4="CCP6"
)

# ------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------

aggregate_dup_symbols <- function(mat, symbols) {
  ok <- !is.na(symbols) & symbols != ""
  mat <- mat[ok, , drop = FALSE]
  symbols <- symbols[ok]

  out <- t(vapply(
    split(seq_len(nrow(mat)), symbols),
    function(ix) colMeans(mat[ix, , drop = FALSE]),
    numeric(ncol(mat))
  ))

  as.matrix(out)
}

make_symbol_map <- function(ensembl_ids) {
  tibble(
    ensembl_gene_id = unique(ensembl_ids),
    external_gene_name = AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = unique(ensembl_ids),
      keytype = "ENSEMBL",
      column = "SYMBOL",
      multiVals = "first"
    )
  ) %>%
    filter(!is.na(external_gene_name), external_gene_name != "")
}

vst_to_symbol_df <- function(vst_matrix, feature_genes, gene_map) {
  ensembl_ids <- sub("\\..*$", "", rownames(vst_matrix))
  symbols <- gene_map$external_gene_name[
    match(ensembl_ids, gene_map$ensembl_gene_id)
  ]

  vst_by_symbol <- aggregate_dup_symbols(vst_matrix, symbols)
  matched_genes <- intersect(feature_genes, rownames(vst_by_symbol))

  t(vst_by_symbol[matched_genes, , drop = FALSE]) %>%
    as.data.frame()
}

# ------------------------------------------------------------------
# Load TCGA data
# ------------------------------------------------------------------

dds_files <- list.files(
  input_dir,
  pattern = "dds_.*_matched_group\\.rds$",
  full.names = TRUE
)

tcga_data <- list()

for (file in dds_files) {
  project <- sub("dds_(.*)_matched_group\\.rds", "\\1", basename(file))
  dds <- readRDS(file)

  vsd <- vst(dds, blind = TRUE)
  vst_mat <- assay(vsd)

  gene_map <- make_symbol_map(sub("\\..*$", "", rownames(vst_mat)))

  expr <- tryCatch(
    vst_to_symbol_df(vst_mat, feature_genes, gene_map),
    error = function(e) NULL
  )

  if (is.null(expr)) next

  expr$sample <- rownames(expr)
  expr$condition <- as.character(colData(dds)[expr$sample, "condition"])
  expr$project <- project

  tcga_data[[length(tcga_data) + 1]] <- expr
}

tcga_df <- bind_rows(tcga_data) %>%
  select(any_of(feature_genes), condition, sample, project) %>%
  mutate(
    condition_bin = factor(
      ifelse(condition == "Primary Tumor", "Tumor", "Normal"),
      levels = c("Normal", "Tumor")
    )
  ) %>%
  filter(if_all(all_of(intersect(feature_genes, colnames(.))), ~ !is.na(.)))

# ------------------------------------------------------------------
# Feature elimination
# ------------------------------------------------------------------

set.seed(123)

current_genes <- feature_genes[feature_genes %in% colnames(tcga_df)]
y <- tcga_df$condition_bin
iterations <- list()

ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

for (i in seq(length(current_genes), 5, by = -1)) {

  x <- tcga_df[, current_genes, drop = FALSE]

  rf_cv <- train(
    x = x,
    y = y,
    method = "rf",
    trControl = ctrl,
    metric = "ROC"
  )

  rf <- randomForest(
    x = x,
    y = y,
    mtry = rf_cv$bestTune$mtry,
    ntree = 500,
    importance = TRUE
  )

  imp <- importance(rf, type = 1) %>%
    as.data.frame() %>%
    rownames_to_column("Gene") %>%
    arrange(MeanDecreaseAccuracy) %>%
    mutate(
      Genes_in_model = length(current_genes),
      Gene_label = gene_labels[Gene],
      Removed_next = row_number() == 1,
      Gene_ranked = paste(Gene_label, Genes_in_model, sep = "_")
    )

  iterations[[length(current_genes)]] <- list(
    genes = current_genes,
    auc = max(rf_cv$results$ROC),
    importance = imp
  )

  current_genes <- setdiff(current_genes, imp$Gene[1])
}

# ------------------------------------------------------------------
# Collect results
# ------------------------------------------------------------------

auc_df <- map_dfr(iterations, ~ tibble(
  Genes = length(.x$genes),
  TCGA_AUC = .x$auc
))

importance_all <- map_dfr(iterations, "importance")

# ------------------------------------------------------------------
# Plot 1: Importance (faceted)
# ------------------------------------------------------------------

p_importance <- importance_all %>%
  mutate(
    Panel = paste0(Genes_in_model, " genes"),
    Gene_ranked = factor(Gene_ranked, levels = Gene_ranked)
  ) %>%
  ggplot(aes(x = Gene_ranked, y = MeanDecreaseAccuracy, fill = Removed_next)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ Panel, scales = "free_y") +
  scale_x_discrete(labels = function(x) sub("_.*", "", x)) +
  scale_fill_manual(values = c("TRUE" = "grey70", "FALSE" = "#009ACD")) +
  theme_classic() +
  labs(y = "Mean decrease in accuracy", x = NULL)

print(p_importance)

# ------------------------------------------------------------------
# Plot 2: AUC + segmented regression
# ------------------------------------------------------------------

tcga_auc_df <- auc_df %>% arrange(desc(Genes))

lm_fit <- lm(TCGA_AUC ~ Genes, data = tcga_auc_df)

seg_fit <- segmented(lm_fit, seg.Z = ~Genes, psi = list(Genes = 9))

breakpoint <- seg_fit$psi[1, "Est."]
tcga_auc_df$fit <- fitted(seg_fit)

p_auc <- ggplot(tcga_auc_df, aes(Genes, TCGA_AUC)) +
  geom_line(color = "#009ACD") +
  geom_point(color = "#009ACD", size = 3) +
  geom_line(aes(y = fit), color = "black") +
  geom_vline(xintercept = breakpoint, linetype = "dashed", color = "red") +
  scale_x_reverse() +
  theme_classic() +
  labs(x = "Number of genes", y = "TCGA-CV AUC")

print(p_auc)

# ------------------------------------------------------------------
# Statistics
# ------------------------------------------------------------------

print(summary(seg_fit))
print(davies.test(lm_fit, seg.Z = ~Genes))
