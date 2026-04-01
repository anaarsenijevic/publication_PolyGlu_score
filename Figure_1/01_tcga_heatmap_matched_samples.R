# TCGA matched-sample heatmap of TTLL/CCP family dysregulation
#

#
# Input:
# - dds_*_matched_group.rds files
#
# Notes:
# - Edit 'input_dir' below before running
# - CCP genes are displayed using CCP1-6 labels

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(pheatmap)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# ------------------------------------------------------------------
# Input path
# ------------------------------------------------------------------

input_dir <- "path-to-your-dds-files"

if (!dir.exists(input_dir)) {
  stop("Please set 'input_dir' to the folder containing dds_*_matched_group.rds files.")
}

# ------------------------------------------------------------------
# Genes of interest
# ------------------------------------------------------------------

feature_genes <- c(
  "TTLL1", "TTLL2", "TTLL4", "TTLL5", "TTLL6", "TTLL7", "TTLL9", "TTLL11", "TTLL13",
  "AGTPBP1", "AGBL2", "AGBL3", "AGBL1", "AGBL5", "AGBL4"
)

symbol_labels <- c(
  "TTLL1", "TTLL2", "TTLL4", "TTLL5", "TTLL6", "TTLL7", "TTLL9", "TTLL11", "TTLL13",
  "CCP1", "CCP2", "CCP3", "CCP4", "CCP5", "CCP6"
)
names(symbol_labels) <- feature_genes

# ------------------------------------------------------------------
# Read DESeq2 objects and collect fold changes
# ------------------------------------------------------------------

rds_files <- list.files(
  input_dir,
  pattern = "dds_.*_matched_group\\.rds$",
  full.names = TRUE
)

fc_matrix <- list()
fdr_matrix <- list()

for (file in rds_files) {
  project <- sub("dds_(.*)_matched_group\\.rds", "\\1", basename(file))
  message("Processing ", project)

  dds <- readRDS(file)

  ens_ids <- rownames(dds)
  ens_ids_clean <- sub("\\.\\d+$", "", ens_ids)

  symbols <- mapIds(
    org.Hs.eg.db,
    keys = ens_ids_clean,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )

  mcols(dds)$symbol <- symbols

  res <- results(dds, contrast = c("condition", "Primary Tumor", "Solid Tissue Normal"))
  res_df <- as.data.frame(res)
  res_df$symbol <- mcols(dds)$symbol

  res_filtered <- res_df %>%
    filter(symbol %in% feature_genes) %>%
    group_by(symbol) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    complete(symbol = feature_genes) %>%
    arrange(factor(symbol, levels = feature_genes))

  fc_vector <- res_filtered$log2FoldChange
  fdr_vector <- res_filtered$padj

  names(fc_vector) <- res_filtered$symbol
  names(fdr_vector) <- res_filtered$symbol

  fc_matrix[[project]] <- fc_vector[feature_genes]
  fdr_matrix[[project]] <- fdr_vector[feature_genes]
}

fc_df <- do.call(cbind, fc_matrix)
fdr_df <- do.call(cbind, fdr_matrix)

rownames(fc_df) <- symbol_labels[feature_genes]
rownames(fdr_df) <- symbol_labels[feature_genes]

# ------------------------------------------------------------------
# Heatmap settings
# ------------------------------------------------------------------

max_abs_fc <- max(abs(fc_df), na.rm = TRUE)
color_breaks <- seq(-max_abs_fc, max_abs_fc, length.out = 101)
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# ------------------------------------------------------------------
# Full heatmap
# ------------------------------------------------------------------

pheatmap(
  fc_df,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = color_palette,
  breaks = color_breaks,
  na_col = "gray90",
  main = "Log2 fold change in primary tumours vs normal tissues",
  cellwidth = 30,
  cellheight = 15,
  fontsize_row = 10,
  fontsize_col = 10
)

# ------------------------------------------------------------------
# Significant values only
# ------------------------------------------------------------------

sig_mask <- fdr_df < 0.05
fc_masked <- fc_df
fc_masked[!sig_mask] <- NA

pheatmap(
  fc_masked,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = color_palette,
  breaks = color_breaks,
  na_col = "gray90",
  main = "Log2 fold change in primary tumours vs normal tissues (FDR < 0.05)",
  cellwidth = 30,
  cellheight = 15,
  fontsize_row = 10,
  fontsize_col = 10
)
