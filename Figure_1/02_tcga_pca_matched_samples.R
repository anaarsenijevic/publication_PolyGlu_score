# TCGA matched-sample PCA of TTLL/CCP expression
#
#
# Input:
# - dds_*_matched_group.rds files
#
# Notes:
# - Edit 'input_dir' below before running
# - PCA is run on VST-transformed expression values
# - The plot shows all samples together

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(biomaRt)
  library(ggplot2)
  library(purrr)
  library(tibble)
})

# ------------------------------------------------------------------
# Input path
# ------------------------------------------------------------------

input_dir <- "path-to-your-dds-files"

if (!dir.exists(input_dir)) {
  stop("Please set 'input_dir' to the folder containing dds_*_matched_group.rds files.")
}

# ------------------------------------------------------------------
# Gene sets
# ------------------------------------------------------------------

ttll_genes <- c("TTLL1", "TTLL2", "TTLL4", "TTLL5", "TTLL6", "TTLL7", "TTLL9", "TTLL11", "TTLL13")
ccp_genes  <- c("AGBL1", "AGBL2", "AGBL3", "AGBL4", "AGBL5", "AGTPBP1")
all_genes  <- unique(c(ttll_genes, ccp_genes))

# ------------------------------------------------------------------
# Ensembl to gene symbol mapping
# ------------------------------------------------------------------

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_map <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  mart = ensembl
) %>%
  filter(external_gene_name != "") %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

# ------------------------------------------------------------------
# Read and process matched TCGA data
# ------------------------------------------------------------------

dds_files <- list.files(
  input_dir,
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
  expr_mat <- assay(vsd)
  rownames(expr_mat) <- symbols

  matched_genes <- intersect(all_genes, rownames(expr_mat))
  if (!length(matched_genes)) next

  expr <- t(expr_mat[matched_genes, , drop = FALSE]) %>% as.data.frame()
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
  select(all_of(all_genes), condition, sample, project) %>%
  filter(if_all(all_of(all_genes), ~ !is.na(.))) %>%
  mutate(
    condition_bin = factor(
      ifelse(condition == "Primary Tumor", "Tumor", "Normal"),
      levels = c("Normal", "Tumor")
    )
  )

col_map <- c("Tumor" = "#FF6A6A", "Normal" = "lightskyblue")

# ------------------------------------------------------------------
# PCA across all matched samples
# ------------------------------------------------------------------

expr_all <- tcga_clean %>%
  select(where(is.numeric)) %>%
  as.matrix()

rownames(expr_all) <- tcga_clean$sample

pca_all <- prcomp(expr_all, center = TRUE, scale. = TRUE)
ve_all <- summary(pca_all)$importance[2, 1:2] * 100

pca_df_all <- data.frame(
  PC1 = pca_all$x[, 1],
  PC2 = pca_all$x[, 2],
  condition_bin = tcga_clean$condition_bin
)

p_pca_all <- ggplot(pca_df_all, aes(PC1, PC2, color = condition_bin)) +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_manual(values = col_map) +
  labs(
    title = "PCA of TCGA matched samples using TTLL and CCP genes",
    x = paste0("PC1 (", round(ve_all[1], 1), "%)"),
    y = paste0("PC2 (", round(ve_all[2], 1), "%)")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p_pca_all)

