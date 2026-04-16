# TCGA matched-sample PCA of 9-gene model expression
#
# This script reproduces the previously deposited TCGA matched-sample PCA
# analysis of TTLL/CCP expression, using the reduced 9-gene panel derived
# from random forest feature selection instead of the original full panel.
#
# Input:
# - dds_*_matched_group.rds files
#
# Notes:
# - Edit 'input_dir' below before running
# - PCA is run on VST-transformed expression values
# - One PCA is generated per TCGA project

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
# 9-gene model
# ------------------------------------------------------------------

feature_genes <- c(
  "TTLL4", "TTLL5", "TTLL6", "TTLL7", "TTLL11",
  "AGBL1", "AGBL3", "AGBL4", "AGBL5"
)

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

  matched_genes <- intersect(feature_genes, rownames(expr_mat))
  if (!length(matched_genes)) next

  expr <- t(expr_mat[matched_genes, , drop = FALSE]) %>% as.data.frame()

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
  select(all_of(feature_genes), condition, sample, project) %>%
  filter(if_all(all_of(feature_genes), ~ !is.na(.))) %>%
  mutate(
    condition_bin = factor(
      ifelse(condition == "Primary Tumor", "Tumor", "Normal"),
      levels = c("Normal", "Tumor")
    )
  )

col_map <- c("Tumor" = "#FF6A6A", "Normal" = "lightskyblue")

# ------------------------------------------------------------------
# PCA within each project
# ------------------------------------------------------------------

split_idx <- split(seq_len(nrow(tcga_clean)), tcga_clean$project)

per_proj_df <- purrr::map_dfr(names(split_idx), function(proj) {
  idx <- split_idx[[proj]]

  Xp <- tcga_clean[idx, ] %>%
    select(all_of(feature_genes)) %>%
    as.matrix()

  if (nrow(Xp) < 3 || ncol(Xp) < 2) {
    return(NULL)
  }

  rownames(Xp) <- tcga_clean$sample[idx]

  pca_p <- prcomp(Xp, center = TRUE, scale. = TRUE)
  ve_p <- summary(pca_p)$importance[2, 1:2] * 100

  tibble(
    PC1 = pca_p$x[, 1],
    PC2 = pca_p$x[, 2],
    condition_bin = tcga_clean$condition_bin[idx],
    project = proj,
    project_lab = sprintf("%s (PC1 %.1f%%; PC2 %.1f%%)", proj, ve_p[1], ve_p[2])
  )
})

per_proj_df$project_lab <- factor(
  per_proj_df$project_lab,
  levels = per_proj_df %>%
    distinct(project, project_lab) %>%
    pull(project_lab)
)

p_pca_per_project <- ggplot(
  per_proj_df,
  aes(x = PC1, y = PC2, color = condition_bin)
) +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_manual(values = col_map) +
  labs(
    title = "PCA of matched TCGA samples by project (9-gene model)",
    x = "PC1",
    y = "PC2"
  ) +
  facet_wrap(~ project_lab, scales = "free") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

print(p_pca_per_project)
