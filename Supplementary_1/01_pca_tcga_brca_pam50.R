# TCGA-BRCA PCA by PAM50 subtype
#
# Input:
# - TCGA-BRCA_RNAseq.rds
#
#This script:
#Loads TCGA-BRCA RNA-seq data (DESeq2 object)
#Performs variance stabilizing transformation (VST)
#Maps Ensembl gene IDs to gene symbols
#Extracts expression of TTLL and CCP genes involved in microtubule polyglutamylation
#Filters for primary tumor samples and aggregates data at the patient level
#Annotates samples with PAM50 molecular subtypes
#Performs principal component analysis (PCA) on the selected gene panel
#Visualizes the PCA colored by PAM50 subtype
#Saves PCA plots and associated data tables
#
#
# Notes:
# - Edit 'input_file' and 'output_dir' before running
# - PCA is run on VST-transformed expression values
# - Only primary tumor samples are included
# - PCA is based on 15 TTLL/CCP polyglutamylation-related genes

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(ggplot2)
})

set.seed(123)

# ------------------------------------------------------------------
# Input and output paths
# ------------------------------------------------------------------

input_file <- "path-to-your-file/TCGA-BRCA_RNAseq.rds"


if (!file.exists(input_file)) {
  stop("Please set 'input_file' to the location of TCGA-BRCA_RNAseq.rds.")
}



# ------------------------------------------------------------------
# Gene set
# ------------------------------------------------------------------

ttll_genes <- c(
  "TTLL1", "TTLL2", "TTLL4", "TTLL5", "TTLL6", "TTLL7",
  "TTLL9", "TTLL11", "TTLL13"
)

ccp_genes <- c(
  "AGBL1", "AGBL2", "AGBL3", "AGBL4", "AGBL5", "AGTPBP1"
)

pca_genes <- unique(c(ttll_genes, ccp_genes))

# ------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------

aggregate_duplicate_symbols <- function(expr_mat, symbols) {
  keep <- !is.na(symbols) & symbols != ""
  expr_mat <- expr_mat[keep, , drop = FALSE]
  symbols <- symbols[keep]

  symbol_split <- split(seq_len(nrow(expr_mat)), symbols)

  aggregated <- vapply(
    symbol_split,
    function(idx) colMeans(expr_mat[idx, , drop = FALSE], na.rm = TRUE),
    numeric(ncol(expr_mat))
  )

  aggregated <- t(aggregated)
  rownames(aggregated) <- names(symbol_split)

  aggregated
}

clean_na_strings <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)

  x[x %in% c(
    "NA", "N/A", "", "Unknown", "Not Reported", "NOT REPORTED",
    "[Not Available]", "[NOT AVAILABLE]", "not reported"
  )] <- NA_character_

  x
}

recode_pam50 <- function(x) {
  y <- clean_na_strings(x)
  y <- stringr::str_trim(y)
  y <- gsub("\\s+", "", y)

  y <- dplyr::case_when(
    stringr::str_to_lower(y) == "luma" ~ "LumA",
    stringr::str_to_lower(y) == "lumb" ~ "LumB",
    stringr::str_to_lower(y) == "basal" ~ "Basal",
    stringr::str_to_lower(y) %in% c("her2", "her2enriched", "her2+") ~ "Her2",
    stringr::str_to_lower(y) %in% c("normal", "normal-like", "normallike") ~ "Normal",
    stringr::str_to_lower(y) %in% c("claudin-low", "claudinlow") ~ "claudin-low",
    TRUE ~ y
  )

  factor(
    y,
    levels = c("LumA", "LumB", "Her2", "Basal", "claudin-low", "Normal")
  )
}

# ------------------------------------------------------------------
# Plot colours
# ------------------------------------------------------------------

pam50_colors <- c(
  "LumA" = "#87CEFA",
  "LumB" = "#CD96CD",
  "Her2" = "#FFB5C5",
  "Basal" = "#CD6090",
  "claudin-low" = "#FFD700",
  "Normal" = "lightskyblue1"
)

# ------------------------------------------------------------------
# Load TCGA-BRCA RNA-seq object
# ------------------------------------------------------------------

brca_se <- readRDS(input_file)
metadata <- as.data.frame(colData(brca_se))

# ------------------------------------------------------------------
# VST transformation
# ------------------------------------------------------------------

dds <- DESeqDataSet(brca_se, design = ~ 1)
vsd <- vst(dds, blind = TRUE)
expr_mat <- assay(vsd)

# ------------------------------------------------------------------
# Map Ensembl IDs to gene symbols
# ------------------------------------------------------------------

ensembl_ids <- sub("\\..*$", "", rownames(expr_mat))

symbols <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = ensembl_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

expr_symbol <- aggregate_duplicate_symbols(expr_mat, symbols)

# ------------------------------------------------------------------
# Extract genes of interest
# ------------------------------------------------------------------

matched_genes <- intersect(pca_genes, rownames(expr_symbol))
missing_genes <- setdiff(pca_genes, matched_genes)

message("Matched genes: ", paste(matched_genes, collapse = ", "))

if (length(missing_genes) > 0) {
  message("Missing genes: ", paste(missing_genes, collapse = ", "))
}

if (length(matched_genes) < 2) {
  stop("Fewer than two PCA genes were found in the expression matrix.")
}

expr_df <- t(expr_symbol[matched_genes, , drop = FALSE]) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample")

# ------------------------------------------------------------------
# Keep primary tumor samples only and collapse to patient level
# ------------------------------------------------------------------

tumor_samples <- rownames(metadata)[metadata$sample_type == "Primary Tumor"]

expr_tumor <- expr_df %>%
  filter(sample %in% tumor_samples) %>%
  mutate(patient = substr(sample, 1, 12))

expr_patient <- expr_tumor %>%
  group_by(patient) %>%
  summarise(
    across(all_of(matched_genes), ~ median(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# ------------------------------------------------------------------
# Fetch PAM50 subtype annotations
# ------------------------------------------------------------------

pam50_table <- TCGAquery_subtype("BRCA") %>%
  mutate(patient = substr(patient, 1, 12)) %>%
  transmute(
    patient,
    PAM50_raw = BRCA_Subtype_PAM50,
    PAM50 = recode_pam50(BRCA_Subtype_PAM50)
  )

# ------------------------------------------------------------------
# Join expression with PAM50 annotation
# ------------------------------------------------------------------

brca_pca_data <- expr_patient %>%
  inner_join(pam50_table, by = "patient") %>%
  filter(!is.na(PAM50)) %>%
  mutate(PAM50 = droplevels(PAM50))

message("PAM50 subtype counts:")
print(table(brca_pca_data$PAM50))

# ------------------------------------------------------------------
# PCA
# ------------------------------------------------------------------

pca_input <- brca_pca_data %>%
  dplyr::select(all_of(matched_genes)) %>%
  as.matrix()

rownames(pca_input) <- brca_pca_data$patient

pca_result <- prcomp(pca_input, center = TRUE, scale. = TRUE)
variance_explained <- summary(pca_result)$importance[2, 1:2] * 100

pca_plot_data <- data.frame(
  patient = brca_pca_data$patient,
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  PAM50 = brca_pca_data$PAM50
)

# ------------------------------------------------------------------
# Plot PCA by PAM50 subtype
# ------------------------------------------------------------------

p_pam50 <- ggplot(
  pca_plot_data,
  aes(x = PC1, y = PC2, color = PAM50)
) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_manual(values = pam50_colors, drop = FALSE) +
  labs(
    title = "TCGA-BRCA PCA by PAM50 subtype",
    subtitle = "Primary tumor samples only; TTLL/CCP polyglutamylation-related genes",
    x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
    color = "PAM50 subtype"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right",
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.grid = element_blank()
  )

print(p_pam50)

# ------------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------------

