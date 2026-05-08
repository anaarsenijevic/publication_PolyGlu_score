# TCGA-BRCA tumor-normal PCA within PAM50 subtypes
#
# This script:
# - Loads a matched TCGA-BRCA DESeq2 object containing tumor and normal samples
# - Performs variance stabilizing transformation (VST)
# - Maps Ensembl gene IDs to gene symbols
# - Extracts TTLL and CCP genes involved in microtubule polyglutamylation
# - Assigns PAM50 subtype annotations to each patient
# - Allows matched normal samples to inherit the PAM50 subtype of the paired tumor
# - Filters for PAM50 subtypes containing both tumor and normal samples
# - Performs PCA separately within each PAM50 subtype
# - Visualizes tumor-normal separation within each subtype
#
# Input:
# - dds_TCGA-BRCA_matched_group.rds
#
# Output:
# - Subtype-specific tumor-normal PCA plot (.pdf, .png)
# - Expression table with PAM50 annotations
# - PCA scores table
#
# Notes:
# - PCA is performed independently within each PAM50 subtype
# - Only subtypes with both tumor and normal samples are retained
# - The analysis uses TTLL/CCP polyglutamylation-related genes
# - Normal samples are assigned the PAM50 subtype of their matched tumor patient

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(TCGAbiolinks)
  library(ggplot2)
})

set.seed(123)

# ------------------------------------------------------------------
# Input and output paths
# ------------------------------------------------------------------

input_file <- "path-to-your-file/dds_TCGA-BRCA_matched_group.rds"

if (!file.exists(input_file)) {
  stop("Please set 'input_file' to the location of dds_TCGA-BRCA_matched_group.rds.")
}


# ------------------------------------------------------------------
# Gene set
# ------------------------------------------------------------------

ttll_genes <- c(
  "TTLL1", "TTLL4", "TTLL5", "TTLL6", "TTLL7",
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
    stringr::str_to_lower(y) %in% c("claudin-low", "claudinlow") ~ "claudin-low",
    stringr::str_to_lower(y) %in% c("normal", "normal-like", "normallike") ~ "Normal",
    stringr::str_to_lower(y) == "nc" ~ "NC",
    TRUE ~ y
  )

  factor(
    y,
    levels = c("LumA", "LumB", "Her2", "Basal", "claudin-low", "Normal", "NC")
  )
}

# ------------------------------------------------------------------
# Plot colours
# ------------------------------------------------------------------

condition_colors <- c(
  "Normal" = "lightskyblue",
  "Tumor" = "#FF6A6A"
)

# ------------------------------------------------------------------
# Load matched TCGA-BRCA DESeq2 object
# ------------------------------------------------------------------

dds_brca <- readRDS(input_file)

message("Condition counts in matched DESeq2 object:")
print(table(colData(dds_brca)$condition, useNA = "ifany"))

# ------------------------------------------------------------------
# VST transformation
# ------------------------------------------------------------------

vsd_brca <- vst(dds_brca, blind = TRUE)
expr_mat <- assay(vsd_brca)

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
  tibble::rownames_to_column("sample") %>%
  mutate(
    condition = as.character(colData(dds_brca)[sample, "condition"]),
    condition_bin = case_when(
      condition == "Primary Tumor" ~ "Tumor",
      condition %in% c("Solid Tissue Normal", "Normal") ~ "Normal",
      TRUE ~ NA_character_
    ),
    patient = substr(sample, 1, 12)
  ) %>%
  filter(!is.na(condition_bin)) %>%
  mutate(
    condition_bin = factor(condition_bin, levels = c("Normal", "Tumor"))
  )

message("Expression table condition counts:")
print(table(expr_df$condition_bin, useNA = "ifany"))

# ------------------------------------------------------------------
# Fetch PAM50 subtype annotations
# ------------------------------------------------------------------

pam50_table <- TCGAquery_subtype("BRCA") %>%
  mutate(patient = substr(patient, 1, 12)) %>%
  transmute(
    patient,
    PAM50_raw = BRCA_Subtype_PAM50,
    PAM50 = recode_pam50(BRCA_Subtype_PAM50)
  ) %>%
  distinct(patient, .keep_all = TRUE)

message("PAM50 counts from TCGAquery_subtype:")
print(table(pam50_table$PAM50, useNA = "ifany"))

# ------------------------------------------------------------------
# Join PAM50 annotations to matched tumor and normal samples
# ------------------------------------------------------------------

expr_pam50 <- expr_df %>%
  inner_join(pam50_table, by = "patient") %>%
  filter(!is.na(PAM50), PAM50 != "NC") %>%
  mutate(PAM50 = droplevels(PAM50))

message("Samples by PAM50 and condition:")
print(table(expr_pam50$PAM50, expr_pam50$condition_bin, useNA = "ifany"))

# ------------------------------------------------------------------
# Keep only PAM50 subtypes with both tumor and normal samples
# ------------------------------------------------------------------

counts_df <- expr_pam50 %>%
  count(PAM50, condition_bin) %>%
  tidyr::complete(PAM50, condition_bin, fill = list(n = 0)) %>%
  tidyr::pivot_wider(
    names_from = condition_bin,
    values_from = n,
    values_fill = 0
  )

if (!"Tumor" %in% colnames(counts_df)) {
  counts_df$Tumor <- 0
}

if (!"Normal" %in% colnames(counts_df)) {
  counts_df$Normal <- 0
}

valid_subtypes <- counts_df %>%
  filter(Tumor > 0, Normal > 0) %>%
  pull(PAM50)

expr_pam50_valid <- expr_pam50 %>%
  filter(PAM50 %in% valid_subtypes) %>%
  mutate(PAM50 = droplevels(PAM50))

message("Valid PAM50 subtypes with both tumor and normal samples:")
print(valid_subtypes)

message("Final samples by PAM50 and condition:")
print(table(expr_pam50_valid$PAM50, expr_pam50_valid$condition_bin))

if (nrow(expr_pam50_valid) == 0) {
  stop("No PAM50 subtype contains both tumor and normal samples.")
}

# ------------------------------------------------------------------
# PCA separately within each PAM50 subtype
# ------------------------------------------------------------------

per_subtype_pca_data <- purrr::map_dfr(
  levels(expr_pam50_valid$PAM50),
  function(subtype) {

    subtype_data <- expr_pam50_valid %>%
      filter(PAM50 == subtype)

    if (nrow(subtype_data) < 4) {
      return(NULL)
    }

    if (length(unique(subtype_data$condition_bin)) < 2) {
      return(NULL)
    }

    pca_input <- subtype_data %>%
      dplyr::select(all_of(matched_genes)) %>%
      as.matrix()

    rownames(pca_input) <- subtype_data$sample

    pca_result <- prcomp(pca_input, center = TRUE, scale. = TRUE)
    variance_explained <- summary(pca_result)$importance[2, 1:2] * 100

    tibble(
      sample = subtype_data$sample,
      patient = subtype_data$patient,
      PC1 = pca_result$x[, 1],
      PC2 = pca_result$x[, 2],
      condition_bin = subtype_data$condition_bin,
      PAM50 = subtype,
      facet_label = paste0(
        subtype,
        "\nPC1: ", round(variance_explained[1], 1),
        "% | PC2: ", round(variance_explained[2], 1), "%"
      )
    )
  }
)

if (nrow(per_subtype_pca_data) == 0) {
  stop("No subtype had enough samples for PCA.")
}

# ------------------------------------------------------------------
# Plot subtype-specific PCA
# ------------------------------------------------------------------

p_per_subtype <- ggplot(
  per_subtype_pca_data,
  aes(x = PC1, y = PC2, color = condition_bin)
) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_manual(values = condition_colors) +
  facet_wrap(~ facet_label, scales = "free") +
  labs(
    title = "TCGA-BRCA subtype-specific PCA",
    subtitle = "Separate PCA within each PAM50 subtype; matched tumor-normal samples",
    x = "PC1",
    y = "PC2",
    color = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom",
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(p_per_subtype)
