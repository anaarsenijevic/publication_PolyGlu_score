#!/usr/bin/env Rscript

# TCGA-BRCA PAM50 association with pan-cancer RF score
#
# This script:
# 1. trains pan-cancer and BRCA random forest classifiers using matched TCGA tumour-normal cohorts
# 2. scores TCGA-BRCA tumour patients using the fixed 9-gene panel
# 3. joins PAM50 subtype annotations from TCGAbiolinks
# 4. compares PAM50 subtype proportions in High vs Low quartiles of the pan-cancer RF score
# 5. performs a Chi-square test and generates a stacked proportion plot
# 6. saves PAM50 subtype percentages for each RF-score group
#
# Input:
# - dds_*_matched_group.rds files
# - TCGA-BRCA_RNAseq.rds
#
# Notes:
# - Edit all input paths below before running
# - No z-scoring is applied
# - The fixed 9-gene panel is:
#   TTLL4, TTLL5, TTLL6, TTLL7, TTLL11, AGBL1, AGBL3, AGBL4, AGBL5
# - High and Low groups are defined as the top and bottom quartiles of rf_prob_pan

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(biomaRt)
  library(caret)
  library(randomForest)
  library(scales)
  library(stringr)
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

set.seed(123)

# ------------------------------------------------------------------
# Input paths
# ------------------------------------------------------------------

tcga_dir <- "path-to-tcga-dds-files"
gdc_dir <- "path-to-gdc-rds-files"
out_dir <- "path-to-output-folder"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(tcga_dir)) {
  stop("Please set 'tcga_dir' to the folder containing dds_*_matched_group.rds files.")
}

if (!dir.exists(gdc_dir)) {
  stop("Please set 'gdc_dir' to the folder containing TCGA-BRCA_RNAseq.rds.")
}

brca_rds <- file.path(gdc_dir, "TCGA-BRCA_RNAseq.rds")

if (!file.exists(brca_rds)) {
  stop("Could not find: ", brca_rds)
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

make_tb_group <- function(df, score_col, q = 0.75) {
  hi <- as.numeric(quantile(df[[score_col]], q, na.rm = TRUE))
  lo <- as.numeric(quantile(df[[score_col]], 1 - q, na.rm = TRUE))

  df %>%
    mutate(
      group = case_when(
        .data[[score_col]] >= hi ~ "High",
        .data[[score_col]] <= lo ~ "Low",
        TRUE ~ NA_character_
      ),
      group = factor(group, levels = c("Low", "High"))
    ) %>%
    filter(!is.na(group))
}

fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 1e-4) return(sprintf("p = %.1e", p))
  sprintf("p = %.4f", p)
}

na_str_to_na <- function(x) {
  y <- str_trim(as.character(x))
  y[y %in% c(
    "NA", "N/A", "", "Unknown", "Not Reported", "NOT REPORTED",
    "[Not Available]", "[NOT AVAILABLE]", "not reported"
  )] <- NA_character_
  y
}

recode_pam50_raw <- function(x) {
  y <- na_str_to_na(x)
  y <- str_trim(y)
  y <- gsub("\\s+", "", y)

  y <- case_when(
    str_to_lower(y) %in% c("luma") ~ "LumA",
    str_to_lower(y) %in% c("lumb") ~ "LumB",
    str_to_lower(y) %in% c("basal") ~ "Basal",
    str_to_lower(y) %in% c("her2", "her2enriched", "her2+") ~ "Her2",
    str_to_lower(y) %in% c("normal", "normal-like", "normallike") ~ "Normal",
    TRUE ~ y
  )

  factor(y, levels = c("Basal", "Her2", "LumB", "LumA", "Normal"))
}

pal_PAM50 <- c(
  "Basal" = "#CD6090",
  "Her2" = "#FFB5C5",
  "LumB" = "#CD96CD",
  "LumA" = "#87CEFA",
  "Normal" = "lightskyblue1"
)

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
# Build training matrix from matched TCGA cohorts
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

all_tcga <- list()

for (file in dds_files) {
  project <- sub("dds_(.*)_matched_group\\.rds$", "\\1", basename(file))
  message("  training: ", project)

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

  cond <- tryCatch(colData(dds)[expr$sample, "condition"], error = function(e) NULL)
  if (is.null(cond)) next

  expr$condition <- as.character(cond)

  all_tcga[[length(all_tcga) + 1]] <- expr
}

tcga_df <- bind_rows(all_tcga)

if (nrow(tcga_df) == 0) {
  stop("No training data could be assembled from matched TCGA cohorts.")
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
# Train pan-cancer and BRCA RF models
# ------------------------------------------------------------------

ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = TRUE
)

rf_pan <- caret::train(
  x = tcga_clean[, feature_genes],
  y = tcga_clean$condition_bin,
  method = "rf",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5
)

tcga_brca_train <- tcga_clean %>%
  filter(grepl("BRCA", project, ignore.case = TRUE))

if (nrow(tcga_brca_train) == 0 || length(unique(tcga_brca_train$condition_bin)) < 2) {
  stop("Could not assemble a valid TCGA-BRCA training set.")
}

rf_brca <- caret::train(
  x = tcga_brca_train[, feature_genes],
  y = tcga_brca_train$condition_bin,
  method = "rf",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5
)

saveRDS(rf_pan, file.path(out_dir, "rf_pan_tcga_9gene_RAWVST.rds"))
saveRDS(rf_brca, file.path(out_dir, "rf_brca_tcga_9gene_RAWVST.rds"))

# ------------------------------------------------------------------
# Score TCGA-BRCA tumour patients
# ------------------------------------------------------------------

message("Scoring TCGA-BRCA tumour patients...")

brca_se <- readRDS(brca_rds)
meta <- as.data.frame(colData(brca_se))

dds_all <- DESeqDataSet(brca_se, design = ~1)
vsd_all <- vst(dds_all, blind = TRUE)
mat_all <- assay(vsd_all)

ens <- sub("\\..*$", "", rownames(mat_all))
sym <- gene_map$external_gene_name[match(ens, gene_map$ensembl_gene_id)]
mat_sym <- aggregate_dup_symbols(mat_all, sym, mode = "mean")

matched <- intersect(feature_genes, rownames(mat_sym))

if (length(matched) < 8) {
  stop("Fewer than 8 panel genes were found in TCGA-BRCA RNA-seq data.")
}

expr <- t(mat_sym[matched, , drop = FALSE]) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample")

tumour_samples <- rownames(meta)[meta$sample_type == "Primary Tumor"]

expr_tum <- expr %>%
  filter(sample %in% tumour_samples) %>%
  mutate(patient = substr(sample, 1, 12))

brca_expr_pat <- expr_tum %>%
  group_by(patient) %>%
  summarise(
    across(all_of(feature_genes), ~ median(.x, na.rm = TRUE)),
    .groups = "drop"
  )

X_all <- brca_expr_pat %>%
  dplyr::select(all_of(feature_genes))

tcga_brca_scored <- brca_expr_pat %>%
  mutate(
    rf_prob_pan = as.numeric(predict(rf_pan, newdata = X_all, type = "prob")[, "Tumor"]),
    rf_prob_brca = as.numeric(predict(rf_brca, newdata = X_all, type = "prob")[, "Tumor"])
  )

write.csv(
  tcga_brca_scored,
  file.path(out_dir, "tcga_brca_scored_patients.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Add PAM50 subtype annotations
# ------------------------------------------------------------------

message("Retrieving TCGA-BRCA PAM50 annotations...")

tcga_pam50 <- TCGAquery_subtype("BRCA") %>%
  mutate(patient = substr(patient, 1, 12)) %>%
  transmute(patient, PAM50_src = BRCA_Subtype_PAM50)

tcga_brca_pam50 <- tcga_brca_scored %>%
  left_join(tcga_pam50, by = "patient") %>%
  mutate(PAM50 = recode_pam50_raw(PAM50_src))

cat("\nTCGA-BRCA samples scored:", nrow(tcga_brca_pam50), "\n")
cat("PAM50 non-missing:", sum(!is.na(tcga_brca_pam50$PAM50)), "\n")

write.csv(
  tcga_brca_pam50,
  file.path(out_dir, "tcga_brca_scored_with_PAM50.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# PAM50 association with High vs Low pan-cancer RF score
# ------------------------------------------------------------------

q <- 0.75

d_pan_pam <- make_tb_group(tcga_brca_pam50, "rf_prob_pan", q = q) %>%
  filter(!is.na(PAM50)) %>%
  mutate(
    PAM50 = factor(
      PAM50,
      levels = c("Basal", "Her2", "LumB", "LumA", "Normal")
    )
  )

if (nrow(d_pan_pam) < 10) {
  stop("Too few TCGA-BRCA samples with PAM50 annotations after High/Low grouping.")
}

tb <- table(d_pan_pam$group, d_pan_pam$PAM50)
chi_p <- suppressWarnings(chisq.test(tb)$p.value)

ann <- tibble(
  x = "High",
  y = 0.98,
  label = paste0(fmt_p(chi_p), "\nChi-square (RxC)\n", "n=", sum(tb))
)

g_pan_pam50 <- ggplot(d_pan_pam, aes(x = group, fill = PAM50)) +
  geom_bar(position = "fill", color = "white", linewidth = 0.25) +
  scale_y_continuous(
    labels = scales::percent,
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  geom_text(
    data = ann,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    vjust = 1,
    hjust = 1,
    size = 3.4
  ) +
  scale_fill_manual(values = pal_PAM50, drop = FALSE, na.translate = FALSE) +
  labs(
    title = "TCGA-BRCA: PAM50 subtype vs pan-cancer RF score",
    x = NULL,
    y = "Proportion",
    fill = "PAM50"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

print(g_pan_pam50)

ggsave(
  file.path(out_dir, "TCGA_BRCA_PAM50_vs_PAN_RF_score_q0.75.pdf"),
  g_pan_pam50,
  width = 6.5,
  height = 4.5
)

# ------------------------------------------------------------------
# PAM50 percentages by RF-score group
# ------------------------------------------------------------------

pam50_perc <- d_pan_pam %>%
  count(group, PAM50) %>%
  group_by(group) %>%
  mutate(
    total = sum(n),
    percent = 100 * n / total
  ) %>%
  ungroup()

print(pam50_perc)

write.csv(
  pam50_perc,
  file.path(out_dir, "TCGA_BRCA_PAM50_percentages_by_PAN_RF_group.csv"),
  row.names = FALSE
)

assoc_result <- tibble(
  score = "PAN",
  variable = "PAM50",
  test = "Chi-square (RxC)",
  p = chi_p,
  n = sum(tb)
)

write.csv(
  assoc_result,
  file.path(out_dir, "TCGA_BRCA_PAM50_vs_PAN_RF_score_chisq.csv"),
  row.names = FALSE
)

message("\nFinished. Outputs written to:\n", out_dir)
