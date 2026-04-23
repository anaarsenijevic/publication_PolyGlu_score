#!/usr/bin/env Rscript

# Pan-cancer aneuploidy score analysis for the 9-gene TTLL/AGBL RF model
#
# This script:
# 1. reads matched TCGA tumour-normal DESeq2 objects
# 2. trains a pan-cancer random forest classifier on the fixed 9-gene panel
# 3. scores tumour samples across TCGA projects using GDC RNA-seq RDS files
# 4. reads aneuploidy scores from the published supplementary Excel table
# 5. compares Aneuploidy Score (AS) in High vs Low quartiles of RF pan-cancer score
# 6. performs a Wilcoxon rank-sum test and generates a pan-cancer boxplot
#
# Input:
# - dds_*_matched_group.rds files
# - TCGA-<PROJECT>_RNAseq.rds files
# - NIHMS958047-supplement-1.xlsx
#
# Notes:
# - Edit all input paths below before running
# - Only tumour samples are used for the AS analysis
# - The fixed 9-gene panel is:
#   TTLL4, TTLL5, TTLL6, TTLL7, TTLL11, AGBL1, AGBL3, AGBL4, AGBL5
# - Aneuploidy scores were obtained directly from Supplementary Table 2
#   in Taylor et al., 2018

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(biomaRt)
  library(caret)
  library(randomForest)
  library(readxl)
  library(stringr)
  library(ggplot2)
})

set.seed(123)

# ------------------------------------------------------------------
# Input paths
# ------------------------------------------------------------------

tcga_dir <- "path-to-tcga-dds-files"
gdc_root <- "path-to-gdc-rds-files"
supp_excel <- "path-to-aneuploidy-excel/NIHMS958047-supplement-1.xlsx"
out_dir <- "path-to-output-folder"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(tcga_dir)) {
  stop("Please set 'tcga_dir' to the folder containing dds_*_matched_group.rds files.")
}
if (!dir.exists(gdc_root)) {
  stop("Please set 'gdc_root' to the folder containing TCGA-<PROJECT>_RNAseq.rds files.")
}
if (!file.exists(supp_excel)) {
  stop("Could not find: ", supp_excel)
}

# ------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------

q <- 0.75

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

find_col <- function(df, patterns) {
  nms <- names(df)
  for (pat in patterns) {
    hit <- which(grepl(pat, nms, ignore.case = TRUE))
    if (length(hit)) return(nms[hit[1]])
  }
  NA_character_
}

tcga_sample_type_code <- function(x) {
  s <- as.character(x)
  code <- str_match(s, "-(\\d\\d)[A-Za-z]?$")[, 2]
  suppressWarnings(as.integer(code))
}

is_tumour_tcga <- function(x) {
  code <- tcga_sample_type_code(x)
  !is.na(code) & code >= 1 & code <= 9
}

make_highlow_group <- function(df, score_col, q = 0.75) {
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

wilcox_highlow <- function(d, value_col = "AS") {
  d <- d %>% filter(!is.na(group), !is.na(.data[[value_col]]))

  low <- d %>% filter(group == "Low") %>% pull(.data[[value_col]])
  high <- d %>% filter(group == "High") %>% pull(.data[[value_col]])
  w <- suppressWarnings(wilcox.test(high, low, exact = FALSE))

  tibble(
    n = nrow(d),
    low_n = length(low),
    high_n = length(high),
    low_med = median(low, na.rm = TRUE),
    high_med = median(high, na.rm = TRUE),
    effect = median(high, na.rm = TRUE) - median(low, na.rm = TRUE),
    p = w$p.value
  )
}

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-4) return(sprintf("%.1e", p))
  sprintf("%.4f", p)
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
# Train pan-cancer RF model from matched TCGA cohorts
# ------------------------------------------------------------------

message("Loading TCGA matched DESeq2 objects...")
dds_files <- list.files(
  tcga_dir,
  pattern = "dds_.*_matched_group\\.rds$",
  full.names = TRUE
)

if (length(dds_files) == 0) {
  stop("No TCGA DESeq2 files found in 'tcga_dir'.")
}

tcga_list <- list()

for (f in dds_files) {
  proj <- sub("dds_(.*)_matched_group\\.rds$", "\\1", basename(f))
  message("  training: ", proj)

  dds <- readRDS(f)
  vsd <- vst(dds, blind = TRUE)
  mat <- assay(vsd)

  ens <- sub("\\..*$", "", rownames(mat))
  syms <- gene_map$external_gene_name[match(ens, gene_map$ensembl_gene_id)]
  mat_sym <- aggregate_dup_symbols(mat, syms)

  matched <- intersect(feature_genes, rownames(mat_sym))
  if (!length(matched)) next

  expr <- as.data.frame(t(mat_sym[matched, , drop = FALSE]))
  expr$sample <- rownames(expr)

  cond <- tryCatch(colData(dds)[expr$sample, "condition"], error = function(e) NULL)
  if (is.null(cond)) next

  expr$condition_bin <- factor(
    ifelse(as.character(cond) == "Primary Tumor", "Tumor", "Normal"),
    levels = c("Normal", "Tumor")
  )

  tcga_list[[proj]] <- expr
}

tcga_train <- bind_rows(tcga_list)

if (nrow(tcga_train) == 0) {
  stop("No training samples were assembled from the TCGA matched cohorts.")
}

ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

rf_pan <- caret::train(
  x = tcga_train[, feature_genes],
  y = tcga_train$condition_bin,
  method = "rf",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5
)

saveRDS(rf_pan, file.path(out_dir, "rf_pan_tcga_9gene_RAWVST.rds"))

# ------------------------------------------------------------------
# Read Aneuploidy Score table
# ------------------------------------------------------------------

message("Reading Aneuploidy Score table...")
preview <- read_excel(
  supp_excel,
  sheet = 1,
  col_names = FALSE,
  .name_repair = "minimal",
  n_max = 250
)

header_row <- which(
  apply(preview, 1, function(r) {
    rr <- tolower(trimws(as.character(unlist(r))))
    any(grepl("^sample$", rr)) &&
      any(grepl("aneuploidy.*score|aneuploidyscore|\\bas\\b", rr))
  })
)[1]

if (is.na(header_row)) {
  stop("Could not detect the header row in the supplementary Excel file.")
}

supp <- read_excel(
  supp_excel,
  sheet = 1,
  skip = header_row - 1,
  col_names = TRUE
) %>% as.data.frame()

sample_col <- find_col(supp, c("^sample$", "sample\\s*id", "tcga"))
as_col <- find_col(supp, c(
  "aneuploidy.*score",
  "aneuploidyscore",
  "aneuploidy\\s*score\\(as\\)",
  "as\\)$",
  "\\bas\\b"
))

if (is.na(sample_col) || is.na(as_col)) {
  stop("Could not detect sample and/or AS columns in the supplementary table.")
}

as_tbl <- supp %>%
  transmute(
    sample = as.character(.data[[sample_col]]),
    AS = suppressWarnings(as.numeric(.data[[as_col]]))
  ) %>%
  filter(!is.na(sample), !is.na(AS)) %>%
  mutate(
    sample = str_replace_all(sample, "\\.", "-"),
    sample = str_replace_all(sample, "_", "-"),
    tumour = is_tumour_tcga(sample),
    patient = substr(sample, 1, 12)
  ) %>%
  filter(tumour) %>%
  select(patient, AS) %>%
  distinct(patient, .keep_all = TRUE)

write.csv(
  as_tbl,
  file.path(out_dir, "AS_table_from_excel_tumourOnly.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Score TCGA tumour samples across projects
# ------------------------------------------------------------------

message("Scoring TCGA tumour samples...")
rds_files <- list.files(
  gdc_root,
  pattern = "^TCGA-[A-Z0-9]+_RNAseq\\.rds$",
  full.names = TRUE
)

if (length(rds_files) == 0) {
  stop("No TCGA project RNA-seq RDS files found in 'gdc_root'.")
}

score_one_project <- function(rds_path) {
  proj <- sub("^TCGA-([A-Z0-9]+)_RNAseq\\.rds$", "\\1", basename(rds_path))
  message("  scoring: ", proj)

  se <- readRDS(rds_path)
  meta <- as.data.frame(colData(se))

  dds <- DESeqDataSet(se, design = ~1)
  vsd <- vst(dds, blind = TRUE)
  mat <- assay(vsd)

  ens <- sub("\\..*$", "", rownames(mat))
  syms <- gene_map$external_gene_name[match(ens, gene_map$ensembl_gene_id)]
  mat_sym <- aggregate_dup_symbols(mat, syms)

  matched <- intersect(feature_genes, rownames(mat_sym))
  if (length(matched) < 8) return(NULL)

  expr <- as.data.frame(t(mat_sym[matched, , drop = FALSE])) %>%
    tibble::rownames_to_column("sample")

  if ("sample_type" %in% names(meta)) {
    tumour_samples <- rownames(meta)[meta$sample_type == "Primary Tumor"]
    if (length(tumour_samples) == 0) {
      tumour_samples <- expr$sample[is_tumour_tcga(expr$sample)]
    }
  } else {
    tumour_samples <- expr$sample[is_tumour_tcga(expr$sample)]
  }

  expr_tum <- expr %>%
    filter(sample %in% tumour_samples) %>%
    mutate(patient = substr(sample, 1, 12)) %>%
    group_by(patient) %>%
    summarise(across(all_of(feature_genes), ~median(.x, na.rm = TRUE)), .groups = "drop")

  if (nrow(expr_tum) < 10) return(NULL)

  X <- expr_tum %>% select(all_of(feature_genes))

  expr_tum %>%
    mutate(
      rf_prob_pan = as.numeric(predict(rf_pan, newdata = X, type = "prob")[, "Tumor"])
    )
}

tcga_scored <- bind_rows(lapply(rds_files, score_one_project))

if (nrow(tcga_scored) == 0) {
  stop("No scored tumour samples were produced from the TCGA project RDS files.")
}

# ------------------------------------------------------------------
# Join RF scores with Aneuploidy Score
# ------------------------------------------------------------------

dat <- tcga_scored %>%
  inner_join(as_tbl, by = "patient")

if (nrow(dat) == 0) {
  stop("No overlapping patients were found between scored tumours and AS table.")
}

write.csv(
  dat,
  file.path(out_dir, "tcga_scored_tumours_with_AS.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# High vs Low RF pan-cancer score comparison
# ------------------------------------------------------------------

dat_hl <- make_highlow_group(dat, "rf_prob_pan", q = q)

stats <- wilcox_highlow(dat_hl, value_col = "AS") %>%
  mutate(
    score_name = "RF_PAN",
    q = q,
    p_adj = p.adjust(p, method = "BH")
  )

write.csv(
  stats,
  file.path(out_dir, "AS_highlow_PANCAN_RFscore_stats.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------

label_txt <- with(
  stats,
  paste0(
    "Wilcoxon\n",
    "p=", fmt_p(p), "\n",
    "BH=", fmt_p(p_adj), "\n",
    "Delta median=", signif(effect, 3), "\n",
    "n=", n, " (Low=", low_n, ", High=", high_n, ")"
  )
)

p <- ggplot(dat_hl, aes(x = group, y = AS, fill = group)) +
  geom_boxplot(outlier.alpha = 0.25, width = 0.65) +
  scale_fill_manual(values = c("Low" = "pink", "High" = "#8B475D")) +
  annotate("text", x = 2, y = Inf, label = label_txt, vjust = 1.15, hjust = 1, size = 3.6) +
  labs(
    title = paste0("Pan-cancer Aneuploidy Score in High vs Low RF score quartiles (q=", q, ")"),
    x = NULL,
    y = "Aneuploidy Score (AS)",
    fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p)

ggsave(
  file.path(out_dir, "AS_highlow_PANCAN_RFscore_boxplot.pdf"),
  p,
  width = 5.5,
  height = 4.5
)
