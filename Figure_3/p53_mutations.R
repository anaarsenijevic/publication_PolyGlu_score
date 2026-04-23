#!/usr/bin/env Rscript

# Pan-cancer TP53 mutation status in High vs Low RF-score patients
#
# This script:
# 1. reads matched TCGA tumour-normal DESeq2 objects
# 2. trains a pan-cancer random forest classifier on the fixed 9-gene panel
# 3. scores tumour patients across TCGA projects using local RNA-seq RDS files
# 4. retrieves TP53 mutation calls per cancer type from GDC MAF files
# 5. compares TP53 mutation status in High vs Low quartiles of pan-cancer RF score
# 6. performs a Fisher's exact test and generates a pan-cancer stacked bar plot
#
# Input:
# - dds_*_matched_group.rds files
# - TCGA-<PROJECT>_RNAseq.rds files
# - GDC MAF files or downloadable GDC mutation data
#
# Notes:
# - Edit all input paths below before running
# - Only tumour samples are used for the pan-cancer scoring analysis
# - The fixed 9-gene panel is:
#   TTLL4, TTLL5, TTLL6, TTLL7, TTLL11, AGBL1, AGBL3, AGBL4, AGBL5

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(biomaRt)
  library(caret)
  library(stringr)
  library(TCGAbiolinks)
  library(scales)
  library(readr)
})

set.seed(123)

# ------------------------------------------------------------------
# Input paths
# ------------------------------------------------------------------

tcga_dir <- "path-to-tcga-dds-files"
gdc_dir <- "path-to-gdc-rnaseq-rds-files"
cache_root <- "path-to-gdc-maf-cache"
out_dir <- "path-to-output-folder"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(tcga_dir)) {
  stop("Please set 'tcga_dir' to the folder containing dds_*_matched_group.rds files.")
}
if (!dir.exists(gdc_dir)) {
  stop("Please set 'gdc_dir' to the folder containing TCGA-<PROJECT>_RNAseq.rds files.")
}

# ------------------------------------------------------------------
# Fixed 9-gene panel
# ------------------------------------------------------------------

feature_genes <- c(
  "TTLL4", "TTLL5", "TTLL6", "TTLL7", "TTLL11",
  "AGBL1", "AGBL3", "AGBL4", "AGBL5"
)

pal_p53 <- c("WT" = "#CAE1FF", "Mutated" = "lightsteelblue4")

# ------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------

aggregate_dup_symbols <- function(mat, symbols, mode = "mean") {
  keep <- !is.na(symbols) & symbols != ""
  mat <- mat[keep, , drop = FALSE]
  symbols <- symbols[keep]

  if (mode == "sum") {
    out <- rowsum(mat, group = symbols)
  } else {
    split_idx <- split(seq_len(nrow(mat)), symbols)
    out <- vapply(
      split_idx,
      function(ix) colMeans(mat[ix, , drop = FALSE]),
      numeric(ncol(mat))
    )
    out <- t(out)
    rownames(out) <- names(split_idx)
  }

  out
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

fisher_2x2 <- function(df, group_col = "group", var_col = "p53_status") {
  d <- df %>% filter(!is.na(.data[[group_col]]), !is.na(.data[[var_col]]))

  if (nrow(d) < 10) {
    return(list(p = NA_real_, n = nrow(d)))
  }

  x <- droplevels(factor(d[[group_col]], levels = c("Low", "High")))
  y <- droplevels(factor(d[[var_col]], levels = c("WT", "Mutated")))

  if (nlevels(x) < 2 || nlevels(y) < 2) {
    return(list(p = NA_real_, n = nrow(d)))
  }

  tb <- table(x, y)
  tb <- tb[c("Low", "High"), c("WT", "Mutated"), drop = FALSE]

  list(
    p = fisher.test(tb)$p.value,
    n = sum(tb),
    tb = tb
  )
}

plot_tp53_hi_lo <- function(df,
                            score_col = "rf_prob_pan",
                            q = 0.75,
                            title = "TCGA-Pan: TP53 status vs RF pan-cancer score",
                            palette = pal_p53,
                            base_size = 13) {
  d <- make_tb_group(df, score_col, q = q) %>%
    mutate(p53_status = factor(p53_status, levels = c("WT", "Mutated"))) %>%
    filter(!is.na(p53_status))

  if (nrow(d) < 10) {
    return(ggplot() + theme_void() + labs(title = "Too few data"))
  }

  ft <- fisher_2x2(d, "group", "p53_status")

  ann <- tibble(
    x = "High",
    y = 0.98,
    label = paste0(fmt_p(ft$p), "\nFisher (2x2)\n", "n = ", ft$n)
  )

  ggplot(d, aes(x = group, fill = p53_status)) +
    geom_bar(position = "fill", color = "white", linewidth = 0.25) +
    scale_y_continuous(
      labels = percent,
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    scale_fill_manual(values = palette, drop = FALSE, na.translate = FALSE) +
    geom_text(
      data = ann,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      vjust = 1,
      size = 3.4
    ) +
    labs(
      title = title,
      x = NULL,
      y = "Proportion",
      fill = "TP53"
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold")
    )
}

get_tp53_patients_one_project <- function(project_id,
                                          cache_root,
                                          verbose = TRUE,
                                          workflow_type = "Aliquot Ensemble Somatic Variant Merging and Masking",
                                          force_download = FALSE) {
  suppressPackageStartupMessages({
    library(TCGAbiolinks)
    library(dplyr)
    library(readr)
  })

  proj_dir <- file.path(cache_root, project_id)
  dir.create(proj_dir, recursive = TRUE, showWarnings = FALSE)

  cache_rds <- file.path(proj_dir, paste0(project_id, "_TP53_mut_patients.rds"))

  if (!force_download && file.exists(cache_rds)) {
    if (verbose) message("TP53 MAF: ", project_id, " | using cached RDS")
    return(readRDS(cache_rds))
  }

  maf_paths <- list.files(
    proj_dir,
    pattern = "\\.maf(\\.gz)?$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (!force_download && length(maf_paths) > 0) {
    if (verbose) message("TP53 MAF: ", project_id, " | MAFs already present, parsing only")
  } else {
    if (verbose) message("TP53 MAF: ", project_id, " | downloading")

    q <- tryCatch(
      GDCquery(
        project = project_id,
        data.category = "Simple Nucleotide Variation",
        data.type = "Masked Somatic Mutation",
        workflow.type = workflow_type
      ),
      error = function(e) NULL
    )

    if (is.null(q) || nrow(q$results[[1]]) == 0) {
      saveRDS(NULL, cache_rds)
      return(NULL)
    }

    oldwd <- getwd()
    on.exit(setwd(oldwd), add = TRUE)
    setwd(proj_dir)

    GDCdownload(q, directory = proj_dir, method = "api")

    tars <- list.files(proj_dir, pattern = "\\.tar\\.gz$", full.names = TRUE)
    if (length(tars) > 0) {
      newest_tar <- tars[which.max(file.info(tars)$mtime)]
      if (verbose) message("  extracting: ", basename(newest_tar))
      utils::untar(newest_tar, exdir = proj_dir)
    }

    maf_paths <- list.files(
      proj_dir,
      pattern = "\\.maf(\\.gz)?$",
      recursive = TRUE,
      full.names = TRUE
    )

    if (!length(maf_paths)) {
      saveRDS(NULL, cache_rds)
      return(NULL)
    }
  }

  read_one_maf <- function(p) {
    x <- suppressWarnings(
      readr::read_tsv(p, show_col_types = FALSE, progress = FALSE, comment = "#")
    )

    if (!all(c("Hugo_Symbol", "Tumor_Sample_Barcode") %in% names(x))) {
      return(NULL)
    }

    x %>%
      filter(Hugo_Symbol == "TP53") %>%
      transmute(
        patient = substr(as.character(Tumor_Sample_Barcode), 1, 12),
        cancer_type = sub("^TCGA-", "", project_id),
        p53_status = "Mutated"
      )
  }

  tp53 <- bind_rows(lapply(maf_paths, read_one_maf)) %>%
    filter(!is.na(patient), patient != "") %>%
    distinct(patient, cancer_type, .keep_all = TRUE)

  if (!nrow(tp53)) {
    tp53 <- NULL
  }

  saveRDS(tp53, cache_rds)
  tp53
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
# Load matched TCGA cohorts and train pan-cancer RF
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

tcga_df <- map_dfr(dds_files, function(f) {
  project <- sub("dds_(.*)_matched_group\\.rds$", "\\1", basename(f))
  message("  training: ", project)

  dds <- readRDS(f)
  vsd <- vst(dds, blind = TRUE)
  mat <- assay(vsd)

  ens <- sub("\\..*$", "", rownames(mat))
  sym <- gene_map$external_gene_name[match(ens, gene_map$ensembl_gene_id)]
  mat_sym <- aggregate_dup_symbols(mat, sym, mode = "mean")

  matched <- intersect(feature_genes, rownames(mat_sym))
  if (!length(matched)) {
    return(tibble())
  }

  expr <- t(mat_sym[matched, , drop = FALSE]) %>% as.data.frame()
  expr$sample <- rownames(expr)
  expr$project <- project

  cond <- tryCatch(colData(dds)[expr$sample, "condition"], error = function(e) NA)
  expr$condition <- as.character(cond)

  as_tibble(expr)
})

if (nrow(tcga_df) == 0) {
  stop("No training data could be assembled from the matched TCGA cohorts.")
}

tcga_clean <- tcga_df %>%
  select(all_of(feature_genes), sample, project, condition) %>%
  filter(if_all(all_of(feature_genes), ~ !is.na(.))) %>%
  mutate(
    condition_bin = factor(
      ifelse(condition == "Primary Tumor", "Tumor", "Normal"),
      levels = c("Normal", "Tumor")
    )
  )

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

saveRDS(rf_pan, file.path(out_dir, "rf_pan_tcga_9gene_RAWVST.rds"))

# ------------------------------------------------------------------
# Score pan-TCGA tumour patients
# ------------------------------------------------------------------

message("Scoring tumour patients across TCGA projects...")
rds_files <- list.files(
  gdc_dir,
  pattern = "^TCGA-[A-Z0-9]+_RNAseq\\.rds$",
  full.names = TRUE
)

if (length(rds_files) == 0) {
  stop("No TCGA RNA-seq RDS files found in 'gdc_dir'.")
}

score_one_project_pan <- function(rds_path) {
  project_id <- sub("_RNAseq\\.rds$", "", basename(rds_path))
  message("  scoring: ", project_id)

  se <- readRDS(rds_path)
  meta <- as.data.frame(colData(se))

  dds <- DESeqDataSet(se, design = ~1)
  vsd <- vst(dds, blind = TRUE)
  mat <- assay(vsd)

  ens <- sub("\\..*$", "", rownames(mat))
  sym <- gene_map$external_gene_name[match(ens, gene_map$ensembl_gene_id)]
  mat_sym <- aggregate_dup_symbols(mat, sym, mode = "mean")

  matched <- intersect(feature_genes, rownames(mat_sym))
  if (length(matched) < 8) {
    return(NULL)
  }

  expr <- t(mat_sym[matched, , drop = FALSE]) %>%
    as.data.frame() %>%
    rownames_to_column("sample")

  tumor_samples <- rownames(meta)[meta$sample_type == "Primary Tumor"]
  if (!length(tumor_samples)) {
    return(NULL)
  }

  pat <- expr %>%
    filter(sample %in% tumor_samples) %>%
    mutate(
      patient = substr(sample, 1, 12),
      cancer_type = sub("^TCGA-", "", project_id)
    ) %>%
    group_by(patient, cancer_type) %>%
    summarise(across(all_of(feature_genes), ~ median(.x, na.rm = TRUE)), .groups = "drop")

  X <- pat %>% select(all_of(feature_genes))

  pat %>%
    mutate(
      rf_prob_pan = as.numeric(predict(rf_pan, newdata = X, type = "prob")[, "Tumor"])
    )
}

tcga_pan_scored <- bind_rows(lapply(rds_files, score_one_project_pan))

if (nrow(tcga_pan_scored) == 0) {
  stop("No tumour patients were scored across TCGA projects.")
}

write.csv(
  tcga_pan_scored,
  file.path(out_dir, "tcga_pan_scored_patients.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Retrieve TP53 mutation calls from GDC MAF files
# ------------------------------------------------------------------

projects <- unique(tcga_pan_scored$cancer_type)
project_ids <- paste0("TCGA-", projects)

tp53_list <- lapply(
  project_ids,
  get_tp53_patients_one_project,
  cache_root = cache_root,
  verbose = TRUE,
  force_download = FALSE
)

tp53_pancan_by_cancer <- bind_rows(tp53_list) %>%
  filter(!is.na(patient), patient != "") %>%
  distinct(patient, cancer_type, .keep_all = TRUE)

tcga_pan2 <- tcga_pan_scored %>%
  left_join(
    tp53_pancan_by_cancer %>%
      transmute(patient, cancer_type, p53_status = "Mutated"),
    by = c("patient", "cancer_type")
  ) %>%
  mutate(
    p53_status = ifelse(is.na(p53_status), "WT", p53_status),
    p53_status = factor(p53_status, levels = c("WT", "Mutated"))
  )

write.csv(
  tcga_pan2,
  file.path(out_dir, "tcga_pan_scored_with_TP53.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Define High and Low RF-score groups
# ------------------------------------------------------------------

tcga_pan_hilo <- make_tb_group(tcga_pan2, "rf_prob_pan", q = 0.75)

write.csv(
  tcga_pan_hilo,
  file.path(out_dir, "tcga_pan_scored_with_TP53_highlow.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Plot overall pan-cancer TP53 status in High vs Low RF-score groups
# ------------------------------------------------------------------

g_pan_overall <- plot_tp53_hi_lo(
  tcga_pan2,
  score_col = "rf_prob_pan",
  q = 0.75,
  title = "TCGA-Pan: TP53 mutation status vs RF pan-cancer score",
  palette = pal_p53
)

print(g_pan_overall)

ggsave(
  file.path(out_dir, "PanTCGA_TP53_vs_PANscore_overall.png"),
  g_pan_overall,
  width = 6,
  height = 4,
  dpi = 300
)
