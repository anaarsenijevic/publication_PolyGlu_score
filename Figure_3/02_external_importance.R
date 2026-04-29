#!/usr/bin/env Rscript

# External drop-one-gene importance for the 9-gene TTLL/AGBL RF model
#
# This script:
# 1. trains a pan-cancer random forest classifier on matched TCGA tumour-normal cohorts
# 2. loads four external validation datasets:
#    - BRCA patient cohort
#    - GSE87410 lung cohort
#    - CPTAC kidney cohort
#    - Korean CRC cohort
# 3. evaluates the full 9-gene model on each external dataset
# 4. retrains the model after dropping one gene at a time
# 5. computes external importance as:
#      AUC(full model) - AUC(drop-one model)
# 6. plots external importance averaged across datasets and separately per dataset
#
# Input:
# - dds_*_matched_group.rds files
# - BRCA Salmon quant.sf files
# - GENCODE GTF annotation file
# - GSE87410_raw_counts_GRCh38.p13_NCBI.tsv.gz
# - GSM_sample_map.csv
# - CPTAC kidney STAR-count files in GDC_kidney_filtered/
# - CMCBSN_expectedcount_342.txt
#
# Notes:
# - Edit all input paths below before running
# - The fixed 9-gene panel is:
#   TTLL4, TTLL5, TTLL6, TTLL7, TTLL11, AGBL1, AGBL3, AGBL4, AGBL5
# - External importance is positive when dropping a gene reduces AUC
# - CCP aliases are used only for plot labels:
#   AGBL1 = CCP4, AGBL3 = CCP3, AGBL4 = CCP6, AGBL5 = CCP5

suppressPackageStartupMessages({
  library(DESeq2)
  library(tximport)
  library(tidyverse)
  library(biomaRt)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GenomicFeatures)
  library(caret)
  library(randomForest)
  library(pROC)
  library(ggplot2)
  library(forcats)
  library(TCGAbiolinks)
})

set.seed(42)

# ------------------------------------------------------------------
# Input paths
# ------------------------------------------------------------------

tcga_dir <- "path-to-tcga-dds-files"

brca_dir <- "path-to-brca-salmon-folder"
gtf_path <- "path-to-gencode.gtf"

lung_dir <- "path-to-GSE87410-folder"
lung_counts_file <- file.path(lung_dir, "GSE87410_raw_counts_GRCh38.p13_NCBI.tsv.gz")
lung_map_file <- file.path(lung_dir, "GSM_sample_map.csv")

kidney_dir <- "path-to-CPTAC-kidney-folder"
kidney_counts_dir <- file.path(kidney_dir, "GDC_kidney_filtered")

crc_dir <- "path-to-korean-crc-folder"
crc_counts_file <- file.path(crc_dir, "CMCBSN_expectedcount_342.txt")

out_dir <- "path-to-output-folder"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(tcga_dir)) stop("Please set 'tcga_dir'.")
if (!dir.exists(brca_dir)) stop("Please set 'brca_dir'.")
if (!file.exists(gtf_path)) stop("Could not find: ", gtf_path)
if (!file.exists(lung_counts_file)) stop("Could not find: ", lung_counts_file)
if (!file.exists(lung_map_file)) stop("Could not find: ", lung_map_file)
if (!dir.exists(kidney_counts_dir)) stop("Please set 'kidney_dir' so GDC_kidney_filtered/ exists.")
if (!file.exists(crc_counts_file)) stop("Could not find: ", crc_counts_file)

# ------------------------------------------------------------------
# Fixed 9-gene panel and labels
# ------------------------------------------------------------------

feature_genes <- c(
  "TTLL4", "TTLL5", "TTLL6", "TTLL7", "TTLL11",
  "AGBL1", "AGBL3", "AGBL4", "AGBL5"
)

ccp_map <- c(
  "AGTPBP1" = "CCP1",
  "AGBL2" = "CCP2",
  "AGBL3" = "CCP3",
  "AGBL1" = "CCP4",
  "AGBL5" = "CCP5",
  "AGBL4" = "CCP6"
)

to_ccp <- function(x) {
  ifelse(x %in% names(ccp_map), ccp_map[x], x)
}

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

fit_rf_panel <- function(train_df, labels, genes, ctrl, tune_len = 5) {
  caret::train(
    x = train_df[, genes, drop = FALSE],
    y = labels,
    method = "rf",
    metric = "ROC",
    trControl = ctrl,
    tuneLength = tune_len
  )
}

eval_on_set <- function(model, genes, ext_set) {
  missing <- setdiff(genes, colnames(ext_set$expr))

  if (length(missing) > 0) {
    stop(
      "External dataset is missing required genes: ",
      paste(missing, collapse = ", ")
    )
  }

  x <- ext_set$expr[, genes, drop = FALSE]
  x[] <- lapply(x, as.numeric)

  prob <- predict(model, newdata = x, type = "prob")[, "Tumor"]
  truth <- factor(ext_set$truth, levels = c("Normal", "Tumor"))

  as.numeric(pROC::auc(truth, prob))
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
# Load matched TCGA cohorts and train RF
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

for (f in dds_files) {
  project <- sub("dds_(.*)_matched_group\\.rds$", "\\1", basename(f))
  message("  training data: ", project)

  dds <- readRDS(f)
  vsd <- vst(dds, blind = TRUE)
  mat <- assay(vsd)

  ens <- sub("\\..*$", "", rownames(mat))
  sym <- gene_map$external_gene_name[match(ens, gene_map$ensembl_gene_id)]
  mat_sym <- aggregate_dup_symbols(mat, sym, mode = "mean")

  matched <- intersect(feature_genes, rownames(mat_sym))
  if (!length(matched)) next

  expr <- as.data.frame(t(mat_sym[matched, , drop = FALSE]))
  expr$sample <- rownames(expr)
  expr$project <- project

  cond <- tryCatch(
    colData(dds)[expr$sample, "condition"],
    error = function(e) NULL
  )

  if (is.null(cond)) next

  expr$condition <- as.character(cond)
  tcga_list[[project]] <- expr
}

tcga_df <- bind_rows(tcga_list)

if (nrow(tcga_df) == 0) {
  stop("No TCGA training data could be assembled.")
}

tcga_clean <- tcga_df %>%
  filter(!is.na(condition)) %>%
  dplyr::select(all_of(feature_genes), condition, sample, project) %>%
  filter(if_all(all_of(feature_genes), ~ !is.na(.))) %>%
  mutate(
    condition_bin = factor(
      ifelse(condition == "Primary Tumor", "Tumor", "Normal"),
      levels = c("Normal", "Tumor")
    )
  )

y_tcga <- tcga_clean$condition_bin
panel_eval_genes <- intersect(feature_genes, colnames(tcga_clean))

if (length(panel_eval_genes) < 5) {
  stop("Fewer than five panel genes are available in TCGA.")
}

ctrl_pan <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

tune_len <- 5

rf_full <- fit_rf_panel(
  train_df = tcga_clean,
  labels = y_tcga,
  genes = panel_eval_genes,
  ctrl = ctrl_pan,
  tune_len = tune_len
)

saveRDS(
  rf_full,
  file.path(out_dir, "rf_pan_tcga_9gene_RAWVST.rds")
)

# ------------------------------------------------------------------
# External dataset 1: BRCA patient cohort
# ------------------------------------------------------------------

message("Loading BRCA patient cohort...")

samples_brca <- tibble::tribble(
  ~sample, ~tissue, ~group, ~path,
  "1NAT", "NAT", "non-invasive", "salmon_1NAT/quant.sf",
  "2PT", "PT", "non-invasive", "salmon_2PT/quant.sf",
  "3NAT", "NAT", "non-invasive", "salmon_3NAT/quant.sf",
  "4PT", "PT", "non-invasive", "salmon_4PT/quant.sf",
  "5NAT", "NAT", "metastatic", "salmon_5NAT/quant.sf",
  "6PT", "PT", "metastatic", "salmon_6PT/quant.sf",
  "7NAT", "NAT", "metastatic", "salmon_7NAT/quant.sf",
  "8PT", "PT", "metastatic", "salmon_8PT/quant.sf",
  "9NAT", "NAT", "metastatic", "salmon_9NAT/quant.sf",
  "10PT", "PT", "metastatic", "salmon_10PT/quant.sf",
  "11NAT", "NAT", "TNBC", "salmon_11NAT/quant.sf",
  "12PT", "PT", "TNBC", "salmon_12PT/quant.sf",
  "13NAT", "NAT", "TNBC", "salmon_13NAT/quant.sf",
  "14PT", "PT", "TNBC", "salmon_14PT/quant.sf",
  "15NAT", "NAT", "TNBC", "salmon_15NAT/quant.sf",
  "16PT", "PT", "TNBC", "salmon_16PT/quant.sf",
  "17NAT", "NAT", "non-invasive", "salmon_17NAT/quant.sf",
  "18PT", "PT", "non-invasive", "salmon_18PT/quant.sf"
)

files_brca <- setNames(
  file.path(brca_dir, samples_brca$path),
  samples_brca$sample
)

missing_brca <- files_brca[!file.exists(files_brca)]
if (length(missing_brca) > 0) {
  stop("Missing BRCA Salmon files: ", paste(missing_brca, collapse = ", "))
}

txdb <- makeTxDbFromGFF(gtf_path, format = "gtf")

tx2gene <- AnnotationDbi::select(
  txdb,
  keys = keys(txdb, keytype = "TXNAME"),
  columns = "GENEID",
  keytype = "TXNAME"
) %>%
  mutate(
    TXNAME = sub("\\..*", "", TXNAME),
    GENEID = sub("\\..*", "", GENEID)
  ) %>%
  distinct()

txi_b <- tximport(
  files_brca,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE
)

dds_b <- DESeqDataSetFromTximport(
  txi_b,
  colData = samples_brca %>%
    dplyr::select(sample, tissue, group) %>%
    column_to_rownames("sample"),
  design = ~ tissue
)

vsd_b <- vst(dds_b, blind = TRUE)
mat_b <- assay(vsd_b)

ens_b <- sub("\\..*", "", rownames(mat_b))
sym_b <- mapIds(
  org.Hs.eg.db,
  keys = ens_b,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

rownames(mat_b) <- ifelse(is.na(sym_b) | sym_b == "", ens_b, sym_b)

brca_expr_full <- as.data.frame(t(mat_b))
brca_truth <- factor(
  ifelse(samples_brca$tissue == "PT", "Tumor", "Normal"),
  levels = c("Normal", "Tumor")
)

# ------------------------------------------------------------------
# External dataset 2: GSE87410 lung cohort
# ------------------------------------------------------------------

message("Loading GSE87410 lung cohort...")

raw_counts_lung <- readr::read_tsv(
  lung_counts_file,
  show_col_types = FALSE
)

counts_lung <- as.data.frame(raw_counts_lung[, -1])
rownames(counts_lung) <- raw_counts_lung[[1]]
counts_lung[] <- lapply(counts_lung, as.numeric)

gsm_map <- readr::read_csv(
  lung_map_file,
  col_names = c("GSM", "Label"),
  show_col_types = FALSE
) %>%
  tidyr::separate(Label, into = c("Project", "Patient", "Tissue"), sep = "_") %>%
  filter(Tissue %in% c("T", "N"))

common_lung <- intersect(gsm_map$GSM, colnames(counts_lung))

counts_lung <- counts_lung[, common_lung, drop = FALSE]

meta_lung <- gsm_map %>%
  filter(GSM %in% common_lung) %>%
  column_to_rownames("GSM")

dds_lung <- DESeqDataSetFromMatrix(
  countData = counts_lung,
  colData = meta_lung,
  design = ~ Tissue
)

vsd_lung <- vst(dds_lung, blind = TRUE)
mat_lung <- assay(vsd_lung)

sym_lung <- mapIds(
  org.Hs.eg.db,
  keys = rownames(mat_lung),
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

rownames(mat_lung) <- ifelse(
  is.na(sym_lung) | sym_lung == "",
  rownames(mat_lung),
  sym_lung
)

lung_expr_full <- as.data.frame(t(mat_lung))
lung_truth <- factor(
  ifelse(meta_lung$Tissue == "T", "Tumor", "Normal"),
  levels = c("Normal", "Tumor")
)

# ------------------------------------------------------------------
# External dataset 3: CPTAC kidney cohort
# ------------------------------------------------------------------

message("Loading CPTAC kidney cohort...")

count_files <- list.files(
  kidney_counts_dir,
  pattern = "\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(count_files) == 0) {
  stop("No kidney STAR-count files found in: ", kidney_counts_dir)
}

expr_list <- list()

for (f in count_files) {
  dat <- readr::read_tsv(f, skip = 1, show_col_types = FALSE) %>%
    dplyr::select(gene_id, gene_name, count = unstranded) %>%
    mutate(
      gene_id = sub("\\..*", "", gene_id),
      count = as.numeric(count)
    ) %>%
    filter(!is.na(gene_name)) %>%
    group_by(gene_name) %>%
    summarise(count = sum(count), .groups = "drop")

  sid <- basename(dirname(f))

  expr_list[[sid]] <- dat %>%
    setNames(c("gene_name", sid))
}

expr_merged_kidney <- Reduce(
  function(a, b) full_join(a, b, by = "gene_name"),
  expr_list
) %>%
  column_to_rownames("gene_name")

query_expr <- GDCquery(
  project = "CPTAC-3",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

metadata_expr <- getResults(
  query_expr,
  cols = c("file_id", "cases.submitter_id", "sample_type")
)

metadata_kidney <- metadata_expr %>%
  filter(
    file_id %in% colnames(expr_merged_kidney),
    sample_type %in% c("Primary Tumor", "Solid Tissue Normal")
  ) %>%
  dplyr::select(sample = file_id, condition = sample_type) %>%
  mutate(
    condition = factor(
      ifelse(grepl("Normal", condition), "Normal", "Tumor"),
      levels = c("Normal", "Tumor")
    )
  ) %>%
  distinct() %>%
  column_to_rownames("sample")

dds_kidney <- DESeqDataSetFromMatrix(
  expr_merged_kidney[, rownames(metadata_kidney), drop = FALSE],
  metadata_kidney,
  design = ~ condition
)

vsd_kidney <- vst(dds_kidney, blind = TRUE)
mat_kidney <- assay(vsd_kidney)

kidney_expr_full <- as.data.frame(t(mat_kidney))
kidney_truth <- metadata_kidney$condition

# ------------------------------------------------------------------
# External dataset 4: Korean CRC cohort
# ------------------------------------------------------------------

message("Loading Korean CRC cohort...")

raw_crc <- read.delim(
  crc_counts_file,
  check.names = FALSE
)

gene_symbols_crc <- raw_crc[[1]]

counts_crc <- as.data.frame(raw_crc[, -1, drop = FALSE])
rownames(counts_crc) <- gene_symbols_crc

counts_crc[] <- lapply(
  counts_crc,
  function(x) as.numeric(as.character(x))
)

counts_sym_crc <- as.data.frame(
  rowsum(
    as.matrix(counts_crc),
    group = rownames(counts_crc),
    reorder = FALSE,
    na.rm = TRUE
  )
)

sids_crc <- colnames(counts_sym_crc)

tissue_crc <- ifelse(
  grepl("-T-", sids_crc),
  "T",
  ifelse(grepl("-N-", sids_crc), "N", NA)
)

sample_info_crc <- data.frame(
  GSM = sids_crc,
  Tissue = tissue_crc,
  row.names = sids_crc,
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(Tissue))

counts_sym_crc <- counts_sym_crc[, rownames(sample_info_crc), drop = FALSE]

dds_crc <- DESeqDataSetFromMatrix(
  countData = round(counts_sym_crc),
  colData = sample_info_crc,
  design = ~ Tissue
)

vsd_crc <- vst(dds_crc, blind = TRUE)
mat_crc <- assay(vsd_crc)

crc_expr_full <- as.data.frame(t(mat_crc))
crc_truth <- factor(
  ifelse(sample_info_crc$Tissue == "T", "Tumor", "Normal"),
  levels = c("Normal", "Tumor")
)

# ------------------------------------------------------------------
# Assemble external sets
# ------------------------------------------------------------------

ext_sets <- list(
  BRCA = list(expr = brca_expr_full, truth = brca_truth),
  Lung = list(expr = lung_expr_full, truth = lung_truth),
  Kidney = list(expr = kidney_expr_full, truth = kidney_truth),
  CRC_KR = list(expr = crc_expr_full, truth = crc_truth)
)

missing_by_set <- purrr::map(
  ext_sets,
  ~ setdiff(panel_eval_genes, colnames(.x$expr))
)

missing_any <- purrr::keep(missing_by_set, ~ length(.x) > 0)

if (length(missing_any) > 0) {
  print(missing_any)
  stop("At least one external dataset is missing required panel genes.")
}

# ------------------------------------------------------------------
# Baseline external AUC with full model
# ------------------------------------------------------------------

baseline_perf <- map_dfr(names(ext_sets), function(nm) {
  auc_full <- eval_on_set(rf_full, panel_eval_genes, ext_sets[[nm]])

  tibble(
    Dataset = nm,
    Gene = "All 9 genes",
    AUC = auc_full
  )
})

# ------------------------------------------------------------------
# Drop-one-gene external importance
# ------------------------------------------------------------------

message("Running drop-one-gene external evaluation...")

drop_perf <- map_dfr(panel_eval_genes, function(gene_drop) {
  genes_left <- setdiff(panel_eval_genes, gene_drop)

  rf_drop <- fit_rf_panel(
    train_df = tcga_clean,
    labels = y_tcga,
    genes = genes_left,
    ctrl = ctrl_pan,
    tune_len = tune_len
  )

  map_dfr(names(ext_sets), function(nm) {
    auc_drop <- eval_on_set(rf_drop, genes_left, ext_sets[[nm]])

    tibble(
      Dataset = nm,
      Gene = gene_drop,
      AUC = auc_drop
    )
  })
})

all_perf <- bind_rows(baseline_perf, drop_perf)

write.csv(
  all_perf,
  file.path(out_dir, "external_drop_one_auc_all_results.csv"),
  row.names = FALSE
)

delta_auc <- all_perf %>%
  group_by(Dataset) %>%
  mutate(
    auc_full = AUC[Gene == "All 9 genes"],
    delta_auc = AUC - auc_full
  ) %>%
  ungroup() %>%
  filter(Gene != "All 9 genes") %>%
  mutate(
    importance = -delta_auc,
    GeneLab = unname(to_ccp(Gene))
  )

write.csv(
  delta_auc,
  file.path(out_dir, "external_drop_one_importance_by_dataset.csv"),
  row.names = FALSE
)

avg_importance <- delta_auc %>%
  group_by(GeneLab) %>%
  summarise(
    mean_importance = mean(importance, na.rm = TRUE),
    sd_importance = sd(importance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_importance))

write.csv(
  avg_importance,
  file.path(out_dir, "external_drop_one_importance_average.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Plot average external importance
# ------------------------------------------------------------------

p_avg <- ggplot(
  avg_importance %>%
    mutate(GeneLab = factor(GeneLab, levels = rev(GeneLab))),
  aes(x = GeneLab, y = mean_importance)
) +
  geom_col(fill = "plum", width = 0.75) +
  geom_errorbar(
    aes(
      ymin = mean_importance - sd_importance,
      ymax = mean_importance + sd_importance
    ),
    width = 0.2
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(
    title = "Average external drop-one-gene importance",
    subtitle = "Importance = AUC(full model) - AUC(drop-one model)",
    x = "Dropped gene",
    y = "Mean importance across external datasets"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  )

print(p_avg)

ggsave(
  file.path(out_dir, "external_importance_average.pdf"),
  p_avg,
  width = 6.5,
  height = 4.8
)

# ------------------------------------------------------------------
# Plot per-dataset external importance
# ------------------------------------------------------------------

global_order <- avg_importance %>%
  arrange(desc(mean_importance)) %>%
  pull(GeneLab)

p_dataset <- ggplot(
  delta_auc %>%
    mutate(GeneLab = factor(GeneLab, levels = rev(global_order))),
  aes(x = GeneLab, y = importance, fill = Dataset)
) +
  geom_col(width = 0.75, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  facet_wrap(~ Dataset, ncol = 2) +
  scale_fill_manual(values = c(
    "BRCA" = "steelblue",
    "CRC_KR" = "grey60",
    "Lung" = "forestgreen",
    "Kidney" = "purple"
  )) +
  labs(
    title = "External drop-one-gene importance by dataset",
    subtitle = "Genes ordered by mean importance across external datasets",
    x = "Dropped gene",
    y = "Importance: AUC(full) - AUC(drop)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  )

print(p_dataset)

ggsave(
  file.path(out_dir, "external_importance_by_dataset.pdf"),
  p_dataset,
  width = 9,
  height = 6.5
)

message("\nFinished. Outputs written to:\n", out_dir)
