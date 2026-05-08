# 15-gene signature PCA in breast patient cohort
#
#
# Input:
# - Salmon quant.sf files for all samples
# - Gencode GTF annotation file
#
# Notes:
# - Edit 'annotation_dir' and 'salmon_dir' below before running
# - CCP genes are displayed internally as AGTPBP1 / AGBL1-5
# - The paired PT-NAT PCA uses tumour-minus-matched-normal expression on the VST scale

suppressPackageStartupMessages({
  library(DESeq2)
  library(tximport)
  library(tidyverse)
  library(GenomicFeatures)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(patchwork)
})

# ------------------------------------------------------------------
# Input paths
# ------------------------------------------------------------------

annotation_dir <- "path-to-annotation-folder"
salmon_dir     <- "path-to-salmon-folder"

gtf_file <- file.path(annotation_dir, "annotations", "gencode.v29.annotation.gtf")

if (!dir.exists(annotation_dir)) {
  stop("Please set 'annotation_dir' to the folder containing the GTF annotation.")
}

if (!dir.exists(salmon_dir)) {
  stop("Please set 'salmon_dir' to the folder containing the Salmon quantification folders.")
}

if (!file.exists(gtf_file)) {
  stop("GTF file not found: ", gtf_file)
}

# ------------------------------------------------------------------
# Sample metadata
# ------------------------------------------------------------------

samples_patient <- tibble::tribble(
  ~sample,   ~tissue, ~group,         ~path,
  "1NAT",    "NAT",   "non-invasive", "salmon_1NAT/quant.sf",
  "2PT",     "PT",    "non-invasive", "salmon_2PT/quant.sf",
  "3NAT",    "NAT",   "non-invasive", "salmon_3NAT/quant.sf",
  "4PT",     "PT",    "non-invasive", "salmon_4PT/quant.sf",
  "5NAT",    "NAT",   "metastatic",   "salmon_5NAT/quant.sf",
  "6PT",     "PT",    "metastatic",   "salmon_6PT/quant.sf",
  "7NAT",    "NAT",   "metastatic",   "salmon_7NAT/quant.sf",
  "8PT",     "PT",    "metastatic",   "salmon_8PT/quant.sf",
  "9NAT",    "NAT",   "metastatic",   "salmon_9NAT/quant.sf",
  "10PT",    "PT",    "metastatic",   "salmon_10PT/quant.sf",
  "11NAT",   "NAT",   "TNBC",         "salmon_11NAT/quant.sf",
  "12PT",    "PT",    "TNBC",         "salmon_12PT/quant.sf",
  "13NAT",   "NAT",   "TNBC",         "salmon_13NAT/quant.sf",
  "14PT",    "PT",    "TNBC",         "salmon_14PT/quant.sf",
  "15NAT",   "NAT",   "TNBC",         "salmon_15NAT/quant.sf",
  "16PT",    "PT",    "TNBC",         "salmon_16PT/quant.sf",
  "17NAT",   "NAT",   "non-invasive", "salmon_17NAT/quant.sf",
  "18PT",    "PT",    "non-invasive", "salmon_18PT/quant.sf"
)

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
# Build tx2gene from GTF
# ------------------------------------------------------------------

message("Building tx2gene from GTF...")
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")

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

# ------------------------------------------------------------------
# Import Salmon quantifications
# ------------------------------------------------------------------

files <- file.path(salmon_dir, samples_patient$path)
names(files) <- samples_patient$sample

missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0) {
  stop("The following Salmon files were not found:\n", paste(missing_files, collapse = "\n"))
}

message("Importing Salmon quantifications...")
txi_patient <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE
)

sample_table <- samples_patient %>%
  column_to_rownames("sample")

dds_patient <- DESeqDataSetFromTximport(
  txi = txi_patient,
  colData = sample_table,
  design = ~ tissue
)

vsd_patient <- vst(dds_patient, blind = TRUE)
vst_mat <- assay(vsd_patient)
colnames(vst_mat) <- rownames(colData(dds_patient))

# ------------------------------------------------------------------
# Map Ensembl IDs to HGNC symbols
# ------------------------------------------------------------------

aggregate_dup_symbols <- function(mat, symbols, mode = "mean") {
  ok <- !is.na(symbols) & symbols != ""
  mat <- mat[ok, , drop = FALSE]
  symbols <- symbols[ok]

  if (mode == "sum") {
    out <- rowsum(mat, group = symbols)
  } else {
    sp <- split(seq_len(nrow(mat)), symbols)
    out <- vapply(
      sp,
      function(ix) colMeans(mat[ix, , drop = FALSE]),
      numeric(ncol(mat))
    )
    out <- t(out)
    rownames(out) <- names(sp)
  }

  out
}

ens_ids <- sub("\\..*", "", rownames(vst_mat))

gene_symbols <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = ens_ids,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

vst_sym <- aggregate_dup_symbols(vst_mat, gene_symbols, mode = "mean")

have <- intersect(feature_genes, rownames(vst_sym))
missing <- setdiff(feature_genes, have)

cat("Present genes:\n")
print(have)
cat("Missing genes:\n")
print(missing)

if (length(have) < 10) {
  stop("Too few of the 15 signature genes were found after annotation.")
}

# ------------------------------------------------------------------
# Build expression matrix for the 15-gene signature
# ------------------------------------------------------------------

patient_expr_df <- t(vst_sym[have, , drop = FALSE]) %>% as.data.frame()

for (g in missing) {
  patient_expr_df[[g]] <- NA_real_
}

patient_expr_df <- patient_expr_df[, feature_genes, drop = FALSE]

for (g in colnames(patient_expr_df)) {
  if (anyNA(patient_expr_df[[g]])) {
    med_g <- median(patient_expr_df[[g]], na.rm = TRUE)
    patient_expr_df[[g]][is.na(patient_expr_df[[g]])] <- if (is.finite(med_g)) med_g else 0
  }
}

# ------------------------------------------------------------------
# Metadata formatting
# ------------------------------------------------------------------

patient_meta <- samples_patient %>%
  mutate(
    tissue = factor(tissue, levels = c("NAT", "PT")),
    group  = factor(group, levels = c("non-invasive", "metastatic", "TNBC"))
  )

patient_expr_df <- patient_expr_df[patient_meta$sample, , drop = FALSE]

# ------------------------------------------------------------------
# PCA helper
# ------------------------------------------------------------------

make_pca_plot <- function(expr_df, meta_df, color_var, title_txt, color_values) {
  sample_ids <- intersect(rownames(expr_df), meta_df$sample)

  expr_df <- expr_df[sample_ids, , drop = FALSE]
  meta_df <- meta_df %>%
    filter(sample %in% sample_ids) %>%
    slice(match(sample_ids, sample))

  keep <- apply(expr_df, 2, var, na.rm = TRUE) > 0
  expr_df <- expr_df[, keep, drop = FALSE]

  if (nrow(expr_df) < 2 || ncol(expr_df) < 2) {
    stop("Not enough samples or variable genes for PCA.")
  }

  pca_res <- prcomp(expr_df, center = TRUE, scale. = TRUE)
  var_exp <- 100 * (pca_res$sdev^2 / sum(pca_res$sdev^2))

  pca_df <- as.data.frame(pca_res$x[, 1:2, drop = FALSE]) %>%
    rownames_to_column("sample") %>%
    left_join(meta_df, by = "sample")

  ggplot(pca_df, aes(x = PC1, y = PC2, color = .data[[color_var]])) +
    geom_point(size = 4.5, shape = 16, alpha = 1) +
    scale_color_manual(values = color_values) +
    labs(
      title = title_txt,
      x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
      y = paste0("PC2 (", round(var_exp[2], 1), "%)")
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "bottom"
    )
}

tissue_colors <- c(
  "NAT" = "lightskyblue",
  "PT"  = "#FF6A6A"
)

group_colors <- c(
  "non-invasive" = "#87CEFA",
  "metastatic"   = "#FFB6C1",
  "TNBC"         = "#8B1A1A"
)

# ------------------------------------------------------------------
# PCA 1: all samples coloured by tissue
# ------------------------------------------------------------------

p_pca_all <- make_pca_plot(
  expr_df = patient_expr_df,
  meta_df = patient_meta,
  color_var = "tissue",
  title_txt = "PCA of Breast Patient Samples (15-Gene Signature)",
  color_values = tissue_colors
)

print(p_pca_all)

# ------------------------------------------------------------------
# PCA 2: each subgroup separately, coloured by tissue
# ------------------------------------------------------------------

plot_pca_by_group <- function(expr_df, meta_df, subgroup_name) {
  meta_sub <- meta_df %>% filter(group == subgroup_name)
  expr_sub <- expr_df[meta_sub$sample, , drop = FALSE]

  make_pca_plot(
    expr_df = expr_sub,
    meta_df = meta_sub,
    color_var = "tissue",
    title_txt = subgroup_name,
    color_values = tissue_colors
  )
}

p_noninv <- plot_pca_by_group(patient_expr_df, patient_meta, "non-invasive")
p_met    <- plot_pca_by_group(patient_expr_df, patient_meta, "metastatic")
p_tnbc   <- plot_pca_by_group(patient_expr_df, patient_meta, "TNBC")

combined_group_pca <- (p_noninv + p_met + p_tnbc) +
  plot_annotation(
    title = "PCA by Breast Cancer Subgroup",
    subtitle = "Each PCA computed separately within subgroup using the 15-gene signature"
  )

print(combined_group_pca)

# ------------------------------------------------------------------
# PCA 3: paired PT-NAT differences per patient
# ------------------------------------------------------------------

patient_meta2 <- patient_meta %>%
  mutate(
    sample_num = as.numeric(stringr::str_extract(sample, "^[0-9]+")),
    patient_id = ceiling(sample_num / 2)
  )

patient_ids <- sort(unique(patient_meta2$patient_id))

logPTNAT_mat <- matrix(
  NA_real_,
  nrow = length(feature_genes),
  ncol = length(patient_ids),
  dimnames = list(feature_genes, paste0("P", patient_ids))
)

patient_group <- patient_meta2 %>%
  distinct(patient_id, group) %>%
  arrange(patient_id)

for (pid in patient_ids) {
  meta_pid <- patient_meta2 %>% filter(patient_id == pid)

  pt_sample  <- meta_pid %>% filter(tissue == "PT")  %>% pull(sample)
  nat_sample <- meta_pid %>% filter(tissue == "NAT") %>% pull(sample)

  if (
    length(pt_sample) == 1 &&
    length(nat_sample) == 1 &&
    pt_sample %in% rownames(patient_expr_df) &&
    nat_sample %in% rownames(patient_expr_df)
  ) {
    logPTNAT_mat[, paste0("P", pid)] <-
      as.numeric(patient_expr_df[pt_sample, feature_genes] -
                   patient_expr_df[nat_sample, feature_genes])
  } else {
    message("Missing PT/NAT pair for patient ", pid)
  }
}

ratio_mat <- t(logPTNAT_mat) %>% as.data.frame()
ratio_mat$Patient <- rownames(ratio_mat)

pca_ratio_meta <- patient_group %>%
  mutate(Patient = paste0("P", patient_id)) %>%
  select(Patient, group)

ratio_expr <- ratio_mat %>%
  select(all_of(feature_genes)) %>%
  as.data.frame()

rownames(ratio_expr) <- ratio_mat$Patient

keep <- apply(ratio_expr, 2, var, na.rm = TRUE) > 0
ratio_expr <- ratio_expr[, keep, drop = FALSE]

pca_ratio_res <- prcomp(ratio_expr, center = TRUE, scale. = TRUE)
var_exp_ratio <- 100 * (pca_ratio_res$sdev^2 / sum(pca_ratio_res$sdev^2))

pca_ratio_df <- as.data.frame(pca_ratio_res$x[, 1:2, drop = FALSE]) %>%
  rownames_to_column("Patient") %>%
  left_join(pca_ratio_meta, by = "Patient")

p_pca_ratio <- ggplot(pca_ratio_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 4.5, shape = 16, alpha = 1) +
  scale_color_manual(values = group_colors) +
  labs(
    title = "PCA of Paired PT-NAT Differences (15-Gene Signature)",
    subtitle = "Each point represents tumour-induced change relative to matched normal",
    x = paste0("PC1 (", round(var_exp_ratio[1], 1), "%)"),
    y = paste0("PC2 (", round(var_exp_ratio[2], 1), "%)"),
    color = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )

print(p_pca_ratio)
