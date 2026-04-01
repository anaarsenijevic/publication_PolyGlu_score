# BRCA patient-group heatmap of TTLL and CCP family dysregulation
#
#
# Inputs:
# - Salmon quant.sf files listed in samples_patient
# - Gencode v29 GTF file
#
# Notes:
# - Edit the paths below before running
# - NAT = normal adjacent tissue
# - PT  = primary tumour

suppressPackageStartupMessages({
  library(DESeq2)
  library(tximport)
  library(tidyverse)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(pheatmap)
  library(GenomicFeatures)
  library(biomaRt)
})

# ------------------------------------------------------------------
# Input paths
# ------------------------------------------------------------------

patient_dir <- "path-to-your-salmon-folder"
gtf_file <- "path-to-your-gencode-v29.gtf"

if (!dir.exists(patient_dir)) {
  stop("Please set 'patient_dir' to the folder containing the Salmon sample folders.")
}

if (!file.exists(gtf_file)) {
  stop("Please set 'gtf_file' to the GTF annotation file.")
}

# ------------------------------------------------------------------
# Sample metadata
# ------------------------------------------------------------------

samples_patient <- tibble::tribble(
  ~sample, ~tissue, ~group,         ~path,
  "1NAT",  "NAT",   "non-invasive", "salmon_1NAT/quant.sf",
  "2PT",   "PT",    "non-invasive", "salmon_2PT/quant.sf",
  "3NAT",  "NAT",   "non-invasive", "salmon_3NAT/quant.sf",
  "4PT",   "PT",    "non-invasive", "salmon_4PT/quant.sf",
  "5NAT",  "NAT",   "metastatic",   "salmon_5NAT/quant.sf",
  "6PT",   "PT",    "metastatic",   "salmon_6PT/quant.sf",
  "7NAT",  "NAT",   "metastatic",   "salmon_7NAT/quant.sf",
  "8PT",   "PT",    "metastatic",   "salmon_8PT/quant.sf",
  "9NAT",  "NAT",   "metastatic",   "salmon_9NAT/quant.sf",
  "10PT",  "PT",    "metastatic",   "salmon_10PT/quant.sf",
  "11NAT", "NAT",   "TNBC",         "salmon_11NAT/quant.sf",
  "12PT",  "PT",    "TNBC",         "salmon_12PT/quant.sf",
  "13NAT", "NAT",   "TNBC",         "salmon_13NAT/quant.sf",
  "14PT",  "PT",    "TNBC",         "salmon_14PT/quant.sf",
  "15NAT", "NAT",   "TNBC",         "salmon_15NAT/quant.sf",
  "16PT",  "PT",    "TNBC",         "salmon_16PT/quant.sf",
  "17NAT", "NAT",   "non-invasive", "salmon_17NAT/quant.sf",
  "18PT",  "PT",    "non-invasive", "salmon_18PT/quant.sf"
) %>%
  mutate(path = file.path(patient_dir, path))

missing_files <- samples_patient$path[!file.exists(samples_patient$path)]
if (length(missing_files) > 0) {
  stop(
    "The following Salmon quantification files were not found:\n",
    paste(missing_files, collapse = "\n")
  )
}

# ------------------------------------------------------------------
# Gene panel and display labels
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
# Helper functions
# ------------------------------------------------------------------

clean_ensembl_ids <- function(ids) {
  sub("\\..*", "", ids)
}

aggregate_dup_symbols <- function(mat, symbols, mode = c("mean", "sum")[1]) {
  ok <- !is.na(symbols) & symbols != ""
  mat <- mat[ok, , drop = FALSE]
  symbols <- symbols[ok]

  if (mode == "sum") {
    rowsum(mat, group = symbols)
  } else {
    sp <- split(seq_len(nrow(mat)), symbols)
    out <- vapply(
      sp,
      function(ix) colMeans(mat[ix, , drop = FALSE]),
      numeric(ncol(mat))
    )
    out <- t(out)
    rownames(out) <- names(sp)
    out
  }
}

# ------------------------------------------------------------------
# Transcript-to-gene mapping from GTF
# ------------------------------------------------------------------

txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")

tx2gene <- AnnotationDbi::select(
  txdb,
  keys = keys(txdb, "TXNAME"),
  columns = "GENEID",
  keytype = "TXNAME"
) %>%
  mutate(
    TXNAME = clean_ensembl_ids(TXNAME),
    GENEID = clean_ensembl_ids(GENEID)
  ) %>%
  distinct()

# ------------------------------------------------------------------
# Import Salmon quantifications
# ------------------------------------------------------------------

files <- samples_patient$path
names(files) <- samples_patient$sample

txi_patient <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE
)

sample_table <- samples_patient %>%
  select(sample, tissue, group) %>%
  column_to_rownames("sample")

dds_patient <- DESeqDataSetFromTximport(
  txi = txi_patient,
  colData = sample_table,
  design = ~ tissue
)

# ------------------------------------------------------------------
# VST expression and gene-symbol mapping
# ------------------------------------------------------------------

vsd_patient <- vst(dds_patient, blind = TRUE)
vst_expr <- assay(vsd_patient)

ens_ids <- clean_ensembl_ids(rownames(vst_expr))

mart <- useEnsembl(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl"
)

gene_map_pat <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = unique(ens_ids),
  mart = mart
)

sym <- gene_map_pat$external_gene_name[
  match(ens_ids, gene_map_pat$ensembl_gene_id)
]

vst_sym <- aggregate_dup_symbols(vst_expr, sym, mode = "mean")

have <- intersect(feature_genes, rownames(vst_sym))
missing <- setdiff(feature_genes, have)

if (length(have) < 7) {
  stop(
    "Too few panel genes in patient data. Present: ",
    paste(have, collapse = ", "),
    " | Missing: ",
    paste(missing, collapse = ", ")
  )
}

patient_expr_df <- t(vst_sym[have, , drop = FALSE]) %>% as.data.frame()

for (g in missing) {
  patient_expr_df[[g]] <- NA_real_
}

patient_expr_df <- patient_expr_df[, feature_genes, drop = FALSE]

# ------------------------------------------------------------------
# Gene mapping for DESeq2 results
# ------------------------------------------------------------------

rownames(txi_patient$counts) <- clean_ensembl_ids(rownames(txi_patient$counts))
ens_ids_clean <- clean_ensembl_ids(rownames(txi_patient$counts))

gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ens_ids_clean,
  mart = mart
) %>%
  filter(hgnc_symbol != "") %>%
  distinct()

panel_ens_ids <- gene_map$ensembl_gene_id

# ------------------------------------------------------------------
# Assign paired patient IDs
# ------------------------------------------------------------------

samples_patient <- samples_patient %>%
  mutate(patient_id = ceiling(as.numeric(gsub("[A-Za-z]", "", sample)) / 2))

# ------------------------------------------------------------------
# DESeq2 within each patient group
# ------------------------------------------------------------------

groups <- unique(samples_patient$group)
log2fc_list <- list()
padj_list <- list()

for (grp in groups) {
  message("Running DESeq2 for group: ", grp)

  sub_meta <- samples_patient %>% filter(group == grp)
  sub_samples <- sub_meta$sample

  txi_sub <- lapply(txi_patient[1:3], function(x) {
    x[, match(sub_samples, colnames(x)), drop = FALSE]
  })
  txi_sub$countsFromAbundance <- txi_patient$countsFromAbundance

  sub_meta <- sub_meta %>%
    column_to_rownames("sample")

  dds <- DESeqDataSetFromTximport(
    txi_sub,
    colData = sub_meta,
    design = ~ patient_id + tissue
  )

  dds$tissue <- relevel(dds$tissue, ref = "NAT")

  keep <- rowSums(counts(dds)) > 10 | rownames(dds) %in% panel_ens_ids
  dds <- dds[keep, ]

  dds <- DESeq(dds)
  res <- results(dds, contrast = c("tissue", "PT", "NAT"))

  res_df <- as.data.frame(res) %>%
    mutate(ensembl_gene_id = rownames(.)) %>%
    left_join(gene_map, by = "ensembl_gene_id") %>%
    filter(hgnc_symbol %in% feature_genes) %>%
    select(hgnc_symbol, log2FoldChange, padj)

  log2fc_list[[grp]] <- res_df %>%
    select(hgnc_symbol, log2FoldChange) %>%
    rename_with(~ grp, .cols = log2FoldChange)

  padj_list[[grp]] <- res_df %>%
    select(hgnc_symbol, padj) %>%
    rename_with(~ grp, .cols = padj)
}

# ------------------------------------------------------------------
# Combine log2 fold changes and adjusted p-values
# ------------------------------------------------------------------

combined_fc <- reduce(log2fc_list, full_join, by = "hgnc_symbol") %>%
  column_to_rownames("hgnc_symbol")

combined_padj <- reduce(padj_list, full_join, by = "hgnc_symbol") %>%
  column_to_rownames("hgnc_symbol")

desired_order <- c("non-invasive", "metastatic", "TNBC")

combined_fc_ordered <- combined_fc[feature_genes, desired_order, drop = FALSE]
combined_padj_ordered <- combined_padj[feature_genes, desired_order, drop = FALSE]

rownames(combined_fc_ordered) <- symbol_labels[feature_genes]
rownames(combined_padj_ordered) <- symbol_labels[feature_genes]

# ------------------------------------------------------------------
# Heatmap settings
# ------------------------------------------------------------------

max_abs_fc <- max(abs(combined_fc_ordered), na.rm = TRUE)
color_breaks <- seq(-max_abs_fc, max_abs_fc, length.out = 101)
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# ------------------------------------------------------------------
# Full heatmap
# ------------------------------------------------------------------

pheatmap(
  combined_fc_ordered,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = color_palette,
  breaks = color_breaks,
  na_col = "gray90",
  main = "Log2 fold change per BRCA patient group — TTLL/CCPs"
)

# ------------------------------------------------------------------
# Significant-only heatmap
# ------------------------------------------------------------------

combined_fc_sig <- combined_fc_ordered
combined_fc_sig[is.na(combined_padj_ordered) | combined_padj_ordered >= 0.05] <- NA

pheatmap(
  combined_fc_sig,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = color_palette,
  breaks = color_breaks,
  na_col = "gray90",
  main = "Log2 fold change per BRCA patient group — TTLL/CCPs (FDR < 0.05)"
)
