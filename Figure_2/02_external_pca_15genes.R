# External validation PCA using the 15-gene panel
#
# This script:
# - Loads four external RNA-seq datasets: BRCA-9, Lung, Kidney, and CRC
# - Performs VST normalization where required
# - Maps gene IDs to gene symbols
# - Extracts expression of 15 TTLL/CCP polyglutamylation-related genes
# - Performs PCA independently for each dataset
# - Visualizes tumor-normal separation for each dataset
#
# Inputs:
# - BRCA-9: Salmon quant.sf files + GTF annotation
# - Lung: GSE87410 raw counts + metadata
# - Kidney: CPTAC-3 RNA-seq counts
# - CRC: Korean CRC count matrix
#
# Output:
# - 15-gene PCA plot for all four datasets (.pdf, .png)
#
# Notes:
# - Samples are colored by Tumor vs Normal
# - PCA is performed independently per dataset
# - Gene mapping uses org.Hs.eg.db, avoiding BioMart dependency

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(tximport)
  library(GenomicFeatures)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(TCGAbiolinks)
  library(patchwork)
})

set.seed(123)

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------

brca_annotation_dir <- "path-to/Seq_analysis"
brca_salmon_dir     <- "path-to/Patient_samples/salmon"

lung_dir   <- "path-to/GSE87410"
kidney_dir <- "path-to/KIRC_CPTAC"
crc_dir    <- "path-to/KoreanCRC"
# ------------------------------------------------------------------
# 15-gene panel
# ------------------------------------------------------------------

panel_15 <- c(
  "TTLL1", "TTLL2", "TTLL4", "TTLL5", "TTLL6", "TTLL7", "TTLL9", "TTLL11", "TTLL13",
  "AGBL1", "AGBL2", "AGBL3", "AGBL4", "AGBL5", "AGTPBP1"
)

condition_colors <- c(
  "Tumor" = "#FF6A6A",
  "Normal" = "lightskyblue"
)

# ------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------

collapse_duplicate_symbols <- function(expr_mat) {
  expr_mat <- expr_mat[
    !is.na(rownames(expr_mat)) & rownames(expr_mat) != "",
    ,
    drop = FALSE
  ]

  as.data.frame(
    rowsum(
      as.matrix(expr_mat),
      group = rownames(expr_mat),
      reorder = FALSE
    )
  )
}

map_ensembl_to_symbol <- function(ensembl_ids) {
  AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = sub("\\..*", "", ensembl_ids),
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
}

map_entrez_to_symbol <- function(entrez_ids) {
  AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = as.character(entrez_ids),
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
  )
}

make_pca_plot <- function(expr_mat, metadata, genes, dataset_name) {

  genes_present <- intersect(genes, rownames(expr_mat))
  missing_genes <- setdiff(genes, genes_present)

  message(dataset_name, " | genes present: ",
          length(genes_present), "/", length(genes))

  if (length(missing_genes) > 0) {
    message(dataset_name, " | missing genes: ",
            paste(missing_genes, collapse = ", "))
  }

  if (length(genes_present) < 2) {
    stop("Fewer than two panel genes found for ", dataset_name)
  }

  common_samples <- intersect(colnames(expr_mat), rownames(metadata))

  expr_mat <- expr_mat[, common_samples, drop = FALSE]
  metadata <- metadata[common_samples, , drop = FALSE]

  pca_input <- t(expr_mat[genes_present, , drop = FALSE])
  pca_input <- as.matrix(pca_input)

  keep <- complete.cases(pca_input)
  pca_input <- pca_input[keep, , drop = FALSE]
  metadata <- metadata[keep, , drop = FALSE]

  pca_result <- prcomp(pca_input, center = TRUE, scale. = TRUE)
  variance_explained <- summary(pca_result)$importance[2, 1:2] * 100

  pca_df <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    condition = metadata$condition
  )

  ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 2.5, alpha = 0.9) +
    scale_color_manual(values = condition_colors) +
    labs(
      title = dataset_name,
      subtitle = "15-gene panel",
      x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
      y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
      color = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom",
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    )
}

# ============================================================
# 1. BRCA-9 dataset
# ============================================================

samples_patient <- tibble::tribble(
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

setwd(brca_annotation_dir)

txdb <- makeTxDbFromGFF(
  "annotations/gencode.v29.annotation.gtf",
  format = "gtf"
)

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

setwd(brca_salmon_dir)

files <- samples_patient$path
names(files) <- samples_patient$sample

txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE
)

sample_table <- samples_patient %>%
  dplyr::select(sample, tissue, group) %>%
  mutate(condition = ifelse(tissue == "PT", "Tumor", "Normal")) %>%
  column_to_rownames("sample")

dds_brca <- DESeqDataSetFromTximport(
  txi,
  colData = sample_table,
  design = ~ tissue
)

vsd_brca <- vst(dds_brca, blind = TRUE)
expr_brca <- assay(vsd_brca)

brca_symbols <- map_ensembl_to_symbol(rownames(expr_brca))
rownames(expr_brca) <- brca_symbols
expr_brca <- collapse_duplicate_symbols(expr_brca)

meta_brca <- sample_table %>%
  dplyr::select(condition)

# ============================================================
# 2. Lung dataset: GSE87410
# ============================================================

setwd(lung_dir)

raw_lung <- readr::read_tsv(
  "GSE87410_raw_counts_GRCh38.p13_NCBI.tsv.gz",
  show_col_types = FALSE
)

counts_lung <- as.data.frame(raw_lung[, -1])
rownames(counts_lung) <- raw_lung[[1]]
counts_lung[] <- lapply(counts_lung, as.numeric)

gsm_map <- readr::read_csv(
  "GSM_sample_map.csv",
  col_names = c("GSM", "Label"),
  show_col_types = FALSE
) %>%
  tidyr::separate(
    Label,
    into = c("Project", "Patient", "Tissue"),
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(Tissue = ifelse(Tissue %in% c("T", "N"), Tissue, NA)) %>%
  filter(!is.na(Tissue))

valid_lung_samples <- intersect(colnames(counts_lung), gsm_map$GSM)
counts_lung <- counts_lung[, valid_lung_samples, drop = FALSE]

meta_lung_all <- gsm_map %>%
  filter(GSM %in% valid_lung_samples) %>%
  dplyr::select(GSM, Project, Tissue) %>%
  column_to_rownames("GSM")

dds_lung <- DESeqDataSetFromMatrix(
  countData = round(counts_lung),
  colData = meta_lung_all,
  design = ~ Project + Tissue
)

vsd_lung <- vst(dds_lung, blind = TRUE)
expr_lung <- assay(vsd_lung)

lung_symbols <- map_entrez_to_symbol(rownames(expr_lung))
rownames(expr_lung) <- ifelse(
  is.na(lung_symbols) | lung_symbols == "",
  rownames(expr_lung),
  lung_symbols
)

expr_lung <- collapse_duplicate_symbols(expr_lung)

keep_lung <- rownames(meta_lung_all)[meta_lung_all$Project %in% c("LUAD", "LUSC")]

expr_lung <- expr_lung[, keep_lung, drop = FALSE]

meta_lung <- meta_lung_all[keep_lung, , drop = FALSE] %>%
  mutate(condition = ifelse(Tissue == "T", "Tumor", "Normal")) %>%
  dplyr::select(condition)

# ============================================================
# 3. Kidney dataset: CPTAC-3
# ============================================================

setwd(kidney_dir)

count_files <- list.files(
  "GDC_kidney_filtered",
  pattern = "\\.rna_seq\\.augmented_star_gene_counts\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

expr_list <- list()

for (file in count_files) {

  dat <- tryCatch({
    readr::read_tsv(file, skip = 1, show_col_types = FALSE) %>%
      dplyr::select(gene_id, gene_name, count = unstranded) %>%
      mutate(
        gene_id = sub("\\..*", "", gene_id),
        count = as.numeric(count)
      ) %>%
      filter(!is.na(gene_name), !is.na(count)) %>%
      group_by(gene_name) %>%
      summarise(count = sum(count), .groups = "drop")
  }, error = function(e) {
    message("Skipping file: ", file)
    NULL
  })

  if (!is.null(dat)) {
    sample_id <- basename(dirname(file))
    colnames(dat)[colnames(dat) == "count"] <- sample_id
    expr_list[[sample_id]] <- dat
  }
}

expr_kidney_merged <- reduce(expr_list, full_join, by = "gene_name") %>%
  column_to_rownames("gene_name")

expr_kidney_merged[is.na(expr_kidney_merged)] <- 0

query_expr <- GDCquery(
  project = "CPTAC-3",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

metadata_expr <- getResults(
  query_expr,
  cols = c("file_id", "cases.submitter_id", "sample_type", "file_name")
)

meta_kidney <- metadata_expr %>%
  filter(
    file_id %in% colnames(expr_kidney_merged),
    sample_type %in% c("Primary Tumor", "Solid Tissue Normal")
  ) %>%
  dplyr::select(sample = file_id, sample_type) %>%
  mutate(
    condition = ifelse(grepl("Normal", sample_type), "Normal", "Tumor")
  ) %>%
  distinct(sample, .keep_all = TRUE) %>%
  column_to_rownames("sample") %>%
  dplyr::select(condition)

expr_kidney_counts <- expr_kidney_merged[, rownames(meta_kidney), drop = FALSE]

dds_kidney <- DESeqDataSetFromMatrix(
  countData = round(expr_kidney_counts),
  colData = meta_kidney,
  design = ~ condition
)

vsd_kidney <- vst(dds_kidney, blind = TRUE)
expr_kidney <- assay(vsd_kidney)
expr_kidney <- collapse_duplicate_symbols(expr_kidney)

# ============================================================
# 4. CRC dataset: Korean CRC
# ============================================================

setwd(crc_dir)

raw_crc <- read.delim(
  "CMCBSN_expectedcount_342.txt",
  check.names = FALSE
)

counts_crc <- as.data.frame(raw_crc[, -1, drop = FALSE])
rownames(counts_crc) <- raw_crc[[1]]
counts_crc[] <- lapply(counts_crc, function(x) as.numeric(as.character(x)))

counts_crc <- collapse_duplicate_symbols(counts_crc)

crc_sample_ids <- colnames(counts_crc)

crc_tissue <- ifelse(
  grepl("-T-", crc_sample_ids),
  "T",
  ifelse(grepl("-N-", crc_sample_ids), "N", NA)
)

meta_crc <- data.frame(
  sample = crc_sample_ids,
  Tissue = crc_tissue,
  row.names = crc_sample_ids,
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(Tissue)) %>%
  mutate(condition = ifelse(Tissue == "T", "Tumor", "Normal")) %>%
  dplyr::select(condition)

counts_crc <- counts_crc[, rownames(meta_crc), drop = FALSE]

dds_crc <- DESeqDataSetFromMatrix(
  countData = round(counts_crc),
  colData = meta_crc,
  design = ~ condition
)

vsd_crc <- vst(dds_crc, blind = TRUE)
expr_crc <- assay(vsd_crc)
expr_crc <- collapse_duplicate_symbols(expr_crc)

# ============================================================
# Make 15-gene PCA plots
# ============================================================

p_brca_15 <- make_pca_plot(
  expr_mat = expr_brca,
  metadata = meta_brca,
  genes = panel_15,
  dataset_name = "BRCA-9"
)

p_lung_15 <- make_pca_plot(
  expr_mat = expr_lung,
  metadata = meta_lung,
  genes = panel_15,
  dataset_name = "Lung"
)

p_kidney_15 <- make_pca_plot(
  expr_mat = expr_kidney,
  metadata = meta_kidney,
  genes = panel_15,
  dataset_name = "Kidney"
)

p_crc_15 <- make_pca_plot(
  expr_mat = expr_crc,
  metadata = meta_crc,
  genes = panel_15,
  dataset_name = "CRC"
)

p_pca_15 <- (
  p_brca_15 | p_lung_15 | p_kidney_15 | p_crc_15
) +
  plot_annotation(
    title = "PCA of external validation datasets using the 15-gene panel",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16)
    )
  )

print(p_pca_15)

