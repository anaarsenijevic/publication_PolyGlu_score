# Iterative random forest feature reduction:
# train on TCGA matched tumour-normal cohorts and test on external BRCA patient samples
#
#
# Input:
# - dds_*_matched_group.rds files
# - Salmon quant.sf files for the patient cohort
# - Gencode GTF annotation file
#
# Notes:
# - Edit the paths below before running
# - Ensembl IDs are mapped to gene symbols using biomaRt
# - Duplicate symbols are not collapsed here; one symbol match is retained per row as in the original workflow

suppressPackageStartupMessages({
  library(DESeq2)
  library(tximport)
  library(tidyverse)
  library(biomaRt)
  library(randomForest)
  library(caret)
  library(ggplot2)
  library(pROC)
  library(AnnotationDbi)
  library(GenomicFeatures)
  library(forcats)
  library(tibble)
  library(purrr)
})

# ------------------------------------------------------------------
# Input paths
# ------------------------------------------------------------------

tcga_dir       <- "path-to-tcga-dds-files"
annotation_dir <- "path-to-annotation-folder"
salmon_dir     <- "path-to-salmon-folder"

gtf_file <- file.path(annotation_dir, "annotations", "gencode.v29.annotation.gtf")

if (!dir.exists(tcga_dir)) {
  stop("Please set 'tcga_dir' to the folder containing dds_*_matched_group.rds files.")
}

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
# Define gene panel
# ------------------------------------------------------------------

ttll_genes <- c(
  "TTLL1", "TTLL2", "TTLL4", "TTLL5", "TTLL6", "TTLL7", "TTLL9", "TTLL11", "TTLL13"
)

ccp_genes <- c(
  "AGBL1", "AGBL2", "AGBL3", "AGBL4", "AGBL5", "AGTPBP1"
)

feature_genes <- unique(c(ttll_genes, ccp_genes))

symbol_labels <- c(
  "TTLL1", "TTLL2", "TTLL4", "TTLL5", "TTLL6", "TTLL7", "TTLL9", "TTLL11", "TTLL13",
  "CCP1", "CCP2", "CCP3", "CCP4", "CCP5", "CCP6"
)
gene_label_map <- setNames(symbol_labels, feature_genes)

# ------------------------------------------------------------------
# Gene symbol mapping using biomaRt
# ------------------------------------------------------------------

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_map <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  mart = ensembl
) %>%
  filter(external_gene_name != "") %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

# ==================================================================
# 1. LOAD AND PROCESS TCGA DESeq2 OBJECTS
# ==================================================================

dds_files <- list.files(
  tcga_dir,
  pattern = "dds_.*_matched_group\\.rds$",
  full.names = TRUE
)

if (length(dds_files) == 0) {
  stop("No TCGA DESeq2 files found in 'tcga_dir'.")
}

all_data <- list()

for (file in dds_files) {
  project <- sub("dds_(.*)_matched_group\\.rds", "\\1", basename(file))
  cat("Checking:", project, "\n")

  dds <- readRDS(file)

  ensembl_ids <- sub("\\..*$", "", rownames(dds))
  symbols <- gene_map$external_gene_name[match(ensembl_ids, gene_map$ensembl_gene_id)]
  names(symbols) <- ensembl_ids

  vsd <- vst(dds, blind = TRUE)
  norm_counts <- assay(vsd)
  rownames(norm_counts) <- symbols

  matched_genes <- intersect(feature_genes, rownames(norm_counts))
  if (length(matched_genes) == 0) {
    cat("No matched genes in", project, "\n")
    next
  }

  expr <- t(norm_counts[matched_genes, , drop = FALSE]) %>% as.data.frame()
  expr$sample <- rownames(expr)

  cond <- tryCatch({
    colData(dds)[expr$sample, "condition"]
  }, error = function(e) {
    cat("Sample mismatch in", project, "\n")
    return(NULL)
  })

  if (is.null(cond)) next

  expr$condition <- as.character(cond)
  expr$project <- project
  all_data[[length(all_data) + 1]] <- expr

  cat("Added:", nrow(expr), "samples from", project, "\n")
}

tcga_df <- bind_rows(all_data)

tcga_df <- tcga_df %>%
  dplyr::select(any_of(feature_genes), condition, sample, project) %>%
  filter(if_all(any_of(feature_genes), ~ !is.na(.))) %>%
  mutate(
    condition_bin = factor(
      ifelse(condition == "Primary Tumor", "Tumor", "Normal"),
      levels = c("Normal", "Tumor")
    )
  )

if (nrow(tcga_df) == 0) {
  stop("No TCGA samples available after preprocessing.")
}

# ==================================================================
# 2. LOAD AND PROCESS BRCA PATIENT SAMPLES
# ==================================================================

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

files <- file.path(salmon_dir, samples_patient$path)
names(files) <- samples_patient$sample

missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0) {
  stop("The following Salmon files were not found:\n", paste(missing_files, collapse = "\n"))
}

txi_patient <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE
)

sample_table <- samples_patient %>%
  dplyr::select(sample, tissue, group) %>%
  column_to_rownames("sample")

dds_patient <- DESeqDataSetFromTximport(
  txi_patient,
  colData = sample_table,
  design = ~ tissue
)

vsd_patient <- vst(dds_patient, blind = TRUE)
vst_expr <- assay(vsd_patient)

ensembl_ids_patient <- sub("\\..*", "", rownames(vst_expr))

gene_map_patient <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = ensembl_ids_patient,
  mart = ensembl
)

symbols_patient <- gene_map_patient$external_gene_name[
  match(ensembl_ids_patient, gene_map_patient$ensembl_gene_id)
]

rownames(vst_expr) <- symbols_patient
vst_expr <- vst_expr[!is.na(rownames(vst_expr)) & rownames(vst_expr) != "", , drop = FALSE]

available_patient_genes <- intersect(feature_genes, rownames(vst_expr))
missing_patient_genes <- setdiff(feature_genes, available_patient_genes)

if (length(available_patient_genes) == 0) {
  stop("No feature genes were found in the patient expression matrix.")
}

patient_expr_df <- t(vst_expr[available_patient_genes, , drop = FALSE]) %>%
  as.data.frame()

for (g in missing_patient_genes) {
  patient_expr_df[[g]] <- NA_real_
}

patient_expr_df <- patient_expr_df[, feature_genes, drop = FALSE]

for (g in colnames(patient_expr_df)) {
  if (anyNA(patient_expr_df[[g]])) {
    med_g <- median(patient_expr_df[[g]], na.rm = TRUE)
    patient_expr_df[[g]][is.na(patient_expr_df[[g]])] <- if (is.finite(med_g)) med_g else 0
  }
}

# ==================================================================
# 3. ITERATIVE MODEL LOOP
# ==================================================================

y <- tcga_df$condition_bin

current_genes <- feature_genes[
  feature_genes %in% intersect(colnames(tcga_df), colnames(patient_expr_df))
]

iterations <- list()

if (length(current_genes) < 6) {
  stop("Fewer than 6 shared genes are available between TCGA and patient datasets.")
}

for (i in seq(length(current_genes), 6, by = -1)) {
  valid_genes <- intersect(current_genes, intersect(colnames(tcga_df), colnames(patient_expr_df)))

  if (length(valid_genes) < 6) break

  cat("Iteration with", length(valid_genes), "genes\n")

  x_tcga <- tcga_df[, valid_genes, drop = FALSE]
  x_patient <- patient_expr_df[, valid_genes, drop = FALSE]

  ctrl <- trainControl(
    method = "cv",
    number = 5,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final"
  )

  set.seed(123)
  rf_cv <- train(
    x = x_tcga,
    y = y,
    method = "rf",
    trControl = ctrl,
    metric = "ROC"
  )

  best_mtry <- rf_cv$bestTune$mtry

  final_model <- randomForest(
    x = x_tcga,
    y = y,
    mtry = best_mtry,
    ntree = 500,
    importance = TRUE
  )

  importance_df <- importance(final_model, type = 1) %>%
    as.data.frame() %>%
    rownames_to_column("Gene") %>%
    arrange(desc(MeanDecreaseAccuracy))

  prob_patient <- predict(final_model, newdata = x_patient, type = "prob")
  pred_patient <- predict(final_model, newdata = x_patient)

  iterations[[paste0("Genes_", length(valid_genes))]] <- list(
    genes = valid_genes,
    auc = max(rf_cv$results$ROC, na.rm = TRUE),
    acc = max(rf_cv$results$Accuracy, na.rm = TRUE),
    importance = importance_df,
    patient_probs = prob_patient,
    patient_preds = pred_patient,
    model = final_model,
    rf_cv = rf_cv
  )

  least_gene <- tail(importance_df$Gene, 1)
  current_genes <- setdiff(valid_genes, least_gene)
}

# ==================================================================
# 4. COMPARE TCGA CV AUC VS PATIENT AUC
# ==================================================================

truth_labels <- factor(
  ifelse(samples_patient$tissue == "PT", "Tumor", "Normal"),
  levels = c("Normal", "Tumor")
)

auc_df <- map_dfr(names(iterations), function(name) {
  iter <- iterations[[name]]
  gene_count <- length(iter$genes)

  tcga_auc <- iter$auc

  if (!is.null(iter$patient_probs)) {
    roc_obj <- tryCatch(
      roc(truth_labels, iter$patient_probs[, "Tumor"], quiet = TRUE),
      error = function(e) NULL
    )
    patient_auc <- if (!is.null(roc_obj)) as.numeric(auc(roc_obj)) else NA_real_
  } else {
    patient_auc <- NA_real_
  }

  tibble(
    Genes = gene_count,
    TCGA_AUC = tcga_auc,
    Patient_AUC = patient_auc
  )
})

auc_df_long <- pivot_longer(
  auc_df,
  cols = c(TCGA_AUC, Patient_AUC),
  names_to = "Dataset",
  values_to = "AUC"
)

ggplot(auc_df_long, aes(x = Genes, y = AUC, color = Dataset)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_x_reverse(breaks = sort(unique(auc_df$Genes), decreasing = TRUE)) +
  scale_color_manual(values = c("TCGA_AUC" = "darkgreen", "Patient_AUC" = "goldenrod")) +
  labs(
    title = "AUC vs Gene Panel Size",
    subtitle = "Cross-validated TCGA vs external patient samples",
    x = "Number of Genes in Model",
    y = "AUC",
    color = "Dataset"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# ==================================================================
# 5. PER-ITERATION IMPORTANCE PLOTS
# ==================================================================

iter_names <- names(iterations)

for (k in seq_along(iter_names)) {
  curr <- iter_names[k]
  next_iter <- if (k < length(iter_names)) iter_names[k + 1] else NA

  imp <- iterations[[curr]]$importance %>%
    mutate(
      Gene = dplyr::recode(Gene, !!!gene_label_map),
      RemovedNext = if (!is.na(next_iter)) {
        dropped <- setdiff(iterations[[curr]]$genes, iterations[[next_iter]]$genes)
        recoded_drop <- dplyr::recode(dropped, !!!gene_label_map)
        Gene %in% recoded_drop
      } else {
        FALSE
      },
      Gene = forcats::fct_reorder(Gene, MeanDecreaseAccuracy)
    )

  p <- ggplot(imp, aes(x = Gene, y = MeanDecreaseAccuracy, fill = RemovedNext)) +
    geom_col(width = 0.7, show.legend = FALSE) +
    coord_flip() +
    labs(
      title = paste0("RF Permutation Importance — ", gsub("_", " ", curr)),
      subtitle = if (!is.na(next_iter)) "Highlighted = gene removed next round" else NULL,
      x = "Gene",
      y = "Mean Decrease in Accuracy (OOB)"
    ) +
    scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "tomato")) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  print(p)
}

# ==================================================================
# 6. FACETED IMPORTANCE PLOT ACROSS ITERATIONS
# ==================================================================

importance_all <- map_dfr(names(iterations), function(nm) {
  imp <- iterations[[nm]]$importance
  imp$Iteration <- nm
  imp
})

importance_all <- importance_all %>%
  mutate(
    Gene = dplyr::recode(Gene, !!!gene_label_map),
    Genes = as.integer(gsub("Genes_", "", Iteration))
  )

importance_all <- importance_all %>%
  group_by(Iteration) %>%
  arrange(desc(MeanDecreaseAccuracy), Gene, .by_group = TRUE) %>%
  mutate(
    RemovedNext = Gene == dplyr::last(Gene),
    Gene_key = paste0(Iteration, "__", Gene)
  ) %>%
  ungroup()

levels_vec <- importance_all %>%
  group_by(Iteration) %>%
  summarize(levels = list(rev(Gene_key)), .groups = "drop") %>%
  tidyr::unnest(levels) %>%
  pull(levels)

importance_all <- importance_all %>%
  mutate(
    Gene_fac = factor(Gene_key, levels = levels_vec),
    Iteration = factor(
      paste0(Genes, " genes"),
      levels = paste0(sort(unique(Genes), decreasing = TRUE), " genes")
    )
  )

ggplot(importance_all, aes(x = Gene_fac, y = MeanDecreaseAccuracy, fill = RemovedNext)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  coord_flip() +
  facet_wrap(~ Iteration, scales = "free_y") +
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "grey70")) +
  scale_x_discrete(labels = function(x) sub("^.*__", "", x)) +
  labs(
    title = "Permutation Importance Across Iterations",
    x = "Gene",
    y = "Mean Decrease in Accuracy (OOB)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold", size = 12)
  )
