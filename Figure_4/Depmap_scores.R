#!/usr/bin/env Rscript

# DepMap PolyGlu-TP score in selected cancer cell lines
#
# This script:
# 1. trains a pan-cancer RF classifier on matched TCGA tumour-normal cohorts
#    using log2(TPM + 1) expression for a fixed 9-gene panel
# 2. applies the trained RF model to DepMap log2(TPM + 1) expression data
# 3. joins DepMap model metadata
# 4. extracts selected cell lines
# 5. plots the PolyGlu-TP score across all selected cell lines
#
# Input:
# - dds_*_matched_group.rds files
# - OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv
# - Model.csv
#
# Notes:
# - Edit all input paths below before running
# - No z-scoring is applied
# - The fixed 9-gene panel is:
#   TTLL4, TTLL5, TTLL6, TTLL7, TTLL11, AGBL1, AGBL3, AGBL4, AGBL5
# - The TCGA training matrix is converted from counts to TPM before RF training
# - DepMap expression input is expected to already be log2(TPM + 1)

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(biomaRt)
  library(caret)
  library(randomForest)
  library(pROC)
  library(readr)
  library(ggplot2)
})

set.seed(123)

# ------------------------------------------------------------------
# Input paths
# ------------------------------------------------------------------

tcga_dir <- "path-to-tcga-dds-files"
depmap_expr_file <- "path-to-DepMap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
depmap_model_file <- "path-to-DepMap/Model.csv"
out_dir <- "path-to-output-folder"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(tcga_dir)) {
  stop("Please set 'tcga_dir' to the folder containing dds_*_matched_group.rds files.")
}
if (!file.exists(depmap_expr_file)) {
  stop("Could not find: ", depmap_expr_file)
}
if (!file.exists(depmap_model_file)) {
  stop("Could not find: ", depmap_model_file)
}

# ------------------------------------------------------------------
# Fixed 9-gene panel
# ------------------------------------------------------------------

panel_genes <- c(
  "TTLL4", "TTLL5", "TTLL6", "TTLL7", "TTLL11",
  "AGBL1", "AGBL3", "AGBL4", "AGBL5"
)

# ------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------

strip_ver <- function(x) sub("\\..*$", "", x)

aggregate_dup_symbols <- function(mat, symbols, mode = "mean") {
  stopifnot(nrow(mat) == length(symbols))

  mat <- as.matrix(mat)
  mode(mat) <- "numeric"

  keep <- !is.na(symbols) & symbols != ""
  mat <- mat[keep, , drop = FALSE]
  symbols <- symbols[keep]

  split_idx <- split(seq_len(nrow(mat)), symbols)

  out <- lapply(split_idx, function(i) {
    submat <- mat[i, , drop = FALSE]
    if (mode == "mean") {
      colMeans(submat, na.rm = TRUE)
    } else {
      colSums(submat, na.rm = TRUE)
    }
  })

  out <- do.call(rbind, out)
  rownames(out) <- names(split_idx)
  out
}

tpm_from_counts <- function(counts_mat, length_kb_vec) {
  stopifnot(all(rownames(counts_mat) %in% names(length_kb_vec)))

  length_kb_vec <- length_kb_vec[rownames(counts_mat)]
  rate <- counts_mat / length_kb_vec
  t(t(rate) / colSums(rate)) * 1e6
}

pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit)) hit[1] else NA_character_
}

# ------------------------------------------------------------------
# Ensembl symbol and gene length mapping
# ------------------------------------------------------------------

message("Fetching Ensembl symbol and gene length mapping...")

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_map <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "transcript_length"),
  mart = ensembl
) %>%
  filter(external_gene_name != "") %>%
  group_by(ensembl_gene_id, external_gene_name) %>%
  summarise(
    gene_length = mean(transcript_length, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

len_by_symbol <- gene_map %>%
  filter(!is.na(gene_length)) %>%
  group_by(external_gene_name) %>%
  summarise(
    length_kb = mean(gene_length) / 1000,
    .groups = "drop"
  )

# ------------------------------------------------------------------
# Build TCGA log2(TPM + 1) training matrix
# ------------------------------------------------------------------

message("Building TCGA log2(TPM + 1) training matrix...")

dds_files <- list.files(
  tcga_dir,
  pattern = "dds_.*_matched_group\\.rds$",
  full.names = TRUE
)

if (length(dds_files) == 0) {
  stop("No matched TCGA DESeq2 files found in 'tcga_dir'.")
}

tcga_list <- list()

for (file in dds_files) {
  project <- sub("dds_(.*)_matched_group\\.rds$", "\\1", basename(file))
  message("  processing: ", project)

  dds <- readRDS(file)

  counts_raw <- counts(dds)
  ens <- strip_ver(rownames(counts_raw))
  sym <- gene_map$external_gene_name[match(ens, gene_map$ensembl_gene_id)]

  counts_sym <- aggregate_dup_symbols(counts_raw, sym, mode = "sum")

  len_vec <- len_by_symbol$length_kb[
    match(rownames(counts_sym), len_by_symbol$external_gene_name)
  ]
  names(len_vec) <- rownames(counts_sym)

  keep_len <- !is.na(len_vec) & len_vec > 0

  counts_sym <- counts_sym[keep_len, , drop = FALSE]
  len_vec <- len_vec[keep_len]

  tpm_sym <- tpm_from_counts(counts_sym, len_vec)
  log_tpm_sym <- log2(tpm_sym + 1)

  matched <- intersect(panel_genes, rownames(log_tpm_sym))
  if (!length(matched)) next

  expr_tpm <- as.data.frame(t(log_tpm_sym[matched, , drop = FALSE]))
  expr_tpm$sample <- rownames(expr_tpm)
  expr_tpm$project <- project

  cond <- tryCatch(
    colData(dds)[expr_tpm$sample, "condition"],
    error = function(e) NULL
  )

  if (is.null(cond)) next

  expr_tpm$condition <- as.character(cond)

  tcga_list[[length(tcga_list) + 1]] <- expr_tpm
}

tcga_tpm_df <- bind_rows(tcga_list)

if (nrow(tcga_tpm_df) == 0) {
  stop("No TCGA TPM training data could be assembled.")
}

missing_cols <- setdiff(panel_genes, colnames(tcga_tpm_df))
if (length(missing_cols)) {
  for (g in missing_cols) tcga_tpm_df[[g]] <- NA_real_
}

tcga_tpm_df <- tcga_tpm_df %>%
  mutate(
    condition_bin = factor(
      ifelse(condition == "Primary Tumor", "Tumor", "Normal"),
      levels = c("Normal", "Tumor")
    )
  )

x_tpm <- tcga_tpm_df[, panel_genes, drop = FALSE]
y_tpm <- droplevels(tcga_tpm_df$condition_bin)

keep <- complete.cases(x_tpm) & !is.na(y_tpm)

x_tpm <- x_tpm[keep, , drop = FALSE]
y_tpm <- y_tpm[keep]

if (nrow(x_tpm) <= 50 || nlevels(y_tpm) != 2) {
  stop("Insufficient complete TCGA TPM training data.")
}

# ------------------------------------------------------------------
# Train RF model in log2(TPM + 1) space
# ------------------------------------------------------------------

message("Training pan-cancer RF model in log2(TPM + 1) space...")

ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

rf_pan_tpm <- caret::train(
  x = x_tpm,
  y = y_tpm,
  method = "rf",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5
)

cv_tpm <- rf_pan_tpm$pred %>%
  filter(mtry == rf_pan_tpm$bestTune$mtry)

roc_tpm <- pROC::roc(
  response = cv_tpm$obs,
  predictor = cv_tpm$Tumor,
  quiet = TRUE
)

cat(sprintf(
  "\nTCGA TPM RF 5-fold CV AUC = %.3f\n",
  as.numeric(pROC::auc(roc_tpm))
))

saveRDS(
  rf_pan_tpm,
  file.path(out_dir, "rf_pan_tcga_9gene_logTPM_model.rds")
)

# ------------------------------------------------------------------
# Load DepMap expression and metadata
# ------------------------------------------------------------------

message("Loading DepMap expression and metadata...")

depmap_expr <- readr::read_csv(depmap_expr_file, show_col_types = FALSE)
depmap_model <- readr::read_csv(depmap_model_file, show_col_types = FALSE)

id_cols <- intersect(
  c(
    "...1", "ProfileID", "SequencingID",
    "ModelID", "IsDefaultEntryForModel",
    "ModelConditionID", "IsDefaultEntryForMC"
  ),
  colnames(depmap_expr)
)

if (length(id_cols) < 1) {
  stop("Could not identify ID columns in DepMap expression file.")
}

dep_ids <- depmap_expr %>%
  dplyr::select(all_of(id_cols))

gene_cols <- setdiff(colnames(depmap_expr), id_cols)
bare_gene <- sub(" \\(.*$", "", gene_cols)

col_index <- sapply(panel_genes, function(g) which(bare_gene == g)[1])
names(col_index) <- panel_genes

missing_dep <- panel_genes[is.na(col_index)]

if (length(missing_dep) > 0) {
  warning(
    "Panel genes missing in DepMap expression file: ",
    paste(missing_dep, collapse = ", ")
  )
}

present_genes <- panel_genes[!is.na(col_index)]

if (length(present_genes) < 3) {
  stop("Fewer than three panel genes are present in DepMap expression data.")
}

expr_dep <- matrix(
  NA_real_,
  nrow = nrow(depmap_expr),
  ncol = length(panel_genes),
  dimnames = list(NULL, panel_genes)
)

for (g in present_genes) {
  expr_dep[, g] <- depmap_expr[[gene_cols[col_index[g]]]]
}

x_dep <- as.data.frame(expr_dep)
x_dep[] <- lapply(x_dep, function(z) suppressWarnings(as.numeric(z)))

prob_dep <- predict(
  rf_pan_tpm,
  newdata = x_dep,
  type = "prob"
)[, "Tumor"]

scores_dep <- bind_cols(
  dep_ids,
  tibble(score_TPM_model = as.numeric(prob_dep))
)

# ------------------------------------------------------------------
# Join DepMap metadata
# ------------------------------------------------------------------

message("Joining DepMap Model.csv metadata...")

expr_modelid_col <- pick_col(scores_dep, c("ModelID"))
meta_modelid_col <- pick_col(depmap_model, c("ModelID", "model_id", "DepMap_ID"))

if (!is.na(expr_modelid_col) && !is.na(meta_modelid_col)) {
  scores_dep <- scores_dep %>%
    left_join(depmap_model, by = setNames(meta_modelid_col, expr_modelid_col))
} else {
  warning("Could not find ModelID join between expression and Model.csv.")
}

cellline_col <- pick_col(scores_dep, c("CellLineName", "cell_line_name", "cell_line"))
stripped_col <- pick_col(scores_dep, c("StrippedCellLineName", "stripped_cell_line_name"))
lineage_col <- pick_col(scores_dep, c("OncotreeLineage", "oncotree_lineage", "lineage"))

scores_dep <- scores_dep %>%
  mutate(
    CellLineName = if (!is.na(cellline_col)) .data[[cellline_col]] else NA_character_,
    StrippedCellLineName = case_when(
      !is.na(stripped_col) ~ as.character(.data[[stripped_col]]),
      !is.na(cellline_col) ~ toupper(gsub("[^A-Za-z0-9]", "", .data[[cellline_col]])),
      TRUE ~ NA_character_
    ),
    OncotreeLineage = if (!is.na(lineage_col)) .data[[lineage_col]] else NA_character_
  )

# ------------------------------------------------------------------
# Keep default DepMap entries and selected cell lines
# ------------------------------------------------------------------

scores_default <- scores_dep

if (all(c("IsDefaultEntryForModel", "IsDefaultEntryForMC") %in% names(scores_default))) {
  scores_default <- scores_default %>%
    filter(IsDefaultEntryForModel == "Yes", IsDefaultEntryForMC == "Yes")
}

scores_default <- scores_default %>%
  arrange(desc(score_TPM_model))

write_csv(
  scores_default,
  file.path(out_dir, "DepMap_PolyGlu_TPM_scores_all_default.csv")
)

target_map <- tibble::tribble(
  ~Stripped, ~cell_line,
  "HELA", "HeLa",
  "MDAMB468", "MDA-MB-468",
  "MDAMB231", "MDA-MB-231",
  "MDAMB453", "MDA-MB-453",
  "U2OS", "U2OS",
  "A549", "A549",
  "SW480", "SW480",
  "SW48", "SW48",
  "BT474", "BT-474",
  "HCT116", "HCT116",
  "HT29", "HT29"
)

dep_scores_all <- scores_default %>%
  mutate(StrippedCellLineName = toupper(StrippedCellLineName)) %>%
  inner_join(target_map, by = c("StrippedCellLineName" = "Stripped")) %>%
  transmute(
    cell_line = cell_line,
    dataset = "DepMap",
    Tumor_score = score_TPM_model
  )

if (nrow(dep_scores_all) == 0) {
  stop("None of the selected cell lines were found in DepMap default entries.")
}

write_csv(
  dep_scores_all,
  file.path(out_dir, "DepMap_selected_cell_lines_PolyGlu_TPM_scores.csv")
)

message("\nSelected cell lines found:")
print(table(dep_scores_all$cell_line))

# ------------------------------------------------------------------
# Plot PolyGlu-TP score across selected cell lines
# ------------------------------------------------------------------

summary_df <- dep_scores_all %>%
  group_by(cell_line) %>%
  summarise(
    mean_score = mean(Tumor_score, na.rm = TRUE),
    sd_score = sd(Tumor_score, na.rm = TRUE),
    se_score = ifelse(n() > 1, sd_score / sqrt(n()), 0),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_score))

cell_order <- rev(summary_df$cell_line)

dep_scores_all <- dep_scores_all %>%
  mutate(cell_line = factor(cell_line, levels = cell_order))

summary_df <- summary_df %>%
  mutate(cell_line = factor(cell_line, levels = cell_order))

p_dep_all <- ggplot(dep_scores_all, aes(x = cell_line, y = Tumor_score)) +
  geom_col(
    data = summary_df,
    aes(x = cell_line, y = mean_score, fill = mean_score),
    width = 0.6,
    color = NA,
    alpha = 0.9,
    inherit.aes = FALSE
  ) +
  geom_jitter(
    width = 0.12,
    size = 2,
    alpha = 0.8,
    color = "black"
  ) +
  geom_errorbar(
    data = summary_df,
    aes(
      x = cell_line,
      ymin = pmax(0, mean_score - se_score),
      ymax = pmin(1, mean_score + se_score)
    ),
    width = 0.25,
    color = "black",
    alpha = 0.7,
    inherit.aes = FALSE
  ) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    color = "grey40"
  ) +
  scale_fill_gradient(
    low = "pink",
    high = "#8B475D",
    name = "Mean score"
  ) +
  labs(
    title = "DepMap: PolyGlu-TP score across selected cell lines",
    x = NULL,
    y = "Predicted tumor probability"
  ) +
  coord_flip() +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title = element_text(color = "black", face = "bold"),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.margin = ggplot2::margin(15, 20, 15, 15)
  )

print(p_dep_all)

ggsave(
  file.path(out_dir, "DepMap_selected_cell_lines_PolyGlu_TP_score.pdf"),
  p_dep_all,
  width = 7,
  height = 5
)

ggsave(
  file.path(out_dir, "DepMap_selected_cell_lines_PolyGlu_TP_score.png"),
  p_dep_all,
  width = 7,
  height = 5,
  dpi = 300
)

message("\nFinished. Outputs written to:\n", out_dir)
