#!/usr/bin/env Rscript

# ==============================================================================
# TCGA pan-cancer PolyGlu-TP association with aneuploidy
#
# Analysis
# - Train the 9-gene pan-cancer random forest classifier on matched TCGA
#   tumour/normal datasets using VST-transformed expression values.
# - Score TCGA primary tumours and summarize multiple primary-tumour samples
#   at the patient level by the median expression of each panel gene.
# - Match PolyGlu-TP scores to published aneuploidy scores.
# - Define Low (<=25th percentile) and High (>=75th percentile) PolyGlu-TP
#   groups within each cancer type.
# - Compare aneuploidy scores using a cancer-type-stratified
#   Wilcoxon-Mann-Whitney (van Elteren) test.
# ==============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(biomaRt)
  library(caret)
  library(randomForest)
  library(readxl)
  library(stringr)
  library(ggplot2)
  library(coin)
})

set.seed(123)


# ------------------------------------------------------------------
# Input paths
# ------------------------------------------------------------------

tcga_dir <- "path-to-tcga-dds-files"
gdc_dir <- "path-to-gdc-rnaseq-rds-files"
aneuploidy_file <- "path-to-aneuploidy-score-file.xlsx"
out_dir <- "path-to-output-folder"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(tcga_dir)) {
  stop("Please set 'tcga_dir' to the folder containing dds_*_matched_group.rds files.")
}

if (!dir.exists(gdc_dir)) {
  stop("Please set 'gdc_dir' to the folder containing TCGA-<PROJECT>_RNAseq.rds files.")
}

if (!file.exists(aneuploidy_file)) {
  stop("Please set 'aneuploidy_file' to the Taylor et al. aneuploidy-score Excel file.")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

feature_genes <- c(
  "TTLL4", "TTLL5", "TTLL6", "TTLL7", "TTLL11",
  "AGBL1", "AGBL3", "AGBL4", "AGBL5"
)

# ==============================================================================
# HELPERS
# ==============================================================================

aggregate_dup_symbols <- function(mat, symbols) {
  ok <- !is.na(symbols) & symbols != ""
  mat <- mat[ok, , drop = FALSE]
  symbols <- symbols[ok]

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

find_col <- function(df, patterns) {
  nms <- names(df)

  for (pat in patterns) {
    ix <- which(grepl(pat, nms, ignore.case = TRUE))
    if (length(ix)) return(nms[ix[1]])
  }

  NA_character_
}

tcga_sample_type_code <- function(x) {
  s <- as.character(x)
  code <- stringr::str_match(s, "-(\\d\\d)[A-Za-z]?$")[, 2]
  suppressWarnings(as.integer(code))
}

is_tumour_tcga <- function(x) {
  code <- tcga_sample_type_code(x)
  !is.na(code) & code >= 1 & code <= 9
}

# ==============================================================================
# 1. TRAIN PAN-CANCER RANDOM FOREST
# ==============================================================================

message("== Training pan-cancer RF model ==")

stopifnot(dir.exists(tcga_dir))

dds_files <- list.files(
  tcga_dir,
  pattern = "dds_.*_matched_group\\.rds$",
  full.names = TRUE
)

stopifnot(length(dds_files) > 0)

ensembl <- biomaRt::useMart(
  "ensembl",
  dataset = "hsapiens_gene_ensembl"
)

gene_map <- biomaRt::getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  mart = ensembl
) %>%
  dplyr::filter(external_gene_name != "") %>%
  dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)

all_tcga <- list()

for (file in dds_files) {
  project <- sub(
    "dds_(.*)_matched_group\\.rds$",
    "\\1",
    basename(file)
  )

  message("  training: ", project)

  dds <- readRDS(file)
  vsd <- DESeq2::vst(dds, blind = TRUE)
  mat <- SummarizedExperiment::assay(vsd)

  ens <- sub("\\..*$", "", rownames(mat))
  sym <- gene_map$external_gene_name[
    match(ens, gene_map$ensembl_gene_id)
  ]

  mat_sym <- aggregate_dup_symbols(mat, sym)

  matched <- intersect(feature_genes, rownames(mat_sym))
  if (!length(matched)) next

  expr <- t(mat_sym[matched, , drop = FALSE]) %>%
    as.data.frame()

  expr$sample <- rownames(expr)

  cond <- tryCatch(
    SummarizedExperiment::colData(dds)[expr$sample, "condition"],
    error = function(e) NULL
  )

  if (is.null(cond)) next

  expr$condition_bin <- factor(
    ifelse(
      as.character(cond) == "Primary Tumor",
      "Tumor",
      "Normal"
    ),
    levels = c("Normal", "Tumor")
  )

  all_tcga[[length(all_tcga) + 1]] <- expr
}

tcga_train <- dplyr::bind_rows(all_tcga)
stopifnot(nrow(tcga_train) > 0)

missing_train_genes <- setdiff(feature_genes, colnames(tcga_train))
if (length(missing_train_genes)) {
  stop(
    "Training matrix is missing required genes: ",
    paste(missing_train_genes, collapse = ", ")
  )
}

tcga_train <- tcga_train %>%
  dplyr::filter(
    dplyr::if_all(
      dplyr::all_of(feature_genes),
      ~ !is.na(.)
    )
  )

ctrl <- caret::trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = caret::twoClassSummary
)

set.seed(123)

rf_pan <- caret::train(
  x = tcga_train[, feature_genes, drop = FALSE],
  y = tcga_train$condition_bin,
  method = "rf",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5
)

saveRDS(
  rf_pan,
  file.path(out_dir, "rf_pan_tcga_9gene_RAWVST.rds")
)

# ==============================================================================
# 2. READ ANEUPLOIDY SCORE TABLE
# ==============================================================================

message("== Reading aneuploidy-score table ==")

preview <- readxl::read_excel(
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

stopifnot(!is.na(header_row))

supp <- readxl::read_excel(
  supp_excel,
  sheet = 1,
  skip = header_row - 1,
  col_names = TRUE
) %>%
  as.data.frame()

sample_col <- find_col(
  supp,
  c("^sample$", "sample\\s*id", "tcga")
)

as_col <- find_col(
  supp,
  c(
    "aneuploidy.*score",
    "aneuploidyscore",
    "aneuploidy\\s*score\\(as\\)",
    "as\\)$",
    "\\bas\\b"
  )
)

stopifnot(
  !is.na(sample_col),
  !is.na(as_col)
)

as_tbl <- supp %>%
  dplyr::transmute(
    sample = as.character(.data[[sample_col]]),
    AS = suppressWarnings(as.numeric(.data[[as_col]]))
  ) %>%
  dplyr::filter(
    !is.na(sample),
    !is.na(AS)
  ) %>%
  dplyr::mutate(
    sample = stringr::str_replace_all(sample, "\\.", "-"),
    sample = stringr::str_replace_all(sample, "_", "-"),
    tumour = is_tumour_tcga(sample),
    patient = substr(sample, 1, 12)
  ) %>%
  dplyr::filter(tumour) %>%
  dplyr::select(patient, AS) %>%
  dplyr::distinct(patient, .keep_all = TRUE)

# ==============================================================================
# 3. SCORE TCGA PRIMARY-TUMOUR PATIENTS
# ==============================================================================

message("== Scoring TCGA primary-tumour patients ==")

stopifnot(dir.exists(gdc_root))

rds_files <- list.files(
  gdc_root,
  pattern = "^TCGA-[A-Z0-9]+_RNAseq\\.rds$",
  full.names = TRUE
)

stopifnot(length(rds_files) > 0)

score_one_project <- function(rds_path) {
  proj <- sub(
    "^TCGA-([A-Z0-9]+)_RNAseq\\.rds$",
    "\\1",
    basename(rds_path)
  )

  message("  scoring: ", proj)

  se <- readRDS(rds_path)
  meta <- as.data.frame(SummarizedExperiment::colData(se))

  dds <- DESeq2::DESeqDataSet(se, design = ~ 1)
  vsd <- DESeq2::vst(dds, blind = TRUE)
  mat <- SummarizedExperiment::assay(vsd)

  ens <- sub("\\..*$", "", rownames(mat))
  sym <- gene_map$external_gene_name[
    match(ens, gene_map$ensembl_gene_id)
  ]

  mat_sym <- aggregate_dup_symbols(mat, sym)

  matched <- intersect(feature_genes, rownames(mat_sym))

  if (length(matched) < length(feature_genes)) {
    message(
      "    skipping ",
      proj,
      ": only ",
      length(matched),
      "/",
      length(feature_genes),
      " panel genes found."
    )
    return(NULL)
  }

  expr <- t(mat_sym[feature_genes, , drop = FALSE]) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample")

  if ("sample_type" %in% names(meta)) {
    tumour_samples <- rownames(meta)[meta$sample_type == "Primary Tumor"]

    if (!length(tumour_samples)) {
      tumour_samples <- expr$sample[
        is_tumour_tcga(expr$sample)
      ]
    }
  } else {
    tumour_samples <- expr$sample[
      is_tumour_tcga(expr$sample)
    ]
  }

  expr_pat <- expr %>%
    dplyr::filter(sample %in% tumour_samples) %>%
    dplyr::mutate(
      patient = substr(sample, 1, 12)
    ) %>%
    dplyr::group_by(patient) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(feature_genes),
        ~ median(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Cancer_Type = paste0("TCGA-", proj)
    )

  X <- expr_pat %>%
    dplyr::select(
      dplyr::all_of(feature_genes)
    )

  expr_pat$rf_prob_pan <- as.numeric(
    predict(
      rf_pan,
      newdata = X,
      type = "prob"
    )[, "Tumor"]
  )

  expr_pat
}

tcga_scored <- dplyr::bind_rows(
  lapply(rds_files, score_one_project)
)

stopifnot(nrow(tcga_scored) > 0)

# ==============================================================================
# 4. MATCH POLYGLU-TP SCORES TO ANEUPLOIDY SCORES
# ==============================================================================

dat <- tcga_scored %>%
  dplyr::inner_join(
    as_tbl,
    by = "patient"
  )

readr::write_csv(
  dat,
  file.path(out_dir, "tcga_scored_tumours_with_AS.csv")
)

# ==============================================================================
# 5. DEFINE WITHIN-CANCER LOW/HIGH GROUPS
# ==============================================================================

dat_hl <- dat %>%
  dplyr::filter(
    !is.na(rf_prob_pan),
    !is.na(AS),
    !is.na(Cancer_Type)
  ) %>%
  dplyr::group_by(Cancer_Type) %>%
  dplyr::mutate(
    q1 = stats::quantile(
      rf_prob_pan,
      0.25,
      na.rm = TRUE
    ),
    q3 = stats::quantile(
      rf_prob_pan,
      0.75,
      na.rm = TRUE
    ),
    group = dplyr::case_when(
      rf_prob_pan <= q1 ~ "Low",
      rf_prob_pan >= q3 ~ "High",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(group)) %>%
  dplyr::mutate(
    group = factor(
      group,
      levels = c("Low", "High")
    ),
    Cancer_Type = factor(Cancer_Type)
  )

# ==============================================================================
# 6. DESCRIPTIVE STATISTICS AND STRATIFIED TEST
# ==============================================================================

stats <- dat_hl %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(
    n = dplyr::n(),
    median_AS = median(AS, na.rm = TRUE),
    IQR_AS = IQR(AS, na.rm = TRUE),
    .groups = "drop"
  )

print(as.data.frame(stats))

low_n <- stats$n[stats$group == "Low"]
high_n <- stats$n[stats$group == "High"]

low_med <- stats$median_AS[stats$group == "Low"]
high_med <- stats$median_AS[stats$group == "High"]

delta_med <- high_med - low_med

ve_test <- coin::wilcox_test(
  AS ~ group | Cancer_Type,
  data = dat_hl,
  distribution = "asymptotic"
)

ve_p <- coin::pvalue(ve_test)

print(ve_test)

cat("\n========================================\n")
cat("ANEUPLOIDY ANALYSIS\n")
cat("Within-cancer PolyGlu-TP quartiles\n")
cat("========================================\n")
cat("Low: n =", low_n, "| median AS =", low_med, "\n")
cat("High: n =", high_n, "| median AS =", high_med, "\n")
cat("Delta median =", delta_med, "\n")
cat(
  "Cancer-type-stratified Wilcoxon p =",
  format(ve_p, scientific = TRUE),
  "\n"
)

cancer_counts <- dat_hl %>%
  dplyr::count(
    Cancer_Type,
    group,
    name = "n"
  ) %>%
  tidyr::pivot_wider(
    names_from = group,
    values_from = n,
    values_fill = 0
  )

print(as.data.frame(cancer_counts))

# ==============================================================================
# 7. SAVE STATISTICS
# ==============================================================================

stats_out <- tibble::tibble(
  analysis = "Within-cancer PolyGlu-TP quartiles",
  low_n = low_n,
  high_n = high_n,
  low_median_AS = low_med,
  high_median_AS = high_med,
  delta_median = delta_med,
  stratified_wilcoxon_p = ve_p
)

readr::write_csv(
  stats_out,
  file.path(
    out_dir,
    "AS_highlow_withinCancer_RFscore_stats.csv"
  )
)

readr::write_csv(
  cancer_counts,
  file.path(
    out_dir,
    "AS_highlow_withinCancer_counts_by_cancer.csv"
  )
)

# ==============================================================================
# 8. PLOT
# ==============================================================================

p_label <- if (ve_p < 2.2e-16) {
  "p < 2.2 \u00d7 10\u207b\u00b9\u2076"
} else {
  paste0(
    "p = ",
    format(
      ve_p,
      scientific = TRUE,
      digits = 2
    )
  )
}

stat_label <- paste0(
  "Stratified Wilcoxon ",
  p_label,
  "\n\u0394median = ",
  delta_med
)

p_as <- ggplot2::ggplot(
  dat_hl,
  ggplot2::aes(
    x = group,
    y = AS,
    fill = group
  )
) +
  ggplot2::geom_boxplot(
    width = 0.62,
    linewidth = 0.7,
    outlier.alpha = 0.18
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Low" = "pink",
      "High" = "#8B475D"
    )
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      "Low" = "Low\n(\u226425th percentile)",
      "High" = "High\n(\u226575th percentile)"
    )
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Aneuploidy score",
    fill = NULL
  ) +
  ggplot2::annotate(
    "text",
    x = 1.5,
    y = max(dat_hl$AS, na.rm = TRUE) * 0.97,
    label = stat_label,
    hjust = 0.5,
    vjust = 1,
    size = 4.2
  ) +
  ggplot2::theme_classic(
    base_size = 14
  ) +
  ggplot2::theme(
    legend.position = "none",
    axis.text = ggplot2::element_text(
      color = "black"
    ),
    axis.title = ggplot2::element_text(
      color = "black",
      face = "bold"
    )
  )

print(p_as)

