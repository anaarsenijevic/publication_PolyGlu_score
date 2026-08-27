#!/usr/bin/env Rscript

# ==============================================================================
# TCGA pan-cancer PolyGlu-TP association with TP53 mutation status
#
# Analysis
# - Train the 9-gene pan-cancer random forest classifier on matched TCGA
#   tumour/normal datasets using VST-transformed expression values.
# - Score TCGA primary tumours and summarize multiple primary-tumour samples
#   at the patient level by the median expression of each panel gene.
# - Obtain GDC Masked Somatic Mutation MAF files and identify patients with
#   at least one TP53 mutation.
# - Define Low (<=25th percentile) and High (>=75th percentile) PolyGlu-TP
#   groups separately within each cancer type.
# - Test the association between PolyGlu-TP group and TP53 mutation status
#   using a Cochran-Mantel-Haenszel test stratified by cancer type.
#
# Main outputs
# - rf_pan_tcga_9gene_RAWVST.rds
# - tcga_pancancer_scored_only.csv
# - TCGA_panCancer_TP53_mutated_patients.rds
# - TCGA_panCancer_TP53_mutated_patients.csv
# - TP53_counts_by_cancer_and_group.csv
# - TP53_withinCancer_CMH_stats.csv
# - TP53_withinCancer_CMH_barplot.pdf
# - TP53_withinCancer_CMH_barplot.png
# ==============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(caret)
  library(randomForest)
  library(scales)
  library(stringr)
  library(SummarizedExperiment)
  library(TCGAbiolinks)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

set.seed(123)

# ==============================================================================
# USER SETTINGS
# ==============================================================================

tcga_dir <- "path-to-tcga-dds-files"
gdc_dir <- "path-to-gdc-rnaseq-rds-files"
maf_dir <- "path-to-gdc-maf-cache"
out_dir <- "path-to-output-folder"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(maf_dir, recursive = TRUE, showWarnings = FALSE)

# GDCdownload created extraction problems on the network share in the working
# environment. A local cache is therefore used for the raw MAF files. The
# processed patient-level TP53 table is saved to out_dir and can be reused.

feature_genes <- c(
  "TTLL4", "TTLL5", "TTLL6", "TTLL7", "TTLL11",
  "AGBL1", "AGBL3", "AGBL4", "AGBL5"
)

cancer_types <- c(
  "TCGA-UCEC", "TCGA-THCA", "TCGA-STAD", "TCGA-READ", "TCGA-PRAD",
  "TCGA-LUSC", "TCGA-LUAD", "TCGA-LIHC", "TCGA-KIRP", "TCGA-KIRC",
  "TCGA-KICH", "TCGA-HNSC", "TCGA-ESCA", "TCGA-COAD", "TCGA-CESC",
  "TCGA-BLCA", "TCGA-BRCA"
)

pal_p53 <- c(
  "WT" = "#CAE1FF",
  "Mutated" = "lightsteelblue4"
)

# ==============================================================================
# HELPERS
# ==============================================================================

aggregate_dup_symbols <- function(
    mat,
    symbols,
    mode = c("mean", "sum")[1]
) {
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

map_ens_to_symbol <- function(ens_ids) {
  AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = ens_ids,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
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

  ens_ids <- sub("\\..*$", "", rownames(mat))
  syms <- unname(map_ens_to_symbol(ens_ids))
  mat_sym <- aggregate_dup_symbols(
    mat,
    syms,
    mode = "mean"
  )

  matched <- intersect(
    feature_genes,
    rownames(mat_sym)
  )

  if (!length(matched)) next

  expr <- t(mat_sym[matched, , drop = FALSE]) %>%
    as.data.frame()

  expr$sample <- rownames(expr)
  expr$project <- project

  cond <- tryCatch(
    SummarizedExperiment::colData(dds)[expr$sample, "condition"],
    error = function(e) NULL
  )

  if (is.null(cond)) next

  expr$condition <- as.character(cond)
  all_tcga[[length(all_tcga) + 1]] <- expr
}

tcga_df <- dplyr::bind_rows(all_tcga)
stopifnot(nrow(tcga_df) > 0)

for (g in setdiff(feature_genes, colnames(tcga_df))) {
  tcga_df[[g]] <- NA_real_
}

tcga_clean <- tcga_df %>%
  dplyr::select(
    dplyr::all_of(feature_genes),
    sample,
    project,
    condition
  ) %>%
  dplyr::filter(
    dplyr::if_all(
      dplyr::all_of(feature_genes),
      ~ !is.na(.)
    )
  ) %>%
  dplyr::mutate(
    condition_bin = factor(
      ifelse(
        condition == "Primary Tumor",
        "Tumor",
        "Normal"
      ),
      levels = c("Normal", "Tumor")
    )
  )

ctrl <- caret::trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = caret::twoClassSummary,
  savePredictions = TRUE
)

set.seed(123)

rf_pan <- caret::train(
  x = tcga_clean[, feature_genes, drop = FALSE],
  y = tcga_clean$condition_bin,
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
# 2. SCORE TCGA PRIMARY-TUMOUR PATIENTS
# ==============================================================================

message("== Scoring TCGA primary-tumour patients ==")

stopifnot(dir.exists(gdc_dir))

expression_list <- list()

for (cancer in cancer_types) {
  file_path <- file.path(
    gdc_dir,
    paste0(cancer, "_RNAseq.rds")
  )

  if (!file.exists(file_path)) {
    message("Skipping missing file: ", file_path)
    next
  }

  message("  scoring: ", cancer)

  se <- readRDS(file_path)
  meta <- as.data.frame(SummarizedExperiment::colData(se))

  assays_avail <- SummarizedExperiment::assayNames(se)

  counts_assay <- if ("unstranded" %in% assays_avail) {
    "unstranded"
  } else {
    assays_avail[1]
  }

  raw_counts <- SummarizedExperiment::assay(
    se,
    counts_assay
  )

  dds_all <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(raw_counts),
    colData = as.data.frame(
      SummarizedExperiment::colData(se)
    ),
    design = ~ 1
  )

  vsd_all <- DESeq2::vst(
    dds_all,
    blind = TRUE
  )

  mat_all <- SummarizedExperiment::assay(vsd_all)

  rd <- as.data.frame(
    SummarizedExperiment::rowData(se)
  )

  gene_name_col <- intersect(
    c(
      "gene_name",
      "hgnc_symbol",
      "external_gene_name"
    ),
    names(rd)
  )

  if (length(gene_name_col) > 0) {
    symbols_for_rows <- rd[[gene_name_col[1]]]
  } else {
    ens <- if ("gene_id" %in% names(rd)) {
      rd$gene_id
    } else {
      rownames(mat_all)
    }

    ens <- sub("\\..*$", "", ens)

    symbols_for_rows <- unname(
      map_ens_to_symbol(ens)
    )
  }

  mat_sym <- aggregate_dup_symbols(
    mat_all,
    symbols_for_rows,
    mode = "mean"
  )

  matched <- intersect(
    feature_genes,
    rownames(mat_sym)
  )

  if (length(matched) < length(feature_genes)) {
    message(
      "    skipping ",
      cancer,
      ": only ",
      length(matched),
      "/",
      length(feature_genes),
      " panel genes found."
    )
    next
  }

  expr <- t(
    mat_sym[feature_genes, , drop = FALSE]
  ) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample")

  if ("sample_type" %in% names(meta)) {
    tumor_samples <- rownames(meta)[
      meta$sample_type == "Primary Tumor"
    ]
  } else {
    tumor_samples <- expr$sample[
      substr(expr$sample, 14, 15) == "01"
    ]
  }

  expr_pat <- expr %>%
    dplyr::filter(
      sample %in% tumor_samples
    ) %>%
    dplyr::mutate(
      patient = substr(sample, 1, 12),
      Cancer_Type = cancer
    ) %>%
    dplyr::group_by(
      patient,
      Cancer_Type
    ) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(feature_genes),
        ~ median(.x, na.rm = TRUE)
      ),
      .groups = "drop"
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

  expression_list[[cancer]] <- expr_pat
}

tcga_pan_scored <- dplyr::bind_rows(
  expression_list
)

stopifnot(nrow(tcga_pan_scored) > 0)

readr::write_csv(
  tcga_pan_scored,
  file.path(
    out_dir,
    "tcga_pancancer_scored_only.csv"
  )
)

# ==============================================================================
# 3. DOWNLOAD/PARSE TP53 MUTATION DATA
# ==============================================================================

get_tp53_status_tcga <- function(
    projects,
    maf_dir
) {
  dir.create(
    maf_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  out <- list()

  for (prj in projects) {
    message("\n========================================")
    message("Processing mutation data: ", prj)
    message("========================================")

    q <- tryCatch(
      TCGAbiolinks::GDCquery(
        project = prj,
        data.category = "Simple Nucleotide Variation",
        data.type = "Masked Somatic Mutation",
        workflow.type =
          "Aliquot Ensemble Somatic Variant Merging and Masking"
      ),
      error = function(e) {
        message(
          "Query failed for ",
          prj,
          ": ",
          e$message
        )
        NULL
      }
    )

    if (is.null(q)) next

    query_results <- tryCatch(
      TCGAbiolinks::getResults(q),
      error = function(e) {
        message(
          "Could not extract query results for ",
          prj,
          ": ",
          e$message
        )
        NULL
      }
    )

    if (is.null(query_results)) next

    expected_n <- nrow(query_results)

    message(
      "GDC reports ",
      expected_n,
      " mutation files for ",
      prj
    )

    all_existing <- list.files(
      maf_dir,
      pattern = "\\.maf(\\.gz)?$",
      recursive = TRUE,
      full.names = TRUE
    )

    existing_mafs <- all_existing[
      grepl(
        paste0("[/\\\\]", prj, "[/\\\\]"),
        all_existing
      )
    ]

    existing_mafs <- existing_mafs[
      file.exists(existing_mafs)
    ]

    cached_n <- length(existing_mafs)

    message(
      "Valid cached MAF files found for ",
      prj,
      ": ",
      cached_n
    )

    if (cached_n < expected_n) {
      message(
        "Cache incomplete for ",
        prj,
        " (",
        cached_n,
        "/",
        expected_n,
        "). Redownloading project."
      )

      project_cache <- file.path(
        maf_dir,
        prj
      )

      if (dir.exists(project_cache)) {
        unlink(
          project_cache,
          recursive = TRUE,
          force = TRUE
        )
      }

      download_ok <- tryCatch(
        {
          TCGAbiolinks::GDCdownload(
            q,
            directory = maf_dir,
            method = "api",
            files.per.chunk = 50
          )

          TRUE
        },
        error = function(e) {
          message(
            "Download failed for ",
            prj,
            ": ",
            e$message
          )
          FALSE
        }
      )

      if (!download_ok) next
    } else {
      message(
        "Using complete cached MAF set for ",
        prj,
        " (",
        cached_n,
        "/",
        expected_n,
        ")"
      )
    }

    all_mafs <- list.files(
      maf_dir,
      pattern = "\\.maf(\\.gz)?$",
      recursive = TRUE,
      full.names = TRUE
    )

    maf_files <- all_mafs[
      grepl(
        paste0("[/\\\\]", prj, "[/\\\\]"),
        all_mafs
      )
    ]

    maf_files <- maf_files[
      file.exists(maf_files)
    ]

    message(
      "Readable MAF files found for ",
      prj,
      ": ",
      length(maf_files)
    )

    if (!length(maf_files)) {
      message(
        "No usable MAF files found for ",
        prj
      )
      next
    }

    if (length(maf_files) < expected_n) {
      message(
        "WARNING: only ",
        length(maf_files),
        " of ",
        expected_n,
        " expected MAF files were found for ",
        prj
      )
    }

    maf_list <- lapply(
      maf_files,
      function(f) {
        tryCatch(
          readr::read_tsv(
            f,
            comment = "#",
            show_col_types = FALSE,
            progress = FALSE,
            col_types = readr::cols(
              .default = readr::col_character()
            )
          ),
          error = function(e) {
            message(
              "Failed to read ",
              basename(f),
              ": ",
              e$message
            )
            NULL
          }
        )
      }
    )

    maf_list <- maf_list[
      !vapply(
        maf_list,
        is.null,
        logical(1)
      )
    ]

    if (!length(maf_list)) {
      message(
        "No readable MAF files for ",
        prj
      )
      next
    }

    maf <- dplyr::bind_rows(maf_list)

    required_cols <- c(
      "Hugo_Symbol",
      "Tumor_Sample_Barcode"
    )

    if (!all(required_cols %in% colnames(maf))) {
      message(
        "Required columns missing for ",
        prj
      )
      next
    }

    tp53_proj <- maf %>%
      dplyr::filter(
        Hugo_Symbol == "TP53"
      ) %>%
      dplyr::transmute(
        patient = substr(
          Tumor_Sample_Barcode,
          1,
          12
        ),
        Cancer_Type = prj,
        p53_status = "Mutated"
      ) %>%
      dplyr::filter(
        !is.na(patient),
        patient != ""
      ) %>%
      dplyr::distinct(
        patient,
        Cancer_Type,
        .keep_all = TRUE
      )

    message(
      "TP53-mutated patients identified for ",
      prj,
      ": ",
      nrow(tp53_proj)
    )

    out[[prj]] <- tp53_proj
  }

  if (!length(out)) {
    return(
      tibble::tibble(
        patient = character(),
        Cancer_Type = character(),
        p53_status = character()
      )
    )
  }

  dplyr::bind_rows(out)
}

tp53_rds <- file.path(
  out_dir,
  "TCGA_panCancer_TP53_mutated_patients.rds"
)

tp53_csv <- file.path(
  out_dir,
  "TCGA_panCancer_TP53_mutated_patients.csv"
)

if (file.exists(tp53_rds)) {
  message(
    "Loading previously processed TP53 mutation table: ",
    tp53_rds
  )

  tp53_pat <- readRDS(
    tp53_rds
  )
} else {
  message(
    "Processed TP53 table not found; downloading/parsing GDC MAFs."
  )

  tp53_pat <- get_tp53_status_tcga(
    projects = cancer_types,
    maf_dir = maf_dir
  )

  saveRDS(
    tp53_pat,
    tp53_rds
  )

  readr::write_csv(
    tp53_pat,
    tp53_csv
  )

  message(
    "Saved processed TP53 mutation table."
  )
}

tp53_by_cancer <- tp53_pat %>%
  dplyr::count(
    Cancer_Type,
    name = "n_TP53_mutated"
  ) %>%
  dplyr::arrange(Cancer_Type)

print(
  as.data.frame(tp53_by_cancer),
  row.names = FALSE
)

# ==============================================================================
# 4. JOIN TP53 STATUS TO SCORED TUMOURS
# ==============================================================================

tcga_tumor_tp53 <- tcga_pan_scored %>%
  dplyr::left_join(
    tp53_pat,
    by = c(
      "patient",
      "Cancer_Type"
    )
  ) %>%
  dplyr::mutate(
    p53_status = ifelse(
      is.na(p53_status),
      "WT",
      "Mutated"
    ),
    p53_status = factor(
      p53_status,
      levels = c(
        "WT",
        "Mutated"
      )
    )
  )

cat("\n========================================\n")
cat("TP53 STATUS IN SCORED TUMOURS\n")
cat("========================================\n")

print(
  table(
    tcga_tumor_tp53$p53_status
  )
)

print(
  prop.table(
    table(
      tcga_tumor_tp53$p53_status
    )
  )
)

# ==============================================================================
# 5. DEFINE WITHIN-CANCER LOW/HIGH POLYGLU-TP GROUPS
# ==============================================================================

tp53_hl <- tcga_tumor_tp53 %>%
  dplyr::filter(
    !is.na(rf_prob_pan),
    !is.na(p53_status),
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
  dplyr::filter(
    !is.na(group)
  ) %>%
  dplyr::mutate(
    group = factor(
      group,
      levels = c(
        "Low",
        "High"
      )
    ),
    Cancer_Type = factor(
      Cancer_Type
    )
  )

# ==============================================================================
# 6. DESCRIPTIVE STATISTICS
# ==============================================================================

tp53_summary <- tp53_hl %>%
  dplyr::count(
    group,
    p53_status,
    name = "n"
  ) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(
    total = sum(n),
    proportion = n / total
  ) %>%
  dplyr::ungroup()

print(
  as.data.frame(tp53_summary),
  row.names = FALSE
)

tp53_counts_by_cancer <- tp53_hl %>%
  dplyr::count(
    Cancer_Type,
    group,
    p53_status,
    name = "n"
  )

readr::write_csv(
  tp53_counts_by_cancer,
  file.path(
    out_dir,
    "TP53_counts_by_cancer_and_group.csv"
  )
)

# ==============================================================================
# 7. COCHRAN-MANTEL-HAENSZEL TEST STRATIFIED BY CANCER TYPE
# ==============================================================================

cmh_array <- xtabs(
  ~ group + p53_status + Cancer_Type,
  data = tp53_hl
)

valid_strata <- vapply(
  seq_len(dim(cmh_array)[3]),
  function(i) {
    x <- cmh_array[, , i]

    all(rowSums(x) > 0) &&
      all(colSums(x) > 0)
  },
  logical(1)
)

cmh_array_valid <- cmh_array[
  ,
  ,
  valid_strata,
  drop = FALSE
]

cat(
  "\nCancer strata included in CMH test: ",
  sum(valid_strata),
  " / ",
  length(valid_strata),
  "\n",
  sep = ""
)

cmh_test <- stats::mantelhaen.test(
  cmh_array_valid
)

print(cmh_test)

cmh_p <- cmh_test$p.value
cmh_or <- unname(cmh_test$estimate)
cmh_ci <- cmh_test$conf.int

# ==============================================================================
# 8. SUMMARY AND OUTPUT TABLE
# ==============================================================================

pooled_table <- table(
  tp53_hl$group,
  tp53_hl$p53_status
)

low_n <- sum(
  pooled_table["Low", ]
)

high_n <- sum(
  pooled_table["High", ]
)

low_mut <- pooled_table[
  "Low",
  "Mutated"
]

high_mut <- pooled_table[
  "High",
  "Mutated"
]

low_mut_pct <- 100 * low_mut / low_n
high_mut_pct <- 100 * high_mut / high_n

cat("\n========================================\n")
cat("TP53 ASSOCIATION ANALYSIS\n")
cat("Within-cancer PolyGlu-TP quartiles\n")
cat("========================================\n")

cat(
  "Low: n = ",
  low_n,
  " | TP53 mutated = ",
  low_mut,
  " (",
  round(low_mut_pct, 1),
  "%)\n",
  sep = ""
)

cat(
  "High: n = ",
  high_n,
  " | TP53 mutated = ",
  high_mut,
  " (",
  round(high_mut_pct, 1),
  "%)\n",
  sep = ""
)

cat(
  "CMH p = ",
  format(cmh_p, scientific = TRUE),
  "\n",
  sep = ""
)

cat(
  "Common odds ratio = ",
  round(cmh_or, 3),
  "\n",
  sep = ""
)

cat(
  "95% CI = ",
  round(cmh_ci[1], 3),
  " to ",
  round(cmh_ci[2], 3),
  "\n",
  sep = ""
)

tp53_stats <- tibble::tibble(
  analysis = "Within-cancer PolyGlu-TP quartiles",
  low_n = low_n,
  high_n = high_n,
  low_TP53_mutated = low_mut,
  high_TP53_mutated = high_mut,
  low_TP53_mutated_percent = low_mut_pct,
  high_TP53_mutated_percent = high_mut_pct,
  CMH_common_OR = cmh_or,
  CMH_CI_low = cmh_ci[1],
  CMH_CI_high = cmh_ci[2],
  CMH_p = cmh_p
)

readr::write_csv(
  tp53_stats,
  file.path(
    out_dir,
    "TP53_withinCancer_CMH_stats.csv"
  )
)

# ==============================================================================
# 9. PLOT
# ==============================================================================

p_label <- if (cmh_p < 2.2e-16) {
  "CMH p < 2.2 \u00d7 10\u207b\u00b9\u2076"
} else {
  paste0(
    "CMH p = ",
    format(
      cmh_p,
      scientific = TRUE,
      digits = 2
    )
  )
}

stat_label <- paste0(
  p_label,
  "\nOR = ",
  round(cmh_or, 2),
  " (95% CI ",
  round(cmh_ci[1], 2),
  "\u2013",
  round(cmh_ci[2], 2),
  ")"
)

p_tp53 <- ggplot2::ggplot(
  tp53_hl,
  ggplot2::aes(
    x = group,
    fill = p53_status
  )
) +
  ggplot2::geom_bar(
    position = "fill",
    color = "white",
    linewidth = 0.3,
    width = 0.72
  ) +
  ggplot2::scale_fill_manual(
    values = pal_p53,
    drop = FALSE
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    ),
    expand = ggplot2::expansion(
      mult = c(0.02, 0.05)
    )
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      "Low" = "Low\n(\u226425th percentile)",
      "High" = "High\n(\u226575th percentile)"
    )
  ) +
  ggplot2::annotate(
    "text",
    x = 1.5,
    y = 0.97,
    label = stat_label,
    hjust = 0.5,
    vjust = 1,
    size = 4
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Proportion",
    fill = "TP53"
  ) +
  ggplot2::theme_classic(
    base_size = 14
  ) +
  ggplot2::theme(
    axis.text = ggplot2::element_text(
      color = "black"
    ),
    axis.title = ggplot2::element_text(
      color = "black",
      face = "bold"
    ),
    legend.position = "right"
  )

print(p_tp53)
